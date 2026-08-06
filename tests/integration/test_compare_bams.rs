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
use std::fs;
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
    assert!(
        stdout.contains("RESULT: BAM groupings are EQUIVALENT"),
        "expected 'RESULT: BAM groupings are EQUIVALENT' in output, got:\n{stdout}"
    );
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
    assert!(
        stdout.contains("RESULT: BAM files are IDENTICAL"),
        "expected 'RESULT: BAM files are IDENTICAL' in output, got:\n{stdout}"
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

    let (code, stdout, _stderr) = run_compare(&bam1, &bam2, "grouping", &["--ignore-order"]);
    assert_eq!(code, Some(0), "Expected success for ignore-order grouping mode, stdout:\n{stdout}");
    assert!(
        stdout.contains("RESULT: BAM groupings are EQUIVALENT"),
        "expected 'RESULT: BAM groupings are EQUIVALENT' in output, got:\n{stdout}"
    );
}

/// Runs `CompareBams::execute()` in-process under `--mode grouping` (plus any `extra_args`,
/// e.g. `--ignore-order`) and returns the error's full display string. Panics if `execute`
/// succeeds, or if it fails with a [`CompareMismatch`] (a genuine record-level DIFFER) rather
/// than the MI-less hard-fail this helper exists to capture (`molecule_runs` rejects the first
/// record with no/unparseable MI) — mirrors `run_compare_content_expect_err`'s
/// success/DIFFER-vs-hard-error distinction for the grouping-mode MI-less rejection.
fn run_compare_grouping_expect_err(bam1: &Path, bam2: &Path, extra_args: &[&str]) -> String {
    let mut args: Vec<&str> =
        vec!["bams", bam1.to_str().unwrap(), bam2.to_str().unwrap(), "--mode", "grouping"];
    args.extend(extra_args);
    let cmd = CompareBams::try_parse_from(args).expect("failed to parse compare bams args");
    match cmd.execute("fgumi compare bams") {
        Ok(()) => panic!("expected compare bams to hard-fail on MI-less input"),
        Err(e) if e.is::<CompareMismatch>() => {
            panic!("expected the MI-less rejection error, got a CompareMismatch (DIFFER): {e:#}")
        }
        Err(e) => format!("{e:#}"),
    }
}

// ---------------------------------------------------------------------------
// Missing MI tag behavior (streaming molecule-join engine)
//
// The old (retired) key-join engine explicitly counted "missing MI" records and always
// DIFFERed when both files entirely lacked MI tags, on the theory that "no grouping has
// been verified" must never read as EQUIVALENT. The streaming molecule-join engine enforces
// the same principle, but *per record and upstream*: `molecule_runs` yields an `Err` on the
// first record whose MI tag is missing or unparseable (`base == None`), and
// `molecule_join_compare` propagates that error via `?`. So `--mode grouping`/`--command
// group` require already-grouped (same-MI-consecutive, every-read-MI-tagged) input, and
// feeding any non-empty input with an MI-less record on either side is a hard error rather
// than silently reporting MATCH or an ordinary molecule-count DIFFER.
//
// Enforcing this per record (rather than with a whole-input "never saw an MI" guard) is what
// closes the *partial*-MI-less soundness hole. Consider a side with some records tagged and
// some not, or an entirely-MI-less side whose single spanning run happens to
// canonical-id/content-match a real molecule on the grouped side: because the join compares
// content *excluding* MI, such a pair could report a false MATCH despite one side never
// having verified any grouping. A whole-input guard that only checked "did this side ever see
// an MI" would pass the mixed case (it *did* see one); the per-record check in `molecule_runs`
// rejects the first MI-less record regardless. See
// `test_grouping_mi_less_side_with_matching_records_on_other_side_now_errors` (the
// canonical-id-collision false MATCH), `test_grouping_partial_mi_less_pair_now_errors` (the
// count-asymmetry case that used to DIFFER), and
// `test_grouping_mixed_mi_within_one_side_errors` (the mixed-within-a-side case a whole-input
// guard would miss).
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
    // fully-MI-less rejection outcome.
    let err = run_compare_grouping_expect_err(&bam1, &bam2, &["--ignore-order"]);
    assert!(err.contains("MI-tagged"), "got: {err}");
}

/// Review fix (soundness hole closed): a partial-MI-less pair — one side grouped, the other
/// not — is now a hard error, not an ordinary molecule-count DIFFER. Before this fix the check
/// only fired when *both* sides were fully MI-less, so this case fell through to the normal
/// molecule-join path and `DIFFERed` via `bam1_molecules` (2) != `bam2_molecules` (1);
/// `molecule_runs`' per-record MI check now rejects the MI-less side up front, regardless of
/// whether the counts would have matched.
#[test]
fn test_grouping_partial_mi_less_pair_now_errors() {
    let tmp = TempDir::new().unwrap();
    let header = create_minimal_header("chr1", 10000);

    // bam1: two distinct MI-tagged molecules (2 runs).
    let grouped = vec![mi_record(b"read1", 100, "1"), mi_record(b"read2", 200, "2")];
    // bam2: same two records, but with NO MI tag at all -> molecule_runs folds them into
    // one run (consecutive base == None).
    let ungrouped = vec![mapped_record(b"read1", 100), mapped_record(b"read2", 200)];

    let bam1 = tmp.path().join("a.bam");
    let bam2 = tmp.path().join("b.bam");
    write_bam(&bam1, &header, &grouped);
    write_bam(&bam2, &header, &ungrouped);

    let err = run_compare_grouping_expect_err(&bam1, &bam2, &[]);
    assert!(err.contains("MI-tagged"), "got: {err}");
}

/// Independent oracle (mixed-within-a-side): the case a whole-input "did this side ever see an
/// MI" guard would MISS. bam2 is *partially* grouped — its first record carries an MI, its
/// second does not — so such a guard would pass it (the side did see an MI). Because the join
/// compares content excluding MI, `read2` losing its MI would otherwise still content-match its
/// tagged counterpart in bam1, a false MATCH for a tool that dropped MI on one molecule. The
/// per-record check in `molecule_runs` rejects bam2's first MI-less record (`read2`) outright.
/// This is a genuinely different failure mode from the fully-MI-less-side tests above, not a
/// self-consistency rephrasing of them.
#[test]
fn test_grouping_mixed_mi_within_one_side_errors() {
    let tmp = TempDir::new().unwrap();
    let header = create_minimal_header("chr1", 10000);

    // bam1: two fully MI-tagged molecules.
    let grouped = vec![mi_record(b"read1", 100, "1"), mi_record(b"read2", 200, "2")];
    // bam2: read1 keeps its MI, read2 has dropped it -> a *mixed* side. A whole-input guard
    // that only asks "did this side ever see an MI" would pass this (it did, on read1); the
    // per-record check rejects read2.
    let mixed = vec![mi_record(b"read1", 100, "1"), mapped_record(b"read2", 200)];

    let bam1 = tmp.path().join("a.bam");
    let bam2 = tmp.path().join("b.bam");
    write_bam(&bam1, &header, &grouped);
    write_bam(&bam2, &header, &mixed);

    let err = run_compare_grouping_expect_err(&bam1, &bam2, &[]);
    assert!(err.contains("MI-tagged"), "got: {err}");
    assert!(err.contains("not grouped"), "got: {err}");
}

/// CRITICAL soundness regression: the false-MATCH hole this fix closes. bam1 is entirely
/// MI-less; absent a per-record check `molecule_runs` would fold its two records into a single
/// spanning run. bam2 has ONE real MI-tagged molecule containing the *same* two records (same
/// canonical id — the lexicographically smallest read name — and the same member content).
/// Before this fix the check only fired when *both* sides were fully MI-less, so this pair fell
/// through to the ordinary molecule-join path: one run per side, matching canonical id,
/// matching members -> a false MATCH despite bam1 never having verified any grouping at all.
/// `molecule_runs` now rejects bam1's first MI-less record outright, before any match is made.
#[test]
fn test_grouping_mi_less_side_with_matching_records_on_other_side_now_errors() {
    let tmp = TempDir::new().unwrap();
    let header = create_minimal_header("chr1", 10000);

    // bam1: entirely MI-less -> one spanning run, canon = min("read1", "read2") = "read1".
    let bam1_records = vec![mapped_record(b"read1", 100), mapped_record(b"read2", 200)];
    // bam2: the SAME two records (content-identical aside from the added MI tag),
    // MI-tagged as one real molecule -> also one run with canon "read1" and the same
    // member records. Absent the fix, this pair would MATCH.
    let bam2_records = vec![mi_record(b"read1", 100, "1"), mi_record(b"read2", 200, "1")];

    let bam1 = tmp.path().join("a.bam");
    let bam2 = tmp.path().join("b.bam");
    write_bam(&bam1, &header, &bam1_records);
    write_bam(&bam2, &header, &bam2_records);

    let err = run_compare_grouping_expect_err(&bam1, &bam2, &[]);
    assert!(err.contains("MI-tagged"), "got: {err}");
}

/// Review fix (empty-vs-empty exemption): two EMPTY grouped BAMs (zero records on both
/// sides) must MATCH, not be rejected as MI-less. There are no records at all, so
/// `molecule_runs` yields no runs and never encounters an MI-less record to reject — and
/// there is also nothing to compare, so this is a vacuous MATCH, not the "misuse" case
/// (a non-empty input with an MI-less record) the per-record check exists to catch.
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

/// Review fix (error-vs-diff precedence): both existing MI-less tests above use *identical*
/// content on both sides, so the MI-less rejection was never distinguished from an ordinary
/// content DIFFER also firing. This test uses genuinely different content (same read
/// names/keys so the records would otherwise pair up and surface a POS content diff — see
/// `RecordKey`'s doc comment on POS not being part of a primary record's identity) on two
/// fully MI-less BAMs, and locks that the MI-less rejection still wins: the compare hard-fails
/// with the "not grouped / MI-tagged" message rather than reporting an ordinary DIFFER.
#[test]
fn test_grouping_mode_mi_missing_and_content_differs_still_errors() {
    let tmp = TempDir::new().unwrap();
    let header = create_minimal_header("chr1", 10000);

    // No MI tag on either side, but genuinely different content (different positions) so
    // that, absent the MI-less rejection, this would be an ordinary content DIFFER, not a MATCH.
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
    assert!(
        stdout.contains("RESULT: BAM groupings are EQUIVALENT"),
        "expected 'RESULT: BAM groupings are EQUIVALENT' in output, got:\n{stdout}"
    );
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
    assert!(
        stdout.contains("RESULT: BAM groupings are EQUIVALENT"),
        "expected 'RESULT: BAM groupings are EQUIVALENT' in output, got:\n{stdout}"
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

/// `positional_compare` is publicly re-exported, so a direct caller bypasses the CLI's
/// `--batch-size` parser guard. A zero batch size must error at the engine boundary (the guard
/// in `start_raw_batch_reader`) rather than hang: `read_raw_batch` would otherwise yield empty,
/// non-EOF batches forever and the consumer's `recv()` would never make progress.
#[test]
fn test_positional_compare_zero_batch_size_errors() {
    let tmp = TempDir::new().unwrap();
    let header = create_minimal_header("chr1", 10000);
    let records = vec![mapped_record(b"read1", 100)];

    let bam1 = tmp.path().join("a.bam");
    let bam2 = tmp.path().join("b.bam");
    write_bam(&bam1, &header, &records);
    write_bam(&bam2, &header, &records);

    let err = positional_compare(&bam1, &bam2, 1, 0, 100, ContentPredicate::Exact, None)
        .expect_err("batch_size 0 must error at the engine boundary, not hang");
    assert!(
        err.to_string().contains("batch size must be at least 1"),
        "expected the batch-size guard message, got: {err}"
    );
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

    let outcome = positional_compare(&bam1, &bam2, 1, 64, 100, ContentPredicate::Exact, None)
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

    let outcome = positional_compare(&bam1, &bam2, 1, 64, 100, ContentPredicate::Exact, None)
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

    let outcome = positional_compare(&bam1, &bam2, 1, 64, 100, ContentPredicate::Exact, None)
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

    let outcome = positional_compare(&bam1, &bam2, 1, 64, 100, ContentPredicate::Exact, None)
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

    let outcome = positional_compare(&bam1, &bam2, 1, 64, 100, ContentPredicate::Exact, None)
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

    // molecule A = {read1,read2}, molecule B = {read3}. bam1 emits A then B (MI 1, 2).
    let records1 = vec![
        mi_record(b"read1", 100, "1"),
        mi_record(b"read2", 200, "1"),
        mi_record(b"read3", 300, "2"),
    ];
    // bam2 emits B then A (reordered). MI is a monotonic emission-order counter, so B (emitted
    // first) gets 1 and A gets 2 — the same molecules, different physical order AND renumbered,
    // exactly as a different tool / tie-break would produce.
    let records2 = vec![
        mi_record(b"read3", 300, "1"),
        mi_record(b"read1", 100, "2"),
        mi_record(b"read2", 200, "2"),
    ];

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

/// `molecule_join_compare` must enforce its own header-compatibility precondition rather
/// than trusting the CLI caller to have run `require_compatible_headers` first: called
/// directly (as this test does, bypassing `CompareBams::execute`) on two BAMs with
/// incompatible `@SQ` dictionaries, it must return an `Err` naming the `@SQ` mismatch
/// rather than silently molecule-joining records against the wrong reference dictionary.
#[test]
fn molecule_join_compare_rejects_incompatible_sq_dictionaries() {
    let tmp = TempDir::new().unwrap();
    let bam1 = write_bam_with_ref(tmp.path(), "a.bam", "chr1", 1000);
    let bam2 = write_bam_with_ref(tmp.path(), "b.bam", "chr1", 2000);

    let err = molecule_join_compare(&bam1, &bam2, 10)
        .expect_err("incompatible @SQ dictionaries must be rejected by molecule_join_compare");
    assert!(format!("{err:#}").contains("@SQ"), "got: {err:#}");
}

/// CRITICAL soundness regression: `--max-diffs 0` must not collapse `is_match()` into a
/// bare molecule-count comparison. `push_diff` (see `engines/mod.rs`) is capped by
/// `max_diffs`, so with `max_diffs == 0` `diff_details` is *always* empty regardless of
/// how many real diffs were found — an `is_match()` keyed off `diff_details.is_empty()`
/// would report EQUIVALENT for this genuine content mismatch (equal molecule counts, one
/// matched pair's SEQ mutated) purely because the cap ate the evidence. `is_match()` must
/// instead key off the `matched` counter, which only increments on a molecule pair
/// `compare_molecule` found fully equivalent — independent of `max_diffs` entirely.
#[test]
fn molecule_join_compare_content_diff_survives_max_diffs_zero() {
    let tmp = TempDir::new().unwrap();
    let header = create_minimal_header("chr1", 10000);

    let records1 = vec![
        mapped_record_with_mi_and_seq(b"read1", 100, "1", b"ACGTACGT"),
        mapped_record_with_mi_and_seq(b"read2", 200, "1", b"ACGTACGT"),
    ];
    // Same grouping (both MI=1 on both sides, so bam1_molecules == bam2_molecules == 1),
    // but read2's SEQ is mutated on bam2 — a real content diff that `max_diffs == 0` would
    // hide entirely from `diff_details`.
    let records2 = vec![
        mapped_record_with_mi_and_seq(b"read1", 100, "1", b"ACGTACGT"),
        mapped_record_with_mi_and_seq(b"read2", 200, "1", b"ACGTACGA"),
    ];

    let bam1 = tmp.path().join("a.bam");
    let bam2 = tmp.path().join("b.bam");
    write_bam(&bam1, &header, &records1);
    write_bam(&bam2, &header, &records2);

    let outcome =
        molecule_join_compare(&bam1, &bam2, 0).expect("molecule_join_compare should succeed");

    assert!(
        outcome.diff_details.is_empty(),
        "max_diffs == 0 must suppress diff_details entirely (this is the trap): {outcome:?}"
    );
    assert_eq!(outcome.bam1_molecules, outcome.bam2_molecules, "molecule counts must be equal");
    assert_eq!(outcome.matched, 0, "the sole molecule pair must fail to match on the SEQ diff");
    assert!(
        !outcome.is_match(),
        "a real content diff must DIFFER even with max_diffs == 0 and equal molecule counts \
         (this is the CRITICAL soundness bug): {outcome:?}"
    );
}

/// Control for the above: a fully-equivalent pair must still MATCH under `max_diffs == 0`
/// — the cap must not turn a genuine match into a spurious DIFFER either.
#[test]
fn molecule_join_compare_equivalent_pair_matches_with_max_diffs_zero() {
    let tmp = TempDir::new().unwrap();
    let header = create_minimal_header("chr1", 10000);

    let records = vec![mi_record(b"read1", 100, "1"), mi_record(b"read2", 200, "1")];

    let bam1 = tmp.path().join("a.bam");
    let bam2 = tmp.path().join("b.bam");
    write_bam(&bam1, &header, &records);
    write_bam(&bam2, &header, &records);

    let outcome =
        molecule_join_compare(&bam1, &bam2, 0).expect("molecule_join_compare should succeed");

    assert_eq!(outcome.matched, outcome.bam1_molecules);
    assert_eq!(outcome.matched, outcome.bam2_molecules);
    assert!(
        outcome.is_match(),
        "a genuinely equivalent pair must MATCH even with max_diffs == 0: {outcome:?}"
    );
}

/// CRITICAL soundness regression: the "balanced-residual" case. bam1 has one molecule
/// (`only1`) that bam2 lacks, and bam2 has a *different* molecule (`only2`) that bam1
/// lacks, alongside one molecule (`shared`) present, identically, on both sides. Each
/// file therefore has 2 molecules total, so `bam1_molecules == bam2_molecules` even
/// though a real presence difference exists on both sides. Combined with `max_diffs ==
/// 0` (which forces `diff_details` empty regardless of how many diffs were found — see
/// `molecule_join_compare_content_diff_survives_max_diffs_zero`), the OLD verdict
/// (`diff_details.is_empty() && bam1_molecules == bam2_molecules`) would report a FALSE
/// MATCH here: both conditions hold despite two molecules genuinely differing in
/// presence. The NEW verdict (`matched == bam1_molecules && matched == bam2_molecules`)
/// correctly DIFFERs because only the one `shared` molecule pair actually matched
/// (`matched == 1 < bam1_molecules == 2`).
#[test]
fn molecule_join_compare_balanced_residual_differs_with_max_diffs_zero() {
    let tmp = TempDir::new().unwrap();
    let header = create_minimal_header("chr1", 10000);

    // bam1: one molecule shared with bam2, plus one molecule ("only1") bam2 lacks.
    let records1 = vec![
        mi_record(b"shared1", 100, "1"),
        mi_record(b"shared2", 100, "1"),
        mi_record(b"only1", 200, "2"),
    ];
    // bam2: the same shared molecule, plus a *different* molecule ("only2") bam1 lacks.
    // Total molecule count on each side is 2, so bam1_molecules == bam2_molecules.
    let records2 = vec![
        mi_record(b"shared1", 100, "1"),
        mi_record(b"shared2", 100, "1"),
        mi_record(b"only2", 300, "3"),
    ];

    let bam1 = tmp.path().join("a.bam");
    let bam2 = tmp.path().join("b.bam");
    write_bam(&bam1, &header, &records1);
    write_bam(&bam2, &header, &records2);

    let outcome =
        molecule_join_compare(&bam1, &bam2, 0).expect("molecule_join_compare should succeed");

    assert_eq!(outcome.bam1_molecules, 2);
    assert_eq!(outcome.bam2_molecules, 2);
    assert_eq!(
        outcome.bam1_molecules, outcome.bam2_molecules,
        "molecule counts must be equal (the balanced-residual setup): {outcome:?}"
    );
    assert!(
        outcome.diff_details.is_empty(),
        "max_diffs == 0 must suppress diff_details entirely (this is the trap): {outcome:?}"
    );
    assert_eq!(
        outcome.matched, 1,
        "only the shared molecule pair should match; the two residuals must not count: \
         {outcome:?}"
    );
    assert!(
        !outcome.is_match(),
        "a balanced residual (equal molecule counts, but one molecule only-in-bam1 and a \
         DIFFERENT molecule only-in-bam2) must DIFFER even with max_diffs == 0 — a verdict \
         keyed off diff_details.is_empty() && equal counts would falsely MATCH here: {outcome:?}"
    );
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

    // Same content, molecules reordered relative to each other, and MI values renumbered —
    // exactly what a cross-tool `group` comparison must tolerate. MI is a monotonic
    // emission-order counter, so the molecule emitted first here ({read3}) gets the lower MI (5)
    // and {read1,read2} the higher (9): reordered AND renumbered, each file still MI-monotonic.
    let records2 = vec![
        mi_record(b"read3", 300, "5"),
        mi_record(b"read1", 100, "9"),
        mi_record(b"read2", 200, "9"),
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
    assert!(
        stdout.contains("RESULT: BAM groupings are EQUIVALENT"),
        "expected 'RESULT: BAM groupings are EQUIVALENT' in output, got:\n{stdout}"
    );
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

/// CRITICAL soundness regression, end to end: the same content bug as the test above, but
/// run through the real CLI with `--max-diffs 0` (a plausible CI invocation — `--max-diffs`
/// is an unvalidated `usize`, so `0` is accepted). `push_diff` suppresses every diff line at
/// this cap, so a pre-fix `is_match()` keyed off `diff_details.is_empty()` would report
/// EQUIVALENT (exit 0) for a genuine content mismatch. Must still DIFFER (nonzero exit).
#[test]
fn test_command_group_content_diff_differs_with_max_diffs_zero() {
    let tmp = TempDir::new().unwrap();
    let header = create_minimal_header("chr1", 10000);

    let records1 = vec![
        mapped_record_with_mi_and_seq(b"read1", 100, "1", b"ACGTACGT"),
        mapped_record_with_mi_and_seq(b"read2", 200, "1", b"ACGTACGT"),
    ];
    let records2 = vec![
        mapped_record_with_mi_and_seq(b"read1", 100, "1", b"ACGTACGT"),
        mapped_record_with_mi_and_seq(b"read2", 200, "1", b"ACGTACGA"),
    ];

    let bam1 = tmp.path().join("a.bam");
    let bam2 = tmp.path().join("b.bam");
    write_bam(&bam1, &header, &records1);
    write_bam(&bam2, &header, &records2);

    let (code, stdout, stderr) = run_compare_command(&bam1, &bam2, "group", &["--max-diffs", "0"]);

    assert_eq!(
        code,
        Some(1),
        "`--max-diffs 0` must report a real content diff as DIFFER (exit 1, not the \
         EQUIVALENT exit 0 and not a signal/parse failure), stdout:\n{stdout}\nstderr:\n{stderr}"
    );
    assert!(stdout.contains("DIFFER"), "Expected DIFFER in output, got:\n{stdout}");
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

    assert_eq!(
        code,
        Some(1),
        "content-identical, fully MI-less BAMs must hard-fail via `--command group` \
         (`molecule_runs` rejects the MI-less input with an `Err`, so `main` exits 1 — \
         not a signal/parse failure), stdout:\n{stdout}\nstderr:\n{stderr}"
    );
    assert!(stderr.contains("MI-tagged"), "expected the MI-less rejection message, got:\n{stderr}");
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

/// The run *multiset* comparison itself: an equal-key run of size 2 that lines up on both
/// sides (same key, same length — so neither the run-boundary/desync path nor a run-length
/// difference can catch it), where one record's content differs. Only comparing the two runs'
/// record multisets detects this. The expectations are hand-computed rather than derived from
/// the engine: exactly one of the two runs differs, so `run_mismatches == 1`; the runs stay
/// aligned so there is no desync; and every record on both sides is read.
#[test]
fn test_sort_verify_content_diff_within_aligned_multi_record_run_differs() {
    let tmp = TempDir::new().unwrap();
    let header = create_coordinate_sorted_header("chr1", 10000);

    // The tied run at pos 100: readA is identical on both sides, readB's sequence differs by
    // its last base in bam2 while keeping the same name and coordinate (so it stays in the
    // same run and both files remain correctly coordinate-sorted).
    let read_a = fragment_record(b"readA", 0, 100, false);
    let read_b = fragment_record(b"readB", 0, 100, false);
    let read_b_mutated = {
        let mut b = SamBuilder::new();
        b.read_name(b"readB")
            .sequence(b"ACGTACGA") // last base differs from `fragment_record`'s ACGTACGT
            .qualities(&[30; 8])
            .ref_id(0)
            .pos(99) // same 1-based pos 100 as `read_b`
            .mapq(60);
        b.build()
    };
    // A second, fully identical run at pos 200 — it must not be counted as a mismatch.
    let read_c = fragment_record(b"readC", 0, 200, false);

    let bam1 = tmp.path().join("a.bam");
    let bam2 = tmp.path().join("b.bam");
    write_bam(&bam1, &header, &[read_a.clone(), read_b, read_c.clone()]);
    write_bam(&bam2, &header, &[read_a, read_b_mutated, read_c]);

    let outcome =
        sort_verify_compare(&bam1, &bam2, 100).expect("sort_verify_compare should succeed");

    assert_eq!(outcome.bam1_violations, 0, "bam1 is itself correctly sorted: {outcome:?}");
    assert_eq!(outcome.bam2_violations, 0, "bam2 is itself correctly sorted: {outcome:?}");
    assert!(!outcome.presence_mismatch, "the runs line up on both sides: {outcome:?}");
    assert_eq!(outcome.bam1_count, 3, "{outcome:?}");
    assert_eq!(outcome.bam2_count, 3, "{outcome:?}");
    assert_eq!(
        outcome.run_mismatches, 1,
        "exactly the pos-100 run's multiset differs; the pos-200 run matches: {outcome:?}"
    );
    assert!(
        !outcome.is_match(),
        "a content difference inside an aligned tied run must DIFFER: {outcome:?}"
    );
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

/// Appends an additional `@SQ` entry to `header`, yielding a reference dictionary one
/// longer than the input's.
fn header_with_extra_reference(mut header: Header, ref_name: &str, ref_len: usize) -> Header {
    use bstr::BString;
    use noodles::sam::header::record::value::{Map, map::ReferenceSequence};
    use std::num::NonZeroUsize;

    header.reference_sequences_mut().insert(
        BString::from(ref_name),
        Map::<ReferenceSequence>::new(
            NonZeroUsize::new(ref_len).expect("reference length must be non-zero"),
        ),
    );
    header
}

/// Companion to the `@SQ` *length* mismatch above, covering the other way a reference
/// dictionary can diverge: a differing `@SQ` *count* (`nref`). bam2 declares a second
/// reference that bam1 does not, while every record still lives on `chr1` and both files
/// remain correctly coordinate-sorted — so neither the per-file order check nor the run
/// comparison can see it, and only the header comparison catches the divergence.
#[test]
fn test_sort_verify_header_sq_count_mismatch_differs() {
    let tmp = TempDir::new().unwrap();
    let header1 = create_coordinate_sorted_header("chr1", 10000);
    let header2 =
        header_with_extra_reference(create_coordinate_sorted_header("chr1", 10000), "chr2", 20000);

    // Pin the divergence the test is about: same first @SQ, different reference counts.
    assert_eq!(header1.reference_sequences().len(), 1, "bam1 declares one @SQ");
    assert_eq!(header2.reference_sequences().len(), 2, "bam2 declares two @SQ");

    let records = vec![fragment_record(b"read1", 0, 100, false)];

    let bam1 = tmp.path().join("a.bam");
    let bam2 = tmp.path().join("b.bam");
    write_bam(&bam1, &header1, &records);
    write_bam(&bam2, &header2, &records);

    let outcome =
        sort_verify_compare(&bam1, &bam2, 100).expect("sort_verify_compare should succeed");

    assert_eq!(outcome.bam1_violations, 0, "bam1 is itself correctly sorted: {outcome:?}");
    assert_eq!(outcome.bam2_violations, 0, "bam2 is itself correctly sorted: {outcome:?}");
    assert_eq!(outcome.bam1_count, 1, "{outcome:?}");
    assert_eq!(outcome.bam2_count, 1, "{outcome:?}");
    assert_eq!(outcome.run_mismatches, 0, "record content itself is identical");
    assert!(outcome.header_mismatch, "a @SQ count mismatch must be flagged: {outcome:?}");
    assert!(
        !outcome.is_match(),
        "an nref divergence must DIFFER even with matching sort order/content: {outcome:?}"
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
    assert!(
        stdout.contains("RESULT: BAM files are IDENTICAL"),
        "expected 'RESULT: BAM files are IDENTICAL' in output, got:\n{stdout}"
    );
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
/// rejected (don't trust the @HD tag), independent of record content.
///
/// Asserts the specific diagnostic rather than merely that the word "order" appears: the
/// message must name the declared order that was violated and the offending input, so that a
/// change to the order check cannot degrade the report while still "containing order".
#[test]
fn content_mode_rejects_records_not_in_declared_coordinate_order() {
    let tmp = TempDir::new().unwrap();
    let bad = write_bam_declaring_coordinate_but_unsorted(tmp.path()); // pos 500 then pos 100
    let good = write_sorted_coordinate_bam(tmp.path());
    let err = run_compare_content_expect_err(&good, &bad);
    assert!(err.contains("violate the declared"), "got: {err}");
    assert!(err.contains("Coordinate"), "must name the declared order violated: {err}");
    let bad_name = bad.file_name().unwrap().to_str().unwrap();
    assert!(err.contains(bad_name), "must name the offending input {bad_name}: {err}");
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

    // bam2: M_b before M_a (reversed file order), MI renumbered. Both tools assign the MI as a
    // monotonically increasing counter in emission order, so the molecule emitted first (M_b)
    // gets the lower MI (7) and M_a the higher (9) — different values from bam1, still monotonic.
    let mut bam2_records = vec![mi_record(b"mol_b1", 200, "7"), mi_record(b"mol_b2", 200, "7")];
    bam2_records.extend([mi_record(b"mol_a1", 100, "9"), mi_record(b"mol_a2", 100, "9")]);

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

/// Deterministic residual diagnostics: the EOF residual drain used to iterate
/// `pending1`/`pending2` (both `AHashMap`) via `.drain()`, so which residual canonical ids
/// got reported (once the count exceeded `max_diffs`) and their relative order varied
/// across runs with the hash map's iteration order. bam1 has five residual-only molecules
/// (no counterpart in bam2) with canonical ids that sort very differently from their
/// insertion order; with `max_diffs` capping the reported residuals to 3, the residual
/// diagnostics must consistently report the *lexicographically first* three ids, in sorted
/// order, run after run.
/// Both residual sides are covered: the residual molecules live entirely on one file, the
/// other is empty, and the sorted, `max_diffs`-capped `only in bamN` diagnostics must be
/// deterministic run-to-run. `pending2` (the `only in bam2` side) is exercised as well as
/// `pending1`, so a regression to hash-order or wrong capping on either side is caught.
#[rstest]
#[case::residuals_only_in_bam1(true, "bam1")]
#[case::residuals_only_in_bam2(false, "bam2")]
fn molecule_join_residual_diagnostics_are_deterministic_and_sorted(
    #[case] residuals_in_bam1: bool,
    #[case] side_label: &str,
) {
    let tmp = TempDir::new().unwrap();
    let header = create_minimal_header("chr1", 10000);

    // Insertion order deliberately not sorted, so a hash-order-dependent drain would
    // report a different subset/order than the lexicographically-sorted one.
    let canon_ids = ["zebra", "mango", "apple", "banana", "cherry"];
    let records: Vec<RawRecord> = canon_ids
        .iter()
        .enumerate()
        .map(|(i, name)| mi_record(name.as_bytes(), 100, &i.to_string()))
        .collect();

    let bam1 = tmp.path().join("residuals1.bam");
    let bam2 = tmp.path().join("residuals2.bam");
    let (records1, records2): (&[RawRecord], &[RawRecord]) =
        if residuals_in_bam1 { (&records, &[]) } else { (&[], &records) };
    write_bam(&bam1, &header, records1);
    write_bam(&bam2, &header, records2);

    // Lexicographically first three of {zebra, mango, apple, banana, cherry}.
    let expected: Vec<String> = ["apple", "banana", "cherry"]
        .iter()
        .map(|id| format!("molecule {id} only in {side_label}"))
        .collect();

    let max_diffs = 3;
    for _ in 0..5 {
        let out = molecule_join_compare(&bam1, &bam2, max_diffs)
            .expect("molecule_join_compare should succeed");
        assert!(!out.is_match());
        assert_eq!(
            out.diff_details, expected,
            "residual diagnostics must be sorted and deterministic across runs: {:?}",
            out.diff_details
        );
    }
}

/// Regression guard for the multiset-membership fix: `RecordKey` is collision-resistant,
/// not collision-free (`record_key.rs:20-24`) — for a primary record it reduces to
/// `(name, segment)`, with no locus discriminator, so two records that share a name and
/// segment but differ in position/content share one `RecordKey`. `index_by_key` used to
/// build a `RecordKey -> &RawRecord` `BTreeMap` (last-wins on a collision), so two records
/// sharing a key inside one molecule silently collapsed to one entry — a genuine 3-vs-2
/// physical-record multiplicity difference reported MATCH. `index_by_key` now maps each
/// `RecordKey` to *all* its member records (a multiset), and `compare_molecule` compares
/// per-key multiplicity, so this same fixture must now DIFFER.
///
/// bam1's molecule holds THREE physical records: `r1`, plus two same-name/same-segment `dup`
/// records at different positions (200 and 300) that collide to one `RecordKey`. bam2's
/// molecule holds only TWO: `r1` and a `dup` matching the position (300) of whichever bam1
/// record used to win the `BTreeMap` collapse.
#[test]
fn molecule_join_duplicate_record_key_within_a_molecule_now_differs() {
    let tmp = TempDir::new().unwrap();
    let header = create_minimal_header("chr1", 10000);

    // bam1: r1, plus two colliding same-name/same-segment "dup" records (pos 200 then 300).
    let bam1_records =
        vec![mi_record(b"r1", 100, "5"), mi_record(b"dup", 200, "5"), mi_record(b"dup", 300, "5")];
    // bam2: r1, plus only one "dup" record (pos 300) -> a 2-vs-1 multiplicity difference
    // for the "dup" RecordKey.
    let bam2_records = vec![mi_record(b"r1", 100, "9"), mi_record(b"dup", 300, "9")];

    let bam1 = tmp.path().join("dup_key1.bam");
    let bam2 = tmp.path().join("dup_key2.bam");
    write_bam(&bam1, &header, &bam1_records);
    write_bam(&bam2, &header, &bam2_records);

    let out =
        molecule_join_compare(&bam1, &bam2, 10).expect("molecule_join_compare should succeed");
    assert!(
        !out.is_match(),
        "multiset membership must catch the 2-vs-1 \"dup\" RecordKey multiplicity \
         difference, not silently collapse it to a MATCH: {out:?}"
    );
    assert!(
        out.diff_details.iter().any(|d| d.contains("dup") && d.contains("multiplicity")),
        "expected a multiplicity diff naming the colliding \"dup\" key: {:?}",
        out.diff_details
    );
}

// ============================================================================
// Task 10: the streaming molecule-join must never spill to disk.
// ============================================================================
//
// Unlike the retired keyjoin engine (deleted in Task 7), `molecule_join_compare` is a
// streaming two-sided hash-join: it holds only the bounded set of currently in-flight
// molecules in memory and never writes a sorted run, scratch file, or spill directory. The
// test below proves that by pointing the process's temp-dir env vars at a *non-directory* for
// a large grouped-pair compare: any attempt to materialize a scratch file there fails with
// `ENOTDIR`, so a spilling implementation could not succeed — and the compare succeeds.

/// Writes two large grouped BAMs encoding `n_molecules` molecules apiece, each molecule
/// comprising two member reads (`m{i}_1`, `m{i}_2`) that share an MI tag.
///
/// bam2's molecules are emitted in adjacent-pair-swapped physical order relative to bam1
/// (`0,1,2,3,...` becomes `1,0,3,2,...`) — a bounded-distance reordering (never more than one
/// molecule of lookahead) that still forces the join to buffer more than a single molecule at a
/// time, without requiring unbounded buffering the way a fully-scrambled order would. The MI in
/// each file is assigned by *file position* (a monotonic counter, as both grouping tools do), so
/// a physical molecule carries a different MI in bam2 than in bam1 — exercising the join's
/// renumber-tolerance and order-independence while keeping each file's MI monotonic.
///
/// Used by [`molecule_join_does_not_spill_to_disk`] to prove the streaming molecule-join
/// engine never spills: peak buffering scales with the (bounded) number of in-flight
/// molecules, not with file size.
fn write_large_grouped_pair(dir: &Path, n_molecules: usize) -> (PathBuf, PathBuf) {
    let header = create_minimal_header("chr1", 10_000);

    // A molecule's *identity* is its read names (`m{physical}_*`); its MI is assigned separately
    // by emission-order file position so each file stays MI-monotonic.
    let molecule_records = |physical: usize, mi: usize| -> Vec<RawRecord> {
        let mi_s = mi.to_string();
        vec![
            mi_record(format!("m{physical}_1").as_bytes(), 100, &mi_s),
            mi_record(format!("m{physical}_2").as_bytes(), 100, &mi_s),
        ]
    };

    // bam1: physical molecule i at position i, MI = i (monotonic).
    let bam1_records: Vec<RawRecord> =
        (0..n_molecules).flat_map(|i| molecule_records(i, i)).collect();

    // bam2: adjacent-pair-swapped physical order, MI = file position (still monotonic).
    let mut swapped_order: Vec<usize> = (0..n_molecules).collect();
    for pair in swapped_order.chunks_mut(2) {
        pair.reverse();
    }
    let bam2_records: Vec<RawRecord> = swapped_order
        .into_iter()
        .enumerate()
        .flat_map(|(pos, physical)| molecule_records(physical, pos))
        .collect();

    let bam1 = dir.join("large_grouped1.bam");
    let bam2 = dir.join("large_grouped2.bam");
    write_bam(&bam1, &header, &bam1_records);
    write_bam(&bam2, &header, &bam2_records);
    (bam1, bam2)
}

/// Grouping comparison never spills: a large input with molecules interleaved at bounded
/// distance compares with bounded buffering, and writes nothing at all to the process's temp
/// location.
///
/// Inspecting a temp *directory* after the run cannot catch a spill that a spilling engine
/// creates and then deletes mid-compare — the directory is empty again by the time the check
/// runs. This test closes that gap by pointing the standard temp-dir env vars
/// (`TMPDIR`/`TEMP`/`TMP`) at a regular *file* rather than a directory: any attempt to
/// materialize a scratch file *under* that path fails with `ENOTDIR`, so a spilling
/// implementation would abort. The compare succeeding therefore proves nothing was spilled —
/// transient or persisted.
///
/// The BAM fixtures live in their own directory, separate from the non-directory temp target,
/// so the fixtures themselves are unaffected.
///
/// This drives the real `fgumi compare bams --mode grouping` subprocess, rather than calling
/// `molecule_join_compare` in-process, so the env vars are scoped to the child process; the
/// test crate denies `unsafe_code`, and `std::env::set_var` requires `unsafe` as of the 2024
/// edition (mutating process-global env state is not thread-safe), so scoping via
/// `Command::env` on a subprocess is the sound way to do this here.
#[test]
fn molecule_join_does_not_spill_to_disk() {
    let fixtures = TempDir::new().unwrap();
    let (bam1, bam2) = write_large_grouped_pair(fixtures.path(), 50_000 /* molecules */);

    // A regular file, not a directory: creating any scratch file *under* this path fails with
    // ENOTDIR, so a spilling implementation cannot succeed.
    let tmp_holder = TempDir::new().unwrap();
    let non_dir_tmp = tmp_holder.path().join("not-a-directory");
    std::fs::write(&non_dir_tmp, b"").expect("failed to create non-directory temp target");

    let output = Command::new(env!("CARGO_BIN_EXE_fgumi"))
        .args(["compare", "bams"])
        .arg(&bam1)
        .arg(&bam2)
        .args(["--mode", "grouping"])
        .env("TMPDIR", &non_dir_tmp)
        .env("TEMP", &non_dir_tmp)
        .env("TMP", &non_dir_tmp)
        .output()
        .expect("failed to run fgumi compare bams");

    assert!(
        output.status.success(),
        "the large grouped pair must compare EQUIVALENT with spilling made impossible (temp \
         dir pointed at a non-directory), got exit {:?}\nstdout: {}\nstderr: {}",
        output.status.code(),
        String::from_utf8_lossy(&output.stdout),
        String::from_utf8_lossy(&output.stderr),
    );
}

// ============================================================================
// Single-pass sort-verify: complete per-file order verification past a desync
// ============================================================================
//
// The single-streaming-pass engine reads each file once, folding the per-file monotonic
// sort-order check into the same read the run-grouped comparison consumes (via an
// order-checking stream adapter). The invariant most at risk under that fold: once the two
// files' runs desynchronize (a run-boundary mismatch or one file exhausting first), the
// comparison stops but each file's *tail* must still be pulled through its order tracker, so
// the per-file violation count and first violation stay complete — identical to the old
// end-to-end `verify_sort_order` pass. These tests pin that behavior for a violation that
// lands in the tail, *past* the desync point (discovered by the drain path, not by the
// comparison).

/// Independently recompute the coordinate-order verification of `path` via the reference
/// `fgumi_sort::verify_sort_order` two-pass primitive (the exact fold the single-pass engine
/// folds inline), returning `(violations, first_violation)`. Used to lock the equivalence.
fn reference_coordinate_violations(path: &Path, nref: u32) -> (u64, Option<(u64, String)>) {
    let mut reader = fgumi_sort::RawBamRecordReader::new(fs::File::open(path).expect("open bam"))
        .expect("reader");
    reader.skip_header().expect("skip header");
    let (_total, violations, first) = fgumi_sort::verify_sort_order(
        reader,
        |bam: &[u8]| fgumi_sort::extract_coordinate_key_inline(bam, nref),
        |key: &u64, prev: &u64| key < prev,
    )
    .expect("verify_sort_order");
    (violations, first)
}

/// A desynced coordinate-sorted pair where *each* file additionally hides a monotonic-order
/// violation deep in its tail — after the run boundaries stop lining up. The comparison
/// desyncs at the second run (bam1's pos-200 group vs bam2's pos-300 group don't align) and
/// never resyncs, so the violating tail records (`readD` @ pos 50 in bam1, `readZ` @ pos 10
/// in bam2) are only ever seen by the drain path. Both must still be counted, with the
/// correct 1-based record number and read name, matching the reference two-pass computation.
#[test]
fn test_sort_verify_counts_order_violation_in_tail_past_desync() {
    let tmp = TempDir::new().unwrap();
    let header = create_coordinate_sorted_header("chr1", 10000);

    // bam1: readA/readB/readC are non-decreasing; readD @ pos 50 violates (50 < 210) and
    // sits past the run boundary the comparison desyncs on.
    let records1 = vec![
        fragment_record(b"readA", 0, 100, false),
        fragment_record(b"readB", 0, 200, false),
        fragment_record(b"readC", 0, 210, false),
        fragment_record(b"readD", 0, 50, false),
    ];
    // bam2: readX/readY are non-decreasing after the shared readA; readZ @ pos 10 violates
    // (10 < 310), likewise in the drained tail.
    let records2 = vec![
        fragment_record(b"readA", 0, 100, false),
        fragment_record(b"readX", 0, 300, false),
        fragment_record(b"readY", 0, 310, false),
        fragment_record(b"readZ", 0, 10, false),
    ];

    let bam1 = tmp.path().join("a.bam");
    let bam2 = tmp.path().join("b.bam");
    write_bam(&bam1, &header, &records1);
    write_bam(&bam2, &header, &records2);

    let outcome =
        sort_verify_compare(&bam1, &bam2, 100).expect("sort_verify_compare should succeed");

    // The streams desynchronized (run boundaries do not align), and no resync happened.
    assert!(outcome.presence_mismatch, "the pair must desync on the mismatched run: {outcome:?}");

    // Every record in each file was still read: counts are complete.
    assert_eq!(outcome.bam1_count, 4, "bam1 fully drained past the desync: {outcome:?}");
    assert_eq!(outcome.bam2_count, 4, "bam2 fully drained past the desync: {outcome:?}");

    // The tail-past-desync order violation in each file is fully counted, with the correct
    // 1-based record number and read name.
    assert_eq!(outcome.bam1_violations, 1, "bam1's tail violation must be counted: {outcome:?}");
    assert_eq!(outcome.bam1_first_violation, Some((4, "readD".to_string())));
    assert_eq!(outcome.bam2_violations, 1, "bam2's tail violation must be counted: {outcome:?}");
    assert_eq!(outcome.bam2_first_violation, Some((4, "readZ".to_string())));

    // And it matches the reference two-pass computation exactly.
    let nref = u32::try_from(header.reference_sequences().len()).expect("nref fits in u32");
    assert_eq!(reference_coordinate_violations(&bam1, nref), (1, Some((4, "readD".to_string()))));
    assert_eq!(reference_coordinate_violations(&bam2, nref), (1, Some((4, "readZ".to_string()))));

    assert!(!outcome.is_match(), "a desynced pair with mis-sorted tails must DIFFER: {outcome:?}");
}

proptest::proptest! {
    #![proptest_config(proptest::prelude::ProptestConfig::with_cases(64))]

    /// Parity/equivalence property: over randomized coordinate-keyed inputs — including
    /// pairs that desynchronize — the single-pass engine's per-file `(violations,
    /// first_violation)` must equal the reference `fgumi_sort::verify_sort_order` two-pass
    /// computation run independently over each file. This locks the single-pass fold to the
    /// primitive it replaced, for both correctly-sorted and mis-sorted inputs.
    #[test]
    fn prop_sort_verify_per_file_violations_match_reference(
        positions1 in proptest::collection::vec(1i32..500, 1..12),
        positions2 in proptest::collection::vec(1i32..500, 1..12),
    ) {
        let tmp = TempDir::new().unwrap();
        let header = create_coordinate_sorted_header("chr1", 10000);
        let nref = u32::try_from(header.reference_sequences().len()).expect("nref fits in u32");

        let to_records = |positions: &[i32]| -> Vec<RawRecord> {
            positions
                .iter()
                .enumerate()
                .map(|(i, &pos)| fragment_record(format!("r{i}").as_bytes(), 0, pos, false))
                .collect()
        };
        let records1 = to_records(&positions1);
        let records2 = to_records(&positions2);

        let bam1 = tmp.path().join("a.bam");
        let bam2 = tmp.path().join("b.bam");
        write_bam(&bam1, &header, &records1);
        write_bam(&bam2, &header, &records2);

        let outcome =
            sort_verify_compare(&bam1, &bam2, 100).expect("sort_verify_compare should succeed");

        let (ref1_violations, ref1_first) = reference_coordinate_violations(&bam1, nref);
        let (ref2_violations, ref2_first) = reference_coordinate_violations(&bam2, nref);

        proptest::prop_assert_eq!(outcome.bam1_violations, ref1_violations);
        proptest::prop_assert_eq!(outcome.bam1_first_violation, ref1_first);
        proptest::prop_assert_eq!(outcome.bam2_violations, ref2_violations);
        proptest::prop_assert_eq!(outcome.bam2_first_violation, ref2_first);

        // Record counts must also be complete regardless of any desync.
        proptest::prop_assert_eq!(outcome.bam1_count, records1.len() as u64);
        proptest::prop_assert_eq!(outcome.bam2_count, records2.len() as u64);
    }
}

// ============================================================================
// Both engines read each input exactly once
// ============================================================================

/// Serves `source`'s bytes over a fresh FIFO in `dir`, returning it and its feeder.
///
/// A FIFO delivers its bytes exactly once and cannot be re-opened for a second
/// read, so it fails any consumer that opens its input twice — which is what these
/// tests are for.
///
/// The write end is opened write-only, which blocks until the consumer attaches its
/// read end. That blocking open is the hand-off: opening `O_RDWR` instead would
/// return immediately, and this thread could then copy the whole (small) file into
/// the pipe buffer and close both ends before the consumer ever opened it, leaving
/// the bytes discarded and the consumer waiting on a writer that no longer exists.
/// A failed copy returns quietly rather than panicking so the consumer's own error
/// is the one that reports.
///
/// The returned handle must only be joined once the consumer is known to have
/// finished — see [`bounded`] — because a consumer that dies before opening the
/// FIFO leaves this thread blocked in `open` forever.
#[cfg(unix)]
fn serve_over_fifo(
    dir: &Path,
    name: &str,
    source: &Path,
) -> (PathBuf, std::thread::JoinHandle<()>) {
    let fifo = dir.join(name);
    let status = Command::new("mkfifo").arg(&fifo).status().expect("failed to run mkfifo");
    assert!(status.success(), "mkfifo failed for {}", fifo.display());

    let source = source.to_path_buf();
    let fifo_path = fifo.clone();
    let feeder = std::thread::spawn(move || {
        // `source` is a regular file this test just wrote, so opening it cannot
        // block: tolerating a failure here would buy no hang-safety and would
        // instead strand the consumer on an empty FIFO, making the timeout below
        // report a second open that never happened. The FIFO open and the copy are
        // the calls that can legitimately block or fail (`EPIPE`) when the consumer
        // dies, so those stay quiet.
        let mut src = fs::File::open(&source).expect("open FIFO source");
        let Ok(mut sink) = fs::File::create(&fifo_path) else { return };
        let _ = std::io::copy(&mut src, &mut sink);
    });

    (fifo, feeder)
}

/// Runs `engine` on its own thread, failing the test if it has not returned within
/// a generous bound.
///
/// A second open of a FIFO **blocks** rather than erroring — there is no writer
/// left to attach — so an engine that regresses to two opens would hang the suite
/// instead of failing it, and CI would report a timeout with no indication of
/// which property broke. This turns that into a named assertion. The stuck thread
/// is left behind deliberately: it cannot be cancelled, and the harness reaps it
/// when the process exits.
#[cfg(unix)]
fn bounded<T, Engine>(engine: Engine) -> anyhow::Result<T>
where
    Engine: FnOnce() -> anyhow::Result<T> + Send + 'static,
    T: Send + 'static,
{
    let (tx, rx) = std::sync::mpsc::channel();
    std::thread::spawn(move || {
        let _ = tx.send(engine());
    });
    rx.recv_timeout(std::time::Duration::from_secs(30)).unwrap_or_else(|_| {
        panic!(
            "engine did not return within 30s on a FIFO input; the likely cause is a \
             second open of the input, which blocks forever on a FIFO because no writer \
             is left to attach"
        )
    })
}

/// Joins a [`serve_over_fifo`] feeder, which is safe only once its consumer has
/// returned: the consumer having drained the FIFO is what lets the feeder's copy
/// finish.
#[cfg(unix)]
fn join_feeder(feeder: std::thread::JoinHandle<()>) {
    feeder.join().expect("FIFO feeder thread panicked");
}

/// `sort_verify_compare` must take a single open per input.
///
/// It needs both the header (to pick a key extractor) and the records, and used to
/// take them from two separate opens of the same path. That quietly required the
/// input to be re-openable, so the engine rejected a FIFO or a process
/// substitution that `fgumi sort` itself reads happily. Both now come from one
/// stream, the header parsed through a tee with its bytes replayed.
#[cfg(unix)]
#[test]
fn test_sort_verify_reads_a_non_reopenable_input() {
    let tmp = TempDir::new().unwrap();
    let header = create_coordinate_sorted_header("chr1", 10000);
    let records =
        vec![fragment_record(b"read1", 0, 100, false), fragment_record(b"read2", 0, 200, false)];

    let bam = tmp.path().join("a.bam");
    write_bam(&bam, &header, &records);

    // Both arguments are FIFOs: serving only the first would leave a regression that
    // reopens just the second indistinguishable from a passing run.
    let (fifo1, feeder1) = serve_over_fifo(tmp.path(), "a.fifo", &bam);
    let (fifo2, feeder2) = serve_over_fifo(tmp.path(), "b.fifo", &bam);
    let outcome = bounded(move || sort_verify_compare(&fifo1, &fifo2, 100))
        .expect("sort_verify_compare must accept inputs it can only open once");
    join_feeder(feeder1);
    join_feeder(feeder2);

    assert!(outcome.is_match(), "the same data over a FIFO must match itself: {outcome:?}");
    assert_eq!(outcome.bam1_count, records.len() as u64, "every record must be read from the FIFO");
    assert_eq!(outcome.bam2_count, records.len() as u64);
    assert_eq!(outcome.bam1_violations, 0);
}

/// `molecule_join_compare` must take a single open per input, for the same reason
/// as [`test_sort_verify_reads_a_non_reopenable_input`]: it reads the header to
/// check the two inputs are compatible, and the records to join molecules.
#[cfg(unix)]
#[test]
fn test_molecule_join_reads_a_non_reopenable_input() {
    let tmp = TempDir::new().unwrap();
    let header = create_minimal_header("chr1", 10000);
    let records = vec![mi_record(b"read1", 100, "1"), mi_record(b"read2", 200, "1")];

    let bam = tmp.path().join("a.bam");
    write_bam(&bam, &header, &records);

    // Both arguments are FIFOs, for the reason given in
    // [`test_sort_verify_reads_a_non_reopenable_input`].
    let (fifo1, feeder1) = serve_over_fifo(tmp.path(), "a.fifo", &bam);
    let (fifo2, feeder2) = serve_over_fifo(tmp.path(), "b.fifo", &bam);
    let outcome = bounded(move || molecule_join_compare(&fifo1, &fifo2, 100))
        .expect("molecule_join_compare must accept inputs it can only open once");
    join_feeder(feeder1);
    join_feeder(feeder2);

    assert!(outcome.is_match(), "the same data over a FIFO must match itself: {outcome:?}");
    assert_eq!(outcome.matched, 1, "the single molecule must be matched");
    assert_eq!(outcome.bam1_molecules, 1);
    assert_eq!(outcome.bam2_molecules, 1);
}

/// Waits for `child`, failing the test if it has not exited within a generous bound.
///
/// A second open of a FIFO **blocks** rather than erroring — there is no writer left
/// to attach — so a command that regresses to two opens would hang this test
/// indefinitely instead of failing it, and CI would report a timeout with no
/// indication of which property broke. The child is killed so the run does not leak
/// a blocked process.
///
/// Both pipes are drained on their own threads for the whole wait, rather than only
/// after the child exits. `stdout` and `stderr` are both `piped()`, so a child whose
/// combined output outgrows the OS pipe buffer (16 KiB on macOS, 64 KiB on Linux)
/// would block in `write` before reaching `exit`, `try_wait` would never return
/// `Some`, and the timeout would fire — blaming a second open that never happened.
/// Today's output is far short of that, but the diagnostic must not be able to lie.
#[cfg(unix)]
fn wait_bounded(mut child: std::process::Child, label: &str) -> std::process::Output {
    use std::io::Read;

    let mut child_stdout = child.stdout.take().expect("stdout was piped");
    let mut child_stderr = child.stderr.take().expect("stderr was piped");
    let drain_stdout = std::thread::spawn(move || {
        let mut buf = Vec::new();
        child_stdout.read_to_end(&mut buf).expect("drain fgumi stdout");
        buf
    });
    let drain_stderr = std::thread::spawn(move || {
        let mut buf = Vec::new();
        child_stderr.read_to_end(&mut buf).expect("drain fgumi stderr");
        buf
    });

    let deadline = std::time::Instant::now() + std::time::Duration::from_secs(30);
    loop {
        match child.try_wait().expect("poll fgumi") {
            Some(status) => {
                return std::process::Output {
                    status,
                    stdout: drain_stdout.join().expect("stdout drain thread"),
                    stderr: drain_stderr.join().expect("stderr drain thread"),
                };
            }
            None if std::time::Instant::now() >= deadline => {
                let _ = child.kill();
                let _ = child.wait();
                panic!(
                    "compare bams --command {label} did not exit within 30s on a FIFO input; \
                     the likely cause is a second open of the input, which blocks forever on a \
                     FIFO because no writer is left to attach"
                );
            }
            None => std::thread::sleep(std::time::Duration::from_millis(20)),
        }
    }
}

/// `compare bams` must accept inputs it can only open once, end to end.
///
/// The engine-level tests above pin the engines; this pins the *command*, which is
/// where the second open used to live: `execute` read both headers for its
/// compatibility precondition and the engine then reopened both paths for their
/// records, so the CLI required a re-openable input even though each half read the
/// stream once. Both single-pass presets are covered — `--command sort` and
/// `--command group` — because they dispatch to different engines.
///
/// Content mode is deliberately absent: it makes two passes over each input (order
/// verification, then the record comparison) and so cannot stream a FIFO at all.
#[cfg(unix)]
#[rstest]
#[case::sort_preset("sort")]
#[case::group_preset("group")]
fn test_compare_bams_cli_reads_a_non_reopenable_input(#[case] preset: &str) {
    let tmp = TempDir::new().unwrap();
    // `create_minimal_header` advertises SO:unsorted GO:query SS:template-coordinate,
    // which satisfies the sort-verify engine and the grouping engine alike.
    let header = create_minimal_header("chr1", 10000);
    let records = vec![mi_record(b"read1", 100, "1"), mi_record(b"read2", 200, "1")];

    let bam = tmp.path().join("a.bam");
    write_bam(&bam, &header, &records);

    // Both arguments are FIFOs. Serving only the first would leave the second
    // re-openable, so a regression that reopened just `bam2` — `execute` opens the two
    // independently — would still pass.
    let (fifo1, feeder1) = serve_over_fifo(tmp.path(), "a.fifo", &bam);
    let (fifo2, feeder2) = serve_over_fifo(tmp.path(), "b.fifo", &bam);

    let child = Command::new(env!("CARGO_BIN_EXE_fgumi"))
        .args(["compare", "bams"])
        .arg(&fifo1)
        .arg(&fifo2)
        .args(["--command", preset])
        .stdout(std::process::Stdio::piped())
        .stderr(std::process::Stdio::piped())
        .spawn()
        .expect("failed to spawn fgumi");
    let output = wait_bounded(child, preset);

    // Assert before joining: a regression leaves a feeder blocked in its `open`
    // forever, so joining first would hang rather than fail. See `serve_over_fifo`.
    assert!(
        output.status.success(),
        "compare bams --command {preset} rejected FIFOs it accepts as regular files:\n{}",
        String::from_utf8_lossy(&output.stderr)
    );
    join_feeder(feeder1);
    join_feeder(feeder2);

    // Each preset has its own result line — `execute_sort_verify` prints "BAM files
    // are IDENTICAL", `execute_grouping` prints "BAM groupings are EQUIVALENT". Accepting
    // either would let a dispatch regression that ran the wrong engine for the preset
    // still pass, so pin the one this preset must print.
    let expected_result_line = match preset {
        "sort" => "RESULT: BAM files are IDENTICAL",
        "group" => "RESULT: BAM groupings are EQUIVALENT",
        other => panic!("unexpected preset {other}: add its RESULT line to this match"),
    };
    let stdout = String::from_utf8_lossy(&output.stdout);
    assert!(
        stdout.contains(expected_result_line),
        "the same data over a FIFO must compare equal to itself: --command {preset} must print \
         {expected_result_line:?}, got:\n{stdout}"
    );
}

// ---------------------------------------------------------------------------
// Content mode: sort-order verification
//
// Content mode pairs records purely by position, so an `@HD`-declared order the
// records do not actually honor would silently corrupt that pairing rather than
// surfacing as an honest diff. The check runs inside the comparison's own record
// pass (`OrderCheck`, folded into `positional_compare`) rather than as a separate
// traversal of each input, and these tests pin the contract that fold has to
// preserve: that a mis-sorted input is rejected, that the reported violation
// *count* reflects the whole file rather than stopping at the first bad record,
// and that bam1 is reported before bam2 when both are mis-sorted.
// ---------------------------------------------------------------------------

/// Writes a coordinate-declaring BAM whose records appear at `positions` in the
/// given order — so a caller can declare `SO:coordinate` while emitting records
/// that violate it.
fn write_coordinate_declared_bam(path: &Path, positions: &[i32]) {
    let header = create_coordinate_sorted_header("chr1", 100_000);
    let records: Vec<RawRecord> = positions
        .iter()
        .enumerate()
        .map(|(i, pos)| mapped_record(format!("read{i:04}").as_bytes(), *pos))
        .collect();
    write_bam(path, &header, &records);
}

/// The reported violation count covers the whole file, not just the first bad
/// record. This is why the check accumulates instead of failing fast, and it is
/// the property a "simplify to fail-fast" refactor would silently break.
#[test]
fn content_mode_reports_every_sort_order_violation_not_just_the_first() {
    let tmp = TempDir::new().unwrap();
    let sorted = tmp.path().join("sorted.bam");
    let mis_sorted = tmp.path().join("mis_sorted.bam");

    write_coordinate_declared_bam(&sorted, &[100, 200, 300, 400, 500, 600]);
    // Two independent backwards steps: 150 after 200, and 350 after 400.
    write_coordinate_declared_bam(&mis_sorted, &[100, 200, 150, 400, 350, 600]);

    let (code, _stdout, stderr) = run_compare(&sorted, &mis_sorted, "content", &[]);
    assert_eq!(code, Some(1), "a mis-sorted input must fail: {stderr}");
    assert!(
        stderr.contains("2 record(s) violate"),
        "violation count must cover the whole file, not stop at the first: {stderr}"
    );
}

/// When both inputs are mis-sorted, bam1 is reported. Pins the evaluation order,
/// which is what keeps the diagnostic stable now that both files are checked in
/// one pass rather than one-then-the-other.
#[test]
fn content_mode_reports_bam1_when_both_inputs_are_mis_sorted() {
    let tmp = TempDir::new().unwrap();
    let first = tmp.path().join("first.bam");
    let second = tmp.path().join("second.bam");

    write_coordinate_declared_bam(&first, &[100, 200, 150, 400]);
    write_coordinate_declared_bam(&second, &[100, 200, 150, 400]);

    let (code, _stdout, stderr) = run_compare(&first, &second, "content", &[]);
    assert_eq!(code, Some(1), "mis-sorted inputs must fail: {stderr}");
    assert!(
        stderr.contains("first.bam") && !stderr.contains("second.bam"),
        "bam1 must be the input reported when both are mis-sorted: {stderr}"
    );
}

/// A correctly-sorted pair that declares its order must still compare normally —
/// the order check must not reject well-formed input, at any thread count. The
/// thread parameter covers the split decompression budget as well as the check.
#[rstest]
fn content_mode_accepts_input_that_honors_its_declared_sort_order(#[values(1, 4)] threads: usize) {
    let tmp = TempDir::new().unwrap();
    let a = tmp.path().join("a.bam");
    let b = tmp.path().join("b.bam");

    write_coordinate_declared_bam(&a, &[100, 200, 300, 400]);
    write_coordinate_declared_bam(&b, &[100, 200, 300, 400]);

    assert!(
        run_compare_in_process(&a, &b, "content", &["-t", &threads.to_string()]),
        "correctly sorted identical BAMs must match at --threads {threads}"
    );
}

/// Once pairing stops because one input ran out, `positional_compare` keeps draining
/// the longer side — and must keep checking its order while doing so. The drained
/// records are the ones no counterpart ever forced it to look at, so an order check
/// that stopped at the last paired record would report a violation count covering
/// only part of the file, and miss a file whose *only* mis-sorted records are in its
/// unpaired tail. The batch size is 1 so the tail genuinely arrives as extra batches
/// through the drain loop rather than inside the last paired batch.
///
/// Both directions are covered because the two sides drain through separate loops.
#[rstest]
#[case::bam1_is_longer(true)]
#[case::bam2_is_longer(false)]
fn positional_compare_checks_sort_order_in_the_unpaired_tail(#[case] bam1_is_longer: bool) {
    let tmp = TempDir::new().unwrap();
    let short = tmp.path().join("short.bam");
    let long = tmp.path().join("long.bam");

    write_coordinate_declared_bam(&short, &[100, 200]);
    // The tail (everything past record 2) steps backwards at 300 -> 250.
    write_coordinate_declared_bam(&long, &[100, 200, 300, 250, 400]);

    let (bam1, bam2) = if bam1_is_longer { (&long, &short) } else { (&short, &long) };
    let err = positional_compare(
        bam1,
        bam2,
        1,
        1,
        100,
        ContentPredicate::Exact,
        Some(fgumi_sort::SortOrder::Coordinate),
    )
    .expect_err("a violation in the unpaired tail must still be reported");

    let msg = err.to_string();
    assert!(
        msg.contains("1 record(s) violate"),
        "the drained tail must be order-checked too: {msg}"
    );
    assert!(msg.contains("long.bam"), "must name the mis-sorted input: {msg}");
}

/// Truncating a BAM mid-BGZF-block makes its reader fail partway through. That must
/// surface as an error naming the offending input, never as a short file: a
/// silently-short read turns a corrupt input into a record-count DIFFER, or a false
/// IDENTICAL when both sides are damaged alike, from the tool whose only job is
/// deciding whether two files agree.
///
/// Both positions are covered because each input is received through its own branch,
/// and both the paired phase and the drain phase are covered because they receive
/// through separate loops. Which phase sees the failure is set by how many records
/// the intact side holds, since the damaged side fails partway through 20,000:
/// one record makes the intact side reach EOF first, so the failure arrives while
/// the damaged side is being drained; 20,000 keeps both sides live, so it arrives
/// while they are still being paired. A small intact count would put every case in
/// the drain loop and leave the pairing loop untested.
#[rstest]
#[case::damaged_bam1_while_pairing(true, 20_000)]
#[case::damaged_bam2_while_pairing(false, 20_000)]
#[case::damaged_bam1_while_draining(true, 1)]
#[case::damaged_bam2_while_draining(false, 1)]
fn positional_compare_reports_a_damaged_input_rather_than_a_short_one(
    #[case] damage_bam1: bool,
    #[case] intact_records: i32,
) {
    let tmp = TempDir::new().unwrap();
    let intact = tmp.path().join("intact.bam");
    let damaged = tmp.path().join("damaged.bam");

    let positions: Vec<i32> = (0..intact_records).map(|i| 100 + i).collect();
    write_coordinate_declared_bam(&intact, &positions);

    // Enough records to span several BGZF blocks, so the cut lands well past the
    // header and the reader delivers batches before it fails.
    let long_positions: Vec<i32> = (0..20_000).map(|i| 100 + i).collect();
    let whole = tmp.path().join("whole.bam");
    write_coordinate_declared_bam(&whole, &long_positions);
    let bytes = fs::read(&whole).expect("read whole BAM");
    fs::write(&damaged, &bytes[..bytes.len() * 3 / 5]).expect("write truncated BAM");

    let (bam1, bam2) = if damage_bam1 { (&damaged, &intact) } else { (&intact, &damaged) };
    let err = positional_compare(
        bam1,
        bam2,
        1,
        1,
        100,
        ContentPredicate::Exact,
        Some(fgumi_sort::SortOrder::Coordinate),
    )
    .expect_err("a truncated input must error, not read as a shorter file");

    let msg = err.to_string();
    assert!(
        msg.contains("damaged.bam"),
        "the error must name the input that failed to read: {msg}"
    );
}
