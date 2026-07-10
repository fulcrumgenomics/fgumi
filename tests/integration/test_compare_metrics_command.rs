//! Integration tests for the `fgumi compare metrics` command's CLI surface.
//!
//! Unit tests for the row key-join engine itself live alongside the
//! implementation in `src/lib/commands/compare/metrics.rs`; these tests instead
//! exercise the command end to end: argument parsing, the `RESULT: metrics
//! files ...` report line, `--quiet`/exit-code behavior, and the multi-column
//! `--key-columns` flag via the real CLI surface.

use clap::Parser;
use fgumi_lib::commands::command::Command as FgumiCommand;
use fgumi_lib::commands::compare::{CompareMetrics, CompareMismatch};
use std::path::Path;
use std::process::Command;
use tempfile::TempDir;

use crate::helpers::write_tsv;

/// Runs `CompareMetrics::execute()` in-process; returns `true` on match, `false`
/// on `CompareMismatch` (DIFFER). Panics on any other error.
fn run_compare_in_process(file1: &Path, file2: &Path, extra_args: &[&str]) -> bool {
    let mut args: Vec<&str> = vec!["metrics", file1.to_str().unwrap(), file2.to_str().unwrap()];
    args.extend(extra_args);
    let cmd = CompareMetrics::try_parse_from(args).expect("failed to parse compare metrics args");
    match cmd.execute("fgumi compare metrics") {
        Ok(()) => true,
        Err(e) if e.is::<CompareMismatch>() => false,
        Err(e) => panic!("compare metrics hit unexpected error: {e:#}"),
    }
}

/// Runs `fgumi compare metrics` via subprocess and returns `(exit_code, stdout, stderr)`.
///
/// Returns the *exact* exit code (not just success/failure) so callers can distinguish a
/// clean DIFFER (`Some(1)`, set explicitly for `CompareMismatch` in `main`) from a usage
/// error (`Some(2)`), a panic/abort (`Some(101)`/`Some(134)`), or a signal (`None`) — a
/// nonzero exit that is *not* 1 must never be mistaken for a difference result. `code()`
/// is `None` only when the process was terminated by a signal without exiting. `stderr` is
/// returned so a caller can additionally rule out a generic runtime error (which `main`
/// also exits 1 for, but via the `anyhow` `Error:` report on stderr) masquerading as a
/// clean DIFFER when stdout is suppressed by `--quiet`.
fn run_compare_subprocess(
    file1: &Path,
    file2: &Path,
    extra_args: &[&str],
) -> (Option<i32>, String, String) {
    let output = Command::new(env!("CARGO_BIN_EXE_fgumi"))
        .args(["compare", "metrics"])
        .arg(file1)
        .arg(file2)
        .args(extra_args)
        .output()
        .expect("failed to run fgumi");
    (
        output.status.code(),
        String::from_utf8_lossy(&output.stdout).to_string(),
        String::from_utf8_lossy(&output.stderr).to_string(),
    )
}

#[test]
fn identical_files_match_in_process() {
    let dir = TempDir::new().expect("tempdir");
    let content = "umi\tcount\nAAA\t3\nCCC\t5\n";
    let f1 = write_tsv(dir.path(), "a.txt", content);
    let f2 = write_tsv(dir.path(), "b.txt", content);
    assert!(run_compare_in_process(&f1, &f2, &[]));
}

#[test]
fn reordered_rows_still_match_in_process() {
    let dir = TempDir::new().expect("tempdir");
    let f1 = write_tsv(dir.path(), "a.txt", "family_size\tcount\n1\t10\n2\t5\n3\t1\n");
    let f2 = write_tsv(dir.path(), "b.txt", "family_size\tcount\n3\t1\n1\t10\n2\t5\n");
    assert!(run_compare_in_process(&f1, &f2, &[]), "row order must not affect the key-join");
}

#[test]
fn multi_column_key_columns_flag_joins_correctly() {
    let dir = TempDir::new().expect("tempdir");
    let f1 = write_tsv(dir.path(), "a.txt", "ab_size\tba_size\tcount\n1\t0\t100\n2\t0\t20\n");
    let f2 = write_tsv(dir.path(), "b.txt", "ab_size\tba_size\tcount\n2\t0\t20\n1\t0\t100\n");
    assert!(run_compare_in_process(&f1, &f2, &["--key-columns", "ab_size,ba_size"]));
}

/// Exercises the CLI's exit code and the `RESULT: metrics files DIFFER` report
/// line, which downstream tooling (per the `compare bams` precedent) greps for.
#[test]
fn subprocess_reports_result_differ_and_exit_code_1() {
    let dir = TempDir::new().expect("tempdir");
    let f1 = write_tsv(dir.path(), "a.txt", "key\tvalue\nA\t1\n");
    let f2 = write_tsv(dir.path(), "b.txt", "key\tvalue\nA\t2\n");
    let (code, stdout, _stderr) = run_compare_subprocess(&f1, &f2, &[]);
    assert_eq!(
        code,
        Some(1),
        "a clean DIFFER must exit exactly 1, not crash/usage-error: {stdout}"
    );
    assert!(stdout.contains("RESULT: metrics files DIFFER"), "stdout was: {stdout}");
    assert!(!stdout.contains("only in"), "no presence-only rows in this fixture: {stdout}");
}

#[test]
fn subprocess_reports_result_identical_with_verbose() {
    let dir = TempDir::new().expect("tempdir");
    let content = "key\tvalue\nA\t1\n";
    let f1 = write_tsv(dir.path(), "a.txt", content);
    let f2 = write_tsv(dir.path(), "b.txt", content);
    let (code, stdout, _stderr) = run_compare_subprocess(&f1, &f2, &["--verbose"]);
    assert_eq!(code, Some(0), "identical metrics files must exit exactly 0: {stdout}");
    assert!(stdout.contains("RESULT: metrics files are IDENTICAL"), "stdout was: {stdout}");
}

#[test]
fn subprocess_quiet_suppresses_output_on_differ() {
    let dir = TempDir::new().expect("tempdir");
    let f1 = write_tsv(dir.path(), "a.txt", "key\tvalue\nA\t1\n");
    let f2 = write_tsv(dir.path(), "b.txt", "key\tvalue\nA\t2\n");
    let (code, stdout, stderr) = run_compare_subprocess(&f1, &f2, &["--quiet"]);
    assert_eq!(
        code,
        Some(1),
        "a quiet DIFFER must exit exactly 1, not crash/usage-error: {stderr}"
    );
    assert!(stdout.is_empty(), "quiet mode must suppress stdout: {stdout}");
    // A generic runtime error also exits 1 but reports via `anyhow`'s `Error:` line on
    // stderr; assert its absence so a real error can't masquerade as a clean quiet DIFFER.
    assert!(!stderr.contains("Error:"), "quiet DIFFER must not emit a runtime error: {stderr}");
}

#[test]
fn localized_differ_names_the_offending_key() {
    let dir = TempDir::new().expect("tempdir");
    let f1 = write_tsv(dir.path(), "a.txt", "family_size\tcount\n1\t10\n2\t5\n3\t1\n");
    let f2 = write_tsv(dir.path(), "b.txt", "family_size\tcount\n1\t10\n2\t5\n");
    let (code, stdout, _stderr) = run_compare_subprocess(&f1, &f2, &[]);
    assert_eq!(code, Some(1), "a clean DIFFER must exit exactly 1: {stdout}");
    assert!(stdout.contains("3 only in file1"), "stdout was: {stdout}");
}
