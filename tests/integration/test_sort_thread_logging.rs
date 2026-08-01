//! Integration tests for the thread counts `fgumi sort` reports.
//!
//! `--sort-threads` and `--merge-threads` each default to `--threads`, so the
//! flag alone does not describe a run that set only the per-phase overrides.
//! These tests pin the reported counts to the ones the phases actually use.

use rstest::rstest;
use std::fmt::Write as _;
use std::path::{Path, PathBuf};
use std::process::Command;
use tempfile::TempDir;

/// Writes a small unsorted SAM and returns its path.
///
/// SAM keeps the fixture readable and avoids a samtools dependency; sort accepts
/// it directly.
fn write_unsorted_sam(dir: &Path) -> PathBuf {
    let mut sam = String::from("@HD\tVN:1.6\tSO:unsorted\n@SQ\tSN:chr1\tLN:10000\n");
    for pos in [500, 100, 400, 200, 300] {
        writeln!(sam, "q{pos}\t0\tchr1\t{pos}\t60\t10M\t*\t0\t0\tACGTACGTAC\tIIIIIIIIII")
            .expect("write SAM record");
    }
    let path = dir.join("unsorted.sam");
    std::fs::write(&path, sam).expect("write SAM fixture");
    path
}

/// Runs `fgumi sort` at info verbosity and returns its stderr.
///
/// Pinned to info rather than `-v` so the engine's debug config dump — which
/// reports the `--threads` default alongside its own phase breakdown — stays out
/// of the captured lines.
fn sort_and_capture_logs(extra_args: &[&str]) -> String {
    let tmp = TempDir::new().expect("tempdir");
    let input = write_unsorted_sam(tmp.path());
    let output = tmp.path().join("sorted.bam");

    let result = Command::new(env!("CARGO_BIN_EXE_fgumi"))
        .env("RUST_LOG", "info")
        .args(["sort", "-i"])
        .arg(&input)
        .arg("-o")
        .arg(&output)
        .args(["--order", "coordinate"])
        .args(extra_args)
        .output()
        .expect("run fgumi sort");

    let stderr = String::from_utf8_lossy(&result.stderr).into_owned();
    assert!(result.status.success(), "fgumi sort failed:\n{stderr}");
    stderr
}

/// Every `Threads:` line the run logged, stripped of its log prefix.
fn threads_lines(stderr: &str) -> Vec<String> {
    stderr
        .lines()
        .filter_map(|line| line.split_once("Threads: "))
        .map(|(_, counts)| counts.trim().to_string())
        .collect()
}

/// Every `Threads:` line a run logs must name the counts its phases actually
/// used.
///
/// * `per_phase_overrides_without_threads` — overrides with `--threads` left
///   alone used to report `Threads: 1` while the phases ran 8 and 16 workers.
/// * `sort_override_below_threads` — a phase override below `--threads` is the
///   documented use (cede cores to an upstream producer, keep the merge wide),
///   so both counts must still show.
/// * `matching_phase_counts` — with no overrides the phases agree, so the report
///   stays a single count rather than repeating it twice.
#[rstest]
#[case::per_phase_overrides_without_threads(
    &["--sort-threads", "8", "--merge-threads", "16"],
    "sort 8, merge 16"
)]
#[case::sort_override_below_threads(&["--threads", "6", "--sort-threads", "2"], "sort 2, merge 6")]
#[case::matching_phase_counts(&["--threads", "4"], "4")]
fn reported_thread_counts_match_the_phases(#[case] extra_args: &[&str], #[case] expected: &str) {
    let stderr = sort_and_capture_logs(extra_args);
    let reported = threads_lines(&stderr);

    assert!(!reported.is_empty(), "no `Threads:` line was logged:\n{stderr}");
    for counts in &reported {
        assert_eq!(counts, expected, "misreported thread counts:\n{stderr}");
    }
}
