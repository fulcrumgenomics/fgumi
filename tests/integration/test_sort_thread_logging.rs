//! Integration tests for what `fgumi sort` reports about its thread counts and
//! the memory budget they scale.
//!
//! `--sort-threads` and `--merge-threads` are caps within `--threads` (each
//! defaults to it). These tests pin the reported counts to the ones the phases
//! actually use, the clamp warning for an override above `--threads`, and the
//! per-thread memory budget to the effective sort-phase count.

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
/// used: overrides above `--threads` are clamped to it, below are honoured.
#[rstest]
#[case::overrides_above_unset_threads_clamp_to_one(
    &["--sort-threads", "8", "--merge-threads", "16"],
    "1"
)]
#[case::sort_override_below_threads(&["--threads", "6", "--sort-threads", "2"], "sort 2, merge 6")]
#[case::merge_override_below_threads(&["--threads", "6", "--merge-threads", "3"], "sort 6, merge 3")]
#[case::matching_phase_counts(&["--threads", "4"], "4")]
#[case::sort_above_threads_clamps(&["--threads", "4", "--sort-threads", "32"], "4")]
fn reported_thread_counts_match_the_phases(#[case] extra_args: &[&str], #[case] expected: &str) {
    let stderr = sort_and_capture_logs(extra_args);
    let reported = threads_lines(&stderr);
    assert!(!reported.is_empty(), "no `Threads:` line was logged:\n{stderr}");
    for counts in &reported {
        assert_eq!(counts, expected, "misreported thread counts:\n{stderr}");
    }
}

/// An override above `--threads` warns exactly once per flag, naming the
/// requested and effective counts; an override within `--threads` is silent.
#[rstest]
#[case::sort_clamped(&["--threads", "4", "--sort-threads", "32"], &[
    "--sort-threads 32 exceeds --threads 4; the sort phase cannot use more workers than the pool has, using 4 (raise --threads to widen the pool)",
])]
#[case::both_clamped(&["--threads", "2", "--sort-threads", "8", "--merge-threads", "8"], &[
    "--sort-threads 8 exceeds --threads 2; the sort phase cannot use more workers than the pool has, using 2 (raise --threads to widen the pool)",
    "--merge-threads 8 exceeds --threads 2; the merge phase cannot use more workers than the pool has, using 2 (raise --threads to widen the pool)",
])]
#[case::within_threads_is_silent(&["--threads", "4", "--sort-threads", "2", "--merge-threads", "4"], &[])]
fn clamped_overrides_warn_once_per_flag(#[case] extra_args: &[&str], #[case] expected: &[&str]) {
    let stderr = sort_and_capture_logs(extra_args);
    let warnings: Vec<&str> = stderr
        .lines()
        .filter(|l| {
            l.contains(" WARN ") && (l.contains("sort-threads") || l.contains("merge-threads"))
        })
        .map(|l| l.split_once("] ").map_or(l, |(_, msg)| msg))
        .collect();
    assert_eq!(warnings, expected, "clamp warnings:\n{stderr}");
}

/// The `Max memory:` line names the multiplier the budget resolved with: the
/// effective sort-phase count, attributed to `--sort-threads` only when it
/// lowered the count. The resolved byte count is pinned independently by
/// `test_memory_budget_follows_effective_phase1`; `expected_line` is asserted
/// separately because the total goes through `bytesize`'s `Display`.
#[rstest]
#[case::sort_threads_above_threads_clamps(
    &["--max-memory", "100M", "--threads", "16", "--sort-threads", "32"],
    "MiB/thread x 16 threads, from --threads)",
    "Max memory: 1.5 GiB (95.4 MiB/thread x 16 threads, from --threads)"
)]
#[case::sort_threads_below_threads_shrinks(
    &["--max-memory", "100M", "--threads", "4", "--sort-threads", "2"],
    "MiB/thread x 2 threads, from --sort-threads)",
    "Max memory: 190.7 MiB (95.4 MiB/thread x 2 threads, from --sort-threads)"
)]
#[case::sort_threads_zero_is_one(
    &["--max-memory", "100M", "--threads", "4", "--sort-threads", "0"],
    "MiB/thread x 1 threads, from --sort-threads)",
    "Max memory: 95.4 MiB (95.4 MiB/thread x 1 threads, from --sort-threads)"
)]
fn reported_memory_budget_scales_by_the_effective_sort_phase(
    #[case] extra_args: &[&str],
    #[case] expected_multiplier: &str,
    #[case] expected_line: &str,
) {
    let stderr = sort_and_capture_logs(extra_args);
    assert!(
        stderr.contains(expected_multiplier),
        "budget did not scale by the effective sort-phase count:\n{stderr}"
    );
    assert!(
        stderr.contains(expected_line),
        "unexpected rendering of the budget (bytesize formatting change?):\n{stderr}"
    );
}
