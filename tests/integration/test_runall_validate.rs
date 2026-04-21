//! Validation tests for `fgumi runall` pre-flight checks.
//!
//! Exercises `Runall::validate()` via the CLI surface: verifies that
//! cross-parameter invariants (e.g. `--filter::min-reads` required for plans
//! that run the consensus / filter stages) are enforced before any expensive
//! work is started.

use std::process::Command;
use tempfile::TempDir;

/// Path to the compiled `fgumi` binary built by this crate's test harness.
fn fgumi_bin() -> &'static str {
    env!("CARGO_BIN_EXE_fgumi")
}

/// Run `fgumi runall` with the given args and return (success, stderr).
fn run_runall(args: &[&str]) -> (bool, String) {
    let output =
        Command::new(fgumi_bin()).arg("runall").args(args).output().expect("spawn fgumi runall");
    let stderr = String::from_utf8_lossy(&output.stderr).to_string();
    (output.status.success(), stderr)
}

/// When the plan runs consensus or filter, `--filter::min-reads` must be set.
/// Absence is a validation error surfaced before the pipeline spins up.
#[test]
fn runall_stop_after_filter_without_min_reads_fails_validation() {
    let tmp = TempDir::new().unwrap();
    let in_bam = tmp.path().join("in.bam");
    std::fs::write(&in_bam, b"not a real bam").unwrap();
    let out_bam = tmp.path().join("out.bam");

    let (ok, stderr) = run_runall(&[
        "--start-from",
        "group",
        "--stop-after",
        "filter",
        "--input",
        in_bam.to_str().unwrap(),
        "--output",
        out_bam.to_str().unwrap(),
        "--group::strategy",
        "adjacency",
    ]);
    assert!(
        !ok,
        "expected validation failure when --filter::min-reads is missing on a filter plan"
    );
    assert!(
        stderr.contains("--filter::min-reads is required"),
        "stderr should explain the missing flag; got:\n{stderr}"
    );
}

/// Same as above but for `--stop-after consensus` (the stage also consumes
/// `min_reads` via `VanillaUmiConsensusOptions`).
#[test]
fn runall_stop_after_consensus_without_min_reads_fails_validation() {
    let tmp = TempDir::new().unwrap();
    let in_bam = tmp.path().join("in.bam");
    std::fs::write(&in_bam, b"not a real bam").unwrap();
    let out_bam = tmp.path().join("out.bam");

    let (ok, stderr) = run_runall(&[
        "--start-from",
        "group",
        "--stop-after",
        "consensus",
        "--input",
        in_bam.to_str().unwrap(),
        "--output",
        out_bam.to_str().unwrap(),
        "--group::strategy",
        "adjacency",
    ]);
    assert!(
        !ok,
        "expected validation failure when --filter::min-reads is missing on a consensus plan"
    );
    assert!(
        stderr.contains("--filter::min-reads is required"),
        "stderr should explain the missing flag; got:\n{stderr}"
    );
}

/// Plans that stop before consensus / filter do NOT require `--filter::min-reads`.
/// Validation should not fail for a missing `min_reads` flag in that case. We
/// don't run the full pipeline to EOF here — we just look for the
/// min-reads validation error specifically, since the BAM is a stub.
#[test]
fn runall_stop_after_group_does_not_require_min_reads() {
    let tmp = TempDir::new().unwrap();
    let in_bam = tmp.path().join("in.bam");
    std::fs::write(&in_bam, b"not a real bam").unwrap();
    let out_bam = tmp.path().join("out.bam");

    let (_ok, stderr) = run_runall(&[
        "--start-from",
        "group",
        "--stop-after",
        "group",
        "--input",
        in_bam.to_str().unwrap(),
        "--output",
        out_bam.to_str().unwrap(),
        "--group::strategy",
        "adjacency",
    ]);
    assert!(
        !stderr.contains("--filter::min-reads is required"),
        "min-reads must not be required for --stop-after group; stderr:\n{stderr}"
    );
}
