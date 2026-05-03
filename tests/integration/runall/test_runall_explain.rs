//! Integration tests for `fgumi runall --explain`.
//!
//! `--explain` short-circuits before validation, file I/O, and signal
//! handler installation. It must exit zero, print a chain plan to
//! stdout, and never create the final or `.tmp` output files.

use std::path::{Path, PathBuf};
use std::process::Command;

use tempfile::TempDir;

fn fgumi_bin() -> &'static str {
    env!("CARGO_BIN_EXE_fgumi")
}

fn assert_no_output_artifacts(output: &Path) {
    let tmp = output.with_extension("bam.tmp");
    assert!(!tmp.exists(), "explain must not create tmp output: {}", tmp.display());
    assert!(!output.exists(), "explain must not create final output: {}", output.display());
}

#[test]
fn explain_prints_chain_plan_for_unimplemented_chain() {
    let dir = TempDir::new().unwrap();
    let output: PathBuf = dir.path().join("out.bam");

    let result = Command::new(fgumi_bin())
        .args(["runall", "--explain", "--start-from", "extract", "--stop-after", "filter"])
        .arg("--input")
        .arg("ignored.fq.gz")
        .arg("--output")
        .arg(&output)
        .output()
        .expect("spawn");

    assert!(
        result.status.success(),
        "--explain must exit zero; stderr={}",
        String::from_utf8_lossy(&result.stderr),
    );
    let stdout = String::from_utf8(result.stdout).expect("utf8 stdout");
    assert!(stdout.contains("fgumi runall"), "missing header in: {stdout}");
    assert!(stdout.contains("Extract"), "missing start stage: {stdout}");
    assert!(stdout.contains("Filter"), "missing stop stage: {stdout}");
    assert!(
        stdout.contains("not yet implemented"),
        "unimplemented chain must be flagged: {stdout}",
    );
    assert!(stdout.contains("extract -> correct"), "missing stage chain: {stdout}");

    assert_no_output_artifacts(&output);
}

#[test]
fn explain_does_not_validate_input_paths() {
    // The input does NOT exist on disk; --explain must succeed anyway.
    let dir = TempDir::new().unwrap();
    let output = dir.path().join("out.bam");

    let result = Command::new(fgumi_bin())
        .args(["runall", "--explain", "--start-from", "group", "--stop-after", "group"])
        .arg("--input")
        .arg("does_not_exist.bam")
        .arg("--output")
        .arg(&output)
        .output()
        .expect("spawn");

    assert!(
        result.status.success(),
        "--explain must skip input-existence checks; stderr={}",
        String::from_utf8_lossy(&result.stderr),
    );
}

#[test]
fn explain_works_with_smaller_chain_shape() {
    let dir = TempDir::new().unwrap();
    let output = dir.path().join("out.bam");

    let result = Command::new(fgumi_bin())
        .args(["runall", "--explain", "--start-from", "extract", "--stop-after", "extract"])
        .arg("--input")
        .arg("r1.fq")
        .arg("--output")
        .arg(&output)
        .output()
        .expect("spawn");

    assert!(result.status.success());
    let stdout = String::from_utf8(result.stdout).expect("utf8 stdout");
    assert!(stdout.contains("Extract"));
    // Only one stage — the chain should not contain "->"
    let chain_line =
        stdout.lines().find(|l| l.contains("stage chain")).expect("must have stage chain line");
    assert!(!chain_line.contains("->"), "single-stage chain must not have arrow: {chain_line}");
}
