//! Integration tests for `fgumi runall` dispatch errors.
//!
//! The PR 1 foundation only registers the `NotImplementedRunner`
//! fallback, so every `(start_from, stop_after)` pair that passes
//! validation must exit non-zero with the documented "not yet
//! implemented" stderr. These tests guard against the dispatch
//! layer silently swallowing the error, or creating partial output
//! files that don't get cleaned up.

use std::path::{Path, PathBuf};
use std::process::Command;

use tempfile::TempDir;

fn fgumi_bin() -> &'static str {
    env!("CARGO_BIN_EXE_fgumi")
}

/// Build a placeholder input BAM that exists on disk so validation
/// passes. The contents are not real BAM bytes; the chain runner
/// fails before attempting to read them.
fn placeholder_input(dir: &Path, name: &str) -> PathBuf {
    let p = dir.join(name);
    std::fs::write(&p, b"PLACEHOLDER").unwrap();
    p
}

#[test]
fn unsupported_chain_exits_with_not_yet_implemented_error() {
    let dir = TempDir::new().unwrap();
    let input = placeholder_input(dir.path(), "input.bam");
    let output = dir.path().join("out.bam");

    let result = Command::new(fgumi_bin())
        .args(["runall", "--start-from", "group", "--stop-after", "group"])
        .arg("--input")
        .arg(&input)
        .arg("--output")
        .arg(&output)
        .output()
        .expect("spawn fgumi");

    assert!(
        !result.status.success(),
        "runall foundation must reject every chain shape with an error",
    );
    let stderr = String::from_utf8_lossy(&result.stderr);
    assert!(
        stderr.contains("not yet implemented"),
        "stderr must contain the documented error; got: {stderr}",
    );

    // Atomic output cleanup: the .tmp must not be left behind, and
    // the final --output must not exist (the chain runner errored
    // before producing anything renameable).
    let tmp = output.with_extension("bam.tmp");
    assert!(!tmp.exists(), "tmp output {} must be cleaned up", tmp.display());
    assert!(!output.exists(), "final output {} must not exist on error", output.display());
}

#[test]
fn invalid_chain_shape_is_rejected_by_validation() {
    let dir = TempDir::new().unwrap();
    let input = placeholder_input(dir.path(), "input.bam");
    let output = dir.path().join("out.bam");

    // Group → Extract is impossible (Extract is "before" Group in the
    // pipeline order). Validation rejects this before any chain
    // dispatch happens.
    let result = Command::new(fgumi_bin())
        .args(["runall", "--start-from", "group", "--stop-after", "extract"])
        .arg("--input")
        .arg(&input)
        .arg("--output")
        .arg(&output)
        .output()
        .expect("spawn fgumi");

    assert!(!result.status.success());
    let stderr = String::from_utf8_lossy(&result.stderr);
    assert!(
        stderr.contains("stop-after must be at or after the start stage"),
        "expected validation error; got: {stderr}",
    );
}

#[test]
fn missing_input_file_is_rejected_by_validation() {
    let dir = TempDir::new().unwrap();
    let output = dir.path().join("out.bam");

    let result = Command::new(fgumi_bin())
        .args(["runall", "--start-from", "group", "--stop-after", "group"])
        .arg("--input")
        .arg(dir.path().join("does_not_exist.bam"))
        .arg("--output")
        .arg(&output)
        .output()
        .expect("spawn fgumi");

    assert!(!result.status.success());
    let stderr = String::from_utf8_lossy(&result.stderr);
    assert!(stderr.contains("Input BAM"), "expected file-existence error; got: {stderr}",);
}
