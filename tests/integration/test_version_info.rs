//! Tests that `fgumi --version` emits an enriched multi-line build banner
//! containing crate version, git commit, build date, rustc version, target
//! triple, enabled features, and allocator.
//!
//! This output is invaluable for bug reports — it pins down the exact binary
//! that produced a given result. The integration contract:
//!   * `-V` prints a single line `fgumi <semver>` (clap default short form).
//!   * `--version` prints the long form with all build metadata labeled.

use std::process::Command;

fn fgumi_bin() -> &'static str {
    env!("CARGO_BIN_EXE_fgumi")
}

#[test]
fn short_dash_v_is_one_liner() {
    let out = Command::new(fgumi_bin()).arg("-V").output().expect("spawn fgumi -V");
    assert!(out.status.success(), "fgumi -V exited non-zero");
    let stdout = String::from_utf8_lossy(&out.stdout);
    // Single logical line: one trailing newline, no embedded ones.
    let trimmed = stdout.trim_end_matches('\n');
    assert!(!trimmed.contains('\n'), "expected single-line -V output, got: {stdout:?}");
    assert!(trimmed.starts_with("fgumi "), "expected 'fgumi ' prefix, got: {trimmed:?}");
}

#[test]
fn long_dash_dash_version_includes_all_metadata() {
    let out = Command::new(fgumi_bin()).arg("--version").output().expect("spawn fgumi --version");
    assert!(out.status.success(), "fgumi --version exited non-zero");
    let stdout = String::from_utf8_lossy(&out.stdout);

    // First line: "fgumi <version>".
    let first = stdout.lines().next().expect("--version produced no output");
    assert!(first.starts_with("fgumi "), "expected 'fgumi ' prefix, got: {first:?}");

    // Metadata fields (labels are defined in src/lib/version.rs::LONG_VERSION).
    for label in ["commit:", "built:", "rustc:", "target:", "profile:", "features:", "allocator:"] {
        assert!(
            stdout.contains(label),
            "--version output missing '{label}' field; full output:\n{stdout}"
        );
    }

    // Allocator is always mimalloc.
    let expected_allocator = "mimalloc";
    assert!(
        stdout.contains(expected_allocator),
        "--version output missing allocator '{expected_allocator}'; full output:\n{stdout}"
    );
}
