//! Progress-bar non-regression tests.
//!
//! `runall` is invoked as a subprocess with stderr redirected (via
//! [`std::process::Command::output`]), which means stderr is NOT a TTY. In
//! that environment [`indicatif::MultiProgress`] must pick a hidden draw
//! target and emit none of the cursor-movement escape sequences that a
//! rendered bar would use. This test guards that contract: if a future
//! change accidentally forces a visible bar in non-TTY mode, the cursor
//! escapes would leak into piped logs and CI output.
//!
//! Note: the `tracing` subscriber emits plain color escapes (`ESC[32m`,
//! `ESC[0m`, ...). Those are benign for log scrapers and are not the
//! target of this test. We specifically look for the cursor-movement
//! sequences indicatif uses to overwrite previous bar frames in place.

use std::process::Command;

use crate::helpers::grouped_bam_fixture::build_small_grouped_bam;

fn fgumi_bin() -> &'static str {
    env!("CARGO_BIN_EXE_fgumi")
}

/// Cursor-movement escape sequences that a rendered progress bar uses to
/// overwrite itself. Tracing color codes do not include these.
///
/// - `\x1b[2K` clears the current line
/// - `\x1b[1A` / `\x1b[NA` move the cursor up N lines
/// - `\r` alone is not sufficient (tracing does not emit it but neither is
///   it indicatif-exclusive), so we require the pair `\r\x1b[2K` which
///   indicatif emits from `InMemoryTerm::clear_line` and the real
///   `Term::clear_line` alike.
const INDICATIF_CLEAR_LINE: &str = "\r\u{1b}[2K";
const INDICATIF_CURSOR_UP_PREFIX: &str = "\u{1b}[";

/// Returns true if `stderr` contains any cursor-movement escape that only
/// a rendered progress bar would produce.
fn contains_progress_bar_escapes(stderr: &str) -> bool {
    if stderr.contains(INDICATIF_CLEAR_LINE) {
        return true;
    }
    // Look for `ESC [ <digit>+ A` (cursor up). Not foolproof but close.
    for (start, _) in stderr.match_indices(INDICATIF_CURSOR_UP_PREFIX) {
        // Peek at the next up-to-4 characters after the prefix.
        let tail = &stderr[start + INDICATIF_CURSOR_UP_PREFIX.len()..];
        let mut digits = 0usize;
        for ch in tail.chars().take(4) {
            if ch.is_ascii_digit() {
                digits += 1;
            } else if ch == 'A' && digits >= 1 {
                return true;
            } else {
                break;
            }
        }
    }
    false
}

#[test]
fn runall_non_tty_does_not_emit_ansi_escape_codes() {
    let tmp = tempfile::tempdir().unwrap();
    let input = tmp.path().join("input.bam");
    build_small_grouped_bam(&input);
    let output = tmp.path().join("output.bam");

    let result = Command::new(fgumi_bin())
        .args(["runall", "--start-from", "group", "--threads", "1"])
        .args(["--filter::min-reads", "1"])
        .arg("--input")
        .arg(&input)
        .arg("--output")
        .arg(&output)
        // Belt-and-suspenders: even if indicatif's TTY detection ever
        // misbehaves (e.g. on a platform where Term features look
        // color-capable without a real TTY), this env var forces the
        // draw target to hidden. The test passes without the env var
        // too under normal subprocess-capture semantics — set here to
        // avoid a flaky CI agent that exposes an interactive-shell
        // descriptor.
        .env("FGUMI_HIDE_PROGRESS", "1")
        .output()
        .expect("spawn fgumi runall");

    assert!(
        result.status.success(),
        "fgumi runall failed: {}\nstderr:\n{}",
        result.status,
        String::from_utf8_lossy(&result.stderr),
    );

    let stderr = String::from_utf8_lossy(&result.stderr);
    assert!(
        !contains_progress_bar_escapes(&stderr),
        "stderr contains progress-bar escape sequences (bars leaked into non-TTY output):\n{stderr}",
    );
}

/// Without the env var, the subprocess pipe should also be detected as
/// non-TTY by indicatif's own TTY check. Keeping this separate exercises
/// indicatif's built-in behavior — the env var in the first test masks
/// any upstream regression there.
#[test]
fn runall_subprocess_pipe_suppresses_bars_by_default() {
    let tmp = tempfile::tempdir().unwrap();
    let input = tmp.path().join("input.bam");
    build_small_grouped_bam(&input);
    let output = tmp.path().join("output.bam");

    let result = Command::new(fgumi_bin())
        .args(["runall", "--start-from", "group", "--threads", "1"])
        .args(["--filter::min-reads", "1"])
        .arg("--input")
        .arg(&input)
        .arg("--output")
        .arg(&output)
        // Ensure CI indicator is NOT set — otherwise we'd be testing the
        // same forced-hide path as the prior test. We want indicatif's
        // own TTY probe to decide.
        .env_remove("CI")
        .env_remove("FGUMI_HIDE_PROGRESS")
        .output()
        .expect("spawn fgumi runall");

    assert!(
        result.status.success(),
        "fgumi runall failed: {}\nstderr:\n{}",
        result.status,
        String::from_utf8_lossy(&result.stderr),
    );

    let stderr = String::from_utf8_lossy(&result.stderr);
    assert!(
        !contains_progress_bar_escapes(&stderr),
        "stderr contains progress-bar escape sequences under default TTY probe:\n{stderr}",
    );
}
