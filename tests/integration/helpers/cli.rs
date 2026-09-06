//! Helpers for spawning the compiled `fgumi` binary and capturing its stderr.
//!
//! A handful of tests need to observe a real logged banner or warning line
//! (not just a command's exit status or output BAM), which requires spawning
//! the actual binary rather than calling `Command::execute()` in-process: the
//! in-process log-capture harness (`fgumi_lib::commands::common::test_log_capture`)
//! is `pub(crate)` inside `fgumi_lib` and unreachable from this external
//! integration-test crate.

use std::path::Path;

/// Spawn `fgumi <subcommand> -i <input> -o <output> <extra_args>...` and
/// return its captured stderr, asserting the run succeeded.
///
/// `RUST_LOG=info` is set so the run's `info!`/`warn!` banners reach stderr.
pub fn run_and_capture_logs(
    subcommand: &str,
    input: &Path,
    output: &Path,
    extra_args: &[&str],
) -> String {
    let result = std::process::Command::new(env!("CARGO_BIN_EXE_fgumi"))
        .env("RUST_LOG", "info")
        .args([subcommand, "-i"])
        .arg(input)
        .arg("-o")
        .arg(output)
        .args(extra_args)
        .output()
        .unwrap_or_else(|e| panic!("failed to spawn `fgumi {subcommand}`: {e}"));
    let stderr = String::from_utf8_lossy(&result.stderr).into_owned();
    assert!(result.status.success(), "fgumi {subcommand} failed:\n{stderr}");
    stderr
}
