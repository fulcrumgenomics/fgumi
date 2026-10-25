//! Helpers for spawning the compiled `fgumi` binary and capturing its stderr.
//!
//! A handful of tests need to observe a real logged banner or warning line
//! (not just a command's exit status or output BAM), which requires spawning
//! the actual binary rather than calling `Command::execute()` in-process: the
//! in-process log-capture harness (`fgumi_lib::commands::common::test_log_capture`)
//! is `pub(crate)` inside `fgumi_lib` and unreachable from this external
//! integration-test crate.

use std::ffi::OsStr;
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
    let mut args: Vec<&OsStr> = vec![subcommand.as_ref(), "-i".as_ref(), input.as_os_str()];
    args.extend([OsStr::new("-o"), output.as_os_str()]);
    args.extend(extra_args.iter().map(OsStr::new));
    run_args_and_capture_logs(&args)
}

/// `path` as a `&str`, for a command-line argument; panics if it is not
/// valid UTF-8.
pub fn utf8_path(path: &std::path::Path) -> &str {
    path.to_str().expect("test path is valid UTF-8")
}

/// Spawn `fgumi <args>...` and return its captured stderr, asserting the run
/// succeeded. For command lines that do not fit [`run_and_capture_logs`]'s
/// `-i`/`-o` shape (e.g. `fgumi sort --verify`, which takes no output).
///
/// `RUST_LOG=info` is set so the run's `info!`/`warn!` banners reach stderr.
pub fn run_args_and_capture_logs(args: &[&OsStr]) -> String {
    let display = args.iter().map(|a| a.to_string_lossy()).collect::<Vec<_>>().join(" ");
    let result = std::process::Command::new(env!("CARGO_BIN_EXE_fgumi"))
        .env("RUST_LOG", "info")
        .args(args)
        .output()
        .unwrap_or_else(|e| panic!("failed to spawn `fgumi {display}`: {e}"));
    let stderr = String::from_utf8_lossy(&result.stderr).into_owned();
    assert!(result.status.success(), "fgumi {display} failed:\n{stderr}");
    stderr
}
