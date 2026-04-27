//! Every subcommand's --help includes an EXAMPLES section.
//!
//! Users copy-paste from --help, so each subcommand must advertise at least
//! one concrete invocation. This smoke test spawns `fgumi <sub> --help` for
//! every non-feature-gated subcommand and asserts that the token
//! "EXAMPLES:" appears in the rendered help output.

use std::process::Command;

/// Path to the `fgumi` binary produced for this integration test.
fn fgumi_bin() -> &'static str {
    env!("CARGO_BIN_EXE_fgumi")
}

/// All non-feature-gated subcommands of `fgumi`. `compare` and `simulate`
/// are excluded because they are gated behind Cargo features and may not
/// be compiled into every build.
const SUBCOMMANDS: &[&str] = &[
    "extract",
    "fastq",
    "correct",
    "group",
    "dedup",
    "simplex",
    "duplex",
    "codec",
    "filter",
    "clip",
    "sort",
    "zipper",
    "merge",
    "downsample",
    "duplex-metrics",
    "simplex-metrics",
    "review",
    "runall",
];

#[test]
fn every_subcommand_has_examples_section() {
    let mut missing: Vec<&str> = Vec::new();
    for sub in SUBCOMMANDS {
        let out = Command::new(fgumi_bin())
            .args([sub, "--help"])
            .output()
            .unwrap_or_else(|e| panic!("failed to spawn `fgumi {sub} --help`: {e}"));
        assert!(
            out.status.success(),
            "`fgumi {sub} --help` exited non-zero: {}",
            String::from_utf8_lossy(&out.stderr)
        );
        let stdout = String::from_utf8_lossy(&out.stdout);
        if !stdout.contains("EXAMPLES:") {
            missing.push(*sub);
        }
    }
    assert!(missing.is_empty(), "subcommands missing EXAMPLES: in --help output: {missing:?}");
}
