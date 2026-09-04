//! Real-aligner-subprocess helpers shared by the `runall` integration tests.
//!
//! `runall`'s Align-stage tests (`test_runall_command.rs`,
//! `test_runall_chain_transitions.rs`) both run end to end only when a real
//! aligner binary is on `PATH`, and validate the assembled spec/reject with a
//! clear message otherwise; this guard used to be duplicated verbatim in both
//! files.

use std::path::Path;
use std::process::{Command, Stdio};

/// Returns the first real aligner binary found on `PATH` (`bwa-mem3`
/// preferred, then classic `bwa`), or `None` if neither is installed.
pub fn aligner_binary() -> Option<&'static str> {
    ["bwa-mem3", "bwa"].into_iter().find(|bin| which::which(bin).is_ok())
}

/// Runs `<binary> index <reference>`, panicking on failure. Only called after
/// the caller has confirmed `binary` is on `PATH` via [`aligner_binary`], so a
/// failure here is a real setup regression, not a skip condition.
pub fn build_aligner_index(reference: &Path, binary: &str) {
    let status = Command::new(binary)
        .args(["index", reference.to_str().expect("reference path is valid UTF-8")])
        .stdout(Stdio::null())
        .stderr(Stdio::null())
        .status()
        .unwrap_or_else(|e| panic!("failed to run `{binary} index`: {e}"));
    assert!(status.success(), "`{binary} index` failed with status {status}");
}
