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

/// Whether `bwa-mem3` is on `PATH` for a real-aligner test. Without it the test
/// skips (returns `false`), unless `FGUMI_BWA_MEM3_REQUIRE_TOOLS` is set (the CI
/// job that installs bwa-mem3 sets it), which makes a missing binary a failure.
pub fn bwa_mem3_or_skip(test: &str) -> bool {
    if which::which("bwa-mem3").is_ok() {
        return true;
    }
    assert!(
        std::env::var_os("FGUMI_BWA_MEM3_REQUIRE_TOOLS").is_none(),
        "FGUMI_BWA_MEM3_REQUIRE_TOOLS is set but no bwa-mem3 is on PATH ({test})"
    );
    eprintln!("no bwa-mem3 on PATH: skipping {test}");
    false
}

/// The bwa-mem3 indexed reference at `FGUMI_BWA_MEM3_TEST_REF` (the FASTA, with
/// its index files and `.dict` beside it) for a real-aligner test, or `None` to
/// skip it. Like [`bwa_mem3_or_skip`], a missing one fails instead when
/// `FGUMI_BWA_MEM3_REQUIRE_TOOLS` is set (the CI job that builds it).
pub fn bwa_mem3_test_ref(test: &str) -> Option<std::path::PathBuf> {
    let reference = std::env::var_os("FGUMI_BWA_MEM3_TEST_REF");
    if reference.is_none() {
        assert!(
            std::env::var_os("FGUMI_BWA_MEM3_REQUIRE_TOOLS").is_none(),
            "FGUMI_BWA_MEM3_REQUIRE_TOOLS is set but FGUMI_BWA_MEM3_TEST_REF (a bwa-mem3 \
             indexed reference) is not ({test})"
        );
        eprintln!("FGUMI_BWA_MEM3_TEST_REF not set: skipping {test}");
    }
    reference.map(std::path::PathBuf::from)
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
