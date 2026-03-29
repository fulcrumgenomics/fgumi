//! Shared helpers for locating reference data and external tools in integration tests.

use std::path::{Path, PathBuf};

/// Path to the tracked `PhiX` reference FASTA, relative to the repo root.
const TRACKED_PHIX_FASTA: &str = "tests/data/references/phix.fasta";

/// Check if `bwa` is on the PATH.
pub fn bwa_available() -> bool {
    which::which("bwa").is_ok()
}

/// Locate the `PhiX` reference FASTA for alignment tests.
///
/// Checks in order:
/// 1. The `FGUMI_TEST_REFERENCE` environment variable (if set and file exists).
/// 2. The tracked copy at `tests/data/references/phix.fasta`.
///
/// Returns `None` if neither is available.
pub fn phix_reference() -> Option<PathBuf> {
    // Environment variable override
    if let Some(path) =
        std::env::var("FGUMI_TEST_REFERENCE").ok().map(PathBuf::from).filter(|p| p.exists())
    {
        return Some(path);
    }

    // Tracked copy — resolve from CARGO_MANIFEST_DIR (repo root)
    let manifest_dir = std::env::var("CARGO_MANIFEST_DIR").ok()?;
    let tracked = PathBuf::from(manifest_dir).join(TRACKED_PHIX_FASTA);
    tracked.exists().then_some(tracked)
}

/// Check that a bwa index exists for the given reference.
pub fn has_bwa_index(ref_path: &Path) -> bool {
    let bwt = PathBuf::from(format!("{}.bwt", ref_path.display()));
    bwt.exists()
}

/// Whether we are running in CI (GitHub Actions sets `CI=true`).
fn in_ci() -> bool {
    std::env::var("CI").is_ok()
}

/// Skip locally, fail in CI.
fn skip_or_fail(msg: &str) -> Option<PathBuf> {
    assert!(!in_ci(), "CI: {msg}");
    eprintln!("SKIP: {msg}");
    None
}

/// Check all prerequisites for alignment tests. Returns the reference path,
/// or `None` (skip) when running locally. In CI, panics on missing prerequisites
/// so the test fails loudly rather than being silently skipped.
pub fn check_alignment_prerequisites() -> Option<PathBuf> {
    if !bwa_available() {
        return skip_or_fail("bwa not found on PATH — install bwa");
    }
    let Some(ref_path) = phix_reference() else {
        return skip_or_fail(
            "PhiX reference not found. Set FGUMI_TEST_REFERENCE env var \
             or ensure tests/data/references/phix.fasta exists.",
        );
    };
    if !has_bwa_index(&ref_path) {
        return skip_or_fail(&format!(
            "bwa index not found for {}. Run `bwa index` first.",
            ref_path.display()
        ));
    }
    Some(ref_path)
}
