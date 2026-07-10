//! Helper utilities for integration tests.
//!
//! Individual helpers are used by different subsets of the test modules
//! (some of which are feature-gated), so some helpers appear unused under
//! certain feature combinations. Silence the resulting warnings wholesale
//! so `--no-default-features` stays clean.

#![allow(dead_code, unused_imports)]

pub mod assertions;
pub mod bam_generator;

pub use assertions::*;
pub use bam_generator::*;

/// Writes `contents` to `dir/name` and returns the path.
///
/// Shared by the `compare` integration tests (`test_compare_metrics_command.rs`,
/// `test_compare_mutation.rs`) for building ad hoc TSV fixtures.
pub fn write_tsv(dir: &std::path::Path, name: &str, contents: &str) -> std::path::PathBuf {
    let path = dir.join(name);
    std::fs::write(&path, contents).expect("failed to write temp TSV");
    path
}
