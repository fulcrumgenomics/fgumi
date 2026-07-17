//! Comparison commands for developer testing and validation.
//!
//! This module provides commands to compare BAM files and metrics files
//! for testing and validation of fgumi output against other tools.

pub mod bams;
pub(crate) mod engines;
pub mod metrics;
pub(crate) mod molecule;
pub(crate) mod raw_compare;
pub(crate) mod record_key;

use crate::commands::command::Command;
use anyhow::Result;
use clap::{Parser, Subcommand};

pub use bams::CompareBams;
pub use metrics::CompareMetrics;

// Narrow facade over the `engines` module tree: `engines` itself stays
// `pub(crate)` (it is an internal implementation detail of `compare`), but the
// separate-crate integration tests for the positional and key-join engines need a
// public entry point. Re-export exactly the items those tests use rather than
// widening `engines` to `pub`, so the rest of `engines`' internals stay out of
// the public API surface.
pub use engines::content::ContentPredicate;
pub use engines::header::compare_headers;
pub use engines::keyjoin::{
    KeyJoinConfig, KeyJoinOutcome, canonicalize_to_queryname, keyjoin_compare,
};
pub use engines::positional::positional_compare;
pub use engines::sort_verify::{SortVerifyOutcome, sort_verify_compare};

/// Marker error returned when a compare subcommand detects a mismatch between
/// the two inputs.
///
/// The CLI maps this to exit code 1 so the shell contract is preserved.
/// In-process callers (e.g. integration tests running under `cargo-llvm-cov`) can
/// downcast an `anyhow::Error` to this type to distinguish "mismatch found" from a
/// genuine I/O or logic error.
#[derive(Debug, thiserror::Error)]
#[error("compare detected a mismatch: {0}")]
pub struct CompareMismatch(pub String);

/// Compare files for testing and validation.
#[derive(Parser, Debug)]
#[command(
    name = "compare",
    about = "\x1b[38;5;166m[UTILITIES]\x1b[0m      \x1b[36mCompare files for testing and validation\x1b[0m"
)]
pub struct Compare {
    #[command(subcommand)]
    pub command: CompareCommand,
}

impl Command for Compare {
    fn execute(&self, command_line: &str) -> Result<()> {
        self.command.execute(command_line)
    }
}

#[derive(Subcommand, Debug)]
pub enum CompareCommand {
    Bams(CompareBams),
    Metrics(CompareMetrics),
}

impl CompareCommand {
    fn execute(&self, command_line: &str) -> Result<()> {
        match self {
            Self::Bams(cmd) => cmd.execute(command_line),
            Self::Metrics(cmd) => cmd.execute(command_line),
        }
    }
}
