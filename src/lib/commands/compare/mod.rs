//! Comparison commands for developer testing and validation.
//!
//! This module provides commands to compare BAM files and metrics files
//! for testing and validation of fgumi output against other tools.

pub mod bam_roundtrip;
pub mod bams;
pub mod metrics;
pub(crate) mod raw_compare;

use crate::commands::command::Command;
use anyhow::Result;
use clap::{Parser, Subcommand};

pub use bam_roundtrip::CompareBamRoundtrip;
pub use bams::CompareBams;
pub use metrics::CompareMetrics;

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
    BamRoundtrip(CompareBamRoundtrip),
    Metrics(CompareMetrics),
}

impl CompareCommand {
    fn execute(&self, command_line: &str) -> Result<()> {
        match self {
            Self::Bams(cmd) => cmd.execute(command_line),
            Self::BamRoundtrip(cmd) => cmd.execute(command_line),
            Self::Metrics(cmd) => cmd.execute(command_line),
        }
    }
}
