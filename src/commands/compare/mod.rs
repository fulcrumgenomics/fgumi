//! Comparison commands for developer testing and validation.
//!
//! This module provides commands to compare BAM files and metrics files
//! for testing and validation of fgumi output against other tools.

pub mod bams;
pub mod metrics;

use crate::commands::command::Command;
use anyhow::Result;
use clap::{Parser, Subcommand};

pub use bams::CompareBams;
pub use metrics::CompareMetrics;

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
