//! Simulation commands for generating synthetic test data.
//!
//! This module provides commands to generate synthetic sequencing data
//! for benchmarking and testing the fgumi pipeline.

pub mod common;
pub mod consensus_reads;
pub mod correct_reads;
pub mod fastq_reads;
pub mod grouped_reads;
pub mod mapped_reads;
pub mod sort;

use crate::commands::command::Command;
use anyhow::Result;
use clap::{Parser, Subcommand};

pub use consensus_reads::ConsensusReads;
pub use correct_reads::CorrectReads;
pub use fastq_reads::FastqReads;
pub use grouped_reads::GroupedReads;
pub use mapped_reads::MappedReads;

/// Generate synthetic test data for benchmarking fgumi.
#[derive(Parser, Debug)]
#[command(
    name = "simulate",
    about = "\x1b[38;5;166m[UTILITIES]\x1b[0m      \x1b[36mGenerate synthetic test data\x1b[0m"
)]
pub struct Simulate {
    #[command(subcommand)]
    pub command: SimulateCommand,
}

impl Command for Simulate {
    fn execute(&self, command_line: &str) -> Result<()> {
        self.command.execute(command_line)
    }
}

#[derive(Subcommand, Debug)]
pub enum SimulateCommand {
    FastqReads(FastqReads),
    MappedReads(MappedReads),
    GroupedReads(GroupedReads),
    ConsensusReads(ConsensusReads),
    CorrectReads(CorrectReads),
}

impl SimulateCommand {
    fn execute(&self, command_line: &str) -> Result<()> {
        match self {
            Self::FastqReads(cmd) => cmd.execute(command_line),
            Self::MappedReads(cmd) => cmd.execute(command_line),
            Self::GroupedReads(cmd) => cmd.execute(command_line),
            Self::ConsensusReads(cmd) => cmd.execute(command_line),
            Self::CorrectReads(cmd) => cmd.execute(command_line),
        }
    }
}
