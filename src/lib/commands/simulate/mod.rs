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

/// Compute the BAM bin for an alignment range using SAM spec §5.3 reg2bin.
///
/// `start_1based` and `end_1based` are 1-based inclusive coordinates. Returns the
/// SAM "unmapped bin" (4680) when either bound is `None`. Used by the simulate
/// commands when constructing mapped BAM records via `fgumi_raw_bam::SamBuilder`,
/// which does not compute the bin field automatically.
///
/// The bin layer constants `((1 << k) - 1) / 7` are written verbatim from the SAM
/// specification reference C code; clippy's `eq_op` lint flags the smallest layer
/// (`((1 << 3) - 1) / 7 == 7 / 7 == 1`) as a tautology, but we keep the form to
/// stay textually identical to the spec.
#[allow(clippy::eq_op, clippy::cast_possible_truncation)]
#[must_use]
pub fn region_to_bin(start_1based: Option<u32>, end_1based: Option<u32>) -> u16 {
    /// SAM spec §4.2.1: `reg2bin(-1, 0)` = 4680.
    const UNMAPPED_BIN: u16 = 4680;

    let (Some(start_1), Some(end_1)) = (start_1based, end_1based) else {
        return UNMAPPED_BIN;
    };
    let start = (start_1.saturating_sub(1)) as usize;
    let end = (end_1.saturating_sub(1)) as usize;

    let bin = if start >> 14 == end >> 14 {
        ((1 << 15) - 1) / 7 + (start >> 14)
    } else if start >> 17 == end >> 17 {
        ((1 << 12) - 1) / 7 + (start >> 17)
    } else if start >> 20 == end >> 20 {
        ((1 << 9) - 1) / 7 + (start >> 20)
    } else if start >> 23 == end >> 23 {
        ((1 << 6) - 1) / 7 + (start >> 23)
    } else if start >> 26 == end >> 26 {
        ((1 << 3) - 1) / 7 + (start >> 26)
    } else {
        0
    };

    // Truncate overflowing bin IDs (matches noodles behaviour).
    bin as u16
}
