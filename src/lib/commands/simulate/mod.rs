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
/// Delegates to the canonical [`fgumi_raw_bam::reg2bin`] so the simulate encoder
/// and the raw clip/zipper pipelines share one implementation of the binning
/// scheme. `reg2bin` takes a 0-based half-open `[beg, end)` interval, so the
/// 1-based inclusive `start`/`end` map to `start - 1` and `end` respectively.
#[must_use]
pub fn region_to_bin(start_1based: Option<u32>, end_1based: Option<u32>) -> u16 {
    let (Some(start_1), Some(end_1)) = (start_1based, end_1based) else {
        return fgumi_raw_bam::UNMAPPED_BIN;
    };
    let beg = i32::try_from(start_1.saturating_sub(1)).unwrap_or(i32::MAX);
    // `reg2bin` takes an exclusive end; clamp to 1 so a degenerate `end_1 == 0`
    // maps to the same bin the pre-delegation code produced (which saturated the
    // 0-based inclusive end to 0, i.e. an exclusive end of 1).
    let end = i32::try_from(end_1.max(1)).unwrap_or(i32::MAX);
    fgumi_raw_bam::reg2bin(beg, end)
}

#[cfg(test)]
mod tests {
    use super::region_to_bin;
    use rstest::rstest;

    #[rstest]
    // 1-based inclusive coordinates -> BAM bin (SAM spec §5.3 reg2bin).
    #[case::single_base(Some(1), Some(1), 4681)]
    #[case::within_first_16kb(Some(101), Some(300), 4681)]
    #[case::second_16kb_window(Some(16417), Some(16500), 4682)]
    #[case::straddles_16kb_boundary(Some(16301), Some(16500), 585)]
    #[case::unmapped_start(None, Some(100), 4680)]
    #[case::unmapped_end(Some(100), None, 4680)]
    #[case::unmapped_both(None, None, 4680)]
    // Degenerate 1-based end of 0 clamps to the first bin, matching pre-refactor behavior.
    #[case::degenerate_zero_end(Some(1), Some(0), 4681)]
    fn region_to_bin_matches_sam_spec(
        #[case] start_1based: Option<u32>,
        #[case] end_1based: Option<u32>,
        #[case] expected: u16,
    ) {
        assert_eq!(region_to_bin(start_1based, end_1based), expected);
    }
}
