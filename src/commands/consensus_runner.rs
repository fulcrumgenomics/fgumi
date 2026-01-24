//! Common infrastructure for consensus calling commands.
//!
//! This module provides shared traits and utilities used by simplex, duplex, and codec
//! consensus calling commands to reduce code duplication.

use anyhow::Result;
use bstr::BString;
use fgumi_lib::bam_io::create_bam_reader;
use fgumi_lib::consensus::codec_caller::CodecConsensusStats;
use fgumi_lib::consensus_caller::{ConsensusCallingStats, make_prefix_from_header};
use fgumi_lib::logging::OperationTimer;
use fgumi_lib::metrics::consensus::ConsensusMetrics;
use fgumi_lib::overlapping_consensus::CorrectionStats;
use log::info;
use noodles::sam::Header;
use noodles::sam::alignment::record_buf::RecordBuf;
use noodles::sam::header::record::value::Map;
use noodles::sam::header::record::value::map::ReadGroup;
use noodles::sam::header::record::value::map::header::tag as header_tag;
use std::path::PathBuf;

use crate::commands::common::ThreadingOptions;

/// Trait for converting command-specific statistics to metrics.
pub trait ConsensusStatsOps: Clone + Default + Send {
    /// Merges another stats instance into this one.
    fn merge(&mut self, other: &Self);

    /// Converts to `ConsensusMetrics` for TSV output.
    fn to_metrics(&self) -> ConsensusMetrics;
}

impl ConsensusStatsOps for ConsensusCallingStats {
    fn merge(&mut self, other: &Self) {
        ConsensusCallingStats::merge(self, other);
    }

    fn to_metrics(&self) -> ConsensusMetrics {
        let mut metrics = ConsensusMetrics::new();
        metrics.total_input_reads = self.total_reads as u64;
        metrics.consensus_reads = self.consensus_reads as u64;
        metrics.filtered_reads = self.filtered_reads as u64;

        // Convert rejection reasons to centralized reasons
        for (reason, count) in &self.rejection_reasons {
            metrics.add_rejection(reason.to_centralized(), *count as u64);
        }

        metrics
    }
}

impl ConsensusStatsOps for CodecConsensusStats {
    fn merge(&mut self, other: &Self) {
        self.total_input_reads += other.total_input_reads;
        self.consensus_reads_generated += other.consensus_reads_generated;
        self.reads_filtered += other.reads_filtered;
        self.consensus_reads_rejected_hdd += other.consensus_reads_rejected_hdd;
        self.consensus_bases_emitted += other.consensus_bases_emitted;
        self.consensus_duplex_bases_emitted += other.consensus_duplex_bases_emitted;
        self.duplex_disagreement_base_count += other.duplex_disagreement_base_count;

        // Merge rejection reasons
        for (reason, count) in &other.rejection_reasons {
            *self.rejection_reasons.entry(*reason).or_insert(0) += count;
        }
    }

    fn to_metrics(&self) -> ConsensusMetrics {
        let mut metrics = ConsensusMetrics::new();
        metrics.total_input_reads = self.total_input_reads;
        metrics.consensus_reads = self.consensus_reads_generated;
        metrics.filtered_reads = self.reads_filtered;

        // Convert rejection reasons to centralized reasons
        for (reason, count) in &self.rejection_reasons {
            metrics.add_rejection(reason.to_centralized(), *count as u64);
        }

        metrics
    }
}

/// Creates an output header for unmapped consensus reads.
///
/// This creates a header with:
/// - A single read group with the specified ID
/// - Sort order set to "unknown" and group order to "query"
/// - A comment indicating the number of input read groups
/// - A @PG record with version and command line
pub fn create_unmapped_consensus_header(
    input_header: &Header,
    read_group_id: &str,
    comment_prefix: &str,
    version: &str,
    command_line: &str,
) -> Result<Header> {
    let mut output_header = Header::builder();

    // Create read group
    let new_rg = Map::<ReadGroup>::builder().build()?;
    output_header = output_header.add_read_group(BString::from(read_group_id), new_rg);

    // Set sort order
    let header_map = Map::<noodles::sam::header::record::value::map::Header>::builder()
        .insert(header_tag::SORT_ORDER, BString::from("unknown"))
        .insert(header_tag::GROUP_ORDER, BString::from("query"))
        .build()?;
    output_header = output_header.set_header(header_map);

    // Add comment
    let rg_count = input_header.read_groups().len();
    output_header = output_header.add_comment(format!(
        "{comment_prefix} {read_group_id} contains consensus reads generated from {rg_count} input read groups."
    ));

    // Add @PG record
    output_header = fgumi_lib::header::add_pg_to_builder(output_header, version, command_line)?;

    Ok(output_header.build())
}

/// Logs overlapping consensus statistics if enabled.
pub fn log_overlapping_stats(stats: &CorrectionStats) {
    info!("Overlapping consensus statistics:");
    info!("  Overlapping bases examined: {}", stats.overlapping_bases);
    info!("  Bases agreeing: {}", stats.bases_agreeing);
    info!("  Bases disagreeing: {}", stats.bases_disagreeing);
    info!("  Bases corrected: {}", stats.bases_corrected);
}

/// Type alias for consensus result tuple: (UMI, consensus reads result, stats, rejects, overlapping stats)
pub type ConsensusResultTuple<S> =
    (String, Result<Vec<RecordBuf>>, S, Vec<RecordBuf>, Option<CorrectionStats>);

/// Configuration for consensus command execution.
#[derive(Clone)]
pub struct ConsensusConfig {
    pub input: PathBuf,
    pub output: PathBuf,
    pub rejects: Option<PathBuf>,
    pub stats_path: Option<PathBuf>,
    pub threading: ThreadingOptions,
    pub overlapping_enabled: bool,
    pub tag: String,
    pub operation_name: String,
}

impl ConsensusConfig {
    /// Returns the number of threads from threading options.
    pub fn num_threads(&self) -> usize {
        self.threading.num_threads()
    }

    /// Returns whether parallel processing is enabled.
    pub fn is_parallel(&self) -> bool {
        self.threading.is_parallel()
    }

    /// Returns whether rejects tracking is enabled.
    pub fn track_rejects(&self) -> bool {
        self.rejects.is_some()
    }
}

/// Helper for executing consensus commands with common setup/teardown.
pub struct ConsensusExecutionContext {
    pub config: ConsensusConfig,
    pub timer: OperationTimer,
    pub header: Header,
    pub output_header: Header,
    pub read_name_prefix: String,
}

impl ConsensusExecutionContext {
    /// Creates a new execution context.
    ///
    /// Opens the input BAM to read the header, creates the output header,
    /// and sets up the read name prefix.
    pub fn new(
        config: ConsensusConfig,
        output_header: Header,
        read_name_prefix_override: Option<String>,
    ) -> Result<Self> {
        let timer = OperationTimer::new(&config.operation_name);

        // Open input to get header
        let (_reader, header) = create_bam_reader(&config.input, 1)?;

        // Determine read name prefix
        let read_name_prefix =
            read_name_prefix_override.unwrap_or_else(|| make_prefix_from_header(&header));

        Ok(Self { config, timer, header, output_header, read_name_prefix })
    }

    /// Creates a new execution context with an unmapped consensus header.
    pub fn new_with_unmapped_header(
        config: ConsensusConfig,
        read_group_id: &str,
        comment_prefix: &str,
        version: &str,
        command_line: &str,
        read_name_prefix_override: Option<String>,
    ) -> Result<Self> {
        let timer = OperationTimer::new(&config.operation_name);

        // Open input to get header
        let (_reader, header) = create_bam_reader(&config.input, 1)?;

        // Create output header for unmapped consensus reads
        let output_header = create_unmapped_consensus_header(
            &header,
            read_group_id,
            comment_prefix,
            version,
            command_line,
        )?;

        // Determine read name prefix
        let read_name_prefix =
            read_name_prefix_override.unwrap_or_else(|| make_prefix_from_header(&header));

        Ok(Self { config, timer, header, output_header, read_name_prefix })
    }

    /// Creates a new execution context with a passthrough header.
    pub fn new_with_passthrough_header(
        config: ConsensusConfig,
        read_name_prefix_override: Option<String>,
    ) -> Result<Self> {
        let timer = OperationTimer::new(&config.operation_name);

        // Open input to get header
        let (_reader, header) = create_bam_reader(&config.input, 1)?;

        // Use input header as output header (passthrough)
        let output_header = header.clone();

        // Determine read name prefix
        let read_name_prefix =
            read_name_prefix_override.unwrap_or_else(|| make_prefix_from_header(&header));

        Ok(Self { config, timer, header, output_header, read_name_prefix })
    }

    /// Logs completion with consensus count.
    pub fn log_completion(&self, consensus_count: u64) {
        self.timer.log_completion(consensus_count);
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_stats_to_metrics_empty() {
        let stats = ConsensusCallingStats::new();
        let metrics = stats.to_metrics();
        assert_eq!(metrics.total_input_reads, 0);
        assert_eq!(metrics.consensus_reads, 0);
        assert_eq!(metrics.filtered_reads, 0);
    }

    #[test]
    fn test_stats_to_metrics_with_values() {
        let mut stats = ConsensusCallingStats::new();
        stats.total_reads = 100;
        stats.consensus_reads = 50;
        stats.filtered_reads = 10;

        let metrics = stats.to_metrics();
        assert_eq!(metrics.total_input_reads, 100);
        assert_eq!(metrics.consensus_reads, 50);
        assert_eq!(metrics.filtered_reads, 10);
    }

    #[test]
    fn test_stats_merge() {
        let mut stats1 = ConsensusCallingStats::new();
        stats1.total_reads = 100;
        stats1.consensus_reads = 50;

        let mut stats2 = ConsensusCallingStats::new();
        stats2.total_reads = 200;
        stats2.consensus_reads = 75;

        stats1.merge(&stats2);

        assert_eq!(stats1.total_reads, 300);
        assert_eq!(stats1.consensus_reads, 125);
    }
}
