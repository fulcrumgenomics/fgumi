//! Enhanced logging utilities for formatted output.
//!
//! This module provides consistent, user-friendly logging utilities for metrics,
//! progress tracking, and operation summaries.
//!
//! `format_count`, `format_duration`, `format_rate`, and `OperationTimer` are
//! re-exported from `fgumi-cli-common` (single source of truth).  `format_percent`
//! and the consensus/UMI summary helpers are umbrella-only and defined here.

pub use fgumi_cli_common::{OperationTimer, format_count, format_duration, format_rate};

use crate::metrics::{ConsensusMetrics, UmiGroupingMetrics};

/// Formats a percentage with specified decimal places.
///
/// # Examples
///
/// ```
/// use fgumi_lib::logging::format_percent;
///
/// assert_eq!(format_percent(0.9543, 2), "95.43%");
/// assert_eq!(format_percent(0.5, 1), "50.0%");
/// assert_eq!(format_percent(1.0, 0), "100%");
/// ```
#[must_use]
pub fn format_percent(value: f64, decimals: usize) -> String {
    format!("{:.decimals$}%", value * 100.0, decimals = decimals)
}

/// Logs a formatted summary of consensus metrics.
///
/// Outputs key metrics including input/output counts, rejection breakdown,
/// and quality statistics.
///
/// # Examples
///
/// ```no_run
/// use fgumi_lib::logging::log_consensus_summary;
/// use fgumi_lib::metrics::ConsensusMetrics;
///
/// let mut metrics = ConsensusMetrics::new();
/// metrics.total_input_reads = 10_000;
/// metrics.consensus_reads = 8_000;
/// metrics.filtered_reads = 2_000;
///
/// log_consensus_summary(&metrics);
/// ```
#[allow(clippy::cast_precision_loss)]
pub fn log_consensus_summary(metrics: &ConsensusMetrics) {
    log::info!("Consensus Calling Summary:");
    log::info!("  Input reads: {}", format_count(metrics.total_input_reads));
    log::info!("  Consensus reads: {}", format_count(metrics.consensus_reads));
    log::info!("  Filtered reads: {}", format_count(metrics.filtered_reads));

    if metrics.total_input_reads > 0 {
        let pass_rate = metrics.consensus_reads as f64 / metrics.total_input_reads as f64;
        log::info!("  Pass rate: {}", format_percent(pass_rate, 2));
    }

    if metrics.consensus_reads > 0 {
        log::info!("  Avg reads per consensus: {:.1}", metrics.avg_input_reads_per_consensus);
        log::info!("  Avg depth: {:.1}", metrics.avg_raw_read_depth);
    }

    if metrics.total_rejections() > 0 {
        log::info!("  Top rejection reasons:");
        let summary = metrics.rejection_summary();
        let mut sorted: Vec<_> = summary.iter().collect();
        sorted.sort_by(|a, b| b.1.cmp(a.1));

        for (reason, count) in sorted.iter().take(5) {
            log::info!("    {}: {}", reason.description(), format_count(**count));
        }
    }
}

/// Logs a formatted summary of UMI grouping metrics.
///
/// Outputs grouping statistics including molecule counts and family sizes.
///
/// # Examples
///
/// ```no_run
/// use fgumi_lib::logging::log_umi_grouping_summary;
/// use fgumi_lib::metrics::UmiGroupingMetrics;
///
/// let mut metrics = UmiGroupingMetrics::new();
/// metrics.total_records = 10_000;
/// metrics.accepted_records = 9_500;
/// metrics.unique_molecule_ids = 1_000;
///
/// log_umi_grouping_summary(&metrics);
/// ```
#[allow(clippy::cast_precision_loss)]
pub fn log_umi_grouping_summary(metrics: &UmiGroupingMetrics) {
    log::info!("UMI Grouping Summary:");
    log::info!("  Total records: {}", format_count(metrics.total_records));
    log::info!("  Accepted records: {}", format_count(metrics.accepted_records));

    if metrics.total_records > 0 {
        let accept_rate = metrics.accepted_records as f64 / metrics.total_records as f64;
        log::info!("  Accept rate: {}", format_percent(accept_rate, 2));
    }

    log::info!("  Unique molecules: {}", format_count(metrics.unique_molecule_ids));
    log::info!("  Total families: {}", format_count(metrics.total_families));

    if metrics.total_families > 0 {
        log::info!("  Avg reads/molecule: {:.1}", metrics.avg_reads_per_molecule);
        log::info!("  Median reads/molecule: {}", metrics.median_reads_per_molecule);
        log::info!("  Min reads/molecule: {}", metrics.min_reads_per_molecule);
        log::info!("  Max reads/molecule: {}", metrics.max_reads_per_molecule);
    }

    // Log filter counts (matching fgbio's logging style)
    if metrics.discarded_non_pf > 0 {
        log::info!("Filtered out {} non-PF records.", format_count(metrics.discarded_non_pf));
    }
    log::info!(
        "Filtered out {} records due to mapping issues.",
        format_count(metrics.discarded_poor_alignment)
    );
    log::info!(
        "Filtered out {} records that contained one or more Ns in their UMIs.",
        format_count(metrics.discarded_ns_in_umi)
    );
    if metrics.discarded_umi_too_short > 0 {
        log::info!(
            "Filtered out {} records that contained UMIs that were too short.",
            format_count(metrics.discarded_umi_too_short)
        );
    }
}

#[cfg(test)]
#[allow(clippy::cast_precision_loss)]
mod tests {
    use super::*;
    use std::time::Duration;

    #[test]
    fn test_format_percent() {
        assert_eq!(format_percent(0.9543, 2), "95.43%");
        assert_eq!(format_percent(0.5, 1), "50.0%");
        assert_eq!(format_percent(1.0, 0), "100%");
        assert_eq!(format_percent(0.0, 2), "0.00%");
    }

    #[test]
    fn test_format_duration() {
        assert_eq!(format_duration(Duration::from_secs(0)), "0s");
        assert_eq!(format_duration(Duration::from_secs(45)), "45s");
        assert_eq!(format_duration(Duration::from_secs(60)), "1m");
        assert_eq!(format_duration(Duration::from_secs(135)), "2m 15s");
        assert_eq!(format_duration(Duration::from_secs(3600)), "1h");
        assert_eq!(format_duration(Duration::from_secs(5400)), "1h 30m");
    }

    #[test]
    fn test_format_rate() {
        assert_eq!(format_rate(1000, Duration::from_secs(1)), "1,000 items/s");
        assert_eq!(format_rate(60, Duration::from_secs(60)), "1 items/s");
        assert_eq!(format_rate(30, Duration::from_secs(60)), "30.0 items/min");
        // Near-zero duration
        assert!(format_rate(1000, Duration::from_nanos(1)).contains("items/s"));
    }

    #[test]
    fn test_operation_timer() {
        let timer = OperationTimer::new("Test");
        timer.log_completion(1000);
    }

    #[test]
    fn test_log_consensus_summary() {
        // Empty metrics
        log_consensus_summary(&ConsensusMetrics::new());

        // With data
        let mut metrics = ConsensusMetrics::new();
        metrics.total_input_reads = 10000;
        metrics.consensus_reads = 8000;
        metrics.avg_input_reads_per_consensus = 1.25;
        metrics.rejected_insufficient_support = 1000;
        log_consensus_summary(&metrics);
    }

    #[test]
    fn test_log_umi_grouping_summary() {
        // Empty metrics
        log_umi_grouping_summary(&UmiGroupingMetrics::new());

        // With data and discards
        let mut metrics = UmiGroupingMetrics::new();
        metrics.total_records = 10000;
        metrics.accepted_records = 8000;
        metrics.total_families = 950;
        metrics.discarded_non_pf = 500;
        metrics.discarded_ns_in_umi = 200;
        log_umi_grouping_summary(&metrics);
    }
}
