//! Enhanced logging utilities for formatted output.
//!
//! This module provides consistent, user-friendly logging utilities for metrics,
//! progress tracking, and operation summaries.

use std::time::{Duration, Instant};

use crate::metrics::{ConsensusMetrics, UmiGroupingMetrics};
use crate::rejection::format_count;

/// Formats a percentage with specified decimal places.
///
/// # Arguments
///
/// * `value` - The fraction (0.0-1.0) to format as percentage
/// * `decimals` - Number of decimal places to include
///
/// # Returns
///
/// A string formatted as "XX.XX%" (e.g., "95.43%")
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

/// Formats a duration in human-readable form.
///
/// # Arguments
///
/// * `duration` - The duration to format
///
/// # Returns
///
/// A human-readable string (e.g., "2m 15s", "1h 30m", "45s")
///
/// # Examples
///
/// ```
/// use fgumi_lib::logging::format_duration;
/// use std::time::Duration;
///
/// assert_eq!(format_duration(Duration::from_secs(45)), "45s");
/// assert_eq!(format_duration(Duration::from_secs(135)), "2m 15s");
/// assert_eq!(format_duration(Duration::from_secs(5400)), "1h 30m");
/// ```
#[must_use]
pub fn format_duration(duration: Duration) -> String {
    let secs = duration.as_secs();
    if secs < 60 {
        format!("{secs}s")
    } else if secs < 3600 {
        let mins = secs / 60;
        let remaining_secs = secs % 60;
        if remaining_secs == 0 { format!("{mins}m") } else { format!("{mins}m {remaining_secs}s") }
    } else {
        let hours = secs / 3600;
        let mins = (secs % 3600) / 60;
        if mins == 0 { format!("{hours}h") } else { format!("{hours}h {mins}m") }
    }
}

/// Formats a rate (items per second) with appropriate units.
///
/// # Arguments
///
/// * `count` - Number of items processed
/// * `duration` - Time taken to process items
///
/// # Returns
///
/// A formatted rate string (e.g., "1,234 reads/s", "50 reads/min")
///
/// # Examples
///
/// ```
/// use fgumi_lib::logging::format_rate;
/// use std::time::Duration;
///
/// assert_eq!(format_rate(1000, Duration::from_secs(1)), "1,000 items/s");
/// assert_eq!(format_rate(600, Duration::from_secs(60)), "10 items/s");
/// ```
#[must_use]
#[allow(clippy::cast_precision_loss, clippy::cast_possible_truncation, clippy::cast_sign_loss)]
pub fn format_rate(count: u64, duration: Duration) -> String {
    let secs = duration.as_secs_f64();
    if secs < 0.001 {
        return format!("{} items/s", format_count(count));
    }

    let rate = count as f64 / secs;
    if rate >= 1.0 {
        format!("{} items/s", format_count(rate as u64))
    } else {
        let items_per_min = count as f64 / (secs / 60.0);
        format!("{items_per_min:.1} items/min")
    }
}

/// Logs a formatted summary of consensus metrics.
///
/// Outputs key metrics including input/output counts, rejection breakdown,
/// and quality statistics.
///
/// # Arguments
///
/// * `metrics` - The consensus metrics to summarize
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
/// # Arguments
///
/// * `metrics` - The UMI grouping metrics to summarize
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
    // Non-PF only logged if > 0
    if metrics.discarded_non_pf > 0 {
        log::info!("Filtered out {} non-PF records.", format_count(metrics.discarded_non_pf));
    }
    // Poor alignment always logged (like fgbio)
    log::info!(
        "Filtered out {} records due to mapping issues.",
        format_count(metrics.discarded_poor_alignment)
    );
    // Ns in UMI always logged (like fgbio)
    log::info!(
        "Filtered out {} records that contained one or more Ns in their UMIs.",
        format_count(metrics.discarded_ns_in_umi)
    );
    // UMI too short only logged if > 0
    if metrics.discarded_umi_too_short > 0 {
        log::info!(
            "Filtered out {} records that contained UMIs that were too short.",
            format_count(metrics.discarded_umi_too_short)
        );
    }
}

/// Operation timing and summary helper.
///
/// Tracks operation timing and provides formatted summary output.
///
/// # Examples
///
/// ```no_run
/// use fgumi_lib::logging::OperationTimer;
///
/// let timer = OperationTimer::new("Processing reads");
///
/// // ... do work ...
///
/// timer.log_completion(10_000); // Log with item count
/// ```
pub struct OperationTimer {
    operation: String,
    start_time: Instant,
}

impl OperationTimer {
    /// Creates a new operation timer and logs the start.
    #[must_use]
    pub fn new(operation: &str) -> Self {
        log::info!("{operation} ...");
        Self { operation: operation.to_string(), start_time: Instant::now() }
    }

    /// Logs the completion with item count and rate.
    pub fn log_completion(&self, count: u64) {
        let duration = self.start_time.elapsed();
        log::info!(
            "{} completed: {} in {} ({})",
            self.operation,
            format_count(count),
            format_duration(duration),
            format_rate(count, duration)
        );
    }
}

#[cfg(test)]
#[allow(clippy::cast_precision_loss)]
mod tests {
    use super::*;

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
