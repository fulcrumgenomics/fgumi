#![deny(unsafe_code)]

//! Structured metric types and TSV writer for fgumi operations.
//!
//! This crate provides:
//! - [`Metric`] and [`ProcessingMetrics`] traits for extensible metric types
//! - Metric structs for consensus, grouping, correction, and duplex operations
//! - [`rejection`] module for rejection reason tracking
//! - [`writer`] module for TSV file output

#[cfg(feature = "clip")]
pub mod clip;
pub mod consensus;
pub mod correct;
pub mod downsampling;
pub mod duplex;
pub mod group;
pub mod inline_collector;
pub mod intervals;
pub mod rejection;
pub mod shared;
pub mod simplex;
pub mod template_info;
pub mod writer;

use serde::{Deserialize, Serialize};

/// Number of decimal places used for float metrics (matches fgbio).
pub const FLOAT_PRECISION: usize = 6;

/// Formats a float value with the standard precision for metrics.
///
/// This ensures consistent float formatting across all metrics output,
/// matching fgbio's 6 decimal place precision.
///
/// # Example
/// ```
/// use fgumi_metrics::format_float;
/// assert_eq!(format_float(0.9), "0.900000");
/// assert_eq!(format_float(0.0), "0.000000");
/// ```
#[must_use]
pub fn format_float(value: f64) -> String {
    format!("{value:.FLOAT_PRECISION$}")
}

/// Computes `numerator / denominator`, returning 0.0 if the denominator is zero.
#[must_use]
#[expect(clippy::cast_precision_loss, reason = "metric counts never exceed 2^53")]
pub fn frac(numerator: usize, denominator: usize) -> f64 {
    if denominator > 0 { numerator as f64 / denominator as f64 } else { 0.0 }
}

/// Computes `numerator / denominator` for `u64` values, returning 0.0 if denominator is zero.
#[must_use]
#[expect(clippy::cast_precision_loss, reason = "metric counts never exceed 2^53")]
pub fn frac_u64(numerator: u64, denominator: u64) -> f64 {
    if denominator > 0 { numerator as f64 / denominator as f64 } else { 0.0 }
}

/// A metric type that can be serialized to TSV files.
///
/// All metric types in fgumi implement this trait, providing a consistent
/// interface for serialization and identification.
pub trait Metric: Serialize + for<'de> Deserialize<'de> + Clone + Default {
    /// Human-readable name for this metric type.
    ///
    /// Used in error messages and logging when writing metrics files.
    fn metric_name() -> &'static str;
}

/// Common interface for metrics that track processing pipeline counts.
///
/// This trait provides a consistent way to access input, output, and filtered
/// counts across different metric types, enabling generic summary output.
pub trait ProcessingMetrics {
    /// Total number of input items (reads, records, etc.) processed.
    fn total_input(&self) -> u64;

    /// Total number of output items (consensus reads, accepted records, etc.) produced.
    fn total_output(&self) -> u64;

    /// Total number of items filtered out or rejected.
    fn total_filtered(&self) -> u64;

    /// Processing efficiency as a percentage (output / input * 100).
    fn efficiency(&self) -> f64 {
        frac_u64(self.total_output(), self.total_input()) * 100.0
    }
}

// Re-export commonly used types
#[cfg(feature = "clip")]
pub use clip::{ClipCounts, ClippingMetrics, ClippingMetricsCollection, ReadType};
pub use consensus::{ConsensusKvMetric, ConsensusMetrics};
pub use correct::UmiCorrectionMetrics;
pub use downsampling::{DOWNSAMPLING_FRACTIONS, compute_hash_fraction};
pub use duplex::{
    DuplexFamilySizeMetric, DuplexMetricsCollector, DuplexUmiMetric, DuplexYieldMetric,
    FamilySizeMetric,
};
pub use group::{FamilySizeMetrics, PositionGroupSizeMetrics, UmiGroupingMetrics};
pub use inline_collector::{InlineCollector, InlineMetricsCollector, MetricsCollector};
pub use intervals::{Interval, overlaps_intervals, parse_intervals};
pub use rejection::{RejectionReason, format_count};
pub use shared::UmiMetric;
pub use simplex::{SimplexFamilySizeMetric, SimplexMetricsCollector, SimplexYieldMetric};
pub use template_info::{ReadInfoKey, TemplateInfo, TemplateMetadata, compute_template_metadata};
pub use writer::{read_metrics, read_metrics_auto, write_metrics};

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_frac_normal() {
        assert!((frac(3, 4) - 0.75).abs() < f64::EPSILON);
    }

    #[test]
    fn test_frac_zero_denominator() {
        assert!((frac(5, 0)).abs() < f64::EPSILON);
    }

    #[test]
    fn test_frac_zero_numerator() {
        assert!((frac(0, 10)).abs() < f64::EPSILON);
    }

    #[test]
    fn test_frac_u64_normal() {
        assert!((frac_u64(3, 4) - 0.75).abs() < f64::EPSILON);
    }

    #[test]
    fn test_frac_u64_zero_denominator() {
        assert!((frac_u64(5, 0)).abs() < f64::EPSILON);
    }

    #[test]
    fn test_frac_u64_zero_numerator() {
        assert!((frac_u64(0, 10)).abs() < f64::EPSILON);
    }

    #[test]
    fn test_processing_metrics_consensus() {
        let metrics = ConsensusMetrics {
            total_input_reads: 1000,
            consensus_reads: 800,
            filtered_reads: 200,
            ..Default::default()
        };

        assert_eq!(metrics.total_input(), 1000);
        assert_eq!(metrics.total_output(), 800);
        assert_eq!(metrics.total_filtered(), 200);
        assert!((metrics.efficiency() - 80.0).abs() < f64::EPSILON);
    }

    #[test]
    fn test_processing_metrics_grouping() {
        let metrics = UmiGroupingMetrics {
            total_records: 1000,
            accepted_records: 900,
            discarded_non_pf: 50,
            discarded_poor_alignment: 30,
            discarded_ns_in_umi: 20,
            ..Default::default()
        };

        assert_eq!(metrics.total_input(), 1000);
        assert_eq!(metrics.total_output(), 900);
        assert_eq!(metrics.total_filtered(), 100);
        assert!((metrics.efficiency() - 90.0).abs() < f64::EPSILON);
    }

    #[test]
    fn test_processing_metrics_zero_input() {
        let metrics = ConsensusMetrics::default();

        assert_eq!(metrics.total_input(), 0);
        assert_eq!(metrics.total_output(), 0);
        assert_eq!(metrics.total_filtered(), 0);
        assert!((metrics.efficiency()).abs() < f64::EPSILON);
    }

    #[test]
    fn test_processing_metrics_generic_usage() {
        fn log_efficiency(m: &impl ProcessingMetrics) -> f64 {
            m.efficiency()
        }

        let consensus =
            ConsensusMetrics { total_input_reads: 100, consensus_reads: 50, ..Default::default() };

        let grouping =
            UmiGroupingMetrics { total_records: 100, accepted_records: 75, ..Default::default() };

        assert!((log_efficiency(&consensus) - 50.0).abs() < f64::EPSILON);
        assert!((log_efficiency(&grouping) - 75.0).abs() < f64::EPSILON);
    }
}
