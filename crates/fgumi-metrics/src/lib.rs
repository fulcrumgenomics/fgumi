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
pub mod duplex;
pub mod group;
pub mod rejection;
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
        if self.total_input() == 0 {
            0.0
        } else {
            #[expect(clippy::cast_precision_loss, reason = "read counts never exceed 2^53")]
            let result = self.total_output() as f64 / self.total_input() as f64 * 100.0;
            result
        }
    }
}

// Re-export commonly used types
#[cfg(feature = "clip")]
pub use clip::{ClippingMetrics, ClippingMetricsCollection, ReadType};
pub use consensus::{ConsensusKvMetric, ConsensusMetrics};
pub use correct::UmiCorrectionMetrics;
pub use duplex::{
    DuplexFamilySizeMetric, DuplexMetricsCollector, DuplexUmiMetric, DuplexYieldMetric,
    FamilySizeMetric, UmiMetric,
};
pub use group::{FamilySizeMetrics, UmiGroupingMetrics};
pub use rejection::{RejectionReason, format_count};
pub use writer::write_metrics;

#[cfg(test)]
mod tests {
    use super::*;

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
