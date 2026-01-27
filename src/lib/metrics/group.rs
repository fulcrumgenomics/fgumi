//! Metrics for the `group` command.
//!
//! This module provides metrics types for tracking UMI grouping operations
//! and family size distributions.

use serde::{Deserialize, Serialize};

use super::Metric;

/// Metrics for UMI grouping operations.
///
/// These metrics track how reads are grouped by UMI and provide insight into
/// data quality and molecule representation.
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct UmiGroupingMetrics {
    /// Total SAM records processed
    pub total_records: u64,

    /// Records accepted for grouping
    pub accepted_records: u64,

    /// Records discarded (not passing filter)
    pub discarded_non_pf: u64,

    /// Records discarded (poor alignment quality)
    pub discarded_poor_alignment: u64,

    /// Records discarded (Ns in UMI)
    pub discarded_ns_in_umi: u64,

    /// Records discarded (UMI too short)
    pub discarded_umi_too_short: u64,

    /// Records discarded (duplicate UMI at same position)
    pub discarded_duplicate_umi: u64,

    /// Number of unique molecule IDs assigned
    pub unique_molecule_ids: u64,

    /// Total number of UMI families/groups
    pub total_families: u64,

    /// Average reads per molecule
    pub avg_reads_per_molecule: f64,

    /// Median reads per molecule
    pub median_reads_per_molecule: u64,

    /// Minimum reads per molecule
    pub min_reads_per_molecule: u64,

    /// Maximum reads per molecule
    pub max_reads_per_molecule: u64,
}

impl UmiGroupingMetrics {
    /// Creates a new UMI grouping metrics struct with all counts initialized to zero.
    #[must_use]
    pub fn new() -> Self {
        Self {
            total_records: 0,
            accepted_records: 0,
            discarded_non_pf: 0,
            discarded_poor_alignment: 0,
            discarded_ns_in_umi: 0,
            discarded_umi_too_short: 0,
            discarded_duplicate_umi: 0,
            unique_molecule_ids: 0,
            total_families: 0,
            avg_reads_per_molecule: 0.0,
            median_reads_per_molecule: 0,
            min_reads_per_molecule: 0,
            max_reads_per_molecule: 0,
        }
    }
}

impl Default for UmiGroupingMetrics {
    fn default() -> Self {
        Self::new()
    }
}

impl Metric for UmiGroupingMetrics {
    fn metric_name() -> &'static str {
        "UMI grouping"
    }
}

impl super::ProcessingMetrics for UmiGroupingMetrics {
    fn total_input(&self) -> u64 {
        self.total_records
    }

    fn total_output(&self) -> u64 {
        self.accepted_records
    }

    fn total_filtered(&self) -> u64 {
        self.discarded_non_pf
            + self.discarded_poor_alignment
            + self.discarded_ns_in_umi
            + self.discarded_umi_too_short
            + self.discarded_duplicate_umi
    }
}

/// Family size distribution metrics.
///
/// Describes the distribution of UMI family sizes in the dataset.
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct FamilySizeMetrics {
    /// Family size (number of reads)
    pub family_size: usize,

    /// Number of families with this size
    pub count: u64,

    /// Fraction of all families with this size
    pub fraction: f64,

    /// Cumulative fraction (families with size >= this value)
    pub fraction_gt_or_eq_family_size: f64,
}

impl FamilySizeMetrics {
    /// Creates a new family size metric.
    #[must_use]
    pub fn new(family_size: usize) -> Self {
        Self { family_size, count: 0, fraction: 0.0, fraction_gt_or_eq_family_size: 0.0 }
    }
}

impl Default for FamilySizeMetrics {
    fn default() -> Self {
        Self::new(0)
    }
}

impl Metric for FamilySizeMetrics {
    fn metric_name() -> &'static str {
        "family size"
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_umi_grouping_metrics_new() {
        let metrics = UmiGroupingMetrics::new();
        assert_eq!(metrics.total_records, 0);
        assert_eq!(metrics.accepted_records, 0);
        assert_eq!(metrics.unique_molecule_ids, 0);
    }

    #[test]
    fn test_umi_grouping_metrics_default() {
        let metrics = UmiGroupingMetrics::default();
        assert_eq!(metrics.total_records, 0);
        assert_eq!(metrics.accepted_records, 0);
    }

    #[test]
    fn test_family_size_metrics_new() {
        let metrics = FamilySizeMetrics::new(5);
        assert_eq!(metrics.family_size, 5);
        assert_eq!(metrics.count, 0);
        assert!(metrics.fraction.abs() < f64::EPSILON);
    }

    #[test]
    fn test_metric_trait_impl() {
        assert_eq!(UmiGroupingMetrics::metric_name(), "UMI grouping");
        assert_eq!(FamilySizeMetrics::metric_name(), "family size");
    }
}
