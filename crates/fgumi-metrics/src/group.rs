//! Metrics for the `group` command.
//!
//! This module provides metrics types for tracking UMI grouping operations
//! and family size distributions.

use serde::{Deserialize, Serialize};

use crate::Metric;

/// Build a size distribution from (size, count) pairs.
///
/// Sorts by ascending size, computes per-entry fractions and reverse cumulative
/// fractions (fraction of entries with size >= this value), then returns a `Vec`
/// in ascending size order. The `ctor` closure maps `(size, count, fraction,
/// cumulative_fraction)` into the caller's metric type.
#[allow(clippy::cast_precision_loss)]
fn build_size_distribution<T>(
    counts: impl IntoIterator<Item = (usize, u64)>,
    ctor: impl Fn(usize, u64, f64, f64) -> T,
) -> Vec<T> {
    let mut sorted: Vec<_> = counts.into_iter().collect();
    sorted.sort_by_key(|(size, _)| *size);

    let total: f64 = sorted.iter().map(|(_, count)| *count as f64).sum();
    if total == 0.0 {
        return Vec::new();
    }

    let mut metrics = Vec::with_capacity(sorted.len());
    let mut cumulative = 0.0;
    for &(size, count) in sorted.iter().rev() {
        let fraction = count as f64 / total;
        cumulative += fraction;
        metrics.push(ctor(size, count, fraction, cumulative));
    }
    metrics.reverse();
    metrics
}

/// Metrics for UMI grouping operations.
///
/// These metrics track how reads are grouped by UMI and provide insight into
/// data quality and molecule representation.
#[derive(Debug, Clone, Default, Serialize, Deserialize)]
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
        Self::default()
    }
}

impl Metric for UmiGroupingMetrics {
    fn metric_name() -> &'static str {
        "UMI grouping"
    }
}

impl crate::ProcessingMetrics for UmiGroupingMetrics {
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
    }
}

/// Family size distribution metrics.
///
/// Describes the distribution of UMI family sizes in the dataset.
#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
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

    /// Build family size metrics from (`family_size`, count) pairs.
    ///
    /// Returns a `Vec` sorted by ascending family size, with cumulative
    /// fractions computed from largest to smallest.
    #[must_use]
    pub fn from_size_counts(counts: impl IntoIterator<Item = (usize, u64)>) -> Vec<Self> {
        build_size_distribution(counts, |family_size, count, fraction, cumulative| Self {
            family_size,
            count,
            fraction,
            fraction_gt_or_eq_family_size: cumulative,
        })
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

/// Position group size distribution metrics.
///
/// Describes how many distinct UMI families share the same genomic position
/// (start coordinate and strand). This provides insight into the complexity
/// of the library at each position.
#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
pub struct PositionGroupSizeMetrics {
    /// Position group size (number of families at a position)
    pub position_group_size: usize,

    /// Number of positions with this group size
    pub count: u64,

    /// Fraction of all positions with this group size
    pub fraction: f64,

    /// Cumulative fraction (positions with group size >= this value)
    pub fraction_gt_or_eq_position_group_size: f64,
}

impl PositionGroupSizeMetrics {
    /// Creates a new position group size metric.
    #[must_use]
    pub fn new(position_group_size: usize) -> Self {
        Self {
            position_group_size,
            count: 0,
            fraction: 0.0,
            fraction_gt_or_eq_position_group_size: 0.0,
        }
    }

    /// Build position group size metrics from (`position_group_size`, count) pairs.
    ///
    /// Returns a `Vec` sorted by ascending position group size, with cumulative
    /// fractions computed from largest to smallest.
    #[must_use]
    pub fn from_size_counts(counts: impl IntoIterator<Item = (usize, u64)>) -> Vec<Self> {
        build_size_distribution(counts, |position_group_size, count, fraction, cumulative| Self {
            position_group_size,
            count,
            fraction,
            fraction_gt_or_eq_position_group_size: cumulative,
        })
    }
}

impl Default for PositionGroupSizeMetrics {
    fn default() -> Self {
        Self::new(0)
    }
}

impl Metric for PositionGroupSizeMetrics {
    fn metric_name() -> &'static str {
        "position group size"
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

    #[test]
    fn test_from_size_counts() {
        let counts = vec![(3, 1u64), (1, 1), (2, 1)];
        let metrics = FamilySizeMetrics::from_size_counts(counts);
        assert_eq!(metrics.len(), 3);
        assert_eq!(metrics[0].family_size, 1);
        assert_eq!(metrics[1].family_size, 2);
        assert_eq!(metrics[2].family_size, 3);
        // Each is 1/3 of total
        assert!((metrics[0].fraction - 1.0 / 3.0).abs() < 1e-10);
        // Cumulative from largest: size 3 = 1/3, size 2 = 2/3, size 1 = 1.0
        assert!((metrics[0].fraction_gt_or_eq_family_size - 1.0).abs() < 1e-10);
        assert!((metrics[2].fraction_gt_or_eq_family_size - 1.0 / 3.0).abs() < 1e-10);
    }

    #[test]
    fn test_from_size_counts_empty() {
        let metrics = FamilySizeMetrics::from_size_counts(std::iter::empty());
        assert!(metrics.is_empty());
    }

    #[test]
    fn test_position_group_size_metrics_new() {
        let metrics = PositionGroupSizeMetrics::new(5);
        assert_eq!(metrics.position_group_size, 5);
        assert_eq!(metrics.count, 0);
        assert!(metrics.fraction.abs() < f64::EPSILON);
    }

    #[test]
    fn test_position_group_size_metric_name() {
        assert_eq!(PositionGroupSizeMetrics::metric_name(), "position group size");
    }

    #[test]
    fn test_position_group_size_from_size_counts() {
        let counts = vec![(3, 1u64), (1, 1), (2, 1)];
        let metrics = PositionGroupSizeMetrics::from_size_counts(counts);
        assert_eq!(metrics.len(), 3);
        assert_eq!(metrics[0].position_group_size, 1);
        assert_eq!(metrics[1].position_group_size, 2);
        assert_eq!(metrics[2].position_group_size, 3);
        assert!((metrics[0].fraction - 1.0 / 3.0).abs() < 1e-10);
        assert!((metrics[0].fraction_gt_or_eq_position_group_size - 1.0).abs() < 1e-10);
        assert!((metrics[2].fraction_gt_or_eq_position_group_size - 1.0 / 3.0).abs() < 1e-10);
    }

    #[test]
    fn test_position_group_size_from_size_counts_empty() {
        let metrics = PositionGroupSizeMetrics::from_size_counts(std::iter::empty());
        assert!(metrics.is_empty());
    }
}
