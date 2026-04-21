//! Metrics for the `simplex-metrics` command.
//!
//! This module provides QC metrics for simplex sequencing experiments,
//! including family size distributions, UMI frequencies, and simplex yield
//! at multiple sampling levels.

use serde::{Deserialize, Serialize};
use std::collections::HashMap;

use crate::inline_collector::MetricsCollector;
use crate::shared::{UmiCountTracker, UmiMetric};
use crate::template_info::TemplateMetadata;
use crate::{Metric, frac};

/// Metrics quantifying the distribution of CS and SS read family sizes.
///
/// Two kinds of families are described:
/// - **CS** (Coordinate & Strand): families grouped by unclipped 5' genomic positions and strands
/// - **SS** (Single Strand): single-strand families using UMIs, not linking opposing strands
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct SimplexFamilySizeMetric {
    /// The family size (number of read pairs grouped together)
    pub family_size: usize,
    /// Count of CS families with this size
    pub cs_count: usize,
    /// Fraction of all CS families with this size
    pub cs_fraction: f64,
    /// Fraction of CS families with size >= `family_size`
    pub cs_fraction_gt_or_eq_size: f64,
    /// Count of SS families with this size
    pub ss_count: usize,
    /// Fraction of all SS families with this size
    pub ss_fraction: f64,
    /// Fraction of SS families with size >= `family_size`
    pub ss_fraction_gt_or_eq_size: f64,
}

impl SimplexFamilySizeMetric {
    /// Creates a new simplex family size metric with all counts and fractions initialized to zero.
    #[must_use]
    pub fn new(family_size: usize) -> Self {
        Self {
            family_size,
            cs_count: 0,
            cs_fraction: 0.0,
            cs_fraction_gt_or_eq_size: 0.0,
            ss_count: 0,
            ss_fraction: 0.0,
            ss_fraction_gt_or_eq_size: 0.0,
        }
    }
}

impl Default for SimplexFamilySizeMetric {
    fn default() -> Self {
        Self::new(0)
    }
}

impl Metric for SimplexFamilySizeMetric {
    fn metric_name() -> &'static str {
        "simplex family size"
    }
}

/// Metrics sampled at various levels of coverage via random downsampling for simplex experiments.
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct SimplexYieldMetric {
    /// Approximate fraction of full dataset used
    pub fraction: f64,
    /// Number of read pairs upon which metrics are based
    pub read_pairs: usize,
    /// Number of CS (Coordinate & Strand) families
    pub cs_families: usize,
    /// Number of SS (Single-Strand by UMI) families
    pub ss_families: usize,
    /// Mean SS family size
    pub mean_ss_family_size: f64,
    /// Number of SS singleton families (size 1)
    pub ss_singletons: usize,
    /// Fraction of SS families that are singletons
    pub ss_singleton_fraction: f64,
    /// Number of SS families with size >= consensus minimum
    pub ss_consensus_families: usize,
}

impl SimplexYieldMetric {
    /// Creates a new simplex yield metric with the given fraction and all other fields zeroed.
    #[must_use]
    pub fn new(fraction: f64) -> Self {
        Self {
            fraction,
            read_pairs: 0,
            cs_families: 0,
            ss_families: 0,
            mean_ss_family_size: 0.0,
            ss_singletons: 0,
            ss_singleton_fraction: 0.0,
            ss_consensus_families: 0,
        }
    }
}

impl Default for SimplexYieldMetric {
    fn default() -> Self {
        Self::new(0.0)
    }
}

impl Metric for SimplexYieldMetric {
    fn metric_name() -> &'static str {
        "simplex yield"
    }
}

/// Collector for simplex sequencing metrics.
///
/// Tracks CS and SS family sizes and UMI frequencies.
pub struct SimplexMetricsCollector {
    /// CS family size counts: `family_size` -> count
    cs_family_sizes: HashMap<usize, usize>,
    /// SS family size counts: `family_size` -> count
    ss_family_sizes: HashMap<usize, usize>,
    /// UMI observation tracking
    umi_counts: UmiCountTracker,
}

impl SimplexMetricsCollector {
    /// Creates an empty metrics collector.
    #[must_use]
    pub fn new() -> Self {
        Self {
            cs_family_sizes: HashMap::new(),
            ss_family_sizes: HashMap::new(),
            umi_counts: UmiCountTracker::new(),
        }
    }

    /// Records a CS (Coordinate+Strand) family of the given size.
    pub fn record_cs_family(&mut self, size: usize) {
        *self.cs_family_sizes.entry(size).or_insert(0) += 1;
    }

    /// Records an SS (Single-Strand) family of the given size.
    pub fn record_ss_family(&mut self, size: usize) {
        *self.ss_family_sizes.entry(size).or_insert(0) += 1;
    }

    /// Records a UMI observation.
    ///
    /// # Arguments
    /// * `umi` - The consensus UMI sequence
    /// * `raw_count` - Number of raw observations of this UMI
    /// * `error_count` - Number of raw observations that had errors (differed from consensus)
    /// * `is_unique` - Whether this is a unique family observation
    pub fn record_umi(&mut self, umi: &str, raw_count: usize, error_count: usize, is_unique: bool) {
        self.umi_counts.record(umi, raw_count, error_count, is_unique);
    }

    /// Merge another collector into this one by summing all histogram counts and UMI observations.
    pub fn merge(&mut self, other: Self) {
        for (size, count) in other.cs_family_sizes {
            *self.cs_family_sizes.entry(size).or_insert(0) += count;
        }
        for (size, count) in other.ss_family_sizes {
            *self.ss_family_sizes.entry(size).or_insert(0) += count;
        }
        self.umi_counts.merge(other.umi_counts);
    }

    /// Generates family size metrics with fractions and cumulative fractions.
    ///
    /// Returns one [`SimplexFamilySizeMetric`] per observed family size (from 1 to the
    /// maximum across CS and SS), with per-size fractions and reverse-cumulative fractions.
    ///
    /// # Panics
    ///
    /// Panics if the internal array-of-two maximum computation fails, which cannot
    /// happen since the array is always non-empty.
    #[must_use]
    pub fn family_size_metrics(&self) -> Vec<SimplexFamilySizeMetric> {
        let max_size = self
            .cs_family_sizes
            .keys()
            .copied()
            .max()
            .unwrap_or(0)
            .max(self.ss_family_sizes.keys().copied().max().unwrap_or(0));

        if max_size == 0 {
            return Vec::new();
        }

        let cs_total: usize = self.cs_family_sizes.values().sum();
        let ss_total: usize = self.ss_family_sizes.values().sum();

        let mut metrics = Vec::new();
        for size in 1..=max_size {
            let mut metric = SimplexFamilySizeMetric::new(size);

            metric.cs_count = *self.cs_family_sizes.get(&size).unwrap_or(&0);
            metric.cs_fraction = frac(metric.cs_count, cs_total);

            metric.ss_count = *self.ss_family_sizes.get(&size).unwrap_or(&0);
            metric.ss_fraction = frac(metric.ss_count, ss_total);

            metrics.push(metric);
        }

        // Calculate cumulative fractions (>= size)
        for i in (0..metrics.len()).rev() {
            let next_coord_strand =
                if i + 1 < metrics.len() { metrics[i + 1].cs_fraction_gt_or_eq_size } else { 0.0 };
            let next_single_strand =
                if i + 1 < metrics.len() { metrics[i + 1].ss_fraction_gt_or_eq_size } else { 0.0 };

            metrics[i].cs_fraction_gt_or_eq_size = metrics[i].cs_fraction + next_coord_strand;
            metrics[i].ss_fraction_gt_or_eq_size = metrics[i].ss_fraction + next_single_strand;
        }

        metrics
    }

    /// Generates UMI metrics by delegating to the internal [`UmiCountTracker`].
    #[must_use]
    pub fn umi_metrics(&self) -> Vec<UmiMetric> {
        self.umi_counts.to_metrics()
    }
}

impl Default for SimplexMetricsCollector {
    fn default() -> Self {
        Self::new()
    }
}

impl MetricsCollector for SimplexMetricsCollector {
    fn record_group(&mut self, templates: &[TemplateMetadata<'_>]) {
        if templates.is_empty() {
            return;
        }

        // CS family: entire coordinate group
        self.record_cs_family(templates.len());

        // SS families: group by MI (base_umi with strand)
        let mut ss_counts: HashMap<&str, usize> = HashMap::new();
        for m in templates {
            *ss_counts.entry(m.template.mi.as_str()).or_insert(0) += 1;
        }
        for count in ss_counts.values() {
            self.record_ss_family(*count);
        }
    }

    fn merge(&mut self, other: Self) {
        // Delegate to the inherent merge method (inherent methods take priority)
        self.merge(other);
    }

    fn write_metrics(&self, _prefix: &std::path::Path) -> anyhow::Result<()> {
        Ok(())
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    // =========================================================================
    // SimplexFamilySizeMetric tests
    // =========================================================================

    #[test]
    fn test_simplex_family_size_metric_new() {
        let metric = SimplexFamilySizeMetric::new(5);
        assert_eq!(metric.family_size, 5);
        assert_eq!(metric.cs_count, 0);
        assert!(metric.cs_fraction.abs() < f64::EPSILON);
        assert!(metric.cs_fraction_gt_or_eq_size.abs() < f64::EPSILON);
        assert_eq!(metric.ss_count, 0);
        assert!(metric.ss_fraction.abs() < f64::EPSILON);
        assert!(metric.ss_fraction_gt_or_eq_size.abs() < f64::EPSILON);
    }

    // =========================================================================
    // SimplexYieldMetric tests
    // =========================================================================

    #[test]
    fn test_simplex_yield_metric_new() {
        let metric = SimplexYieldMetric::new(0.5);
        assert!((metric.fraction - 0.5).abs() < f64::EPSILON);
        assert_eq!(metric.read_pairs, 0);
        assert_eq!(metric.cs_families, 0);
        assert_eq!(metric.ss_families, 0);
        assert!(metric.mean_ss_family_size.abs() < f64::EPSILON);
        assert_eq!(metric.ss_singletons, 0);
        assert!(metric.ss_singleton_fraction.abs() < f64::EPSILON);
        assert_eq!(metric.ss_consensus_families, 0);
    }

    // =========================================================================
    // Metric trait tests
    // =========================================================================

    #[test]
    fn test_metric_trait_impl() {
        assert_eq!(SimplexFamilySizeMetric::metric_name(), "simplex family size");
        assert_eq!(SimplexYieldMetric::metric_name(), "simplex yield");
    }

    // =========================================================================
    // SimplexMetricsCollector tests
    // =========================================================================

    #[test]
    fn test_record_cs_family() {
        let mut collector = SimplexMetricsCollector::new();
        collector.record_cs_family(5);
        collector.record_cs_family(5);
        collector.record_cs_family(10);

        let metrics = collector.family_size_metrics();
        let size_5 = metrics
            .iter()
            .find(|m| m.family_size == 5)
            .expect("family_size 5 metric should be present");
        assert_eq!(size_5.cs_count, 2);
        let size_10 = metrics
            .iter()
            .find(|m| m.family_size == 10)
            .expect("family_size 10 metric should be present");
        assert_eq!(size_10.cs_count, 1);
    }

    #[test]
    fn test_record_ss_family() {
        let mut collector = SimplexMetricsCollector::new();
        collector.record_ss_family(3);
        collector.record_ss_family(3);
        collector.record_ss_family(3);

        let metrics = collector.family_size_metrics();
        let size_3 = metrics
            .iter()
            .find(|m| m.family_size == 3)
            .expect("family_size 3 metric should be present");
        assert_eq!(size_3.ss_count, 3);
    }

    #[test]
    fn test_family_size_metrics_fractions() {
        let mut collector = SimplexMetricsCollector::new();
        // 4 SS families: 2x size 1, 1x size 2, 1x size 3
        collector.record_ss_family(1);
        collector.record_ss_family(1);
        collector.record_ss_family(2);
        collector.record_ss_family(3);

        let metrics = collector.family_size_metrics();

        let size_1 = metrics
            .iter()
            .find(|m| m.family_size == 1)
            .expect("family_size 1 metric should be present");
        assert_eq!(size_1.ss_count, 2);
        assert!((size_1.ss_fraction - 0.5).abs() < 0.001); // 2/4 = 0.5
        assert!((size_1.ss_fraction_gt_or_eq_size - 1.0).abs() < 0.001); // All >= 1

        let size_2 = metrics
            .iter()
            .find(|m| m.family_size == 2)
            .expect("family_size 2 metric should be present");
        assert_eq!(size_2.ss_count, 1);
        assert!((size_2.ss_fraction - 0.25).abs() < 0.001); // 1/4 = 0.25
        assert!((size_2.ss_fraction_gt_or_eq_size - 0.5).abs() < 0.001); // size 2 + size 3

        let size_3 = metrics
            .iter()
            .find(|m| m.family_size == 3)
            .expect("family_size 3 metric should be present");
        assert_eq!(size_3.ss_count, 1);
        assert!((size_3.ss_fraction - 0.25).abs() < 0.001); // 1/4 = 0.25
        assert!((size_3.ss_fraction_gt_or_eq_size - 0.25).abs() < 0.001); // Only size 3
    }

    #[test]
    fn test_empty_collector() {
        let collector = SimplexMetricsCollector::new();

        let family_metrics = collector.family_size_metrics();
        assert!(family_metrics.is_empty());

        let umi_metrics = collector.umi_metrics();
        assert!(umi_metrics.is_empty());
    }

    #[test]
    fn test_record_umi() {
        let mut collector = SimplexMetricsCollector::new();
        collector.record_umi("AAAA", 10, 2, true);
        collector.record_umi("AAAA", 5, 1, false);
        collector.record_umi("CCCC", 8, 0, true);

        let metrics = collector.umi_metrics();
        assert_eq!(metrics.len(), 2);

        let aaaa =
            metrics.iter().find(|m| m.umi == "AAAA").expect("AAAA UMI metric should be present");
        assert_eq!(aaaa.raw_observations, 15); // 10 + 5
        assert_eq!(aaaa.raw_observations_with_errors, 3); // 2 + 1
        assert_eq!(aaaa.unique_observations, 1); // Only one unique

        let cccc =
            metrics.iter().find(|m| m.umi == "CCCC").expect("CCCC UMI metric should be present");
        assert_eq!(cccc.raw_observations, 8);
        assert_eq!(cccc.unique_observations, 1);
    }

    // =========================================================================
    // SimplexMetricsCollector::merge tests
    // =========================================================================

    #[test]
    fn test_simplex_metrics_collector_merge() {
        let mut a = SimplexMetricsCollector::new();
        a.record_cs_family(3);
        a.record_cs_family(3);
        a.record_ss_family(2);
        a.record_umi("AAAA", 5, 1, true);

        let mut b = SimplexMetricsCollector::new();
        b.record_cs_family(3);
        b.record_cs_family(5);
        b.record_ss_family(2);
        b.record_ss_family(4);
        b.record_umi("AAAA", 3, 0, true);
        b.record_umi("CCCC", 2, 0, true);

        a.merge(b);

        let metrics = a.family_size_metrics();
        let cs3 = metrics.iter().find(|m| m.family_size == 3).unwrap();
        assert_eq!(cs3.cs_count, 3);
        let cs5 = metrics.iter().find(|m| m.family_size == 5).unwrap();
        assert_eq!(cs5.cs_count, 1);
        let ss2 = metrics.iter().find(|m| m.family_size == 2).unwrap();
        assert_eq!(ss2.ss_count, 2);
        let ss4 = metrics.iter().find(|m| m.family_size == 4).unwrap();
        assert_eq!(ss4.ss_count, 1);

        let umi_metrics = a.umi_metrics();
        assert_eq!(umi_metrics.len(), 2);
    }

    #[test]
    fn test_simplex_implements_metrics_collector() {
        use crate::inline_collector::MetricsCollector;
        use crate::template_info::{TemplateInfo, TemplateMetadata};

        let mut collector = SimplexMetricsCollector::new();
        let t1 = TemplateInfo {
            mi: "1/A".to_string(),
            rx: "AAAA".to_string(),
            ref_name: Some("chr1".to_string()),
            position: Some(100),
            end_position: Some(200),
            hash_fraction: 0.5,
            read_info_key: crate::template_info::ReadInfoKey::default(),
        };
        let t2 = TemplateInfo {
            mi: "1/A".to_string(),
            rx: "AAAA".to_string(),
            ref_name: Some("chr1".to_string()),
            position: Some(100),
            end_position: Some(200),
            hash_fraction: 0.6,
            read_info_key: crate::template_info::ReadInfoKey::default(),
        };
        let t3 = TemplateInfo {
            mi: "2/A".to_string(),
            rx: "CCCC".to_string(),
            ref_name: Some("chr1".to_string()),
            position: Some(100),
            end_position: Some(200),
            hash_fraction: 0.7,
            read_info_key: crate::template_info::ReadInfoKey::default(),
        };

        let metadata: Vec<TemplateMetadata<'_>> = vec![
            TemplateMetadata {
                template: &t1,
                base_umi: "1",
                is_a_strand: true,
                is_b_strand: false,
            },
            TemplateMetadata {
                template: &t2,
                base_umi: "1",
                is_a_strand: true,
                is_b_strand: false,
            },
            TemplateMetadata {
                template: &t3,
                base_umi: "2",
                is_a_strand: true,
                is_b_strand: false,
            },
        ];

        collector.record_group(&metadata);

        let metrics = collector.family_size_metrics();
        // CS family: entire group = 3 templates
        let cs3 = metrics.iter().find(|m| m.family_size == 3).unwrap();
        assert_eq!(cs3.cs_count, 1);
        // SS families: MI "1/A" has 2 templates, MI "2/A" has 1 template
        let ss2 = metrics.iter().find(|m| m.family_size == 2).unwrap();
        assert_eq!(ss2.ss_count, 1);
        let ss1 = metrics.iter().find(|m| m.family_size == 1).unwrap();
        assert_eq!(ss1.ss_count, 1);
    }

    #[test]
    fn test_simplex_metrics_collector_merge_empty() {
        let mut a = SimplexMetricsCollector::new();
        a.record_cs_family(3);
        let b = SimplexMetricsCollector::new();
        a.merge(b);
        let metrics = a.family_size_metrics();
        // family_size_metrics returns entries for sizes 1..=max, so 3 entries
        assert_eq!(metrics.len(), 3);
        let cs3 = metrics.iter().find(|m| m.family_size == 3).unwrap();
        assert_eq!(cs3.cs_count, 1);
    }
}
