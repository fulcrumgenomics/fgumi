//! Metrics for the `duplex_metrics` command.
//!
//! This module provides comprehensive QC metrics for duplex sequencing experiments,
//! including family size distributions, UMI frequencies, and duplex yield at multiple
//! sampling levels.

use serde::{Deserialize, Serialize};
use std::collections::HashMap;

use crate::inline_collector::MetricsCollector;
use crate::shared::{UmiCountTracker, UmiMetric};
use crate::template_info::TemplateMetadata;
use crate::{Metric, frac};

/// Metrics quantifying the distribution of different kinds of read family sizes.
///
/// Three kinds of families are described:
/// - **CS** (Coordinate & Strand): families grouped by unclipped 5' genomic positions and strands
/// - **SS** (Single Strand): single-strand families using UMIs, not linking opposing strands
/// - **DS** (Double Strand): families combining single-strand families from opposite strands
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct FamilySizeMetric {
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
    /// Count of DS families with this size
    pub ds_count: usize,
    /// Fraction of all DS families with this size
    pub ds_fraction: f64,
    /// Fraction of DS families with size >= `family_size`
    pub ds_fraction_gt_or_eq_size: f64,
}

impl FamilySizeMetric {
    /// Creates a new family size metric
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
            ds_count: 0,
            ds_fraction: 0.0,
            ds_fraction_gt_or_eq_size: 0.0,
        }
    }
}

impl Default for FamilySizeMetric {
    fn default() -> Self {
        Self::new(0)
    }
}

impl Metric for FamilySizeMetric {
    fn metric_name() -> &'static str {
        "duplex family size"
    }
}

/// Metrics describing double-stranded (duplex) tag families by AB and BA strand sizes.
///
/// For a given tag family, `ab` is the larger sub-family and `ba` is the smaller one.
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct DuplexFamilySizeMetric {
    /// Number of reads in the AB sub-family (larger)
    pub ab_size: usize,
    /// Number of reads in the BA sub-family (smaller)
    pub ba_size: usize,
    /// Count of families with these AB/BA sizes
    pub count: usize,
    /// Fraction of all duplex families with these sizes
    pub fraction: f64,
    /// Fraction of duplex families with AB >= `ab_size` and BA >= `ba_size`
    pub fraction_gt_or_eq_size: f64,
}

impl DuplexFamilySizeMetric {
    /// Creates a new duplex family size metric
    #[must_use]
    pub fn new(ab_size: usize, ba_size: usize) -> Self {
        Self { ab_size, ba_size, count: 0, fraction: 0.0, fraction_gt_or_eq_size: 0.0 }
    }
}

impl Default for DuplexFamilySizeMetric {
    fn default() -> Self {
        Self::new(0, 0)
    }
}

impl Metric for DuplexFamilySizeMetric {
    fn metric_name() -> &'static str {
        "duplex AB/BA family size"
    }
}

impl Ord for DuplexFamilySizeMetric {
    fn cmp(&self, other: &Self) -> std::cmp::Ordering {
        self.ab_size.cmp(&other.ab_size).then_with(|| self.ba_size.cmp(&other.ba_size))
    }
}

impl PartialOrd for DuplexFamilySizeMetric {
    fn partial_cmp(&self, other: &Self) -> Option<std::cmp::Ordering> {
        Some(self.cmp(other))
    }
}

impl Eq for DuplexFamilySizeMetric {}

impl PartialEq for DuplexFamilySizeMetric {
    fn eq(&self, other: &Self) -> bool {
        self.ab_size == other.ab_size && self.ba_size == other.ba_size
    }
}

/// Metrics sampled at various levels of coverage via random downsampling.
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct DuplexYieldMetric {
    /// Approximate fraction of full dataset used
    pub fraction: f64,
    /// Number of read pairs upon which metrics are based
    pub read_pairs: usize,
    /// Number of CS (Coordinate & Strand) families
    pub cs_families: usize,
    /// Number of SS (Single-Strand by UMI) families
    pub ss_families: usize,
    /// Number of DS (Double-Strand by UMI) families
    pub ds_families: usize,
    /// Number of DS families that are duplexes (min reads on both strands)
    pub ds_duplexes: usize,
    /// Fraction of DS families that are duplexes
    pub ds_fraction_duplexes: f64,
    /// Expected fraction of DS families that should be duplexes under ideal model
    pub ds_fraction_duplexes_ideal: f64,
}

impl DuplexYieldMetric {
    /// Creates a new yield metric
    #[must_use]
    pub fn new(fraction: f64) -> Self {
        Self {
            fraction,
            read_pairs: 0,
            cs_families: 0,
            ss_families: 0,
            ds_families: 0,
            ds_duplexes: 0,
            ds_fraction_duplexes: 0.0,
            ds_fraction_duplexes_ideal: 0.0,
        }
    }
}

impl Default for DuplexYieldMetric {
    fn default() -> Self {
        Self::new(0.0)
    }
}

impl Metric for DuplexYieldMetric {
    fn metric_name() -> &'static str {
        "duplex yield"
    }
}

/// Metrics describing observed duplex UMI sequences and their frequencies.
///
/// Duplex UMIs are normalized to F1R2 orientation (positive strand first).
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct DuplexUmiMetric {
    /// The duplex UMI sequence (possibly corrected, F1R2 normalized)
    pub umi: String,
    /// Number of read pairs observing this duplex UMI
    pub raw_observations: usize,
    /// Subset of raw observations that underwent correction
    pub raw_observations_with_errors: usize,
    /// Number of double-stranded tag families observing this duplex UMI
    pub unique_observations: usize,
    /// Fraction of all raw observations
    pub fraction_raw_observations: f64,
    /// Fraction of all unique observations
    pub fraction_unique_observations: f64,
    /// Expected fraction based on individual UMI frequencies
    pub fraction_unique_observations_expected: f64,
}

impl DuplexUmiMetric {
    /// Creates a new duplex UMI metric
    #[must_use]
    pub fn new(umi: String) -> Self {
        Self {
            umi,
            raw_observations: 0,
            raw_observations_with_errors: 0,
            unique_observations: 0,
            fraction_raw_observations: 0.0,
            fraction_unique_observations: 0.0,
            fraction_unique_observations_expected: 0.0,
        }
    }
}

impl Default for DuplexUmiMetric {
    fn default() -> Self {
        Self::new(String::new())
    }
}

impl Metric for DuplexUmiMetric {
    fn metric_name() -> &'static str {
        "duplex UMI"
    }
}

/// Collector for duplex sequencing metrics.
///
/// Tracks family sizes, UMI frequencies, and yields at multiple sampling levels.
pub struct DuplexMetricsCollector {
    /// Whether to collect duplex UMI counts (memory intensive)
    collect_duplex_umi_counts: bool,

    // Family size tracking
    cs_family_sizes: HashMap<usize, usize>,
    ss_family_sizes: HashMap<usize, usize>,
    ds_family_sizes: HashMap<usize, usize>,
    duplex_family_sizes: HashMap<(usize, usize), usize>,

    // UMI tracking
    umi_counts: UmiCountTracker,
    duplex_umi_counts: UmiCountTracker,
}

impl DuplexMetricsCollector {
    /// Creates a new metrics collector
    #[must_use]
    pub fn new(collect_duplex_umi_counts: bool) -> Self {
        Self {
            collect_duplex_umi_counts,
            cs_family_sizes: HashMap::new(),
            ss_family_sizes: HashMap::new(),
            ds_family_sizes: HashMap::new(),
            duplex_family_sizes: HashMap::new(),
            umi_counts: UmiCountTracker::new(),
            duplex_umi_counts: UmiCountTracker::new(),
        }
    }

    /// Merge another collector into this one by summing all histogram counts and UMI observations.
    pub fn merge(&mut self, other: Self) {
        for (size, count) in other.cs_family_sizes {
            *self.cs_family_sizes.entry(size).or_insert(0) += count;
        }
        for (size, count) in other.ss_family_sizes {
            *self.ss_family_sizes.entry(size).or_insert(0) += count;
        }
        for (size, count) in other.ds_family_sizes {
            *self.ds_family_sizes.entry(size).or_insert(0) += count;
        }
        for (key, count) in other.duplex_family_sizes {
            *self.duplex_family_sizes.entry(key).or_insert(0) += count;
        }
        self.umi_counts.merge(other.umi_counts);
        self.duplex_umi_counts.merge(other.duplex_umi_counts);
    }

    /// Records a CS (Coordinate+Strand) family
    pub fn record_cs_family(&mut self, size: usize) {
        *self.cs_family_sizes.entry(size).or_insert(0) += 1;
    }

    /// Records an SS (Single-Strand) family
    pub fn record_ss_family(&mut self, size: usize) {
        *self.ss_family_sizes.entry(size).or_insert(0) += 1;
    }

    /// Records a DS (Double-Strand) family
    pub fn record_ds_family(&mut self, size: usize) {
        *self.ds_family_sizes.entry(size).or_insert(0) += 1;
    }

    /// Records a duplex family with AB and BA sizes
    pub fn record_duplex_family(&mut self, ab_size: usize, ba_size: usize) {
        // Ensure ab >= ba by convention
        let (ab, ba) = if ab_size >= ba_size { (ab_size, ba_size) } else { (ba_size, ab_size) };
        *self.duplex_family_sizes.entry((ab, ba)).or_insert(0) += 1;
    }

    /// Records a UMI observation
    ///
    /// # Arguments
    /// * `umi` - The consensus UMI sequence
    /// * `raw_count` - Number of raw observations of this UMI
    /// * `error_count` - Number of raw observations that had errors (differed from consensus)
    /// * `is_unique` - Whether this is a unique family observation
    pub fn record_umi(&mut self, umi: &str, raw_count: usize, error_count: usize, is_unique: bool) {
        self.umi_counts.record(umi, raw_count, error_count, is_unique);
    }

    /// Records a duplex UMI observation
    ///
    /// # Arguments
    /// * `umi` - The duplex consensus UMI sequence
    /// * `raw_count` - Number of raw observations of this duplex UMI
    /// * `error_count` - Number of raw observations that had errors
    /// * `is_unique` - Whether this is a unique family observation
    pub fn record_duplex_umi(
        &mut self,
        umi: &str,
        raw_count: usize,
        error_count: usize,
        is_unique: bool,
    ) {
        if !self.collect_duplex_umi_counts {
            return;
        }
        self.duplex_umi_counts.record(umi, raw_count, error_count, is_unique);
    }

    /// Generates family size metrics
    ///
    /// # Panics
    ///
    /// Panics if the internal array-of-three maximum computation fails, which cannot
    /// happen since the array is always non-empty.
    #[must_use]
    pub fn family_size_metrics(&self) -> Vec<FamilySizeMetric> {
        // Find max family size across all types
        let max_size = *[
            self.cs_family_sizes.keys().max().unwrap_or(&0),
            self.ss_family_sizes.keys().max().unwrap_or(&0),
            self.ds_family_sizes.keys().max().unwrap_or(&0),
        ]
        .iter()
        .max()
        .expect("array of three elements always has a maximum");

        let coord_strand_total: usize = self.cs_family_sizes.values().sum();
        let single_strand_total: usize = self.ss_family_sizes.values().sum();
        let double_strand_total: usize = self.ds_family_sizes.values().sum();

        let mut metrics = Vec::new();
        for size in 1..=*max_size {
            let mut metric = FamilySizeMetric::new(size);

            metric.cs_count = *self.cs_family_sizes.get(&size).unwrap_or(&0);
            metric.cs_fraction = frac(metric.cs_count, coord_strand_total);

            metric.ss_count = *self.ss_family_sizes.get(&size).unwrap_or(&0);
            metric.ss_fraction = frac(metric.ss_count, single_strand_total);

            metric.ds_count = *self.ds_family_sizes.get(&size).unwrap_or(&0);
            metric.ds_fraction = frac(metric.ds_count, double_strand_total);

            metrics.push(metric);
        }

        // Calculate cumulative fractions (>= size)
        for i in (0..metrics.len()).rev() {
            let next_coord_strand =
                if i + 1 < metrics.len() { metrics[i + 1].cs_fraction_gt_or_eq_size } else { 0.0 };
            let next_single_strand =
                if i + 1 < metrics.len() { metrics[i + 1].ss_fraction_gt_or_eq_size } else { 0.0 };
            let next_double_strand =
                if i + 1 < metrics.len() { metrics[i + 1].ds_fraction_gt_or_eq_size } else { 0.0 };

            metrics[i].cs_fraction_gt_or_eq_size = metrics[i].cs_fraction + next_coord_strand;
            metrics[i].ss_fraction_gt_or_eq_size = metrics[i].ss_fraction + next_single_strand;
            metrics[i].ds_fraction_gt_or_eq_size = metrics[i].ds_fraction + next_double_strand;
        }

        metrics
    }

    /// Generates duplex family size metrics
    #[must_use]
    pub fn duplex_family_size_metrics(&self) -> Vec<DuplexFamilySizeMetric> {
        let total: usize = self.duplex_family_sizes.values().sum();

        let mut metrics: Vec<_> = self
            .duplex_family_sizes
            .iter()
            .map(|((ab, ba), count)| {
                let mut metric = DuplexFamilySizeMetric::new(*ab, *ba);
                metric.count = *count;
                metric.fraction = frac(*count, total);
                metric
            })
            .collect();

        metrics.sort();

        // Calculate 2D cumulative fractions: fraction of families with AB >= ab AND BA >= ba
        // This matches fgbio's definition in DuplexFamilySizeMetric
        //
        // Build a 2D suffix sum grid to avoid O(n²) per-metric iteration.
        if total > 0 {
            let max_ab = self.duplex_family_sizes.keys().map(|(a, _)| *a).max().unwrap_or(0);
            let max_ba = self.duplex_family_sizes.keys().map(|(_, b)| *b).max().unwrap_or(0);
            let cols = max_ba + 1;
            let mut grid = vec![0usize; (max_ab + 1) * cols];
            for (&(a, b), &count) in &self.duplex_family_sizes {
                grid[a * cols + b] = count;
            }
            // Accumulate suffix sums: first along ba (columns right-to-left),
            // then along ab (rows bottom-to-top)
            for a in 0..=max_ab {
                for b in (0..max_ba).rev() {
                    grid[a * cols + b] += grid[a * cols + b + 1];
                }
            }
            for b in 0..=max_ba {
                for a in (0..max_ab).rev() {
                    grid[a * cols + b] += grid[(a + 1) * cols + b];
                }
            }
            for metric in &mut metrics {
                let cumulative_count = grid[metric.ab_size * cols + metric.ba_size];
                metric.fraction_gt_or_eq_size = frac(cumulative_count, total);
            }
        }

        metrics
    }

    /// Generates UMI metrics
    #[must_use]
    pub fn umi_metrics(&self) -> Vec<UmiMetric> {
        self.umi_counts.to_metrics()
    }

    /// Generates duplex UMI metrics
    ///
    /// # Arguments
    /// * `umi_metrics` - Individual UMI metrics used to calculate expected frequencies
    #[must_use]
    pub fn duplex_umi_metrics(&self, umi_metrics: &[UmiMetric]) -> Vec<DuplexUmiMetric> {
        if !self.collect_duplex_umi_counts {
            return Vec::new();
        }

        // Build a map of individual UMI -> fraction_unique_observations for lookup
        let single_umi_fractions: HashMap<&str, f64> =
            umi_metrics.iter().map(|m| (m.umi.as_str(), m.fraction_unique_observations)).collect();

        let total_raw = self.duplex_umi_counts.total_raw();
        let total_unique = self.duplex_umi_counts.total_unique();

        let mut metrics: Vec<_> = self
            .duplex_umi_counts
            .iter()
            .map(|(umi, raw, errors, unique)| {
                let mut metric = DuplexUmiMetric::new(umi.to_string());
                metric.raw_observations = raw;
                metric.raw_observations_with_errors = errors;
                metric.unique_observations = unique;
                metric.fraction_raw_observations = frac(raw, total_raw);
                metric.fraction_unique_observations = frac(unique, total_unique);

                // Calculate expected fraction based on individual UMI frequencies
                // Expected frequency = freq(umi1) * freq(umi2) assuming independence
                metric.fraction_unique_observations_expected =
                    if let Some((umi1, umi2)) = umi.split_once('-') {
                        let freq1 = single_umi_fractions.get(umi1).copied().unwrap_or(0.0);
                        let freq2 = single_umi_fractions.get(umi2).copied().unwrap_or(0.0);
                        freq1 * freq2
                    } else {
                        0.0
                    };

                metric
            })
            .collect();

        // Sort by unique observations descending (matching Scala's sort order)
        metrics.sort_by(|a, b| b.unique_observations.cmp(&a.unique_observations));
        metrics
    }
}

impl Default for DuplexMetricsCollector {
    fn default() -> Self {
        Self::new(false)
    }
}

impl MetricsCollector for DuplexMetricsCollector {
    fn record_group(&mut self, templates: &[TemplateMetadata<'_>]) {
        if templates.is_empty() {
            return;
        }

        // CS family: entire coordinate group
        self.record_cs_family(templates.len());

        // SS families: group by full MI (including strand suffix)
        let mut ss_counts: HashMap<&str, usize> = HashMap::new();
        for m in templates {
            *ss_counts.entry(m.template.mi.as_str()).or_insert(0) += 1;
        }
        for count in ss_counts.values() {
            self.record_ss_family(*count);
        }

        // DS families: group by base_umi, count A and B strands
        let mut ds_groups: HashMap<&str, (usize, usize)> = HashMap::new();
        for m in templates {
            let entry = ds_groups.entry(m.base_umi).or_insert((0, 0));
            if m.is_a_strand {
                entry.0 += 1;
            }
            if m.is_b_strand {
                entry.1 += 1;
            }
        }
        for (ab_count, ba_count) in ds_groups.values() {
            let total = ab_count + ba_count;
            self.record_ds_family(total);
            self.record_duplex_family(*ab_count, *ba_count);
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
    // FamilySizeMetric tests
    // =========================================================================

    #[test]
    fn test_family_size_metric_new() {
        let metric = FamilySizeMetric::new(5);
        assert_eq!(metric.family_size, 5);
        assert_eq!(metric.cs_count, 0);
        assert!(metric.cs_fraction.abs() < f64::EPSILON);
        assert_eq!(metric.ss_count, 0);
        assert_eq!(metric.ds_count, 0);
    }

    // =========================================================================
    // DuplexFamilySizeMetric tests
    // =========================================================================

    #[test]
    fn test_duplex_family_size_metric_new() {
        let metric = DuplexFamilySizeMetric::new(10, 5);
        assert_eq!(metric.ab_size, 10);
        assert_eq!(metric.ba_size, 5);
        assert_eq!(metric.count, 0);
        assert!(metric.fraction.abs() < f64::EPSILON);
    }

    #[test]
    fn test_duplex_family_size_metric_ordering() {
        let m1 = DuplexFamilySizeMetric::new(5, 3);
        let m2 = DuplexFamilySizeMetric::new(5, 4);
        let m3 = DuplexFamilySizeMetric::new(6, 2);

        // m1 < m2 (same ab, smaller ba)
        assert!(m1 < m2);
        // m1 < m3 (smaller ab)
        assert!(m1 < m3);
        // m2 < m3 (smaller ab)
        assert!(m2 < m3);
    }

    #[test]
    fn test_duplex_family_size_metric_equality() {
        let m1 = DuplexFamilySizeMetric::new(5, 3);
        let m2 = DuplexFamilySizeMetric::new(5, 3);
        let m3 = DuplexFamilySizeMetric::new(5, 4);

        assert_eq!(m1, m2);
        assert_ne!(m1, m3);
    }

    // =========================================================================
    // DuplexYieldMetric tests
    // =========================================================================

    #[test]
    fn test_duplex_yield_metric_new() {
        let metric = DuplexYieldMetric::new(0.5);
        assert!((metric.fraction - 0.5).abs() < f64::EPSILON);
        assert_eq!(metric.read_pairs, 0);
        assert_eq!(metric.cs_families, 0);
        assert_eq!(metric.ds_duplexes, 0);
    }

    // =========================================================================
    // DuplexUmiMetric tests
    // =========================================================================

    #[test]
    fn test_duplex_umi_metric_new() {
        let metric = DuplexUmiMetric::new("ACGT-TGCA".to_string());
        assert_eq!(metric.umi, "ACGT-TGCA");
        assert_eq!(metric.raw_observations, 0);
        assert!(metric.fraction_unique_observations_expected.abs() < f64::EPSILON);
    }

    // =========================================================================
    // DuplexMetricsCollector tests
    // =========================================================================

    #[test]
    fn test_record_cs_family() {
        let mut collector = DuplexMetricsCollector::new(false);
        collector.record_cs_family(5);
        collector.record_cs_family(5);
        collector.record_cs_family(10);

        let metrics = collector.family_size_metrics();
        // Size 5 should have count 2
        let size_5 = metrics
            .iter()
            .find(|m| m.family_size == 5)
            .expect("family_size 5 metric should be present");
        assert_eq!(size_5.cs_count, 2);
        // Size 10 should have count 1
        let size_10 = metrics
            .iter()
            .find(|m| m.family_size == 10)
            .expect("family_size 10 metric should be present");
        assert_eq!(size_10.cs_count, 1);
    }

    #[test]
    fn test_record_ss_family() {
        let mut collector = DuplexMetricsCollector::new(false);
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
    fn test_record_ds_family() {
        let mut collector = DuplexMetricsCollector::new(false);
        collector.record_ds_family(2);

        let metrics = collector.family_size_metrics();
        let size_2 = metrics
            .iter()
            .find(|m| m.family_size == 2)
            .expect("family_size 2 metric should be present");
        assert_eq!(size_2.ds_count, 1);
    }

    #[test]
    fn test_record_duplex_family_normalization() {
        let mut collector = DuplexMetricsCollector::new(false);
        // Record with ab < ba - should be normalized
        collector.record_duplex_family(3, 5);
        // Record with ab > ba - already normalized
        collector.record_duplex_family(5, 3);

        let metrics = collector.duplex_family_size_metrics();
        // Both should map to (5, 3)
        assert_eq!(metrics.len(), 1);
        assert_eq!(metrics[0].ab_size, 5);
        assert_eq!(metrics[0].ba_size, 3);
        assert_eq!(metrics[0].count, 2);
    }

    #[test]
    fn test_record_umi() {
        let mut collector = DuplexMetricsCollector::new(false);
        collector.record_umi("AAAA", 10, 2, true);
        collector.record_umi("AAAA", 5, 1, false); // Not unique
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

    #[test]
    fn test_record_duplex_umi_disabled() {
        let mut collector = DuplexMetricsCollector::new(false); // Disabled
        collector.record_duplex_umi("AAAA-TTTT", 10, 0, true);

        let umi_metrics = collector.umi_metrics();
        let duplex_metrics = collector.duplex_umi_metrics(&umi_metrics);
        assert!(duplex_metrics.is_empty());
    }

    #[test]
    fn test_record_duplex_umi_enabled() {
        let mut collector = DuplexMetricsCollector::new(true); // Enabled
        collector.record_duplex_umi("AAAA-TTTT", 10, 2, true);

        // Need to also record individual UMIs for expected calculation
        collector.record_umi("AAAA", 5, 0, true);
        collector.record_umi("TTTT", 5, 0, true);

        let umi_metrics = collector.umi_metrics();
        let duplex_metrics = collector.duplex_umi_metrics(&umi_metrics);
        assert_eq!(duplex_metrics.len(), 1);
        assert_eq!(duplex_metrics[0].umi, "AAAA-TTTT");
        assert_eq!(duplex_metrics[0].raw_observations, 10);
    }

    #[test]
    fn test_family_size_metrics_fractions() {
        let mut collector = DuplexMetricsCollector::new(false);
        // 4 families total: 2 size-1, 1 size-2, 1 size-3
        collector.record_cs_family(1);
        collector.record_cs_family(1);
        collector.record_cs_family(2);
        collector.record_cs_family(3);

        let metrics = collector.family_size_metrics();

        let size_1 = metrics
            .iter()
            .find(|m| m.family_size == 1)
            .expect("family_size 1 metric should be present");
        assert_eq!(size_1.cs_count, 2);
        assert!((size_1.cs_fraction - 0.5).abs() < 0.001); // 2/4 = 0.5
        assert!((size_1.cs_fraction_gt_or_eq_size - 1.0).abs() < 0.001); // All >= 1

        let size_3 = metrics
            .iter()
            .find(|m| m.family_size == 3)
            .expect("family_size 3 metric should be present");
        assert_eq!(size_3.cs_count, 1);
        assert!((size_3.cs_fraction - 0.25).abs() < 0.001); // 1/4 = 0.25
        assert!((size_3.cs_fraction_gt_or_eq_size - 0.25).abs() < 0.001); // Only size 3
    }

    #[test]
    fn test_duplex_family_size_metrics_sorting() {
        let mut collector = DuplexMetricsCollector::new(false);
        collector.record_duplex_family(5, 3);
        collector.record_duplex_family(2, 1);
        collector.record_duplex_family(5, 2);

        let metrics = collector.duplex_family_size_metrics();
        // Should be sorted by (ab, ba)
        assert_eq!(metrics[0].ab_size, 2);
        assert_eq!(metrics[0].ba_size, 1);
        assert_eq!(metrics[1].ab_size, 5);
        assert_eq!(metrics[1].ba_size, 2);
        assert_eq!(metrics[2].ab_size, 5);
        assert_eq!(metrics[2].ba_size, 3);
    }

    #[test]
    fn test_duplex_umi_expected_frequency() {
        let mut collector = DuplexMetricsCollector::new(true);

        // Individual UMIs with known frequencies
        collector.record_umi("AAAA", 10, 0, true);
        collector.record_umi("TTTT", 10, 0, true);
        // Total unique = 2, so each has fraction 0.5

        // Duplex UMI
        collector.record_duplex_umi("AAAA-TTTT", 5, 0, true);

        let umi_metrics = collector.umi_metrics();
        let duplex_metrics = collector.duplex_umi_metrics(&umi_metrics);

        assert_eq!(duplex_metrics.len(), 1);
        // Expected = 0.5 * 0.5 = 0.25
        assert!((duplex_metrics[0].fraction_unique_observations_expected - 0.25).abs() < 0.001);
    }

    #[test]
    fn test_empty_collector() {
        let collector = DuplexMetricsCollector::new(false);

        let family_metrics = collector.family_size_metrics();
        assert!(family_metrics.is_empty());

        let duplex_metrics = collector.duplex_family_size_metrics();
        assert!(duplex_metrics.is_empty());

        let umi_metrics = collector.umi_metrics();
        assert!(umi_metrics.is_empty());
    }

    #[test]
    fn test_metric_trait_impl() {
        assert_eq!(FamilySizeMetric::metric_name(), "duplex family size");
        assert_eq!(DuplexFamilySizeMetric::metric_name(), "duplex AB/BA family size");
        assert_eq!(DuplexYieldMetric::metric_name(), "duplex yield");
        assert_eq!(UmiMetric::metric_name(), "UMI");
        assert_eq!(DuplexUmiMetric::metric_name(), "duplex UMI");
    }

    #[test]
    fn test_duplex_implements_metrics_collector() {
        use crate::inline_collector::MetricsCollector;
        use crate::template_info::{TemplateInfo, TemplateMetadata};

        let mut collector = DuplexMetricsCollector::new(false);
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
            mi: "1/B".to_string(),
            rx: "AAAA".to_string(),
            ref_name: Some("chr1".to_string()),
            position: Some(100),
            end_position: Some(200),
            hash_fraction: 0.6,
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
                is_a_strand: false,
                is_b_strand: true,
            },
        ];

        collector.record_group(&metadata);

        let metrics = collector.family_size_metrics();
        // CS family: entire group = 2 templates
        let cs2 = metrics.iter().find(|m| m.family_size == 2).unwrap();
        assert_eq!(cs2.cs_count, 1);
        // SS families: MI "1/A" -> 1, MI "1/B" -> 1
        let ss1 = metrics.iter().find(|m| m.family_size == 1).unwrap();
        assert_eq!(ss1.ss_count, 2);

        // DS family: base_umi "1" has 1 A-strand and 1 B-strand -> DS size 2
        let ds2 = metrics.iter().find(|m| m.family_size == 2).unwrap();
        assert_eq!(ds2.ds_count, 1);

        // Duplex family: ab=1, ba=1
        let duplex = collector.duplex_family_size_metrics();
        let d11 = duplex.iter().find(|m| m.ab_size == 1 && m.ba_size == 1).unwrap();
        assert_eq!(d11.count, 1);
    }

    // =========================================================================
    // DuplexMetricsCollector::merge tests
    // =========================================================================

    #[test]
    fn test_duplex_metrics_collector_merge() {
        let mut a = DuplexMetricsCollector::new(false);
        a.record_cs_family(3);
        a.record_ss_family(2);
        a.record_ds_family(1);
        a.record_duplex_family(2, 3);
        a.record_umi("AAAA", 5, 1, true);

        let mut b = DuplexMetricsCollector::new(false);
        b.record_cs_family(3);
        b.record_ss_family(4);
        b.record_ds_family(1);
        b.record_duplex_family(2, 3);
        b.record_duplex_family(1, 1);
        b.record_umi("AAAA", 3, 0, true);
        b.record_umi("CCCC", 2, 0, true);

        a.merge(b);

        let metrics = a.family_size_metrics();
        let cs3 = metrics.iter().find(|m| m.family_size == 3).unwrap();
        assert_eq!(cs3.cs_count, 2);
        let ss2 = metrics.iter().find(|m| m.family_size == 2).unwrap();
        assert_eq!(ss2.ss_count, 1);
        let ss4 = metrics.iter().find(|m| m.family_size == 4).unwrap();
        assert_eq!(ss4.ss_count, 1);
        let ds1 = metrics.iter().find(|m| m.family_size == 1).unwrap();
        assert_eq!(ds1.ds_count, 2);

        let duplex_metrics = a.duplex_family_size_metrics();
        // record_duplex_family(2, 3) normalizes to (3, 2) since ab >= ba
        let d32 = duplex_metrics.iter().find(|m| m.ab_size == 3 && m.ba_size == 2).unwrap();
        assert_eq!(d32.count, 2);
        let d11 = duplex_metrics.iter().find(|m| m.ab_size == 1 && m.ba_size == 1).unwrap();
        assert_eq!(d11.count, 1);

        let umi_metrics = a.umi_metrics();
        assert_eq!(umi_metrics.len(), 2);
    }
}
