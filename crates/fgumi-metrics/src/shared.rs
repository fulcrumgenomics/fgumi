//! Shared UMI tracking types used across simplex and duplex metrics.
//!
//! This module provides [`UmiCountTracker`] for accumulating raw, error, and unique
//! observation counts per UMI sequence, and [`UmiMetric`] for the corresponding
//! serializable metric output.

use std::collections::HashMap;

use serde::{Deserialize, Serialize};

use crate::{Metric, frac};

/// Metrics describing observed UMI sequences and their observation frequencies.
///
/// UMI sequences may be corrected using information within a double-stranded tag family.
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct UmiMetric {
    /// The UMI sequence (possibly corrected)
    pub umi: String,
    /// Number of read pairs observing this UMI (after correction)
    pub raw_observations: usize,
    /// Subset of raw observations that underwent correction
    pub raw_observations_with_errors: usize,
    /// Number of double-stranded tag families observing this UMI
    pub unique_observations: usize,
    /// Fraction of all raw observations
    pub fraction_raw_observations: f64,
    /// Fraction of all unique observations
    pub fraction_unique_observations: f64,
}

impl UmiMetric {
    /// Creates a new UMI metric
    #[must_use]
    pub fn new(umi: String) -> Self {
        Self {
            umi,
            raw_observations: 0,
            raw_observations_with_errors: 0,
            unique_observations: 0,
            fraction_raw_observations: 0.0,
            fraction_unique_observations: 0.0,
        }
    }
}

impl Default for UmiMetric {
    fn default() -> Self {
        Self::new(String::new())
    }
}

impl Metric for UmiMetric {
    fn metric_name() -> &'static str {
        "UMI"
    }
}

/// Tracks raw, error, and unique UMI observation counts.
///
/// Each UMI string maps to a tuple of `(raw_count, error_count, unique_count)`.
/// Use [`record`](Self::record) to accumulate observations and [`to_metrics`](Self::to_metrics)
/// to produce sorted [`UmiMetric`] output.
pub struct UmiCountTracker {
    /// Maps UMI string to `(raw_count, error_count, unique_count)`.
    counts: HashMap<String, (usize, usize, usize)>,
}

impl UmiCountTracker {
    /// Creates an empty tracker.
    #[must_use]
    pub fn new() -> Self {
        Self { counts: HashMap::new() }
    }

    /// Records an observation for a UMI.
    pub fn record(&mut self, umi: &str, raw_count: usize, error_count: usize, is_unique: bool) {
        if let Some(entry) = self.counts.get_mut(umi) {
            entry.0 += raw_count;
            entry.1 += error_count;
            if is_unique {
                entry.2 += 1;
            }
        } else {
            self.counts.insert(umi.to_string(), (raw_count, error_count, usize::from(is_unique)));
        }
    }

    /// Total raw observations across all UMIs.
    #[must_use]
    pub fn total_raw(&self) -> usize {
        self.counts.values().map(|(raw, _, _)| raw).sum()
    }

    /// Total unique observations across all UMIs.
    #[must_use]
    pub fn total_unique(&self) -> usize {
        self.counts.values().map(|(_, _, unique)| unique).sum()
    }

    /// Iterates over all tracked UMIs, yielding `(umi, raw_count, error_count, unique_count)`.
    pub(crate) fn iter(&self) -> impl Iterator<Item = (&str, usize, usize, usize)> {
        self.counts.iter().map(|(umi, &(raw, errors, unique))| (umi.as_str(), raw, errors, unique))
    }

    /// Generates [`UmiMetric`] entries sorted alphabetically by UMI sequence.
    ///
    /// Computes fractional observations relative to the totals tracked by this instance.
    #[must_use]
    pub fn to_metrics(&self) -> Vec<UmiMetric> {
        let total_raw = self.total_raw();
        let total_unique = self.total_unique();

        let mut metrics: Vec<_> = self
            .iter()
            .map(|(umi, raw, errors, unique)| UmiMetric {
                umi: umi.to_string(),
                raw_observations: raw,
                raw_observations_with_errors: errors,
                unique_observations: unique,
                fraction_raw_observations: frac(raw, total_raw),
                fraction_unique_observations: frac(unique, total_unique),
            })
            .collect();

        // Sort by UMI string (matching fgbio's sortBy(_.umi))
        metrics.sort_by(|a, b| a.umi.cmp(&b.umi));
        metrics
    }
}

impl Default for UmiCountTracker {
    fn default() -> Self {
        Self::new()
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    // =========================================================================
    // UmiCountTracker tests
    // =========================================================================

    #[test]
    fn test_umi_count_tracker_empty() {
        let tracker = UmiCountTracker::new();
        assert_eq!(tracker.total_raw(), 0);
        assert_eq!(tracker.total_unique(), 0);
        assert_eq!(tracker.iter().count(), 0);
    }

    #[test]
    fn test_umi_count_tracker_default() {
        let tracker = UmiCountTracker::default();
        assert_eq!(tracker.total_raw(), 0);
        assert_eq!(tracker.total_unique(), 0);
        assert_eq!(tracker.iter().count(), 0);
    }

    #[test]
    fn test_umi_count_tracker_record_and_iter() {
        let mut tracker = UmiCountTracker::new();
        tracker.record("AAAA", 10, 2, true);
        tracker.record("AAAA", 5, 1, false);
        tracker.record("CCCC", 8, 0, true);

        assert_eq!(tracker.total_raw(), 23); // 15 + 8
        assert_eq!(tracker.total_unique(), 2); // 1 + 1

        let mut items: Vec<_> = tracker.iter().collect();
        items.sort_by(|a, b| a.0.cmp(b.0));
        assert_eq!(items.len(), 2);
        assert_eq!(items[0], ("AAAA", 15, 3, 1));
        assert_eq!(items[1], ("CCCC", 8, 0, 1));
    }

    // =========================================================================
    // UmiMetric tests
    // =========================================================================

    #[test]
    fn test_umi_metric_new() {
        let metric = UmiMetric::new("ACGT".to_string());
        assert_eq!(metric.umi, "ACGT");
        assert_eq!(metric.raw_observations, 0);
        assert_eq!(metric.unique_observations, 0);
    }

    // =========================================================================
    // UmiCountTracker::to_metrics tests
    // =========================================================================

    #[test]
    fn test_to_metrics_sorting() {
        let mut tracker = UmiCountTracker::new();
        tracker.record("ZZZZ", 1, 0, true);
        tracker.record("AAAA", 1, 0, true);
        tracker.record("MMMM", 1, 0, true);

        let metrics = tracker.to_metrics();
        // Should be sorted alphabetically
        assert_eq!(metrics[0].umi, "AAAA");
        assert_eq!(metrics[1].umi, "MMMM");
        assert_eq!(metrics[2].umi, "ZZZZ");
    }

    #[test]
    fn test_to_metrics_fractions() {
        let mut tracker = UmiCountTracker::new();
        tracker.record("AAAA", 10, 2, true);
        tracker.record("AAAA", 5, 1, false);
        tracker.record("CCCC", 8, 0, true);

        let metrics = tracker.to_metrics();
        assert_eq!(metrics.len(), 2);

        let aaaa =
            metrics.iter().find(|m| m.umi == "AAAA").expect("AAAA UMI metric should be present");
        assert_eq!(aaaa.raw_observations, 15);
        assert_eq!(aaaa.raw_observations_with_errors, 3);
        assert_eq!(aaaa.unique_observations, 1);
        assert!((aaaa.fraction_raw_observations - 15.0 / 23.0).abs() < f64::EPSILON);
        assert!((aaaa.fraction_unique_observations - 0.5).abs() < f64::EPSILON);

        let cccc =
            metrics.iter().find(|m| m.umi == "CCCC").expect("CCCC UMI metric should be present");
        assert_eq!(cccc.raw_observations, 8);
        assert_eq!(cccc.unique_observations, 1);
        assert!((cccc.fraction_raw_observations - 8.0 / 23.0).abs() < f64::EPSILON);
        assert!((cccc.fraction_unique_observations - 0.5).abs() < f64::EPSILON);
    }

    #[test]
    fn test_to_metrics_empty() {
        let tracker = UmiCountTracker::new();
        let metrics = tracker.to_metrics();
        assert!(metrics.is_empty());
    }
}
