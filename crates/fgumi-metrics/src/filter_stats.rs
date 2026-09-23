//! Filtering statistics metric emitted by `fgumi filter --stats`.

use serde::{Deserialize, Serialize};

use crate::Metric;

/// Summary statistics for a `fgumi filter` run: one row describing how many
/// reads were seen, kept, and rejected, plus the pass rate.
#[derive(Debug, Clone, PartialEq, Serialize, Deserialize, Default)]
pub struct FilterStatsMetrics {
    /// Total reads examined.
    pub total_reads: u64,

    /// Reads that passed all filters.
    pub passed_reads: u64,

    /// Reads that were rejected by at least one filter.
    pub failed_reads: u64,

    /// Fraction of `total_reads` that passed (`0.0` when `total_reads == 0`).
    #[serde(with = "crate::float")]
    pub pass_rate: f64,
}

impl FilterStatsMetrics {
    /// Builds the metric from raw counts, computing the pass rate safely.
    #[must_use]
    pub fn from_counts(total: u64, passed: u64, failed: u64) -> Self {
        #[allow(clippy::cast_precision_loss)]
        let pass_rate = if total > 0 { passed as f64 / total as f64 } else { 0.0 };
        Self { total_reads: total, passed_reads: passed, failed_reads: failed, pass_rate }
    }
}

impl Metric for FilterStatsMetrics {
    fn metric_name() -> &'static str {
        "filter stats"
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::writer::{read_metrics_auto, write_metrics_auto};
    use tempfile::NamedTempFile;

    #[test]
    fn from_counts_computes_pass_rate() {
        let m = FilterStatsMetrics::from_counts(10, 7, 3);
        assert_eq!(m.total_reads, 10);
        assert_eq!(m.passed_reads, 7);
        assert_eq!(m.failed_reads, 3);
        assert!((m.pass_rate - 0.7).abs() < 1e-12);
    }

    #[test]
    fn from_counts_zero_total_is_zero_pass_rate() {
        let m = FilterStatsMetrics::from_counts(0, 0, 0);
        assert!(m.pass_rate.abs() < f64::EPSILON);
        assert!(m.pass_rate.is_finite());
    }

    #[test]
    fn header_is_declared_column_order() {
        let tmp = NamedTempFile::new().unwrap();
        write_metrics_auto(tmp.path(), &[FilterStatsMetrics::from_counts(4, 1, 3)]).unwrap();
        let content = std::fs::read_to_string(tmp.path()).unwrap();
        let header = content.lines().next().unwrap();
        assert_eq!(header, "total_reads\tpassed_reads\tfailed_reads\tpass_rate");
    }

    #[test]
    fn round_trips_through_metric_tsv() {
        let tmp = NamedTempFile::new().unwrap();
        let original = FilterStatsMetrics::from_counts(1000, 950, 50);
        write_metrics_auto(tmp.path(), std::slice::from_ref(&original)).unwrap();
        let back: Vec<FilterStatsMetrics> = read_metrics_auto(tmp.path()).unwrap();
        assert_eq!(back, vec![original]);
    }

    #[test]
    fn empty_write_is_header_only_and_reparses() {
        let tmp = NamedTempFile::new().unwrap();
        let empty: Vec<FilterStatsMetrics> = vec![];
        write_metrics_auto(tmp.path(), &empty).unwrap();
        let content = std::fs::read_to_string(tmp.path()).unwrap();
        assert_eq!(content.lines().count(), 1);
        let back: Vec<FilterStatsMetrics> = read_metrics_auto(tmp.path()).unwrap();
        assert!(back.is_empty());
    }
}
