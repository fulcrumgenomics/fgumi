//! Metrics for paired UMI half-match detection.
//!
//! This module provides metrics types for tracking when UMI families share only
//! one half of their paired UMI (left or right adapter) but not the full UMI.
//! This helps identify systematic issues in paired UMI protocols.

use serde::{Deserialize, Serialize};

use super::Metric;

/// Output row for half-match metrics TSV.
///
/// Each row represents a single metric with values for left-half matches,
/// right-half matches, and combined (both) counts.
#[derive(Debug, Clone, Serialize, Deserialize, Default)]
pub struct PairedHalfMatchMetric {
    /// Name of the metric
    pub metric: String,
    /// Value for left-half UMI matches
    pub left_half: f64,
    /// Value for right-half UMI matches
    pub right_half: f64,
    /// Combined value (both left and right)
    pub both: f64,
    /// Human-readable description of the metric
    pub description: String,
}

impl Metric for PairedHalfMatchMetric {
    fn metric_name() -> &'static str {
        "paired half-match"
    }
}

/// Thread-local collector for half-match statistics.
///
/// Collects statistics about UMI families that share one half of their
/// paired UMI but not both halves. This can indicate systematic issues
/// in paired UMI protocols.
#[derive(Debug, Clone, Default)]
pub struct HalfMatchCollector {
    /// Number of position groups analyzed (groups with 2+ families)
    pub position_groups_analyzed: u64,
    /// Total families in multi-family groups
    pub families_analyzed: u64,
    /// Families with a left-half-only match to another family
    pub families_with_left_half_match: u64,
    /// Families with a right-half-only match to another family
    pub families_with_right_half_match: u64,
    /// Families with any half-match (left or right, deduplicated)
    pub families_with_any_half_match: u64,
    /// Number of family pairs sharing only left half
    pub family_pairs_left_collision: u64,
    /// Number of family pairs sharing only right half
    pub family_pairs_right_collision: u64,
    /// Position groups where ALL families share the same left half
    pub position_groups_all_same_left: u64,
    /// Position groups where ALL families share the same right half
    pub position_groups_all_same_right: u64,
    /// Count of families involved in left half-matches from forward strand only (`PairedA` only)
    pub half_matches_forward_strand_left: u64,
    /// Count of families involved in right half-matches from forward strand only (`PairedA` only)
    pub half_matches_forward_strand_right: u64,
    /// Count of families involved in left half-matches from reverse strand only (`PairedB` only)
    pub half_matches_reverse_strand_left: u64,
    /// Count of families involved in right half-matches from reverse strand only (`PairedB` only)
    pub half_matches_reverse_strand_right: u64,
    /// Count of families involved in left half-matches with both strands (duplex)
    pub half_matches_duplex_strand_left: u64,
    /// Count of families involved in right half-matches with both strands (duplex)
    pub half_matches_duplex_strand_right: u64,
}

impl HalfMatchCollector {
    /// Merge another collector into this one.
    pub fn merge(&mut self, other: &Self) {
        self.position_groups_analyzed += other.position_groups_analyzed;
        self.families_analyzed += other.families_analyzed;
        self.families_with_left_half_match += other.families_with_left_half_match;
        self.families_with_right_half_match += other.families_with_right_half_match;
        self.families_with_any_half_match += other.families_with_any_half_match;
        self.family_pairs_left_collision += other.family_pairs_left_collision;
        self.family_pairs_right_collision += other.family_pairs_right_collision;
        self.position_groups_all_same_left += other.position_groups_all_same_left;
        self.position_groups_all_same_right += other.position_groups_all_same_right;
        self.half_matches_forward_strand_left += other.half_matches_forward_strand_left;
        self.half_matches_forward_strand_right += other.half_matches_forward_strand_right;
        self.half_matches_reverse_strand_left += other.half_matches_reverse_strand_left;
        self.half_matches_reverse_strand_right += other.half_matches_reverse_strand_right;
        self.half_matches_duplex_strand_left += other.half_matches_duplex_strand_left;
        self.half_matches_duplex_strand_right += other.half_matches_duplex_strand_right;
    }

    /// Safe division that returns 0.0 when denominator is 0.
    #[inline]
    fn safe_divide(numerator: u64, denominator: u64) -> f64 {
        if denominator == 0 {
            0.0
        } else {
            #[allow(clippy::cast_precision_loss)]
            {
                numerator as f64 / denominator as f64
            }
        }
    }

    /// Convert collected statistics to TSV output rows.
    #[must_use]
    #[allow(clippy::cast_precision_loss)]
    pub fn to_metrics(&self) -> Vec<PairedHalfMatchMetric> {
        // Use the deduplicated count for "both" column to avoid double-counting
        // families that have both left and right half-matches
        let families_with_half_match_both = self.families_with_any_half_match;
        let family_pairs_collision_both =
            self.family_pairs_left_collision + self.family_pairs_right_collision;
        let position_groups_all_same_both =
            self.position_groups_all_same_left + self.position_groups_all_same_right;
        let half_matches_forward_both =
            self.half_matches_forward_strand_left + self.half_matches_forward_strand_right;
        let half_matches_reverse_both =
            self.half_matches_reverse_strand_left + self.half_matches_reverse_strand_right;
        let half_matches_duplex_both =
            self.half_matches_duplex_strand_left + self.half_matches_duplex_strand_right;

        vec![
            PairedHalfMatchMetric {
                metric: "position_groups_analyzed".to_string(),
                left_half: self.position_groups_analyzed as f64,
                right_half: self.position_groups_analyzed as f64,
                both: self.position_groups_analyzed as f64,
                description: "Position groups with 2+ families".to_string(),
            },
            PairedHalfMatchMetric {
                metric: "families_analyzed".to_string(),
                left_half: self.families_analyzed as f64,
                right_half: self.families_analyzed as f64,
                both: self.families_analyzed as f64,
                description: "Total families in multi-family groups".to_string(),
            },
            PairedHalfMatchMetric {
                metric: "families_with_half_match_only".to_string(),
                left_half: self.families_with_left_half_match as f64,
                right_half: self.families_with_right_half_match as f64,
                both: families_with_half_match_both as f64,
                description: "Families sharing half-UMI but not full UMI".to_string(),
            },
            PairedHalfMatchMetric {
                metric: "fraction_families_half_match".to_string(),
                left_half: Self::safe_divide(
                    self.families_with_left_half_match,
                    self.families_analyzed,
                ),
                right_half: Self::safe_divide(
                    self.families_with_right_half_match,
                    self.families_analyzed,
                ),
                both: Self::safe_divide(families_with_half_match_both, self.families_analyzed),
                description: "Fraction of families with half-match".to_string(),
            },
            PairedHalfMatchMetric {
                metric: "family_pairs_with_collision".to_string(),
                left_half: self.family_pairs_left_collision as f64,
                right_half: self.family_pairs_right_collision as f64,
                both: family_pairs_collision_both as f64,
                description: "Family pairs sharing only one half-UMI".to_string(),
            },
            PairedHalfMatchMetric {
                metric: "position_groups_all_same_half".to_string(),
                left_half: self.position_groups_all_same_left as f64,
                right_half: self.position_groups_all_same_right as f64,
                both: position_groups_all_same_both as f64,
                description: "Groups where ALL families share same half (potential adapter issue)"
                    .to_string(),
            },
            PairedHalfMatchMetric {
                metric: "fraction_groups_all_same_half".to_string(),
                left_half: Self::safe_divide(
                    self.position_groups_all_same_left,
                    self.position_groups_analyzed,
                ),
                right_half: Self::safe_divide(
                    self.position_groups_all_same_right,
                    self.position_groups_analyzed,
                ),
                both: Self::safe_divide(
                    position_groups_all_same_both,
                    self.position_groups_analyzed,
                ),
                description: "Fraction of groups with systematic pattern".to_string(),
            },
            PairedHalfMatchMetric {
                metric: "half_matches_forward_strand".to_string(),
                left_half: self.half_matches_forward_strand_left as f64,
                right_half: self.half_matches_forward_strand_right as f64,
                both: half_matches_forward_both as f64,
                description: "Families in half-match pairs from forward strand only".to_string(),
            },
            PairedHalfMatchMetric {
                metric: "half_matches_reverse_strand".to_string(),
                left_half: self.half_matches_reverse_strand_left as f64,
                right_half: self.half_matches_reverse_strand_right as f64,
                both: half_matches_reverse_both as f64,
                description: "Families in half-match pairs from reverse strand only".to_string(),
            },
            PairedHalfMatchMetric {
                metric: "half_matches_duplex_strand".to_string(),
                left_half: self.half_matches_duplex_strand_left as f64,
                right_half: self.half_matches_duplex_strand_right as f64,
                both: half_matches_duplex_both as f64,
                description: "Families in half-match pairs with both strands (duplex)".to_string(),
            },
        ]
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_half_match_collector_default() {
        let collector = HalfMatchCollector::default();
        assert_eq!(collector.position_groups_analyzed, 0);
        assert_eq!(collector.families_analyzed, 0);
        assert_eq!(collector.families_with_left_half_match, 0);
        assert_eq!(collector.families_with_right_half_match, 0);
    }

    #[test]
    fn test_half_match_collector_merge() {
        let mut c1 = HalfMatchCollector {
            position_groups_analyzed: 10,
            families_analyzed: 50,
            families_with_left_half_match: 5,
            families_with_right_half_match: 3,
            families_with_any_half_match: 6, // Some overlap
            family_pairs_left_collision: 6,
            family_pairs_right_collision: 4,
            position_groups_all_same_left: 1,
            position_groups_all_same_right: 2,
            half_matches_forward_strand_left: 10,
            half_matches_forward_strand_right: 8,
            half_matches_reverse_strand_left: 12,
            half_matches_reverse_strand_right: 9,
            half_matches_duplex_strand_left: 3,
            half_matches_duplex_strand_right: 2,
        };

        let c2 = HalfMatchCollector {
            position_groups_analyzed: 20,
            families_analyzed: 100,
            families_with_left_half_match: 10,
            families_with_right_half_match: 7,
            families_with_any_half_match: 14, // Some overlap
            family_pairs_left_collision: 12,
            family_pairs_right_collision: 8,
            position_groups_all_same_left: 3,
            position_groups_all_same_right: 4,
            half_matches_forward_strand_left: 20,
            half_matches_forward_strand_right: 16,
            half_matches_reverse_strand_left: 24,
            half_matches_reverse_strand_right: 18,
            half_matches_duplex_strand_left: 5,
            half_matches_duplex_strand_right: 4,
        };

        c1.merge(&c2);

        assert_eq!(c1.position_groups_analyzed, 30);
        assert_eq!(c1.families_analyzed, 150);
        assert_eq!(c1.families_with_left_half_match, 15);
        assert_eq!(c1.families_with_right_half_match, 10);
        assert_eq!(c1.families_with_any_half_match, 20);
        assert_eq!(c1.family_pairs_left_collision, 18);
        assert_eq!(c1.family_pairs_right_collision, 12);
        assert_eq!(c1.position_groups_all_same_left, 4);
        assert_eq!(c1.position_groups_all_same_right, 6);
        assert_eq!(c1.half_matches_forward_strand_left, 30);
        assert_eq!(c1.half_matches_forward_strand_right, 24);
        assert_eq!(c1.half_matches_reverse_strand_left, 36);
        assert_eq!(c1.half_matches_reverse_strand_right, 27);
        assert_eq!(c1.half_matches_duplex_strand_left, 8);
        assert_eq!(c1.half_matches_duplex_strand_right, 6);
    }

    #[test]
    fn test_to_metrics_format() {
        let collector = HalfMatchCollector {
            position_groups_analyzed: 100,
            families_analyzed: 500,
            families_with_left_half_match: 10,
            families_with_right_half_match: 8,
            families_with_any_half_match: 15, // Less than 10+8 due to overlap
            family_pairs_left_collision: 12,
            family_pairs_right_collision: 9,
            position_groups_all_same_left: 2,
            position_groups_all_same_right: 1,
            half_matches_forward_strand_left: 5,
            half_matches_forward_strand_right: 4,
            half_matches_reverse_strand_left: 6,
            half_matches_reverse_strand_right: 5,
            half_matches_duplex_strand_left: 2,
            half_matches_duplex_strand_right: 1,
        };

        let metrics = collector.to_metrics();

        assert_eq!(metrics.len(), 10);

        // Check first metric
        assert_eq!(metrics[0].metric, "position_groups_analyzed");
        assert!((metrics[0].left_half - 100.0).abs() < f64::EPSILON);
        assert!((metrics[0].both - 100.0).abs() < f64::EPSILON);

        // Check families_analyzed
        assert_eq!(metrics[1].metric, "families_analyzed");
        assert!((metrics[1].left_half - 500.0).abs() < f64::EPSILON);

        // Check families_with_half_match_only - "both" uses deduplicated count
        assert_eq!(metrics[2].metric, "families_with_half_match_only");
        assert!((metrics[2].left_half - 10.0).abs() < f64::EPSILON);
        assert!((metrics[2].right_half - 8.0).abs() < f64::EPSILON);
        assert!((metrics[2].both - 15.0).abs() < f64::EPSILON); // Uses families_with_any_half_match

        // Check fraction calculation - uses deduplicated count for "both"
        assert_eq!(metrics[3].metric, "fraction_families_half_match");
        assert!((metrics[3].left_half - 0.02).abs() < 1e-10);
        assert!((metrics[3].right_half - 0.016).abs() < 1e-10);
        assert!((metrics[3].both - 0.03).abs() < 1e-10); // 15/500 = 0.03

        // Check duplex strand metric
        assert_eq!(metrics[9].metric, "half_matches_duplex_strand");
        assert!((metrics[9].left_half - 2.0).abs() < f64::EPSILON);
        assert!((metrics[9].right_half - 1.0).abs() < f64::EPSILON);
        assert!((metrics[9].both - 3.0).abs() < f64::EPSILON);
    }

    #[test]
    fn test_safe_divide_zero() {
        let collector = HalfMatchCollector::default();
        let metrics = collector.to_metrics();

        // All fraction metrics should be 0.0 when denominator is 0
        assert_eq!(metrics[3].metric, "fraction_families_half_match");
        assert!((metrics[3].left_half - 0.0).abs() < f64::EPSILON);
        assert!((metrics[3].right_half - 0.0).abs() < f64::EPSILON);
        assert!((metrics[3].both - 0.0).abs() < f64::EPSILON);

        assert_eq!(metrics[6].metric, "fraction_groups_all_same_half");
        assert!((metrics[6].left_half - 0.0).abs() < f64::EPSILON);
    }

    #[test]
    fn test_metric_trait_impl() {
        assert_eq!(PairedHalfMatchMetric::metric_name(), "paired half-match");
    }
}
