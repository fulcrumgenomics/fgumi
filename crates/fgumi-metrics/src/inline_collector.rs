//! Inline metrics collection with downsampling and interval filtering.
//!
//! Provides the [`InlineMetricsCollector`] that wraps simplex or duplex collectors
//! with deterministic downsampling, optional interval filtering, and position group
//! size tracking. This is the core abstraction for inline metrics collection.

use std::collections::HashMap;
use std::path::{Path, PathBuf};

use anyhow::Result;

use crate::downsampling::DOWNSAMPLING_FRACTIONS;
use crate::group::PositionGroupSizeMetrics;
use crate::intervals::{Interval, overlaps_intervals};
use crate::template_info::{TemplateInfo, TemplateMetadata, compute_template_metadata};
use crate::writer::write_metrics;

/// Trait for the underlying simplex/duplex collectors to support generic dispatch.
pub trait MetricsCollector: Default + Send {
    /// Record a single coordinate+strand group of downsampled templates.
    fn record_group(&mut self, templates: &[TemplateMetadata<'_>]);
    /// Merge another collector's counts into this one.
    fn merge(&mut self, other: Self);
    /// Write collector-specific metrics files with the given prefix.
    ///
    /// # Errors
    ///
    /// Returns an error if any metrics file cannot be written.
    fn write_metrics(&self, prefix: &Path) -> Result<()>;
}

/// Inline metrics collector that wraps a [`MetricsCollector`] with downsampling,
/// interval filtering, and position group size tracking.
///
/// One collector is maintained per downsampling fraction (20 total). Position group
/// sizes are always collected from the full (undownsampled) data.
pub struct InlineMetricsCollector<C: MetricsCollector> {
    /// One collector per downsampling fraction (20 fractions = 20 collectors).
    collectors: Vec<C>,
    /// Optional interval filter (empty = no filtering).
    intervals: Vec<Interval>,
    /// Position group size counter: maps `position_group_size` -> count.
    /// Not downsampled -- always collected from the full data.
    position_group_sizes: HashMap<usize, u64>,
}

/// Builds an output path by appending a dot-separated suffix to the prefix.
///
/// For example, `with_extension(Path::new("/out/prefix"), "position_group_sizes.txt")`
/// produces `/out/prefix.position_group_sizes.txt`.
fn with_extension(prefix: &Path, suffix: &str) -> PathBuf {
    let mut path = prefix.as_os_str().to_owned();
    path.push(".");
    path.push(suffix);
    PathBuf::from(path)
}

impl<C: MetricsCollector> InlineMetricsCollector<C> {
    /// Creates a new inline metrics collector with the given interval filter.
    ///
    /// Creates one collector per downsampling fraction (20 total). If `intervals`
    /// is empty, no interval filtering is applied.
    #[must_use]
    pub fn new(intervals: Vec<Interval>) -> Self {
        let collectors = (0..DOWNSAMPLING_FRACTIONS.len()).map(|_| C::default()).collect();
        Self { collectors, intervals, position_group_sizes: HashMap::new() }
    }

    /// Records a coordinate group of templates for metrics collection.
    ///
    /// This method:
    /// 1. Skips empty template slices
    /// 2. Applies interval filtering if configured
    /// 3. Sub-groups templates by [`ReadInfoKey`] (coordinate+strand geometry)
    /// 4. For each sub-group: counts distinct MI families for position group size
    /// 5. Pre-computes [`TemplateMetadata`] and records per downsampling fraction
    pub fn record_coordinate_group(&mut self, templates: &[TemplateInfo]) {
        if templates.is_empty() {
            return;
        }

        // Apply interval filtering
        let filtered: Vec<&TemplateInfo> = if self.intervals.is_empty() {
            templates.iter().collect()
        } else {
            templates
                .iter()
                .filter(|t| {
                    overlaps_intervals(
                        t.ref_name.as_deref(),
                        t.position,
                        t.end_position,
                        &self.intervals,
                    )
                })
                .collect()
        };

        if filtered.is_empty() {
            return;
        }

        // Sub-group filtered templates by ReadInfoKey so that templates at the same
        // coordinate position but with different strand/mate geometry are counted as
        // separate position groups (matching the standalone metrics path).
        let mut sub_groups: HashMap<&crate::template_info::ReadInfoKey, Vec<&TemplateInfo>> =
            HashMap::new();
        for t in &filtered {
            sub_groups.entry(&t.read_info_key).or_default().push(t);
        }

        for sub_group in sub_groups.values() {
            // Count distinct MI families for position group size (from full filtered data)
            let position_group_size = {
                let mut mis: Vec<&str> = sub_group.iter().map(|t| t.mi.as_str()).collect();
                mis.sort_unstable();
                mis.dedup();
                mis.len()
            };
            *self.position_group_sizes.entry(position_group_size).or_insert(0) += 1;

            // Clone sub-group into an owned Vec for compute_template_metadata.
            let sub_templates: Vec<TemplateInfo> = sub_group.iter().copied().cloned().collect();
            let metadata = compute_template_metadata(&sub_templates);

            // For each downsampling fraction, filter by hash_fraction and record
            for (i, &fraction) in DOWNSAMPLING_FRACTIONS.iter().enumerate() {
                let downsampled: Vec<TemplateMetadata<'_>> = metadata
                    .iter()
                    .filter(|m| m.template.hash_fraction <= fraction)
                    .map(|m| TemplateMetadata {
                        template: m.template,
                        base_umi: m.base_umi,
                        is_a_strand: m.is_a_strand,
                        is_b_strand: m.is_b_strand,
                    })
                    .collect();

                if !downsampled.is_empty() {
                    self.collectors[i].record_group(&downsampled);
                }
            }
        }
    }

    /// Merges another collector into this one.
    ///
    /// Combines position group size counts and merges each per-fraction collector.
    pub fn merge(&mut self, other: Self) {
        for (size, count) in other.position_group_sizes {
            *self.position_group_sizes.entry(size).or_insert(0) += count;
        }
        for (mine, theirs) in self.collectors.iter_mut().zip(other.collectors) {
            mine.merge(theirs);
        }
    }

    /// Writes all metrics files with the given prefix.
    ///
    /// Writes position group size distribution and delegates to the last collector
    /// (100% fraction) for collector-specific metrics.
    ///
    /// # Errors
    ///
    /// Returns an error if any metrics file cannot be written.
    pub fn write_metrics(&self, prefix: &Path) -> Result<()> {
        // Write position group sizes
        let pgs_path = with_extension(prefix, "position_group_sizes.txt");
        let pgs_metrics = PositionGroupSizeMetrics::from_size_counts(
            self.position_group_sizes.iter().map(|(&size, &count)| (size, count)),
        );
        write_metrics(&pgs_path, &pgs_metrics, "position group size")?;

        // Delegate to the last collector (100% fraction)
        self.collectors[DOWNSAMPLING_FRACTIONS.len() - 1].write_metrics(prefix)?;

        Ok(())
    }

    /// Returns a reference to the per-fraction collectors.
    #[must_use]
    pub fn collectors(&self) -> &[C] {
        &self.collectors
    }

    /// Returns a reference to the position group size counts.
    #[must_use]
    pub fn position_group_sizes(&self) -> &HashMap<usize, u64> {
        &self.position_group_sizes
    }
}

/// Type-erased inline metrics collector (simplex or duplex).
///
/// This enum wraps the generic `InlineMetricsCollector<C>` to allow both simplex
/// and duplex collectors to be used interchangeably without generics at call sites.
pub enum InlineCollector {
    /// Simplex consensus metrics collector.
    Simplex(InlineMetricsCollector<crate::simplex::SimplexMetricsCollector>),
    /// Duplex consensus metrics collector.
    Duplex(InlineMetricsCollector<crate::duplex::DuplexMetricsCollector>),
}

impl InlineCollector {
    /// Record a coordinate group of templates for metrics collection.
    pub fn record_coordinate_group(&mut self, templates: &[TemplateInfo]) {
        match self {
            Self::Simplex(c) => c.record_coordinate_group(templates),
            Self::Duplex(c) => c.record_coordinate_group(templates),
        }
    }

    /// Merge another collector into this one.
    ///
    /// # Panics
    ///
    /// Panics if `self` and `other` are different variants.
    pub fn merge(&mut self, other: Self) {
        match (self, other) {
            (Self::Simplex(a), Self::Simplex(b)) => a.merge(b),
            (Self::Duplex(a), Self::Duplex(b)) => a.merge(b),
            _ => panic!("cannot merge InlineCollector variants of different types"),
        }
    }

    /// Write all metrics files with the given prefix.
    ///
    /// # Errors
    ///
    /// Returns an error if any metrics file cannot be written.
    pub fn write_metrics(&self, prefix: &std::path::Path) -> Result<()> {
        match self {
            Self::Simplex(c) => c.write_metrics(prefix),
            Self::Duplex(c) => c.write_metrics(prefix),
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[derive(Default)]
    struct MockCollector {
        groups_recorded: usize,
    }

    impl MetricsCollector for MockCollector {
        fn record_group(&mut self, templates: &[TemplateMetadata<'_>]) {
            if !templates.is_empty() {
                self.groups_recorded += 1;
            }
        }

        fn merge(&mut self, other: Self) {
            self.groups_recorded += other.groups_recorded;
        }

        fn write_metrics(&self, _prefix: &Path) -> Result<()> {
            Ok(())
        }
    }

    fn make_template(mi: &str, hash_fraction: f64) -> TemplateInfo {
        use crate::template_info::ReadInfoKey;
        TemplateInfo {
            mi: mi.to_string(),
            rx: "AAA-TTT".to_string(),
            ref_name: Some("chr1".to_string()),
            position: Some(100),
            end_position: Some(200),
            hash_fraction,
            read_info_key: ReadInfoKey::default(),
        }
    }

    #[test]
    fn test_new_creates_correct_number_of_collectors() {
        let collector: InlineMetricsCollector<MockCollector> =
            InlineMetricsCollector::new(Vec::new());
        assert_eq!(collector.collectors().len(), 20);
    }

    #[test]
    fn test_record_coordinate_group_tracks_position_group_size() {
        let mut collector: InlineMetricsCollector<MockCollector> =
            InlineMetricsCollector::new(Vec::new());

        // Three templates with three distinct MI values
        let templates =
            vec![make_template("1/A", 0.5), make_template("1/B", 0.5), make_template("2/A", 0.5)];
        collector.record_coordinate_group(&templates);

        // Three distinct MI values: "1/A", "1/B", "2/A"
        assert_eq!(collector.position_group_sizes().get(&3), Some(&1));
    }

    #[test]
    fn test_record_coordinate_group_empty_skipped() {
        let mut collector: InlineMetricsCollector<MockCollector> =
            InlineMetricsCollector::new(Vec::new());

        collector.record_coordinate_group(&[]);

        assert!(collector.position_group_sizes().is_empty());
        // No groups should have been recorded on any collector
        for c in collector.collectors() {
            assert_eq!(c.groups_recorded, 0);
        }
    }

    #[test]
    fn test_interval_filtering() {
        let intervals = vec![Interval { ref_name: "chr1".to_string(), start: 50, end: 150 }];
        let mut collector: InlineMetricsCollector<MockCollector> =
            InlineMetricsCollector::new(intervals);

        // Template that overlaps the interval (chr1:100-200 overlaps chr1:50-150)
        let overlapping = TemplateInfo {
            mi: "1/A".to_string(),
            rx: "AAA-TTT".to_string(),
            ref_name: Some("chr1".to_string()),
            position: Some(100),
            end_position: Some(200),
            hash_fraction: 0.5,
            read_info_key: crate::template_info::ReadInfoKey::default(),
        };

        // Template that does NOT overlap (chr2 != chr1)
        let non_overlapping = TemplateInfo {
            mi: "2/A".to_string(),
            rx: "AAA-TTT".to_string(),
            ref_name: Some("chr2".to_string()),
            position: Some(100),
            end_position: Some(200),
            hash_fraction: 0.5,
            read_info_key: crate::template_info::ReadInfoKey::default(),
        };

        collector.record_coordinate_group(&[overlapping, non_overlapping]);

        // Only one template passes filtering, so position group size = 1
        assert_eq!(collector.position_group_sizes().get(&1), Some(&1));
        // The last collector (100%, index 19) should have recorded one group
        assert_eq!(collector.collectors()[19].groups_recorded, 1);
    }

    #[test]
    fn test_interval_filtering_all_filtered() {
        let intervals = vec![Interval { ref_name: "chrX".to_string(), start: 0, end: 100 }];
        let mut collector: InlineMetricsCollector<MockCollector> =
            InlineMetricsCollector::new(intervals);

        // Neither template overlaps chrX:0-100
        let templates = vec![make_template("1/A", 0.5), make_template("2/A", 0.5)];
        collector.record_coordinate_group(&templates);

        // Everything filtered out, so no position group sizes recorded
        assert!(collector.position_group_sizes().is_empty());
        for c in collector.collectors() {
            assert_eq!(c.groups_recorded, 0);
        }
    }

    #[test]
    fn test_merge_combines_collectors() {
        let mut collector1: InlineMetricsCollector<MockCollector> =
            InlineMetricsCollector::new(Vec::new());
        let mut collector2: InlineMetricsCollector<MockCollector> =
            InlineMetricsCollector::new(Vec::new());

        collector1.record_coordinate_group(&[make_template("1/A", 0.5)]);
        collector2.record_coordinate_group(&[make_template("2/A", 0.5)]);

        collector1.merge(collector2);

        // Both recorded position group size 1, so count should be 2
        assert_eq!(collector1.position_group_sizes().get(&1), Some(&2));
        // Last collector should have merged groups from both
        assert_eq!(collector1.collectors()[19].groups_recorded, 2);
    }

    #[test]
    fn test_downsampling_filters_by_hash_fraction() {
        let mut collector: InlineMetricsCollector<MockCollector> =
            InlineMetricsCollector::new(Vec::new());

        // Template with very high hash_fraction (0.99) — only in the 100% bucket
        let templates = vec![make_template("1/A", 0.99)];
        collector.record_coordinate_group(&templates);

        // First fraction is 0.05: hash_fraction 0.99 > 0.05, so no templates pass
        assert_eq!(collector.collectors()[0].groups_recorded, 0);

        // Last fraction is 1.00: hash_fraction 0.99 <= 1.00, so template passes
        assert_eq!(collector.collectors()[19].groups_recorded, 1);

        // Fraction index 18 is 0.95: hash_fraction 0.99 > 0.95, so no templates pass
        assert_eq!(collector.collectors()[18].groups_recorded, 0);
    }

    #[test]
    fn test_with_extension() {
        let path = with_extension(Path::new("/out/prefix"), "position_group_sizes.txt");
        assert_eq!(path, PathBuf::from("/out/prefix.position_group_sizes.txt"));
    }

    fn make_template_with_key(
        mi: &str,
        hash_fraction: f64,
        key: crate::template_info::ReadInfoKey,
    ) -> TemplateInfo {
        TemplateInfo {
            mi: mi.to_string(),
            rx: "AAA-TTT".to_string(),
            ref_name: Some("chr1".to_string()),
            position: Some(100),
            end_position: Some(200),
            hash_fraction,
            read_info_key: key,
        }
    }

    #[test]
    fn test_record_coordinate_group_sub_groups_by_read_info_key() {
        use crate::template_info::ReadInfoKey;

        let mut collector: InlineMetricsCollector<MockCollector> =
            InlineMetricsCollector::new(Vec::new());

        // Two templates at the same position but with different ReadInfoKeys (different strands).
        // These should be counted as two separate position groups, each of size 1.
        let key_fwd = ReadInfoKey {
            ref_index1: Some(0),
            start1: 100,
            strand1: false,
            ref_index2: Some(0),
            start2: 200,
            strand2: true,
        };
        let key_rev = ReadInfoKey {
            ref_index1: Some(0),
            start1: 100,
            strand1: true,
            ref_index2: Some(0),
            start2: 200,
            strand2: false,
        };

        let templates = vec![
            make_template_with_key("1/A", 0.5, key_fwd),
            make_template_with_key("2/A", 0.5, key_rev),
        ];
        collector.record_coordinate_group(&templates);

        // With strand-aware sub-grouping, we should see two groups of size 1 (not one group of size 2).
        let pgs = collector.position_group_sizes();
        assert_eq!(
            pgs.get(&1),
            Some(&2),
            "Should have two position groups of size 1, got: {pgs:?}"
        );
        assert_eq!(pgs.get(&2), None, "Should NOT have a position group of size 2, got: {pgs:?}");
    }

    #[test]
    fn test_record_coordinate_group_same_key_stays_grouped() {
        use crate::template_info::ReadInfoKey;

        let mut collector: InlineMetricsCollector<MockCollector> =
            InlineMetricsCollector::new(Vec::new());

        // Two templates with the SAME ReadInfoKey should remain in one position group of size 2.
        let key = ReadInfoKey {
            ref_index1: Some(0),
            start1: 100,
            strand1: false,
            ref_index2: Some(0),
            start2: 200,
            strand2: true,
        };

        let templates = vec![
            make_template_with_key("1/A", 0.5, key.clone()),
            make_template_with_key("2/A", 0.5, key),
        ];
        collector.record_coordinate_group(&templates);

        let pgs = collector.position_group_sizes();
        assert_eq!(pgs.get(&2), Some(&1), "Should have one position group of size 2, got: {pgs:?}");
    }

    #[test]
    fn test_write_metrics_creates_files() {
        let dir = tempfile::tempdir().unwrap();
        let prefix = dir.path().join("test_prefix");

        let mut collector: InlineMetricsCollector<MockCollector> =
            InlineMetricsCollector::new(Vec::new());
        collector.record_coordinate_group(&[make_template("1/A", 0.5)]);

        collector.write_metrics(&prefix).unwrap();

        let pgs_path = with_extension(&prefix, "position_group_sizes.txt");
        assert!(pgs_path.exists(), "position_group_sizes.txt should exist");
    }
}
