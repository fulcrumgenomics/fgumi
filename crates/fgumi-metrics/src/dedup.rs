//! Metrics for the `dedup` command.

use serde::{Deserialize, Serialize};

use crate::Metric;

/// Metrics written by `fgumi dedup --metrics`.
///
/// Three sections, in column order:
///
/// 1. **Filtering** — `filtered_templates` through `passthrough_templates`.
///    `dedup` is a read *filter* as well as a duplicate marker: templates failing
///    the pre-grouping filter are dropped and never reach the output. The
///    `filtered_*` columns give the count per reason in **template** units — note
///    that `group`'s fgbio-shaped `.grouping_metrics.txt` counts primary
///    *records*, so the two files' discard counts are not comparable.
/// 2. **Deduplication** — `total_templates` through `duplicate_reads`, covering
///    only what survived filtering.
/// 3. **Diagnostics** — `secondary_reads`, `supplementary_reads`,
///    `missing_tc_tag`.
///
/// Templates read from the input are `filtered_templates + total_templates`; that
/// sum is not emitted as its own column because it is derivable, matching
/// [`crate::group::UmiGroupingMetrics`]' treatment of its derivable fields.
#[derive(Debug, Clone, Default, Serialize, Deserialize)]
pub struct DeduplicationMetrics {
    /// Templates dropped by the filter, for any reason.
    pub filtered_templates: u64,
    /// Dropped because a record was shorter than the minimum BAM record length
    /// (corrupt input).
    pub filtered_malformed_record: u64,
    /// Dropped because the template had neither a primary R1 nor a primary R2.
    pub filtered_no_primary_reads: u64,
    /// Dropped because no primary read was mapped. Always zero with
    /// `--include-unmapped`, which routes these to `passthrough_templates`.
    pub filtered_unmapped: u64,
    /// Dropped because a primary read carried the QC-fail flag.
    pub filtered_not_passing_filter: u64,
    /// Dropped because a mapped primary read's mapping quality was below
    /// `--min-map-q`.
    pub filtered_low_mapping_quality: u64,
    /// Dropped because a primary read's `MQ` tag was below `--min-map-q`.
    pub filtered_low_mate_mapping_quality: u64,
    /// Dropped because a primary read had no UMI tag.
    pub filtered_missing_umi: u64,
    /// Dropped because the UMI contained an `N` base.
    pub filtered_ns_in_umi: u64,
    /// Dropped because the UMI was shorter than `--min-umi-length`.
    pub filtered_umi_too_short: u64,
    /// Templates emitted untouched by `--include-unmapped`: never filtered, never
    /// marked, never MI-tagged. Counted as unique in the columns below.
    pub passthrough_templates: u64,
    /// Templates written to the output (accepted by the filter, plus
    /// pass-throughs) — not the number read from the input.
    pub total_templates: u64,
    /// Templates kept as the representative of their UMI family.
    pub unique_templates: u64,
    /// Templates marked as duplicates.
    pub duplicate_templates: u64,
    /// `duplicate_templates / total_templates`, or 0 when nothing was output.
    #[serde(with = "crate::float")]
    pub duplicate_rate: f64,
    /// Reads written to the output.
    pub total_reads: u64,
    /// Reads not marked as duplicates.
    pub unique_reads: u64,
    /// Reads marked as duplicates.
    pub duplicate_reads: u64,
    /// Secondary alignments seen in the output.
    pub secondary_reads: u64,
    /// Supplementary alignments seen in the output.
    pub supplementary_reads: u64,
    /// Secondary/supplementary reads lacking the `tc` tag that `fgumi zipper`
    /// adds and template-coordinate ordering requires.
    pub missing_tc_tag: u64,
}

impl DeduplicationMetrics {
    /// Creates a metrics struct with all counts initialized to zero.
    #[must_use]
    pub fn new() -> Self {
        Self::default()
    }
}

impl Metric for DeduplicationMetrics {
    fn metric_name() -> &'static str {
        "deduplication"
    }
}

impl crate::ProcessingMetrics for DeduplicationMetrics {
    fn total_input(&self) -> u64 {
        self.filtered_templates + self.total_templates
    }

    fn total_output(&self) -> u64 {
        self.total_templates
    }

    fn total_filtered(&self) -> u64 {
        self.filtered_templates
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::template_filter::TemplateFilterReason;
    use fgoxide::io::DelimFile;
    use tempfile::NamedTempFile;

    fn header_of(metrics: DeduplicationMetrics) -> String {
        let file = NamedTempFile::new().expect("temp file");
        DelimFile::default().write_tsv(file.path(), [metrics]).expect("write");
        std::fs::read_to_string(file.path())
            .expect("read")
            .lines()
            .next()
            .expect("header")
            .to_string()
    }

    /// Every filter reason must have a column. This is the test that would have
    /// caught fgumi#739: the counts existed and were merged correctly, then were
    /// dropped at the serialization boundary because no column held them.
    #[test]
    fn every_filter_reason_has_a_column() {
        let header = header_of(DeduplicationMetrics::default());
        let columns: Vec<&str> = header.split('\t').collect();
        for reason in TemplateFilterReason::ALL {
            let column = format!("filtered_{}", reason.column_suffix());
            assert!(columns.contains(&column.as_str()), "missing {column} for {reason:?}");
        }
    }

    /// Column order is the file's contract; pin it whole so a field reorder is a
    /// deliberate, reviewed change rather than a silent schema break.
    #[test]
    fn serializes_expected_columns_in_order() {
        assert_eq!(
            header_of(DeduplicationMetrics::default()),
            "filtered_templates\tfiltered_malformed_record\tfiltered_no_primary_reads\t\
             filtered_unmapped\tfiltered_not_passing_filter\tfiltered_low_mapping_quality\t\
             filtered_low_mate_mapping_quality\tfiltered_missing_umi\tfiltered_ns_in_umi\t\
             filtered_umi_too_short\tpassthrough_templates\ttotal_templates\tunique_templates\t\
             duplicate_templates\tduplicate_rate\ttotal_reads\tunique_reads\tduplicate_reads\t\
             secondary_reads\tsupplementary_reads\tmissing_tc_tag"
        );
    }

    #[test]
    fn processing_metrics_reconciles_input_from_the_columns() {
        use crate::ProcessingMetrics;
        let metrics = DeduplicationMetrics {
            filtered_templates: 10,
            total_templates: 90,
            ..DeduplicationMetrics::default()
        };
        assert_eq!(metrics.total_input(), 100);
        assert_eq!(metrics.total_output(), 90);
        assert_eq!(metrics.total_filtered(), 10);
    }
}
