//! Metrics for the `dedup` command.
//!
//! [`DeduplicationMetrics`] is the serializable row type written by `fgumi dedup
//! --metrics`. One row is written per library observed in the input, plus a
//! final aggregate row summing across all libraries (see the `dedup` command
//! for the exact library-splitting and total-row naming).

use serde::{Deserialize, Serialize};

use crate::Metric;
use crate::library_size::estimate_library_size;

/// Metrics written by `fgumi dedup --metrics`, one row per library plus an
/// aggregate total row.
///
/// Columns are grouped, in order:
///
/// 0. **Identity** — `library`: the library the row summarizes, or a sentinel
///    ("Unknown Library" for reads with no resolvable `@RG`/`LB`, "All Reads"
///    for the aggregate row) — see the `dedup` command for the exact naming.
/// 1. **Filtering** — `filtered_templates` through `passthrough_templates`.
///    `dedup` is a read *filter* as well as a duplicate marker: templates failing
///    the pre-grouping filter are dropped and never reach the output. The
///    `filtered_*` columns give the count per reason in **template** units — note
///    that `group`'s fgbio-shaped `.grouping_metrics.txt` counts primary
///    *records*, so the two files' discard counts are not comparable.
/// 2. **Deduplication** — `total_templates` through `duplicate_reads`, covering
///    only what survived filtering. These count what *passed the filter*, not what
///    reaches the output file: under `--remove-duplicates` the duplicates are
///    dropped at write time, after they have been counted here.
/// 3. **Diagnostics** — `secondary_reads`, `supplementary_reads`,
///    `missing_tc_tag`.
/// 4. **Complexity** — `percent_duplication` (read-level, Picard-parity) and the
///    Lander-Waterman `estimated_library_size`.
///
/// Templates read from the input are `filtered_templates + total_templates`; that
/// sum is not emitted as its own column because it is derivable, matching
/// `UmiGroupingMetrics`' treatment of its derivable fields.
// NB: plain code spans above, never intra-doc links — this doc is rendered verbatim
// into docs/src/metrics/deduplication-metrics.md, where a rustdoc link would come
// out as literal markdown. Keep field docs to a single short line: they become the
// Description column of that page's table.
#[derive(Debug, Clone, Default, PartialEq, Serialize, Deserialize)]
pub struct DeduplicationMetrics {
    /// Library name, "Unknown Library", or "All Reads" (the aggregate row)
    pub library: String,
    /// Templates dropped by the filter (any reason)
    pub filtered_templates: u64,
    /// Templates dropped (record shorter than the minimum BAM record length)
    pub filtered_malformed_record: u64,
    /// Templates dropped (no primary R1 or R2)
    pub filtered_no_primary_reads: u64,
    /// Templates dropped (no mapped primary read)
    pub filtered_unmapped: u64,
    /// Templates dropped (QC-fail flag)
    pub filtered_not_passing_filter: u64,
    /// Templates dropped (mapping quality below `--min-map-q`)
    pub filtered_low_mapping_quality: u64,
    /// Templates dropped (mate `MQ` below `--min-map-q`)
    pub filtered_low_mate_mapping_quality: u64,
    /// Templates dropped (no UMI tag)
    pub filtered_missing_umi: u64,
    /// Templates dropped (N base in the UMI)
    pub filtered_ns_in_umi: u64,
    /// Templates dropped (UMI shorter than `--min-umi-length`)
    pub filtered_umi_too_short: u64,
    /// Templates emitted untouched by `--include-unmapped`
    pub passthrough_templates: u64,
    /// Templates that passed the filter
    pub total_templates: u64,
    /// Templates kept as their family's representative
    pub unique_templates: u64,
    /// Templates marked as duplicates
    pub duplicate_templates: u64,
    /// `duplicate_templates` / `total_templates`
    #[serde(with = "crate::float")]
    pub duplicate_rate: f64,
    /// Reads in templates that passed the filter
    pub total_reads: u64,
    /// Reads not marked as duplicates
    pub unique_reads: u64,
    /// Reads marked as duplicates
    pub duplicate_reads: u64,
    /// Secondary alignments that passed the filter
    pub secondary_reads: u64,
    /// Supplementary alignments that passed the filter
    pub supplementary_reads: u64,
    /// Secondary/supplementary reads lacking the `tc` tag
    pub missing_tc_tag: u64,
    /// Read-level duplicate fraction (`duplicate_reads` / `total_reads`).
    ///
    /// Distinct from the template-level `duplicate_rate` above; this is the
    /// Picard `PERCENT_DUPLICATION`-parity column.
    #[serde(with = "crate::float")]
    pub percent_duplication: f64,
    /// Lander-Waterman estimate of the complete library size, or empty when the
    /// estimate is undefined (no templates, or no duplicates observed).
    ///
    /// fgumi counts *templates*, not read pairs directly; this treats
    /// `total_templates`/`unique_templates` as `n_pairs`/`n_unique` (each
    /// template stands in for one read pair in paired-end data). This is the
    /// same approximation `duplicate_rate` above already makes by reporting a
    /// template-level rather than a read-pair-level rate.
    pub estimated_library_size: Option<u64>,
}

/// Raw per-library (or aggregate) counts used to build a [`DeduplicationMetrics`]
/// row.
///
/// Everything except the three derived columns — `duplicate_rate`,
/// `percent_duplication`, and `estimated_library_size` — is copied through 1:1
/// by [`DeduplicationMetrics::from_counts`], which computes those three. Using a
/// struct rather than a 20-plus-argument constructor keeps the call sites (one
/// per library, plus the total) readable and self-labeling.
#[derive(Debug, Clone, Default, PartialEq)]
pub struct DeduplicationCounts {
    /// Library name, "Unknown Library", or "All Reads".
    pub library: String,
    /// Templates dropped by the filter (any reason).
    pub filtered_templates: u64,
    /// Templates dropped (record shorter than the minimum BAM record length).
    pub filtered_malformed_record: u64,
    /// Templates dropped (no primary R1 or R2).
    pub filtered_no_primary_reads: u64,
    /// Templates dropped (no mapped primary read).
    pub filtered_unmapped: u64,
    /// Templates dropped (QC-fail flag).
    pub filtered_not_passing_filter: u64,
    /// Templates dropped (mapping quality below `--min-map-q`).
    pub filtered_low_mapping_quality: u64,
    /// Templates dropped (mate `MQ` below `--min-map-q`).
    pub filtered_low_mate_mapping_quality: u64,
    /// Templates dropped (no UMI tag).
    pub filtered_missing_umi: u64,
    /// Templates dropped (N base in the UMI).
    pub filtered_ns_in_umi: u64,
    /// Templates dropped (UMI shorter than `--min-umi-length`).
    pub filtered_umi_too_short: u64,
    /// Templates emitted untouched by `--include-unmapped`.
    pub passthrough_templates: u64,
    /// Templates that passed the filter.
    pub total_templates: u64,
    /// Templates kept as their family's representative.
    pub unique_templates: u64,
    /// Templates marked as duplicates.
    pub duplicate_templates: u64,
    /// Reads in templates that passed the filter.
    pub total_reads: u64,
    /// Reads not marked as duplicates.
    pub unique_reads: u64,
    /// Reads marked as duplicates.
    pub duplicate_reads: u64,
    /// Secondary alignments that passed the filter.
    pub secondary_reads: u64,
    /// Supplementary alignments that passed the filter.
    pub supplementary_reads: u64,
    /// Secondary/supplementary reads lacking the `tc` tag.
    pub missing_tc_tag: u64,
}

impl DeduplicationMetrics {
    /// Creates a metrics struct with all counts initialized to zero.
    #[must_use]
    pub fn new() -> Self {
        Self::default()
    }

    /// Builds a metrics row from summed counts for one library (or the aggregate
    /// total across all libraries), computing the derived `duplicate_rate`,
    /// `percent_duplication`, and `estimated_library_size` columns.
    #[must_use]
    pub fn from_counts(counts: DeduplicationCounts) -> Self {
        let DeduplicationCounts {
            library,
            filtered_templates,
            filtered_malformed_record,
            filtered_no_primary_reads,
            filtered_unmapped,
            filtered_not_passing_filter,
            filtered_low_mapping_quality,
            filtered_low_mate_mapping_quality,
            filtered_missing_umi,
            filtered_ns_in_umi,
            filtered_umi_too_short,
            passthrough_templates,
            total_templates,
            unique_templates,
            duplicate_templates,
            total_reads,
            unique_reads,
            duplicate_reads,
            secondary_reads,
            supplementary_reads,
            missing_tc_tag,
        } = counts;
        Self {
            library,
            filtered_templates,
            filtered_malformed_record,
            filtered_no_primary_reads,
            filtered_unmapped,
            filtered_not_passing_filter,
            filtered_low_mapping_quality,
            filtered_low_mate_mapping_quality,
            filtered_missing_umi,
            filtered_ns_in_umi,
            filtered_umi_too_short,
            passthrough_templates,
            total_templates,
            unique_templates,
            duplicate_templates,
            duplicate_rate: crate::frac_u64(duplicate_templates, total_templates),
            total_reads,
            unique_reads,
            duplicate_reads,
            secondary_reads,
            supplementary_reads,
            missing_tc_tag,
            percent_duplication: crate::frac_u64(duplicate_reads, total_reads),
            estimated_library_size: estimate_library_size(total_templates, unique_templates),
        }
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

/// Serializable row for the `--duplication-ladder` sampled duplication
/// saturation curve: "after `templates_seen` templates processed (in
/// coordinate order) for this library, what cumulative fraction were
/// duplicates".
///
/// One row per (library, snapshot) — a library gets a row each time its
/// cumulative `templates_seen` crosses a multiple of `--ladder-interval`,
/// plus a final row at its true total. See the `dedup` command for the
/// serial-order accumulation this depends on.
#[derive(Debug, Clone, Default, PartialEq, Serialize, Deserialize)]
pub struct DuplicationLadderMetrics {
    /// Library name this row summarizes (or "Unknown Library"; never an
    /// aggregate "All Reads" row — the ladder is inherently per-library).
    pub library: String,
    /// Cumulative templates seen for this library at this snapshot.
    pub templates_seen: u64,
    /// Cumulative duplicate fraction (`duplicate_templates / templates_seen`)
    /// at this snapshot.
    #[serde(with = "crate::float")]
    pub duplicate_fraction: f64,
}

impl DuplicationLadderMetrics {
    /// Builds one snapshot row from cumulative per-library counts.
    #[must_use]
    pub fn new(library: String, templates_seen: u64, duplicate_templates: u64) -> Self {
        Self {
            library,
            templates_seen,
            duplicate_fraction: crate::frac_u64(duplicate_templates, templates_seen),
        }
    }
}

impl Metric for DuplicationLadderMetrics {
    fn metric_name() -> &'static str {
        "duplication_ladder"
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::template_filter::TemplateFilterReason;
    use fgoxide::io::DelimFile;
    use rstest::rstest;
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
    /// deliberate, reviewed change rather than a silent schema break. `library`
    /// leads; the filtering/dedup/diagnostics block matches the pre-per-library
    /// schema; `percent_duplication` and `estimated_library_size` trail.
    #[test]
    fn serializes_expected_columns_in_order() {
        assert_eq!(
            header_of(DeduplicationMetrics::default()),
            "library\tfiltered_templates\tfiltered_malformed_record\tfiltered_no_primary_reads\t\
             filtered_unmapped\tfiltered_not_passing_filter\tfiltered_low_mapping_quality\t\
             filtered_low_mate_mapping_quality\tfiltered_missing_umi\tfiltered_ns_in_umi\t\
             filtered_umi_too_short\tpassthrough_templates\ttotal_templates\tunique_templates\t\
             duplicate_templates\tduplicate_rate\ttotal_reads\tunique_reads\tduplicate_reads\t\
             secondary_reads\tsupplementary_reads\tmissing_tc_tag\tpercent_duplication\t\
             estimated_library_size"
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

    #[test]
    fn metric_name_is_deduplication() {
        assert_eq!(DeduplicationMetrics::metric_name(), "deduplication");
    }

    /// `(total, unique, duplicate)` counts — shared shape for both the
    /// template and read triples below.
    type Counts = (u64, u64, u64);

    #[rstest]
    #[case::with_duplicates("lib1".to_string(), (100, 10, 90), (200, 20, 180), 0.9, 0.9, true)]
    #[case::no_duplicates("lib2".to_string(), (50, 50, 0), (100, 100, 0), 0.0, 0.0, false)]
    #[case::empty_library("lib3".to_string(), (0, 0, 0), (0, 0, 0), 0.0, 0.0, false)]
    fn from_counts_computes_rates_and_library_size(
        #[case] library: String,
        #[case] templates: Counts,
        #[case] reads: Counts,
        #[case] expected_duplicate_rate: f64,
        #[case] expected_percent_duplication: f64,
        #[case] expect_library_size_defined: bool,
    ) {
        let (total_templates, unique_templates, duplicate_templates) = templates;
        let (total_reads, unique_reads, duplicate_reads) = reads;
        let metrics = DeduplicationMetrics::from_counts(DeduplicationCounts {
            library: library.clone(),
            total_templates,
            unique_templates,
            duplicate_templates,
            total_reads,
            unique_reads,
            duplicate_reads,
            ..DeduplicationCounts::default()
        });

        assert_eq!(metrics.library, library);
        assert_eq!(metrics.total_templates, total_templates);
        assert_eq!(metrics.unique_templates, unique_templates);
        assert_eq!(metrics.duplicate_templates, duplicate_templates);
        assert!((metrics.duplicate_rate - expected_duplicate_rate).abs() < 1e-9);
        assert_eq!(metrics.total_reads, total_reads);
        assert_eq!(metrics.unique_reads, unique_reads);
        assert_eq!(metrics.duplicate_reads, duplicate_reads);
        assert!((metrics.percent_duplication - expected_percent_duplication).abs() < 1e-9);
        assert_eq!(
            metrics.estimated_library_size.is_some(),
            expect_library_size_defined,
            "estimated_library_size definedness mismatch: got {:?}",
            metrics.estimated_library_size
        );
    }

    /// `from_counts` copies the per-reason filter counts through unchanged, and
    /// `ProcessingMetrics` reconciles input from them.
    #[test]
    fn from_counts_carries_filter_counts_through() {
        use crate::ProcessingMetrics;
        let metrics = DeduplicationMetrics::from_counts(DeduplicationCounts {
            library: "lib1".to_string(),
            filtered_templates: 7,
            filtered_unmapped: 4,
            filtered_missing_umi: 3,
            passthrough_templates: 2,
            total_templates: 90,
            unique_templates: 80,
            duplicate_templates: 10,
            ..DeduplicationCounts::default()
        });
        assert_eq!(metrics.filtered_templates, 7);
        assert_eq!(metrics.filtered_unmapped, 4);
        assert_eq!(metrics.filtered_missing_umi, 3);
        assert_eq!(metrics.passthrough_templates, 2);
        assert_eq!(metrics.total_input(), 97);
        assert_eq!(metrics.total_filtered(), 7);
    }

    #[test]
    fn serializes_with_library_column_first_and_empty_cell_for_none() {
        // lib1 has duplicates (estimate defined); lib2 has none (estimate undefined).
        let metrics = vec![
            DeduplicationMetrics::from_counts(DeduplicationCounts {
                library: "lib1".to_string(),
                total_templates: 100,
                unique_templates: 10,
                duplicate_templates: 90,
                total_reads: 200,
                unique_reads: 20,
                duplicate_reads: 180,
                ..DeduplicationCounts::default()
            }),
            DeduplicationMetrics::from_counts(DeduplicationCounts {
                library: "lib2".to_string(),
                total_templates: 50,
                unique_templates: 50,
                total_reads: 100,
                unique_reads: 100,
                ..DeduplicationCounts::default()
            }),
        ];

        let file = NamedTempFile::new().expect("temp file");
        DelimFile::default().write_tsv(file.path(), &metrics).expect("write");
        let content = std::fs::read_to_string(file.path()).expect("read");
        let mut lines = content.lines();

        let header: Vec<&str> = lines.next().expect("header line").split('\t').collect();
        assert_eq!(header.first(), Some(&"library"), "library column must be first: {header:?}");
        let size_col = header
            .iter()
            .position(|&h| h == "estimated_library_size")
            .expect("estimated_library_size column present");

        let lib1_row: Vec<&str> = lines.next().expect("lib1 row").split('\t').collect();
        assert_eq!(lib1_row.first(), Some(&"lib1"));
        assert!(!lib1_row[size_col].is_empty(), "lib1 has duplicates: {lib1_row:?}");

        let lib2_row: Vec<&str> = lines.next().expect("lib2 row").split('\t').collect();
        assert_eq!(lib2_row.first(), Some(&"lib2"));
        assert!(
            lib2_row[size_col].is_empty(),
            "lib2's None estimate must render empty: {lib2_row:?}"
        );
    }

    #[test]
    fn ladder_metric_computes_cumulative_fraction() {
        let row = DuplicationLadderMetrics::new("lib1".to_string(), 100, 25);
        assert_eq!(row.library, "lib1");
        assert_eq!(row.templates_seen, 100);
        assert!((row.duplicate_fraction - 0.25).abs() < 1e-9);
        // Div-by-zero guard: zero templates yields a defined 0.0 fraction.
        let empty = DuplicationLadderMetrics::new("lib2".to_string(), 0, 0);
        assert!(empty.duplicate_fraction.abs() < 1e-9);
    }
}
