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
/// 4. **Pair/orphan breakdown** — the Picard `DuplicationMetrics` split. Because
///    `dedup` dedups each mate in its own position group, reads are classified by
///    their SAM flags: `mapped_pairs`/`duplicate_pairs` count mapped read *pairs*
///    (both mates mapped; Picard `READ_PAIRS_EXAMINED`/`READ_PAIR_DUPLICATES`,
///    already halved from the two mates to whole pairs), `mapped_orphans`/
///    `duplicate_orphans` count mapped reads with no mapped mate (Picard
///    `UNPAIRED_READS_EXAMINED`/`UNPAIRED_READ_DUPLICATES`), and `unmapped_pairs`,
///    `unmapped_orphans`, `unmated_templates` are read-unit diagnostics. Every
///    counted primary read falls in exactly one mapping bucket, so (ignoring the
///    at-most-one truncated half per library from integer halving)
///    `2*mapped_pairs + mapped_orphans + unmapped_pairs + unmapped_orphans ==
///    total_templates` and `2*duplicate_pairs + duplicate_orphans ==
///    duplicate_templates` — but only for templates holding a single primary
///    read, which is every deduplicated template (each mate is dedup'd in its
///    own position group). The one exception is a `--include-unmapped`
///    pass-through: a fully unmapped pair is emitted verbatim as one template
///    holding both mates, so it adds 1 to `total_templates` but 2 to
///    `unmapped_pairs`, and the left-hand side then exceeds `total_templates`.
/// 5. **Complexity** — `percent_duplication` (Picard-parity pair/orphan units) and
///    the Lander-Waterman `estimated_library_size`.
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
    /// Sample name — the comma-joined unique `@RG SM:` values from the header,
    /// or the `--sample` override. Constant across every row of a run; leads
    /// the schema to match dupblaster/Picard `--stats` output.
    pub sample: String,
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
    /// Mapped read pairs examined, both mates mapped (Picard `READ_PAIRS_EXAMINED`)
    pub mapped_pairs: u64,
    /// Mapped read pairs marked duplicate (Picard `READ_PAIR_DUPLICATES`)
    pub duplicate_pairs: u64,
    /// Mapped reads with no mapped mate: single-end or mate unmapped (Picard `UNPAIRED_READS_EXAMINED`)
    pub mapped_orphans: u64,
    /// Single-mapped-end reads marked duplicate (Picard `UNPAIRED_READ_DUPLICATES`)
    pub duplicate_orphans: u64,
    /// Unmapped primary reads carrying the paired flag (diagnostic)
    pub unmapped_pairs: u64,
    /// Unmapped primary reads not carrying the paired flag (diagnostic)
    pub unmapped_orphans: u64,
    /// Primary reads with no paired flag, i.e. single-end (diagnostic)
    pub unmated_templates: u64,
    /// Picard `PERCENT_DUPLICATION`-parity duplicate fraction, in primary-read
    /// units: `(2*duplicate_pairs + duplicate_orphans) / (2*mapped_pairs +
    /// mapped_orphans)`.
    ///
    /// Distinct from the template-level `duplicate_rate` above. Counts each
    /// mapped pair as its two examined reads and excludes unmapped, secondary,
    /// and supplementary records, matching Picard `MarkDuplicates`. Zero when the
    /// denominator (no mapped primary reads) is zero.
    #[serde(with = "crate::float")]
    pub percent_duplication: f64,
    /// Lander-Waterman estimate of the complete library size, or empty when the
    /// estimate is undefined (no mapped pairs, or no duplicate pairs observed).
    ///
    /// Estimated from complete mapped pairs only (`mapped_pairs` as `n_pairs`,
    /// `mapped_pairs - duplicate_pairs` as `n_unique`), matching Picard
    /// `MarkDuplicates`, which excludes orphans from the library-size estimate.
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
    /// Sample name (comma-joined `@RG SM:` values or the `--sample` override).
    pub sample: String,
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
    /// Mapped read pairs examined, both mates mapped (already halved to pairs).
    pub mapped_pairs: u64,
    /// Mapped read pairs marked duplicate (already halved to pairs).
    pub duplicate_pairs: u64,
    /// Mapped reads with no mapped mate (single-end or mate unmapped).
    pub mapped_orphans: u64,
    /// Single-mapped-end reads marked duplicate.
    pub duplicate_orphans: u64,
    /// Unmapped primary reads carrying the paired flag.
    pub unmapped_pairs: u64,
    /// Unmapped primary reads not carrying the paired flag.
    pub unmapped_orphans: u64,
    /// Primary reads with no paired flag (single-end).
    pub unmated_templates: u64,
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
            sample,
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
            mapped_pairs,
            duplicate_pairs,
            mapped_orphans,
            duplicate_orphans,
            unmapped_pairs,
            unmapped_orphans,
            unmated_templates,
        } = counts;
        Self {
            sample,
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
            mapped_pairs,
            duplicate_pairs,
            mapped_orphans,
            duplicate_orphans,
            unmapped_pairs,
            unmapped_orphans,
            unmated_templates,
            // Picard `PERCENT_DUPLICATION`: each mapped pair counts as its two
            // examined reads, orphans as one, excluding unmapped/secondary/
            // supplementary. Zero when no mapped primary reads were examined.
            percent_duplication: picard_percent_duplication(
                mapped_pairs,
                duplicate_pairs,
                mapped_orphans,
                duplicate_orphans,
            ),
            // Lander-Waterman estimate from complete mapped pairs only (Picard
            // excludes orphans). `saturating_sub` guards a future bug from
            // underflowing rather than panicking.
            estimated_library_size: estimate_library_size(
                mapped_pairs,
                mapped_pairs.saturating_sub(duplicate_pairs),
            ),
        }
    }
}

/// Picard `PERCENT_DUPLICATION` in primary-read units:
/// `(2*duplicate_pairs + duplicate_orphans) / (2*mapped_pairs + mapped_orphans)`.
///
/// Each mapped pair contributes its two examined reads; orphans contribute one.
/// Unmapped, secondary, and supplementary records are excluded (they are not in
/// these counts). Returns `0.0` when the denominator is zero (no mapped primary
/// reads examined), matching `frac_u64`'s zero-denominator convention.
#[expect(clippy::cast_precision_loss, reason = "metric counts never exceed 2^53")]
fn picard_percent_duplication(
    mapped_pairs: u64,
    duplicate_pairs: u64,
    mapped_orphans: u64,
    duplicate_orphans: u64,
) -> f64 {
    let examined = 2 * mapped_pairs + mapped_orphans;
    if examined == 0 {
        0.0
    } else {
        (2 * duplicate_pairs + duplicate_orphans) as f64 / examined as f64
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
    /// Templates added since the previous snapshot (this window's depth). Equals
    /// `templates_seen` on the first snapshot.
    pub window_templates: u64,
    /// Marginal duplicate fraction over just this window's templates
    /// (`window_duplicate_templates / window_templates`). Often the more legible
    /// view of the saturation curve than the cumulative `duplicate_fraction`,
    /// since it isolates each depth band instead of averaging over all prior
    /// ones. Mirrors dupblaster's per-window complexity columns.
    #[serde(with = "crate::float")]
    pub window_duplicate_fraction: f64,
}

impl DuplicationLadderMetrics {
    /// Builds one snapshot row from cumulative and per-window per-library counts.
    ///
    /// `templates_seen`/`duplicate_templates` are the cumulative totals at this
    /// snapshot; `window_templates`/`window_duplicate_templates` cover only the
    /// templates added since the previous snapshot.
    #[must_use]
    pub fn new(
        library: String,
        templates_seen: u64,
        duplicate_templates: u64,
        window_templates: u64,
        window_duplicate_templates: u64,
    ) -> Self {
        Self {
            library,
            templates_seen,
            duplicate_fraction: crate::frac_u64(duplicate_templates, templates_seen),
            window_templates,
            window_duplicate_fraction: crate::frac_u64(
                window_duplicate_templates,
                window_templates,
            ),
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
            "sample\tlibrary\tfiltered_templates\tfiltered_malformed_record\tfiltered_no_primary_reads\t\
             filtered_unmapped\tfiltered_not_passing_filter\tfiltered_low_mapping_quality\t\
             filtered_low_mate_mapping_quality\tfiltered_missing_umi\tfiltered_ns_in_umi\t\
             filtered_umi_too_short\tpassthrough_templates\ttotal_templates\tunique_templates\t\
             duplicate_templates\tduplicate_rate\ttotal_reads\tunique_reads\tduplicate_reads\t\
             secondary_reads\tsupplementary_reads\tmissing_tc_tag\t\
             mapped_pairs\tduplicate_pairs\tmapped_orphans\tduplicate_orphans\t\
             unmapped_pairs\tunmapped_orphans\tunmated_templates\t\
             percent_duplication\testimated_library_size"
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

    /// `(mapped_pairs, duplicate_pairs, mapped_orphans, duplicate_orphans)` — the
    /// pair/orphan units that drive `percent_duplication` (Picard `(2*dup_pairs +
    /// dup_orphans) / (2*mapped_pairs + mapped_orphans)`) and, from mapped pairs
    /// only, `estimated_library_size`.
    type PairOrphan = (u64, u64, u64, u64);

    #[rstest]
    // duplicate_rate stays template-level; percent_duplication and library size
    // are now driven by the pair/orphan columns, not the read totals.
    #[case::paired_with_duplicates(
        "lib1".to_string(), (100, 10, 90), (100, 90, 0, 0), 0.9, 0.9, true)]
    #[case::paired_no_duplicates(
        "lib2".to_string(), (50, 50, 0), (50, 0, 0, 0), 0.0, 0.0, false)]
    #[case::empty_library("lib3".to_string(), (0, 0, 0), (0, 0, 0, 0), 0.0, 0.0, false)]
    // Orphan-only library: pct = 10/30; no mapped pairs, so library size undefined.
    #[case::orphans_only(
        "lib4".to_string(), (30, 20, 10), (0, 0, 30, 10), 1.0 / 3.0, 1.0 / 3.0, false)]
    fn from_counts_computes_rates_and_library_size(
        #[case] library: String,
        #[case] templates: Counts,
        #[case] pair_orphan: PairOrphan,
        #[case] expected_duplicate_rate: f64,
        #[case] expected_percent_duplication: f64,
        #[case] expect_library_size_defined: bool,
    ) {
        let (total_templates, unique_templates, duplicate_templates) = templates;
        let (mapped_pairs, duplicate_pairs, mapped_orphans, duplicate_orphans) = pair_orphan;
        let metrics = DeduplicationMetrics::from_counts(DeduplicationCounts {
            library: library.clone(),
            total_templates,
            unique_templates,
            duplicate_templates,
            mapped_pairs,
            duplicate_pairs,
            mapped_orphans,
            duplicate_orphans,
            ..DeduplicationCounts::default()
        });

        assert_eq!(metrics.library, library);
        assert_eq!(metrics.total_templates, total_templates);
        assert_eq!(metrics.unique_templates, unique_templates);
        assert_eq!(metrics.duplicate_templates, duplicate_templates);
        assert!((metrics.duplicate_rate - expected_duplicate_rate).abs() < 1e-9);
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

    /// The pair/orphan breakdown counts are copied through `from_counts`
    /// unchanged, and the documented reconciliation invariants hold on the
    /// resulting row.
    #[test]
    fn from_counts_carries_pair_orphan_breakdown_through() {
        // Pairs are in whole-pair units (already halved); orphan/unmapped counts
        // are read units. Values chosen so the documented read-unit invariants
        // hold: 2*35 + 20 + 6 + 4 = 100, and 2*12 + 6 = 30.
        let metrics = DeduplicationMetrics::from_counts(DeduplicationCounts {
            library: "lib1".to_string(),
            total_templates: 100,
            duplicate_templates: 30,
            mapped_pairs: 35,
            duplicate_pairs: 12,
            mapped_orphans: 20,
            duplicate_orphans: 6,
            unmapped_pairs: 6,
            unmapped_orphans: 4,
            unmated_templates: 24,
            ..DeduplicationCounts::default()
        });
        assert_eq!(metrics.mapped_pairs, 35);
        assert_eq!(metrics.duplicate_pairs, 12);
        assert_eq!(metrics.mapped_orphans, 20);
        assert_eq!(metrics.duplicate_orphans, 6);
        assert_eq!(metrics.unmapped_pairs, 6);
        assert_eq!(metrics.unmapped_orphans, 4);
        assert_eq!(metrics.unmated_templates, 24);
        // Documented read-unit invariants: each pair is two counted primary
        // reads, so pairs count double against the template total.
        assert_eq!(
            2 * metrics.mapped_pairs
                + metrics.mapped_orphans
                + metrics.unmapped_pairs
                + metrics.unmapped_orphans,
            metrics.total_templates
        );
        assert_eq!(
            2 * metrics.duplicate_pairs + metrics.duplicate_orphans,
            metrics.duplicate_templates
        );
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
                // Mapped pairs drive the library-size estimate: 90 of 100 pairs
                // duplicate, so the Lander-Waterman estimate is defined.
                mapped_pairs: 100,
                duplicate_pairs: 90,
                ..DeduplicationCounts::default()
            }),
            DeduplicationMetrics::from_counts(DeduplicationCounts {
                library: "lib2".to_string(),
                total_templates: 50,
                unique_templates: 50,
                total_reads: 100,
                unique_reads: 100,
                // All unique pairs, so the estimate is undefined (empty cell).
                mapped_pairs: 50,
                duplicate_pairs: 0,
                ..DeduplicationCounts::default()
            }),
        ];

        let file = NamedTempFile::new().expect("temp file");
        DelimFile::default().write_tsv(file.path(), &metrics).expect("write");
        let content = std::fs::read_to_string(file.path()).expect("read");
        let mut lines = content.lines();

        let header: Vec<&str> = lines.next().expect("header line").split('\t').collect();
        assert_eq!(header.first(), Some(&"sample"), "sample column must be first: {header:?}");
        assert_eq!(header.get(1), Some(&"library"), "library column must be second: {header:?}");
        let size_col = header
            .iter()
            .position(|&h| h == "estimated_library_size")
            .expect("estimated_library_size column present");

        let lib1_row: Vec<&str> = lines.next().expect("lib1 row").split('\t').collect();
        assert_eq!(lib1_row.get(1), Some(&"lib1"), "library is the second column: {lib1_row:?}");
        assert!(!lib1_row[size_col].is_empty(), "lib1 has duplicates: {lib1_row:?}");

        let lib2_row: Vec<&str> = lines.next().expect("lib2 row").split('\t').collect();
        assert_eq!(lib2_row.get(1), Some(&"lib2"), "library is the second column: {lib2_row:?}");
        assert!(
            lib2_row[size_col].is_empty(),
            "lib2's None estimate must render empty: {lib2_row:?}"
        );
    }

    #[test]
    fn ladder_metric_computes_cumulative_and_window_fractions() {
        // Cumulative 25/100 = 0.25; this window 20/40 = 0.5 (a denser band).
        let row = DuplicationLadderMetrics::new("lib1".to_string(), 100, 25, 40, 20);
        assert_eq!(row.library, "lib1");
        assert_eq!(row.templates_seen, 100);
        assert!((row.duplicate_fraction - 0.25).abs() < 1e-9);
        assert_eq!(row.window_templates, 40);
        assert!((row.window_duplicate_fraction - 0.5).abs() < 1e-9);
        // Div-by-zero guard: zero templates yields defined 0.0 fractions.
        let empty = DuplicationLadderMetrics::new("lib2".to_string(), 0, 0, 0, 0);
        assert!(empty.duplicate_fraction.abs() < 1e-9);
        assert!(empty.window_duplicate_fraction.abs() < 1e-9);
    }
}
