//! Common infrastructure for consensus calling commands.
//!
//! This module provides shared traits and utilities used by simplex, duplex, and codec
//! consensus calling commands to reduce code duplication.

#[cfg(feature = "codec")]
use crate::consensus::codec_caller::CodecConsensusStats;
use crate::consensus_caller::ConsensusCallingStats;
use crate::metrics::consensus::ConsensusMetrics;
use crate::overlapping_consensus::CorrectionStats;
use anyhow::Result;
use bstr::BString;
use noodles::sam::Header;
use noodles::sam::header::record::value::Map;
use noodles::sam::header::record::value::map::ReadGroup;
use noodles::sam::header::record::value::map::header::tag as header_tag;
use noodles::sam::header::record::value::map::read_group::tag as rg_tag;
use noodles::sam::header::record::value::map::tag::Other;
use tracing::info;

/// Trait for converting command-specific statistics to metrics.
pub trait ConsensusStatsOps: Clone + Default + Send {
    /// Merges another stats instance into this one.
    fn merge(&mut self, other: &Self);

    /// Converts to `ConsensusMetrics` for TSV output.
    fn to_metrics(&self) -> ConsensusMetrics;
}

impl ConsensusStatsOps for ConsensusCallingStats {
    fn merge(&mut self, other: &Self) {
        ConsensusCallingStats::merge(self, other);
    }

    fn to_metrics(&self) -> ConsensusMetrics {
        let mut metrics = ConsensusMetrics::new();
        metrics.total_input_reads = self.total_reads as u64;
        metrics.consensus_reads = self.consensus_reads as u64;
        metrics.filtered_reads = self.filtered_reads as u64;

        // Convert rejection reasons to centralized reasons
        for (reason, count) in &self.rejection_reasons {
            metrics.add_rejection(reason.to_centralized(), *count as u64);
        }

        metrics
    }
}

#[cfg(feature = "codec")]
impl ConsensusStatsOps for CodecConsensusStats {
    fn merge(&mut self, other: &Self) {
        self.total_input_reads += other.total_input_reads;
        self.consensus_reads_generated += other.consensus_reads_generated;
        self.reads_filtered += other.reads_filtered;
        self.consensus_reads_rejected_hdd += other.consensus_reads_rejected_hdd;
        self.consensus_bases_emitted += other.consensus_bases_emitted;
        self.consensus_duplex_bases_emitted += other.consensus_duplex_bases_emitted;
        self.duplex_disagreement_base_count += other.duplex_disagreement_base_count;

        // Merge rejection reasons
        for (reason, count) in &other.rejection_reasons {
            *self.rejection_reasons.entry(*reason).or_insert(0) += count;
        }
    }

    fn to_metrics(&self) -> ConsensusMetrics {
        let mut metrics = ConsensusMetrics::new();
        metrics.total_input_reads = self.total_input_reads;
        metrics.consensus_reads = self.consensus_reads_generated;
        metrics.filtered_reads = self.reads_filtered;

        // Convert rejection reasons to centralized reasons
        for (reason, count) in &self.rejection_reasons {
            metrics.add_rejection(reason.to_centralized(), *count as u64);
        }

        metrics
    }
}

/// Collapses a read group tag across all input read groups.
///
/// Iterates all input read groups, extracts the given tag, deduplicates values
/// (preserving insertion order), and returns them comma-joined. Returns `None`
/// if no read groups have the tag set.
///
/// For the `PL` (platform) tag, values are uppercased per the SAM spec before
/// deduplication.
fn collapse_read_group_tag(input_header: &Header, tag: Other<rg_tag::Standard>) -> Option<String> {
    let is_platform = tag == rg_tag::PLATFORM;
    let mut seen = Vec::new();
    let mut seen_set = std::collections::HashSet::new();

    for (_id, rg) in input_header.read_groups() {
        if let Some(value) = rg.other_fields().get(&tag) {
            let value = value.to_string();
            if value.is_empty() {
                continue;
            }
            let value = if is_platform { value.to_uppercase() } else { value };
            if seen_set.insert(value.clone()) {
                seen.push(value);
            }
        }
    }

    if seen.is_empty() { None } else { Some(seen.join(",")) }
}

/// Creates an output header for unmapped consensus reads.
///
/// This creates a header with:
/// - A single read group with the specified ID and attributes collapsed from all
///   input read groups (SM, LB, PL, PU, CN, DS tags are deduplicated and
///   comma-joined when multiple distinct values exist)
/// - Sort order set to "unknown" and group order to "query"
/// - A comment indicating the number of input read groups
/// - A @PG record with version and command line
pub fn create_unmapped_consensus_header(
    input_header: &Header,
    read_group_id: &str,
    comment_prefix: &str,
    command_line: &str,
) -> Result<Header> {
    let mut output_header = Header::builder();

    // Create read group with collapsed attributes from all input read groups
    let mut rg_builder = Map::<ReadGroup>::builder();
    let tags_to_collapse = [
        rg_tag::SAMPLE,
        rg_tag::LIBRARY,
        rg_tag::PLATFORM,
        rg_tag::PLATFORM_UNIT,
        rg_tag::SEQUENCING_CENTER,
        rg_tag::DESCRIPTION,
    ];
    for tag in tags_to_collapse {
        if let Some(value) = collapse_read_group_tag(input_header, tag) {
            rg_builder = rg_builder.insert(tag, value);
        }
    }
    let new_rg = rg_builder.build()?;
    output_header = output_header.add_read_group(BString::from(read_group_id), new_rg);

    // Set sort order
    let header_map = Map::<noodles::sam::header::record::value::map::Header>::builder()
        .insert(header_tag::SORT_ORDER, BString::from("unknown"))
        .insert(header_tag::GROUP_ORDER, BString::from("query"))
        .build()?;
    output_header = output_header.set_header(header_map);

    // Add comment
    let rg_count = input_header.read_groups().len();
    output_header = output_header.add_comment(format!(
        "{comment_prefix} {read_group_id} contains consensus reads generated from {rg_count} input read groups."
    ));

    // Add @PG record
    output_header = crate::commands::common::add_pg_to_builder(output_header, command_line)?;

    Ok(output_header.build())
}

/// Logs overlapping consensus statistics if enabled.
pub fn log_overlapping_stats(stats: &CorrectionStats) {
    info!("Overlapping consensus statistics:");
    info!("  Overlapping bases examined: {}", stats.overlapping_bases);
    info!("  Bases agreeing: {}", stats.bases_agreeing);
    info!("  Bases disagreeing: {}", stats.bases_disagreeing);
    info!("  Bases corrected: {}", stats.bases_corrected);
}

#[cfg(test)]
mod tests {
    use super::*;

    /// Builds a header with the given read groups for testing.
    #[allow(clippy::type_complexity)]
    fn build_test_header(rgs: Vec<(&str, Vec<(Other<rg_tag::Standard>, &str)>)>) -> Header {
        let mut builder = Header::builder();
        for (id, tags) in rgs {
            let mut rg_builder = Map::<ReadGroup>::builder();
            for (tag, value) in tags {
                rg_builder = rg_builder.insert(tag, value.to_string());
            }
            builder = builder.add_read_group(
                BString::from(id),
                rg_builder.build().expect("failed to build read group"),
            );
        }
        builder.build()
    }

    #[test]
    fn test_collapse_read_group_tag_no_read_groups() {
        let header = Header::builder().build();
        assert_eq!(collapse_read_group_tag(&header, rg_tag::SAMPLE), None);
    }

    #[test]
    fn test_collapse_read_group_tag_single_value() {
        let header = build_test_header(vec![("RG1", vec![(rg_tag::SAMPLE, "SampleA")])]);
        assert_eq!(collapse_read_group_tag(&header, rg_tag::SAMPLE), Some("SampleA".to_string()),);
    }

    #[test]
    fn test_collapse_read_group_tag_duplicate_values() {
        let header = build_test_header(vec![
            ("RG1", vec![(rg_tag::SAMPLE, "SampleA")]),
            ("RG2", vec![(rg_tag::SAMPLE, "SampleA")]),
        ]);
        assert_eq!(collapse_read_group_tag(&header, rg_tag::SAMPLE), Some("SampleA".to_string()),);
    }

    #[test]
    fn test_collapse_read_group_tag_multiple_distinct_values() {
        let header = build_test_header(vec![
            ("RG1", vec![(rg_tag::SAMPLE, "SampleA")]),
            ("RG2", vec![(rg_tag::SAMPLE, "SampleB")]),
        ]);
        assert_eq!(
            collapse_read_group_tag(&header, rg_tag::SAMPLE),
            Some("SampleA,SampleB".to_string()),
        );
    }

    #[test]
    fn test_collapse_read_group_tag_missing_from_some_rgs() {
        let header =
            build_test_header(vec![("RG1", vec![(rg_tag::SAMPLE, "SampleA")]), ("RG2", vec![])]);
        assert_eq!(collapse_read_group_tag(&header, rg_tag::SAMPLE), Some("SampleA".to_string()),);
    }

    #[test]
    fn test_collapse_read_group_tag_all_missing() {
        let header = build_test_header(vec![("RG1", vec![]), ("RG2", vec![])]);
        assert_eq!(collapse_read_group_tag(&header, rg_tag::SAMPLE), None);
    }

    #[test]
    fn test_collapse_read_group_tag_platform_uppercased() {
        let header = build_test_header(vec![
            ("RG1", vec![(rg_tag::PLATFORM, "illumina")]),
            ("RG2", vec![(rg_tag::PLATFORM, "Illumina")]),
        ]);
        // Both should uppercase to "ILLUMINA" and deduplicate
        assert_eq!(
            collapse_read_group_tag(&header, rg_tag::PLATFORM),
            Some("ILLUMINA".to_string()),
        );
    }

    #[test]
    fn test_collapse_read_group_tag_platform_multiple_distinct() {
        let header = build_test_header(vec![
            ("RG1", vec![(rg_tag::PLATFORM, "illumina")]),
            ("RG2", vec![(rg_tag::PLATFORM, "ONT")]),
        ]);
        assert_eq!(
            collapse_read_group_tag(&header, rg_tag::PLATFORM),
            Some("ILLUMINA,ONT".to_string()),
        );
    }

    #[test]
    fn test_create_unmapped_consensus_header_collapses_tags() {
        let input_header = build_test_header(vec![
            (
                "RG1",
                vec![
                    (rg_tag::SAMPLE, "SampleA"),
                    (rg_tag::LIBRARY, "LibA"),
                    (rg_tag::PLATFORM, "illumina"),
                    (rg_tag::PLATFORM_UNIT, "FlowcellA.1"),
                    (rg_tag::SEQUENCING_CENTER, "CenterX"),
                    (rg_tag::DESCRIPTION, "Run 1"),
                ],
            ),
            (
                "RG2",
                vec![
                    (rg_tag::SAMPLE, "SampleA"),
                    (rg_tag::LIBRARY, "LibB"),
                    (rg_tag::PLATFORM, "ILLUMINA"),
                    (rg_tag::PLATFORM_UNIT, "FlowcellB.2"),
                    (rg_tag::SEQUENCING_CENTER, "CenterX"),
                    (rg_tag::DESCRIPTION, "Run 2"),
                ],
            ),
        ]);

        let output_header = create_unmapped_consensus_header(
            &input_header,
            "consensus",
            "Read group",
            "fgumi simplex",
        )
        .expect("create_unmapped_consensus_header should succeed");

        let read_groups = output_header.read_groups();
        assert_eq!(read_groups.len(), 1);
        let rg = read_groups.get(&BString::from("consensus")).expect("consensus RG not found");

        // SM: both have SampleA -> deduplicated
        assert_eq!(
            rg.other_fields().get(&rg_tag::SAMPLE).expect("SM tag not found").to_string(),
            "SampleA"
        );
        // LB: LibA and LibB -> comma-joined
        assert_eq!(
            rg.other_fields().get(&rg_tag::LIBRARY).expect("LB tag not found").to_string(),
            "LibA,LibB"
        );
        // PL: illumina and ILLUMINA -> uppercased and deduped
        assert_eq!(
            rg.other_fields().get(&rg_tag::PLATFORM).expect("PL tag not found").to_string(),
            "ILLUMINA"
        );
        // PU: two distinct values
        assert_eq!(
            rg.other_fields().get(&rg_tag::PLATFORM_UNIT).expect("PU tag not found").to_string(),
            "FlowcellA.1,FlowcellB.2",
        );
        // CN: both CenterX -> deduplicated
        assert_eq!(
            rg.other_fields()
                .get(&rg_tag::SEQUENCING_CENTER)
                .expect("CN tag not found")
                .to_string(),
            "CenterX",
        );
        // DS: two distinct descriptions
        assert_eq!(
            rg.other_fields().get(&rg_tag::DESCRIPTION).expect("DS tag not found").to_string(),
            "Run 1,Run 2",
        );
    }

    #[test]
    fn test_create_unmapped_consensus_header_no_input_rgs() {
        let input_header = Header::builder().build();

        let output_header =
            create_unmapped_consensus_header(&input_header, "A", "Read group", "fgumi simplex")
                .expect("create_unmapped_consensus_header should succeed");

        let read_groups = output_header.read_groups();
        assert_eq!(read_groups.len(), 1);
        let rg = read_groups.get(&BString::from("A")).expect("RG A not found");
        // No tags should be set when there are no input read groups
        assert!(rg.other_fields().get(&rg_tag::SAMPLE).is_none());
        assert!(rg.other_fields().get(&rg_tag::LIBRARY).is_none());
        assert!(rg.other_fields().get(&rg_tag::PLATFORM).is_none());
    }

    #[test]
    fn test_stats_to_metrics_empty() {
        let stats = ConsensusCallingStats::new();
        let metrics = stats.to_metrics();
        assert_eq!(metrics.total_input_reads, 0);
        assert_eq!(metrics.consensus_reads, 0);
        assert_eq!(metrics.filtered_reads, 0);
    }

    #[test]
    fn test_stats_to_metrics_with_values() {
        let mut stats = ConsensusCallingStats::new();
        stats.total_reads = 100;
        stats.consensus_reads = 50;
        stats.filtered_reads = 10;

        let metrics = stats.to_metrics();
        assert_eq!(metrics.total_input_reads, 100);
        assert_eq!(metrics.consensus_reads, 50);
        assert_eq!(metrics.filtered_reads, 10);
    }

    #[test]
    fn test_stats_merge() {
        let mut stats1 = ConsensusCallingStats::new();
        stats1.total_reads = 100;
        stats1.consensus_reads = 50;

        let mut stats2 = ConsensusCallingStats::new();
        stats2.total_reads = 200;
        stats2.consensus_reads = 75;

        stats1.merge(&stats2);

        assert_eq!(stats1.total_reads, 300);
        assert_eq!(stats1.consensus_reads, 125);
    }
}
