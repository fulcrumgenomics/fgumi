//! `SimplexMetrics` command for collecting QC metrics from simplex sequencing data.
//!
//! This module implements single-pass metric collection with deterministic downsampling
//! at 20 levels (5%, 10%, ..., 100%). It produces:
//! - Family size distributions (CS and SS families)
//! - Yield curves at each downsampling fraction
//! - UMI observation frequencies
//! - Optional PDF plots via an embedded R script

use crate::logging::OperationTimer;
use crate::metrics::simplex::{SimplexMetricsCollector, SimplexYieldMetric};
use crate::simple_umi_consensus::SimpleUmiConsensusCaller;
use crate::validation::validate_file_exists;
use anyhow::{Context, Result};
use clap::Parser;
use fgoxide::io::DelimFile;
use log::info;
use std::path::PathBuf;

use super::command::Command;
use super::shared_metrics::{
    DOWNSAMPLING_FRACTIONS, TemplateInfo, compute_template_metadata, execute_r_script,
    is_r_available, parse_intervals, process_templates_from_bam, validate_not_consensus_bam,
};

/// Embedded R script for PDF plot generation (bundled with binary).
const R_SCRIPT: &str = include_str!("../../../resources/CollectSimplexSeqMetrics.R");

/// Collects comprehensive QC metrics for simplex sequencing experiments.
#[derive(Parser, Debug)]
#[command(
    name = "simplex-metrics",
    author,
    version,
    about = "\x1b[38;5;173m[POST-CONSENSUS]\x1b[0m \x1b[36mCollect QC metrics for simplex sequencing data\x1b[0m",
    long_about = r#"
Collects a suite of metrics to QC simplex sequencing data.

## Inputs

The input to this tool must be a BAM file that is either:

1. The exact BAM output by the `group` tool (in the sort-order it was produced in)
2. A BAM file that has MI tags present on all reads (usually set by `group` and has been
   sorted into template-coordinate order

Calculation of metrics may be restricted to a set of regions using the `--intervals` parameter.
This can significantly affect results as off-target reads often have very different properties
than on-target reads due to the lack of enrichment.

## Outputs

The following output files are produced:

1. **<output>.family_sizes.txt**: metrics on the frequency of CS and SS families of different sizes
2. **<output>.simplex_yield_metrics.txt**: summary QC metrics produced using 5%, 10%, 15%...100% of the data
3. **<output>.umi_counts.txt**: metrics on the frequency of observations of UMIs within reads and tag families
4. **<output>.simplex_qc.pdf**: (optional) a series of plots generated from the preceding metrics files for
                               visualization. This file is only produced if R is available with the required
                               packages (ggplot2 and scales). Use `--description` to customize plot titles.

Within the metrics files the prefixes `CS` and `SS` are used to mean:

* **CS**: tag families where membership is defined solely on matching genome coordinates and strand
* **SS**: single-stranded tag families where membership is defined by genome coordinates, strand and UMI
"#
)]
pub struct SimplexMetrics {
    /// Input BAM file (UMI-grouped, from `group`).
    #[arg(short = 'i', long = "input")]
    pub input: PathBuf,

    /// Output prefix for metrics files.
    #[arg(short = 'o', long = "output")]
    pub output: PathBuf,

    /// Minimum reads per SS family to count as a consensus family in yield metrics.
    #[arg(long = "min-reads", default_value = "1")]
    pub min_reads: usize,

    /// Optional intervals file to restrict analysis (BED or Picard interval list format).
    #[arg(short = 'l', long = "intervals")]
    pub intervals: Option<PathBuf>,

    /// Optional sample name or description for PDF plot titles.
    #[arg(long = "description")]
    pub description: Option<String>,
}

impl Command for SimplexMetrics {
    fn execute(&self, _command_line: &str) -> Result<()> {
        info!("SimplexMetrics");
        info!("  Input: {}", self.input.display());
        info!("  Output prefix: {}", self.output.display());
        info!("  Min reads: {}", self.min_reads);

        let timer = OperationTimer::new("Computing simplex metrics");

        // Validate inputs
        validate_file_exists(&self.input, "input BAM file")?;

        // Check that input is not a consensus BAM
        validate_not_consensus_bam(&self.input)?;

        // Load intervals if provided
        let intervals = if let Some(intervals_path) = &self.intervals {
            info!("  Loading intervals from: {}", intervals_path.display());
            let intervals = parse_intervals(intervals_path)?;
            info!("  Loaded {} intervals", intervals.len());
            intervals
        } else {
            Vec::new()
        };

        let fractions = &DOWNSAMPLING_FRACTIONS;

        // Create collectors for each fraction (single-pass architecture)
        let mut collectors: Vec<SimplexMetricsCollector> =
            fractions.iter().map(|_| SimplexMetricsCollector::new()).collect();

        let mut umi_consensus_caller = SimpleUmiConsensusCaller::default();

        info!("Processing templates in single pass at {} sampling fractions...", fractions.len());

        // Single pass: process templates and update all applicable collectors
        let (total_template_count, fraction_template_counts) = process_templates_from_bam(
            &self.input,
            &intervals,
            fractions.len(),
            |group, fraction_counts| {
                Self::process_coordinate_group(
                    group,
                    fractions,
                    &mut collectors,
                    &mut umi_consensus_caller,
                    fraction_counts,
                );
            },
        )?;

        info!("Processed {total_template_count} templates");

        // Generate yield metrics from each collector
        let mut yield_metrics = Vec::new();
        for ((&fraction, collector), &read_pairs) in
            fractions.iter().zip(collectors.iter()).zip(fraction_template_counts.iter())
        {
            let yield_metric =
                Self::generate_yield_metric(collector, fraction, read_pairs, self.min_reads);
            yield_metrics.push(yield_metric);
        }

        // Use the 100% fraction collector for main metrics
        let main_collector =
            collectors.pop().expect("collectors is non-empty (always includes 100% fraction)");

        // Generate and write metrics
        info!("Writing metrics...");

        // Family size metrics
        let family_size_metrics = main_collector.family_size_metrics();
        let family_size_path = format!("{}.family_sizes.txt", self.output.display());
        DelimFile::default()
            .write_tsv(&family_size_path, family_size_metrics)
            .with_context(|| format!("Failed to write family size metrics: {family_size_path}"))?;
        info!("Wrote family size metrics to {family_size_path}");

        // UMI metrics
        let umi_metrics = main_collector.umi_metrics();
        let umi_path = format!("{}.umi_counts.txt", self.output.display());
        DelimFile::default()
            .write_tsv(&umi_path, umi_metrics)
            .with_context(|| format!("Failed to write UMI metrics: {umi_path}"))?;
        info!("Wrote UMI metrics to {umi_path}");

        // Yield metrics
        let yield_path = format!("{}.simplex_yield_metrics.txt", self.output.display());
        DelimFile::default()
            .write_tsv(&yield_path, yield_metrics)
            .with_context(|| format!("Failed to write yield metrics: {yield_path}"))?;
        info!("Wrote yield metrics to {yield_path}");

        // Generate PDF plots using R script (optional)
        let pdf_path = format!("{}.simplex_qc.pdf", self.output.display());
        if is_r_available() {
            let description = self.description.as_deref().unwrap_or("Sample");
            match execute_r_script(
                R_SCRIPT,
                &[&family_size_path, &yield_path, &umi_path, &pdf_path, description],
                "fgumi_CollectSimplexSeqMetrics.R",
            ) {
                Ok(()) => info!("Generated PDF plots: {pdf_path}"),
                Err(e) => {
                    log::warn!("Failed to generate PDF plots: {e}. Continuing without plots.");
                    log::warn!(
                        "To enable PDF generation, ensure R is installed with ggplot2 and scales packages:"
                    );
                    log::warn!("  install.packages(c(\"ggplot2\", \"scales\"))");
                }
            }
        } else {
            log::warn!(
                "R or required packages (ggplot2, scales) not available. Skipping PDF generation."
            );
            log::warn!("To enable PDF generation, install R and required packages:");
            log::warn!("  install.packages(c(\"ggplot2\", \"scales\"))");
        }

        info!("Done!");
        timer.log_completion(total_template_count as u64);
        Ok(())
    }
}

impl SimplexMetrics {
    /// Processes a single coordinate group for all downsampling fractions.
    ///
    /// For each fraction, filters templates by hash, records CS family size (the entire
    /// group), groups by MI tag for SS families, and (at 100% only) collects UMI
    /// observations via consensus calling per UMI position.
    fn process_coordinate_group(
        group: &[TemplateInfo],
        fractions: &[f64],
        collectors: &mut [SimplexMetricsCollector],
        umi_consensus_caller: &mut SimpleUmiConsensusCaller,
        fraction_template_counts: &mut [usize],
    ) {
        use std::collections::HashMap;

        if group.is_empty() {
            return;
        }

        // Pre-compute metadata once for the entire group
        let metadata = compute_template_metadata(group);
        let last_fraction_idx = fractions.len() - 1;

        let mut ss_groups: HashMap<&str, usize> = HashMap::new();

        for (idx, &fraction) in fractions.iter().enumerate() {
            // Filter once per fraction
            let downsampled: Vec<_> =
                metadata.iter().filter(|m| m.template.hash_fraction <= fraction).collect();

            if downsampled.is_empty() {
                continue;
            }

            // CS family size
            fraction_template_counts[idx] += downsampled.len();
            collectors[idx].record_cs_family(downsampled.len());

            // Group by MI tag for SS families
            ss_groups.clear();
            for m in &downsampled {
                *ss_groups.entry(m.template.mi.as_str()).or_default() += 1;
            }
            for &ss_size in ss_groups.values() {
                collectors[idx].record_ss_family(ss_size);
            }

            // UMI metrics only at the 100% fraction (last index)
            if idx == last_fraction_idx {
                // Group by base_umi (MI without strand suffix) and collect RX tags
                let mut umi_groups: HashMap<&str, Vec<&str>> = HashMap::new();
                for m in &downsampled {
                    umi_groups.entry(m.base_umi).or_default().push(m.template.rx.as_str());
                }

                for rx_tags in umi_groups.values() {
                    // For simplex: no strand-swapping. Split each RX by '-' for
                    // multi-component UMIs and call consensus per position.
                    let split_rx: Vec<Vec<&str>> =
                        rx_tags.iter().map(|rx| rx.split('-').collect()).collect();
                    let num_components = split_rx.first().map_or(0, Vec::len);

                    for pos in 0..num_components {
                        let umis_at_pos: Vec<String> = split_rx
                            .iter()
                            .filter_map(|parts| parts.get(pos).map(|s| (*s).to_string()))
                            .filter(|s| !s.is_empty())
                            .collect();

                        if umis_at_pos.is_empty() {
                            continue;
                        }

                        let (consensus, _had_errors) = umi_consensus_caller.consensus(&umis_at_pos);
                        let raw_count = umis_at_pos.len();
                        let error_count = umis_at_pos.iter().filter(|u| **u != consensus).count();
                        collectors[idx].record_umi(&consensus, raw_count, error_count, true);
                    }
                }
            }
        }
    }

    /// Generates a yield metric from a collector at a specific downsampling fraction.
    ///
    /// Computes summary statistics including CS and SS family counts, mean SS family size,
    /// singleton fraction, and number of SS families meeting the minimum read threshold.
    fn generate_yield_metric(
        collector: &SimplexMetricsCollector,
        fraction: f64,
        read_pairs: usize,
        min_reads: usize,
    ) -> SimplexYieldMetric {
        let family_size_metrics = collector.family_size_metrics();

        let cs_families: usize = family_size_metrics.iter().map(|m| m.cs_count).sum();
        let ss_families: usize = family_size_metrics.iter().map(|m| m.ss_count).sum();

        // Total reads in SS families = sum(family_size * ss_count)
        let total_ss_reads: usize =
            family_size_metrics.iter().map(|m| m.family_size * m.ss_count).sum();
        let mean_ss_family_size =
            if ss_families > 0 { total_ss_reads as f64 / ss_families as f64 } else { 0.0 };

        let ss_singletons: usize =
            family_size_metrics.iter().find(|m| m.family_size == 1).map_or(0, |m| m.ss_count);
        let ss_singleton_fraction =
            if ss_families > 0 { ss_singletons as f64 / ss_families as f64 } else { 0.0 };

        let ss_consensus_families: usize = family_size_metrics
            .iter()
            .filter(|m| m.family_size >= min_reads)
            .map(|m| m.ss_count)
            .sum();

        SimplexYieldMetric {
            fraction,
            read_pairs,
            cs_families,
            ss_families,
            mean_ss_family_size,
            ss_singletons,
            ss_singleton_fraction,
            ss_consensus_families,
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::metrics::simplex::{SimplexFamilySizeMetric, SimplexMetricsCollector};
    use anyhow::Result;
    use fgoxide::io::DelimFile;
    use fgumi_raw_bam::{
        SamBuilder as RawSamBuilder, flags, raw_record_to_record_buf, testutil::encode_op,
    };
    use noodles::bam;
    use noodles::sam;
    use noodles::sam::alignment::io::Write;
    use std::num::NonZeroUsize;
    use tempfile::{NamedTempFile, TempDir};

    fn create_test_header() -> sam::Header {
        use noodles::sam::header::record::value::Map;
        use noodles::sam::header::record::value::map::ReferenceSequence;

        let mut builder = sam::Header::builder();
        builder = builder.add_reference_sequence(
            bstr::BString::from("chr1"),
            Map::<ReferenceSequence>::new(
                NonZeroUsize::new(248_956_422).expect("non-zero chromosome length"),
            ),
        );
        builder = builder.add_reference_sequence(
            bstr::BString::from("chr2"),
            Map::<ReferenceSequence>::new(
                NonZeroUsize::new(242_193_529).expect("non-zero chromosome length"),
            ),
        );
        builder.build()
    }

    fn to_record_buf(raw: fgumi_raw_bam::RawRecord) -> sam::alignment::RecordBuf {
        raw_record_to_record_buf(&raw, &sam::Header::default())
            .expect("raw_record_to_record_buf failed in test")
    }

    #[allow(clippy::cast_sign_loss)]
    fn build_test_pair(
        name: &str,
        ref_id: usize,
        pos1: i32,
        pos2: i32,
        rx_umi: &str,
        mi_tag: &str,
    ) -> (sam::alignment::RecordBuf, sam::alignment::RecordBuf) {
        let seq = vec![b'A'; 100];
        let quals = vec![30u8; 100];
        let cigar = encode_op(0, 100); // 100M

        let mut b1 = RawSamBuilder::new();
        b1.read_name(name.as_bytes())
            .flags(flags::PAIRED | flags::FIRST_SEGMENT | flags::MATE_REVERSE)
            .ref_id(ref_id as i32)
            .pos(pos1 - 1)
            .mapq(60)
            .cigar_ops(&[cigar])
            .sequence(&seq)
            .qualities(&quals)
            .mate_ref_id(ref_id as i32)
            .mate_pos(pos2 - 1);
        b1.add_string_tag(b"RX", rx_umi.as_bytes());
        b1.add_string_tag(b"MI", mi_tag.as_bytes());
        let r1 = to_record_buf(b1.build());

        let mut b2 = RawSamBuilder::new();
        b2.read_name(name.as_bytes())
            .flags(flags::PAIRED | flags::LAST_SEGMENT | flags::REVERSE)
            .ref_id(ref_id as i32)
            .pos(pos2 - 1)
            .mapq(60)
            .cigar_ops(&[cigar])
            .sequence(&seq)
            .qualities(&quals)
            .mate_ref_id(ref_id as i32)
            .mate_pos(pos1 - 1);
        b2.add_string_tag(b"RX", rx_umi.as_bytes());
        b2.add_string_tag(b"MI", mi_tag.as_bytes());
        let r2 = to_record_buf(b2.build());

        (r1, r2)
    }

    fn create_test_bam(records: Vec<sam::alignment::RecordBuf>) -> Result<NamedTempFile> {
        let temp_file = NamedTempFile::new()?;
        let header = create_test_header();
        let mut writer = bam::io::writer::Builder.build_from_path(temp_file.path())?;
        writer.write_header(&header)?;
        for record in records {
            writer.write_alignment_record(&header, &record)?;
        }
        drop(writer);
        Ok(temp_file)
    }

    #[test]
    fn test_generate_yield_metric_basic() {
        let mut collector = SimplexMetricsCollector::new();
        collector.record_ss_family(1);
        collector.record_ss_family(2);
        collector.record_ss_family(3);
        collector.record_cs_family(6);

        let metric = SimplexMetrics::generate_yield_metric(&collector, 1.0, 6, 2);

        assert_eq!(metric.cs_families, 1);
        assert_eq!(metric.ss_families, 3);
        assert!((metric.mean_ss_family_size - 2.0).abs() < 0.001);
        assert_eq!(metric.ss_singletons, 1);
        assert!((metric.ss_singleton_fraction - 1.0 / 3.0).abs() < 0.001);
        assert_eq!(metric.ss_consensus_families, 2); // sizes 2 and 3 meet min_reads=2
    }

    #[test]
    fn test_generate_yield_metric_empty() {
        let collector = SimplexMetricsCollector::new();
        let metric = SimplexMetrics::generate_yield_metric(&collector, 0.5, 0, 1);

        assert_eq!(metric.cs_families, 0);
        assert_eq!(metric.ss_families, 0);
        assert!(metric.mean_ss_family_size.abs() < f64::EPSILON);
        assert_eq!(metric.ss_singletons, 0);
        assert_eq!(metric.ss_consensus_families, 0);
    }

    #[test]
    fn test_generate_yield_metric_all_singletons() {
        let mut collector = SimplexMetricsCollector::new();
        collector.record_ss_family(1);
        collector.record_ss_family(1);
        collector.record_ss_family(1);

        let metric = SimplexMetrics::generate_yield_metric(&collector, 1.0, 3, 2);

        assert_eq!(metric.ss_families, 3);
        assert_eq!(metric.ss_singletons, 3);
        assert!((metric.ss_singleton_fraction - 1.0).abs() < f64::EPSILON);
        assert_eq!(metric.ss_consensus_families, 0);
    }

    #[test]
    fn test_generate_accurate_family_size_counts() -> Result<()> {
        let mut records = Vec::new();

        // Coordinate group 1 at pos 100-200: MI "1/A" x2, MI "2/A" x1
        // → 1 CS family size 3, SS families: size 2 (MI 1/A) and size 1 (MI 2/A)
        let (r1, r2) = build_test_pair("grp1_mi1_0", 0, 100, 200, "AAA", "1/A");
        records.push(r1);
        records.push(r2);
        let (r1, r2) = build_test_pair("grp1_mi1_1", 0, 100, 200, "AAA", "1/A");
        records.push(r1);
        records.push(r2);
        let (r1, r2) = build_test_pair("grp1_mi2_0", 0, 100, 200, "CCC", "2/A");
        records.push(r1);
        records.push(r2);

        // Coordinate group 2 at pos 300-400: MI "3/A" x1
        // → 1 CS family size 1, 1 SS family size 1
        let (r1, r2) = build_test_pair("grp2_mi3_0", 0, 300, 400, "GGG", "3/A");
        records.push(r1);
        records.push(r2);

        let input = create_test_bam(records)?;
        let output_dir = TempDir::new()?;
        let output = output_dir.path().join("output");

        let cmd = SimplexMetrics {
            input: input.path().to_path_buf(),
            output: output.clone(),
            min_reads: 1,
            intervals: None,
            description: None,
        };
        cmd.execute("test")?;

        let family_size_path = format!("{}.family_sizes.txt", output.display());
        let metrics: Vec<SimplexFamilySizeMetric> =
            DelimFile::default().read_tsv(&family_size_path)?;

        // Size 1: 1 CS family (group 2), 2 SS families (MI 2/A, MI 3/A)
        let size_1 = metrics
            .iter()
            .find(|m| m.family_size == 1)
            .expect("expected family size 1 metric not found");
        assert_eq!(size_1.cs_count, 1);
        assert_eq!(size_1.ss_count, 2);

        // Size 2: 0 CS, 1 SS (MI 1/A from group 1)
        let size_2 = metrics
            .iter()
            .find(|m| m.family_size == 2)
            .expect("expected family size 2 metric not found");
        assert_eq!(size_2.cs_count, 0);
        assert_eq!(size_2.ss_count, 1);

        // Size 3: 1 CS (group 1), 0 SS
        let size_3 = metrics
            .iter()
            .find(|m| m.family_size == 3)
            .expect("expected family size 3 metric not found");
        assert_eq!(size_3.cs_count, 1);
        assert_eq!(size_3.ss_count, 0);

        Ok(())
    }

    #[test]
    fn test_downsampling_generates_multiple_fractions() -> Result<()> {
        let mut records = Vec::new();

        // Create 50 templates at the same position with unique MIs
        for i in 0..50 {
            let (r1, r2) =
                build_test_pair(&format!("read_{i}"), 0, 100, 200, "AAA", &format!("{i}/A"));
            records.push(r1);
            records.push(r2);
        }

        let input = create_test_bam(records)?;
        let output_dir = TempDir::new()?;
        let output = output_dir.path().join("output");

        let cmd = SimplexMetrics {
            input: input.path().to_path_buf(),
            output: output.clone(),
            min_reads: 1,
            intervals: None,
            description: None,
        };
        cmd.execute("test")?;

        let yield_path = format!("{}.simplex_yield_metrics.txt", output.display());
        let metrics: Vec<SimplexYieldMetric> = DelimFile::default().read_tsv(&yield_path)?;

        assert_eq!(metrics.len(), 20);

        for i in 1..metrics.len() {
            assert!(metrics[i].fraction > metrics[i - 1].fraction);
            assert!(metrics[i].ss_families >= metrics[i - 1].ss_families);
        }

        Ok(())
    }

    #[test]
    fn test_default_simplex_metrics_parameters() {
        let cmd = SimplexMetrics {
            input: PathBuf::from("test.bam"),
            output: PathBuf::from("output"),
            min_reads: 1,
            intervals: None,
            description: None,
        };
        assert_eq!(cmd.min_reads, 1);
    }

    #[test]
    fn test_reject_consensus_bam() -> Result<()> {
        let temp_file = NamedTempFile::new()?;
        let header = create_test_header();
        let mut writer = bam::io::writer::Builder.build_from_path(temp_file.path())?;
        writer.write_header(&header)?;

        let mut b = RawSamBuilder::new();
        b.read_name(b"read1")
            .flags(flags::PAIRED | flags::FIRST_SEGMENT)
            .ref_id(0)
            .pos(99)
            .mapq(60)
            .cigar_ops(&[encode_op(0, 100)])
            .sequence(&[b'A'; 100])
            .qualities(&[30u8; 100]);
        b.add_string_tag(b"cD", b"5");
        let record = to_record_buf(b.build());
        writer.write_alignment_record(&header, &record)?;
        drop(writer);

        let result = validate_not_consensus_bam(temp_file.path());
        assert!(result.is_err());
        assert!(result.unwrap_err().to_string().contains("consensus BAM"));

        Ok(())
    }
}
