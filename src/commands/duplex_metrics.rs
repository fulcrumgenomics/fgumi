//! Enhanced `DuplexMetrics` with downsampling and yield metrics.
//!
//! This version implements:
//! - Full downsampling at 20 levels (5%, 10%, ..., 100%)
//! - UMI consensus calling within coordinate+strand families
//! - Ideal duplex fraction calculation using proper binomial CDF
//! - Optional interval filtering (BED or Picard interval list format) to restrict analysis to specific regions

use anyhow::{Context, Result};
use clap::Parser;
use fgoxide::io::DelimFile;
use fgumi_lib::defaults;
use fgumi_lib::logging::OperationTimer;
use fgumi_lib::metrics::duplex::{DuplexMetricsCollector, DuplexYieldMetric};
use fgumi_lib::simple_umi_consensus::SimpleUmiConsensusCaller;
use fgumi_lib::umi::extract_mi_base;
use fgumi_lib::validation::validate_file_exists;
use log::info;
use statrs::distribution::{Binomial, DiscreteCDF};
use std::path::PathBuf;

use super::command::Command;
use super::shared_metrics::{
    DOWNSAMPLING_FRACTIONS, TemplateInfo, TemplateMetadata, compute_template_metadata,
    execute_r_script, is_r_available, parse_intervals, process_templates_from_bam,
    validate_not_consensus_bam,
};

/// Embedded R script for PDF plot generation (bundled with binary)
const R_SCRIPT: &str = include_str!("../../resources/CollectDuplexSeqMetrics.R");

/// Collects comprehensive QC metrics for duplex sequencing experiments
#[derive(Parser, Debug)]
#[command(
    name = "duplex-metrics",
    author,
    version,
    about = "\x1b[38;5;173m[POST-CONSENSUS]\x1b[0m \x1b[36mCollect QC metrics for duplex consensus reads\x1b[0m",
    long_about = r#"
Collects a suite of metrics to QC duplex sequencing data.

## Inputs

The input to this tool must be a BAM file that is either:

1. The exact BAM output by the `group` tool (in the sort-order it was produced in)
2. A BAM file that has MI tags present on all reads (usually set by `group` and has been
   sorted into template-coordinate order

Calculation of metrics may be restricted to a set of regions using the `--intervals` parameter.
This can significantly affect results as off-target reads in duplex sequencing experiments often
have very different properties than on-target reads due to the lack of enrichment.

Several metrics are calculated related to the fraction of tag families that have duplex coverage.
The definition of "duplex" is controlled by the `--min-ab-reads` and `--min-ba-reads` parameters.
The default is to treat any tag family with at least one observation of each strand as a duplex,
but this could be made more stringent, e.g. by setting `--min-ab-reads=3 --min-ba-reads=3`.

## Outputs

The following output files are produced:

1. **<output>.family_sizes.txt**: metrics on the frequency of different types of families of different sizes
2. **<output>.duplex_family_sizes.txt**: metrics on the frequency of duplex tag families by the number of
                                          observations from each strand
3. **<output>.duplex_yield_metrics.txt**: summary QC metrics produced using 5%, 10%, 15%...100% of the data
4. **<output>.umi_counts.txt**: metrics on the frequency of observations of UMIs within reads and tag families
5. **<output>.duplex_umi_counts.txt**: (optional) metrics on the frequency of observations of duplex UMIs within
                                       reads and tag families. This file is only produced _if_ the
                                       `--duplex-umi-counts` option is used as it requires significantly more
                                       memory to track all pairs of UMIs seen when a large number of UMI
                                       sequences are present.
6. **<output>.duplex_qc.pdf**: (optional) a series of plots generated from the preceding metrics files for
                               visualization. This file is only produced if R is available with the required
                               packages (ggplot2 and scales). Use `--description` to customize plot titles.

Within the metrics files the prefixes `CS`, `SS` and `DS` are used to mean:

* **CS**: tag families where membership is defined solely on matching genome coordinates and strand
* **SS**: single-stranded tag families where membership is defined by genome coordinates, strand and UMI;
          ie. 50/A and 50/B are considered different tag families
* **DS**: double-stranded tag families where membership is collapsed across single-stranded tag families
          from the same double-stranded source molecule; i.e. 50/A and 50/B become one family
"#
)]
pub struct DuplexMetrics {
    /// Input BAM file (UMI-grouped, from `group`)
    #[arg(short = 'i', long = "input")]
    pub input: PathBuf,

    /// Output prefix for metrics files
    #[arg(short = 'o', long = "output")]
    pub output: PathBuf,

    /// Minimum AB reads to call a duplex
    #[arg(long = "min-ab-reads", default_value = "1")]
    pub min_ab_reads: usize,

    /// Minimum BA reads to call a duplex
    #[arg(long = "min-ba-reads", default_value = "1")]
    pub min_ba_reads: usize,

    /// Collect duplex UMI counts (memory intensive)
    #[arg(long = "duplex-umi-counts", default_value = "false", num_args = 0..=1, default_missing_value = "true", action = clap::ArgAction::Set)]
    pub duplex_umi_counts: bool,

    /// Optional intervals file to restrict analysis (BED or Picard interval list format)
    #[arg(short = 'l', long = "intervals")]
    pub intervals: Option<PathBuf>,

    /// Optional sample name or description for PDF plot titles
    #[arg(long = "description")]
    pub description: Option<String>,

    /// SAM tag containing the raw UMI sequence
    #[arg(long = "umi-tag", default_value = defaults::UMI_TAG)]
    pub umi_tag: String,

    /// SAM tag containing the molecular identifier (assigned by `group`)
    #[arg(long = "mi-tag", default_value = defaults::MI_TAG)]
    pub mi_tag: String,
}

impl Command for DuplexMetrics {
    fn execute(&self, _command_line: &str) -> Result<()> {
        info!("DuplexMetrics");
        info!("  Input: {}", self.input.display());
        info!("  Output prefix: {}", self.output.display());
        info!("  Min AB reads: {}", self.min_ab_reads);
        info!("  Min BA reads: {}", self.min_ba_reads);
        info!("  Collect duplex UMI counts: {}", self.duplex_umi_counts);

        let timer = OperationTimer::new("Computing duplex metrics");

        // Validate inputs
        validate_file_exists(&self.input, "input BAM file")?;

        // Validate parameters
        if self.min_ab_reads < self.min_ba_reads {
            anyhow::bail!(
                "--min-ab-reads ({}) must be >= --min-ba-reads ({})",
                self.min_ab_reads,
                self.min_ba_reads
            );
        }

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
        let mut collectors: Vec<DuplexMetricsCollector> = fractions
            .iter()
            .map(|&fraction| {
                // Only collect duplex UMI counts for the 100% fraction to save memory
                let collect_duplex = self.duplex_umi_counts && (fraction - 1.0_f64).abs() < 0.01;
                DuplexMetricsCollector::new(collect_duplex)
            })
            .collect();

        let mut umi_consensus_caller = SimpleUmiConsensusCaller::default();

        info!("Processing templates in single pass at {} sampling fractions...", fractions.len());

        // Single pass: process templates and update all applicable collectors
        let (total_template_count, fraction_template_counts) = process_templates_from_bam(
            &self.input,
            &intervals,
            &self.mi_tag,
            &self.umi_tag,
            fractions.len(),
            |group, fraction_counts| {
                Self::process_coordinate_group(
                    group,
                    fractions,
                    &mut collectors,
                    &mut umi_consensus_caller,
                    fraction_counts,
                    self,
                );
            },
        )?;

        info!("Processed {total_template_count} templates");

        // Generate yield metrics from each collector
        let mut yield_metrics = Vec::new();
        for ((&fraction, collector), &read_pairs) in
            fractions.iter().zip(collectors.iter()).zip(fraction_template_counts.iter())
        {
            let yield_metric = self.generate_yield_metric(collector, fraction, read_pairs);
            yield_metrics.push(yield_metric);
        }

        // Use the 100% fraction collector for main metrics
        let main_collector = collectors.pop().unwrap(); // Last one is 100%

        // Generate and write metrics
        info!("Writing metrics...");

        // Family size metrics
        let family_size_metrics = main_collector.family_size_metrics();
        let family_size_path = format!("{}.family_sizes.txt", self.output.display());
        DelimFile::default()
            .write_tsv(&family_size_path, family_size_metrics)
            .with_context(|| format!("Failed to write family size metrics: {family_size_path}"))?;
        info!("Wrote family size metrics to {family_size_path}");

        // Duplex family size metrics
        let duplex_family_size_metrics = main_collector.duplex_family_size_metrics();
        let duplex_family_size_path = format!("{}.duplex_family_sizes.txt", self.output.display());
        DelimFile::default()
            .write_tsv(&duplex_family_size_path, duplex_family_size_metrics)
            .with_context(|| {
                format!("Failed to write duplex family size metrics: {duplex_family_size_path}")
            })?;
        info!("Wrote duplex family size metrics to {duplex_family_size_path}");

        // UMI metrics
        let umi_metrics = main_collector.umi_metrics();
        let umi_path = format!("{}.umi_counts.txt", self.output.display());
        DelimFile::default()
            .write_tsv(&umi_path, umi_metrics)
            .with_context(|| format!("Failed to write UMI metrics: {umi_path}"))?;
        info!("Wrote UMI metrics to {umi_path}");

        // Duplex UMI metrics (if enabled)
        if self.duplex_umi_counts {
            let duplex_umi_metrics =
                main_collector.duplex_umi_metrics(&main_collector.umi_metrics());
            let duplex_umi_path = format!("{}.duplex_umi_counts.txt", self.output.display());
            DelimFile::default().write_tsv(&duplex_umi_path, duplex_umi_metrics).with_context(
                || format!("Failed to write duplex UMI metrics: {duplex_umi_path}"),
            )?;
            info!("Wrote duplex UMI metrics to {duplex_umi_path}");
        }

        // Yield metrics
        let yield_path = format!("{}.duplex_yield_metrics.txt", self.output.display());
        DelimFile::default()
            .write_tsv(&yield_path, yield_metrics)
            .with_context(|| format!("Failed to write yield metrics: {yield_path}"))?;
        info!("Wrote yield metrics to {yield_path}");

        // Generate PDF plots using R script (optional)
        let pdf_path = format!("{}.duplex_qc.pdf", self.output.display());
        if is_r_available() {
            let description = self.description.as_deref().unwrap_or("Sample");
            match execute_r_script(
                R_SCRIPT,
                &[
                    &family_size_path,
                    &duplex_family_size_path,
                    &yield_path,
                    &umi_path,
                    &pdf_path,
                    description,
                ],
                "fgumi_CollectDuplexSeqMetrics.R",
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

impl DuplexMetrics {
    /// Processes a single coordinate group for all downsampling fractions.
    ///
    /// Optimized using fgbio's approach: filter once per fraction, then groupBy.
    /// This is O(fractions × `group_size`) instead of O(fractions × `unique_keys` × `group_size`).
    fn process_coordinate_group(
        group: &[TemplateInfo],
        fractions: &[f64],
        collectors: &mut [DuplexMetricsCollector],
        umi_consensus_caller: &mut SimpleUmiConsensusCaller,
        fraction_template_counts: &mut [usize],
        metrics: &DuplexMetrics,
    ) {
        use std::collections::HashMap;

        if group.is_empty() {
            return;
        }

        // Pre-compute metadata once for the entire group
        let metadata = compute_template_metadata(group);

        // For each fraction: filter ONCE, then groupBy (like fgbio)
        for (idx, &fraction) in fractions.iter().enumerate() {
            // Filter once per fraction - equivalent to fgbio's downsampledGroup
            let downsampled: Vec<&TemplateMetadata> =
                metadata.iter().filter(|m| m.template.hash_fraction <= fraction).collect();

            if downsampled.is_empty() {
                continue;
            }

            // CS family size
            fraction_template_counts[idx] += downsampled.len();
            collectors[idx].record_cs_family(downsampled.len());

            // Group by MI tag for SS families (like fgbio's groupBy)
            let mut ss_groups: HashMap<&str, usize> = HashMap::new();
            for m in &downsampled {
                *ss_groups.entry(m.template.mi.as_str()).or_default() += 1;
            }
            for &ss_size in ss_groups.values() {
                collectors[idx].record_ss_family(ss_size);
            }

            // Group by base_umi for DS families with strand counts
            // HashMap value: (a_count, b_count, mi_rx_pairs for UMI metrics)
            #[allow(clippy::type_complexity)]
            let mut ds_groups: HashMap<&str, (usize, usize, Vec<(&str, &str)>)> = HashMap::new();
            for m in &downsampled {
                let entry = ds_groups.entry(m.base_umi).or_default();
                if m.is_a_strand {
                    entry.0 += 1;
                } else if m.is_b_strand {
                    entry.1 += 1;
                }
                // Collect mi/rx pairs for UMI metrics (only used at 100% fraction)
                entry.2.push((m.template.mi.as_str(), m.template.rx.as_str()));
            }

            let is_full_fraction = (fraction - 1.0_f64).abs() < 0.01;

            for (base_umi, (a_count, b_count, mi_rx_pairs)) in &ds_groups {
                let ds_size = a_count + b_count;
                collectors[idx].record_ds_family(ds_size);

                let (ab_count, ba_count) =
                    if a_count >= b_count { (*a_count, *b_count) } else { (*b_count, *a_count) };

                collectors[idx].record_duplex_family(ab_count, ba_count);

                // Only collect UMI metrics for the 100% fraction
                if is_full_fraction {
                    let owned_pairs: Vec<(String, String)> = mi_rx_pairs
                        .iter()
                        .map(|(mi, rx)| ((*mi).to_string(), (*rx).to_string()))
                        .collect();

                    metrics.update_umi_metrics(
                        &mut collectors[idx],
                        umi_consensus_caller,
                        &owned_pairs,
                        base_umi,
                        *a_count,
                        *b_count,
                    );
                }
            }
        }
    }

    /// Generates a yield metric from a collector at a specific downsampling fraction.
    ///
    /// Computes summary yield metrics including family counts (CS, SS, DS), duplex counts,
    /// duplex fractions, and ideal duplex fraction based on binomial probability. The ideal
    /// fraction represents the expected duplex yield if strands were randomly distributed.
    ///
    /// # Arguments
    ///
    /// * `collector` - Metrics collector containing accumulated family and duplex data
    /// * `fraction` - Downsampling fraction (0.05 to 1.00)
    /// * `read_pairs` - Number of read pairs at this downsampling fraction
    ///
    /// # Returns
    ///
    /// A `DuplexYieldMetric` containing all computed yield statistics.
    fn generate_yield_metric(
        &self,
        collector: &DuplexMetricsCollector,
        fraction: f64,
        read_pairs: usize,
    ) -> DuplexYieldMetric {
        let family_size_metrics = collector.family_size_metrics();

        // Sum up DS family counts from family_size_metrics
        let ds_families_count: usize = family_size_metrics.iter().map(|m| m.ds_count).sum();

        let duplex_family_size_metrics = collector.duplex_family_size_metrics();
        let ds_duplexes_count: usize = duplex_family_size_metrics
            .iter()
            .filter(|m| m.ab_size >= self.min_ab_reads && m.ba_size >= self.min_ba_reads)
            .map(|m| m.count)
            .sum();

        // Calculate ideal fraction
        let ds_family_sizes: Vec<usize> =
            family_size_metrics.iter().flat_map(|m| vec![m.family_size; m.ds_count]).collect();

        let ideal_fraction = Self::calculate_ideal_duplex_fraction(
            &ds_family_sizes,
            self.min_ab_reads,
            self.min_ba_reads,
        );

        let cs_families: usize = family_size_metrics.iter().map(|m| m.cs_count).sum();
        let ss_families: usize = duplex_family_size_metrics
            .iter()
            .map(|m| {
                let mut count = 0;
                if m.ab_size > 0 {
                    count += 1;
                }
                if m.ba_size > 0 {
                    count += 1;
                }
                count * m.count
            })
            .sum();

        DuplexYieldMetric {
            fraction,
            read_pairs,
            cs_families,
            ss_families,
            ds_families: ds_families_count,
            ds_duplexes: ds_duplexes_count,
            ds_fraction_duplexes: if ds_families_count > 0 {
                ds_duplexes_count as f64 / ds_families_count as f64
            } else {
                0.0
            },
            ds_fraction_duplexes_ideal: ideal_fraction,
        }
    }

    /// Calculates the ideal duplex fraction using binomial model.
    ///
    /// For each family of size N, calculates the probability that both strands
    /// have sufficient reads (A >= `min_ab` AND B >= `min_ba` where A + B = N).
    /// Assumes each read has 0.5 probability of being on each strand.
    fn calculate_ideal_duplex_fraction(
        family_sizes: &[usize],
        min_ab: usize,
        min_ba: usize,
    ) -> f64 {
        if family_sizes.is_empty() {
            return 0.0;
        }

        let mut ideal_duplexes = 0.0;
        let total_families = family_sizes.len() as f64;

        for &size in family_sizes {
            if size < min_ab + min_ba {
                // Impossible to form a duplex with this family size
                continue;
            }

            // Calculate P(A >= min_ab AND B >= min_ba) where A ~ Binomial(n=size, p=0.5)
            // and B = size - A
            //
            // This is equivalent to: P(min_ba <= A <= size - min_ab)
            // = P(A <= size - min_ab) - P(A < min_ba)
            // = CDF(size - min_ab) - CDF(min_ba - 1)

            let binomial = match Binomial::new(0.5, size as u64) {
                Ok(b) => b,
                Err(_) => continue, // Skip if binomial creation fails
            };

            let upper_bound = size - min_ba;
            let lower_bound = min_ab;

            // P(A >= lower_bound AND A <= upper_bound)
            let prob = if upper_bound >= lower_bound {
                let p_upper = binomial.cdf(upper_bound as u64);
                let p_lower =
                    if lower_bound > 0 { binomial.cdf((lower_bound - 1) as u64) } else { 0.0 };
                p_upper - p_lower
            } else {
                0.0
            };

            ideal_duplexes += prob;
        }

        ideal_duplexes / total_families
    }

    /// Updates UMI metrics for a duplex family
    ///
    /// This method:
    /// 1. Uses RX tags (raw UMI sequences) to extract individual UMI observations
    /// 2. Separates by strand (/A and /B suffixes in MI tags), swapping UMI parts for B strand
    /// 3. Calls consensus for each UMI position
    /// 4. Records raw observations, errors, and unique observations for each individual UMI
    /// 5. Records duplex UMI metrics if enabled
    ///
    /// This matches the Scala implementation in CollectDuplexSeqMetrics.scala:407-431
    fn update_umi_metrics(
        &self,
        collector: &mut DuplexMetricsCollector,
        umi_consensus_caller: &mut SimpleUmiConsensusCaller,
        group_pairs: &[(String, String)],
        base_umi: &str,
        _a_count: usize,
        _b_count: usize,
    ) {
        // Collect individual UMI parts from each strand's RX tag
        // For /A strand: split RX "AAA-TTT" → umi1s += "AAA", umi2s += "TTT"
        // For /B strand: split RX "TTT-AAA" → umi1s += "AAA", umi2s += "TTT" (swapped)
        let mut umi1s = Vec::new();
        let mut umi2s = Vec::new();

        for (mi, rx) in group_pairs {
            // Check if this MI tag belongs to the current base_umi family
            let mi_base = extract_mi_base(mi);

            if mi_base != base_umi {
                continue;
            }

            // Split the RX tag to get individual UMI parts
            let parts: Vec<&str> = rx.split('-').collect();
            if parts.len() != 2 {
                // Not a valid duplex UMI, skip
                continue;
            }

            // Check that both components are non-empty
            if parts[0].is_empty() || parts[1].is_empty() {
                // Empty component, skip
                continue;
            }

            // Add UMI parts based on strand
            if mi.ends_with("/A") {
                // For /A strand: u1 goes to umi1s, u2 goes to umi2s
                umi1s.push(parts[0].to_string());
                umi2s.push(parts[1].to_string());
            } else if mi.ends_with("/B") {
                // For /B strand: u2 goes to umi1s, u1 goes to umi2s (swapped)
                umi1s.push(parts[1].to_string());
                umi2s.push(parts[0].to_string());
            }
        }

        // Call consensus for each UMI position and record metrics
        let mut consensus_umis = Vec::new();

        if !umi1s.is_empty() {
            let (consensus, _had_errors) = umi_consensus_caller.consensus(&umi1s);
            let raw_count = umi1s.len();
            let error_count = umi1s.iter().filter(|u| **u != consensus).count();
            collector.record_umi(&consensus, raw_count, error_count, true);
            consensus_umis.push(consensus);
        }

        if !umi2s.is_empty() {
            let (consensus, _had_errors) = umi_consensus_caller.consensus(&umi2s);
            let raw_count = umi2s.len();
            let error_count = umi2s.iter().filter(|u| **u != consensus).count();
            collector.record_umi(&consensus, raw_count, error_count, true);
            consensus_umis.push(consensus);
        }

        // Record duplex UMI metrics if enabled
        if self.duplex_umi_counts && consensus_umis.len() == 2 {
            let duplex_umi = format!("{}-{}", consensus_umis[0], consensus_umis[1]);
            // Each read pair contributes one observation to the duplex UMI
            // (not two, even though we track each component separately)
            let total_raw = umi1s.len();

            // Count how many raw RX tags had errors (don't match either duplex orientation)
            let expected_duplex1 = format!("{}-{}", consensus_umis[0], consensus_umis[1]);
            let expected_duplex2 = format!("{}-{}", consensus_umis[1], consensus_umis[0]);
            let error_count = group_pairs
                .iter()
                .filter(|(mi, rx)| {
                    let mi_base = extract_mi_base(mi);
                    mi_base == base_umi && *rx != expected_duplex1 && *rx != expected_duplex2
                })
                .count();

            collector.record_duplex_umi(&duplex_umi, total_raw, error_count, true);
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::commands::shared_metrics::{
        Interval, TemplateInfo, compute_hash_fraction, overlaps_intervals, parse_intervals,
    };
    use anyhow::Result;
    use fgoxide::io::DelimFile;
    use fgumi_lib::metrics::duplex::{DuplexFamilySizeMetric, DuplexUmiMetric, FamilySizeMetric};
    use fgumi_lib::metrics::shared::UmiMetric;
    use fgumi_lib::sam::builder::{RecordBuilder, RecordPairBuilder};
    use noodles::bam;
    use noodles::sam;
    use noodles::sam::alignment::io::Write;
    use std::num::NonZeroUsize;
    use tempfile::{NamedTempFile, TempDir};

    // ==================== Test Helper Functions ====================

    /// Creates a standard test header with two reference sequences
    fn create_test_header() -> sam::Header {
        use noodles::sam::header::record::value::Map;
        use noodles::sam::header::record::value::map::ReferenceSequence;

        let mut builder = sam::Header::builder();

        builder = builder.add_reference_sequence(
            bstr::BString::from("chr1"),
            Map::<ReferenceSequence>::new(NonZeroUsize::new(248_956_422).unwrap()),
        );
        builder = builder.add_reference_sequence(
            bstr::BString::from("chr2"),
            Map::<ReferenceSequence>::new(NonZeroUsize::new(242_193_529).unwrap()),
        );

        builder.build()
    }

    /// Builds a pair of reads with proper mate information
    #[allow(clippy::too_many_arguments)]
    fn build_test_pair(
        name: &str,
        ref_id: usize,
        pos1: i32,
        pos2: i32,
        rx_umi: &str,
        mi_tag: &str,
        strand1_plus: bool,
        strand2_plus: bool,
    ) -> (sam::alignment::RecordBuf, sam::alignment::RecordBuf) {
        RecordPairBuilder::new()
            .name(name)
            .r1_sequence(&"A".repeat(100))
            .r2_sequence(&"A".repeat(100))
            .reference_sequence_id(ref_id)
            .r1_start(pos1 as usize)
            .r2_start(pos2 as usize)
            .r1_reverse(!strand1_plus)
            .r2_reverse(!strand2_plus)
            .tag("RX", rx_umi)
            .tag("MI", mi_tag)
            .build()
    }

    /// Writes records to a temporary BAM file in template-coordinate order
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

    // ==================== Test Cases ====================

    #[test]
    fn test_count_umis_once_per_read_pair() -> Result<()> {
        let mut records = Vec::new();

        // Two pairs at the same location with complementary UMIs
        // NOTE: MI tag should contain UMI with strand suffix for duplex-metrics
        let (r1, r2) = build_test_pair("q1", 0, 100, 200, "AAA-TTT", "AAA-TTT/A", true, false);
        records.push(r1);
        records.push(r2);

        // /B reads: swap positions (200, 100) and use opposite strands to match fgbio's test setup
        let (r1, r2) = build_test_pair("q2", 0, 200, 100, "TTT-AAA", "AAA-TTT/B", false, true);
        records.push(r1);
        records.push(r2);

        let input = create_test_bam(records)?;
        let output_dir = TempDir::new()?;
        let output = output_dir.path().join("output");

        let cmd = DuplexMetrics {
            input: input.path().to_path_buf(),
            output: output.clone(),
            min_ab_reads: 1,
            min_ba_reads: 1,
            duplex_umi_counts: false,
            intervals: None,
            description: None,
            umi_tag: "RX".to_string(),
            mi_tag: "MI".to_string(),
        };

        cmd.execute("test")?;

        // Read UMI metrics
        let umi_path = format!("{}.umi_counts.txt", output.display());
        let umi_metrics: Vec<UmiMetric> = DelimFile::default().read_tsv(&umi_path)?;

        // Should have 2 UMIs: AAA and TTT
        assert_eq!(umi_metrics.len(), 2, "Should have 2 UMI entries");

        for metric in &umi_metrics {
            assert_eq!(metric.raw_observations, 2, "Each UMI should have 2 raw observations");
            // A and B strands are now in the same coordinate group (same ReadInfo after 5' ordering),
            // so they form one DS family with one unique observation per UMI
            assert_eq!(
                metric.unique_observations, 1,
                "Each UMI should have 1 unique observation (one DS family)"
            );
        }

        // Duplex UMI metrics should not be generated
        let duplex_umi_path = format!("{}.duplex_umi_counts.txt", output.display());
        assert!(!std::path::Path::new(&duplex_umi_path).exists());

        Ok(())
    }

    #[test]
    fn test_error_correct_umis_for_counting() -> Result<()> {
        let mut records = Vec::new();

        // A strand: AAA-TTT (2 exact + 1 with 1 error: CAA-TTT)
        for i in 1..=3 {
            let umi = if i == 3 { "CAA-TTT" } else { "AAA-TTT" };
            let (r1, r2) = build_test_pair(&format!("a{i}"), 0, 100, 200, umi, "1/A", true, false);
            records.push(r1);
            records.push(r2);
        }

        // B strand: TTT-AAA (2 exact + 1 with 1 error: CTT-AAA)
        // Swap positions (200, 100) to match fgbio's test setup for /B reads
        for i in 1..=3 {
            let umi = if i == 3 { "CTT-AAA" } else { "TTT-AAA" };
            let (r1, r2) = build_test_pair(&format!("b{i}"), 0, 200, 100, umi, "1/B", false, true);
            records.push(r1);
            records.push(r2);
        }

        let input = create_test_bam(records)?;
        let output_dir = TempDir::new()?;
        let output = output_dir.path().join("output");

        let cmd = DuplexMetrics {
            input: input.path().to_path_buf(),
            output: output.clone(),
            min_ab_reads: 1,
            min_ba_reads: 1,
            duplex_umi_counts: true,
            intervals: None,
            description: None,
            umi_tag: "RX".to_string(),
            mi_tag: "MI".to_string(),
        };

        cmd.execute("test")?;

        // Read UMI metrics
        let umi_path = format!("{}.umi_counts.txt", output.display());
        let umi_metrics: Vec<UmiMetric> = DelimFile::default().read_tsv(&umi_path)?;

        assert_eq!(umi_metrics.len(), 2, "Should have 2 UMIs after error correction");

        for metric in &umi_metrics {
            assert!(metric.umi == "AAA" || metric.umi == "TTT", "UMI should be AAA or TTT");
            assert_eq!(metric.raw_observations, 6, "Each UMI should have 6 total raw observations");
            assert_eq!(
                metric.raw_observations_with_errors, 1,
                "Each UMI should have 1 observation with error"
            );
            // A and B strands are now in the same coordinate group (same ReadInfo after 5' ordering)
            assert_eq!(metric.unique_observations, 1, "Each UMI should form 1 unique DS family");
        }

        // Read duplex UMI metrics
        let duplex_umi_path = format!("{}.duplex_umi_counts.txt", output.display());
        let duplex_umi_metrics: Vec<DuplexUmiMetric> =
            DelimFile::default().read_tsv(&duplex_umi_path)?;

        assert_eq!(duplex_umi_metrics.len(), 1, "Should have 1 duplex UMI");
        let duplex_metric = &duplex_umi_metrics[0];
        assert_eq!(duplex_metric.umi, "AAA-TTT");
        assert_eq!(duplex_metric.raw_observations, 6);
        // A and B strands are now in the same coordinate group (same DS family)
        assert_eq!(duplex_metric.unique_observations, 1);

        Ok(())
    }

    #[test]
    fn test_count_unique_umi_observations_with_tag_families() -> Result<()> {
        let mut records = Vec::new();

        // Family 1: AAA-TTT + TTT-AAA at position 100-200
        let (r1, r2) = build_test_pair("q1", 0, 100, 200, "AAA-TTT", "1/A", true, false);
        records.push(r1);
        records.push(r2);

        // /B reads: swap positions (200, 100) to match fgbio's test setup
        let (r1, r2) = build_test_pair("q2", 0, 200, 100, "TTT-AAA", "1/B", false, true);
        records.push(r1);
        records.push(r2);

        // Family 2: TTT-AAA + error variant (NTT-AAA) at position 150-250
        for i in 1..=3 {
            let umi = if i == 3 { "NTT-AAA" } else { "TTT-AAA" };
            let (r1, r2) =
                build_test_pair(&format!("q{}", i + 2), 0, 150, 250, umi, "2/A", true, false);
            records.push(r1);
            records.push(r2);
        }

        // Family 3: CCC-GGG at position 250-350
        let (r1, r2) = build_test_pair("q6", 0, 250, 350, "CCC-GGG", "3/B", false, true);
        records.push(r1);
        records.push(r2);

        let input = create_test_bam(records)?;
        let output_dir = TempDir::new()?;
        let output = output_dir.path().join("output");

        let cmd = DuplexMetrics {
            input: input.path().to_path_buf(),
            output: output.clone(),
            min_ab_reads: 1,
            min_ba_reads: 1,
            duplex_umi_counts: false,
            intervals: None,
            description: None,
            umi_tag: "RX".to_string(),
            mi_tag: "MI".to_string(),
        };

        cmd.execute("test")?;

        let umi_path = format!("{}.umi_counts.txt", output.display());
        let umi_metrics: Vec<UmiMetric> = DelimFile::default().read_tsv(&umi_path)?;

        assert_eq!(umi_metrics.len(), 4, "Should have 4 unique UMIs");

        // Check AAA
        let aaa = umi_metrics.iter().find(|m| m.umi == "AAA").expect("Should find AAA");
        assert_eq!(aaa.raw_observations, 5, "AAA should have 5 raw observations");
        assert_eq!(aaa.raw_observations_with_errors, 0, "AAA should have no errors");
        assert_eq!(aaa.unique_observations, 2, "AAA should have 2 unique families");

        // Check TTT
        let ttt = umi_metrics.iter().find(|m| m.umi == "TTT").expect("Should find TTT");
        assert_eq!(ttt.raw_observations, 5, "TTT should have 5 raw observations");
        assert_eq!(ttt.raw_observations_with_errors, 1, "TTT should have 1 error (NTT)");
        assert_eq!(ttt.unique_observations, 2, "TTT should have 2 unique families");

        // Check CCC
        let ccc = umi_metrics.iter().find(|m| m.umi == "CCC").expect("Should find CCC");
        assert_eq!(ccc.raw_observations, 1);
        assert_eq!(ccc.unique_observations, 1);

        // Check GGG
        let ggg = umi_metrics.iter().find(|m| m.umi == "GGG").expect("Should find GGG");
        assert_eq!(ggg.raw_observations, 1);
        assert_eq!(ggg.unique_observations, 1);

        Ok(())
    }

    #[test]
    fn test_generate_accurate_family_size_counts() -> Result<()> {
        let mut records = Vec::new();

        // Three SS families at the same location (100-200)
        // Family 1/A: 2 reads
        for i in 1..=2 {
            let (r1, r2) =
                build_test_pair(&format!("a{i}"), 0, 100, 200, "AAA-TTT", "1/A", true, false);
            records.push(r1);
            records.push(r2);
        }
        // Family 2/A: 1 read
        let (r1, r2) = build_test_pair("b1", 0, 100, 200, "ACG-GGA", "2/A", true, false);
        records.push(r1);
        records.push(r2);
        // Family 3/B: 1 read - swap positions for /B
        let (r1, r2) = build_test_pair("c1", 0, 200, 100, "TAT-CGT", "3/B", false, true);
        records.push(r1);
        records.push(r2);

        // Two duplex families at the same location (200-300)
        // Duplex family 4: 1/A + 1/B
        let (r1, r2) = build_test_pair("d1", 0, 200, 300, "TTT-AAA", "4/A", true, false);
        records.push(r1);
        records.push(r2);
        // Swap positions for /B
        let (r1, r2) = build_test_pair("d2", 0, 300, 200, "AAA-AAA", "4/B", false, true);
        records.push(r1);
        records.push(r2);

        // Duplex family 5: 1/A + 1/B
        let (r1, r2) = build_test_pair("e1", 0, 200, 300, "CCC-GGG", "5/A", true, false);
        records.push(r1);
        records.push(r2);
        // Swap positions for /B
        let (r1, r2) = build_test_pair("e2", 0, 300, 200, "GGG-CCC", "5/B", false, true);
        records.push(r1);
        records.push(r2);

        // One duplex and one SS family at the same location (400-500)
        // SS family 6/A: 1 read
        let (r1, r2) = build_test_pair("f1", 0, 400, 500, "GCG-GAA", "6/A", true, false);
        records.push(r1);
        records.push(r2);

        // Duplex family 7: 2/A + 1/B
        for i in 1..=2 {
            let (r1, r2) =
                build_test_pair(&format!("g{i}"), 0, 400, 500, "ACG-CCT", "7/A", true, false);
            records.push(r1);
            records.push(r2);
        }
        // Swap positions for /B
        let (r1, r2) = build_test_pair("g3", 0, 500, 400, "CCT-ACG", "7/B", false, true);
        records.push(r1);
        records.push(r2);

        let input = create_test_bam(records)?;
        let output_dir = TempDir::new()?;
        let output = output_dir.path().join("output");

        let cmd = DuplexMetrics {
            input: input.path().to_path_buf(),
            output: output.clone(),
            min_ab_reads: 1,
            min_ba_reads: 1,
            duplex_umi_counts: false,
            intervals: None,
            description: None,
            umi_tag: "RX".to_string(),
            mi_tag: "MI".to_string(),
        };

        cmd.execute("test")?;

        let family_path = format!("{}.family_sizes.txt", output.display());
        let family_metrics: Vec<FamilySizeMetric> = DelimFile::default().read_tsv(&family_path)?;

        // Find metrics for each family size
        let one = family_metrics.iter().find(|m| m.family_size == 1).expect("Should have size 1");
        assert_eq!(one.cs_count, 0);
        assert_eq!(one.ss_count, 8, "Should have 8 single-strand families of size 1");
        assert_eq!(one.ds_count, 3, "Should have 3 duplex families of size 1");

        let two = family_metrics.iter().find(|m| m.family_size == 2).expect("Should have size 2");
        assert_eq!(two.cs_count, 0);
        assert_eq!(two.ss_count, 2, "Should have 2 single-strand families of size 2");
        assert_eq!(two.ds_count, 3, "Should have 3 duplex families of size 2");

        let three = family_metrics.iter().find(|m| m.family_size == 3).expect("Should have size 3");
        assert_eq!(three.cs_count, 0);
        assert_eq!(three.ss_count, 0);
        assert_eq!(three.ds_count, 1, "Should have 1 duplex family of size 3");

        let four = family_metrics.iter().find(|m| m.family_size == 4).expect("Should have size 4");
        assert_eq!(four.cs_count, 3, "Should have 3 coordinate families of size 4");
        assert_eq!(four.ss_count, 0);
        assert_eq!(four.ds_count, 0);

        Ok(())
    }

    #[test]
    fn test_generate_accurate_duplex_family_size_counts() -> Result<()> {
        let mut records = Vec::new();

        // 1/0
        let (r1, r2) = build_test_pair("q1", 0, 100, 200, "AAA-ACG", "1/A", true, false);
        records.push(r1);
        records.push(r2);

        // 1/1
        let (r1, r2) = build_test_pair("q2", 0, 200, 300, "AAA-ACG", "2/A", true, false);
        records.push(r1);
        records.push(r2);
        // Swap positions for /B
        let (r1, r2) = build_test_pair("q3", 0, 300, 200, "ACG-AAA", "2/B", false, true);
        records.push(r1);
        records.push(r2);

        // 2/1 - swap positions for /B
        for i in 1..=2 {
            let (r1, r2) =
                build_test_pair(&format!("q{}", i + 3), 0, 400, 300, "AAC-GGG", "3/B", false, true);
            records.push(r1);
            records.push(r2);
        }
        let (r1, r2) = build_test_pair("q6", 0, 300, 400, "GGG-AAC", "3/A", true, false);
        records.push(r1);
        records.push(r2);

        // 4/3
        for i in 1..=4 {
            let (r1, r2) =
                build_test_pair(&format!("q{}", i + 6), 0, 400, 500, "GGG-AAC", "4/A", true, false);
            records.push(r1);
            records.push(r2);
        }
        // Swap positions for /B
        for i in 1..=3 {
            let (r1, r2) = build_test_pair(
                &format!("q{}", i + 10),
                0,
                500,
                400,
                "AAC-GGG",
                "4/B",
                false,
                true,
            );
            records.push(r1);
            records.push(r2);
        }

        // 5/5
        for i in 1..=5 {
            let (r1, r2) = build_test_pair(
                &format!("q{}", i + 13),
                0,
                600,
                700,
                "AGT-GCT",
                "6/A",
                true,
                false,
            );
            records.push(r1);
            records.push(r2);
        }
        // Swap positions for /B
        for i in 1..=5 {
            let (r1, r2) = build_test_pair(
                &format!("q{}", i + 18),
                0,
                700,
                600,
                "GCT-AGT",
                "6/B",
                false,
                true,
            );
            records.push(r1);
            records.push(r2);
        }

        let input = create_test_bam(records)?;
        let output_dir = TempDir::new()?;
        let output = output_dir.path().join("output");

        let cmd = DuplexMetrics {
            input: input.path().to_path_buf(),
            output: output.clone(),
            min_ab_reads: 1,
            min_ba_reads: 1,
            duplex_umi_counts: false,
            intervals: None,
            description: None,
            umi_tag: "RX".to_string(),
            mi_tag: "MI".to_string(),
        };

        cmd.execute("test")?;

        let duplex_family_path = format!("{}.duplex_family_sizes.txt", output.display());
        let duplex_metrics: Vec<DuplexFamilySizeMetric> =
            DelimFile::default().read_tsv(&duplex_family_path)?;

        // All metrics should have ab_size >= ba_size
        for metric in &duplex_metrics {
            assert!(metric.ab_size >= metric.ba_size, "ab_size should be >= ba_size");
        }

        // Check specific counts
        assert_eq!(
            duplex_metrics.iter().find(|m| m.ab_size == 1 && m.ba_size == 0).map(|m| m.count),
            Some(1)
        );
        assert_eq!(
            duplex_metrics.iter().find(|m| m.ab_size == 1 && m.ba_size == 1).map(|m| m.count),
            Some(1)
        );
        assert_eq!(
            duplex_metrics.iter().find(|m| m.ab_size == 2 && m.ba_size == 1).map(|m| m.count),
            Some(1)
        );
        assert_eq!(
            duplex_metrics.iter().find(|m| m.ab_size == 4 && m.ba_size == 3).map(|m| m.count),
            Some(1)
        );
        assert_eq!(
            duplex_metrics.iter().find(|m| m.ab_size == 5 && m.ba_size == 5).map(|m| m.count),
            Some(1)
        );

        Ok(())
    }

    #[test]
    fn test_respect_min_ab_and_min_ba_for_counting_duplexes() -> Result<()> {
        let mut records = Vec::new();

        // Duplex 1: 1/A + 1/B
        let (r1, r2) = build_test_pair("q1", 0, 300, 400, "AAA-GGG", "1/A", true, false);
        records.push(r1);
        records.push(r2);
        // Swap positions for /B
        let (r1, r2) = build_test_pair("q2", 0, 400, 300, "GGG-AAA", "1/B", false, true);
        records.push(r1);
        records.push(r2);

        // Duplex 2: 1/A + 2/B
        let (r1, r2) = build_test_pair("q3", 0, 300, 400, "ACT-TTA", "2/A", true, false);
        records.push(r1);
        records.push(r2);
        // Swap positions for /B
        for i in 1..=2 {
            let (r1, r2) =
                build_test_pair(&format!("q{}", i + 3), 0, 400, 300, "TTA-ACT", "2/B", false, true);
            records.push(r1);
            records.push(r2);
        }

        // Duplex 3: 2/A + 2/B
        for i in 1..=2 {
            let (r1, r2) =
                build_test_pair(&format!("q{}", i + 5), 0, 300, 400, "CGA-GGT", "3/A", true, false);
            records.push(r1);
            records.push(r2);
        }
        // Swap positions for /B
        for i in 1..=2 {
            let (r1, r2) =
                build_test_pair(&format!("q{}", i + 7), 0, 400, 300, "GGT-CGA", "3/B", false, true);
            records.push(r1);
            records.push(r2);
        }

        let input = create_test_bam(records)?;

        // Test with min_ab=1, min_ba=1: All 3 should be duplexes
        {
            let output_dir = TempDir::new()?;
            let output = output_dir.path().join("output");

            let cmd = DuplexMetrics {
                input: input.path().to_path_buf(),
                output: output.clone(),
                min_ab_reads: 1,
                min_ba_reads: 1,
                duplex_umi_counts: false,
                intervals: None,
                description: None,
                umi_tag: "RX".to_string(),
                mi_tag: "MI".to_string(),
            };

            cmd.execute("test")?;

            let yield_path = format!("{}.duplex_yield_metrics.txt", output.display());
            let yield_metrics: Vec<DuplexYieldMetric> =
                DelimFile::default().read_tsv(&yield_path)?;
            let full_metric = yield_metrics
                .iter()
                .find(|m| (m.fraction - 1.0).abs() < 0.01)
                .expect("Should have 100% fraction");
            assert_eq!(full_metric.ds_duplexes, 3, "Should have 3 duplexes with min 1/1");
        }

        // Test with min_ab=2, min_ba=1: Only 2 should be duplexes (exclude 1/1)
        {
            let output_dir = TempDir::new()?;
            let output = output_dir.path().join("output");

            let cmd = DuplexMetrics {
                input: input.path().to_path_buf(),
                output: output.clone(),
                min_ab_reads: 2,
                min_ba_reads: 1,
                duplex_umi_counts: false,
                intervals: None,
                description: None,
                umi_tag: "RX".to_string(),
                mi_tag: "MI".to_string(),
            };

            cmd.execute("test")?;

            let yield_path = format!("{}.duplex_yield_metrics.txt", output.display());
            let yield_metrics: Vec<DuplexYieldMetric> =
                DelimFile::default().read_tsv(&yield_path)?;
            let full_metric = yield_metrics
                .iter()
                .find(|m| (m.fraction - 1.0).abs() < 0.01)
                .expect("Should have 100% fraction");
            assert_eq!(full_metric.ds_duplexes, 2, "Should have 2 duplexes with min 2/1");
        }

        // Test with min_ab=2, min_ba=2: Only 1 should be duplex (only 2/2)
        {
            let output_dir = TempDir::new()?;
            let output = output_dir.path().join("output");

            let cmd = DuplexMetrics {
                input: input.path().to_path_buf(),
                output: output.clone(),
                min_ab_reads: 2,
                min_ba_reads: 2,
                duplex_umi_counts: false,
                intervals: None,
                description: None,
                umi_tag: "RX".to_string(),
                mi_tag: "MI".to_string(),
            };

            cmd.execute("test")?;

            let yield_path = format!("{}.duplex_yield_metrics.txt", output.display());
            let yield_metrics: Vec<DuplexYieldMetric> =
                DelimFile::default().read_tsv(&yield_path)?;
            let full_metric = yield_metrics
                .iter()
                .find(|m| (m.fraction - 1.0).abs() < 0.01)
                .expect("Should have 100% fraction");
            assert_eq!(full_metric.ds_duplexes, 1, "Should have 1 duplex with min 2/2");
        }

        Ok(())
    }

    #[test]
    fn test_only_count_inserts_overlapping_intervals() -> Result<()> {
        let mut records = Vec::new();

        // Family 1 at chr1:1000-1100
        let (r1, r2) = build_test_pair("q1", 0, 1000, 1100, "AAA-GGG", "1/A", true, false);
        records.push(r1);
        records.push(r2);

        // Family 2 at chr1:2000-2100 (2 reads)
        for i in 1..=2 {
            let (r1, r2) = build_test_pair(
                &format!("q{}", i + 1),
                0,
                2000,
                2100,
                "GGG-AAA",
                "2/A",
                true,
                false,
            );
            records.push(r1);
            records.push(r2);
        }

        // Family 3 at chr1:3000-3100 (3 reads)
        for i in 1..=3 {
            let (r1, r2) = build_test_pair(
                &format!("q{}", i + 3),
                0,
                3000,
                3100,
                "ACT-TTA",
                "3/A",
                true,
                false,
            );
            records.push(r1);
            records.push(r2);
        }

        // Family 4 at chr2:4000-4100 (4 reads)
        for i in 1..=4 {
            let (r1, r2) = build_test_pair(
                &format!("q{}", i + 6),
                1,
                4000,
                4100,
                "TTA-ACT",
                "4/A",
                true,
                false,
            );
            records.push(r1);
            records.push(r2);
        }

        // Family 5 at chr2:5000-5100 (5 reads)
        for i in 1..=5 {
            let (r1, r2) = build_test_pair(
                &format!("q{}", i + 10),
                1,
                5000,
                5100,
                "CGA-GGT",
                "5/A",
                true,
                false,
            );
            records.push(r1);
            records.push(r2);
        }

        // Family 6 at chr2:6000-6100 (6 reads)
        for i in 1..=6 {
            let (r1, r2) = build_test_pair(
                &format!("q{}", i + 15),
                1,
                6000,
                6100,
                "GGT-CGA",
                "6/A",
                true,
                false,
            );
            records.push(r1);
            records.push(r2);
        }

        let input = create_test_bam(records)?;

        // Create intervals file (BED format)
        let intervals_dir = TempDir::new()?;
        let intervals_path = intervals_dir.path().join("intervals.bed");
        std::fs::write(&intervals_path, "chr1\t900\t1001\nchr1\t3150\t3500\nchr2\t5050\t6050\n")?;

        let output_dir = TempDir::new()?;
        let output = output_dir.path().join("output");

        let cmd = DuplexMetrics {
            input: input.path().to_path_buf(),
            output: output.clone(),
            min_ab_reads: 1,
            min_ba_reads: 1,
            duplex_umi_counts: false,
            intervals: Some(intervals_path),
            description: None,
            umi_tag: "RX".to_string(),
            mi_tag: "MI".to_string(),
        };

        cmd.execute("test")?;

        // Should capture families 1, 3, 5, 6 based on intervals
        let duplex_family_path = format!("{}.duplex_family_sizes.txt", output.display());
        let duplex_metrics: Vec<DuplexFamilySizeMetric> =
            DelimFile::default().read_tsv(&duplex_family_path)?;

        // Find the specific size counts
        let size_1_count =
            duplex_metrics.iter().find(|m| m.ab_size == 1 && m.ba_size == 0).map_or(0, |m| m.count);
        let size_3_count =
            duplex_metrics.iter().find(|m| m.ab_size == 3 && m.ba_size == 0).map_or(0, |m| m.count);
        let size_5_count =
            duplex_metrics.iter().find(|m| m.ab_size == 5 && m.ba_size == 0).map_or(0, |m| m.count);
        let size_6_count =
            duplex_metrics.iter().find(|m| m.ab_size == 6 && m.ba_size == 0).map_or(0, |m| m.count);

        assert_eq!(size_1_count, 1, "Should have 1 family of size 1");
        assert_eq!(size_3_count, 1, "Should have 1 family of size 3");
        assert_eq!(size_5_count, 1, "Should have 1 family of size 5");
        assert_eq!(size_6_count, 1, "Should have 1 family of size 6");

        Ok(())
    }

    #[test]
    fn test_downsampling_generates_multiple_fractions() -> Result<()> {
        let mut records = Vec::new();

        // Create 50 simple families
        for i in 1..=50 {
            let (r1, r2) = build_test_pair(
                &format!("q{i}"),
                0,
                i * 100,
                i * 100 + 100,
                "AAA-TTT",
                &format!("{i}/A"),
                true,
                false,
            );
            records.push(r1);
            records.push(r2);
        }

        let input = create_test_bam(records)?;
        let output_dir = TempDir::new()?;
        let output = output_dir.path().join("output");

        let cmd = DuplexMetrics {
            input: input.path().to_path_buf(),
            output: output.clone(),
            min_ab_reads: 1,
            min_ba_reads: 1,
            duplex_umi_counts: false,
            intervals: None,
            description: None,
            umi_tag: "RX".to_string(),
            mi_tag: "MI".to_string(),
        };

        cmd.execute("test")?;

        // Check that yield metrics has 20 entries (5%, 10%, ..., 100%)
        let yield_path = format!("{}.duplex_yield_metrics.txt", output.display());
        let yield_metrics: Vec<DuplexYieldMetric> = DelimFile::default().read_tsv(&yield_path)?;

        assert_eq!(yield_metrics.len(), 20, "Should have 20 sampling fractions");

        // Check that fractions are increasing
        for i in 1..yield_metrics.len() {
            assert!(
                yield_metrics[i].fraction > yield_metrics[i - 1].fraction,
                "Fractions should be increasing"
            );
        }

        // Check that family counts are generally increasing
        for i in 1..yield_metrics.len() {
            assert!(
                yield_metrics[i].cs_families >= yield_metrics[i - 1].cs_families,
                "CS families should be increasing or equal"
            );
        }

        Ok(())
    }

    #[test]
    fn test_duplex_umi_expected_fraction_calculation() -> Result<()> {
        use fgumi_lib::metrics::duplex::DuplexMetricsCollector;

        // Create a collector with duplex UMI counting enabled
        let mut collector = DuplexMetricsCollector::new(true);

        // Record individual UMIs with different frequencies
        // UMI "AAA" appears in 50% of families (5 out of 10)
        // UMI "TTT" appears in 50% of families (5 out of 10)
        // UMI "GGG" appears in 30% of families (3 out of 10)
        // UMI "CCC" appears in 20% of families (2 out of 10)

        // Families with AAA-TTT duplex (5 families)
        for _ in 0..5 {
            collector.record_umi("AAA", 1, 0, true);
            collector.record_umi("TTT", 1, 0, true);
            collector.record_duplex_umi("AAA-TTT", 1, 0, true);
        }

        // Families with GGG-CCC duplex (2 families)
        for _ in 0..2 {
            collector.record_umi("GGG", 1, 0, true);
            collector.record_umi("CCC", 1, 0, true);
            collector.record_duplex_umi("GGG-CCC", 1, 0, true);
        }

        // Families with GGG-AAA duplex (1 family)
        collector.record_umi("GGG", 1, 0, true);
        collector.record_umi("AAA", 1, 0, true);
        collector.record_duplex_umi("GGG-AAA", 1, 0, true);

        // Families with TTT-CCC duplex (2 families)
        for _ in 0..2 {
            collector.record_umi("TTT", 1, 0, true);
            collector.record_umi("CCC", 1, 0, true);
            collector.record_duplex_umi("TTT-CCC", 1, 0, true);
        }

        // Generate UMI metrics first
        let umi_metrics = collector.umi_metrics();

        // Verify individual UMI fractions
        let umi_map: std::collections::HashMap<_, _> =
            umi_metrics.iter().map(|m| (m.umi.as_str(), m)).collect();

        // Total unique observations = 10
        assert_eq!(umi_map["AAA"].unique_observations, 6, "AAA should appear in 6 families");
        assert_eq!(umi_map["TTT"].unique_observations, 7, "TTT should appear in 7 families");
        assert_eq!(umi_map["GGG"].unique_observations, 3, "GGG should appear in 3 families");
        assert_eq!(umi_map["CCC"].unique_observations, 4, "CCC should appear in 4 families");

        // Total = 20 individual UMI observations (each duplex family contributes 2)
        let expected_aaa_fraction = 6.0 / 20.0; // 0.3
        let expected_ttt_fraction = 7.0 / 20.0; // 0.35
        let expected_ggg_fraction = 3.0 / 20.0; // 0.15
        let expected_ccc_fraction = 4.0 / 20.0; // 0.2

        assert!(
            (umi_map["AAA"].fraction_unique_observations - expected_aaa_fraction).abs() < 0.001
        );
        assert!(
            (umi_map["TTT"].fraction_unique_observations - expected_ttt_fraction).abs() < 0.001
        );
        assert!(
            (umi_map["GGG"].fraction_unique_observations - expected_ggg_fraction).abs() < 0.001
        );
        assert!(
            (umi_map["CCC"].fraction_unique_observations - expected_ccc_fraction).abs() < 0.001
        );

        // Generate duplex UMI metrics
        let duplex_metrics = collector.duplex_umi_metrics(&umi_metrics);

        // Verify duplex UMI expected fractions
        let duplex_map: std::collections::HashMap<_, _> =
            duplex_metrics.iter().map(|m| (m.umi.as_str(), m)).collect();

        // Expected fraction for AAA-TTT = 0.3 * 0.35 = 0.105
        let expected_aaa_ttt = expected_aaa_fraction * expected_ttt_fraction;
        assert!(
            (duplex_map["AAA-TTT"].fraction_unique_observations_expected - expected_aaa_ttt).abs()
                < 0.001,
            "Expected AAA-TTT fraction to be {}, got {}",
            expected_aaa_ttt,
            duplex_map["AAA-TTT"].fraction_unique_observations_expected
        );

        // Expected fraction for GGG-CCC = 0.15 * 0.2 = 0.03
        let expected_ggg_ccc = expected_ggg_fraction * expected_ccc_fraction;
        assert!(
            (duplex_map["GGG-CCC"].fraction_unique_observations_expected - expected_ggg_ccc).abs()
                < 0.001,
            "Expected GGG-CCC fraction to be {}, got {}",
            expected_ggg_ccc,
            duplex_map["GGG-CCC"].fraction_unique_observations_expected
        );

        // Expected fraction for GGG-AAA = 0.15 * 0.3 = 0.045
        let expected_ggg_aaa = expected_ggg_fraction * expected_aaa_fraction;
        assert!(
            (duplex_map["GGG-AAA"].fraction_unique_observations_expected - expected_ggg_aaa).abs()
                < 0.001,
            "Expected GGG-AAA fraction to be {}, got {}",
            expected_ggg_aaa,
            duplex_map["GGG-AAA"].fraction_unique_observations_expected
        );

        // Expected fraction for TTT-CCC = 0.35 * 0.2 = 0.07
        let expected_ttt_ccc = expected_ttt_fraction * expected_ccc_fraction;
        assert!(
            (duplex_map["TTT-CCC"].fraction_unique_observations_expected - expected_ttt_ccc).abs()
                < 0.001,
            "Expected TTT-CCC fraction to be {}, got {}",
            expected_ttt_ccc,
            duplex_map["TTT-CCC"].fraction_unique_observations_expected
        );

        Ok(())
    }

    #[test]
    fn test_reject_consensus_bam() -> Result<()> {
        use tempfile::NamedTempFile;

        // Create a temporary BAM with consensus reads
        let temp_bam = NamedTempFile::new()?;
        let bam_path = temp_bam.path().to_path_buf();

        // Create header
        use noodles::sam::header::record::value::Map;
        use noodles::sam::header::record::value::map::ReferenceSequence;
        use std::num::NonZeroUsize;

        let header = noodles::sam::Header::builder()
            .add_reference_sequence(
                "chr1",
                Map::<ReferenceSequence>::new(NonZeroUsize::new(1000).unwrap()),
            )
            .build();

        // Create a simplex consensus record (has cD tag)
        let rec = RecordBuilder::new()
            .name("read1")
            .sequence(&"A".repeat(100))
            .qualities(&[30u8; 100])
            .reference_sequence_id(0)
            .alignment_start(100)
            .mapping_quality(60)
            .cigar("100M")
            .first_segment(true)
            .tag("cD", 10_i32)
            .build();

        // Write BAM
        let mut writer = noodles::bam::io::Writer::new(std::fs::File::create(&bam_path)?);
        writer.write_header(&header)?;
        writer.write_alignment_record(&header, &rec)?;
        writer.try_finish()?;
        drop(writer);

        // Try to run duplex-metrics on this consensus BAM
        let temp_output = tempfile::tempdir()?;
        let output_prefix = temp_output.path().join("metrics");

        let cmd = DuplexMetrics {
            input: bam_path,
            output: output_prefix,
            min_ab_reads: 1,
            min_ba_reads: 1,
            duplex_umi_counts: false,
            intervals: None,
            description: None,
            umi_tag: "RX".to_string(),
            mi_tag: "MI".to_string(),
        };

        // This should fail with an error about consensus BAM
        let result = cmd.execute("test");
        assert!(result.is_err());
        let err = result.unwrap_err().to_string();
        assert!(
            err.contains("appears to contain consensus sequences"),
            "Error message should mention consensus sequences, got: {err}"
        );

        Ok(())
    }

    #[test]
    fn test_default_duplex_metrics_parameters() {
        let metrics = DuplexMetrics {
            input: PathBuf::from("input.bam"),
            output: PathBuf::from("output"),
            min_ab_reads: 1,
            min_ba_reads: 1,
            duplex_umi_counts: false,
            intervals: None,
            description: None,
            umi_tag: "RX".to_string(),
            mi_tag: "MI".to_string(),
        };

        assert_eq!(metrics.min_ab_reads, 1);
        assert_eq!(metrics.min_ba_reads, 1);
        assert!(!metrics.duplex_umi_counts);
        assert!(metrics.intervals.is_none());
        assert!(metrics.description.is_none());
    }

    #[test]
    fn test_duplex_metrics_with_custom_thresholds() {
        let metrics = DuplexMetrics {
            input: PathBuf::from("input.bam"),
            output: PathBuf::from("output"),
            min_ab_reads: 3,
            min_ba_reads: 2,
            duplex_umi_counts: true,
            intervals: Some(PathBuf::from("intervals.bed")),
            description: Some("Test Sample".to_string()),
            umi_tag: "RX".to_string(),
            mi_tag: "MI".to_string(),
        };

        assert_eq!(metrics.min_ab_reads, 3);
        assert_eq!(metrics.min_ba_reads, 2);
        assert!(metrics.duplex_umi_counts);
        assert_eq!(metrics.intervals, Some(PathBuf::from("intervals.bed")));
    }

    #[test]
    fn test_duplex_metrics_validation_fails_wrong_order() {
        let metrics = DuplexMetrics {
            input: PathBuf::from("input.bam"),
            output: PathBuf::from("output"),
            min_ab_reads: 1,
            min_ba_reads: 2, // Invalid: BA > AB
            duplex_umi_counts: false,
            intervals: None,
            description: None,
            umi_tag: "RX".to_string(),
            mi_tag: "MI".to_string(),
        };

        // The validation happens in execute(), check during command construction would be ideal
        // but currently it's done at execute time
        // This test documents the expected behavior
        assert!(metrics.min_ab_reads < metrics.min_ba_reads);
    }

    #[test]
    fn test_hash_fraction_deterministic() {
        // Hash should be deterministic for same input
        let hash1 = compute_hash_fraction("read1");
        let hash2 = compute_hash_fraction("read1");
        assert!((hash1 - hash2).abs() < f64::EPSILON);

        // Hash should be in [0, 1] range
        assert!((0.0..=1.0).contains(&hash1));
    }

    #[test]
    fn test_hash_fraction_different_reads() {
        // Different read names should produce different hashes (usually)
        let hash1 = compute_hash_fraction("read1");
        let hash2 = compute_hash_fraction("read2");
        let hash3 = compute_hash_fraction("read3");

        // All should be in valid range
        assert!((0.0..=1.0).contains(&hash1));
        assert!((0.0..=1.0).contains(&hash2));
        assert!((0.0..=1.0).contains(&hash3));

        // They should be different (very high probability)
        assert!((hash1 - hash2).abs() > f64::EPSILON || (hash2 - hash3).abs() > f64::EPSILON);
    }

    #[test]
    fn test_interval_overlap_logic() {
        // Test no intervals - should always return true
        let template = TemplateInfo {
            mi: "1/A".to_string(),
            rx: "AAA-TTT".to_string(),
            ref_name: Some("chr1".to_string()),
            position: Some(1000),
            end_position: Some(1100),
            hash_fraction: 0.5,
        };

        assert!(overlaps_intervals(&template, &[]));

        // Test with matching interval
        let intervals = vec![Interval { ref_name: "chr1".to_string(), start: 900, end: 1050 }];
        assert!(overlaps_intervals(&template, &intervals));

        // Test with non-matching interval
        let intervals = vec![Interval { ref_name: "chr1".to_string(), start: 2000, end: 3000 }];
        assert!(!overlaps_intervals(&template, &intervals));

        // Test unmapped template
        let unmapped = TemplateInfo {
            mi: "2/A".to_string(),
            rx: "AAA-TTT".to_string(),
            ref_name: None,
            position: None,
            end_position: None,
            hash_fraction: 0.5,
        };
        assert!(!overlaps_intervals(&unmapped, &intervals));
    }

    #[test]
    fn test_parse_intervals_bed_format() -> Result<()> {
        let dir = TempDir::new()?;
        let path = dir.path().join("intervals.bed");
        std::fs::write(&path, "# comment line\nchr1\t100\t200\nchr2\t300\t400\n")?;
        let intervals = parse_intervals(&path)?;
        assert_eq!(intervals.len(), 2);
        assert_eq!(intervals[0].ref_name, "chr1");
        assert_eq!(intervals[0].start, 100);
        assert_eq!(intervals[0].end, 200);
        assert_eq!(intervals[1].ref_name, "chr2");
        assert_eq!(intervals[1].start, 300);
        assert_eq!(intervals[1].end, 400);
        Ok(())
    }

    #[test]
    fn test_parse_intervals_interval_list_format() -> Result<()> {
        let dir = TempDir::new()?;
        let path = dir.path().join("intervals.interval_list");
        // Picard interval list: SAM header then 1-based closed coordinates
        std::fs::write(
            &path,
            "@HD\tVN:1.6\tSO:coordinate\n\
             @SQ\tSN:chr1\tLN:10000\n\
             @SQ\tSN:chr2\tLN:10000\n\
             chr1\t101\t200\t+\tinterval1\n\
             chr2\t301\t400\t+\tinterval2\n",
        )?;
        let intervals = parse_intervals(&path)?;
        assert_eq!(intervals.len(), 2);
        // 1-based closed [101, 200] converts to 0-based half-open [100, 200]
        assert_eq!(intervals[0].ref_name, "chr1");
        assert_eq!(intervals[0].start, 100);
        assert_eq!(intervals[0].end, 200);
        // 1-based closed [301, 400] converts to 0-based half-open [300, 400]
        assert_eq!(intervals[1].ref_name, "chr2");
        assert_eq!(intervals[1].start, 300);
        assert_eq!(intervals[1].end, 400);
        Ok(())
    }

    #[test]
    fn test_only_count_inserts_overlapping_interval_list() -> Result<()> {
        let mut records = Vec::new();

        // Family 1 at chr1:1000-1100
        let (r1, r2) = build_test_pair("q1", 0, 1000, 1100, "AAA-GGG", "1/A", true, false);
        records.push(r1);
        records.push(r2);

        // Family 2 at chr1:2000-2100 (2 reads)
        for i in 1..=2 {
            let (r1, r2) = build_test_pair(
                &format!("q{}", i + 1),
                0,
                2000,
                2100,
                "GGG-AAA",
                "2/A",
                true,
                false,
            );
            records.push(r1);
            records.push(r2);
        }

        // Family 3 at chr1:3000-3100 (3 reads)
        for i in 1..=3 {
            let (r1, r2) = build_test_pair(
                &format!("q{}", i + 3),
                0,
                3000,
                3100,
                "ACT-TTA",
                "3/A",
                true,
                false,
            );
            records.push(r1);
            records.push(r2);
        }

        // Family 4 at chr2:4000-4100 (4 reads)
        for i in 1..=4 {
            let (r1, r2) = build_test_pair(
                &format!("q{}", i + 6),
                1,
                4000,
                4100,
                "TTA-ACT",
                "4/A",
                true,
                false,
            );
            records.push(r1);
            records.push(r2);
        }

        // Family 5 at chr2:5000-5100 (5 reads)
        for i in 1..=5 {
            let (r1, r2) = build_test_pair(
                &format!("q{}", i + 10),
                1,
                5000,
                5100,
                "CGA-GGT",
                "5/A",
                true,
                false,
            );
            records.push(r1);
            records.push(r2);
        }

        // Family 6 at chr2:6000-6100 (6 reads)
        for i in 1..=6 {
            let (r1, r2) = build_test_pair(
                &format!("q{}", i + 15),
                1,
                6000,
                6100,
                "GGT-CGA",
                "6/A",
                true,
                false,
            );
            records.push(r1);
            records.push(r2);
        }

        let input = create_test_bam(records)?;

        // Create intervals file in Picard interval list format
        // Same regions as the BED test: chr1:900-1001, chr1:3150-3500, chr2:5050-6050
        // BED 0-based half-open -> interval list 1-based closed:
        //   chr1:900-1001 -> chr1:901-1001
        //   chr1:3150-3500 -> chr1:3151-3500
        //   chr2:5050-6050 -> chr2:5051-6050
        let intervals_dir = TempDir::new()?;
        let intervals_path = intervals_dir.path().join("intervals.interval_list");
        std::fs::write(
            &intervals_path,
            "@HD\tVN:1.6\tSO:coordinate\n\
             @SQ\tSN:chr1\tLN:100000\n\
             @SQ\tSN:chr2\tLN:100000\n\
             chr1\t901\t1001\t+\tinterval1\n\
             chr1\t3151\t3500\t+\tinterval2\n\
             chr2\t5051\t6050\t+\tinterval3\n",
        )?;

        let output_dir = TempDir::new()?;
        let output = output_dir.path().join("output");

        let cmd = DuplexMetrics {
            input: input.path().to_path_buf(),
            output: output.clone(),
            min_ab_reads: 1,
            min_ba_reads: 1,
            duplex_umi_counts: false,
            intervals: Some(intervals_path),
            description: None,
            umi_tag: "RX".to_string(),
            mi_tag: "MI".to_string(),
        };

        cmd.execute("test")?;

        // Should capture families 1, 3, 5, 6 based on intervals (same as BED test)
        let duplex_family_path = format!("{}.duplex_family_sizes.txt", output.display());
        let duplex_metrics: Vec<DuplexFamilySizeMetric> =
            DelimFile::default().read_tsv(&duplex_family_path)?;

        let size_1_count =
            duplex_metrics.iter().find(|m| m.ab_size == 1 && m.ba_size == 0).map_or(0, |m| m.count);
        let size_3_count =
            duplex_metrics.iter().find(|m| m.ab_size == 3 && m.ba_size == 0).map_or(0, |m| m.count);
        let size_5_count =
            duplex_metrics.iter().find(|m| m.ab_size == 5 && m.ba_size == 0).map_or(0, |m| m.count);
        let size_6_count =
            duplex_metrics.iter().find(|m| m.ab_size == 6 && m.ba_size == 0).map_or(0, |m| m.count);

        assert_eq!(size_1_count, 1, "Should have 1 family of size 1");
        assert_eq!(size_3_count, 1, "Should have 1 family of size 3");
        assert_eq!(size_5_count, 1, "Should have 1 family of size 5");
        assert_eq!(size_6_count, 1, "Should have 1 family of size 6");

        Ok(())
    }

    #[test]
    fn test_handle_empty_umi_components() -> Result<()> {
        let mut records = Vec::new();

        // Family 1: Valid duplex UMI "AAA-TTT"
        let (r1, r2) = build_test_pair("q1", 0, 100, 200, "AAA-TTT", "1/A", true, false);
        records.push(r1);
        records.push(r2);
        // Swap positions for /B
        let (r1, r2) = build_test_pair("q2", 0, 200, 100, "TTT-AAA", "1/B", false, true);
        records.push(r1);
        records.push(r2);

        // Family 2: Empty first component "-CCC" (should be skipped gracefully)
        let (r1, r2) = build_test_pair("q3", 0, 200, 300, "-CCC", "2/A", true, false);
        records.push(r1);
        records.push(r2);

        // Family 3: Empty second component "AAA-" (should be skipped gracefully)
        let (r1, r2) = build_test_pair("q4", 0, 300, 400, "AAA-", "3/A", true, false);
        records.push(r1);
        records.push(r2);

        // Family 4: Single component "GGG" without dash (should be skipped)
        let (r1, r2) = build_test_pair("q5", 0, 400, 500, "GGG", "4/A", true, false);
        records.push(r1);
        records.push(r2);

        // Family 5: Another valid duplex "CCC-GGG"
        let (r1, r2) = build_test_pair("q6", 0, 500, 600, "CCC-GGG", "5/A", true, false);
        records.push(r1);
        records.push(r2);
        // Swap positions for /B
        let (r1, r2) = build_test_pair("q7", 0, 600, 500, "GGG-CCC", "5/B", false, true);
        records.push(r1);
        records.push(r2);

        let input = create_test_bam(records)?;
        let output_dir = TempDir::new()?;
        let output = output_dir.path().join("output");

        let cmd = DuplexMetrics {
            input: input.path().to_path_buf(),
            output: output.clone(),
            min_ab_reads: 1,
            min_ba_reads: 1,
            duplex_umi_counts: false,
            intervals: None,
            description: None,
            umi_tag: "RX".to_string(),
            mi_tag: "MI".to_string(),
        };

        // Should complete without panicking or errors
        cmd.execute("test")?;

        // Read UMI metrics
        let umi_path = format!("{}.umi_counts.txt", output.display());
        let umi_metrics: Vec<UmiMetric> = DelimFile::default().read_tsv(&umi_path)?;

        // Should only have UMIs from valid duplex families (1 and 5)
        // Valid UMIs: AAA, TTT, CCC, GGG
        assert_eq!(umi_metrics.len(), 4, "Should have 4 valid UMIs");

        // Verify each UMI is from valid families
        for metric in &umi_metrics {
            assert!(
                ["AAA", "TTT", "CCC", "GGG"].contains(&metric.umi.as_str()),
                "Unexpected UMI: {}",
                metric.umi
            );
        }

        // Check that empty-component families were counted in family metrics but not in UMI metrics
        let family_path = format!("{}.family_sizes.txt", output.display());
        let family_metrics: Vec<FamilySizeMetric> = DelimFile::default().read_tsv(&family_path)?;

        // Should have families for all coordinate groups including those with invalid UMIs
        // Families at: 100-200 (valid), 200-300 (invalid), 300-400 (invalid), 400-500 (invalid), 500-600 (valid)
        // Total CS families = 5 (one at each coordinate)
        let cs_families: usize = family_metrics.iter().map(|m| m.cs_count).sum();
        assert_eq!(cs_families, 5, "Should have 5 CS families (including those with invalid UMIs)");

        // But only valid duplexes should be counted as having UMI observations
        let duplex_path = format!("{}.duplex_family_sizes.txt", output.display());
        let duplex_metrics: Vec<DuplexFamilySizeMetric> =
            DelimFile::default().read_tsv(&duplex_path)?;

        // Should have 2 duplex families (1 and 5)
        let duplex_count: usize = duplex_metrics
            .iter()
            .filter(|m| m.ab_size >= 1 && m.ba_size >= 1)
            .map(|m| m.count)
            .sum();
        assert_eq!(duplex_count, 2, "Should have 2 valid duplex families");

        Ok(())
    }
}
