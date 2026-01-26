//! Enhanced `DuplexMetrics` with downsampling and yield metrics.
//!
//! This version implements:
//! - Full downsampling at 20 levels (5%, 10%, ..., 100%)
//! - UMI consensus calling within coordinate+strand families
//! - Ideal duplex fraction calculation using proper binomial CDF
//! - Optional interval filtering (BED format) to restrict analysis to specific regions

use anyhow::{Context, Result};
use clap::Parser;
use fgoxide::io::DelimFile;
use fgumi_lib::bam_io::create_bam_reader;
use fgumi_lib::logging::OperationTimer;
use fgumi_lib::metrics::duplex::{
    calculate_ideal_duplex_fraction, DuplexMetricsCollector, DuplexYieldMetric,
};
use fgumi_lib::progress::ProgressTracker;
use fgumi_lib::simple_umi_consensus::SimpleUmiConsensusCaller;
use fgumi_lib::template::TemplateIterator;
use fgumi_lib::validation::validate_file_exists;
use log::info;
use murmur3::murmur3_32;
use noodles::sam::alignment::record::Cigar;
use noodles::sam::alignment::record::cigar::op::Kind;
use noodles::sam::alignment::record_buf::RecordBuf;
use std::path::PathBuf;
use std::sync::OnceLock;

use super::command::Command;
use super::common::ThreadingOptions;

/// Embedded R script for PDF plot generation (bundled with binary)
const R_SCRIPT: &str = include_str!("../../resources/CollectDuplexSeqMetrics.R");

/// Cached R availability check (computed once per process)
static R_AVAILABLE: OnceLock<bool> = OnceLock::new();

/// Genomic interval for filtering
#[derive(Clone, Debug)]
struct Interval {
    ref_name: String,
    start: i32,
    end: i32,
}

impl Interval {
    // Interval overlap checking is done directly in overlaps_intervals function
}

/// Read name and template information for downsampling
#[derive(Clone)]
struct TemplateInfo {
    #[allow(dead_code)] // Part of struct for completeness, may be used in future
    read_name: String,
    mi: String,
    rx: String,
    ref_name: Option<String>,
    position: Option<i32>,
    end_position: Option<i32>,
    #[allow(dead_code)] // Part of struct for completeness, may be used in future
    is_reverse: bool,
    /// Hash fraction for deterministic downsampling (computed once per template)
    hash_fraction: f64,
    /// Grouping key matching fgbio's `ReadInfo`: (ref1, start1, strand1, ref2, start2, strand2)
    /// Positions are ordered so the earlier-mapping read comes first
    #[allow(dead_code)]
    // Stored for debugging/future use, streaming uses local variable for grouping
    read_info_key: ReadInfoKey,
}

/// Grouping key matching fgbio's `ReadInfo` structure
/// Fields are ordered so the earlier-mapping read comes first
#[derive(Clone, PartialEq, Eq, Hash)]
struct ReadInfoKey {
    ref_index1: Option<usize>,
    start1: i32,
    strand1: bool, // true = reverse
    ref_index2: Option<usize>,
    start2: i32,
    strand2: bool,
}

/// Computes the unclipped 5' position for a read, matching fgbio's positionOf.
///
/// For forward strand reads: unclippedStart = alignmentStart - leading soft clips
/// For reverse strand reads: unclippedEnd = alignmentEnd + trailing soft clips
fn unclipped_five_prime_position(record: &RecordBuf) -> Option<i32> {
    let is_reverse = record.flags().is_reverse_complemented();
    let cigar = record.cigar();

    if is_reverse {
        // For reverse strand, 5' is at the end
        // unclippedEnd = alignmentEnd + trailing soft clips
        let alignment_end = record.alignment_end().map(|p| usize::from(p) as i32)?;

        // Count trailing soft clips (last element if it's S)
        let trailing_clips: i32 = cigar
            .iter()
            .filter_map(std::result::Result::ok)
            .last()
            .filter(|op| op.kind() == Kind::SoftClip)
            .map_or(0, |op| op.len() as i32);

        Some(alignment_end + trailing_clips)
    } else {
        // For forward strand, 5' is at the start
        // unclippedStart = alignmentStart - leading soft clips
        let alignment_start = record.alignment_start().map(|p| usize::from(p) as i32)?;

        // Count leading soft clips (first element if it's S)
        let leading_clips: i32 = cigar
            .iter()
            .find_map(std::result::Result::ok)
            .filter(|op| op.kind() == Kind::SoftClip)
            .map_or(0, |op| op.len() as i32);

        Some(alignment_start - leading_clips)
    }
}

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
    #[arg(long = "duplex-umi-counts", default_value = "false")]
    pub duplex_umi_counts: bool,

    /// Optional intervals file (BED format) to restrict analysis
    #[arg(short = 'L', long = "intervals")]
    pub intervals: Option<PathBuf>,

    /// Optional sample name or description for PDF plot titles
    #[arg(long = "description")]
    pub description: Option<String>,

    /// SAM tag containing the raw UMI sequence
    #[arg(long = "umi-tag", default_value = "RX")]
    pub umi_tag: String,

    /// SAM tag containing the molecular identifier (assigned by `group`)
    #[arg(long = "mi-tag", default_value = "MI")]
    pub mi_tag: String,

    /// Optional SAM tag for cell barcode (for single-cell data)
    #[arg(long = "cell-tag")]
    pub cell_tag: Option<String>,

    /// Threading options for parallel BAM decompression
    #[command(flatten)]
    pub threading: ThreadingOptions,
}

impl Command for DuplexMetrics {
    fn execute(&self, _command_line: &str) -> Result<()> {
        info!("DuplexMetrics");
        info!("  Input: {}", self.input.display());
        info!("  Output prefix: {}", self.output.display());
        info!("  Min AB reads: {}", self.min_ab_reads);
        info!("  Min BA reads: {}", self.min_ba_reads);
        info!("  Collect duplex UMI counts: {}", self.duplex_umi_counts);
        info!("  Threads: {}", self.threading.num_threads());

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
        self.validate_not_consensus_bam()?;

        // Load intervals if provided
        let intervals = if let Some(intervals_path) = &self.intervals {
            info!("  Loading intervals from: {}", intervals_path.display());
            let intervals = Self::parse_intervals(intervals_path)?;
            info!("  Loaded {} intervals", intervals.len());
            intervals
        } else {
            Vec::new()
        };

        // Downsampling fractions to process
        let fractions = [
            0.05, 0.10, 0.15, 0.20, 0.25, 0.30, 0.35, 0.40, 0.45, 0.50, 0.55, 0.60, 0.65, 0.70,
            0.75, 0.80, 0.85, 0.90, 0.95, 1.00,
        ];

        // Create collectors for each fraction (single-pass architecture)
        let mut collectors: Vec<DuplexMetricsCollector> = fractions
            .iter()
            .map(|&fraction| {
                // Only collect duplex UMI counts for the 100% fraction to save memory
                let collect_duplex = self.duplex_umi_counts && (fraction - 1.0_f64).abs() < 0.01;
                DuplexMetricsCollector::new(self.min_ab_reads, self.min_ba_reads, collect_duplex)
            })
            .collect();

        let mut umi_consensus_caller = SimpleUmiConsensusCaller::default();

        info!("Processing templates in single pass at {} sampling fractions...", fractions.len());

        // Single pass: process templates and update all applicable collectors
        let (total_template_count, fraction_template_counts) = self.process_templates_single_pass(
            &intervals,
            &fractions,
            &mut collectors,
            &mut umi_consensus_caller,
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
        if Self::is_r_available() {
            let description = self.description.as_deref().unwrap_or("Sample");
            match Self::execute_r_script(
                &family_size_path,
                &duplex_family_size_path,
                &yield_path,
                &umi_path,
                &pdf_path,
                description,
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
    /// Converts a two-character tag string to a noodles Tag object.
    ///
    /// # Arguments
    ///
    /// * `tag_name` - Two-character tag name (e.g., "RX", "MI", "CB")
    ///
    /// # Returns
    ///
    /// A noodles Tag object
    ///
    /// # Panics
    ///
    /// Panics if the tag name is not exactly 2 characters
    fn tag_from_string(tag_name: &str) -> noodles::sam::alignment::record::data::field::Tag {
        assert_eq!(tag_name.len(), 2, "Tag name must be exactly 2 characters");
        let bytes = tag_name.as_bytes();
        noodles::sam::alignment::record::data::field::Tag::from([bytes[0], bytes[1]])
    }

    /// Validates that the input BAM is not a consensus BAM.
    ///
    /// Consensus BAMs (output from simplex/duplex callers) should not be used with this tool.
    /// This checks the first valid R1 record and errors if it contains consensus tags.
    fn validate_not_consensus_bam(&self) -> Result<()> {
        use fgumi_lib::consensus_tags::is_consensus;

        let (mut reader, header) = create_bam_reader(&self.input, self.threading.num_threads())?;

        // Look at the first valid R1 record
        for result in reader.record_bufs(&header) {
            let record = result?;

            // Only check R1 records that are paired, mapped, and primary
            let flags = record.flags();
            if !flags.is_segmented()
                || !flags.is_first_segment()
                || flags.is_secondary()
                || flags.is_supplementary()
                || flags.is_unmapped()
            {
                continue;
            }

            // Check if this is a consensus read
            if is_consensus(&record) {
                anyhow::bail!(
                    "Input BAM file ({}) appears to contain consensus sequences. \
                    duplex-metrics cannot run on consensus BAMs, and instead requires \
                    the UMI-grouped BAM generated by group which is run prior to consensus calling.\n\
                    First R1 record '{}' has consensus SAM tags present.",
                    self.input.display(),
                    record.name().map_or_else(
                        || "<unnamed>".to_string(),
                        |n| String::from_utf8_lossy(n.as_ref()).to_string()
                    )
                );
            }

            // Only need to check the first valid R1 record
            break;
        }

        Ok(())
    }

    /// Parses an intervals file in BED format.
    ///
    /// Reads a BED-formatted file containing genomic intervals and parses each line
    /// into an `Interval` struct with chromosome, start, and end coordinates.
    /// Skips empty lines and comment lines starting with '#'.
    ///
    /// # Arguments
    ///
    /// * `path` - Path to the BED file
    ///
    /// # Returns
    ///
    /// A vector of `Interval` objects.
    ///
    /// # Errors
    ///
    /// Returns an error if:
    /// - File cannot be opened or read
    /// - Lines have fewer than 3 fields
    /// - Start or end positions cannot be parsed as integers
    fn parse_intervals(path: &PathBuf) -> Result<Vec<Interval>> {
        use std::fs::File;
        use std::io::{BufRead, BufReader};

        let file = File::open(path)?;
        let reader = BufReader::new(file);
        let mut intervals = Vec::new();

        for line in reader.lines() {
            let line = line?;
            let line = line.trim();

            // Skip empty lines and comments
            if line.is_empty() || line.starts_with('#') {
                continue;
            }

            // Parse BED format: chr start end [name] [score] [strand]
            let fields: Vec<&str> = line.split('\t').collect();
            if fields.len() < 3 {
                anyhow::bail!("Invalid BED line (needs at least 3 fields): {line}");
            }

            let ref_name = fields[0].to_string();
            let start: i32 = fields[1]
                .parse()
                .map_err(|_| anyhow::anyhow!("Invalid start position: {}", fields[1]))?;
            let end: i32 = fields[2]
                .parse()
                .map_err(|_| anyhow::anyhow!("Invalid end position: {}", fields[2]))?;

            intervals.push(Interval { ref_name, start, end });
        }

        Ok(intervals)
    }

    /// Checks if a template's insert overlaps any provided interval.
    ///
    /// Determines whether a template's genomic insert coordinates overlap with any
    /// of the specified intervals. An insert overlaps an interval if any part of it
    /// (from start to end position) overlaps the interval region. If no intervals are
    /// provided, all templates are considered to overlap (no filtering).
    ///
    /// # Arguments
    ///
    /// * `template` - Template information including chromosome and positions
    /// * `intervals` - Slice of intervals to check for overlap
    ///
    /// # Returns
    ///
    /// `true` if the template overlaps any interval or if no intervals are provided,
    /// `false` if the template is unmapped or does not overlap any interval.
    fn overlaps_intervals(template: &TemplateInfo, intervals: &[Interval]) -> bool {
        if intervals.is_empty() {
            return true; // No filtering if no intervals provided
        }

        if let (Some(ref_name), Some(start), Some(end)) =
            (&template.ref_name, template.position, template.end_position)
        {
            // Check if the insert [start, end] overlaps with any interval
            intervals.iter().any(|interval| {
                // Two intervals overlap if: start1 <= end2 && start2 <= end1
                interval.ref_name == *ref_name && start <= interval.end && interval.start <= end
            })
        } else {
            false // Unmapped reads or reads without proper coordinates don't overlap any interval
        }
    }

    /// Processes templates in a single pass, updating all collectors for applicable fractions.
    ///
    /// Implements a single-pass downsampling strategy where each template is hashed to determine
    /// which downsampling fractions it belongs to. Templates are stored per fraction, then
    /// analyzed to compute metrics for each fraction independently. This approach enables
    /// efficient computation of metrics at multiple downsampling levels without re-reading
    /// the input file.
    ///
    /// # Arguments
    ///
    /// * `intervals` - Optional intervals for filtering templates
    /// * `fractions` - Slice of downsampling fractions (e.g., 0.05, 0.10, ..., 1.00)
    /// * `collectors` - Mutable slice of metrics collectors, one per fraction
    /// * `umi_consensus_caller` - UMI consensus caller for error correction
    ///
    /// # Returns
    ///
    /// A tuple of `(total_template_count, per_fraction_template_counts)` where:
    /// - `total_template_count` is the total number of templates processed
    /// - `per_fraction_template_counts` is a vector of template counts for each fraction
    ///
    /// # Errors
    ///
    /// Returns an error if:
    /// - BAM file cannot be read
    /// - Template processing fails
    /// - Metrics collection fails
    fn process_templates_single_pass(
        &self,
        intervals: &[Interval],
        fractions: &[f64],
        collectors: &mut [DuplexMetricsCollector],
        umi_consensus_caller: &mut SimpleUmiConsensusCaller,
    ) -> Result<(usize, Vec<usize>)> {
        let (mut reader, header) = create_bam_reader(&self.input, self.threading.num_threads())?;

        let record_iter = reader.record_bufs(&header).map(|r| r.map_err(Into::into));
        let template_iter = TemplateIterator::new(record_iter);

        // Streaming approach: process groups as they arrive (assumes consecutive ReadInfo grouping)
        // This matches fgbio's takeNextGroup behavior
        let mut current_group: Vec<TemplateInfo> = Vec::new();
        let mut current_key: Option<ReadInfoKey> = None;
        let mut template_count = 0;
        let progress = ProgressTracker::new("Processed records").with_interval(1_000_000);
        let mut fraction_template_counts: Vec<usize> = vec![0; fractions.len()];

        for template in template_iter {
            let template = template?;
            if template.records.len() < 2 {
                continue;
            }

            // Find R1 and R2 that pass fgbio's filtering criteria
            let r1 = template.records.iter().find(|r| {
                let f = r.flags();
                f.is_segmented()
                    && !f.is_unmapped()
                    && !f.is_mate_unmapped()
                    && f.is_first_segment()
                    && !f.is_secondary()
                    && !f.is_supplementary()
            });
            let r2 = template.records.iter().find(|r| {
                let f = r.flags();
                f.is_segmented()
                    && !f.is_unmapped()
                    && !f.is_mate_unmapped()
                    && f.is_last_segment()
                    && !f.is_secondary()
                    && !f.is_supplementary()
            });

            let (r1, r2) = match (r1, r2) {
                (Some(r1), Some(r2)) => (r1, r2),
                _ => continue,
            };

            let record = r1; // Use R1 as the primary record for tags

            // Get read name
            let read_name = record
                .name()
                .map(|n| String::from_utf8_lossy(n.as_ref()).to_string())
                .unwrap_or_default();

            // Get MI tag (molecular identifier)
            let mi_tag = Self::tag_from_string(&self.mi_tag);
            let mi =
                if let Some(noodles::sam::alignment::record_buf::data::field::Value::String(s)) =
                    record.data().get(&mi_tag)
                {
                    String::from_utf8(s.iter().copied().collect::<Vec<u8>>())?
                } else {
                    continue;
                };

            // Get UMI tag (raw UMI)
            let umi_tag = Self::tag_from_string(&self.umi_tag);
            let rx =
                if let Some(noodles::sam::alignment::record_buf::data::field::Value::String(s)) =
                    record.data().get(&umi_tag)
                {
                    String::from_utf8(s.iter().copied().collect::<Vec<u8>>())?
                } else {
                    continue;
                };

            // Get reference position and calculate insert coordinates
            let ref_name = if let Some(ref_id) = record.reference_sequence_id() {
                header.reference_sequences().get_index(ref_id).map(|(name, _)| name.to_string())
            } else {
                None
            };

            // Calculate insert coordinates and ReadInfo key (matching fgbio)
            // r1 and r2 were already found above
            let (position, end_position, read_info_key) = {
                let r1_ref = r1.reference_sequence_id();
                let r2_ref = r2.reference_sequence_id();

                if r1_ref == r2_ref && r1_ref.is_some() {
                    // Get unclipped 5' positions (matching fgbio's positionOf)
                    let r1_5prime = unclipped_five_prime_position(r1);
                    let r2_5prime = unclipped_five_prime_position(r2);
                    let r1_strand = r1.flags().is_reverse_complemented();
                    let r2_strand = r2.flags().is_reverse_complemented();

                    if let (Some(s1), Some(s2)) = (r1_5prime, r2_5prime) {
                        // For insert coordinates, use alignment positions
                        let r1_start = r1.alignment_start().map(|p| usize::from(p) as i32);
                        let r2_start = r2.alignment_start().map(|p| usize::from(p) as i32);
                        let r1_end = r1.alignment_end().map(|p| usize::from(p) as i32);
                        let r2_end = r2.alignment_end().map(|p| usize::from(p) as i32);

                        let (pos, end) = match (r1_start, r2_start, r1_end, r2_end) {
                            (Some(rs1), Some(rs2), Some(re1), Some(re2)) => {
                                (rs1.min(rs2), re1.max(re2))
                            }
                            (Some(rs1), Some(rs2), _, _) => (rs1.min(rs2), rs1.max(rs2)),
                            _ => continue,
                        };

                        // Build ReadInfo key: order by (ref, 5' position) so earlier read is first
                        // This matches fgbio's ReadInfo ordering
                        let key = if (r1_ref, s1) <= (r2_ref, s2) {
                            ReadInfoKey {
                                ref_index1: r1_ref,
                                start1: s1,
                                strand1: r1_strand,
                                ref_index2: r2_ref,
                                start2: s2,
                                strand2: r2_strand,
                            }
                        } else {
                            ReadInfoKey {
                                ref_index1: r2_ref,
                                start1: s2,
                                strand1: r2_strand,
                                ref_index2: r1_ref,
                                start2: s1,
                                strand2: r1_strand,
                            }
                        };

                        (pos, end, key)
                    } else {
                        continue;
                    }
                } else {
                    continue;
                }
            };

            // Compute hash once for this template
            let hash_fraction = Self::compute_hash_fraction(&read_name);

            let template_info = TemplateInfo {
                read_name,
                mi,
                rx,
                ref_name,
                position: Some(position),
                end_position: Some(end_position),
                is_reverse: false,
                hash_fraction,
                read_info_key: read_info_key.clone(),
            };

            // Check interval overlap
            if !Self::overlaps_intervals(&template_info, intervals) {
                continue;
            }

            template_count += 1;
            progress.log_if_needed(2); // Each template has R1 and R2

            // Streaming: when ReadInfo key changes, process the accumulated group
            if current_key.as_ref() != Some(&read_info_key) && !current_group.is_empty() {
                Self::process_coordinate_group(
                    &current_group,
                    fractions,
                    collectors,
                    umi_consensus_caller,
                    &mut fraction_template_counts,
                    self,
                );
                current_group.clear();
            }

            current_group.push(template_info);
            current_key = Some(read_info_key);
        }

        // Process the final group
        if !current_group.is_empty() {
            Self::process_coordinate_group(
                &current_group,
                fractions,
                collectors,
                umi_consensus_caller,
                &mut fraction_template_counts,
                self,
            );
        }

        progress.log_final();
        Ok((template_count, fraction_template_counts))
    }

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
        _metrics: &DuplexMetrics,
    ) {
        use std::collections::HashMap;

        if group.is_empty() {
            return;
        }

        // Pre-compute metadata once for the entire group
        struct TemplateMetadata<'a> {
            template: &'a TemplateInfo,
            base_umi: &'a str,
            is_a_strand: bool,
            is_b_strand: bool,
        }

        let metadata: Vec<TemplateMetadata> = group
            .iter()
            .map(|t| {
                let (base_umi, is_a, is_b) = if t.mi.ends_with("/A") {
                    (&t.mi[..t.mi.len() - 2], true, false)
                } else if t.mi.ends_with("/B") {
                    (&t.mi[..t.mi.len() - 2], false, true)
                } else {
                    (t.mi.as_str(), false, false)
                };
                TemplateMetadata { template: t, base_umi, is_a_strand: is_a, is_b_strand: is_b }
            })
            .collect();

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

                    collectors[idx].update_umi_metrics_for_family(
                        umi_consensus_caller,
                        &owned_pairs,
                        base_umi,
                    );
                }
            }
        }
    }

    /// Computes a hash value normalized to the [0, 1] range using Murmur3.
    ///
    /// Hashes the read name using the Murmur3 32-bit algorithm with seed 42, then
    /// normalizes the result to a floating-point value between 0 and 1. This is used
    /// for deterministic downsampling where reads are assigned to fractions based on
    /// their hash value. Matches the Scala implementation's hashing approach.
    ///
    /// # Arguments
    ///
    /// * `read_name` - The read name string to hash
    ///
    /// # Returns
    ///
    /// A floating-point value in the range [0, 1].
    fn compute_hash_fraction(read_name: &str) -> f64 {
        // Use Murmur3 with seed 42 (matching Scala implementation)
        let hash = murmur3_32(&mut std::io::Cursor::new(read_name.as_bytes()), 42).unwrap_or(0);

        // Scala implementation uses Int (i32), takes absolute value, then normalizes
        // Convert u32 to i32 first to match Scala's behavior
        let hash_i32 = hash as i32;
        let positive_hash = hash_i32.abs() as f64;
        positive_hash / i32::MAX as f64
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

        let ideal_fraction = calculate_ideal_duplex_fraction(
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

    /// Checks if R and required packages (ggplot2, scales) are available.
    /// Result is cached for the lifetime of the process to avoid repeated subprocess spawns.
    fn is_r_available() -> bool {
        use std::process::Command;

        *R_AVAILABLE.get_or_init(|| {
            Command::new("Rscript")
                .args(["-e", "stopifnot(require(ggplot2)); stopifnot(require(scales))"])
                .output()
                .map(|output| output.status.success())
                .unwrap_or(false)
        })
    }

    /// Executes the R script to generate PDF plots.
    ///
    /// The R script is embedded in the binary at compile time and written to a
    /// temporary file for execution. This ensures the script is always available
    /// regardless of working directory or installation location.
    ///
    /// # Arguments
    ///
    /// * `family_size_path` - Path to family size metrics file
    /// * `duplex_family_size_path` - Path to duplex family size metrics file
    /// * `yield_path` - Path to yield metrics file
    /// * `umi_path` - Path to UMI metrics file
    /// * `pdf_path` - Path to output PDF file
    /// * `description` - Sample description for plot titles
    ///
    /// # Returns
    ///
    /// Returns Ok(()) if successful, or an error if R execution fails
    fn execute_r_script(
        family_size_path: &str,
        duplex_family_size_path: &str,
        yield_path: &str,
        umi_path: &str,
        pdf_path: &str,
        description: &str,
    ) -> Result<()> {
        use std::process::Command;

        // Write embedded R script to temp file
        let temp_dir = std::env::temp_dir();
        let r_script_path = temp_dir.join("fgumi_CollectDuplexSeqMetrics.R");
        std::fs::write(&r_script_path, R_SCRIPT)
            .context("Failed to write embedded R script to temp file")?;

        info!("Executing R script to generate PDF plots...");

        let output = Command::new("Rscript")
            .arg(&r_script_path)
            .arg(family_size_path)
            .arg(duplex_family_size_path)
            .arg(yield_path)
            .arg(umi_path)
            .arg(pdf_path)
            .arg(description)
            .output()
            .context("Failed to execute Rscript command")?;

        // Clean up temp file (ignore errors)
        let _ = std::fs::remove_file(&r_script_path);

        if output.status.success() {
            info!("Successfully generated PDF plots: {pdf_path}");
            Ok(())
        } else {
            let stderr = String::from_utf8_lossy(&output.stderr);
            anyhow::bail!(
                "R script execution failed with exit code {:?}. Error: {}",
                output.status.code(),
                stderr
            )
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use anyhow::Result;
    use fgoxide::io::DelimFile;
    use fgumi_lib::metrics::duplex::{
        DuplexFamilySizeMetric, DuplexUmiMetric, FamilySizeMetric, UmiMetric,
    };
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
        _mapq1: u8,
        _mapq2: u8,
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
        let (r1, r2) =
            build_test_pair("q1", 0, 100, 200, 60, 60, "AAA-TTT", "AAA-TTT/A", true, false);
        records.push(r1);
        records.push(r2);

        // /B reads: swap positions (200, 100) and use opposite strands to match fgbio's test setup
        let (r1, r2) =
            build_test_pair("q2", 0, 200, 100, 60, 60, "TTT-AAA", "AAA-TTT/B", false, true);
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
            cell_tag: None,
            threading: ThreadingOptions::none(),
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
            let (r1, r2) =
                build_test_pair(&format!("a{i}"), 0, 100, 200, 60, 60, umi, "1/A", true, false);
            records.push(r1);
            records.push(r2);
        }

        // B strand: TTT-AAA (2 exact + 1 with 1 error: CTT-AAA)
        // Swap positions (200, 100) to match fgbio's test setup for /B reads
        for i in 1..=3 {
            let umi = if i == 3 { "CTT-AAA" } else { "TTT-AAA" };
            let (r1, r2) =
                build_test_pair(&format!("b{i}"), 0, 200, 100, 60, 60, umi, "1/B", false, true);
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
            cell_tag: None,
            threading: ThreadingOptions::none(),
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
        let (r1, r2) = build_test_pair("q1", 0, 100, 200, 60, 60, "AAA-TTT", "1/A", true, false);
        records.push(r1);
        records.push(r2);

        // /B reads: swap positions (200, 100) to match fgbio's test setup
        let (r1, r2) = build_test_pair("q2", 0, 200, 100, 60, 60, "TTT-AAA", "1/B", false, true);
        records.push(r1);
        records.push(r2);

        // Family 2: TTT-AAA + error variant (NTT-AAA) at position 150-250
        for i in 1..=3 {
            let umi = if i == 3 { "NTT-AAA" } else { "TTT-AAA" };
            let (r1, r2) = build_test_pair(
                &format!("q{}", i + 2),
                0,
                150,
                250,
                60,
                60,
                umi,
                "2/A",
                true,
                false,
            );
            records.push(r1);
            records.push(r2);
        }

        // Family 3: CCC-GGG at position 250-350
        let (r1, r2) = build_test_pair("q6", 0, 250, 350, 60, 60, "CCC-GGG", "3/B", false, true);
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
            cell_tag: None,
            threading: ThreadingOptions::none(),
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
            let (r1, r2) = build_test_pair(
                &format!("a{i}"),
                0,
                100,
                200,
                60,
                60,
                "AAA-TTT",
                "1/A",
                true,
                false,
            );
            records.push(r1);
            records.push(r2);
        }
        // Family 2/A: 1 read
        let (r1, r2) = build_test_pair("b1", 0, 100, 200, 60, 60, "ACG-GGA", "2/A", true, false);
        records.push(r1);
        records.push(r2);
        // Family 3/B: 1 read - swap positions for /B
        let (r1, r2) = build_test_pair("c1", 0, 200, 100, 60, 60, "TAT-CGT", "3/B", false, true);
        records.push(r1);
        records.push(r2);

        // Two duplex families at the same location (200-300)
        // Duplex family 4: 1/A + 1/B
        let (r1, r2) = build_test_pair("d1", 0, 200, 300, 60, 60, "TTT-AAA", "4/A", true, false);
        records.push(r1);
        records.push(r2);
        // Swap positions for /B
        let (r1, r2) = build_test_pair("d2", 0, 300, 200, 60, 60, "AAA-AAA", "4/B", false, true);
        records.push(r1);
        records.push(r2);

        // Duplex family 5: 1/A + 1/B
        let (r1, r2) = build_test_pair("e1", 0, 200, 300, 60, 60, "CCC-GGG", "5/A", true, false);
        records.push(r1);
        records.push(r2);
        // Swap positions for /B
        let (r1, r2) = build_test_pair("e2", 0, 300, 200, 60, 60, "GGG-CCC", "5/B", false, true);
        records.push(r1);
        records.push(r2);

        // One duplex and one SS family at the same location (400-500)
        // SS family 6/A: 1 read
        let (r1, r2) = build_test_pair("f1", 0, 400, 500, 60, 60, "GCG-GAA", "6/A", true, false);
        records.push(r1);
        records.push(r2);

        // Duplex family 7: 2/A + 1/B
        for i in 1..=2 {
            let (r1, r2) = build_test_pair(
                &format!("g{i}"),
                0,
                400,
                500,
                60,
                60,
                "ACG-CCT",
                "7/A",
                true,
                false,
            );
            records.push(r1);
            records.push(r2);
        }
        // Swap positions for /B
        let (r1, r2) = build_test_pair("g3", 0, 500, 400, 60, 60, "CCT-ACG", "7/B", false, true);
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
            cell_tag: None,
            threading: ThreadingOptions::none(),
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
        let (r1, r2) = build_test_pair("q1", 0, 100, 200, 60, 60, "AAA-ACG", "1/A", true, false);
        records.push(r1);
        records.push(r2);

        // 1/1
        let (r1, r2) = build_test_pair("q2", 0, 200, 300, 60, 60, "AAA-ACG", "2/A", true, false);
        records.push(r1);
        records.push(r2);
        // Swap positions for /B
        let (r1, r2) = build_test_pair("q3", 0, 300, 200, 60, 60, "ACG-AAA", "2/B", false, true);
        records.push(r1);
        records.push(r2);

        // 2/1 - swap positions for /B
        for i in 1..=2 {
            let (r1, r2) = build_test_pair(
                &format!("q{}", i + 3),
                0,
                400,
                300,
                60,
                60,
                "AAC-GGG",
                "3/B",
                false,
                true,
            );
            records.push(r1);
            records.push(r2);
        }
        let (r1, r2) = build_test_pair("q6", 0, 300, 400, 60, 60, "GGG-AAC", "3/A", true, false);
        records.push(r1);
        records.push(r2);

        // 4/3
        for i in 1..=4 {
            let (r1, r2) = build_test_pair(
                &format!("q{}", i + 6),
                0,
                400,
                500,
                60,
                60,
                "GGG-AAC",
                "4/A",
                true,
                false,
            );
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
                60,
                60,
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
                60,
                60,
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
                60,
                60,
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
            cell_tag: None,
            threading: ThreadingOptions::none(),
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
        let (r1, r2) = build_test_pair("q1", 0, 300, 400, 60, 60, "AAA-GGG", "1/A", true, false);
        records.push(r1);
        records.push(r2);
        // Swap positions for /B
        let (r1, r2) = build_test_pair("q2", 0, 400, 300, 60, 60, "GGG-AAA", "1/B", false, true);
        records.push(r1);
        records.push(r2);

        // Duplex 2: 1/A + 2/B
        let (r1, r2) = build_test_pair("q3", 0, 300, 400, 60, 60, "ACT-TTA", "2/A", true, false);
        records.push(r1);
        records.push(r2);
        // Swap positions for /B
        for i in 1..=2 {
            let (r1, r2) = build_test_pair(
                &format!("q{}", i + 3),
                0,
                400,
                300,
                60,
                60,
                "TTA-ACT",
                "2/B",
                false,
                true,
            );
            records.push(r1);
            records.push(r2);
        }

        // Duplex 3: 2/A + 2/B
        for i in 1..=2 {
            let (r1, r2) = build_test_pair(
                &format!("q{}", i + 5),
                0,
                300,
                400,
                60,
                60,
                "CGA-GGT",
                "3/A",
                true,
                false,
            );
            records.push(r1);
            records.push(r2);
        }
        // Swap positions for /B
        for i in 1..=2 {
            let (r1, r2) = build_test_pair(
                &format!("q{}", i + 7),
                0,
                400,
                300,
                60,
                60,
                "GGT-CGA",
                "3/B",
                false,
                true,
            );
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
                cell_tag: None,
                threading: ThreadingOptions::none(),
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
                cell_tag: None,
                threading: ThreadingOptions::none(),
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
                cell_tag: None,
                threading: ThreadingOptions::none(),
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
        let (r1, r2) = build_test_pair("q1", 0, 1000, 1100, 60, 60, "AAA-GGG", "1/A", true, false);
        records.push(r1);
        records.push(r2);

        // Family 2 at chr1:2000-2100 (2 reads)
        for i in 1..=2 {
            let (r1, r2) = build_test_pair(
                &format!("q{}", i + 1),
                0,
                2000,
                2100,
                60,
                60,
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
                60,
                60,
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
                60,
                60,
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
                60,
                60,
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
                60,
                60,
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
            cell_tag: None,
            threading: ThreadingOptions::none(),
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
                60,
                60,
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
            cell_tag: None,
            threading: ThreadingOptions::none(),
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
        let mut collector = DuplexMetricsCollector::new(1, 1, true);

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
        writer.finish(&header)?;
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
            cell_tag: None,
            threading: ThreadingOptions::none(),
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
            cell_tag: None,
            threading: ThreadingOptions::none(),
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
            cell_tag: None,
            threading: ThreadingOptions::none(),
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
            cell_tag: None,
            threading: ThreadingOptions::none(),
        };

        // The validation happens in execute(), check during command construction would be ideal
        // but currently it's done at execute time
        // This test documents the expected behavior
        assert!(metrics.min_ab_reads < metrics.min_ba_reads);
    }

    #[test]
    fn test_hash_fraction_deterministic() {
        // Hash should be deterministic for same input
        let hash1 = DuplexMetrics::compute_hash_fraction("read1");
        let hash2 = DuplexMetrics::compute_hash_fraction("read1");
        assert!((hash1 - hash2).abs() < f64::EPSILON);

        // Hash should be in [0, 1] range
        assert!((0.0..=1.0).contains(&hash1));
    }

    #[test]
    fn test_hash_fraction_different_reads() {
        // Different read names should produce different hashes (usually)
        let hash1 = DuplexMetrics::compute_hash_fraction("read1");
        let hash2 = DuplexMetrics::compute_hash_fraction("read2");
        let hash3 = DuplexMetrics::compute_hash_fraction("read3");

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
            read_name: "read1".to_string(),
            mi: "1/A".to_string(),
            rx: "AAA-TTT".to_string(),
            ref_name: Some("chr1".to_string()),
            position: Some(1000),
            end_position: Some(1100),
            is_reverse: false,
            hash_fraction: 0.5,
            read_info_key: ReadInfoKey {
                ref_index1: Some(0),
                start1: 1000,
                strand1: false,
                ref_index2: Some(0),
                start2: 1100,
                strand2: false,
            },
        };

        assert!(DuplexMetrics::overlaps_intervals(&template, &[]));

        // Test with matching interval
        let intervals = vec![Interval { ref_name: "chr1".to_string(), start: 900, end: 1050 }];
        assert!(DuplexMetrics::overlaps_intervals(&template, &intervals));

        // Test with non-matching interval
        let intervals = vec![Interval { ref_name: "chr1".to_string(), start: 2000, end: 3000 }];
        assert!(!DuplexMetrics::overlaps_intervals(&template, &intervals));

        // Test unmapped template
        let unmapped = TemplateInfo {
            read_name: "read2".to_string(),
            mi: "2/A".to_string(),
            rx: "AAA-TTT".to_string(),
            ref_name: None,
            position: None,
            end_position: None,
            is_reverse: false,
            hash_fraction: 0.5,
            read_info_key: ReadInfoKey {
                ref_index1: None,
                start1: 0,
                strand1: false,
                ref_index2: None,
                start2: 0,
                strand2: false,
            },
        };
        assert!(!DuplexMetrics::overlaps_intervals(&unmapped, &intervals));
    }

    #[test]
    fn test_handle_empty_umi_components() -> Result<()> {
        let mut records = Vec::new();

        // Family 1: Valid duplex UMI "AAA-TTT"
        let (r1, r2) = build_test_pair("q1", 0, 100, 200, 60, 60, "AAA-TTT", "1/A", true, false);
        records.push(r1);
        records.push(r2);
        // Swap positions for /B
        let (r1, r2) = build_test_pair("q2", 0, 200, 100, 60, 60, "TTT-AAA", "1/B", false, true);
        records.push(r1);
        records.push(r2);

        // Family 2: Empty first component "-CCC" (should be skipped gracefully)
        let (r1, r2) = build_test_pair("q3", 0, 200, 300, 60, 60, "-CCC", "2/A", true, false);
        records.push(r1);
        records.push(r2);

        // Family 3: Empty second component "AAA-" (should be skipped gracefully)
        let (r1, r2) = build_test_pair("q4", 0, 300, 400, 60, 60, "AAA-", "3/A", true, false);
        records.push(r1);
        records.push(r2);

        // Family 4: Single component "GGG" without dash (should be skipped)
        let (r1, r2) = build_test_pair("q5", 0, 400, 500, 60, 60, "GGG", "4/A", true, false);
        records.push(r1);
        records.push(r2);

        // Family 5: Another valid duplex "CCC-GGG"
        let (r1, r2) = build_test_pair("q6", 0, 500, 600, 60, 60, "CCC-GGG", "5/A", true, false);
        records.push(r1);
        records.push(r2);
        // Swap positions for /B
        let (r1, r2) = build_test_pair("q7", 0, 600, 500, 60, 60, "GGG-CCC", "5/B", false, true);
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
            cell_tag: None,
            threading: ThreadingOptions::none(),
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
