//! `SimplexMetrics` command for collecting QC metrics from simplex sequencing data.
//!
//! This module implements single-pass metric collection with deterministic downsampling
//! at 20 levels (5%, 10%, ..., 100%). It produces:
//! - Family size distributions (CS and SS families)
//! - Yield curves at each downsampling fraction
//! - UMI observation frequencies
//! - Optional PDF plots via an embedded R script

use crate::commands::common::{
    CompressionOptions, QueueMemoryOptions, SchedulerOptions, ThreadingOptions,
};
use crate::logging::OperationTimer;
use crate::metrics::simplex::SimplexMetricsCollector;
use crate::simple_umi_consensus::SimpleUmiConsensusCaller;
use crate::validation::validate_input_exists;
use anyhow::Result;
use clap::Parser;
use log::info;
use std::path::PathBuf;

use super::command::Command;
use super::shared_metrics::{
    DOWNSAMPLING_FRACTIONS, execute_r_script, is_r_available, parse_intervals,
    process_templates_from_bam, record_simplex_coordinate_group,
};

/// Projected options the chain builder needs for the simplex-metrics stage.
///
/// A free-standing struct (mirrors `RetagOptions`/`FilterOptions`) carrying the
/// three fields `add_metrics` reads: the output prefix, the `--min-reads`
/// threshold, and the optional intervals path. Threading/compression come from
/// the [`crate::pipeline::chains::SingleStageContext`], so they are not carried
/// here.
#[derive(Debug, Clone)]
pub struct SimplexMetricsOptions {
    /// Output prefix for the metrics TSV files.
    pub output: PathBuf,
    /// Minimum reads per SS family to count as a consensus family in yield metrics.
    pub min_reads: usize,
    /// Optional intervals file (BED or Picard interval-list) restricting analysis.
    pub intervals: Option<PathBuf>,
}

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
                               NOTE: the PDF (and `--description`) are not produced on the parallel
                               `--threads` path, which writes only the metrics TSVs; omit `--threads` for the PDF.

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

    /// Optional intervals file to restrict analysis. Both BED and Picard interval-list formats are
    /// auto-detected (an fgumi superset; fgbio's duplex analog accepts only the Picard interval list).
    #[arg(short = 'l', long = "intervals")]
    pub intervals: Option<PathBuf>,

    /// Optional sample name or description for PDF plot titles. When omitted, fgumi uses the
    /// literal "Sample" (fgbio instead derives the sample/library name from the BAM `@RG` header,
    /// so plot titles differ unless this is set).
    #[arg(long = "description")]
    pub description: Option<String>,

    /// Threading options. When `--threads N` is set, metrics collection runs on
    /// the typed-step pipeline (parallel BGZF decode + MI-grouping + parallel
    /// per-thread metric accumulation). Absent `--threads`, the original
    /// single-pass streaming collector runs, producing byte-identical output.
    #[command(flatten)]
    pub threading: ThreadingOptions,

    /// Scheduler and pipeline stats options (chain path only).
    #[command(flatten)]
    pub scheduler_opts: SchedulerOptions,

    /// Pipeline queue memory options (chain path only).
    #[command(flatten)]
    pub queue_memory: QueueMemoryOptions,
}

impl SimplexMetrics {
    /// Project the clap fields into [`SimplexMetricsOptions`] for the chain builder.
    #[must_use]
    pub fn to_simplex_metrics_options(&self) -> SimplexMetricsOptions {
        SimplexMetricsOptions {
            output: self.output.clone(),
            min_reads: self.min_reads,
            intervals: self.intervals.clone(),
        }
    }

    /// Run simplex-metrics on the declarative chain builder: a single
    /// `Stage::Metrics` chain with the simplex slot filled and a `SinkSpec::None`
    /// (no BAM output). Selected when `--threads` is set; the metrics files are
    /// written byte-identically to the serial path by the chain's finalize hook.
    fn execute_chain(&self, command_line: &str) -> Result<()> {
        use crate::commands::common::BamIoOptions;
        use crate::pipeline::chains::{ChainSpec, SingleStageContext, StageOptionsBag, build_for};

        // The metrics chain reads only `io.input` (+ CRC/async-reader policy);
        // `io.output` is ignored (metrics carry their own `--output` prefix in
        // the stage options and the sink is `None`), so a dummy output keeps the
        // shared `BamIoOptions`/`SingleStageContext` shape without a BAM write.
        let io = BamIoOptions {
            input: self.input.clone(),
            output: PathBuf::from("/dev/null"),
            ..Default::default()
        };
        let stage_opts = StageOptionsBag {
            simplex_metrics: Some(self.to_simplex_metrics_options()),
            ..Default::default()
        };
        // Metrics write no BAM, so compression is never used; a default keeps the
        // uniform `SingleStageContext` shape without exposing an inert
        // `--compression-level` flag on this command.
        let compression = CompressionOptions::default();
        let ctx = SingleStageContext {
            io: &io,
            threading: &self.threading,
            compression: &compression,
            scheduler: &self.scheduler_opts,
            queue_memory: &self.queue_memory,
            command_line,
        };
        let spec = ChainSpec::single_stage_metrics(stage_opts, &ctx);
        build_for(spec)?.run()
    }
}

impl Command for SimplexMetrics {
    fn execute(&self, command_line: &str) -> Result<()> {
        info!("SimplexMetrics");
        info!("  Input: {}", self.input.display());
        info!("  Output prefix: {}", self.output.display());
        info!("  Min reads: {}", self.min_reads);

        // Validate inputs
        validate_input_exists(&self.input, "input BAM")?;

        // fgbio's analog validates minReads >= 1; --min-reads 0 would label every
        // family (size >= 0) a consensus family, so reject it (SIMM3-02).
        if self.min_reads == 0 {
            anyhow::bail!("--min-reads must be >= 1 (got {})", self.min_reads);
        }

        // With `--threads`, run the parallel typed-step chain (parallel decode +
        // MI-grouping + per-thread metric accumulation). It writes the SAME TSV
        // files (family_sizes/umi_counts/simplex_yield_metrics) byte-identically
        // via the chain's finalize hook. The R/PDF step is chain-independent and
        // not reproduced there yet, so the chain path emits the three metrics
        // TSVs only (no PDF); the default (serial) path below keeps the PDF.
        if self.threading.threads.is_some() {
            // The parallel path writes the metrics TSVs only — it does not run the
            // R/PDF plot step. Surface that (and the resulting `--description`
            // no-op) so the flags read as inert rather than silently dropped,
            // mirroring the `warn_unwired_pipeline_flags` diagnostics.
            log::warn!(
                "--threads runs metrics on the parallel pipeline, which writes the metrics \
                 TSVs only and does not generate the *.simplex_qc.pdf plot; run without \
                 --threads to produce the PDF"
            );
            if self.description.is_some() {
                log::warn!(
                    "--description customizes the PDF plot title only, so it is ignored on \
                     the --threads path (no PDF is generated)"
                );
            }
            return self.execute_chain(command_line);
        }

        let timer = OperationTimer::new("Computing simplex metrics");

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
                record_simplex_coordinate_group(
                    group,
                    fractions,
                    &mut collectors,
                    &mut umi_consensus_caller,
                    fraction_counts,
                )?;
                // simplex UMIs are N-component by design, so there is no
                // wrong-segment-count guard (cf. DXM-03 for duplex).
                Ok(())
            },
        )?;

        info!("Processed {total_template_count} templates");

        // Generate yield metrics from each collector
        let mut yield_metrics = Vec::new();
        for ((&fraction, collector), &read_pairs) in
            fractions.iter().zip(collectors.iter()).zip(fraction_template_counts.iter())
        {
            let yield_metric = collector.to_yield_metric(fraction, read_pairs, self.min_reads);
            yield_metrics.push(yield_metric);
        }

        // Use the 100% fraction collector for main metrics
        let main_collector =
            collectors.pop().expect("collectors is non-empty (always includes 100% fraction)");

        // Generate and write metrics
        info!("Writing metrics...");

        // Family size metrics. Use `write_metrics_auto` (not a bare `DelimFile::write_tsv`)
        // so an empty distribution still produces a header-only, fgbio-readable file.
        let family_size_metrics = main_collector.family_size_metrics();
        let family_size_path = format!("{}.family_sizes.txt", self.output.display());
        crate::metrics::writer::write_metrics_auto(&family_size_path, &family_size_metrics)?;
        info!("Wrote family size metrics to {family_size_path}");

        // UMI metrics
        let umi_metrics = main_collector.umi_metrics();
        let umi_path = format!("{}.umi_counts.txt", self.output.display());
        crate::metrics::writer::write_metrics_auto(&umi_path, &umi_metrics)?;
        info!("Wrote UMI metrics to {umi_path}");

        // Yield metrics
        let yield_path = format!("{}.simplex_yield_metrics.txt", self.output.display());
        crate::metrics::writer::write_metrics_auto(&yield_path, &yield_metrics)?;
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

#[cfg(test)]
mod tests {
    use super::*;
    use crate::commands::shared_metrics::{
        consensus_guard_record, ensure_not_consensus_record, is_consensus_guard_record,
    };
    use crate::metrics::shared::UmiMetric;
    use crate::metrics::simplex::{
        SimplexFamilySizeMetric, SimplexMetricsCollector, SimplexYieldMetric,
    };
    use crate::sam::SamTag;
    use anyhow::Result;
    use fgoxide::io::DelimFile;
    use fgumi_raw_bam::{
        SamBuilder as RawSamBuilder, flags, raw_record_to_record_buf, testutil::encode_op,
    };
    use noodles::bam;
    use noodles::sam;
    use noodles::sam::alignment::io::Write;
    use rstest::rstest;
    use std::num::NonZeroUsize;
    use std::path::Path;
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
        b1.add_string_tag(SamTag::RX, rx_umi.as_bytes());
        b1.add_string_tag(SamTag::MI, mi_tag.as_bytes());
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
        b2.add_string_tag(SamTag::RX, rx_umi.as_bytes());
        b2.add_string_tag(SamTag::MI, mi_tag.as_bytes());
        let r2 = to_record_buf(b2.build());

        (r1, r2)
    }

    /// Like [`build_test_pair`] but with an explicit R1 strand, for constructing realistic
    /// duplex geometry (the BA strand has R1 reverse and its mate positions swapped).
    fn build_test_pair_stranded(
        name: &str,
        ref_id: usize,
        pos1: i32,
        pos2: i32,
        rx_umi: &str,
        mi_tag: &str,
        r1_reverse: bool,
    ) -> (sam::alignment::RecordBuf, sam::alignment::RecordBuf) {
        let seq = vec![b'A'; 100];
        let quals = vec![30u8; 100];
        let cigar = encode_op(0, 100); // 100M
        let (r1_rev, r1_mate_rev) =
            if r1_reverse { (flags::REVERSE, 0) } else { (0, flags::MATE_REVERSE) };
        let (r2_rev, r2_mate_rev) =
            if r1_reverse { (0, flags::MATE_REVERSE) } else { (flags::REVERSE, 0) };

        let mut b1 = RawSamBuilder::new();
        b1.read_name(name.as_bytes())
            .flags(flags::PAIRED | flags::FIRST_SEGMENT | r1_rev | r1_mate_rev)
            .ref_id(ref_id as i32)
            .pos(pos1 - 1)
            .mapq(60)
            .cigar_ops(&[cigar])
            .sequence(&seq)
            .qualities(&quals)
            .mate_ref_id(ref_id as i32)
            .mate_pos(pos2 - 1);
        b1.add_string_tag(SamTag::RX, rx_umi.as_bytes());
        b1.add_string_tag(SamTag::MI, mi_tag.as_bytes());
        let r1 = to_record_buf(b1.build());

        let mut b2 = RawSamBuilder::new();
        b2.read_name(name.as_bytes())
            .flags(flags::PAIRED | flags::LAST_SEGMENT | r2_rev | r2_mate_rev)
            .ref_id(ref_id as i32)
            .pos(pos2 - 1)
            .mapq(60)
            .cigar_ops(&[cigar])
            .sequence(&seq)
            .qualities(&quals)
            .mate_ref_id(ref_id as i32)
            .mate_pos(pos1 - 1);
        b2.add_string_tag(SamTag::RX, rx_umi.as_bytes());
        b2.add_string_tag(SamTag::MI, mi_tag.as_bytes());
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

    /// SIMM3-02: `--min-reads 0` is rejected (fgbio validates `minReads >= 1`);
    /// otherwise every family of size >= 0 is silently labeled a consensus family.
    #[test]
    fn test_simplex_metrics_rejects_min_reads_zero() -> Result<()> {
        let (r1, r2) = build_test_pair("q1", 0, 100, 200, "AAA", "1/A");
        let input = create_test_bam(vec![r1, r2])?;
        let output_dir = TempDir::new()?;
        let cmd = SimplexMetrics {
            input: input.path().to_path_buf(),
            output: output_dir.path().join("output"),
            min_reads: 0,
            intervals: None,
            description: None,
            threading: crate::commands::common::ThreadingOptions { threads: None },
            scheduler_opts: crate::commands::common::SchedulerOptions::default(),
            queue_memory: crate::commands::common::QueueMemoryOptions::default(),
        };
        let err = cmd.execute("test").expect_err("must reject --min-reads 0");
        assert!(err.to_string().contains("min-reads must be >= 1"), "unexpected: {err}");
        Ok(())
    }

    /// The typed-step chain path (`--threads N`) must produce byte-identical
    /// metrics TSVs to the original single-pass serial path (`--threads` absent).
    ///
    /// Builds a BAM spanning many coordinate groups (distinct positions) with a
    /// few MI families each, so a 4-worker run shards groups across workers and
    /// exercises the `BoundaryReorder` cross-batch closing. Runs the command
    /// serially and with `--threads 4` into separate output prefixes, then diffs
    /// every emitted TSV byte-for-byte — the end-to-end guard that the parallel
    /// per-thread accumulation + merge reassembles the metrics losslessly
    /// regardless of how groups are partitioned across workers.
    #[test]
    fn chain_threads_output_matches_serial() -> Result<()> {
        // Many small coordinate groups (distinct positions) with 1-3 templates
        // across 1-2 MI families each, so a 4-worker run distributes groups. RX
        // is drawn from a small pool so a coordinate group carries a MIX of
        // distinct and tied UMIs — otherwise a constant RX makes umi_counts.txt
        // parity blind to any within-family read reordering the chain path might
        // introduce (which `pair_records_by_read_name` is written to preserve).
        const RX_POOL: [&str; 4] = ["ACGT", "TGCA", "GGCC", "AATT"];
        let mut records = Vec::new();
        for g in 0..300i32 {
            let pos1 = 100 + g * 10;
            let pos2 = pos1 + 100;
            for k in 0..=(g % 3) {
                let mi = format!("{}/A", g * 2 + (k % 2));
                #[allow(clippy::cast_sign_loss)]
                let rx = RX_POOL[((g + k) as usize) % RX_POOL.len()];
                let (r1, r2) = build_test_pair(&format!("g{g}_k{k}"), 0, pos1, pos2, rx, &mi);
                records.push(r1);
                records.push(r2);
            }
        }
        let input = create_test_bam(records)?;
        let dir = TempDir::new()?;

        let run = |threads: Option<usize>, prefix: &str| -> Result<std::path::PathBuf> {
            let out = dir.path().join(prefix);
            let cmd = SimplexMetrics {
                input: input.path().to_path_buf(),
                output: out.clone(),
                min_reads: 1,
                intervals: None,
                description: None,
                threading: crate::commands::common::ThreadingOptions { threads },
                scheduler_opts: crate::commands::common::SchedulerOptions::default(),
                queue_memory: crate::commands::common::QueueMemoryOptions::default(),
            };
            cmd.execute("test")?;
            Ok(out)
        };

        let serial = run(None, "serial")?;
        let chain = run(Some(4), "chain")?;

        for suffix in ["family_sizes.txt", "umi_counts.txt", "simplex_yield_metrics.txt"] {
            let s = std::fs::read(format!("{}.{suffix}", serial.display()))?;
            let c = std::fs::read(format!("{}.{suffix}", chain.display()))?;
            assert_eq!(s, c, "{suffix} differs between the serial and --threads 4 chain paths");
        }
        Ok(())
    }

    /// `--intervals` must be honored identically on the chain (`--threads`) path
    /// and the serial path. Builds groups spanning a wide position range, writes a
    /// BED covering only part of it (so some groups are excluded), and diffs the
    /// TSVs between serial+intervals and `--threads 4`+intervals — the guard that
    /// the intervals wiring on the parallel path filters the same groups, so a
    /// broken chain-path wiring cannot silently produce unfiltered counts.
    #[test]
    fn chain_threads_intervals_match_serial() -> Result<()> {
        let mut records = Vec::new();
        for g in 0..120i32 {
            let pos1 = 100 + g * 20;
            let pos2 = pos1 + 100;
            let mi = format!("{g}/A");
            let (r1, r2) = build_test_pair(&format!("g{g}"), 0, pos1, pos2, "ACGT", &mi);
            records.push(r1);
            records.push(r2);
        }
        let input = create_test_bam(records)?;
        let dir = TempDir::new()?;
        // BED (0-based half-open) covering only chr1:0-1200, so groups whose
        // coordinate is past ~1200 are excluded on both paths — the intervals
        // must actually filter, not pass everything through.
        let bed = dir.path().join("regions.bed");
        std::fs::write(&bed, "chr1\t0\t1200\n")?;

        let run = |threads: Option<usize>, prefix: &str| -> Result<std::path::PathBuf> {
            let out = dir.path().join(prefix);
            let cmd = SimplexMetrics {
                input: input.path().to_path_buf(),
                output: out.clone(),
                min_reads: 1,
                intervals: Some(bed.clone()),
                description: None,
                threading: crate::commands::common::ThreadingOptions { threads },
                scheduler_opts: crate::commands::common::SchedulerOptions::default(),
                queue_memory: crate::commands::common::QueueMemoryOptions::default(),
            };
            cmd.execute("test")?;
            Ok(out)
        };

        let serial = run(None, "serial")?;
        let chain = run(Some(4), "chain")?;

        for suffix in ["family_sizes.txt", "umi_counts.txt", "simplex_yield_metrics.txt"] {
            let s = std::fs::read(format!("{}.{suffix}", serial.display()))?;
            let c = std::fs::read(format!("{}.{suffix}", chain.display()))?;
            assert_eq!(
                s, c,
                "{suffix} differs with --intervals between the serial and --threads 4 paths"
            );
        }
        Ok(())
    }

    /// SIMM3-01: simplex-metrics must reject duplex-UMI input (a base UMI observed on
    /// both the /A and /B strands at one coordinate), rather than silently producing
    /// garbage UMI counts from the mixed-orientation RX consensus.
    #[test]
    fn test_simplex_metrics_rejects_duplex_input() {
        // Real duplex geometry for base UMI "1": the AB strand (R1 forward @100, R2
        // reverse @200) and the BA strand (R1 reverse @200, R2 forward @100, UMI halves
        // swapped) canonicalize to the same coordinate/strand key and so co-group.
        let mut records = Vec::new();
        let (r1, r2) = build_test_pair_stranded("q1", 0, 100, 200, "AAA-TTT", "1/A", false);
        records.push(r1);
        records.push(r2);
        let (r1, r2) = build_test_pair_stranded("q2", 0, 200, 100, "TTT-AAA", "1/B", true);
        records.push(r1);
        records.push(r2);

        let input = create_test_bam(records).expect("write test bam");
        let output_dir = TempDir::new().expect("tempdir");
        let output = output_dir.path().join("output");
        let cmd = SimplexMetrics {
            input: input.path().to_path_buf(),
            output,
            min_reads: 1,
            intervals: None,
            description: None,
            threading: crate::commands::common::ThreadingOptions { threads: None },
            scheduler_opts: crate::commands::common::SchedulerOptions::default(),
            queue_memory: crate::commands::common::QueueMemoryOptions::default(),
        };
        let message =
            cmd.execute("test").expect_err("duplex-UMI input must be rejected").to_string();
        // Pin the full diagnostic contract, not just the word "duplex": the message must both
        // name the offending input ("duplex-UMI data") and give the actionable next step (run
        // "duplex-metrics"), so a regression that drops either half is caught.
        assert!(
            message.contains("duplex-UMI data"),
            "error should name duplex-UMI data: {message}"
        );
        assert!(
            message.contains("duplex-metrics"),
            "error should point at duplex-metrics: {message}"
        );
    }

    /// SIMM3-01, chain path: the `--threads` typed-step path must reject
    /// duplex-UMI input just like the serial path, and the duplex-UMI diagnostic
    /// must survive the pipeline's `io::Error` wrapping. The rejection lives in
    /// the shared `record_coordinate_group`, but the chain records inside a
    /// worker and surfaces failures through `io::Error::other`, so this pins that
    /// the error still propagates to the caller with the actionable wording.
    #[test]
    fn chain_path_rejects_duplex_input() {
        let mut records = Vec::new();
        let (r1, r2) = build_test_pair_stranded("q1", 0, 100, 200, "AAA-TTT", "1/A", false);
        records.push(r1);
        records.push(r2);
        let (r1, r2) = build_test_pair_stranded("q2", 0, 200, 100, "TTT-AAA", "1/B", true);
        records.push(r1);
        records.push(r2);

        let input = create_test_bam(records).expect("write test bam");
        let output_dir = TempDir::new().expect("tempdir");
        let cmd = SimplexMetrics {
            input: input.path().to_path_buf(),
            output: output_dir.path().join("output"),
            min_reads: 1,
            intervals: None,
            description: None,
            threading: crate::commands::common::ThreadingOptions { threads: Some(4) },
            scheduler_opts: crate::commands::common::SchedulerOptions::default(),
            queue_memory: crate::commands::common::QueueMemoryOptions::default(),
        };
        let message =
            cmd.execute("test").expect_err("duplex-UMI input must be rejected").to_string();
        assert!(
            message.contains("duplex-UMI data"),
            "chain error should name duplex-UMI data: {message}"
        );
    }

    #[test]
    fn test_generate_yield_metric_basic() {
        let mut collector = SimplexMetricsCollector::new();
        collector.record_ss_family(1);
        collector.record_ss_family(2);
        collector.record_ss_family(3);
        collector.record_cs_family(6);

        let metric = collector.to_yield_metric(1.0, 6, 2);

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
        let metric = collector.to_yield_metric(0.5, 0, 1);

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

        let metric = collector.to_yield_metric(1.0, 3, 2);

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
            threading: crate::commands::common::ThreadingOptions { threads: None },
            scheduler_opts: crate::commands::common::SchedulerOptions::default(),
            queue_memory: crate::commands::common::QueueMemoryOptions::default(),
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

    /// SIM-01: `simplex-metrics` must count UMIs even when one component of a
    /// multi-part UMI is an empty molecule-end half (e.g. `CCC-`, `-GGG`), recording
    /// the empty half as an empty-string (`""`) UMI. This mirrors the DXM-01 fix in
    /// `duplex-metrics` and fgbio's single-strand UMI counting
    /// (`CollectDuplexSeqMetrics` records both halves via `split("-", -1)`), keeping
    /// the two metrics tools internally consistent. Previously the per-position
    /// collection dropped empty components, undercounting single-index designs.
    #[test]
    fn test_count_umis_with_empty_molecule_end_half() -> Result<()> {
        let mut records = Vec::new();

        // Family 1: normal two-part UMI "AAA-TTT" → records "AAA" and "TTT".
        let (r1, r2) = build_test_pair("q1", 0, 100, 200, "AAA-TTT", "1/A");
        records.push(r1);
        records.push(r2);
        // Family 2: empty SECOND half "CCC-" → records "CCC" and "".
        let (r1, r2) = build_test_pair("q2", 0, 1000, 1100, "CCC-", "2/A");
        records.push(r1);
        records.push(r2);
        // Family 3: empty FIRST half "-GGG" → records "" and "GGG".
        let (r1, r2) = build_test_pair("q3", 0, 2000, 2100, "-GGG", "3/A");
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
            threading: crate::commands::common::ThreadingOptions { threads: None },
            scheduler_opts: crate::commands::common::SchedulerOptions::default(),
            queue_memory: crate::commands::common::QueueMemoryOptions::default(),
        };
        cmd.execute("test")?;

        let umi_path = format!("{}.umi_counts.txt", output.display());
        let umi_metrics: Vec<UmiMetric> = DelimFile::default().read_tsv(&umi_path)?;
        let by_umi: std::collections::HashMap<&str, &UmiMetric> =
            umi_metrics.iter().map(|m| (m.umi.as_str(), m)).collect();

        // Present halves ("CCC", "GGG") are counted, and each empty molecule-end half
        // is recorded as an empty-string UMI (CCC- second half + -GGG first half → 2).
        assert_eq!(
            umi_metrics.len(),
            5,
            "expected 5 UMIs (present halves + empty UMI): {:?}",
            umi_metrics.iter().map(|m| m.umi.as_str()).collect::<Vec<_>>()
        );
        for umi in ["", "AAA", "TTT", "CCC", "GGG"] {
            assert!(by_umi.contains_key(umi), "missing UMI {umi:?}");
        }
        assert_eq!(by_umi[""].unique_observations, 2, "empty halves counted as empty UMI");
        assert_eq!(by_umi["CCC"].unique_observations, 1);
        assert_eq!(by_umi["GGG"].unique_observations, 1);

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
            threading: crate::commands::common::ThreadingOptions { threads: None },
            scheduler_opts: crate::commands::common::SchedulerOptions::default(),
            queue_memory: crate::commands::common::QueueMemoryOptions::default(),
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
            threading: crate::commands::common::ThreadingOptions { threads: None },
            scheduler_opts: crate::commands::common::SchedulerOptions::default(),
            queue_memory: crate::commands::common::QueueMemoryOptions::default(),
        };
        assert_eq!(cmd.min_reads, 1);
    }

    /// Both consensus flavours must be rejected — simplex (`cD` without the
    /// `aD`+`bD` pair) and duplex (`aD`+`bD`, with or without `cD`) — while a
    /// plain grouped read must not be. The check runs against a single record —
    /// `process_templates_from_bam` applies it to the first qualifying record of
    /// its own pass — so the test builds records rather than a BAM file.
    #[rstest]
    #[case::simplex_consensus(&[(SamTag::CD, 5)], true)]
    #[case::duplex_consensus(&[(SamTag::AD, 3), (SamTag::BD, 2)], true)]
    #[case::duplex_consensus_with_depth(&[(SamTag::CD, 5), (SamTag::AD, 3), (SamTag::BD, 2)], true)]
    #[case::grouped_not_consensus(&[], false)]
    fn test_reject_consensus_reads(#[case] int_tags: &[(SamTag, i32)], #[case] is_consensus: bool) {
        let mut b = RawSamBuilder::new();
        b.read_name(b"read1")
            .flags(flags::PAIRED | flags::FIRST_SEGMENT)
            .ref_id(0)
            .pos(99)
            .mapq(60)
            .cigar_ops(&[encode_op(0, 100)])
            .sequence(&[b'A'; 100])
            .qualities(&[30u8; 100]);
        for (tag, value) in int_tags {
            b.add_int_tag(*tag, *value);
        }
        let record = b.build();

        assert!(
            is_consensus_guard_record(&record),
            "a paired primary R1 must be the record the guard inspects"
        );

        let result = ensure_not_consensus_record(&record, Path::new("in.bam"));
        assert_eq!(result.is_err(), is_consensus);
        if let Err(e) = result {
            assert!(e.to_string().contains("consensus"), "unexpected error: {e}");
        }
    }

    /// The `--threads` chain path must reject a consensus BAM just like the serial
    /// path — the guard is folded into `MetricsSink` (not a pre-flight re-open, to
    /// preserve stdin/pipe inputs). Builds a mapped, paired-primary pair carrying
    /// the simplex consensus `cD` depth tag (so it survives the pregroup filter and
    /// reaches the sink guard), runs `--threads 4`, and asserts the consensus-BAM
    /// error surfaces through the pipeline's `io::Error` wrapping.
    #[test]
    fn chain_path_rejects_consensus_bam() -> Result<()> {
        // A mapped, paired-primary pair with RX/MI and the simplex consensus `cD`
        // tag. `to_record_buf` round-trips it through the BAM writer/reader so the
        // chain sees the tag on the raw record.
        fn consensus_pair(name: &str) -> (sam::alignment::RecordBuf, sam::alignment::RecordBuf) {
            let seq = vec![b'A'; 100];
            let quals = vec![30u8; 100];
            let cigar = encode_op(0, 100);
            let mut b1 = RawSamBuilder::new();
            b1.read_name(name.as_bytes())
                .flags(flags::PAIRED | flags::FIRST_SEGMENT)
                .ref_id(0)
                .pos(99)
                .mapq(60)
                .cigar_ops(&[cigar])
                .sequence(&seq)
                .qualities(&quals)
                .mate_ref_id(0)
                .mate_pos(199);
            b1.add_string_tag(SamTag::RX, b"ACGT");
            b1.add_string_tag(SamTag::MI, b"1");
            b1.add_int_tag(SamTag::CD, 5);
            let mut b2 = RawSamBuilder::new();
            b2.read_name(name.as_bytes())
                .flags(flags::PAIRED | flags::LAST_SEGMENT)
                .ref_id(0)
                .pos(199)
                .mapq(60)
                .cigar_ops(&[cigar])
                .sequence(&seq)
                .qualities(&quals)
                .mate_ref_id(0)
                .mate_pos(99);
            b2.add_string_tag(SamTag::RX, b"ACGT");
            b2.add_string_tag(SamTag::MI, b"1");
            b2.add_int_tag(SamTag::CD, 5);
            (to_record_buf(b1.build()), to_record_buf(b2.build()))
        }

        let (r1, r2) = consensus_pair("c1");
        let input = create_test_bam(vec![r1, r2])?;
        let output_dir = TempDir::new()?;
        let cmd = SimplexMetrics {
            input: input.path().to_path_buf(),
            output: output_dir.path().join("output"),
            min_reads: 1,
            intervals: None,
            description: None,
            threading: crate::commands::common::ThreadingOptions { threads: Some(4) },
            scheduler_opts: crate::commands::common::SchedulerOptions::default(),
            queue_memory: crate::commands::common::QueueMemoryOptions::default(),
        };
        let message = cmd
            .execute("test")
            .expect_err("consensus BAM must be rejected on the chain path")
            .to_string();
        assert!(
            message.contains("appears to contain consensus sequences"),
            "chain error should reject the consensus BAM: {message}"
        );
        Ok(())
    }

    /// Builds one record with `record_flags` and no consensus tags.
    fn guard_candidate(name: &[u8], record_flags: u16) -> fgumi_raw_bam::RawRecord {
        let mut b = RawSamBuilder::new();
        b.read_name(name)
            .flags(record_flags)
            .ref_id(0)
            .pos(99)
            .mapq(60)
            .cigar_ops(&[encode_op(0, 100)])
            .sequence(&[b'A'; 100])
            .qualities(&[30u8; 100]);
        b.build()
    }

    /// Which record of a template the guard inspects.
    ///
    /// A paired primary R1 is preferred, but a template without one must still
    /// yield a record: a fragment or single-end consensus BAM has no R1 anywhere,
    /// and returning `None` for it is what let it through the guard entirely. Only
    /// a template of nothing but secondary/supplementary records has no answer,
    /// which defers the check to the next template rather than skipping it.
    #[rstest]
    #[case::prefers_paired_r1(
        &[
            (&b"r2"[..], flags::PAIRED | flags::LAST_SEGMENT),
            (&b"r1"[..], flags::PAIRED | flags::FIRST_SEGMENT),
        ],
        Some(&b"r1"[..])
    )]
    #[case::falls_back_to_fragment(&[(&b"frag"[..], 0)], Some(&b"frag"[..]))]
    #[case::falls_back_to_r2_when_r1_absent(
        &[(&b"r2"[..], flags::PAIRED | flags::LAST_SEGMENT)],
        Some(&b"r2"[..])
    )]
    #[case::skips_secondary_and_supplementary(
        &[(&b"sec"[..], flags::SECONDARY), (&b"sup"[..], flags::SUPPLEMENTARY)],
        None
    )]
    #[case::falls_back_past_secondary(
        &[(&b"sec"[..], flags::SECONDARY), (&b"frag"[..], 0)],
        Some(&b"frag"[..])
    )]
    fn test_consensus_guard_record_selection(
        #[case] records: &[(&[u8], u16)],
        #[case] expected_name: Option<&[u8]>,
    ) {
        let records: Vec<_> = records
            .iter()
            .map(|(name, record_flags)| guard_candidate(name, *record_flags))
            .collect();

        let selected = consensus_guard_record(&records)
            .map(|raw| fgumi_raw_bam::read_name(raw.as_ref()).to_vec());

        assert_eq!(selected.as_deref(), expected_name);
    }
}
