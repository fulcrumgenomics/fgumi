//! CLI for the `fgumi runall` command.
//!
//! Combines multiple processing stages into a single command, streaming data between
//! stages in memory to eliminate intermediate BAM I/O. Supports multiple entry points:
//!
//! - `--start-from group`: reads MI-tagged input, runs consensus + filter.
//! - `--start-from sort`: reads mapped BAM, sorts by template-coordinate, groups by
//!   position, assigns UMIs, calls consensus, filters, and writes output (serial).
//! - `--start-from zipper`: reads unmapped BAM + mapped SAM/BAM, merges (zipper),
//!   then sorts, groups, calls consensus, filters, and writes output.
//! - `--start-from extract`: reads FASTQ files, extracts UMIs, aligns, zippers,
//!   sorts, groups, calls consensus, filters, and writes output.

mod aligner_helpers;
mod consensus_helpers;
mod metrics;
pub mod options;
mod parallel;
mod processing_sink;
mod run_extract;
mod run_group;
mod run_sort;
pub mod run_zipper;
mod validate;

use std::path::{Path, PathBuf};
use std::sync::Arc;
use std::sync::atomic::{AtomicBool, Ordering};

use anyhow::{Context, Result};
use clap::Args;
use fgumi_lib::bam_io::create_raw_bam_writer;
use fgumi_lib::logging::log_consensus_summary;
use fgumi_lib::sort::SortStats;
use fgumi_lib::sort::sink::SortedRecordSink;

use crate::commands::command::Command;
use crate::commands::consensus_runner::{ConsensusStatsOps, create_unmapped_consensus_header};
use metrics::{AtomicOutputGuard, PipelineMetrics};
use options::{
    ConsensusMode, MultiAlignerOptions, MultiCodecOptions, MultiConsensusOptions,
    MultiCorrectOptions, MultiDuplexOptions, MultiExtractOptions, MultiFilterOptions,
    MultiGroupOptions, MultiSortOptions, MultiZipperOptions, StartFrom, StopAfter,
};

/// Run all processing stages in a single command.
///
/// Eliminates intermediate BAM I/O by streaming data between stages
/// (sort, group, consensus, filter) in memory.
/// Use `--start-from` to select the entry point based on your input data:
///
/// - `group`: Input is a template-coordinate sorted BAM with MI tags already assigned.
///   Runs consensus calling and filtering only.
/// - `sort`: Input is a mapped BAM (e.g., output of `fgumi zipper`). Runs the full
///   sort -> group -> assign UMIs -> consensus -> filter chain (serial implementation).
/// - `zipper`: Input is a mapped SAM/BAM (`--input`) plus an unmapped BAM (`--unmapped`).
///   Runs zipper merge -> sort -> group -> consensus -> filter.
/// - `extract`: Input is FASTQ files. Runs extract -> fastq -> aligner -> zipper -> sort ->
///   group -> consensus -> filter (full pipeline).
#[derive(Args, Debug)]
#[command(
    about = "\x1b[38;5;82m[END-TO-END]\x1b[0m     \x1b[36mRun any or all of the consensus-calling workflow in one command\x1b[0m"
)]
pub struct Runall {
    /// Input BAM/SAM file path (or FASTQ files for --start-from extract).
    #[arg(short = 'i', long, num_args = 1..)]
    pub input: Vec<PathBuf>,

    /// Output BAM file path.
    #[arg(short = 'o', long)]
    pub output: PathBuf,

    /// Reference FASTA file (required for --start-from sort, zipper, and extract).
    #[arg(short = 'r', long)]
    pub reference: Option<PathBuf>,

    /// Unmapped BAM file with UMI tags (required for --start-from zipper).
    #[arg(long)]
    pub unmapped: Option<PathBuf>,

    /// Entry point stage for the runall pipeline (default: extract for full pipeline).
    #[arg(long, default_value = "extract", value_enum)]
    pub start_from: StartFrom,

    /// Stop the pipeline after this stage (default: filter).
    ///
    /// Must be at or after the `--start-from` stage.  For example,
    /// `--start-from sort --stop-after sort` writes a sorted BAM and exits.
    #[arg(long, default_value = "filter", value_enum)]
    pub stop_after: StopAfter,

    /// Consensus calling algorithm.
    #[arg(long = "consensus", default_value = "simplex", value_enum)]
    pub consensus_mode: ConsensusMode,

    /// Number of threads for BAM I/O (default: 1 for serial mode).
    #[arg(short = 't', long, default_value = "1")]
    pub threads: usize,

    /// Compression level for output BAM (0-9).
    #[arg(long, default_value = "6")]
    pub compression_level: u32,

    /// Path to write pipeline stage metrics TSV.
    ///
    /// Writes a tab-separated file with columns: `stage`, `wall_time_secs`,
    /// `records_in`, `records_out`.
    #[arg(long)]
    pub metrics: Option<PathBuf>,

    // --- Scoped stage options ---
    /// Sort stage options.
    #[command(flatten)]
    pub sort_opts: MultiSortOptions,

    /// Group stage options.
    #[command(flatten)]
    pub group_opts: MultiGroupOptions,

    /// Filter stage options.
    #[command(flatten)]
    pub filter_opts: MultiFilterOptions,

    /// Consensus calling options.
    #[command(flatten)]
    pub consensus_opts: MultiConsensusOptions,

    /// Extract stage options.
    #[command(flatten)]
    pub extract_opts: MultiExtractOptions,

    /// Aligner subprocess options.
    #[command(flatten)]
    pub aligner_opts: MultiAlignerOptions,

    /// Zipper merge options.
    #[command(flatten)]
    pub zipper_opts: MultiZipperOptions,

    /// UMI correction options.
    #[command(flatten)]
    pub correct_opts: MultiCorrectOptions,

    /// Duplex consensus options.
    #[command(flatten)]
    pub duplex_opts: MultiDuplexOptions,

    /// CODEC consensus options.
    #[command(flatten)]
    pub codec_opts: MultiCodecOptions,
}

/// Convert a 2-character SAM tag string to a `[u8; 2]` array.
///
/// # Panics
///
/// Panics if `tag` is not exactly 2 bytes long. Callers must validate tag
/// lengths before calling this (typically in [`Runall::validate`]).
fn tag_bytes(tag: &str) -> [u8; 2] {
    let bytes = tag.as_bytes();
    [bytes[0], bytes[1]]
}

impl Runall {
    /// Run the serial pipeline path: create a [`ProcessingSink`], invoke the caller's
    /// sort function, finish the sink, and log stats + metrics.
    ///
    /// Both `run_from_sort` (serial fallback) and `sort_and_process_records` (serial
    /// fallback) delegate here, passing their specific sort invocation as a closure.
    ///
    /// # Arguments
    ///
    /// * `header` - Input BAM header (used to build consensus/output header).
    /// * `output` - Path for the output BAM file.
    /// * `sort_fn` - Closure that sorts records into the [`ProcessingSink`] and returns
    ///   [`SortStats`]. Callers pass either `sort_template_to_sink` or `sort_from_iter`.
    /// * `cancel` - Cancellation flag (SIGINT/SIGTERM).
    /// * `label` - Human-readable label for log messages (e.g. "from sort").
    /// * `command_line` - Full command line for PG header record.
    fn run_serial_pipeline(
        &self,
        header: &noodles::sam::Header,
        output: &Path,
        sort_fn: impl FnOnce(&mut processing_sink::ProcessingSink) -> Result<SortStats>,
        cancel: &Arc<AtomicBool>,
        label: &str,
        command_line: &str,
    ) -> Result<()> {
        let read_name_prefix = self
            .consensus_opts
            .consensus_read_name_prefix
            .clone()
            .unwrap_or_else(|| fgumi_lib::consensus_caller::make_prefix_from_header(header));

        // For --stop-after group, use the input header (mapped records need reference
        // sequences). For full pipeline (through consensus), consensus reads are unmapped.
        let output_header = if self.stop_after == StopAfter::Group {
            header.clone()
        } else {
            create_unmapped_consensus_header(
                header,
                &self.consensus_opts.consensus_read_group_id,
                &format!("Pipeline consensus ({label})"),
                command_line,
            )?
        };

        let writer =
            create_raw_bam_writer(output, &output_header, self.threads, self.compression_level)?;

        let assigner = self.group_opts.group_strategy.new_assigner(self.group_opts.group_max_edits);
        let caller = self.build_consensus_caller(read_name_prefix)?;
        let umi_tag_bytes = tag_bytes(&self.group_opts.group_umi_tag);

        let mut sink = processing_sink::ProcessingSink::new(
            writer,
            caller,
            assigner,
            umi_tag_bytes,
            self.build_filter_config(),
            self.filter_opts.filter_require_ss_agreement,
            self.build_group_filter_config(),
            self.stop_after,
            cancel.clone(),
        );

        let pipeline_start = std::time::Instant::now();
        let sort_stats = sort_fn(&mut sink)?;
        sink.finish().context("Failed to finish processing sink")?;
        let pipeline_elapsed = pipeline_start.elapsed().as_secs_f64();

        log::info!(
            "Sort stats: {} records read, {} output",
            sort_stats.total_records,
            sort_stats.output_records
        );

        let stats = sink.caller.statistics();
        let consensus_metrics = stats.to_metrics();
        log_consensus_summary(&consensus_metrics);

        log::info!(
            "Total consensus reads written: {}, position groups: {}",
            sink.total_consensus_reads,
            sink.position_groups_processed
        );
        log::info!(
            "Pipeline complete ({label}): {} input records -> {} consensus reads",
            sink.total_input_records,
            sink.total_consensus_reads
        );

        if let Some(ref metrics_path) = self.metrics {
            let mut pm = PipelineMetrics::new();
            pm.add(
                &format!("sort+group+consensus ({label})"),
                pipeline_elapsed,
                sort_stats.total_records,
                sink.total_consensus_reads as u64,
            );
            pm.write_to_file(metrics_path)?;
            log::info!("Pipeline metrics written to {}", metrics_path.display());
        }

        Ok(())
    }
}

impl Command for Runall {
    fn execute(&self, command_line: &str) -> Result<()> {
        self.validate()?;

        // Set up signal handling: SIGINT/SIGTERM set the cancel flag
        let cancel = Arc::new(AtomicBool::new(false));
        let cancel_for_handler = cancel.clone();
        ctrlc::set_handler(move || {
            cancel_for_handler.store(true, Ordering::SeqCst);
        })
        .context("Failed to set signal handler")?;

        // Atomic output: write to a temp file, rename on success
        let tmp_output = self.output.with_extension("bam.tmp");
        let guard = AtomicOutputGuard { tmp_path: tmp_output.clone() };

        let result = match self.start_from {
            StartFrom::Group => self.run_from_group(command_line, &tmp_output, &cancel),
            StartFrom::Sort => self.run_from_sort(command_line, &tmp_output, &cancel),
            StartFrom::Zipper => self.run_from_zipper(command_line, &tmp_output, &cancel),
            StartFrom::Correct => self.run_from_correct(command_line, &tmp_output, &cancel),
            StartFrom::Fastq => self.run_from_fastq(command_line, &tmp_output, &cancel),
            StartFrom::Align => self.run_from_align(command_line, &tmp_output, &cancel),
            StartFrom::Extract => self.run_from_extract(command_line, &tmp_output, &cancel),
        };

        match result {
            Ok(()) => {
                // Rename temp file to final output
                std::fs::rename(&tmp_output, &self.output).with_context(|| {
                    format!(
                        "Failed to rename {} to {}",
                        tmp_output.display(),
                        self.output.display()
                    )
                })?;
                // Disarm the guard so it doesn't delete the (now-renamed) temp file
                guard.disarm();
                Ok(())
            }
            Err(e) => {
                // Guard's Drop will clean up the temp file
                Err(e)
            }
        }
    }
}
