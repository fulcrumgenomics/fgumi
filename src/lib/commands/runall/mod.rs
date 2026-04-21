//! CLI for the `fgumi runall` command.
//!
//! Combines multiple processing stages into a single command, streaming data between
//! stages in memory to eliminate intermediate BAM I/O. Supports multiple entry points:
//!
//! - `--start-from group`: reads MI-tagged input, runs consensus + filter.
//! - `--start-from sort`: reads mapped BAM, sorts by template-coordinate, groups by
//!   position, assigns UMIs, calls consensus, filters, and writes output (serial).
//! - `--start-from fastq`: reads unmapped BAM (UMIs already extracted), runs
//!   fastq -> aligner -> zipper -> sort -> group -> consensus -> filter.
//! - `--start-from correct`: reads unmapped BAM with UMIs, runs
//!   correct -> fastq -> aligner -> zipper -> sort -> group -> consensus -> filter.
//! - `--start-from extract`: reads FASTQ files, extracts UMIs, aligns, zippers,
//!   sorts, groups, calls consensus, filters, and writes output.

mod consensus_helpers;
mod extract_helpers;
mod metrics;
pub mod options;
mod parallel;
pub mod plan;
pub mod planner;
pub mod runner;
mod validate;

use std::path::PathBuf;
use std::sync::Arc;
use std::sync::atomic::{AtomicBool, Ordering};

use anyhow::{Context, Result};
use clap::Args;

use crate::commands::command::Command;
use metrics::AtomicOutputGuard;
use options::{
    ConsensusMode, MultiAlignerOptions, MultiCodecOptions, MultiConsensusOptions,
    MultiCorrectOptions, MultiDuplexOptions, MultiExtractOptions, MultiFilterOptions,
    MultiGroupOptions, MultiSortOptions, MultiZipperOptions, StartFrom, StopAfter,
};

/// Progress rendering mode for the runall command. Maps to
/// [`crate::progress::Mode`] via the `From` impl below.
#[derive(clap::ValueEnum, Debug, Clone, Copy, PartialEq, Eq)]
pub enum ProgressMode {
    /// Detect TTY vs non-TTY and render accordingly. Default.
    Auto,
    /// Force the interactive multi-bar dashboard.
    Dashboard,
    /// Force periodic logfmt heartbeat lines suitable for log tailing.
    Heartbeat,
    /// Disable the tracker entirely. Every producer call is a no-op.
    None,
}

impl From<ProgressMode> for crate::progress::Mode {
    fn from(m: ProgressMode) -> Self {
        match m {
            ProgressMode::Auto => crate::progress::Mode::Auto,
            ProgressMode::Dashboard => crate::progress::Mode::Dashboard,
            ProgressMode::Heartbeat => crate::progress::Mode::Heartbeat,
            ProgressMode::None => crate::progress::Mode::None,
        }
    }
}

/// Install a single process-wide SIGINT/SIGTERM handler that sets the shared
/// cancellation flag. `ctrlc::set_handler` errors if called more than once
/// per process, so we guard with a `Once` to make repeated runall invocations
/// (e.g. in a single test process) safe.
///
/// The signal handler stores into `GLOBAL_CANCEL_FLAG`; each runall call takes
/// a clone of this `Arc<AtomicBool>` so cancellation is observed by whichever
/// pipeline is currently running. This is safe because fgumi runs one pipeline
/// at a time per process.
static SIGNAL_HANDLER_INIT: std::sync::Once = std::sync::Once::new();
static GLOBAL_CANCEL_FLAG: std::sync::OnceLock<Arc<AtomicBool>> = std::sync::OnceLock::new();

fn install_signal_handler_once() -> Result<Arc<AtomicBool>> {
    let flag = GLOBAL_CANCEL_FLAG.get_or_init(|| Arc::new(AtomicBool::new(false))).clone();
    let mut install_err: Option<ctrlc::Error> = None;
    SIGNAL_HANDLER_INIT.call_once(|| {
        let handler_flag = flag.clone();
        if let Err(e) = ctrlc::set_handler(move || {
            tracing::warn!("received interrupt signal; cancelling pipeline");
            handler_flag.store(true, Ordering::SeqCst);
        }) {
            install_err = Some(e);
        }
    });
    if let Some(e) = install_err {
        return Err(e).context("Failed to set signal handler");
    }
    // Clear any stale cancel state from a prior invocation in the same process
    // (only relevant in tests that drive runall multiple times). A stale `true`
    // would otherwise make the next pipeline insta-cancel.
    flag.store(false, Ordering::SeqCst);
    Ok(flag)
}

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
/// - `fastq`: Input is an unmapped BAM (UMIs already extracted). Runs fastq -> aligner
///   -> zipper -> sort -> group -> consensus -> filter.
/// - `correct`: Input is an unmapped BAM with UMIs. Runs correct -> fastq -> aligner
///   -> zipper -> sort -> group -> consensus -> filter.
/// - `extract`: Input is FASTQ files. Runs extract -> fastq -> aligner -> zipper -> sort ->
///   group -> consensus -> filter (full pipeline).
const RUNALL_EXAMPLES: &str = r"EXAMPLES:
    # Full pipeline from FASTQ to filtered consensus BAM
    fgumi runall \
        --start-from extract \
        --input r1.fastq.gz r2.fastq.gz \
        --extract::sample SAMPLE1 --extract::library LIB1 \
        --extract::read-structures 8M+T 8M+T \
        --reference ref.fa \
        --output final.bam

    # Skip extraction and alignment; start from a zippered, unsorted mapped BAM
    fgumi runall \
        --start-from sort \
        --input mapped.bam \
        --reference ref.fa \
        --output final.bam

    # Consensus-only (input is already template-coordinate sorted and grouped)
    fgumi runall \
        --start-from group \
        --input grouped.bam \
        --output consensus.bam \
        --consensus simplex

    # Print the resolved plan without running
    fgumi runall --explain \
        --start-from extract \
        --input r1.fq.gz r2.fq.gz \
        --extract::sample SAMPLE1 --extract::library LIB1 \
        --extract::read-structures 8M+T 8M+T \
        --reference ref.fa \
        --output final.bam
";

#[derive(Args, Debug)]
#[command(
    about = "\x1b[38;5;82m[END-TO-END]\x1b[0m     \x1b[36mRun any or all of the consensus-calling workflow in one command\x1b[0m",
    after_help = RUNALL_EXAMPLES,
)]
pub struct Runall {
    /// Input BAM/SAM file path (or FASTQ files for --start-from extract).
    #[arg(short = 'i', long, num_args = 1..)]
    pub input: Vec<PathBuf>,

    /// Output BAM file path.
    #[arg(short = 'o', long)]
    pub output: PathBuf,

    /// Reference FASTA file (required for --start-from sort, fastq, correct, and extract).
    #[arg(short = 'r', long)]
    pub reference: Option<PathBuf>,

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

    /// Progress rendering mode. Default `auto` selects a TTY dashboard when
    /// stderr is interactive and periodic heartbeat lines otherwise.
    #[arg(long, value_enum, default_value_t = ProgressMode::Auto)]
    pub progress: ProgressMode,

    /// Compression level for output BAM (0-9).
    ///
    /// Level 1 is fastest with larger files. Level 9 produces smallest files
    /// but is slowest.
    #[arg(long, default_value = "1")]
    pub compression_level: u32,

    /// Path to write pipeline stage metrics TSV.
    ///
    /// Writes a tab-separated file with columns: `stage`, `wall_time_secs`,
    /// `records_in`, `records_out`.
    #[arg(long)]
    pub metrics: Option<PathBuf>,

    /// Print the resolved pipeline plan as a human-readable report and exit
    /// without executing. Useful for understanding what runall will do with
    /// a given set of inputs and flags before committing to a long run.
    #[arg(long, default_value_t = false)]
    pub explain: bool,

    /// Output prefix for consensus metrics files (family sizes, position group
    /// sizes, UMI counts, yield projections).
    #[arg(long = "consensus-metrics")]
    pub consensus_metrics: Option<PathBuf>,

    /// Intervals file (BED or Picard interval list) for restricting consensus
    /// metrics to target regions.
    #[arg(long = "intervals")]
    pub intervals: Option<PathBuf>,

    /// Methylation mode for consensus calling (em-seq or taps).
    #[arg(long = "methylation-mode", value_enum)]
    pub methylation_mode: Option<crate::commands::common::MethylationModeArg>,

    /// Restore unconverted bases during zipper merge (for bisulfite/EM-seq data).
    #[arg(long = "restore-unconverted-bases", default_value = "false", num_args = 0..=1,
          default_missing_value = "true", action = clap::ArgAction::Set,
          value_parser = crate::commands::common::parse_bool)]
    pub restore_unconverted_bases: bool,

    /// Minimum depth for methylation calls (comma-separated).
    #[arg(long = "min-methylation-depth", value_delimiter = ',')]
    pub min_methylation_depth: Vec<usize>,

    /// Require strand agreement for methylation calls.
    #[arg(long = "require-strand-methylation-agreement", default_value = "false", num_args = 0..=1,
          default_missing_value = "true", action = clap::ArgAction::Set,
          value_parser = crate::commands::common::parse_bool)]
    pub require_strand_methylation_agreement: bool,

    /// Minimum conversion fraction for methylation filtering.
    #[arg(long = "min-conversion-fraction")]
    pub min_conversion_fraction: Option<f64>,

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

impl Command for Runall {
    fn execute(&self, command_line: &str) -> Result<()> {
        // --explain short-circuits: dump the resolved plan without validation,
        // file I/O, atomic-output setup, or signal handling.
        if self.explain {
            let tmp_output = self.output.with_extension("bam.tmp");
            let plan_with_finalize = self.build_plan_explain(command_line, &tmp_output)?;
            println!("{}", plan_with_finalize.plan.explain(command_line));
            return Ok(());
        }

        self.validate()?;

        // Spin up (or not) the progress tracker. The guard is held for the
        // duration of this function so the tracker thread drains after
        // `pipeline_finished` fires below.
        let _progress_guard = crate::progress::init(self.progress.into());
        crate::progress::pipeline_started("runall", None);

        // Set up signal handling: SIGINT/SIGTERM set the cancel flag. Installed
        // exactly once per process so that repeated runall invocations (tests,
        // REPL-style usage) don't fail with ctrlc's "already installed" error.
        let cancel = install_signal_handler_once()?;

        // Atomic output: write to a temp file, rename on success
        let tmp_output = self.output.with_extension("bam.tmp");
        let guard = AtomicOutputGuard { tmp_path: tmp_output.clone() };

        let result = (|| {
            let plan = self.build_plan(command_line, &tmp_output)?;
            let token = crate::runall::engine::cancel::CancelToken::from_arc(cancel.clone());
            plan.run(&token)
        })();

        // A cancelled pipeline returns `Ok(())` because workers/sinks treat
        // cancellation as a clean shutdown, not an error. Detect the cancel
        // flag here and convert to an error so the AtomicOutputGuard discards
        // the `.tmp` file rather than promoting it to the final output path.
        let cancelled = cancel.load(Ordering::SeqCst);
        let result = if cancelled {
            Err(result.err().unwrap_or_else(|| {
                anyhow::anyhow!("pipeline cancelled by signal; output discarded")
            }))
        } else {
            result
        };

        let exec_result = match result {
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
        };

        // Announce pipeline completion (ok / not-ok) before the tracker
        // guard drops at end-of-scope, so the summary reflects the real
        // outcome.
        crate::progress::pipeline_finished(exec_result.is_ok());
        exec_result
    }
}
