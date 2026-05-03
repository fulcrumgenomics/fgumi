//! `fgumi runall` — chain-runner-based end-to-end command (foundation).
//!
//! `fgumi runall` is the umbrella command that runs any consecutive
//! slice of the consensus-calling workflow in a single process. PR 1
//! ships the foundation: the full CLI surface (every per-stage option
//! group), validation, atomic-output handling, signal handling, the
//! [`dispatch::ChainRunner`] trait, the [`dispatch::ChainRegistry`]
//! that picks a runner for each `(start_from, stop_after)` pair, and
//! a single [`chains::NotImplementedRunner`] that bails cleanly for
//! every shape with the documented "not yet implemented" error.
//!
//! Real chain runners land in PR 2 (`(Extract, Extract)` on
//! `unified_pipeline::fastq`) and PR 3 (`(Correct, Correct)` and
//! `(Extract, Correct)` on `unified_pipeline::bam` / `::fastq`).
//! Each new runner is pushed onto the front of the registry built
//! by [`Runall::build_chain_registry`]; the
//! [`chains::NotImplementedRunner`] keeps handling unimplemented
//! shapes from the back.

pub mod chains;
pub mod dispatch;
pub mod infra;
pub mod options;
pub mod validate;

use std::path::PathBuf;
use std::sync::atomic::Ordering;

use anyhow::{Context, Result};
use clap::Args;

use crate::commands::command::Command;
use dispatch::{CancelToken, ChainContext, ChainRegistry, DispatchContext};
use infra::atomic_output::AtomicOutputGuard;
use infra::explain::{NOT_IMPLEMENTED_RUNNER_NAME, render_explain};
use infra::progress::{self, ProgressMode};
use infra::signal::install_signal_handler_once;
use options::{
    ConsensusMode, MultiAlignerOptions, MultiCodecOptions, MultiConsensusOptions,
    MultiCorrectOptions, MultiDuplexOptions, MultiExtractOptions, MultiFilterOptions,
    MultiGroupOptions, MultiSortOptions, MultiZipperOptions, StartFrom, StopAfter,
};

/// Examples shown in `fgumi runall --help`.
const RUNALL_EXAMPLES: &str = r"EXAMPLES:
    # Full pipeline from FASTQ to filtered consensus BAM
    fgumi runall \
        --start-from extract --stop-after filter \
        --input r1.fastq.gz r2.fastq.gz \
        --extract::sample SAMPLE1 --extract::library LIB1 \
        --extract::read-structures 8M+T 8M+T \
        --reference ref.fa \
        --filter::min-reads 1 \
        --output final.bam

    # Skip extraction and alignment; start from a zippered, unsorted mapped BAM
    fgumi runall \
        --start-from sort --stop-after filter \
        --input mapped.bam \
        --reference ref.fa \
        --filter::min-reads 1 \
        --output final.bam

    # Print the resolved chain plan without running anything
    fgumi runall --explain \
        --start-from extract --stop-after filter \
        --input r1.fq.gz r2.fq.gz \
        --extract::sample SAMPLE1 --extract::library LIB1 \
        --extract::read-structures 8M+T 8M+T \
        --reference ref.fa \
        --output final.bam
";

/// `fgumi runall` — run any or all of the consensus-calling workflow in
/// one command.
///
/// See the module docs for the architecture; see [`Runall::execute`]
/// for the lifecycle.
#[derive(Args, Debug)]
#[command(
    about = "[END-TO-END]    Run any or all of the consensus-calling workflow in one command",
    after_help = RUNALL_EXAMPLES,
)]
pub struct Runall {
    /// Input BAM/SAM file path (or FASTQ files for --start-from extract).
    #[arg(short = 'i', long, num_args = 1..)]
    pub input: Vec<PathBuf>,

    /// Output BAM file path.
    #[arg(short = 'o', long)]
    pub output: PathBuf,

    /// Reference FASTA file (required for --start-from sort, fastq, correct, and extract
    /// when the chain reaches the alignment stage).
    #[arg(short = 'r', long)]
    pub reference: Option<PathBuf>,

    /// Entry point stage for the runall pipeline.
    #[arg(long, default_value = "extract", value_enum)]
    pub start_from: StartFrom,

    /// Stop the pipeline after this stage. Must be at or after the
    /// `--start-from` stage.
    #[arg(long, default_value = "filter", value_enum)]
    pub stop_after: StopAfter,

    /// Consensus calling algorithm.
    #[arg(long = "consensus", default_value = "simplex", value_enum)]
    pub consensus_mode: ConsensusMode,

    /// Number of threads for BAM I/O (default: 1 for serial mode).
    #[arg(short = 't', long, default_value = "1")]
    pub threads: usize,

    /// Progress rendering mode.
    #[arg(long, value_enum, default_value_t = ProgressMode::Auto)]
    pub progress: ProgressMode,

    /// Compression level for output BAM (0-9).
    #[arg(long, default_value = "1")]
    pub compression_level: u32,

    /// Pipeline queue memory tuning. Mirrors standalone commands'
    /// `--queue-memory` / `--queue-memory-per-thread` flags. Chain
    /// runners propagate these into the underlying pipeline config
    /// (`FastqPipelineConfig` / `BamPipelineConfig`) so users can tune
    /// memory pressure on the same axis they're used to from the
    /// standalone commands.
    #[command(flatten)]
    pub queue_memory: crate::commands::common::QueueMemoryOptions,

    /// Path to write pipeline stage metrics TSV. Per-chain metrics
    /// support is feature-gated by individual chain runners; the
    /// PR 1 fallback runner ignores this flag.
    #[arg(long)]
    pub metrics: Option<PathBuf>,

    /// Print the resolved pipeline plan as a human-readable report
    /// and exit without executing.
    #[arg(long, default_value_t = false)]
    pub explain: bool,

    /// Output prefix for consensus metrics files.
    #[arg(long = "consensus-metrics")]
    pub consensus_metrics: Option<PathBuf>,

    /// Intervals file (BED or Picard interval list) for restricting
    /// consensus metrics to target regions.
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
    #[arg(long = "require-strand-methylation-agreement", default_value = "false",
          num_args = 0..=1, default_missing_value = "true",
          action = clap::ArgAction::Set,
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

impl Runall {
    /// Build the chain registry.
    ///
    /// PR 2 registers [`chains::ExtractChainRunner`] for the
    /// `(Extract, Extract)` shape; PR 3 will add the correct and
    /// extract→correct runners ahead of it. Real runners are pushed
    /// onto the front so their `supports()` checks fire before the
    /// catch-all [`chains::NotImplementedRunner`] fallback.
    fn build_chain_registry(&self) -> ChainRegistry {
        ChainRegistry::new()
            .register(Box::new(chains::ExtractChainRunner::from_runall(self)))
            .register(Box::new(chains::NotImplementedRunner))
    }

    /// Build a [`DispatchContext`] borrow that captures the CLI
    /// flags relevant to runner selection. Borrows from `self` so
    /// callers control the lifetime.
    fn dispatch_context(&self) -> DispatchContext<'_> {
        DispatchContext {
            start_from: self.start_from,
            stop_after: self.stop_after,
            metrics_path: self.metrics.as_deref(),
        }
    }

    /// Path runall writes its output to before atomically renaming
    /// onto `--output` on success. Sibling-of-output to keep the
    /// rename within the same filesystem.
    fn tmp_output_path(&self) -> PathBuf {
        self.output.with_extension("bam.tmp")
    }

    /// Run the `--explain` short-circuit. No validation, no I/O,
    /// no signal handler — just print and exit `Ok(())`.
    fn run_explain(&self, command_line: &str) -> Result<()> {
        let dispatch = self.dispatch_context();
        let registry = self.build_chain_registry();
        let runner_name = registry
            .select(&dispatch)
            .map_or(NOT_IMPLEMENTED_RUNNER_NAME, dispatch::ChainRunner::name);
        let text = render_explain(self.start_from, self.stop_after, runner_name, command_line);
        println!("{text}");
        Ok(())
    }

    /// Promote the `.tmp` output to the final `--output` path.
    /// Disarms the guard on success; on failure the guard's `Drop`
    /// will clean up the `.tmp` file.
    fn finalize_atomic_output(
        &self,
        result: Result<()>,
        cancelled: bool,
        guard: AtomicOutputGuard,
    ) -> Result<()> {
        // A cancelled pipeline returns `Ok(())` because workers/sinks
        // treat cancellation as a clean shutdown. Detect the cancel
        // flag here and convert to an error so the AtomicOutputGuard
        // discards the `.tmp` file rather than promoting it.
        let result = if cancelled {
            Err(result.err().unwrap_or_else(|| {
                anyhow::anyhow!("pipeline cancelled by signal; output discarded")
            }))
        } else {
            result
        };

        match result {
            Ok(()) => {
                let tmp = guard.tmp_path().to_path_buf();
                std::fs::rename(&tmp, &self.output).with_context(|| {
                    format!("Failed to rename {} to {}", tmp.display(), self.output.display())
                })?;
                guard.disarm();
                Ok(())
            }
            Err(e) => Err(e),
        }
    }
}

impl Runall {
    /// Run the command with a caller-supplied registry. Used by tests
    /// to inject a `TestChainRunner` ahead of the production
    /// `NotImplementedRunner` fallback. Production code calls
    /// [`Runall::execute`], which builds the production registry via
    /// [`Runall::build_chain_registry`].
    fn execute_with_registry(&self, command_line: &str, registry: ChainRegistry) -> Result<()> {
        // --explain short-circuits before any side effects.
        if self.explain {
            return self.run_explain(command_line);
        }

        self.validate()?;

        let _progress_guard = progress::init(self.progress);
        progress::pipeline_started("runall");

        // Install (idempotently) the SIGINT/SIGTERM handler. The
        // returned flag is shared across the chain runner via
        // CancelToken.
        let cancel_flag = install_signal_handler_once()?;
        let cancel = CancelToken::from_flag(cancel_flag.clone());

        // Atomic output: write to a `.tmp` sibling, rename on success.
        let tmp_output = self.tmp_output_path();
        let guard = AtomicOutputGuard::new(tmp_output.clone());

        let dispatch = self.dispatch_context();
        let chain_ctx =
            ChainContext { dispatch, tmp_output: tmp_output.clone(), cancel, command_line };

        let result = registry.dispatch(chain_ctx);
        let cancelled = cancel_flag.load(Ordering::SeqCst);
        let final_result = self.finalize_atomic_output(result, cancelled, guard);

        progress::pipeline_finished(final_result.is_ok());
        final_result
    }
}

impl Command for Runall {
    fn execute(&self, command_line: &str) -> Result<()> {
        let registry = self.build_chain_registry();
        self.execute_with_registry(command_line, registry)
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::commands::runall::dispatch::tests::TestChainRunner;

    /// Builder helper for tests that need a registry with a
    /// `TestChainRunner` ahead of the `NotImplementedRunner` fallback.
    fn registry_with_test_runner_for(
        start: StartFrom,
        stop: StopAfter,
    ) -> (ChainRegistry, std::sync::Arc<TestChainRunner>) {
        let runner = std::sync::Arc::new(TestChainRunner::new("Test", start, stop));
        // Wrap the Arc in a thin shim so the registry can own a Box.
        struct Forwarder(std::sync::Arc<TestChainRunner>);
        impl dispatch::ChainRunner for Forwarder {
            fn name(&self) -> &'static str {
                self.0.name()
            }
            fn supports(&self, ctx: &DispatchContext<'_>) -> bool {
                self.0.supports(ctx)
            }
            fn run(&self, ctx: ChainContext<'_>) -> Result<()> {
                self.0.run(ctx)
            }
        }
        let registry = ChainRegistry::new()
            .register(Box::new(Forwarder(runner.clone())))
            .register(Box::new(chains::NotImplementedRunner));
        (registry, runner)
    }

    #[test]
    fn registry_routes_to_test_runner_when_present() {
        let (registry, runner) =
            registry_with_test_runner_for(StartFrom::Extract, StopAfter::Extract);
        let dir = tempfile::TempDir::new().unwrap();
        let tmp = dir.path().join("out.bam.tmp");
        let cancel =
            CancelToken::from_flag(std::sync::Arc::new(std::sync::atomic::AtomicBool::new(false)));
        let ctx = ChainContext {
            dispatch: DispatchContext {
                start_from: StartFrom::Extract,
                stop_after: StopAfter::Extract,
                metrics_path: None,
            },
            tmp_output: tmp.clone(),
            cancel,
            command_line: "fgumi runall",
        };

        registry.dispatch(ctx).expect("dispatch ok");
        assert_eq!(runner.call_count(), 1);
        assert!(tmp.exists());
    }

    #[test]
    fn registry_falls_back_to_not_implemented_for_unsupported_chains() {
        let (registry, _runner) =
            registry_with_test_runner_for(StartFrom::Extract, StopAfter::Extract);
        let dir = tempfile::TempDir::new().unwrap();
        let tmp = dir.path().join("out.bam.tmp");
        let cancel =
            CancelToken::from_flag(std::sync::Arc::new(std::sync::atomic::AtomicBool::new(false)));
        let ctx = ChainContext {
            dispatch: DispatchContext {
                start_from: StartFrom::Sort,
                stop_after: StopAfter::Filter,
                metrics_path: None,
            },
            tmp_output: tmp,
            cancel,
            command_line: "",
        };

        let err = registry.dispatch(ctx).expect_err("must bail");
        assert!(format!("{err:#}").contains("not yet implemented"));
    }

    /// Build a minimal `Runall` that passes `validate()` for the
    /// `(Group, Group)` chain. Callers need only set `output` (and
    /// optionally other flags) before invoking
    /// [`Runall::execute_with_registry`].
    fn make_runall_for_group_group(input: PathBuf, output: PathBuf) -> Runall {
        use crate::commands::runall::options::{
            AlignerOptions, CodecOptions, ConsensusOptions, CorrectOptions, DuplexOptions,
            ExtractOptions, FilterOptions, GroupOptions, MultiAlignerOptions, MultiCodecOptions,
            MultiConsensusOptions, MultiCorrectOptions, MultiDuplexOptions, MultiExtractOptions,
            MultiFilterOptions, MultiGroupOptions, MultiSortOptions, MultiZipperOptions,
            SortOptions, ZipperOptions,
        };
        Runall {
            input: vec![input],
            output,
            reference: None,
            start_from: StartFrom::Group,
            stop_after: StopAfter::Group,
            consensus_mode: ConsensusMode::Simplex,
            threads: 1,
            progress: ProgressMode::None,
            compression_level: 1,
            metrics: None,
            explain: false,
            consensus_metrics: None,
            intervals: None,
            methylation_mode: None,
            restore_unconverted_bases: false,
            min_methylation_depth: vec![],
            require_strand_methylation_agreement: false,
            min_conversion_fraction: None,
            sort_opts: MultiSortOptions::from(SortOptions::default()),
            group_opts: MultiGroupOptions::from(GroupOptions::default()),
            filter_opts: MultiFilterOptions::from(FilterOptions::default()),
            consensus_opts: MultiConsensusOptions::from(ConsensusOptions::default()),
            extract_opts: MultiExtractOptions::from(ExtractOptions::default()),
            aligner_opts: MultiAlignerOptions::from(AlignerOptions::default()),
            zipper_opts: MultiZipperOptions::from(ZipperOptions::default()),
            correct_opts: MultiCorrectOptions::from(CorrectOptions::default()),
            duplex_opts: MultiDuplexOptions::from(DuplexOptions::default()),
            codec_opts: MultiCodecOptions::from(CodecOptions::default()),
            queue_memory: crate::commands::common::QueueMemoryOptions::default(),
        }
    }

    #[test]
    fn execute_with_test_runner_renames_tmp_to_final_output() {
        let dir = tempfile::TempDir::new().unwrap();
        let input = dir.path().join("input.bam");
        std::fs::write(&input, b"placeholder").unwrap();
        let output = dir.path().join("out.bam");

        let runall = make_runall_for_group_group(input, output.clone());
        let (registry, runner) = registry_with_test_runner_for(StartFrom::Group, StopAfter::Group);

        runall.execute_with_registry("fgumi runall", registry).expect("execute ok");

        assert_eq!(runner.call_count(), 1);
        assert!(output.exists(), "final output must exist after rename");
        // The .tmp must have been disarmed-and-renamed, not left behind.
        assert!(
            !output.with_extension("bam.tmp").exists(),
            "tmp output must be gone after successful rename",
        );
    }

    #[test]
    fn execute_with_unsupported_chain_returns_error_and_cleans_tmp() {
        let dir = tempfile::TempDir::new().unwrap();
        let input = dir.path().join("input.bam");
        std::fs::write(&input, b"placeholder").unwrap();
        let output = dir.path().join("out.bam");

        let runall = make_runall_for_group_group(input, output.clone());
        // Empty registry → falls through to NotImplementedRunner via the
        // `bail!` in ChainRegistry::dispatch() OR via a fallback runner;
        // here we use the production registry which has the
        // NotImplementedRunner fallback.
        let registry = runall.build_chain_registry();

        let err = runall
            .execute_with_registry("fgumi runall", registry)
            .expect_err("not implemented yet");
        assert!(format!("{err:#}").contains("not yet implemented"));
        assert!(!output.exists(), "no final output on error");
        assert!(
            !output.with_extension("bam.tmp").exists(),
            "tmp output must be cleaned up on error",
        );
    }

    /// Test-only `Parser` wrapper so we can exercise `Runall`'s clap
    /// parsing end-to-end. `Runall` itself derives `Args` because it is
    /// composed into the top-level `Subcommand` enum in `main.rs`; for
    /// unit tests we need a real `Parser` to call `try_parse_from`.
    #[derive(Debug, clap::Parser)]
    #[command(name = "runall")]
    struct RunallCli {
        #[command(flatten)]
        runall: Runall,
    }

    /// Common arg prefix that satisfies `Runall::validate` for the
    /// `(Extract, Extract)` chain — used by the queue-memory parsing
    /// tests below so they can assert the parsed values without each
    /// test having to re-spell the entire required argument set.
    const EXTRACT_EXTRACT_BASE_ARGS: &[&str] = &[
        "runall",
        "--start-from",
        "extract",
        "--stop-after",
        "extract",
        "--input",
        "r1.fq.gz",
        "--input",
        "r2.fq.gz",
        "--output",
        "out.bam",
        "--extract::sample",
        "S",
        "--extract::library",
        "L",
    ];

    #[test]
    fn parses_queue_memory_flags() {
        use clap::Parser;
        let mut args: Vec<&str> = EXTRACT_EXTRACT_BASE_ARGS.to_vec();
        args.extend_from_slice(&["--queue-memory", "1.5GB", "--queue-memory-per-thread", "false"]);
        let cli = RunallCli::try_parse_from(args).expect("valid CLI args should parse");
        assert_eq!(cli.runall.queue_memory.queue_memory, "1.5GB");
        assert!(!cli.runall.queue_memory.queue_memory_per_thread);
    }

    #[test]
    fn queue_memory_defaults_match_standalone() {
        use clap::Parser;
        let cli = RunallCli::try_parse_from(EXTRACT_EXTRACT_BASE_ARGS.to_vec())
            .expect("valid CLI args should parse");
        assert_eq!(cli.runall.queue_memory.queue_memory, "768");
        assert!(cli.runall.queue_memory.queue_memory_per_thread);
    }
}
