#![deny(unsafe_code)]

use anyhow::Result;
use clap::Parser;
use clap::builder::styling::{AnsiColor, Effects, Styles};

/// Custom styles for CLI help output
const STYLES: Styles = Styles::styled()
    .header(AnsiColor::Green.on_default().effects(Effects::BOLD))
    .usage(AnsiColor::Green.on_default().effects(Effects::BOLD))
    .literal(AnsiColor::Cyan.on_default().effects(Effects::BOLD))
    .placeholder(AnsiColor::Cyan.on_default());
use env_logger::Env;
use fgumi_lib::commands::clip::Clip;
#[cfg(feature = "consensus")]
use fgumi_lib::commands::codec::Codec;
use fgumi_lib::commands::command::Command;
#[cfg(feature = "compare")]
use fgumi_lib::commands::compare::Compare;
#[cfg(feature = "compare")]
use fgumi_lib::commands::compare::CompareMismatch;
use fgumi_lib::commands::correct::CorrectUmis;
use fgumi_lib::commands::dedup::MarkDuplicates;
use fgumi_lib::commands::downsample::Downsample;
#[cfg(feature = "consensus")]
use fgumi_lib::commands::duplex::Duplex;
#[cfg(feature = "consensus")]
use fgumi_lib::commands::duplex_metrics::DuplexMetrics;
use fgumi_lib::commands::extract::Extract;
use fgumi_lib::commands::fastq::Fastq;
use fgumi_lib::commands::filter::Filter;
use fgumi_lib::commands::group::GroupReadsByUmi;
use fgumi_lib::commands::merge::Merge;
use fgumi_lib::commands::review::Review;
#[cfg(feature = "consensus")]
use fgumi_lib::commands::runall::RunAll;
#[cfg(feature = "consensus")]
use fgumi_lib::commands::simplex::Simplex;
#[cfg(feature = "consensus")]
use fgumi_lib::commands::simplex_metrics::SimplexMetrics;
#[cfg(feature = "simulate")]
use fgumi_lib::commands::simulate::Simulate;
use fgumi_lib::commands::sort::Sort;
use fgumi_lib::commands::zipper::Zipper;
use log::info;

/// Commands that require feature flags to be enabled.
/// Format: (`command_name`, `feature_name`)
const FEATURE_GATED_COMMANDS: &[(&str, &str)] = &[
    #[cfg(not(feature = "consensus"))]
    ("simplex", "consensus"),
    #[cfg(not(feature = "consensus"))]
    ("simplex-metrics", "consensus"),
    #[cfg(not(feature = "consensus"))]
    ("duplex", "consensus"),
    #[cfg(not(feature = "consensus"))]
    ("duplex-metrics", "consensus"),
    #[cfg(not(feature = "consensus"))]
    ("codec", "consensus"),
    // `runall` fuses the consensus stages, so it too requires the `consensus`
    // feature; it was previously omitted from this hint list (S5d-009).
    #[cfg(not(feature = "consensus"))]
    ("runall", "consensus"),
    #[cfg(not(feature = "compare"))]
    ("compare", "compare"),
    #[cfg(not(feature = "simulate"))]
    ("simulate", "simulate"),
];

#[cfg(feature = "dhat-heap")]
#[global_allocator]
static GLOBAL: dhat::Alloc = dhat::Alloc;

#[cfg(not(feature = "dhat-heap"))]
#[global_allocator]
static GLOBAL: mimalloc::MiMalloc = mimalloc::MiMalloc;

#[derive(Parser, Debug)]
#[command(styles = STYLES)]
struct Args {
    /// Enable verbose (debug-level) logging. Equivalent to setting `RUST_LOG=debug`.
    #[arg(short, long, global = true)]
    verbose: bool,
    #[clap(subcommand)]
    subcommand: Subcommand,
}

#[derive(Parser, Debug)]
#[command(version)]
#[allow(clippy::large_enum_variant)]
enum Subcommand {
    // Grouping
    #[command(display_order = 1)]
    Extract(Extract),
    #[command(display_order = 2)]
    Correct(CorrectUmis),

    // Alignment
    #[command(display_order = 3)]
    Fastq(Fastq),
    #[command(display_order = 4)]
    Zipper(Zipper),
    #[command(display_order = 5)]
    Sort(Sort),
    #[command(display_order = 6)]
    Merge(Merge),

    // Group
    #[command(display_order = 7)]
    Group(GroupReadsByUmi),

    // Deduplication
    #[command(display_order = 8)]
    Dedup(MarkDuplicates),

    // Consensus Calling
    #[cfg(feature = "consensus")]
    #[command(display_order = 9)]
    Simplex(Simplex),
    #[cfg(feature = "consensus")]
    #[command(display_order = 10)]
    Duplex(Duplex),
    #[cfg(feature = "consensus")]
    #[command(display_order = 11)]
    Codec(Codec),
    #[cfg(feature = "consensus")]
    #[command(display_order = 12)]
    Runall(RunAll),

    // Post-consensus
    #[command(display_order = 13)]
    Filter(Filter),
    #[command(display_order = 14)]
    Clip(Clip),
    #[cfg(feature = "consensus")]
    #[command(display_order = 15)]
    DuplexMetrics(DuplexMetrics),
    #[cfg(feature = "consensus")]
    #[command(display_order = 16)]
    SimplexMetrics(SimplexMetrics),
    #[command(display_order = 17)]
    Review(Review),

    // Utilities
    #[command(display_order = 18)]
    Downsample(Downsample),
    #[cfg(feature = "compare")]
    #[command(display_order = 19)]
    Compare(Compare),
    #[cfg(feature = "simulate")]
    #[command(display_order = 20)]
    Simulate(Simulate),
}

impl Subcommand {
    fn execute(&self, command_line: &str) -> Result<()> {
        match self {
            Self::Extract(cmd) => cmd.execute(command_line),
            Self::Correct(cmd) => cmd.execute(command_line),
            Self::Fastq(cmd) => cmd.execute(command_line),
            Self::Zipper(cmd) => cmd.execute(command_line),
            Self::Sort(cmd) => fgumi_lib::commands::sort::execute_sort_command(cmd, command_line),
            Self::Merge(cmd) => cmd.execute(command_line),
            Self::Group(cmd) => cmd.execute(command_line),
            Self::Dedup(cmd) => cmd.execute(command_line),
            #[cfg(feature = "consensus")]
            Self::Simplex(cmd) => cmd.execute(command_line),
            #[cfg(feature = "consensus")]
            Self::Duplex(cmd) => cmd.execute(command_line),
            #[cfg(feature = "consensus")]
            Self::Codec(cmd) => cmd.execute(command_line),
            #[cfg(feature = "consensus")]
            Self::Runall(cmd) => cmd.execute(command_line),
            Self::Filter(cmd) => cmd.execute(command_line),
            Self::Clip(cmd) => cmd.execute(command_line),
            #[cfg(feature = "consensus")]
            Self::DuplexMetrics(cmd) => cmd.execute(command_line),
            #[cfg(feature = "consensus")]
            Self::SimplexMetrics(cmd) => cmd.execute(command_line),
            Self::Review(cmd) => cmd.execute(command_line),
            Self::Downsample(cmd) => cmd.execute(command_line),
            #[cfg(feature = "compare")]
            Self::Compare(cmd) => cmd.execute(command_line),
            #[cfg(feature = "simulate")]
            Self::Simulate(cmd) => cmd.execute(command_line),
        }
    }
}

fn main() -> Result<()> {
    #[cfg(feature = "dhat-heap")]
    let _profiler = dhat::Profiler::new_heap();

    // Capture full command line BEFORE clap parsing for @PG records
    let command_line = std::env::args().collect::<Vec<_>>().join(" ");

    // Install the git-augmented version string so `fgumi sort`'s @PG record
    // (built in the framework-light `fgumi-sort-cli` crate, which cannot see
    // the umbrella's `built.rs`) matches the rest of the toolchain. Set once
    // at startup; idempotent — the boolean return is intentionally ignored.
    let _ = fgumi_sort_cli::version::set_version_override(fgumi_lib::version::VERSION.clone());

    let args = match Args::try_parse() {
        Ok(args) => args,
        Err(e) => {
            // Check if this is an unrecognized subcommand that's behind a feature flag
            if let clap::error::ErrorKind::InvalidSubcommand = e.kind() {
                let err_str = e.to_string();
                for (cmd, feature) in FEATURE_GATED_COMMANDS {
                    if err_str.contains(&format!("'{cmd}'")) {
                        eprintln!(
                            "error: The '{cmd}' command requires the '{feature}' feature.\n\n\
                             Rebuild with: cargo build --release --features {feature}\n"
                        );
                        std::process::exit(2);
                    }
                }
            }
            e.exit();
        }
    };

    let default_level = if args.verbose { "debug" } else { "info" };
    env_logger::Builder::from_env(Env::default().default_filter_or(default_level)).init();

    info!("Running fgumi version {}", fgumi_lib::version::VERSION.as_str());

    let result = args.subcommand.execute(&command_line);

    #[cfg(feature = "compare")]
    if let Err(ref e) = result {
        if e.downcast_ref::<CompareMismatch>().is_some() {
            std::process::exit(1);
        }
    }

    result
}
