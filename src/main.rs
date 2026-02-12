#![deny(unsafe_code)]
pub mod commands;
mod version;

use anyhow::Result;
use clap::Parser;
use clap::builder::styling::{AnsiColor, Effects, Styles};

/// Custom styles for CLI help output
const STYLES: Styles = Styles::styled()
    .header(AnsiColor::Green.on_default().effects(Effects::BOLD))
    .usage(AnsiColor::Green.on_default().effects(Effects::BOLD))
    .literal(AnsiColor::Cyan.on_default().effects(Effects::BOLD))
    .placeholder(AnsiColor::Cyan.on_default());
use commands::clip::Clip;
use commands::codec::Codec;
use commands::command::Command;
#[cfg(feature = "compare")]
use commands::compare::Compare;
use commands::correct::CorrectUmis;
use commands::dedup::MarkDuplicates;
use commands::downsample::Downsample;
use commands::duplex::Duplex;
use commands::duplex_metrics::DuplexMetrics;
use commands::extract::Extract;
use commands::fastq::Fastq;
use commands::filter::Filter;
use commands::group::GroupReadsByUmi;
use commands::review::Review;
use commands::simplex::Simplex;
#[cfg(feature = "simulate")]
use commands::simulate::Simulate;
use commands::sort::Sort;
use commands::zipper::Zipper;
use enum_dispatch::enum_dispatch;
use env_logger::Env;
use log::info;

/// Commands that require feature flags to be enabled.
/// Format: (`command_name`, `feature_name`)
const FEATURE_GATED_COMMANDS: &[(&str, &str)] = &[
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
    #[clap(subcommand)]
    subcommand: Subcommand,
}

#[enum_dispatch(Command)]
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

    // Group
    #[command(display_order = 6)]
    Group(GroupReadsByUmi),

    // Deduplication
    #[command(display_order = 7)]
    Dedup(MarkDuplicates),

    // Consensus Calling
    #[command(display_order = 8)]
    Simplex(Simplex),
    #[command(display_order = 9)]
    Duplex(Duplex),
    #[command(display_order = 10)]
    Codec(Codec),

    // Post-consensus
    #[command(display_order = 11)]
    Filter(Filter),
    #[command(display_order = 12)]
    Clip(Clip),
    #[command(display_order = 13)]
    DuplexMetrics(DuplexMetrics),
    #[command(display_order = 14)]
    Review(Review),

    // Utilities
    #[command(display_order = 15)]
    Downsample(Downsample),
    #[cfg(feature = "compare")]
    #[command(display_order = 16)]
    Compare(Compare),
    #[cfg(feature = "simulate")]
    #[command(display_order = 17)]
    Simulate(Simulate),
}

fn main() -> Result<()> {
    #[cfg(feature = "dhat-heap")]
    let _profiler = dhat::Profiler::new_heap();

    env_logger::Builder::from_env(Env::default().default_filter_or("info")).init();

    // Capture full command line BEFORE clap parsing for @PG records
    let command_line = std::env::args().collect::<Vec<_>>().join(" ");

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

    info!("Running fgumi version {}", version::VERSION.as_str());
    args.subcommand.execute(&command_line)
}
