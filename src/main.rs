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
use fgumi_lib::commands::copy_umi::CopyUmi;
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
use fgumi_lib::commands::retag::Retag;
use fgumi_lib::commands::review::Review;
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

/// Names the full set of spellings boolean flags accept.
///
/// Per-flag help advertises only `[<true|false>]`: `BoolishValueParser` takes
/// twelve literals, and naming all of them on every boolean flag would be
/// unreadable, so the complete set is stated once here instead.
const BOOL_FLAG_HELP: &str = "Boolean flags shown as [<true|false>] take an optional value, and \
                              mean true when given with none. They accept true, t, yes, y, on, 1 \
                              and false, f, no, n, off, 0, in any capitalization.";

#[derive(Parser, Debug)]
#[command(styles = STYLES, after_help = BOOL_FLAG_HELP)]
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

    // Post-consensus
    #[command(display_order = 12)]
    Filter(Filter),
    #[command(display_order = 13)]
    Clip(Clip),
    #[cfg(feature = "consensus")]
    #[command(display_order = 14)]
    DuplexMetrics(DuplexMetrics),
    #[cfg(feature = "consensus")]
    #[command(display_order = 15)]
    SimplexMetrics(SimplexMetrics),
    #[command(display_order = 16)]
    Review(Review),

    // Utilities
    #[command(display_order = 17)]
    Downsample(Downsample),
    #[command(display_order = 18)]
    CopyUmi(CopyUmi),
    #[command(display_order = 19)]
    Retag(Retag),
    #[cfg(feature = "compare")]
    #[command(display_order = 20)]
    Compare(Compare),
    #[cfg(feature = "simulate")]
    #[command(display_order = 21)]
    Simulate(Simulate),
}

impl Subcommand {
    fn execute(&self, command_line: &str) -> Result<()> {
        match self {
            Self::Extract(cmd) => cmd.execute(command_line),
            Self::Correct(cmd) => cmd.execute(command_line),
            Self::Fastq(cmd) => cmd.execute(command_line),
            Self::Zipper(cmd) => cmd.execute(command_line),
            Self::Sort(cmd) => cmd.execute(command_line),
            Self::Merge(cmd) => cmd.execute(command_line),
            Self::Group(cmd) => cmd.execute(command_line),
            Self::Dedup(cmd) => cmd.execute(command_line),
            #[cfg(feature = "consensus")]
            Self::Simplex(cmd) => cmd.execute(command_line),
            #[cfg(feature = "consensus")]
            Self::Duplex(cmd) => cmd.execute(command_line),
            #[cfg(feature = "consensus")]
            Self::Codec(cmd) => cmd.execute(command_line),
            Self::Filter(cmd) => cmd.execute(command_line),
            Self::Clip(cmd) => cmd.execute(command_line),
            #[cfg(feature = "consensus")]
            Self::DuplexMetrics(cmd) => cmd.execute(command_line),
            #[cfg(feature = "consensus")]
            Self::SimplexMetrics(cmd) => cmd.execute(command_line),
            Self::Review(cmd) => cmd.execute(command_line),
            Self::Downsample(cmd) => cmd.execute(command_line),
            Self::CopyUmi(cmd) => cmd.execute(command_line),
            Self::Retag(cmd) => cmd.execute(command_line),
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
    if let Err(ref e) = result
        && e.downcast_ref::<CompareMismatch>().is_some()
    {
        std::process::exit(1);
    }

    result
}

#[cfg(test)]
mod tests {
    use super::Args;
    use clap::CommandFactory;

    /// Every literal `BoolishValueParser` accepts, lower-cased.
    ///
    /// clap marks ten of the twelve hidden, but `Arg::get_possible_values`
    /// returns all of them — which is what lets `is_boolish` decide exactly
    /// whether a flag is a boolean instead of guessing from its shape.
    const BOOLISH_LITERALS: [&str; 12] =
        ["true", "t", "yes", "y", "on", "1", "false", "f", "no", "n", "off", "0"];

    /// The placeholder every boolean flag must render.
    const BOOL_VALUE_NAME: &str = "true|false";

    /// Commands the walk must reach, so a refactor that stops recursing into
    /// subcommands fails loudly instead of passing vacuously.
    const EXPECTED_COMMANDS: [&str; 6] = ["sort", "filter", "group", "dedup", "clip", "correct"];

    /// A boolean flag located during the walk, carrying everything the
    /// declaration assertions need so they never re-walk the tree.
    ///
    /// `command` is a path rather than the `clap::Command` itself because the
    /// two tests that need one — the rendered-help and bare-parse tests —
    /// mutate their copy, so they resolve a fresh one through
    /// [`resolve_command`].
    struct BoolFlag {
        /// Full command path, e.g. `fgumi compare bams`.
        command: String,
        /// The flag's id, e.g. `memory_per_thread`.
        id: String,
        /// The flag's long spelling, e.g. `memory-per-thread`.
        long: Option<String>,
        /// `value_name`s the flag declares, if any.
        value_names: Option<Vec<String>>,
        /// Whether the flag suppresses clap's `[possible values: ...]` line.
        hides_possible_values: bool,
        /// Whether the flag is hidden from help, and so renders no placeholder.
        hidden: bool,
    }

    /// Whether `arg` is parsed by `BoolishValueParser`.
    ///
    /// Keyed on the parser's own accepted set, not on a structural proxy such as
    /// `num_args = 0..=1`: an optional-value flag that is not a boolean — an
    /// optional path, an optional count — has that same shape, and holding it to
    /// the boolean placeholder would be wrong.
    fn is_boolish(arg: &clap::Arg) -> bool {
        let accepted: Vec<String> =
            arg.get_possible_values().iter().map(|v| v.get_name().to_ascii_lowercase()).collect();
        accepted.len() == BOOLISH_LITERALS.len()
            && accepted.iter().all(|value| BOOLISH_LITERALS.contains(&value.as_str()))
    }

    /// Collect every boolean flag in the CLI, recursing through subcommands so
    /// nested commands (`compare bams`, `simulate fastq`) are covered too.
    fn collect_bool_flags(command: &clap::Command, path: &str, out: &mut Vec<BoolFlag>) {
        for arg in command.get_arguments().filter(|arg| is_boolish(arg)) {
            out.push(BoolFlag {
                command: path.to_string(),
                id: arg.get_id().to_string(),
                long: arg.get_long().map(ToString::to_string),
                value_names: arg
                    .get_value_names()
                    .map(|names| names.iter().map(ToString::to_string).collect()),
                hides_possible_values: arg.is_hide_possible_values_set(),
                hidden: arg.is_hide_set(),
            });
        }
        for sub in command.get_subcommands() {
            collect_bool_flags(sub, &format!("{path} {}", sub.get_name()), out);
        }
    }

    /// Re-resolve a `BoolFlag::command` path back to the `clap::Command` it
    /// names, e.g. `fgumi compare bams` to the `bams` subcommand.
    fn resolve_command(path: &str) -> clap::Command {
        let mut command = Args::command();
        for name in path.split_whitespace().skip(1) {
            let sub = command
                .get_subcommands()
                .find(|sub| sub.get_name() == name)
                .unwrap_or_else(|| panic!("no subcommand `{name}` under `{path}`"))
                .clone();
            command = sub;
        }
        command
    }

    /// The long spelling of `flag`, which every boolean flag is required to
    /// have: a bare boolean is only usable as `--flag`, and a short-only
    /// boolean could not render the placeholder this module asserts.
    fn long_of(flag: &BoolFlag) -> &str {
        flag.long.as_deref().unwrap_or_else(|| {
            panic!("`{} --{}` is a boolean flag with no long spelling", flag.command, flag.id)
        })
    }

    /// Walk the CLI once and assert it was actually walked.
    fn bool_flags() -> Vec<BoolFlag> {
        let mut flags = Vec::new();
        collect_bool_flags(&Args::command(), "fgumi", &mut flags);

        // Guard against a vacuous pass: if `is_boolish` or the subcommand
        // recursion stops matching, every loop over this would check nothing.
        assert!(
            flags.len() >= 40,
            "expected the CLI to expose its boolean flags, got {}",
            flags.len()
        );
        for expected in EXPECTED_COMMANDS {
            assert!(
                flags.iter().any(|flag| flag.command.ends_with(expected)),
                "walk never reached `{expected}`",
            );
        }
        flags
    }

    /// Every boolean flag must declare `value_name = "true|false"`.
    ///
    /// Without it clap derives the placeholder from the field name, so
    /// `--memory-per-thread` renders as `--memory-per-thread
    /// [<MEMORY_PER_THREAD>]` — which reads as a request for a per-thread memory
    /// size and invites `--memory-per-thread 10G`. `QueueMemoryOptions` is
    /// `#[command(flatten)]`-ed into most commands, so the placeholder has to
    /// carry this on its own in each of them.
    #[test]
    fn test_bool_args_name_their_accepted_values() {
        for flag in bool_flags() {
            assert_eq!(
                flag.value_names.as_deref(),
                Some([BOOL_VALUE_NAME.to_string()].as_slice()),
                "`{} --{}` must declare `value_name = \"{BOOL_VALUE_NAME}\"`; without it clap \
                 renders the upper-cased field name, which reads as a value to supply rather \
                 than a bool",
                flag.command,
                flag.id,
            );
        }
    }

    /// Every boolean flag must suppress clap's `[possible values: ...]` line.
    ///
    /// `BoolishValueParser` hides ten of its twelve literals, so clap would
    /// advertise `true, false` alone — naming two of the twelve accepted
    /// spellings is worse than naming none, which is why the placeholder carries
    /// the contract instead.
    #[test]
    fn test_bool_args_hide_their_partial_possible_values() {
        for flag in bool_flags() {
            assert!(
                flag.hides_possible_values,
                "`{} --{}` must set `hide_possible_values = true`; clap hides ten of the twelve \
                 boolish literals, so help would advertise only `true, false`",
                flag.command, flag.id,
            );
        }
    }

    /// Since per-flag help advertises only two of the twelve accepted
    /// spellings, the top-level help has to name the rest — and has to keep
    /// naming them if the parser's set ever changes.
    #[test]
    fn test_top_level_help_names_every_accepted_spelling() {
        let help = super::BOOL_FLAG_HELP.to_ascii_lowercase();
        for literal in BOOLISH_LITERALS {
            assert!(
                help.split(|c: char| !c.is_ascii_alphanumeric()).any(|word| word == literal),
                "top-level help does not name the accepted spelling `{literal}`; help was: {help}",
            );
        }
        assert!(
            Args::command().render_help().to_string().contains("Boolean flags shown as"),
            "`fgumi --help` does not carry the boolean-flag note",
        );
    }

    /// The declared `value_name` is only half the contract — assert what a user
    /// actually sees, for every flag rather than once per command.
    ///
    /// The rendered placeholder also depends on `num_args`: `0..=1` renders
    /// `[<true|false>]`, while a required `1` renders `<true|false>` and would
    /// silently break bare `--verify` while still satisfying the declaration
    /// test above. Asserting per flag rather than per command matters because
    /// most commands own several booleans, and a command-level assertion passes
    /// as long as *one* of them renders the placeholder. Flags marked
    /// `hide = true` are exempt: they render nothing at all, so there is no
    /// placeholder to hold them to — the declaration tests above still cover
    /// them, as does the bare-parse test below.
    #[test]
    fn test_bool_args_render_an_optional_boolean_placeholder() {
        let mut checked = 0_usize;
        for flag in bool_flags().into_iter().filter(|flag| !flag.hidden) {
            checked += 1;
            let help = resolve_command(&flag.command).render_help().to_string();
            let expected = format!("--{} [<{BOOL_VALUE_NAME}>]", long_of(&flag));
            assert!(
                help.contains(&expected),
                "`{} --help` renders no `{expected}`; a boolean flag must render its value as \
                 optional, since a required `<{BOOL_VALUE_NAME}>` breaks the bare `--{}` form; \
                 help was:\n{help}",
                flag.command,
                long_of(&flag),
            );
        }
        // Guard against a vacuous pass: if every boolean flag were hidden, or
        // the walk stopped finding them, the loop above would assert nothing.
        assert!(checked >= 40, "expected the CLI to expose visible boolean flags, got {checked}");
    }

    /// Passing a boolean flag with no value must mean `true`.
    ///
    /// `num_args = 0..=1` alone only makes the value optional — without
    /// `default_missing_value = "true"` clap errors on the bare form, so the
    /// documented "mean true when given with none" contract needs a parse to
    /// hold it. Required arguments are relaxed first so each command can be
    /// probed with its boolean flag alone.
    #[test]
    fn test_bool_args_parse_their_bare_long_option_as_true() {
        for flag in bool_flags() {
            let long = long_of(&flag).to_string();
            let command = resolve_command(&flag.command).mut_args(|arg| arg.required(false));
            let matches = command
                .try_get_matches_from(["probe", &format!("--{long}")])
                .unwrap_or_else(|error| {
                    panic!("`{} --{long}` (no value) failed to parse: {error}", flag.command)
                });
            assert_eq!(
                matches.get_one::<bool>(&flag.id),
                Some(&true),
                "`{} --{long}` (no value) must parse as `true`; set \
                 `default_missing_value = \"true\"`",
                flag.command,
            );
        }
    }
}
