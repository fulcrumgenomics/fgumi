//! Parity tests over fixtures shaped like the fgumi options structs
//! `multi_options` exists to re-expose.
//!
//! The macro's classification rules only matter insofar as they hold for the
//! attribute *combinations* `runall` will actually annotate, so the fixtures
//! below are modelled on the two commands with the widest `#[arg(...)]` surface:
//! `Sort` (`src/lib/commands/sort.rs`) and `GroupReadsByUmi`
//! (`src/lib/commands/group.rs`). Between them they cover a string-form
//! `default_value` on a `Display`-less type, a `long` override paired with a
//! `short`, a repeatable `Vec`, a ranged `value_parser`, a hidden expert flag,
//! `#[arg(skip)]` slots with and without a declared value, a required
//! `value_enum`, and the `num_args = 0..=1` boolean form every fgumi flag uses.
//!
//! These are **fixtures, not mirrors**. `fgumi-cli-macros` sits below the fgumi
//! commands in the dependency graph — they will annotate their options structs
//! with `multi_options`, not the other way round — so its tests deliberately do
//! not reach into `fgumi_lib` for the production definitions: the value types are
//! local stand-ins, and the fields and defaults are not kept in lockstep (the
//! production `Sort`, for one, spells its memory default `"768M"`, not
//! `"768MiB"`). What is under test is the macro's contract, not fgumi's option
//! set — a production option change is expected to leave these tests untouched
//! and is covered by that command's own tests. What warrants a fixture here is a
//! genuinely new *attribute shape* reaching a struct the macro must re-expose.
//!
//! The central assertion is parity: parsing the standalone command and parsing
//! the prefixed `Multi*` companion with the same inputs must produce identical
//! option structs. That covers every classification rule at once — a dropped
//! default, a misclassified `bool`, a leaked alias or a lost `value_parser` all
//! surface as a field that disagrees.

use std::path::PathBuf;

use clap::{Args, CommandFactory, Parser, ValueEnum};
use fgumi_cli_macros::multi_options;

// ─────────────────────────────────────────────────────────────────────────────
// Local stand-ins for the fgumi value types
// ─────────────────────────────────────────────────────────────────────────────

/// Stand-in for `fgumi_lib`'s `MemoryLimit`, which parses "768MiB"-style values
/// and deliberately does **not** implement `Display` — the reason the macro must
/// preserve the string-form `default_value` rather than rewriting it to
/// `default_value_t`.
#[derive(Debug, Clone, PartialEq, Eq, Default)]
pub struct MemoryLimit(pub u64);

fn parse_memory(value: &str) -> Result<MemoryLimit, String> {
    let (digits, multiplier) = match value.strip_suffix("MiB") {
        Some(digits) => (digits, 1024 * 1024),
        None => match value.strip_suffix("GiB") {
            Some(digits) => (digits, 1024 * 1024 * 1024),
            None => (value, 1),
        },
    };
    let n = digits.parse::<u64>().map_err(|e| format!("invalid memory value {value:?}: {e}"))?;
    // `checked_mul` so an oversized `GiB` value maps to a parser error rather than
    // panicking on multiplication overflow in a debug build.
    let bytes = n
        .checked_mul(multiplier)
        .ok_or_else(|| format!("invalid memory value {value:?}: overflows u64"))?;
    Ok(MemoryLimit(bytes))
}

/// Stand-in for the `parse_bool` value parser shared by every fgumi boolean flag.
fn parse_bool(value: &str) -> Result<bool, String> {
    match value {
        "true" | "yes" | "t" | "y" => Ok(true),
        "false" | "no" | "f" | "n" => Ok(false),
        other => Err(format!("invalid boolean {other:?}")),
    }
}

#[derive(Debug, Clone, Copy, PartialEq, Eq, Default, ValueEnum)]
pub enum SpillCodec {
    #[default]
    Zstd,
    Bgzf,
}

#[derive(Debug, Clone, Copy, PartialEq, Eq, Default, ValueEnum)]
pub enum Strategy {
    #[default]
    Identity,
    Edit,
    Adjacency,
    Paired,
}

#[derive(Debug, Clone, Copy, PartialEq, Eq, Default)]
pub struct SortOrderArg;

// ─────────────────────────────────────────────────────────────────────────────
// `SortOptions` — string-form defaults, a `Display`-less value type, a long
// override with a short flag, a repeatable Vec, hidden expert flags, a skip slot
// ─────────────────────────────────────────────────────────────────────────────

#[multi_options("sort", "Sort Options")]
#[derive(Args, Debug, Clone, PartialEq)]
pub struct SortOptions {
    /// Maximum memory for in-memory sorting.
    ///
    /// Default is "768MiB" per thread (matching samtools' 768 MiB). Explicit
    /// values like "512MiB", "1GiB", "4GiB" are per-thread when
    /// --memory-per-thread is enabled (default).
    ///
    /// When the limit is reached, sorted chunks spill to temporary files.
    #[arg(short = 'm', long = "max-memory", default_value = "768MiB", value_parser = parse_memory)]
    pub max_memory: MemoryLimit,

    /// Scale memory limit by thread count (samtools behavior).
    #[arg(long = "memory-per-thread", default_value = "true", num_args = 0..=1, default_missing_value = "true", action = clap::ArgAction::Set, value_parser = parse_bool)]
    pub memory_per_thread: bool,

    /// Temporary directory for intermediate files. Repeatable.
    #[arg(short = 'T', long = "tmp-dir", action = clap::ArgAction::Append)]
    pub tmp_dirs: Vec<PathBuf>,

    /// Compression level for temporary chunk files (0-9).
    #[arg(long = "temp-compression", default_value = "1", value_parser = clap::value_parser!(u32).range(0..=9))]
    pub temp_compression: u32,

    /// Codec for temporary spill files: `zstd` (default) or `bgzf`.
    #[arg(long = "temp-codec", default_value = "zstd")]
    pub temp_codec: SpillCodec,

    /// Worker threads for the accumulation/sort/spill phase (Phase 1).
    #[arg(long = "sort-threads")]
    pub sort_threads: Option<usize>,

    /// Phase-2 spill decompression granularity (expert tuning).
    #[arg(long = "file-granularity", default_value = "false", num_args = 0..=1, default_missing_value = "true", action = clap::ArgAction::Set, value_parser = parse_bool, hide = true)]
    pub file_granularity: bool,

    /// Sort order (chain-builder slot), populated by the command, never by clap.
    #[arg(skip)]
    pub order: SortOrderArg,
}

impl Default for SortOptions {
    fn default() -> Self {
        Self {
            max_memory: MemoryLimit(768 * 1024 * 1024),
            memory_per_thread: true,
            tmp_dirs: Vec::new(),
            temp_compression: 1,
            temp_codec: SpillCodec::Zstd,
            sort_threads: None,
            file_granularity: false,
            order: SortOrderArg,
        }
    }
}

#[derive(Parser, Debug)]
struct StandaloneSort {
    #[command(flatten)]
    opts: SortOptions,
}

#[derive(Parser, Debug)]
struct PrefixedSort {
    #[command(flatten)]
    opts: MultiSortOptions,
}

#[test]
fn sort_options_defaults_match_between_standalone_and_prefixed() {
    let standalone = StandaloneSort::try_parse_from(["fgumi"]).expect("standalone defaults");
    let prefixed = PrefixedSort::try_parse_from(["fgumi"]).expect("prefixed defaults");

    assert_eq!(
        prefixed.opts.validate().expect("validate"),
        standalone.opts,
        "every re-exposed default must equal the standalone command's default"
    );
}

#[test]
fn sort_options_supplied_values_match_between_standalone_and_prefixed() {
    let standalone = StandaloneSort::try_parse_from([
        "fgumi",
        "--max-memory",
        "4GiB",
        "--memory-per-thread",
        "false",
        "--tmp-dir",
        "/scratch/a",
        "--tmp-dir",
        "/scratch/b",
        "--temp-compression",
        "9",
        "--temp-codec",
        "bgzf",
        "--sort-threads",
        "8",
    ])
    .expect("standalone parse");

    let prefixed = PrefixedSort::try_parse_from([
        "fgumi",
        "--sort::max-memory",
        "4GiB",
        "--sort::memory-per-thread",
        "false",
        "--sort::tmp-dir",
        "/scratch/a",
        "--sort::tmp-dir",
        "/scratch/b",
        "--sort::temp-compression",
        "9",
        "--sort::temp-codec",
        "bgzf",
        "--sort::sort-threads",
        "8",
    ])
    .expect("prefixed parse");

    assert_eq!(prefixed.opts.validate().expect("validate"), standalone.opts);
}

/// Assert a parse failed for the stated reason, not merely that it failed. A bare
/// `is_err()` would be satisfied by any unrelated failure a later fixture change
/// introduced.
#[track_caller]
fn assert_parse_error_kind(
    result: Result<PrefixedSort, clap::Error>,
    expected: clap::error::ErrorKind,
    what: &str,
) {
    let err = result.map(|_| ()).expect_err(&format!("{what} must be rejected"));
    assert_eq!(err.kind(), expected, "{what}: expected {expected:?}, got {:?}: {err}", err.kind());
}

#[test]
fn sort_options_value_parsers_are_still_enforced_on_the_prefixed_side() {
    // `temp_compression`'s ranged value_parser and `max_memory`'s custom parser
    // both arrive as `#[arg(...)]` metas the macro copies verbatim. The failure
    // must come from the value parser, not from clap failing to know the flag.
    assert_parse_error_kind(
        PrefixedSort::try_parse_from(["fgumi", "--sort::temp-compression", "10"]),
        clap::error::ErrorKind::ValueValidation,
        "the ranged value_parser must reject 10",
    );
    assert_parse_error_kind(
        PrefixedSort::try_parse_from(["fgumi", "--sort::max-memory", "not-a-size"]),
        clap::error::ErrorKind::ValueValidation,
        "the custom value_parser must reject a malformed size",
    );
}

#[test]
fn sort_options_hidden_flag_stays_hidden_when_re_exposed() {
    // The contract is the `hide` setting on the registered argument, not the
    // absence of a substring from rendered text.
    let command = PrefixedSort::command();
    let arg = command
        .get_arguments()
        .find(|arg| arg.get_long() == Some("sort::file-granularity"))
        .expect("the prefixed flag should still be registered");
    assert!(arg.is_hide_set(), "hide = true must be preserved onto the prefixed flag");
}

#[test]
fn sort_options_short_flags_are_not_propagated() {
    // `-m` and `-T` belong to the standalone command; on runall they would
    // collide with every other stage's short flags.
    assert!(StandaloneSort::try_parse_from(["fgumi", "-m", "1GiB"]).is_ok());
    assert_parse_error_kind(
        PrefixedSort::try_parse_from(["fgumi", "-m", "1GiB"]),
        clap::error::ErrorKind::UnknownArgument,
        "-m",
    );
    assert_parse_error_kind(
        PrefixedSort::try_parse_from(["fgumi", "-T", "/scratch"]),
        clap::error::ErrorKind::UnknownArgument,
        "-T",
    );
}

#[test]
fn every_standalone_long_flag_has_a_prefixed_counterpart() {
    let standalone_flags: Vec<String> = StandaloneSort::command()
        .get_arguments()
        .filter_map(|arg| arg.get_long().map(ToString::to_string))
        .filter(|long| long != "help")
        .collect();
    assert!(!standalone_flags.is_empty(), "fixture should declare long flags");

    let prefixed_flags: Vec<String> = PrefixedSort::command()
        .get_arguments()
        .filter_map(|arg| arg.get_long().map(ToString::to_string))
        .collect();

    // Compared as sets, not one-directionally: a containment check passes even
    // when the prefixed command exposes an extra flag the standalone one never
    // declared.
    let mut expected: Vec<String> =
        standalone_flags.iter().map(|long| format!("sort::{long}")).collect();
    let mut actual: Vec<String> =
        prefixed_flags.iter().filter(|long| *long != "help").cloned().collect();
    expected.sort();
    actual.sort();

    assert_eq!(
        actual, expected,
        "the prefixed command must expose exactly the prefixed counterparts, no more and no less"
    );
}

#[test]
fn skip_slot_is_not_exposed_but_survives_the_round_trip() {
    // The contract is that no `--sort::order` argument is registered — not that
    // the word "order" is absent from the help text, which any future doc
    // comment could break while the contract still holds.
    assert!(
        PrefixedSort::command().get_arguments().all(|arg| arg.get_long() != Some("sort::order")),
        "the skip slot must not become a flag"
    );

    let original = SortOptions { order: SortOrderArg, ..SortOptions::default() };
    let multi: MultiSortOptions = original.clone().into();
    assert_eq!(multi.validate().expect("validate"), original);
}

// ─────────────────────────────────────────────────────────────────────────────
// `GroupOptions` — a required `value_enum`, several skip slots
// ─────────────────────────────────────────────────────────────────────────────

#[multi_options("group", "Group Options")]
#[derive(Args, Debug, Clone, PartialEq)]
pub struct GroupOptions {
    /// Minimum mapping quality.
    #[arg(short = 'm', long = "min-map-q", default_value = "1")]
    pub min_map_q: u8,

    /// Include non-PF reads.
    #[arg(short = 'n', long = "include-non-pf-reads", default_value = "false", num_args = 0..=1, default_missing_value = "true", action = clap::ArgAction::Set, value_parser = parse_bool)]
    pub include_non_pf_reads: bool,

    /// UMI assignment strategy. Required on the standalone command.
    #[arg(short = 's', long = "strategy", value_enum)]
    pub strategy: Strategy,

    /// Minimum UMI length.
    #[arg(short = 'l', long = "min-umi-length")]
    pub min_umi_length: Option<usize>,

    /// Strategy actually used, resolved by the command.
    #[arg(skip = Strategy::Identity)]
    pub effective_strategy: Strategy,

    /// Edits actually used, resolved by the command.
    #[arg(skip)]
    pub effective_edits: u32,
}

impl Default for GroupOptions {
    fn default() -> Self {
        Self {
            min_map_q: 1,
            include_non_pf_reads: false,
            strategy: Strategy::Identity,
            min_umi_length: None,
            effective_strategy: Strategy::Identity,
            effective_edits: 0,
        }
    }
}

#[derive(Parser, Debug)]
struct StandaloneGroup {
    #[command(flatten)]
    opts: GroupOptions,
}

#[derive(Parser, Debug)]
struct PrefixedGroup {
    #[command(flatten)]
    opts: MultiGroupOptions,
}

#[test]
fn group_options_match_when_the_required_strategy_is_supplied() {
    let standalone = StandaloneGroup::try_parse_from(["fgumi", "--strategy", "adjacency"])
        .expect("standalone parse");
    let prefixed = PrefixedGroup::try_parse_from(["fgumi", "--group::strategy", "adjacency"])
        .expect("prefixed parse");

    assert_eq!(prefixed.opts.validate().expect("validate"), standalone.opts);
}

#[test]
fn group_required_field_is_staged_not_enforced_by_clap() {
    // The standalone command refuses to parse without --strategy. On the runall
    // side clap must accept the parse so `validate()` can name the stage and the
    // flag that is missing.
    let err = StandaloneGroup::try_parse_from(["fgumi"])
        .map(|_| ())
        .expect_err("the standalone command requires --strategy");
    assert_eq!(
        err.kind(),
        clap::error::ErrorKind::MissingRequiredArgument,
        "expected clap itself to demand --strategy, got {:?}: {err}",
        err.kind()
    );

    let prefixed = PrefixedGroup::try_parse_from(["fgumi"]).expect("parse must be staged");
    let err = prefixed.opts.validate().expect_err("validate must reject the missing strategy");
    let msg = format!("{err:#}");
    assert!(msg.contains("--group::strategy"), "error should name the prefixed flag: {msg}");
    assert!(msg.contains("required when group is selected"), "error should name the stage: {msg}");
}

// ─────────────────────────────────────────────────────────────────────────────
// Two stages flattened side by side — the shape `runall` actually builds
// ─────────────────────────────────────────────────────────────────────────────

/// `runall` flattens every stage into one command. This is where an un-prefixed
/// flag, a propagated short, or a leaked `next_help_heading` would surface as a
/// clap panic or a mis-filed argument.
#[derive(Parser, Debug)]
struct RunAllLike {
    #[command(flatten)]
    sort: MultiSortOptions,

    #[command(flatten)]
    group: MultiGroupOptions,

    /// A flag belonging to runall itself, declared after both stages.
    #[arg(long)]
    threads: Option<usize>,
}

#[test]
fn two_stages_and_a_parent_flag_coexist_in_one_command() {
    let parsed = RunAllLike::try_parse_from([
        "fgumi",
        "--sort::max-memory",
        "2GiB",
        "--group::strategy",
        "paired",
        "--threads",
        "4",
    ])
    .expect("a two-stage command must build and parse");

    assert_eq!(
        parsed.sort.validate().expect("sort").max_memory,
        MemoryLimit(2 * 1024 * 1024 * 1024)
    );
    assert_eq!(parsed.group.validate().expect("group").strategy, Strategy::Paired);
    assert_eq!(parsed.threads, Some(4));
}

#[test]
fn each_stage_gets_its_own_help_heading_and_the_parent_keeps_its_own() {
    // Asserted on each argument's own `help_heading` rather than on the byte
    // offsets of headings within `render_long_help()`: rendered help is clap's
    // layout concern — it depends on argument ordering and, with the
    // `wrap_help` feature active anywhere in the workspace feature graph, on
    // terminal width — so offset comparisons can hold while the headings are
    // wrong, and break while they are right.
    let command = RunAllLike::command();
    let heading_of = |long: &str| {
        command
            .get_arguments()
            .find(|arg| arg.get_long() == Some(long))
            .unwrap_or_else(|| panic!("--{long} should be registered"))
            .get_help_heading()
            .map(ToString::to_string)
    };

    assert_eq!(
        heading_of("sort::max-memory"),
        Some("Sort Options".to_string()),
        "a sort field must be filed under the sort stage's heading"
    );
    assert_eq!(
        heading_of("group::strategy"),
        Some("Group Options".to_string()),
        "a group field must be filed under the group stage's heading"
    );
    assert_eq!(
        heading_of("threads"),
        None,
        "the parent's own flag must not inherit a stage heading"
    );
}

#[test]
fn stages_sharing_a_flag_name_do_not_collide() {
    // Both stages declare a `-m` short and a min/max flag on the standalone side.
    // Prefixing is what keeps them apart; a leak would make clap panic while
    // building the command above, so reaching this assertion is the test.
    let flags: Vec<String> = RunAllLike::command()
        .get_arguments()
        .filter_map(|arg| arg.get_long().map(ToString::to_string))
        .collect();
    let mut sorted = flags.clone();
    sorted.sort();
    sorted.dedup();
    assert_eq!(sorted.len(), flags.len(), "no flag name may be declared twice: {flags:?}");
}
