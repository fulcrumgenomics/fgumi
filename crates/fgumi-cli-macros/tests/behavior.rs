//! Behavioral contract tests for `multi_options`.
//!
//! Where `smoke.rs` covers the three field-kind branches, this file pins the
//! properties that make the re-exposed flags *faithful* to the standalone
//! command: the declared clap default is the one that applies, help headings
//! do not leak onto the parent command, short/long help stay split, `cfg`
//! gating survives, aliases are namespaced, and the
//! `From` + `validate()` pair is a genuine round trip.
//!
//! Each fixture is built around an attribute shape a real fgumi options struct
//! actually uses; `real_world.rs` assembles those shapes into larger fixtures
//! modelled on the `Sort` and `GroupReadsByUmi` commands.

use std::path::PathBuf;

use clap::{Args, CommandFactory, Parser};
use fgumi_cli_macros::multi_options;

// ─────────────────────────────────────────────────────────────────────────────
// Declared clap defaults are authoritative
// ─────────────────────────────────────────────────────────────────────────────

/// The `Default` impl below deliberately DISAGREES with every `default_value*`
/// attribute. A macro that re-derives defaults from `Struct::default()` would
/// hand the prefixed flag the `Default` value; preserving the original
/// attribute verbatim keeps the standalone and prefixed flags identical.
#[multi_options("drift", "Drift Options")]
#[derive(Args, Debug, Clone, PartialEq)]
pub struct DriftOptions {
    /// String-form default.
    #[arg(long, default_value = "42")]
    pub string_defaulted: u32,

    /// Typed default.
    #[arg(long, default_value_t = 7)]
    pub typed_defaulted: u32,

    /// OS-string default — `default_value_os` must count as "has a default",
    /// not be misclassified as a required field.
    #[arg(long, default_value_os = "/tmp/os-default")]
    pub os_defaulted: PathBuf,
}

impl Default for DriftOptions {
    fn default() -> Self {
        // Every value here disagrees with the clap default above.
        Self {
            string_defaulted: 999,
            typed_defaulted: 999,
            os_defaulted: PathBuf::from("/tmp/wrong"),
        }
    }
}

#[derive(Parser, Debug)]
struct DriftWrapper {
    #[command(flatten)]
    opts: MultiDriftOptions,
}

#[test]
fn declared_clap_defaults_win_over_a_disagreeing_default_impl() {
    let parsed = DriftWrapper::try_parse_from(["test-prog"]).expect("parse with all defaults");
    let opts = parsed.opts.validate().expect("validate");

    assert_eq!(opts.string_defaulted, 42, "string-form default_value must be preserved verbatim");
    assert_eq!(opts.typed_defaulted, 7, "default_value_t must be preserved verbatim");
    assert_eq!(
        opts.os_defaulted,
        PathBuf::from("/tmp/os-default"),
        "default_value_os must be preserved verbatim"
    );
}

#[test]
fn default_value_os_field_is_not_treated_as_required() {
    // A field misclassified as required would be wrapped in `Option<T>` and
    // rejected by `validate()` when omitted; it also would not carry the
    // default. Parsing with no flags at all proves it is optional.
    let parsed = DriftWrapper::try_parse_from(["test-prog"]).expect("parse");
    assert!(parsed.opts.validate().is_ok(), "default_value_os field must not be required");
}

/// Asserted per argument rather than against the rendered help text: a bare
/// `help.contains("42")` also passes when the value appears anywhere else in
/// the help, including on a different flag.
#[rstest::rstest]
#[case::string_form("drift::string-defaulted", "42")]
#[case::typed("drift::typed-defaulted", "7")]
#[case::os_string("drift::os-defaulted", "/tmp/os-default")]
fn prefixed_defaults_are_advertised_in_help(#[case] long: &str, #[case] expected: &str) {
    let command = DriftWrapper::command();
    let arg = command
        .get_arguments()
        .find(|arg| arg.get_long() == Some(long))
        .unwrap_or_else(|| panic!("--{long} should be registered"));
    let defaults: Vec<String> =
        arg.get_default_values().iter().map(|value| value.to_string_lossy().into_owned()).collect();
    assert_eq!(
        defaults,
        vec![expected.to_string()],
        "--{long} must advertise its declared default"
    );
}

// ─────────────────────────────────────────────────────────────────────────────
// Help headings do not leak onto the parent command
// ─────────────────────────────────────────────────────────────────────────────

#[multi_options("stage", "Stage Options")]
#[derive(Args, Debug, Clone, PartialEq)]
pub struct StageOptions {
    /// A stage knob.
    #[arg(long, default_value_t = 1)]
    pub knob: u32,
}

impl Default for StageOptions {
    fn default() -> Self {
        Self { knob: 1 }
    }
}

#[test]
fn parent_args_declared_after_a_flattened_multi_keep_their_own_heading() {
    /// A parent command that declares its own argument *after* flattening the
    /// generated struct. With a struct-level `next_help_heading`, clap applies
    /// the heading to every argument registered after it — silently filing
    /// `--after-the-flatten` under "Stage Options".
    #[derive(Parser, Debug)]
    struct Wrapper {
        #[command(flatten)]
        stage: MultiStageOptions,

        /// Declared after the flatten; belongs to the parent, not the stage.
        #[arg(long)]
        after_the_flatten: Option<u32>,
    }

    // Asserted on each argument's heading rather than on where it lands in the
    // rendered help: help ordering is clap's to change, and a parent argument
    // rendered before the heading may still carry it.
    let command = Wrapper::command();
    let heading_of = |long: &str| {
        command
            .get_arguments()
            .find(|arg| arg.get_long() == Some(long))
            .unwrap_or_else(|| panic!("--{long} should be registered"))
            .get_help_heading()
            .map(ToString::to_string)
    };

    assert_eq!(
        heading_of("stage::knob"),
        Some("Stage Options".to_string()),
        "the stage's own argument must carry the stage heading"
    );
    assert_eq!(
        heading_of("after-the-flatten"),
        None,
        "--after-the-flatten belongs to the parent and must carry no stage heading"
    );
}

// A parent that declares no doc comment of its own. The comment below is
// deliberately NOT a rustdoc comment: clap adopts a parent's own doc as its
// description, which would mask the very leak this test looks for.
//
// clap adopts a flattened `Args` struct's doc comment as the parent command's
// description whenever the parent has none. The generated companion carries
// rustdoc (crates that deny(missing_docs) require it), so the macro must reset
// `about` explicitly — otherwise an undocumented parent would describe itself
// with the companion's boilerplate, and with several stages flattened clap
// would arbitrarily pick whichever came first.
#[derive(Parser, Debug)]
struct UndocumentedParent {
    #[command(flatten)]
    opts: MultiStageOptions,
}

#[test]
fn the_companions_own_docs_never_become_the_parents_description() {
    // Asserted on the command's own metadata rather than on rendered help, so
    // the test cannot be satisfied by clap merely laying the description out
    // somewhere the substring check does not look.
    let command = UndocumentedParent::command();
    let about = command.get_about().map(ToString::to_string);
    let long_about = command.get_long_about().map(ToString::to_string);

    // The contract is that an undocumented parent stays undocumented, not
    // merely that the companion's *current* boilerplate wording is absent — a
    // substring check would start passing the moment that wording changed,
    // leak and all.
    for (label, description) in [("about", &about), ("long_about", &long_about)] {
        assert_eq!(
            *description, None,
            "an undocumented parent must keep an unset {label}, got {description:?}"
        );
    }
}

// ─────────────────────────────────────────────────────────────────────────────
// Field-level help metadata
// ─────────────────────────────────────────────────────────────────────────────

/// A field-level `help_heading` cannot be namespaced per stage. Emitted after
/// the companion's own heading it would win, filing the field outside its
/// stage — and two stages declaring the same heading would merge into one
/// section naming neither. `help` carries no such ambiguity: it is the field's
/// own documentation and reads identically on both commands.
#[multi_options("meta", "Meta Options")]
#[derive(Args, Debug, Clone, PartialEq)]
pub struct MetaOptions {
    /// Doc-derived help.
    #[arg(long, help_heading = "Escape Hatch")]
    pub own_heading: Option<u32>,

    /// Doc-derived help that clap must not prefer.
    #[arg(long, help = "explicit help wins")]
    pub own_help: Option<u32>,
}

#[derive(Parser, Debug)]
struct MetaWrapper {
    #[command(flatten)]
    opts: MultiMetaOptions,
}

#[test]
fn a_field_level_help_heading_cannot_escape_the_stages_heading() {
    let command = MetaWrapper::command();
    let arg = command
        .get_arguments()
        .find(|arg| arg.get_long() == Some("meta::own-heading"))
        .expect("--meta::own-heading should be registered");

    assert_eq!(
        arg.get_help_heading().map(ToString::to_string),
        Some("Meta Options".to_string()),
        "the stage's heading must win over the field's own"
    );
}

#[test]
fn an_explicit_help_is_preserved_verbatim() {
    let command = MetaWrapper::command();
    let arg = command
        .get_arguments()
        .find(|arg| arg.get_long() == Some("meta::own-help"))
        .expect("--meta::own-help should be registered");

    assert_eq!(
        arg.get_help().map(ToString::to_string),
        Some("explicit help wins".to_string()),
        "an explicit help is the field's own documentation and must survive re-exposure"
    );
}

// ─────────────────────────────────────────────────────────────────────────────
// Short/long help split survives re-exposure
// ─────────────────────────────────────────────────────────────────────────────

#[multi_options("docs", "Docs Options")]
#[derive(Args, Debug, Clone, PartialEq)]
pub struct DocsOptions {
    /// Short summary line.
    ///
    /// Longer explanation that clap shows only for `--help`, never for `-h`.
    /// Real fgumi options structs carry several such paragraphs per field.
    #[arg(long, default_value_t = 1)]
    pub documented: u32,
}

impl Default for DocsOptions {
    fn default() -> Self {
        Self { documented: 1 }
    }
}

#[derive(Parser, Debug)]
struct DocsWrapper {
    #[command(flatten)]
    opts: MultiDocsOptions,
}

/// Returns the `--docs::documented` argument of the prefixed companion.
///
/// The short/long split is asserted on `Arg::get_help()`/`get_long_help()`
/// rather than on `render_help()`: rendered help is clap's layout concern and
/// wraps to the terminal width once the `wrap_help` feature is anywhere in the
/// workspace feature graph, so a substring check there can fail while the
/// contract holds — and pass while it does not, since the text of one argument
/// is indistinguishable from any other's in the rendered block.
fn documented_arg(command: &clap::Command) -> &clap::Arg {
    command
        .get_arguments()
        .find(|arg| arg.get_long() == Some("docs::documented"))
        .expect("--docs::documented should be registered")
}

#[test]
fn short_help_shows_only_the_first_doc_paragraph() {
    let command = DocsWrapper::command();
    // clap trims the trailing period when it derives short help from a doc
    // comment, so match the sentence without it.
    let short = documented_arg(&command).get_help().map(ToString::to_string);
    let short = short.expect("the documented field must carry short help");
    assert!(short.contains("Short summary line"), "short help should show the summary: {short}");
    assert!(
        !short.contains("Longer explanation"),
        "short help must not carry the long explanation: {short}"
    );
}

#[test]
fn long_help_shows_every_doc_paragraph() {
    let command = DocsWrapper::command();
    let long = documented_arg(&command).get_long_help().map(ToString::to_string);
    let long = long.expect("the documented field must carry long help");
    assert!(long.contains("Short summary line."), "long help should show the summary: {long}");
    assert!(
        long.contains("Longer explanation"),
        "long help should show the full explanation: {long}"
    );
}

// ─────────────────────────────────────────────────────────────────────────────
// cfg-gated fields
// ─────────────────────────────────────────────────────────────────────────────

/// `#[cfg(any())]` is never satisfied and `#[cfg(all())]` always is, so this
/// fixture exercises both sides of cfg forwarding without needing a feature
/// flag. A macro that drops `#[cfg]` emits a Multi field (and `validate()` /
/// `From` arms) for a field that does not exist on the original struct, which
/// does not compile — so merely building this file is most of the test.
#[multi_options("gated", "Gated Options")]
#[derive(Args, Debug, Clone, PartialEq)]
pub struct GatedOptions {
    /// Always present.
    #[arg(long, default_value_t = 1)]
    pub always: u32,

    /// Compiled out.
    #[cfg(any())]
    #[arg(long, default_value_t = 2)]
    pub never: u32,
}

impl Default for GatedOptions {
    fn default() -> Self {
        Self { always: 1 }
    }
}

#[derive(Parser, Debug)]
struct GatedWrapper {
    #[command(flatten)]
    opts: MultiGatedOptions,
}

#[test]
fn cfg_gated_out_field_is_absent_from_the_multi_struct() {
    // Constructing the Multi struct with only the surviving field proves the
    // gated-out field was not emitted.
    let multi = MultiGatedOptions { gated_always: 5 };
    let opts = multi.validate().expect("validate");
    assert_eq!(opts.always, 5);

    assert!(
        GatedWrapper::command().get_arguments().all(|arg| arg.get_long() != Some("gated::never")),
        "the cfg-ed out field must not be registered as an argument"
    );
}

// ─────────────────────────────────────────────────────────────────────────────
// Aliases are namespaced; short flags are never propagated
// ─────────────────────────────────────────────────────────────────────────────

/// Mirrors `fgumi clip`'s `--ref` / `-r` shape: a long override, a short flag
/// and a long alias on one required field. `multi_aliased` additionally covers
/// the list spellings, which take a different code path to the single-valued
/// ones and would otherwise only be checked for "does not error".
#[multi_options("alias", "Alias Options")]
#[derive(Args, Debug, Clone, PartialEq)]
pub struct AliasOptions {
    /// Reference fasta, with a short flag and a long alias.
    #[arg(long = "reference", short = 'r', alias = "ref")]
    pub reference: PathBuf,

    /// Carries every long-alias spelling at once.
    #[arg(
        long = "multi-aliased",
        aliases = ["one", "two"],
        visible_alias = "vis",
        visible_aliases = ["vis-one", "vis-two"],
        short_alias = 's',
        default_value_t = 0
    )]
    pub multi_aliased: u32,
}

#[derive(Parser, Debug)]
struct AliasWrapper {
    #[command(flatten)]
    opts: MultiAliasOptions,
}

#[test]
fn long_alias_is_namespaced_under_the_prefix() {
    let parsed = AliasWrapper::try_parse_from(["test-prog", "--alias::ref", "/tmp/a.fa"])
        .expect("the prefixed alias should parse");
    let opts = parsed.opts.validate().expect("validate");
    assert_eq!(opts.reference, PathBuf::from("/tmp/a.fa"));
}

/// Assert a parse failed *because clap does not know the flag* — the contract
/// under test — rather than for any reason at all.
///
/// A bare `is_err()` would pass by accident: every companion field is defaulted,
/// absent-able, or staged as `Option<T>`, so clap never fails a parse for a
/// missing value today, and a later fixture change could satisfy the assertion
/// with an unrelated failure.
#[track_caller]
fn assert_unknown_argument(result: Result<AliasWrapper, clap::Error>, flag: &str) {
    let err = result.map(|_| ()).expect_err(&format!("{flag} must be rejected"));
    assert_eq!(
        err.kind(),
        clap::error::ErrorKind::UnknownArgument,
        "{flag} must be rejected as an unknown argument, got {:?}: {err}",
        err.kind()
    );
}

#[test]
fn unprefixed_alias_does_not_leak_onto_the_parent_command() {
    assert_unknown_argument(
        AliasWrapper::try_parse_from(["test-prog", "--ref", "/tmp/a.fa"]),
        "--ref",
    );
}

#[test]
fn short_flag_is_not_propagated() {
    assert_unknown_argument(AliasWrapper::try_parse_from(["test-prog", "-r", "/tmp/a.fa"]), "-r");
    assert_unknown_argument(
        AliasWrapper::try_parse_from(["test-prog", "--alias::reference", "/tmp/a.fa", "-s", "1"]),
        "-s",
    );
}

/// Every long-alias spelling — `aliases`, `visible_alias`, `visible_aliases` —
/// must reach clap re-prefixed, not merely survive classification.
#[rstest::rstest]
#[case::alias_list_first("--alias::one")]
#[case::alias_list_second("--alias::two")]
#[case::visible_alias("--alias::vis")]
#[case::visible_alias_list_first("--alias::vis-one")]
#[case::visible_alias_list_second("--alias::vis-two")]
fn every_long_alias_spelling_is_namespaced(#[case] flag: &str) {
    let parsed =
        AliasWrapper::try_parse_from(["test-prog", "--alias::reference", "/tmp/a.fa", flag, "7"])
            .unwrap_or_else(|e| panic!("{flag} should parse: {e}"));
    assert_eq!(parsed.opts.validate().expect("validate").multi_aliased, 7);
}

/// The un-prefixed spellings must not reach the parent command, where two
/// stages re-exposing the same options struct would collide on them.
#[rstest::rstest]
#[case::alias_list_first("--one")]
#[case::alias_list_second("--two")]
#[case::visible_alias("--vis")]
#[case::visible_alias_list_first("--vis-one")]
#[case::visible_alias_list_second("--vis-two")]
fn no_long_alias_spelling_leaks_unprefixed(#[case] flag: &str) {
    assert_unknown_argument(
        AliasWrapper::try_parse_from(["test-prog", "--alias::reference", "/tmp/a.fa", flag, "7"]),
        flag,
    );
}

// ─────────────────────────────────────────────────────────────────────────────
// Bare `bool` stays a valueless flag
// ─────────────────────────────────────────────────────────────────────────────

#[multi_options("flag", "Flag Options")]
#[derive(Args, Debug, Clone, PartialEq)]
pub struct FlagOptions {
    /// Bare bool — clap's `SetTrue` idiom, absent means `false`.
    #[arg(long)]
    pub enabled: bool,
}

#[derive(Parser, Debug)]
struct FlagWrapper {
    #[command(flatten)]
    opts: MultiFlagOptions,
}

#[test]
fn bare_bool_defaults_to_false_without_being_required() {
    let parsed = FlagWrapper::try_parse_from(["test-prog"]).expect("parse with no flags");
    let opts = parsed.opts.validate().expect("a bare bool must not be required");
    assert!(!opts.enabled);
}

#[test]
fn bare_bool_takes_no_value() {
    let parsed =
        FlagWrapper::try_parse_from(["test-prog", "--flag::enabled"]).expect("valueless flag");
    assert!(parsed.opts.validate().expect("validate").enabled);
}

// ─────────────────────────────────────────────────────────────────────────────
// `required` on absent-able types is enforced by validate()
// ─────────────────────────────────────────────────────────────────────────────

#[multi_options("req", "Req Options")]
#[derive(Args, Debug, Clone, PartialEq)]
pub struct RequiredOptions {
    /// Required despite being `Option<T>` on the standalone command.
    #[arg(long, required = true)]
    pub needed: Option<u32>,

    /// Required despite being `Vec<T>` on the standalone command.
    #[arg(long, required = true, action = clap::ArgAction::Append)]
    pub needed_many: Vec<u32>,
}

#[derive(Parser, Debug)]
struct RequiredWrapper {
    #[command(flatten)]
    opts: MultiRequiredOptions,
}

#[test]
fn required_option_field_is_enforced_by_validate() {
    // clap must not enforce it during parsing (staged validation owns that),
    // but validate() must still refuse the missing value.
    let parsed = RequiredWrapper::try_parse_from(["test-prog", "--req::needed-many", "1"])
        .expect("parse should succeed; required-ness is staged");
    let err = parsed.opts.validate().expect_err("validate must reject the missing Option field");
    let msg = format!("{err:#}");
    // `--req::needed` is a prefix of `--req::needed-many`, so match the whole
    // "<flag> is required" phrase: a bare substring check would pass even if
    // validate() named the wrong field.
    assert!(
        msg.contains("--req::needed is required"),
        "error should name the Option field, not --req::needed-many: {msg}"
    );
}

#[test]
fn required_vec_field_is_enforced_by_validate() {
    let parsed = RequiredWrapper::try_parse_from(["test-prog", "--req::needed", "1"])
        .expect("parse should succeed; required-ness is staged");
    let err = parsed.opts.validate().expect_err("validate must reject the empty Vec field");
    let msg = format!("{err:#}");
    assert!(
        msg.contains("--req::needed-many is required"),
        "error should name the Vec field: {msg}"
    );
}

#[test]
fn required_fields_validate_once_supplied() {
    let parsed = RequiredWrapper::try_parse_from([
        "test-prog",
        "--req::needed",
        "3",
        "--req::needed-many",
        "4",
    ])
    .expect("parse");
    let opts = parsed.opts.validate().expect("validate");
    assert_eq!(opts.needed, Some(3));
    assert_eq!(opts.needed_many, vec![4]);
}

// ─────────────────────────────────────────────────────────────────────────────
// Vec fields pass through; the struct `Default` is not consulted
// ─────────────────────────────────────────────────────────────────────────────

/// `Default` returns a non-empty Vec, but clap never consults a struct's
/// `Default` — the standalone command yields an empty Vec when the flag is
/// omitted, so the prefixed flag must too.
#[multi_options("vecs", "Vec Options")]
#[derive(Args, Debug, Clone, PartialEq)]
pub struct VecOptions {
    /// Repeatable values.
    #[arg(long, action = clap::ArgAction::Append)]
    pub values: Vec<u32>,
}

impl Default for VecOptions {
    fn default() -> Self {
        Self { values: vec![7, 8, 9] }
    }
}

#[test]
fn omitted_vec_matches_the_standalone_command_not_the_struct_default() {
    #[derive(Parser, Debug)]
    struct Wrapper {
        #[command(flatten)]
        opts: MultiVecOptions,
    }
    #[derive(Parser, Debug)]
    struct Standalone {
        #[command(flatten)]
        opts: VecOptions,
    }

    let standalone = Standalone::try_parse_from(["test-prog"]).expect("parse standalone");
    let prefixed = Wrapper::try_parse_from(["test-prog"]).expect("parse prefixed");

    assert!(standalone.opts.values.is_empty(), "standalone yields an empty Vec");
    assert_eq!(
        prefixed.opts.validate().expect("validate").values,
        standalone.opts.values,
        "the prefixed flag must agree with the standalone command"
    );
}

// ─────────────────────────────────────────────────────────────────────────────
// From + validate() is a lossless round trip
// ─────────────────────────────────────────────────────────────────────────────

/// Mirrors `GroupOptions`' shape: CLI fields plus `#[arg(skip)]` slots that the
/// command fills in itself.
#[multi_options("trip", "Trip Options")]
#[derive(Args, Debug, Clone, PartialEq)]
pub struct TripOptions {
    /// A normal flag.
    #[arg(long, default_value_t = 1)]
    pub knob: u32,

    /// Skipped slot with an explicit expression.
    #[arg(skip = 5u32)]
    pub skipped_with_expr: u32,

    /// Bare skipped slot.
    #[arg(skip)]
    pub bare_skipped: u32,
}

impl Default for TripOptions {
    fn default() -> Self {
        Self { knob: 1, skipped_with_expr: 5, bare_skipped: 11 }
    }
}

#[test]
fn round_trip_preserves_skip_field_values() {
    // Values that match neither the skip expression nor the struct Default, so
    // a re-derived value is distinguishable from a preserved one.
    let original = TripOptions { knob: 2, skipped_with_expr: 100, bare_skipped: 200 };
    let multi: MultiTripOptions = original.clone().into();
    let back = multi.validate().expect("validate");
    assert_eq!(back, original, "From + validate must preserve skip-field values");
}

#[test]
fn parsed_skip_fields_fall_back_to_their_declared_value() {
    #[derive(Parser, Debug)]
    struct Wrapper {
        #[command(flatten)]
        opts: MultiTripOptions,
    }

    let parsed = Wrapper::try_parse_from(["test-prog"]).expect("parse");
    let opts = parsed.opts.validate().expect("validate");
    assert_eq!(opts.skipped_with_expr, 5, "skip = expr supplies the parse-time value");
    assert_eq!(opts.bare_skipped, 11, "bare skip falls back to the struct Default");
}

#[test]
fn skip_fields_are_not_exposed_as_flags() {
    #[derive(Parser, Debug)]
    struct Wrapper {
        #[command(flatten)]
        opts: MultiTripOptions,
    }
    // Assert on the registered arguments, which is the actual contract; a help
    // substring check would break on any doc comment that happens to use the word.
    let command = Wrapper::command();
    let longs: Vec<&str> = command.get_arguments().filter_map(clap::Arg::get_long).collect();
    assert!(!longs.contains(&"trip::skipped-with-expr"), "skip fields must not become flags");
    assert!(!longs.contains(&"trip::bare-skipped"), "skip fields must not become flags");
    assert!(longs.contains(&"trip::knob"), "the non-skip field should still be registered");
}

#[test]
fn try_from_is_the_canonical_conversion() {
    let original = TripOptions { knob: 3, skipped_with_expr: 4, bare_skipped: 5 };
    let multi = MultiTripOptions::from(original.clone());
    let back = TripOptions::try_from(multi).expect("TryFrom should succeed");
    assert_eq!(back, original);
}
