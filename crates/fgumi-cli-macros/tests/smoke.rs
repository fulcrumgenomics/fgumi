//! Smoke test for the `multi_options` attribute macro.
//!
//! Exercises the three field-kind branches the macro generates code
//! for: `Option<T>`, defaulted (`default_value_t`), and required
//! (non-`Option`, no default). Verifies that:
//!
//!   * The original struct compiles unchanged with its bare flags.
//!   * The generated `Multi<Name>` struct exposes the same fields
//!     under `--<prefix>::<flag>` names.
//!   * `Multi<Name>::validate()` succeeds when required fields are
//!     supplied and fails (with a helpful error) when they are not.
//!   * The round-trip `From<<Name>>` impl preserves field values.

use clap::{Args, Parser};
use fgumi_cli_macros::multi_options;

/// Standalone options struct annotated with `multi_options`. The
/// generated `MultiFooOptions` is what the smoke test pokes at.
///
/// Note: `Default` is hand-rolled (not derived) so the
/// `default_value_t` the macro reads from `FooOptions::default()`
/// agrees with the `default_value_t` on the original `#[arg(...)]`
/// attribute. The `#[derive(Default)]` shortcut would give `u32 = 0`
/// for `defaulted_knob`, which would disagree with the standalone
/// CLI's `default_value_t = 7`.
#[multi_options("foo", "Foo Options")]
#[derive(Args, Debug, Clone, PartialEq)]
pub struct FooOptions {
    /// Optional knob (None when not passed).
    #[arg(long)]
    pub optional_knob: Option<u32>,

    /// Defaulted knob — keeps its `default_value_t` on both sides.
    #[arg(long, default_value_t = 7)]
    pub defaulted_knob: u32,

    /// Required knob — missing on the prefixed side fails `validate()`.
    #[arg(long)]
    pub required_knob: u32,
}

impl Default for FooOptions {
    fn default() -> Self {
        Self { optional_knob: None, defaulted_knob: 7, required_knob: 0 }
    }
}

#[test]
fn multi_struct_validate_round_trips_when_all_fields_supplied() {
    // Build a MultiFooOptions by hand (avoiding clap parsing so the
    // test doesn't depend on a full clap Command harness).
    let multi = MultiFooOptions {
        foo_optional_knob: Some(42),
        foo_defaulted_knob: 3,
        foo_required_knob: Some(100),
    };

    let opts = multi.validate().expect("validate should succeed with all fields set");
    assert_eq!(opts.optional_knob, Some(42));
    assert_eq!(opts.defaulted_knob, 3);
    assert_eq!(opts.required_knob, 100);
}

#[test]
fn multi_struct_validate_fails_when_required_missing() {
    let multi =
        MultiFooOptions { foo_optional_knob: None, foo_defaulted_knob: 7, foo_required_knob: None };

    let err = multi.validate().expect_err("validate should fail when required field is None");
    let msg = format!("{err:#}");
    assert!(
        msg.contains("--foo::required-knob"),
        "expected '--foo::required-knob' in error, got: {msg}"
    );
    assert!(msg.contains("required when foo is selected"), "got: {msg}");
}

#[test]
fn from_original_round_trips_field_values() {
    let original = FooOptions { optional_knob: Some(11), defaulted_knob: 22, required_knob: 33 };
    let multi: MultiFooOptions = original.clone().into();
    assert_eq!(multi.foo_optional_knob, Some(11));
    assert_eq!(multi.foo_defaulted_knob, 22);
    assert_eq!(multi.foo_required_knob, Some(33));

    let back = multi.validate().expect("round-trip validate should succeed");
    assert_eq!(back, original);
}

#[test]
fn multi_struct_parses_prefixed_flags_via_clap() {
    // Use a wrapper struct with `#[command(flatten)]` because
    // `Args`-derived structs can only be parsed via a parent
    // `Parser`-derived command.
    #[derive(Parser, Debug)]
    struct Wrapper {
        #[command(flatten)]
        foo_opts: MultiFooOptions,
    }

    let wrapper = Wrapper::try_parse_from([
        "test-prog",
        "--foo::required-knob",
        "5",
        "--foo::optional-knob",
        "9",
    ])
    .expect("parse should succeed");
    let opts = wrapper.foo_opts.validate().expect("validate after parse");
    assert_eq!(opts.required_knob, 5);
    assert_eq!(opts.optional_knob, Some(9));
    assert_eq!(opts.defaulted_knob, 7); // unchanged
}

/// Coverage for the `default_value` (string-form) and `Vec<T>`
/// (pass-through) handling extensions. Mirrors the structure of
/// `FooOptions` but exercises the two paths separately:
///
///   * `string_defaulted` uses `default_value = "..."` (the string
///     form parsed via `value_parser`); the macro must treat this as
///     "defaulted" not "required".
///   * `repeated_paths` is a `Vec<u32>` — the macro must treat Vec<T>
///     like Option<T> (pass-through, empty Vec default, no
///     `default_value_t` rewrite).
#[multi_options("bar", "Bar Options")]
#[derive(Args, Debug, Clone, PartialEq)]
pub struct BarOptions {
    /// Defaulted via the string form of `default_value`.
    #[arg(long, default_value = "42")]
    pub string_defaulted: u32,

    /// Repeatable Vec<T> field.
    #[arg(long, action = clap::ArgAction::Append)]
    pub repeated_paths: Vec<u32>,
}

impl Default for BarOptions {
    fn default() -> Self {
        Self { string_defaulted: 42, repeated_paths: Vec::new() }
    }
}

#[test]
fn macro_treats_default_value_string_form_as_defaulted_not_required() {
    // No string_defaulted on the cli → uses the Default.
    #[derive(Parser, Debug)]
    struct Wrapper {
        #[command(flatten)]
        bar_opts: MultiBarOptions,
    }
    let wrapper = Wrapper::try_parse_from(["test-prog"]).expect("parse with all defaults");
    let opts = wrapper.bar_opts.validate().expect("validate");
    assert_eq!(opts.string_defaulted, 42, "string-form default should apply");
    assert!(opts.repeated_paths.is_empty(), "Vec<T> default is empty");
}

#[test]
fn macro_treats_vec_t_as_repeating_pass_through() {
    #[derive(Parser, Debug)]
    struct Wrapper {
        #[command(flatten)]
        bar_opts: MultiBarOptions,
    }
    let wrapper = Wrapper::try_parse_from([
        "test-prog",
        "--bar::repeated-paths",
        "1",
        "--bar::repeated-paths",
        "2",
        "--bar::repeated-paths",
        "3",
    ])
    .expect("parse with three --repeated-paths");
    let opts = wrapper.bar_opts.validate().expect("validate");
    assert_eq!(opts.repeated_paths, vec![1, 2, 3]);
}

/// Coverage for four attribute-handling fixes:
///
///   * `skip_with_expr` uses `#[arg(skip = expr)]`. The macro must use
///     the provided expression verbatim in `validate()` rather than
///     hard-wiring `BazOptions::default().skip_with_expr`.
///   * `bare_skipped` uses a bare `#[arg(skip)]` (no `= expr`). The
///     macro must fall back to `BazOptions::default().bare_skipped`
///     (the struct default), NOT `u32::default()`.
///   * `conditionally_required` carries `#[arg(required = true)]`. The
///     macro must strip the `required = ...` name-value form so the
///     generated Multi field (an `Option<T>`) is not forced by clap
///     during parsing — required-ness is handled by `validate()`.
///   * `tmp_dirs` carries `#[arg(long = "tmp-dir")]`, a `long` override
///     that differs from the kebab of the field name. The macro must
///     honor the override so the Multi flag is `--baz::tmp-dir`, NOT
///     `--baz::tmp-dirs`.
#[multi_options("baz", "Baz Options")]
#[derive(Args, Debug, Clone, PartialEq)]
pub struct BazOptions {
    /// Skipped field with an explicit expression value. Not exposed on
    /// the CLI; `validate()` must use this expression, not `Default`.
    #[arg(skip = 99u32)]
    pub skip_with_expr: u32,

    /// Bare-skipped field. Not exposed on the CLI; `validate()` must
    /// pull its value from `BazOptions::default()`, not `u32::default()`.
    #[arg(skip)]
    pub bare_skipped: u32,

    /// Required (non-default) field that also carries `required = true`.
    /// The Multi side must wrap this as `Option<T>` and not propagate
    /// `required = true` to clap.
    #[arg(long, required = true)]
    pub conditionally_required: u32,

    /// Repeatable field with a `long` override (`--tmp-dir`, not
    /// `--tmp-dirs`). The Multi side must derive `--baz::tmp-dir`.
    #[arg(long = "tmp-dir", action = clap::ArgAction::Append)]
    pub tmp_dirs: Vec<u32>,
}

impl Default for BazOptions {
    fn default() -> Self {
        // `skip_with_expr`'s Default is deliberately non-zero (and different
        // from the `skip = 99u32` expression) so the test proves the macro
        // emits the skip expression, not this Default value. `bare_skipped`'s
        // Default is deliberately non-zero (and different from `u32::default()`
        // == 0) so the bare-skip test proves the macro pulls from this struct
        // Default rather than the field type's Default.
        Self {
            skip_with_expr: 1,
            bare_skipped: 77,
            conditionally_required: 0,
            tmp_dirs: Vec::new(),
        }
    }
}

#[test]
fn skip_with_expr_uses_provided_expression_not_default() {
    // The skip expression is `99`, while `BazOptions::default().skip_with_expr`
    // is `0`. validate() must yield `99`, proving the expression is honored.
    let multi = MultiBazOptions { baz_conditionally_required: Some(5), baz_tmp_dirs: Vec::new() };
    let opts = multi.validate().expect("validate should succeed");
    assert_eq!(opts.skip_with_expr, 99, "skip = 99 expression must be used verbatim");
    assert_eq!(opts.conditionally_required, 5);
}

#[test]
fn bare_skip_falls_back_to_struct_default_not_type_default() {
    // `bare_skipped` is `#[arg(skip)]` (no `= expr`), so the macro must
    // populate it from `BazOptions::default().bare_skipped` (== 77), NOT
    // `u32::default()` (== 0). It is not a CLI flag, so MultiBazOptions does
    // not carry it.
    let multi = MultiBazOptions { baz_conditionally_required: Some(1), baz_tmp_dirs: Vec::new() };
    let opts = multi.validate().expect("validate should succeed");
    assert_eq!(
        opts.bare_skipped, 77,
        "bare #[arg(skip)] must use the struct Default (77), not u32::default() (0)"
    );
}

#[test]
fn long_override_renames_multi_flag() {
    // `tmp_dirs` carries `#[arg(long = "tmp-dir")]`, so the Multi flag must be
    // `--baz::tmp-dir` (honoring the override), and `--baz::tmp-dirs` (the
    // kebab of the field name) must be rejected.
    #[derive(Parser, Debug)]
    struct Wrapper {
        #[command(flatten)]
        baz_opts: MultiBazOptions,
    }

    // The overridden flag name parses.
    let wrapper = Wrapper::try_parse_from([
        "test-prog",
        "--baz::conditionally-required",
        "1",
        "--baz::tmp-dir",
        "3",
        "--baz::tmp-dir",
        "4",
    ])
    .expect("--baz::tmp-dir should parse");
    let opts = wrapper.baz_opts.validate().expect("validate after parse");
    assert_eq!(opts.tmp_dirs, vec![3, 4]);

    // The field-name kebab (`--baz::tmp-dirs`) must NOT be a valid flag,
    // proving the macro honored the `long = "tmp-dir"` override.
    let err = Wrapper::try_parse_from([
        "test-prog",
        "--baz::conditionally-required",
        "1",
        "--baz::tmp-dirs",
        "3",
    ]);
    assert!(
        err.is_err(),
        "--baz::tmp-dirs must be rejected (override renamed it to --baz::tmp-dir)"
    );
}

#[test]
fn required_name_value_does_not_force_clap_requirement() {
    // `required = true` on the original field must be stripped: parsing
    // the Multi side WITHOUT the prefixed flag must succeed (clap must not
    // enforce it), and the conditional requirement is surfaced by validate().
    #[derive(Parser, Debug)]
    struct Wrapper {
        #[command(flatten)]
        baz_opts: MultiBazOptions,
    }

    // No --baz::conditionally-required supplied: parse must succeed because
    // the Multi field is Option<T> and `required = true` was stripped.
    let wrapper =
        Wrapper::try_parse_from(["test-prog"]).expect("parse should succeed without required flag");
    let err = wrapper
        .baz_opts
        .validate()
        .expect_err("validate should fail when conditionally-required is missing");
    let msg = format!("{err:#}");
    assert!(
        msg.contains("--baz::conditionally-required"),
        "expected '--baz::conditionally-required' in error, got: {msg}"
    );

    // Supplying it parses and validates.
    let wrapper = Wrapper::try_parse_from(["test-prog", "--baz::conditionally-required", "7"])
        .expect("parse should succeed when flag supplied");
    let opts = wrapper.baz_opts.validate().expect("validate after parse");
    assert_eq!(opts.conditionally_required, 7);
    assert_eq!(opts.skip_with_expr, 99);
}
