#![deny(unsafe_code)]

//! Proc-macro support for fgumi CLI option re-exposure.
//!
//! Provides [`multi_options`] — an attribute macro that pairs a standalone
//! command's `clap::Args` options struct with a sibling `Multi<Name>` struct
//! whose CLI flags are prefixed and grouped under a help heading. Used by the
//! `runall` command to re-expose every per-stage option of `fgumi sort` /
//! `fgumi group` / `fgumi simplex` / `fgumi duplex` / `fgumi codec` without
//! hand-maintaining a parallel option set on `RunAll`.
//!
//! # Flag naming
//!
//! A re-exposed flag is named `--<prefix>::<flag>` (e.g. `--sort::max-memory`).
//! The `::` separator is deliberate: `runall` flattens several stages into one
//! command, and every stage independently owns flags like `--threads` or
//! `--max-memory`. A plain kebab prefix (`--sort-max-memory`) would be
//! ambiguous with a stage that genuinely has a `--max` flag taking a `memory`
//! value, and would read as one flag name rather than as a stage-qualified one.
//! `::` cannot appear in a clap-derived flag name, so it can never collide with
//! a real flag, and it makes the stage qualifier obvious in `--help`.
//!
//! # Field-kind handling
//!
//! * `#[arg(skip)]` / `#[arg(skip = expr)]` fields are invisible to the CLI.
//!   They are carried onto the Multi struct as skip fields too, so their values
//!   survive a `From` + `validate()` round trip; when the Multi struct is
//!   parsed from the command line they take the skip expression's value, or the
//!   original struct's `Default` for a bare `skip`.
//! * `Option<T>`, `Vec<T>` and bare `bool` fields keep their type: clap already
//!   treats each as absent-able (`None`, empty, `false`).
//! * Fields with any `default_value*` attribute keep their type and their
//!   default — the attribute is copied verbatim, so the prefixed flag advertises
//!   and applies exactly the default the standalone command does.
//! * Every other field is required: it becomes `Option<T>` on the Multi struct,
//!   and the generated `validate()` returns an error when the corresponding
//!   `--<prefix>::<flag>` is missing. This keeps required-ness staged — clap
//!   never refuses to parse, so `runall` can report which stage is missing what.
//!
//! # What is and is not carried over
//!
//! Copied verbatim onto the generated field: `#[doc]` comments (each attribute
//! separately, so clap's short/long help split survives), `#[cfg]` and
//! `#[cfg_attr]` gates, every `default_value*`, and every other `#[arg(...)]`
//! key the macro does not classify (`value_parser`, `value_name`, `action`,
//! `num_args`, `hide`, `env`, …). A `#[cfg]` additionally gates the generated
//! conversion arms, so a field compiled out is absent from all three sites;
//! `#[cfg_attr]` never removes a field, so it lands on the field alone.
//!
//! Rewritten: `long` (re-prefixed) and long-form aliases (`alias`,
//! `aliases`, `visible_alias`, `visible_aliases` — each re-prefixed the same
//! way, so no un-namespaced flag ever reaches the parent command).
//!
//! Dropped: `short` and every short alias — one letter cannot be namespaced per
//! stage. A field-level `help_heading` is dropped for the same reason: the
//! companion files every field under its stage's heading, and a field-level one
//! is emitted after it and wins, so the field would escape its stage — and two
//! stages declaring the same heading would merge into one section naming
//! neither. `required` is also dropped from the clap side and enforced by
//! `validate()` instead, so a missing value is reported by `runall`'s staged
//! validation rather than by clap's parser.
//!
//! `help` and `long_help` are preserved verbatim, like any other unclassified
//! key: they are the field's own documentation and read identically on both
//! commands. One consequence is that an explicit `help` on a required field
//! replaces the generated "Required when `<stage>` is selected." line, since
//! clap prefers `help` over `#[doc]`; the staged `validate()` error still names
//! both the flag and the stage.
//!
//! Rejected at expansion time, with a spanned build error: cross-field
//! reference keys (`requires`, `conflicts_with`, …) whose arg ids would dangle
//! once fields are prefixed; `id` / `name` overrides; positional arguments,
//! whether declared by `index` or by declaring neither `long` nor `short`
//! (several stages' positionals would be mutually ambiguous, and clap panics
//! when a positional carries a `long`); clap's `key(value)` call form for any
//! key the macro classifies; a clap attribute hidden behind `#[cfg_attr]`,
//! which the macro cannot classify through; a `required` that could never be
//! enforced, because the field always holds a value (defaulted — including a
//! defaulted `Option<T>` or `Vec<T>` — skipped, or a bare `bool`);
//! `#[command(...)]` on a field; struct-level `#[command(...)]` /
//! `#[group(...)]`; generic structs, since the companion and its conversions
//! are emitted without generic parameters; and the legacy `#[clap(...)]` /
//! `#[structopt(...)]` spellings, which every classifier would otherwise
//! silently ignore.
//!
//! Note that `env` is preserved verbatim: two stages re-exposing the same
//! options struct therefore read the same environment variable, exactly as the
//! two standalone commands would.

use proc_macro::TokenStream;
use proc_macro2::TokenStream as TokenStream2;
use quote::{format_ident, quote};
use syn::{Lit, Meta};

/// Attribute macro that passes the original options struct through unchanged and
/// generates a `Multi<OriginalName>` companion whose flags are named
/// `--<prefix>::<flag>` and filed under `<heading>` in `--help`, plus a
/// `validate()` method, a `TryFrom<Multi<OriginalName>>` impl and a
/// `From<OriginalName>` impl.
///
/// Usage: `#[multi_options("prefix", "Help Heading")]`
///
/// The generated items inherit the annotated struct's visibility.
///
/// See the crate-level docs for the field-kind classification rules and for the
/// full list of what is carried over, rewritten, dropped and rejected.
///
/// # Errors
///
/// Emits a spanned compile error — pointing at the offending field, attribute or
/// literal — rather than panicking. See the crate docs' "What is and is not
/// carried over" section for the rejected forms.
#[proc_macro_attribute]
pub fn multi_options(attr: TokenStream, item: TokenStream) -> TokenStream {
    let item = TokenStream2::from(item);
    match expand(TokenStream2::from(attr), &item) {
        Ok(expanded) => expanded.into(),
        Err(error) => {
            // Re-emit the annotated item next to the diagnostic. Without it, every
            // downstream reference to the struct raises its own "cannot find type"
            // error and buries the one that actually explains the problem.
            let compile_error = error.into_compile_error();
            quote! { #item #compile_error }.into()
        }
    }
}

/// Expand `#[multi_options(...)]`, returning either the generated tokens or the
/// accumulated diagnostics.
fn expand(attr: TokenStream2, item: &TokenStream2) -> syn::Result<TokenStream2> {
    let args = syn::parse2::<MultiOptionsArgs>(attr)?;
    let input = parse_annotated_struct(item)?;
    let fields = named_fields(&input)?;

    let context = ExpandContext {
        struct_name: &input.ident,
        vis: &input.vis,
        prefix: &args.prefix,
        heading: &args.heading,
    };
    let generated = generate_fields(fields, &context)?;
    Ok(render(&input, &generated, &context))
}

/// Parse the annotated item, rejecting anything that is not a struct.
fn parse_annotated_struct(item: &TokenStream2) -> syn::Result<syn::ItemStruct> {
    let parsed = syn::parse2::<syn::Item>(item.clone())?;
    let syn::Item::Struct(input) = parsed else {
        return Err(syn::Error::new_spanned(
            parsed,
            "multi_options only supports structs with named fields",
        ));
    };
    reject_unsupported_struct_attrs(&input.attrs)?;
    // The companion struct and both conversion impls are emitted without generic
    // parameters, so a generic options struct expands into code that cannot
    // compile — and every resulting error names the type parameter rather than
    // this macro.
    if !input.generics.params.is_empty() {
        return Err(syn::Error::new_spanned(&input.generics, GENERICS_MSG));
    }
    // `Generics`' `ToTokens` emits only the `<...>` params, so spanning on it
    // when they are empty would collapse the span to the call site and lose the
    // field-accurate diagnostics every other rejection here produces.
    if let Some(where_clause) = &input.generics.where_clause {
        return Err(syn::Error::new_spanned(where_clause, GENERICS_MSG));
    }
    Ok(input)
}

/// Borrow the struct's named fields, rejecting tuple and unit structs.
fn named_fields(
    input: &syn::ItemStruct,
) -> syn::Result<&syn::punctuated::Punctuated<syn::Field, syn::Token![,]>> {
    match &input.fields {
        syn::Fields::Named(named) => Ok(&named.named),
        other => Err(syn::Error::new_spanned(
            other,
            "multi_options only supports structs with named fields",
        )),
    }
}

/// Classify and generate every field, reporting all bad fields at once rather
/// than making the author fix them one build at a time.
fn generate_fields(
    fields: &syn::punctuated::Punctuated<syn::Field, syn::Token![,]>,
    context: &ExpandContext<'_>,
) -> syn::Result<Vec<GeneratedField>> {
    let mut errors: Option<syn::Error> = None;
    let mut generated = Vec::with_capacity(fields.len());
    for field in fields {
        match ParsedField::from_field(field, context.prefix) {
            Ok(parsed) => generated.push(parsed.generate(context)),
            Err(error) => match &mut errors {
                Some(accumulated) => accumulated.combine(error),
                None => errors = Some(error),
            },
        }
    }
    match errors {
        Some(error) => Err(error),
        None => Ok(generated),
    }
}

/// Assemble the original struct, the generated companion, and the conversions.
fn render(
    input: &syn::ItemStruct,
    generated: &[GeneratedField],
    context: &ExpandContext<'_>,
) -> TokenStream2 {
    let struct_name = context.struct_name;
    let multi_struct_name = format_ident!("Multi{}", struct_name);
    let multi_fields = generated.iter().map(|g| &g.multi_field);
    let validate_arms = generated.iter().map(|g| &g.validate_arm);
    let from_arms = generated.iter().map(|g| &g.from_arm);
    let vis = context.vis;

    quote! {
        #input

        /// Prefixed options struct generated by `#[multi_options]` for the
        /// runall command. Carries the same fields as the original options
        /// struct but exposes each via `--<prefix>::<flag>`, filed under the
        /// stage's help heading.
        #[derive(::clap::Args, Debug, Clone)]
        // clap adopts a flattened struct's doc comment as the parent command's
        // description when the parent declares none. This struct's rustdoc is
        // for docs.rs (and for crates that deny(missing_docs)) — it would be
        // nonsense as a command description, and with several stages flattened
        // clap would arbitrarily pick whichever came first.
        #[command(about = None, long_about = None)]
        #vis struct #multi_struct_name {
            #(#multi_fields)*
        }

        impl #multi_struct_name {
            /// Validate required fields and convert to the original options
            /// struct. Returns `Err` naming the missing `--<prefix>::<flag>`
            /// when a field the standalone command requires was not supplied.
            #vis fn validate(self) -> ::anyhow::Result<#struct_name> {
                <#struct_name as ::core::convert::TryFrom<Self>>::try_from(self)
            }
        }

        impl ::core::convert::TryFrom<#multi_struct_name> for #struct_name {
            type Error = ::anyhow::Error;

            fn try_from(opts: #multi_struct_name) -> ::anyhow::Result<Self> {
                Ok(Self {
                    #(#validate_arms)*
                })
            }
        }

        impl ::core::convert::From<#struct_name> for #multi_struct_name {
            fn from(opts: #struct_name) -> Self {
                Self {
                    #(#from_arms)*
                }
            }
        }
    }
}

/// Everything the per-field generator needs from the annotated struct.
struct ExpandContext<'a> {
    struct_name: &'a syn::Ident,
    vis: &'a syn::Visibility,
    prefix: &'a str,
    heading: &'a str,
}

// ─────────────────────────────────────────────────────────────────────────────
// Macro arguments
// ─────────────────────────────────────────────────────────────────────────────

/// Parsed arguments for `#[multi_options("prefix", "heading")]`.
struct MultiOptionsArgs {
    prefix: String,
    heading: String,
}

/// The one-line usage reminder attached to every argument-shape diagnostic.
const USAGE: &str = "multi_options requires two string literal arguments: \
                     #[multi_options(\"prefix\", \"heading\")]";

/// Rejection message for a generic annotated struct, shared by the type-parameter
/// and `where`-clause arms so the two cannot drift apart.
const GENERICS_MSG: &str = "multi_options does not support generic structs: the generated Multi struct and its \
     conversions are emitted without generic parameters, so the expansion would not compile. \
     Use a concrete options struct.";

impl syn::parse::Parse for MultiOptionsArgs {
    fn parse(input: syn::parse::ParseStream) -> syn::Result<Self> {
        let prefix_lit: syn::LitStr =
            input.parse().map_err(|e| syn::Error::new(e.span(), USAGE))?;
        input.parse::<syn::Token![,]>().map_err(|e| syn::Error::new(e.span(), USAGE))?;
        let heading_lit: syn::LitStr =
            input.parse().map_err(|e| syn::Error::new(e.span(), USAGE))?;
        if !input.is_empty() {
            return Err(syn::Error::new(input.span(), USAGE));
        }

        let prefix = prefix_lit.value();
        validate_prefix(&prefix).map_err(|msg| syn::Error::new(prefix_lit.span(), msg))?;
        let heading = heading_lit.value();
        if heading.is_empty() {
            return Err(syn::Error::new(
                heading_lit.span(),
                "multi_options: the help heading must not be empty",
            ));
        }

        Ok(Self { prefix, heading })
    }
}

/// Check that a prefix can be spliced into both a flag name and an identifier.
///
/// The prefix appears in the flag (`--<prefix>::<flag>`) and in the generated
/// field identifier (`<prefix>_<field>`), so an empty or non-identifier prefix
/// would otherwise surface as `--::flag` or as an opaque `format_ident!` panic
/// with no mention of `multi_options`.
fn validate_prefix(prefix: &str) -> Result<(), String> {
    let Some(first) = prefix.chars().next() else {
        return Err("multi_options: the prefix must not be empty (it would generate flags named \
             `--::<flag>`)"
            .to_string());
    };
    if !first.is_ascii_alphabetic() {
        return Err(format!(
            "multi_options: the prefix must start with an ASCII letter, got `{prefix}` — it also \
             becomes the leading segment of the generated field identifier `<prefix>_<field>`"
        ));
    }
    if let Some(bad) =
        prefix.chars().find(|c| !(c.is_ascii_alphanumeric() || *c == '_' || *c == '-'))
    {
        return Err(format!(
            "multi_options: the prefix may only contain ASCII letters, digits, `_` and `-`, but \
             `{prefix}` contains `{bad}`"
        ));
    }
    Ok(())
}

/// Reject struct-level clap configuration the macro cannot faithfully reproduce.
///
/// The generated struct is built from the fields alone, so any struct-level clap
/// setting would apply to the standalone command and silently not to the
/// re-exposed one.
fn reject_unsupported_struct_attrs(attrs: &[syn::Attribute]) -> syn::Result<()> {
    for attr in attrs {
        let path = attr.path();
        if path.is_ident("clap") || path.is_ident("structopt") {
            return Err(legacy_spelling_error(attr));
        }
        if path.is_ident("group") {
            return Err(syn::Error::new_spanned(
                attr,
                "multi_options does not support a struct-level #[group(...)]: the group names its \
                 members by their unprefixed arg ids, which do not exist on the generated Multi \
                 struct. Enforce the grouping in the command's validate()/resolve() instead.",
            ));
        }
        if path.is_ident("command") {
            return Err(syn::Error::new_spanned(
                attr,
                "multi_options does not carry a struct-level #[command(...)] onto the generated \
                 Multi struct, so the standalone and re-exposed commands would silently diverge. \
                 Move the setting onto the individual #[arg(...)] attributes, or drop it.",
            ));
        }
    }
    Ok(())
}

/// Reject a `#[cfg_attr(...)]` that hides a clap attribute behind a condition.
///
/// Every classifier keys on a literal `#[arg(...)]` / `#[command(...)]`, so a
/// `#[cfg_attr(unix, arg(long, default_value_t = 3))]` is invisible: the field is
/// classified as required, wrapped in `Option<T>`, and *then* the forwarded
/// attribute applies `default_value_t` to the wrapped type — a wall of type
/// errors that never mentions `multi_options`. Classifying through the condition
/// is not possible either, since the macro cannot evaluate `cfg` predicates.
fn reject_conditional_clap_attr(attr: &syn::Attribute, field: &syn::Field) -> syn::Result<()> {
    let Ok(metas) =
        attr.parse_args_with(syn::punctuated::Punctuated::<Meta, syn::Token![,]>::parse_terminated)
    else {
        // Not a shape we can inspect; leave it to rustc.
        return Ok(());
    };
    // The first meta is the `cfg` predicate; the rest are the attributes it gates.
    for meta in metas.iter().skip(1) {
        let path = meta.path();
        if path.is_ident("arg") || path.is_ident("clap") || path.is_ident("command") {
            let key = path.get_ident().map_or_else(|| "arg".to_string(), ToString::to_string);
            return Err(syn::Error::new_spanned(
                attr,
                format!(
                    "multi_options: field `{}` hides a clap attribute behind #[cfg_attr(…, \
                     {key}(…))]. The macro classifies fields from their literal #[arg(...)] \
                     attributes and cannot evaluate a cfg predicate, so this one would be ignored \
                     and the field misclassified. Apply #[cfg] to the field and write the \
                     #[arg(...)] unconditionally.",
                    field_name(field)
                ),
            ));
        }
    }
    Ok(())
}

/// Build the diagnostic for clap's legacy attribute spellings.
///
/// Every classifier in this macro keys on `#[arg(...)]` / `#[command(...)]`, so a
/// `#[clap(skip)]` would be invisible — the field would be exposed as a required
/// CLI flag instead of being skipped. Reject rather than silently misclassify.
fn legacy_spelling_error(attr: &syn::Attribute) -> syn::Error {
    syn::Error::new_spanned(
        attr,
        "multi_options does not support the legacy #[clap(...)] / #[structopt(...)] spelling: \
         every classifier keys on #[arg(...)] and #[command(...)], so this attribute would be \
         silently ignored and the field misclassified. Use the #[arg(...)] / #[command(...)] \
         spelling.",
    )
}

// ─────────────────────────────────────────────────────────────────────────────
// clap `#[arg(...)]` key tables
// ─────────────────────────────────────────────────────────────────────────────

/// clap `#[arg(...)]` keys that reference *another argument by its string id*.
///
/// The Multi struct renames every field to `<prefix>_<field>` and never emits a
/// matching arg alias, so any of these ids would dangle on the Multi side and
/// clap would panic ("arg id `x` not defined") when it builds the `runall`
/// command. Reject them at macro-expansion time instead (D1).
const CROSS_REFERENCE_ARG_KEYS: &[&str] = &[
    "requires",
    "requires_all",
    "requires_if",
    "requires_ifs",
    "conflicts_with",
    "conflicts_with_all",
    "overrides_with",
    "overrides_with_all",
    "required_if_eq",
    "required_if_eq_all",
    "required_if_eq_any",
    "required_unless_present",
    "required_unless_present_any",
    "required_unless_present_all",
    "default_value_if",
    "default_value_ifs",
    "default_values_if",
    "default_values_ifs",
    "group",
    "groups",
];

/// clap `#[arg(...)]` keys that rename the argument itself.
///
/// The Multi struct derives every arg id from its own prefixed field name, so an
/// explicit id would reintroduce the unprefixed name and collide across stages.
const RENAMING_ARG_KEYS: &[&str] = &["id", "name"];

/// clap `#[arg(...)]` keys that declare the argument to be positional.
///
/// A positional argument has no flag name to prefix, and clap panics outright
/// when one is given a `long` — which is exactly what the Multi struct emits.
const POSITIONAL_ARG_KEYS: &[&str] = &["index"];

/// clap `#[arg(...)]` keys that declare a default value, in every spelling.
///
/// A field carrying any of these is optional on the standalone command, so it
/// must not be classified as required. Deliberately excludes
/// `default_missing_value`, which is the value used when a flag is passed
/// *without* one and says nothing about whether the flag may be omitted.
const DEFAULT_VALUE_ARG_KEYS: &[&str] = &[
    "default_value",
    "default_value_t",
    "default_value_os",
    "default_value_os_t",
    "default_values",
    "default_values_t",
    "default_values_os",
    "default_values_os_t",
];

/// clap `#[arg(...)]` keys this macro classifies via `Meta::NameValue` /
/// `Meta::Path` (to strip, rewrite, or read them).
///
/// clap also accepts the equivalent `key(value)` call form, which arrives as a
/// `Meta::List` and would slip past every classifier — silently changing the
/// flag name, the required-ness, the default, or leaking an un-prefixed alias.
/// Reject the call form for these keys (D2).
const CALL_FORM_SENSITIVE_ARG_KEYS: &[&str] = &[
    "long",
    "short",
    "required",
    "skip",
    "alias",
    "aliases",
    "visible_alias",
    "visible_aliases",
    "short_alias",
    "short_aliases",
    "visible_short_alias",
    "visible_short_aliases",
    "help_heading",
    "default_value",
    "default_value_t",
    "default_value_os",
    "default_value_os_t",
    "default_values",
    "default_values_t",
    "default_values_os",
    "default_values_os_t",
];

/// Long-form alias keys, which are re-prefixed exactly like `long`.
const LONG_ALIAS_ARG_KEYS: &[&str] = &["alias", "visible_alias"];

/// Long-form alias keys taking a list of aliases; every element is re-prefixed.
const LONG_ALIAS_LIST_ARG_KEYS: &[&str] = &["aliases", "visible_aliases"];

/// Short-form alias keys, dropped for the same reason as `short`.
const SHORT_ALIAS_ARG_KEYS: &[&str] =
    &["short_alias", "short_aliases", "visible_short_alias", "visible_short_aliases"];

// ─────────────────────────────────────────────────────────────────────────────
// Per-field parsing
// ─────────────────────────────────────────────────────────────────────────────

/// How a `#[arg(skip)]` field obtains its value.
enum Skip {
    /// Bare `#[arg(skip)]` — falls back to the original struct's `Default`.
    Bare,
    /// `#[arg(skip = expr)]` — uses the expression verbatim.
    Expr(syn::Expr),
}

/// One field of the annotated struct, classified in a single pass over its
/// attributes.
///
/// Every `#[arg(...)]` meta is visited exactly once and routed to the piece of
/// state it affects, so the key tables above are the single source of truth for
/// classification — earlier revisions re-walked `field.attrs` in six independent
/// helpers, which let the three key lists drift apart.
struct ParsedField<'a> {
    ident: &'a syn::Ident,
    ty: &'a syn::Type,
    /// `#[doc]`, `#[cfg]` and `#[cfg_attr]` attributes, copied onto the
    /// generated field verbatim.
    forwarded: Vec<&'a syn::Attribute>,
    /// The `#[cfg]` subset, which must additionally gate the `TryFrom` and
    /// `From` arms so a gated-out field is absent from all three sites.
    ///
    /// `#[cfg_attr]` is deliberately excluded: it never removes a field, so the
    /// arms do not need it, and only `#[cfg]` is valid on a struct-expression
    /// field — forwarding `#[cfg_attr(unix, serde(skip))]` onto one would expand
    /// to `#[serde(skip)]` on an expression, which does not compile.
    cfgs: Vec<&'a syn::Attribute>,
    /// `#[arg(...)]` entries carried onto the generated field, already rendered
    /// as `, key = value` continuations ready to splice into a new `#[arg(...)]`.
    preserved: Vec<TokenStream2>,
    skip: Option<Skip>,
    /// The value of an explicit `#[arg(long = "...")]` override.
    long_override: Option<String>,
    has_default: bool,
    /// Whether the field carries `#[arg(required)]` / `#[arg(required = true)]`.
    required: bool,
    /// Whether the field declares a flag name — `long` or `short`, in any
    /// spelling. A field declaring neither is positional to clap.
    has_flag_name: bool,
}

impl<'a> ParsedField<'a> {
    /// Classify one field, or return every diagnostic it earns.
    fn from_field(field: &'a syn::Field, prefix: &str) -> syn::Result<Self> {
        let ident = field.ident.as_ref().ok_or_else(|| {
            syn::Error::new_spanned(field, "multi_options only supports structs with named fields")
        })?;
        let mut parsed = Self {
            ident,
            ty: &field.ty,
            forwarded: Vec::new(),
            cfgs: Vec::new(),
            preserved: Vec::new(),
            skip: None,
            long_override: None,
            has_default: false,
            required: false,
            has_flag_name: false,
        };

        for attr in &field.attrs {
            let path = attr.path();
            if path.is_ident("doc") {
                parsed.forwarded.push(attr);
            } else if path.is_ident("cfg") {
                parsed.forwarded.push(attr);
                parsed.cfgs.push(attr);
            } else if path.is_ident("cfg_attr") {
                reject_conditional_clap_attr(attr, field)?;
                parsed.forwarded.push(attr);
            } else if path.is_ident("clap") || path.is_ident("structopt") {
                return Err(legacy_spelling_error(attr));
            } else if path.is_ident("command") {
                return Err(field_command_error(attr, field));
            } else if path.is_ident("arg") {
                for meta in parse_arg_metas(attr, field)? {
                    parsed.absorb_arg_meta(&meta, field, prefix)?;
                }
            }
        }

        parsed.check_required_is_enforceable(field)?;
        parsed.check_is_not_positional(field)?;
        Ok(parsed)
    }

    /// Reject a field clap would treat as a positional argument.
    ///
    /// In clap's derive a field with neither `long` nor `short` is positional.
    /// The Multi struct always emits a `long`, which would silently convert
    /// `fgumi sort <input>` into `--sort::input <input>` — and, when the field
    /// carries an explicit `index`, makes clap panic outright while building the
    /// runall command. A positional cannot be namespaced per stage anyway: with
    /// several stages flattened into one command, their positionals would be
    /// mutually ambiguous.
    fn check_is_not_positional(&self, field: &syn::Field) -> syn::Result<()> {
        if self.skip.is_some() || self.has_flag_name {
            return Ok(());
        }
        Err(syn::Error::new_spanned(
            field,
            format!(
                "multi_options: field `{}` declares neither `long` nor `short`, so clap treats it \
                 as a positional argument. A positional has no flag name to prefix, and several \
                 stages flattened into one runall command would have mutually ambiguous \
                 positionals. Add an explicit `#[arg(long)]` (or `#[arg(long = \"...\")]`) to the \
                 field.",
                field_name(field)
            ),
        ))
    }

    /// Route one `#[arg(...)]` meta to the state it affects.
    fn absorb_arg_meta(
        &mut self,
        meta: &Meta,
        field: &syn::Field,
        prefix: &str,
    ) -> syn::Result<()> {
        let Some(key) = meta.path().get_ident().map(ToString::to_string) else {
            self.preserved.push(quote! { , #meta });
            return Ok(());
        };
        let key = key.as_str();
        let name = field_name(field);
        reject_unsupported_arg_key(meta, field, &name, key)?;

        match key {
            "skip" => {
                self.skip = Some(match meta {
                    Meta::NameValue(nv) => Skip::Expr(nv.value.clone()),
                    _ => Skip::Bare,
                });
            }
            // The Multi field declares its own prefixed `long`; a bare `long`
            // just means "kebab of the field name", which is the fallback.
            "long" => {
                self.has_flag_name = true;
                if let Meta::NameValue(nv) = meta {
                    self.long_override = Some(string_literal(&nv.value).ok_or_else(|| {
                        non_literal_error(field, &name, "long", "a string literal")
                    })?);
                }
            }
            // One letter cannot be namespaced per stage.
            "short" => self.has_flag_name = true,
            _ if SHORT_ALIAS_ARG_KEYS.contains(&key) => {}
            // Neither can a help heading. The companion files every field under
            // the stage's heading; a field-level one is emitted after it and
            // wins, so the field escapes its stage — and two stages declaring
            // the same heading would merge into one section that names neither.
            "help_heading" => {}
            // Required-ness is enforced by validate(), not by clap, so runall
            // can name the stage that is missing a value.
            "required" => {
                self.required = match meta {
                    Meta::Path(_) => true,
                    Meta::NameValue(nv) => bool_literal(&nv.value).ok_or_else(|| {
                        non_literal_error(field, &name, "required", "`true` or `false`")
                    })?,
                    Meta::List(_) => unreachable!("call form rejected above"),
                };
            }
            _ if LONG_ALIAS_ARG_KEYS.contains(&key) => {
                self.preserved.push(prefixed_alias(meta, field, &name, key, prefix)?);
            }
            _ if LONG_ALIAS_LIST_ARG_KEYS.contains(&key) => {
                self.preserved.push(prefixed_alias_list(meta, field, &name, key, prefix)?);
            }
            _ => {
                if DEFAULT_VALUE_ARG_KEYS.contains(&key) {
                    self.has_default = true;
                }
                self.preserved.push(quote! { , #meta });
            }
        }
        Ok(())
    }

    /// Reject `#[arg(required)]` on a field whose generated form always holds a
    /// value.
    ///
    /// `validate()` enforces required-ness by observing absence — a `None`, an
    /// empty `Vec`, or the `Option<T>` wrapper the macro adds. A skipped field, a
    /// defaulted field and a bare `bool` are never absent, so the requirement
    /// would be silently unenforceable.
    ///
    /// `Option<T>` and `Vec<T>` are absent-able only while they carry no default:
    /// clap applies a declared `default_value` to those types too (an omitted
    /// `Option<u32>` with `default_value = "42"` parses as `Some(42)`, not
    /// `None`), so the generated `is_none()` / `is_empty()` check would never
    /// fire and the requirement would be lost in silence.
    fn check_required_is_enforceable(&self, field: &syn::Field) -> syn::Result<()> {
        if !self.required {
            return Ok(());
        }
        let name = field_name(field);
        if self.skip.is_some() {
            return Err(syn::Error::new_spanned(
                field,
                format!(
                    "multi_options: field `{name}` combines #[arg(required …)] with #[arg(skip)]. \
                     A skipped field is never supplied on the command line, so the requirement \
                     could never be satisfied."
                ),
            ));
        }
        let absent_able = (is_option_type(self.ty) || is_vec_type(self.ty)) && !self.has_default;
        if self.is_required_kind() || absent_able {
            return Ok(());
        }
        Err(syn::Error::new_spanned(
            field,
            format!(
                "multi_options: field `{name}` combines #[arg(required …)] with a default value or \
                 a bare `bool`. The generated field always holds a value — clap applies a declared \
                 default to `Option<T>` and `Vec<T>` as well — so the Multi side cannot distinguish \
                 \"not supplied\" from \"supplied the default\" and the requirement would be \
                 silently unenforceable. Drop `required`, or drop the default."
            ),
        ))
    }

    /// Whether the field must be wrapped in `Option<T>` on the Multi struct and
    /// enforced by `validate()`.
    ///
    /// `Option<T>`, `Vec<T>` and bare `bool` are all absent-able as-is, and a
    /// field with a default is never missing.
    fn is_required_kind(&self) -> bool {
        !self.has_default
            && !is_option_type(self.ty)
            && !is_vec_type(self.ty)
            && !is_bool_type(self.ty)
    }

    /// Emit the Multi-struct field plus the matching `TryFrom` and `From` arms.
    fn generate(&self, ctx: &ExpandContext<'_>) -> GeneratedField {
        let field_ident = self.ident;
        let field_type = self.ty;
        let prefixed_ident = format_ident!("{}_{}", ctx.prefix.replace('-', "_"), field_ident);
        let forwarded = &self.forwarded;
        let cfgs = &self.cfgs;
        let preserved = &self.preserved;
        let vis = ctx.vis;

        // Skipped fields are invisible to the CLI on both sides. Carrying them
        // onto the Multi struct as skip fields keeps `From` + `validate()` a
        // genuine round trip instead of re-deriving the value from `Default`.
        if let Some(skip) = &self.skip {
            let value = match skip {
                Skip::Expr(expr) => quote! { #expr },
                Skip::Bare => {
                    let struct_name = ctx.struct_name;
                    quote! { #struct_name::default().#field_ident }
                }
            };
            return GeneratedField {
                multi_field: quote! {
                    #(#forwarded)*
                    #[arg(skip = #value)]
                    #vis #prefixed_ident: #field_type,
                },
                validate_arm: quote! { #(#cfgs)* #field_ident: opts.#prefixed_ident, },
                from_arm: quote! { #(#cfgs)* #prefixed_ident: opts.#field_ident, },
            };
        }

        // Honor an explicit `#[arg(long = "...")]` so the prefixed flag tracks the
        // standalone command's flag name (`tmp_dirs` with `long = "tmp-dir"`
        // becomes `--<prefix>::tmp-dir`, not `--<prefix>::tmp-dirs`).
        let base_name =
            self.long_override.clone().unwrap_or_else(|| field_ident.to_string().replace('_', "-"));
        let long_name = format!("{}::{}", ctx.prefix, base_name);
        let heading = ctx.heading;
        let prefix = ctx.prefix;
        let error_msg = format!("--{long_name} is required when {prefix} is selected");

        if self.is_required_kind() {
            let required_doc = format!("Required when {prefix} is selected.");
            return GeneratedField {
                multi_field: quote! {
                    #(#forwarded)*
                    #[doc = #required_doc]
                    #[arg(long = #long_name, help_heading = #heading #(#preserved)*)]
                    #vis #prefixed_ident: Option<#field_type>,
                },
                validate_arm: quote! {
                    #(#cfgs)*
                    #field_ident: opts.#prefixed_ident
                        .ok_or_else(|| ::anyhow::anyhow!(#error_msg))?,
                },
                from_arm: quote! { #(#cfgs)* #prefixed_ident: Some(opts.#field_ident), },
            };
        }

        // Absent-able or defaulted: the field keeps its type and its declared
        // default. `required` on such a field is dropped from the clap side, so
        // re-assert it here — otherwise a flag the standalone command demands
        // would be quietly optional on the Multi side.
        let value = if self.required && is_option_type(field_type) {
            quote! {{
                let value = opts.#prefixed_ident;
                if value.is_none() {
                    ::anyhow::bail!(#error_msg);
                }
                value
            }}
        } else if self.required && is_vec_type(field_type) {
            quote! {{
                let value = opts.#prefixed_ident;
                if value.is_empty() {
                    ::anyhow::bail!(#error_msg);
                }
                value
            }}
        } else {
            quote! { opts.#prefixed_ident }
        };

        GeneratedField {
            multi_field: quote! {
                #(#forwarded)*
                #[arg(long = #long_name, help_heading = #heading #(#preserved)*)]
                #vis #prefixed_ident: #field_type,
            },
            validate_arm: quote! { #(#cfgs)* #field_ident: #value, },
            from_arm: quote! { #(#cfgs)* #prefixed_ident: opts.#field_ident, },
        }
    }
}

/// The three token streams one source field contributes to the expansion.
struct GeneratedField {
    multi_field: TokenStream2,
    validate_arm: TokenStream2,
    from_arm: TokenStream2,
}

// ─────────────────────────────────────────────────────────────────────────────
// Attribute helpers
// ─────────────────────────────────────────────────────────────────────────────

/// Render a field's name for a diagnostic.
///
/// `multi_options` only accepts named-field structs, so the `None` arm is
/// unreachable in practice; it exists so a diagnostic never panics while
/// reporting another diagnostic.
fn field_name(field: &syn::Field) -> String {
    field.ident.as_ref().map_or_else(|| "<unnamed>".to_string(), ToString::to_string)
}

/// Parse the comma-separated metas out of one `#[arg(...)]` attribute.
///
/// An unparseable `#[arg(...)]` is a hard build error: silently treating it as
/// "no attributes" would misclassify the field — a dropped `skip` exposes a
/// field that should have no CLI flag, a dropped `long` renames it, a dropped
/// `default_value` makes it required.
fn parse_arg_metas(
    attr: &syn::Attribute,
    field: &syn::Field,
) -> syn::Result<syn::punctuated::Punctuated<Meta, syn::Token![,]>> {
    attr.parse_args_with(syn::punctuated::Punctuated::<Meta, syn::Token![,]>::parse_terminated)
        .map_err(|e| {
            syn::Error::new_spanned(
                attr,
                format!(
                    "multi_options: failed to parse #[arg(...)] on field `{}`: {e}",
                    field_name(field)
                ),
            )
        })
}

/// Build the diagnostic for a field-level `#[command(...)]`.
///
/// `flatten` and `subcommand` nest another struct whose fields the macro cannot
/// reach to prefix; any other field-level `#[command(...)]` key would simply not
/// be carried over.
fn field_command_error(attr: &syn::Attribute, field: &syn::Field) -> syn::Error {
    let name = field_name(field);
    let metas = attr
        .parse_args_with(syn::punctuated::Punctuated::<Meta, syn::Token![,]>::parse_terminated)
        .ok();
    let nests = metas.is_some_and(|metas| {
        metas
            .iter()
            .any(|meta| meta.path().is_ident("flatten") || meta.path().is_ident("subcommand"))
    });
    if nests {
        syn::Error::new_spanned(
            attr,
            format!(
                "multi_options does not support #[command(flatten)] / #[command(subcommand)] on \
                 field `{name}`. The nested struct's fields cannot be reached to prefix them; \
                 inline them directly."
            ),
        )
    } else {
        syn::Error::new_spanned(
            attr,
            format!(
                "multi_options does not carry a field-level #[command(...)] onto the generated \
                 Multi struct, so field `{name}` would behave differently on the two commands. \
                 Move the setting onto #[arg(...)], or drop it."
            ),
        )
    }
}

/// Reject an `#[arg(...)]` key the macro cannot faithfully re-expose.
///
/// Covers three latent traps: cross-field reference keys whose arg-id strings
/// dangle once the field is prefixed (D1), explicit id overrides that would
/// reintroduce the un-prefixed name, and clap's `key(value)` call form for any
/// key the macro only recognizes in `key = value` / bare-`key` form (D2).
fn reject_unsupported_arg_key(
    meta: &Meta,
    field: &syn::Field,
    name: &str,
    key: &str,
) -> syn::Result<()> {
    if CROSS_REFERENCE_ARG_KEYS.contains(&key) {
        return Err(syn::Error::new_spanned(
            field,
            format!(
                "multi_options: field `{name}` uses #[arg({key} …)], which references another \
                 argument by its unprefixed id. The Multi struct renames fields to \
                 `<prefix>_<field>`, so that id would dangle and clap would panic when it builds \
                 the runall command. Enforce this coupling in the command's validate()/resolve() \
                 instead (see AlignerOptions::resolve)."
            ),
        ));
    }
    if RENAMING_ARG_KEYS.contains(&key) {
        return Err(syn::Error::new_spanned(
            field,
            format!(
                "multi_options: field `{name}` uses #[arg({key} …)], which overrides the \
                 argument's id. The Multi struct derives every id from its own prefixed field \
                 name, so an explicit id would reintroduce the un-prefixed name and collide \
                 between stages. Drop it."
            ),
        ));
    }
    if POSITIONAL_ARG_KEYS.contains(&key) {
        return Err(syn::Error::new_spanned(
            field,
            format!(
                "multi_options: field `{name}` uses #[arg({key} …)], which makes it a positional \
                 argument. The Multi struct gives every field a prefixed `long`, and clap panics \
                 when a positional has one (\"is a positional argument and can't have short or \
                 long name versions\"). A positional cannot be namespaced per stage — expose it \
                 as a flag with `#[arg(long)]` instead."
            ),
        ));
    }
    if matches!(meta, Meta::List(_)) && CALL_FORM_SENSITIVE_ARG_KEYS.contains(&key) {
        return Err(syn::Error::new_spanned(
            field,
            format!(
                "multi_options: field `{name}` uses the call form #[arg({key}(…))]. Use the \
                 `{key} = …` name-value form (or bare `{key}`) — the macro only classifies those \
                 spellings and would mishandle the call form."
            ),
        ));
    }
    Ok(())
}

/// Re-prefix a single long alias so it can never reach the parent command
/// un-namespaced.
fn prefixed_alias(
    meta: &Meta,
    field: &syn::Field,
    name: &str,
    key: &str,
    prefix: &str,
) -> syn::Result<TokenStream2> {
    let Meta::NameValue(nv) = meta else {
        return Err(non_literal_error(field, name, key, "a string literal"));
    };
    let alias = string_literal(&nv.value)
        .ok_or_else(|| non_literal_error(field, name, key, "a string literal"))?;
    let prefixed = format!("{prefix}::{alias}");
    let key = format_ident!("{}", key);
    Ok(quote! { , #key = #prefixed })
}

/// Re-prefix every element of an alias list (`aliases = ["a", "b"]`).
fn prefixed_alias_list(
    meta: &Meta,
    field: &syn::Field,
    name: &str,
    key: &str,
    prefix: &str,
) -> syn::Result<TokenStream2> {
    const EXPECTED: &str = "an array of string literals";
    let Meta::NameValue(nv) = meta else {
        return Err(non_literal_error(field, name, key, EXPECTED));
    };
    let syn::Expr::Array(array) = &nv.value else {
        return Err(non_literal_error(field, name, key, EXPECTED));
    };
    let prefixed = array
        .elems
        .iter()
        .map(|elem| {
            string_literal(elem)
                .map(|alias| format!("{prefix}::{alias}"))
                .ok_or_else(|| non_literal_error(field, name, key, EXPECTED))
        })
        .collect::<syn::Result<Vec<_>>>()?;
    let key = format_ident!("{}", key);
    Ok(quote! { , #key = [#(#prefixed),*] })
}

/// Build the diagnostic for a classified key whose value is not the literal form
/// the macro can read.
fn non_literal_error(field: &syn::Field, name: &str, key: &str, expected: &str) -> syn::Error {
    syn::Error::new_spanned(
        field,
        format!(
            "multi_options: field `{name}` uses #[arg({key} = …)] with a value that is not \
             {expected}. The macro reads this key at expansion time — to rewrite the flag name, \
             re-prefix the alias, or classify required-ness — so it cannot be an arbitrary \
             expression."
        ),
    )
}

/// Read a string literal out of an attribute value.
fn string_literal(expr: &syn::Expr) -> Option<String> {
    if let syn::Expr::Lit(expr_lit) = expr
        && let Lit::Str(lit) = &expr_lit.lit
    {
        return Some(lit.value());
    }
    None
}

/// Read a boolean literal out of an attribute value.
fn bool_literal(expr: &syn::Expr) -> Option<bool> {
    if let syn::Expr::Lit(expr_lit) = expr
        && let Lit::Bool(lit) = &expr_lit.lit
    {
        return Some(lit.value);
    }
    None
}

// ─────────────────────────────────────────────────────────────────────────────
// Type predicates
// ─────────────────────────────────────────────────────────────────────────────

/// Check whether a type is `Vec<T>`.
///
/// clap collects a `Vec<T>` naturally (empty when the flag is absent), so it
/// needs no default to be absent-able.
fn is_vec_type(ty: &syn::Type) -> bool {
    if let syn::Type::Path(type_path) = ty
        && let Some(segment) = type_path.path.segments.last()
    {
        return segment.ident == "Vec";
    }
    false
}

/// Check whether a type is a bare `bool`.
///
/// clap gives a bare `bool` field `ArgAction::SetTrue` — a valueless flag
/// defaulting to `false` — so it needs no `default_value` to be absent-able.
fn is_bool_type(ty: &syn::Type) -> bool {
    if let syn::Type::Path(type_path) = ty
        && let Some(segment) = type_path.path.segments.last()
    {
        return segment.ident == "bool" && segment.arguments.is_none();
    }
    false
}

/// Check whether a type is `Option<T>`.
fn is_option_type(ty: &syn::Type) -> bool {
    if let syn::Type::Path(type_path) = ty
        && let Some(segment) = type_path.path.segments.last()
    {
        return segment.ident == "Option";
    }
    false
}

#[cfg(test)]
mod tests {
    use super::*;
    use rstest::rstest;
    use syn::parse::Parser;

    /// Parse a single named struct field from tokens (e.g.
    /// `#[arg(long, requires = "x")] pub y: u32`).
    fn named_field(tokens: TokenStream2) -> syn::Field {
        syn::Field::parse_named.parse2(tokens).expect("parse named field")
    }

    /// Classify a field the way `expand` does, returning the diagnostic message
    /// on rejection.
    fn classify(tokens: TokenStream2) -> Result<(), String> {
        let field = named_field(tokens);
        ParsedField::from_field(&field, "p").map(|_| ()).map_err(|e| e.to_string())
    }

    /// Every `#[arg(...)]` spelling the macro classifies must be accepted in its
    /// bare and name-value forms.
    #[rstest]
    #[case::short_and_string_default(quote! { #[arg(long, short = 'x', default_value = "7")] pub a: u32 })]
    #[case::long_override_and_typed_default(
        quote! { #[arg(long = "max-memory", default_value_t = 5)] pub b: usize }
    )]
    #[case::bare_long_on_option(quote! { #[arg(long)] pub c: Option<u32> })]
    #[case::vec_with_required(
        quote! { #[arg(long, value_delimiter = ',', required = true)] pub d: Vec<usize> }
    )]
    #[case::short_only(quote! { #[arg(short = 'x')] pub e: u32 })]
    #[case::bare_skip(quote! { #[arg(skip)] pub f: u32 })]
    #[case::skip_with_expr(quote! { #[arg(skip = 7u32)] pub g: u32 })]
    #[case::long_alias(quote! { #[arg(long, alias = "ref")] pub h: u32 })]
    #[case::alias_list(quote! { #[arg(long, aliases = ["ref", "fasta"])] pub i: u32 })]
    #[case::short_alias(quote! { #[arg(long, short_alias = 'r')] pub j: u32 })]
    #[case::os_default(quote! { #[arg(long, default_value_os = "/tmp/x")] pub k: PathBuf })]
    #[case::doc_and_cfg(quote! { /// docs
        #[cfg(unix)]
        #[arg(long)] pub l: Option<u32> })]
    fn accepts_supported_arg_forms(#[case] tokens: TokenStream2) {
        assert_eq!(classify(tokens), Ok(()), "expected supported form to be accepted");
    }

    /// D1: every key that names another argument by id would dangle once the
    /// field is prefixed. The macro classifies on the key alone, so the uniform
    /// `key = "other"` spelling exercises each one.
    #[rstest]
    fn rejects_every_cross_reference_key(
        #[values(
            "requires",
            "requires_all",
            "requires_if",
            "requires_ifs",
            "conflicts_with",
            "conflicts_with_all",
            "overrides_with",
            "overrides_with_all",
            "required_if_eq",
            "required_if_eq_all",
            "required_if_eq_any",
            "required_unless_present",
            "required_unless_present_any",
            "required_unless_present_all",
            "default_value_if",
            "default_value_ifs",
            "default_values_if",
            "default_values_ifs",
            "group",
            "groups"
        )]
        key: &str,
    ) {
        let key_ident = format_ident!("{}", key);
        let err = classify(quote! { #[arg(long, #key_ident = "other")] pub a: u32 })
            .expect_err("cross-reference key must be rejected");
        assert!(err.contains(key), "message should name the key {key:?}: {err}");
        assert!(err.contains("field `a`"), "message should name the field: {err}");
    }

    /// Exhaustive by construction: drives the constant itself, so a key added to
    /// `CROSS_REFERENCE_ARG_KEYS` is covered without touching the table above.
    #[test]
    fn every_cross_reference_key_in_the_table_is_rejected() {
        for key in CROSS_REFERENCE_ARG_KEYS {
            let key_ident = format_ident!("{}", key);
            let err = classify(quote! { #[arg(long, #key_ident = "other")] pub a: u32 })
                .expect_err("every cross-reference key must be rejected");
            assert!(err.contains(key), "message should name the key {key:?}: {err}");
        }
    }

    /// Every key the classifier reads or rewrites must also be rejected in
    /// clap's `key(value)` call form — that form arrives as a `Meta::List` and
    /// matches no classifier arm, so a key missing from
    /// `CALL_FORM_SENSITIVE_ARG_KEYS` is silently preserved verbatim and the
    /// rewrite is skipped. Adding a spelling to any table below without adding
    /// it here reopens the D2 hole for that key.
    #[test]
    fn every_classified_key_is_call_form_sensitive() {
        let classified = ["skip", "long", "short", "required", "help_heading"]
            .iter()
            .copied()
            .chain(DEFAULT_VALUE_ARG_KEYS.iter().copied())
            .chain(LONG_ALIAS_ARG_KEYS.iter().copied())
            .chain(LONG_ALIAS_LIST_ARG_KEYS.iter().copied())
            .chain(SHORT_ALIAS_ARG_KEYS.iter().copied());
        for key in classified {
            assert!(
                CALL_FORM_SENSITIVE_ARG_KEYS.contains(&key),
                "classified key {key:?} is missing from CALL_FORM_SENSITIVE_ARG_KEYS, so its \
                 `{key}(…)` call form would slip past the classifier"
            );
        }
    }

    /// Every `default_value*` spelling must mark the field as defaulted; a
    /// missed one is silently reclassified as required and gains a spurious
    /// `Option<T>` wrapper plus a false "Required when …" help line.
    #[test]
    fn every_default_value_spelling_makes_a_field_non_required() {
        for key in DEFAULT_VALUE_ARG_KEYS {
            let key_ident = format_ident!("{}", key);
            let field = named_field(quote! { #[arg(long, #key_ident = "x")] pub a: u32 });
            let parsed = ParsedField::from_field(&field, "p")
                .unwrap_or_else(|e| panic!("{key} should classify: {e}"));
            assert!(parsed.has_default, "{key} must count as a default");
            assert!(!parsed.is_required_kind(), "{key} must not leave the field required");
        }
    }

    /// A positional argument has no flag name to prefix, and clap panics when one
    /// carries a `long` — which the companion always emits.
    #[rstest]
    #[case::no_attributes(quote! { pub a: u32 }, "positional")]
    #[case::arg_without_long_or_short(quote! { #[arg(value_name = "FILE")] pub b: u32 }, "positional")]
    #[case::explicit_index(quote! { #[arg(index = 1)] pub c: u32 }, "index")]
    fn rejects_positional_fields(#[case] tokens: TokenStream2, #[case] needle: &str) {
        let err = classify(tokens).expect_err("a positional field must be rejected");
        assert!(err.contains(needle), "message should mention {needle:?}: {err}");
    }

    /// A skipped field is never a CLI argument, so it needs no flag name.
    #[test]
    fn skipped_fields_are_exempt_from_the_positional_check() {
        assert_eq!(classify(quote! { #[arg(skip)] pub a: u32 }), Ok(()));
    }

    /// The macro cannot evaluate a `cfg` predicate, so a clap attribute behind
    /// one would be ignored and the field misclassified.
    #[rstest]
    #[case::arg(quote! { #[cfg_attr(unix, arg(long, default_value_t = 3))] pub a: u32 })]
    #[case::command(quote! { #[cfg_attr(unix, command(flatten))] pub b: u32 })]
    #[case::legacy_clap(quote! { #[cfg_attr(unix, clap(long))] pub c: u32 })]
    fn rejects_clap_attributes_hidden_behind_cfg_attr(#[case] tokens: TokenStream2) {
        let err = classify(tokens).expect_err("a conditional clap attribute must be rejected");
        assert!(err.contains("cfg_attr"), "message should name the attribute: {err}");
    }

    /// A `#[cfg_attr]` carrying no clap attribute is forwarded untouched.
    #[test]
    fn allows_cfg_attr_without_a_clap_attribute() {
        assert_eq!(
            classify(quote! { #[cfg_attr(unix, doc = "unix only")] #[arg(long)] pub a: u32 }),
            Ok(())
        );
    }

    /// D2: the `key(value)` call form arrives as `Meta::List` and would slip past
    /// the `Meta::NameValue` / `Meta::Path` classifiers.
    #[rstest]
    #[case::long(quote! { #[arg(long("x"))] pub a: u32 }, "long")]
    #[case::short(quote! { #[arg(long, short('x'))] pub b: u32 }, "short")]
    #[case::default_value_t(quote! { #[arg(default_value_t(4))] pub c: u32 }, "default_value_t")]
    #[case::required(quote! { #[arg(long, required(true))] pub d: u32 }, "required")]
    #[case::alias(quote! { #[arg(long, alias("ref"))] pub e: u32 }, "alias")]
    #[case::visible_alias(
        quote! { #[arg(long, visible_alias("ref"))] pub f: u32 },
        "visible_alias"
    )]
    #[case::default_value_os(
        quote! { #[arg(long, default_value_os("/tmp/x"))] pub g: PathBuf },
        "default_value_os"
    )]
    fn rejects_call_form_for_classified_keys(#[case] tokens: TokenStream2, #[case] needle: &str) {
        let err = classify(tokens).expect_err("call form must be rejected");
        assert!(err.contains(needle), "message should name the key {needle:?}: {err}");
    }

    #[test]
    fn allows_call_form_for_unclassified_keys() {
        // `value_parser(...)` in call form is preserved verbatim and is not one of
        // the classified keys, so it must not be rejected.
        assert_eq!(
            classify(quote! { #[arg(long, value_parser(clap::value_parser!(u32)))] pub a: u32 }),
            Ok(())
        );
    }

    /// Keys the macro reads at expansion time cannot be arbitrary expressions.
    #[rstest]
    #[case::long(quote! { #[arg(long = SOME_CONST)] pub a: u32 }, "long")]
    #[case::alias(quote! { #[arg(long, alias = SOME_CONST)] pub b: u32 }, "alias")]
    #[case::aliases(quote! { #[arg(long, aliases = SOME_CONST)] pub c: u32 }, "aliases")]
    #[case::aliases_of_non_literals(
        quote! { #[arg(long, aliases = [SOME_CONST])] pub d: u32 },
        "aliases"
    )]
    #[case::required(quote! { #[arg(long, required = SOME_CONST)] pub e: u32 }, "required")]
    fn rejects_non_literal_values_for_classified_keys(
        #[case] tokens: TokenStream2,
        #[case] needle: &str,
    ) {
        let err = classify(tokens).expect_err("non-literal value must be rejected");
        assert!(err.contains(needle), "message should name the key {needle:?}: {err}");
    }

    /// An explicit arg id would reintroduce the un-prefixed name.
    #[rstest]
    #[case::id(quote! { #[arg(long, id = "other")] pub a: u32 }, "id")]
    #[case::name(quote! { #[arg(long, name = "other")] pub b: u32 }, "name")]
    fn rejects_renaming_keys(#[case] tokens: TokenStream2, #[case] needle: &str) {
        let err = classify(tokens).expect_err("renaming key must be rejected");
        assert!(err.contains(needle), "message should name the key {needle:?}: {err}");
    }

    /// The legacy spellings bypass every classifier, so they must not be ignored.
    #[rstest]
    #[case::clap_skip(quote! { #[clap(skip)] pub a: u32 })]
    #[case::clap_long(quote! { #[clap(long, default_value = "3")] pub b: u32 })]
    #[case::structopt(quote! { #[structopt(long)] pub c: u32 })]
    fn rejects_legacy_attribute_spellings(#[case] tokens: TokenStream2) {
        let err = classify(tokens).expect_err("legacy spelling must be rejected");
        assert!(err.contains("#[clap(...)]"), "message should name the spelling: {err}");
    }

    /// A field-level `#[command(...)]` either nests a struct the macro cannot
    /// reach, or would silently not be carried over.
    #[rstest]
    #[case::flatten(quote! { #[command(flatten)] pub a: Inner }, "flatten")]
    #[case::subcommand(quote! { #[command(subcommand)] pub b: Inner }, "subcommand")]
    #[case::other(quote! { #[command(next_help_heading = "x")] pub c: u32 }, "field-level")]
    fn rejects_field_level_command_attrs(#[case] tokens: TokenStream2, #[case] needle: &str) {
        let err = classify(tokens).expect_err("field-level #[command(...)] must be rejected");
        assert!(err.contains(needle), "message should mention {needle:?}: {err}");
    }

    /// `#[command(flatten)]` is matched on the parsed meta, not on a substring of
    /// the attribute's tokens, so a key that merely *contains* "flatten" is
    /// reported as an ordinary unsupported field-level attribute.
    #[test]
    fn flatten_guard_does_not_substring_match() {
        let err = classify(quote! { #[command(help_heading = "flatten me")] pub a: u32 })
            .expect_err("field-level #[command(...)] is rejected");
        assert!(
            err.contains("field-level"),
            "a value merely containing \"flatten\" must not be reported as a flatten: {err}"
        );
    }

    /// `required` is only enforceable where absence is observable.
    ///
    /// The `Option<T>` / `Vec<T>` cases are the subtle ones: clap applies a
    /// declared default even to those types, so the generated `is_none()` /
    /// `is_empty()` check never fires and the requirement is silently lost.
    #[rstest]
    #[case::defaulted(quote! { #[arg(long, default_value_t = 3, required = true)] pub a: u32 })]
    #[case::bare_bool(quote! { #[arg(long, required = true)] pub b: bool })]
    #[case::skipped(quote! { #[arg(skip, required = true)] pub c: u32 })]
    #[case::defaulted_option(
        quote! { #[arg(long, default_value = "3", required = true)] pub d: Option<u32> }
    )]
    #[case::defaulted_vec(
        quote! { #[arg(long, default_values_t = [1u32], required = true)] pub e: Vec<u32> }
    )]
    fn rejects_unenforceable_required(#[case] tokens: TokenStream2) {
        let err = classify(tokens).expect_err("unenforceable `required` must be rejected");
        assert!(err.contains("required"), "message should name the key: {err}");
    }

    /// Without a default, absence IS observable on both types, so `required` is
    /// enforceable and must be accepted.
    #[rstest]
    #[case::option(quote! { #[arg(long, required = true)] pub a: Option<u32> })]
    #[case::vec(quote! { #[arg(long, required = true)] pub b: Vec<u32> })]
    fn accepts_enforceable_required_on_absent_able_types(#[case] tokens: TokenStream2) {
        assert_eq!(classify(tokens), Ok(()));
    }

    /// A generic struct expands into a companion emitted without generics, so the
    /// expansion cannot compile — and the resulting errors name `T`, never
    /// `multi_options`.
    #[rstest]
    #[case::type_param(quote! { pub struct Opts<T> { #[arg(long)] pub a: T } })]
    #[case::lifetime(quote! { pub struct Opts<'a> { #[arg(long)] pub a: &'a str } })]
    #[case::where_clause(
        quote! { pub struct Opts where u32: Clone { #[arg(long)] pub a: u32 } }
    )]
    fn rejects_generic_structs(#[case] tokens: TokenStream2) {
        let err = parse_annotated_struct(&tokens)
            .expect_err("a generic struct must be rejected")
            .to_string();
        assert!(err.contains("generic"), "message should mention generics: {err}");
    }

    /// The ordinary non-generic case must still pass the same guard.
    #[test]
    fn accepts_a_plain_struct() {
        assert!(parse_annotated_struct(&quote! { pub struct Opts { pub a: u32 } }).is_ok());
    }

    /// `required = false` is a no-op and must not trip the enforceability check.
    #[test]
    fn accepts_required_false_on_a_defaulted_field() {
        assert_eq!(
            classify(quote! { #[arg(long, default_value_t = 3, required = false)] pub a: u32 }),
            Ok(())
        );
    }

    /// `required` has a bare form as well as the name-value one, and only the
    /// name-value `false` disables it.
    #[rstest]
    #[case::bare(quote! { #[arg(long, required)] pub a: Option<u32> }, true)]
    #[case::explicit_true(quote! { #[arg(long, required = true)] pub a: Option<u32> }, true)]
    #[case::explicit_false(quote! { #[arg(long, required = false)] pub a: Option<u32> }, false)]
    #[case::absent(quote! { #[arg(long)] pub a: Option<u32> }, false)]
    fn recognizes_every_required_spelling(#[case] tokens: TokenStream2, #[case] expected: bool) {
        let field = named_field(tokens);
        let parsed = ParsedField::from_field(&field, "p").expect("classify");
        assert_eq!(parsed.required, expected);
    }

    /// A bare alias key carries no value to re-prefix, so it cannot be honored.
    #[rstest]
    #[case::bare_alias(quote! { #[arg(long, alias)] pub a: u32 }, "alias")]
    #[case::bare_aliases(quote! { #[arg(long, aliases)] pub b: u32 }, "aliases")]
    fn rejects_valueless_alias_keys(#[case] tokens: TokenStream2, #[case] needle: &str) {
        let err = classify(tokens).expect_err("a valueless alias key must be rejected");
        assert!(err.contains(needle), "message should name the key {needle:?}: {err}");
    }

    /// Struct-level clap configuration is not reproduced on the companion, so it
    /// must be rejected rather than silently applied to one command only.
    #[rstest]
    #[case::legacy_clap(quote! { #[clap(next_help_heading = "x")] }, "#[clap(...)]")]
    #[case::legacy_structopt(quote! { #[structopt(name = "x")] }, "#[clap(...)]")]
    #[case::group(quote! { #[group(required = true)] }, "#[group(...)]")]
    #[case::command(quote! { #[command(next_help_heading = "x")] }, "#[command(...)]")]
    fn rejects_unsupported_struct_level_attrs(
        #[case] attr_tokens: TokenStream2,
        #[case] needle: &str,
    ) {
        let item: syn::ItemStruct = syn::parse2(quote! {
            #attr_tokens
            pub struct Opts { pub a: u32 }
        })
        .expect("parse struct");
        let err = reject_unsupported_struct_attrs(&item.attrs)
            .expect_err("struct-level attr must be rejected")
            .to_string();
        assert!(err.contains(needle), "message should name {needle:?}: {err}");
    }

    /// `#[derive(...)]` and doc comments ride along on every annotated struct and
    /// must not be mistaken for clap configuration.
    #[test]
    fn accepts_derive_and_doc_attrs_on_the_struct() {
        let item: syn::ItemStruct = syn::parse2(quote! {
            /// Docs.
            #[derive(clap::Args, Debug, Clone)]
            pub struct Opts { pub a: u32 }
        })
        .expect("parse struct");
        assert!(reject_unsupported_struct_attrs(&item.attrs).is_ok());
    }

    /// The attribute takes exactly two non-empty string literals.
    #[rstest]
    #[case::well_formed(quote! { "sort", "Sort Options" }, true)]
    #[case::missing_heading(quote! { "sort" }, false)]
    #[case::missing_comma(quote! { "sort" "Sort Options" }, false)]
    #[case::extra_argument(quote! { "sort", "Sort Options", "extra" }, false)]
    #[case::empty(quote! {}, false)]
    #[case::empty_heading(quote! { "sort", "" }, false)]
    #[case::non_literal(quote! { sort, "Sort Options" }, false)]
    fn parses_only_two_string_literal_arguments(
        #[case] tokens: TokenStream2,
        #[case] expected_ok: bool,
    ) {
        assert_eq!(syn::parse2::<MultiOptionsArgs>(tokens).is_ok(), expected_ok);
    }

    /// The prefix is spliced into both a flag name and an identifier.
    #[rstest]
    #[case::simple("sort", true)]
    #[case::with_digit("codec2", true)]
    #[case::with_underscore("read_group", true)]
    #[case::with_dash("read-group", true)]
    #[case::empty("", false)]
    #[case::leading_digit("2fast", false)]
    #[case::leading_dash("-sort", false)]
    #[case::colon("so::rt", false)]
    #[case::space("so rt", false)]
    #[case::non_ascii("sørt", false)]
    fn validate_prefix_accepts_only_identifier_fragments(
        #[case] prefix: &str,
        #[case] expected_ok: bool,
    ) {
        assert_eq!(validate_prefix(prefix).is_ok(), expected_ok, "prefix {prefix:?}");
    }

    fn parse_type(tokens: TokenStream2) -> syn::Type {
        syn::parse2(tokens).expect("parse type")
    }

    /// `is_vec_type` drives whether a field is absent-able, so a non-path type —
    /// reference, array, tuple, slice — must answer `false` rather than panic or
    /// match on the last path segment of something that has none.
    #[rstest]
    #[case::vec(quote! { Vec<u32> }, true)]
    #[case::qualified_vec(quote! { std::vec::Vec<u32> }, true)]
    #[case::option_not_vec(quote! { Option<u32> }, false)]
    #[case::plain(quote! { u32 }, false)]
    #[case::reference(quote! { &str }, false)]
    #[case::array(quote! { [u8; 4] }, false)]
    #[case::tuple(quote! { (u32, u32) }, false)]
    #[case::unit(quote! { () }, false)]
    fn is_vec_type_only_matches_path_types_named_vec(
        #[case] tokens: TokenStream2,
        #[case] expected: bool,
    ) {
        assert_eq!(is_vec_type(&parse_type(tokens)), expected);
    }

    /// Same for `is_option_type`, which decides whether a field is treated as
    /// optional.
    #[rstest]
    #[case::option(quote! { Option<u32> }, true)]
    #[case::qualified_option(quote! { std::option::Option<u32> }, true)]
    #[case::vec_not_option(quote! { Vec<u32> }, false)]
    #[case::plain(quote! { u32 }, false)]
    #[case::reference(quote! { &str }, false)]
    #[case::array(quote! { [u8; 4] }, false)]
    #[case::tuple(quote! { (u32, u32) }, false)]
    #[case::unit(quote! { () }, false)]
    fn is_option_type_only_matches_path_types_named_option(
        #[case] tokens: TokenStream2,
        #[case] expected: bool,
    ) {
        assert_eq!(is_option_type(&parse_type(tokens)), expected);
    }

    /// `is_bool_type` keeps clap's valueless `SetTrue` idiom out of the required
    /// classification. `Option<bool>` and `Vec<bool>` are *not* bare bools.
    #[rstest]
    #[case::bare_bool(quote! { bool }, true)]
    #[case::qualified_bool(quote! { core::primitive::bool }, true)]
    #[case::option_bool(quote! { Option<bool> }, false)]
    #[case::vec_bool(quote! { Vec<bool> }, false)]
    #[case::plain(quote! { u32 }, false)]
    #[case::reference(quote! { &bool }, false)]
    #[case::tuple(quote! { (bool, bool) }, false)]
    #[case::unit(quote! { () }, false)]
    fn is_bool_type_only_matches_the_bare_bool_path(
        #[case] tokens: TokenStream2,
        #[case] expected: bool,
    ) {
        assert_eq!(is_bool_type(&parse_type(tokens)), expected);
    }
}
