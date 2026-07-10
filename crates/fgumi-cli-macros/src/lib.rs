#![deny(unsafe_code)]

//! Proc-macro support for fgumi CLI option re-exposure.
//!
//! Provides [`multi_options`] — an attribute macro that pairs a
//! standalone command's `clap::Args` options struct with a sibling
//! `Multi<Name>` struct whose CLI flags are prefixed (`--<prefix>::<flag>`)
//! and grouped under a help heading. Used by the `runall` command to
//! re-expose every per-stage option of `fgumi sort` / `fgumi group` /
//! `fgumi simplex` / `fgumi duplex` / `fgumi codec` without
//! hand-maintaining a parallel option set on `RunAll`.
//!
//! Field-kind handling:
//!
//!   * `#[arg(skip)]` / `#[arg(skip = expr)]` fields are invisible to the
//!     CLI and are not emitted into the Multi struct. `validate()` fills
//!     them from the original struct's `Default::default()`.
//!   * `Option<T>` fields stay `Option<T>` in the Multi struct.
//!   * Non-`Option<T>` fields with `default_value_t` keep their type
//!     and default — `default_value_t` is rewritten to read the
//!     default from the original struct's `Default::default()`.
//!   * Non-`Option<T>` fields without `default_value_t` are required:
//!     they become `Option<T>` in the Multi struct, and the generated
//!     `validate()` method returns an error when the corresponding
//!     `--<prefix>::<flag>` is missing.
//!
//! The original options struct is emitted unchanged. The generated
//! `Multi<Name>` struct has a `validate(self) -> anyhow::Result<<Name>>`
//! method to convert back to the original struct, and a
//! `From<<Name>> for Multi<Name>` impl for the reverse direction
//! (useful in tests).

use proc_macro::TokenStream;
use quote::{ToTokens, format_ident, quote};
use syn::{Lit, Meta};

/// Attribute macro that passes through the original options struct unchanged
/// and directly generates a `Multi<OriginalName>` struct with
/// `--<prefix>::<field-name>` CLI flags, plus a `validate()` method and a
/// `From<<OriginalName>>` conversion.
///
/// Usage: `#[multi_options("prefix", "Help Heading")]`
///
/// See the crate-level docs for the field-kind classification rules.
///
/// # Panics
///
/// Panics (with a clear build-time message) if applied to:
/// * a non-struct item, or a struct without named fields;
/// * a struct containing a `#[command(flatten)]` field (the struct must be
///   flat — re-exposing nested clap-flattened state is out of scope);
/// * a field with a cross-field clap reference attribute — `requires`,
///   `conflicts_with`, `required_if_eq`, `required_unless_present`, … — whose
///   arg-id string would dangle once the field is prefixed. Enforce such
///   couplings in the command's `validate()`/`resolve()` instead;
/// * a field using clap's `key(value)` call form for `long`, `short`,
///   `default_value`, `default_value_t`, or `required` — use the `key = value`
///   (or bare `key`) spelling, which is the only form the macro classifies.
#[proc_macro_attribute]
pub fn multi_options(attr: TokenStream, item: TokenStream) -> TokenStream {
    let args = syn::parse::<MultiOptionsArgs>(attr).expect(
        "multi_options requires two string literal arguments: \
         #[multi_options(\"prefix\", \"heading\")]",
    );
    let prefix = &args.prefix;
    let heading = &args.heading;

    let input = syn::parse::<syn::ItemStruct>(item).expect("multi_options requires a struct");
    let struct_name = &input.ident;
    let multi_struct_name = format_ident!("Multi{}", struct_name);

    let fields = match &input.fields {
        syn::Fields::Named(named) => &named.named,
        _ => panic!("multi_options only supports structs with named fields"),
    };

    // Reject any #[command(flatten)] inside the original struct —
    // re-prefixing fields of a nested clap struct is out of scope.
    for field in fields {
        for attr in &field.attrs {
            if attr.path().is_ident("command") {
                let tokens = attr.meta.to_token_stream().to_string();
                assert!(
                    !tokens.contains("flatten"),
                    "multi_options does not support #[command(flatten)] on field `{}`. \
                     Inline the flattened struct's fields directly.",
                    field.ident.as_ref().unwrap()
                );
            }
        }
        // Reject #[arg(...)] forms the macro cannot faithfully re-expose on the
        // Multi struct (dangling cross-references, and the `key(value)` call
        // form) — fail loud here rather than panic in clap or silently change
        // behavior at runall parse time.
        if let Some(err) = unsupported_arg_attr_error(field) {
            panic!("{err}");
        }
    }

    let mut multi_field_tokens = Vec::new();
    let mut validate_field_tokens = Vec::new();
    let mut from_original_field_tokens = Vec::new();

    for field in fields {
        generate_field_tokens(
            field,
            struct_name,
            prefix,
            &mut multi_field_tokens,
            &mut validate_field_tokens,
            &mut from_original_field_tokens,
        );
    }

    let original = input.to_token_stream();

    let expanded = quote! {
        #original

        /// Prefixed options struct generated by `#[multi_options]` for
        /// the runall command. Carries the same fields as the original
        /// options struct but exposes each via `--<prefix>::<flag>` and
        /// groups them under a help heading.
        #[derive(::clap::Args, Debug, Clone)]
        #[command(next_help_heading = #heading)]
        pub struct #multi_struct_name {
            #(#multi_field_tokens)*
        }

        impl #multi_struct_name {
            /// Validate required fields and convert to the original
            /// options struct. Returns `Err` if any field that was
            /// required (no default + non-`Option`) was not supplied
            /// on the command line.
            pub fn validate(self) -> ::anyhow::Result<#struct_name> {
                let opts = self;
                Ok(#struct_name {
                    #(#validate_field_tokens)*
                })
            }
        }

        impl From<#struct_name> for #multi_struct_name {
            fn from(opts: #struct_name) -> Self {
                Self {
                    #(#from_original_field_tokens)*
                }
            }
        }
    };

    TokenStream::from(expanded)
}

/// Returns `true` if the field has `#[arg(skip)]` or `#[arg(skip = ...)]`.
///
/// Skip fields are invisible to the CLI; the original struct's `Default`
/// impl (or the skip expression) provides their value. The Multi struct
/// does not expose them as CLI flags and `validate()` initialises them
/// from `Default::default()`.
fn is_skip_field(field: &syn::Field) -> bool {
    for attr in &field.attrs {
        if !attr.path().is_ident("arg") {
            continue;
        }
        if let Ok(nested) = attr
            .parse_args_with(syn::punctuated::Punctuated::<Meta, syn::Token![,]>::parse_terminated)
        {
            for meta in &nested {
                match meta {
                    Meta::Path(path) if path.is_ident("skip") => return true,
                    Meta::NameValue(nv) if nv.path.is_ident("skip") => return true,
                    _ => {}
                }
            }
        }
    }
    false
}

/// Extract the expression from `#[arg(skip = expr)]`, if present.
///
/// Returns `Some(expr)` for the name-value form `#[arg(skip = some_expr)]`
/// and `None` for the bare-path form `#[arg(skip)]` (or no skip at all).
/// The returned expression is used verbatim in `validate()` so the skip
/// value is preserved instead of being overwritten by the struct's
/// `Default`.
fn extract_skip_expr(field: &syn::Field) -> Option<syn::Expr> {
    for attr in &field.attrs {
        if !attr.path().is_ident("arg") {
            continue;
        }
        if let Ok(nested) = attr
            .parse_args_with(syn::punctuated::Punctuated::<Meta, syn::Token![,]>::parse_terminated)
        {
            for meta in &nested {
                if let Meta::NameValue(nv) = meta {
                    if nv.path.is_ident("skip") {
                        return Some(nv.value.clone());
                    }
                }
            }
        }
    }
    None
}

/// Generate Multi-struct field, validate assignment, and From<Original> assignment
/// for a single field of the Options struct.
fn generate_field_tokens(
    field: &syn::Field,
    struct_name: &syn::Ident,
    prefix: &str,
    multi_field_tokens: &mut Vec<proc_macro2::TokenStream>,
    validate_field_tokens: &mut Vec<proc_macro2::TokenStream>,
    from_original_field_tokens: &mut Vec<proc_macro2::TokenStream>,
) {
    let field_ident = field.ident.as_ref().expect("named field");
    let field_type = &field.ty;

    // Skip fields are not exposed as CLI flags. Use the original struct's
    // Default to populate them in validate(). The Multi struct does not
    // carry skip fields at all, so no multi_field_tokens or
    // from_original_field_tokens entry is needed.
    if is_skip_field(field) {
        // Honor an explicit `#[arg(skip = expr)]`: use the provided
        // expression verbatim rather than `#struct_name::default()`, which
        // would discard the intended value and force a `Default` impl on the
        // whole struct. A bare `#[arg(skip)]` falls back to the original
        // struct's `Default`, matching clap's own bare-skip semantics.
        let value = if let Some(expr) = extract_skip_expr(field) {
            quote! { #expr }
        } else {
            quote! { #struct_name::default().#field_ident }
        };
        validate_field_tokens.push(quote! {
            #field_ident: #value,
        });
        // No multi_field_tokens (not in Multi struct).
        // No from_original_field_tokens (not in Multi struct).
        return;
    }

    // Honor any explicit `#[arg(long = "...")]` override on the
    // original field; fall back to kebab(field_ident) when no
    // override is set. This keeps the Multi-side flag name aligned
    // with the standalone command's flag name (e.g. `tmp_dirs` with
    // `long = "tmp-dir"` becomes `--<prefix>::tmp-dir`, not
    // `--<prefix>::tmp-dirs`).
    let base_name =
        extract_long_override(field).unwrap_or_else(|| field_ident.to_string().replace('_', "-"));
    let doc_comment = extract_doc_string(&field.attrs);
    let is_option = is_option_type(field_type);
    let is_vec = is_vec_type(field_type);
    let has_default = has_default_value_t(field);
    let preserved = extract_preserved_arg_attrs(field);

    let long_name = format!("{prefix}::{base_name}");
    let prefixed_ident = format_ident!("{}_{}", prefix.replace('-', "_"), field_ident);

    // A field is "required" if it's not Option<T>, not Vec<T>, and has
    // no `default_value` / `default_value_t`. `Vec<T>` collects via
    // `ArgAction::Append` so the natural default is an empty Vec; we
    // pass it through without rewriting `default_value_t`.
    let is_required = !is_option && !is_vec && !has_default;

    if is_required {
        // Avoid a leading space when the source field has no doc comment
        // (`doc_comment` is then the empty string).
        let required_doc = if doc_comment.is_empty() {
            format!("Required when {prefix} is selected.")
        } else {
            format!("{doc_comment} Required when {prefix} is selected.")
        };
        let error_msg = format!("--{long_name} is required when {prefix} is selected");

        multi_field_tokens.push(quote! {
            #[doc = #required_doc]
            #[arg(long = #long_name #(#preserved)*)]
            pub #prefixed_ident: Option<#field_type>,
        });
        validate_field_tokens.push(quote! {
            #field_ident: opts.#prefixed_ident
                .ok_or_else(|| ::anyhow::anyhow!(#error_msg))?,
        });
        from_original_field_tokens.push(quote! {
            #prefixed_ident: Some(opts.#field_ident),
        });
    } else if is_option || is_vec {
        // `Option<T>` and `Vec<T>` both pass through on the Multi side:
        // clap handles missing-on-cli as `None` / empty Vec respectively,
        // and no `default_value_t` rewrite is possible (clap rejects
        // `default_value_t` on `Vec<T>`, and the macro strips any string
        // `default_value`).
        let doc_attr = doc_attr_tokens(&doc_comment);
        multi_field_tokens.push(quote! {
            #doc_attr
            #[arg(long = #long_name #(#preserved)*)]
            pub #prefixed_ident: #field_type,
        });
        if is_vec {
            // Because clap can't carry a Vec default through to the Multi
            // side, an omitted flag arrives as an empty Vec. Backfill the
            // original struct's `Default` value when empty so the runall
            // surface matches the standalone command without a per-call
            // backfill in the caller. A genuinely-required Vec (e.g.
            // filter's `min_reads`, whose `Default` is also empty) stays
            // empty here, so the caller's own validation surfaces the
            // "required" error exactly as the standalone command does.
            validate_field_tokens.push(quote! {
                #field_ident: if opts.#prefixed_ident.is_empty() {
                    #struct_name::default().#field_ident
                } else {
                    opts.#prefixed_ident
                },
            });
        } else {
            validate_field_tokens.push(quote! {
                #field_ident: opts.#prefixed_ident,
            });
        }
        from_original_field_tokens.push(quote! {
            #prefixed_ident: opts.#field_ident,
        });
    } else {
        // Defaulted field. The Multi-side reads its `default_value_t`
        // from the original struct's `Default::default()`, so the
        // type MUST implement both `Default` and `Display` for the
        // generated code to compile and match the original's CLI
        // default.
        let default_attr = quote! { , default_value_t = #struct_name::default().#field_ident };
        let doc_attr = doc_attr_tokens(&doc_comment);
        multi_field_tokens.push(quote! {
            #doc_attr
            #[arg(long = #long_name #default_attr #(#preserved)*)]
            pub #prefixed_ident: #field_type,
        });
        validate_field_tokens.push(quote! {
            #field_ident: opts.#prefixed_ident,
        });
        from_original_field_tokens.push(quote! {
            #prefixed_ident: opts.#field_ident,
        });
    }
}

/// Parsed arguments for `#[multi_options("prefix", "heading")]`.
struct MultiOptionsArgs {
    prefix: String,
    heading: String,
}

impl syn::parse::Parse for MultiOptionsArgs {
    fn parse(input: syn::parse::ParseStream) -> syn::Result<Self> {
        let prefix_lit: syn::LitStr = input.parse()?;
        input.parse::<syn::Token![,]>()?;
        let heading_lit: syn::LitStr = input.parse()?;
        Ok(Self { prefix: prefix_lit.value(), heading: heading_lit.value() })
    }
}

/// Check whether a field has a default — either `default_value_t = …`
/// (Rust-typed default) or `default_value = "…"` (string default that
/// clap parses via the `value_parser`). Both indicate the field is
/// optional on the original CLI; the Multi-side rewrites it to
/// `default_value_t = Self::default().field` regardless.
fn has_default_value_t(field: &syn::Field) -> bool {
    for attr in &field.attrs {
        if !attr.path().is_ident("arg") {
            continue;
        }
        if let Ok(nested) = attr
            .parse_args_with(syn::punctuated::Punctuated::<Meta, syn::Token![,]>::parse_terminated)
        {
            for meta in &nested {
                match meta {
                    Meta::NameValue(nv)
                        if nv.path.is_ident("default_value_t")
                            || nv.path.is_ident("default_value") =>
                    {
                        return true;
                    }
                    Meta::Path(path)
                        if path.is_ident("default_value_t") || path.is_ident("default_value") =>
                    {
                        return true;
                    }
                    _ => {}
                }
            }
        }
    }
    false
}

/// clap `#[arg(...)]` keys that reference *another argument by its string id*.
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

/// clap `#[arg(...)]` keys this macro classifies via `Meta::NameValue` / `Meta::Path`
/// (to strip, rewrite, or read them). clap also accepts the equivalent `key(value)`
/// call form, which arrives as `Meta::List` and would slip past every classifier —
/// silently changing the flag name, the required-ness, or emitting a duplicate
/// `long`. Reject the call form for these keys (D2).
const CALL_FORM_SENSITIVE_ARG_KEYS: &[&str] =
    &["long", "short", "default_value", "default_value_t", "required"];

/// Returns an error message if a field carries an `#[arg(...)]` attribute that
/// `multi_options` cannot faithfully re-expose on the Multi struct, or `None`
/// when the field is safe.
///
/// Covers two latent traps: (D1) cross-field reference attrs whose arg-id
/// strings dangle after the field is prefixed, and (D2) clap's `key(value)`
/// call form for keys the macro only recognizes in `key = value` / bare-`key`
/// form.
fn unsupported_arg_attr_error(field: &syn::Field) -> Option<String> {
    let field_name =
        field.ident.as_ref().map_or_else(|| "<unnamed>".to_string(), ToString::to_string);
    for attr in &field.attrs {
        if !attr.path().is_ident("arg") {
            continue;
        }
        let nested = attr
            .parse_args_with(syn::punctuated::Punctuated::<Meta, syn::Token![,]>::parse_terminated)
            .unwrap_or_else(|e| {
                panic!("multi_options: failed to parse #[arg(...)] on field `{field_name}`: {e}")
            });
        for meta in &nested {
            let (ident, is_call_form) = match meta {
                Meta::NameValue(nv) => (nv.path.get_ident().map(ToString::to_string), false),
                Meta::Path(p) => (p.get_ident().map(ToString::to_string), false),
                Meta::List(l) => (l.path.get_ident().map(ToString::to_string), true),
            };
            let Some(ident) = ident else { continue };
            if CROSS_REFERENCE_ARG_KEYS.contains(&ident.as_str()) {
                return Some(format!(
                    "multi_options: field `{field_name}` uses #[arg({ident} …)], which references \
                     another argument by its unprefixed id. The Multi struct renames fields to \
                     `<prefix>_<field>`, so that id would dangle and clap would panic when it \
                     builds the runall command. Enforce this coupling in the command's \
                     validate()/resolve() instead (see AlignerOptions::resolve)."
                ));
            }
            if is_call_form && CALL_FORM_SENSITIVE_ARG_KEYS.contains(&ident.as_str()) {
                return Some(format!(
                    "multi_options: field `{field_name}` uses the call form #[arg({ident}(…))]. \
                     Use the `{ident} = …` name-value form (or bare `{ident}`) — the macro only \
                     classifies those spellings and would mishandle the call form."
                ));
            }
        }
    }
    None
}

/// Check whether a type is `Vec<T>`. The macro treats `Vec<T>` like
/// `Option<T>`: pass-through with no `default_value_t` rewrite, since
/// clap collects empty Vecs naturally and rejects `default_value_t`
/// on Vec fields.
fn is_vec_type(ty: &syn::Type) -> bool {
    if let syn::Type::Path(type_path) = ty {
        if let Some(segment) = type_path.path.segments.last() {
            return segment.ident == "Vec";
        }
    }
    false
}

/// Check whether a type is `Option<T>`.
fn is_option_type(ty: &syn::Type) -> bool {
    if let syn::Type::Path(type_path) = ty {
        if let Some(segment) = type_path.path.segments.last() {
            return segment.ident == "Option";
        }
    }
    false
}

/// Extract `#[arg(...)]` attributes that should be preserved on the
/// generated Multi-struct field (everything except `long`, `short`,
/// and `default_value_t`). Returned as a `Vec<TokenStream>` ready to
/// be spliced into the `#[arg(...)]` attribute of the new field as
/// `, key = value` continuations.
fn extract_preserved_arg_attrs(field: &syn::Field) -> Vec<proc_macro2::TokenStream> {
    let mut preserved = Vec::new();
    for attr in &field.attrs {
        if !attr.path().is_ident("arg") {
            continue;
        }
        if let Ok(nested) = attr
            .parse_args_with(syn::punctuated::Punctuated::<Meta, syn::Token![,]>::parse_terminated)
        {
            for meta in &nested {
                match meta {
                    Meta::NameValue(nv) => {
                        let name = nv.path.get_ident().map(ToString::to_string);
                        // Strip `required = ...` here for the same reason the
                        // bare `required` path is stripped below: the Multi-side
                        // handles required-ness via the `Option<T>` wrap +
                        // `validate()` path, so propagating `required = true`
                        // would make clap enforce the requirement during
                        // parsing and bypass the staged validation.
                        if !matches!(
                            name.as_deref(),
                            Some(
                                "long" | "short" | "default_value_t" | "default_value" | "required"
                            )
                        ) {
                            preserved.push(quote! { , #nv });
                        }
                    }
                    Meta::Path(path) => {
                        let name = path.get_ident().map(ToString::to_string);
                        // Strip `long`/`short` (we set our own `long`
                        // and never propagate short flags), and strip
                        // bare `required` — the Multi-side handles
                        // required-ness via the `Option<T>` wrap +
                        // `validate()` path, so a passthrough
                        // `required` would conflict.
                        if !matches!(name.as_deref(), Some("long" | "short" | "required")) {
                            preserved.push(quote! { , #path });
                        }
                    }
                    Meta::List(_) => {
                        preserved.push(quote! { , #meta });
                    }
                }
            }
        }
    }
    preserved
}

/// Read the `#[arg(long = "...")]` override on a field, if any.
///
/// Returns the string value when the field has an explicit
/// `#[arg(long = "tmp-dir")]` (or similar) attribute. Returns `None`
/// for the bare-path form (`#[arg(long)]`, meaning "use the
/// kebab of the field name") or when no `long` is set at all.
fn extract_long_override(field: &syn::Field) -> Option<String> {
    for attr in &field.attrs {
        if !attr.path().is_ident("arg") {
            continue;
        }
        if let Ok(nested) = attr
            .parse_args_with(syn::punctuated::Punctuated::<Meta, syn::Token![,]>::parse_terminated)
        {
            for meta in &nested {
                if let Meta::NameValue(nv) = meta {
                    if nv.path.is_ident("long") {
                        if let syn::Expr::Lit(expr_lit) = &nv.value {
                            if let Lit::Str(lit) = &expr_lit.lit {
                                return Some(lit.value());
                            }
                        }
                    }
                }
            }
        }
    }
    None
}

/// Build a `#[doc = "..."]` attribute token for a Multi-struct field, or an
/// empty token stream when the source field had no doc comment.
///
/// Emitting `#[doc = ""]` for an undocumented field is harmless but produces a
/// stray empty doc attribute in the generated code; returning nothing keeps the
/// generated output clean and avoids a leading-space artifact in concatenated
/// doc strings.
fn doc_attr_tokens(doc_comment: &str) -> proc_macro2::TokenStream {
    if doc_comment.is_empty() {
        quote! {}
    } else {
        quote! { #[doc = #doc_comment] }
    }
}

/// Join all `#[doc = "..."]` attributes into a single trimmed string.
fn extract_doc_string(attrs: &[syn::Attribute]) -> String {
    let lines: Vec<String> = attrs
        .iter()
        .filter_map(|attr| {
            if !attr.path().is_ident("doc") {
                return None;
            }
            if let Meta::NameValue(nv) = &attr.meta {
                if let syn::Expr::Lit(expr_lit) = &nv.value {
                    if let Lit::Str(lit) = &expr_lit.lit {
                        return Some(lit.value());
                    }
                }
            }
            None
        })
        .collect();

    lines.iter().map(|l| l.trim()).collect::<Vec<_>>().join(" ").trim().to_string()
}

#[cfg(test)]
mod tests {
    use super::*;
    use rstest::rstest;
    use syn::parse::Parser;

    /// Parse a single named struct field from tokens (e.g.
    /// `#[arg(long, requires = "x")] pub y: u32`).
    fn named_field(tokens: proc_macro2::TokenStream) -> syn::Field {
        syn::Field::parse_named.parse2(tokens).expect("parse named field")
    }

    #[test]
    fn accepts_supported_arg_forms() {
        // Bare + name-value forms the macro classifies are all fine.
        for tokens in [
            quote! { #[arg(long, short = 'x', default_value = "7")] pub a: u32 },
            quote! { #[arg(long = "max-memory", default_value_t = 5)] pub b: usize },
            quote! { #[arg(long)] pub c: Option<u32> },
            quote! { #[arg(long, value_delimiter = ',', required = true)] pub d: Vec<usize> },
            quote! { pub e: u32 },
        ] {
            let field = named_field(tokens);
            assert_eq!(
                unsupported_arg_attr_error(&field),
                None,
                "expected supported form to be accepted"
            );
        }
    }

    // D1: attrs that reference another arg by id would dangle after prefixing.
    // Each row is its own `#[case]` so a failure names the scenario. The
    // `needle` asserts the message names the offending attr key; `field_ident`
    // asserts it names the actual field (the prior tautological `contains('a')
    // || …` check passed for almost any string).
    #[rstest]
    #[case::requires(
        quote! { #[arg(long, requires = "other")] pub a: u32 },
        "a",
        "requires"
    )]
    #[case::conflicts_with(
        quote! { #[arg(long, conflicts_with = "other")] pub b: u32 },
        "b",
        "conflicts_with"
    )]
    #[case::required_if_eq(
        quote! { #[arg(long, required_if_eq("mode", "x"))] pub c: u32 },
        "c",
        "required_if_eq"
    )]
    #[case::required_unless_present(
        quote! { #[arg(long, required_unless_present = "other")] pub d: u32 },
        "d",
        "required_unless_present"
    )]
    #[case::default_values_if(
        quote! { #[arg(long, default_values_if("mode", "x", Some("v")))] pub e: u32 },
        "e",
        "default_values_if"
    )]
    #[case::default_values_ifs(
        quote! { #[arg(long, default_values_ifs([("mode", "x", Some("v"))]))] pub f: u32 },
        "f",
        "default_values_ifs"
    )]
    fn rejects_cross_reference_attrs(
        #[case] tokens: proc_macro2::TokenStream,
        #[case] field_ident: &str,
        #[case] needle: &str,
    ) {
        let field = named_field(tokens);
        let err =
            unsupported_arg_attr_error(&field).expect("cross-reference attr must be rejected");
        assert!(err.contains(needle), "message should name the attr {needle:?}: {err}");
        assert!(
            err.contains(&format!("field `{field_ident}`")),
            "message should name the field {field_ident:?}: {err}"
        );
    }

    #[test]
    fn rejects_call_form_for_classified_keys() {
        // D2: the `key(value)` call form slips past the NameValue/Path classifiers.
        for (tokens, needle) in [
            (quote! { #[arg(long("x"))] pub a: u32 }, "long"),
            (quote! { #[arg(long, short('x'))] pub b: u32 }, "short"),
            (quote! { #[arg(default_value_t(4))] pub c: u32 }, "default_value_t"),
            (quote! { #[arg(long, required(true))] pub d: u32 }, "required"),
        ] {
            let field = named_field(tokens);
            let err = unsupported_arg_attr_error(&field).expect("call form must be rejected");
            assert!(err.contains(needle), "message should name the key {needle:?}: {err}");
        }
    }

    #[test]
    fn allows_call_form_for_unclassified_keys() {
        // `value_parser(...)` in call form is preserved verbatim and is not one of
        // the classified keys, so it must not be rejected.
        let field =
            named_field(quote! { #[arg(long, value_parser(clap::value_parser!(u32)))] pub a: u32 });
        assert_eq!(unsupported_arg_attr_error(&field), None);
    }
}
