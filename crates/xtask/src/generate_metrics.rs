use anyhow::{Context, Result};
use std::fmt::Write;
use std::fs;
use std::path::Path;

/// A parsed metric page with its metadata.
pub struct MetricPage {
    /// Struct name (e.g., `UmiGroupingMetrics`).
    pub name: String,
    /// Relative path within `docs/src` (e.g., "metrics/umi-grouping-metrics.md").
    pub path: String,
}

/// A metric struct parsed from source code.
struct ParsedMetric {
    /// Struct name.
    name: String,
    /// Struct-level doc comment.
    doc: String,
    /// Fields: (name, type, doc comment).
    fields: Vec<(String, String, String)>,
}

/// Generate metric reference markdown pages by parsing source files with `syn`.
pub fn generate(docs_src: &Path) -> Result<Vec<MetricPage>> {
    let metrics_dir = docs_src.join("metrics");
    if metrics_dir.exists() {
        fs::remove_dir_all(&metrics_dir)?;
    }
    fs::create_dir_all(&metrics_dir)?;

    let metrics_src = Path::new("crates/fgumi-metrics/src");
    let mut all_metrics = Vec::new();

    for entry in fs::read_dir(metrics_src)? {
        let entry = entry?;
        let path = entry.path();
        if path.extension().is_some_and(|e| e == "rs") {
            let filename = path
                .file_name()
                .expect("path from read_dir should have a file name")
                .to_string_lossy()
                .to_string();
            if filename == "lib.rs" || filename == "writer.rs" || filename == "rejection.rs" {
                continue;
            }
            let source =
                fs::read_to_string(&path).with_context(|| format!("reading {}", path.display()))?;
            let parsed = parse_metric_structs(&source)?;
            all_metrics.extend(parsed);
        }
    }

    let mut pages = Vec::new();

    for metric in &all_metrics {
        let markdown = render_metric_page(metric);
        let slug = name_to_slug(&metric.name);
        let filename = format!("{slug}.md");
        fs::write(metrics_dir.join(&filename), &markdown)?;

        pages.push(MetricPage { name: metric.name.clone(), path: format!("metrics/{filename}") });
    }

    let index = render_metrics_index(&all_metrics);
    fs::write(metrics_dir.join("README.md"), index)?;

    Ok(pages)
}

/// Parse a Rust source file to find structs that implement `Metric`.
fn parse_metric_structs(source: &str) -> Result<Vec<ParsedMetric>> {
    let file = syn::parse_file(source)?;

    // Find struct names that have `impl Metric for StructName`
    let metric_impl_names: Vec<String> = file
        .items
        .iter()
        .filter_map(|item| {
            if let syn::Item::Impl(impl_item) = item {
                if let Some((_, trait_path, _)) = &impl_item.trait_ {
                    let trait_name = trait_path.segments.last().map(|s| s.ident.to_string())?;
                    if trait_name == "Metric" {
                        if let syn::Type::Path(type_path) = impl_item.self_ty.as_ref() {
                            return type_path.path.segments.last().map(|s| s.ident.to_string());
                        }
                    }
                }
            }
            None
        })
        .collect();

    let mut metrics = Vec::new();

    for item in &file.items {
        if let syn::Item::Struct(struct_item) = item {
            let name = struct_item.ident.to_string();
            if !metric_impl_names.contains(&name) {
                continue;
            }

            let doc = extract_doc_comments(&struct_item.attrs);
            let fields = extract_fields(struct_item);
            metrics.push(ParsedMetric { name, doc, fields });
        }
    }

    Ok(metrics)
}

/// Extract doc comments from attributes.
fn extract_doc_comments(attrs: &[syn::Attribute]) -> String {
    let lines: Vec<String> = attrs
        .iter()
        .filter_map(|attr| {
            if attr.path().is_ident("doc") {
                if let syn::Meta::NameValue(nv) = &attr.meta {
                    if let syn::Expr::Lit(lit) = &nv.value {
                        if let syn::Lit::Str(s) = &lit.lit {
                            return Some(s.value());
                        }
                    }
                }
            }
            None
        })
        .collect();

    lines
        .iter()
        .map(|l| l.strip_prefix(' ').unwrap_or(l))
        .collect::<Vec<_>>()
        .join("\n")
        .trim()
        .to_string()
}

/// Extract field names, types, and doc comments from a struct.
fn extract_fields(struct_item: &syn::ItemStruct) -> Vec<(String, String, String)> {
    match &struct_item.fields {
        syn::Fields::Named(fields) => fields
            .named
            .iter()
            .map(|f| {
                let name = f.ident.as_ref().map(ToString::to_string).unwrap_or_default();
                let ty = type_to_string(&f.ty);
                let doc = extract_doc_comments(&f.attrs);
                (name, ty, doc)
            })
            .collect(),
        _ => Vec::new(),
    }
}

/// Convert a `syn::Type` to a readable string, preserving generic arguments.
fn type_to_string(ty: &syn::Type) -> String {
    match ty {
        syn::Type::Path(p) => {
            let segments: Vec<_> = p
                .path
                .segments
                .iter()
                .map(|s| {
                    let ident = s.ident.to_string();
                    match &s.arguments {
                        syn::PathArguments::AngleBracketed(args) => {
                            let inner: Vec<_> = args
                                .args
                                .iter()
                                .map(|arg| match arg {
                                    syn::GenericArgument::Type(t) => type_to_string(t),
                                    _ => quote::quote!(#arg).to_string(),
                                })
                                .collect();
                            format!("{}<{}>", ident, inner.join(", "))
                        }
                        _ => ident,
                    }
                })
                .collect();
            segments.join("::")
        }
        _ => quote::quote!(#ty).to_string(),
    }
}

/// Render a single metric's markdown page.
fn render_metric_page(metric: &ParsedMetric) -> String {
    let mut md = String::new();

    let _ = writeln!(md, "# {}\n", metric.name);

    if !metric.doc.is_empty() {
        let _ = writeln!(md, "{}\n", metric.doc);
    }

    if !metric.fields.is_empty() {
        md.push_str("## Fields\n\n");
        md.push_str("| Column | Type | Description |\n");
        md.push_str("|--------|------|-------------|\n");

        for (name, ty, doc) in &metric.fields {
            let doc_oneline = doc.replace('\n', " ").replace('|', "\\|");
            let _ = writeln!(md, "| `{name}` | `{ty}` | {doc_oneline} |");
        }
        md.push('\n');
    }

    md
}

/// Render the metrics index page.
fn render_metrics_index(metrics: &[ParsedMetric]) -> String {
    let mut md = String::new();
    md.push_str("# Metrics Reference\n\n");
    md.push_str("Auto-generated from fgumi metric struct definitions.\n\n");
    md.push_str("| Metric | Description |\n");
    md.push_str("|--------|-------------|\n");

    for metric in metrics {
        let slug = name_to_slug(&metric.name);
        let first_line = metric.doc.lines().next().unwrap_or("");
        let _ = writeln!(md, "| [{}]({slug}.md) | {first_line} |", metric.name);
    }

    md.push('\n');
    md
}

/// Convert a `CamelCase` struct name to a kebab-case slug.
fn name_to_slug(name: &str) -> String {
    let mut slug = String::new();
    for (i, c) in name.chars().enumerate() {
        if c.is_uppercase() && i > 0 {
            slug.push('-');
        }
        slug.push(c.to_ascii_lowercase());
    }
    slug
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_name_to_slug() {
        assert_eq!(name_to_slug("UmiGroupingMetrics"), "umi-grouping-metrics");
        assert_eq!(name_to_slug("FamilySizeMetric"), "family-size-metric");
        assert_eq!(name_to_slug("ConsensusMetrics"), "consensus-metrics");
    }

    #[test]
    fn test_parse_metric_struct() {
        let source = r#"
            use serde::{Serialize, Deserialize};

            /// Metrics about UMI grouping.
            ///
            /// Detailed description here.
            #[derive(Debug, Clone, Serialize, Deserialize)]
            pub struct MyMetric {
                /// Total number of records processed.
                pub total_records: u64,
                /// Fraction of records accepted.
                pub fraction_accepted: f64,
            }

            impl Metric for MyMetric {
                fn metric_name() -> &'static str { "my metric" }
            }
        "#;

        let metrics = parse_metric_structs(source).unwrap();
        assert_eq!(metrics.len(), 1);
        assert_eq!(metrics[0].name, "MyMetric");
        assert!(metrics[0].doc.contains("Metrics about UMI grouping"));
        assert_eq!(metrics[0].fields.len(), 2);
        assert_eq!(metrics[0].fields[0].0, "total_records");
        assert_eq!(metrics[0].fields[0].1, "u64");
        assert!(metrics[0].fields[0].2.contains("Total number"));
        assert_eq!(metrics[0].fields[1].0, "fraction_accepted");
        assert_eq!(metrics[0].fields[1].1, "f64");
    }

    #[test]
    fn test_parse_ignores_non_metric_structs() {
        let source = r"
            #[derive(Debug)]
            pub struct NotAMetric {
                pub x: i32,
            }
        ";

        let metrics = parse_metric_structs(source).unwrap();
        assert!(metrics.is_empty());
    }

    #[test]
    fn test_render_metric_page() {
        let metric = ParsedMetric {
            name: "TestMetric".to_string(),
            doc: "A test metric.".to_string(),
            fields: vec![
                ("count".to_string(), "u64".to_string(), "The count".to_string()),
                ("rate".to_string(), "f64".to_string(), "The rate".to_string()),
            ],
        };
        let md = render_metric_page(&metric);
        assert!(md.contains("# TestMetric"));
        assert!(md.contains("A test metric."));
        assert!(md.contains("| `count` | `u64` | The count |"));
    }

    #[test]
    fn test_type_to_string_preserves_generics() {
        let source = r#"
            pub struct GenericMetric {
                pub opt_count: Option<u64>,
                pub simple: f64,
            }

            impl Metric for GenericMetric {
                fn metric_name() -> &'static str { "generic" }
            }
        "#;

        let metrics = parse_metric_structs(source).unwrap();
        assert_eq!(metrics.len(), 1);
        assert_eq!(metrics[0].fields[0].1, "Option<u64>");
        assert_eq!(metrics[0].fields[1].1, "f64");
    }

    #[test]
    fn test_metrics_index_uses_relative_links() {
        let metrics = vec![ParsedMetric {
            name: "UmiGroupingMetrics".to_string(),
            doc: "Grouping stats.".to_string(),
            fields: vec![],
        }];
        let index = render_metrics_index(&metrics);
        // Link should be relative (no "metrics/" prefix), since README.md lives in metrics/
        assert!(index.contains("(umi-grouping-metrics.md)"));
        assert!(!index.contains("(metrics/umi-grouping-metrics.md)"));
    }
}
