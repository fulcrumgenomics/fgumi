//! The published metric-file column contract.
//!
//! `crates/fgumi-metrics/metric_columns.json` records the ordered columns of every
//! metric TSV fgumi emits. Downstream report tooling pins its schemas to it. Keys are
//! `<namespace>.<file>`, where the namespace names the producing command family:
//!
//! - `group`, `dedup`, `correct`, `clip`, `filter`, `copy_umi`, `retag`, `downsample`,
//!   `review` — the command of that name.
//! - `simplex.*` / `duplex.*` — the files written by `simplex-metrics` / `duplex-metrics`
//!   (and by the consensus callers' inline `--metrics`).
//! - `consensus.stats` — the `--stats` file of `simplex`, `duplex` and `codec`.
//!
//! This module guards three properties of that contract:
//!
//! 1. **Drift** — the committed manifest matches the live struct headers.
//! 2. **Completeness** — every struct in the workspace that derives serde's `Serialize`
//!    is either listed as an emitted metric file or explicitly allowlisted with a
//!    reason, so a new serialized output cannot silently go unpublished. (This keys on
//!    serialized structs, not on call sites: an emitter that hand-formats its TSV with
//!    no struct is invisible to it.)
//! 3. **Float encoding** — every serialized `f64` field carries `crate::float`, so
//!    non-finite values are written as fgbio's `Infinity`/`NaN` tokens (Rust's default
//!    `inf` is unparseable by fgbio's `Metric.read`).
//!
//! Sources are parsed with `syn`, so the scanners are independent of formatting.

use std::collections::{BTreeMap, BTreeSet};
use std::path::{Path, PathBuf};

use fgoxide::io::DelimFile;
use fgumi_lib::commands::copy_umi::CopyUmiMetric;
use fgumi_lib::commands::downsample::DownsampleHistogramMetric;
use fgumi_lib::commands::retag::RetagMetric;
use fgumi_lib::variant_review::ConsensusVariantReviewInfo;
use fgumi_metrics::Metric;
use fgumi_metrics::clip::ClippingMetrics;
use fgumi_metrics::consensus::ConsensusKvMetric;
use fgumi_metrics::correct::UmiCorrectionMetrics;
use fgumi_metrics::dedup::{DeduplicationMetrics, DuplicationLadderMetrics};
use fgumi_metrics::duplex::{
    DuplexFamilySizeMetric, DuplexUmiMetric, DuplexYieldMetric, FamilySizeMetric,
};
use fgumi_metrics::filter_stats::FilterStatsMetrics;
use fgumi_metrics::group::{FamilySizeMetrics, PositionGroupSizeMetrics, UmiGroupingMetrics};
use fgumi_metrics::shared::UmiMetric;
use fgumi_metrics::simplex::{SimplexFamilySizeMetric, SimplexYieldMetric};
use fgumi_metrics::writer::write_metrics_auto;
use serde::Serialize;
use tempfile::NamedTempFile;

const MANIFEST: &str =
    concat!(env!("CARGO_MANIFEST_DIR"), "/crates/fgumi-metrics/metric_columns.json");

/// Serialized structs that are intentionally absent from the manifest, with the reason.
const NOT_EMITTED_AS_A_FILE: &[(&str, &str)] = &[(
    "ConsensusMetrics",
    "counter struct projected into the `consensus.stats` key/value rows (ConsensusKvMetric)",
)];

// ---------------------------------------------------------------------------
// Manifest construction
// ---------------------------------------------------------------------------

fn header_of(path: &Path) -> Vec<String> {
    let content = std::fs::read_to_string(path).expect("read metrics file");
    content.lines().next().expect("header line").split('\t').map(str::to_string).collect()
}

/// Columns emitted for a `Metric` type via the standard metrics writer.
fn columns_of<T: Metric>() -> Vec<String> {
    let tmp = NamedTempFile::new().expect("temp file");
    write_metrics_auto(tmp.path(), &[T::default()]).expect("write metrics");
    header_of(tmp.path())
}

/// Columns emitted for a row type written directly via `DelimFile`.
fn columns_via_delim<T: Serialize>(row: T) -> Vec<String> {
    let tmp = NamedTempFile::new().expect("temp file");
    DelimFile::default().write_tsv(tmp.path(), [row]).expect("write tsv");
    header_of(tmp.path())
}

/// `(manifest key, struct name, columns)` for every metric file fgumi emits.
fn emitted_metric_files() -> Vec<(&'static str, &'static str, Vec<String>)> {
    let kv = ConsensusKvMetric::new("key", "value".to_string(), "description");
    let review_columns: Vec<String> =
        ConsensusVariantReviewInfo::tsv_header().split('\t').map(str::to_string).collect();
    let histogram = || columns_via_delim(DownsampleHistogramMetric::default());
    vec![
        ("clip.metrics", "ClippingMetrics", columns_of::<ClippingMetrics>()),
        ("consensus.stats", "ConsensusKvMetric", columns_via_delim(kv)),
        ("copy_umi.metrics", "CopyUmiMetric", columns_via_delim(CopyUmiMetric::default())),
        ("correct.metrics", "UmiCorrectionMetrics", columns_of::<UmiCorrectionMetrics>()),
        (
            "dedup.duplication_ladder",
            "DuplicationLadderMetrics",
            columns_of::<DuplicationLadderMetrics>(),
        ),
        ("dedup.family_sizes", "FamilySizeMetrics", columns_of::<FamilySizeMetrics>()),
        ("dedup.metrics", "DeduplicationMetrics", columns_of::<DeduplicationMetrics>()),
        ("downsample.histogram_kept", "DownsampleHistogramMetric", histogram()),
        ("downsample.histogram_rejected", "DownsampleHistogramMetric", histogram()),
        (
            "duplex.duplex_family_sizes",
            "DuplexFamilySizeMetric",
            columns_of::<DuplexFamilySizeMetric>(),
        ),
        ("duplex.duplex_umi_counts", "DuplexUmiMetric", columns_of::<DuplexUmiMetric>()),
        ("duplex.duplex_yield_metrics", "DuplexYieldMetric", columns_of::<DuplexYieldMetric>()),
        ("duplex.family_sizes", "FamilySizeMetric", columns_of::<FamilySizeMetric>()),
        ("duplex.umi_counts", "UmiMetric", columns_of::<UmiMetric>()),
        ("filter.stats", "FilterStatsMetrics", columns_of::<FilterStatsMetrics>()),
        ("group.family_sizes", "FamilySizeMetrics", columns_of::<FamilySizeMetrics>()),
        ("group.grouping_metrics", "UmiGroupingMetrics", columns_of::<UmiGroupingMetrics>()),
        (
            "group.position_group_sizes",
            "PositionGroupSizeMetrics",
            columns_of::<PositionGroupSizeMetrics>(),
        ),
        ("retag.metrics", "RetagMetric", columns_via_delim(RetagMetric::default())),
        ("review.details", "ConsensusVariantReviewInfo", review_columns),
        (
            "simplex.family_sizes",
            "SimplexFamilySizeMetric",
            columns_of::<SimplexFamilySizeMetric>(),
        ),
        ("simplex.simplex_yield_metrics", "SimplexYieldMetric", columns_of::<SimplexYieldMetric>()),
        ("simplex.umi_counts", "UmiMetric", columns_of::<UmiMetric>()),
    ]
}

fn expected_manifest() -> BTreeMap<String, Vec<String>> {
    emitted_metric_files().into_iter().map(|(key, _, columns)| (key.to_string(), columns)).collect()
}

// ---------------------------------------------------------------------------
// Source scanning
// ---------------------------------------------------------------------------

/// A named field of a scanned struct.
#[derive(Debug, PartialEq)]
struct ScannedField {
    name: String,
    /// The field is `f64` or `Option<f64>`.
    is_float: bool,
    /// A `#[serde(with = "...::float")]` attribute is present.
    has_float_attr: bool,
    /// `#[serde(skip)]` / `#[serde(skip_serializing)]`: never written, so never checked.
    is_skipped: bool,
}

/// A struct that derives serde's `Serialize`.
#[derive(Debug)]
struct ScannedStruct {
    name: String,
    fields: Vec<ScannedField>,
}

/// The identifier-like words of an attribute's tokens, e.g.
/// `#[serde(skip_serializing_if = "f")]` yields `serde`, `skip_serializing_if`, `f`.
fn attr_words(attr: &syn::Attribute) -> Vec<String> {
    quote::quote!(#attr)
        .to_string()
        .split(|c: char| !(c.is_alphanumeric() || c == '_'))
        .filter(|w| !w.is_empty())
        .map(str::to_string)
        .collect()
}

/// True for `#[derive(.. Serialize ..)]` and `#[cfg_attr(.., derive(.. Serialize ..))]`.
fn derives_serialize(attrs: &[syn::Attribute]) -> bool {
    attrs.iter().any(|attr| {
        (attr.path().is_ident("derive") || attr.path().is_ident("cfg_attr")) && {
            let words = attr_words(attr);
            words.iter().any(|w| w == "derive") && words.iter().any(|w| w == "Serialize")
        }
    })
}

fn is_cfg_test(attrs: &[syn::Attribute]) -> bool {
    attrs
        .iter()
        .any(|attr| attr.path().is_ident("cfg") && attr_words(attr).iter().any(|w| w == "test"))
}

fn serde_attrs(attrs: &[syn::Attribute]) -> impl Iterator<Item = &syn::Attribute> {
    attrs.iter().filter(|attr| attr.path().is_ident("serde"))
}

/// True for `f64` and `Option<f64>`.
fn is_float_type(ty: &syn::Type) -> bool {
    let syn::Type::Path(type_path) = ty else { return false };
    if type_path.path.is_ident("f64") {
        return true;
    }
    let Some(last) = type_path.path.segments.last() else { return false };
    let syn::PathArguments::AngleBracketed(args) = &last.arguments else { return false };
    last.ident == "Option"
        && args
            .args
            .iter()
            .any(|arg| matches!(arg, syn::GenericArgument::Type(t) if is_float_type(t)))
}

fn scan_struct(item: &syn::ItemStruct) -> ScannedStruct {
    let fields = match &item.fields {
        syn::Fields::Named(named) => named
            .named
            .iter()
            .map(|field| ScannedField {
                name: field.ident.as_ref().map(ToString::to_string).unwrap_or_default(),
                is_float: is_float_type(&field.ty),
                has_float_attr: serde_attrs(&field.attrs)
                    .any(|attr| quote::quote!(#attr).to_string().contains("::float\"")),
                is_skipped: serde_attrs(&field.attrs).any(|attr| {
                    attr_words(attr).iter().any(|w| w == "skip" || w == "skip_serializing")
                }),
            })
            .collect(),
        _ => Vec::new(),
    };
    ScannedStruct { name: item.ident.to_string(), fields }
}

fn collect_serialized_structs(items: &[syn::Item], out: &mut Vec<ScannedStruct>) {
    for item in items {
        match item {
            syn::Item::Struct(item) if derives_serialize(&item.attrs) => {
                out.push(scan_struct(item));
            }
            syn::Item::Mod(module) if !is_cfg_test(&module.attrs) => {
                if let Some((_, items)) = &module.content {
                    collect_serialized_structs(items, out);
                }
            }
            _ => {}
        }
    }
}

/// Every struct in `src` that derives serde's `Serialize`, in source order, skipping
/// `#[cfg(test)]` modules.
fn scan_serialized_structs(src: &str) -> Vec<ScannedStruct> {
    let file = syn::parse_file(src).expect("parse Rust source");
    let mut out = Vec::new();
    collect_serialized_structs(&file.items, &mut out);
    out
}

/// `<struct>.<field>` for every written `f64` / `Option<f64>` field of a serialized
/// struct in `src` that lacks the fgbio float encoding.
fn f64_fields_missing_float(src: &str) -> Vec<String> {
    scan_serialized_structs(src)
        .into_iter()
        .flat_map(|s| {
            s.fields
                .into_iter()
                .filter(|f| f.is_float && !f.has_float_attr && !f.is_skipped)
                .map(move |f| format!("{}.{}", s.name, f.name))
        })
        .collect()
}

/// All `.rs` files under the root crate's `src/` and every workspace crate's `src/`,
/// excluding the `xtask` developer tool (its fixtures are not fgumi outputs).
fn workspace_sources() -> Vec<PathBuf> {
    fn collect(dir: &Path, out: &mut Vec<PathBuf>) {
        for entry in std::fs::read_dir(dir).expect("read dir") {
            let path = entry.expect("dir entry").path();
            if path.is_dir() {
                collect(&path, out);
            } else if path.extension().is_some_and(|e| e == "rs") {
                out.push(path);
            }
        }
    }
    let root = Path::new(env!("CARGO_MANIFEST_DIR"));
    let mut files = Vec::new();
    collect(&root.join("src"), &mut files);
    for entry in std::fs::read_dir(root.join("crates")).expect("read crates") {
        let crate_dir = entry.expect("crate entry").path();
        if crate_dir.file_name().is_some_and(|n| n == "xtask") {
            continue;
        }
        let src = crate_dir.join("src");
        if src.is_dir() {
            collect(&src, &mut files);
        }
    }
    files
}

fn read(path: &Path) -> String {
    std::fs::read_to_string(path).expect("read source file")
}

// ---------------------------------------------------------------------------
// Tests
// ---------------------------------------------------------------------------

#[test]
fn manifest_keys_are_unique() {
    let mut seen = BTreeSet::new();
    let duplicates: Vec<&str> = emitted_metric_files()
        .iter()
        .map(|(key, _, _)| *key)
        .filter(|key| !seen.insert(*key))
        .collect();
    assert!(duplicates.is_empty(), "duplicate manifest keys: {duplicates:?}");
}

#[test]
fn committed_manifest_matches_structs() {
    let expected = expected_manifest();
    let text = std::fs::read_to_string(MANIFEST).expect("read metric_columns.json");
    let committed: BTreeMap<String, Vec<String>> =
        serde_json::from_str(&text).expect("parse metric_columns.json");
    assert_eq!(
        committed,
        expected,
        "metric_columns.json is stale; regenerate it with:\n{}",
        serde_json::to_string_pretty(&expected).expect("serialize manifest"),
    );
}

#[test]
fn every_serialized_struct_is_in_manifest_or_allowlisted() {
    let covered: BTreeSet<&str> = emitted_metric_files().iter().map(|(_, name, _)| *name).collect();
    let allowlisted: BTreeSet<&str> = NOT_EMITTED_AS_A_FILE.iter().map(|(name, _)| *name).collect();
    let unaccounted: Vec<String> = workspace_sources()
        .iter()
        .flat_map(|path| {
            scan_serialized_structs(&read(path))
                .into_iter()
                .filter(|s| {
                    !covered.contains(s.name.as_str()) && !allowlisted.contains(s.name.as_str())
                })
                .map(|s| format!("{} ({})", s.name, path.display()))
                .collect::<Vec<_>>()
        })
        .collect();
    assert!(
        unaccounted.is_empty(),
        "serialized structs missing from metric_columns.json (add to emitted_metric_files, \
         or to NOT_EMITTED_AS_A_FILE with a reason):\n{}",
        unaccounted.join("\n")
    );
}

#[test]
fn scanner_flags_f64_fields_without_float_encoding() {
    let src = r#"
#[derive(Debug, Serialize, Deserialize)]
pub struct ProbeMetric {
    /// Encoded.
    #[serde(with = "crate::float")]
    pub good: f64,
    /// Not encoded.
    pub bad: f64,
    pub maybe: Option<f64>,
    /// Never serialized, so never checked.
    #[serde(skip)]
    pub skipped: f64,
    pub count: u64,
}

#[derive(Debug, Clone)]
pub struct NotSerializedMetric {
    pub ignored: f64,
}
"#;
    assert_eq!(f64_fields_missing_float(src), vec!["ProbeMetric.bad", "ProbeMetric.maybe"]);
}

#[test]
fn scanner_sees_every_serialized_struct_layout() {
    let src = r#"
#[derive(
    Debug,
    Clone,
    Default,
    PartialEq,
    serde::Serialize,
    serde::Deserialize,
)]
pub(crate) struct WrappedDerive {
    pub(crate) rate: f64,
}

#[cfg_attr(feature = "x", derive(serde::Serialize))]
pub struct CfgAttrDerive {
    rate: f64,
}

#[derive(Serialize)]
pub struct NotNamedLikeAMetric {
    #[serde(skip_serializing_if = "Option::is_none")]
    pub maybe: Option<u64>,
    #[serde(skip_deserializing)]
    pub written: f64,
    #[serde(skip)]
    pub never_written: f64,
}

#[derive(Debug, Clone)]
pub struct SerializeState {
    pub rate: f64,
}

#[cfg(test)]
mod tests {
    #[derive(Serialize)]
    struct TestOnly {
        rate: f64,
    }
}
"#;
    let names: Vec<String> = scan_serialized_structs(src).into_iter().map(|s| s.name).collect();
    assert_eq!(names, vec!["WrappedDerive", "CfgAttrDerive", "NotNamedLikeAMetric"]);
    assert_eq!(
        f64_fields_missing_float(src),
        vec!["WrappedDerive.rate", "CfgAttrDerive.rate", "NotNamedLikeAMetric.written"]
    );
}

#[test]
fn every_serialized_f64_field_uses_fgbio_float_encoding() {
    let violations: Vec<String> = workspace_sources()
        .iter()
        .flat_map(|path| {
            f64_fields_missing_float(&read(path))
                .into_iter()
                .map(|v| format!("{v} ({})", path.display()))
                .collect::<Vec<_>>()
        })
        .collect();
    assert!(
        violations.is_empty(),
        "serialized f64 fields must carry #[serde(with = \"crate::float\")] so non-finite \
         values serialize as fgbio's Infinity/NaN (an Option<f64> field needs an \
         option-aware float helper):\n{}",
        violations.join("\n")
    );
}
