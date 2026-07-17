//! Compare two TSV metrics files by row key-join, with optional float precision
//! tolerance.
//!
//! fgbio and fgumi may legitimately emit metrics rows in different **order** or
//! (per fgbio-parity bugs tracked in fgumi#498: dense-vs-sparse family sizes,
//! missing zero-seeded rows, nondeterministic `duplex_umi_counts` ties) different
//! **sets**. A naive positional line-by-line diff cascades on the first such
//! divergence: every subsequent line looks different even though only one row
//! actually disagrees. This module instead outer-joins the two files' data rows on
//! a key column set (default: the first column; override with `--key-columns` for
//! multi-column keys like `duplex_family_sizes`'s `ab_size,ba_size`), so a
//! set/order difference is *localized* to the specific key(s) involved instead of
//! poisoning the entire comparison — while still failing hard (`DIFFER`) on any
//! unmatched or mismatched row, since those set/order differences are themselves
//! real fgbio-parity bugs the comparator must catch.
//!
//! The only tolerated divergence is float **representation** (fgbio's
//! `DecimalFormat("0.######")` vs fgumi's full precision): floats are compared
//! numerically within `--rel-tol`/`--abs-tol` (optionally pre-rounded via
//! `--precision`); integers and strings are compared exactly. Column *names* must
//! match exactly between the two files (a column-set difference is itself a
//! `DIFFER` — e.g. fgumi#498's `UmiGroupingMetrics` 12-vs-5 columns).
//!
//! **Key columns are also numeric-representation-normalized before joining** (see
//! [`canonicalize_key_field`]): a key value that parses as a number is canonicalized
//! (int/float parse, float optionally rounded to `--precision`) before being used to
//! join rows, so e.g. `duplex_yield_metrics`'s `fraction` key of `1.0` in one file and
//! `1` in the other still join as the same row instead of producing a spurious
//! presence-only `DIFFER`. Non-numeric keys are unaffected and still join by exact
//! string equality.

use crate::commands::command::Command;
use crate::commands::common::parse_bool;
use crate::commands::compare::engines::push_diff;
use crate::logging::OperationTimer;
use crate::validation::validate_file_exists;
use ahash::AHashMap;
use anyhow::{Context, Result, anyhow, bail, ensure};
use clap::Parser;
use log::info;
use std::cmp::Ordering;
use std::collections::BTreeMap;
use std::collections::BTreeSet;
use std::fs::File;
use std::io::{BufRead, BufReader};
use std::path::{Path, PathBuf};

/// Compare two TSV metrics files by row key-join.
///
/// Rows are matched between the two files by a key column set (default: the
/// first column), not by line position, so rows may legitimately appear in
/// different order between fgbio and fgumi without cascading a diff. Any row
/// whose key is present in only one file, or whose non-key values disagree, is
/// reported individually.
#[derive(Debug, Parser)]
#[command(
    name = "metrics",
    about = "Compare two TSV metrics files",
    long_about = r#"
Compare two TSV metrics files by row key-join.

The comparator parses both files' header (column names) and data rows, verifies
the two files declare the same column set, then outer-joins the data rows on a
key column set:
- Integers are compared exactly
- Floats are optionally rounded to a specified precision, then compared with
  relative/absolute tolerance (the only tolerated divergence -- this accounts
  for fgbio's `DecimalFormat("0.######")` vs fgumi's full-precision output)
- Strings are compared exactly
- A row whose key is present in only one file is reported as a localized
  difference naming that key, instead of cascading into every following line
- Rows may appear in different order between the two files: the key-join
  matches them regardless of position

By default the join key is the first column (works for metrics files keyed by
`umi`, `family_size`, `fraction`, or `key`). Pass `--key-columns` for a
multi-column key, e.g. `duplex_family_sizes`'s `ab_size,ba_size`.

Example usage:
  fgumi compare metrics file1.txt file2.txt
  fgumi compare metrics file1.txt file2.txt --precision 6
  fgumi compare metrics file1.txt file2.txt --rel-tol 1e-6 --abs-tol 1e-9
  fgumi compare metrics file1.txt file2.txt --key-columns ab_size,ba_size
  fgumi compare metrics file1.txt file2.txt --quiet
"#
)]
pub struct CompareMetrics {
    /// First TSV file
    #[arg(index = 1)]
    pub file1: PathBuf,

    /// Second TSV file
    #[arg(index = 2)]
    pub file2: PathBuf,

    /// Round floats to this many decimal places (-1 to disable rounding)
    #[arg(long = "precision", default_value = "6")]
    pub precision: i32,

    /// Relative tolerance for float comparison
    #[arg(long = "rel-tol", default_value = "1e-9")]
    pub rel_tol: f64,

    /// Absolute tolerance for float comparison
    #[arg(long = "abs-tol", default_value = "1e-9")]
    pub abs_tol: f64,

    /// Column(s) to join rows on, comma-separated: each entry is either a column
    /// name (from the header) or a 0-based column index. Defaults to the first
    /// column. Use a multi-entry key for metrics keyed by more than one column,
    /// e.g. `duplex_family_sizes`'s `--key-columns ab_size,ba_size`.
    #[arg(long = "key-columns", value_delimiter = ',')]
    pub key_columns: Option<Vec<String>>,

    /// Maximum number of differences to display
    #[arg(short = 'm', long = "max-diffs", default_value = "20")]
    pub max_diffs: usize,

    /// Quiet mode - the exit code is the only result signal (0=equal, 1=different):
    /// suppresses the stdout report and this command's own informational logging. Global
    /// startup logging stays under `RUST_LOG` control (like fgbio's `--log-level`).
    #[arg(short = 'q', long = "quiet", default_value = "false", num_args = 0..=1, default_missing_value = "true", action = clap::ArgAction::Set, value_parser = parse_bool)]
    pub quiet: bool,

    /// Verbose mode - print success message when files match
    #[arg(short = 'v', long = "verbose", default_value = "false", num_args = 0..=1, default_missing_value = "true", action = clap::ArgAction::Set, value_parser = parse_bool)]
    pub verbose: bool,
}

/// Parsed value from a TSV field.
#[derive(Debug, Clone)]
enum ParsedValue {
    Int(i64),
    Float(f64),
    String(String),
}

impl std::fmt::Display for ParsedValue {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        match self {
            ParsedValue::Int(i) => write!(f, "{i}"),
            ParsedValue::Float(fl) => write!(f, "{fl}"),
            ParsedValue::String(s) => write!(f, "{s}"),
        }
    }
}

/// Parse a string value as int, float (optionally rounded), or keep as string.
fn parse_value(value: &str, precision: Option<u32>) -> ParsedValue {
    // Try to parse as integer first
    if let Ok(i) = value.parse::<i64>() {
        return ParsedValue::Int(i);
    }

    // Try to parse as float
    if let Ok(f) = value.parse::<f64>() {
        let f = if let Some(p) = precision {
            let multiplier = 10_f64.powi(p as i32);
            (f * multiplier).round() / multiplier
        } else {
            f
        };
        return ParsedValue::Float(f);
    }

    // Keep as string
    ParsedValue::String(value.to_string())
}

/// Compare two values, using tolerance for floats.
fn values_equal(val1: &ParsedValue, val2: &ParsedValue, rel_tol: f64, abs_tol: f64) -> bool {
    match (val1, val2) {
        (ParsedValue::Int(i1), ParsedValue::Int(i2)) => i1 == i2,
        (ParsedValue::Float(f1), ParsedValue::Float(f2)) => {
            // Treat NaN == NaN (unlike IEEE 754 semantics where NaN != NaN)
            if f1.is_nan() && f2.is_nan() {
                return true;
            }
            // Handle infinities: +inf == +inf, -inf == -inf, but +inf != -inf
            if f1.is_infinite() && f2.is_infinite() {
                return f1.signum() == f2.signum();
            }
            // If only one is infinite, they are not equal
            if f1.is_infinite() || f2.is_infinite() {
                return false;
            }
            let diff = (f1 - f2).abs();
            let max_val = f1.abs().max(f2.abs());
            diff <= (rel_tol * max_val).max(abs_tol)
        }
        // Try to compare as floats if one is int and one is float
        (ParsedValue::Int(i), ParsedValue::Float(f))
        | (ParsedValue::Float(f), ParsedValue::Int(i)) => {
            let f1 = *i as f64;
            let f2 = *f;
            let diff = (f1 - f2).abs();
            let max_val = f1.abs().max(f2.abs());
            diff <= (rel_tol * max_val).max(abs_tol)
        }
        (ParsedValue::String(s1), ParsedValue::String(s2)) => s1 == s2,
        _ => false,
    }
}

/// A parsed TSV file: header column names (in on-disk order) plus data rows.
///
/// Blank lines are skipped (including a trailing newline's empty final "line");
/// the first non-blank line is the header, everything after is data.
#[derive(Debug)]
struct Tsv {
    /// Column names, in on-disk order.
    columns: Vec<String>,
    /// Column name -> index into `columns` (and thus into each row of `rows`).
    index_by_name: AHashMap<String, usize>,
    /// Data rows; each row has exactly `columns.len()` fields.
    rows: Vec<Vec<String>>,
}

impl Tsv {
    /// Parse `path` as a header + data-rows TSV.
    ///
    /// # Errors
    ///
    /// Returns an error if `path` cannot be read, or if a data row's field count
    /// does not match the header's (a malformed file, not a cross-file
    /// difference -- those are reported by [`compare_metrics`] instead).
    fn parse(path: &Path) -> Result<Self> {
        let file = File::open(path).with_context(|| format!("opening {}", path.display()))?;
        let mut lines = BufReader::new(file).lines();

        let mut columns: Vec<String> = Vec::new();
        for line in lines.by_ref() {
            let line = line.with_context(|| format!("reading {}", path.display()))?;
            let trimmed = line.trim_end_matches(['\n', '\r']);
            if trimmed.is_empty() {
                continue;
            }
            columns = trimmed.split('\t').map(str::to_owned).collect();
            break;
        }

        let index_by_name: AHashMap<String, usize> =
            columns.iter().enumerate().map(|(i, name)| (name.clone(), i)).collect();
        // A duplicate column name would collapse in `index_by_name` (last wins) and in the
        // downstream `BTreeSet` column-set check, silently dropping one of the same-named
        // columns from every row comparison — a potential false MATCH. Reject it up front.
        ensure!(
            index_by_name.len() == columns.len(),
            "{}: header has duplicate column name(s) (columns: {:?})",
            path.display(),
            columns
        );

        let mut rows = Vec::new();
        for line in lines {
            let line = line.with_context(|| format!("reading {}", path.display()))?;
            let trimmed = line.trim_end_matches(['\n', '\r']);
            if trimmed.is_empty() {
                continue;
            }
            let fields: Vec<String> = trimmed.split('\t').map(str::to_owned).collect();
            ensure!(
                fields.len() == columns.len(),
                "{}: data row has {} field(s), expected {} (header: {:?})",
                path.display(),
                fields.len(),
                columns.len(),
                columns
            );
            rows.push(fields);
        }

        Ok(Self { columns, index_by_name, rows })
    }
}

/// One `--key-columns` token, not yet resolved against a specific file's header.
#[derive(Debug, Clone)]
enum KeyColumnSpec {
    /// A literal column name.
    Name(String),
    /// A 0-based column index, resolved positionally against the *reference*
    /// header (always file1's) to obtain a canonical name -- see
    /// [`resolve_key_names`].
    Index(usize),
}

/// Parse each `--key-columns` token as a column index (if it parses as a plain
/// `usize`) or a literal column name otherwise.
fn parse_key_column_specs(tokens: &[String]) -> Vec<KeyColumnSpec> {
    tokens
        .iter()
        .map(|token| {
            let token = token.trim();
            match token.parse::<usize>() {
                Ok(index) => KeyColumnSpec::Index(index),
                Err(_) => KeyColumnSpec::Name(token.to_owned()),
            }
        })
        .collect()
}

/// Resolve `specs` to canonical column names using `reference_columns` (always
/// file1's header, so an index is unambiguous): a `Name` spec is used verbatim
/// (after checking it exists), an `Index` spec is resolved positionally.
///
/// The resulting names are then looked up **by name** in each file's own header
/// independently (see [`key_indices`]), so a key column resolves correctly even
/// if the two files declare their columns in a different order.
///
/// # Errors
///
/// Returns an error if a named key column does not exist in `reference_columns`,
/// or an index is out of range.
fn resolve_key_names(specs: &[KeyColumnSpec], reference_columns: &[String]) -> Result<Vec<String>> {
    specs
        .iter()
        .map(|spec| match spec {
            KeyColumnSpec::Name(name) => {
                ensure!(
                    reference_columns.iter().any(|c| c == name),
                    "key column {name:?} not found in header: {reference_columns:?}"
                );
                Ok(name.clone())
            }
            KeyColumnSpec::Index(index) => reference_columns.get(*index).cloned().ok_or_else(|| {
                anyhow!(
                    "key column index {index} out of range (header has {} column(s): {reference_columns:?})",
                    reference_columns.len()
                )
            }),
        })
        .collect()
}

/// Look up `key_names`' positions within `tsv`'s own header (by name -- see
/// [`resolve_key_names`] for why this is not simply reused across files).
///
/// # Errors
///
/// Returns an error if a key name is not present in `tsv`'s header. This should
/// not happen for the shared key names computed from the column-set check in
/// [`compare_metrics`], but is checked defensively.
fn key_indices(tsv: &Tsv, key_names: &[String]) -> Result<Vec<usize>> {
    key_names
        .iter()
        .map(|name| {
            tsv.index_by_name.get(name).copied().ok_or_else(|| {
                anyhow!("key column {name:?} not present in header: {:?}", tsv.columns)
            })
        })
        .collect()
}

/// Canonicalize a single key-column field for the join.
///
/// Parses `raw` the same way [`parse_value`] parses a data value (int, then float
/// optionally rounded to `precision`, else left as a string) and formats the result back
/// via [`ParsedValue`]'s `Display`. This collapses any purely-representational numeric
/// difference -- `"1"` vs `"1.0"`, `"0.500"` vs `"0.5"`, a leading `"+"` sign, or a
/// `--precision`-losing rounding difference -- to an identical join key (e.g. seen in
/// `duplex_yield_metrics`'s `fraction` key column: fgbio emits `1.0`, fgumi emits `1`). A
/// non-numeric value is returned verbatim, so non-numeric keys still join by exact string
/// equality -- normalization only ever *merges* keys that already denote the same number,
/// it never changes which rows are genuinely distinct.
fn canonicalize_key_field(raw: &str, precision: Option<u32>) -> String {
    parse_value(raw, precision).to_string()
}

/// Group `tsv`'s data rows by the (numeric-normalized, see
/// [`canonicalize_key_field`]) value of `key_idxs` -- exactly one row per key.
///
/// A `BTreeMap` (rather than a hash map) so the outer join in [`compare_metrics`]
/// can walk both files' groups as two sorted streams, merge-join style -- the same
/// discipline `engines::sort_verify::sort_verify_compare` uses for BAM records.
///
/// # Errors
///
/// Returns an error if the key column(s) do not uniquely identify rows (a key value
/// appears on more than one row). Every well-formed metric row is uniquely identified
/// by its dimension column(s) -- each fgbio/fgumi metric is emitted from a map or
/// histogram keyed by exactly those dimensions -- so a duplicate key means the chosen
/// key set is under-specified, not that the rows legitimately collide. Silently
/// reconciling such rows would let the oracle mask a real difference; instead we fail
/// loud and direct the caller to name the row's full identity via `--key-columns`.
fn group_by_key(
    tsv: &Tsv,
    key_idxs: &[usize],
    precision: Option<u32>,
    file: &Path,
) -> Result<BTreeMap<Vec<String>, usize>> {
    let mut groups: BTreeMap<Vec<String>, usize> = BTreeMap::new();
    for (row_idx, row) in tsv.rows.iter().enumerate() {
        let key: Vec<String> =
            key_idxs.iter().map(|&i| canonicalize_key_field(&row[i], precision)).collect();
        if groups.insert(key.clone(), row_idx).is_some() {
            bail!(
                "{}: key {} appears on multiple rows -- the chosen key column(s) do not \
                 uniquely identify rows; pass --key-columns to name the row's full identity",
                file.display(),
                format_key(&key),
            );
        }
    }
    Ok(groups)
}

/// Compare one matched `(file1 row, file2 row)` pair's non-key columns.
///
/// Columns are aligned **by name** (via each file's own `index_by_name`), not by
/// position: the two files may declare columns in different orders as long as
/// they share the same column set (already asserted by [`compare_metrics`]
/// before this is ever called).
fn diff_row_values(
    tsv1: &Tsv,
    row1_idx: usize,
    tsv2: &Tsv,
    row2_idx: usize,
    non_key_columns: &[String],
    cfg: &MetricsCompareConfig,
) -> Vec<String> {
    let row1 = &tsv1.rows[row1_idx];
    let row2 = &tsv2.rows[row2_idx];
    non_key_columns
        .iter()
        .filter_map(|name| {
            // Column-set equality (checked before grouping) guarantees these indices exist.
            let raw1 = &row1[tsv1.index_by_name[name]];
            let raw2 = &row2[tsv2.index_by_name[name]];
            let val1 = parse_value(raw1, cfg.precision);
            let val2 = parse_value(raw2, cfg.precision);
            (!values_equal(&val1, &val2, cfg.rel_tol, cfg.abs_tol)).then(|| {
                format!("column {name:?}: {raw1:?} != {raw2:?} (parsed: {val1} != {val2})")
            })
        })
        .collect()
}

/// Format a (possibly multi-column) key for a diagnostic message.
fn format_key(key: &[String]) -> String {
    if key.len() == 1 { key[0].clone() } else { format!("({})", key.join(", ")) }
}

/// Configuration for [`compare_metrics`].
pub struct MetricsCompareConfig {
    /// Decimal places to round floats to before comparing (`None` disables
    /// rounding).
    pub precision: Option<u32>,
    /// Relative tolerance for float comparison.
    pub rel_tol: f64,
    /// Absolute tolerance for float comparison.
    pub abs_tol: f64,
    /// Column(s) (name or 0-based index, see [`KeyColumnSpec`]) to join rows on.
    /// `None` defaults to the first file's first column.
    pub key_columns: Option<Vec<String>>,
    /// Maximum number of entries collected in
    /// [`MetricsCompareOutcome::diff_details`].
    pub max_diffs: usize,
}

/// Outcome of a [`compare_metrics`] run.
#[derive(Debug, Default)]
pub struct MetricsCompareOutcome {
    /// Number of data rows read from file1.
    pub rows1: usize,
    /// Number of data rows read from file2.
    pub rows2: usize,
    /// Number of keys present in both files (regardless of whether their values
    /// matched).
    pub matched_keys: usize,
    /// Number of keys present only in file1.
    pub only_in_file1: usize,
    /// Number of keys present only in file2.
    pub only_in_file2: usize,
    /// Number of matched keys whose rows disagreed on at least one non-key
    /// column beyond tolerance.
    pub value_diffs: usize,
    /// `true` if the two files declared different column *sets* (a hard
    /// `DIFFER`; row-level comparison is not attempted since column alignment
    /// would be undefined).
    pub column_set_mismatch: bool,
    /// Human-readable diff strings (column-set, presence, and value mismatches),
    /// capped at the caller-supplied `max_diffs`.
    pub diff_details: Vec<String>,
}

impl MetricsCompareOutcome {
    /// Returns `true` iff the two files matched: same column set, no key present
    /// on only one side, and every matched key's row(s) agreed (within float
    /// tolerance).
    #[must_use]
    pub fn is_match(&self) -> bool {
        !self.column_set_mismatch
            && self.only_in_file1 == 0
            && self.only_in_file2 == 0
            && self.value_diffs == 0
    }
}

/// Compare two metrics TSV files by outer-joining their data rows on a key
/// column set (see this module's doc comment).
///
/// # Errors
///
/// Returns an error if either file cannot be parsed as a TSV, or if
/// `cfg.key_columns` names/indexes a column that does not exist.
pub fn compare_metrics(
    file1: &Path,
    file2: &Path,
    cfg: &MetricsCompareConfig,
) -> Result<MetricsCompareOutcome> {
    // Reject nonsensical tolerances up front: a negative tolerance makes two identical finite
    // floats "differ" (`diff <= negative` is never true), and a non-finite (NaN/±inf)
    // tolerance makes arbitrary values "match" (`diff <= inf` is always true), either way
    // silently corrupting the comparison. Fail loud instead.
    ensure!(
        cfg.rel_tol.is_finite() && cfg.rel_tol >= 0.0,
        "relative tolerance must be finite and non-negative (got {})",
        cfg.rel_tol
    );
    ensure!(
        cfg.abs_tol.is_finite() && cfg.abs_tol >= 0.0,
        "absolute tolerance must be finite and non-negative (got {})",
        cfg.abs_tol
    );

    let tsv1 = Tsv::parse(file1)?;
    let tsv2 = Tsv::parse(file2)?;

    // Resolve any explicit `--key-columns` up front, *before* the column-set and
    // empty-file early returns below. Otherwise `--key-columns <missing>` against
    // two empty files (or files whose column sets differ) would bypass key
    // resolution entirely and spuriously report MATCH — the user named a key that
    // does not exist and must get the resolution error regardless of file
    // content. An explicitly empty key list is itself rejected. The default
    // (first-column) key only makes sense once both files are known to share a
    // non-empty column set, so it stays deferred until after those checks.
    let explicit_key_names = match &cfg.key_columns {
        Some(tokens) => {
            ensure!(!tokens.is_empty(), "at least one key column must be specified");
            Some(resolve_key_names(&parse_key_column_specs(tokens), &tsv1.columns)?)
        }
        None => None,
    };

    let mut outcome = MetricsCompareOutcome {
        rows1: tsv1.rows.len(),
        rows2: tsv2.rows.len(),
        ..Default::default()
    };

    // Column *set* must match exactly -- a set difference (e.g. fgumi#498's
    // UmiGroupingMetrics 12-vs-5 columns) is itself a DIFFER, and column
    // alignment for row comparison is undefined without it.
    let set1: BTreeSet<&str> = tsv1.columns.iter().map(String::as_str).collect();
    let set2: BTreeSet<&str> = tsv2.columns.iter().map(String::as_str).collect();
    if set1 != set2 {
        outcome.column_set_mismatch = true;
        let only1: Vec<&str> = set1.difference(&set2).copied().collect();
        let only2: Vec<&str> = set2.difference(&set1).copied().collect();
        if !only1.is_empty() {
            push_diff(&mut outcome.diff_details, cfg.max_diffs, || {
                format!("columns only in file1: {only1:?}")
            });
        }
        if !only2.is_empty() {
            push_diff(&mut outcome.diff_details, cfg.max_diffs, || {
                format!("columns only in file2: {only2:?}")
            });
        }
        return Ok(outcome);
    }

    // Both files have the same (possibly empty) column set. An empty header
    // means an empty file with no data rows possible either -- trivially a
    // match, and there is no column to default the join key to.
    if tsv1.columns.is_empty() {
        return Ok(outcome);
    }

    // Reached only once both files share a non-empty column set (checks above),
    // so `tsv1.columns[0]` is a valid default when no explicit key was given.
    let key_names = explicit_key_names.unwrap_or_else(|| vec![tsv1.columns[0].clone()]);
    let key_idx1 = key_indices(&tsv1, &key_names)?;
    let key_idx2 = key_indices(&tsv2, &key_names)?;

    let non_key_columns: Vec<String> =
        tsv1.columns.iter().filter(|c| !key_names.contains(c)).cloned().collect();

    let groups1 = group_by_key(&tsv1, &key_idx1, cfg.precision, file1)?;
    let groups2 = group_by_key(&tsv2, &key_idx2, cfg.precision, file2)?;

    // Merge-join the two sorted key streams (mirrors
    // `engines::sort_verify::sort_verify_compare`'s no-resync discipline): always compare
    // the two current keys and advance whichever side is behind (or both, on a
    // match), so a presence-only key is reported once and never causes the join
    // to desync.
    let mut iter1 = groups1.into_iter();
    let mut iter2 = groups2.into_iter();
    let mut cur1 = iter1.next();
    let mut cur2 = iter2.next();

    loop {
        match (cur1.take(), cur2.take()) {
            (None, None) => break,
            (Some((key, _row)), None) => {
                outcome.only_in_file1 += 1;
                push_diff(&mut outcome.diff_details, cfg.max_diffs, || {
                    format!("row {} only in file1", format_key(&key))
                });
                cur1 = iter1.next();
            }
            (None, Some((key, _row))) => {
                outcome.only_in_file2 += 1;
                push_diff(&mut outcome.diff_details, cfg.max_diffs, || {
                    format!("row {} only in file2", format_key(&key))
                });
                cur2 = iter2.next();
            }
            (Some((k1, r1)), Some((k2, r2))) => match k1.cmp(&k2) {
                Ordering::Less => {
                    outcome.only_in_file1 += 1;
                    push_diff(&mut outcome.diff_details, cfg.max_diffs, || {
                        format!("row {} only in file1", format_key(&k1))
                    });
                    cur1 = iter1.next();
                    cur2 = Some((k2, r2));
                }
                Ordering::Greater => {
                    outcome.only_in_file2 += 1;
                    push_diff(&mut outcome.diff_details, cfg.max_diffs, || {
                        format!("row {} only in file2", format_key(&k2))
                    });
                    cur1 = Some((k1, r1));
                    cur2 = iter2.next();
                }
                Ordering::Equal => {
                    outcome.matched_keys += 1;
                    let diffs = diff_row_values(&tsv1, r1, &tsv2, r2, &non_key_columns, cfg);
                    if !diffs.is_empty() {
                        outcome.value_diffs += 1;
                        push_diff(&mut outcome.diff_details, cfg.max_diffs, || {
                            format!("row {}: {}", format_key(&k1), diffs.join("; "))
                        });
                    }
                    cur1 = iter1.next();
                    cur2 = iter2.next();
                }
            },
        }
    }

    Ok(outcome)
}

impl Command for CompareMetrics {
    fn execute(&self, _command_line: &str) -> Result<()> {
        validate_file_exists(&self.file1, "First file")?;
        validate_file_exists(&self.file2, "Second file")?;

        // `--quiet` means "only the exit code communicates the comparison result": suppress
        // not just the stdout report but also this command's own informational stderr logging
        // (the timer line and the identical/differ `info!`s). The process-wide startup banner
        // is emitted by `main` before dispatch and stays under `RUST_LOG` control — the same
        // orthogonal-logging model as fgbio's global `--log-level`.
        let timer = (!self.quiet).then(|| OperationTimer::new("Comparing metrics"));

        let precision = if self.precision >= 0 { Some(self.precision as u32) } else { None };
        let cfg = MetricsCompareConfig {
            precision,
            rel_tol: self.rel_tol,
            abs_tol: self.abs_tol,
            key_columns: self.key_columns.clone(),
            max_diffs: self.max_diffs,
        };

        let outcome = compare_metrics(&self.file1, &self.file2, &cfg)?;
        let is_equal = outcome.is_match();

        if !self.quiet {
            if !is_equal {
                if outcome.column_set_mismatch {
                    println!("Column set mismatch between files:");
                    for diff in &outcome.diff_details {
                        println!("  {diff}");
                    }
                } else {
                    let total = outcome.only_in_file1 + outcome.only_in_file2 + outcome.value_diffs;
                    println!(
                        "Found {total} difference(s) between files ({} row key(s) matched):",
                        outcome.matched_keys
                    );
                    for diff in &outcome.diff_details {
                        println!("  {diff}");
                    }
                    if total > outcome.diff_details.len() {
                        println!("  ... and {} more", total - outcome.diff_details.len());
                    }
                }
                println!("RESULT: metrics files DIFFER");
            } else if self.verbose {
                println!(
                    "RESULT: metrics files are IDENTICAL ({} row key(s) matched)",
                    outcome.matched_keys
                );
            }
        }

        if is_equal {
            if let Some(timer) = timer {
                info!("Metrics files are identical");
                timer.log_completion(outcome.matched_keys as u64);
            }
            Ok(())
        } else {
            if !self.quiet {
                info!("Metrics files differ");
            }
            Err(super::CompareMismatch("metrics files differ".to_owned()).into())
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use rstest::rstest;
    use std::io::Write as _;
    use tempfile::tempdir;

    /// Writes `contents` (already TSV-formatted, tab-separated) to `dir/name`.
    fn write_tsv(dir: &Path, name: &str, contents: &str) -> PathBuf {
        let path = dir.join(name);
        let mut f = File::create(&path).expect("create temp tsv");
        f.write_all(contents.as_bytes()).expect("write temp tsv");
        path
    }

    fn default_cfg() -> MetricsCompareConfig {
        MetricsCompareConfig {
            precision: Some(6),
            rel_tol: 1e-9,
            abs_tol: 1e-9,
            key_columns: None,
            max_diffs: 20,
        }
    }

    /// A negative or non-finite tolerance is rejected at the `compare_metrics` boundary
    /// before any comparison runs — a negative tolerance would make identical finite floats
    /// "differ", and a non-finite one would make arbitrary values "match".
    #[rstest]
    #[case::rel_negative(-1e-9, 1e-9, "relative tolerance")]
    #[case::abs_negative(1e-9, -1e-9, "absolute tolerance")]
    #[case::rel_nan(f64::NAN, 1e-9, "relative tolerance")]
    #[case::abs_inf(1e-9, f64::INFINITY, "absolute tolerance")]
    #[case::rel_neg_inf(f64::NEG_INFINITY, 1e-9, "relative tolerance")]
    fn compare_metrics_rejects_invalid_tolerance(
        #[case] rel_tol: f64,
        #[case] abs_tol: f64,
        #[case] expected_msg: &str,
    ) {
        let dir = tempdir().expect("tempdir");
        let f1 = write_tsv(dir.path(), "a.txt", "key\tvalue\nA\t1\n");
        let f2 = write_tsv(dir.path(), "b.txt", "key\tvalue\nA\t1\n");
        let cfg = MetricsCompareConfig { rel_tol, abs_tol, ..default_cfg() };
        let err = compare_metrics(&f1, &f2, &cfg)
            .expect_err("an invalid tolerance must be a hard error, not a silent comparison");
        assert!(
            err.to_string().contains(expected_msg),
            "error must name the offending tolerance, got: {err}"
        );
    }

    /// An explicit `--key-columns` naming a column that does not exist must be rejected even
    /// when both files are empty (no columns): the column-set and empty-file early returns
    /// must not bypass key resolution and let `--key-columns missing` spuriously report MATCH.
    #[test]
    fn compare_metrics_rejects_missing_explicit_key_even_for_empty_files() {
        let dir = tempdir().expect("tempdir");
        let f1 = write_tsv(dir.path(), "a.txt", "");
        let f2 = write_tsv(dir.path(), "b.txt", "");
        let cfg = MetricsCompareConfig {
            key_columns: Some(vec!["missing".to_string()]),
            ..default_cfg()
        };
        let err = compare_metrics(&f1, &f2, &cfg).expect_err(
            "an explicit key column that does not exist must error, even for empty files",
        );
        assert!(
            err.to_string().contains("missing"),
            "error must name the missing key column, got: {err}"
        );
    }

    /// An explicitly empty `--key-columns` list is meaningless and must be rejected up front
    /// rather than silently falling back to the default first-column key.
    #[test]
    fn compare_metrics_rejects_empty_key_column_list() {
        let dir = tempdir().expect("tempdir");
        let f1 = write_tsv(dir.path(), "a.txt", "key\tvalue\nA\t1\n");
        let f2 = write_tsv(dir.path(), "b.txt", "key\tvalue\nA\t1\n");
        let cfg = MetricsCompareConfig { key_columns: Some(vec![]), ..default_cfg() };
        let err = compare_metrics(&f1, &f2, &cfg)
            .expect_err("an explicitly empty key-column list must be rejected");
        assert!(
            err.to_string().contains("at least one key column"),
            "error must explain the empty key list, got: {err}"
        );
    }

    /// The preserved happy path: two empty files with *no* explicit key trivially match —
    /// there are no rows to compare and no column to default the join key to.
    #[test]
    fn compare_metrics_empty_files_without_explicit_key_match() {
        let dir = tempdir().expect("tempdir");
        let f1 = write_tsv(dir.path(), "a.txt", "");
        let f2 = write_tsv(dir.path(), "b.txt", "");
        let outcome = compare_metrics(&f1, &f2, &default_cfg())
            .expect("two empty files with no explicit key must compare cleanly");
        assert!(outcome.is_match(), "two empty files trivially match: {outcome:?}");
    }

    #[test]
    fn test_nan_equality() {
        // NaN values should be considered equal
        let nan1 = ParsedValue::Float(f64::NAN);
        let nan2 = ParsedValue::Float(f64::NAN);
        assert!(values_equal(&nan1, &nan2, 1e-9, 1e-9));
    }

    #[test]
    fn test_positive_infinity_equality() {
        // +Infinity values should be considered equal
        let inf1 = ParsedValue::Float(f64::INFINITY);
        let inf2 = ParsedValue::Float(f64::INFINITY);
        assert!(values_equal(&inf1, &inf2, 1e-9, 1e-9));
    }

    #[test]
    fn test_negative_infinity_equality() {
        // -Infinity values should be considered equal
        let neg_inf1 = ParsedValue::Float(f64::NEG_INFINITY);
        let neg_inf2 = ParsedValue::Float(f64::NEG_INFINITY);
        assert!(values_equal(&neg_inf1, &neg_inf2, 1e-9, 1e-9));
    }

    #[test]
    fn test_mixed_infinity_not_equal() {
        // +Infinity and -Infinity should NOT be equal
        let pos_inf = ParsedValue::Float(f64::INFINITY);
        let neg_inf = ParsedValue::Float(f64::NEG_INFINITY);
        assert!(!values_equal(&pos_inf, &neg_inf, 1e-9, 1e-9));
    }

    #[test]
    fn test_nan_not_equal_to_number() {
        // NaN and a regular number should NOT be equal
        let nan = ParsedValue::Float(f64::NAN);
        let num = ParsedValue::Float(1.0);
        assert!(!values_equal(&nan, &num, 1e-9, 1e-9));
    }

    #[test]
    fn test_parse_nan_string() {
        // "NaN" string should parse to NaN float
        let parsed = parse_value("NaN", Some(6));
        let ParsedValue::Float(f) = parsed else {
            unreachable!("Expected Float variant for NaN");
        };
        assert!(f.is_nan());
    }

    #[test]
    fn test_parse_infinity_string() {
        // "Infinity" and "inf" should parse to infinity
        let parsed = parse_value("Infinity", Some(6));
        let ParsedValue::Float(f) = parsed else {
            unreachable!("Expected Float variant for Infinity");
        };
        assert!(f.is_infinite() && f.is_sign_positive());

        let parsed = parse_value("inf", Some(6));
        let ParsedValue::Float(f) = parsed else {
            unreachable!("Expected Float variant for inf");
        };
        assert!(f.is_infinite() && f.is_sign_positive());
    }

    // ---- key-join engine -------------------------------------------------

    #[test]
    fn identical_files_match() {
        let dir = tempdir().expect("tempdir");
        let content = "umi\tcount\nAAA\t3\nCCC\t5\n";
        let f1 = write_tsv(dir.path(), "a.txt", content);
        let f2 = write_tsv(dir.path(), "b.txt", content);
        let outcome = compare_metrics(&f1, &f2, &default_cfg()).expect("compare should succeed");
        assert!(outcome.is_match(), "{:?}", outcome.diff_details);
        assert_eq!(outcome.matched_keys, 2);
    }

    #[test]
    fn same_rows_different_order_match() {
        // The fix under test: a positional diff would cascade here, but the
        // key-join must tolerate row reordering.
        let dir = tempdir().expect("tempdir");
        let f1 = write_tsv(dir.path(), "a.txt", "umi\tcount\nAAA\t3\nCCC\t5\nGGG\t1\n");
        let f2 = write_tsv(dir.path(), "b.txt", "umi\tcount\nGGG\t1\nCCC\t5\nAAA\t3\n");
        let outcome = compare_metrics(&f1, &f2, &default_cfg()).expect("compare should succeed");
        assert!(
            outcome.is_match(),
            "row order must not affect the key-join: {:?}",
            outcome.diff_details
        );
        assert_eq!(outcome.matched_keys, 3);
    }

    #[test]
    fn row_only_in_file1_is_localized_differ() {
        // Simulates fgumi#498's dense-vs-sparse family_size rows: file2 is
        // missing the zero-seeded/sparse family_size=3 row entirely.
        let dir = tempdir().expect("tempdir");
        let f1 = write_tsv(dir.path(), "a.txt", "family_size\tcount\n1\t10\n2\t5\n3\t1\n");
        let f2 = write_tsv(dir.path(), "b.txt", "family_size\tcount\n1\t10\n2\t5\n");
        let outcome = compare_metrics(&f1, &f2, &default_cfg()).expect("compare should succeed");
        assert!(!outcome.is_match());
        assert_eq!(outcome.only_in_file1, 1);
        assert_eq!(outcome.only_in_file2, 0);
        // Unaffected rows must NOT show up as diffs -- this is the localization
        // guarantee the key-join exists to provide.
        assert_eq!(outcome.value_diffs, 0);
        assert_eq!(outcome.matched_keys, 2);
        assert!(
            outcome.diff_details.iter().any(|d| d.contains("3 only in file1")),
            "{:?}",
            outcome.diff_details
        );
    }

    #[test]
    fn row_only_in_file2_is_localized_differ() {
        let dir = tempdir().expect("tempdir");
        let f1 = write_tsv(dir.path(), "a.txt", "family_size\tcount\n1\t10\n");
        let f2 = write_tsv(dir.path(), "b.txt", "family_size\tcount\n1\t10\n2\t5\n");
        let outcome = compare_metrics(&f1, &f2, &default_cfg()).expect("compare should succeed");
        assert!(!outcome.is_match());
        assert_eq!(outcome.only_in_file2, 1);
        assert!(
            outcome.diff_details.iter().any(|d| d.contains("2 only in file2")),
            "{:?}",
            outcome.diff_details
        );
    }

    #[test]
    fn value_beyond_epsilon_differs() {
        let dir = tempdir().expect("tempdir");
        let f1 = write_tsv(dir.path(), "a.txt", "key\tvalue\nA\t1.0\n");
        let f2 = write_tsv(dir.path(), "b.txt", "key\tvalue\nA\t1.1\n");
        let outcome = compare_metrics(&f1, &f2, &default_cfg()).expect("compare should succeed");
        assert!(!outcome.is_match());
        assert_eq!(outcome.value_diffs, 1);
    }

    #[test]
    fn float_representation_only_matches_within_epsilon() {
        // The one accepted divergence: fgbio's rounded `DecimalFormat` output vs
        // fgumi's full-precision float must still MATCH.
        let dir = tempdir().expect("tempdir");
        let f1 = write_tsv(dir.path(), "a.txt", "key\tvalue\nA\t0.333333\n");
        let f2 = write_tsv(dir.path(), "b.txt", "key\tvalue\nA\t0.3333333333333333\n");
        let outcome = compare_metrics(&f1, &f2, &default_cfg()).expect("compare should succeed");
        assert!(outcome.is_match(), "{:?}", outcome.diff_details);
    }

    #[test]
    fn different_column_sets_differ() {
        let dir = tempdir().expect("tempdir");
        let f1 = write_tsv(dir.path(), "a.txt", "umi\tcount\tfraction\nAAA\t3\t0.5\n");
        let f2 = write_tsv(dir.path(), "b.txt", "umi\tcount\nAAA\t3\n");
        let outcome = compare_metrics(&f1, &f2, &default_cfg()).expect("compare should succeed");
        assert!(!outcome.is_match());
        assert!(outcome.column_set_mismatch);
        assert!(
            outcome.diff_details.iter().any(|d| d.contains("fraction")),
            "{:?}",
            outcome.diff_details
        );
    }

    #[test]
    fn multi_column_key_join() {
        // Mirrors duplex_family_sizes, keyed on (ab_size, ba_size).
        let dir = tempdir().expect("tempdir");
        let f1 = write_tsv(
            dir.path(),
            "a.txt",
            "ab_size\tba_size\tcount\n1\t0\t100\n1\t1\t50\n2\t0\t20\n",
        );
        // Row order differs and each file has one key the other lacks, but the
        // two shared (ab_size, ba_size) keys must still join and compare equal.
        let f2 = write_tsv(
            dir.path(),
            "b.txt",
            "ab_size\tba_size\tcount\n1\t1\t50\n1\t0\t100\n3\t0\t5\n",
        );
        let cfg = MetricsCompareConfig {
            key_columns: Some(vec!["ab_size".to_owned(), "ba_size".to_owned()]),
            ..default_cfg()
        };
        let outcome = compare_metrics(&f1, &f2, &cfg).expect("compare should succeed");
        assert!(!outcome.is_match());
        assert_eq!(outcome.matched_keys, 2, "{:?}", outcome.diff_details);
        assert_eq!(outcome.only_in_file1, 1, "(2, 0) only in file1");
        assert_eq!(outcome.only_in_file2, 1, "(3, 0) only in file2");
        assert_eq!(outcome.value_diffs, 0, "{:?}", outcome.diff_details);
    }

    #[test]
    fn key_columns_accepts_indices() {
        let dir = tempdir().expect("tempdir");
        let f1 = write_tsv(dir.path(), "a.txt", "ab_size\tba_size\tcount\n1\t0\t100\n2\t1\t9\n");
        let f2 = write_tsv(dir.path(), "b.txt", "ab_size\tba_size\tcount\n2\t1\t9\n1\t0\t100\n");
        let cfg = MetricsCompareConfig {
            key_columns: Some(vec!["0".to_owned(), "1".to_owned()]),
            ..default_cfg()
        };
        let outcome = compare_metrics(&f1, &f2, &cfg).expect("compare should succeed");
        assert!(outcome.is_match(), "{:?}", outcome.diff_details);
        assert_eq!(outcome.matched_keys, 2);
    }

    /// A key that does not uniquely identify rows (here `key=A` twice, differing only
    /// in a non-key column) is an under-specified key, not a legitimate collision:
    /// every fgbio/fgumi metric is emitted from a map/histogram keyed by its
    /// dimensions, so a duplicate key can only mean the wrong key was chosen. The
    /// oracle must fail loud (never silently reconcile, which could mask a real diff)
    /// and point at `--key-columns`.
    #[test]
    fn duplicate_key_under_current_key_set_is_rejected() {
        let dir = tempdir().expect("tempdir");
        let f1 = write_tsv(dir.path(), "a.txt", "key\tvalue\nA\t1\nA\t2\n");
        let f2 = write_tsv(dir.path(), "b.txt", "key\tvalue\nA\t2\nA\t1\n");
        let err = compare_metrics(&f1, &f2, &default_cfg())
            .expect_err("a duplicate key must be rejected, not multiset-matched");
        let msg = err.to_string();
        assert!(
            msg.contains("appears on multiple rows") && msg.contains("--key-columns"),
            "error must name the duplicate key and direct to --key-columns: {msg}"
        );
    }

    /// The remedy the error directs to: naming the row's full identity via
    /// `--key-columns` makes the *same* data compare cleanly, because each row is now
    /// uniquely keyed. (Order still differs between the files; the join handles it.)
    #[test]
    fn unique_key_columns_resolve_the_duplicate() {
        let dir = tempdir().expect("tempdir");
        let f1 = write_tsv(dir.path(), "a.txt", "key\tvalue\nA\t1\nA\t2\n");
        let f2 = write_tsv(dir.path(), "b.txt", "key\tvalue\nA\t2\nA\t1\n");
        let cfg = MetricsCompareConfig {
            key_columns: Some(vec!["key".to_owned(), "value".to_owned()]),
            ..default_cfg()
        };
        let outcome = compare_metrics(&f1, &f2, &cfg).expect("compare should succeed");
        assert!(outcome.is_match(), "{:?}", outcome.diff_details);
        assert_eq!(outcome.matched_keys, 2);
    }

    #[rstest]
    // fgbio/fgumi may format the same logical key differently (`1` vs `1.0` for a
    // `fraction`-style key column, as observed in `duplex_yield_metrics`) -- the join
    // must normalize numeric representation and still treat these as the same row.
    #[case::int_vs_float_dot_zero("1", "1.0")]
    #[case::trailing_zeros("0.500", "0.5")]
    #[case::leading_plus_and_int("+2", "2")]
    fn numeric_key_join_normalizes_representation(#[case] key1: &str, #[case] key2: &str) {
        let dir = tempdir().expect("tempdir");
        let f1 = write_tsv(dir.path(), "a.txt", &format!("fraction\tcount\n{key1}\t100\n"));
        let f2 = write_tsv(dir.path(), "b.txt", &format!("fraction\tcount\n{key2}\t100\n"));
        let outcome = compare_metrics(&f1, &f2, &default_cfg()).expect("compare should succeed");
        assert!(
            outcome.is_match(),
            "{key1:?} and {key2:?} must join as the same numeric key: {:?}",
            outcome.diff_details
        );
        assert_eq!(outcome.matched_keys, 1);
        assert_eq!(outcome.only_in_file1, 0);
        assert_eq!(outcome.only_in_file2, 0);
    }

    #[test]
    fn numeric_key_join_still_separates_genuinely_different_keys() {
        // A real key difference (not just a representation difference) must NOT join --
        // normalization must not paper over an actual presence-only divergence.
        let dir = tempdir().expect("tempdir");
        let f1 = write_tsv(dir.path(), "a.txt", "fraction\tcount\n1\t100\n");
        let f2 = write_tsv(dir.path(), "b.txt", "fraction\tcount\n2\t100\n");
        let outcome = compare_metrics(&f1, &f2, &default_cfg()).expect("compare should succeed");
        assert!(!outcome.is_match());
        assert_eq!(outcome.matched_keys, 0);
        assert_eq!(outcome.only_in_file1, 1);
        assert_eq!(outcome.only_in_file2, 1);
    }

    #[test]
    fn numeric_key_join_normalizes_multi_column_key() {
        // The normalization must apply per-key-column, not just to a single-column key.
        let dir = tempdir().expect("tempdir");
        let f1 =
            write_tsv(dir.path(), "a.txt", "ab_size\tba_size\tcount\n1\t0.0\t100\n2.0\t1\t50\n");
        let f2 =
            write_tsv(dir.path(), "b.txt", "ab_size\tba_size\tcount\n1.0\t0\t100\n2\t1.0\t50\n");
        let cfg = MetricsCompareConfig {
            key_columns: Some(vec!["ab_size".to_owned(), "ba_size".to_owned()]),
            ..default_cfg()
        };
        let outcome = compare_metrics(&f1, &f2, &cfg).expect("compare should succeed");
        assert!(outcome.is_match(), "{:?}", outcome.diff_details);
        assert_eq!(outcome.matched_keys, 2);
    }

    #[rstest]
    #[case::exact_match("A\t1\n", "A\t1\n", true)]
    #[case::value_differs("A\t1\n", "A\t2\n", false)]
    #[case::int_vs_float_within_tolerance("A\t1\n", "A\t1.0000000001\n", true)]
    fn simple_value_cases(#[case] row1: &str, #[case] row2: &str, #[case] expect_match: bool) {
        let dir = tempdir().expect("tempdir");
        let f1 = write_tsv(dir.path(), "a.txt", &format!("key\tvalue\n{row1}"));
        let f2 = write_tsv(dir.path(), "b.txt", &format!("key\tvalue\n{row2}"));
        let outcome = compare_metrics(&f1, &f2, &default_cfg()).expect("compare should succeed");
        assert_eq!(outcome.is_match(), expect_match, "{:?}", outcome.diff_details);
    }

    #[test]
    fn empty_files_match() {
        let dir = tempdir().expect("tempdir");
        let f1 = write_tsv(dir.path(), "a.txt", "");
        let f2 = write_tsv(dir.path(), "b.txt", "");
        let outcome = compare_metrics(&f1, &f2, &default_cfg()).expect("compare should succeed");
        assert!(outcome.is_match(), "{:?}", outcome.diff_details);
    }

    #[test]
    fn max_diffs_caps_diff_details_but_not_counts() {
        let dir = tempdir().expect("tempdir");
        use std::fmt::Write as _;
        let mut c1 = String::from("key\tvalue\n");
        let mut c2 = String::from("key\tvalue\n");
        for i in 0..10 {
            writeln!(c1, "k{i}\t{i}").expect("write to String never fails");
            writeln!(c2, "k{i}\t{}", i + 100).expect("write to String never fails");
        }
        let f1 = write_tsv(dir.path(), "a.txt", &c1);
        let f2 = write_tsv(dir.path(), "b.txt", &c2);
        let cfg = MetricsCompareConfig { max_diffs: 3, ..default_cfg() };
        let outcome = compare_metrics(&f1, &f2, &cfg).expect("compare should succeed");
        assert_eq!(outcome.value_diffs, 10);
        assert_eq!(outcome.diff_details.len(), 3);
    }
}
