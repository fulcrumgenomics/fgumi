//! Compare two BAM files for equality of core SAM fields and tag values.
//!
//! Tags are compared by value regardless of their order in the record.
//! This is useful for verifying that two BAM files are functionally equivalent
//! even if they were produced by different tools.
//!
//! # Performance
//!
//! This module is optimized for large BAM files with:
//! - Multi-threaded BGZF decompression (`--threads`)
//! - Batch processing with double buffering
//! - Lazy tag comparison (avoids allocations when tags match)
//! - Zero-copy field comparisons where possible
//! - Fast hashing with ahash

use crate::logging::OperationTimer;
use crate::sam::SamTag;
use crate::validation::validate_file_exists;
use ahash::{AHashMap, AHashSet};
use anyhow::{Result, anyhow, bail};
use clap::{Parser, ValueEnum};
use crossbeam_channel::{Receiver, bounded};
use fgumi_bam_io::ProgressTracker;
use fgumi_bam_io::{RawBamReaderAuto, create_raw_bam_reader};
use fgumi_raw_bam::fields as raw_fields;
use fgumi_raw_bam::{RawRecord, find_int_tag, find_string_tag};
use itertools::Itertools;
use log::info;
use noodles::sam::Header;
use noodles::sam::alignment::record::data::field::Tag;
use noodles::sam::alignment::record_buf::RecordBuf;
use noodles::sam::alignment::record_buf::data::field::Value;
use noodles::sam::alignment::record_buf::data::field::value::Array;
use rayon::prelude::*;
use std::collections::BTreeMap;
use std::path::PathBuf;
use std::thread;

use crate::commands::command::Command;
use crate::commands::common::{MemoryReserve, parse_bool, resolve_memory_budget};
use crate::commands::sort::TMP_DIRS_ENV;

use super::engines::content::ContentPredicate;
use super::engines::keyjoin::{self, KeyJoinConfig};
use super::engines::positional::positional_compare;
use super::raw_compare::{raw_compare_structured, raw_records_byte_equal};
use super::record_key::{self, RecordKey};

/// Comparison mode for BAM files
#[derive(Debug, Clone, Copy, Default, ValueEnum)]
pub enum CompareMode {
    /// Full comparison: verify MI groupings AND compare all BAM content
    /// Both files must be in the same order (e.g., query-name sorted).
    /// This combines grouping equivalence check with full content comparison.
    #[default]
    Full,
    /// Content comparison: all fields and tags must match exactly
    /// This is a pure record-by-record comparison without MI grouping analysis.
    Content,
    /// Grouping comparison: verify MI groupings are equivalent (for grouped BAMs)
    /// Both files must be in the same order (e.g., query-name sorted).
    /// Validates read names and R1/R2 flags match, then verifies that reads
    /// with the same MI in one file have the same MI in the other.
    /// Does NOT compare other BAM content (sequence, quality, other tags).
    Grouping,
}

/// Preset comparison settings for a specific fgumi pipeline stage.
///
/// Each variant encodes canonical `--mode` and `--ignore-order` defaults for
/// comparing BAM output from that stage, including cases (e.g. cross-tool
/// `group` comparison against fgbio) where MI values or record order may
/// legitimately differ. Explicit `--mode` or `--ignore-order` flags override
/// the preset.
#[derive(Debug, Clone, Copy, ValueEnum)]
pub enum CommandPreset {
    /// Extract output: no MI tags; exact content comparison.
    Extract,
    /// Zipper output: preserves MI tags unchanged; exact content comparison.
    Zipper,
    /// Sort output: order *is* the payload, so this bypasses `--mode`/`ContentPredicate`
    /// entirely and routes to a dedicated engine
    /// ([`super::engines::sort_verify::sort_verify_compare`]): the shared sort order is
    /// detected from both inputs' `@HD` header, each file is verified to be independently
    /// correctly ordered, and the two files are compared as a multiset grouped by maximal
    /// equal-core-sort-key run — tolerating intra-run reordering (coordinate ties, and the
    /// documented template-coordinate name-hash-vs-lexical `SORT-01` residue) while still
    /// catching a genuine mis-sort, a missing/extra record, or any content difference. An
    /// explicit `--mode` or `--ignore-order` alongside `--command sort` is rejected, since
    /// neither concept applies to this engine.
    Sort,
    /// Correct output: modifies RX tag only; exact content comparison.
    Correct,
    /// Dedup output: deterministic; exact content comparison.
    Dedup,
    /// Group output: MI values and record order may differ between tools or
    /// runs. Compares via the key-join engine
    /// ([`super::engines::keyjoin::keyjoin_compare`]) under
    /// [`ContentPredicate::ExactMinusMi`] plus a separate fgumi-MI/fgbio-MI
    /// bijection check; see [`CommandPreset::resolve`] for the full rationale.
    Group,
    /// Simplex consensus output: positional, exact content comparison
    /// (see [`ContentPredicate::ExactConsensus`](super::engines::content::ContentPredicate::ExactConsensus)).
    Simplex,
    /// Duplex consensus output: positional, exact content comparison
    /// (see [`ContentPredicate::ExactConsensus`](super::engines::content::ContentPredicate::ExactConsensus)).
    Duplex,
    /// CODEC consensus output: positional, exact content comparison
    /// (see [`ContentPredicate::ExactConsensus`](super::engines::content::ContentPredicate::ExactConsensus)).
    Codec,
    /// Filter output: consensus reads with reads dropped, but the same depth tags
    /// (`cD`/`cM`/`cE`, and duplex `aD`/`aM`/`bD`/`bM`/`aE`/`bE`) as the consensus command
    /// that produced them, passed through unchanged. Compares via the same predicate as
    /// consensus output
    /// ([`ContentPredicate::ExactConsensus`](super::engines::content::ContentPredicate::ExactConsensus));
    /// see [`CommandPreset::resolve`] for the full rationale.
    Filter,
}

impl CommandPreset {
    /// Canonical `(mode, ignore_order, content_predicate)` resolution for this preset — the
    /// single source of truth for how a preset resolves to comparison behavior.
    /// `defaults()` and `content_predicate()` below are thin wrappers over this method, and
    /// `effective_settings()` and `execute()`'s predicate selection consume their results in
    /// turn: there is exactly one place to update when a preset's resolved behavior changes.
    ///
    /// Groups and rationale:
    ///
    /// - `Extract | Zipper | Sort | Correct | Dedup` → `(Content, false, Exact)`: plain
    ///   record-by-record comparison, no accepted divergence. (`Sort`'s mapping here is
    ///   vestigial: `CompareBams::execute` special-cases `CommandPreset::Sort` and returns
    ///   before ever consulting this method, routing to
    ///   [`super::engines::sort_verify::sort_verify_compare`] instead, which has no notion of
    ///   `CompareMode`/`ContentPredicate` at all. This arm exists only so the match stays
    ///   exhaustive without a wildcard.)
    /// - `Filter | Simplex | Duplex | Codec` → `(Content, false, ExactConsensus)`: all four
    ///   deal in consensus reads carrying the depth tags (`cD`/`cM`/`cE`, and duplex
    ///   `aD`/`aM`/`bD`/`bM`/`aE`/`bE`). `Simplex`/`Duplex`/`Codec` are the commands that
    ///   write those tags; `filter` only drops reads afterward, it never rewrites them. All
    ///   four compare positionally under
    ///   [`ContentPredicate::ExactConsensus`](super::engines::content::ContentPredicate::ExactConsensus),
    ///   which (now that fgumi clamps the depth tags to fgbio's `Short` ceiling) is an exact
    ///   content comparison — retained as a distinct preset predicate for these consensus
    ///   commands but behaviourally equivalent to `Exact`.
    /// - `Group` → `(Grouping, true, ExactMinusMi)`: MI values and record order may
    ///   legitimately differ between tools or runs, so `Group` is the only preset that
    ///   verifies grouping equivalence instead of routing to `execute_content`. It compares
    ///   via the key-join engine ([`super::engines::keyjoin::keyjoin_compare`]): records are
    ///   paired by [`RecordKey`](super::record_key::RecordKey) after canonicalizing both
    ///   inputs to queryname order, content is compared under
    ///   [`ContentPredicate::ExactMinusMi`](super::engines::content::ContentPredicate::ExactMinusMi)
    ///   (everything except the MI tag), and the fgumi-MI/fgbio-MI mapping observed across
    ///   matched pairs must separately be a consistent bijection — the predicate excludes MI
    ///   precisely because that bijection check, not the content predicate, is what verifies
    ///   MI equivalence.
    ///
    /// Exhaustive over every [`CommandPreset`] variant so that adding a new preset forces a
    /// conscious choice here rather than silently defaulting via a wildcard arm.
    fn resolve(self) -> (CompareMode, bool, ContentPredicate) {
        match self {
            Self::Extract | Self::Zipper | Self::Sort | Self::Correct | Self::Dedup => {
                (CompareMode::Content, false, ContentPredicate::Exact)
            }
            Self::Filter | Self::Simplex | Self::Duplex | Self::Codec => {
                (CompareMode::Content, false, ContentPredicate::ExactConsensus)
            }
            Self::Group => (CompareMode::Grouping, true, ContentPredicate::ExactMinusMi),
        }
    }

    /// Canonical `(mode, ignore_order)` defaults for this preset. See [`Self::resolve`] for
    /// the full resolution table and rationale.
    fn defaults(self) -> (CompareMode, bool) {
        let (mode, ignore_order, _) = self.resolve();
        (mode, ignore_order)
    }

    /// The [`ContentPredicate`] to use for this preset's content comparison. See
    /// [`Self::resolve`] for the full resolution table and rationale. Note that for `Group`
    /// (resolved mode `CompareMode::Grouping`), this value is *not* consulted by the
    /// key-join engine — [`engines::keyjoin::keyjoin_compare`] hardcodes
    /// [`ContentPredicate::ExactMinusMi`] internally and is not configurable via
    /// `--command`/`--mode`. `Group` never reaches `execute_content`, but this method is
    /// still live for it (its resolved value is exercised only by the
    /// `group_preset_content_predicate_is_exact_minus_mi` test below, documenting the
    /// intended predicate even though it is not threaded anywhere at runtime).
    fn content_predicate(self) -> ContentPredicate {
        self.resolve().2
    }
}

/// Compare two BAM files for equality.
///
/// Compares core SAM fields (QNAME, FLAG, RNAME, POS, MAPQ, CIGAR, RNEXT, PNEXT, TLEN, SEQ, QUAL)
/// and tag values. Tags are compared by value regardless of order.
#[derive(Debug, Parser)]
#[command(
    name = "bams",
    about = "Compare two BAM files for equality",
    long_about = r#"
Compare two BAM files for equality of core SAM fields and tag values.

This tool compares BAM files record-by-record, checking:
- Core SAM fields: QNAME, FLAG, RNAME, POS, MAPQ, CIGAR, RNEXT, PNEXT, TLEN, SEQ, QUAL
- Tag values (order-independent comparison)

Tags are compared by value only - the order of tags within a record does not matter.
This allows comparing BAM files produced by different tools that may serialize tags
in different orders.

MODES:

  full (default):
    Combines grouping equivalence check with full content comparison.
    Both files must be in the same order (e.g., query-name sorted).
    Verifies MI groupings are equivalent AND all other fields match.

  content:
    Pure record-by-record comparison of all fields and tags.
    Does not analyze MI groupings - just compares raw content.

  grouping:
    For comparing grouped BAM files where MI assignment order may differ.
    Without --ignore-order, both files MUST be in the same order (e.g.,
    query-name sorted with `fgumi sort --order queryname`); validates that
    read names/R1/R2 flags match and that reads sharing an MI in one file
    share an MI in the other, but does NOT compare other BAM content.
    With --ignore-order (the `group` preset's default), each input is
    internally canonicalized to queryname order and merge-joined by record
    identity: every matched pair is also compared under EXACT-MI (every
    field except the MI tag), so a non-MI content difference now DIFFERs
    even when the MI grouping itself is untouched, in addition to the MI
    bijection check.

COMMAND PRESETS (--command):

  Use `--command <stage>` to apply canonical `--mode` and `--ignore-order`
  defaults for comparing output from a specific fgumi pipeline stage. This is
  especially useful for cross-tool comparisons (e.g. fgumi vs. fgbio) where
  MI values or record order may legitimately differ. Explicit `--mode` or
  `--ignore-order` flags override the preset.

  Command         --mode      --ignore-order   Notes
  ─────────────────────────────────────────────────────────────────────────
  extract         content     false            No MI tags; deterministic
  zipper          content     false            Preserves MI tags unchanged
  sort            (dedicated sort-verify engine; --mode/--ignore-order rejected)
  correct         content     false            Modifies RX tag only, not MI
  dedup           content     false            Deterministic
  filter          content     false            Passes through MI/depth tags unchanged; exact (like consensus)
  group           grouping    true             Key-join: EXACT-MI content + MI bijection (cross-tool)
  simplex         content     false            Saturation-aware exact (cD/cM/cE carve-out)
  duplex          content     false            Saturation-aware exact (cD/cM/cE carve-out)
  codec           content     false            Saturation-aware exact (cD/cM/cE carve-out)

  `sort` verifies order instead of comparing content positionally: it detects the
  shared sort order from both inputs' @HD header, checks each file is itself
  correctly ordered, and compares the two files as a multiset grouped by maximal
  equal-core-sort-key run — tolerating intra-run tie reordering (coordinate ties,
  and fgumi's template-coordinate name-hash-vs-lexical tie residue) while still
  catching a mis-sort, a missing/extra record, or any content difference.

  Examples:

    # Preset equivalents of the above:
    fgumi compare bams --command extract a.bam b.bam
    fgumi compare bams --command sort     a.bam b.bam
    fgumi compare bams --command group    a.bam b.bam
    fgumi compare bams --command simplex  a.bam b.bam

    # Preset + explicit override (e.g. same-tool group comparison):
    fgumi compare bams --command group --mode full --ignore-order=false a.bam b.bam

Example usage:
  fgumi compare bams bam1.bam bam2.bam                    # full mode (default)
  fgumi compare bams bam1.bam bam2.bam --mode content    # content only
  fgumi compare bams bam1.bam bam2.bam --mode grouping   # grouping only
  fgumi compare bams bam1.bam bam2.bam --mode grouping --ignore-order  # consensus output
  fgumi compare bams bam1.bam bam2.bam --command simplex  # preset for simplex
  fgumi compare bams bam1.bam bam2.bam --max-diffs 20
"#
)]
pub struct CompareBams {
    /// First BAM file
    #[arg(index = 1)]
    pub bam1: PathBuf,

    /// Second BAM file
    #[arg(index = 2)]
    pub bam2: PathBuf,

    /// Use preset comparison settings for a specific fgumi pipeline stage.
    /// Sets `--mode` and `--ignore-order` to the canonical defaults for
    /// comparing BAM output from that stage (see COMMAND PRESETS above).
    /// Explicit `--mode` or `--ignore-order` flags override the preset.
    #[arg(long = "command", short = 'c')]
    pub command: Option<CommandPreset>,

    /// Comparison mode: 'full' (MI grouping + content, for group output),
    /// 'content' (all fields, for extract/zipper/sort/correct/dedup/filter output),
    /// 'grouping' (MI equivalence only, for simplex/duplex/codec output).
    /// Overrides `--command` preset if both given. Defaults to 'full' when
    /// neither `--mode` nor `--command` is set.
    #[arg(long = "mode")]
    pub mode: Option<CompareMode>,

    /// Maximum number of differences to report in detail
    #[arg(short = 'm', long = "max-diffs", default_value = "10")]
    pub max_diffs: usize,

    /// Quiet mode - only exit code indicates result (0=equal, 1=different)
    #[arg(short = 'q', long = "quiet", default_value = "false", num_args = 0..=1, default_missing_value = "true", action = clap::ArgAction::Set, value_parser = parse_bool)]
    pub quiet: bool,

    /// Ignore record order when comparing in grouping mode.
    /// Required for comparing output from consensus commands (simplex/duplex/codec)
    /// when run with --threads, as parallel processing causes non-deterministic ordering.
    /// Only valid with --mode grouping.
    /// Overrides `--command` preset if both given. Defaults to false when
    /// neither `--ignore-order` nor `--command` is set.
    #[arg(long = "ignore-order", num_args = 0..=1, default_missing_value = "true", action = clap::ArgAction::Set, value_parser = parse_bool)]
    pub ignore_order: Option<bool>,

    /// Initial buffer size for --ignore-order mode (number of records)
    #[arg(long = "buffer-size", default_value = "1000")]
    pub buffer_size: usize,

    /// Number of threads for BGZF decompression and parallel comparison.
    /// Using more threads significantly speeds up processing of large BAM files.
    #[arg(short = 't', long = "threads", default_value = "1")]
    pub threads: usize,

    /// Batch size for parallel processing (number of records per batch).
    /// Larger batches reduce synchronization overhead but use more memory.
    #[arg(long = "batch-size", default_value = "10000")]
    pub batch_size: usize,

    /// Total memory budget for the internal queryname-canonicalization sort used by
    /// `--command group`'s key-join engine. Ignored by every other mode/preset.
    /// This is a total budget, not a per-thread budget: it is not multiplied by
    /// `--threads`.
    #[arg(long = "sort-memory", default_value = "512M", value_parser = crate::commands::common::parse_memory)]
    pub sort_memory: crate::commands::common::MemoryLimit,

    /// Temporary directory for the internal queryname-canonicalization sort's spill
    /// files (`--command group` only). Repeatable, same semantics as `fgumi sort
    /// -T`/`--tmp-dir`. Falls back to `FGUMI_TMP_DIRS` (see `fgumi sort --help`),
    /// then a disk-backed default (never a bare system temp directory, which may be
    /// tmpfs on some hosts).
    #[arg(long = "sort-tmp-dir", action = clap::ArgAction::Append)]
    pub sort_tmp_dirs: Vec<PathBuf>,
}

/// Statistics from comparing two BAM files.
#[derive(Debug, Default)]
struct CompareStats {
    bam1_count: u64,
    bam2_count: u64,
    core_matches: u64,
    core_diffs: u64,
    tag_matches: u64,
    tag_diffs: u64,
    tag_order_diffs: u64,
    /// Number of paired records in BAM1 that are missing an MI tag (full mode only).
    missing_mi_bam1: u64,
    /// Number of paired records in BAM2 that are missing an MI tag (full mode only).
    missing_mi_bam2: u64,
    diff_details: Vec<DiffDetail>,
}

/// Details about a single difference found.
#[derive(Debug)]
struct DiffDetail {
    record_num: u64,
    qname: String,
    flags: String,
    diff_type: DiffType,
    diffs: Vec<String>,
}

#[derive(Debug)]
enum DiffType {
    CountMismatch,
    CoreDiff,
    TagDiff,
    ReadNameMismatch,
    FlagMismatch,
}

impl std::fmt::Display for DiffType {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        match self {
            DiffType::CountMismatch => write!(f, "count_mismatch"),
            DiffType::CoreDiff => write!(f, "core_diff"),
            DiffType::TagDiff => write!(f, "tag_diff"),
            DiffType::ReadNameMismatch => write!(f, "read_name_mismatch"),
            DiffType::FlagMismatch => write!(f, "flag_mismatch"),
        }
    }
}

/// Core SAM field names for reporting.
///
/// `pub(crate)` so the content-predicate engine (`engines::content`) can reuse it when
/// rendering per-field diff strings, avoiding a duplicate field-name list.
pub(crate) const FIELD_NAMES: [&str; 11] =
    ["QNAME", "FLAG", "RNAME", "POS", "MAPQ", "CIGAR", "RNEXT", "PNEXT", "TLEN", "SEQ", "QUAL"];

/// Format CIGAR from raw BAM record bytes as a SAM-style string.
fn format_cigar_raw(bam: &[u8]) -> String {
    let s = fgumi_raw_bam::cigar_to_string_from_raw(bam);
    if s.is_empty() { "*".to_string() } else { s }
}

/// Format sequence from raw BAM record bytes as an ASCII string.
fn format_sequence_raw(bam: &[u8]) -> String {
    let view = fgumi_raw_bam::RawRecordView::new(bam);
    if view.l_seq() == 0 {
        "*".to_string()
    } else {
        // extract_sequence returns uppercase ASCII bases (A/C/G/T/N/=…)
        String::from_utf8_lossy(&fgumi_raw_bam::extract_sequence(bam)).into_owned()
    }
}

/// Extract core SAM fields as comparable strings directly from raw BAM bytes.
///
/// RNAME and RNEXT are resolved through `header` (requires typed API).
/// All other fields are read from raw bytes via `RawRecordView`.
///
/// `pub(crate)` (and taking `&[u8]` rather than `&RawRecord`) so the content-predicate
/// engine (`engines::content`) can reuse it for diff-string rendering; `&RawRecord`
/// call sites still work unchanged via `Deref<Target = [u8]>` coercion.
pub(crate) fn get_core_fields_raw(raw: &[u8], header: &noodles::sam::Header) -> [String; 11] {
    let view = fgumi_raw_bam::RawRecordView::new(raw);

    let qname = String::from_utf8_lossy(view.read_name()).into_owned();
    let flag = view.flags().to_string();

    // ref_id is 0-based index into header reference sequences; -1 = unmapped
    let ref_id = view.ref_id();
    let rname = if ref_id < 0 {
        "*".to_string()
    } else {
        header
            .reference_sequences()
            .get_index(ref_id as usize)
            .map(|(name, _)| name.to_string())
            .unwrap_or_else(|| "*".to_string())
    };

    // pos is 0-based; SAM/display convention is 1-based (0 when unmapped)
    let pos_raw = view.pos();
    let pos = if pos_raw < 0 { "0".to_string() } else { (pos_raw + 1).to_string() };

    let mapq_raw = view.mapq();
    let mapq = if mapq_raw == 255 { "255".to_string() } else { mapq_raw.to_string() };

    let cigar = format_cigar_raw(raw);

    let mate_ref_id = view.mate_ref_id();
    let rnext = if mate_ref_id < 0 {
        "*".to_string()
    } else {
        header
            .reference_sequences()
            .get_index(mate_ref_id as usize)
            .map(|(name, _)| name.to_string())
            .unwrap_or_else(|| "*".to_string())
    };

    let mate_pos_raw = view.mate_pos();
    let pnext = if mate_pos_raw < 0 { "0".to_string() } else { (mate_pos_raw + 1).to_string() };

    let tlen = view.template_length().to_string();

    let seq = format_sequence_raw(raw);

    // Quality scores in raw BAM are 0-based Phred. Per SAM/BAM spec, absent QUAL is
    // signaled by ALL bytes being 0xFF — match that exactly to avoid mis-classifying
    // malformed records with a stray 0xFF. Saturate the Phred→ASCII conversion at '~'
    // (Phred 93) to avoid u8 wraparound for any remaining out-of-range bytes.
    let qual_bytes = fgumi_raw_bam::quality_scores_slice(raw);
    let qual = if qual_bytes.is_empty() || qual_bytes.iter().all(|&q| q == 0xFF) {
        "*".to_string()
    } else {
        qual_bytes.iter().map(|&q| (q.saturating_add(33).min(b'~')) as char).collect()
    };

    [qname, flag, rname, pos, mapq, cigar, rnext, pnext, tlen, seq, qual]
}

/// Format a tag as a two-character string.
fn format_tag(tag: Tag) -> String {
    let bytes: [u8; 2] = tag.into();
    String::from_utf8_lossy(&bytes).to_string()
}

/// Extract tags as a sorted map for order-independent comparison.
fn get_tags_map(record: &RecordBuf) -> BTreeMap<String, String> {
    record.data().iter().map(|(tag, value)| (format_tag(tag), format_tag_value(value))).collect()
}

/// Format a tag value for comparison/display.
fn format_tag_value(value: &Value) -> String {
    match value {
        Value::Character(c) => format!("{c}"),
        Value::Int8(i) => format!("{i}"),
        Value::UInt8(i) => format!("{i}"),
        Value::Int16(i) => format!("{i}"),
        Value::UInt16(i) => format!("{i}"),
        Value::Int32(i) => format!("{i}"),
        Value::UInt32(i) => format!("{i}"),
        Value::Float(f) => format!("{f}"),
        Value::String(s) => s.to_string(),
        Value::Hex(h) => format!("{h:?}"),
        Value::Array(arr) => format_array(arr),
    }
}

/// Format an array value for display, showing both type and values.
fn format_array(arr: &Array) -> String {
    match arr {
        Array::Int8(v) => format!("B:c,{}", v.iter().map(ToString::to_string).join(",")),
        Array::UInt8(v) => format!("B:C,{}", v.iter().map(ToString::to_string).join(",")),
        Array::Int16(v) => format!("B:s,{}", v.iter().map(ToString::to_string).join(",")),
        Array::UInt16(v) => format!("B:S,{}", v.iter().map(ToString::to_string).join(",")),
        Array::Int32(v) => format!("B:i,{}", v.iter().map(ToString::to_string).join(",")),
        Array::UInt32(v) => format!("B:I,{}", v.iter().map(ToString::to_string).join(",")),
        Array::Float(v) => format!("B:f,{}", v.iter().map(ToString::to_string).join(",")),
    }
}

/// Convert a record's name to a String, returning "*" if missing.
fn record_name_to_string(record: &RecordBuf) -> String {
    record.name().map_or_else(|| "*".to_string(), std::string::ToString::to_string)
}

/// Check if the first-in-template (R1) flag is set in raw BAM record bytes.
#[inline]
fn is_first_segment_raw(raw: &RawRecord) -> bool {
    fgumi_raw_bam::RawRecordView::new(raw.as_ref()).is_first_segment()
}

// ============================================================================
// Batch processing types
// ============================================================================

/// Result of comparing a single pair of records
#[derive(Debug)]
struct RecordCompareResult {
    core_match: bool,
    tags_match: bool,
    tag_order_match: bool,
    diff_detail: Option<DiffDetail>,
}

/// Result of comparing a single pair of records in grouping mode
#[derive(Debug)]
struct GroupingCompareResult {
    record_num: u64,
    record_key: RecordKey,
    /// Read name as String, only populated when needed for error reporting
    read_name_for_display: Option<String>,
    mi1: Option<MiKey>,
    mi2: Option<MiKey>,
    name_match: bool,
    flag_match: bool,
    diff_detail: Option<DiffDetail>,
}

// ============================================================================
// Types and MI map helpers (using ahash)
// ============================================================================

/// Key used to group records by molecular identifier during comparison.
///
/// Paired-UMI grouping strategies (fgumi and fgbio) emit MI as a Z-type string
/// encoded `<id>/<A|B>`, where the suffix distinguishes the two strand
/// orientations of the same double-stranded molecule. Treating a `PairedA(n)`
/// and a `PairedB(n)` as the same group would mask a real disagreement, so the
/// comparator must keep the suffix distinct from the integer id.
#[derive(Copy, Clone, Debug, PartialEq, Eq, Hash)]
pub(crate) enum MiKey {
    /// Plain integer MI (single-strand assigners, or an `MI:i:<int>` record).
    Int(i64),
    /// Paired-strand MI: `base` is the molecule id; `strand` is `b'A'` or `b'B'`.
    Strand { base: i64, strand: u8 },
}

impl MiKey {
    /// The molecule-level id, ignoring any duplex strand suffix.
    ///
    /// For duplex output the two strands of one molecule share a `base` but
    /// differ in `strand` (`/A` vs `/B`). Grouping-equivalence checks that key
    /// on the full [`MiKey`] treat the strands as independent groups, so they
    /// cannot detect a *strand-pairing* difference — one file bundling `X/A` +
    /// `X/B` into a single molecule while another splits them into `X/A` +
    /// `Y/A`. Comparing at the `base` level catches that: reads that share a
    /// molecule in one file must share a molecule in the other.
    pub(crate) fn base(&self) -> i64 {
        match self {
            MiKey::Int(v) => *v,
            MiKey::Strand { base, .. } => *base,
        }
    }
}

impl std::fmt::Display for MiKey {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        match self {
            MiKey::Int(v) => write!(f, "{v}"),
            MiKey::Strand { base, strand } => write!(f, "{base}/{}", *strand as char),
        }
    }
}

/// Count molecule-level (strand-suffix-stripped) grouping mismatches in both
/// directions: the number of `base` groups in either file whose reads map to
/// more than one `base` in the other file.
///
/// This complements the full-[`MiKey`] grouping check (which keeps `/A`/`/B`
/// distinct to catch strand *disagreements*) by catching strand-*pairing*
/// differences the full-key check is structurally blind to: if BAM1 pairs
/// `X/A` + `X/B` into molecule `X` and BAM2 splits those same reads into
/// `X/A` + `Y/A`, the full-key mapping is still a consistent bijection
/// (`X/A↔X/A`, `X/B↔Y/A`) with zero mismatches, yet the duplex molecules
/// genuinely differ. For single-strand (`Int`) MIs `base()` is the id itself,
/// so this reduces to the same partition check and never adds false mismatches.
fn count_base_pairing_mismatches(
    mi_map1: &AHashMap<RecordKey, MiKey>,
    mi_map2: &AHashMap<RecordKey, MiKey>,
) -> u64 {
    fn count_one_direction(
        from: &AHashMap<RecordKey, MiKey>,
        to: &AHashMap<RecordKey, MiKey>,
    ) -> u64 {
        let mut base_to_base: AHashMap<i64, (i64, bool)> = AHashMap::new();
        for (record_key, mi_from) in from {
            if let Some(mi_to) = to.get(record_key) {
                base_to_base
                    .entry(mi_from.base())
                    .and_modify(|(first, has_mismatch)| {
                        if *first != mi_to.base() {
                            *has_mismatch = true;
                        }
                    })
                    .or_insert((mi_to.base(), false));
            }
        }
        base_to_base.values().filter(|(_, has_mismatch)| *has_mismatch).count() as u64
    }

    count_one_direction(mi_map1, mi_map2) + count_one_direction(mi_map2, mi_map1)
}

/// Build a map from MI value to set of record identity keys.
fn build_mi_groups_compact(
    mi_map: &AHashMap<RecordKey, MiKey>,
) -> AHashMap<MiKey, AHashSet<RecordKey>> {
    let mut groups: AHashMap<MiKey, AHashSet<RecordKey>> = AHashMap::new();
    for (record_key, mi) in mi_map {
        groups.entry(*mi).or_default().insert(record_key.clone());
    }
    groups
}

/// Extract the MI tag from raw BAM record bytes.
///
/// Accepts three forms:
/// - Integer-typed MI (`MI:i:<int>`) → `MiKey::Int`.
/// - String-typed integer MI (`MI:Z:<int>`) → `MiKey::Int`.
/// - String-typed paired MI (`MI:Z:<int>/A` or `/B`) → `MiKey::Strand`.
///
/// Any other string payload (e.g. non-numeric prefix, unknown strand suffix)
/// yields `None`, matching the "missing MI" treatment used by the caller.
pub(crate) fn get_mi_tag_raw(raw: &RawRecord) -> Option<MiKey> {
    let aux = raw_fields::aux_data_slice(raw.as_ref());
    if let Some(v) = find_int_tag(aux, SamTag::MI) {
        return Some(MiKey::Int(v));
    }
    let bytes = find_string_tag(aux, SamTag::MI)?;
    let s = std::str::from_utf8(bytes).ok()?;
    if let Some((base_str, strand_str)) = s.rsplit_once('/') {
        let base = base_str.parse::<i64>().ok()?;
        let strand = match strand_str.as_bytes() {
            b"A" => b'A',
            b"B" => b'B',
            _ => return None,
        };
        Some(MiKey::Strand { base, strand })
    } else {
        s.parse::<i64>().ok().map(MiKey::Int)
    }
}

/// Statistics for grouping comparison mode.
#[derive(Debug, Default)]
struct GroupingStats {
    total_records: u64,
    order_mismatches: u64,
    flag_mismatches: u64,
    missing_mi_bam1: u64,
    missing_mi_bam2: u64,
    grouping_mismatches: u64,
    unique_groups_bam1: usize,
    unique_groups_bam2: usize,
    count_mismatch: bool,
}

// ============================================================================
// Double-buffered batch reading
// ============================================================================

/// Message type for the double-buffered raw reader channel.
///
/// `pub(crate)` so the positional engine (`engines::positional`) can drive the
/// same double-buffered reader as `execute_content`.
pub(crate) enum RawBatchMessage {
    /// A batch of raw records.
    Batch(Vec<RawRecord>),
    /// End of file reached.
    Eof,
    /// Error occurred during reading.
    Error(String),
}

/// Reads a single batch of raw records from a BAM reader.
fn read_raw_batch(
    reader: &mut RawBamReaderAuto,
    batch_size: usize,
) -> std::io::Result<(Vec<RawRecord>, bool)> {
    let mut batch = Vec::with_capacity(batch_size);
    let mut record = RawRecord::new();

    for _ in 0..batch_size {
        if reader.read_record(&mut record)? == 0 {
            return Ok((batch, true)); // EOF
        }
        batch.push(std::mem::take(&mut record));
    }
    Ok((batch, false))
}

/// Starts a background reader thread that sends raw record batches through a channel.
/// Returns a receiver for the batches and the BAM header.
///
/// `pub(crate)` so the positional engine (`engines::positional`) can reuse the same
/// double-buffered reader as `execute_content`.
pub(crate) fn start_raw_batch_reader(
    path: PathBuf,
    threads: usize,
    batch_size: usize,
) -> Result<(Receiver<RawBatchMessage>, Header)> {
    // Open the reader on the main thread to get the header
    let (mut reader, header) = create_raw_bam_reader(&path, threads)?;

    // Create a bounded channel (double buffering = 2 slots)
    let (tx, rx) = bounded::<RawBatchMessage>(2);

    // Spawn the reader thread
    thread::spawn(move || {
        loop {
            match read_raw_batch(&mut reader, batch_size) {
                Ok((batch, eof)) => {
                    if !batch.is_empty() && tx.send(RawBatchMessage::Batch(batch)).is_err() {
                        break; // Receiver dropped
                    }
                    if eof {
                        let _ = tx.send(RawBatchMessage::Eof);
                        break;
                    }
                }
                Err(e) => {
                    let _ = tx.send(RawBatchMessage::Error(e.to_string()));
                    break;
                }
            }
        }
    });

    Ok((rx, header))
}

/// Deserialize raw BAM record bytes into a noodles `RecordBuf`.
///
/// **Intentional non-raw decode.** Other production hot paths avoid
/// `raw_records_to_record_bufs` because it builds a synthetic BAM stream per call;
/// here it is acceptable because this function only runs *off* the comparison hot
/// path — after the byte/structured tiers have already flagged a pair as
/// mismatched and we need typed field values to render a human-readable diff.
///
/// The raw bytes are the BAM record body WITHOUT the 4-byte length prefix.
///
/// # Errors
///
/// Returns an error if the raw bytes cannot be parsed as a valid BAM record.
fn deserialize_raw_record(raw: &RawRecord) -> Result<RecordBuf> {
    let bufs = fgumi_raw_bam::raw_records_to_record_bufs(&[raw.as_ref().to_vec()])?;
    bufs.into_iter().next().ok_or_else(|| anyhow!("failed to deserialize raw BAM record"))
}

/// Compares a batch of raw record pairs in parallel using a three-tier strategy.
///
/// **Tier 1:** Full byte memcmp — if records are byte-identical, return immediately.
/// **Tier 2:** Structured raw comparison — compare core fields and tags at the byte
///   level without full deserialization. If core and tags match (possibly with different
///   tag order), return without deserializing.
/// **Tier 3:** Targeted deserialization for diff reporting — deserialize to `RecordBuf`
///   only for tag-diff cases that need typed tag display; render core-field diffs
///   directly from raw bytes.
///
/// Returns `(results, core_matches, core_diffs, tag_matches, tag_diffs, tag_order_diffs)`.
fn compare_raw_batch_parallel(
    batch1: &[RawRecord],
    batch2: &[RawRecord],
    header1: &Header,
    header2: &Header,
    start_index: u64,
) -> (Vec<RecordCompareResult>, usize, usize, usize, usize, usize) {
    let results: Vec<_> = batch1
        .par_iter()
        .zip(batch2.par_iter())
        .enumerate()
        .map(|(i, (r1, r2))| {
            let record_num = start_index + i as u64 + 1;

            // Tier 1: Full byte memcmp — handles the common case (identical BAMs)
            if raw_records_byte_equal(r1, r2) {
                return RecordCompareResult {
                    core_match: true,
                    tags_match: true,
                    tag_order_match: true,
                    diff_detail: None,
                };
            }

            // Tier 2: Structured raw comparison — avoids full deserialization
            // Note: RawCompareResult.tags_match = byte-identical tags,
            //       RawCompareResult.tag_order_match = semantically equal (order-independent)
            // RecordCompareResult.tags_match = semantic match, .tag_order_match = byte-identical
            let raw_result = raw_compare_structured(r1, r2);
            if raw_result.core_match && raw_result.tag_order_match {
                return RecordCompareResult {
                    core_match: true,
                    tags_match: true,
                    tag_order_match: raw_result.tags_match,
                    diff_detail: None,
                };
            }

            // Tier 3: diff reporting — two sub-cases:
            //   (a) core matches but tags differ: deserialize for typed tag value display
            //   (b) core fields differ: use raw bytes directly (no deserialization needed)
            let diff_detail = if raw_result.core_match {
                // (a) Core matches but tags differ — deserialize for typed tag value display
                let rec1 = match deserialize_raw_record(r1) {
                    Ok(r) => r,
                    Err(e) => {
                        return RecordCompareResult {
                            core_match: false,
                            tags_match: false,
                            tag_order_match: false,
                            diff_detail: Some(DiffDetail {
                                record_num,
                                qname: format!("<deserialization error: {e}>"),
                                flags: String::new(),
                                diff_type: DiffType::TagDiff,
                                diffs: vec![format!("Failed to deserialize record from BAM1: {e}")],
                            }),
                        };
                    }
                };
                let rec2 = match deserialize_raw_record(r2) {
                    Ok(r) => r,
                    Err(e) => {
                        return RecordCompareResult {
                            core_match: false,
                            tags_match: false,
                            tag_order_match: false,
                            diff_detail: Some(DiffDetail {
                                record_num,
                                qname: record_name_to_string(&rec1),
                                flags: rec1.flags().bits().to_string(),
                                diff_type: DiffType::TagDiff,
                                diffs: vec![format!("Failed to deserialize record from BAM2: {e}")],
                            }),
                        };
                    }
                };
                let tags1 = get_tags_map(&rec1);
                let tags2 = get_tags_map(&rec2);
                let mut all_tags: Vec<&String> = tags1.keys().chain(tags2.keys()).collect();
                all_tags.sort();
                all_tags.dedup();
                let diffs: Vec<String> = all_tags
                    .iter()
                    .filter_map(|tag| {
                        let v1 = tags1.get(*tag);
                        let v2 = tags2.get(*tag);
                        if v1 == v2 {
                            None
                        } else {
                            Some(format!(
                                "{tag}:\n{}\n",
                                CompareBams::format_diff(
                                    format!("{v1:?}"),
                                    format!("{v2:?}"),
                                    "      "
                                )
                            ))
                        }
                    })
                    .collect();
                let qname = record_name_to_string(&rec1);
                Some(DiffDetail {
                    record_num,
                    qname,
                    flags: rec1.flags().bits().to_string(),
                    diff_type: DiffType::TagDiff,
                    diffs,
                })
            } else {
                // (b) Core fields differ — extract from raw bytes via RawRecordView
                let core1 = get_core_fields_raw(r1, header1);
                let core2 = get_core_fields_raw(r2, header2);
                let diffs: Vec<String> = core1
                    .iter()
                    .zip(core2.iter())
                    .zip(FIELD_NAMES.iter())
                    .filter(|((v1, v2), _)| v1 != v2)
                    .map(|((v1, v2), name)| {
                        format!(
                            "{name}:\n{}\n",
                            CompareBams::format_diff(
                                format!("{v1:?}"),
                                format!("{v2:?}"),
                                "      "
                            )
                        )
                    })
                    .collect();
                Some(DiffDetail {
                    record_num,
                    qname: core1[0].clone(),
                    flags: core1[1].clone(),
                    diff_type: DiffType::CoreDiff,
                    diffs,
                })
            };

            RecordCompareResult {
                core_match: raw_result.core_match,
                tags_match: raw_result.tag_order_match,
                tag_order_match: raw_result.tags_match,
                diff_detail,
            }
        })
        .collect();

    // Aggregate stats
    let mut core_matches = 0usize;
    let mut core_diffs = 0usize;
    let mut tag_matches = 0usize;
    let mut tag_diffs = 0usize;
    let mut tag_order_diffs = 0usize;

    for r in &results {
        if r.core_match {
            core_matches += 1;
        } else {
            core_diffs += 1;
        }
        if r.tags_match {
            tag_matches += 1;
            if !r.tag_order_match {
                tag_order_diffs += 1;
            }
        } else {
            tag_diffs += 1;
        }
    }

    (results, core_matches, core_diffs, tag_matches, tag_diffs, tag_order_diffs)
}

/// Compare a batch of raw records for grouping data (read name, R1/R2 flag, MI tag).
///
/// Uses zero-copy field accessors on raw BAM bytes instead of decoded `RecordBuf`.
/// Read names are compared as raw bytes; String conversion only happens for error reporting.
///
/// Always populates `diff_detail` (and `read_name_for_display` for missing-MI errors)
/// when a mismatch is found; the caller is responsible for truncating to `max_diffs`
/// once results are ordered, since the parallel batch has no way to know how many
/// diffs precede any given record within the same batch.
fn compare_raw_batch_grouping_parallel(
    batch1: &[RawRecord],
    batch2: &[RawRecord],
    start_index: u64,
) -> Vec<GroupingCompareResult> {
    batch1
        .par_iter()
        .zip(batch2.par_iter())
        .enumerate()
        .map(|(i, (r1, r2))| {
            let record_num = start_index + i as u64 + 1;
            let name1_bytes = fgumi_raw_bam::RawRecordView::new(r1.as_ref()).read_name();
            let name2_bytes = fgumi_raw_bam::RawRecordView::new(r2.as_ref()).read_name();
            let is_read1_r1 = is_first_segment_raw(r1);
            let is_read1_r2 = is_first_segment_raw(r2);

            let name_match = name1_bytes == name2_bytes;
            let flag_match = is_read1_r1 == is_read1_r2;

            let record_key = record_key::record_key(r1);
            let mi1 = get_mi_tag_raw(r1);
            let mi2 = get_mi_tag_raw(r2);

            let diff_detail = if !name_match || !flag_match {
                // Only allocate Strings when there is an actual mismatch to report.
                let qname1 = String::from_utf8_lossy(name1_bytes).into_owned();
                if name_match {
                    Some(DiffDetail {
                        record_num,
                        qname: qname1,
                        flags: fgumi_raw_bam::RawRecordView::new(r1.as_ref()).flags().to_string(),
                        diff_type: DiffType::FlagMismatch,
                        diffs: vec![format!(
                            "R1/R2 flags differ: is_read1={} vs is_read1={}",
                            is_read1_r1, is_read1_r2
                        )],
                    })
                } else {
                    let qname2 = String::from_utf8_lossy(name2_bytes).into_owned();
                    Some(DiffDetail {
                        record_num,
                        qname: qname1,
                        flags: fgumi_raw_bam::RawRecordView::new(r1.as_ref()).flags().to_string(),
                        diff_type: DiffType::ReadNameMismatch,
                        diffs: vec![format!(
                            "Read names differ: '{}' vs '{}'",
                            String::from_utf8_lossy(name1_bytes),
                            qname2
                        )],
                    })
                }
            } else {
                None
            };

            // Populate read name for display for missing-MI error reporting
            let read_name_for_display = if mi1.is_none() ^ mi2.is_none() {
                Some(String::from_utf8_lossy(name1_bytes).into_owned())
            } else {
                None
            };

            GroupingCompareResult {
                record_num,
                record_key,
                read_name_for_display,
                mi1,
                mi2,
                name_match,
                flag_match,
                diff_detail,
            }
        })
        .collect()
}

impl Command for CompareBams {
    fn execute(&self, _command_line: &str) -> Result<()> {
        validate_file_exists(&self.bam1, "First BAM")?;
        validate_file_exists(&self.bam2, "Second BAM")?;

        // `sort` has no notion of `CompareMode`/`ContentPredicate` — it verifies sort
        // order and compares by sort-key run instead of pairing records positionally or
        // by key-join (see `CommandPreset::Sort`'s doc comment). Intercept before
        // `effective_settings()`/predicate resolution so those stay meaningful for every
        // other preset.
        if matches!(self.command, Some(CommandPreset::Sort)) {
            if self.mode.is_some() {
                anyhow::bail!(
                    "--mode is not valid with --command sort; sort verification uses its \
                     own dedicated engine (see `fgumi compare bams --help`)"
                );
            }
            if self.ignore_order.is_some() {
                anyhow::bail!(
                    "--ignore-order is not valid with --command sort; sort verification \
                     uses its own dedicated engine (see `fgumi compare bams --help`)"
                );
            }
            let timer = OperationTimer::new("Comparing BAMs");
            let total_records = self.execute_sort_verify()?;
            timer.log_completion(total_records);
            return Ok(());
        }

        let (mode, ignore_order) = self.effective_settings();

        // Only bail when the user explicitly passed `--ignore-order=true` with a
        // non-grouping mode. Preset-inherited `ignore_order` is silently dropped
        // by `effective_settings()` when the resolved mode isn't `Grouping`, so
        // combinations like `--command simplex --mode content` are not user
        // errors and must not produce an `--ignore-order` error message.
        if self.ignore_order == Some(true) && !matches!(mode, CompareMode::Grouping) {
            anyhow::bail!("--ignore-order is only valid with --mode grouping");
        }

        if let Some(preset) = self.command {
            let overridden = self.mode.is_some() || self.ignore_order.is_some();
            let preset_name = preset
                .to_possible_value()
                .expect("CommandPreset is not skipped")
                .get_name()
                .to_string();
            let mode_name = mode
                .to_possible_value()
                .expect("CompareMode is not skipped")
                .get_name()
                .to_string();
            info!(
                "Using --command {} preset: mode={}, ignore-order={}{}",
                preset_name,
                mode_name,
                ignore_order,
                if overridden { " (with explicit overrides)" } else { "" }
            );
        }

        let timer = OperationTimer::new("Comparing BAMs");

        // When a `--command` preset is given, its content predicate comes from
        // `CommandPreset::resolve` (see that method's doc comment for the full per-preset
        // rationale, e.g. why `Simplex`/`Duplex`/`Codec`/`Filter` use the
        // `ExactConsensus`). This predicate is only ever consumed below by the `Content`
        // mode branch (`self.execute_content(predicate)`); `CompareMode::Grouping` routes to
        // `execute_grouping_with`, which takes no predicate at all — the ordered
        // `execute_grouping` path has no content check to gate, and the `--ignore-order`
        // key-join path (`execute_grouping_unordered`) hardcodes
        // `ContentPredicate::ExactMinusMi` inside `engines::keyjoin::keyjoin_compare` itself,
        // not configurable via `--command`/`--mode`. The `None if matches!(mode,
        // CompareMode::Grouping)` arm below therefore only matters for the (rare, preset-less)
        // case where `--mode content` isn't in play; it does not affect the key-join path.
        let predicate = match self.command {
            Some(preset) => preset.content_predicate(),
            None if matches!(mode, CompareMode::Grouping) => ContentPredicate::ExactMinusMi,
            None => ContentPredicate::Exact,
        };

        let total_records = match mode {
            CompareMode::Full => self.execute_full()?,
            CompareMode::Content => self.execute_content(predicate)?,
            CompareMode::Grouping => self.execute_grouping_with(ignore_order)?,
        };

        timer.log_completion(total_records);
        Ok(())
    }
}

impl CompareBams {
    /// Resolve the effective `(mode, ignore_order)` pair.
    ///
    /// Precedence: explicit flag > `--command` preset default > built-in default
    /// (`full`, `false`). `ignore_order` is silently coerced to `false` when the
    /// resolved mode isn't `Grouping`, so preset-inherited `ignore_order=true`
    /// is dropped cleanly when an explicit `--mode` narrows to a non-grouping
    /// comparison (e.g. `--command simplex --mode content`). An *explicit*
    /// `--ignore-order=true` with a non-grouping mode is still rejected in
    /// `execute()`.
    fn effective_settings(&self) -> (CompareMode, bool) {
        let preset = self.command.map(|p| p.defaults());
        let mode = self.mode.or(preset.map(|(m, _)| m)).unwrap_or_default();
        let ignore_order = self.ignore_order.or(preset.map(|(_, io)| io)).unwrap_or(false);
        let ignore_order = ignore_order && matches!(mode, CompareMode::Grouping);
        (mode, ignore_order)
    }

    fn format_diff(left: String, right: String, leading: &str) -> String {
        let left_vec: Vec<char> = left.chars().collect();
        let right_vec: Vec<char> = right.chars().collect();
        let mut diffs: Vec<usize> = Vec::new();

        let min_length = left.len().min(right.len());
        let max_length = left.len().max(right.len());

        let mut left_str = String::new();
        let mut right_str = String::new();
        let mut aln_str = String::new();
        let mut i = 0;
        while i < min_length {
            left_str.push(left_vec[i]);
            right_str.push(right_vec[i]);
            if left_vec[i] == right_vec[i] {
                aln_str.push(' ');
            } else {
                aln_str.push('X');
                diffs.push(i);
            }
            i += 1;
        }
        while i < max_length {
            if i < left_vec.len() {
                left_str.push(left_vec[i]);
            } else {
                left_str.push('-');
            }
            if i < right_vec.len() {
                right_str.push(right_vec[i]);
            } else {
                right_str.push('-');
            }
            aln_str.push('-');
            diffs.push(i);
            i += 1;
        }
        let diff = diffs.iter().map(|i| format!("{i}")).join(", ");
        format!("{leading}{left_str}\n{leading}{aln_str}\n{leading}{right_str}\n{leading}{diff}")
    }

    /// Execute content comparison mode.
    ///
    /// Delegates pairing and equality entirely to the positional engine
    /// ([`positional_compare`]): records are paired purely by index, a
    /// [`RecordKey`](super::record_key::RecordKey) mismatch stops pairing immediately
    /// (never resyncing), and remaining pairs are compared under `predicate` (plain
    /// `Exact`, or the `ExactConsensus` for the consensus and `filter`
    /// presets — see `execute()`).
    ///
    /// This preserves the external report contract other tooling depends on: the
    /// `RESULT: BAM files are IDENTICAL` / `RESULT: BAM files DIFFER` line (the
    /// benchmark suite greps `RESULT:.*DIFFER`) and the exit-1-on-mismatch behavior via
    /// [`super::CompareMismatch`]. The detailed per-stat breakdown *is not* preserved
    /// byte-for-byte: the old `core_matches`/`tag_matches`/`tag_diffs`/`tag_order_diffs`
    /// counters and the "tags in different order" note are specific to the retired
    /// batch-parallel comparator and have no equivalent in
    /// [`PositionalOutcome`](super::engines::positional::PositionalOutcome) — this
    /// report now shows record counts, a content-diff count, and (if pairing
    /// desynced) the first `RecordKey` mismatch index instead.
    fn execute_content(&self, predicate: ContentPredicate) -> Result<u64> {
        info!(
            "Starting content comparison with {} threads, batch size {}",
            self.threads, self.batch_size
        );

        let outcome = positional_compare(
            &self.bam1,
            &self.bam2,
            self.threads,
            self.batch_size,
            self.max_diffs,
            predicate,
        )?;

        let is_equal = outcome.is_match();

        if !self.quiet {
            println!("=== BAM Comparison Results (content mode) ===");
            println!("BAM1: {}", self.bam1.display());
            println!("BAM2: {}", self.bam2.display());
            println!();
            println!("Record counts: {} vs {}", outcome.bam1_count, outcome.bam2_count);
            println!("Content diffs: {}", outcome.content_diffs);
            if let Some(index) = outcome.key_mismatch_at {
                println!("First RecordKey mismatch: record {index} (pairing stopped)");
            }
            println!();

            if is_equal {
                println!("RESULT: BAM files are IDENTICAL (core fields and tag values match)");
            } else {
                println!("RESULT: BAM files DIFFER");
                if !outcome.diff_details.is_empty() {
                    println!("\nFirst {} differences:", outcome.diff_details.len());
                    for detail in &outcome.diff_details {
                        println!("  {detail}");
                    }
                }
            }
        }

        if is_equal {
            info!("BAM files are identical");
            Ok(outcome.bam1_count)
        } else {
            info!("BAM files differ");
            Err(super::CompareMismatch("BAM files differ".to_owned()).into())
        }
    }

    /// Execute sort-order verification for the `--command sort` preset.
    ///
    /// Delegates entirely to
    /// [`sort_verify_compare`](super::engines::sort_verify::sort_verify_compare): detects
    /// the shared sort order from both inputs' `@HD` header, verifies each file is itself
    /// correctly ordered, and compares the two files as a multiset grouped by maximal
    /// equal-core-sort-key run (see that function's doc comment, and
    /// [`CommandPreset::Sort`]'s). Preserves the same external report contract as
    /// `execute_content`: the `RESULT: BAM files are IDENTICAL` / `RESULT: BAM files
    /// DIFFER` line and exit-1-on-mismatch behavior via [`super::CompareMismatch`].
    fn execute_sort_verify(&self) -> Result<u64> {
        info!("Using --command sort preset: sort-key-run verification (no positional pairing)");

        let outcome = super::engines::sort_verify::sort_verify_compare(
            &self.bam1,
            &self.bam2,
            self.max_diffs,
        )?;

        let is_equal = outcome.is_match();

        if !self.quiet {
            println!("=== BAM Comparison Results (sort mode) ===");
            println!("BAM1: {}", self.bam1.display());
            println!("BAM2: {}", self.bam2.display());
            println!();
            println!("Detected sort order: {:?}", outcome.sort_order);
            println!("Record counts: {} vs {}", outcome.bam1_count, outcome.bam2_count);
            print!("bam1 sort-order violations: {}", outcome.bam1_violations);
            if let Some((record_num, name)) = &outcome.bam1_first_violation {
                print!(" (first at record {record_num}: {name})");
            }
            println!();
            print!("bam2 sort-order violations: {}", outcome.bam2_violations);
            if let Some((record_num, name)) = &outcome.bam2_first_violation {
                print!(" (first at record {record_num}: {name})");
            }
            println!();
            println!("Sort-key-run multiset mismatches: {}", outcome.run_mismatches);
            println!();

            if is_equal {
                println!(
                    "RESULT: BAM files are IDENTICAL (sort order verified; run multisets match)"
                );
            } else {
                println!("RESULT: BAM files DIFFER");
                if !outcome.diff_details.is_empty() {
                    println!("\nFirst {} differences:", outcome.diff_details.len());
                    for detail in &outcome.diff_details {
                        println!("  {detail}");
                    }
                }
            }
        }

        if is_equal {
            info!("BAM files are identical");
            Ok(outcome.bam1_count)
        } else {
            info!("BAM files differ");
            Err(super::CompareMismatch("BAM files differ".to_owned()).into())
        }
    }

    /// Execute full comparison mode
    ///
    /// This mode combines grouping equivalence check with full content comparison.
    /// Both files must be in the same order (e.g., query-name sorted).
    /// First verifies MI groupings are equivalent, then compares all other fields.
    /// Uses parallel batch processing with double buffering for performance.
    ///
    /// Unlike `execute_content` (which delegates to the positional engine), this path is
    /// its own index-paired comparator and does not go through `RecordKey`-gated pairing —
    /// it remains reorder-unsound. Retained pending the compare-hardening Phase 3+ rework
    /// (see `docs/superpowers/plans/2026-07-08-compare-hardening-phase2.md`).
    fn execute_full(&self) -> Result<u64> {
        let mut stats = CompareStats::default();
        let mut grouping_stats = GroupingStats::default();
        let mut grouping_errors: Vec<String> = Vec::new();
        let batch_size = self.batch_size;

        info!("Starting full comparison with {} threads, batch size {}", self.threads, batch_size);

        // Maps: record_key -> MI value for each BAM (compact MiKey representation)
        let mut mi_map1: AHashMap<RecordKey, MiKey> = AHashMap::new();
        let mut mi_map2: AHashMap<RecordKey, MiKey> = AHashMap::new();

        // Start double-buffered readers for both BAM files
        let (rx1, header1) = start_raw_batch_reader(self.bam1.clone(), self.threads, batch_size)?;
        let (rx2, header2) = start_raw_batch_reader(self.bam2.clone(), self.threads, batch_size)?;

        // Progress tracking
        let progress = ProgressTracker::new("Processed records").with_interval(1_000_000);

        // Process batches
        let mut bam1_eof = false;
        let mut bam2_eof = false;
        let mut pending_batch1: Option<Vec<RawRecord>> = None;
        let mut pending_batch2: Option<Vec<RawRecord>> = None;
        let mut current_index = 0u64;

        loop {
            // Get next batch from BAM1 if needed
            if pending_batch1.is_none() && !bam1_eof {
                match rx1.recv() {
                    Ok(RawBatchMessage::Batch(batch)) => {
                        stats.bam1_count += batch.len() as u64;
                        pending_batch1 = Some(batch);
                    }
                    Ok(RawBatchMessage::Eof) => bam1_eof = true,
                    Ok(RawBatchMessage::Error(e)) => bail!("Error reading BAM1: {e}"),
                    Err(_) => bam1_eof = true,
                }
            }

            // Get next batch from BAM2 if needed
            if pending_batch2.is_none() && !bam2_eof {
                match rx2.recv() {
                    Ok(RawBatchMessage::Batch(batch)) => {
                        stats.bam2_count += batch.len() as u64;
                        pending_batch2 = Some(batch);
                    }
                    Ok(RawBatchMessage::Eof) => bam2_eof = true,
                    Ok(RawBatchMessage::Error(e)) => bail!("Error reading BAM2: {e}"),
                    Err(_) => bam2_eof = true,
                }
            }

            // Check for completion
            match (&pending_batch1, &pending_batch2) {
                (None, None) => break,
                (Some(_), None) | (None, Some(_)) => {
                    if stats.diff_details.len() < self.max_diffs {
                        stats.diff_details.push(DiffDetail {
                            record_num: current_index,
                            qname: "N/A".to_string(),
                            flags: "N/A".to_string(),
                            diff_type: DiffType::CountMismatch,
                            diffs: vec!["BAM files have different number of records".to_string()],
                        });
                    }
                    // Drain remaining batches for accurate counts
                    if pending_batch1.is_some() {
                        while let Ok(msg) = rx1.recv() {
                            if let RawBatchMessage::Batch(batch) = msg {
                                stats.bam1_count += batch.len() as u64;
                            }
                        }
                    }
                    if pending_batch2.is_some() {
                        while let Ok(msg) = rx2.recv() {
                            if let RawBatchMessage::Batch(batch) = msg {
                                stats.bam2_count += batch.len() as u64;
                            }
                        }
                    }
                    break;
                }
                (Some(_), Some(_)) => {}
            }

            let batch1 = pending_batch1.take().expect("guarded by (Some, Some) match above");
            let batch2 = pending_batch2.take().expect("guarded by (Some, Some) match above");

            let min_len = batch1.len().min(batch2.len());
            let (cmp_batch1, remainder1) = batch1.split_at(min_len);
            let (cmp_batch2, remainder2) = batch2.split_at(min_len);

            // Compare batches in parallel and collect MI data
            let (results, core_m, core_d, tag_m, tag_d, tag_ord) = compare_raw_batch_parallel(
                cmp_batch1,
                cmp_batch2,
                &header1,
                &header2,
                current_index,
            );

            // Process grouping data from the batch (sequential for MI map building)
            for (r1, r2) in cmp_batch1.iter().zip(cmp_batch2.iter()) {
                grouping_stats.total_records += 1;

                let name1_bytes = fgumi_raw_bam::RawRecordView::new(r1.as_ref()).read_name();
                let name2_bytes = fgumi_raw_bam::RawRecordView::new(r2.as_ref()).read_name();

                if name1_bytes != name2_bytes {
                    grouping_stats.order_mismatches += 1;
                    continue;
                }

                let record_key = record_key::record_key(r1);
                // Track missing MI tags explicitly so that two ungrouped BAMs don't
                // appear equivalent simply because neither inserts into its map.
                match (get_mi_tag_raw(r1), get_mi_tag_raw(r2)) {
                    (Some(mi1), Some(mi2)) => {
                        mi_map1.insert(record_key.clone(), mi1);
                        mi_map2.insert(record_key, mi2);
                    }
                    (None, Some(mi2)) => {
                        stats.missing_mi_bam1 += 1;
                        mi_map2.insert(record_key, mi2);
                    }
                    (Some(mi1), None) => {
                        stats.missing_mi_bam2 += 1;
                        mi_map1.insert(record_key, mi1);
                    }
                    (None, None) => {
                        stats.missing_mi_bam1 += 1;
                        stats.missing_mi_bam2 += 1;
                    }
                }
            }

            stats.core_matches += core_m as u64;
            stats.core_diffs += core_d as u64;
            stats.tag_matches += tag_m as u64;
            stats.tag_diffs += tag_d as u64;
            stats.tag_order_diffs += tag_ord as u64;

            for r in results {
                if let Some(detail) = r.diff_detail
                    && stats.diff_details.len() < self.max_diffs
                {
                    stats.diff_details.push(detail);
                }
            }

            current_index += min_len as u64;
            progress.log_if_needed(min_len as u64);

            if !remainder1.is_empty() {
                pending_batch1 = Some(remainder1.to_vec());
            }
            if !remainder2.is_empty() {
                pending_batch2 = Some(remainder2.to_vec());
            }
        }

        progress.log_final();
        info!("Phase 2: Verifying grouping equivalence...");

        // Phase 2: Verify grouping equivalence (using compact representation)
        let mi_to_reads1 = build_mi_groups_compact(&mi_map1);
        let mi_to_reads2 = build_mi_groups_compact(&mi_map2);

        let unique_mi1_count = mi_to_reads1.len();
        let unique_mi2_count = mi_to_reads2.len();

        // For each MI group in BAM1, verify all reads have the same MI in BAM2
        for (mi1, record_keys) in &mi_to_reads1 {
            let mi2_values: AHashSet<MiKey> =
                record_keys.iter().filter_map(|k| mi_map2.get(k).copied()).collect();

            if mi2_values.len() > 1 {
                grouping_stats.grouping_mismatches += 1;
                if grouping_errors.len() < self.max_diffs {
                    grouping_errors.push(format!(
                        "MI group '{}' in BAM1 ({} reads) maps to {} different MIs in BAM2: [{}]",
                        mi1,
                        record_keys.len(),
                        mi2_values.len(),
                        mi2_values.iter().take(5).map(MiKey::to_string).join(", ")
                    ));
                }
            }
        }

        // Verify the reverse: each MI group in BAM2 maps to single MI in BAM1
        for (mi2, record_keys) in &mi_to_reads2 {
            let mi1_values: AHashSet<MiKey> =
                record_keys.iter().filter_map(|k| mi_map1.get(k).copied()).collect();

            if mi1_values.len() > 1 {
                grouping_stats.grouping_mismatches += 1;
                if grouping_errors.len() < self.max_diffs {
                    grouping_errors.push(format!(
                        "MI group '{}' in BAM2 ({} reads) maps to {} different MIs in BAM1: [{}]",
                        mi2,
                        record_keys.len(),
                        mi1_values.len(),
                        mi1_values.iter().take(5).map(MiKey::to_string).join(", ")
                    ));
                }
            }
        }

        let content_equal =
            stats.bam1_count == stats.bam2_count && stats.core_diffs == 0 && stats.tag_diffs == 0;
        let grouping_equal = grouping_stats.order_mismatches == 0
            && grouping_stats.grouping_mismatches == 0
            && stats.missing_mi_bam1 == 0
            && stats.missing_mi_bam2 == 0;
        let is_equal = content_equal && grouping_equal;

        if !self.quiet {
            println!("=== BAM Comparison Results (full mode) ===");
            println!("BAM1: {}", self.bam1.display());
            println!("BAM2: {}", self.bam2.display());
            println!();
            println!("--- Content Comparison ---");
            println!("Record counts: {} vs {}", stats.bam1_count, stats.bam2_count);
            println!("Core field matches: {}", stats.core_matches);
            println!("Core field diffs: {}", stats.core_diffs);
            println!("Tag value matches: {}", stats.tag_matches);
            println!("Tag value diffs: {}", stats.tag_diffs);
            println!("Tag order diffs (values match): {}", stats.tag_order_diffs);
            println!();
            println!("--- Grouping Comparison ---");
            println!("Total records compared: {}", grouping_stats.total_records);
            println!("Unique MI values: {} vs {}", unique_mi1_count, unique_mi2_count);
            println!("Order mismatches: {}", grouping_stats.order_mismatches);
            println!("Missing MI in BAM1: {}", stats.missing_mi_bam1);
            println!("Missing MI in BAM2: {}", stats.missing_mi_bam2);
            println!("Grouping mismatches: {}", grouping_stats.grouping_mismatches);
            println!();

            if is_equal {
                println!("RESULT: BAM files are IDENTICAL (content and groupings match)");
                if stats.tag_order_diffs > 0 {
                    println!(
                        "  Note: {} records have tags in different order",
                        stats.tag_order_diffs
                    );
                }
                if unique_mi1_count != unique_mi2_count {
                    println!(
                        "  Note: Different number of unique MI values ({} vs {}), but groupings match",
                        unique_mi1_count, unique_mi2_count
                    );
                }
            } else {
                println!("RESULT: BAM files DIFFER");
                if !stats.diff_details.is_empty() {
                    println!("\nContent differences (first {}):", stats.diff_details.len());
                    for detail in &stats.diff_details {
                        println!("  Record {}: {}", detail.record_num, detail.qname);
                        println!("    Flag: {}", detail.flags);
                        println!("    Type: {}", detail.diff_type);
                        for d in &detail.diffs {
                            println!("      {d}");
                        }
                    }
                }
                if !grouping_errors.is_empty() {
                    println!("\nGrouping mismatches (first {}):", grouping_errors.len());
                    for err in &grouping_errors {
                        println!("  {err}");
                    }
                }
            }
        }

        if is_equal {
            info!("BAM files are identical (content and groupings)");
            Ok(stats.bam1_count)
        } else {
            info!("BAM files differ");
            Err(super::CompareMismatch("BAM files differ (content and groupings)".to_owned())
                .into())
        }
    }

    /// Execute grouping comparison mode
    ///
    /// This mode compares grouped BAM files where MI assignment order may differ.
    /// Both files must be in the same order (e.g., query-name sorted), unless
    /// --ignore-order is specified.
    /// We validate read names and R1/R2 flags match, then verify that reads
    /// with the same MI in one file have the same MI in the other.
    /// Uses parallel batch processing with double buffering for performance.
    /// Dispatches between ordered and unordered grouping comparison based on
    /// the resolved `ignore_order` flag (which may come from `--ignore-order`
    /// directly or from a `--command` preset).
    ///
    /// Like `execute_full`, `execute_grouping` (the *ordered* path, no `--ignore-order`)
    /// remains reorder-unsound (no `RecordKey`-gated pairing) and MI-equivalence-only; it is
    /// out of scope for the compare-hardening Phase 3 key-join rework (see
    /// `docs/superpowers/plans/2026-07-08-compare-hardening-phase3.md`'s Task 3.4
    /// "Self-review" for the scope boundary). Neither branch consumes a content predicate:
    /// `execute_grouping` (ordered) never took one, and `execute_grouping_unordered` (the
    /// `--ignore-order` path, which is what the `group` preset uses) delegates to the
    /// key-join engine (`engines::keyjoin::keyjoin_compare`), which always compares content
    /// under a hardcoded [`ContentPredicate::ExactMinusMi`] (not configurable via
    /// `--command`/`--mode`) in addition to the MI bijection check.
    fn execute_grouping_with(&self, ignore_order: bool) -> Result<u64> {
        if ignore_order { self.execute_grouping_unordered() } else { self.execute_grouping() }
    }

    fn execute_grouping(&self) -> Result<u64> {
        let mut stats = GroupingStats::default();
        let mut diff_details: Vec<DiffDetail> = Vec::new();
        let batch_size = self.batch_size;

        info!(
            "Starting grouping comparison with {} threads, batch size {}",
            self.threads, batch_size
        );

        // Maps: record_key -> MI value for each BAM (compact MiKey representation)
        let mut mi_map1: AHashMap<RecordKey, MiKey> = AHashMap::new();
        let mut mi_map2: AHashMap<RecordKey, MiKey> = AHashMap::new();

        // Start double-buffered raw readers for both BAM files
        let (rx1, _header1) = start_raw_batch_reader(self.bam1.clone(), self.threads, batch_size)?;
        let (rx2, _header2) = start_raw_batch_reader(self.bam2.clone(), self.threads, batch_size)?;

        // Progress tracking
        let progress = ProgressTracker::new("Processed records").with_interval(1_000_000);

        // Process batches
        let mut bam1_eof = false;
        let mut bam2_eof = false;
        let mut pending_batch1: Option<Vec<RawRecord>> = None;
        let mut pending_batch2: Option<Vec<RawRecord>> = None;
        let mut current_index = 0u64;

        loop {
            // Get next batch from BAM1 if needed
            if pending_batch1.is_none() && !bam1_eof {
                match rx1.recv() {
                    Ok(RawBatchMessage::Batch(batch)) => {
                        pending_batch1 = Some(batch);
                    }
                    Ok(RawBatchMessage::Eof) => bam1_eof = true,
                    Ok(RawBatchMessage::Error(e)) => bail!("Error reading BAM1: {e}"),
                    Err(_) => bam1_eof = true,
                }
            }

            // Get next batch from BAM2 if needed
            if pending_batch2.is_none() && !bam2_eof {
                match rx2.recv() {
                    Ok(RawBatchMessage::Batch(batch)) => {
                        pending_batch2 = Some(batch);
                    }
                    Ok(RawBatchMessage::Eof) => bam2_eof = true,
                    Ok(RawBatchMessage::Error(e)) => bail!("Error reading BAM2: {e}"),
                    Err(_) => bam2_eof = true,
                }
            }

            // Check for completion
            match (&pending_batch1, &pending_batch2) {
                (None, None) => break,
                (Some(_), None) | (None, Some(_)) => {
                    stats.count_mismatch = true;
                    if diff_details.len() < self.max_diffs {
                        diff_details.push(DiffDetail {
                            record_num: current_index,
                            qname: "N/A".to_string(),
                            flags: "N/A".to_string(),
                            diff_type: DiffType::CountMismatch,
                            diffs: vec!["BAM files have different number of records".to_string()],
                        });
                    }
                    // Drain remaining batches so downstream readers finish cleanly.
                    if pending_batch1.is_some() {
                        while rx1.recv().is_ok() {}
                    }
                    if pending_batch2.is_some() {
                        while rx2.recv().is_ok() {}
                    }
                    break;
                }
                (Some(_), Some(_)) => {}
            }

            let batch1 = pending_batch1.take().expect("guarded by (Some, Some) match above");
            let batch2 = pending_batch2.take().expect("guarded by (Some, Some) match above");

            let min_len = batch1.len().min(batch2.len());
            let (cmp_batch1, remainder1) = batch1.split_at(min_len);
            let (cmp_batch2, remainder2) = batch2.split_at(min_len);

            // Compare batches in parallel for grouping data
            let results =
                compare_raw_batch_grouping_parallel(cmp_batch1, cmp_batch2, current_index);

            // Process results and build MI maps
            for r in results {
                stats.total_records += 1;

                if !r.name_match {
                    stats.order_mismatches += 1;
                    if let Some(detail) = r.diff_detail
                        && diff_details.len() < self.max_diffs
                    {
                        diff_details.push(detail);
                    }
                    continue;
                }

                if !r.flag_match {
                    stats.flag_mismatches += 1;
                    if let Some(detail) = r.diff_detail
                        && diff_details.len() < self.max_diffs
                    {
                        diff_details.push(detail);
                    }
                    continue;
                }

                match (r.mi1, r.mi2) {
                    (Some(mi1_val), Some(mi2_val)) => {
                        mi_map1.insert(r.record_key.clone(), mi1_val);
                        mi_map2.insert(r.record_key, mi2_val);
                    }
                    (None, Some(_)) => {
                        stats.missing_mi_bam1 += 1;
                        if diff_details.len() < self.max_diffs {
                            diff_details.push(DiffDetail {
                                record_num: r.record_num,
                                qname: r.read_name_for_display.unwrap_or_else(|| "?".to_string()),
                                flags: "N/A".to_string(),
                                diff_type: DiffType::TagDiff,
                                diffs: vec!["MI tag missing in BAM1".to_string()],
                            });
                        }
                    }
                    (Some(_), None) => {
                        stats.missing_mi_bam2 += 1;
                        if diff_details.len() < self.max_diffs {
                            diff_details.push(DiffDetail {
                                record_num: r.record_num,
                                qname: r.read_name_for_display.unwrap_or_else(|| "?".to_string()),
                                flags: "N/A".to_string(),
                                diff_type: DiffType::TagDiff,
                                diffs: vec!["MI tag missing in BAM2".to_string()],
                            });
                        }
                    }
                    (None, None) => {
                        // Count missing on both sides so two un-grouped BAMs cannot
                        // pass the grouping check by virtue of silently matching.
                        stats.missing_mi_bam1 += 1;
                        stats.missing_mi_bam2 += 1;
                        if diff_details.len() < self.max_diffs {
                            diff_details.push(DiffDetail {
                                record_num: r.record_num,
                                qname: r.read_name_for_display.unwrap_or_else(|| "?".to_string()),
                                flags: "N/A".to_string(),
                                diff_type: DiffType::TagDiff,
                                diffs: vec!["MI tag missing in both BAMs".to_string()],
                            });
                        }
                    }
                }
            }

            current_index += min_len as u64;
            progress.log_if_needed(min_len as u64);

            if !remainder1.is_empty() {
                pending_batch1 = Some(remainder1.to_vec());
            }
            if !remainder2.is_empty() {
                pending_batch2 = Some(remainder2.to_vec());
            }
        }

        progress.log_final();
        info!("Phase 2: Verifying grouping equivalence...");

        // Phase 2: Verify grouping equivalence (using compact representation)
        let bam1_groups = build_mi_groups_compact(&mi_map1);
        let bam2_groups = build_mi_groups_compact(&mi_map2);

        stats.unique_groups_bam1 = bam1_groups.len();
        stats.unique_groups_bam2 = bam2_groups.len();

        // Check BAM1 groups -> BAM2
        let mut grouping_errors: Vec<String> = Vec::new();
        for (mi1, record_keys) in &bam1_groups {
            let mi2_values: AHashSet<MiKey> = record_keys
                .iter()
                .filter_map(|record_key| mi_map2.get(record_key).copied())
                .collect();

            if mi2_values.len() > 1 {
                stats.grouping_mismatches += 1;
                if grouping_errors.len() < self.max_diffs {
                    grouping_errors.push(format!(
                        "MI group '{}' in BAM1 ({} reads) maps to {} different MIs in BAM2: [{}]",
                        mi1,
                        record_keys.len(),
                        mi2_values.len(),
                        mi2_values.iter().take(5).map(MiKey::to_string).join(", ")
                    ));
                }
            }
        }

        // Check BAM2 groups -> BAM1
        for (mi2, record_keys) in &bam2_groups {
            let mi1_values: AHashSet<MiKey> = record_keys
                .iter()
                .filter_map(|record_key| mi_map1.get(record_key).copied())
                .collect();

            if mi1_values.len() > 1 {
                stats.grouping_mismatches += 1;
                if grouping_errors.len() < self.max_diffs {
                    grouping_errors.push(format!(
                        "MI group '{}' in BAM2 ({} reads) maps to {} different MIs in BAM1: [{}]",
                        mi2,
                        record_keys.len(),
                        mi1_values.len(),
                        mi1_values.iter().take(5).map(MiKey::to_string).join(", ")
                    ));
                }
            }
        }

        // Verify molecule-level (duplex strand-pairing) equivalence — the
        // per-MiKey checks above keep `/A`/`/B` distinct and so cannot detect a
        // strand-pairing split (see `count_base_pairing_mismatches`). Fold any
        // base-level mismatches into the total so a split makes groupings DIFFER.
        let base_pairing_mismatches = count_base_pairing_mismatches(&mi_map1, &mi_map2);
        if base_pairing_mismatches > 0 {
            stats.grouping_mismatches += base_pairing_mismatches;
            if grouping_errors.len() < self.max_diffs {
                grouping_errors.push(format!(
                    "{base_pairing_mismatches} molecule(s) differ in duplex strand pairing \
                     (reads sharing a molecule in one BAM are split across molecules in the other)"
                ));
            }
        }

        // Determine result
        let is_equivalent = !stats.count_mismatch
            && stats.order_mismatches == 0
            && stats.flag_mismatches == 0
            && stats.missing_mi_bam1 == 0
            && stats.missing_mi_bam2 == 0
            && stats.grouping_mismatches == 0;

        if !self.quiet {
            println!("=== BAM Comparison Results (grouping mode) ===");
            println!("BAM1: {}", self.bam1.display());
            println!("BAM2: {}", self.bam2.display());
            println!();
            println!("Total records compared: {}", stats.total_records);
            if stats.count_mismatch {
                println!("Record count mismatch: BAM files have different number of records");
            }
            println!("Order/name mismatches: {}", stats.order_mismatches);
            println!("R1/R2 flag mismatches: {}", stats.flag_mismatches);
            println!("Missing MI in BAM1: {}", stats.missing_mi_bam1);
            println!("Missing MI in BAM2: {}", stats.missing_mi_bam2);
            println!("Unique MI groups in BAM1: {}", stats.unique_groups_bam1);
            println!("Unique MI groups in BAM2: {}", stats.unique_groups_bam2);
            println!("Grouping mismatches: {}", stats.grouping_mismatches);
            println!();

            if is_equivalent {
                println!("RESULT: BAM groupings are EQUIVALENT");
                println!("  Reads with the same MI in one file have the same MI in the other.");
                if stats.unique_groups_bam1 != stats.unique_groups_bam2 {
                    println!(
                        "  Note: Different number of unique MI values ({} vs {}), but groupings match.",
                        stats.unique_groups_bam1, stats.unique_groups_bam2
                    );
                }
            } else {
                println!("RESULT: BAM groupings DIFFER");

                if !diff_details.is_empty() {
                    println!("\nOrder/tag differences (first {}):", diff_details.len());
                    for detail in &diff_details {
                        println!("  Record {}: {}", detail.record_num, detail.qname);
                        println!("    Type: {}", detail.diff_type);
                        for d in &detail.diffs {
                            println!("      {d}");
                        }
                    }
                }

                if !grouping_errors.is_empty() {
                    println!("\nGrouping mismatches (first {}):", grouping_errors.len());
                    for err in &grouping_errors {
                        println!("  {err}");
                    }
                }
            }
        }

        if is_equivalent {
            info!("BAM groupings are equivalent");
            Ok(stats.total_records)
        } else {
            info!("BAM groupings differ");
            Err(super::CompareMismatch("BAM groupings differ".to_owned()).into())
        }
    }

    /// Execute grouping comparison in order-independent mode via the key-join engine.
    ///
    /// This is what `--command group` uses (its preset sets `ignore_order = true`): both
    /// inputs are canonicalized to queryname order and merge-joined by
    /// [`RecordKey`](super::record_key::RecordKey) (see [`keyjoin::keyjoin_compare`]), so —
    /// unlike the ordered `execute_grouping` path — this also catches a non-MI content
    /// difference on a matched pair, in addition to verifying the fgumi-MI/fgbio-MI mapping
    /// is a consistent bijection. The content check always uses
    /// [`ContentPredicate::ExactMinusMi`], hardcoded inside `keyjoin::keyjoin_compare` itself
    /// — this engine takes no predicate argument and is not configurable via
    /// `--command`/`--mode`.
    ///
    /// # Errors
    ///
    /// Returns an error if either input cannot be canonicalized or read (see
    /// [`keyjoin::keyjoin_compare`]), or [`super::CompareMismatch`] if the two BAMs are
    /// found to differ (non-zero exit via the `Command` trait).
    fn execute_grouping_unordered(&self) -> Result<u64> {
        info!(
            "Starting key-join grouping comparison with {} threads, sort memory {:?}",
            self.threads, self.sort_memory
        );

        let cfg = self.keyjoin_config()?;
        let outcome = keyjoin::keyjoin_compare(&self.bam1, &self.bam2, &cfg)?;
        let is_equivalent = outcome.is_match();

        if !self.quiet {
            println!("=== BAM Comparison Results (grouping mode, order-independent) ===");
            println!("BAM1: {}", self.bam1.display());
            println!("BAM2: {}", self.bam2.display());
            println!();
            println!("Total records in BAM1: {}", outcome.bam1_count);
            println!("Total records in BAM2: {}", outcome.bam2_count);
            println!("Records matched: {}", outcome.matched);
            println!("Records only in BAM1: {}", outcome.only_in_bam1);
            println!("Records only in BAM2: {}", outcome.only_in_bam2);
            println!("Missing MI in BAM1: {}", outcome.missing_mi_bam1);
            println!("Missing MI in BAM2: {}", outcome.missing_mi_bam2);
            println!("Content diffs (excluding MI): {}", outcome.content_diffs);
            println!("MI bijection mismatches: {}", outcome.mi_bijection_mismatches);
            println!();

            if is_equivalent {
                println!("RESULT: BAM groupings are EQUIVALENT");
                println!(
                    "  Content matches (excluding MI) and reads with the same MI in one file \
                     have the same MI in the other."
                );
            } else {
                println!("RESULT: BAM groupings DIFFER");

                if !outcome.diff_details.is_empty() {
                    println!("\nDifferences (first {}):", outcome.diff_details.len());
                    for err in &outcome.diff_details {
                        println!("  {err}");
                    }
                }
            }
        }

        if is_equivalent {
            info!("BAM groupings are equivalent (key-join)");
            Ok(outcome.bam1_count + outcome.bam2_count)
        } else {
            info!("BAM groupings differ (key-join)");
            Err(super::CompareMismatch("BAM groupings differ (order-independent)".to_owned())
                .into())
        }
    }

    /// Resolve this run's [`KeyJoinConfig`] from `--threads`/`--sort-memory`/
    /// `--sort-tmp-dir`/`--max-diffs`.
    ///
    /// `--sort-memory` is a **total** budget, not per-thread (`resolve_memory_budget`'s
    /// `per_thread = false`), matching the flag's documented semantics. `sort_tmp_dirs` is
    /// resolved via [`keyjoin::resolve_sort_tmp_dirs`], which never returns empty (falling
    /// back to a disk-backed default rather than tmpfs) — so both the canonicalization
    /// sort's spill chunks and the temp canonicalized BAMs it writes always land on the
    /// same resolved, disk-backed directory.
    ///
    /// # Errors
    ///
    /// Returns an error if the memory budget cannot be resolved (e.g. `--threads 0`).
    fn keyjoin_config(&self) -> Result<KeyJoinConfig> {
        let sort_memory =
            resolve_memory_budget(self.sort_memory, MemoryReserve::Auto, self.threads, false)?;
        let sort_tmp_dirs = keyjoin::resolve_sort_tmp_dirs(
            &self.sort_tmp_dirs,
            std::env::var(TMP_DIRS_ENV).ok().as_deref(),
        );
        Ok(KeyJoinConfig {
            threads: self.threads,
            sort_memory,
            sort_tmp_dirs,
            max_diffs: self.max_diffs,
        })
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use clap::Parser;
    use rstest::rstest;

    fn parse(args: &[&str]) -> CompareBams {
        let mut argv = vec!["bams", "a.bam", "b.bam"];
        argv.extend_from_slice(args);
        CompareBams::try_parse_from(argv).expect("parse")
    }

    // ---- count_base_pairing_mismatches (duplex strand-pairing check) --------

    use super::record_key::Segment;

    /// Build a distinct [`RecordKey`] for use as a test map key; only the
    /// `name` differs between ids, which is sufficient for uniqueness here.
    fn test_key(id: u8) -> RecordKey {
        RecordKey {
            name: vec![id],
            segment: Segment::Fragment,
            secondary: false,
            supplementary: false,
            multimap_locus: None,
        }
    }

    /// A strand-pairing split must be flagged: BAM1 pairs both strands into
    /// molecule `0` (`0/A` + `0/B`); BAM2 splits those same reads into two
    /// molecules (`0/A` + `1/A`). The full-MiKey check sees a consistent
    /// bijection (0 mismatches), so this base-level check is what catches it.
    #[test]
    fn count_base_pairing_mismatches_flags_strand_split() {
        let mut m1: AHashMap<RecordKey, MiKey> = AHashMap::new();
        m1.insert(test_key(1), MiKey::Strand { base: 0, strand: b'A' });
        m1.insert(test_key(2), MiKey::Strand { base: 0, strand: b'B' });
        let mut m2: AHashMap<RecordKey, MiKey> = AHashMap::new();
        m2.insert(test_key(1), MiKey::Strand { base: 0, strand: b'A' });
        m2.insert(test_key(2), MiKey::Strand { base: 1, strand: b'A' });
        assert!(
            count_base_pairing_mismatches(&m1, &m2) > 0,
            "a duplex molecule split across two bases must be flagged"
        );
    }

    /// Pure molecule renumbering (same pairing, different base id) is
    /// equivalent and must not be flagged.
    #[test]
    fn count_base_pairing_mismatches_zero_for_molecule_relabel() {
        let mut m1: AHashMap<RecordKey, MiKey> = AHashMap::new();
        m1.insert(test_key(1), MiKey::Strand { base: 0, strand: b'A' });
        m1.insert(test_key(2), MiKey::Strand { base: 0, strand: b'B' });
        let mut m2: AHashMap<RecordKey, MiKey> = AHashMap::new();
        m2.insert(test_key(1), MiKey::Strand { base: 9, strand: b'A' });
        m2.insert(test_key(2), MiKey::Strand { base: 9, strand: b'B' });
        assert_eq!(count_base_pairing_mismatches(&m1, &m2), 0);
    }

    /// Swapping the `/A` and `/B` labels of a molecule is a naming difference,
    /// not a pairing difference — equivalent at the base level.
    #[test]
    fn count_base_pairing_mismatches_zero_for_ab_strand_swap() {
        let mut m1: AHashMap<RecordKey, MiKey> = AHashMap::new();
        m1.insert(test_key(1), MiKey::Strand { base: 0, strand: b'A' });
        m1.insert(test_key(2), MiKey::Strand { base: 0, strand: b'B' });
        let mut m2: AHashMap<RecordKey, MiKey> = AHashMap::new();
        m2.insert(test_key(1), MiKey::Strand { base: 0, strand: b'B' });
        m2.insert(test_key(2), MiKey::Strand { base: 0, strand: b'A' });
        assert_eq!(count_base_pairing_mismatches(&m1, &m2), 0);
    }

    /// For single-strand (`Int`) MIs `base()` is the id itself, so a consistent
    /// relabel reduces to the ordinary partition check and adds no false
    /// mismatches.
    #[test]
    fn count_base_pairing_mismatches_zero_for_simplex_relabel() {
        let mut m1: AHashMap<RecordKey, MiKey> = AHashMap::new();
        m1.insert(test_key(1), MiKey::Int(5));
        m1.insert(test_key(2), MiKey::Int(5));
        let mut m2: AHashMap<RecordKey, MiKey> = AHashMap::new();
        m2.insert(test_key(1), MiKey::Int(7));
        m2.insert(test_key(2), MiKey::Int(7));
        assert_eq!(count_base_pairing_mismatches(&m1, &m2), 0);
    }

    #[rstest]
    #[case(CommandPreset::Extract)]
    #[case(CommandPreset::Zipper)]
    #[case(CommandPreset::Sort)]
    #[case(CommandPreset::Correct)]
    #[case(CommandPreset::Dedup)]
    #[case(CommandPreset::Filter)]
    fn preset_defaults_content_stages_map_to_content_no_ignore_order(#[case] stage: CommandPreset) {
        let (mode, ignore) = stage.defaults();
        assert!(matches!(mode, CompareMode::Content), "{stage:?} → {mode:?}");
        assert!(!ignore, "{stage:?} → ignore_order {ignore}");
    }

    /// Consensus presets (simplex/duplex/codec) reroute to positional `Content`
    /// comparison under the `ExactConsensus` predicate (see
    /// `execute()`'s predicate selection) rather than MI-grouping equivalence — this
    /// is a breaking strictness change (see the commit message and the
    /// compare-hardening design spec's §"Accepted divergences").
    #[rstest]
    #[case(CommandPreset::Simplex)]
    #[case(CommandPreset::Duplex)]
    #[case(CommandPreset::Codec)]
    fn consensus_presets_are_content_exact(#[case] stage: CommandPreset) {
        let (mode, ignore) = stage.defaults();
        assert!(matches!(mode, CompareMode::Content), "{stage:?} → {mode:?}, must be positional");
        assert!(!ignore, "{stage:?} → ignore_order {ignore}");
    }

    #[test]
    fn preset_defaults_grouping_stage_maps_to_grouping_with_ignore_order() {
        let (mode, ignore) = CommandPreset::Group.defaults();
        assert!(matches!(mode, CompareMode::Grouping), "Group → {mode:?}");
        assert!(ignore, "Group → ignore_order {ignore}");
    }

    /// `Group`'s content predicate must be `ExactMinusMi` (Task 3.4): the key-join engine
    /// checks everything except the MI tag as content, and checks the MI tag separately via
    /// the fgumi-MI/fgbio-MI bijection.
    #[test]
    fn group_preset_content_predicate_is_exact_minus_mi() {
        assert_eq!(CommandPreset::Group.content_predicate(), ContentPredicate::ExactMinusMi);
    }

    /// `filter` output is consensus reads that still carry the depth tags (`cD`/`cM`/`cE`,
    /// and duplex `aD`/`aM`/`bD`/`bM`/`aE`/`bE`) as the consensus command that produced them
    /// -- `filter` only drops reads, it never rewrites these tags. So `Filter` routes through
    /// the same `ExactConsensus` predicate as `Simplex`/`Duplex`/`Codec`, keeping the four
    /// consensus commands on one predicate.
    #[test]
    fn filter_preset_content_predicate_is_exact_consensus() {
        assert_eq!(CommandPreset::Filter.content_predicate(), ContentPredicate::ExactConsensus);
    }

    #[test]
    fn effective_settings_fall_back_to_built_in_defaults() {
        let args = parse(&[]);
        let (mode, ignore) = args.effective_settings();
        assert!(matches!(mode, CompareMode::Full));
        assert!(!ignore);
    }

    #[test]
    fn effective_settings_use_preset_when_only_command_given() {
        let args = parse(&["--command", "group"]);
        let (mode, ignore) = args.effective_settings();
        assert!(matches!(mode, CompareMode::Grouping));
        assert!(ignore);

        let args = parse(&["--command", "extract"]);
        let (mode, ignore) = args.effective_settings();
        assert!(matches!(mode, CompareMode::Content));
        assert!(!ignore);
    }

    #[test]
    fn effective_settings_explicit_mode_overrides_preset() {
        // --mode full overrides --command group's Grouping default
        let args = parse(&["--command", "group", "--mode", "full"]);
        let (mode, ignore) = args.effective_settings();
        assert!(matches!(mode, CompareMode::Full));
        // Preset-inherited ignore_order=true is silently dropped because the
        // resolved mode isn't Grouping; without this, the later bail would
        // blame the user for an --ignore-order flag they never passed.
        assert!(!ignore);
    }

    #[test]
    fn effective_settings_drops_preset_ignore_order_when_explicit_mode_not_grouping() {
        // --command simplex now presets (Content, false) itself, so this exercises the
        // more general case via --command group (Grouping, true) narrowed by an
        // explicit --mode content, proving ignore_order=true never leaks through when
        // the resolved mode isn't Grouping.
        let args = parse(&["--command", "group", "--mode", "content"]);
        let (mode, ignore) = args.effective_settings();
        assert!(matches!(mode, CompareMode::Content));
        assert!(!ignore);
    }

    #[test]
    fn effective_settings_explicit_ignore_order_overrides_preset() {
        // --ignore-order=false overrides --command group's true default
        let args = parse(&["--command", "group", "--ignore-order", "false"]);
        let (mode, ignore) = args.effective_settings();
        // mode still inherits from preset (Grouping)
        assert!(matches!(mode, CompareMode::Grouping));
        assert!(!ignore);
    }

    // ==================== MI tag extraction tests ====================

    /// Build a minimal unmapped `RawRecord` carrying only the supplied aux bytes.
    /// A 3-byte read name keeps the name+NUL length word-aligned for CIGAR.
    fn raw_record_with_aux(aux: &[u8]) -> RawRecord {
        let bytes = fgumi_raw_bam::make_bam_bytes(-1, -1, 4, b"rea", &[], 0, -1, -1, aux);
        RawRecord::from(bytes)
    }

    fn aux_mi_string(payload: &[u8]) -> Vec<u8> {
        let mut aux = vec![b'M', b'I', b'Z'];
        aux.extend_from_slice(payload);
        aux.push(0);
        aux
    }

    fn aux_mi_i32(v: i32) -> Vec<u8> {
        let mut aux = vec![b'M', b'I', b'i'];
        aux.extend_from_slice(&v.to_le_bytes());
        aux
    }

    #[test]
    fn get_mi_tag_raw_parses_integer_string_form() {
        let rec = raw_record_with_aux(&aux_mi_string(b"42"));
        assert_eq!(get_mi_tag_raw(&rec), Some(MiKey::Int(42)));
    }

    #[test]
    fn get_mi_tag_raw_parses_integer_aux_form() {
        let rec = raw_record_with_aux(&aux_mi_i32(42));
        assert_eq!(get_mi_tag_raw(&rec), Some(MiKey::Int(42)));
    }

    #[test]
    fn get_mi_tag_raw_parses_paired_a_suffix() {
        let rec = raw_record_with_aux(&aux_mi_string(b"0/A"));
        assert_eq!(get_mi_tag_raw(&rec), Some(MiKey::Strand { base: 0, strand: b'A' }));
    }

    #[test]
    fn get_mi_tag_raw_parses_paired_b_suffix() {
        let rec = raw_record_with_aux(&aux_mi_string(b"7/B"));
        assert_eq!(get_mi_tag_raw(&rec), Some(MiKey::Strand { base: 7, strand: b'B' }));
    }

    #[test]
    fn get_mi_tag_raw_distinguishes_paired_a_from_b_with_same_base() {
        // The essential invariant: a /A and /B grouping of the same molecule
        // must be different keys so downstream duplex strand checks work.
        let a = raw_record_with_aux(&aux_mi_string(b"13/A"));
        let b = raw_record_with_aux(&aux_mi_string(b"13/B"));
        assert_ne!(get_mi_tag_raw(&a), get_mi_tag_raw(&b));
    }

    #[test]
    fn get_mi_tag_raw_missing_tag_returns_none() {
        let rec = raw_record_with_aux(&[]);
        assert_eq!(get_mi_tag_raw(&rec), None);
    }

    #[test]
    fn get_mi_tag_raw_rejects_non_integer_non_strand_string() {
        let rec = raw_record_with_aux(&aux_mi_string(b"nope"));
        assert_eq!(get_mi_tag_raw(&rec), None);
    }

    #[test]
    fn get_mi_tag_raw_rejects_unknown_strand_suffix() {
        let rec = raw_record_with_aux(&aux_mi_string(b"5/C"));
        assert_eq!(get_mi_tag_raw(&rec), None);
    }

    #[test]
    fn mikey_display_matches_bam_encoding() {
        assert_eq!(MiKey::Int(42).to_string(), "42");
        assert_eq!(MiKey::Strand { base: 3, strand: b'A' }.to_string(), "3/A");
        assert_eq!(MiKey::Strand { base: 3, strand: b'B' }.to_string(), "3/B");
    }
}
