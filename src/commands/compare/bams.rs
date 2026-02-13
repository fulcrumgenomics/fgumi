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

use ahash::{AHashMap, AHashSet};
use anyhow::{Result, anyhow, bail};
use clap::{Parser, ValueEnum};
use crossbeam_channel::{Receiver, bounded};
use fgumi_lib::bam_io::{BamReaderAuto, create_bam_reader};
use fgumi_lib::logging::OperationTimer;
use fgumi_lib::progress::ProgressTracker;
use fgumi_lib::validation::validate_file_exists;
use itertools::Itertools;
use log::info;
use noodles::sam::Header;
use noodles::sam::alignment::record::Cigar as CigarTrait;
use noodles::sam::alignment::record::cigar::op::Kind;
use noodles::sam::alignment::record::data::field::Tag;
use noodles::sam::alignment::record_buf::RecordBuf;
use noodles::sam::alignment::record_buf::data::field::Value;
use noodles::sam::alignment::record_buf::data::field::value::Array;
use rayon::prelude::*;
use std::collections::{BTreeMap, HashSet};
use std::path::{Path, PathBuf};
use std::thread;

use crate::commands::command::Command;

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
    Both files MUST be in the same order (e.g., query-name sorted with `fgumi sort --order queryname`).
    Validates that:
    1. Read names and R1/R2 flags match between files
    2. Reads with the same MI in file 1 have the same MI in file 2 (and vice versa)
    Does NOT compare other BAM content (sequence, quality, other tags).
    This proves the grouping is semantically equivalent even if MI values differ.

RECOMMENDED SETTINGS BY COMMAND:

  When comparing output from different fgumi commands, use these settings:

  Command         --mode      --ignore-order   Notes
  ─────────────────────────────────────────────────────────────────────────
  extract         content     false            No MI tags; deterministic
  zipper          content     false            Preserves MI tags unchanged
  group           full        false            Has MI tags; deterministic
  simplex         grouping    true             Non-deterministic with --threads
  duplex          grouping    true             Non-deterministic with --threads
  codec           grouping    true             Non-deterministic with --threads
  filter          content     false            Passes through MI tags unchanged
  clip            content     false            Does not modify MI tags
  correct         content     false            Modifies RX tag only, not MI
  downsample      content     false            Deterministic with seed
  review          content     false            Preserves MI tags

  For simulate subcommands:
  mapped-reads    content     false            Template-coordinate sorted
  grouped-reads   full        false            Has MI tags; deterministic
  consensus-reads content     false            Unmapped with consensus tags

  Examples for each command:

    # extract output (no MI tags)
    fgumi compare bams extracted1.bam extracted2.bam --mode content

    # group output (has MI tags)
    fgumi compare bams grouped1.bam grouped2.bam --mode full

    # simplex/duplex/codec output (MI grouping, non-deterministic order)
    fgumi compare bams consensus1.bam consensus2.bam --mode grouping --ignore-order

    # filter/clip/correct/downsample output (preserves content)
    fgumi compare bams filtered1.bam filtered2.bam --mode content

Example usage:
  fgumi compare bams bam1.bam bam2.bam                    # full mode (default)
  fgumi compare bams bam1.bam bam2.bam --mode content    # content only
  fgumi compare bams bam1.bam bam2.bam --mode grouping   # grouping only
  fgumi compare bams bam1.bam bam2.bam --mode grouping --ignore-order  # consensus output
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

    /// Comparison mode: 'full' (MI grouping + content, for group output),
    /// 'content' (all fields, for extract/filter/clip/correct/downsample output),
    /// 'grouping' (MI equivalence only, for simplex/duplex/codec output)
    #[arg(long = "mode", default_value = "full")]
    pub mode: CompareMode,

    /// Maximum number of differences to report in detail
    #[arg(short = 'm', long = "max-diffs", default_value = "10")]
    pub max_diffs: usize,

    /// Quiet mode - only exit code indicates result (0=equal, 1=different)
    #[arg(short = 'q', long = "quiet")]
    pub quiet: bool,

    /// Ignore record order when comparing in grouping mode.
    /// Required for comparing output from consensus commands (simplex/duplex/codec)
    /// when run with --threads, as parallel processing causes non-deterministic ordering.
    /// Only valid with --mode grouping.
    #[arg(long = "ignore-order", default_value = "false")]
    pub ignore_order: bool,

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
const FIELD_NAMES: [&str; 11] =
    ["QNAME", "FLAG", "RNAME", "POS", "MAPQ", "CIGAR", "RNEXT", "PNEXT", "TLEN", "SEQ", "QUAL"];

/// Convert CIGAR op kind to its SAM character representation.
fn cigar_kind_to_char(kind: Kind) -> char {
    match kind {
        Kind::Match => 'M',
        Kind::Insertion => 'I',
        Kind::Deletion => 'D',
        Kind::Skip => 'N',
        Kind::SoftClip => 'S',
        Kind::HardClip => 'H',
        Kind::Pad => 'P',
        Kind::SequenceMatch => '=',
        Kind::SequenceMismatch => 'X',
    }
}

/// Format CIGAR as a string.
fn format_cigar(record: &RecordBuf) -> String {
    let cigar = record.cigar();
    let ops: Vec<String> = cigar
        .iter()
        .flatten()
        .map(|op| format!("{}{}", op.len(), cigar_kind_to_char(op.kind())))
        .collect();
    if ops.is_empty() { "*".to_string() } else { ops.join("") }
}

/// Format sequence as a string.
fn format_sequence(record: &RecordBuf) -> String {
    let seq = record.sequence();
    if seq.is_empty() {
        "*".to_string()
    } else {
        seq.as_ref()
            .iter()
            .map(|&b| match b {
                b'A' | b'a' => 'A',
                b'C' | b'c' => 'C',
                b'G' | b'g' => 'G',
                b'T' | b't' => 'T',
                _ => 'N',
            })
            .collect()
    }
}

/// Extract core SAM fields as comparable strings.
fn get_core_fields(record: &RecordBuf, header: &noodles::sam::Header) -> [String; 11] {
    let qname = record_name_to_string(record);
    let flag = record.flags().bits().to_string();

    let rname = record
        .reference_sequence_id()
        .and_then(|id| header.reference_sequences().get_index(id).map(|(name, _)| name.to_string()))
        .unwrap_or_else(|| "*".to_string());

    // Position is 0-based internally, report as 1-based
    let pos = record.alignment_start().map_or_else(|| "0".to_string(), |p| p.get().to_string());

    let mapq = record.mapping_quality().map_or_else(|| "255".to_string(), |q| q.get().to_string());

    let cigar = format_cigar(record);

    let rnext = record
        .mate_reference_sequence_id()
        .and_then(|id| header.reference_sequences().get_index(id).map(|(name, _)| name.to_string()))
        .unwrap_or_else(|| "*".to_string());

    let pnext =
        record.mate_alignment_start().map_or_else(|| "0".to_string(), |p| p.get().to_string());

    let tlen = record.template_length().to_string();

    let seq = format_sequence(record);

    let qual = if record.quality_scores().is_empty() {
        "*".to_string()
    } else {
        record.quality_scores().as_ref().iter().map(|q| (q + 33) as char).collect()
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

/// Get the MI tag value from a record.
fn get_mi_tag(record: &RecordBuf) -> Option<String> {
    let mi_tag = Tag::from([b'M', b'I']);
    record.data().get(&mi_tag).map(format_tag_value)
}

/// Convert a record's name to a String, returning "*" if missing.
fn record_name_to_string(record: &RecordBuf) -> String {
    record.name().map_or_else(|| "*".to_string(), std::string::ToString::to_string)
}

// ============================================================================
// Zero-copy comparison helpers
// ============================================================================

/// Extract the numeric value from an integer Value as i64 for semantic comparison.
fn value_to_i64(v: &Value) -> Option<i64> {
    match v {
        Value::Int8(i) => Some(i64::from(*i)),
        Value::UInt8(i) => Some(i64::from(*i)),
        Value::Int16(i) => Some(i64::from(*i)),
        Value::UInt16(i) => Some(i64::from(*i)),
        Value::Int32(i) => Some(i64::from(*i)),
        Value::UInt32(i) => Some(i64::from(*i)),
        _ => None,
    }
}

/// Compare two tag values without allocating strings (unless they differ).
/// For floats, we use bit-exact comparison since BAM comparison requires exact equality.
/// For integers, we compare semantically (by numeric value, not storage type).
#[allow(clippy::float_cmp)]
fn tag_values_equal(v1: &Value, v2: &Value) -> bool {
    match (v1, v2) {
        (Value::Character(a), Value::Character(b)) => a == b,
        // For integers, compare by numeric value regardless of storage type
        (
            Value::Int8(_)
            | Value::UInt8(_)
            | Value::Int16(_)
            | Value::UInt16(_)
            | Value::Int32(_)
            | Value::UInt32(_),
            Value::Int8(_)
            | Value::UInt8(_)
            | Value::Int16(_)
            | Value::UInt16(_)
            | Value::Int32(_)
            | Value::UInt32(_),
        ) => value_to_i64(v1) == value_to_i64(v2),
        (Value::Float(a), Value::Float(b)) => a == b,
        (Value::String(a), Value::String(b)) => a == b,
        (Value::Hex(a), Value::Hex(b)) => a == b,
        (Value::Array(a), Value::Array(b)) => arrays_equal(a, b),
        _ => false, // Different types (e.g., string vs integer)
    }
}

/// Extract array elements as i64 for semantic comparison.
fn array_to_i64_vec(a: &Array) -> Option<Vec<i64>> {
    match a {
        Array::Int8(v) => Some(v.iter().map(|&i| i64::from(i)).collect()),
        Array::UInt8(v) => Some(v.iter().map(|&i| i64::from(i)).collect()),
        Array::Int16(v) => Some(v.iter().map(|&i| i64::from(i)).collect()),
        Array::UInt16(v) => Some(v.iter().map(|&i| i64::from(i)).collect()),
        Array::Int32(v) => Some(v.iter().map(|&i| i64::from(i)).collect()),
        Array::UInt32(v) => Some(v.iter().map(|&i| i64::from(i)).collect()),
        Array::Float(_) => None,
    }
}

/// Check if an array is an integer array type.
fn is_integer_array(a: &Array) -> bool {
    matches!(
        a,
        Array::Int8(_)
            | Array::UInt8(_)
            | Array::Int16(_)
            | Array::UInt16(_)
            | Array::Int32(_)
            | Array::UInt32(_)
    )
}

/// Compare two arrays for equality.
/// For integer arrays, we compare semantically (by numeric values, not storage type).
/// For float arrays, we require exact type match.
#[allow(clippy::float_cmp)]
fn arrays_equal(a: &Array, b: &Array) -> bool {
    // For integer arrays, compare semantically
    if is_integer_array(a) && is_integer_array(b) {
        return array_to_i64_vec(a) == array_to_i64_vec(b);
    }
    // For float arrays, require exact type match
    match (a, b) {
        (Array::Float(va), Array::Float(vb)) => va == vb,
        _ => false, // Different types or mixed int/float
    }
}

/// Lazy tag comparison: returns true if tags are equal (values only, order independent).
/// Avoids `BTreeMap` allocation when tags match.
fn tags_equal_lazy(record1: &RecordBuf, record2: &RecordBuf) -> bool {
    let data1 = record1.data();
    let data2 = record2.data();

    // Quick length check first
    if data1.len() != data2.len() {
        return false;
    }

    // Compare each tag from record1 against record2
    for (tag, val1) in data1.iter() {
        match data2.get(&tag) {
            Some(val2) if tag_values_equal(val1, val2) => continue,
            _ => return false,
        }
    }
    true
}

/// Check if tag order matches (assuming values already match).
fn tags_order_matches(record1: &RecordBuf, record2: &RecordBuf) -> bool {
    let iter1 = record1.data().iter();
    let iter2 = record2.data().iter();

    for ((tag1, _), (tag2, _)) in iter1.zip(iter2) {
        if tag1 != tag2 {
            return false;
        }
    }
    true
}

/// Zero-copy comparison of core fields. Returns true if all core fields match.
fn core_fields_equal(record1: &RecordBuf, record2: &RecordBuf) -> bool {
    // Compare in order of likelihood to differ / cost to compare
    // Flags are cheap and often differ
    if record1.flags() != record2.flags() {
        return false;
    }

    // Names - compare Option<&BStr> directly
    if record1.name() != record2.name() {
        return false;
    }

    // Reference sequence ID
    if record1.reference_sequence_id() != record2.reference_sequence_id() {
        return false;
    }

    // Position
    if record1.alignment_start() != record2.alignment_start() {
        return false;
    }

    // Mapping quality
    if record1.mapping_quality() != record2.mapping_quality() {
        return false;
    }

    // Mate reference sequence ID
    if record1.mate_reference_sequence_id() != record2.mate_reference_sequence_id() {
        return false;
    }

    // Mate position
    if record1.mate_alignment_start() != record2.mate_alignment_start() {
        return false;
    }

    // Template length
    if record1.template_length() != record2.template_length() {
        return false;
    }

    // Sequence - compare bytes directly
    if record1.sequence().as_ref() != record2.sequence().as_ref() {
        return false;
    }

    // Quality scores
    if record1.quality_scores().as_ref() != record2.quality_scores().as_ref() {
        return false;
    }

    // CIGAR - compare ops
    let cigar1 = record1.cigar();
    let cigar2 = record2.cigar();
    let ops1: Vec<_> = cigar1.iter().flatten().collect();
    let ops2: Vec<_> = cigar2.iter().flatten().collect();
    ops1 == ops2
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
    read_key: ReadKey,
    mi1: Option<String>,
    mi2: Option<String>,
    name_match: bool,
    flag_match: bool,
    diff_detail: Option<DiffDetail>,
}

// ============================================================================
// Types and MI map helpers (using ahash)
// ============================================================================

/// A unique identifier for a read: (`read_name`, `is_read1`)
/// For display purposes only - the actual key is a hash
type ReadKey = (String, bool);

/// Compact read key using hash - saves ~70 bytes per entry vs (String, bool)
type ReadKeyHash = u64;

/// Compute a hash for a read key (`read_name`, `is_read1`)
#[inline]
fn hash_read_key(name: &str, is_read1: bool) -> ReadKeyHash {
    use std::hash::{Hash, Hasher};
    let mut hasher = ahash::AHasher::default();
    name.hash(&mut hasher);
    is_read1.hash(&mut hasher);
    hasher.finish()
}

/// Build a map from MI value to set of read key hashes.
fn build_mi_groups_compact(
    mi_map: &AHashMap<ReadKeyHash, i64>,
) -> AHashMap<i64, AHashSet<ReadKeyHash>> {
    let mut groups: AHashMap<i64, AHashSet<ReadKeyHash>> = AHashMap::new();
    for (read_key_hash, mi) in mi_map {
        groups.entry(*mi).or_default().insert(*read_key_hash);
    }
    groups
}

/// Build a map from MI value to set of read keys that have that MI.
/// (Used by full mode which needs full keys for error reporting)
fn build_mi_groups<'a>(
    mi_map: &'a AHashMap<ReadKey, String>,
) -> AHashMap<&'a String, AHashSet<&'a ReadKey>> {
    let mut groups: AHashMap<&'a String, AHashSet<&'a ReadKey>> = AHashMap::new();
    for (read_key, mi) in mi_map {
        groups.entry(mi).or_default().insert(read_key);
    }
    groups
}

/// Result from parallel MI extraction for a single record
struct MiExtractResult {
    key_hash: ReadKeyHash,
    mi: Option<i64>,
}

/// Build an MI map from a BAM file using parallel batch processing.
/// Returns a map from read key hash to MI value, along with total record count.
fn build_mi_map_parallel(
    path: &Path,
    threads: usize,
    batch_size: usize,
) -> Result<(AHashMap<ReadKeyHash, i64>, u64)> {
    let (rx, _header) = start_batch_reader(path.to_path_buf(), threads, batch_size)?;

    let mut mi_map: AHashMap<ReadKeyHash, i64> = AHashMap::new();
    let mut total_records: u64 = 0;

    loop {
        match rx.recv() {
            Ok(BatchMessage::Batch(batch)) => {
                // Extract MI values in parallel using rayon
                let results: Vec<MiExtractResult> = batch
                    .par_iter()
                    .map(|record| {
                        let name = record_name_to_string(record);
                        let is_read1 = record.flags().is_first_segment();
                        let key_hash = hash_read_key(&name, is_read1);

                        // Parse MI tag as i64 directly
                        let mi = get_mi_tag(record).and_then(|s| s.parse::<i64>().ok());

                        MiExtractResult { key_hash, mi }
                    })
                    .collect();

                // Insert into map (sequential, but fast)
                for r in results {
                    total_records += 1;
                    if let Some(mi_val) = r.mi {
                        mi_map.insert(r.key_hash, mi_val);
                    }
                }
            }
            Ok(BatchMessage::Eof) => break,
            Ok(BatchMessage::Error(e)) => bail!("Error reading BAM: {e}"),
            Err(_) => break, // Channel closed
        }
    }

    Ok((mi_map, total_records))
}

/// Statistics for grouping comparison mode.
#[derive(Debug, Default)]
struct GroupingStats {
    total_records: u64,
    order_mismatches: u64,
    missing_mi_bam1: u64,
    missing_mi_bam2: u64,
    grouping_mismatches: u64,
    unique_groups_bam1: usize,
    unique_groups_bam2: usize,
}

/// Statistics for unordered grouping comparison mode.
#[derive(Debug, Default)]
struct UnorderedGroupingStats {
    total_bam1: u64,
    total_bam2: u64,
    matched: u64,
    only_in_bam1: u64,
    only_in_bam2: u64,
    grouping_mismatches: u64,
    unique_groups_bam1: usize,
    unique_groups_bam2: usize,
}

// ============================================================================
// Double-buffered batch reading
// ============================================================================

/// Message type for the double-buffered reader channel
enum BatchMessage {
    /// A batch of records
    Batch(Vec<RecordBuf>),
    /// End of file reached
    Eof,
    /// Error occurred during reading
    Error(String),
}

/// Reads a single batch of records from a BAM reader.
fn read_batch(
    reader: &mut BamReaderAuto,
    header: &Header,
    batch_size: usize,
) -> std::io::Result<(Vec<RecordBuf>, bool)> {
    let mut batch = Vec::with_capacity(batch_size);
    let mut record = RecordBuf::default();

    for _ in 0..batch_size {
        if reader.read_record_buf(header, &mut record)? == 0 {
            return Ok((batch, true)); // EOF
        }
        batch.push(std::mem::take(&mut record));
    }
    Ok((batch, false))
}

/// Starts a background reader thread that sends batches through a channel.
/// Returns a receiver for the batches.
fn start_batch_reader(
    path: PathBuf,
    threads: usize,
    batch_size: usize,
) -> Result<(Receiver<BatchMessage>, Header)> {
    // Open the reader on the main thread to get the header
    let (mut reader, header) = create_bam_reader(&path, threads)?;
    let header_clone = header.clone();

    // Create a bounded channel (double buffering = 2 slots)
    let (tx, rx) = bounded::<BatchMessage>(2);

    // Spawn the reader thread
    thread::spawn(move || {
        loop {
            match read_batch(&mut reader, &header_clone, batch_size) {
                Ok((batch, eof)) => {
                    if !batch.is_empty() && tx.send(BatchMessage::Batch(batch)).is_err() {
                        break; // Receiver dropped
                    }
                    if eof {
                        let _ = tx.send(BatchMessage::Eof);
                        break;
                    }
                }
                Err(e) => {
                    let _ = tx.send(BatchMessage::Error(e.to_string()));
                    break;
                }
            }
        }
    });

    Ok((rx, header))
}

/// Compares a batch of record pairs in parallel using rayon.
fn compare_batch_parallel(
    batch1: &[RecordBuf],
    batch2: &[RecordBuf],
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

            // Use zero-copy comparison first
            let core_match = core_fields_equal(r1, r2);

            // Use lazy tag comparison
            let tags_match = tags_equal_lazy(r1, r2);
            let tag_order_match = if tags_match { tags_order_matches(r1, r2) } else { false };

            // Collect detailed diff if there's a mismatch.
            // Note: We collect all diffs here (regardless of max_diffs) because the batch
            // index `i` doesn't reflect actual diff count. We truncate to max_diffs later
            // when aggregating results, as diffs are rare and this has minimal overhead.
            let diff_detail = if !core_match {
                let core1 = get_core_fields(r1, header1);
                let core2 = get_core_fields(r2, header2);
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
            } else if !tags_match {
                let tags1 = get_tags_map(r1);
                let tags2 = get_tags_map(r2);
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
                let qname = record_name_to_string(r1);
                Some(DiffDetail {
                    record_num,
                    qname,
                    flags: r1.flags().bits().to_string(),
                    diff_type: DiffType::TagDiff,
                    diffs,
                })
            } else {
                None
            };

            RecordCompareResult { core_match, tags_match, tag_order_match, diff_detail }
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

/// Compares a batch of record pairs for grouping mode in parallel.
fn compare_batch_grouping_parallel(
    batch1: &[RecordBuf],
    batch2: &[RecordBuf],
    start_index: u64,
    max_diffs: usize,
    existing_diffs: usize,
) -> Vec<GroupingCompareResult> {
    batch1
        .par_iter()
        .zip(batch2.par_iter())
        .enumerate()
        .map(|(i, (r1, r2))| {
            let record_num = start_index + i as u64 + 1;
            let qname1 = record_name_to_string(r1);
            let qname2 = record_name_to_string(r2);
            let is_read1_r1 = r1.flags().is_first_segment();
            let is_read1_r2 = r2.flags().is_first_segment();

            let name_match = qname1 == qname2;
            let flag_match = is_read1_r1 == is_read1_r2;

            let read_key = (qname1.clone(), is_read1_r1);
            let mi1 = get_mi_tag(r1);
            let mi2 = get_mi_tag(r2);

            let diff_detail = if !name_match && existing_diffs + i < max_diffs {
                Some(DiffDetail {
                    record_num,
                    qname: qname1.clone(),
                    flags: r1.flags().bits().to_string(),
                    diff_type: DiffType::ReadNameMismatch,
                    diffs: vec![format!("Read names differ: '{}' vs '{}'", qname1, qname2)],
                })
            } else if name_match && !flag_match && existing_diffs + i < max_diffs {
                Some(DiffDetail {
                    record_num,
                    qname: qname1.clone(),
                    flags: r1.flags().bits().to_string(),
                    diff_type: DiffType::FlagMismatch,
                    diffs: vec![format!(
                        "R1/R2 flags differ: is_read1={} vs is_read1={}",
                        is_read1_r1, is_read1_r2
                    )],
                })
            } else {
                None
            };

            GroupingCompareResult {
                record_num,
                read_key,
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

        let timer = OperationTimer::new("Comparing BAMs");

        let total_records = match self.mode {
            CompareMode::Full => self.execute_full()?,
            CompareMode::Content => self.execute_content()?,
            CompareMode::Grouping => self.execute_grouping()?,
        };

        timer.log_completion(total_records);
        Ok(())
    }
}

impl CompareBams {
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

    /// Execute content comparison mode
    /// Compares all BAM fields record-by-record without MI grouping analysis.
    /// Uses parallel batch processing with double buffering for performance.
    fn execute_content(&self) -> Result<u64> {
        let mut stats = CompareStats::default();
        let batch_size = self.batch_size;

        info!(
            "Starting content comparison with {} threads, batch size {}",
            self.threads, batch_size
        );

        // Start double-buffered readers for both BAM files
        let (rx1, header1) = start_batch_reader(self.bam1.clone(), self.threads, batch_size)?;
        let (rx2, header2) = start_batch_reader(self.bam2.clone(), self.threads, batch_size)?;

        // Progress tracking
        let progress = ProgressTracker::new("Processed records").with_interval(1_000_000);

        // Process batches
        let mut bam1_eof = false;
        let mut bam2_eof = false;
        let mut pending_batch1: Option<Vec<RecordBuf>> = None;
        let mut pending_batch2: Option<Vec<RecordBuf>> = None;
        let mut current_index = 0u64;

        loop {
            // Get next batch from BAM1 if needed
            if pending_batch1.is_none() && !bam1_eof {
                match rx1.recv() {
                    Ok(BatchMessage::Batch(batch)) => {
                        stats.bam1_count += batch.len() as u64;
                        pending_batch1 = Some(batch);
                    }
                    Ok(BatchMessage::Eof) => bam1_eof = true,
                    Ok(BatchMessage::Error(e)) => bail!("Error reading BAM1: {e}"),
                    Err(_) => bam1_eof = true,
                }
            }

            // Get next batch from BAM2 if needed
            if pending_batch2.is_none() && !bam2_eof {
                match rx2.recv() {
                    Ok(BatchMessage::Batch(batch)) => {
                        stats.bam2_count += batch.len() as u64;
                        pending_batch2 = Some(batch);
                    }
                    Ok(BatchMessage::Eof) => bam2_eof = true,
                    Ok(BatchMessage::Error(e)) => bail!("Error reading BAM2: {e}"),
                    Err(_) => bam2_eof = true,
                }
            }

            // Check for completion
            match (&pending_batch1, &pending_batch2) {
                (None, None) => break,
                (Some(_), None) | (None, Some(_)) => {
                    // One file exhausted before the other
                    if stats.diff_details.len() < self.max_diffs {
                        stats.diff_details.push(DiffDetail {
                            record_num: current_index,
                            qname: "N/A".to_string(),
                            flags: "N/A".to_string(),
                            diff_type: DiffType::CountMismatch,
                            diffs: vec!["BAM files have different number of records".to_string()],
                        });
                    }
                    // Drain remaining batches to get accurate counts
                    if pending_batch1.is_some() {
                        while let Ok(msg) = rx1.recv() {
                            if let BatchMessage::Batch(batch) = msg {
                                stats.bam1_count += batch.len() as u64;
                            }
                        }
                    }
                    if pending_batch2.is_some() {
                        while let Ok(msg) = rx2.recv() {
                            if let BatchMessage::Batch(batch) = msg {
                                stats.bam2_count += batch.len() as u64;
                            }
                        }
                    }
                    break;
                }
                (Some(_), Some(_)) => {}
            }

            // Compare batches in parallel
            let batch1 = pending_batch1.take().unwrap();
            let batch2 = pending_batch2.take().unwrap();

            // Handle unequal batch sizes
            let min_len = batch1.len().min(batch2.len());
            let (cmp_batch1, remainder1) = batch1.split_at(min_len);
            let (cmp_batch2, remainder2) = batch2.split_at(min_len);

            // Compare the aligned portions in parallel
            let (results, core_m, core_d, tag_m, tag_d, tag_ord) =
                compare_batch_parallel(cmp_batch1, cmp_batch2, &header1, &header2, current_index);

            stats.core_matches += core_m as u64;
            stats.core_diffs += core_d as u64;
            stats.tag_matches += tag_m as u64;
            stats.tag_diffs += tag_d as u64;
            stats.tag_order_diffs += tag_ord as u64;

            // Collect diff details (limited by max_diffs)
            for r in results {
                if let Some(detail) = r.diff_detail {
                    if stats.diff_details.len() < self.max_diffs {
                        stats.diff_details.push(detail);
                    }
                }
            }

            current_index += min_len as u64;
            progress.log_if_needed(min_len as u64);

            // Handle remainders - put them back as pending
            if !remainder1.is_empty() {
                pending_batch1 = Some(remainder1.to_vec());
            }
            if !remainder2.is_empty() {
                pending_batch2 = Some(remainder2.to_vec());
            }
        }

        progress.log_final();

        let is_equal =
            stats.bam1_count == stats.bam2_count && stats.core_diffs == 0 && stats.tag_diffs == 0;

        if !self.quiet {
            println!("=== BAM Comparison Results (content mode) ===");
            println!("BAM1: {}", self.bam1.display());
            println!("BAM2: {}", self.bam2.display());
            println!();
            println!("Record counts: {} vs {}", stats.bam1_count, stats.bam2_count);
            println!("Core field matches: {}", stats.core_matches);
            println!("Core field diffs: {}", stats.core_diffs);
            println!("Tag value matches: {}", stats.tag_matches);
            println!("Tag value diffs: {}", stats.tag_diffs);
            println!("Tag order diffs (values match): {}", stats.tag_order_diffs);
            println!();

            if is_equal {
                println!("RESULT: BAM files are IDENTICAL (core fields and tag values match)");
                if stats.tag_order_diffs > 0 {
                    println!(
                        "  Note: {} records have tags in different order",
                        stats.tag_order_diffs
                    );
                }
            } else {
                println!("RESULT: BAM files DIFFER");
                if !stats.diff_details.is_empty() {
                    println!("\nFirst {} differences:", stats.diff_details.len());
                    for detail in &stats.diff_details {
                        println!("  Record {}: {}", detail.record_num, detail.qname);
                        println!("    Flag: {}", detail.flags);
                        println!("    Type: {}", detail.diff_type);
                        for d in &detail.diffs {
                            println!("      {d}");
                        }
                    }
                }
            }
        }

        if is_equal {
            info!("BAM files are identical");
            Ok(stats.bam1_count)
        } else {
            info!("BAM files differ");
            std::process::exit(1);
        }
    }

    /// Execute full comparison mode
    ///
    /// This mode combines grouping equivalence check with full content comparison.
    /// Both files must be in the same order (e.g., query-name sorted).
    /// First verifies MI groupings are equivalent, then compares all other fields.
    /// Uses parallel batch processing with double buffering for performance.
    fn execute_full(&self) -> Result<u64> {
        let mut stats = CompareStats::default();
        let mut grouping_stats = GroupingStats::default();
        let mut grouping_errors: Vec<String> = Vec::new();
        let batch_size = self.batch_size;

        info!("Starting full comparison with {} threads, batch size {}", self.threads, batch_size);

        // Maps: read_key -> MI value for each BAM (using ahash)
        let mut mi_map1: AHashMap<ReadKey, String> = AHashMap::new();
        let mut mi_map2: AHashMap<ReadKey, String> = AHashMap::new();

        // Start double-buffered readers for both BAM files
        let (rx1, header1) = start_batch_reader(self.bam1.clone(), self.threads, batch_size)?;
        let (rx2, header2) = start_batch_reader(self.bam2.clone(), self.threads, batch_size)?;

        // Progress tracking
        let progress = ProgressTracker::new("Processed records").with_interval(1_000_000);

        // Process batches
        let mut bam1_eof = false;
        let mut bam2_eof = false;
        let mut pending_batch1: Option<Vec<RecordBuf>> = None;
        let mut pending_batch2: Option<Vec<RecordBuf>> = None;
        let mut current_index = 0u64;

        loop {
            // Get next batch from BAM1 if needed
            if pending_batch1.is_none() && !bam1_eof {
                match rx1.recv() {
                    Ok(BatchMessage::Batch(batch)) => {
                        stats.bam1_count += batch.len() as u64;
                        pending_batch1 = Some(batch);
                    }
                    Ok(BatchMessage::Eof) => bam1_eof = true,
                    Ok(BatchMessage::Error(e)) => bail!("Error reading BAM1: {e}"),
                    Err(_) => bam1_eof = true,
                }
            }

            // Get next batch from BAM2 if needed
            if pending_batch2.is_none() && !bam2_eof {
                match rx2.recv() {
                    Ok(BatchMessage::Batch(batch)) => {
                        stats.bam2_count += batch.len() as u64;
                        pending_batch2 = Some(batch);
                    }
                    Ok(BatchMessage::Eof) => bam2_eof = true,
                    Ok(BatchMessage::Error(e)) => bail!("Error reading BAM2: {e}"),
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
                            if let BatchMessage::Batch(batch) = msg {
                                stats.bam1_count += batch.len() as u64;
                            }
                        }
                    }
                    if pending_batch2.is_some() {
                        while let Ok(msg) = rx2.recv() {
                            if let BatchMessage::Batch(batch) = msg {
                                stats.bam2_count += batch.len() as u64;
                            }
                        }
                    }
                    break;
                }
                (Some(_), Some(_)) => {}
            }

            let batch1 = pending_batch1.take().unwrap();
            let batch2 = pending_batch2.take().unwrap();

            let min_len = batch1.len().min(batch2.len());
            let (cmp_batch1, remainder1) = batch1.split_at(min_len);
            let (cmp_batch2, remainder2) = batch2.split_at(min_len);

            // Compare batches in parallel and collect MI data
            let (results, core_m, core_d, tag_m, tag_d, tag_ord) =
                compare_batch_parallel(cmp_batch1, cmp_batch2, &header1, &header2, current_index);

            // Process grouping data from the batch (sequential for MI map building)
            for (r1, r2) in cmp_batch1.iter().zip(cmp_batch2.iter()) {
                grouping_stats.total_records += 1;

                let qname1 = record_name_to_string(r1);
                let qname2 = record_name_to_string(r2);
                let is_read1 = r1.flags().is_first_segment();

                if qname1 != qname2 {
                    grouping_stats.order_mismatches += 1;
                    continue;
                }

                let read_key: ReadKey = (qname1, is_read1);
                if let Some(mi1) = get_mi_tag(r1) {
                    mi_map1.insert(read_key.clone(), mi1);
                }
                if let Some(mi2) = get_mi_tag(r2) {
                    mi_map2.insert(read_key, mi2);
                }
            }

            stats.core_matches += core_m as u64;
            stats.core_diffs += core_d as u64;
            stats.tag_matches += tag_m as u64;
            stats.tag_diffs += tag_d as u64;
            stats.tag_order_diffs += tag_ord as u64;

            for r in results {
                if let Some(detail) = r.diff_detail {
                    if stats.diff_details.len() < self.max_diffs {
                        stats.diff_details.push(detail);
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

        // Phase 2: Verify grouping equivalence
        let unique_mi1: AHashSet<&String> = mi_map1.values().collect();
        let unique_mi2: AHashSet<&String> = mi_map2.values().collect();

        // Build reverse maps: MI -> set of read keys
        let mi_to_reads1 = build_mi_groups(&mi_map1);
        let mi_to_reads2 = build_mi_groups(&mi_map2);

        // For each MI group in BAM1, verify all reads have the same MI in BAM2
        for (mi1, reads1) in &mi_to_reads1 {
            let mi2_values: HashSet<Option<&String>> =
                reads1.iter().map(|k| mi_map2.get(*k)).collect();

            if mi2_values.len() > 1 || mi2_values.contains(&None) {
                grouping_stats.grouping_mismatches += 1;
                if grouping_errors.len() < self.max_diffs {
                    let mi2_strs: Vec<String> =
                        mi2_values.iter().map(|m| format!("{m:?}")).collect();
                    grouping_errors.push(format!(
                        "MI group '{mi1}' in BAM1 maps to multiple MIs in BAM2: {mi2_strs:?}"
                    ));
                }
            }
        }

        // Verify the reverse: each MI group in BAM2 maps to single MI in BAM1
        for (mi2, reads2) in &mi_to_reads2 {
            let mi1_values: HashSet<Option<&String>> =
                reads2.iter().map(|k| mi_map1.get(*k)).collect();

            if mi1_values.len() > 1 || mi1_values.contains(&None) {
                grouping_stats.grouping_mismatches += 1;
                if grouping_errors.len() < self.max_diffs {
                    let mi1_strs: Vec<String> =
                        mi1_values.iter().map(|m| format!("{m:?}")).collect();
                    grouping_errors.push(format!(
                        "MI group '{mi2}' in BAM2 maps to multiple MIs in BAM1: {mi1_strs:?}"
                    ));
                }
            }
        }

        let content_equal =
            stats.bam1_count == stats.bam2_count && stats.core_diffs == 0 && stats.tag_diffs == 0;
        let grouping_equal =
            grouping_stats.order_mismatches == 0 && grouping_stats.grouping_mismatches == 0;
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
            println!("Unique MI values: {} vs {}", unique_mi1.len(), unique_mi2.len());
            println!("Order mismatches: {}", grouping_stats.order_mismatches);
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
                if unique_mi1.len() != unique_mi2.len() {
                    println!(
                        "  Note: Different number of unique MI values ({} vs {}), but groupings match",
                        unique_mi1.len(),
                        unique_mi2.len()
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
            std::process::exit(1);
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
    fn execute_grouping(&self) -> Result<u64> {
        if self.ignore_order {
            return self.execute_grouping_unordered();
        }

        let mut stats = GroupingStats::default();
        let mut diff_details: Vec<DiffDetail> = Vec::new();
        let batch_size = self.batch_size;

        info!(
            "Starting grouping comparison with {} threads, batch size {}",
            self.threads, batch_size
        );

        // Maps: read_key_hash -> MI value for each BAM (compact representation)
        // Using u64 hash for keys and i64 for MI values saves ~80% memory
        let mut mi_map1: AHashMap<ReadKeyHash, i64> = AHashMap::new();
        let mut mi_map2: AHashMap<ReadKeyHash, i64> = AHashMap::new();

        // Start double-buffered readers for both BAM files
        let (rx1, _header1) = start_batch_reader(self.bam1.clone(), self.threads, batch_size)?;
        let (rx2, _header2) = start_batch_reader(self.bam2.clone(), self.threads, batch_size)?;

        // Progress tracking
        let progress = ProgressTracker::new("Processed records").with_interval(1_000_000);

        // Process batches
        let mut bam1_eof = false;
        let mut bam2_eof = false;
        let mut pending_batch1: Option<Vec<RecordBuf>> = None;
        let mut pending_batch2: Option<Vec<RecordBuf>> = None;
        let mut current_index = 0u64;

        loop {
            // Get next batch from BAM1 if needed
            if pending_batch1.is_none() && !bam1_eof {
                match rx1.recv() {
                    Ok(BatchMessage::Batch(batch)) => {
                        pending_batch1 = Some(batch);
                    }
                    Ok(BatchMessage::Eof) => bam1_eof = true,
                    Ok(BatchMessage::Error(e)) => bail!("Error reading BAM1: {e}"),
                    Err(_) => bam1_eof = true,
                }
            }

            // Get next batch from BAM2 if needed
            if pending_batch2.is_none() && !bam2_eof {
                match rx2.recv() {
                    Ok(BatchMessage::Batch(batch)) => {
                        pending_batch2 = Some(batch);
                    }
                    Ok(BatchMessage::Eof) => bam2_eof = true,
                    Ok(BatchMessage::Error(e)) => bail!("Error reading BAM2: {e}"),
                    Err(_) => bam2_eof = true,
                }
            }

            // Check for completion
            match (&pending_batch1, &pending_batch2) {
                (None, None) => break,
                (Some(_), None) | (None, Some(_)) => {
                    if diff_details.len() < self.max_diffs {
                        diff_details.push(DiffDetail {
                            record_num: current_index,
                            qname: "N/A".to_string(),
                            flags: "N/A".to_string(),
                            diff_type: DiffType::CountMismatch,
                            diffs: vec!["BAM files have different number of records".to_string()],
                        });
                    }
                    bail!("BAM files have different record counts");
                }
                (Some(_), Some(_)) => {}
            }

            let batch1 = pending_batch1.take().unwrap();
            let batch2 = pending_batch2.take().unwrap();

            let min_len = batch1.len().min(batch2.len());
            let (cmp_batch1, remainder1) = batch1.split_at(min_len);
            let (cmp_batch2, remainder2) = batch2.split_at(min_len);

            // Compare batches in parallel for grouping data
            let results = compare_batch_grouping_parallel(
                cmp_batch1,
                cmp_batch2,
                current_index,
                self.max_diffs,
                diff_details.len(),
            );

            // Process results and build MI maps
            for r in results {
                stats.total_records += 1;

                if !r.name_match {
                    stats.order_mismatches += 1;
                    if let Some(detail) = r.diff_detail {
                        if diff_details.len() < self.max_diffs {
                            diff_details.push(detail);
                        }
                    }
                    continue;
                }

                if !r.flag_match {
                    stats.order_mismatches += 1;
                    if let Some(detail) = r.diff_detail {
                        if diff_details.len() < self.max_diffs {
                            diff_details.push(detail);
                        }
                    }
                    continue;
                }

                // Use compact representation: hash the read key and parse MI as i64
                let key_hash = hash_read_key(&r.read_key.0, r.read_key.1);

                match (&r.mi1, &r.mi2) {
                    (Some(m1), Some(m2)) => {
                        // Parse MI strings as i64 for compact storage
                        if let (Ok(mi1_val), Ok(mi2_val)) = (m1.parse::<i64>(), m2.parse::<i64>()) {
                            mi_map1.insert(key_hash, mi1_val);
                            mi_map2.insert(key_hash, mi2_val);
                        }
                    }
                    (None, Some(_)) => {
                        stats.missing_mi_bam1 += 1;
                        if diff_details.len() < self.max_diffs {
                            diff_details.push(DiffDetail {
                                record_num: r.record_num,
                                qname: r.read_key.0,
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
                                qname: r.read_key.0,
                                flags: "N/A".to_string(),
                                diff_type: DiffType::TagDiff,
                                diffs: vec!["MI tag missing in BAM2".to_string()],
                            });
                        }
                    }
                    (None, None) => {}
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
        for (mi1, read_hashes) in &bam1_groups {
            let mi2_values: AHashSet<i64> =
                read_hashes.iter().filter_map(|key_hash| mi_map2.get(key_hash).copied()).collect();

            if mi2_values.len() > 1 {
                stats.grouping_mismatches += 1;
                if grouping_errors.len() < self.max_diffs {
                    grouping_errors.push(format!(
                        "MI group '{}' in BAM1 ({} reads) maps to {} different MIs in BAM2: {:?}",
                        mi1,
                        read_hashes.len(),
                        mi2_values.len(),
                        mi2_values.iter().take(5).collect::<Vec<_>>()
                    ));
                }
            }
        }

        // Check BAM2 groups -> BAM1
        for (mi2, read_hashes) in &bam2_groups {
            let mi1_values: AHashSet<i64> =
                read_hashes.iter().filter_map(|key_hash| mi_map1.get(key_hash).copied()).collect();

            if mi1_values.len() > 1 {
                stats.grouping_mismatches += 1;
                if grouping_errors.len() < self.max_diffs {
                    grouping_errors.push(format!(
                        "MI group '{}' in BAM2 ({} reads) maps to {} different MIs in BAM1: {:?}",
                        mi2,
                        read_hashes.len(),
                        mi1_values.len(),
                        mi1_values.iter().take(5).collect::<Vec<_>>()
                    ));
                }
            }
        }

        // Determine result
        let is_equivalent = stats.order_mismatches == 0
            && stats.missing_mi_bam1 == 0
            && stats.missing_mi_bam2 == 0
            && stats.grouping_mismatches == 0;

        if !self.quiet {
            println!("=== BAM Comparison Results (grouping mode) ===");
            println!("BAM1: {}", self.bam1.display());
            println!("BAM2: {}", self.bam2.display());
            println!();
            println!("Total records compared: {}", stats.total_records);
            println!("Order/name mismatches: {}", stats.order_mismatches);
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
            std::process::exit(1);
        }
    }

    /// Execute grouping comparison in order-independent mode
    ///
    /// Uses parallel batch processing to build MI maps for both BAM files,
    /// then compares them for set membership and MI equivalence.
    /// This approach enables full use of multi-threaded BGZF decompression
    /// and rayon parallel processing.
    fn execute_grouping_unordered(&self) -> Result<u64> {
        let mut stats = UnorderedGroupingStats::default();
        let mut grouping_errors: Vec<String> = Vec::new();

        info!(
            "Starting order-independent grouping comparison with {} threads, batch size {}",
            self.threads, self.batch_size
        );

        // Create a thread pool with the specified number of threads
        // This controls BOTH rayon parallelism and BGZF decompression
        let pool = rayon::ThreadPoolBuilder::new()
            .num_threads(self.threads)
            .build()
            .map_err(|e| anyhow!("Failed to create thread pool: {e}"))?;

        // Run all parallel work inside the controlled thread pool
        let result = pool.install(|| -> Result<()> {
            // Phase 1: Build MI map for BAM1
            info!("Phase 1: Building MI map for BAM1...");
            let (mi_map1, total_bam1) =
                build_mi_map_parallel(&self.bam1, self.threads, self.batch_size)?;
            stats.total_bam1 = total_bam1;
            info!("BAM1: {total_bam1} records");

            // Phase 2: Build MI map for BAM2
            info!("Phase 2: Building MI map for BAM2...");
            let (mi_map2, total_bam2) =
                build_mi_map_parallel(&self.bam2, self.threads, self.batch_size)?;
            stats.total_bam2 = total_bam2;
            info!("BAM2: {total_bam2} records");

            // Phase 3: Compare set membership in parallel
            info!("Phase 3: Comparing set membership (parallel)...");
            let ((only_in_bam1, matched), only_in_bam2) = rayon::join(
                || {
                    // Count keys only in BAM1 and matched keys in one pass
                    mi_map1
                        .par_iter()
                        .fold(
                            || (0u64, 0u64),
                            |(only, matched), (k, _)| {
                                if mi_map2.contains_key(k) {
                                    (only, matched + 1)
                                } else {
                                    (only + 1, matched)
                                }
                            },
                        )
                        .reduce(|| (0, 0), |(a1, a2), (b1, b2)| (a1 + b1, a2 + b2))
                },
                || {
                    // Count keys only in BAM2
                    mi_map2.par_iter().filter(|(k, _)| !mi_map1.contains_key(k)).count() as u64
                },
            );
            stats.only_in_bam1 = only_in_bam1;
            stats.matched = matched;
            stats.only_in_bam2 = only_in_bam2;

            // Phase 4: Verify MI grouping equivalence (memory-efficient)
            // Instead of building full group maps, we verify consistency in a single pass
            // For each MI in BAM1, track the first MI seen in BAM2 - any deviation is a mismatch
            info!("Phase 4: Verifying grouping equivalence...");
            let max_diffs = self.max_diffs;

            // Check BAM1 groups -> BAM2: for each mi1, all reads should map to same mi2
            // Use a map: mi1 -> (first_mi2_seen, count, has_mismatch)
            let mut mi1_to_mi2: AHashMap<i64, (i64, u64, bool)> = AHashMap::new();
            for (key_hash, mi1) in &mi_map1 {
                if let Some(&mi2) = mi_map2.get(key_hash) {
                    mi1_to_mi2
                        .entry(*mi1)
                        .and_modify(|(first_mi2, count, has_mismatch)| {
                            *count += 1;
                            if *first_mi2 != mi2 {
                                *has_mismatch = true;
                            }
                        })
                        .or_insert((mi2, 1, false));
                }
            }

            stats.unique_groups_bam1 = mi1_to_mi2.len();
            let mismatches1: Vec<_> = mi1_to_mi2
                .iter()
                .filter(|(_, (_, _, has_mismatch))| *has_mismatch)
                .map(|(mi1, (_, count, _))| {
                    format!("MI group '{mi1}' in BAM1 ({count} reads) maps to multiple MIs in BAM2")
                })
                .collect();

            // Check BAM2 groups -> BAM1: for each mi2, all reads should map to same mi1
            let mut mi2_to_mi1: AHashMap<i64, (i64, u64, bool)> = AHashMap::new();
            for (key_hash, mi2) in &mi_map2 {
                if let Some(&mi1) = mi_map1.get(key_hash) {
                    mi2_to_mi1
                        .entry(*mi2)
                        .and_modify(|(first_mi1, count, has_mismatch)| {
                            *count += 1;
                            if *first_mi1 != mi1 {
                                *has_mismatch = true;
                            }
                        })
                        .or_insert((mi1, 1, false));
                }
            }

            stats.unique_groups_bam2 = mi2_to_mi1.len();
            let mismatches2: Vec<_> = mi2_to_mi1
                .iter()
                .filter(|(_, (_, _, has_mismatch))| *has_mismatch)
                .map(|(mi2, (_, count, _))| {
                    format!("MI group '{mi2}' in BAM2 ({count} reads) maps to multiple MIs in BAM1")
                })
                .collect();

            stats.grouping_mismatches = (mismatches1.len() + mismatches2.len()) as u64;
            grouping_errors.extend(mismatches1.into_iter().take(max_diffs));
            grouping_errors.extend(
                mismatches2.into_iter().take(max_diffs.saturating_sub(grouping_errors.len())),
            );

            Ok(())
        });

        result?;

        // Determine result
        let is_equivalent =
            stats.only_in_bam1 == 0 && stats.only_in_bam2 == 0 && stats.grouping_mismatches == 0;

        if !self.quiet {
            println!("=== BAM Comparison Results (grouping mode, order-independent) ===");
            println!("BAM1: {}", self.bam1.display());
            println!("BAM2: {}", self.bam2.display());
            println!();
            println!("Total records in BAM1: {}", stats.total_bam1);
            println!("Total records in BAM2: {}", stats.total_bam2);
            println!("Records matched: {}", stats.matched);
            println!("Records only in BAM1: {}", stats.only_in_bam1);
            println!("Records only in BAM2: {}", stats.only_in_bam2);
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

                if !grouping_errors.is_empty() {
                    println!("\nDifferences (first {}):", grouping_errors.len());
                    for err in &grouping_errors {
                        println!("  {err}");
                    }
                }
            }
        }

        if is_equivalent {
            info!("BAM groupings are equivalent (order-independent)");
            Ok(stats.total_bam1 + stats.total_bam2)
        } else {
            info!("BAM groupings differ");
            std::process::exit(1);
        }
    }
}
