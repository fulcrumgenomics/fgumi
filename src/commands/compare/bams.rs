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
use fgumi_lib::bam_io::{RawBamReaderAuto, create_raw_bam_reader};
use fgumi_lib::logging::OperationTimer;
use fgumi_lib::progress::ProgressTracker;
use fgumi_lib::validation::validate_file_exists;
use fgumi_raw_bam::fields as raw_fields;
use fgumi_raw_bam::{RawRecord, find_int_tag, find_string_tag};
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
use std::collections::BTreeMap;
use std::path::{Path, PathBuf};
use std::thread;

use crate::commands::command::Command;
use crate::commands::common::parse_bool;

use super::raw_compare::{raw_compare_structured, raw_records_byte_equal};

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
    #[arg(short = 'q', long = "quiet", default_value = "false", num_args = 0..=1, default_missing_value = "true", action = clap::ArgAction::Set, value_parser = parse_bool)]
    pub quiet: bool,

    /// Ignore record order when comparing in grouping mode.
    /// Required for comparing output from consensus commands (simplex/duplex/codec)
    /// when run with --threads, as parallel processing causes non-deterministic ordering.
    /// Only valid with --mode grouping.
    #[arg(long = "ignore-order", default_value = "false", num_args = 0..=1, default_missing_value = "true", action = clap::ArgAction::Set, value_parser = parse_bool)]
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

/// Convert a record's name to a String, returning "*" if missing.
fn record_name_to_string(record: &RecordBuf) -> String {
    record.name().map_or_else(|| "*".to_string(), std::string::ToString::to_string)
}

/// Check if the first-in-template (R1) flag is set in raw BAM record bytes.
fn is_first_segment_raw(raw: &RawRecord) -> bool {
    let flags = raw_fields::flags(raw.as_ref());
    flags & raw_fields::flags::FIRST_SEGMENT != 0
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
    key_hash: ReadKeyHash,
    /// Read name as String, only populated when needed for error reporting
    read_name_for_display: Option<String>,
    mi1: Option<i64>,
    mi2: Option<i64>,
    name_match: bool,
    flag_match: bool,
    diff_detail: Option<DiffDetail>,
}

// ============================================================================
// Types and MI map helpers (using ahash)
// ============================================================================

/// Compact read key using hash - saves ~70 bytes per entry vs (String, bool)
type ReadKeyHash = u64;

/// Compute a hash for a read key from raw bytes (avoids String allocation).
#[inline]
fn hash_read_key_raw(name: &[u8], is_read1: bool) -> ReadKeyHash {
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

/// Result from parallel MI extraction for a single record
struct MiExtractResult {
    key_hash: ReadKeyHash,
    mi: Option<i64>,
}

/// Extract MI tag from raw BAM record bytes directly as i64.
/// Tries integer-type MI first, then falls back to string-type MI parsed as i64.
fn get_mi_tag_raw_i64(raw: &RawRecord) -> Option<i64> {
    let aux = raw_fields::aux_data_slice(raw.as_ref());
    if let Some(v) = find_int_tag(aux, b"MI") {
        return Some(v);
    }
    find_string_tag(aux, b"MI").and_then(|bytes| std::str::from_utf8(bytes).ok()?.parse().ok())
}

/// Build an MI map from a BAM file using parallel batch processing.
/// Returns a map from read key hash to MI value, along with total record count.
fn build_mi_map_parallel(
    path: &Path,
    threads: usize,
    batch_size: usize,
) -> Result<(AHashMap<ReadKeyHash, i64>, u64)> {
    let (rx, _header) = start_raw_batch_reader(path.to_path_buf(), threads, batch_size)?;

    let mut mi_map: AHashMap<ReadKeyHash, i64> = AHashMap::new();
    let mut total_records: u64 = 0;

    loop {
        match rx.recv() {
            Ok(RawBatchMessage::Batch(batch)) => {
                // Extract MI values in parallel using rayon
                let results: Vec<MiExtractResult> = batch
                    .par_iter()
                    .map(|raw| {
                        let name_bytes = raw_fields::read_name(raw.as_ref());
                        let is_read1 = is_first_segment_raw(raw);
                        let key_hash = hash_read_key_raw(name_bytes, is_read1);

                        let mi = get_mi_tag_raw_i64(raw);

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
            Ok(RawBatchMessage::Eof) => break,
            Ok(RawBatchMessage::Error(e)) => bail!("Error reading BAM: {e}"),
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

/// Message type for the double-buffered raw reader channel.
enum RawBatchMessage {
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
fn start_raw_batch_reader(
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
/// This is used only when records differ at the raw byte level and we need
/// human-readable field values for diff reporting. The raw bytes are the BAM
/// record body WITHOUT the 4-byte length prefix.
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
/// **Tier 3:** Full deserialization — only when records actually differ, deserialize
///   to `RecordBuf` for human-readable diff reporting.
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

            // Tier 3: Full deserialization for diff reporting
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
                            diff_type: DiffType::CoreDiff,
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
                            diff_type: DiffType::CoreDiff,
                            diffs: vec![format!("Failed to deserialize record from BAM2: {e}")],
                        }),
                    };
                }
            };

            let diff_detail = if raw_result.core_match {
                // Core matches but tags differ
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
                // Core fields differ
                let core1 = get_core_fields(&rec1, header1);
                let core2 = get_core_fields(&rec2, header2);
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
fn compare_raw_batch_grouping_parallel(
    batch1: &[RawRecord],
    batch2: &[RawRecord],
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
            let name1_bytes = raw_fields::read_name(r1.as_ref());
            let name2_bytes = raw_fields::read_name(r2.as_ref());
            let is_read1_r1 = is_first_segment_raw(r1);
            let is_read1_r2 = is_first_segment_raw(r2);

            let name_match = name1_bytes == name2_bytes;
            let flag_match = is_read1_r1 == is_read1_r2;

            let key_hash = hash_read_key_raw(name1_bytes, is_read1_r1);
            let mi1 = get_mi_tag_raw_i64(r1);
            let mi2 = get_mi_tag_raw_i64(r2);

            let within_diff_budget = existing_diffs + i < max_diffs;
            let needs_detail = within_diff_budget && (!name_match || !flag_match);
            let diff_detail = if needs_detail {
                // Only allocate Strings when we actually need them for error reporting
                let qname1 = String::from_utf8_lossy(name1_bytes).into_owned();
                if name_match {
                    Some(DiffDetail {
                        record_num,
                        qname: qname1,
                        flags: raw_fields::flags(r1.as_ref()).to_string(),
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
                        flags: raw_fields::flags(r1.as_ref()).to_string(),
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

            // Populate read name for display when we need it for missing-MI error reporting
            let needs_name = within_diff_budget && (mi1.is_none() ^ mi2.is_none());
            let read_name_for_display = if needs_name {
                Some(String::from_utf8_lossy(name1_bytes).into_owned())
            } else {
                None
            };

            GroupingCompareResult {
                record_num,
                key_hash,
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

            // Compare batches in parallel
            let batch1 = pending_batch1.take().expect("guarded by (Some, Some) match above");
            let batch2 = pending_batch2.take().expect("guarded by (Some, Some) match above");

            // Handle unequal batch sizes
            let min_len = batch1.len().min(batch2.len());
            let (cmp_batch1, remainder1) = batch1.split_at(min_len);
            let (cmp_batch2, remainder2) = batch2.split_at(min_len);

            // Compare the aligned portions in parallel
            let (results, core_m, core_d, tag_m, tag_d, tag_ord) = compare_raw_batch_parallel(
                cmp_batch1,
                cmp_batch2,
                &header1,
                &header2,
                current_index,
            );

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

        // Maps: read_key_hash -> MI value for each BAM (compact i64 representation)
        let mut mi_map1: AHashMap<ReadKeyHash, i64> = AHashMap::new();
        let mut mi_map2: AHashMap<ReadKeyHash, i64> = AHashMap::new();

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

                let name1_bytes = raw_fields::read_name(r1.as_ref());
                let name2_bytes = raw_fields::read_name(r2.as_ref());
                let is_read1 = is_first_segment_raw(r1);

                if name1_bytes != name2_bytes {
                    grouping_stats.order_mismatches += 1;
                    continue;
                }

                let key_hash = hash_read_key_raw(name1_bytes, is_read1);
                if let Some(mi1) = get_mi_tag_raw_i64(r1) {
                    mi_map1.insert(key_hash, mi1);
                }
                if let Some(mi2) = get_mi_tag_raw_i64(r2) {
                    mi_map2.insert(key_hash, mi2);
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

        // Phase 2: Verify grouping equivalence (using compact representation)
        let mi_to_reads1 = build_mi_groups_compact(&mi_map1);
        let mi_to_reads2 = build_mi_groups_compact(&mi_map2);

        let unique_mi1_count = mi_to_reads1.len();
        let unique_mi2_count = mi_to_reads2.len();

        // For each MI group in BAM1, verify all reads have the same MI in BAM2
        for (mi1, read_hashes) in &mi_to_reads1 {
            let mi2_values: AHashSet<i64> =
                read_hashes.iter().filter_map(|k| mi_map2.get(k).copied()).collect();

            if mi2_values.len() > 1 {
                grouping_stats.grouping_mismatches += 1;
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

        // Verify the reverse: each MI group in BAM2 maps to single MI in BAM1
        for (mi2, read_hashes) in &mi_to_reads2 {
            let mi1_values: AHashSet<i64> =
                read_hashes.iter().filter_map(|k| mi_map1.get(k).copied()).collect();

            if mi1_values.len() > 1 {
                grouping_stats.grouping_mismatches += 1;
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
            println!("Unique MI values: {} vs {}", unique_mi1_count, unique_mi2_count);
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

            let batch1 = pending_batch1.take().expect("guarded by (Some, Some) match above");
            let batch2 = pending_batch2.take().expect("guarded by (Some, Some) match above");

            let min_len = batch1.len().min(batch2.len());
            let (cmp_batch1, remainder1) = batch1.split_at(min_len);
            let (cmp_batch2, remainder2) = batch2.split_at(min_len);

            // Compare batches in parallel for grouping data
            let results = compare_raw_batch_grouping_parallel(
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

                match (r.mi1, r.mi2) {
                    (Some(mi1_val), Some(mi2_val)) => {
                        mi_map1.insert(r.key_hash, mi1_val);
                        mi_map2.insert(r.key_hash, mi2_val);
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
