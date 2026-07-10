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
use anyhow::Result;
use clap::{Parser, ValueEnum};
use crossbeam_channel::{Receiver, bounded};
use fgumi_bam_io::{RawBamReaderAuto, create_raw_bam_reader};
use fgumi_raw_bam::fields as raw_fields;
use fgumi_raw_bam::{RawRecord, find_int_tag, find_string_tag};
use log::info;
use noodles::sam::Header;
use std::path::PathBuf;
use std::thread;

use crate::commands::command::Command;
use crate::commands::common::{MemoryReserve, parse_bool, resolve_memory_budget};
use crate::commands::sort::TMP_DIRS_ENV;

use super::engines::content::ContentPredicate;
use super::engines::keyjoin::{self, KeyJoinConfig};
use super::engines::positional::positional_compare;

/// Comparison mode for BAM files
#[derive(Debug, Clone, Copy, Default, ValueEnum)]
pub enum CompareMode {
    /// Content comparison (default): all fields and tags must match exactly.
    ///
    /// A pure record-by-record comparison. Records are paired by their
    /// [`RecordKey`](super::record_key::RecordKey) identity (see
    /// [`super::engines::positional`]), so a difference in record order is
    /// reported honestly as a difference rather than silently masked — this is
    /// the *sound* default.
    #[default]
    Content,
    /// Grouping comparison: verify MI groupings are equivalent (for grouped BAMs).
    ///
    /// Canonicalizes both inputs to queryname order and merge-joins by
    /// [`RecordKey`](super::record_key::RecordKey), then verifies the
    /// fgumi-MI/fgbio-MI mapping is a consistent
    /// bijection (order-independent; see [`super::engines::keyjoin`]). Also
    /// checks non-MI content on each matched pair under
    /// [`ContentPredicate::ExactMinusMi`](super::engines::content::ContentPredicate::ExactMinusMi).
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
    /// (see [`ContentPredicate::Exact`](super::engines::content::ContentPredicate::Exact)).
    Simplex,
    /// Duplex consensus output: positional, exact content comparison
    /// (see [`ContentPredicate::Exact`](super::engines::content::ContentPredicate::Exact)).
    Duplex,
    /// CODEC consensus output: positional, exact content comparison
    /// (see [`ContentPredicate::Exact`](super::engines::content::ContentPredicate::Exact)).
    Codec,
    /// Filter output: consensus reads with reads dropped, but the same depth tags
    /// (`cD`/`cM`/`cE`, and duplex `aD`/`aM`/`bD`/`bM`/`aE`/`bE`) as the consensus command
    /// that produced them, passed through unchanged. Compares via the same predicate as
    /// consensus output
    /// ([`ContentPredicate::Exact`](super::engines::content::ContentPredicate::Exact));
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
    /// - `Filter | Simplex | Duplex | Codec` → `(Content, false, Exact)`: all four
    ///   deal in consensus reads carrying the depth tags (`cD`/`cM`/`cE`, and duplex
    ///   `aD`/`aM`/`bD`/`bM`/`aE`/`bE`). `Simplex`/`Duplex`/`Codec` are the commands that
    ///   write those tags; `filter` only drops reads afterward, it never rewrites them. All
    ///   four compare positionally under
    ///   [`ContentPredicate::Exact`](super::engines::content::ContentPredicate::Exact):
    ///   because fgumi clamps the depth tags to fgbio's `Short` ceiling at the source, the
    ///   tags are bit-identical and an exact comparison is both sound and complete — no
    ///   consensus-specific predicate is needed.
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
                (CompareMode::Content, false, ContentPredicate::Exact)
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

  content (default):
    Pure record-by-record comparison of all fields and tags. Records are
    paired by their RecordKey identity, so a difference in record order is
    reported honestly as a difference rather than silently masked (the sound
    default). Does not analyze MI groupings - just compares raw content.

  grouping:
    For comparing grouped BAM files where MI assignment order may differ.
    Order-independent: each input is internally canonicalized to queryname
    order and merge-joined by record identity (so --ignore-order is a no-op
    here, retained only for backward compatibility). Every matched pair is
    compared under EXACT-MI (every field except the MI tag), so a non-MI
    content difference DIFFERs even when the MI grouping itself is untouched,
    in addition to verifying the fgumi-MI/fgbio-MI mapping is a consistent
    bijection (including the base-level duplex strand-pairing check).

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
  simplex         content     false            Exact (cD/cM/cE compared exactly; clamped at source)
  duplex          content     false            Exact (cD/cM/cE + aD/aM/bD/bM/aE/bE, exact)
  codec           content     false            Exact (same as simplex/duplex)

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

    # Preset + explicit override (e.g. same-tool group comparison via content):
    fgumi compare bams --command group --mode content a.bam b.bam

Example usage:
  fgumi compare bams bam1.bam bam2.bam                    # content mode (default)
  fgumi compare bams bam1.bam bam2.bam --mode content    # content (explicit)
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

    /// Comparison mode: 'content' (all fields and tags, order-sound positional
    /// pairing — for extract/zipper/correct/dedup/simplex/duplex/codec/filter
    /// output), or 'grouping' (MI equivalence via the order-independent key-join
    /// engine — for group output). Overrides `--command` preset if both given.
    /// Defaults to 'content' when neither `--mode` nor `--command` is set.
    #[arg(long = "mode")]
    pub mode: Option<CompareMode>,

    /// Maximum number of differences to report in detail
    #[arg(short = 'm', long = "max-diffs", default_value = "10")]
    pub max_diffs: usize,

    /// Quiet mode - only exit code indicates result (0=equal, 1=different)
    #[arg(short = 'q', long = "quiet", default_value = "false", num_args = 0..=1, default_missing_value = "true", action = clap::ArgAction::Set, value_parser = parse_bool)]
    pub quiet: bool,

    /// Ignore record order when comparing. Only valid with `--mode grouping`
    /// (rejected with any other mode). Note that grouping mode is *already*
    /// order-independent — it canonicalizes both inputs to queryname order via
    /// the key-join engine — so this flag is effectively a no-op there, retained
    /// for backward compatibility.
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
            let timer = (!self.quiet).then(|| OperationTimer::new("Comparing BAMs"));
            let total_records = self.execute_sort_verify()?;
            if let Some(timer) = timer {
                timer.log_completion(total_records);
            }
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

        // In `--quiet` mode only the exit code communicates the result, so suppress this
        // command's own informational stderr logging and timer (matching `compare metrics`).
        if let Some(preset) = self.command.filter(|_| !self.quiet) {
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

        let timer = (!self.quiet).then(|| OperationTimer::new("Comparing BAMs"));

        // This predicate is only ever consumed below by the `Content` mode branch
        // (`self.execute_content(predicate)`); `CompareMode::Grouping` routes to the
        // order-independent key-join engine (`execute_grouping`), which takes no predicate at
        // all — it hardcodes `ContentPredicate::ExactMinusMi` inside
        // `engines::keyjoin::keyjoin_compare` itself, not configurable via `--command`/`--mode`
        // (so the `Grouping` arm's value is never actually used).
        //
        // An *explicit* `--mode content` (`self.mode == Some(Content)`) always means exact
        // content comparison — it overrides a preset's own predicate too, so
        // `--command group --mode content` compares MI exactly (a preset supplies its
        // relaxed predicate only when the user did *not* pin the mode themselves). A preset's
        // `content_predicate` (see `CommandPreset::resolve`, e.g. why
        // `Simplex`/`Duplex`/`Codec`/`Filter` use plain `Exact` while `group` uses
        // `ExactMinusMi`) is therefore consulted only in the preset-with-implicit-mode case.
        let predicate = match (mode, self.mode, self.command) {
            (CompareMode::Content, Some(_), _) => ContentPredicate::Exact,
            (CompareMode::Content, None, Some(preset)) => preset.content_predicate(),
            (CompareMode::Content, None, None) => ContentPredicate::Exact,
            (CompareMode::Grouping, _, _) => ContentPredicate::ExactMinusMi,
        };

        let total_records = match mode {
            CompareMode::Content => self.execute_content(predicate)?,
            CompareMode::Grouping => self.execute_grouping()?,
        };

        if let Some(timer) = timer {
            timer.log_completion(total_records);
        }
        Ok(())
    }
}

impl CompareBams {
    /// Resolve the effective `(mode, ignore_order)` pair.
    ///
    /// Precedence: explicit flag > `--command` preset default > built-in default
    /// (`content`, `false`). `ignore_order` is silently coerced to `false` when the
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

    /// Execute content comparison mode.
    ///
    /// Delegates pairing and equality entirely to the positional engine
    /// ([`positional_compare`]): records are paired purely by index, a
    /// [`RecordKey`](super::record_key::RecordKey) mismatch stops pairing immediately
    /// (never resyncing), and remaining pairs are compared under `predicate` (`Exact`
    /// for the default and the consensus/`filter` presets, or `ExactMinusMi` when a
    /// grouping preset supplies it — see `execute()`).
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
        if !self.quiet {
            info!(
                "Starting content comparison with {} threads, batch size {}",
                self.threads, self.batch_size
            );
        }

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
            if !self.quiet {
                info!("BAM files are identical");
            }
            Ok(outcome.bam1_count)
        } else {
            if !self.quiet {
                info!("BAM files differ");
            }
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
        if !self.quiet {
            info!("Using --command sort preset: sort-key-run verification (no positional pairing)");
        }

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
            if !self.quiet {
                info!("BAM files are identical");
            }
            Ok(outcome.bam1_count)
        } else {
            if !self.quiet {
                info!("BAM files differ");
            }
            Err(super::CompareMismatch("BAM files differ".to_owned()).into())
        }
    }

    /// Execute grouping comparison via the order-independent key-join engine.
    ///
    /// This is the sole `CompareMode::Grouping` path (and what `--command group` uses).
    /// Both inputs are canonicalized to queryname order and merge-joined by
    /// [`RecordKey`](super::record_key::RecordKey) (see [`keyjoin::keyjoin_compare`]), so it
    /// is inherently order-independent — `--ignore-order` is therefore a no-op under
    /// `--mode grouping`. Besides verifying the fgumi-MI/fgbio-MI mapping is a consistent
    /// bijection (including the `base()`-level duplex strand-pairing check), it also catches
    /// a non-MI content difference on any matched pair. The content check always uses
    /// [`ContentPredicate::ExactMinusMi`], hardcoded inside `keyjoin::keyjoin_compare` itself
    /// — this engine takes no predicate argument and is not configurable via
    /// `--command`/`--mode`.
    ///
    /// # Errors
    ///
    /// Returns an error if either input cannot be canonicalized or read (see
    /// [`keyjoin::keyjoin_compare`]), or [`super::CompareMismatch`] if the two BAMs are
    /// found to differ (non-zero exit via the `Command` trait).
    fn execute_grouping(&self) -> Result<u64> {
        if !self.quiet {
            info!(
                "Starting key-join grouping comparison with {} threads, sort memory {:?}",
                self.threads, self.sort_memory
            );
        }

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
            if !self.quiet {
                info!("BAM groupings are equivalent (key-join)");
            }
            Ok(outcome.bam1_count + outcome.bam2_count)
        } else {
            if !self.quiet {
                info!("BAM groupings differ (key-join)");
            }
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
    /// comparison under the `Exact` predicate (see `execute()`'s predicate selection)
    /// rather than MI-grouping equivalence — this is a breaking strictness change (see
    /// the commit message and the compare-hardening design spec's §"Accepted divergences").
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
    /// the same `Exact` predicate as `Simplex`/`Duplex`/`Codec`, keeping the four consensus
    /// commands on one predicate.
    #[test]
    fn filter_preset_content_predicate_is_exact() {
        assert_eq!(CommandPreset::Filter.content_predicate(), ContentPredicate::Exact);
    }

    #[test]
    fn effective_settings_fall_back_to_built_in_defaults() {
        let args = parse(&[]);
        let (mode, ignore) = args.effective_settings();
        assert!(matches!(mode, CompareMode::Content));
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
        // --mode grouping overrides --command extract's Content default.
        let args = parse(&["--command", "extract", "--mode", "grouping"]);
        let (mode, ignore) = args.effective_settings();
        assert!(matches!(mode, CompareMode::Grouping));
        // extract presets ignore_order=false and none was passed explicitly.
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
