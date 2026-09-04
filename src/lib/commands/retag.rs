//! Rewrite SAM aux tags: copy, move, or delete a tag by name.
//!
//! `fgumi` is intentionally opinionated about tag names — it consumes and produces
//! a fixed set (`RX`, `MI`, ...) rather than exposing per-tool options to override
//! them. `retag` is the sanctioned escape hatch: a pure per-record tag rewriter so
//! `fgumi`'s core stays opinionated while users convert to/from vendor- or
//! lab-specific tags at the boundary.
//!
//! # Grammar
//!
//! Operations are positional and applied left-to-right, per record:
//!
//! - `SRC::copy::DST` — copy `SRC`'s aux value **verbatim** (type byte + value
//!   bytes) into `DST`; `SRC` stays; an existing `DST` is overwritten.
//! - `SRC::move::DST` — sugar for `SRC::copy::DST` then `SRC::delete`.
//! - `SRC::delete` — drop `SRC`.
//!
//! Because operations apply in order they compose, e.g. `RX::copy::BX RX::copy::CB
//! RX::delete` fans one tag out to two and drops the original. A record missing
//! `SRC` simply skips that operation.
//!
//! This tool is strictly tag↔tag on the records: it never touches read names or
//! sequence. It does synthesize an `@HD` line when the input lacks one and always
//! chains an `@PG` provenance record onto the header. Read-name UMI handling stays
//! in `fgumi extract`.

use std::fmt;
use std::path::PathBuf;
use std::str::FromStr;
use std::sync::Arc;
use std::sync::atomic::AtomicU64;

use anyhow::{Result, bail, ensure};
use clap::Parser;
use fgumi_raw_bam::{RawRecord, append_raw_tag, remove_tag};
use serde::{Deserialize, Serialize};

use crate::commands::command::Command;
use crate::commands::common::{
    BamIoOptions, CompressionOptions, QueueMemoryOptions, SchedulerOptions, ThreadingOptions,
    reject_output_collisions,
};
use crate::per_thread_accumulator::PerThreadAccumulator;
use crate::sam::SamTag;

/// A single tag-rewrite operation, parsed from one `SRC::op::DST` positional argument.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum RetagOp {
    /// Copy `src`'s aux value verbatim into `dst`, leaving `src` in place.
    Copy { src: SamTag, dst: SamTag },
    /// Copy `src` into `dst`, then delete `src`.
    Move { src: SamTag, dst: SamTag },
    /// Delete `src`.
    Delete { src: SamTag },
}

impl RetagOp {
    /// The operation keyword (`"copy"`, `"move"`, or `"delete"`), for metrics and logging.
    #[must_use]
    pub fn kind(&self) -> &'static str {
        match self {
            RetagOp::Copy { .. } => "copy",
            RetagOp::Move { .. } => "move",
            RetagOp::Delete { .. } => "delete",
        }
    }

    /// The source tag the operation reads (and, for `delete`, removes).
    #[must_use]
    pub fn src(&self) -> SamTag {
        match self {
            RetagOp::Copy { src, .. } | RetagOp::Move { src, .. } | RetagOp::Delete { src } => *src,
        }
    }
}

impl fmt::Display for RetagOp {
    /// Render the operation back to its canonical `SRC::op::DST` / `SRC::delete` form.
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            RetagOp::Copy { src, dst } => write!(f, "{src}::copy::{dst}"),
            RetagOp::Move { src, dst } => write!(f, "{src}::move::{dst}"),
            RetagOp::Delete { src } => write!(f, "{src}::delete"),
        }
    }
}

/// Parse a two-byte SAM tag, wrapping the validation error with the offending field.
fn parse_tag(field: &str) -> Result<SamTag> {
    SamTag::from_str(field).map_err(|e| anyhow::anyhow!("invalid retag tag {field:?}: {e}"))
}

impl FromStr for RetagOp {
    type Err = anyhow::Error;

    /// Parse one operation from `SRC::copy::DST`, `SRC::move::DST`, or `SRC::delete`.
    ///
    /// The `::` separator splits the fields; each `SRC`/`DST` is validated against
    /// the SAM aux-tag pattern via [`SamTag`]. Self-referential `copy`/`move`
    /// (`X::copy::X`, `X::move::X`) are rejected: the first is a no-op typo, the
    /// second silently deletes `X`.
    fn from_str(s: &str) -> Result<Self> {
        let fields: Vec<&str> = s.split("::").collect();
        match fields.as_slice() {
            [src, "delete"] => Ok(RetagOp::Delete { src: parse_tag(src)? }),
            [src, "copy", dst] => {
                let (src, dst) = (parse_tag(src)?, parse_tag(dst)?);
                ensure!(
                    src != dst,
                    "'{src}::copy::{dst}' is a no-op: source and destination are the same tag"
                );
                Ok(RetagOp::Copy { src, dst })
            }
            [src, "move", dst] => {
                let (src, dst) = (parse_tag(src)?, parse_tag(dst)?);
                ensure!(
                    src != dst,
                    "'{src}::move::{dst}' would silently delete '{src}': \
                     source and destination are the same tag"
                );
                Ok(RetagOp::Move { src, dst })
            }
            [_, "delete", _] => {
                bail!("'delete' takes no destination tag; use 'SRC::delete', got: {s:?}")
            }
            [_, op, _] => bail!(
                "unknown operation {op:?}; expected 'copy' or 'move' in 'SRC::op::DST', got: {s:?}"
            ),
            _ => bail!(
                "could not parse retag operation {s:?}; \
                 expected 'SRC::copy::DST', 'SRC::move::DST', or 'SRC::delete'"
            ),
        }
    }
}

/// Per-operation application counts, accumulated across all records.
///
/// One [`OpCounts`] tracks a single [`RetagOp`]. `records_applied` counts records
/// where the source tag was present (so the operation did something);
/// `src_missing` counts records where it was absent (the operation was skipped);
/// `dst_overwritten` counts `copy`/`move` records whose destination already held a
/// value that was replaced (always `0` for `delete`).
#[derive(Debug, Clone, Copy, Default, PartialEq, Eq)]
pub(crate) struct OpCounts {
    /// Records where the source tag was present and the operation was applied.
    pub records_applied: u64,
    /// Records (copy/move) where an existing destination value was overwritten.
    pub dst_overwritten: u64,
    /// Records where the source tag was absent, so the operation was skipped.
    pub src_missing: u64,
}

/// Copy `src`'s aux value verbatim into `dst`, overwriting any existing `dst`.
///
/// Returns `true` when `src` was present (and thus copied). The source value bytes
/// are copied to an owned buffer before the record is mutated, so the immutable
/// borrow is released before [`remove_tag`]/[`append_raw_tag`] run.
fn copy_tag(record: &mut RawRecord, src: SamTag, dst: SamTag, counts: &mut OpCounts) -> bool {
    // Grab the source entry's on-disk (type byte, value bytes) verbatim. `value_bytes`
    // is the exact encoding (a Z value keeps its NUL, a B array its subtype+count+elems),
    // so the destination entry is byte-identical apart from the tag identifier.
    let src_entry = record
        .tags()
        .iter()
        .find(|entry| entry.tag == *src)
        .map(|entry| (entry.type_byte, entry.value_bytes.to_vec()));

    let Some((type_byte, value)) = src_entry else {
        counts.src_missing += 1;
        return false;
    };

    if record.tags().contains(dst) {
        remove_tag(record.as_mut_vec(), dst);
        counts.dst_overwritten += 1;
    }
    append_raw_tag(record.as_mut_vec(), dst, type_byte, &value);
    counts.records_applied += 1;
    true
}

/// Apply a single [`RetagOp`] to one record, updating its [`OpCounts`] in place.
///
/// Operations are pure per-record edits: `copy` adds `dst` and keeps `src`, `move`
/// additionally drops `src`, and `delete` drops `src`. A record missing the source
/// tag is left unchanged and counted in `src_missing`.
pub(crate) fn apply_op(record: &mut RawRecord, op: RetagOp, counts: &mut OpCounts) {
    match op {
        RetagOp::Copy { src, dst } => {
            copy_tag(record, src, dst, counts);
        }
        RetagOp::Move { src, dst } => {
            // The parser rejects `src == dst`, so removing `src` never touches the
            // just-written `dst`.
            if copy_tag(record, src, dst, counts) {
                remove_tag(record.as_mut_vec(), src);
            }
        }
        RetagOp::Delete { src } => {
            if record.tags().contains(src) {
                remove_tag(record.as_mut_vec(), src);
                counts.records_applied += 1;
            } else {
                counts.src_missing += 1;
            }
        }
    }
}

/// One metrics row per operation, written to the optional `--metrics` TSV.
#[derive(Debug, Clone, Default, Serialize, Deserialize)]
pub struct RetagMetric {
    /// The operation as written on the command line, e.g. `RX::copy::BX`.
    pub operation: String,
    /// The operation keyword: `copy`, `move`, or `delete`.
    pub kind: String,
    /// Records where the source tag was present and the operation was applied.
    pub records_applied: u64,
    /// Records (copy/move) where an existing destination value was overwritten.
    pub dst_overwritten: u64,
    /// Records where the source tag was absent, so the operation was skipped.
    pub src_missing: u64,
}

impl RetagMetric {
    /// Build a metrics row from an operation and its accumulated counts.
    pub(crate) fn from_counts(op: RetagOp, counts: &OpCounts) -> Self {
        Self {
            operation: op.to_string(),
            kind: op.kind().to_string(),
            records_applied: counts.records_applied,
            dst_overwritten: counts.dst_overwritten,
            src_missing: counts.src_missing,
        }
    }
}

impl fgumi_metrics::Metric for RetagMetric {
    fn metric_name() -> &'static str {
        "retag"
    }
}

/// Projected options the chain builder needs for the retag stage.
///
/// A free-standing struct (mirrors filter/group) rather than `Some(self.clone())`,
/// because `Retag` does not derive `Clone`. Compression/threading come from the
/// [`crate::pipeline::chains::SingleStageContext`], so they are not carried here.
#[derive(Debug, Clone)]
pub struct RetagOptions {
    /// Tag-rewrite operations, applied left-to-right per record.
    pub operations: Vec<RetagOp>,
    /// Optional per-operation metrics TSV.
    pub metrics: Option<PathBuf>,
}

/// Shared per-run state for a retag chain, built by [`RetagOptions::setup_pipeline`].
pub(crate) struct RetagPipelineSetup {
    /// Per-thread per-operation counters, positionally indexed by operation.
    pub(crate) collected_metrics: Arc<PerThreadAccumulator<Vec<OpCounts>>>,
    /// Global processed-record counter — the finalize hook's `record_count`.
    pub(crate) progress_counter: Arc<AtomicU64>,
}

/// Captures the retag chain process closure needs. retag's transform ignores the
/// header, so — unlike filter's `process_captures` — this takes no header argument.
pub(crate) struct RetagProcessCaptures {
    /// The operations to apply per record, in order (shared, read-only).
    pub(crate) operations: Arc<Vec<RetagOp>>,
    /// Global processed-record counter, incremented once per record.
    pub(crate) progress: Arc<AtomicU64>,
}

impl RetagOptions {
    /// Build the shared per-thread accumulator + progress counter. Nothing here
    /// is fallible (contrast filter, which loads a reference), so no `Result`.
    pub(crate) fn setup_pipeline(&self, num_threads: usize) -> RetagPipelineSetup {
        RetagPipelineSetup {
            collected_metrics: PerThreadAccumulator::<Vec<OpCounts>>::new(num_threads),
            progress_counter: Arc::new(AtomicU64::new(0)),
        }
    }

    /// Build the process-closure captures from these options and the setup.
    pub(crate) fn process_captures(&self, setup: &RetagPipelineSetup) -> RetagProcessCaptures {
        RetagProcessCaptures {
            operations: Arc::new(self.operations.clone()),
            progress: Arc::clone(&setup.progress_counter),
        }
    }
}

/// Rewrite SAM aux tags by copying, moving, or deleting them by name.
#[derive(Debug, Parser)]
#[command(
    name = "retag",
    about = "\x1b[38;5;166m[UTILITIES]\x1b[0m      \x1b[36mRewrite SAM tags (copy/move/delete)\x1b[0m",
    long_about = r#"
Rewrite SAM auxiliary tags from one name to another.

fgumi is intentionally opinionated about tag names (RX, MI, ...). retag is the sanctioned
escape hatch for interoperating with data that uses vendor- or lab-specific tags, or for
handing fgumi output to downstream tools that expect different names. It is strictly
tag-to-tag on the records: read names and sequence are never touched. It does synthesize an
@HD line when the input lacks one and always adds an @PG provenance record to the header.

Operations are positional and applied left-to-right, per record:

  SRC::copy::DST   copy SRC's aux value verbatim (type byte + bytes) into DST; SRC stays;
                   an existing DST is overwritten
  SRC::move::DST   sugar for SRC::copy::DST then SRC::delete
  SRC::delete      drop SRC

Because operations apply in order, they compose:

  # fan one tag out to two, then drop the original
  fgumi retag -i in.bam -o out.bam RX::copy::BX RX::copy::CB RX::delete

  # chain through
  fgumi retag -i in.bam -o out.bam RX::move::BX BX::move::CB

Each SRC/DST is validated against the SAM aux-tag pattern [A-Za-z][A-Za-z0-9]. A record
missing SRC simply skips that operation. Output is written in the same order as the input.

With --metrics, one row per operation is written (operation, kind, records_applied,
dst_overwritten, src_missing). A warning is logged for any operation that matched zero
records, which catches tag typos.
"#
)]
pub struct Retag {
    /// Input and output BAM files.
    #[command(flatten)]
    pub io: BamIoOptions,

    /// Tag-rewrite operations, applied left-to-right per record:
    /// `SRC::copy::DST`, `SRC::move::DST`, or `SRC::delete`.
    #[arg(value_name = "SRC::op::DST", required = true, num_args = 1..)]
    pub operations: Vec<RetagOp>,

    /// Optional TSV file for per-operation metrics.
    ///
    /// The per-operation counts are summed across the pipeline's worker threads,
    /// so the metrics are identical regardless of `--threads`. (This differs from
    /// `clip`, whose per-read metrics are not summable and so are rejected under
    /// `--threads`.)
    #[arg(short = 'M', long = "metrics")]
    pub metrics: Option<PathBuf>,

    /// Compression options for the output BAM.
    #[command(flatten)]
    pub compression: CompressionOptions,

    /// Threading options. retag always runs on the declarative chain builder (the
    /// typed-step pipeline) via `execute_chain`; `--threads N` only sets the
    /// worker count. Absent, the chain runs at a single worker.
    #[command(flatten)]
    pub threading: ThreadingOptions,

    /// Scheduler and pipeline stats options.
    #[command(flatten)]
    pub scheduler_opts: SchedulerOptions,

    /// Pipeline queue memory options.
    #[command(flatten)]
    pub queue_memory: QueueMemoryOptions,
}

impl Retag {
    /// Project the clap fields into [`RetagOptions`] for the chain builder.
    #[must_use]
    pub fn to_retag_options(&self) -> RetagOptions {
        RetagOptions { operations: self.operations.clone(), metrics: self.metrics.clone() }
    }

    /// Reject a `--output`/`--metrics` path that resolves to the same file as
    /// `--input`, which would clobber the BAM being read.
    ///
    /// Paths are compared after [`std::fs::canonicalize`] resolves symlinks and
    /// `.`/`..`. A write target that does not exist yet cannot alias an existing
    /// input, so an unresolvable target is left for the writer to create.
    fn reject_write_aliasing_input(&self) -> Result<()> {
        let Ok(input_canon) = std::fs::canonicalize(&self.io.input) else {
            return Ok(());
        };
        let write_targets = std::iter::once((self.io.output.as_path(), "--output"))
            .chain(self.metrics.as_deref().map(|m| (m, "--metrics")));
        for (path, flag) in write_targets {
            if std::fs::canonicalize(path).is_ok_and(|canon| canon == input_canon) {
                bail!(
                    "{flag} '{}' is the same file as --input '{}'; choose a different path",
                    path.display(),
                    self.io.input.display()
                );
            }
        }
        Ok(())
    }
}

impl Retag {
    /// Run the retag on the declarative chain builder
    /// (`ChainSpec::single_stage(Stage::Retag, …)` → `build_for(spec)?.run()`).
    ///
    /// The chain is the only execution path: `execute` always dispatches here,
    /// with or without `--threads` (absent `--threads` runs the chain at a single
    /// worker). All user-facing diagnostics — the `Starting Retag` banner, the
    /// `OperationTimer`, the Input/Output/Operation lines, and the
    /// summary/warn/metrics finalize hooks — are emitted inside
    /// `ChainBuilder::add_retag`, matching the `add_filter`/`add_dedup`
    /// convention, so `execute_chain` itself only surfaces the CRC policy and
    /// runs the pipeline (mirrors `dedup`/`group` `execute_chain`).
    fn execute_chain(&self, command_line: &str) -> Result<()> {
        use crate::pipeline::chains::{
            ChainSpec, SingleStageContext, Stage, StageOptionsBag, build_for,
        };
        self.io.log_effective_check_crc();
        let stage_opts =
            StageOptionsBag { retag: Some(self.to_retag_options()), ..Default::default() };
        let ctx = SingleStageContext {
            io: &self.io,
            threading: &self.threading,
            compression: &self.compression,
            scheduler: &self.scheduler_opts,
            queue_memory: &self.queue_memory,
            command_line,
        };
        let spec = ChainSpec::single_stage(Stage::Retag, stage_opts, &ctx);
        build_for(spec)?.run()
    }
}

/// Sum the per-thread accumulator slots into a single, positionally-indexed
/// counts vector (one entry per operation).
///
/// Each worker merges into its own slot, so after the pipeline drains the final
/// per-operation totals are the element-wise sum across every populated slot; an
/// untouched slot is empty and contributes nothing (the `zip` stops early). The
/// counters are `u64`, so the sum is exact and order-independent — the result is
/// identical regardless of how the scheduler spread records across slots. Pulled
/// out into a shared free function — used by the chain finalize hooks
/// (`RetagFinalizeHook` and `RetagMetricsFinalizeHook`) to reduce the per-thread
/// accumulator — so the cross-slot merge is unit-testable with several
/// deliberately-populated slots, rather than left to scheduler luck.
pub(crate) fn sum_slot_counts(
    collected: &PerThreadAccumulator<Vec<OpCounts>>,
    n_ops: usize,
) -> Vec<OpCounts> {
    let mut counts = vec![OpCounts::default(); n_ops];
    for slot in collected.slots() {
        let slot = slot.lock();
        for (agg, part) in counts.iter_mut().zip(slot.iter()) {
            agg.records_applied += part.records_applied;
            agg.dst_overwritten += part.dst_overwritten;
            agg.src_missing += part.src_missing;
        }
    }
    counts
}

impl Command for Retag {
    fn execute(&self, command_line: &str) -> Result<()> {
        // stdin paths (`-` / `/dev/stdin`) are exempt from the existence check.
        self.io.validate()?;
        // Reject two outputs resolving to one destination before any writer opens.
        let mut outputs: Vec<(&std::path::Path, &str)> =
            vec![(self.io.output.as_path(), "--output")];
        if let Some(path) = &self.metrics {
            outputs.push((path.as_path(), "--metrics"));
        }
        reject_output_collisions(&outputs)?;
        // A write target that aliases the input would clobber the file being read.
        // Reject it up front, resolving symlinks/`.`/`..` via canonicalize, before any
        // reader or writer opens (matching `merge`'s output-vs-input guard).
        self.reject_write_aliasing_input()?;
        // The declarative chain is the only execution path. It emits its own CRC
        // log, `Starting Retag` banner + `OperationTimer`, and the
        // summary/warn/metrics finalize hooks inside `ChainBuilder::add_retag`
        // (running any of those here first would double-log and pre-consume
        // stdin), so `execute` does only the pre-flight validation above and then
        // dispatches. Absent `--threads` runs the chain at a single worker.
        self.execute_chain(command_line)
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use fgumi_bam_io::{create_raw_bam_reader_with_opts, create_raw_bam_writer};
    use fgumi_raw_bam::{RawTagsView, SamBuilder, aux_data_slice};
    use noodles::sam::Header;
    use rstest::rstest;

    fn tag(s: &str) -> SamTag {
        s.parse().expect("valid tag")
    }

    // ── projector parity (to_retag_options) ────────────────────────────────

    #[test]
    fn to_retag_options_carries_every_tuning_flag() {
        let cmd = Retag::try_parse_from([
            "retag",
            "-i",
            "in.bam",
            "-o",
            "out.bam",
            "-M",
            "m.tsv",
            "RX::copy::BX",
            "BX::move::CB",
        ])
        .expect("parse");
        let o = cmd.to_retag_options();
        assert_eq!(o.operations, cmd.operations);
        assert_eq!(o.metrics.as_deref(), Some(std::path::Path::new("m.tsv")));
    }

    #[test]
    fn to_retag_options_carries_defaults() {
        let cmd = Retag::try_parse_from(["retag", "-i", "in.bam", "-o", "out.bam", "RX::delete"])
            .expect("parse");
        let o = cmd.to_retag_options();
        assert_eq!(o.operations, cmd.operations);
        assert_eq!(o.metrics, None);
    }

    /// Read a tag's exact on-disk `(type_byte, value_bytes)`, for verbatim-copy assertions.
    fn raw_entry(record: &RawRecord, t: SamTag) -> Option<(u8, Vec<u8>)> {
        RawTagsView::new(aux_data_slice(record.as_ref()))
            .iter()
            .find(|e| e.tag == *t)
            .map(|e| (e.type_byte, e.value_bytes.to_vec()))
    }

    /// Apply one op to a record starting from zeroed counts; return the counts.
    fn apply_one(record: &mut RawRecord, op: RetagOp) -> OpCounts {
        let mut counts = OpCounts::default();
        apply_op(record, op, &mut counts);
        counts
    }

    #[rstest]
    #[case::copy("RX::copy::BX", RetagOp::Copy { src: tag("RX"), dst: tag("BX") })]
    #[case::move_("RX::move::BX", RetagOp::Move { src: tag("RX"), dst: tag("BX") })]
    #[case::delete("RX::delete", RetagOp::Delete { src: tag("RX") })]
    #[case::lowercase_local_tag("ob::move::RX", RetagOp::Move { src: tag("ob"), dst: tag("RX") })]
    fn parses_valid_operations(#[case] input: &str, #[case] expected: RetagOp) {
        assert_eq!(input.parse::<RetagOp>().expect("should parse"), expected);
    }

    #[rstest]
    #[case::self_copy("RX::copy::RX")]
    #[case::self_move("RX::move::RX")]
    #[case::unknown_op("RX::rename::BX")]
    #[case::delete_with_dst("RX::delete::BX")]
    #[case::copy_missing_dst("RX::copy")]
    #[case::bare_tag("RX")]
    #[case::empty("")]
    #[case::bad_src_tag("R::copy::BX")]
    #[case::bad_dst_tag("RX::copy::B")]
    #[case::three_char_tag("RXX::copy::BX")]
    fn rejects_invalid_operations(#[case] input: &str) {
        assert!(input.parse::<RetagOp>().is_err(), "expected {input:?} to be rejected");
    }

    #[test]
    fn kind_and_src_report_the_operation() {
        assert_eq!("RX::copy::BX".parse::<RetagOp>().unwrap().kind(), "copy");
        assert_eq!("RX::move::BX".parse::<RetagOp>().unwrap().kind(), "move");
        assert_eq!("RX::delete".parse::<RetagOp>().unwrap().kind(), "delete");
        assert_eq!("RX::delete".parse::<RetagOp>().unwrap().src(), tag("RX"));
    }

    // ── apply engine ───────────────────────────────────────────────────────

    /// Build a single-record BAM record carrying one string tag.
    fn record_with(t: &str, value: &[u8]) -> RawRecord {
        let mut b = SamBuilder::new();
        b.read_name(b"r1").add_string_tag(tag(t), value);
        b.build()
    }

    fn record_with_rx(value: &[u8]) -> RawRecord {
        record_with("RX", value)
    }

    /// Apply a sequence of ops to a record, returning the per-op counts.
    fn apply_all(record: &mut RawRecord, ops: &[RetagOp]) -> Vec<OpCounts> {
        let mut counts = vec![OpCounts::default(); ops.len()];
        for (op, c) in ops.iter().zip(counts.iter_mut()) {
            apply_op(record, *op, c);
        }
        counts
    }

    #[test]
    fn copy_adds_destination_and_keeps_source() {
        let mut record = record_with_rx(b"ACGT");
        let counts = apply_one(&mut record, "RX::copy::BX".parse().unwrap());

        assert_eq!(record.tags().find_string(tag("RX")), Some(b"ACGT".as_ref()), "source kept");
        assert_eq!(record.tags().find_string(tag("BX")), Some(b"ACGT".as_ref()), "dest written");
        assert_eq!(counts, OpCounts { records_applied: 1, dst_overwritten: 0, src_missing: 0 });
    }

    #[test]
    fn copy_preserves_an_integer_value_verbatim() {
        // Source NM is a small int; `add_int_tag` encodes it with the smallest type.
        let mut b = SamBuilder::new();
        b.read_name(b"r1").add_int_tag(tag("NM"), 2);
        let mut record = b.build();
        let before = raw_entry(&record, tag("NM")).expect("NM present");

        apply_one(&mut record, "NM::copy::XN".parse().unwrap());

        assert_eq!(
            raw_entry(&record, tag("XN")),
            Some(before),
            "type byte + value copied verbatim"
        );
    }

    #[test]
    fn copy_preserves_a_b_array_value_verbatim() {
        let mut b = SamBuilder::new();
        b.read_name(b"r1").add_array_i32(tag("id"), &[-100_000, 0, 100_000]);
        let mut record = b.build();
        let before = raw_entry(&record, tag("id")).expect("id present");

        apply_one(&mut record, "id::copy::XA".parse().unwrap());

        assert_eq!(raw_entry(&record, tag("XA")), Some(before), "B array copied verbatim");
    }

    #[test]
    fn copy_overwrites_an_existing_destination() {
        let mut b = SamBuilder::new();
        b.read_name(b"r1").add_string_tag(tag("RX"), b"ACGT").add_string_tag(tag("BX"), b"OLD");
        let mut record = b.build();

        let counts = apply_one(&mut record, "RX::copy::BX".parse().unwrap());

        assert_eq!(record.tags().find_string(tag("BX")), Some(b"ACGT".as_ref()), "dest replaced");
        assert_eq!(counts, OpCounts { records_applied: 1, dst_overwritten: 1, src_missing: 0 });
    }

    #[test]
    fn copy_skips_a_record_missing_the_source() {
        let mut record = record_with("MI", b"7");
        let counts = apply_one(&mut record, "RX::copy::BX".parse().unwrap());

        assert!(!record.tags().contains(tag("BX")), "no destination written");
        assert_eq!(counts, OpCounts { records_applied: 0, dst_overwritten: 0, src_missing: 1 });
    }

    #[test]
    fn move_writes_destination_and_drops_source() {
        let mut record = record_with_rx(b"ACGT");
        let counts = apply_one(&mut record, "RX::move::BX".parse().unwrap());

        assert!(!record.tags().contains(tag("RX")), "source removed");
        assert_eq!(record.tags().find_string(tag("BX")), Some(b"ACGT".as_ref()), "dest written");
        assert_eq!(counts, OpCounts { records_applied: 1, dst_overwritten: 0, src_missing: 0 });
    }

    #[test]
    fn move_skips_a_record_missing_the_source() {
        let mut record = record_with_rx(b"ACGT");
        let counts = apply_one(&mut record, "MI::move::BX".parse().unwrap());

        assert!(record.tags().contains(tag("RX")), "untouched");
        assert!(!record.tags().contains(tag("BX")));
        assert_eq!(counts, OpCounts { records_applied: 0, dst_overwritten: 0, src_missing: 1 });
    }

    #[test]
    fn delete_removes_a_present_source() {
        let mut record = record_with_rx(b"ACGT");
        let counts = apply_one(&mut record, "RX::delete".parse().unwrap());

        assert!(!record.tags().contains(tag("RX")));
        assert_eq!(counts, OpCounts { records_applied: 1, dst_overwritten: 0, src_missing: 0 });
    }

    #[test]
    fn delete_skips_a_missing_source() {
        let mut record = record_with("MI", b"7");
        let counts = apply_one(&mut record, "RX::delete".parse().unwrap());

        assert!(record.tags().contains(tag("MI")), "untouched");
        assert_eq!(counts, OpCounts { records_applied: 0, dst_overwritten: 0, src_missing: 1 });
    }

    #[test]
    fn ordered_ops_fan_a_tag_out_then_drop_the_original() {
        // RX::copy::BX RX::copy::CB RX::delete — fan one tag to two, drop the source.
        let ops: Vec<RetagOp> = ["RX::copy::BX", "RX::copy::CB", "RX::delete"]
            .iter()
            .map(|s| s.parse().unwrap())
            .collect();
        let mut record = record_with_rx(b"ACGT");
        apply_all(&mut record, &ops);

        assert!(!record.tags().contains(tag("RX")), "original dropped");
        assert_eq!(record.tags().find_string(tag("BX")), Some(b"ACGT".as_ref()));
        assert_eq!(record.tags().find_string(tag("CB")), Some(b"ACGT".as_ref()));
    }

    #[test]
    fn ordered_moves_chain_through_intermediate_tags() {
        // RX::move::BX BX::move::CB — chain the value through and leave only CB.
        let ops: Vec<RetagOp> =
            ["RX::move::BX", "BX::move::CB"].iter().map(|s| s.parse().unwrap()).collect();
        let mut record = record_with_rx(b"ACGT");
        apply_all(&mut record, &ops);

        assert!(!record.tags().contains(tag("RX")));
        assert!(!record.tags().contains(tag("BX")));
        assert_eq!(record.tags().find_string(tag("CB")), Some(b"ACGT".as_ref()));
    }

    // ── end-to-end: real BAM round-trip through the Retag command ────────────

    /// Write raw records to a BAM at `path` under a default (reference-less) header.
    fn write_bam(path: &std::path::Path, records: &[RawRecord]) {
        let header = noodles::sam::Header::default();
        let mut writer = create_raw_bam_writer(path, &header, 1, 1).expect("open writer");
        for r in records {
            writer.write_raw_record(r.as_ref()).expect("write record");
        }
        writer.finish().expect("finish writer");
    }

    /// Read all records back from a BAM at `path`.
    fn read_bam(path: &std::path::Path) -> Vec<RawRecord> {
        let (mut reader, _header) =
            create_raw_bam_reader_with_opts(path, 1, fgumi_bam_io::PipelineReaderOpts::default())
                .expect("open reader");
        let mut out = Vec::new();
        let mut rec = RawRecord::new();
        while reader.read_record(&mut rec).expect("read record") != 0 {
            out.push(rec.clone());
        }
        out
    }

    #[test]
    fn end_to_end_rewrites_tags_and_writes_metrics() {
        let dir = tempfile::TempDir::new().expect("temp dir");
        let input = dir.path().join("in.bam");
        let output = dir.path().join("out.bam");
        let metrics = dir.path().join("retag.tsv");

        // r1 carries RX + MI; r2 carries only RX (MI absent).
        let r1 = {
            let mut b = SamBuilder::new();
            b.read_name(b"r1").add_string_tag(tag("RX"), b"ACGT").add_string_tag(tag("MI"), b"1");
            b.build()
        };
        let r2 = {
            let mut b = SamBuilder::new();
            b.read_name(b"r2").add_string_tag(tag("RX"), b"TTTT");
            b.build()
        };
        write_bam(&input, &[r1, r2]);

        // Copy RX->BX (matches both), delete MI (matches one), delete a never-present tag.
        let cmd = Retag {
            io: BamIoOptions {
                input,
                output: output.clone(),
                async_reader: false,
                check_crc: false,
                no_check_crc: false,
            },
            operations: ["RX::copy::BX", "MI::delete", "ZZ::delete"]
                .iter()
                .map(|s| s.parse().unwrap())
                .collect(),
            metrics: Some(metrics.clone()),
            compression: CompressionOptions { compression_level: 1 },
            threading: ThreadingOptions::none(),
            scheduler_opts: SchedulerOptions::default(),
            queue_memory: QueueMemoryOptions::default(),
        };
        cmd.execute("fgumi retag").expect("retag should succeed");

        // Output preserves record order and rewrites tags per record.
        let out = read_bam(&output);
        assert_eq!(out.len(), 2);
        assert_eq!(out[0].tags().find_string(tag("RX")), Some(b"ACGT".as_ref()), "r1 RX kept");
        assert_eq!(out[0].tags().find_string(tag("BX")), Some(b"ACGT".as_ref()), "r1 BX copied");
        assert!(!out[0].tags().contains(tag("MI")), "r1 MI deleted");
        assert_eq!(out[1].tags().find_string(tag("RX")), Some(b"TTTT".as_ref()), "r2 RX kept");
        assert_eq!(out[1].tags().find_string(tag("BX")), Some(b"TTTT".as_ref()), "r2 BX copied");

        // Metrics: one row per op, in order, with the expected counts.
        let tsv = std::fs::read_to_string(&metrics).expect("read metrics");
        let lines: Vec<&str> = tsv.lines().collect();
        assert_eq!(lines[0], "operation\tkind\trecords_applied\tdst_overwritten\tsrc_missing");
        assert_eq!(lines[1], "RX::copy::BX\tcopy\t2\t0\t0");
        assert_eq!(lines[2], "MI::delete\tdelete\t1\t0\t1");
        assert_eq!(lines[3], "ZZ::delete\tdelete\t0\t0\t2");
    }

    /// Build a `Retag` command reading `input` with the given write targets.
    fn retag_cmd(input: PathBuf, output: PathBuf, metrics: Option<PathBuf>) -> Retag {
        Retag {
            io: BamIoOptions {
                input,
                output,
                async_reader: false,
                check_crc: false,
                no_check_crc: false,
            },
            operations: vec!["RX::copy::BX".parse().unwrap()],
            metrics,
            compression: CompressionOptions { compression_level: 1 },
            threading: ThreadingOptions::none(),
            scheduler_opts: SchedulerOptions::default(),
            queue_memory: QueueMemoryOptions::default(),
        }
    }

    /// Build a `Retag` command with explicit operations and threading, for the
    /// worker-count invariance tests (single worker vs multi-worker chain).
    fn retag_cmd_with(
        input: PathBuf,
        output: PathBuf,
        metrics: Option<PathBuf>,
        operations: &[&str],
        threading: ThreadingOptions,
    ) -> Retag {
        Retag {
            io: BamIoOptions {
                input,
                output,
                async_reader: false,
                check_crc: false,
                no_check_crc: false,
            },
            operations: operations.iter().map(|s| s.parse().unwrap()).collect(),
            metrics,
            compression: CompressionOptions { compression_level: 1 },
            threading,
            scheduler_opts: SchedulerOptions::default(),
            queue_memory: QueueMemoryOptions::default(),
        }
    }

    #[rstest]
    #[case::metrics_aliases_input(true)]
    #[case::output_aliases_input(false)]
    fn rejects_a_write_target_that_aliases_the_input(#[case] alias_metrics: bool) {
        let dir = tempfile::TempDir::new().expect("temp dir");
        let input = dir.path().join("in.bam");
        write_bam(&input, &[record_with_rx(b"ACGT")]);

        // Point either --metrics or --output at the existing input file.
        let cmd = if alias_metrics {
            retag_cmd(input.clone(), dir.path().join("out.bam"), Some(input.clone()))
        } else {
            retag_cmd(input.clone(), input.clone(), None)
        };

        let err = cmd.execute("fgumi retag").expect_err("aliasing the input must be rejected");
        assert!(err.to_string().contains("same file as --input"), "unexpected error: {err}");

        // The guard runs before any writer opens, so the input BAM is untouched.
        let records = read_bam(&input);
        assert_eq!(records.len(), 1);
        assert_eq!(records[0].tags().find_string(tag("RX")), Some(b"ACGT".as_ref()));
    }

    // ── worker-count invariance: output independent of --threads ─────────────

    /// A header with one `@SQ` (`chr1`), so the mapped records in
    /// [`build_varied_records`] carry a valid reference id.
    fn parity_header() -> Header {
        use noodles::sam::header::record::value::Map;
        use noodles::sam::header::record::value::map::ReferenceSequence;
        use std::num::NonZeroUsize;
        Header::builder()
            .add_reference_sequence(
                b"chr1",
                Map::<ReferenceSequence>::new(NonZeroUsize::new(1000).expect("1000 is non-zero")),
            )
            .build()
    }

    /// Write raw records to a BAM at `path` under the given header.
    fn write_bam_with_header(path: &std::path::Path, header: &Header, records: &[RawRecord]) {
        let mut writer = create_raw_bam_writer(path, header, 1, 1).expect("open writer");
        for r in records {
            writer.write_raw_record(r.as_ref()).expect("write record");
        }
        writer.finish().expect("finish writer");
    }

    /// Read just the header back from a BAM at `path`.
    fn read_bam_header(path: &std::path::Path) -> Header {
        let (_reader, header) =
            create_raw_bam_reader_with_opts(path, 1, fgumi_bam_io::PipelineReaderOpts::default())
                .expect("open reader");
        header
    }

    /// Build 50 varied records spanning the record classes where a
    /// worker-count-dependent decode divergence could hide. Every record carries a
    /// distinct `RX`; read names encode the input index (`r0000`..`r0049`) so order
    /// can be asserted after a multi-worker run. `ZZ` is never present.
    ///
    /// - even `i` also carry `MI`;
    /// - every 5th record pre-carries `BX`, so `RX::copy::BX` overwrites it and the
    ///   `dst_overwritten` counter is exercised (and its cross-worker aggregation checked);
    /// - every 3rd record is a mapped, paired read (real CIGAR, position, `MC`/`RG`
    ///   tags, mate info) — some flagged secondary — so the per-record decode work
    ///   (CIGAR walk, aux/mate parsing) runs on real input, not just trivial
    ///   unmapped primaries.
    fn build_varied_records() -> Vec<RawRecord> {
        use fgumi_raw_bam::flags;
        (0..50u32)
            .map(|i| {
                let mut b = SamBuilder::new();
                let name = format!("r{i:04}");
                let rx = format!("ACGT{i:04}");
                b.read_name(name.as_bytes()).add_string_tag(tag("RX"), rx.as_bytes());
                if i % 2 == 0 {
                    b.add_string_tag(tag("MI"), i.to_string().as_bytes());
                }
                if i % 5 == 0 {
                    // Pre-existing destination for `RX::copy::BX` to overwrite.
                    b.add_string_tag(tag("BX"), b"OLD");
                }
                if i % 3 == 0 {
                    // Mapped, paired read: exercises the CIGAR walk + mate/aux parsing
                    // the per-record decode performs on real input.
                    let mut f = flags::PAIRED;
                    if i % 6 == 0 {
                        f |= flags::SECONDARY;
                    }
                    b.ref_id(0)
                        .pos(i as i32 * 7 + 1)
                        .flags(f)
                        .mate_ref_id(0)
                        .mate_pos(i as i32 * 7 + 50)
                        .cigar_ops(&[8 << 4]) // 8M
                        .sequence(b"ACGTACGT")
                        .qualities(&[30u8; 8])
                        .add_string_tag(tag("MC"), b"8M")
                        .add_string_tag(tag("RG"), b"A");
                } else {
                    // Unmapped primary.
                    b.sequence(b"ACGT").qualities(&[30u8; 4]);
                }
                b.build()
            })
            .collect()
    }

    #[rstest]
    #[case::threads_1(1)]
    #[case::threads_4(4)]
    fn output_is_independent_of_worker_count(#[case] threads: usize) {
        let dir = tempfile::TempDir::new().expect("temp dir");
        let input = dir.path().join("in.bam");
        write_bam_with_header(&input, &parity_header(), &build_varied_records());

        // copy (matches all, overwrites pre-existing BX on every 5th), move (matches
        // evens), delete a never-present tag.
        let ops = ["RX::copy::BX", "MI::move::CB", "ZZ::delete"];

        // No `--threads`: the chain at a single worker.
        let single_worker_out = dir.path().join("single_worker.bam");
        let single_worker_tsv = dir.path().join("single_worker.tsv");
        retag_cmd_with(
            input.clone(),
            single_worker_out.clone(),
            Some(single_worker_tsv.clone()),
            &ops,
            ThreadingOptions::none(),
        )
        .execute("fgumi retag")
        .expect("single-worker retag");

        let threaded_out = dir.path().join("threaded.bam");
        let threaded_tsv = dir.path().join("threaded.tsv");
        retag_cmd_with(
            input.clone(),
            threaded_out.clone(),
            Some(threaded_tsv.clone()),
            &ops,
            ThreadingOptions::new(threads),
        )
        .execute("fgumi retag")
        .expect("threaded retag");

        // Decoded records are identical in order and content. (Compressed BAM bytes
        // may differ — the pipeline writer chunks BGZF blocks differently — so we
        // compare decoded records, not file bytes.)
        let single_worker_recs = read_bam(&single_worker_out);
        let threaded_recs = read_bam(&threaded_out);
        assert_eq!(single_worker_recs.len(), threaded_recs.len(), "record count differs");
        for (i, (a, b)) in single_worker_recs.iter().zip(&threaded_recs).enumerate() {
            assert_eq!(a.as_ref(), b.as_ref(), "record {i} differs between worker counts");
        }

        // The synthesized @HD and chained @PG header must also match across worker
        // counts. Both runs pass the same command line, so the headers are byte-identical.
        assert_eq!(
            read_bam_header(&single_worker_out),
            read_bam_header(&threaded_out),
            "output header differs between worker counts"
        );

        // Per-operation metrics TSV is byte-identical across worker counts (u64 counts
        // sum commutatively, rows emitted in operation order).
        let single_worker_metrics =
            std::fs::read_to_string(&single_worker_tsv).expect("single-worker tsv");
        let threaded_metrics = std::fs::read_to_string(&threaded_tsv).expect("threaded tsv");
        assert_eq!(
            single_worker_metrics, threaded_metrics,
            "metrics TSV differ between worker counts"
        );

        // Guard that the corpus actually exercises the counters whose cross-worker
        // aggregation this test is here to check: `dst_overwritten > 0` (10 records
        // pre-carry BX) and the zero-match `delete` (drives the warning).
        let rx_row = threaded_metrics.lines().nth(1).expect("RX row");
        assert_eq!(rx_row, "RX::copy::BX\tcopy\t50\t10\t0", "dst_overwritten not aggregated");
        let last = threaded_metrics.lines().last().expect("metrics rows");
        assert_eq!(last, "ZZ::delete\tdelete\t0\t0\t50");
    }

    #[test]
    fn threaded_preserves_input_record_order() {
        let dir = tempfile::TempDir::new().expect("temp dir");
        let input = dir.path().join("in.bam");
        write_bam_with_header(&input, &parity_header(), &build_varied_records());

        let output = dir.path().join("out.bam");
        retag_cmd_with(input, output.clone(), None, &["RX::copy::BX"], ThreadingOptions::new(4))
            .execute("fgumi retag")
            .expect("threaded retag");

        // Read names encode the input index; they must come back in input order.
        let out = read_bam(&output);
        let got: Vec<Vec<u8>> = out.iter().map(|r| r.read_name().to_vec()).collect();
        let want: Vec<Vec<u8>> = (0..50u32).map(|i| format!("r{i:04}").into_bytes()).collect();
        assert_eq!(got, want, "records emitted out of input order under --threads");
    }

    #[test]
    fn threads_flag_coexists_with_positional_operations() {
        use clap::Parser;
        let cmd = Retag::try_parse_from([
            "retag",
            "-i",
            "in.bam",
            "-o",
            "out.bam",
            "--threads",
            "4",
            "--compression-level",
            "2",
            "RX::copy::BX",
            "MI::delete",
        ])
        .expect("--threads should coexist with positional SRC::op::DST args");
        assert_eq!(cmd.threading.threads, Some(4));
        assert_eq!(cmd.compression.compression_level, 2);
        assert_eq!(cmd.operations.len(), 2);
    }

    #[test]
    fn sum_slot_counts_merges_multiple_populated_slots() {
        // The end-to-end parity test proves output identity but cannot guarantee the
        // scheduler spread its records across more than one accumulator slot, so it
        // only probabilistically exercises the cross-slot merge. Pin that merge here
        // directly: populate two slots with different per-operation counts and leave a
        // third empty (writing slots by index, bypassing the thread-local slot
        // assignment), then assert the sum is element-wise across both populated slots
        // and that the empty slot contributes nothing.
        let acc = PerThreadAccumulator::<Vec<OpCounts>>::new(3);
        {
            let mut s0 = acc.slots()[0].lock();
            *s0 = vec![
                OpCounts { records_applied: 3, dst_overwritten: 1, src_missing: 2 },
                OpCounts { records_applied: 5, dst_overwritten: 0, src_missing: 4 },
            ];
        }
        {
            let mut s1 = acc.slots()[1].lock();
            *s1 = vec![
                OpCounts { records_applied: 7, dst_overwritten: 2, src_missing: 1 },
                OpCounts { records_applied: 0, dst_overwritten: 0, src_missing: 9 },
            ];
        }
        // acc.slots()[2] stays empty and must contribute nothing.

        let summed = sum_slot_counts(&acc, 2);
        assert_eq!(
            summed,
            vec![
                OpCounts { records_applied: 10, dst_overwritten: 3, src_missing: 3 },
                OpCounts { records_applied: 5, dst_overwritten: 0, src_missing: 13 },
            ],
            "cross-slot per-operation counts must sum element-wise"
        );
    }
}
