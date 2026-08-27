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
use std::io;
use std::path::PathBuf;
use std::str::FromStr;
use std::sync::Arc;

use anyhow::{Result, bail, ensure};
use clap::Parser;
use fgumi_bam_io::{
    ProgressTracker, create_bam_reader_for_pipeline_with_opts, create_raw_bam_reader_with_opts,
    create_raw_bam_writer,
};
use fgumi_raw_bam::{RawRecord, append_raw_tag, remove_tag};
use log::{info, warn};
use noodles::sam::Header;
use serde::{Deserialize, Serialize};

use crate::commands::command::Command;
use crate::commands::common::{
    BamIoOptions, CompressionOptions, QueueMemoryOptions, SchedulerOptions, ThreadingOptions,
    add_pg_record, build_pipeline_config, ensure_hd_record, reject_output_collisions,
    serialize_raw_bam_records,
};
use crate::grouper::SingleRawRecordGrouper;
use crate::logging::OperationTimer;
use crate::per_thread_accumulator::PerThreadAccumulator;
use crate::read_info::LibraryIndex;
use crate::sam::SamTag;
use crate::unified_pipeline::{GroupKeyConfig, Grouper, run_bam_pipeline_from_reader};

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
pub struct OpCounts {
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
pub fn apply_op(record: &mut RawRecord, op: &RetagOp, counts: &mut OpCounts) {
    match op {
        RetagOp::Copy { src, dst } => {
            copy_tag(record, *src, *dst, counts);
        }
        RetagOp::Move { src, dst } => {
            // The parser rejects `src == dst`, so removing `src` never touches the
            // just-written `dst`.
            if copy_tag(record, *src, *dst, counts) {
                remove_tag(record.as_mut_vec(), *src);
            }
        }
        RetagOp::Delete { src } => {
            if record.tags().contains(*src) {
                remove_tag(record.as_mut_vec(), *src);
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
    fn from_counts(op: RetagOp, counts: &OpCounts) -> Self {
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
    /// Works with `--threads`: the per-operation counts are summed across the
    /// pipeline's worker threads, so the metrics are identical to a
    /// single-threaded run. (This differs from `clip`, whose per-read metrics
    /// are not summable and so are rejected under `--threads`.)
    #[arg(short = 'M', long = "metrics")]
    pub metrics: Option<PathBuf>,

    /// Compression options for the output BAM.
    #[command(flatten)]
    pub compression: CompressionOptions,

    /// Threading options: `--threads N` routes through the unified pipeline.
    #[command(flatten)]
    pub threading: ThreadingOptions,

    /// Scheduler and pipeline stats options (used only in `--threads` mode).
    #[command(flatten)]
    pub scheduler_opts: SchedulerOptions,

    /// Pipeline queue memory options (used only in `--threads` mode).
    #[command(flatten)]
    pub queue_memory: QueueMemoryOptions,
}

impl Retag {
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
    /// Single-threaded fast path: a serial read → rewrite → write loop.
    ///
    /// Reached when `--threads` is absent. Opens its own record reader and
    /// writer (both single-thread), synthesizes the `@HD`/`@PG` header, and
    /// returns the record count plus the per-operation counts (positionally
    /// indexed against `self.operations`) for shared summary/metrics reporting.
    fn run_single_threaded(&self, command_line: &str) -> Result<(u64, Vec<OpCounts>)> {
        let (mut reader, header) =
            create_raw_bam_reader_with_opts(&self.io.input, 1, self.io.pipeline_reader_opts())?;
        // Synthesize @HD VN:1.6 SO:unsorted when absent, then chain a @PG for provenance.
        let header = ensure_hd_record(header)?;
        let header = add_pg_record(header, command_line)?;

        let mut writer =
            create_raw_bam_writer(&self.io.output, &header, 1, self.compression.compression_level)?;

        let mut counts = vec![OpCounts::default(); self.operations.len()];
        let mut record_count: u64 = 0;
        let progress = ProgressTracker::new("Processed records").with_interval(1_000_000);

        let mut record = RawRecord::new();
        while reader.read_record(&mut record)? != 0 {
            for (op, op_counts) in self.operations.iter().zip(counts.iter_mut()) {
                apply_op(&mut record, op, op_counts);
            }
            writer.write_raw_record(record.as_ref())?;
            record_count += 1;
            progress.log_if_needed(1);
        }
        progress.log_final();

        // Flush before summarizing so any write error surfaces here, not silently on drop.
        writer.finish()?;
        Ok((record_count, counts))
    }

    /// Multi-threaded path: route records through the unified 7-step pipeline.
    ///
    /// retag is a pure per-record transform, so it uses [`SingleRawRecordGrouper`]
    /// (one record per unit — no template batching) and applies the ops inside the
    /// pipeline's `process_fn`. The pipeline preserves input order, so output order
    /// matches the single-threaded path.
    ///
    /// Per-operation counts are accumulated into a [`PerThreadAccumulator`] and
    /// summed after the pipeline drains. Because [`OpCounts`] fields are `u64` and
    /// addition is commutative, the aggregated counts — and therefore the `--metrics`
    /// TSV — are identical to the single-threaded path regardless of thread count.
    /// This is a deliberate divergence from `clip`, which rejects `--metrics` with
    /// `--threads` because its metrics are non-summable per-read collections; retag's
    /// three summable counters aggregate cleanly, so `--metrics` stays available here.
    fn run_threaded(&self, num_threads: usize, command_line: &str) -> Result<(u64, Vec<OpCounts>)> {
        // The pipeline manages its own I/O lifecycle, so open the streaming reader
        // (header pre-consumed) here rather than the serial record reader.
        let (reader, header) = create_bam_reader_for_pipeline_with_opts(
            &self.io.input,
            self.io.pipeline_reader_opts(),
        )?;
        let header = ensure_hd_record(header)?;
        let header = add_pg_record(header, command_line)?;

        let mut pipeline_config = build_pipeline_config(
            &self.scheduler_opts,
            &self.compression,
            &self.queue_memory,
            &self.io,
            num_threads,
        )?;
        // Raw-byte mode: hand the grouper RawRecord items and skip the per-record
        // cell-barcode scan the default config would otherwise perform. retag never
        // groups by key (one record per group), so the library index the key would
        // carry is irrelevant — pass the default and skip the header scan, matching
        // clip's raw-passthrough setup.
        pipeline_config.group_key_config =
            Some(GroupKeyConfig::new_raw_no_cell(LibraryIndex::default()));

        let n_ops = self.operations.len();
        let ops = self.operations.clone();
        let collected = PerThreadAccumulator::<Vec<OpCounts>>::new(num_threads);
        let collected_for_process = Arc::clone(&collected);

        // One record per group: retag is a pure per-record transform, so no template
        // batching is needed.
        let grouper_fn = |_header: &Header| {
            Box::new(SingleRawRecordGrouper::new()) as Box<dyn Grouper<Group = RawRecord> + Send>
        };
        // Apply every op to the record and merge the resulting counts into this
        // worker's accumulator slot.
        let process_fn = move |mut record: RawRecord| -> io::Result<RawRecord> {
            collected_for_process.with_slot(|slot| {
                // Lazily size the slot to the operation count on first use.
                if slot.is_empty() {
                    slot.resize(n_ops, OpCounts::default());
                }
                for (op, op_counts) in ops.iter().zip(slot.iter_mut()) {
                    apply_op(&mut record, op, op_counts);
                }
            });
            Ok(record)
        };

        // Serialize: write the rewritten record. The pipeline owns the
        // "Processed records" heartbeat (its output stage logs it at each
        // million-record boundary), so this closure must not log progress or it
        // would double the heartbeat line.
        let serialize_fn =
            move |record: RawRecord, _header: &Header, output: &mut Vec<u8>| -> io::Result<u64> {
                serialize_raw_bam_records(std::slice::from_ref(&record), output)
            };

        let record_count = run_bam_pipeline_from_reader(
            pipeline_config,
            reader,
            header,
            &self.io.output,
            None, // reuse the input header for output
            grouper_fn,
            process_fn,
            serialize_fn,
        )?;

        Ok((record_count, sum_slot_counts(&collected, n_ops)))
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
/// out of [`Retag::run_threaded`] so the cross-slot merge is unit-testable with
/// several deliberately-populated slots, rather than left to scheduler luck.
fn sum_slot_counts(collected: &PerThreadAccumulator<Vec<OpCounts>>, n_ops: usize) -> Vec<OpCounts> {
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
        // Surface the effective CRC-verification policy once, at run start, honoring
        // the `--check-crc` doc's promise of a per-run `CRC verify:` line. Both modes
        // take their policy from `pipeline_reader_opts()` / `build_pipeline_config`.
        self.io.log_effective_check_crc();

        let timer = OperationTimer::new("Rewriting tags");
        info!("Starting Retag");
        info!("Input: {}", self.io.input.display());
        info!("Output: {}", self.io.output.display());
        for op in &self.operations {
            info!("Operation: {op}");
        }

        // Route on --threads: absent → serial fast path; present → unified pipeline.
        // Both return the record count and per-operation counts for shared reporting.
        let (record_count, counts) = match self.threading.threads {
            None => self.run_single_threaded(command_line)?,
            Some(threads) => self.run_threaded(threads, command_line)?,
        };

        // Warn on operations that never matched — the usual sign of a mistyped source tag.
        for (op, op_counts) in self.operations.iter().zip(&counts) {
            if op_counts.records_applied == 0 {
                warn!(
                    "operation '{op}' matched zero records: no record carried the source tag '{}'",
                    op.src()
                );
            }
        }

        if let Some(path) = &self.metrics {
            let rows: Vec<RetagMetric> = self
                .operations
                .iter()
                .zip(&counts)
                .map(|(op, c)| RetagMetric::from_counts(*op, c))
                .collect();
            fgumi_metrics::write_metrics(path, &rows, "retag")?;
            info!("Wrote metrics to: {}", path.display());
        }

        info!("=== Summary ===");
        info!("Records processed: {record_count}");
        for (op, op_counts) in self.operations.iter().zip(&counts) {
            info!(
                "{op}: applied={} overwritten={} missing={}",
                op_counts.records_applied, op_counts.dst_overwritten, op_counts.src_missing
            );
        }

        timer.log_completion(record_count);
        Ok(())
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use fgumi_raw_bam::{RawTagsView, SamBuilder, aux_data_slice};
    use rstest::rstest;

    fn tag(s: &str) -> SamTag {
        s.parse().expect("valid tag")
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
        apply_op(record, &op, &mut counts);
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
            apply_op(record, op, c);
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
    /// serial-vs-pipeline parity tests.
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

    // ── threaded pipeline: parity with the single-threaded path ──────────────

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

    /// Build 50 varied records spanning the record classes where a serial-vs-threaded
    /// decode divergence could hide. Every record carries a distinct `RX`; read names
    /// encode the input index (`r0000`..`r0049`) so order can be asserted after a
    /// parallel run. `ZZ` is never present.
    ///
    /// - even `i` also carry `MI`;
    /// - every 5th record pre-carries `BX`, so `RX::copy::BX` overwrites it and the
    ///   `dst_overwritten` counter is exercised (and its threaded aggregation checked);
    /// - every 3rd record is a mapped, paired read (real CIGAR, position, `MC`/`RG`
    ///   tags, mate info) — some flagged secondary — so the threaded path's extra
    ///   per-record decode work (CIGAR walk, aux/mate parsing) runs on real input,
    ///   not just trivial unmapped primaries.
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
                    // the threaded decode performs but the serial loop skips.
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
    fn threaded_output_matches_single_threaded(#[case] threads: usize) {
        let dir = tempfile::TempDir::new().expect("temp dir");
        let input = dir.path().join("in.bam");
        write_bam_with_header(&input, &parity_header(), &build_varied_records());

        // copy (matches all, overwrites pre-existing BX on every 5th), move (matches
        // evens), delete a never-present tag.
        let ops = ["RX::copy::BX", "MI::move::CB", "ZZ::delete"];

        let serial_out = dir.path().join("serial.bam");
        let serial_tsv = dir.path().join("serial.tsv");
        retag_cmd_with(
            input.clone(),
            serial_out.clone(),
            Some(serial_tsv.clone()),
            &ops,
            ThreadingOptions::none(),
        )
        .execute("fgumi retag")
        .expect("serial retag");

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
        let serial_recs = read_bam(&serial_out);
        let threaded_recs = read_bam(&threaded_out);
        assert_eq!(serial_recs.len(), threaded_recs.len(), "record count differs");
        for (i, (a, b)) in serial_recs.iter().zip(&threaded_recs).enumerate() {
            assert_eq!(a.as_ref(), b.as_ref(), "record {i} differs between modes");
        }

        // The synthesized @HD and chained @PG header must also match between modes.
        // Both runs pass the same command line, so the headers are byte-identical.
        assert_eq!(
            read_bam_header(&serial_out),
            read_bam_header(&threaded_out),
            "output header differs between modes"
        );

        // Per-operation metrics TSV is byte-identical between modes (u64 counts sum
        // commutatively, rows emitted in operation order).
        let serial_metrics = std::fs::read_to_string(&serial_tsv).expect("serial tsv");
        let threaded_metrics = std::fs::read_to_string(&threaded_tsv).expect("threaded tsv");
        assert_eq!(serial_metrics, threaded_metrics, "metrics TSV differ between modes");

        // Guard that the corpus actually exercises the counters whose threaded
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
