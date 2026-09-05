//! `copy-umi` — copy the UMI embedded in a BAM read name into the `RX` tag.
//!
//! Ports fgbio's `CopyUmiFromReadName`: the UMI is the last field of the read
//! name (split on `--field-delimiter`, default `:`), normalized to the canonical
//! SAM form (dual UMIs `-`-joined, `r`-prefixed segments reverse-complemented by
//! default), and written to the `RX` tag. It is a non-destructive, streaming
//! BAM→BAM pass: alignments, `SEQ`, and other tags pass through untouched, and
//! only the read name (optionally, with `--remove-umi`) and the `RX` tag change.
//!
//! Unlike `extract` (which parses FASTQ read names during unmapped-BAM creation
//! and uses fgbio's strict ≥8-field rule), `copy-umi` uses fgbio
//! `CopyUmiFromReadName`'s lenient rule: the UMI is simply the last delimited
//! field, with no field-count requirement. The two commands share the
//! normalization step (`crate::umi::read_name::normalize_read_name_umi`) but keep
//! their own field-location logic.

use std::path::{Path, PathBuf};
use std::sync::Arc;
use std::sync::atomic::AtomicU64;

use anyhow::{Context, Result, bail, ensure};
use clap::Parser;
use log::{info, warn};

use crate::commands::command::Command;
use crate::commands::common::{
    BamIoOptions, CompressionOptions, QueueMemoryOptions, SchedulerOptions, ThreadingOptions,
    reject_output_collisions,
};
use crate::per_thread_accumulator::PerThreadAccumulator;
use crate::sam::SamTag;
use crate::umi::read_name::normalize_read_name_umi;
use fgumi_raw_bam::{RawRecord, update_string_tag};

/// Copy the UMI from a BAM read name into the RX tag.
#[derive(Debug, Parser)]
#[command(
    name = "copy-umi",
    about = "\x1b[38;5;166m[UTILITIES]\x1b[0m      \x1b[36mCopy the UMI from the read name into the RX tag\x1b[0m",
    long_about = r#"
Copy the UMI embedded in each read name into the RX SAM tag.

Migrating from fgbio's CopyUmiFromReadName: for a BAM whose reads carry the UMI as
the last field of the read name (e.g. the Illumina 8th field, or a UMI appended by
an upstream tool), this copies that UMI into the RX tag so downstream fgumi tools
(which read UMIs only from RX) can consume it.

It is non-destructive and streams BAM in to BAM out: alignments, sequence, and all
other tags pass through untouched. Only the RX tag is added and, with --remove-umi,
the UMI is trimmed off the read name. This composes in a pipe, e.g. between
`samtools fixmate` and `fgumi sort`.

## UMI location

The UMI is the last `--field-delimiter` (default `:`) field of the read name. There
is no field-count requirement (fgbio CopyUmiFromReadName's lenient rule); a read
name with fewer than two fields, or an empty last field, is an error.

## UMI normalization

The extracted UMI is normalized to the canonical SAM form, matching fgbio:
- multiple UMIs in the field (delimited by `+`, fixed) are joined with `-`;
- segments prefixed with `r` (marking reverse-complemented UMIs, e.g. from BCL
  Convert) are reverse-complemented back to the forward strand by default; pass
  --reverse-complement-r-umis false to instead strip the `r` and keep the stated
  orientation (matching UMI-tools/fgpyo);
- the result is upper-cased and rejected if it contains any character outside
  `ACGTN-`; a field that normalizes to nothing (e.g. a bare `r`) is an error.

An error on any record fails the run (fail-fast), naming the read name.

## Performance

The work is dominated by BGZF (de)compression, not the UMI copy itself, and
output compression is the largest single cost. In a pipe to another tool, pass a
low `--compression-level` (e.g. 1, or 0 for uncompressed BAM straight into
`fgumi sort`) to roughly halve CPU; use `--threads` to parallelize.
"#
)]
#[command(verbatim_doc_comment)]
#[allow(clippy::struct_excessive_bools)]
pub struct CopyUmi {
    /// Input/output BAM paths and reader options.
    #[command(flatten)]
    pub io: BamIoOptions,

    /// Read-name field delimiter; the UMI is the last field. Must be ASCII.
    #[arg(long = "field-delimiter", default_value_t = ':')]
    pub field_delimiter: char,

    /// Remove the UMI (and its delimiter) from the read name after copying it.
    #[arg(long = "remove-umi", value_name = "true|false", default_value = "false", num_args = 0..=1, default_missing_value = "true", action = clap::ArgAction::Set, value_parser = clap::builder::BoolishValueParser::new(), hide_possible_values = true)]
    pub remove_umi: bool,

    /// Reverse-complement `r`-prefixed UMI segments (fgbio default). Set to
    /// `false` to strip the `r` marker without reverse-complementing, matching
    /// UMI-tools/fgpyo.
    #[arg(long = "reverse-complement-r-umis", value_name = "true|false", default_value = "true", num_args = 0..=1, default_missing_value = "true", action = clap::ArgAction::Set, value_parser = clap::builder::BoolishValueParser::new(), hide_possible_values = true)]
    pub reverse_complement_r_umis: bool,

    /// Error if a record already carries an RX tag, instead of overwriting it.
    #[arg(long = "fail-if-tag-present", value_name = "true|false", default_value = "false", num_args = 0..=1, default_missing_value = "true", action = clap::ArgAction::Set, value_parser = clap::builder::BoolishValueParser::new(), hide_possible_values = true)]
    pub fail_if_tag_present: bool,

    /// Optional TSV file for run metrics.
    #[arg(short = 'M', long = "metrics")]
    pub metrics: Option<PathBuf>,

    /// Threading options for the parallel pipeline.
    #[command(flatten)]
    pub threading: ThreadingOptions,

    /// Compression options for the output BAM.
    #[command(flatten)]
    pub compression: CompressionOptions,

    /// Scheduler and pipeline statistics options.
    #[command(flatten)]
    pub scheduler_opts: SchedulerOptions,

    /// Queue memory options.
    #[command(flatten)]
    pub queue_memory: QueueMemoryOptions,
}

/// Per-record outcome flags, aggregated into [`CollectedCopyUmiMetrics`].
///
/// `pub(crate)` so the chain step module can fold the outcome into its batch
/// counters.
#[derive(Debug, Clone, Copy, Default)]
pub(crate) struct RecordOutcome {
    /// An existing RX tag was overwritten with the newly extracted UMI.
    pub(crate) overwrote_rx: bool,
    /// The UMI was trimmed off the read name (`--remove-umi`).
    pub(crate) trimmed_name: bool,
}

/// Aggregated run counters, one instance per pipeline slot.
///
/// Reduced by the chain finalize hooks (`pub(crate)` fields for the chain step
/// module).
#[derive(Debug, Clone, Default)]
pub(crate) struct CollectedCopyUmiMetrics {
    /// Records processed (all successful; a failure aborts the run).
    pub(crate) total_records: u64,
    /// Records where an existing RX tag was overwritten.
    pub(crate) rx_overwritten: u64,
    /// Records whose read name had the UMI trimmed off.
    pub(crate) names_trimmed: u64,
}

/// One metrics row summarizing the run, written to the optional `--metrics` TSV.
#[derive(Debug, Clone, Default, serde::Serialize, serde::Deserialize)]
pub struct CopyUmiMetric {
    /// Total records processed.
    pub total_records: u64,
    /// Records where the RX tag was written (equals `total_records`).
    pub rx_written: u64,
    /// Records where a pre-existing RX tag was overwritten.
    pub rx_overwritten: u64,
    /// Records whose read name had the UMI trimmed off (`--remove-umi`).
    pub names_trimmed: u64,
}

impl fgumi_metrics::Metric for CopyUmiMetric {
    fn metric_name() -> &'static str {
        "copy_umi"
    }
}

/// Copy the read-name UMI into the RX tag of a single record, in place.
///
/// Locates the UMI as the last `field_delimiter` field of the read name,
/// normalizes it via [`normalize_read_name_umi`], writes it to `RX` (overwriting
/// or, with `fail_if_tag_present`, erroring on an existing tag), and optionally
/// trims it off the read name.
///
/// # Errors
/// Fails if the read name has no delimiter, has an empty last field, yields an
/// invalid UMI, or already carries an RX tag when `fail_if_tag_present` is set.
pub(crate) fn copy_umi_into_record(
    record: &mut RawRecord,
    field_delimiter: u8,
    reverse_complement_prefixed: bool,
    remove_umi: bool,
    fail_if_tag_present: bool,
) -> Result<RecordOutcome> {
    // Locate + normalize the UMI while only borrowing the record immutably, then
    // drop the borrow before mutating tags/name.
    let (delim_idx, umi) = {
        let name = record.read_name();
        let Some(idx) = name.iter().rposition(|&b| b == field_delimiter) else {
            bail!(
                "read name '{}' has no '{}'-delimited UMI field",
                String::from_utf8_lossy(name),
                field_delimiter as char
            );
        };
        let last = &name[idx + 1..];
        // fgbio would write an empty RX here; we fail instead — an empty UMI only
        // defers the failure to `group` (mixed UMI lengths) with a worse message.
        if last.is_empty() {
            bail!(
                "read name '{}' has an empty UMI field after the final '{}'",
                String::from_utf8_lossy(name),
                field_delimiter as char
            );
        }
        // `normalize_read_name_umi` rejects a field that normalizes to empty (e.g. a
        // bare `r`), so an empty RX can never reach the output.
        let umi =
            normalize_read_name_umi(last, reverse_complement_prefixed).with_context(|| {
                format!("extracting UMI from read name '{}'", String::from_utf8_lossy(name))
            })?;
        (idx, umi)
    };

    // Write RX, overwriting (or erroring on) a pre-existing tag. `update_string_tag`
    // replaces an existing RX in place and appends the SAM-spec NUL itself; the
    // `contains` probe is still needed for the overwrite counter and fail policy.
    let overwrote_rx = record.tags().contains(SamTag::RX);
    ensure!(
        !(overwrote_rx && fail_if_tag_present),
        "record '{}' already has an RX tag (drop --fail-if-tag-present to overwrite)",
        String::from_utf8_lossy(record.read_name())
    );
    update_string_tag(record.as_mut_vec(), SamTag::RX, &umi);

    // Trim the UMI (and its leading delimiter) off the read name, matching fgbio's
    // `name.substring(0, lastIndexOf(delim))`.
    let trimmed_name = if remove_umi {
        // `delim_idx == 0` (a leading-delimiter name like `:ACGT`) would trim to an
        // empty QNAME, which no downstream reader accepts; reject it instead.
        ensure!(
            delim_idx > 0,
            "removing the UMI from read name '{}' would leave an empty read name",
            String::from_utf8_lossy(record.read_name())
        );
        let new_name = record.read_name()[..delim_idx].to_vec();
        record.set_read_name(&new_name);
        true
    } else {
        false
    };

    Ok(RecordOutcome { overwrote_rx, trimmed_name })
}

impl CopyUmi {
    /// Reject a write target that resolves to the same file as the input, which
    /// would clobber the BAM being read (mirrors `retag`/`merge`).
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

/// Validate that the read-name field delimiter is a single ASCII byte.
///
/// The sole delimiter validator, called by both `execute()`'s reader-free
/// pre-flight and `CopyUmiOptions::process_captures` (so a chain built directly,
/// bypassing `execute()`, still validates) — no divergence between the paths.
pub(crate) fn validate_field_delimiter(c: char) -> Result<u8> {
    u8::try_from(c)
        .ok()
        .filter(u8::is_ascii)
        .with_context(|| format!("--field-delimiter must be a single ASCII character, got '{c}'"))
}

/// Write the single-row `--metrics` TSV from the reduced counters. Called by
/// the chain's success-only metrics finalize hook.
pub(crate) fn write_copy_umi_metrics(path: &Path, totals: &CollectedCopyUmiMetrics) -> Result<()> {
    let row = CopyUmiMetric {
        total_records: totals.total_records,
        rx_written: totals.total_records,
        rx_overwritten: totals.rx_overwritten,
        names_trimmed: totals.names_trimmed,
    };
    fgumi_metrics::write_metrics(path, std::slice::from_ref(&row), "copy_umi")?;
    info!("Wrote metrics to: {}", path.display());
    Ok(())
}

/// Emit the overwrite warning (when any) and the `=== Summary ===` block.
/// Called by the chain's summary finalize hook.
pub(crate) fn warn_and_log_copy_umi_summary(totals: &CollectedCopyUmiMetrics) {
    if totals.rx_overwritten > 0 {
        warn!("overwrote a pre-existing RX tag on {} record(s)", totals.rx_overwritten);
    }
    info!("=== Summary ===");
    info!("Records processed: {}", totals.total_records);
    info!("RX tags written: {}", totals.total_records);
    info!("RX tags overwritten: {}", totals.rx_overwritten);
    info!("Read names trimmed: {}", totals.names_trimmed);
}

/// Copy-umi tuning projected off the parsed CLI flags, independent of how they
/// were supplied. Mirrors `FilterOptions`: the chain builder reads only these
/// values, so they are projected into a plain struct rather than the bag holding
/// the whole `CopyUmi` command.
#[derive(Debug, Clone)]
#[allow(clippy::struct_excessive_bools)]
pub struct CopyUmiOptions {
    /// Read-name field delimiter (the UMI is the last field). Validated to a
    /// single ASCII byte in `CopyUmiOptions::process_captures`.
    pub field_delimiter: char,
    /// Remove the UMI (and its delimiter) from the read name after copying it.
    pub remove_umi: bool,
    /// Reverse-complement `r`-prefixed UMI segments (fgbio default).
    pub reverse_complement_r_umis: bool,
    /// Error if a record already carries an RX tag, instead of overwriting it.
    pub fail_if_tag_present: bool,
    /// Optional TSV file for run metrics.
    pub metrics: Option<PathBuf>,
}

impl CopyUmi {
    /// Project the parsed CLI flags into [`CopyUmiOptions`].
    #[must_use]
    pub fn to_copy_umi_options(&self) -> CopyUmiOptions {
        CopyUmiOptions {
            field_delimiter: self.field_delimiter,
            remove_umi: self.remove_umi,
            reverse_complement_r_umis: self.reverse_complement_r_umis,
            fail_if_tag_present: self.fail_if_tag_present,
            metrics: self.metrics.clone(),
        }
    }
}

/// Shared state for a copy-umi chain run, built by
/// [`CopyUmiOptions::setup_pipeline`].
pub(crate) struct CopyUmiPipelineSetup {
    /// Per-thread run counters, reduced by the finalize hooks.
    pub(crate) collected_metrics: Arc<PerThreadAccumulator<CollectedCopyUmiMetrics>>,
    /// Cross-thread processed-record counter driving the heartbeat.
    pub(crate) progress_counter: Arc<AtomicU64>,
}

/// Captures the chain process closure needs, extracted from
/// [`CopyUmiOptions`] and [`CopyUmiPipelineSetup`]. Copy-umi's transform ignores
/// the header, so — unlike `FilterProcessCaptures` — this carries no header.
pub(crate) struct CopyUmiProcessCaptures {
    /// The validated single-byte read-name field delimiter.
    pub(crate) field_delimiter: u8,
    /// Reverse-complement `r`-prefixed UMI segments (fgbio default).
    pub(crate) reverse_complement_prefixed: bool,
    /// Trim the UMI off the read name after copying it.
    pub(crate) remove_umi: bool,
    /// Error on a pre-existing RX tag instead of overwriting.
    pub(crate) fail_if_tag_present: bool,
    /// Cross-thread processed-record counter for the heartbeat.
    pub(crate) progress: Arc<AtomicU64>,
}

impl CopyUmiOptions {
    /// Build the per-thread metrics accumulators and progress counter shared by
    /// the chain builder's process step and finalize hooks. The
    /// `BamPipelineConfig` is not built here — the chain builder derives it from
    /// its `ChainSpec`.
    pub(crate) fn setup_pipeline(&self, num_threads: usize) -> CopyUmiPipelineSetup {
        CopyUmiPipelineSetup {
            collected_metrics: PerThreadAccumulator::<CollectedCopyUmiMetrics>::new(num_threads),
            progress_counter: Arc::new(AtomicU64::new(0)),
        }
    }

    /// Build the process-closure captures. Validates the field delimiter here too
    /// (not only in `execute()`'s pre-flight), so a chain built directly still
    /// rejects a non-ASCII delimiter.
    pub(crate) fn process_captures(
        &self,
        setup: &CopyUmiPipelineSetup,
    ) -> Result<CopyUmiProcessCaptures> {
        Ok(CopyUmiProcessCaptures {
            field_delimiter: validate_field_delimiter(self.field_delimiter)?,
            reverse_complement_prefixed: self.reverse_complement_r_umis,
            remove_umi: self.remove_umi,
            fail_if_tag_present: self.fail_if_tag_present,
            progress: Arc::clone(&setup.progress_counter),
        })
    }
}

impl Command for CopyUmi {
    fn execute(&self, command_line: &str) -> Result<()> {
        self.io.validate()?;
        // Fail fast on a non-ASCII delimiter before any reader opens; the chain's
        // own `CopyUmiOptions::process_captures` re-derives the validated byte too.
        validate_field_delimiter(self.field_delimiter)?;

        // Reject two writers to one destination, and a write target aliasing input.
        let mut outputs: Vec<(&Path, &str)> = vec![(self.io.output.as_path(), "--output")];
        if let Some(path) = &self.metrics {
            outputs.push((path.as_path(), "--metrics"));
        }
        reject_output_collisions(&outputs)?;
        self.reject_write_aliasing_input()?;

        self.execute_chain(command_line)
    }
}

impl CopyUmi {
    /// Run the copy-umi stage on the declarative chain builder. `add_copy_umi`
    /// re-emits the timer/banner/threading logs and opens its own source, so —
    /// like dedup/group — this method only builds and runs the chain. Copy-umi's
    /// `execute()` never emitted `log_effective_check_crc`, and neither does
    /// `add_copy_umi`.
    fn execute_chain(&self, command_line: &str) -> Result<()> {
        use crate::pipeline::chains::{
            ChainSpec, SingleStageContext, Stage, StageOptionsBag, build_for,
        };

        let stage_opts =
            StageOptionsBag { copy_umi: Some(self.to_copy_umi_options()), ..Default::default() };
        let ctx = SingleStageContext {
            io: &self.io,
            threading: &self.threading,
            compression: &self.compression,
            scheduler: &self.scheduler_opts,
            queue_memory: &self.queue_memory,
            command_line,
        };
        let spec = ChainSpec::single_stage(Stage::CopyUmi, stage_opts, &ctx);
        build_for(spec)?.run()
    }
}
