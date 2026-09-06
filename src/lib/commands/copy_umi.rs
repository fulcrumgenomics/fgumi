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
    add_pg_record, ensure_hd_record, reject_output_collisions,
};
use crate::logging::OperationTimer;
use crate::per_thread_accumulator::PerThreadAccumulator;
use crate::sam::SamTag;
use crate::umi::read_name::normalize_read_name_umi;
use fgumi_bam_io::{ProgressTracker, create_raw_bam_reader_with_opts, create_raw_bam_writer};
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
/// Shared by the serial oracle and the chain finalize hooks, so the two paths
/// reduce the same counters (`pub(crate)` fields for the chain step module).
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

/// Write the single-row `--metrics` TSV from the reduced counters. Shared by the
/// serial oracle and the chain's success-only metrics finalize hook so the two
/// paths cannot drift.
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

/// Emit the overwrite warning (when any) and the `=== Summary ===` block. Shared
/// by the serial oracle and the chain's summary finalize hook so the log content
/// is byte-identical across the two paths.
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
        // Fail fast on a non-ASCII delimiter before any reader opens; each path
        // (serial oracle / chain) re-derives the validated byte itself.
        validate_field_delimiter(self.field_delimiter)?;

        // Reject two writers to one destination, and a write target aliasing input.
        let mut outputs: Vec<(&Path, &str)> = vec![(self.io.output.as_path(), "--output")];
        if let Some(path) = &self.metrics {
            outputs.push((path.as_path(), "--metrics"));
        }
        reject_output_collisions(&outputs)?;
        self.reject_write_aliasing_input()?;

        if self.threading.threads.is_some() {
            return self.execute_chain(command_line);
        }

        // No-`--threads` serial path: the in-process parity oracle.
        let timer = OperationTimer::new("Copying UMIs from read names");
        info!("Starting copy-umi");
        info!("Input: {}", self.io.input.display());
        info!("Output: {}", self.io.output.display());

        // Without this, a scheduler/pipeline flag set here (e.g.
        // `--deadlock-recover`) is silently ignored on the serial path with no
        // diagnostic at all. This is NOT `common::warn_unwired_pipeline_flags`:
        // that helper is chain-specific (its wording says "the chain engine ...",
        // and it deliberately skips `--max-memory`/`--memory-reserve`/
        // `--memory-per-thread` because the chain DOES honor them) — on this
        // serial path there is no chain engine at all, so those queue-memory
        // flags are inert here too and must be warned, with wording that is true
        // for the single-threaded context.
        self.warn_unwired_flags_on_serial_path();

        let totals = self.run_single_threaded(command_line)?;
        warn_and_log_copy_umi_summary(&totals);
        if let Some(path) = &self.metrics {
            write_copy_umi_metrics(path, &totals)?;
        }
        timer.log_completion(totals.total_records);
        Ok(())
    }
}

impl CopyUmi {
    /// Warn about scheduler/pipeline-queue flags that have no effect on THIS
    /// serial (no-`--threads`) path.
    ///
    /// Deliberately NOT `common::warn_unwired_pipeline_flags`: that helper is
    /// chain-specific -- its wording ("the chain engine ...", "the typed-step
    /// pipeline") describes a running chain, and it deliberately skips
    /// `--max-memory`/`--memory-reserve`/`--memory-per-thread` because the chain
    /// DOES honor them (they seed `PipelineConfig::queue_memory_total`). On this
    /// serial path there is no chain and no queue to budget, so reusing that
    /// helper would both assert a false rationale and fail to warn about the
    /// queue-memory flags, which are genuinely inert here. Each flag is warned
    /// only when set away from its default, mirroring how the chain-path helper
    /// gates its own warnings.
    fn warn_unwired_flags_on_serial_path(&self) {
        let requested_scheduler = self.scheduler_opts.strategy();
        if requested_scheduler != crate::unified_pipeline::scheduler::SchedulerStrategy::default() {
            warn!(
                "--scheduler has no effect on the single-threaded path: the serial \
                 read -> copy-umi -> write loop does not use a pluggable scheduler \
                 strategy, so the requested strategy ({requested_scheduler:?}) is \
                 ignored (run with --threads to use the pipeline)"
            );
        }
        if self.scheduler_opts.deadlock_recover_enabled() {
            warn!(
                "--deadlock-recover has no effect on the single-threaded path: there \
                 is no pipeline to deadlock-detect or recover on, only a serial \
                 read -> copy-umi -> write loop (run with --threads to use the \
                 pipeline)"
            );
        }
        let default_queue_memory = QueueMemoryOptions::default();
        if self.queue_memory.max_memory != default_queue_memory.max_memory
            || self.queue_memory.memory_reserve != default_queue_memory.memory_reserve
            || self.queue_memory.memory_per_thread != default_queue_memory.memory_per_thread
        {
            warn!(
                "--max-memory/--memory-reserve/--memory-per-thread have no effect on \
                 the single-threaded path: there are no pipeline queues to budget on \
                 a serial read -> copy-umi -> write loop (run with --threads to use \
                 the pipeline)"
            );
        }
    }

    /// Serial no-`--threads` path: a read → copy-umi → write loop, and the
    /// in-process parity oracle for the chain path. Uses `pipeline_reader_opts()`
    /// so the CRC-check policy matches the chain, and `ensure_hd_record` +
    /// `add_pg_record` so the header matches. A bad/empty UMI (or an existing RX
    /// under `--fail-if-tag-present`) aborts the run via `?`, naming the read name.
    ///
    /// Parity: chain and oracle produce identical output on well-formed BAM, and
    /// BOTH reject a framing-consistent-but-malformed record (`l_read_name` past
    /// the record end) with a clean `InvalidData` error. This serial loop
    /// validates each record via `validate_record_for_decode` before touching
    /// read-name bytes, exactly as the chain path does through `DecodeRecords` —
    /// without that guard the loop would panic in `read_name` on such a record,
    /// a regression from the pre-cutover pipeline, which decoded and thus
    /// validated.
    fn run_single_threaded(&self, command_line: &str) -> Result<CollectedCopyUmiMetrics> {
        let field_delimiter = validate_field_delimiter(self.field_delimiter)?;
        let (mut reader, header) =
            create_raw_bam_reader_with_opts(&self.io.input, 1, self.io.pipeline_reader_opts())?;
        let header = ensure_hd_record(header)?;
        let header = add_pg_record(header, command_line)?;
        let mut writer =
            create_raw_bam_writer(&self.io.output, &header, 1, self.compression.compression_level)?;

        let mut totals = CollectedCopyUmiMetrics::default();
        let progress = ProgressTracker::new("Processed records").with_interval(1_000_000);
        let mut record = RawRecord::new();
        while reader.read_record(&mut record)? != 0 {
            // Validate record framing (record length + l_read_name bound) before
            // touching read-name bytes, exactly as the `--threads` chain path does
            // via `DecodeRecords`. Without this, a framing-consistent-but-malformed
            // record (l_read_name past the record end) would panic here in
            // `read_name` instead of failing with a clean error — a regression from
            // the pre-cutover pipeline, which validated. Both paths now reject a
            // malformed record with the same InvalidData error.
            crate::pipeline::steps::parse::decode::validate_record_for_decode(record.as_ref())?;
            let outcome = copy_umi_into_record(
                &mut record,
                field_delimiter,
                self.reverse_complement_r_umis,
                self.remove_umi,
                self.fail_if_tag_present,
            )?;
            totals.total_records += 1;
            if outcome.overwrote_rx {
                totals.rx_overwritten += 1;
            }
            if outcome.trimmed_name {
                totals.names_trimmed += 1;
            }
            writer.write_raw_record(record.as_ref())?;
            progress.log_if_needed(1);
        }
        progress.log_final();
        // Flush before returning so any write error surfaces here, not on drop.
        writer.finish()?;
        Ok(totals)
    }

    /// Run the copy-umi stage on the declarative chain builder (the `--threads N`
    /// path). The no-`--threads` `CopyUmi::run_single_threaded` loop is the
    /// in-process parity oracle. `add_copy_umi` re-emits the timer/banner/threading
    /// logs and opens its own source, so — like dedup/group — this method only
    /// builds and runs the chain. Copy-umi's `execute()` never emitted
    /// `log_effective_check_crc`, so this does NOT call it either (adding it would
    /// make the `--threads` path log a line the serial path never does).
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
