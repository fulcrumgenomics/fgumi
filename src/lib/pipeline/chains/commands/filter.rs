//! Chain builder for `Stage::Filter`.
//!
//! Phase 2 (T2.14) held the full ~500-LOC chain construction here.
//! Phase 3 (T3a.4) lifts that logic into
//! [`crate::pipeline::chains::builder::ChainBuilder`]; this module
//! now holds the filter-specific types and step factories that the builder
//! imports: `FilterFinalizeHook` and the step-factory functions used by
//! `ChainBuilder::add_filter`.
//!
//! ## Four chain shapes
//!
//! Filter dispatches on `(filter_by_template, track_rejects)`:
//!
//! - `(false, false)`: `ProcessOrdered<DecodedRecordBatch, DecompressedBlock>`
//! - `(false, true)`: `Process2Ordered<DecodedRecordBatch, DecompressedBlock, DecompressedBlock>`
//! - `(true, false)`: `GroupByQueryname` + `ProcessOrdered<BamTemplateBatch, DecompressedBlock>`
//! - `(true, true)`: `GroupByQueryname` + `Process2Ordered<BamTemplateBatch, DecompressedBlock, DecompressedBlock>`
//!
//! For shapes with rejects (`track_rejects = true`), branch 0 is kept (flows
//! to `add_sink`) and branch 1 is the rejects output (wired to its own
//! compress+write inside `add_filter`).

use std::io;
use std::sync::Arc;
use std::sync::atomic::Ordering;

use anyhow::Result;
use fgumi_raw_bam::RawRecord;
use log::{info, warn};

use crate::commands::filter::{CollectedFilterMetrics, Filter, FilterProcessCaptures};
use crate::consensus_filter::{is_primary_read, retained_primary_masked_bases};
use crate::logging::OperationTimer;
use crate::per_thread_accumulator::PerThreadAccumulator;
use crate::pipeline::chains::FinalizeHook;
use crate::pipeline::chains::commands::framed_len;
use crate::pipeline::steps::process::{
    Process2Ordered, ProcessOrdered, process_ordered, process2_ordered,
};
use crate::pipeline::steps::types::{
    BamTemplateBatch, DecodedRecordBatch, DecompressedBlock, RecordBatch, fit_for_queue,
};

// ─────────────────────────────────────────────────────────────────────────────
// FilterFinalizeHook
// ─────────────────────────────────────────────────────────────────────────────

/// Post-pipeline finalize hook for filter. Reduces per-thread metrics and logs
/// the summary banner. Registered on the **always-run** `finalize` list so the
/// summary is reported even after a failure; the stats *file* is written by the
/// success-only [`FilterStatsFinalizeHook`] instead, so a partial run never
/// publishes counts.
///
/// `pub(crate)` so `ChainBuilder` can construct and register it in
/// `add_filter`.
pub(crate) struct FilterFinalizeHook {
    pub(crate) accumulators: Arc<PerThreadAccumulator<CollectedFilterMetrics>>,
    pub(crate) has_rejects: bool,
    pub(crate) timer: OperationTimer,
    /// Records the methylation filters could not evaluate.
    pub(crate) methylation_skips: Arc<crate::commands::filter::MethylationFilterSkips>,
}

impl FinalizeHook for FilterFinalizeHook {
    fn finalize(self: Box<Self>) -> Result<()> {
        let FilterFinalizeHook { accumulators, has_rejects, timer, methylation_skips } = *self;

        let mut total_reads = 0u64;
        let mut passed_reads = 0u64;
        let mut failed_reads = 0u64;
        let mut total_bases_masked = 0u64;
        let mut evaluated_records = 0u64;
        for slot in accumulators.slots() {
            let m = slot.lock();
            total_reads += m.total_records;
            passed_reads += m.passed_records;
            failed_reads += m.failed_records;
            total_bases_masked += m.total_bases_masked;
            evaluated_records += m.evaluated_records;
        }

        info!("Processed {total_reads} reads; kept {passed_reads} and rejected {failed_reads}");
        if has_rejects && failed_reads > 0 {
            info!("Wrote {failed_reads} rejected records to rejects file");
        }
        info!("Total bases masked: {total_bases_masked}");
        for warning in methylation_skip_warnings(evaluated_records, &methylation_skips) {
            warn!("{warning}");
        }

        timer.log_completion(total_reads);

        Ok(())
    }
}

/// The end-of-run warnings for records the methylation filters could not evaluate, or whose
/// modification tags they removed. `evaluated` is the number of records the filter evaluated
/// ([`CollectedFilterMetrics::evaluated_records`]), which the every-record warnings are measured
/// against.
fn methylation_skip_warnings(
    evaluated: u64,
    skips: &crate::commands::filter::MethylationFilterSkips,
) -> Vec<String> {
    let mut warnings = Vec::new();
    let unaligned = skips.unaligned.load(Ordering::Relaxed);
    let no_counts = skips.no_counts.load(Ordering::Relaxed);
    let length_mismatch = skips.length_mismatch.load(Ordering::Relaxed);
    // The skip causes are exclusive per record: together they may cover every record without
    // any one of them doing so, and then the methylation filters also checked nothing.
    let skipped = unaligned + no_counts + length_mismatch;
    let single_cause = [unaligned, no_counts, length_mismatch].contains(&evaluated);
    if skipped > 0 && skipped == evaluated && !single_cause {
        warnings.push(format!(
            "none of the {evaluated} records were checked by the methylation filters (see the \
             reasons below)"
        ));
    }
    // A skip that hit every record means the methylation filters checked nothing at all.
    let mut skip_warning = |count: u64, every: &str, some: String| {
        if count > 0 && count == evaluated {
            warnings.push(format!(
                "none of the {count} records were checked by the methylation filters: {every}"
            ));
        } else if count > 0 {
            warnings.push(some);
        }
    };
    skip_warning(
        unaligned,
        "single-strand records need records aligned to the reference, and every record was \
         unmapped; filter after alignment",
        format!(
            "{unaligned} unmapped single-strand records were not checked by the methylation \
             filters"
        ),
    );
    skip_warning(
        no_counts,
        "no record carries cu/ct methylation counts; call consensus with --methylation-mode",
        format!(
            "{no_counts} records without cu/ct methylation counts were not checked by the \
             methylation filters"
        ),
    );
    skip_warning(
        length_mismatch,
        "no record's cu/ct match its SEQ in length (SEQ hard-clipped after consensus calling?)",
        format!(
            "{length_mismatch} records whose cu/ct do not match SEQ in length (SEQ hard-clipped \
             after consensus calling?) were not checked by the methylation filters"
        ),
    );
    let simplex_agreement = skips.simplex_agreement.load(Ordering::Relaxed);
    if simplex_agreement > 0 {
        warnings.push(format!(
            "--require-strand-methylation-agreement applies to duplex consensus records only; it \
             was not applied to {simplex_agreement} single-strand records"
        ));
    }
    let unreversed = skips.unreversed_reverse_strand.load(Ordering::Relaxed);
    if unreversed > 0 {
        warnings.push(format!(
            "the methylation filters read cu/ct by SEQ position in reference orientation, but \
             --reverse-per-base-tags is not set: unless cu/ct were already reversed (for example \
             by zipper --tags-to-reverse Consensus), they were evaluated at the wrong positions on \
             {unreversed} reverse-mapped records"
        ));
    }
    let dropped_calls = skips.dropped_calls.load(Ordering::Relaxed);
    if dropped_calls > 0 {
        let records = skips.records_with_dropped_calls.load(Ordering::Relaxed);
        warnings.push(format!(
            "the methylation filters dropped {dropped_calls} methylation calls from MM/ML on \
             {records} duplex records (SEQ unchanged)"
        ));
    }
    let tags_removed = skips.tags_removed.load(Ordering::Relaxed);
    if tags_removed > 0 {
        warnings.push(format!(
            "MM/ML/am/bm were removed from {tags_removed} records whose modification tags could \
             not be updated to match SEQ"
        ));
    }
    warnings
}

/// Success-only finalize hook that writes the `--filter::stats` file. Registered
/// on `finalize_on_success` (not the always-run `finalize`) so a failed or
/// partial run never publishes filter counts. Shares the accumulators with
/// [`FilterFinalizeHook`] via `Arc`.
pub(crate) struct FilterStatsFinalizeHook {
    pub(crate) accumulators: Arc<PerThreadAccumulator<CollectedFilterMetrics>>,
    pub(crate) stats_path: std::path::PathBuf,
}

impl FinalizeHook for FilterStatsFinalizeHook {
    fn finalize(self: Box<Self>) -> Result<()> {
        let FilterStatsFinalizeHook { accumulators, stats_path } = *self;

        let mut total_reads = 0u64;
        let mut passed_reads = 0u64;
        let mut failed_reads = 0u64;
        for slot in accumulators.slots() {
            let m = slot.lock();
            total_reads += m.total_records;
            passed_reads += m.passed_records;
            failed_reads += m.failed_records;
        }

        write_filter_stats(&stats_path, total_reads, passed_reads, failed_reads)
    }
}

// ─────────────────────────────────────────────────────────────────────────────
// Shared helper
// ─────────────────────────────────────────────────────────────────────────────

/// Write filtering statistics as a one-row fgbio-`Metric` TSV
/// (`total_reads<TAB>passed_reads<TAB>failed_reads<TAB>pass_rate` + one data row).
pub(crate) fn write_filter_stats(
    path: &std::path::Path,
    total: u64,
    passed: u64,
    failed: u64,
) -> Result<()> {
    use crate::metrics::{FilterStatsMetrics, write_metrics_auto};

    let metric = FilterStatsMetrics::from_counts(total, passed, failed);
    write_metrics_auto(path, std::slice::from_ref(&metric))
}

/// Thin wrapper that calls `Filter::process_record_raw` with captures,
/// so the four pipeline-mode closures share identical call sites.
pub(crate) fn process_record_raw_call(
    record: &mut fgumi_raw_bam::RawRecord,
    captures: &FilterProcessCaptures,
) -> anyhow::Result<(u64, bool)> {
    Filter::process_record_raw(
        record,
        &captures.config,
        captures.reference.as_deref(),
        &captures.header,
        captures.should_reverse_tags,
        captures.min_base_quality,
        captures.require_single_strand_agreement,
        captures.min_mean_base_quality,
        captures.max_no_call_fraction,
        captures.methylation_depth_thresholds.as_ref(),
        captures.require_strand_methylation_agreement,
        captures.min_conversion_fraction,
        captures.methylation_mode,
        &captures.ref_names,
        &captures.methylation_skips,
    )
}

// ─────────────────────────────────────────────────────────────────────────────
// Step factories — extracted from build_filter_chain (T3a.4).
//
// Each factory receives its captured state as plain arguments and returns the
// concrete step type (via the constructor's return type). Returning
// `impl Step<...>` is blocked here because the closure types embed in
// opaque-return position and cannot name themselves, so we return the concrete
// `Process*` structs directly (the same pattern used by the dedup factories
// in `chains::commands::dedup`).
// ─────────────────────────────────────────────────────────────────────────────

/// Emit the per-batch progress milestone and fold this batch's counts into the
/// per-thread accumulator. Shared by all four `build_filter_step_*` factories so
/// the metrics contract lives in exactly one place (four copies of it is the
/// sibling-divergence pattern that has shipped bugs in this module).
/// `evaluated_count` is the records the filter evaluated, which only the template
/// steps can leave below `total_records`.
fn record_batch_metrics(
    captures: &FilterProcessCaptures,
    accumulators: &PerThreadAccumulator<CollectedFilterMetrics>,
    total_records: u64,
    passed_count: u64,
    bases_masked_total: u64,
    evaluated_count: u64,
) {
    let prev = captures.progress.fetch_add(total_records, Ordering::Relaxed);
    if (prev + total_records) / 1_000_000 > prev / 1_000_000 {
        info!("Processed {} records", prev + total_records);
    }
    accumulators.with_slot(|m| {
        m.total_records += total_records;
        m.passed_records += passed_count;
        m.failed_records += total_records - passed_count;
        m.total_bases_masked += bases_masked_total;
        m.evaluated_records += evaluated_count;
    });
}

/// Apply the single-record filter transform to one record and route it.
///
/// Runs [`process_record_raw_call`], writes the record framed to `kept` if it
/// passes, else to `rejected` (or drops it when `rejected` is `None`), and
/// returns `(pass, masked_contribution)` for the caller's batch tallies
/// (`masked_contribution` is the record's fgbio "Total bases masked" share — 0
/// unless it is a retained primary read).
///
/// This is the one per-record body shared by all four single-record filter step
/// builders (decoded/raw × no-rejects/with-rejects), so the transform, the
/// masked-bases accounting, and the keep/reject routing cannot drift between the
/// owned `DecodedRecordBatch` path and the borrowed `RecordBatch` fast path.
fn filter_one_record(
    record: &mut RawRecord,
    captures: &FilterProcessCaptures,
    kept: &mut Vec<u8>,
    rejected: Option<&mut Vec<u8>>,
) -> io::Result<(bool, u64)> {
    let (bases_masked, pass) =
        process_record_raw_call(record, captures).map_err(io::Error::other)?;
    // Match fgbio's "Total bases masked": count masked bases only in a retained
    // primary read (0 for a rejected / secondary / supplementary read).
    let masked_contribution =
        retained_primary_masked_bases(std::slice::from_ref(&*record), &[bases_masked], pass);
    if pass {
        fgumi_raw_bam::write_framed_record(kept, record.as_ref())?;
    } else if let Some(rejected) = rejected {
        fgumi_raw_bam::write_framed_record(rejected, record.as_ref())?;
    }
    Ok((pass, masked_contribution))
}

/// The result of filtering one template: whether each record is kept, the template's share of
/// fgbio's "Total bases masked" (the masked bases of its primary reads when it is kept, else 0),
/// and how many of its records were evaluated.
#[derive(Debug, PartialEq, Eq)]
struct TemplateFilterOutcome {
    keep: Vec<bool>,
    masked_bases: u64,
    evaluated: u64,
}

/// Apply the filter to one template the way fgbio does (`FilterConsensusReads.scala:189-219`,
/// fgbio `origin/main` e51a661).
///
/// The primary reads ([`is_primary_read`]) are evaluated in template order (R1 before R2, as
/// [`Template`](crate::template::Template) orders them) and evaluation stops at the first that
/// fails; the secondary and supplementary records are evaluated only when every primary read
/// passed, and each is then kept on its own result.
///
/// A record left unevaluated is rejected as it was read, so a rejected template never errors on
/// a later record that lacks the consensus tags. It gets only [`Filter::ensure_base_qualities`],
/// the check every record gets (a record with bases but no base qualities is an error whether or
/// not it is evaluated), and [`Filter::orient_per_base_tags`], the change every record gets, as
/// fgbio reverses every read of the template before filtering (`:195-200`). Leaving it
/// unevaluated also skips the rest of [`process_record_raw_call`]: an unevaluated record is not
/// checked against the "--ref is required for mapped reads" guard, and is neither counted in the
/// methylation-filter skip warnings nor in the records they are out of
/// ([`TemplateFilterOutcome::evaluated`]).
///
/// # Errors
///
/// A template with no primary read (only secondary or supplementary records, for example after
/// its primary reads were removed upstream) is an error naming the template, as in fgbio, which
/// throws `"<name> had no R1."` (`:190`). So is any error evaluating a record.
///
/// Shared by both template step builders so the short-circuit cannot drift between them.
fn filter_template_records(
    records: &mut [RawRecord],
    captures: &FilterProcessCaptures,
) -> io::Result<TemplateFilterOutcome> {
    if !records.iter().any(is_primary_read) {
        let name = records
            .first()
            .map(|r| String::from_utf8_lossy(fgumi_raw_bam::RawRecordView::new(r).read_name()))
            .unwrap_or_default();
        return Err(io::Error::other(format!(
            "template {name} has no primary read (only secondary or supplementary records); \
             filter needs each template's primary reads, as fgbio does"
        )));
    }

    let mut keep = vec![false; records.len()];
    let mut masked_bases: u64 = 0;
    let mut evaluated: u64 = 0;
    // The index of the primary read that failed, if one did; the primary reads up to and
    // including it are the only records evaluated.
    let mut failed_at: Option<usize> = None;

    for (idx, record) in records.iter_mut().enumerate().filter(|(_, r)| is_primary_read(r)) {
        evaluated += 1;
        let (masked, pass) = process_record_raw_call(record, captures).map_err(io::Error::other)?;
        if !pass {
            failed_at = Some(idx);
            break;
        }
        keep[idx] = true;
        // Only primary reads reach this tally, so it is `retained_primary_masked_bases`'s rule.
        masked_bases += masked;
    }
    if let Some(failed) = failed_at {
        for (idx, record) in records.iter_mut().enumerate() {
            let was_evaluated = is_primary_read(record) && idx <= failed;
            if !was_evaluated {
                Filter::ensure_base_qualities(record).map_err(io::Error::other)?;
                Filter::orient_per_base_tags(record, captures.should_reverse_tags)
                    .map_err(io::Error::other)?;
            }
        }
        keep.fill(false);
        return Ok(TemplateFilterOutcome { keep, masked_bases: 0, evaluated });
    }

    for (idx, record) in records.iter_mut().enumerate().filter(|(_, r)| !is_primary_read(r)) {
        evaluated += 1;
        let (_, pass) = process_record_raw_call(record, captures).map_err(io::Error::other)?;
        keep[idx] = pass;
    }
    Ok(TemplateFilterOutcome { keep, masked_bases, evaluated })
}

/// Kept and (when `rejects`) rejected output buffers for a batch of `framed`
/// input bytes (see [`framed_len`]).
///
/// Filter emits each input record at most once, framed as `4 + len`, so
/// `framed` bounds the kept output. Reserving it up front replaces the `Vec`
/// doubling (realloc + copy) of serializing into an empty buffer; the
/// reservation is capped at `ceiling`, the step's per-queue byte budget
/// (config), and [`filter_output_blocks`] shrinks a buffer that a heavy filter
/// left mostly empty before it is queued (queues are charged by capacity).
/// Rejects are usually a small fraction of a batch, so their buffer is not
/// reserved: it grows from empty, which keeps its capacity within twice what it
/// holds.
fn filter_output_buffers(
    framed: usize,
    ceiling: usize,
    rejects: bool,
) -> (Vec<u8>, Option<Vec<u8>>) {
    (Vec::with_capacity(framed.min(ceiling)), rejects.then(Vec::new))
}

/// Wrap a batch's kept and rejected bytes as `DecompressedBlock`s, each fitted
/// for its queue (capacity at most `2 * len + 64 KiB`, see [`fit_for_queue`]).
fn filter_output_blocks(
    batch_serial: u64,
    mut kept: Vec<u8>,
    rejected: Option<Vec<u8>>,
) -> (DecompressedBlock, Option<DecompressedBlock>) {
    fit_for_queue(&mut kept);
    (
        DecompressedBlock { batch_serial, bytes: kept },
        rejected.map(|mut bytes| {
            fit_for_queue(&mut bytes);
            DecompressedBlock { batch_serial, bytes }
        }),
    )
}

/// The rejects branch of a with-rejects step's output. A with-rejects step
/// always asks its batch body for rejects, so `rejected` is `Some`; a `None`
/// here is a wiring bug, and an empty rejects block would silently write an
/// empty rejects file, so it panics instead.
///
/// # Panics
///
/// Panics if `rejected` is `None`.
fn rejects_block(
    kept: &DecompressedBlock,
    rejected: Option<DecompressedBlock>,
) -> DecompressedBlock {
    rejected.unwrap_or_else(|| {
        unreachable!(
            "with-rejects filter step got no rejects block for batch {}",
            kept.batch_serial
        )
    })
}

/// A step's per-queue byte budget as a reservation ceiling.
fn reserve_ceiling(limit_bytes: u64) -> usize {
    usize::try_from(limit_bytes).unwrap_or(usize::MAX)
}

/// Per-batch body of the single-read filter steps on the decoded path
/// (`DecodedRecordBatch`). Returns the kept records and, when `rejects` is set,
/// the rejected ones, serialized framed into the buffers
/// [`filter_output_buffers`] sizes (kept reserved once, capped at `ceiling`).
fn filter_single_batch(
    item: DecodedRecordBatch,
    captures: &FilterProcessCaptures,
    accumulators: &PerThreadAccumulator<CollectedFilterMetrics>,
    ceiling: usize,
    rejects: bool,
) -> io::Result<(DecompressedBlock, Option<DecompressedBlock>)> {
    let batch_serial = item.batch_serial();
    let framed = framed_len(item.records().iter().map(fgumi_bam_io::DecodedRecord::raw_bytes));
    let records = item.into_records();
    let records_count = records.len() as u64;
    let (mut kept_bytes, mut rejected_bytes) = filter_output_buffers(framed, ceiling, rejects);
    let mut passed_count: u64 = 0;
    let mut bases_masked_total: u64 = 0;

    for decoded in records {
        let mut record = decoded.into_raw_bytes();
        let (pass, masked) =
            filter_one_record(&mut record, captures, &mut kept_bytes, rejected_bytes.as_mut())?;
        if pass {
            passed_count += 1;
        }
        bases_masked_total += masked;
    }

    record_batch_metrics(
        captures,
        accumulators,
        records_count,
        passed_count,
        bases_masked_total,
        records_count, // every record is evaluated
    );

    Ok(filter_output_blocks(batch_serial, kept_bytes, rejected_bytes))
}

/// Per-batch body of the single-read filter steps on the decode-free fast path
/// (`RecordBatch`, the default `filter` chain). Same contract as
/// [`filter_single_batch`]; iterates borrowed record byte-ranges with one
/// reused scratch `RawRecord`.
fn filter_single_raw_batch(
    item: &RecordBatch,
    captures: &FilterProcessCaptures,
    accumulators: &PerThreadAccumulator<CollectedFilterMetrics>,
    ceiling: usize,
    rejects: bool,
) -> io::Result<(DecompressedBlock, Option<DecompressedBlock>)> {
    let batch_serial = item.batch_serial();
    let records_count = item.len() as u64;
    let framed = framed_len(item.iter_record_bytes());
    let (mut kept_bytes, mut rejected_bytes) = filter_output_buffers(framed, ceiling, rejects);
    let mut passed_count: u64 = 0;
    let mut bases_masked_total: u64 = 0;
    // One scratch buffer, reused across every record in this batch.
    let mut scratch = RawRecord::new();

    crate::pipeline::chains::commands::for_each_raw_record(item, &mut scratch, |record| {
        let (pass, masked) =
            filter_one_record(record, captures, &mut kept_bytes, rejected_bytes.as_mut())?;
        if pass {
            passed_count += 1;
        }
        bases_masked_total += masked;
        Ok(())
    })?;

    record_batch_metrics(
        captures,
        accumulators,
        records_count,
        passed_count,
        bases_masked_total,
        records_count, // every record is evaluated
    );

    Ok(filter_output_blocks(batch_serial, kept_bytes, rejected_bytes))
}

/// Per-batch body of the template-aware filter steps (`BamTemplateBatch`).
/// Same output contract as [`filter_single_batch`]; templates are filtered
/// whole by [`filter_template_records`].
fn filter_template_batch(
    item: BamTemplateBatch,
    captures: &FilterProcessCaptures,
    accumulators: &PerThreadAccumulator<CollectedFilterMetrics>,
    ceiling: usize,
    rejects: bool,
) -> io::Result<(DecompressedBlock, Option<DecompressedBlock>)> {
    let framed = framed_len(
        item.templates().iter().flat_map(crate::template::Template::records).map(AsRef::as_ref),
    );
    let (batch_serial, templates) = item.into_parts();
    let (mut kept_bytes, mut rejected_bytes) = filter_output_buffers(framed, ceiling, rejects);
    let mut total_records: u64 = 0;
    let mut passed_count: u64 = 0;
    let mut bases_masked_total: u64 = 0;
    let mut evaluated_count: u64 = 0;

    for template in templates {
        let mut template_records: Vec<RawRecord> = template.into_records();
        total_records += template_records.len() as u64;
        let outcome = filter_template_records(&mut template_records, captures)?;
        bases_masked_total += outcome.masked_bases;
        evaluated_count += outcome.evaluated;

        for (record, keep) in template_records.into_iter().zip(outcome.keep) {
            let target = if keep {
                passed_count += 1;
                &mut kept_bytes
            } else if let Some(rejected) = rejected_bytes.as_mut() {
                rejected
            } else {
                continue;
            };
            fgumi_raw_bam::write_framed_record(target, record.as_ref())?;
        }
    }

    record_batch_metrics(
        captures,
        accumulators,
        total_records,
        passed_count,
        bases_masked_total,
        evaluated_count,
    );

    Ok(filter_output_blocks(batch_serial, kept_bytes, rejected_bytes))
}

/// Build the single-read, no-rejects filter step.
///
/// `DecodedRecordBatch → DecompressedBlock`. Parallel, `ByItemOrdinal`.
/// Rejected records are dropped; kept records are serialised to raw BAM bytes.
///
/// `pub(crate)` — consumed only by `ChainBuilder::add_filter`.
#[allow(clippy::type_complexity)]
pub(crate) fn build_filter_step_single_no_rejects(
    limit_bytes: u64,
    captures: FilterProcessCaptures,
    accumulators: Arc<PerThreadAccumulator<CollectedFilterMetrics>>,
) -> ProcessOrdered<
    DecodedRecordBatch,
    DecompressedBlock,
    impl Fn(DecodedRecordBatch) -> io::Result<DecompressedBlock> + Send + Sync + 'static,
> {
    process_ordered::<DecodedRecordBatch, DecompressedBlock, _>(
        "FilterProcess",
        limit_bytes,
        move |item: DecodedRecordBatch| -> io::Result<DecompressedBlock> {
            // No rejects: rejected records are dropped.
            Ok(filter_single_batch(
                item,
                &captures,
                &accumulators,
                reserve_ceiling(limit_bytes),
                false,
            )?
            .0)
        },
    )
}

/// Build the single-read, with-rejects filter step.
///
/// `DecodedRecordBatch → (DecompressedBlock kept, DecompressedBlock rejects)`.
/// Parallel, `ByItemOrdinal`. Branch 0 = kept records; branch 1 = rejected.
///
/// `pub(crate)` — consumed only by `ChainBuilder::add_filter`.
#[allow(clippy::type_complexity)]
pub(crate) fn build_filter_step_single_with_rejects(
    limit_bytes: u64,
    captures: FilterProcessCaptures,
    accumulators: Arc<PerThreadAccumulator<CollectedFilterMetrics>>,
) -> Process2Ordered<
    DecodedRecordBatch,
    DecompressedBlock,
    DecompressedBlock,
    impl Fn(
        DecodedRecordBatch,
    ) -> io::Result<
        crate::pipeline::steps::process::Process2Output<DecompressedBlock, DecompressedBlock>,
    > + Send
    + Sync
    + 'static,
> {
    use crate::pipeline::steps::process::Process2Output;

    const _: fn() = || {
        fn assert_heap<T: crate::pipeline::core::item::HeapSize>() {}
        assert_heap::<DecompressedBlock>();
    };

    process2_ordered::<DecodedRecordBatch, DecompressedBlock, DecompressedBlock, _>(
        "FilterProcess",
        limit_bytes,
        limit_bytes,
        move |item: DecodedRecordBatch|
              -> io::Result<Process2Output<DecompressedBlock, DecompressedBlock>> {
            let (kept, rejected) = filter_single_batch(item, &captures, &accumulators, reserve_ceiling(limit_bytes), true)?;
            let rejects = rejects_block(&kept, rejected);
            Ok(Process2Output::both(kept, rejects))
        },
    )
}

/// Build the single-read, no-rejects filter step on the decode-free fast path.
///
/// `RecordBatch → DecompressedBlock`. Identical behavior to
/// [`build_filter_step_single_no_rejects`] — the same `process_record_raw_call`
/// transform and metric tally — but it iterates borrowed record byte-ranges via
/// [`for_each_raw_record`](crate::pipeline::chains::commands::for_each_raw_record)
/// with one reused scratch `RawRecord`, skipping the per-record `DecodedRecord`
/// allocation and the dead `GroupKey` of the `DecodedRecordBatch` decode path.
///
/// `pub(crate)` — consumed only by `ChainBuilder::add_filter` on a BAM source.
pub(crate) fn build_filter_step_single_no_rejects_raw(
    limit_bytes: u64,
    captures: FilterProcessCaptures,
    accumulators: Arc<PerThreadAccumulator<CollectedFilterMetrics>>,
) -> ProcessOrdered<
    RecordBatch,
    DecompressedBlock,
    impl Fn(RecordBatch) -> io::Result<DecompressedBlock> + Send + Sync + 'static,
> {
    process_ordered::<RecordBatch, DecompressedBlock, _>(
        "FilterProcess",
        limit_bytes,
        move |item: RecordBatch| -> io::Result<DecompressedBlock> {
            Ok(filter_single_raw_batch(
                &item,
                &captures,
                &accumulators,
                reserve_ceiling(limit_bytes),
                false,
            )?
            .0)
        },
    )
}

/// Build the single-read, with-rejects filter step on the decode-free fast path.
///
/// `RecordBatch → (DecompressedBlock kept, DecompressedBlock rejects)`. Behavior
/// identical to [`build_filter_step_single_with_rejects`]; iterates borrowed
/// record byte-ranges with one reused scratch. Branch 0 = kept; branch 1 =
/// rejected.
///
/// `pub(crate)` — consumed only by `ChainBuilder::add_filter` on a BAM source.
#[allow(clippy::type_complexity)]
pub(crate) fn build_filter_step_single_with_rejects_raw(
    limit_bytes: u64,
    captures: FilterProcessCaptures,
    accumulators: Arc<PerThreadAccumulator<CollectedFilterMetrics>>,
) -> Process2Ordered<
    RecordBatch,
    DecompressedBlock,
    DecompressedBlock,
    impl Fn(
        RecordBatch,
    ) -> io::Result<
        crate::pipeline::steps::process::Process2Output<DecompressedBlock, DecompressedBlock>,
    > + Send
    + Sync
    + 'static,
> {
    use crate::pipeline::steps::process::Process2Output;

    process2_ordered::<RecordBatch, DecompressedBlock, DecompressedBlock, _>(
        "FilterProcess",
        limit_bytes,
        limit_bytes,
        move |item: RecordBatch|
              -> io::Result<Process2Output<DecompressedBlock, DecompressedBlock>> {
            let (kept, rejected) = filter_single_raw_batch(&item, &captures, &accumulators, reserve_ceiling(limit_bytes), true)?;
            let rejects = rejects_block(&kept, rejected);
            Ok(Process2Output::both(kept, rejects))
        },
    )
}

/// Build the template-aware, no-rejects filter step.
///
/// `BamTemplateBatch → DecompressedBlock`. Parallel, `ByItemOrdinal`.
/// Templates failing the filter are dropped; kept records are serialised.
///
/// `pub(crate)` — consumed only by `ChainBuilder::add_filter`.
#[allow(clippy::type_complexity)]
pub(crate) fn build_filter_step_template_no_rejects(
    limit_bytes: u64,
    captures: FilterProcessCaptures,
    accumulators: Arc<PerThreadAccumulator<CollectedFilterMetrics>>,
) -> ProcessOrdered<
    BamTemplateBatch,
    DecompressedBlock,
    impl Fn(BamTemplateBatch) -> io::Result<DecompressedBlock> + Send + Sync + 'static,
> {
    process_ordered::<BamTemplateBatch, DecompressedBlock, _>(
        "FilterProcess",
        limit_bytes,
        move |item: BamTemplateBatch| -> io::Result<DecompressedBlock> {
            Ok(filter_template_batch(
                item,
                &captures,
                &accumulators,
                reserve_ceiling(limit_bytes),
                false,
            )?
            .0)
        },
    )
}

/// Build the template-aware, with-rejects filter step.
///
/// `BamTemplateBatch → (DecompressedBlock kept, DecompressedBlock rejects)`.
/// Parallel, `ByItemOrdinal`. Branch 0 = kept records; branch 1 = rejected.
///
/// `pub(crate)` — consumed only by `ChainBuilder::add_filter`.
#[allow(clippy::type_complexity)]
pub(crate) fn build_filter_step_template_with_rejects(
    limit_bytes: u64,
    captures: FilterProcessCaptures,
    accumulators: Arc<PerThreadAccumulator<CollectedFilterMetrics>>,
) -> Process2Ordered<
    BamTemplateBatch,
    DecompressedBlock,
    DecompressedBlock,
    impl Fn(
        BamTemplateBatch,
    ) -> io::Result<
        crate::pipeline::steps::process::Process2Output<DecompressedBlock, DecompressedBlock>,
    > + Send
    + Sync
    + 'static,
> {
    use crate::pipeline::steps::process::Process2Output;

    process2_ordered::<BamTemplateBatch, DecompressedBlock, DecompressedBlock, _>(
        "FilterProcess",
        limit_bytes,
        limit_bytes,
        move |item: BamTemplateBatch|
              -> io::Result<Process2Output<DecompressedBlock, DecompressedBlock>> {
            let (kept, rejected) = filter_template_batch(item, &captures, &accumulators, reserve_ceiling(limit_bytes), true)?;
            let rejects = rejects_block(&kept, rejected);
            Ok(Process2Output::both(kept, rejects))
        },
    )
}

#[cfg(test)]
mod stats_tests {
    use super::write_filter_stats;

    #[test]
    fn write_filter_stats_emits_headered_metric_tsv() {
        let tmp = tempfile::NamedTempFile::new().unwrap();
        write_filter_stats(tmp.path(), 1000, 950, 50).unwrap();
        let content = std::fs::read_to_string(tmp.path()).unwrap();
        let mut lines = content.lines();
        assert_eq!(lines.next().unwrap(), "total_reads\tpassed_reads\tfailed_reads\tpass_rate");
        assert_eq!(lines.next().unwrap(), "1000\t950\t50\t0.95");
        assert!(lines.next().is_none(), "expected exactly a header + one data row");
    }

    #[test]
    fn write_filter_stats_zero_total_is_zero_pass_rate() {
        let tmp = tempfile::NamedTempFile::new().unwrap();
        write_filter_stats(tmp.path(), 0, 0, 0).unwrap();
        let content = std::fs::read_to_string(tmp.path()).unwrap();
        let data = content.lines().nth(1).unwrap();
        assert_eq!(data, "0\t0\t0\t0");
    }
}

#[cfg(test)]
mod tests {
    use std::sync::atomic::Ordering;

    /// The end-of-run warnings name what was skipped (or removed), and say so plainly when every
    /// record was skipped: that run checked nothing.
    #[rstest::rstest]
    #[case::none(10, [0, 0, 0, 0, 0, 0, 0, 0], &[])]
    #[case::some_unaligned(10, [3, 0, 0, 0, 0, 0, 0, 0], &["3 unmapped single-strand records"])]
    #[case::all_unaligned(10, [10, 0, 0, 0, 0, 0, 0, 0], &["none of the 10 records were checked"])]
    #[case::simplex_agreement(10, [0, 4, 0, 0, 0, 0, 0, 0], &["not applied to 4 single-strand records"])]
    #[case::some_no_counts(10, [0, 0, 2, 0, 0, 0, 0, 0], &["2 records without cu/ct"])]
    #[case::all_no_counts(10, [0, 0, 10, 0, 0, 0, 0, 0], &["no record carries cu/ct"])]
    #[case::length_mismatch(10, [0, 0, 0, 5, 0, 0, 0, 0], &["5 records whose cu/ct do not match SEQ"])]
    #[case::tags_removed(10, [0, 0, 0, 0, 7, 0, 0, 0], &["removed from 7 records"])]
    #[case::all_skipped_split_causes(
        10,
        [4, 0, 3, 3, 0, 0, 0, 0],
        &[
            "none of the 10 records were checked by the methylation filters (see the reasons below)",
            "4 unmapped single-strand records",
            "3 records without cu/ct",
            "3 records whose cu/ct do not match SEQ",
        ]
    )]
    #[case::unreversed(10, [0, 0, 0, 0, 0, 6, 0, 0], &["wrong positions on 6 reverse-mapped records"])]
    #[case::dropped_calls(10, [0, 0, 0, 0, 0, 0, 5, 2], &["dropped 5 methylation calls from MM/ML on 2 duplex records"])]
    fn test_methylation_skip_warnings(
        #[case] total: u64,
        #[case] counts: [u64; 8],
        #[case] expected: &[&str],
    ) {
        let skips = crate::commands::filter::MethylationFilterSkips::default();
        let [
            unaligned,
            simplex_agreement,
            no_counts,
            length_mismatch,
            tags_removed,
            unreversed,
            dropped_calls,
            records_with_dropped_calls,
        ] = counts;
        skips.unaligned.store(unaligned, Ordering::Relaxed);
        skips.simplex_agreement.store(simplex_agreement, Ordering::Relaxed);
        skips.no_counts.store(no_counts, Ordering::Relaxed);
        skips.length_mismatch.store(length_mismatch, Ordering::Relaxed);
        skips.tags_removed.store(tags_removed, Ordering::Relaxed);
        skips.unreversed_reverse_strand.store(unreversed, Ordering::Relaxed);
        skips.dropped_calls.store(dropped_calls, Ordering::Relaxed);
        skips.records_with_dropped_calls.store(records_with_dropped_calls, Ordering::Relaxed);
        let warnings = super::methylation_skip_warnings(total, &skips);
        assert_eq!(warnings.len(), expected.len(), "{warnings:?}");
        for (warning, needle) in warnings.iter().zip(expected) {
            assert!(warning.contains(needle), "{warning:?} lacks {needle:?}");
        }
    }
}

#[cfg(test)]
mod template_tests {
    use std::sync::Arc;
    use std::sync::atomic::{AtomicU64, Ordering};

    use fgumi_raw_bam::{RawRecord, SamBuilder as RawSamBuilder, flags};

    use crate::commands::filter::{FilterProcessCaptures, MethylationFilterSkips};
    use crate::consensus_filter::{FilterConfig, MethylationDepthThresholds};
    use crate::sam::SamTag;

    /// An unmapped paired read of template `t` with simplex consensus depth `depth` and no
    /// `cu`/`ct` methylation counts.
    fn unmapped_member(flag: u16, depth: i32) -> RawRecord {
        let mut b = RawSamBuilder::new();
        b.read_name(b"t").flags(flag | flags::UNMAPPED).sequence(b"AAAA").qualities(&[30; 4]);
        b.add_int_tag(SamTag::CD, depth).add_float_tag(SamTag::CE, 0.0_f32);
        b.build()
    }

    /// The methylation-filter warnings are measured against the records the filter evaluated,
    /// not every record read: a template whose R1 fails leaves R2 unevaluated and unclassified,
    /// so a run where no evaluated record carries `cu`/`ct` still says that none did.
    #[test]
    fn test_methylation_warnings_count_only_evaluated_records() {
        let skips = Arc::new(MethylationFilterSkips::default());
        let captures = FilterProcessCaptures {
            config: Arc::new(FilterConfig::new(&[5], &[1.0], &[1.0], None, None, 1.0)),
            reference: None,
            min_base_quality: None,
            should_reverse_tags: false,
            min_mean_base_quality: None,
            max_no_call_fraction: 1.0,
            require_single_strand_agreement: false,
            methylation_depth_thresholds: Some(MethylationDepthThresholds::from_values(&[1])),
            require_strand_methylation_agreement: false,
            min_conversion_fraction: None,
            methylation_mode: fgumi_consensus::MethylationMode::EmSeq,
            ref_names: Arc::new(Vec::new()),
            progress: Arc::new(AtomicU64::new(0)),
            header: noodles::sam::Header::default(),
            methylation_skips: Arc::clone(&skips),
        };
        let mut records = vec![
            unmapped_member(flags::PAIRED | flags::FIRST_SEGMENT, 2),
            unmapped_member(flags::PAIRED | flags::LAST_SEGMENT, 10),
        ];

        let outcome =
            super::filter_template_records(&mut records, &captures).expect("filter the template");

        assert_eq!(
            (outcome, skips.no_counts.load(Ordering::Relaxed)),
            (
                super::TemplateFilterOutcome {
                    keep: vec![false, false],
                    masked_bases: 0,
                    evaluated: 1
                },
                1
            ),
            "only R1 is evaluated (and classified)"
        );
        let warnings = super::methylation_skip_warnings(1, &skips);
        assert_eq!(
            warnings,
            vec![
                "none of the 1 records were checked by the methylation filters: no record carries \
                 cu/ct methylation counts; call consensus with --methylation-mode"
                    .to_string()
            ]
        );
    }
}

#[cfg(test)]
mod reserve_tests {
    use std::sync::Arc;
    use std::sync::atomic::AtomicU64;

    use fgumi_bam_io::DecodedRecord;
    use fgumi_raw_bam::{RawRecord, SamBuilder as RawSamBuilder, flags};

    use crate::commands::filter::{
        CollectedFilterMetrics, FilterProcessCaptures, MethylationFilterSkips,
    };
    use crate::consensus_filter::FilterConfig;
    use crate::per_thread_accumulator::PerThreadAccumulator;
    use crate::pipeline::steps::types::{
        BamTemplateBatch, DecodedRecordBatch, DecompressedBlock, RecordBatch,
    };
    use crate::sam::SamTag;
    use crate::template::Template;

    /// Captures for `--min-reads 3`, no reference, no masking.
    fn captures() -> FilterProcessCaptures {
        FilterProcessCaptures {
            config: Arc::new(FilterConfig::new(&[3], &[1.0], &[1.0], None, None, 1.0)),
            reference: None,
            min_base_quality: None,
            should_reverse_tags: false,
            min_mean_base_quality: None,
            max_no_call_fraction: 1.0,
            require_single_strand_agreement: false,
            methylation_depth_thresholds: None,
            require_strand_methylation_agreement: false,
            min_conversion_fraction: None,
            methylation_mode: fgumi_consensus::MethylationMode::Disabled,
            ref_names: Arc::new(Vec::new()),
            progress: Arc::new(AtomicU64::new(0)),
            header: noodles::sam::Header::default(),
            methylation_skips: Arc::new(MethylationFilterSkips::default()),
        }
    }

    fn accumulators() -> Arc<PerThreadAccumulator<CollectedFilterMetrics>> {
        PerThreadAccumulator::new(1)
    }

    /// `n` unmapped simplex-consensus records in templates of two; one template
    /// in `reject_every` has depth 1 and fails `--min-reads 3`, the rest have
    /// depth 5 and pass. `flag` is OR-ed in (paired flags for the template
    /// fixture).
    fn consensus_records(
        n: usize,
        reject_every: usize,
        flag: impl Fn(usize) -> u16,
    ) -> Vec<RawRecord> {
        (0..n)
            .map(|i| {
                let mut b = RawSamBuilder::new();
                let name = format!("r{:06}", i / 2);
                b.read_name(name.as_bytes())
                    .flags(flags::UNMAPPED | flag(i))
                    .sequence(&[b'A'; 50])
                    .qualities(&[30; 50]);
                let depth = if (i / 2) % reject_every == 0 { 1 } else { 5 };
                b.add_int_tag(SamTag::CD, depth).add_float_tag(SamTag::CE, 0.0_f32);
                b.build()
            })
            .collect()
    }

    fn framed(records: &[RawRecord]) -> usize {
        records.iter().map(|r| 4 + r.len()).sum()
    }

    /// Records of one batch: enough that a full-size rejects reservation would
    /// exceed the queue-fit bound of what one in ten records rejects.
    const BATCH_RECORDS: usize = 2_000;

    /// The output-buffer contract for a batch that keeps most records:
    /// - kept is allocated once at the framed input size (an empty `Vec` regrows
    ///   by doubling to a power-of-two-ish capacity, never exactly `framed`);
    /// - rejects are not reserved: a buffer grown from empty stays within twice
    ///   its length, while a full-size reservation would be charged `framed`;
    /// - every record lands in exactly one of the two.
    fn assert_output_contract(
        kept: &DecompressedBlock,
        rejected: Option<&DecompressedBlock>,
        framed: usize,
        rejects: bool,
    ) {
        assert_eq!(kept.bytes.capacity(), framed, "kept buffer regrew or was mis-sized");
        assert!(kept.bytes.len() < framed, "fixture must drop records");
        if rejects {
            let rejected = rejected.expect("rejects requested");
            assert!(!rejected.bytes.is_empty(), "fixture must reject records");
            assert!(
                rejected.bytes.capacity() <= 2 * rejected.bytes.len(),
                "rejects buffer capacity {} for {} bytes: it was reserved, not grown",
                rejected.bytes.capacity(),
                rejected.bytes.len()
            );
            assert_eq!(
                kept.bytes.len() + rejected.bytes.len(),
                framed,
                "every record is kept or rejected exactly once"
            );
        } else {
            assert!(rejected.is_none());
        }
    }

    /// The default filter path (`RecordBatch`, no `--filter-by-template`).
    #[rstest::rstest]
    fn filter_single_raw_output_honours_the_reservation_contract(
        #[values(false, true)] rejects: bool,
    ) {
        let records = consensus_records(BATCH_RECORDS, 10, |_| 0);
        let framed = framed(&records);
        let batch = RecordBatch::new(0, &records);
        let (kept, rejected) = super::filter_single_raw_batch(
            &batch,
            &captures(),
            &accumulators(),
            usize::MAX,
            rejects,
        )
        .expect("filter batch");
        assert_output_contract(&kept, rejected.as_ref(), framed, rejects);
    }

    /// The decoded path (`DecodedRecordBatch`).
    #[rstest::rstest]
    fn filter_single_output_honours_the_reservation_contract(#[values(false, true)] rejects: bool) {
        let records = consensus_records(BATCH_RECORDS, 10, |_| 0);
        let framed = framed(&records);
        let decoded = records
            .into_iter()
            .map(|r| DecodedRecord::from_raw_bytes(r, fgumi_bam_io::GroupKey::default()))
            .collect();
        let batch = DecodedRecordBatch::new(0, decoded);
        let (kept, rejected) =
            super::filter_single_batch(batch, &captures(), &accumulators(), usize::MAX, rejects)
                .expect("filter batch");
        assert_output_contract(&kept, rejected.as_ref(), framed, rejects);
    }

    /// The template path (`--filter-by-template`): paired R1/R2 per template.
    #[rstest::rstest]
    fn filter_template_output_honours_the_reservation_contract(
        #[values(false, true)] rejects: bool,
    ) {
        let records = consensus_records(BATCH_RECORDS, 10, |i| {
            flags::PAIRED | if i % 2 == 0 { flags::FIRST_SEGMENT } else { flags::LAST_SEGMENT }
        });
        let framed = framed(&records);
        let templates = records
            .chunks(2)
            .map(|pair| Template::from_records(pair.to_vec()).expect("template"))
            .collect();
        let batch = BamTemplateBatch::new(0, templates);
        let (kept, rejected) =
            super::filter_template_batch(batch, &captures(), &accumulators(), usize::MAX, rejects)
                .expect("filter batch");
        assert_output_contract(&kept, rejected.as_ref(), framed, rejects);
    }

    /// A filter that drops nearly everything does not hand its queue the
    /// full-size kept reservation: the kept block is fitted to what it holds.
    #[test]
    fn a_heavy_filter_fits_its_kept_block_for_the_queue() {
        let records = consensus_records(BATCH_RECORDS, 1, |_| 0)
            .into_iter()
            .chain(consensus_records(4, 2, |_| 0).into_iter().skip(2))
            .collect::<Vec<_>>();
        let framed = framed(&records);
        let batch = RecordBatch::new(0, &records);
        let (kept, _) =
            super::filter_single_raw_batch(&batch, &captures(), &accumulators(), usize::MAX, false)
                .expect("filter batch");
        assert!(!kept.bytes.is_empty(), "fixture must keep a record");
        assert!(framed > 2 * kept.bytes.len() + 64 * 1024, "fixture must drop most bytes");
        assert!(
            kept.bytes.capacity() <= 2 * kept.bytes.len() + 64 * 1024,
            "kept capacity {} charged for {} bytes",
            kept.bytes.capacity(),
            kept.bytes.len()
        );
    }

    /// The kept reservation never exceeds the step's per-queue byte budget,
    /// and rejects are not reserved at all.
    #[test]
    fn output_buffers_cap_kept_and_leave_rejects_unreserved() {
        let (kept, rejected) = super::filter_output_buffers(1 << 30, 4 << 20, true);
        assert_eq!(kept.capacity(), 4 << 20);
        assert_eq!(rejected.expect("rejects requested").capacity(), 0);
    }

    /// A with-rejects step whose batch body returned no rejects block is a
    /// wiring bug: it fails loudly instead of writing an empty rejects file.
    #[test]
    #[should_panic(expected = "got no rejects block for batch 7")]
    fn a_missing_rejects_block_panics() {
        let kept = DecompressedBlock { batch_serial: 7, bytes: Vec::new() };
        let _ = super::rejects_block(&kept, None);
    }
}
