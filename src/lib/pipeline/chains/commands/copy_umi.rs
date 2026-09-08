//! Chain builder support for `Stage::CopyUmi`.
//!
//! Holds the copy-umi-specific types and step factory the builder imports:
//! `CopyUmiFinalizeHook`, `CopyUmiMetricsFinalizeHook`, and
//! `build_copy_umi_process_step`. Copy-umi is a pure per-record transform
//! (`DecodedRecordBatch → DecompressedBlock`, no template grouping, no rejects),
//! so it mirrors filter's `build_filter_step_single_no_rejects` shape.

use std::io;
use std::sync::Arc;
use std::sync::atomic::Ordering;

use anyhow::Result;
use log::info;

use crate::commands::copy_umi::{
    CollectedCopyUmiMetrics, CopyUmiProcessCaptures, copy_umi_into_record,
    warn_and_log_copy_umi_summary, write_copy_umi_metrics,
};
use crate::logging::OperationTimer;
use crate::per_thread_accumulator::PerThreadAccumulator;
use crate::pipeline::chains::FinalizeHook;
use crate::pipeline::steps::process::{ProcessOrdered, process_ordered};
use crate::pipeline::steps::types::{DecodedRecordBatch, DecompressedBlock};

/// Reduce the per-thread accumulators into a single set of totals. The pipeline
/// is fully drained before any finalize hook runs, so every slot holds its final
/// state.
fn reduce(accumulators: &PerThreadAccumulator<CollectedCopyUmiMetrics>) -> CollectedCopyUmiMetrics {
    let mut totals = CollectedCopyUmiMetrics::default();
    for slot in accumulators.slots() {
        let m = slot.lock();
        totals.total_records += m.total_records;
        totals.rx_overwritten += m.rx_overwritten;
        totals.names_trimmed += m.names_trimmed;
    }
    totals
}

/// Success-only finalize hook: emits the overwrite warning + `=== Summary ===`
/// block (via the shared [`warn_and_log_copy_umi_summary`]) and logs completion.
///
/// Registered on `finalize_on_success` (NOT the always-run `finalize`): a bad
/// record aborts the run via `map_err` before the pipeline drains, so an
/// always-run summary would log a partial summary on a fail-fast abort.
pub(crate) struct CopyUmiFinalizeHook {
    pub(crate) accumulators: Arc<PerThreadAccumulator<CollectedCopyUmiMetrics>>,
    pub(crate) timer: OperationTimer,
}

impl FinalizeHook for CopyUmiFinalizeHook {
    fn finalize(self: Box<Self>) -> Result<()> {
        let CopyUmiFinalizeHook { accumulators, timer } = *self;
        let totals = reduce(&accumulators);
        warn_and_log_copy_umi_summary(&totals);
        timer.log_completion(totals.total_records);
        Ok(())
    }
}

/// Success-only finalize hook that writes the `--metrics` TSV (via the shared
/// [`write_copy_umi_metrics`]), so a failed/partial run never publishes counts.
pub(crate) struct CopyUmiMetricsFinalizeHook {
    pub(crate) accumulators: Arc<PerThreadAccumulator<CollectedCopyUmiMetrics>>,
    pub(crate) metrics_path: std::path::PathBuf,
}

impl FinalizeHook for CopyUmiMetricsFinalizeHook {
    fn finalize(self: Box<Self>) -> Result<()> {
        let CopyUmiMetricsFinalizeHook { accumulators, metrics_path } = *self;
        let totals = reduce(&accumulators);
        write_copy_umi_metrics(&metrics_path, &totals)
    }
}

/// Emit the per-batch progress milestone and fold this batch's counts into the
/// per-thread accumulator (mirrors filter's `record_batch_metrics`).
fn record_batch_metrics(
    captures: &CopyUmiProcessCaptures,
    accumulators: &PerThreadAccumulator<CollectedCopyUmiMetrics>,
    total_records: u64,
    rx_overwritten: u64,
    names_trimmed: u64,
) {
    // Cross-thread heartbeat: workers run this concurrently, so the milestone
    // counter is a shared `AtomicU64` rather than `fgumi_bam_io::ProgressTracker`,
    // which is built for a single reader thread, not concurrent workers. The
    // finalize hooks read the record total from the accumulator, not from this
    // counter — it drives only the periodic log.
    let prev = captures.progress.fetch_add(total_records, Ordering::Relaxed);
    if (prev + total_records) / 1_000_000 > prev / 1_000_000 {
        info!("Processed {} records", prev + total_records);
    }
    accumulators.with_slot(|m| {
        m.total_records += total_records;
        m.rx_overwritten += rx_overwritten;
        m.names_trimmed += names_trimmed;
    });
}

/// Build the copy-umi process step (`DecodedRecordBatch → DecompressedBlock`).
///
/// Parallel, `ByItemOrdinal`. Every record is kept (no filtering, no rejects);
/// the read-name UMI is copied into `RX` in place. A bad/empty UMI (or an
/// existing RX under `--fail-if-tag-present`) aborts the run, matching the
/// pre-cutover pipeline `process_fn`. The error crosses an `io::Error`
/// boundary (`map_err`) with its full `anyhow` cause chain flattened into the
/// message via `{e:#}` first, so the pipeline's step-failure reconstruction
/// (which reads only the `io::Error`'s top-level Display) still surfaces the
/// inner cause, not just the outer "extracting UMI from read name" context.
///
/// Returns the concrete `ProcessOrdered` (not `impl Step`), because the closure
/// type embeds in opaque-return position and cannot name itself — the same
/// pattern as the filter/dedup step factories.
#[allow(clippy::type_complexity)]
pub(crate) fn build_copy_umi_process_step(
    limit_bytes: u64,
    captures: CopyUmiProcessCaptures,
    accumulators: Arc<PerThreadAccumulator<CollectedCopyUmiMetrics>>,
) -> ProcessOrdered<
    DecodedRecordBatch,
    DecompressedBlock,
    impl Fn(DecodedRecordBatch) -> io::Result<DecompressedBlock> + Send + Sync + 'static,
> {
    process_ordered::<DecodedRecordBatch, DecompressedBlock, _>(
        "CopyUmiProcess",
        limit_bytes,
        move |item: DecodedRecordBatch| -> io::Result<DecompressedBlock> {
            let batch_serial = item.batch_serial();
            let records = item.into_records();
            let total_records = records.len() as u64;
            let mut bytes: Vec<u8> = Vec::new();
            let mut rx_overwritten: u64 = 0;
            let mut names_trimmed: u64 = 0;

            for decoded in records {
                let mut record = decoded.into_raw_bytes();
                let outcome = copy_umi_into_record(
                    &mut record,
                    captures.field_delimiter,
                    captures.reverse_complement_prefixed,
                    captures.remove_umi,
                    captures.fail_if_tag_present,
                )
                // `{e:#}` (anyhow's alternate Display) joins the full cause chain into
                // one string before crossing the `io::Error` boundary: `io::Error`'s own
                // Display only ever shows its wrapped error's top-level Display, so a
                // bare `io::Error::other(e)` would silently drop the inner cause (e.g.
                // "Invalid UMI ... illegal character ...") once the pipeline reconstructs
                // a step failure from this `io::Error` — exactly the diagnostic the
                // no-`--threads` path used to show before the chain became the only path.
                .map_err(|e| io::Error::other(format!("{e:#}")))?;
                if outcome.overwrote_rx {
                    rx_overwritten += 1;
                }
                if outcome.trimmed_name {
                    names_trimmed += 1;
                }
                fgumi_raw_bam::write_framed_record(&mut bytes, record.as_ref())?;
            }

            record_batch_metrics(
                &captures,
                &accumulators,
                total_records,
                rx_overwritten,
                names_trimmed,
            );

            Ok(DecompressedBlock { batch_serial, bytes })
        },
    )
}
