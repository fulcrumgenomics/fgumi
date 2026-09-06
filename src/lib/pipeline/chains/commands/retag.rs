//! Chain builder support for `Stage::Retag`.
//!
//! Holds the retag-specific finalize hooks and the process-step factory that
//! `ChainBuilder::add_retag` imports. retag
//! is a pure per-record transform (`DecodedRecordBatch → DecompressedBlock`, no
//! grouping preamble, no rejects), so it mirrors filter's
//! `build_filter_step_single_no_rejects` shape.
//!
//! The per-operation counters (`Vec<OpCounts>`, positionally indexed by op), the
//! per-record apply engine (`apply_op`), and the row-building
//! (`RetagMetric::from_counts`) all live in `commands::retag` and are reused
//! verbatim here, so the metrics rows are defined in exactly one place.
//! `sum_slot_counts` reduces the chain's per-thread accumulator into that shared
//! `Vec<OpCounts>`.

use std::io;
use std::path::PathBuf;
use std::sync::Arc;
use std::sync::atomic::{AtomicU64, Ordering};

use anyhow::Result;
use log::{info, warn};

use crate::commands::retag::{
    OpCounts, RetagMetric, RetagOp, RetagProcessCaptures, apply_op, sum_slot_counts,
};
use crate::logging::OperationTimer;
use crate::per_thread_accumulator::PerThreadAccumulator;
use crate::pipeline::chains::FinalizeHook;
use crate::pipeline::steps::process::{ProcessOrdered, process_ordered};
use crate::pipeline::steps::types::{DecodedRecordBatch, DecompressedBlock};

// ─────────────────────────────────────────────────────────────────────────────
// Finalize hooks
// ─────────────────────────────────────────────────────────────────────────────

/// Always-run finalize hook: logs the `=== Summary ===` banner and completes the
/// timer. Registered on the always-run `finalize` list (matching filter/dedup),
/// so it reports even on a failed/partial run — the summary counts what was
/// processed before the failure rather than being suppressed by an early `?`.
///
/// `record_count` comes from the shared `progress` counter (the process step
/// increments it once per record), NOT from summing `records_applied` across ops
/// — that would undercount whenever an op's source tag is absent on some records.
pub(crate) struct RetagFinalizeHook {
    pub(crate) accumulators: Arc<PerThreadAccumulator<Vec<OpCounts>>>,
    pub(crate) operations: Vec<RetagOp>,
    pub(crate) progress: Arc<AtomicU64>,
    pub(crate) timer: OperationTimer,
}

impl FinalizeHook for RetagFinalizeHook {
    fn finalize(self: Box<Self>) -> Result<()> {
        let RetagFinalizeHook { accumulators, operations, progress, timer } = *self;
        let record_count = progress.load(Ordering::Relaxed);
        let counts = sum_slot_counts(&accumulators, operations.len());

        info!("=== Summary ===");
        info!("Records processed: {record_count}");
        for (op, op_counts) in operations.iter().zip(&counts) {
            info!(
                "{op}: applied={} overwritten={} missing={}",
                op_counts.records_applied, op_counts.dst_overwritten, op_counts.src_missing
            );
        }

        timer.log_completion(record_count);
        Ok(())
    }
}

/// Success-only finalize hook: the warn-on-zero-match loop and the `--metrics`
/// TSV write. Registered on `finalize_on_success` (not the always-run list) so a
/// failed/partial run publishes neither the warning nor the TSV — a zero-match
/// warning or a metrics file are only meaningful once every record was processed.
pub(crate) struct RetagMetricsFinalizeHook {
    pub(crate) accumulators: Arc<PerThreadAccumulator<Vec<OpCounts>>>,
    pub(crate) operations: Vec<RetagOp>,
    pub(crate) metrics: Option<PathBuf>,
}

impl FinalizeHook for RetagMetricsFinalizeHook {
    fn finalize(self: Box<Self>) -> Result<()> {
        let RetagMetricsFinalizeHook { accumulators, operations, metrics } = *self;
        let counts = sum_slot_counts(&accumulators, operations.len());

        // Warn on operations that never matched — the usual sign of a mistyped
        // source tag.
        for (op, op_counts) in operations.iter().zip(&counts) {
            if op_counts.records_applied == 0 {
                warn!(
                    "operation '{op}' matched zero records: no record carried the source tag '{}'",
                    op.src()
                );
            }
        }

        if let Some(path) = &metrics {
            let rows: Vec<RetagMetric> = operations
                .iter()
                .zip(&counts)
                .map(|(op, c)| RetagMetric::from_counts(*op, c))
                .collect();
            fgumi_metrics::write_metrics(path, &rows, "retag")?;
            info!("Wrote metrics to: {}", path.display());
        }

        Ok(())
    }
}

// ─────────────────────────────────────────────────────────────────────────────
// Process-step factory
// ─────────────────────────────────────────────────────────────────────────────

/// Build the retag process step: `DecodedRecordBatch → DecompressedBlock`,
/// parallel, `ByItemOrdinal` (input order preserved). Applies every op to each
/// record via the shared `apply_op`, folding per-op counts into this worker's
/// accumulator slot, then serialises every record (retag keeps all records — no
/// rejects, no filtering).
///
/// Returns the concrete `ProcessOrdered` (not `impl Step`) — the closure type is
/// opaque-return and cannot name itself, matching the filter/dedup factories.
///
/// `pub(crate)` — consumed only by `ChainBuilder::add_retag`.
#[allow(clippy::type_complexity)]
pub(crate) fn build_retag_process_step(
    limit_bytes: u64,
    captures: RetagProcessCaptures,
    accumulators: Arc<PerThreadAccumulator<Vec<OpCounts>>>,
) -> ProcessOrdered<
    DecodedRecordBatch,
    DecompressedBlock,
    impl Fn(DecodedRecordBatch) -> io::Result<DecompressedBlock> + Send + Sync + 'static,
> {
    let n_ops = captures.operations.len();
    process_ordered::<DecodedRecordBatch, DecompressedBlock, _>(
        "RetagProcess",
        limit_bytes,
        move |item: DecodedRecordBatch| -> io::Result<DecompressedBlock> {
            let batch_serial = item.batch_serial();
            let records = item.into_records();
            let total_records = records.len() as u64;
            let mut bytes: Vec<u8> = Vec::new();

            // Accumulate per-op counts into a batch-local vector, then merge into
            // the shared per-thread slot with a SINGLE lock at batch end (mirrors
            // filter's `record_batch_metrics`), rather than locking per record.
            let mut batch_counts = vec![OpCounts::default(); n_ops];
            for decoded in records {
                let mut record = decoded.into_raw_bytes();
                for (op, op_counts) in captures.operations.iter().zip(batch_counts.iter_mut()) {
                    apply_op(&mut record, *op, op_counts);
                }
                fgumi_raw_bam::write_framed_record(&mut bytes, record.as_ref())?;
            }
            accumulators.with_slot(|slot: &mut Vec<OpCounts>| {
                // Lazily size the slot to the operation count on first use
                // (mirrors the pre-cutover run_threaded process closure).
                if slot.is_empty() {
                    slot.resize(n_ops, OpCounts::default());
                }
                for (agg, part) in slot.iter_mut().zip(batch_counts.iter()) {
                    agg.records_applied += part.records_applied;
                    agg.dst_overwritten += part.dst_overwritten;
                    agg.src_missing += part.src_missing;
                }
            });

            // Heartbeat + record-count source for RetagFinalizeHook. A shared
            // `Arc<AtomicU64>` rather than `ProgressTracker`: it doubles as the
            // aggregate record count the finalize hook reads at drain, and is
            // trivially shared across the parallel workers.
            let prev = captures.progress.fetch_add(total_records, Ordering::Relaxed);
            if (prev + total_records) / 1_000_000 > prev / 1_000_000 {
                info!("Processed {} records", prev + total_records);
            }

            Ok(DecompressedBlock { batch_serial, bytes })
        },
    )
}
