//! Chain builder for `Stage::Filter`.
//!
//! Phase 2 (T2.14) held the full ~500-LOC chain construction here.
//! Phase 3 (T3a.4) lifts that logic into
//! [`crate::pipeline::chains::builder::ChainBuilder`]; this module
//! now holds the filter-specific types and step factories that the builder
//! imports: [`FilterFinalizeHook`] and the step-factory functions used by
//! `ChainBuilder::add_filter`.
//!
//! [`build_filter_chain`] shrinks to a ~10-line delegate.
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

use ahash::AHashMap;
use anyhow::Result;
use fgumi_raw_bam::RawRecord;
use log::info;

use crate::commands::filter::{CollectedFilterMetrics, Filter, FilterProcessCaptures};
use crate::logging::OperationTimer;
use crate::per_thread_accumulator::PerThreadAccumulator;
use crate::pipeline::chains::builder::ChainBuilder;
use crate::pipeline::chains::{BuiltPipeline, ChainSpec, FinalizeHook};
use crate::pipeline::steps::process::{
    Process2Ordered, ProcessOrdered, process_ordered, process2_ordered,
};
use crate::pipeline::steps::types::{BamTemplateBatch, DecodedRecordBatch, DecompressedBlock};

// ─────────────────────────────────────────────────────────────────────────────
// FilterFinalizeHook
// ─────────────────────────────────────────────────────────────────────────────

/// Post-pipeline finalize hook for filter. Reduces per-thread metrics,
/// writes the optional stats file, and logs the summary banner.
///
/// `pub(crate)` so [`ChainBuilder`] can construct and register it in
/// `add_filter`.
pub(crate) struct FilterFinalizeHook {
    pub(crate) accumulators: Arc<PerThreadAccumulator<CollectedFilterMetrics>>,
    pub(crate) stats_path: Option<std::path::PathBuf>,
    pub(crate) has_rejects: bool,
    pub(crate) timer: OperationTimer,
}

impl FinalizeHook for FilterFinalizeHook {
    fn finalize(self: Box<Self>) -> Result<()> {
        let FilterFinalizeHook { accumulators, stats_path, has_rejects, timer } = *self;

        let mut total_reads = 0u64;
        let mut passed_reads = 0u64;
        let mut failed_reads = 0u64;
        let mut total_bases_masked = 0u64;
        for slot in accumulators.slots() {
            let m = slot.lock();
            total_reads += m.total_records;
            passed_reads += m.passed_records;
            failed_reads += m.failed_records;
            total_bases_masked += m.total_bases_masked;
        }

        if let Some(ref path) = stats_path {
            write_filter_stats(path, total_reads, passed_reads, failed_reads)?;
        }

        info!("Processed {total_reads} reads; kept {passed_reads} and rejected {failed_reads}");
        if has_rejects && failed_reads > 0 {
            info!("Wrote {failed_reads} rejected records to rejects file");
        }
        info!("Total bases masked: {total_bases_masked}");

        timer.log_completion(total_reads);

        Ok(())
    }
}

// ─────────────────────────────────────────────────────────────────────────────
// Shared helper
// ─────────────────────────────────────────────────────────────────────────────

/// Write filtering statistics to a file.
pub(crate) fn write_filter_stats(
    path: &std::path::Path,
    total: u64,
    passed: u64,
    failed: u64,
) -> Result<()> {
    use std::fs::File;
    use std::io::Write;

    let mut file = File::create(path)?;
    writeln!(file, "total_reads\t{total}")?;
    writeln!(file, "passed_reads\t{passed}")?;
    writeln!(file, "failed_reads\t{failed}")?;
    #[allow(clippy::cast_precision_loss)]
    let pass_rate = if total > 0 { passed as f64 / total as f64 } else { 0.0 };
    writeln!(file, "pass_rate\t{pass_rate:.4}")?;
    Ok(())
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

/// Build the single-read, no-rejects filter step.
///
/// `DecodedRecordBatch → DecompressedBlock`. Parallel, `ByItemOrdinal`.
/// Rejected records are dropped; kept records are serialised to raw BAM bytes.
///
/// `pub(crate)` — consumed only by [`ChainBuilder::add_filter`].
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
            let DecodedRecordBatch { batch_serial, records, .. } = item;
            let records_count = records.len() as u64;
            let mut kept_bytes: Vec<u8> = Vec::new();
            let mut passed_count: u64 = 0;
            let mut bases_masked_total: u64 = 0;

            for decoded in records {
                let mut record = decoded.into_raw_bytes();
                let (bases_masked, pass) =
                    process_record_raw_call(&mut record, &captures).map_err(io::Error::other)?;
                bases_masked_total += bases_masked;

                if pass {
                    passed_count += 1;
                    let block_size = u32::try_from(record.len()).map_err(|_| {
                        io::Error::new(
                            io::ErrorKind::InvalidData,
                            format!("BAM record too large ({} bytes)", record.len()),
                        )
                    })?;
                    kept_bytes.extend_from_slice(&block_size.to_le_bytes());
                    kept_bytes.extend_from_slice(record.as_ref());
                }
                // No rejects: rejected records are simply dropped.
            }

            let prev = captures.progress.fetch_add(records_count, Ordering::Relaxed);
            if (prev + records_count) / 1_000_000 > prev / 1_000_000 {
                info!("Processed {} records", prev + records_count);
            }
            accumulators.with_slot(|m| {
                m.total_records += records_count;
                m.passed_records += passed_count;
                m.failed_records += records_count - passed_count;
                m.total_bases_masked += bases_masked_total;
            });

            Ok(DecompressedBlock { batch_serial, bytes: kept_bytes })
        },
    )
}

/// Build the single-read, with-rejects filter step.
///
/// `DecodedRecordBatch → (DecompressedBlock kept, DecompressedBlock rejects)`.
/// Parallel, `ByItemOrdinal`. Branch 0 = kept records; branch 1 = rejected.
///
/// `pub(crate)` — consumed only by [`ChainBuilder::add_filter`].
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
            let DecodedRecordBatch { batch_serial, records, .. } = item;
            let records_count = records.len() as u64;
            let mut kept_bytes: Vec<u8> = Vec::new();
            let mut rejected_bytes: Vec<u8> = Vec::new();
            let mut passed_count: u64 = 0;
            let mut bases_masked_total: u64 = 0;

            for decoded in records {
                let mut record = decoded.into_raw_bytes();
                let (bases_masked, pass) = process_record_raw_call(&mut record, &captures)
                    .map_err(io::Error::other)?;
                bases_masked_total += bases_masked;

                let target = if pass {
                    passed_count += 1;
                    &mut kept_bytes
                } else {
                    &mut rejected_bytes
                };
                let block_size = u32::try_from(record.len()).map_err(|_| {
                    io::Error::new(
                        io::ErrorKind::InvalidData,
                        format!("BAM record too large ({} bytes)", record.len()),
                    )
                })?;
                target.extend_from_slice(&block_size.to_le_bytes());
                target.extend_from_slice(record.as_ref());
            }

            let prev = captures.progress.fetch_add(records_count, Ordering::Relaxed);
            if (prev + records_count) / 1_000_000 > prev / 1_000_000 {
                info!("Processed {} records", prev + records_count);
            }
            accumulators.with_slot(|m| {
                m.total_records += records_count;
                m.passed_records += passed_count;
                m.failed_records += records_count - passed_count;
                m.total_bases_masked += bases_masked_total;
            });

            Ok(Process2Output::both(
                DecompressedBlock { batch_serial, bytes: kept_bytes },
                DecompressedBlock { batch_serial, bytes: rejected_bytes },
            ))
        },
    )
}

/// Build the template-aware, no-rejects filter step.
///
/// `BamTemplateBatch → DecompressedBlock`. Parallel, `ByItemOrdinal`.
/// Templates failing the filter are dropped; kept records are serialised.
///
/// `pub(crate)` — consumed only by [`ChainBuilder::add_filter`].
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
    use crate::consensus_filter::template_passes;

    process_ordered::<BamTemplateBatch, DecompressedBlock, _>(
        "FilterProcess",
        limit_bytes,
        move |item: BamTemplateBatch| -> io::Result<DecompressedBlock> {
            let BamTemplateBatch { batch_serial, templates, .. } = item;
            let mut kept_bytes: Vec<u8> = Vec::new();
            let mut total_records: u64 = 0;
            let mut passed_count: u64 = 0;
            let mut bases_masked_total: u64 = 0;

            for template in templates {
                let mut template_records: Vec<RawRecord> = template.into_records();
                let mut pass_map: AHashMap<usize, bool> = AHashMap::new();

                for (idx, record) in template_records.iter_mut().enumerate() {
                    total_records += 1;
                    let (masked, pass) =
                        process_record_raw_call(record, &captures).map_err(io::Error::other)?;
                    bases_masked_total += masked;
                    pass_map.insert(idx, pass);
                }

                let template_pass = template_passes(&template_records, &pass_map);

                for (idx, record) in template_records.into_iter().enumerate() {
                    let flags = fgumi_raw_bam::RawRecordView::new(&record).flags();
                    let is_primary = (flags & fgumi_raw_bam::flags::SECONDARY) == 0
                        && (flags & fgumi_raw_bam::flags::SUPPLEMENTARY) == 0;
                    let keep = if is_primary {
                        template_pass
                    } else {
                        template_pass && pass_map.get(&idx).copied().unwrap_or(false)
                    };
                    if keep {
                        passed_count += 1;
                        let block_size = u32::try_from(record.len()).map_err(|_| {
                            io::Error::new(
                                io::ErrorKind::InvalidData,
                                format!("BAM record too large ({} bytes)", record.len()),
                            )
                        })?;
                        kept_bytes.extend_from_slice(&block_size.to_le_bytes());
                        kept_bytes.extend_from_slice(record.as_ref());
                    }
                }
            }

            let prev = captures.progress.fetch_add(total_records, Ordering::Relaxed);
            if (prev + total_records) / 1_000_000 > prev / 1_000_000 {
                info!("Processed {} records", prev + total_records);
            }
            accumulators.with_slot(|m| {
                m.total_records += total_records;
                m.passed_records += passed_count;
                m.failed_records += total_records - passed_count;
                m.total_bases_masked += bases_masked_total;
            });

            Ok(DecompressedBlock { batch_serial, bytes: kept_bytes })
        },
    )
}

/// Build the template-aware, with-rejects filter step.
///
/// `BamTemplateBatch → (DecompressedBlock kept, DecompressedBlock rejects)`.
/// Parallel, `ByItemOrdinal`. Branch 0 = kept records; branch 1 = rejected.
///
/// `pub(crate)` — consumed only by [`ChainBuilder::add_filter`].
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
    use crate::consensus_filter::template_passes;
    use crate::pipeline::steps::process::Process2Output;

    process2_ordered::<BamTemplateBatch, DecompressedBlock, DecompressedBlock, _>(
        "FilterProcess",
        limit_bytes,
        limit_bytes,
        move |item: BamTemplateBatch|
              -> io::Result<Process2Output<DecompressedBlock, DecompressedBlock>> {
            let BamTemplateBatch { batch_serial, templates, .. } = item;
            let mut kept_bytes: Vec<u8> = Vec::new();
            let mut rejected_bytes: Vec<u8> = Vec::new();
            let mut total_records: u64 = 0;
            let mut passed_count: u64 = 0;
            let mut bases_masked_total: u64 = 0;

            for template in templates {
                let mut template_records: Vec<RawRecord> = template.into_records();
                let mut pass_map: AHashMap<usize, bool> = AHashMap::new();

                for (idx, record) in template_records.iter_mut().enumerate() {
                    total_records += 1;
                    let (masked, pass) = process_record_raw_call(record, &captures)
                        .map_err(io::Error::other)?;
                    bases_masked_total += masked;
                    pass_map.insert(idx, pass);
                }

                let template_pass = template_passes(&template_records, &pass_map);

                for (idx, record) in template_records.into_iter().enumerate() {
                    let flags = fgumi_raw_bam::RawRecordView::new(&record).flags();
                    let is_primary = (flags & fgumi_raw_bam::flags::SECONDARY) == 0
                        && (flags & fgumi_raw_bam::flags::SUPPLEMENTARY) == 0;

                    let keep = if is_primary {
                        template_pass
                    } else {
                        template_pass && pass_map.get(&idx).copied().unwrap_or(false)
                    };
                    let target = if keep {
                        passed_count += 1;
                        &mut kept_bytes
                    } else {
                        &mut rejected_bytes
                    };
                    let block_size = u32::try_from(record.len()).map_err(|_| {
                        io::Error::new(
                            io::ErrorKind::InvalidData,
                            format!("BAM record too large ({} bytes)", record.len()),
                        )
                    })?;
                    target.extend_from_slice(&block_size.to_le_bytes());
                    target.extend_from_slice(record.as_ref());
                }
            }

            let prev = captures.progress.fetch_add(total_records, Ordering::Relaxed);
            if (prev + total_records) / 1_000_000 > prev / 1_000_000 {
                info!("Processed {} records", prev + total_records);
            }
            accumulators.with_slot(|m| {
                m.total_records += total_records;
                m.passed_records += passed_count;
                m.failed_records += total_records - passed_count;
                m.total_bases_masked += bases_masked_total;
            });

            Ok(Process2Output::both(
                DecompressedBlock { batch_serial, bytes: kept_bytes },
                DecompressedBlock { batch_serial, bytes: rejected_bytes },
            ))
        },
    )
}

// ─────────────────────────────────────────────────────────────────────────────
// build_filter_chain — thin delegate
// ─────────────────────────────────────────────────────────────────────────────

/// Build a [`BuiltPipeline`] for a filter-only chain.
///
/// Delegates to [`ChainBuilder`]. The full step sequence lives in
/// [`ChainBuilder`]'s `add_filter` method in
/// [`crate::pipeline::chains::builder`].
///
/// # Errors
///
/// Returns input-validation errors or any underlying pipeline construction error.
#[allow(clippy::needless_pass_by_value)]
pub fn build_filter_chain(spec: ChainSpec) -> Result<BuiltPipeline> {
    use crate::pipeline::chains::{Stage, builder::StagePosition};
    let mut chain = ChainBuilder::new(&spec)?;
    chain.add_source()?;
    chain.add_stage(Stage::Filter, StagePosition::Terminal)?;
    chain.add_sink()?;
    chain.build()
}
