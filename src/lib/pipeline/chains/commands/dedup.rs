//! Chain builder for `Stage::Dedup`.
//!
//! Phase 2 (T2.12) held the full ~500-LOC chain construction here.
//! Phase 3 (T3a.3) lifts that logic into
//! [`crate::pipeline::chains::builder::ChainBuilder`]; this module
//! now holds the dedup-specific types and step factories that the builder
//! imports: [`DedupSerializeState`], [`DedupFinalizeHook`], and the three
//! step-factory functions used by `ChainBuilder::add_dedup`.
//!
//! [`build_dedup_chain`] shrinks to a ~10-line delegate.

use std::sync::Arc;
use std::sync::atomic::{AtomicU64, Ordering};

use ahash::AHashMap;
use anyhow::{Result, bail};

use crate::commands::dedup::{
    BatchedProcessedDedupGroups, CollectedDedupMetrics, DedupFilterConfig, DedupMetrics,
    ProcessedDedupGroup, process_position_group, write_dedup_metrics, write_family_size_histogram,
};
use crate::logging::OperationTimer;
use crate::per_thread_accumulator::PerThreadAccumulator;
use crate::pipeline::chains::builder::ChainBuilder;
use crate::pipeline::chains::{BuiltPipeline, ChainSpec, FinalizeHook};
use crate::pipeline::steps::group::position::BatchedRawPositionGroups;
use crate::pipeline::steps::process::{
    MiAssign, ProcessOrdered, ProcessWithWorkerState, mi_assign, process_ordered,
    process_with_worker_state,
};
use crate::pipeline::steps::types::DecompressedBlock;
use crate::sam::SamTag;
use fgumi_raw_bam::RawRecordView;

/// Serialize worker state for the dedup pipeline's output step.
///
/// Defined at module level (not inside a function body) to avoid the
/// `clippy::items_after_statements` lint that fires when a struct is defined
/// after executable statements in a function body.
///
/// `pub(crate)` so [`ChainBuilder`] can reference it when constructing the
/// dedup step sequence in `add_dedup`.
pub(crate) struct DedupSerializeState {
    scratch: Vec<u8>,
    mi_buf: String,
}

impl DedupSerializeState {
    pub(crate) fn new() -> Self {
        Self { scratch: Vec::with_capacity(512), mi_buf: String::with_capacity(16) }
    }

    pub(crate) fn mi_buf_mut(&mut self) -> &mut String {
        &mut self.mi_buf
    }

    /// Load `raw` into scratch, update the MI tag in-place, and return the
    /// resulting bytes. Encapsulates the split-borrow needed to update both
    /// `scratch` and use `mi_buf` as the new tag value.
    pub(crate) fn apply_mi_tag(&mut self, raw: &[u8], assign_tag_bytes: [u8; 2]) -> &[u8] {
        self.scratch.clear();
        self.scratch.extend_from_slice(raw);
        fgumi_raw_bam::update_string_tag(
            &mut self.scratch,
            assign_tag_bytes,
            self.mi_buf.as_bytes(),
        );
        &self.scratch
    }
}

impl crate::pipeline::core::item::HeapSize for DedupSerializeState {}

/// Post-pipeline finalize hook for dedup. Reduces per-thread metrics,
/// writes the optional metrics file and family-size histogram, logs the
/// summary banner, checks for missing `tc` tags, and calls
/// `timer.log_completion`.
///
/// `pub(crate)` so [`ChainBuilder`] can construct and register it in
/// `add_dedup`.
pub(crate) struct DedupFinalizeHook {
    pub(crate) accumulators: Arc<PerThreadAccumulator<CollectedDedupMetrics>>,
    pub(crate) metrics_path: Option<std::path::PathBuf>,
    pub(crate) family_size_histogram_path: Option<std::path::PathBuf>,
    pub(crate) timer: OperationTimer,
}

impl FinalizeHook for DedupFinalizeHook {
    fn finalize(self: Box<Self>) -> Result<()> {
        let DedupFinalizeHook { accumulators, metrics_path, family_size_histogram_path, timer } =
            *self;

        // Reduce per-thread accumulators.
        let mut final_metrics = DedupMetrics::default();
        let mut final_family_sizes: AHashMap<usize, u64> = AHashMap::new();
        for slot in accumulators.slots() {
            let acc = slot.lock();
            final_metrics.merge(&acc.dedup_metrics);
            for (&size, &count) in &acc.family_sizes {
                *final_family_sizes.entry(size).or_insert(0) += count;
            }
        }

        // Write metrics file if requested.
        if let Some(ref path) = metrics_path {
            write_dedup_metrics(&final_metrics, path)?;
        }

        // Write family-size histogram if requested.
        if let Some(ref path) = family_size_histogram_path {
            write_family_size_histogram(&final_family_sizes, path)?;
        }

        // Log summary banner (verbatim from the original execute).
        log::info!(
            "Deduplication complete: {} templates ({} unique, {} duplicates, {:.2}% duplicate rate)",
            final_metrics.total_templates,
            final_metrics.unique_templates,
            final_metrics.duplicate_templates,
            final_metrics.duplicate_rate() * 100.0
        );

        // Check for missing tc tags (bail verbatim from the original execute).
        if final_metrics.missing_tc_tag > 0 {
            bail!(
                "{} secondary/supplementary reads are missing the `tc` tag.\n\n\
                The `tc` tag is required for correct UMI-aware deduplication of \
                secondary and supplementary alignments. This tag is added by \
                `fgumi zipper` during the merge of unmapped and mapped BAMs.\n\n\
                To fix this, re-run your pipeline starting from `fgumi zipper`:\n  \
                fgumi zipper -i aligned.bam --unmapped unmapped.bam -r reference.fa -o merged.bam\n  \
                fgumi sort -i merged.bam -o sorted.bam --order template-coordinate\n  \
                fgumi dedup -i sorted.bam -o deduped.bam",
                final_metrics.missing_tc_tag
            );
        }

        timer.log_completion(final_metrics.total_reads);

        Ok(())
    }
}

// ─────────────────────────────────────────────────────────────────────────────
// Step factories — extracted from ChainBuilder::add_dedup (T3a.3 review M2).
//
// Each factory receives its captured state as plain arguments and returns the
// concrete step type (via the constructor's return type). Returning
// `impl Step<...>` is blocked here because the closure types embed in
// opaque-return position and cannot name themselves, so we return the concrete
// `Process*` structs directly (the same pattern used by the correct-step
// factories in `steps::correct`).
// ─────────────────────────────────────────────────────────────────────────────

/// Build the `DedupProcess` step: parallel, `ByItemOrdinal`, processes each
/// position group through [`process_position_group`] and accumulates per-thread
/// metrics.
///
/// `pub(crate)` — consumed only by [`ChainBuilder::add_dedup`].
#[allow(clippy::type_complexity, clippy::too_many_arguments)]
pub(crate) fn build_process_step(
    limit_bytes: u64,
    filter_config: DedupFilterConfig,
    effective_strategy: crate::assigner::Strategy,
    effective_edits: u32,
    index_threshold: usize,
    raw_tag: SamTag,
    min_umi_length: Option<usize>,
    no_umi: bool,
    accumulators: Arc<PerThreadAccumulator<CollectedDedupMetrics>>,
) -> ProcessOrdered<
    BatchedRawPositionGroups,
    BatchedProcessedDedupGroups,
    impl Fn(BatchedRawPositionGroups) -> std::io::Result<BatchedProcessedDedupGroups>
    + Send
    + Sync
    + 'static,
> {
    process_ordered::<BatchedRawPositionGroups, BatchedProcessedDedupGroups, _>(
        "DedupProcess",
        limit_bytes,
        move |item: BatchedRawPositionGroups| -> std::io::Result<BatchedProcessedDedupGroups> {
            let BatchedRawPositionGroups { batch_serial, groups } = item;
            let mut processed_batch: Vec<ProcessedDedupGroup> = Vec::with_capacity(groups.len());
            let assigner =
                effective_strategy.new_assigner_full(effective_edits, 1, index_threshold);
            for group in groups {
                assigner.reset();
                let processed = process_position_group(
                    group,
                    &filter_config,
                    assigner.as_ref(),
                    raw_tag,
                    min_umi_length,
                    no_umi,
                )?;
                accumulators.with_slot(|acc| {
                    acc.dedup_metrics.merge(&processed.dedup_metrics);
                    for (size, count) in &processed.family_sizes {
                        *acc.family_sizes.entry(*size).or_insert(0) += count;
                    }
                });
                processed_batch.push(processed);
            }
            Ok(BatchedProcessedDedupGroups { batch_serial, groups: processed_batch })
        },
    )
}

/// Build the `MiAssignDedup` step: serial, `ByItemOrdinal`, assigns
/// monotonically increasing MI offsets to each batch.
///
/// `pub(crate)` — consumed only by [`ChainBuilder::add_dedup`].
#[allow(clippy::type_complexity)]
pub(crate) fn build_mi_assign_step(
    limit_bytes: u64,
) -> MiAssign<
    BatchedProcessedDedupGroups,
    BatchedProcessedDedupGroups,
    impl Fn(&mut u64, BatchedProcessedDedupGroups) -> std::io::Result<BatchedProcessedDedupGroups>
    + Send
    + Sync
    + 'static,
> {
    mi_assign::<BatchedProcessedDedupGroups, BatchedProcessedDedupGroups, _>(
        "MiAssignDedup",
        limit_bytes,
        move |next_mi: &mut u64, mut item: BatchedProcessedDedupGroups| {
            for processed in &mut item.groups {
                let count = processed.distinct_mi_count;
                let base = *next_mi;
                *next_mi = next_mi.checked_add(count).ok_or_else(|| {
                    std::io::Error::other(
                        "MoleculeId offset overflow: cumulative MI counter \
                                 exceeded u64::MAX",
                    )
                })?;
                for template in &mut processed.templates {
                    template.mi = template.mi.with_offset(base);
                }
            }
            Ok(item)
        },
    )
}

/// Build the `DedupSerialize` step: parallel, `ByItemOrdinal`, serializes
/// each processed batch to raw BAM bytes with MI tags applied. Only included
/// in the chain when the stage is [`StagePosition::Terminal`]; for
/// [`StagePosition::Intermediate`] the chain tail stays as
/// [`BatchedProcessedDedupGroups`] for the next stage's input.
///
/// `pub(crate)` — consumed only by [`ChainBuilder::add_dedup`].
#[allow(clippy::type_complexity)]
pub(crate) fn build_serialize_step(
    limit_bytes: u64,
    remove_duplicates: bool,
    assign_tag_bytes: [u8; 2],
    progress: Arc<AtomicU64>,
) -> ProcessWithWorkerState<
    BatchedProcessedDedupGroups,
    DecompressedBlock,
    impl Fn(&mut DedupSerializeState, BatchedProcessedDedupGroups) -> std::io::Result<DecompressedBlock>
    + Send
    + Sync
    + 'static,
    DedupSerializeState,
    fn() -> DedupSerializeState,
> {
    use crate::commands::dedup::DUPLICATE_FLAG;

    process_with_worker_state::<
        BatchedProcessedDedupGroups,
        DecompressedBlock,
        _,
        DedupSerializeState,
        _,
    >(
        "DedupSerialize",
        limit_bytes,
        DedupSerializeState::new,
        move |state: &mut DedupSerializeState,
              item: BatchedProcessedDedupGroups|
              -> std::io::Result<DecompressedBlock> {
            let BatchedProcessedDedupGroups { batch_serial, groups } = item;
            let total_records: usize = groups
                .iter()
                .map(|g| g.templates.iter().map(|t| t.records().len()).sum::<usize>())
                .sum();
            let mut output = Vec::with_capacity(total_records * 800);
            let mut total_input_records: u64 = 0;
            for processed in &groups {
                total_input_records += processed.input_record_count;
                for template in &processed.templates {
                    let mi = template.mi;
                    let has_mi = mi.is_assigned();
                    if has_mi {
                        mi.write_with_offset(0, state.mi_buf_mut());
                    }
                    for raw in template.records() {
                        if remove_duplicates
                            && (RawRecordView::new(raw).flags() & DUPLICATE_FLAG) != 0
                        {
                            continue;
                        }
                        if has_mi {
                            // `apply_mi_tag` clears scratch, writes raw,
                            // and applies the MI tag — encapsulates the
                            // split-borrow needed for scratch + mi_buf.
                            let tagged = state.apply_mi_tag(raw, assign_tag_bytes);
                            let block_size = u32::try_from(tagged.len()).map_err(|_| {
                                std::io::Error::new(
                                    std::io::ErrorKind::InvalidData,
                                    format!("BAM record too large ({} bytes)", tagged.len()),
                                )
                            })?;
                            output.extend_from_slice(&block_size.to_le_bytes());
                            output.extend_from_slice(tagged);
                        } else {
                            let block_size = u32::try_from(raw.len()).map_err(|_| {
                                std::io::Error::new(
                                    std::io::ErrorKind::InvalidData,
                                    format!("BAM record too large ({} bytes)", raw.len()),
                                )
                            })?;
                            output.extend_from_slice(&block_size.to_le_bytes());
                            output.extend_from_slice(raw);
                        }
                    }
                }
            }

            let prev = progress.fetch_add(total_input_records, Ordering::Relaxed);
            if (prev + total_input_records) / 1_000_000 > prev / 1_000_000 {
                log::info!("Processed {} records", prev + total_input_records);
            }

            Ok(DecompressedBlock { batch_serial, bytes: output })
        },
    )
}

/// Build a [`BuiltPipeline`] for a dedup-only chain.
///
/// Delegates to [`ChainBuilder`]. The full step sequence lives in
/// [`ChainBuilder`]'s `add_dedup` method in
/// [`crate::pipeline::chains::builder`].
///
/// # Errors
///
/// Returns input-validation errors (BAM not template-coordinate sorted,
/// etc.), or any underlying pipeline construction error.
#[allow(clippy::needless_pass_by_value)]
pub fn build_dedup_chain(spec: ChainSpec) -> Result<BuiltPipeline> {
    use crate::pipeline::chains::{Stage, builder::StagePosition};
    let mut chain = ChainBuilder::new(&spec)?;
    chain.add_source()?;
    chain.add_stage(Stage::Dedup, StagePosition::Terminal)?;
    chain.add_sink()?;
    chain.build()
}
