//! Chain builder for `Stage::Dedup`.
//!
//! Phase 2 (T2.12) held the full ~500-LOC chain construction here.
//! Phase 3 (T3a.3) lifts that logic into
//! [`crate::pipeline::chains::builder::ChainBuilder`]; this module
//! now holds the dedup-specific types and step factories that the builder
//! imports: `DedupSerializeState`, `DedupFinalizeHook`, and the three
//! step-factory functions used by `ChainBuilder::add_dedup`.

use std::sync::Arc;
use std::sync::atomic::{AtomicU64, Ordering};

use ahash::AHashMap;
use anyhow::{Result, bail};

use crate::commands::dedup::{
    BatchedProcessedDedupGroups, CollectedDedupMetrics, DedupCounts, DuplicationLadderRecorder,
    ProcessedDedupGroup, process_position_group, write_dedup_metrics, write_duplication_ladder,
    write_family_size_histogram,
};
use crate::logging::OperationTimer;
use crate::per_thread_accumulator::PerThreadAccumulator;
use crate::pipeline::chains::FinalizeHook;
use crate::pipeline::steps::group::position::BatchedRawPositionGroups;
use crate::pipeline::steps::process::{
    MiAssign, ProcessOrdered, ProcessWithWorkerState, mi_assign, process_ordered,
    process_with_worker_state,
};
use crate::pipeline::steps::types::DecompressedBlock;
use crate::sam::SamTag;
use crate::template_filter::TemplateFilterConfig;
use fgumi_bam_io::LibraryIndex;
use fgumi_raw_bam::RawRecordView;

/// Serialize worker state for the dedup pipeline's output step.
///
/// Defined at module level (not inside a function body) to avoid the
/// `clippy::items_after_statements` lint that fires when a struct is defined
/// after executable statements in a function body.
///
/// `pub(crate)` so `ChainBuilder` can reference it when constructing the
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
/// `pub(crate)` so `ChainBuilder` can construct and register it in
/// `add_dedup`.
pub(crate) struct DedupFinalizeHook {
    pub(crate) accumulators: Arc<PerThreadAccumulator<CollectedDedupMetrics>>,
    pub(crate) metrics_path: Option<std::path::PathBuf>,
    pub(crate) family_size_histogram_path: Option<std::path::PathBuf>,
    /// `--duplication-ladder` output path paired with the recorder accumulated
    /// by the serial MI-assign step; `Some` exactly when the flag is set. A
    /// single `Option<(path, recorder)>` (rather than two coupled `Option`s)
    /// makes the "both or neither" invariant unrepresentable-when-broken, so
    /// finalize needs no runtime assert. The recorder is held behind an `Arc`
    /// shared with the MI-assign step, so the ladder is read here **through the
    /// lock** rather than by `Arc::try_unwrap` — the step's clone may still be
    /// alive when finalize hooks run.
    pub(crate) duplication_ladder:
        Option<(std::path::PathBuf, Arc<parking_lot::Mutex<DuplicationLadderRecorder>>)>,
    /// Library index resolved from the input header, needed for the per-library
    /// metrics rows.
    pub(crate) library_index: LibraryIndex,
    /// Metrics `sample` value (resolved from `--sample` or the read groups).
    pub(crate) sample: String,
    pub(crate) timer: OperationTimer,
}

impl FinalizeHook for DedupFinalizeHook {
    fn finalize(self: Box<Self>) -> Result<()> {
        let DedupFinalizeHook {
            accumulators,
            metrics_path,
            family_size_histogram_path,
            duplication_ladder,
            library_index,
            sample,
            timer,
        } = *self;

        // Reduce per-thread accumulators into per-library counts + a grand total.
        let mut final_counts_by_library: AHashMap<u16, DedupCounts> = AHashMap::new();
        let mut final_family_sizes: AHashMap<usize, u64> = AHashMap::new();
        for slot in accumulators.slots() {
            let acc = slot.lock();
            for (&library_idx, counts) in &acc.dedup_counts_by_library {
                final_counts_by_library.entry(library_idx).or_default().merge(counts);
            }
            for (&size, &count) in &acc.family_sizes {
                *final_family_sizes.entry(size).or_insert(0) += count;
            }
        }

        // Total across all libraries: the summary log, the `tc`-tag check, and
        // the metrics file's aggregate "All Reads" row.
        let mut final_counts = DedupCounts::default();
        for library_counts in final_counts_by_library.values() {
            final_counts.merge(library_counts);
        }

        // Write metrics file if requested (per-library rows + aggregate).
        if let Some(ref path) = metrics_path {
            write_dedup_metrics(
                &final_counts_by_library,
                &final_counts,
                &library_index,
                &sample,
                path,
            )?;
        }

        // Write family-size histogram if requested.
        if let Some(ref path) = family_size_histogram_path {
            write_family_size_histogram(&final_family_sizes, path)?;
        }

        // Write the duplication-saturation ladder if requested. Accessed through
        // the lock (not `Arc::try_unwrap`): the serial MI-assign step holds a
        // clone of this `Arc` that may still be alive here.
        if let Some((path, recorder)) = duplication_ladder {
            let mut guard = recorder.lock();
            guard.finish();
            write_duplication_ladder(&guard, &library_index, &path)?;
        }

        // Log summary banner (verbatim from the original execute).
        log::info!(
            "Deduplication complete: {} templates ({} unique, {} duplicates, {:.2}% duplicate rate)",
            final_counts.total_templates,
            final_counts.unique_templates,
            final_counts.duplicate_templates,
            final_counts.duplicate_rate() * 100.0
        );

        // Check for missing tc tags (bail verbatim from the original execute).
        if final_counts.missing_tc_tag > 0 {
            bail!(
                "{} secondary/supplementary reads are missing the `tc` tag.\n\n\
                The `tc` tag is required for correct UMI-aware deduplication of \
                secondary and supplementary alignments. This tag is added by \
                `fgumi zipper` during the merge of unmapped and mapped BAMs.\n\n\
                To fix this, re-run your pipeline starting from `fgumi zipper`:\n  \
                fgumi zipper -i aligned.bam --unmapped unmapped.bam -r reference.fa -o merged.bam\n  \
                fgumi sort -i merged.bam -o sorted.bam --order template-coordinate\n  \
                fgumi dedup -i sorted.bam -o deduped.bam",
                final_counts.missing_tc_tag
            );
        }

        // Report filtered templates identically to the non-chain path (shared
        // helper, so the two paths cannot drift).
        crate::commands::dedup::log_filtered_templates(&final_counts.filter_counts);

        timer.log_completion(final_counts.total_reads);

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
/// `pub(crate)` — consumed only by `ChainBuilder::add_dedup`.
#[allow(clippy::type_complexity, clippy::too_many_arguments)]
pub(crate) fn build_process_step(
    limit_bytes: u64,
    filter_config: TemplateFilterConfig,
    effective_strategy: crate::assigner::Strategy,
    effective_edits: u32,
    index_threshold: fgumi_umi::IndexThreshold,
    raw_tag: SamTag,
    min_umi_length: Option<usize>,
    no_umi: bool,
    include_unmapped: bool,
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
                    include_unmapped,
                )?;
                accumulators.with_slot(|acc| {
                    // A position group is always single-library, so its counts
                    // belong entirely to `processed.library_idx`.
                    acc.dedup_counts_by_library
                        .entry(processed.library_idx)
                        .or_default()
                        .merge(&processed.dedup_counts);
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
/// When `ladder_recorder` is `Some`, this step also records the
/// `--duplication-ladder` saturation curve — per position group, in the same
/// serial/coordinate-order seam the non-chain path records it (`mi_assign_fn`),
/// **not** in the parallel serialize step. `MiAssign` is `Serial` +
/// `ByItemOrdinal`, so batches reach this closure in input-record order and the
/// groups within a batch are in coordinate order, exactly reproducing the
/// non-chain per-group stream — see the ordering note on
/// [`crate::commands::dedup::DuplicationLadderRecorder`].
///
/// `pub(crate)` — consumed only by `ChainBuilder::add_dedup`.
#[allow(clippy::type_complexity)]
pub(crate) fn build_mi_assign_step(
    limit_bytes: u64,
    ladder_recorder: Option<Arc<parking_lot::Mutex<DuplicationLadderRecorder>>>,
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
            assign_mi_offsets(next_mi, &mut item)?;
            if let Some(recorder) = &ladder_recorder {
                let mut guard = recorder.lock();
                for processed in &item.groups {
                    guard.record(processed.library_idx, &processed.dedup_counts);
                }
            }
            Ok(item)
        },
    )
}

/// Assign global `MoleculeId` offsets to every template in `item`, advancing the
/// running `next_mi` counter by each group's `distinct_mi_count`.
///
/// Each group reserves a contiguous block `[base, base + distinct_mi_count)`
/// from the global counter and shifts its templates' local MI ids by `base`.
/// Mirrors `chains::commands::group::assign_mi_offsets` (which operates on
/// `BatchedProcessedPositionGroups`); the two are kept in step by design.
///
/// # Errors
///
/// Returns an `io::Error` if the cumulative MI counter would exceed
/// [`u64::MAX`].
fn assign_mi_offsets(
    next_mi: &mut u64,
    item: &mut BatchedProcessedDedupGroups,
) -> std::io::Result<()> {
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
    Ok(())
}

/// Build the `DedupSerialize` step: parallel, `ByItemOrdinal`, serializes
/// each processed batch to raw BAM bytes with MI tags applied. Only included
/// in the chain when the stage is `StagePosition::Terminal`; for
/// `StagePosition::Intermediate` the chain tail stays as
/// [`BatchedProcessedDedupGroups`] for the next stage's input.
///
/// `pub(crate)` — consumed only by `ChainBuilder::add_dedup`.
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
                            fgumi_raw_bam::write_framed_record(&mut output, tagged)?;
                        } else {
                            fgumi_raw_bam::write_framed_record(&mut output, raw)?;
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

#[cfg(test)]
mod tests {
    use super::*;
    use crate::template::Template;
    use fgumi_umi::MoleculeId;

    /// Build a [`ProcessedDedupGroup`] carrying `distinct_mi_count` and the
    /// given local-id templates (each as a `MoleculeId::Single`).
    fn make_group(distinct_mi_count: u64, local_ids: &[u64]) -> ProcessedDedupGroup {
        let templates = local_ids
            .iter()
            .map(|&id| {
                let mut template = Template::new(Vec::new());
                template.mi = MoleculeId::Single(id);
                template
            })
            .collect();
        ProcessedDedupGroup {
            templates,
            family_sizes: AHashMap::default(),
            dedup_counts: DedupCounts::default(),
            library_idx: 0,
            input_record_count: 0,
            distinct_mi_count,
        }
    }

    #[test]
    fn assign_mi_offsets_advances_counter_and_shifts_local_ids() {
        // Two groups; the second must be offset past the first's block.
        let mut item = BatchedProcessedDedupGroups {
            batch_serial: 0,
            groups: vec![make_group(2, &[0, 1]), make_group(3, &[0, 1, 2])],
        };
        let mut next_mi = 10;

        assign_mi_offsets(&mut next_mi, &mut item).expect("no overflow");

        // Counter advanced by 2 + 3 over the starting value of 10.
        assert_eq!(next_mi, 15);
        // First group offset by base 10.
        assert_eq!(item.groups[0].templates[0].mi, MoleculeId::Single(10));
        assert_eq!(item.groups[0].templates[1].mi, MoleculeId::Single(11));
        // Second group offset by base 12 (10 + first group's count of 2).
        assert_eq!(item.groups[1].templates[0].mi, MoleculeId::Single(12));
        assert_eq!(item.groups[1].templates[2].mi, MoleculeId::Single(14));
    }

    #[test]
    fn assign_mi_offsets_errors_on_counter_overflow() {
        // A non-zero starting counter plus a u64::MAX block overflows.
        let mut item = BatchedProcessedDedupGroups {
            batch_serial: 0,
            groups: vec![make_group(u64::MAX, &[])],
        };
        let mut next_mi = 1;

        let err = assign_mi_offsets(&mut next_mi, &mut item)
            .expect_err("cumulative MI counter should overflow");
        assert_eq!(err.kind(), std::io::ErrorKind::Other);
        assert!(err.to_string().contains("MoleculeId offset overflow"));
    }
}
