//! Chain builder for `Stage::Group`.
//!
//! Phase 2 (T2.17) held the full ~500-LOC chain construction here.
//! Phase 3 (T3a.7) lifts that logic into
//! [`crate::pipeline::chains::builder::ChainBuilder`]; this module
//! now holds the group-specific types and step factories that the builder
//! imports: [`GroupFinalizeHook`] and the three step-factory functions used
//! by `ChainBuilder::add_group`.
//!
//! [`build_group_chain`] shrinks to a ~10-line delegate.
//!
//! Mirrors the dedup (6423252), filter (8eead93), clip (8e95428) reference
//! migrations.

use std::io;
use std::sync::Arc;

use ahash::AHashMap;
use anyhow::Result;
use log::warn;

use crate::assigner::Strategy;
use crate::commands::group::{
    GroupFilterConfig, GroupMetricsAccumulator, assign_umi_groups_impl, build_grouping_metrics,
    filter_template_raw, write_metrics_for_chain,
};
use crate::grouper::{FilterMetrics, ProcessedPositionGroup, build_templates_from_records};
use crate::logging::OperationTimer;
use crate::per_thread_accumulator::PerThreadAccumulator;
use crate::pipeline::chains::builder::ChainBuilder;
use crate::pipeline::chains::{BuiltPipeline, ChainSpec, FinalizeHook};
use crate::pipeline::steps::group::position::{
    BatchedProcessedPositionGroups, BatchedRawPositionGroups,
};
use crate::pipeline::steps::process::{
    MiAssign, ProcessOrdered, ProcessWithWorkerState, mi_assign, process_ordered,
};
use crate::pipeline::steps::serialize_processed::{
    SerializeState, build_serialize_processed_groups_step,
};
use crate::pipeline::steps::types::DecompressedBlock;
use crate::template::Template;
use crate::umi::parallel_assigner::{
    ParallelAdjacencyAssigner, ParallelEditAssigner, ParallelIdentityAssigner,
    ParallelPairedAssigner,
};

// ─────────────────────────────────────────────────────────────────────────────
// GroupFinalizeHook
// ─────────────────────────────────────────────────────────────────────────────

/// Post-pipeline finalize hook for group. Reduces per-thread metric
/// accumulators, writes the optional metrics files (family-size histogram,
/// grouping-metrics, position-group-size histogram, and `--metrics PREFIX`
/// outputs), logs the UMI grouping summary banner, and calls
/// `timer.log_completion`.
///
/// `pub(crate)` so [`ChainBuilder`] can construct and register it in
/// `add_group`.
pub(crate) struct GroupFinalizeHook {
    pub(crate) accumulators: Arc<PerThreadAccumulator<GroupMetricsAccumulator>>,
    pub(crate) output_path: std::path::PathBuf,
    pub(crate) family_size_histogram: Option<std::path::PathBuf>,
    pub(crate) grouping_metrics: Option<std::path::PathBuf>,
    pub(crate) metrics_prefix: Option<std::path::PathBuf>,
    pub(crate) timer: OperationTimer,
}

impl FinalizeHook for GroupFinalizeHook {
    fn finalize(self: Box<Self>) -> Result<()> {
        use crate::logging::log_umi_grouping_summary;
        use log::info;

        let GroupFinalizeHook {
            accumulators,
            output_path,
            family_size_histogram,
            grouping_metrics,
            metrics_prefix,
            timer,
        } = *self;

        // Reduce per-thread accumulators.
        let mut family_size_counter: AHashMap<usize, u64> = AHashMap::new();
        let mut position_group_size_counter: AHashMap<usize, u64> = AHashMap::new();
        let mut total_filter_metrics = FilterMetrics::new();
        for slot in accumulators.slots() {
            let acc = slot.lock();
            for (&size, &count) in &acc.family_sizes {
                *family_size_counter.entry(size).or_insert(0) += count;
            }
            for (&size, &count) in &acc.position_group_sizes {
                *position_group_size_counter.entry(size).or_insert(0) += count;
            }
            total_filter_metrics.merge(&acc.filter_metrics);
        }

        let metrics = build_grouping_metrics(&total_filter_metrics, &family_size_counter);

        // Log summary banner (verbatim from the original execute).
        log_umi_grouping_summary(&metrics);

        // Write individual flag outputs and --metrics prefix outputs when any
        // output path was requested.
        if family_size_histogram.is_some() || grouping_metrics.is_some() || metrics_prefix.is_some()
        {
            write_metrics_for_chain(
                &metrics,
                &family_size_counter,
                &position_group_size_counter,
                family_size_histogram.as_deref(),
                grouping_metrics.as_deref(),
                metrics_prefix.as_deref(),
            )?;
        }

        timer.log_completion(metrics.accepted_records);

        info!("Wrote output to {}", output_path.display());
        info!("group completed successfully");
        info!("Records accepted by group: {}", metrics.accepted_records);

        Ok(())
    }
}

// ─────────────────────────────────────────────────────────────────────────────
// Step factories — extracted from build_group_chain (T3a.7).
//
// Each factory receives its captured state as plain arguments and returns the
// concrete step type. Returning `impl Step<...>` is blocked here because the
// closure types embed in opaque-return position and cannot name themselves, so
// we return the concrete `Process*` structs directly (the same pattern used by
// the dedup factories in `chains::commands::dedup`).
// ─────────────────────────────────────────────────────────────────────────────

/// Build the `GroupProcess` step: parallel, `ByItemOrdinal`. Filters
/// templates, assigns UMI groups, tracks per-thread family-size and filter
/// metrics, and emits [`BatchedProcessedPositionGroups`].
///
/// `pub(crate)` — consumed only by [`ChainBuilder::add_group`].
#[allow(clippy::type_complexity, clippy::too_many_arguments, clippy::too_many_lines)]
pub(crate) fn build_group_process_step(
    limit_bytes: u64,
    strategy: Strategy,
    effective_edits: u32,
    index_threshold: usize,
    raw_tag: [u8; 2],
    no_umi: bool,
    allow_unmapped: bool,
    num_threads: usize,
    filter_config: GroupFilterConfig,
    accumulators: Arc<PerThreadAccumulator<GroupMetricsAccumulator>>,
) -> ProcessOrdered<
    BatchedRawPositionGroups,
    BatchedProcessedPositionGroups,
    impl Fn(BatchedRawPositionGroups) -> io::Result<BatchedProcessedPositionGroups>
    + Send
    + Sync
    + 'static,
> {
    process_ordered::<BatchedRawPositionGroups, BatchedProcessedPositionGroups, _>(
        "GroupProcess",
        limit_bytes,
        move |item: BatchedRawPositionGroups| -> io::Result<BatchedProcessedPositionGroups> {
            let BatchedRawPositionGroups { batch_serial, groups } = item;
            let mut processed_batch: Vec<ProcessedPositionGroup> = Vec::with_capacity(groups.len());

            let assigner: Box<dyn crate::assigner::UmiAssigner> = if allow_unmapped {
                match strategy {
                    Strategy::Identity => Box::new(ParallelIdentityAssigner::new(num_threads)),
                    Strategy::Edit => {
                        Box::new(ParallelEditAssigner::new(effective_edits, num_threads))
                    }
                    Strategy::Adjacency => {
                        Box::new(ParallelAdjacencyAssigner::new(effective_edits, num_threads))
                    }
                    Strategy::Paired => {
                        Box::new(ParallelPairedAssigner::new(effective_edits, num_threads))
                    }
                }
            } else {
                strategy.new_assigner_full(effective_edits, 1, index_threshold)
            };

            for group in groups {
                let mut filter_metrics = FilterMetrics::new();
                let all_templates =
                    build_templates_from_records(group.records).map_err(io::Error::other)?;
                let input_record_count: u64 =
                    all_templates.iter().map(|t| t.read_count() as u64).sum();

                let filtered_templates: Vec<Template> = all_templates
                    .into_iter()
                    .filter(|t| filter_template_raw(t, &filter_config, &mut filter_metrics))
                    .collect();

                if filtered_templates.is_empty() {
                    accumulators.with_slot(|acc| {
                        acc.record_group(AHashMap::new(), &filter_metrics);
                    });
                    processed_batch.push(ProcessedPositionGroup {
                        templates: Vec::new(),
                        family_sizes: AHashMap::new(),
                        filter_metrics,
                        input_record_count,
                        distinct_mi_count: 0,
                    });
                    continue;
                }

                assigner.reset();

                let mut templates = filtered_templates;
                if let Err(e) = assign_umi_groups_impl(
                    &mut templates,
                    assigner.as_ref(),
                    raw_tag,
                    filter_config.min_umi_length,
                    no_umi,
                ) {
                    warn!("UMI assignment failed, returning empty group: {e}");
                    // The whole group is dropped from the BAM, so none of its
                    // filter-passing templates are actually accepted. Zero the
                    // accepted counter before recording so the summary/metrics
                    // match the (empty) output rather than over-counting.
                    let mut dropped_metrics = filter_metrics.clone();
                    dropped_metrics.accepted_templates = 0;
                    accumulators.with_slot(|acc| {
                        acc.record_group(AHashMap::new(), &dropped_metrics);
                    });
                    processed_batch.push(ProcessedPositionGroup {
                        templates: Vec::new(),
                        family_sizes: AHashMap::new(),
                        filter_metrics,
                        input_record_count,
                        distinct_mi_count: 0,
                    });
                    continue;
                }

                let distinct_mi_count: u64 =
                    templates.iter().filter_map(|t| t.mi.id()).max().map_or(0, |max_id| max_id + 1);

                templates.sort_by(|a, b| {
                    let a_idx = a.mi.to_vec_index();
                    let b_idx = b.mi.to_vec_index();
                    a_idx.cmp(&b_idx).then_with(|| a.name.cmp(&b.name))
                });

                let mut family_sizes: AHashMap<usize, u64> = AHashMap::with_capacity(50);
                if !templates.is_empty() {
                    let mut current_mi = templates[0].mi.to_vec_index();
                    let mut current_count = 1usize;
                    for template in templates.iter().skip(1) {
                        let mi = template.mi.to_vec_index();
                        if mi == current_mi {
                            current_count += 1;
                        } else {
                            if current_mi.is_some() {
                                *family_sizes.entry(current_count).or_insert(0) += 1;
                            }
                            current_mi = mi;
                            current_count = 1;
                        }
                    }
                    if current_mi.is_some() {
                        *family_sizes.entry(current_count).or_insert(0) += 1;
                    }
                }
                accumulators.with_slot(|acc| {
                    acc.record_group(family_sizes.clone(), &filter_metrics);
                });

                processed_batch.push(ProcessedPositionGroup {
                    templates,
                    family_sizes,
                    filter_metrics,
                    input_record_count,
                    distinct_mi_count,
                });
            }

            Ok(BatchedProcessedPositionGroups::new(batch_serial, processed_batch))
        },
    )
}

/// Build the `MiAssignGroups` step: serial, `ByItemOrdinal`. Assigns
/// monotonically increasing MI offsets to each batch of processed position
/// groups.
///
/// Group uses [`BatchedProcessedPositionGroups`] (not the dedup-specific
/// `BatchedProcessedDedupGroups`), so this factory is distinct from
/// `chains::commands::dedup::build_mi_assign_step`.
///
/// `pub(crate)` — consumed only by [`ChainBuilder::add_group`].
#[allow(clippy::type_complexity)]
pub(crate) fn build_group_mi_assign_step(
    limit_bytes: u64,
) -> MiAssign<
    BatchedProcessedPositionGroups,
    BatchedProcessedPositionGroups,
    impl Fn(&mut u64, BatchedProcessedPositionGroups) -> io::Result<BatchedProcessedPositionGroups>
    + Send
    + Sync
    + 'static,
> {
    mi_assign::<BatchedProcessedPositionGroups, BatchedProcessedPositionGroups, _>(
        "MiAssignGroups",
        limit_bytes,
        move |next_mi: &mut u64, mut item: BatchedProcessedPositionGroups| {
            assign_mi_offsets(next_mi, &mut item)?;
            Ok(item)
        },
    )
}

/// Assign global `MoleculeId` offsets to every template in `item`, advancing the
/// running `next_mi` counter by each group's `distinct_mi_count`.
///
/// Each group reserves a contiguous block `[base, base + distinct_mi_count)`
/// from the global counter and shifts its templates' local MI ids by `base`.
///
/// # Errors
///
/// Returns an [`io::Error`] if the cumulative MI counter would exceed
/// [`u64::MAX`].
fn assign_mi_offsets(
    next_mi: &mut u64,
    item: &mut BatchedProcessedPositionGroups,
) -> io::Result<()> {
    for processed in &mut item.groups {
        let count = processed.distinct_mi_count;
        let base = *next_mi;
        *next_mi = next_mi.checked_add(count).ok_or_else(|| {
            io::Error::other(
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

/// Build the `SerializeGroups` step for group: parallel, `ByItemOrdinal`.
/// Serializes each [`BatchedProcessedPositionGroups`] to raw BAM bytes
/// ([`DecompressedBlock`]).
///
/// Delegates to the shared
/// [`build_serialize_processed_groups_step`][`crate::pipeline::steps::serialize_processed::build_serialize_processed_groups_step`]
/// used by both the standalone `fgumi group` chain and the runall
/// `--stop-after group` fused path.
///
/// Only included in the chain when the stage is
/// [`StagePosition::Terminal`][`crate::pipeline::chains::builder::StagePosition`];
/// for `Intermediate` the chain tail stays as [`BatchedProcessedPositionGroups`]
/// for the next stage's input.
///
/// `pub(crate)` — consumed only by [`ChainBuilder::add_group`].
#[allow(clippy::type_complexity)]
pub(crate) fn build_group_serialize_step(
    limit_bytes: u64,
    assign_tag_bytes: [u8; 2],
    progress_counter: Arc<std::sync::atomic::AtomicU64>,
) -> ProcessWithWorkerState<
    BatchedProcessedPositionGroups,
    DecompressedBlock,
    impl Fn(&mut SerializeState, BatchedProcessedPositionGroups) -> io::Result<DecompressedBlock>
    + Send
    + Sync
    + 'static,
    SerializeState,
    impl Fn() -> SerializeState + Send + Sync + 'static,
> {
    build_serialize_processed_groups_step(limit_bytes, assign_tag_bytes, progress_counter)
}

// ─────────────────────────────────────────────────────────────────────────────
// build_group_chain — thin delegate (Phase 3 T3a.7)
// ─────────────────────────────────────────────────────────────────────────────

/// Build a [`BuiltPipeline`] for a group-only chain.
///
/// Delegates to [`ChainBuilder`]. The full step sequence lives in
/// [`ChainBuilder`]'s `add_group` method in
/// [`crate::pipeline::chains::builder`].
///
/// The single-threaded fast path (`--threads` unset, BAM-only) is **not**
/// migrated here — it stays inline in `GroupReadsByUmi::execute`.
///
/// # Errors
///
/// Returns input-validation errors (sort-order mismatch, strategy/UMI
/// incompatibility) or any underlying pipeline construction error.
#[allow(clippy::needless_pass_by_value)]
pub fn build_group_chain(spec: ChainSpec) -> Result<BuiltPipeline> {
    use crate::pipeline::chains::{Stage, builder::StagePosition};
    let mut chain = ChainBuilder::new(&spec)?;
    chain.add_source()?;
    chain.add_stage(Stage::Group, StagePosition::Terminal)?;
    chain.add_sink()?;
    chain.build()
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::grouper::ProcessedPositionGroup;
    use fgumi_umi::MoleculeId;

    /// Build a [`ProcessedPositionGroup`] carrying `distinct_mi_count` and the
    /// given local-id templates (each as a `MoleculeId::Single`).
    fn make_group(distinct_mi_count: u64, local_ids: &[u64]) -> ProcessedPositionGroup {
        let templates = local_ids
            .iter()
            .map(|&id| {
                let mut template = Template::new(Vec::new());
                template.mi = MoleculeId::Single(id);
                template
            })
            .collect();
        ProcessedPositionGroup {
            templates,
            family_sizes: AHashMap::default(),
            filter_metrics: FilterMetrics::default(),
            input_record_count: 0,
            distinct_mi_count,
        }
    }

    #[test]
    fn assign_mi_offsets_advances_counter_and_shifts_local_ids() {
        // Two groups; the second must be offset past the first's block.
        let mut item = BatchedProcessedPositionGroups::new(
            0,
            vec![make_group(2, &[0, 1]), make_group(3, &[0, 1, 2])],
        );
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
        let mut item = BatchedProcessedPositionGroups::new(0, vec![make_group(u64::MAX, &[])]);
        let mut next_mi = 1;

        let err = assign_mi_offsets(&mut next_mi, &mut item)
            .expect_err("cumulative MI counter should overflow");
        assert_eq!(err.kind(), io::ErrorKind::Other);
        assert!(err.to_string().contains("MoleculeId offset overflow"));
    }
}
