//! Chain builder for `Stage::Group`.
//!
//! Phase 2 (T2.17) held the full ~500-LOC chain construction here.
//! Phase 3 (T3a.7) lifts that logic into
//! [`crate::pipeline::chains::builder::ChainBuilder`]; this module
//! now holds the group-specific types and step factories that the builder
//! imports: `GroupFinalizeHook` and the three step-factory functions used
//! by `ChainBuilder::add_group`.
//!
//! Mirrors the dedup (6423252), filter (8eead93), clip (8e95428) reference
//! migrations.

use std::io;
use std::sync::Arc;

use ahash::AHashMap;
use anyhow::Result;

use crate::assigner::Strategy;
use crate::commands::group::{
    GroupMetricsAccumulator, ParallelMinTemplates, assign_umi_groups_impl, build_grouping_metrics,
    create_umi_assigner, should_use_parallel, write_metrics_for_chain,
};
use crate::grouper::{ProcessedPositionGroup, RawPositionGroup, build_templates_from_records};
use crate::half_mapped::copies::CopyStore;
use crate::half_mapped::pair_ends::{MoleculeEnd, PairEndBuffer, PairEnds};
use crate::half_mapped::{EndReader, attach, is_half_mapped, split_off_placed};
use crate::logging::OperationTimer;
use crate::per_thread_accumulator::PerThreadAccumulator;
use crate::pipeline::chains::FinalizeHook;
use crate::pipeline::core::item::{HeapSize, Ordered};
use crate::pipeline::steps::group::position::{
    BatchedProcessedPositionGroups, BatchedRawPositionGroups,
};
use crate::pipeline::steps::process::{
    MiAssign, ProcessOrdered, ProcessWithWorkerState, mi_assign, process_ordered,
};
use crate::pipeline::steps::serialize_processed::{
    SerializeState, build_serialize_processed_groups_step, emit_templates_raw_with_mi,
};
use crate::pipeline::steps::types::DecompressedBlock;
use crate::template::Template;
use crate::template_filter::{TemplateFilterConfig, filter_template};
use fgumi_bam_io::{GroupKey, MemoryEstimate};
use fgumi_metrics::TemplateFilterCounts;
use fgumi_raw_bam::RawRecord;
use fgumi_umi::{IndexThreshold, MoleculeId};

// ─────────────────────────────────────────────────────────────────────────────
// GroupFinalizeHook
// ─────────────────────────────────────────────────────────────────────────────

/// Post-pipeline finalize hook for group. Reduces per-thread metric
/// accumulators, writes the optional metrics files (family-size histogram,
/// grouping-metrics, position-group-size histogram, and `--metrics PREFIX`
/// outputs), logs the UMI grouping summary banner, and calls
/// `timer.log_completion`.
///
/// `pub(crate)` so `ChainBuilder` can construct and register it in
/// `add_group`.
pub(crate) struct GroupFinalizeHook {
    pub(crate) accumulators: Arc<PerThreadAccumulator<GroupMetricsAccumulator>>,
    pub(crate) output_path: std::path::PathBuf,
    pub(crate) family_size_histogram: Option<std::path::PathBuf>,
    pub(crate) grouping_metrics: Option<std::path::PathBuf>,
    pub(crate) metrics_prefix: Option<std::path::PathBuf>,
    /// The half-mapped templates put into the molecules of fully mapped pairs.
    pub(crate) copies: Arc<CopyStore>,
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
            copies,
            timer,
        } = *self;

        // Reduce per-thread accumulators.
        let mut family_size_counter: AHashMap<usize, u64> = AHashMap::new();
        let mut position_group_size_counter: AHashMap<usize, u64> = AHashMap::new();
        let mut total_filter_counts = TemplateFilterCounts::new();
        for slot in accumulators.slots() {
            let acc = slot.lock();
            for (&size, &count) in &acc.family_sizes {
                *family_size_counter.entry(size).or_insert(0) += count;
            }
            for (&size, &count) in &acc.position_group_sizes {
                *position_group_size_counter.entry(size).or_insert(0) += count;
            }
            total_filter_counts.merge(&acc.filter_counts);
        }
        for change in copies.take_size_changes() {
            if let Some(count) = family_size_counter.get_mut(&change.from) {
                *count -= 1;
            }
            *family_size_counter.entry(change.to).or_insert(0) += 1;
        }
        family_size_counter.retain(|_, count| *count > 0);
        if !copies.is_empty() {
            info!("Half-mapped templates in the molecule of a fully mapped pair: {}", copies.len());
        }

        let metrics = build_grouping_metrics(&total_filter_counts, &family_size_counter);

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

/// Filters position groups and assigns their UMI groups, tracking per-thread family-size and
/// filter metrics.
///
/// `pub(crate)` — consumed only by `ChainBuilder::add_group`.
pub(crate) struct GroupProcessor {
    pub(crate) strategy: Strategy,
    pub(crate) effective_edits: u32,
    pub(crate) index_threshold: IndexThreshold,
    pub(crate) raw_tag: [u8; 2],
    pub(crate) parallel_group_min_templates: Option<ParallelMinTemplates>,
    pub(crate) num_threads: usize,
    pub(crate) filter_config: TemplateFilterConfig,
    pub(crate) accumulators: Arc<PerThreadAccumulator<GroupMetricsAccumulator>>,
    pub(crate) consensus_metrics:
        Option<Arc<crate::inline_metrics_collector::ConsensusMetricsCaptures>>,
    /// `Some` exactly when `consensus_metrics` is `Some` — both are built together
    /// in `add_group` only when a downstream consensus stage requested metrics.
    /// `None` on the default (metrics-off) path, so nothing is cloned or retained.
    pub(crate) header: Option<Arc<noodles::sam::Header>>,
    pub(crate) library_index: Option<Arc<fgumi_bam_io::LibraryIndex>>,
    /// `Some` exactly when `--realign-window` is set: the contig names of segment windows.
    /// `None` leaves every code path untouched.
    pub(crate) realign_contigs: Option<Arc<crate::realign::window_group::ContigNames>>,
}

/// A position group's templates that passed the filters.
#[derive(Debug)]
pub(crate) struct FilteredPositionGroup {
    key: GroupKey,
    templates: Vec<Template>,
    filter_counts: TemplateFilterCounts,
    input_record_count: u64,
}

#[derive(Debug)]
pub(crate) enum PendingPositionGroup {
    Processed(ProcessedPositionGroup),
    /// Its half-mapped templates may be copies of earlier fully mapped pairs, so its UMI groups
    /// have not been assigned.
    HalfMapped(FilteredPositionGroup),
}

/// A batch of position groups and the ends of its fully mapped pairs.
pub(crate) struct BatchedPendingPositionGroups {
    batch_serial: u64,
    groups: Vec<PendingPositionGroup>,
    pair_ends: PairEndBuffer<MoleculeEnd>,
}

impl HeapSize for BatchedPendingPositionGroups {
    fn heap_size(&self) -> usize {
        let group_size = |group: &PendingPositionGroup| match group {
            PendingPositionGroup::Processed(group) => group.estimate_heap_size(),
            PendingPositionGroup::HalfMapped(group) => {
                group.templates.iter().map(MemoryEstimate::estimate_heap_size).sum()
            }
        };
        self.groups.iter().map(group_size).sum::<usize>()
            + self.groups.capacity() * std::mem::size_of::<PendingPositionGroup>()
            + self.pair_ends.heap_size()
    }
}

impl Ordered for BatchedPendingPositionGroups {
    fn ordinal(&self) -> u64 {
        self.batch_serial
    }
}

impl GroupProcessor {
    /// Uniform cases build one assigner for the whole batch and reuse it
    /// across every group (preserving the per-batch reuse):
    /// `--allow-unmapped` always goes parallel, and the no-threshold
    /// default always stays sequential. The per-group
    /// `--parallel-group-min-templates` threshold must decide per group
    /// (template counts vary), so it builds in [`Self::process`]. See
    /// `should_use_parallel` / `create_umi_assigner` in `commands::group`.
    fn batch_assigner(&self) -> Option<Box<dyn crate::assigner::UmiAssigner>> {
        let allow_unmapped = self.filter_config.allow_unmapped;
        (allow_unmapped || self.parallel_group_min_templates.is_none()).then(|| {
            create_umi_assigner(
                self.strategy,
                self.effective_edits,
                self.index_threshold,
                self.num_threads,
                allow_unmapped,
            )
        })
    }

    fn end_reader(&self) -> EndReader {
        let by_segment =
            self.strategy.new_assigner(self.effective_edits).split_templates_by_pair_orientation();
        let TemplateFilterConfig { no_umi, min_umi_length, .. } = self.filter_config;
        EndReader::new(self.raw_tag, no_umi, by_segment, min_umi_length)
    }

    fn filter(&self, group: RawPositionGroup) -> io::Result<FilteredPositionGroup> {
        let key = group.group_key;
        let mut filter_counts = TemplateFilterCounts::new();
        let all_templates =
            build_templates_from_records(group.records).map_err(io::Error::other)?;
        let input_record_count = all_templates.iter().map(|t| t.read_count() as u64).sum();
        let templates = all_templates
            .into_iter()
            .filter(|t| filter_template(t, &self.filter_config, &mut filter_counts))
            .collect();
        Ok(FilteredPositionGroup { key, templates, filter_counts, input_record_count })
    }

    /// Assigns a filtered group's UMI groups.
    #[allow(clippy::too_many_lines)]
    fn process(
        &self,
        group: FilteredPositionGroup,
        batch_assigner: Option<&dyn crate::assigner::UmiAssigner>,
        end_reader: EndReader,
    ) -> io::Result<ProcessedPositionGroup> {
        let FilteredPositionGroup { key, mut templates, filter_counts, input_record_count } = group;

        if templates.is_empty() {
            self.accumulators.with_slot(|acc| {
                acc.record_group(AHashMap::new(), &filter_counts);
            });
            return Ok(ProcessedPositionGroup {
                templates: Vec::new(),
                family_sizes: AHashMap::new(),
                filter_counts,
                input_record_count,
                distinct_mi_count: 0,
                windows: None,
            });
        }

        // Pick the assigner for this group: reuse the batch-wide one for
        // the uniform cases, or build a per-group one when the threshold
        // flag is set (the decision depends on this group's size).
        let per_group_assigner;
        let assigner: &dyn crate::assigner::UmiAssigner = if let Some(batch) = batch_assigner {
            batch
        } else {
            let use_parallel = should_use_parallel(
                self.filter_config.allow_unmapped,
                templates.len(),
                self.strategy,
                self.parallel_group_min_templates.as_ref(),
            );
            per_group_assigner = create_umi_assigner(
                self.strategy,
                self.effective_edits,
                self.index_threshold,
                self.num_threads,
                use_parallel,
            );
            per_group_assigner.as_ref()
        };

        assigner.reset();

        self.assign(&mut templates, &key, assigner, end_reader).map_err(|e| {
            // A failed assignment means this position group's reads cannot be
            // grouped at all — an unparseable UMI for the chosen strategy, a
            // missing `RX` tag, mixed UMI lengths. Substituting an empty group
            // would silently drop every read in it and still exit 0, so fail
            // the run and let the user fix the input or the strategy. (Matches
            // the non-chain `Group::execute` path.)
            crate::umi::assign_umi_groups_error(&e)
        })?;

        let distinct_mi_count: u64 =
            templates.iter().filter_map(|t| t.mi.id()).max().map_or(0, |max_id| max_id + 1);

        templates.sort_by(|a, b| {
            let a_idx = a.mi.to_vec_index();
            let b_idx = b.mi.to_vec_index();
            a_idx.cmp(&b_idx).then_with(|| a.name().cmp(b.name()))
        });

        let windows = self
            .realign_contigs
            .as_ref()
            .map(|c| crate::realign::window_group::compute_segment_windows(&templates, c))
            .transpose()
            .map_err(io::Error::other)?;

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
        self.accumulators.with_slot(|acc| {
            acc.record_group(family_sizes.clone(), &filter_counts);
        });

        // T1 fused inline-metrics tap: when a downstream consensus
        // stage in this chain requested `--metrics`, record this
        // position group's templates into the shared accumulator now,
        // reusing the already-decoded `Template`s (no second BAM read).
        // `header`/`library_index` are `Some` exactly when
        // `consensus_metrics` is (built together in `add_group`), so the
        // three-way match both drives the tap and unwraps the metrics
        // header/library index without a separate expect.
        if let (Some(consensus_metrics), Some(header), Some(library_index)) =
            (&self.consensus_metrics, &self.header, &self.library_index)
        {
            let infos = crate::inline_metrics_collector::coordinate_group_from_processed_position(
                &templates,
                header,
                library_index,
            )
            .map_err(io::Error::other)?;
            consensus_metrics
                .accumulator
                .with_slot(|slot| {
                    slot.acc.record_coordinate_group(&infos, &consensus_metrics.intervals)
                })
                .map_err(io::Error::other)?;
        }

        Ok(ProcessedPositionGroup {
            templates,
            family_sizes,
            filter_counts,
            input_record_count,
            distinct_mi_count,
            windows,
        })
    }

    /// Assigns UMI groups among the templates other than the half-mapped ones sorted with a pair
    /// (`tc`), which then join the closest of them, or are grouped among themselves.
    fn assign(
        &self,
        templates: &mut Vec<Template>,
        key: &GroupKey,
        assigner: &dyn crate::assigner::UmiAssigner,
        end_reader: EndReader,
    ) -> Result<()> {
        let TemplateFilterConfig { no_umi, min_umi_length, .. } = self.filter_config;
        let assign = |templates: &mut [Template]| {
            assign_umi_groups_impl(templates, assigner, self.raw_tag, min_umi_length, no_umi)
        };
        let mut placed = split_off_placed(templates, key);
        assign(templates)?;
        if placed.is_empty() {
            return Ok(());
        }
        let placements = attach(templates, &mut placed, key, end_reader, self.effective_edits);
        let mut attached = placements.into_iter().map(|placement| placement.attached);
        let (mut joined, mut alone): (Vec<_>, Vec<_>) =
            placed.into_iter().partition(|_| attached.next().unwrap_or(false));
        assign(&mut alone)?;
        templates.append(&mut joined);
        templates.append(&mut alone);
        Ok(())
    }

    /// Takes the half-mapped templates with a molecule in `known` out of `group`, with the
    /// molecule.
    fn take_copies(
        &self,
        group: &mut FilteredPositionGroup,
        known: &mut PairEnds<MoleculeEnd>,
        end_reader: EndReader,
    ) -> Vec<Template> {
        let key = group.key;
        group
            .templates
            .extract_if(.., |template| {
                let closest = end_reader
                    .read_half(template, 0, &key)
                    .and_then(|half| known.closest(&half, self.effective_edits));
                if let Some((id, second)) = closest {
                    let is_r2 = template.r1().is_none_or(RawRecord::is_unmapped);
                    template.mi = match (self.strategy, is_r2 == second) {
                        (Strategy::Paired, true) => MoleculeId::PairedA(id),
                        (Strategy::Paired, false) => MoleculeId::PairedB(id),
                        _ => MoleculeId::Single(id),
                    };
                }
                closest.is_some()
            })
            .collect()
    }
}

/// Build the `GroupProcess` step: parallel, `ByItemOrdinal`. Filters
/// templates, assigns UMI groups, and tracks per-thread family-size and filter
/// metrics, except in position groups of half-mapped templates.
///
/// `pub(crate)` — consumed only by `ChainBuilder::add_group`.
#[allow(clippy::type_complexity)]
pub(crate) fn build_group_process_step(
    limit_bytes: u64,
    processor: Arc<GroupProcessor>,
) -> ProcessOrdered<
    BatchedRawPositionGroups,
    BatchedPendingPositionGroups,
    impl Fn(BatchedRawPositionGroups) -> io::Result<BatchedPendingPositionGroups>
    + Send
    + Sync
    + 'static,
> {
    let end_reader = processor.end_reader();
    process_ordered::<BatchedRawPositionGroups, BatchedPendingPositionGroups, _>(
        "GroupProcess",
        limit_bytes,
        move |item: BatchedRawPositionGroups| -> io::Result<BatchedPendingPositionGroups> {
            let BatchedRawPositionGroups { batch_serial, groups } = item;
            let batch_assigner = processor.batch_assigner();
            let mut pending = Vec::with_capacity(groups.len());
            let mut pair_ends = PairEndBuffer::default();
            for (index, group) in groups.into_iter().enumerate() {
                let group = processor.filter(group)?;
                let key = group.key;
                if !key.has_mate_position() && group.templates.iter().any(is_half_mapped) {
                    pending.push(PendingPositionGroup::HalfMapped(group));
                    continue;
                }
                let group = processor.process(group, batch_assigner.as_deref(), end_reader)?;
                pair_ends.add_group(end_reader, &key, index, &group.templates);
                pending.push(PendingPositionGroup::Processed(group));
            }
            Ok(BatchedPendingPositionGroups { batch_serial, groups: pending, pair_ends })
        },
    )
}

/// Build the `MiAssignGroups` step: serial, `ByItemOrdinal`. Assigns
/// monotonically increasing MI offsets to each batch of processed position
/// groups, and takes out the half-mapped templates that share a molecule with
/// fully mapped pairs.
///
/// Group uses [`BatchedProcessedPositionGroups`] (not the dedup-specific
/// `BatchedProcessedDedupGroups`), so this factory is distinct from
/// `chains::commands::dedup::build_mi_assign_step`.
///
/// `pub(crate)` — consumed only by `ChainBuilder::add_group`.
#[allow(clippy::type_complexity)]
pub(crate) fn build_group_mi_assign_step(
    limit_bytes: u64,
    processor: Arc<GroupProcessor>,
    copies: Arc<CopyStore>,
    assign_tag_bytes: [u8; 2],
) -> MiAssign<
    BatchedPendingPositionGroups,
    BatchedProcessedPositionGroups,
    impl Fn(&mut u64, BatchedPendingPositionGroups) -> io::Result<BatchedProcessedPositionGroups>
    + Send
    + Sync
    + 'static,
> {
    let end_reader = processor.end_reader();
    let known = parking_lot::Mutex::new(PairEnds::default());
    mi_assign::<BatchedPendingPositionGroups, BatchedProcessedPositionGroups, _>(
        "MiAssignGroups",
        limit_bytes,
        move |next_mi: &mut u64, item: BatchedPendingPositionGroups| {
            let BatchedPendingPositionGroups { batch_serial, groups, pair_ends } = item;
            let mut known = known.lock();
            let mut ends = pair_ends.reader();
            let mut assigner = None;
            let mut first_ids = Vec::with_capacity(groups.len());
            let mut processed = Vec::with_capacity(groups.len());
            for (index, group) in groups.into_iter().enumerate() {
                first_ids.push(*next_mi);
                let mut group = match group {
                    PendingPositionGroup::Processed(group) => group,
                    PendingPositionGroup::HalfMapped(mut group) => {
                        known.add(&mut ends, index, |index| first_ids[index]);
                        let found = processor.take_copies(&mut group, &mut known, end_reader);
                        store_copies(&copies, &found, assign_tag_bytes)?;
                        let assigner = assigner.get_or_insert_with(|| processor.batch_assigner());
                        processor.process(group, assigner.as_deref(), end_reader)?
                    }
                };
                assign_mi_offset(next_mi, &mut group)?;
                processed.push(group);
            }
            known.add(&mut ends, usize::MAX, |index| first_ids[index]);
            known.release(&pair_ends);
            copies.end_batch(*next_mi);
            Ok(BatchedProcessedPositionGroups::new(batch_serial, processed))
        },
    )
}

/// Adds the half-mapped templates `found` to `copies`, as they will be written.
fn store_copies(
    copies: &CopyStore,
    found: &[Template],
    assign_tag_bytes: [u8; 2],
) -> io::Result<()> {
    let (mut records, mut scratch, mut mi_buf) = (Vec::new(), Vec::new(), String::new());
    for template in found {
        records.clear();
        emit_templates_raw_with_mi(
            std::slice::from_ref(template),
            0,
            None,
            assign_tag_bytes,
            &mut scratch,
            &mut mi_buf,
            |bytes| {
                records.extend_from_slice(bytes);
                Ok(())
            },
        )?;
        if let Some(molecule) = template.mi.id() {
            copies.add(molecule, &records)?;
        }
    }
    Ok(())
}

/// Assign global `MoleculeId` offsets to every template in `processed`, advancing
/// the running `next_mi` counter by the group's `distinct_mi_count`.
///
/// The group reserves a contiguous block `[base, base + distinct_mi_count)`
/// from the global counter and shifts its templates' local MI ids by `base`.
///
/// # Errors
///
/// Returns an [`io::Error`] if the cumulative MI counter would exceed
/// [`u64::MAX`].
fn assign_mi_offset(next_mi: &mut u64, processed: &mut ProcessedPositionGroup) -> io::Result<()> {
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
    if let Some(windows) = &mut processed.windows {
        windows.rebase(base);
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
/// `pub(crate)` — consumed only by `ChainBuilder::add_group`.
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

#[cfg(test)]
mod tests {
    use super::*;
    use crate::grouper::ProcessedPositionGroup;

    /// Build a [`ProcessedPositionGroup`] carrying `distinct_mi_count` and the
    /// given local-id templates (each as a `MoleculeId::Single`).
    fn make_group(distinct_mi_count: u64, local_ids: &[u64]) -> ProcessedPositionGroup {
        let templates = local_ids
            .iter()
            .map(|&id| {
                let mut template = Template::new();
                template.mi = MoleculeId::Single(id);
                template
            })
            .collect();
        ProcessedPositionGroup {
            templates,
            family_sizes: AHashMap::default(),
            filter_counts: TemplateFilterCounts::default(),
            input_record_count: 0,
            distinct_mi_count,
            windows: None,
        }
    }

    #[test]
    fn assign_mi_offset_advances_counter_and_shifts_local_ids() {
        // Two groups; the second must be offset past the first's block.
        let mut item = BatchedProcessedPositionGroups::new(
            0,
            vec![make_group(2, &[0, 1]), make_group(3, &[0, 1, 2])],
        );
        let mut next_mi = 10;

        for group in &mut item.groups {
            assign_mi_offset(&mut next_mi, group).expect("no overflow");
        }

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
    fn assign_mi_offset_errors_on_counter_overflow() {
        // A non-zero starting counter plus a u64::MAX block overflows.
        let mut next_mi = 1;

        let err = assign_mi_offset(&mut next_mi, &mut make_group(u64::MAX, &[]))
            .expect_err("cumulative MI counter should overflow");
        assert_eq!(err.kind(), io::ErrorKind::Other);
        assert!(err.to_string().contains("MoleculeId offset overflow"));
    }
}
