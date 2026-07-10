//! Chain builder for `Stage::Duplex`.
//!
//! Phase 2 (T2.19) held the full ~550-LOC chain construction here.
//! Phase 3 (T3a.9) lifts that logic into
//! [`crate::pipeline::chains::builder::ChainBuilder`]; this module
//! now holds the duplex-specific types and step factory that the builder
//! imports: [`DuplexFinalizeHook`] and the `build_duplex_consensus_step_with_rejects` / `build_duplex_consensus_step_kept_only` factories, used
//! by `ChainBuilder::add_duplex`.
//!
//! [`build_duplex_chain`] shrinks to a ~10-line delegate.
//!
//! Mirrors simplex (5cce86c). Duplex's chain is structurally the same with a
//! different consensus caller, a record filter, and an MI-tag transform that
//! strips `/A`/`/B` suffixes.
//!
//! This is the sole standalone execution path — `Duplex::execute` always routes
//! here, and `--threads` unset resolves to a one-worker chain (no separate
//! single-threaded path).

use std::io;
use std::sync::Arc;
use std::sync::atomic::{AtomicU64, Ordering};

use anyhow::Result;
use log::info;

use crate::commands::common::MethylationRef;
use crate::commands::consensus_runner::{ConsensusStatsOps, log_overlapping_stats};
use crate::consensus_caller::{ConsensusCaller, ConsensusCallingStats, ConsensusOutput};
use crate::duplex_consensus_caller::DuplexConsensusCaller;
use crate::logging::OperationTimer;
use crate::mi_group::MiGroup;
use crate::overlapping_consensus::{
    AgreementStrategy, CorrectionStats, DisagreementStrategy, OverlappingBasesConsensusCaller,
    apply_overlapping_consensus,
};
use crate::per_thread_accumulator::PerThreadAccumulator;
use crate::pipeline::chains::builder::ChainBuilder;
use crate::pipeline::chains::{BuiltPipeline, ChainSpec, FinalizeHook};
use crate::pipeline::core::outputs::OrderedBytesTuple2;
use crate::pipeline::core::step::Step;
use crate::pipeline::steps::group::mi::BatchedMiGroups;
use crate::pipeline::steps::process::{
    Process2Output, ProcessWithWorkerState, process_with_worker_state, process2_with_worker_state,
};
use crate::pipeline::steps::types::DecompressedBlock;

// ─────────────────────────────────────────────────────────────────────────────
// CollectedDuplexMetrics
// ─────────────────────────────────────────────────────────────────────────────

/// Per-thread accumulator for duplex consensus metrics (mirrors
/// `commands::duplex::CollectedDuplexMetrics`).
///
/// Merged into final aggregates after the pipeline completes; one instance
/// per worker slot (see [`PerThreadAccumulator`]).
///
/// `pub(crate)` so [`ChainBuilder`] can construct it in `add_duplex`.
#[derive(Default)]
pub(crate) struct CollectedDuplexMetrics {
    /// Consensus calling statistics.
    pub(crate) stats: ConsensusCallingStats,
    /// Overlapping consensus stats (if enabled).
    pub(crate) overlapping_stats: Option<CorrectionStats>,
    /// Number of MI groups processed.
    pub(crate) groups_processed: u64,
}

// ─────────────────────────────────────────────────────────────────────────────
// DuplexState
// ─────────────────────────────────────────────────────────────────────────────

/// Per-worker state for the duplex consensus step.
///
/// Defined at module level (not inside a function body) to avoid the
/// `clippy::items_after_statements` lint that fires when a struct is defined
/// after executable statements in a function body.
pub(crate) struct DuplexState {
    pub(crate) caller: DuplexConsensusCaller,
    pub(crate) overlapping: Option<OverlappingBasesConsensusCaller>,
}

impl crate::pipeline::core::item::HeapSize for DuplexState {}

// ─────────────────────────────────────────────────────────────────────────────
// DuplexFinalizeHook
// ─────────────────────────────────────────────────────────────────────────────

/// Post-pipeline finalize hook for duplex. Reduces per-thread metrics,
/// writes the optional stats file, logs the overlapping-consensus stats
/// (if enabled), logs the summary banner, finalizes the rejects writer,
/// and calls `timer.log_completion`.
///
/// `pub(crate)` so [`ChainBuilder`] can construct and register it in
/// `add_duplex`.
pub(crate) struct DuplexFinalizeHook {
    pub(crate) accumulators: Arc<PerThreadAccumulator<CollectedDuplexMetrics>>,
    pub(crate) stats_path: Option<std::path::PathBuf>,
    pub(crate) overlapping_enabled: bool,
    pub(crate) timer: OperationTimer,
}

impl FinalizeHook for DuplexFinalizeHook {
    fn finalize(self: Box<Self>) -> Result<()> {
        let DuplexFinalizeHook { accumulators, stats_path, overlapping_enabled, timer } = *self;

        // Reduce per-thread accumulators.
        let mut total_groups = 0u64;
        let mut merged_stats = ConsensusCallingStats::new();
        let mut merged_overlapping_stats = CorrectionStats::new();

        for slot in accumulators.slots() {
            let m = slot.lock();
            total_groups += m.groups_processed;
            merged_stats.merge(&m.stats);
            if let Some(ref ocs) = m.overlapping_stats {
                merged_overlapping_stats.merge(ocs);
            }
        }

        if overlapping_enabled {
            log_overlapping_stats(&merged_overlapping_stats);
        }

        info!("Duplex consensus calling complete");
        info!("Total MI groups processed: {total_groups}");

        let metrics = merged_stats.to_metrics();
        let consensus_count = metrics.consensus_reads;
        crate::logging::log_consensus_summary(&metrics);

        if let Some(ref stats_path) = stats_path {
            use fgoxide::io::DelimFile;
            let kv_metrics = metrics.to_kv_metrics();
            DelimFile::default().write_tsv(stats_path, kv_metrics).map_err(|e| {
                anyhow::anyhow!("Failed to write statistics: {}: {e}", stats_path.display())
            })?;
            info!("Wrote statistics to: {}", stats_path.display());
        }

        info!("Wrote {consensus_count} duplex consensus reads");

        // Rejects are emitted on the consensus step's second output branch
        // (fan-out) and finalized by the rejects branch's WriteBgzfFile sink.

        timer.log_completion(consensus_count);

        Ok(())
    }
}

// ─────────────────────────────────────────────────────────────────────────────
// Step factory — extracted from build_duplex_chain (T3a.9).
//
// Returns the concrete `ProcessWithWorkerState` type directly (the same
// pattern used by the other step factories in this module family). Returning
// `impl Step<...>` is blocked here because the closure types embed in
// opaque-return position and cannot name themselves.
// ─────────────────────────────────────────────────────────────────────────────

/// Captures passed into [`build_duplex_consensus_step_with_rejects`] / [`build_duplex_consensus_step_kept_only`] from `add_duplex`.
///
/// Bundles all the cloned scalars and Arcs the closure needs so `add_duplex`
/// can prepare them once and hand them off cleanly.
///
/// `pub(crate)` — consumed only by [`ChainBuilder::add_duplex`] and
/// [`build_duplex_consensus_step_with_rejects`] / [`build_duplex_consensus_step_kept_only`].
#[allow(clippy::struct_excessive_bools)]
pub(crate) struct DuplexConsensusCaptures {
    pub(crate) track_rejects: bool,
    pub(crate) overlapping_enabled: bool,
    pub(crate) methylation_ref: MethylationRef,
    pub(crate) methylation_mode: fgumi_consensus::MethylationMode,
    pub(crate) read_name_prefix: String,
    pub(crate) read_group_id: String,
    pub(crate) min_reads: Vec<usize>,
    pub(crate) min_input_base_quality: u8,
    pub(crate) output_per_base_tags: bool,
    pub(crate) trim: bool,
    pub(crate) max_reads_per_strand: Option<usize>,
    pub(crate) error_rate_pre_umi: u8,
    pub(crate) error_rate_post_umi: u8,
    pub(crate) cell_tag: noodles::sam::alignment::record::data::field::Tag,
    pub(crate) accumulators: Arc<PerThreadAccumulator<CollectedDuplexMetrics>>,
    pub(crate) progress: Arc<AtomicU64>,
}

/// Per-worker init: build the `DuplexConsensusCaller` (+ optional overlapping
/// caller) once, reused across batches. Shared by both step variants.
///
/// # Panics
///
/// Panics if `DuplexConsensusCaller::new` fails (only on invalid `min_reads`,
/// validated before the pipeline is built) — the init closure type
/// `Fn() -> DuplexState` cannot propagate a `Result`.
#[allow(clippy::too_many_arguments, clippy::fn_params_excessive_bools)]
fn make_duplex_consensus_init(
    read_name_prefix: String,
    read_group_id: String,
    min_reads: Vec<usize>,
    min_input_base_quality: u8,
    output_per_base_tags: bool,
    trim: bool,
    max_reads_per_strand: Option<usize>,
    cell_tag: noodles::sam::alignment::record::data::field::Tag,
    track_rejects: bool,
    error_rate_pre_umi: u8,
    error_rate_post_umi: u8,
    methylation_ref: MethylationRef,
    methylation_mode: fgumi_consensus::MethylationMode,
    overlapping_enabled: bool,
) -> impl Fn() -> DuplexState + Send + Sync + 'static {
    move || {
        let mut caller = DuplexConsensusCaller::new(
            read_name_prefix.clone(),
            read_group_id.clone(),
            min_reads.clone(),
            min_input_base_quality,
            output_per_base_tags,
            trim,
            max_reads_per_strand,
            Some(cell_tag),
            track_rejects,
            error_rate_pre_umi,
            error_rate_post_umi,
        )
        .expect("DuplexConsensusCaller::new failed during worker init");
        if let Some((ref reference, ref ref_names)) = methylation_ref {
            caller.set_reference(Arc::clone(reference), Arc::clone(ref_names), methylation_mode);
        }
        let overlapping = if overlapping_enabled {
            Some(OverlappingBasesConsensusCaller::new(
                AgreementStrategy::Consensus,
                DisagreementStrategy::Consensus,
            ))
        } else {
            None
        };
        DuplexState { caller, overlapping }
    }
}

/// Per-batch duplex consensus body, shared by both step variants.
///
/// Returns the consensus `DecompressedBlock` (branch 0) and, when
/// `track_rejects` is set and any record was rejected, a rejects
/// `DecompressedBlock` (branch 1) of `[len][record]`-framed raw-input records
/// in input order — the PR #332 fan-out contract (rejects flow through the
/// ordered serialize/compress stages, not a mutex side-channel; the rejects
/// writer is configured with the input header by `add_duplex`).
fn run_duplex_consensus_batch(
    state: &mut DuplexState,
    item: BatchedMiGroups,
    track_rejects: bool,
    overlapping_enabled: bool,
    accumulators: &Arc<PerThreadAccumulator<CollectedDuplexMetrics>>,
    progress: &Arc<AtomicU64>,
) -> io::Result<(DecompressedBlock, Option<DecompressedBlock>)> {
    let BatchedMiGroups { batch_serial, groups } = item;
    let groups_count = groups.len() as u64;

    let mut all_output = ConsensusOutput::default();
    let mut batch_stats = ConsensusCallingStats::new();
    let mut batch_overlapping = CorrectionStats::new();
    let mut rejects_bytes: Vec<u8> = Vec::new();

    let mut total_input_records: u64 = 0;
    for MiGroup { mi, records: mut group_reads } in groups {
        state.caller.clear();
        total_input_records += group_reads.len() as u64;

        if let Some(ref mut oc) = state.overlapping {
            if crate::commands::duplex::has_both_strands_raw(&group_reads) {
                oc.reset_stats();
                apply_overlapping_consensus(&mut group_reads, oc).map_err(|e| {
                    io::Error::other(format!("Overlapping consensus error for MI {mi}: {e}"))
                })?;
                batch_overlapping.merge(oc.stats());
            }
        }

        let group_output = state
            .caller
            .consensus_reads(group_reads)
            .map_err(|e| io::Error::other(format!("Duplex consensus error for MI {mi}: {e}")))?;
        all_output.merge(group_output);
        batch_stats.merge(&state.caller.statistics());
        if track_rejects {
            for raw in &state.caller.take_rejected_reads() {
                super::append_framed_bytes(&mut rejects_bytes, raw)?;
            }
        }
    }

    // Merge per-batch metrics into this worker's slot.
    accumulators.with_slot(|m| {
        m.stats.merge(&batch_stats);
        if overlapping_enabled {
            m.overlapping_stats.get_or_insert_with(CorrectionStats::new).merge(&batch_overlapping);
        }
        m.groups_processed += groups_count;
    });

    // Progress logging at million-record boundaries.
    let prev = progress.fetch_add(total_input_records, Ordering::Relaxed);
    if (prev + total_input_records) / 1_000_000 > prev / 1_000_000 {
        info!("Processed {} records", prev + total_input_records);
    }

    let consensus = DecompressedBlock { batch_serial, bytes: all_output.data };
    let rejects = if rejects_bytes.is_empty() {
        None
    } else {
        Some(DecompressedBlock { batch_serial, bytes: rejects_bytes })
    };
    Ok((consensus, rejects))
}

/// Build the 2-output `DuplexConsensus` step (used when `--rejects` is set):
/// branch 0 = consensus, branch 1 = rejects.
///
/// `pub(crate)` — consumed only by [`ChainBuilder::add_duplex`].
pub(crate) fn build_duplex_consensus_step_with_rejects(
    limit_bytes: u64,
    cap: DuplexConsensusCaptures,
) -> impl Step<Input = BatchedMiGroups, Outputs = OrderedBytesTuple2<DecompressedBlock, DecompressedBlock>>
{
    let DuplexConsensusCaptures {
        track_rejects,
        overlapping_enabled,
        methylation_ref,
        methylation_mode,
        read_name_prefix,
        read_group_id,
        min_reads,
        min_input_base_quality,
        output_per_base_tags,
        trim,
        max_reads_per_strand,
        error_rate_pre_umi,
        error_rate_post_umi,
        cell_tag,
        accumulators,
        progress,
    } = cap;

    let init = make_duplex_consensus_init(
        read_name_prefix,
        read_group_id,
        min_reads,
        min_input_base_quality,
        output_per_base_tags,
        trim,
        max_reads_per_strand,
        cell_tag,
        track_rejects,
        error_rate_pre_umi,
        error_rate_post_umi,
        methylation_ref,
        methylation_mode,
        overlapping_enabled,
    );
    let body = move |state: &mut DuplexState,
                     item: BatchedMiGroups|
          -> io::Result<Process2Output<DecompressedBlock, DecompressedBlock>> {
        let (consensus, rejects) = run_duplex_consensus_batch(
            state,
            item,
            track_rejects,
            overlapping_enabled,
            &accumulators,
            &progress,
        )?;
        if let Some(r) = rejects {
            Ok(Process2Output::both(consensus, r))
        } else {
            // X5-001: an all-clean batch must still emit a (zero-byte) rejects
            // block so the rejects branch's `ByItemOrdinal` reorder stage sees a
            // dense serial sequence; pushing nothing (`only_a`) leaves a gap at
            // this serial that wedges the rejects sink. The empty `Vec` does not
            // allocate and produces no physical BGZF block.
            let batch_serial = consensus.batch_serial;
            Ok(Process2Output::both(
                consensus,
                DecompressedBlock { batch_serial, bytes: Vec::new() },
            ))
        }
    };

    process2_with_worker_state::<
        BatchedMiGroups,
        DecompressedBlock,
        DecompressedBlock,
        DuplexState,
        _,
        _,
    >("DuplexConsensus", limit_bytes, limit_bytes, init, body)
}

/// Build the 1-output kept-only `DuplexConsensus` step (used when `--rejects`
/// is unset). The framework has no public discard sink, so omitting the rejects
/// branch requires this single-output variant.
///
/// `pub(crate)` — consumed only by [`ChainBuilder::add_duplex`].
#[allow(clippy::type_complexity)]
pub(crate) fn build_duplex_consensus_step_kept_only(
    limit_bytes: u64,
    cap: DuplexConsensusCaptures,
) -> ProcessWithWorkerState<
    BatchedMiGroups,
    DecompressedBlock,
    impl Fn(&mut DuplexState, BatchedMiGroups) -> io::Result<DecompressedBlock> + Send + Sync + 'static,
    DuplexState,
    impl Fn() -> DuplexState + Send + Sync + 'static,
> {
    let DuplexConsensusCaptures {
        track_rejects,
        overlapping_enabled,
        methylation_ref,
        methylation_mode,
        read_name_prefix,
        read_group_id,
        min_reads,
        min_input_base_quality,
        output_per_base_tags,
        trim,
        max_reads_per_strand,
        error_rate_pre_umi,
        error_rate_post_umi,
        cell_tag,
        accumulators,
        progress,
    } = cap;

    let init = make_duplex_consensus_init(
        read_name_prefix,
        read_group_id,
        min_reads,
        min_input_base_quality,
        output_per_base_tags,
        trim,
        max_reads_per_strand,
        cell_tag,
        track_rejects,
        error_rate_pre_umi,
        error_rate_post_umi,
        methylation_ref,
        methylation_mode,
        overlapping_enabled,
    );
    let body =
        move |state: &mut DuplexState, item: BatchedMiGroups| -> io::Result<DecompressedBlock> {
            let (consensus, _rejects) = run_duplex_consensus_batch(
                state,
                item,
                track_rejects,
                overlapping_enabled,
                &accumulators,
                &progress,
            )?;
            Ok(consensus)
        };

    process_with_worker_state::<BatchedMiGroups, DecompressedBlock, _, DuplexState, _>(
        "DuplexConsensus",
        limit_bytes,
        init,
        body,
    )
}

// ─────────────────────────────────────────────────────────────────────────────
// build_duplex_chain — thin delegate (Phase 3 T3a.9)
// ─────────────────────────────────────────────────────────────────────────────

/// Build a [`BuiltPipeline`] for a duplex-only chain.
///
/// Delegates to [`ChainBuilder`]. The full step sequence lives in
/// [`ChainBuilder`]'s `add_duplex` method in
/// [`crate::pipeline::chains::builder`].
///
/// This is the sole standalone execution path — `Duplex::execute` always routes
/// here, and `--threads` unset resolves to a one-worker chain.
///
/// # Errors
///
/// Returns input-validation errors (missing reference for methylation mode,
/// etc.) or any underlying pipeline construction error.
#[allow(clippy::needless_pass_by_value)]
pub fn build_duplex_chain(spec: ChainSpec) -> Result<BuiltPipeline> {
    use crate::pipeline::chains::{Stage, builder::StagePosition};
    let mut chain = ChainBuilder::new(&spec)?;
    chain.add_source()?;
    chain.add_stage(Stage::Duplex, StagePosition::Terminal)?;
    chain.add_sink()?;
    chain.build()
}
