//! Chain builder for `Stage::Simplex`.
//!
//! Phase 2 (T2.18) held the full ~500-LOC chain construction here.
//! Phase 3 (T3a.8) lifts that logic into
//! [`crate::pipeline::chains::builder::ChainBuilder`]; this module
//! now holds the simplex-specific types and step factory that the builder
//! imports: `SimplexFinalizeHook` and the `build_simplex_consensus_step_with_rejects` / `build_simplex_consensus_step_kept_only` factories, used
//! by `ChainBuilder::add_simplex`.
//!
//! Mirrors the dedup (6423252), filter (8eead93), clip (8e95428), group
//! (d371d84) reference migrations.
//!
//! This module supplies the chain-builder pieces for the simplex stage; the
//! chain is constructed via `ChainBuilder` /
//! [`crate::pipeline::chains::build::build_for`]. It is the sole path for
//! `Simplex::execute` (via `execute_chain`), with or without `--threads` —
//! absent `--threads` runs the chain at a single worker.

use std::io;
use std::sync::Arc;
use std::sync::atomic::{AtomicU64, Ordering};

use anyhow::Result;
use log::info;

use crate::commands::common::MethylationRef;
use crate::commands::consensus_runner::{ConsensusStatsOps, log_overlapping_stats};
use crate::consensus_caller::{
    ConsensusCaller, ConsensusCallingStats, ConsensusOutput, RejectionReason,
};
use crate::inline_metrics_collector::push_mi_group_entries;
use crate::logging::OperationTimer;
use crate::mi_group::MiGroup;
use crate::overlapping_consensus::{
    AgreementStrategy, CorrectionStats, DisagreementStrategy, OverlappingBasesConsensusCaller,
    apply_overlapping_consensus,
};
use crate::per_thread_accumulator::PerThreadAccumulator;
use crate::pipeline::chains::FinalizeHook;
use crate::pipeline::core::outputs::OrderedBytesTuple2;
use crate::pipeline::core::step::Step;
use crate::pipeline::steps::group::mi::BatchedMiGroups;
use crate::pipeline::steps::process::{
    Process2Output, ProcessWithWorkerState, process_with_worker_state, process2_with_worker_state,
};
use crate::pipeline::steps::types::DecompressedBlock;
use crate::vanilla_consensus_caller::{VanillaUmiConsensusCaller, VanillaUmiConsensusOptions};

// ─────────────────────────────────────────────────────────────────────────────
// CollectedSimplexMetrics
// ─────────────────────────────────────────────────────────────────────────────

/// Per-thread accumulator for simplex consensus metrics.
///
/// Merged into final aggregates after the pipeline completes; one instance
/// per worker slot (see [`PerThreadAccumulator`]).
///
/// `pub(crate)` so `ChainBuilder` can construct it in `add_simplex`.
#[derive(Default)]
pub(crate) struct CollectedSimplexMetrics {
    /// Consensus calling statistics.
    pub(crate) stats: ConsensusCallingStats,
    /// Overlapping consensus stats (if enabled).
    pub(crate) overlapping_stats: Option<CorrectionStats>,
    /// Number of MI groups processed.
    pub(crate) groups_processed: u64,
}

// ─────────────────────────────────────────────────────────────────────────────
// ConsensusState
// ─────────────────────────────────────────────────────────────────────────────

/// Per-worker state for the simplex consensus step.
///
/// Defined at module level (not inside a function body) to avoid the
/// `clippy::items_after_statements` lint that fires when a struct is defined
/// after executable statements in a function body.
pub(crate) struct ConsensusState {
    pub(crate) caller: VanillaUmiConsensusCaller,
    pub(crate) overlapping: Option<OverlappingBasesConsensusCaller>,
    /// Input header + library index, threaded in for the inline-metrics T2
    /// path so the metrics-on batch body can build `TemplateInfo`s via
    /// `push_mi_group_entries` (Task 11). Present on every worker regardless of
    /// whether metrics are on — only the metrics-on batch body reads them.
    pub(crate) header: Arc<noodles::sam::Header>,
    pub(crate) library_index: Arc<fgumi_bam_io::LibraryIndex>,
}

impl crate::pipeline::core::item::HeapSize for ConsensusState {}

// ─────────────────────────────────────────────────────────────────────────────
// SimplexFinalizeHook
// ─────────────────────────────────────────────────────────────────────────────

/// Post-pipeline finalize hook for simplex. Reduces per-thread metrics,
/// writes the optional stats file, logs the overlapping-consensus stats
/// (if enabled), logs the summary banner, finalizes the rejects writer,
/// and calls `timer.log_completion`.
///
/// `pub(crate)` so `ChainBuilder` can construct and register it in
/// `add_simplex`.
pub(crate) struct SimplexFinalizeHook {
    pub(crate) accumulators: Arc<PerThreadAccumulator<CollectedSimplexMetrics>>,
    pub(crate) stats_path: Option<std::path::PathBuf>,
    pub(crate) overlapping_enabled: bool,
    pub(crate) timer: OperationTimer,
}

impl FinalizeHook for SimplexFinalizeHook {
    fn finalize(self: Box<Self>) -> Result<()> {
        let SimplexFinalizeHook { accumulators, stats_path, overlapping_enabled, timer } = *self;

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

        info!("Consensus calling complete");
        info!("Total MI groups processed: {total_groups}");

        let metrics = merged_stats.to_metrics();
        let consensus_count = metrics.consensus_reads;
        crate::logging::log_consensus_summary(&metrics);

        if let Some(ref stats_path) = stats_path {
            use fgoxide::io::DelimFile;
            let kv_metrics = metrics.to_kv_metrics(fgumi_metrics::ConsensusCallerKind::Vanilla);
            DelimFile::default().write_tsv(stats_path, kv_metrics).map_err(|e| {
                anyhow::anyhow!("Failed to write statistics: {}: {e}", stats_path.display())
            })?;
            info!("Wrote statistics to: {}", stats_path.display());
        }

        info!("Wrote {consensus_count} consensus reads");

        // Rejects are now emitted on the consensus step's second output branch
        // (fan-out) and finalized by the rejects branch's own WriteBgzfFile sink,
        // so there is no rejects writer to finish here.

        timer.log_completion(consensus_count);

        Ok(())
    }
}

// ─────────────────────────────────────────────────────────────────────────────
// Step factory — extracted from build_simplex_chain (T3a.8).
//
// Returns the concrete `ProcessWithWorkerState` type directly (the same
// pattern used by the other step factories in this module family). Returning
// `impl Step<...>` is blocked here because the closure types embed in
// opaque-return position and cannot name themselves.
// ─────────────────────────────────────────────────────────────────────────────

/// Maps `SimplexOptions` (the resolved options-bag entry `add_simplex`
/// receives) into the `VanillaUmiConsensusOptions` the per-worker consensus
/// caller is built from.
///
/// Extracted out of `ChainBuilder::add_simplex` so this mapping is
/// unit-testable on its own: the `test_simplex_chain_matches_single_threaded`
/// parity tests compare two chain runs (single- vs multi-worker) that BOTH go
/// through this same mapping, so a field dropped or swapped inside it would
/// make the two sides agree with each other and still pass. See
/// `simplex_consensus_options_carries_every_tuning_flag` below, which mirrors
/// `to_simplex_options_carries_every_tuning_flag`'s non-default-everywhere
/// discipline for the other half of the CLI-args -> `SimplexOptions` ->
/// `VanillaUmiConsensusOptions` pipeline.
pub(crate) fn simplex_consensus_options(
    simplex: &crate::commands::simplex::SimplexOptions,
    cell_tag: noodles::sam::alignment::record::data::field::Tag,
) -> VanillaUmiConsensusOptions {
    let consensus = simplex.consensus();
    VanillaUmiConsensusOptions {
        tag: "MI".to_string(),
        error_rate_pre_umi: consensus.error_rate_pre_umi,
        error_rate_post_umi: consensus.error_rate_post_umi,
        min_input_base_quality: consensus.min_input_base_quality,
        min_reads: simplex.min_reads,
        max_reads: simplex.max_reads,
        produce_per_base_tags: consensus.output_per_base_tags,
        trim: consensus.trim,
        min_consensus_base_quality: consensus.min_consensus_base_quality,
        cell_tag: Some(cell_tag),
        methylation_mode: simplex.methylation_mode,
        tie_rule: simplex.tie_rule,
    }
}

/// Captures passed into [`build_simplex_consensus_step_with_rejects`] / [`build_simplex_consensus_step_kept_only`] from `add_simplex`.
///
/// Bundles all the cloned scalars and Arcs the closure needs so `add_simplex`
/// can prepare them once and hand them off cleanly.
///
/// `pub(crate)` — consumed only by `ChainBuilder::add_simplex` and
/// [`build_simplex_consensus_step_with_rejects`] / [`build_simplex_consensus_step_kept_only`].
pub(crate) struct SimplexConsensusCaptures {
    pub(crate) track_rejects: bool,
    pub(crate) overlapping_enabled: bool,
    pub(crate) consensus_options: VanillaUmiConsensusOptions,
    pub(crate) read_name_prefix: String,
    pub(crate) read_group_id: String,
    pub(crate) methylation_ref: MethylationRef,
    pub(crate) accumulators: Arc<PerThreadAccumulator<CollectedSimplexMetrics>>,
    pub(crate) min_reads: usize,
    pub(crate) progress: Arc<AtomicU64>,
    /// Threaded onto every `ConsensusState` (see that struct); the metrics-on
    /// step variants read them, the metrics-off variants ignore them.
    pub(crate) header: Arc<noodles::sam::Header>,
    pub(crate) library_index: Arc<fgumi_bam_io::LibraryIndex>,
    /// QC captures for the metrics-ON T2 path; `None` on metrics-off/T1
    /// builds.
    pub(crate) qc_metrics: Option<Arc<crate::inline_metrics_collector::ConsensusMetricsCaptures>>,
}

/// Per-worker init: build the `VanillaUmiConsensusCaller` (+ optional
/// overlapping caller) once, reused across batches. Shared by both step
/// variants.
#[allow(clippy::too_many_arguments)]
fn make_simplex_consensus_init(
    read_name_prefix: String,
    read_group_id: String,
    consensus_options: VanillaUmiConsensusOptions,
    methylation_ref: MethylationRef,
    track_rejects: bool,
    overlapping_enabled: bool,
    header: Arc<noodles::sam::Header>,
    library_index: Arc<fgumi_bam_io::LibraryIndex>,
) -> impl Fn() -> ConsensusState + Send + Sync + 'static {
    move || {
        let mut caller = VanillaUmiConsensusCaller::new_with_rejects_tracking(
            read_name_prefix.clone(),
            read_group_id.clone(),
            consensus_options.clone(),
            track_rejects,
        );
        if let Some((ref reference, ref ref_names)) = methylation_ref {
            caller.set_reference(Arc::clone(reference), Arc::clone(ref_names));
        }
        let overlapping = if overlapping_enabled {
            Some(OverlappingBasesConsensusCaller::new(
                AgreementStrategy::Consensus,
                DisagreementStrategy::Consensus,
            ))
        } else {
            None
        };
        ConsensusState {
            caller,
            overlapping,
            header: Arc::clone(&header),
            library_index: Arc::clone(&library_index),
        }
    }
}

/// Per-batch simplex consensus body, shared by both step variants.
///
/// Returns the consensus `DecompressedBlock` (branch 0) and, when
/// `track_rejects` is set and at least one record was rejected, a rejects
/// `DecompressedBlock` (branch 1) of `[len][record]`-framed raw-input records
/// in input order. Per the PR #332 contract the rejects flow out through the
/// ordered serialize/compress stages (input order) rather than a mutex
/// side-channel (mutex-acquisition order); the rejects writer is configured
/// with the input header by `add_simplex`.
fn run_simplex_consensus_batch(
    state: &mut ConsensusState,
    item: BatchedMiGroups,
    track_rejects: bool,
    overlapping_enabled: bool,
    min_reads: usize,
    accumulators: &Arc<PerThreadAccumulator<CollectedSimplexMetrics>>,
    progress: &Arc<AtomicU64>,
) -> io::Result<(DecompressedBlock, Option<DecompressedBlock>)> {
    let BatchedMiGroups { batch_serial, groups } = item;
    let groups_count = groups.len() as u64;

    let mut all_output = ConsensusOutput::default();
    let mut batch_stats = ConsensusCallingStats::new();
    let mut batch_overlapping = CorrectionStats::new();
    let mut rejects_bytes: Vec<u8> = Vec::new();

    let mut total_input_records: u64 = 0;
    for MiGroup { mi, records: mut raw_records } in groups {
        state.caller.clear();
        total_input_records += raw_records.len() as u64;

        if raw_records.len() < min_reads {
            batch_stats.record_input(raw_records.len());
            batch_stats.record_rejection(RejectionReason::InsufficientReads, raw_records.len());
            if track_rejects {
                for raw in &raw_records {
                    super::append_framed_bytes(&mut rejects_bytes, raw.as_ref())?;
                }
            }
            continue;
        }

        if let Some(ref mut oc) = state.overlapping {
            oc.reset_stats();
            // A failure here must be fatal so the absent-`--threads`
            // (single-worker) and `--threads N` chain runs behave identically.
            // Downgrading to a `RejectionReason::Other` reject would make the
            // same input fail at one worker count but succeed at another.
            apply_overlapping_consensus(&mut raw_records, oc).map_err(|e| {
                io::Error::other(format!("Overlapping consensus error for MI {mi}: {e}"))
            })?;
            batch_overlapping.merge(oc.stats());
        }

        let group_output = state
            .caller
            .consensus_reads(raw_records)
            .map_err(|e| io::Error::other(format!("Consensus error for MI {mi}: {e}")))?;
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

    // Progress logging at million-record boundaries (mirrors legacy's
    // `ProgressTracker::log_if_needed`).
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

/// Metrics-on variant of the per-batch simplex body. Collects one
/// `(TemplateInfo, ReadInfoKey)` entry per qualifying template across **every**
/// family in the batch — done before `run_simplex_consensus_batch` consumes the
/// groups, so families the `min_reads` threshold later rejects are still counted
/// (matching the separate-pass `simplex-metrics` command) — then splits those
/// entries into maximal same-key runs (`split_into_runs`) and classifies each
/// run as an interior (whole, complete-within-this-batch) coordinate group or a
/// batch-boundary run needing cross-batch reassembly (`classify_batch_runs`).
/// Interior runs are recorded directly into this worker's `ConsensusMetricsSlot`
/// under one `with_slot` lock acquisition per batch; boundary runs are
/// submitted to the shared `BoundaryReorder`, which closes coordinate groups
/// incrementally as the batch-serial prefix becomes contiguous, and whichever
/// groups that submission closes are recorded into this worker's slot under a
/// second `with_slot` acquisition, taken after the reorder's own mutex is
/// released (H3 design §6.3 — no lock-order cycle). This replaces the
/// old 3-branch `CoordinateGroupFragment`/serial `MetricsCollectorStep` design:
/// metrics recording now happens inline in the consensus worker body, so the
/// step's Outputs shape is unchanged from the metrics-off variant. Delegates
/// the unchanged consensus calling to `run_simplex_consensus_batch`. Selected
/// at chain-build time only when this stage's `metrics` field is set
/// (`add_simplex`); the metrics-off build calls `run_simplex_consensus_batch`
/// directly, so there is no runtime metrics work on that path (spec §7.1).
#[allow(clippy::too_many_arguments)]
fn run_simplex_consensus_batch_with_metrics(
    state: &mut ConsensusState,
    item: BatchedMiGroups,
    track_rejects: bool,
    overlapping_enabled: bool,
    min_reads: usize,
    accumulators: &Arc<PerThreadAccumulator<CollectedSimplexMetrics>>,
    progress: &Arc<AtomicU64>,
    qc: &Arc<crate::inline_metrics_collector::ConsensusMetricsCaptures>,
) -> io::Result<(DecompressedBlock, Option<DecompressedBlock>)> {
    let batch_serial = item.batch_serial;
    let mut metrics_entries = Vec::new();
    for group in &item.groups {
        push_mi_group_entries(group, &state.header, &state.library_index, &mut metrics_entries)
            .map_err(|e| io::Error::other(format!("metrics conversion error: {e:#}")))?;
    }
    let runs = crate::inline_metrics_collector::split_into_runs(metrics_entries);
    let (interior, boundary) =
        crate::inline_metrics_collector::classify_batch_runs(batch_serial, runs);
    qc.accumulator
        .with_slot(|slot| -> anyhow::Result<()> {
            for (_key, templates) in interior {
                slot.acc.record_coordinate_group(&templates, &qc.intervals)?;
            }
            Ok(())
        })
        .map_err(io::Error::other)?;
    // `submit` takes/releases the reorder mutex here, before the second
    // `with_slot` acquisition below — never nested — so there is no
    // lock-order cycle. Submitted unconditionally, even when `boundary` is
    // empty: batch serials are contiguous, and a batch with no metrics
    // entries still occupies its serial slot (H3 design §6.3).
    let closed = qc.reorder.submit(batch_serial, boundary);
    qc.accumulator
        .with_slot(|slot| -> anyhow::Result<()> {
            for group in &closed {
                slot.acc.record_coordinate_group(group, &qc.intervals)?;
            }
            Ok(())
        })
        .map_err(io::Error::other)?;

    run_simplex_consensus_batch(
        state,
        item,
        track_rejects,
        overlapping_enabled,
        min_reads,
        accumulators,
        progress,
    )
}

/// Build the 2-output `SimplexConsensus` step (used when BOTH `--rejects` and
/// `--metrics` are set): branch 0 = consensus, branch 1 = rejects. Same
/// Outputs shape as [`build_simplex_consensus_step_with_rejects`] — metrics
/// are recorded inline into `cap.qc_metrics`'s per-thread accumulator rather
/// than via an extra output branch. Parallel, `ByItemOrdinal`.
///
/// `pub(crate)` — consumed only by `ChainBuilder::add_simplex`.
///
/// # Panics
///
/// Panics if `cap.qc_metrics` is `None`. Unreachable in practice: `ChainBuilder::add_simplex`
/// selects this metrics-ON builder only when `--metrics` is set, in which case `qc_metrics`
/// is always `Some`.
pub(crate) fn build_simplex_consensus_step_with_rejects_and_metrics(
    limit_bytes: u64,
    cap: SimplexConsensusCaptures,
) -> impl Step<Input = BatchedMiGroups, Outputs = OrderedBytesTuple2<DecompressedBlock, DecompressedBlock>>
{
    let SimplexConsensusCaptures {
        track_rejects,
        overlapping_enabled,
        consensus_options,
        read_name_prefix,
        read_group_id,
        methylation_ref,
        accumulators,
        min_reads,
        progress,
        header,
        library_index,
        qc_metrics,
    } = cap;
    let qc = qc_metrics
        .expect("build_simplex_consensus_step_with_rejects_and_metrics requires qc_metrics");

    let init = make_simplex_consensus_init(
        read_name_prefix,
        read_group_id,
        consensus_options,
        methylation_ref,
        track_rejects,
        overlapping_enabled,
        header,
        library_index,
    );
    let body = move |state: &mut ConsensusState,
                     item: BatchedMiGroups|
          -> io::Result<Process2Output<DecompressedBlock, DecompressedBlock>> {
        let (consensus, rejects) = run_simplex_consensus_batch_with_metrics(
            state,
            item,
            track_rejects,
            overlapping_enabled,
            min_reads,
            &accumulators,
            &progress,
            &qc,
        )?;
        // X5-001 dense-serial rule: emit a (zero-byte) rejects block on an
        // all-clean batch so the rejects branch's reorder stage sees a dense
        // serial sequence, exactly as the metrics-off rejects builder does.
        let rejects = rejects.unwrap_or(DecompressedBlock {
            batch_serial: consensus.batch_serial,
            bytes: Vec::new(),
        });
        Ok(Process2Output::both(consensus, rejects))
    };

    process2_with_worker_state::<
        BatchedMiGroups,
        DecompressedBlock,
        DecompressedBlock,
        ConsensusState,
        _,
        _,
    >("SimplexConsensus", limit_bytes, limit_bytes, init, body)
}

/// Build the 1-output (kept-only) `SimplexConsensus` step (used when
/// `--metrics` is set but `--rejects` is not): same Outputs shape as
/// [`build_simplex_consensus_step_kept_only`] — metrics are recorded inline
/// into `cap.qc_metrics`'s per-thread accumulator rather than via an extra
/// output branch.
///
/// `pub(crate)` — consumed only by `ChainBuilder::add_simplex`.
///
/// # Panics
///
/// Panics if `cap.qc_metrics` is `None`. Unreachable in practice: `ChainBuilder::add_simplex`
/// selects this metrics-ON builder only when `--metrics` is set, in which case `qc_metrics`
/// is always `Some`.
#[allow(clippy::type_complexity)]
pub(crate) fn build_simplex_consensus_step_metrics(
    limit_bytes: u64,
    cap: SimplexConsensusCaptures,
) -> ProcessWithWorkerState<
    BatchedMiGroups,
    DecompressedBlock,
    impl Fn(&mut ConsensusState, BatchedMiGroups) -> io::Result<DecompressedBlock>
    + Send
    + Sync
    + 'static,
    ConsensusState,
    impl Fn() -> ConsensusState + Send + Sync + 'static,
> {
    let SimplexConsensusCaptures {
        track_rejects,
        overlapping_enabled,
        consensus_options,
        read_name_prefix,
        read_group_id,
        methylation_ref,
        accumulators,
        min_reads,
        progress,
        header,
        library_index,
        qc_metrics,
    } = cap;
    let qc = qc_metrics.expect("build_simplex_consensus_step_metrics requires qc_metrics");

    let init = make_simplex_consensus_init(
        read_name_prefix,
        read_group_id,
        consensus_options,
        methylation_ref,
        track_rejects,
        overlapping_enabled,
        header,
        library_index,
    );
    let body =
        move |state: &mut ConsensusState, item: BatchedMiGroups| -> io::Result<DecompressedBlock> {
            let (consensus, _rejects) = run_simplex_consensus_batch_with_metrics(
                state,
                item,
                track_rejects,
                overlapping_enabled,
                min_reads,
                &accumulators,
                &progress,
                &qc,
            )?;
            Ok(consensus)
        };

    process_with_worker_state::<BatchedMiGroups, DecompressedBlock, _, ConsensusState, _>(
        "SimplexConsensus",
        limit_bytes,
        init,
        body,
    )
}

/// Build the 2-output `SimplexConsensus` step (used when `--rejects` is set):
/// branch 0 carries the consensus `DecompressedBlock`, branch 1 carries the
/// rejects `DecompressedBlock`. Parallel, `ByItemOrdinal`.
///
/// `pub(crate)` — consumed only by `ChainBuilder::add_simplex`.
pub(crate) fn build_simplex_consensus_step_with_rejects(
    limit_bytes: u64,
    cap: SimplexConsensusCaptures,
) -> impl Step<Input = BatchedMiGroups, Outputs = OrderedBytesTuple2<DecompressedBlock, DecompressedBlock>>
{
    let SimplexConsensusCaptures {
        track_rejects,
        overlapping_enabled,
        consensus_options,
        read_name_prefix,
        read_group_id,
        methylation_ref,
        accumulators,
        min_reads,
        progress,
        header,
        library_index,
        qc_metrics: _,
    } = cap;

    let init = make_simplex_consensus_init(
        read_name_prefix,
        read_group_id,
        consensus_options,
        methylation_ref,
        track_rejects,
        overlapping_enabled,
        header,
        library_index,
    );
    let body = move |state: &mut ConsensusState,
                     item: BatchedMiGroups|
          -> io::Result<Process2Output<DecompressedBlock, DecompressedBlock>> {
        let (consensus, rejects) = run_simplex_consensus_batch(
            state,
            item,
            track_rejects,
            overlapping_enabled,
            min_reads,
            &accumulators,
            &progress,
        )?;
        if let Some(r) = rejects {
            Ok(Process2Output::both(consensus, r))
        } else {
            // X5-001: an all-clean batch must still emit a (zero-byte) rejects
            // block so the rejects branch's `ByItemOrdinal` reorder stage sees a
            // dense serial sequence. Pushing nothing here (`only_a`) leaves a
            // permanent gap at this serial that wedges `try_pop_in_order`, so
            // the rejects sink never drains. The empty `Vec` does not allocate
            // and produces no physical BGZF block (compress/write of `&[]` are
            // no-ops; the single EOF is appended once at drain).
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
        ConsensusState,
        _,
        _,
    >("SimplexConsensus", limit_bytes, limit_bytes, init, body)
}

/// Build the 1-output (kept-only) `SimplexConsensus` step (used when
/// `--rejects` is unset). `track_rejects` is `false`, so no rejects are
/// produced; the framework has no public discard sink, so omitting the rejects
/// branch requires this single-output variant (mirrors `correct_step_kept_only`).
///
/// `pub(crate)` — consumed only by `ChainBuilder::add_simplex`.
#[allow(clippy::type_complexity)]
pub(crate) fn build_simplex_consensus_step_kept_only(
    limit_bytes: u64,
    cap: SimplexConsensusCaptures,
) -> ProcessWithWorkerState<
    BatchedMiGroups,
    DecompressedBlock,
    impl Fn(&mut ConsensusState, BatchedMiGroups) -> io::Result<DecompressedBlock>
    + Send
    + Sync
    + 'static,
    ConsensusState,
    impl Fn() -> ConsensusState + Send + Sync + 'static,
> {
    let SimplexConsensusCaptures {
        track_rejects,
        overlapping_enabled,
        consensus_options,
        read_name_prefix,
        read_group_id,
        methylation_ref,
        accumulators,
        min_reads,
        progress,
        header,
        library_index,
        qc_metrics: _,
    } = cap;

    let init = make_simplex_consensus_init(
        read_name_prefix,
        read_group_id,
        consensus_options,
        methylation_ref,
        track_rejects,
        overlapping_enabled,
        header,
        library_index,
    );
    let body =
        move |state: &mut ConsensusState, item: BatchedMiGroups| -> io::Result<DecompressedBlock> {
            let (consensus, _rejects) = run_simplex_consensus_batch(
                state,
                item,
                track_rejects,
                overlapping_enabled,
                min_reads,
                &accumulators,
                &progress,
            )?;
            Ok(consensus)
        };

    process_with_worker_state::<BatchedMiGroups, DecompressedBlock, _, ConsensusState, _>(
        "SimplexConsensus",
        limit_bytes,
        init,
        body,
    )
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::commands::simplex::Simplex;
    use clap::Parser;

    /// The `SimplexOptions -> VanillaUmiConsensusOptions` mapping inside
    /// `add_simplex` must carry every tuning flag. The
    /// `test_simplex_chain_matches_single_threaded` parity tests compare two
    /// chain runs that BOTH go through this mapping, so they cannot catch a
    /// field dropped or swapped inside it; this test exercises the mapping
    /// directly, driven through `try_parse_from` + `to_simplex_options()` (not
    /// a hand-built `SimplexOptions` literal) with every value non-default, so
    /// a field read from the wrong source fails rather than coincidentally
    /// matching a default. Mirrors
    /// `to_simplex_options_carries_every_tuning_flag` in
    /// `commands::simplex`, which covers the other half of the pipeline
    /// (CLI args -> `SimplexOptions`).
    #[test]
    fn simplex_consensus_options_carries_every_tuning_flag() {
        let cmd = Simplex::try_parse_from([
            "simplex",
            "-i",
            "in.bam",
            "-o",
            "out.bam",
            "--error-rate-pre-umi",
            "40",
            "--error-rate-post-umi",
            "35",
            "--min-input-base-quality",
            "17",
            "--output-per-base-tags=false",
            "--trim=true",
            "--min-consensus-base-quality",
            "19",
            "--tie-rule",
            "ulp-relative",
            "--min-reads",
            "3",
            "--max-reads",
            "77",
            "--methylation-mode",
            "em-seq",
            "--ref",
            "ref.fa",
        ])
        .expect("parses");
        let simplex_options = cmd.to_simplex_options();

        let cell_tag =
            noodles::sam::alignment::record::data::field::Tag::from(crate::sam::SamTag::CB);
        let opts = simplex_consensus_options(&simplex_options, cell_tag);

        assert_eq!(opts.tag, "MI");
        assert_eq!(opts.error_rate_pre_umi, 40);
        assert_eq!(opts.error_rate_post_umi, 35);
        assert_eq!(opts.min_input_base_quality, 17);
        assert!(!opts.produce_per_base_tags, "an explicit false must not be lost");
        assert!(opts.trim);
        assert_eq!(opts.min_consensus_base_quality, 19);
        assert_eq!(
            opts.tie_rule,
            fgumi_consensus::TieRule::UlpRelative,
            "--tie-rule must reach the mapping"
        );
        assert_eq!(opts.min_reads, 3);
        assert_eq!(opts.max_reads, Some(77));
        assert_eq!(opts.cell_tag, Some(cell_tag));
        assert_eq!(
            opts.methylation_mode,
            fgumi_consensus::MethylationMode::EmSeq,
            "--methylation-mode must reach the mapping"
        );
    }
}
