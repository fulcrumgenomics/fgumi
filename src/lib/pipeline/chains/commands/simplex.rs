//! Chain builder for `Stage::Simplex`.
//!
//! Phase 2 (T2.18) held the full ~500-LOC chain construction here.
//! Phase 3 (T3a.8) lifts that logic into
//! [`crate::pipeline::chains::builder::ChainBuilder`]; this module
//! now holds the simplex-specific types and step factory that the builder
//! imports: [`SimplexFinalizeHook`] and [`build_simplex_consensus_step`] used
//! by `ChainBuilder::add_simplex`.
//!
//! [`build_simplex_chain`] shrinks to a ~10-line delegate.
//!
//! Mirrors the dedup (6423252), filter (8eead93), clip (8e95428), group
//! (d371d84) reference migrations.
//!
//! The single-threaded fast path (`--threads` unset, BAM-only) is **not**
//! migrated here — it stays inline in `Simplex::execute`.

use std::io;
use std::sync::Arc;
use std::sync::atomic::{AtomicU64, Ordering};

use anyhow::Result;
use log::info;
use parking_lot::Mutex;

use crate::commands::common::MethylationRef;
use crate::commands::consensus_runner::{ConsensusStatsOps, log_overlapping_stats};
use crate::consensus_caller::{
    ConsensusCaller, ConsensusCallingStats, ConsensusOutput, RejectionReason,
};
use crate::logging::OperationTimer;
use crate::mi_group::MiGroup;
use crate::overlapping_consensus::{
    AgreementStrategy, CorrectionStats, DisagreementStrategy, OverlappingBasesConsensusCaller,
    apply_overlapping_consensus,
};
use crate::per_thread_accumulator::PerThreadAccumulator;
use crate::pipeline::chains::builder::ChainBuilder;
use crate::pipeline::chains::{BuiltPipeline, ChainSpec, FinalizeHook};
use crate::pipeline::steps::group::mi::BatchedMiGroups;
use crate::pipeline::steps::process::{ProcessWithWorkerState, process_with_worker_state};
use crate::pipeline::steps::types::DecompressedBlock;
use crate::vanilla_consensus_caller::{VanillaUmiConsensusCaller, VanillaUmiConsensusOptions};
use fgumi_bam_io::RawBamWriter;
use fgumi_raw_bam::RawRecord;

// ─────────────────────────────────────────────────────────────────────────────
// CollectedSimplexMetrics
// ─────────────────────────────────────────────────────────────────────────────

/// Per-thread accumulator for simplex consensus metrics (mirrors
/// `commands::simplex::CollectedSimplexMetrics`).
///
/// Merged into final aggregates after the pipeline completes; one instance
/// per worker slot (see [`PerThreadAccumulator`]).
///
/// `pub(crate)` so [`ChainBuilder`] can construct it in `add_simplex`.
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
/// `pub(crate)` so [`ChainBuilder`] can construct and register it in
/// `add_simplex`.
pub(crate) struct SimplexFinalizeHook {
    pub(crate) accumulators: Arc<PerThreadAccumulator<CollectedSimplexMetrics>>,
    pub(crate) stats_path: Option<std::path::PathBuf>,
    pub(crate) overlapping_enabled: bool,
    pub(crate) rejects_writer: Option<Arc<Mutex<Option<RawBamWriter>>>>,
    pub(crate) timer: OperationTimer,
}

impl FinalizeHook for SimplexFinalizeHook {
    fn finalize(self: Box<Self>) -> Result<()> {
        let SimplexFinalizeHook {
            accumulators,
            stats_path,
            overlapping_enabled,
            rejects_writer,
            timer,
        } = *self;

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
            let kv_metrics = metrics.to_kv_metrics();
            DelimFile::default().write_tsv(stats_path, kv_metrics).map_err(|e| {
                anyhow::anyhow!("Failed to write statistics: {}: {e}", stats_path.display())
            })?;
            info!("Wrote statistics to: {}", stats_path.display());
        }

        info!("Wrote {consensus_count} consensus reads");

        // Finalize the rejects writer (always finalize even on success to flush BGZF EOF).
        if let Some(rw_arc) = rejects_writer {
            let maybe_writer = rw_arc.lock().take();
            if let Some(writer) = maybe_writer {
                writer
                    .finish()
                    .map_err(|e| anyhow::anyhow!("Failed to finish rejects file: {e}"))?;
                info!("Rejected reads streamed to rejects file during processing");
            }
        }

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

/// Captures passed into [`build_simplex_consensus_step`] from `add_simplex`.
///
/// Bundles all the cloned scalars and Arcs the closure needs so `add_simplex`
/// can prepare them once and hand them off cleanly.
///
/// `pub(crate)` — consumed only by [`ChainBuilder::add_simplex`] and
/// [`build_simplex_consensus_step`].
pub(crate) struct SimplexConsensusCaptures {
    pub(crate) track_rejects: bool,
    pub(crate) rejects_writer: Option<Arc<Mutex<Option<RawBamWriter>>>>,
    pub(crate) overlapping_enabled: bool,
    pub(crate) consensus_options: VanillaUmiConsensusOptions,
    pub(crate) read_name_prefix: String,
    pub(crate) read_group_id: String,
    pub(crate) methylation_ref: MethylationRef,
    pub(crate) accumulators: Arc<PerThreadAccumulator<CollectedSimplexMetrics>>,
    pub(crate) min_reads: usize,
    pub(crate) progress: Arc<AtomicU64>,
}

/// Build the `SimplexConsensus` step: parallel, `ByItemOrdinal`. Builds a
/// per-worker [`VanillaUmiConsensusCaller`] and optional
/// [`OverlappingBasesConsensusCaller`] once, reuses across batches. Streams
/// rejects inline through `cap.rejects_writer`.
///
/// Output: [`DecompressedBlock`] containing concatenated pre-serialized
/// consensus records, ready for `BgzfCompress`.
///
/// `pub(crate)` — consumed only by [`ChainBuilder::add_simplex`].
#[allow(clippy::type_complexity, clippy::too_many_lines)]
pub(crate) fn build_simplex_consensus_step(
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
        rejects_writer,
        overlapping_enabled,
        consensus_options,
        read_name_prefix,
        read_group_id,
        methylation_ref,
        accumulators,
        min_reads,
        progress,
    } = cap;

    process_with_worker_state::<BatchedMiGroups, DecompressedBlock, _, ConsensusState, _>(
        "SimplexConsensus",
        limit_bytes,
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
            ConsensusState { caller, overlapping }
        },
        move |state: &mut ConsensusState, item: BatchedMiGroups| -> io::Result<DecompressedBlock> {
            let BatchedMiGroups { batch_serial, groups } = item;
            let groups_count = groups.len() as u64;

            let mut all_output = ConsensusOutput::default();
            let mut batch_stats = ConsensusCallingStats::new();
            let mut batch_overlapping = CorrectionStats::new();

            // Stream per-MI-group rejects straight to disk so they
            // don't accumulate in a batch-level Vec. The mutex only
            // serializes the raw-byte append; BGZF compression runs
            // on the writer's own thread pool.
            let flush_raw_records = |recs: &[RawRecord]| -> io::Result<()> {
                if let Some(ref rw_arc) = rejects_writer {
                    if !recs.is_empty() {
                        let mut guard = rw_arc.lock();
                        if let Some(w) = guard.as_mut() {
                            for raw in recs {
                                w.write_raw_record(raw.as_ref())?;
                            }
                        }
                    }
                }
                Ok(())
            };
            let flush_byte_records = |recs: &[Vec<u8>]| -> io::Result<()> {
                if let Some(ref rw_arc) = rejects_writer {
                    if !recs.is_empty() {
                        let mut guard = rw_arc.lock();
                        if let Some(w) = guard.as_mut() {
                            for raw in recs {
                                w.write_raw_record(raw)?;
                            }
                        }
                    }
                }
                Ok(())
            };

            let mut total_input_records: u64 = 0;
            for MiGroup { mi, records: mut raw_records } in groups {
                state.caller.clear();
                total_input_records += raw_records.len() as u64;

                if raw_records.len() < min_reads {
                    batch_stats.record_input(raw_records.len());
                    batch_stats
                        .record_rejection(RejectionReason::InsufficientReads, raw_records.len());
                    if track_rejects {
                        flush_raw_records(&raw_records)?;
                    }
                    continue;
                }

                if let Some(ref mut oc) = state.overlapping {
                    oc.reset_stats();
                    // A failure here must be fatal to match the single-thread
                    // fast path in `src/lib/commands/simplex.rs`, which
                    // propagates `apply_overlapping_consensus` errors with `?`.
                    // Downgrading to a `RejectionReason::Other` reject would make
                    // the same input fail without `--threads` but succeed with it.
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
                    flush_byte_records(&state.caller.take_rejected_reads())?;
                }
            }

            // Merge per-batch metrics into this worker's slot.
            accumulators.with_slot(|m| {
                m.stats.merge(&batch_stats);
                if overlapping_enabled {
                    m.overlapping_stats
                        .get_or_insert_with(CorrectionStats::new)
                        .merge(&batch_overlapping);
                }
                m.groups_processed += groups_count;
            });

            // Progress logging at million-record boundaries (mirrors
            // legacy's `ProgressTracker::log_if_needed`).
            let prev = progress.fetch_add(total_input_records, Ordering::Relaxed);
            if (prev + total_input_records) / 1_000_000 > prev / 1_000_000 {
                info!("Processed {} records", prev + total_input_records);
            }

            Ok(DecompressedBlock { batch_serial, bytes: all_output.data })
        },
    )
}

// ─────────────────────────────────────────────────────────────────────────────
// build_simplex_chain — thin delegate (Phase 3 T3a.8)
// ─────────────────────────────────────────────────────────────────────────────

/// Build a [`BuiltPipeline`] for a simplex-only chain.
///
/// Delegates to [`ChainBuilder`]. The full step sequence lives in
/// [`ChainBuilder`]'s `add_simplex` method in
/// [`crate::pipeline::chains::builder`].
///
/// The single-threaded fast path (`--threads` unset, BAM-only) is **not**
/// migrated here — it stays inline in `Simplex::execute`.
///
/// # Errors
///
/// Returns input-validation errors (missing reference for methylation mode,
/// etc.) or any underlying pipeline construction error.
#[allow(clippy::needless_pass_by_value)]
pub fn build_simplex_chain(spec: ChainSpec) -> Result<BuiltPipeline> {
    use crate::pipeline::chains::{Stage, builder::StagePosition};
    let mut chain = ChainBuilder::new(&spec)?;
    chain.add_source()?;
    chain.add_stage(Stage::Simplex, StagePosition::Terminal)?;
    chain.add_sink()?;
    chain.build()
}
