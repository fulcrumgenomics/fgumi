//! Chain builder for `Stage::Codec`.
//!
//! Phase 2 (T2.20) held the full ~490-LOC chain construction here.
//! Phase 3 (T3a.10) lifts that logic into
//! [`crate::pipeline::chains::builder::ChainBuilder`]; this module
//! now holds the codec-specific types and step factory that the builder
//! imports: [`CodecFinalizeHook`] and [`build_codec_consensus_step`] used
//! by `ChainBuilder::add_codec`.
//!
//! [`build_codec_chain`] shrinks to a ~10-line delegate.
//!
//! Mirrors duplex (2dc8bd8). Codec's chain is structurally the same with a
//! different consensus caller and no MI-tag transform, overlapping consensus,
//! or methylation mode — matching fgbio's `CallCodecConsensusReads` behaviour.
//!
//! The single-threaded fast path (`--threads` unset, BAM-only) is **not**
//! migrated here — it stays inline in `Codec::execute`.

use std::io;
use std::sync::Arc;
use std::sync::atomic::{AtomicU64, Ordering};

use anyhow::Result;
use log::info;
use parking_lot::Mutex;

use crate::commands::consensus_runner::ConsensusStatsOps;
use crate::consensus::codec_caller::{
    CodecConsensusCaller, CodecConsensusOptions, CodecConsensusStats,
};
use crate::consensus_caller::{ConsensusCaller, ConsensusOutput};
use crate::logging::OperationTimer;
use crate::mi_group::MiGroup;
use crate::per_thread_accumulator::PerThreadAccumulator;
use crate::pipeline::chains::builder::ChainBuilder;
use crate::pipeline::chains::{BuiltPipeline, ChainSpec, FinalizeHook};
use crate::pipeline::steps::group::mi::BatchedMiGroups;
use crate::pipeline::steps::process::{ProcessWithWorkerState, process_with_worker_state};
use crate::pipeline::steps::types::DecompressedBlock;
use fgumi_bam_io::RawBamWriter;

// ─────────────────────────────────────────────────────────────────────────────
// CollectedCodecMetrics
// ─────────────────────────────────────────────────────────────────────────────

/// Per-thread accumulator for codec consensus metrics (mirrors
/// `commands::codec::CollectedCodecMetrics`).
///
/// Merged into final aggregates after the pipeline completes; one instance
/// per worker slot (see [`PerThreadAccumulator`]).
///
/// `pub(crate)` so [`ChainBuilder`] can construct it in `add_codec`.
#[derive(Default)]
pub(crate) struct CollectedCodecMetrics {
    /// CODEC consensus calling statistics.
    pub(crate) stats: CodecConsensusStats,
    /// Number of MI groups processed.
    pub(crate) groups_processed: u64,
}

// ─────────────────────────────────────────────────────────────────────────────
// CodecState
// ─────────────────────────────────────────────────────────────────────────────

/// Per-worker state for the codec consensus step.
///
/// Defined at module level (not inside a function body) to avoid the
/// `clippy::items_after_statements` lint that fires when a struct is defined
/// after executable statements in a function body.
pub(crate) struct CodecState {
    pub(crate) caller: CodecConsensusCaller,
}

impl crate::pipeline::core::item::HeapSize for CodecState {}

// ─────────────────────────────────────────────────────────────────────────────
// CodecFinalizeHook
// ─────────────────────────────────────────────────────────────────────────────

/// Post-pipeline finalize hook for codec. Reduces per-thread metrics,
/// writes the optional stats file, logs the summary banner, finalizes
/// the rejects writer, and calls `timer.log_completion`.
///
/// `pub(crate)` so [`ChainBuilder`] can construct and register it in
/// `add_codec`.
pub(crate) struct CodecFinalizeHook {
    pub(crate) accumulators: Arc<PerThreadAccumulator<CollectedCodecMetrics>>,
    pub(crate) stats_path: Option<std::path::PathBuf>,
    pub(crate) rejects_writer: Option<Arc<Mutex<Option<RawBamWriter>>>>,
    pub(crate) timer: OperationTimer,
}

impl FinalizeHook for CodecFinalizeHook {
    fn finalize(self: Box<Self>) -> Result<()> {
        let CodecFinalizeHook { accumulators, stats_path, rejects_writer, timer } = *self;

        // Reduce per-thread accumulators.
        let mut total_groups = 0u64;
        let mut merged_stats = CodecConsensusStats::default();

        for slot in accumulators.slots() {
            let m = slot.lock();
            total_groups += m.groups_processed;
            merged_stats.merge(&m.stats);
        }

        info!("CODEC consensus calling complete");
        info!("Total MI groups processed: {total_groups}");

        let metrics = merged_stats.to_metrics();
        let consensus_count = metrics.consensus_reads;
        crate::logging::log_consensus_summary(&metrics);

        if let Some(ref stats_path) = stats_path {
            use fgoxide::io::DelimFile;
            DelimFile::default().write_tsv(stats_path, [metrics]).map_err(|e| {
                anyhow::anyhow!("Failed to write statistics: {}: {e}", stats_path.display())
            })?;
            info!("Wrote statistics to: {}", stats_path.display());
        }

        info!("Wrote {consensus_count} CODEC consensus reads");

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
// Step factory — extracted from build_codec_chain (T3a.10).
//
// Returns the concrete `ProcessWithWorkerState` type directly (the same
// pattern used by duplex and simplex in this module family). Returning
// `impl Step<...>` is blocked here because the closure types embed in
// opaque-return position and cannot name themselves.
// ─────────────────────────────────────────────────────────────────────────────

/// Captures passed into [`build_codec_consensus_step`] from `add_codec`.
///
/// Bundles all the cloned scalars and Arcs the closure needs so `add_codec`
/// can prepare them once and hand them off cleanly.
///
/// `pub(crate)` — consumed only by [`ChainBuilder::add_codec`] and
/// [`build_codec_consensus_step`].
pub(crate) struct CodecConsensusCaptures {
    pub(crate) track_rejects: bool,
    pub(crate) rejects_writer: Option<Arc<Mutex<Option<RawBamWriter>>>>,
    pub(crate) read_name_prefix: String,
    pub(crate) read_group_id: String,
    pub(crate) consensus_options: CodecConsensusOptions,
    pub(crate) accumulators: Arc<PerThreadAccumulator<CollectedCodecMetrics>>,
    pub(crate) progress: Arc<AtomicU64>,
}

/// Build the `CodecConsensus` step: parallel, `ByItemOrdinal`. Builds a
/// per-worker [`CodecConsensusCaller`] once, reuses across batches. Streams
/// rejects inline through `cap.rejects_writer`.
///
/// Output: [`DecompressedBlock`] containing concatenated pre-serialized
/// consensus records, ready for `BgzfCompress`.
///
/// Unlike duplex, codec does not apply a record filter, MI-tag transform,
/// overlapping consensus, or methylation mode — matching fgbio's
/// `CallCodecConsensusReads` behaviour.
///
/// `pub(crate)` — consumed only by [`ChainBuilder::add_codec`].
///
/// # Panics
///
/// Does not panic during construction. Per-worker init failures are
/// surfaced immediately because `CodecConsensusCaller::new_with_rejects_tracking`
/// is infallible.
#[allow(clippy::type_complexity, clippy::too_many_lines)]
pub(crate) fn build_codec_consensus_step(
    limit_bytes: u64,
    cap: CodecConsensusCaptures,
) -> ProcessWithWorkerState<
    BatchedMiGroups,
    DecompressedBlock,
    impl Fn(&mut CodecState, BatchedMiGroups) -> io::Result<DecompressedBlock> + Send + Sync + 'static,
    CodecState,
    impl Fn() -> CodecState + Send + Sync + 'static,
> {
    let CodecConsensusCaptures {
        track_rejects,
        rejects_writer,
        read_name_prefix,
        read_group_id,
        consensus_options,
        accumulators,
        progress,
    } = cap;

    process_with_worker_state::<BatchedMiGroups, DecompressedBlock, _, CodecState, _>(
        "CodecConsensus",
        limit_bytes,
        move || {
            let caller = CodecConsensusCaller::new_with_rejects_tracking(
                read_name_prefix.clone(),
                read_group_id.clone(),
                consensus_options.clone(),
                track_rejects,
            );
            CodecState { caller }
        },
        move |state: &mut CodecState, item: BatchedMiGroups| -> io::Result<DecompressedBlock> {
            let BatchedMiGroups { batch_serial, groups } = item;
            let groups_count = groups.len() as u64;

            let mut all_output = ConsensusOutput::default();
            let mut batch_stats = CodecConsensusStats::default();

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
            for MiGroup { mi, records } in groups {
                state.caller.clear();
                total_input_records += records.len() as u64;

                let result: anyhow::Result<ConsensusOutput> = state.caller.consensus_reads(records);
                match result {
                    Ok(group_output) => {
                        all_output.merge(group_output);
                        batch_stats.merge(state.caller.statistics());
                        if track_rejects {
                            flush_byte_records(&state.caller.take_rejected_reads())?;
                        }
                    }
                    Err(e) => {
                        // Mirrors legacy: a "duplex disagreement" error
                        // isn't fatal — preserve stats + rejects and move
                        // on. Anything else is a hard error.
                        if e.to_string().contains("duplex disagreement") {
                            batch_stats.merge(state.caller.statistics());
                            if track_rejects {
                                flush_byte_records(&state.caller.take_rejected_reads())?;
                            }
                        } else {
                            return Err(io::Error::other(format!(
                                "CODEC consensus error for MI {mi}: {e}"
                            )));
                        }
                    }
                }
            }

            // Merge per-batch metrics into this worker's slot.
            accumulators.with_slot(|m| {
                m.stats.merge(&batch_stats);
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
// build_codec_chain — thin delegate (Phase 3 T3a.10)
// ─────────────────────────────────────────────────────────────────────────────

/// Build a [`BuiltPipeline`] for a codec-only chain.
///
/// Delegates to [`ChainBuilder`]. The full step sequence lives in
/// [`ChainBuilder`]'s `add_codec` method in
/// [`crate::pipeline::chains::builder`].
///
/// The single-threaded fast path (`--threads` unset, BAM-only) is **not**
/// migrated here — it stays inline in `Codec::execute`.
///
/// # Errors
///
/// Returns input-validation errors (missing file, etc.) or any underlying
/// pipeline construction error.
#[allow(clippy::needless_pass_by_value)]
pub fn build_codec_chain(spec: ChainSpec) -> Result<BuiltPipeline> {
    use crate::pipeline::chains::{Stage, builder::StagePosition};
    let mut chain = ChainBuilder::new(&spec)?;
    chain.add_source()?;
    chain.add_stage(Stage::Codec, StagePosition::Terminal)?;
    chain.add_sink()?;
    chain.build()
}
