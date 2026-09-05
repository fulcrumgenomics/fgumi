//! Chain builder for `Stage::Codec`.
//!
//! Phase 2 (T2.20) held the full ~490-LOC chain construction here.
//! Phase 3 (T3a.10) lifts that logic into
//! [`crate::pipeline::chains::builder::ChainBuilder`]; this module
//! now holds the codec-specific types and step factory that the builder
//! imports: `CodecFinalizeHook` and the `build_codec_consensus_step_with_rejects` / `build_codec_consensus_step_kept_only` factories, used
//! by `ChainBuilder::add_codec`.
//!
//! Mirrors duplex (2dc8bd8). Codec's chain is structurally the same with a
//! different consensus caller and no MI-tag transform, overlapping consensus,
//! or methylation mode — matching fgbio's `CallCodecConsensusReads` behaviour.
//!
//! This is the chain-builder construction path for the codec stage, consumed via
//! `ChainBuilder` / [`crate::pipeline::chains::build::build_for`]. It is the
//! sole path for `Codec::execute` (via `execute_chain`), with or without
//! `--threads` — absent `--threads` runs the chain at a single worker.

use std::io;
use std::sync::Arc;
use std::sync::atomic::{AtomicU64, Ordering};

use anyhow::Result;
use log::info;

use crate::commands::consensus_runner::ConsensusStatsOps;
use crate::consensus::codec_caller::{
    CodecConsensusCaller, CodecConsensusOptions, CodecConsensusStats,
};
use crate::consensus_caller::ConsensusOutput;
use crate::logging::OperationTimer;
use crate::mi_group::MiGroup;
use crate::per_thread_accumulator::PerThreadAccumulator;
use crate::pipeline::chains::FinalizeHook;
use crate::pipeline::core::outputs::OrderedBytesTuple2;
use crate::pipeline::core::step::Step;
use crate::pipeline::steps::group::mi::BatchedMiGroups;
use crate::pipeline::steps::process::{
    Process2Output, ProcessWithWorkerState, process_with_worker_state, process2_with_worker_state,
};
use crate::pipeline::steps::types::DecompressedBlock;

// ─────────────────────────────────────────────────────────────────────────────
// CollectedCodecMetrics
// ─────────────────────────────────────────────────────────────────────────────

/// Per-thread accumulator for codec consensus metrics.
///
/// Merged into final aggregates after the pipeline completes; one instance
/// per worker slot (see [`PerThreadAccumulator`]).
///
/// `pub(crate)` so `ChainBuilder` can construct it in `add_codec`.
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
/// `pub(crate)` so `ChainBuilder` can construct and register it in
/// `add_codec`.
pub(crate) struct CodecFinalizeHook {
    pub(crate) accumulators: Arc<PerThreadAccumulator<CollectedCodecMetrics>>,
    pub(crate) stats_path: Option<std::path::PathBuf>,
    pub(crate) timer: OperationTimer,
}

impl FinalizeHook for CodecFinalizeHook {
    fn finalize(self: Box<Self>) -> Result<()> {
        let CodecFinalizeHook { accumulators, stats_path, timer } = *self;

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
            // Write the key-value stats format, matching the sibling
            // simplex/duplex chain hooks; the wide-table `[metrics]` form
            // diverged from both.
            let kv_metrics = metrics.to_kv_metrics(fgumi_metrics::ConsensusCallerKind::Codec);
            DelimFile::default().write_tsv(stats_path, kv_metrics).map_err(|e| {
                anyhow::anyhow!("Failed to write statistics: {}: {e}", stats_path.display())
            })?;
            info!("Wrote statistics to: {}", stats_path.display());
        }

        info!("Wrote {consensus_count} CODEC consensus reads");

        // Rejects are emitted on the consensus step's second output branch
        // (fan-out) and finalized by the rejects branch's WriteBgzfFile sink.

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

/// Maps `CodecOptions` (the resolved options-bag entry `add_codec` receives)
/// into the `CodecConsensusOptions` the per-worker consensus caller is built
/// from.
///
/// Extracted out of `ChainBuilder::add_codec` so this mapping is
/// unit-testable on its own: the `test_codec_chain_matches_single_threaded`
/// parity tests compare two chain runs (single- vs multi-worker) that BOTH go
/// through this same mapping, so a field dropped or swapped inside it would
/// make the two sides agree with each other and still pass. See
/// `codec_consensus_options_carries_every_tuning_flag` below, which mirrors
/// `to_codec_options_carries_every_tuning_flag`'s non-default-everywhere
/// discipline for the other half of the CLI-args -> `CodecOptions` ->
/// `CodecConsensusOptions` pipeline.
pub(crate) fn codec_consensus_options(
    codec: &crate::commands::codec::CodecOptions,
    cell_tag: noodles::sam::alignment::record::data::field::Tag,
) -> CodecConsensusOptions {
    let consensus = codec.consensus();
    CodecConsensusOptions {
        min_input_base_quality: consensus.min_input_base_quality,
        error_rate_pre_umi: consensus.error_rate_pre_umi,
        error_rate_post_umi: consensus.error_rate_post_umi,
        min_reads_per_strand: codec.min_reads,
        max_reads_per_strand: codec.max_reads,
        min_duplex_length: codec.min_duplex_length,
        single_strand_qual: codec.single_strand_qual,
        outer_bases_qual: codec.outer_bases_qual,
        outer_bases_length: codec.outer_bases_length,
        max_duplex_disagreements: codec.max_duplex_disagreements.unwrap_or(usize::MAX),
        max_duplex_disagreement_rate: codec.max_duplex_disagreement_rate,
        legacy_overlap_window: codec.legacy_overlap_window,
        cell_tag: Some(cell_tag),
        produce_per_base_tags: consensus.output_per_base_tags,
        trim: consensus.trim,
        min_consensus_base_quality: consensus.min_consensus_base_quality,
        // `codec.tie_rule` is already the resolved `TieRule` on `CodecOptions`.
        tie_rule: codec.tie_rule,
    }
}

/// Captures passed into [`build_codec_consensus_step_with_rejects`] / [`build_codec_consensus_step_kept_only`] from `add_codec`.
///
/// Bundles all the cloned scalars and Arcs the closure needs so `add_codec`
/// can prepare them once and hand them off cleanly.
///
/// `pub(crate)` — consumed only by `ChainBuilder::add_codec` and
/// [`build_codec_consensus_step_with_rejects`] / [`build_codec_consensus_step_kept_only`].
pub(crate) struct CodecConsensusCaptures {
    pub(crate) track_rejects: bool,
    pub(crate) read_name_prefix: String,
    pub(crate) read_group_id: String,
    pub(crate) consensus_options: CodecConsensusOptions,
    pub(crate) accumulators: Arc<PerThreadAccumulator<CollectedCodecMetrics>>,
    pub(crate) progress: Arc<AtomicU64>,
}

/// Per-worker init: build the `CodecConsensusCaller` once, reused across
/// batches. Shared by both step variants.
fn make_codec_consensus_init(
    read_name_prefix: String,
    read_group_id: String,
    consensus_options: CodecConsensusOptions,
    track_rejects: bool,
) -> impl Fn() -> CodecState + Send + Sync + 'static {
    move || {
        let caller = CodecConsensusCaller::new_with_rejects_tracking(
            read_name_prefix.clone(),
            read_group_id.clone(),
            consensus_options.clone(),
            track_rejects,
        );
        CodecState { caller }
    }
}

/// Per-batch CODEC consensus body, shared by both step variants.
///
/// Returns the consensus `DecompressedBlock` (branch 0) and, when
/// `track_rejects` is set and any record was rejected, a rejects
/// `DecompressedBlock` (branch 1) of `[len][record]`-framed raw-input records
/// in input order — the PR #332 fan-out contract (rejects flow through the
/// ordered serialize/compress stages, not a mutex side-channel; the rejects
/// writer is configured with the input header by `add_codec`).
fn run_codec_consensus_batch(
    state: &mut CodecState,
    item: BatchedMiGroups,
    track_rejects: bool,
    accumulators: &Arc<PerThreadAccumulator<CollectedCodecMetrics>>,
    progress: &Arc<AtomicU64>,
) -> io::Result<(DecompressedBlock, Option<DecompressedBlock>)> {
    let BatchedMiGroups { batch_serial, groups } = item;
    let groups_count = groups.len() as u64;

    let mut all_output = ConsensusOutput::default();
    let mut batch_stats = CodecConsensusStats::default();
    let mut rejects_bytes: Vec<u8> = Vec::new();

    let mut total_input_records: u64 = 0;
    for MiGroup { mi, records } in groups {
        state.caller.clear();
        total_input_records += records.len() as u64;

        // Use the typed entry point so recoverable "duplex disagreement" errors
        // are classified by the `CodecConsensusError::is_duplex_disagreement()`
        // predicate. The generic `consensus_reads` erases the typed variant
        // into `anyhow::Error`, which would force a fragile substring match on
        // the Display text.
        let result = state.caller.consensus_reads_typed(records);
        match result {
            Ok(group_output) => {
                all_output.merge(group_output);
                batch_stats.merge(state.caller.statistics());
                if track_rejects {
                    for raw in &state.caller.take_rejected_reads() {
                        super::append_framed_bytes(&mut rejects_bytes, raw)?;
                    }
                }
            }
            Err(e) => {
                // A "duplex disagreement" error isn't fatal — preserve stats +
                // rejects and move on. Anything else is hard.
                if e.is_duplex_disagreement() {
                    batch_stats.merge(state.caller.statistics());
                    if track_rejects {
                        for raw in &state.caller.take_rejected_reads() {
                            super::append_framed_bytes(&mut rejects_bytes, raw)?;
                        }
                    }
                } else {
                    // `{:#}` on the anyhow view preserves the nested cause chain,
                    // matching the context group.rs keeps for an equivalent failure.
                    let err = anyhow::Error::from(e);
                    return Err(io::Error::other(format!(
                        "CODEC consensus error for MI {mi}: {err:#}"
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

/// Build the 2-output `CodecConsensus` step (used when `--rejects` is set):
/// branch 0 = consensus, branch 1 = rejects.
///
/// `pub(crate)` — consumed only by `ChainBuilder::add_codec`.
pub(crate) fn build_codec_consensus_step_with_rejects(
    limit_bytes: u64,
    cap: CodecConsensusCaptures,
) -> impl Step<Input = BatchedMiGroups, Outputs = OrderedBytesTuple2<DecompressedBlock, DecompressedBlock>>
{
    let CodecConsensusCaptures {
        track_rejects,
        read_name_prefix,
        read_group_id,
        consensus_options,
        accumulators,
        progress,
    } = cap;

    let init = make_codec_consensus_init(
        read_name_prefix,
        read_group_id,
        consensus_options,
        track_rejects,
    );
    let body = move |state: &mut CodecState,
                     item: BatchedMiGroups|
          -> io::Result<Process2Output<DecompressedBlock, DecompressedBlock>> {
        let (consensus, rejects) =
            run_codec_consensus_batch(state, item, track_rejects, &accumulators, &progress)?;
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
        CodecState,
        _,
        _,
    >("CodecConsensus", limit_bytes, limit_bytes, init, body)
}

/// Build the 1-output kept-only `CodecConsensus` step (used when `--rejects`
/// is unset). The framework has no public discard sink, so omitting the rejects
/// branch requires this single-output variant.
///
/// `pub(crate)` — consumed only by `ChainBuilder::add_codec`.
#[allow(clippy::type_complexity)]
pub(crate) fn build_codec_consensus_step_kept_only(
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
        read_name_prefix,
        read_group_id,
        consensus_options,
        accumulators,
        progress,
    } = cap;

    let init = make_codec_consensus_init(
        read_name_prefix,
        read_group_id,
        consensus_options,
        track_rejects,
    );
    let body =
        move |state: &mut CodecState, item: BatchedMiGroups| -> io::Result<DecompressedBlock> {
            let (consensus, _rejects) =
                run_codec_consensus_batch(state, item, track_rejects, &accumulators, &progress)?;
            Ok(consensus)
        };

    process_with_worker_state::<BatchedMiGroups, DecompressedBlock, _, CodecState, _>(
        "CodecConsensus",
        limit_bytes,
        init,
        body,
    )
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::commands::codec::Codec;
    use clap::Parser;

    /// The `CodecOptions -> CodecConsensusOptions` mapping inside `add_codec`
    /// must carry every tuning flag. The
    /// `test_codec_chain_matches_single_threaded` parity tests compare two
    /// chain runs that BOTH go through this mapping, so they cannot catch a
    /// field dropped or swapped inside it; this test exercises the mapping
    /// directly, driven through `try_parse_from` + `to_codec_options()` (not
    /// a hand-built `CodecOptions` literal) with every value non-default, so
    /// a field read from the wrong source fails rather than coincidentally
    /// matching a default. Mirrors `to_codec_options_carries_every_tuning_flag`
    /// in `commands::codec`, which covers the other half of the pipeline
    /// (CLI args -> `CodecOptions`).
    #[test]
    fn codec_consensus_options_carries_every_tuning_flag() {
        let cmd = Codec::try_parse_from([
            "codec",
            "-i",
            "in.bam",
            "-o",
            "out.bam",
            "--error-rate-pre-umi",
            "42",
            "--error-rate-post-umi",
            "37",
            "--min-input-base-quality",
            "16",
            "--output-per-base-tags=false",
            "--trim=true",
            "--min-consensus-base-quality",
            "23",
            "--tie-rule",
            "ulp-relative",
            "--min-reads",
            "4",
            "--max-reads",
            "88",
            "--min-duplex-length",
            "9",
            "--legacy-overlap-window",
            "--single-strand-qual",
            "12",
            "--outer-bases-qual",
            "15",
            "--outer-bases-length",
            "7",
            "--max-duplex-disagreement-rate",
            "0.25",
            "--max-duplex-disagreements",
            "6",
        ])
        .expect("parses");
        let codec_options = cmd.to_codec_options();

        let cell_tag =
            noodles::sam::alignment::record::data::field::Tag::from(crate::sam::SamTag::CB);
        let opts = codec_consensus_options(&codec_options, cell_tag);

        assert_eq!(opts.min_input_base_quality, 16);
        assert_eq!(opts.error_rate_pre_umi, 42);
        assert_eq!(opts.error_rate_post_umi, 37);
        assert_eq!(opts.min_reads_per_strand, 4);
        assert_eq!(opts.max_reads_per_strand, Some(88));
        assert_eq!(opts.min_duplex_length, 9);
        assert_eq!(opts.single_strand_qual, Some(12));
        assert_eq!(opts.outer_bases_qual, Some(15));
        assert_eq!(opts.outer_bases_length, 7);
        assert_eq!(opts.max_duplex_disagreements, 6);
        assert!((opts.max_duplex_disagreement_rate - 0.25).abs() < f64::EPSILON);
        assert!(opts.legacy_overlap_window, "--legacy-overlap-window must reach the mapping");
        assert_eq!(opts.cell_tag, Some(cell_tag));
        assert!(!opts.produce_per_base_tags, "an explicit false must not be lost");
        assert!(opts.trim);
        assert_eq!(opts.min_consensus_base_quality, 23);
        assert_eq!(
            opts.tie_rule,
            fgumi_consensus::TieRule::UlpRelative,
            "--tie-rule must reach the mapping"
        );
    }
}
