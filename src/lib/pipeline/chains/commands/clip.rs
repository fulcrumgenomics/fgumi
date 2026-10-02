//! Chain builder for `Stage::Clip`.
//!
//! Phase 2 (T2.15) held the full ~400-LOC chain construction here.
//! Phase 3 (T3a.5) lifts that logic into
//! [`crate::pipeline::chains::builder::ChainBuilder`]; this module
//! now holds the clip-specific types and step factories that the builder
//! imports: `ClipAtomicMetrics`, `ClipFinalizeHook`, and the two
//! step-factory functions used by `ChainBuilder::add_clip`.

use std::io;
use std::path::PathBuf;
use std::sync::Arc;
use std::sync::atomic::{AtomicU64, Ordering};

use anyhow::Result;
use log::{info, warn};

use crate::clipper::RawRecordClipper;
use crate::commands::clip::ClipParams;
use crate::logging::OperationTimer;
use crate::metrics::clip::ClippingMetricsCollection;
use crate::per_thread_accumulator::PerThreadAccumulator;
use crate::pipeline::chains::FinalizeHook;
use crate::pipeline::steps::process::{ProcessOrdered, process_ordered};
use crate::pipeline::steps::serialize::SerializeBamRecords;
use crate::pipeline::steps::types::BamTemplateBatch;
use crate::reference::ReferenceReader;
use fgumi_raw_bam::RawRecord;

// ─────────────────────────────────────────────────────────────────────────────
// ClipAtomicMetrics
// ─────────────────────────────────────────────────────────────────────────────

/// Atomic counters for clip metrics. Shared across all worker clones of
/// the `ClipTemplates` step; each closure call adds its per-batch counts
/// via `fetch_add(_, Relaxed)`.
///
/// `pub(crate)` so `ChainBuilder` can construct and pass it to the step
/// factories and the finalize hook in `add_clip`.
pub(crate) struct ClipAtomicMetrics {
    pub(crate) total_templates: AtomicU64,
    pub(crate) overlap_clipped: AtomicU64,
    pub(crate) extend_clipped: AtomicU64,
    /// Records whose modification tags clipping could not keep in step and removed.
    pub(crate) modification_tags_removed: AtomicU64,
}

impl Default for ClipAtomicMetrics {
    fn default() -> Self {
        Self {
            total_templates: AtomicU64::new(0),
            overlap_clipped: AtomicU64::new(0),
            extend_clipped: AtomicU64::new(0),
            modification_tags_removed: AtomicU64::new(0),
        }
    }
}

// ─────────────────────────────────────────────────────────────────────────────
// ClipFinalizeHook
// ─────────────────────────────────────────────────────────────────────────────

/// Post-pipeline finalize hook for clip. Reads atomic counters, logs
/// the summary banner, and calls `timer.log_completion`.
///
/// `pub(crate)` so `ChainBuilder` can construct and register it in
/// `add_clip`.
pub(crate) struct ClipFinalizeHook {
    pub(crate) metrics: Arc<ClipAtomicMetrics>,
    pub(crate) progress_counter: Arc<AtomicU64>,
    pub(crate) timer: OperationTimer,
}

impl FinalizeHook for ClipFinalizeHook {
    fn finalize(self: Box<Self>) -> Result<()> {
        let ClipFinalizeHook { metrics, progress_counter, timer } = *self;

        let total_templates = metrics.total_templates.load(Ordering::Relaxed);
        let total_overlap_clipped = metrics.overlap_clipped.load(Ordering::Relaxed);
        let total_extend_clipped = metrics.extend_clipped.load(Ordering::Relaxed);
        let modification_tags_removed = metrics.modification_tags_removed.load(Ordering::Relaxed);
        let records_written = progress_counter.load(Ordering::Relaxed);

        info!("Total templates processed: {total_templates}");
        info!("Templates with overlap clipping: {total_overlap_clipped}");
        info!("Templates with mate extension clipping: {total_extend_clipped}");
        if modification_tags_removed > 0 {
            warn!(
                "MM/ML/am/bm were removed from {modification_tags_removed} records whose \
                 modification tags could not be updated to match the clipped SEQ"
            );
        }
        info!("Done!");

        timer.log_completion(records_written);

        Ok(())
    }
}

// ─────────────────────────────────────────────────────────────────────────────
// ClipMetricsFinalizeHook
// ─────────────────────────────────────────────────────────────────────────────

/// Success-only finalize hook: reduces the per-thread `--metrics` accumulator
/// and writes the detailed clipping-metrics TSV. Registered on
/// `finalize_on_success` (not the always-run `finalize` list) so a
/// failed/partial run publishes no stale metrics file — a metrics TSV is only
/// meaningful once every record was processed.
///
/// Each worker's slot holds raw (unfinalized) `fragment`/`read_one`/`read_two`
/// counters (see [`ClippingMetricsCollection::merge`]), which are summed
/// across all slots into one collection and then finalized/written via
/// [`ClippingMetricsCollection::finalize_and_write`], with `finalize` running
/// exactly once, after every slot is merged (accumulate, then finalize once at
/// the end). Because the reduction is pure integer addition, that makes the
/// `pair`/`all` aggregates the TSV reports byte-for-byte reproducible regardless
/// of how work was sharded across threads.
pub(crate) struct ClipMetricsFinalizeHook {
    pub(crate) accumulator: Arc<PerThreadAccumulator<ClippingMetricsCollection>>,
    pub(crate) metrics_path: PathBuf,
}

impl FinalizeHook for ClipMetricsFinalizeHook {
    fn finalize(self: Box<Self>) -> Result<()> {
        let ClipMetricsFinalizeHook { accumulator, metrics_path } = *self;

        let mut merged = ClippingMetricsCollection::new();
        for slot in accumulator.slots() {
            merged.merge(&slot.lock());
        }
        merged.finalize_and_write(&metrics_path)?;
        info!("Wrote metrics to: {}", metrics_path.display());

        Ok(())
    }
}

// ─────────────────────────────────────────────────────────────────────────────
// Step factories — extracted from build_clip_chain (T3a.5).
//
// Each factory receives its captured state as plain arguments and returns the
// concrete step type. Returning `impl Step<...>` is blocked here because the
// closure types embed in opaque-return position and cannot name themselves, so
// we return the concrete `Process*` structs directly (the same pattern used by
// the dedup factories in `chains::commands::dedup`).
// ─────────────────────────────────────────────────────────────────────────────

/// Captures passed into [`build_clip_process_step`] from `add_clip`.
///
/// Bundles all the cloned scalars and Arcs the closure needs so `add_clip`
/// can prepare them once and hand them off cleanly.
///
/// `pub(crate)` — consumed only by `ChainBuilder::add_clip` and
/// [`build_clip_process_step`].
///
/// The clipping control flags are bundled into a single `params: ClipParams`;
/// the hot-path closure delegates the whole per-template decision to
/// `params.clip_template(...)` rather than branching on individual flags.
pub(crate) struct ClipProcessCaptures {
    pub(crate) clipping_mode: crate::clipper::ClippingMode,
    pub(crate) auto_clip_attributes: bool,
    /// The per-template clip configuration built from the parsed `Clip` command.
    pub(crate) params: ClipParams,
    pub(crate) header: noodles::sam::Header,
    pub(crate) reference: Arc<ReferenceReader>,
    pub(crate) metrics: Arc<ClipAtomicMetrics>,
    pub(crate) progress: Arc<AtomicU64>,
    /// Per-thread detailed `--metrics` accumulator. `None` when `--metrics` was
    /// not requested, in which case the hot-path closure passes `None` into
    /// `clip_template` and skips detailed collection entirely (the fast path is
    /// unchanged from before this field existed).
    pub(crate) metrics_accumulator: Option<Arc<PerThreadAccumulator<ClippingMetricsCollection>>>,
}

/// Build the `ClipTemplates` step: parallel, `ByItemOrdinal`. Operates in
/// place on the records of each [`BamTemplateBatch`], regenerating tags and
/// emitting a fresh [`BamTemplateBatch`] carrying the same `batch_serial`
/// so downstream order is preserved.
///
/// `pub(crate)` — consumed only by `ChainBuilder::add_clip`.
#[allow(clippy::type_complexity)]
pub(crate) fn build_clip_process_step(
    limit_bytes: u64,
    cap: ClipProcessCaptures,
) -> ProcessOrdered<
    BamTemplateBatch,
    BamTemplateBatch,
    impl Fn(BamTemplateBatch) -> io::Result<BamTemplateBatch> + Send + Sync + 'static,
> {
    process_ordered::<BamTemplateBatch, BamTemplateBatch, _>(
        "ClipTemplates",
        limit_bytes,
        move |batch: BamTemplateBatch| -> io::Result<BamTemplateBatch> {
            // Per-worker clipper: cheap to construct (no large state).
            let clipper = if cap.auto_clip_attributes {
                RawRecordClipper::with_auto_clip(cap.clipping_mode, true)
            } else {
                RawRecordClipper::new(cap.clipping_mode)
            };

            let (batch_serial, mut templates) = batch.into_parts();

            // When `--metrics` is requested, run the batch under this worker's
            // `PerThreadAccumulator` slot so `clip_template` collects detailed
            // per-read base-clip counts; otherwise pass `None` and stay on the
            // fast, collection-free path.
            let BatchCounts {
                templates: local_templates,
                overlap_clipped: local_overlap_clipped,
                extend_clipped: local_extend_clipped,
                modification_tags_removed: local_modification_tags_removed,
                records: local_record_count,
            } = if let Some(accumulator) = &cap.metrics_accumulator {
                accumulator.with_slot(|slot| {
                    clip_templates_in_batch(
                        &mut templates,
                        &cap.params,
                        &clipper,
                        &cap.header,
                        &cap.reference,
                        Some(slot),
                    )
                })?
            } else {
                clip_templates_in_batch(
                    &mut templates,
                    &cap.params,
                    &clipper,
                    &cap.header,
                    &cap.reference,
                    None,
                )?
            };

            // Aggregate metrics (relaxed atomics, lock-free).
            cap.metrics.total_templates.fetch_add(local_templates, Ordering::Relaxed);
            cap.metrics.overlap_clipped.fetch_add(local_overlap_clipped, Ordering::Relaxed);
            cap.metrics.extend_clipped.fetch_add(local_extend_clipped, Ordering::Relaxed);
            cap.metrics
                .modification_tags_removed
                .fetch_add(local_modification_tags_removed, Ordering::Relaxed);

            // Progress logging (record granularity matches legacy).
            let prev = cap.progress.fetch_add(local_record_count, Ordering::Relaxed);
            if (prev + local_record_count) / 1_000_000 > prev / 1_000_000 {
                info!("Processed {} records", prev + local_record_count);
            }

            // Recompute total_bytes since clipping changed record sizes;
            // BamTemplateBatch::new sums Template::heap_size for us.
            Ok(BamTemplateBatch::new(batch_serial, templates))
        },
    )
}

/// Clips every template in a batch in place, optionally accumulating detailed
/// per-read base-clip counts into `metrics`.
///
/// Shared by both the `--metrics`-enabled and fast (`None`) call sites in
/// [`build_clip_process_step`] so the per-template loop — clip, then
/// regenerate alignment tags — exists exactly once. `metrics` is reborrowed
/// on each iteration via `as_deref_mut` so the same `&mut ClippingMetricsCollection`
/// (the calling worker's `PerThreadAccumulator` slot) accumulates across the
/// whole batch.
///
/// Returns the batch's [`BatchCounts`], for the caller to fold into the chain's
/// atomic summary counters and progress tracker.
fn clip_templates_in_batch(
    templates: &mut [crate::template::Template],
    params: &ClipParams,
    clipper: &RawRecordClipper,
    header: &noodles::sam::Header,
    reference: &ReferenceReader,
    mut metrics: Option<&mut ClippingMetricsCollection>,
) -> io::Result<BatchCounts> {
    use crate::alignment_tags::regenerate_alignment_tags_raw_with_scoring;

    let mut counts = BatchCounts::default();

    for template in templates.iter_mut() {
        counts.templates += 1;
        // Mutate the template's records in place. The earlier
        // version of this closure cloned `template.name` and
        // re-allocated a fresh `Template` per template — that
        // showed up as ~6% extra mimalloc CPU in profiling
        // (mi_page_free_list_extend, mi_free) and pushed the
        // new pipeline ~7% behind legacy at threads=4. Mutating
        // in place keeps allocation count near-equal to legacy.
        let records: &mut Vec<RawRecord> = &mut template.records;

        // Delegate to the canonical per-template clip implementation
        // (`ClipParams::clip_template`). It finds the primary pair by SAM flag (so
        // secondary/supplementary reads are handled, not just records[0]/[1]),
        // applies fixed clipping with "ensure at least N including existing
        // clipping" semantics (not "clip N more"), and repairs mate info. When the
        // caller passed a `--metrics` accumulator slot, detailed per-read base-clip
        // counts are collected here too; otherwise `metrics` is `None` and only the
        // returned per-template flags feed the atomic counters.
        let outcome = params
            .clip_template(records, clipper, metrics.as_deref_mut())
            .map_err(io::Error::other)?;
        counts.overlap_clipped += u64::from(outcome.overlap_clipped);
        counts.extend_clipped += u64::from(outcome.extend_clipped);
        counts.modification_tags_removed += outcome.modification_tags_removed;

        // Regenerate alignment tags for every record (always done to match fgbio). NM/UQ follow
        // the record's SEQ convention, as in `filter`: a simplex methylation consensus keeps the
        // converted bases, so their conversions are not counted.
        for record in records.iter_mut() {
            let scoring = fgumi_consensus::filter::conversion_scoring_for_record(record);
            regenerate_alignment_tags_raw_with_scoring(
                record.as_mut_vec(),
                header,
                reference,
                scoring,
            )
            .map_err(io::Error::other)?;
        }

        counts.records += records.len() as u64;
    }

    Ok(counts)
}

/// Per-batch counts from [`clip_templates_in_batch`].
#[derive(Debug, Default)]
struct BatchCounts {
    templates: u64,
    overlap_clipped: u64,
    extend_clipped: u64,
    modification_tags_removed: u64,
    records: u64,
}

/// Build the `SerializeBamRecords` step for clip: parallel, `ByItemOrdinal`.
/// Serializes each [`BamTemplateBatch`] to raw BAM bytes
/// ([`crate::pipeline::steps::types::DecompressedBlock`]).
///
/// Only included in the chain when the stage is
/// [`StagePosition::Terminal`][`crate::pipeline::chains::builder::StagePosition`];
/// for `Intermediate` the chain tail stays as [`BamTemplateBatch`] for the
/// next stage's input.
///
/// `pub(crate)` — consumed only by `ChainBuilder::add_clip`.
pub(crate) fn build_clip_serialize_step(limit_bytes: u64) -> SerializeBamRecords {
    SerializeBamRecords::new(limit_bytes)
}
