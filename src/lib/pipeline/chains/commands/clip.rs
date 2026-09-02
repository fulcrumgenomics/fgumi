//! Chain builder for `Stage::Clip`.
//!
//! Phase 2 (T2.15) held the full ~400-LOC chain construction here.
//! Phase 3 (T3a.5) lifts that logic into
//! [`crate::pipeline::chains::builder::ChainBuilder`]; this module
//! now holds the clip-specific types and step factories that the builder
//! imports: `ClipAtomicMetrics`, `ClipFinalizeHook`, and the two
//! step-factory functions used by `ChainBuilder::add_clip`.

use std::io;
use std::sync::Arc;
use std::sync::atomic::{AtomicU64, Ordering};

use anyhow::Result;
use log::info;

use crate::clipper::RawRecordClipper;
use crate::logging::OperationTimer;
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
}

impl Default for ClipAtomicMetrics {
    fn default() -> Self {
        Self {
            total_templates: AtomicU64::new(0),
            overlap_clipped: AtomicU64::new(0),
            extend_clipped: AtomicU64::new(0),
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
        let records_written = progress_counter.load(Ordering::Relaxed);

        info!("Total templates processed: {total_templates}");
        info!("Templates with overlap clipping: {total_overlap_clipped}");
        info!("Templates with mate extension clipping: {total_extend_clipped}");
        info!("Done!");

        timer.log_completion(records_written);

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
/// The clipping control flags are bundled into a single `params: ClipParams`,
/// shared verbatim with the single-threaded oracle (`Clip::execute_single_threaded`);
/// the hot-path closure delegates the whole per-template decision to
/// `params.clip_template(...)` rather than branching on individual flags.
pub(crate) struct ClipProcessCaptures {
    pub(crate) clipping_mode: crate::clipper::ClippingMode,
    pub(crate) auto_clip_attributes: bool,
    /// The per-template clip configuration, shared verbatim with the
    /// single-threaded `Clip::execute` path (its in-process parity oracle).
    pub(crate) params: crate::commands::clip::ClipParams,
    pub(crate) header: noodles::sam::Header,
    pub(crate) reference: Arc<ReferenceReader>,
    pub(crate) metrics: Arc<ClipAtomicMetrics>,
    pub(crate) progress: Arc<AtomicU64>,
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
            use crate::alignment_tags::regenerate_alignment_tags_raw;

            // Per-worker clipper: cheap to construct (no large state).
            let clipper = if cap.auto_clip_attributes {
                RawRecordClipper::with_auto_clip(cap.clipping_mode, true)
            } else {
                RawRecordClipper::new(cap.clipping_mode)
            };

            let (batch_serial, mut templates) = batch.into_parts();
            let mut local_templates: u64 = 0;
            let mut local_overlap_clipped: u64 = 0;
            let mut local_extend_clipped: u64 = 0;
            let mut local_record_count: u64 = 0;

            for template in &mut templates {
                local_templates += 1;
                // Mutate the template's records in place. The earlier
                // version of this closure cloned `template.name` and
                // re-allocated a fresh `Template` per template — that
                // showed up as ~6% extra mimalloc CPU in profiling
                // (mi_page_free_list_extend, mi_free) and pushed the
                // new pipeline ~7% behind legacy at threads=4. Mutating
                // in place keeps allocation count near-equal to legacy.
                let records: &mut Vec<RawRecord> = &mut template.records;

                // Delegate to the canonical per-template clip implementation —
                // the exact code the single-threaded `Clip::execute` oracle runs
                // (`ClipParams::clip_template`). It finds the primary pair by SAM
                // flag (so secondary/supplementary reads are handled, not just
                // records[0]/[1]), applies fixed clipping with "ensure at least N
                // including existing clipping" semantics (not "clip N more"), and
                // repairs mate info. The `--metrics` detailed collection is never
                // produced under `--threads` (rejected up front), so pass `None`
                // and use the returned per-template flags for the atomic counters.
                let (overlap_clipped, extend_clipped) =
                    cap.params.clip_template(records, &clipper, None).map_err(io::Error::other)?;
                if overlap_clipped {
                    local_overlap_clipped += 1;
                }
                if extend_clipped {
                    local_extend_clipped += 1;
                }

                // Regenerate alignment tags for every record (matches the oracle).
                for record in records.iter_mut() {
                    regenerate_alignment_tags_raw(record.as_mut_vec(), &cap.header, &cap.reference)
                        .map_err(io::Error::other)?;
                }

                local_record_count += records.len() as u64;
            }

            // Aggregate metrics (relaxed atomics, lock-free).
            cap.metrics.total_templates.fetch_add(local_templates, Ordering::Relaxed);
            cap.metrics.overlap_clipped.fetch_add(local_overlap_clipped, Ordering::Relaxed);
            cap.metrics.extend_clipped.fetch_add(local_extend_clipped, Ordering::Relaxed);

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
