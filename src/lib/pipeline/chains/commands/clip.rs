//! Chain builder for `Stage::Clip`.
//!
//! Phase 2 (T2.15) held the full ~400-LOC chain construction here.
//! Phase 3 (T3a.5) lifts that logic into
//! [`crate::pipeline::chains::builder::ChainBuilder`]; this module
//! now holds the clip-specific types and step factories that the builder
//! imports: [`ClipAtomicMetrics`], [`ClipFinalizeHook`], and the two
//! step-factory functions used by `ChainBuilder::add_clip`.
//!
//! [`build_clip_chain`] shrinks to a ~10-line delegate.

use std::io;
use std::sync::Arc;
use std::sync::atomic::{AtomicU64, Ordering};

use anyhow::Result;
use log::info;

use crate::clipper::RawRecordClipper;
use crate::logging::OperationTimer;
use crate::pipeline::chains::builder::ChainBuilder;
use crate::pipeline::chains::{BuiltPipeline, ChainSpec, FinalizeHook};
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
/// `pub(crate)` so [`ChainBuilder`] can construct and pass it to the step
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
/// `pub(crate)` so [`ChainBuilder`] can construct and register it in
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
/// `pub(crate)` — consumed only by [`ChainBuilder::add_clip`] and
/// [`build_clip_process_step`].
///
/// The five boolean fields reflect the five independent clipping control
/// flags from the `Clip` command struct. Refactoring into a state-machine
/// enum is not warranted here — each flag is checked individually in the
/// hot-path closure.
#[allow(clippy::struct_excessive_bools)]
pub(crate) struct ClipProcessCaptures {
    pub(crate) clipping_mode: crate::clipper::ClippingMode,
    pub(crate) auto_clip_attributes: bool,
    pub(crate) upgrade_clipping: bool,
    pub(crate) clip_overlapping_reads: bool,
    pub(crate) clip_extending_past_mate: bool,
    pub(crate) read_one_five_prime: usize,
    pub(crate) read_one_three_prime: usize,
    pub(crate) read_two_five_prime: usize,
    pub(crate) read_two_three_prime: usize,
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
/// `pub(crate)` — consumed only by [`ChainBuilder::add_clip`].
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
            use crate::commands::clip::update_mate_info_raw;

            // Per-worker clipper: cheap to construct (no large state).
            let clipper = if cap.auto_clip_attributes {
                RawRecordClipper::with_auto_clip(cap.clipping_mode, true)
            } else {
                RawRecordClipper::new(cap.clipping_mode)
            };

            let BamTemplateBatch { batch_serial, mut templates, .. } = batch;
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

                #[allow(clippy::len_zero)]
                if records.len() == 1 {
                    let record = &mut records[0];
                    if cap.upgrade_clipping {
                        clipper.upgrade_all_clipping_raw(record).map_err(io::Error::other)?;
                    }
                    if cap.read_one_five_prime > 0 {
                        clipper.clip_5_prime_end_of_alignment(record, cap.read_one_five_prime);
                    }
                    if cap.read_one_three_prime > 0 {
                        clipper.clip_3_prime_end_of_alignment(record, cap.read_one_three_prime);
                    }
                } else if records.len() == 2 {
                    let (r1_slice, r2_slice) = records.split_at_mut(1);
                    let r1 = &mut r1_slice[0];
                    let r2 = &mut r2_slice[0];

                    if cap.upgrade_clipping {
                        clipper.upgrade_all_clipping_raw(r1).map_err(io::Error::other)?;
                        clipper.upgrade_all_clipping_raw(r2).map_err(io::Error::other)?;
                    }

                    let is_r1_first = r1.is_first_segment();
                    let is_r2_last = r2.is_last_segment();

                    if is_r1_first && cap.read_one_five_prime > 0 {
                        clipper.clip_5_prime_end_of_alignment(r1, cap.read_one_five_prime);
                    } else if !is_r1_first && cap.read_two_five_prime > 0 {
                        clipper.clip_5_prime_end_of_alignment(r1, cap.read_two_five_prime);
                    }
                    if is_r1_first && cap.read_one_three_prime > 0 {
                        clipper.clip_3_prime_end_of_alignment(r1, cap.read_one_three_prime);
                    } else if !is_r1_first && cap.read_two_three_prime > 0 {
                        clipper.clip_3_prime_end_of_alignment(r1, cap.read_two_three_prime);
                    }
                    if is_r2_last && cap.read_two_five_prime > 0 {
                        clipper.clip_5_prime_end_of_alignment(r2, cap.read_two_five_prime);
                    } else if !is_r2_last && cap.read_one_five_prime > 0 {
                        clipper.clip_5_prime_end_of_alignment(r2, cap.read_one_five_prime);
                    }
                    if is_r2_last && cap.read_two_three_prime > 0 {
                        clipper.clip_3_prime_end_of_alignment(r2, cap.read_two_three_prime);
                    } else if !is_r2_last && cap.read_one_three_prime > 0 {
                        clipper.clip_3_prime_end_of_alignment(r2, cap.read_one_three_prime);
                    }

                    if cap.clip_overlapping_reads {
                        let (n1, n2) = clipper.clip_overlapping_reads(r1, r2);
                        if n1 > 0 || n2 > 0 {
                            local_overlap_clipped += 1;
                        }
                    }
                    if cap.clip_extending_past_mate {
                        let (n1, n2) = clipper.clip_extending_past_mate_ends(r1, r2);
                        if n1 > 0 || n2 > 0 {
                            local_extend_clipped += 1;
                        }
                    }

                    update_mate_info_raw(r1, r2);
                    update_mate_info_raw(r2, r1);
                }

                // Regenerate alignment tags for every record.
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
/// `pub(crate)` — consumed only by [`ChainBuilder::add_clip`].
pub(crate) fn build_clip_serialize_step(limit_bytes: u64) -> SerializeBamRecords {
    SerializeBamRecords::new(limit_bytes)
}

// ─────────────────────────────────────────────────────────────────────────────
// build_clip_chain — thin delegate (Phase 3 T3a.5)
// ─────────────────────────────────────────────────────────────────────────────

/// Build a [`BuiltPipeline`] for a clip-only chain.
///
/// Delegates to [`ChainBuilder`]. The full step sequence lives in
/// [`ChainBuilder`]'s `add_clip` method in
/// [`crate::pipeline::chains::builder`].
///
/// # Errors
///
/// Returns input-validation errors (reference FASTA missing, no clipping
/// operation requested, `--metrics` + `--threads` conflict, etc.), or any
/// underlying pipeline construction error.
#[allow(clippy::needless_pass_by_value)]
pub fn build_clip_chain(spec: ChainSpec) -> Result<BuiltPipeline> {
    use crate::pipeline::chains::{Stage, builder::StagePosition};
    let mut chain = ChainBuilder::new(&spec)?;
    chain.add_source()?;
    chain.add_stage(Stage::Clip, StagePosition::Terminal)?;
    chain.add_sink()?;
    chain.build()
}
