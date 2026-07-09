//! `CorrectStep` — typed `Step` wrapping UMI correction.
//!
//! Input:  `BamTemplateBatch` (from `GroupByQueryname`).
//!
//! Two factory functions cover the two `--rejects` cases:
//!
//! - `correct_step_with_rejects(cfg)` — `Step` with
//!   `Outputs = OrderedBytesTuple2<BamTemplateBatch, DecompressedBlock>`.
//!   `multi.b0` is the corrected templates (typed), `multi.b1` is the
//!   pre-serialized rejects bytes. Used when `track_rejects = true`.
//!
//! - `correct_step_kept_only(cfg)` — `Step` with
//!   `Outputs = OrderedBytesSingle<BamTemplateBatch>`. Single-output
//!   variant when no rejects file is requested. The framework has no
//!   public `DiscardSink`, and `Pipeline::build()` rejects an unwired
//!   branch, so this is the only way to omit the rejects branch.
//!
//! Per-worker `LruCache<Vec<u8>, UmiMatch>` lives in the framework's
//! `Process2WithWorkerState` (or `ProcessWithWorkerState`) lazy-init
//! slot. Wrapped in `CorrectWorkerState` to satisfy the `S: HeapSize`
//! bound (`LruCache` doesn't impl `HeapSize`; existing precedents like
//! `ConsensusState` use the same `impl HeapSize for State {}` pattern).

use std::io;
use std::num::NonZero;
use std::sync::Arc;
use std::sync::atomic::{AtomicU64, Ordering};

use lru::LruCache;

use crate::commands::correct::{
    self, CollectedCorrectMetrics, CorrectOptions, EncodedUmiSet, RejectionReason, UmiMatch,
};
use crate::per_thread_accumulator::PerThreadAccumulator;
use crate::pipeline::core::item::HeapSize;
use crate::pipeline::core::outputs::{OrderedBytesSingle, OrderedBytesTuple2};
use crate::pipeline::core::step::Step;
use crate::pipeline::steps::process::{
    Process2Output, process_with_worker_state, process2_with_worker_state,
};
use crate::pipeline::steps::types::{BamTemplateBatch, DecompressedBlock};
use crate::template::Template;
use fgumi_raw_bam::RawRecord;

#[cfg(test)]
mod tests;

/// Per-worker scratch state for `CorrectStep`. Wraps `Option<LruCache>`
/// in a newtype so we can `impl HeapSize` for it (the bound on
/// `process2_with_worker_state`'s `S` type parameter). `LruCache` itself
/// can't impl `HeapSize` (orphan rule). Reports `heap_size = 0` because
/// the cache size is bounded by the user-set `cache_size` option, not
/// dynamically grown — same convention as `ConsensusState` /
/// `CodecState` / `DuplexState` in `commands::runall`.
pub(crate) struct CorrectWorkerState {
    cache: Option<LruCache<Vec<u8>, UmiMatch>>,
}

impl HeapSize for CorrectWorkerState {
    fn heap_size(&self) -> usize {
        0
    }
}

/// Immutable configuration shared by `CorrectStep`. All `Arc`-wrapped
/// fields are cheap to clone across workers; per-worker mutable state
/// (the cache) lives in the framework's lazy-init slot, not here.
///
/// `pub(crate)` because the type references `CollectedCorrectMetrics`,
/// which is `pub(crate)`; the new pipeline is wired entirely from
/// `commands::correct::execute_new_pipeline`, so this never needs to
/// cross the crate boundary.
pub(crate) struct CorrectStepConfig {
    pub(crate) encoded_umi_set: Arc<EncodedUmiSet>,
    pub(crate) umi_length: usize,
    pub(crate) umi_tag: [u8; 2],
    pub(crate) opts: CorrectOptions,
    pub(crate) metrics: Arc<PerThreadAccumulator<CollectedCorrectMetrics>>,
    pub(crate) records_emitted: Arc<AtomicU64>,
    pub(crate) output_byte_limit: u64,
    pub(crate) unmatched_umi: String,
}

/// Build the 2-output `CorrectStep`. Used when `--rejects` is set; kept
/// branch carries `BamTemplateBatch`, rejects branch carries
/// pre-serialized BAM bytes inside a `DecompressedBlock`.
pub(crate) fn correct_step_with_rejects(
    cfg: CorrectStepConfig,
) -> impl Step<Input = BamTemplateBatch, Outputs = OrderedBytesTuple2<BamTemplateBatch, DecompressedBlock>>
{
    let cfg = Arc::new(cfg);
    let cfg_init = Arc::clone(&cfg);
    let cfg_run = Arc::clone(&cfg);

    let init = move || CorrectWorkerState {
        cache: if cfg_init.opts.cache_size > 0 {
            Some(LruCache::new(
                NonZero::new(cfg_init.opts.cache_size).expect("cache_size > 0 checked above"),
            ))
        } else {
            None
        },
    };

    let body = move |state: &mut CorrectWorkerState,
                     batch: BamTemplateBatch|
          -> io::Result<Process2Output<BamTemplateBatch, DecompressedBlock>> {
        let (kept, rejects_opt) = run_batch_with_rejects(state, &cfg_run, batch)?;
        if let Some(rejects) = rejects_opt {
            Ok(Process2Output::both(kept, rejects))
        } else {
            // X5-001: an all-clean batch must still emit a (zero-byte) rejects
            // block so the rejects branch's `ByItemOrdinal` reorder stage sees a
            // dense serial sequence; pushing nothing (`only_a`) leaves a gap at
            // this serial that wedges the rejects sink. The empty `Vec` does not
            // allocate and produces no physical BGZF block.
            let batch_serial = kept.batch_serial();
            Ok(Process2Output::both(kept, DecompressedBlock { batch_serial, bytes: Vec::new() }))
        }
    };

    process2_with_worker_state::<
        BamTemplateBatch,
        BamTemplateBatch,
        DecompressedBlock,
        CorrectWorkerState,
        _,
        _,
    >("correct", cfg.output_byte_limit, cfg.output_byte_limit, init, body)
}

/// Build the 1-output (kept-only) `CorrectStep`. Used when `--rejects`
/// is unset. The framework has no public `DiscardSink`, so the only way
/// to omit the rejects branch is to omit it from the chain entirely via
/// this single-output variant.
pub(crate) fn correct_step_kept_only(
    cfg: CorrectStepConfig,
) -> impl Step<Input = BamTemplateBatch, Outputs = OrderedBytesSingle<BamTemplateBatch>> {
    let cfg = Arc::new(cfg);
    let cfg_init = Arc::clone(&cfg);
    let cfg_run = Arc::clone(&cfg);

    let init = move || CorrectWorkerState {
        cache: if cfg_init.opts.cache_size > 0 {
            Some(LruCache::new(
                NonZero::new(cfg_init.opts.cache_size).expect("cache_size > 0 checked above"),
            ))
        } else {
            None
        },
    };

    let body = move |state: &mut CorrectWorkerState,
                     batch: BamTemplateBatch|
          -> io::Result<BamTemplateBatch> {
        let (kept, _no_rejects_collected) = run_batch_kept_only(state, &cfg_run, batch)?;
        Ok(kept)
    };

    process_with_worker_state::<BamTemplateBatch, BamTemplateBatch, _, CorrectWorkerState, _>(
        "correct",
        cfg.output_byte_limit,
        init,
        body,
    )
}

/// Per-batch correction body — 2-output variant. Returns the corrected
/// `BamTemplateBatch` (kept) and an optional rejects `DecompressedBlock`.
///
/// NOTE: `templates_processed` is incremented for ALL templates
/// (matched + unmatched + missing-RX), matching legacy
/// `templates_count = batch.len()` semantics. `merge_match` bumps for
/// matched; the unmatched/missing arms bump explicitly here.
fn run_batch_with_rejects(
    state: &mut CorrectWorkerState,
    cfg: &Arc<CorrectStepConfig>,
    batch: BamTemplateBatch,
) -> io::Result<(BamTemplateBatch, Option<DecompressedBlock>)> {
    let (batch_serial, templates) = batch.into_parts();

    let mut kept_templates: Vec<Template> = Vec::with_capacity(templates.len());
    let mut rejects_bytes: Vec<u8> = Vec::new();
    let mut local_metrics = CollectedCorrectMetrics::default();
    let mut kept_record_count: u64 = 0;

    for template in templates {
        let mut records = template.into_records();
        let num_records = records.len() as u64;

        let umi_opt =
            correct::CorrectUmis::extract_and_validate_template_umi_raw(&records, cfg.umi_tag)
                .map_err(io::Error::other)?;

        match umi_opt {
            None => {
                local_metrics.templates_processed += 1;
                local_metrics.missing_umis += num_records;
                local_metrics.merge_unmatched(&cfg.unmatched_umi, num_records);
                for rec in &records {
                    append_framed_raw_record(&mut rejects_bytes, rec);
                }
            }
            Some(umi) => {
                let correction = correct::CorrectUmis::compute_template_correction(
                    &umi,
                    cfg.umi_length,
                    cfg.opts.revcomp,
                    cfg.opts.max_mismatches,
                    cfg.opts.min_distance_diff,
                    &cfg.encoded_umi_set,
                    &mut state.cache,
                );
                if correction.matched {
                    local_metrics.merge_match(&correction, num_records);
                    for raw in &mut records {
                        correct::CorrectUmis::apply_correction_to_raw(
                            raw,
                            &correction,
                            cfg.umi_tag,
                            cfg.opts.dont_store_original_umis,
                        );
                    }
                    let corrected = Template::from_records(records)
                        .map_err(|e| io::Error::other(format!("Template::from_records: {e}")))?;
                    kept_record_count += num_records;
                    kept_templates.push(corrected);
                } else {
                    local_metrics.templates_processed += 1;
                    match correction.rejection_reason {
                        RejectionReason::WrongLength => local_metrics.wrong_length += num_records,
                        RejectionReason::Mismatched => local_metrics.mismatched += num_records,
                        RejectionReason::None => {}
                    }
                    local_metrics.merge_unmatched(&cfg.unmatched_umi, num_records);
                    for rec in &records {
                        append_framed_raw_record(&mut rejects_bytes, rec);
                    }
                }
            }
        }
    }

    cfg.metrics.with_slot(|m| m.merge_into(&mut local_metrics));
    cfg.records_emitted.fetch_add(kept_record_count, Ordering::Relaxed);

    let kept = BamTemplateBatch::new(batch_serial, kept_templates);
    let rejects = if rejects_bytes.is_empty() {
        None
    } else {
        Some(DecompressedBlock { batch_serial, bytes: rejects_bytes })
    };
    Ok((kept, rejects))
}

/// Per-batch correction body — kept-only variant. Same logic as
/// `run_batch_with_rejects` but doesn't accumulate rejects bytes.
/// Metrics still accumulate.
fn run_batch_kept_only(
    state: &mut CorrectWorkerState,
    cfg: &Arc<CorrectStepConfig>,
    batch: BamTemplateBatch,
) -> io::Result<(BamTemplateBatch, ())> {
    let (batch_serial, templates) = batch.into_parts();

    let mut kept_templates: Vec<Template> = Vec::with_capacity(templates.len());
    let mut local_metrics = CollectedCorrectMetrics::default();
    let mut kept_record_count: u64 = 0;

    for template in templates {
        let mut records = template.into_records();
        let num_records = records.len() as u64;

        let umi_opt =
            correct::CorrectUmis::extract_and_validate_template_umi_raw(&records, cfg.umi_tag)
                .map_err(io::Error::other)?;

        match umi_opt {
            None => {
                local_metrics.templates_processed += 1;
                local_metrics.missing_umis += num_records;
                local_metrics.merge_unmatched(&cfg.unmatched_umi, num_records);
            }
            Some(umi) => {
                let correction = correct::CorrectUmis::compute_template_correction(
                    &umi,
                    cfg.umi_length,
                    cfg.opts.revcomp,
                    cfg.opts.max_mismatches,
                    cfg.opts.min_distance_diff,
                    &cfg.encoded_umi_set,
                    &mut state.cache,
                );
                if correction.matched {
                    local_metrics.merge_match(&correction, num_records);
                    for raw in &mut records {
                        correct::CorrectUmis::apply_correction_to_raw(
                            raw,
                            &correction,
                            cfg.umi_tag,
                            cfg.opts.dont_store_original_umis,
                        );
                    }
                    let corrected = Template::from_records(records)
                        .map_err(|e| io::Error::other(format!("Template::from_records: {e}")))?;
                    kept_record_count += num_records;
                    kept_templates.push(corrected);
                } else {
                    local_metrics.templates_processed += 1;
                    match correction.rejection_reason {
                        RejectionReason::WrongLength => local_metrics.wrong_length += num_records,
                        RejectionReason::Mismatched => local_metrics.mismatched += num_records,
                        RejectionReason::None => {}
                    }
                    local_metrics.merge_unmatched(&cfg.unmatched_umi, num_records);
                }
            }
        }
    }

    cfg.metrics.with_slot(|m| m.merge_into(&mut local_metrics));
    cfg.records_emitted.fetch_add(kept_record_count, Ordering::Relaxed);

    let kept = BamTemplateBatch::new(batch_serial, kept_templates);
    Ok((kept, ()))
}

/// Append one raw BAM record to `dst` using the standard BAM framing:
/// 4-byte LE `block_size` followed by the record body. Matches the
/// legacy rejects writer's framing (the `serialize_fn` body that frames
/// kept records) and what `BgzfCompress` expects in a `DecompressedBlock`.
fn append_framed_raw_record(dst: &mut Vec<u8>, rec: &RawRecord) {
    // BAM record body size is u32-bounded per the spec; the underlying
    // RawRecord buffer was sized to fit a single record, so the cast
    // cannot truncate in practice. Matches the legacy serialize_fn cast.
    #[allow(clippy::cast_possible_truncation)]
    let block_size = rec.len() as u32;
    dst.extend_from_slice(&block_size.to_le_bytes());
    dst.extend_from_slice(rec);
}
