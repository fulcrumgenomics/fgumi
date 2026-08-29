//! `CorrectStep` — typed `Step` wrapping UMI correction.
//!
//! **Ported ahead of its caller.** Every item here is `pub(crate)` -- the
//! config type names `CollectedCorrectMetrics`, which is `pub(crate)`, so this
//! surface cannot be `pub` without leaking command internals. Until the
//! command rewiring points `commands::correct` at this step, nothing outside
//! the module's own tests constructs it, so the whole module is dead code to
//! the compiler. The `allow` below is scoped to this module and should be
//! deleted by the PR that adds the caller.
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

#![allow(dead_code)]

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
    /// Tag the pre-correction sequence is stashed in when original storage
    /// is enabled. Derived from [`crate::commands::correct::Target::original_tag`]
    /// by the caller, alongside `umi_tag`.
    pub(crate) original_tag: [u8; 2],
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
        let mut rejects_bytes: Vec<u8> = Vec::new();
        let kept = run_batch(state, &cfg_run, batch, Some(&mut rejects_bytes))?;
        // X5-001: an all-clean batch must still emit a (zero-byte) rejects block
        // so the rejects branch's `ByItemOrdinal` reorder stage sees a dense
        // serial sequence; pushing nothing (`only_a`) leaves a gap at this
        // serial that wedges the rejects sink. An empty `Vec` does not allocate
        // and produces no physical BGZF block, so the clean and rejecting cases
        // are the same construction.
        let batch_serial = kept.batch_serial();
        Ok(Process2Output::both(kept, DecompressedBlock { batch_serial, bytes: rejects_bytes }))
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

    let body =
        move |state: &mut CorrectWorkerState,
              batch: BamTemplateBatch|
              -> io::Result<BamTemplateBatch> { run_batch(state, &cfg_run, batch, None) };

    process_with_worker_state::<BamTemplateBatch, BamTemplateBatch, _, CorrectWorkerState, _>(
        "correct",
        cfg.output_byte_limit,
        init,
        body,
    )
}

/// Per-batch correction body, shared by both output shapes.
///
/// `rejects` is `Some` only for the 2-output (`--rejects`) variant; the
/// kept-only variant passes `None` and rejected records are simply dropped.
/// One body rather than two so the metrics accounting below has a single
/// definition -- it previously existed in two copies that could drift.
///
/// Metrics follow the legacy `execute` path exactly, which in turn follows
/// fgbio (`CorrectUmis.scala`):
///
/// - `templates_processed` counts every template once, matching legacy's
///   `templates_count = batch.len()`.
/// - Per-UMI crediting is delegated to `CorrectUmis::credit_umi_metrics` and
///   happens for every template that reached matching, *before* the
///   keep/reject decision -- so a template rejected because one segment
///   failed still credits its matched segments, and its unmatched segments
///   credit the all-`N` bucket.
/// - Missing-UMI and wrong-length templates credit **no** per-UMI bucket:
///   fgbio counts them only in `missingUmisRecords` / wrong-length
///   (`CorrectUmis.scala:199-202`), and `compute_template_correction`
///   returns an empty `matches` for `WrongLength`.
fn run_batch(
    state: &mut CorrectWorkerState,
    cfg: &Arc<CorrectStepConfig>,
    batch: BamTemplateBatch,
    mut rejects: Option<&mut Vec<u8>>,
) -> io::Result<BamTemplateBatch> {
    let (batch_serial, templates) = batch.into_parts();

    let mut kept_templates: Vec<Template> = Vec::with_capacity(templates.len());
    let mut local_metrics = CollectedCorrectMetrics::default();
    let mut kept_record_count: u64 = 0;

    for template in templates {
        let mut records = template.into_records();
        let num_records = records.len() as u64;

        // Every template counts once, whatever its outcome.
        local_metrics.templates_processed += 1;

        let umi_opt =
            correct::CorrectUmis::extract_and_validate_template_umi_raw(&records, cfg.umi_tag)
                .map_err(io::Error::other)?;

        let Some(umi) = umi_opt else {
            // No per-UMI credit here: fgbio counts a missing-UMI read only in
            // `missingUmisRecords`.
            local_metrics.missing_umis += num_records;
            if let Some(dst) = rejects.as_deref_mut() {
                for rec in &records {
                    append_framed_raw_record(dst, rec);
                }
            }
            continue;
        };

        let correction = correct::CorrectUmis::compute_template_correction(
            &umi,
            cfg.umi_length,
            cfg.opts.revcomp,
            cfg.opts.max_mismatches,
            cfg.opts.min_distance_diff,
            &cfg.encoded_umi_set,
            &mut state.cache,
        );

        // Credit segments before the keep/reject decision, matching legacy.
        // `matches` is empty for `WrongLength`, so those credit nothing.
        correct::CorrectUmis::credit_umi_metrics(
            &correction.matches,
            num_records,
            &cfg.unmatched_umi,
            &mut local_metrics.umi_matches,
        );

        if correction.matched {
            for raw in &mut records {
                correct::CorrectUmis::apply_correction_to_raw(
                    raw,
                    &correction,
                    cfg.umi_tag,
                    cfg.original_tag,
                    cfg.opts.dont_store_original_umis,
                );
            }
            let corrected = Template::from_records(records)
                .map_err(|e| io::Error::other(format!("Template::from_records: {e}")))?;
            kept_record_count += num_records;
            kept_templates.push(corrected);
        } else {
            match correction.rejection_reason {
                RejectionReason::WrongLength => local_metrics.wrong_length += num_records,
                RejectionReason::Mismatched => local_metrics.mismatched += num_records,
                RejectionReason::None => {}
            }
            if let Some(dst) = rejects.as_deref_mut() {
                for rec in &records {
                    append_framed_raw_record(dst, rec);
                }
            }
        }
    }

    cfg.metrics.with_slot(|m| m.merge_into(&mut local_metrics));
    cfg.records_emitted.fetch_add(kept_record_count, Ordering::Relaxed);

    Ok(BamTemplateBatch::new(batch_serial, kept_templates))
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
