use std::collections::VecDeque;
use std::sync::Mutex;

use super::*;
use crate::commands::correct::{CollectedCorrectMetrics, CorrectOptions, Target};
use crate::pipeline::core::Unpushed;
use crate::pipeline::core::builder::{Pipeline, PipelineConfig};
use crate::pipeline::core::held::HeldSlot;
use crate::pipeline::core::item::Ordered;
use crate::pipeline::core::queues::QueueSpec;
use crate::pipeline::core::reorder::BranchOrdering;
use crate::pipeline::core::step::{StepCtx, StepKind, StepOutcome, StepProfile};
use crate::sam::SamTag;
use crate::template::Template;
use fgumi_raw_bam::SamBuilder;
use rstest::rstest;

/// `CorrectOptions` mirroring the CLI's documented defaults.
///
/// Spelled out rather than derived: a `Default` impl on `CorrectOptions` would
/// hand back `max_mismatches: 0`, `min_distance_diff: 0` and `cache_size: 0`,
/// none of which match the flags' `default_value`s, so it would be a trap for
/// any non-test caller that reached for it.
fn make_default_opts() -> CorrectOptions {
    CorrectOptions {
        metrics: None,
        target: Target::Umi,
        max_mismatches: 2,
        min_distance_diff: 2,
        umis: Vec::new(),
        umi_files: Vec::new(),
        dont_store_original_umis: false,
        cache_size: 100_000,
        min_corrected: None,
        revcomp: false,
        rejects_path: None,
    }
}

fn make_default_cfg() -> CorrectStepConfig {
    let opts = make_default_opts();
    CorrectStepConfig {
        encoded_umi_set: Arc::new(EncodedUmiSet::new(&["AAAA".to_string()])),
        umi_length: 4,
        umi_tag: opts.target.sequence_tag().into(),
        original_tag: opts.target.original_tag().into(),
        opts,
        metrics: PerThreadAccumulator::new(1),
        records_emitted: Arc::new(AtomicU64::new(0)),
        output_byte_limit: 1024,
        unmatched_umi: "NNNN".to_string(),
    }
}

/// Both factories must present the same step identity and kind, and differ
/// only in how many output queues they declare: the rejects variant fans out
/// to two branches, the kept-only variant to one. Every queue is byte-bounded
/// in both -- correction emits variable-size batches, so a count-bounded queue
/// would not bound memory.
#[rstest]
#[case::with_rejects(correct_step_with_rejects(make_default_cfg()).profile(), 2)]
#[case::kept_only(correct_step_kept_only(make_default_cfg()).profile(), 1)]
fn correct_step_profile(#[case] profile: StepProfile, #[case] expected_output_queues: usize) {
    assert_eq!(profile.name, "correct");
    assert!(matches!(profile.kind, StepKind::Parallel));
    assert_eq!(profile.output_queues.len(), expected_output_queues);
    // Every output queue must be byte-bounded at exactly the config's
    // `output_byte_limit`. Accepting any `ByteBounded` value would let a factory
    // that ignored `output_byte_limit` (or wired a different bound) pass.
    let expected_limit = make_default_cfg().output_byte_limit;
    for q in &profile.output_queues {
        let QueueSpec::ByteBounded { limit_bytes } = q else {
            panic!("correct output queues must be byte-bounded");
        };
        assert_eq!(*limit_bytes, expected_limit, "queue byte bound must equal output_byte_limit");
    }
}

// --- Batch-routing tests for `run_batch` (both output shapes) ---

/// Mirror of the 2-output factory's body: run a batch collecting rejects and
/// frame them exactly as the step does.
///
/// The step emits a `DecompressedBlock` unconditionally — an all-clean batch
/// still yields a zero-byte block so the rejects branch's `ByItemOrdinal`
/// reorder stage sees a dense serial sequence (X5-001). This helper does the
/// same, so a test cannot assert a shape the step never produces.
fn run_batch_collecting_rejects(
    state: &mut CorrectWorkerState,
    cfg: &Arc<CorrectStepConfig>,
    batch: BamTemplateBatch,
) -> (BamTemplateBatch, DecompressedBlock) {
    let mut bytes: Vec<u8> = Vec::new();
    let kept = run_batch(state, cfg, batch, Some(&mut bytes)).expect("run_batch with rejects");
    let batch_serial = kept.batch_serial();
    (kept, DecompressedBlock { batch_serial, bytes })
}

/// Build a fresh per-worker state matching how the step's `init` closure
/// constructs it for the default (cache-enabled) `CorrectOptions`.
fn make_state() -> CorrectWorkerState {
    CorrectWorkerState {
        cache: Some(LruCache::new(
            NonZero::new(make_default_opts().cache_size).expect("default cache_size > 0"),
        )),
    }
}

/// Build a single-record (unpaired) template whose `RX` tag holds `umi`.
///
/// A 4-base sequence/quality keeps the record comfortably above the 32-byte
/// minimum that `extract_and_validate_template_umi_raw` / `Template::from_records`
/// require, and a distinct `qname` keeps each template independent.
fn make_template_with_rx(qname: &[u8], umi: &str) -> Template {
    let mut b = SamBuilder::new();
    b.read_name(qname)
        .flags(0) // unmapped, unpaired -> single-record template (R1)
        .sequence(b"ACGT")
        .qualities(b"IIII")
        .add_string_tag(SamTag::RX, umi.as_bytes());
    Template::from_records(vec![b.build()]).expect("single-record template")
}

/// Read the `RX` tag value back from the first record of a template.
fn template_rx(template: &Template) -> String {
    let rx = fgumi_raw_bam::find_string_tag_in_record(&template.records[0], SamTag::RX)
        .expect("RX tag present on kept record");
    String::from_utf8(rx.to_vec()).expect("RX is valid UTF-8")
}

/// Read an arbitrary string tag from the first record of a template, or `None`
/// when the tag is absent — used to assert both presence and absence of the
/// original-UMI tag (`OX`) after correction.
fn template_tag(template: &Template, tag: SamTag) -> Option<String> {
    fgumi_raw_bam::find_string_tag_in_record(&template.records[0], tag)
        .map(|bytes| String::from_utf8(bytes.to_vec()).expect("tag is valid UTF-8"))
}

/// Read the `QNAME` (read name) back from the first record of a template. Used
/// to pin *which* templates were routed to a branch by identity — the corrected
/// `RX` collapses both kept templates to `AAAA`, so an RX-only check cannot tell
/// the intended survivor from a regression that kept the wrong record.
fn template_qname(template: &Template) -> String {
    let name = fgumi_raw_bam::read_name(template.records[0].as_ref());
    String::from_utf8(name.to_vec()).expect("qname is valid UTF-8")
}

/// Decode a BAM-framed reject block — `4-byte LE block_size` + record body,
/// repeated — into the `RX` tag of each rejected record. This is an
/// independent inverse of the production `append_framed_raw_record` writer, so
/// a routing regression that serialized the wrong record can't pass a mere
/// non-empty-bytes check.
fn reject_block_rx_tags(bytes: &[u8]) -> Vec<String> {
    let mut tags = Vec::new();
    let mut offset = 0;
    while offset < bytes.len() {
        let block_size = u32::from_le_bytes(
            bytes[offset..offset + 4].try_into().expect("4-byte block_size prefix"),
        ) as usize;
        offset += 4;
        let body = &bytes[offset..offset + block_size];
        let rx = fgumi_raw_bam::find_string_tag_in_record(body, SamTag::RX)
            .expect("RX tag present on rejected record");
        tags.push(String::from_utf8(rx.to_vec()).expect("RX is valid UTF-8"));
        offset += block_size;
    }
    tags
}

/// With the default config (whitelist `["AAAA"]`, `max_mismatches = 2`):
/// - `RX = "AAAA"` matches exactly -> KEPT, unchanged.
/// - `RX = "AAAT"` is 1 mismatch from `AAAA` (<= 2) -> KEPT, corrected to `AAAA`.
/// - `RX = "TTTT"` is 4 mismatches from `AAAA` (> 2) -> REJECTED.
#[test]
fn run_batch_with_rejects_splits_kept_and_rejected() {
    let cfg = Arc::new(make_default_cfg());
    let mut state = make_state();

    let templates = vec![
        make_template_with_rx(b"exact", "AAAA"),
        make_template_with_rx(b"correctable", "AAAT"),
        make_template_with_rx(b"offlist", "TTTT"),
    ];
    let batch = BamTemplateBatch::new(7, templates);

    let (kept, rejects) = run_batch_collecting_rejects(&mut state, &cfg, batch);

    // Kept branch keeps the on-whitelist + correctable templates only.
    assert_eq!(kept.batch_serial(), 7);
    assert_eq!(kept.templates().len(), 2, "exact + correctable kept");
    let mut kept_umis: Vec<String> = kept.templates().iter().map(template_rx).collect();
    kept_umis.sort();
    // Both kept records carry the corrected/whitelisted UMI "AAAA".
    assert_eq!(kept_umis, vec!["AAAA".to_string(), "AAAA".to_string()]);
    // Pin kept *identity* by QNAME, not just the corrected RX: both kept
    // templates normalize to "AAAA", so a regression that kept the wrong record
    // (e.g. dropped `correctable` and kept `offlist` rewritten to AAAA) would
    // still satisfy the RX check above.
    let mut kept_qnames: Vec<String> = kept.templates().iter().map(template_qname).collect();
    kept_qnames.sort();
    assert_eq!(kept_qnames, vec!["correctable".to_string(), "exact".to_string()]);

    // Rejects branch carries exactly the off-whitelist template. Decode the
    // framed block and assert the rejected record's identity (its untouched
    // `RX = "TTTT"`, 4 mismatches from the `AAAA` whitelist) — a non-empty
    // bytes check alone would pass even if the wrong record were routed here.
    assert_eq!(rejects.batch_serial, 7);
    assert_eq!(
        reject_block_rx_tags(&rejects.bytes),
        vec!["TTTT".to_string()],
        "only the off-whitelist TTTT record is routed to rejects, unchanged"
    );

    // records_emitted counts only the two kept single-record templates.
    assert_eq!(cfg.records_emitted.load(Ordering::Relaxed), 2);

    // Metrics: one rejected (mismatched) template, two kept.
    let slot = &cfg.metrics.slots()[0];
    let metrics = slot.lock();
    assert_eq!(metrics.mismatched, 1, "one off-whitelist template rejected");
}

/// `run_batch` forwards `cfg.original_tag` and `cfg.opts.dont_store_original_umis`
/// to `apply_correction_to_raw`, which stores the pre-correction UMI under the
/// original tag (`OX` for `Target::Umi`) whenever a template is actually
/// corrected — unless suppression is requested. Drive a corrected template
/// (`AAAT` -> `AAAA`, one mismatch) through `run_batch` and assert the original
/// `AAAT` lands under `OX` when storage is enabled and is absent when
/// `dont_store_original_umis` is set. The corrected `RX` is `AAAA` either way.
#[rstest]
#[case::store_original(false, Some("AAAT"))]
#[case::suppress_original(true, None)]
fn run_batch_stores_or_suppresses_original_umi(
    #[case] dont_store_original_umis: bool,
    #[case] expected_ox: Option<&str>,
) {
    let mut cfg = make_default_cfg();
    cfg.opts.dont_store_original_umis = dont_store_original_umis;
    let cfg = Arc::new(cfg);
    let mut state = make_state();

    let batch = BamTemplateBatch::new(0, vec![make_template_with_rx(b"correctable", "AAAT")]);
    let kept = run_batch(&mut state, &cfg, batch, None).expect("run_batch");

    assert_eq!(kept.templates().len(), 1, "the correctable template is kept");
    let template = &kept.templates()[0];
    // The corrected UMI is written to RX regardless of original-UMI storage.
    assert_eq!(template_rx(template), "AAAA", "RX is corrected to the whitelist UMI");
    // The pre-correction UMI is stored under OX only when storage is enabled.
    assert_eq!(
        template_tag(template, SamTag::OX),
        expected_ox.map(str::to_string),
        "OX carries the original UMI iff dont_store_original_umis is false",
    );
}

/// An all-on-whitelist batch produces no rejects (`None`) from the
/// 2-output path, so the framework emits a zero-byte rejects block instead.
/// This drives `run_batch` directly; its end-to-end companion,
/// `correct_step_with_rejects_pipeline_emits_empty_rejects_for_clean_input`,
/// proves the same invariant holds when the real pipeline drives the factory.
#[test]
fn run_batch_with_rejects_no_rejects_when_all_kept() {
    let cfg = Arc::new(make_default_cfg());
    let mut state = make_state();

    let batch = BamTemplateBatch::new(3, vec![make_template_with_rx(b"exact", "AAAA")]);
    let (kept, rejects) = run_batch_collecting_rejects(&mut state, &cfg, batch);

    assert_eq!(kept.templates().len(), 1);
    // X5-001: the block is still emitted, just empty — the rejects branch needs a
    // dense serial sequence, so "no rejects" is a zero-byte block, not a gap.
    assert_eq!(rejects.batch_serial, 3);
    assert!(rejects.bytes.is_empty(), "no rejected records -> zero-byte rejects block");
    assert_eq!(cfg.records_emitted.load(Ordering::Relaxed), 1);
}

/// Build a single-record template with **no** `RX` tag — the missing-UMI case.
fn make_template_without_rx(qname: &[u8]) -> Template {
    let mut b = SamBuilder::new();
    b.read_name(qname).flags(0).sequence(b"ACGT").qualities(b"IIII");
    Template::from_records(vec![b.build()]).expect("single-record template")
}

/// Total credited to one per-UMI bucket, or 0 when the bucket is absent.
fn bucket_total(metrics: &CollectedCorrectMetrics, umi: &str) -> u64 {
    metrics.umi_matches.get(umi).map_or(0, |m| m.total_matches)
}

/// Per-UMI crediting must match fgbio, which the legacy `execute` path
/// implements via `CorrectUmis::credit_umi_metrics`. These are the three cases
/// where a whole-template "unmatched" credit diverges from fgbio's per-segment
/// accounting, so each is pinned independently.
///
/// - A **missing-UMI** template is counted only in `missing_umis`; fgbio never
///   credits the all-`N` bucket for it (`CorrectUmis.scala:199-202`).
/// - A **wrong-length** UMI never reaches matching, so it credits no bucket.
/// - A **mismatched** template still credits its matched segments, and credits
///   the all-`N` bucket only for the segments that actually failed. The
///   `AAAA-TTTT` case is the discriminating one: the template is rejected on
///   its second segment, yet the first still credits `AAAA`. Crediting the
///   whole template to the all-`N` bucket would leave `AAAA` at zero.
#[rstest]
#[case::missing_umi_credits_no_bucket(None, 1, 0, 0, 0, 0)]
#[case::wrong_length_credits_no_bucket(Some("AA"), 0, 1, 0, 0, 0)]
#[case::offlist_credits_only_the_n_bucket(Some("TTTT"), 0, 0, 1, 0, 1)]
#[case::partial_match_credits_both_segments(Some("AAAA-TTTT"), 0, 0, 1, 1, 1)]
#[case::onlist_credits_its_own_bucket(Some("AAAA"), 0, 0, 0, 1, 0)]
fn per_umi_crediting_matches_fgbio(
    #[case] rx: Option<&str>,
    #[case] expected_missing: u64,
    #[case] expected_wrong_length: u64,
    #[case] expected_mismatched: u64,
    #[case] expected_aaaa_bucket: u64,
    #[case] expected_n_bucket: u64,
) {
    let cfg = Arc::new(make_default_cfg());
    let mut state = make_state();

    let template = match rx {
        Some(umi) => make_template_with_rx(b"t", umi),
        None => make_template_without_rx(b"t"),
    };
    let batch = BamTemplateBatch::new(0, vec![template]);
    let _kept = run_batch(&mut state, &cfg, batch, None).expect("run_batch");

    let slot = &cfg.metrics.slots()[0];
    let metrics = slot.lock();
    assert_eq!(metrics.templates_processed, 1, "every template counts exactly once");
    assert_eq!(metrics.missing_umis, expected_missing, "missing_umis");
    assert_eq!(metrics.wrong_length, expected_wrong_length, "wrong_length");
    assert_eq!(metrics.mismatched, expected_mismatched, "mismatched");
    assert_eq!(bucket_total(&metrics, "AAAA"), expected_aaaa_bucket, "AAAA bucket");
    assert_eq!(bucket_total(&metrics, "NNNN"), expected_n_bucket, "all-N bucket");
}

/// The kept-only path keeps exactly the same templates as the with-rejects
/// path's kept branch, but silently DROPS rejected templates (there is no
/// rejects branch to inspect).
#[test]
fn run_batch_kept_only_drops_rejected() {
    let cfg = Arc::new(make_default_cfg());
    let mut state = make_state();

    let templates = vec![
        make_template_with_rx(b"exact", "AAAA"),
        make_template_with_rx(b"correctable", "AAAT"),
        make_template_with_rx(b"offlist", "TTTT"),
    ];
    let batch = BamTemplateBatch::new(11, templates);

    let kept = run_batch(&mut state, &cfg, batch, None).expect("run_batch kept-only");

    // Only the AAAA + correctable templates survive; TTTT is dropped silently.
    assert_eq!(kept.batch_serial(), 11);
    assert_eq!(kept.templates().len(), 2, "TTTT dropped, no rejects branch");
    let mut kept_umis: Vec<String> = kept.templates().iter().map(template_rx).collect();
    kept_umis.sort();
    assert_eq!(kept_umis, vec!["AAAA".to_string(), "AAAA".to_string()]);
    // Pin kept identity by QNAME: the on-whitelist `exact` and the `correctable`
    // template survive; `offlist` (TTTT) is the one dropped.
    let mut kept_qnames: Vec<String> = kept.templates().iter().map(template_qname).collect();
    kept_qnames.sort();
    assert_eq!(kept_qnames, vec!["correctable".to_string(), "exact".to_string()]);

    assert_eq!(cfg.records_emitted.load(Ordering::Relaxed), 2);

    let slot = &cfg.metrics.slots()[0];
    let metrics = slot.lock();
    assert_eq!(metrics.mismatched, 1, "rejected template still counted in metrics");
}

/// Decode a BAM-framed reject block into the `QNAME` of each record. The
/// missing-UMI case carries no `RX` to key on, so identity is pinned by read
/// name instead -- an independent inverse of `append_framed_raw_record`.
fn reject_block_qnames(bytes: &[u8]) -> Vec<String> {
    let mut names = Vec::new();
    let mut offset = 0;
    while offset < bytes.len() {
        let block_size = u32::from_le_bytes(
            bytes[offset..offset + 4].try_into().expect("4-byte block_size prefix"),
        ) as usize;
        offset += 4;
        let body = &bytes[offset..offset + block_size];
        let name = fgumi_raw_bam::read_name(body);
        names.push(String::from_utf8(name.to_vec()).expect("qname is valid UTF-8"));
        offset += block_size;
    }
    names
}

/// A template with no `RX` tag is a missing-UMI reject. On the 2-output
/// (`--rejects`) path its record must be framed into the rejects block
/// unchanged -- the `run_batch` branch a kept-only run never takes -- while it
/// credits only `missing_umis` and no per-UMI bucket.
#[test]
fn run_batch_with_rejects_routes_missing_umi_to_rejects() {
    let cfg = Arc::new(make_default_cfg());
    let mut state = make_state();

    let batch = BamTemplateBatch::new(
        5,
        vec![make_template_with_rx(b"kept", "AAAA"), make_template_without_rx(b"no_umi")],
    );
    let (kept, rejects) = run_batch_collecting_rejects(&mut state, &cfg, batch);

    // Only the on-whitelist template survives; the missing-UMI one is rejected.
    assert_eq!(kept.templates().len(), 1);
    assert_eq!(template_qname(&kept.templates()[0]), "kept");

    // The missing-UMI record reaches the rejects block, keyed by QNAME since it
    // has no RX. A non-empty-bytes check alone would pass on the wrong record.
    assert_eq!(rejects.batch_serial, 5);
    assert_eq!(
        reject_block_qnames(&rejects.bytes),
        vec!["no_umi".to_string()],
        "the missing-UMI record is framed into rejects unchanged",
    );
    assert_eq!(cfg.records_emitted.load(Ordering::Relaxed), 1, "only the kept record is emitted");

    let slot = &cfg.metrics.slots()[0];
    let metrics = slot.lock();
    assert_eq!(metrics.missing_umis, 1, "the RX-less template credits missing_umis");
    assert_eq!(metrics.mismatched, 0, "a missing UMI is not a mismatch");
}

// ============================================================================
// End-to-end pipeline tests
//
// The batch-routing tests above call `run_batch` directly, which leaves the
// factory functions' `init` (per-worker cache construction) and `body`
// (the closure the framework invokes per batch) closures unexercised — they
// only run when the framework drives the step. These tests assemble a real
// `Pipeline` around each factory so those closures execute, at 1 and 4 threads
// and with the cache both enabled and disabled (the `Some`/`None` arms of the
// lazy-init closure).
// ============================================================================

/// Build a `CorrectStepConfig` with an explicit `cache_size`. `cache_size = 0`
/// drives the factory's cache-disabled (`None`) init arm; any positive value
/// drives the `Some(LruCache)` arm.
fn make_cfg_with_cache(cache_size: usize) -> CorrectStepConfig {
    let mut cfg = make_default_cfg();
    cfg.opts.cache_size = cache_size;
    cfg
}

/// The standard mixed batch used by the end-to-end tests: one on-whitelist
/// template (kept, unchanged), one correctable (kept, rewritten to `AAAA`),
/// one off-whitelist (rejected: mismatched), and one missing-UMI (rejected:
/// missing). `serial` distinguishes the batches so a routing regression that
/// dropped or duplicated a batch is visible.
fn mixed_batch(serial: usize) -> BamTemplateBatch {
    let templates = vec![
        make_template_with_rx(b"exact", "AAAA"),
        make_template_with_rx(b"correctable", "AAAT"),
        make_template_with_rx(b"offlist", "TTTT"),
        make_template_without_rx(b"no_umi"),
    ];
    BamTemplateBatch::new(serial as u64, templates)
}

/// An all-clean batch: two templates that both survive correction (`exact`
/// stays `AAAA`, `correctable` is rewritten to `AAAA`). Nothing is rejected, so
/// the factory's rejects branch must still emit one *empty* block for this
/// serial to keep the `ByItemOrdinal` sequence dense (X5-001).
fn clean_batch(serial: usize) -> BamTemplateBatch {
    let templates = vec![
        make_template_with_rx(b"exact", "AAAA"),
        make_template_with_rx(b"correctable", "AAAT"),
    ];
    BamTemplateBatch::new(serial as u64, templates)
}

/// Test-local source that replays pre-built items into the chain. `Exclusive`
/// so it owns its cursor and is never cloned per worker. Held-item retry keeps
/// the ordinal sequence dense for the downstream `ByItemOrdinal` reorder.
struct ReplaySource<T: Send + HeapSize + Ordered + 'static> {
    items: VecDeque<T>,
    held: HeldSlot<Unpushed<T>>,
}

impl<T: Send + HeapSize + Ordered + 'static> ReplaySource<T> {
    fn new(items: Vec<T>) -> Self {
        Self { items: items.into(), held: HeldSlot::new() }
    }
}

impl<T: Send + HeapSize + Ordered + 'static> Step for ReplaySource<T> {
    type Input = ();
    type Outputs = OrderedBytesSingle<T>;

    fn profile(&self) -> StepProfile {
        StepProfile {
            name: "ReplaySource",
            kind: StepKind::Exclusive,
            sticky: false,
            output_queues: vec![QueueSpec::ByteBounded { limit_bytes: 4096 }],
            branch_ordering: vec![BranchOrdering::ByItemOrdinal],
        }
    }

    fn try_run(&mut self, ctx: &mut StepCtx<'_, Self>) -> io::Result<StepOutcome> {
        if let Some(unpushed) = self.held.take() {
            match ctx.outputs.retry(unpushed) {
                Ok(()) => {}
                Err(again) => {
                    self.held.put(again);
                    return Ok(StepOutcome::Contention);
                }
            }
        }
        let Some(item) = self.items.pop_front() else {
            return Ok(StepOutcome::Finished);
        };
        if let Err(unpushed) = ctx.outputs.push(item) {
            self.held.put(unpushed);
        }
        Ok(StepOutcome::Progress)
    }
}

/// Terminal sink that accumulates every item it receives. `Exclusive` so
/// arrival is the chain's output order with no sink-side interleaving.
struct CollectSink<T: Send + HeapSize + 'static> {
    collected: Arc<Mutex<Vec<T>>>,
}

impl<T: Send + HeapSize + 'static> Step for CollectSink<T> {
    type Input = T;
    type Outputs = ();

    fn profile(&self) -> StepProfile {
        StepProfile {
            name: "CollectSink",
            kind: StepKind::Exclusive,
            sticky: false,
            output_queues: vec![],
            branch_ordering: vec![],
        }
    }

    fn try_run(&mut self, ctx: &mut StepCtx<'_, Self>) -> io::Result<StepOutcome> {
        match ctx.input.pop() {
            Some(item) => {
                self.collected.lock().expect("sink mutex not poisoned").push(item);
                Ok(StepOutcome::Progress)
            }
            None if ctx.input.is_drained() => Ok(StepOutcome::Finished),
            None => Ok(StepOutcome::NoProgress),
        }
    }
}

/// Sum a metric across every accumulator slot. Workers claim slots by thread
/// index, so at >1 thread the counts can land in different slots.
fn sum_metrics<F: Fn(&CollectedCorrectMetrics) -> u64>(
    metrics: &PerThreadAccumulator<CollectedCorrectMetrics>,
    field: F,
) -> u64 {
    metrics.slots().iter().map(|slot| field(&slot.lock())).sum()
}

/// Sort the kept QNAMEs collected from every kept batch.
fn sorted_kept_qnames(batches: &Arc<Mutex<Vec<BamTemplateBatch>>>) -> Vec<String> {
    let guard = batches.lock().expect("kept mutex not poisoned");
    let mut names: Vec<String> =
        guard.iter().flat_map(|b| b.templates().iter().map(template_qname)).collect();
    names.sort();
    names
}

/// The batch serials of every kept batch, in the order they arrived at the
/// sink. Both output branches are `ByItemOrdinal`-ordered, so an intact run
/// yields the dense sequence `0..n_batches`; a dropped-then-duplicated batch
/// (which identical records would hide from a QNAME check) shows up here as a
/// repeated or missing serial.
fn kept_batch_serials(batches: &Arc<Mutex<Vec<BamTemplateBatch>>>) -> Vec<u64> {
    batches
        .lock()
        .expect("kept mutex not poisoned")
        .iter()
        .map(BamTemplateBatch::batch_serial)
        .collect()
}

/// The batch serials of every rejects block, in sink-arrival order. The factory
/// emits exactly one block per input batch (empty when nothing was rejected),
/// so these must also form the dense sequence `0..n_batches`.
fn reject_block_serials(blocks: &Arc<Mutex<Vec<DecompressedBlock>>>) -> Vec<u64> {
    blocks.lock().expect("rejects mutex not poisoned").iter().map(|b| b.batch_serial).collect()
}

/// The 2-output factory, driven through a real `Pipeline`, must split every
/// batch into a kept branch (on-whitelist + correctable templates) and a
/// rejects branch (framed off-whitelist + missing-UMI records), preserving the
/// per-batch routing across `n_batches` batches and both thread counts. Running
/// it end-to-end is what exercises the factory's `init` and `body` closures.
#[rstest]
fn correct_step_with_rejects_pipeline_splits_every_batch(
    #[values(1, 4)] threads: usize,
    #[values(0, 100_000)] cache_size: usize,
) {
    let n_batches: usize = 3;
    let cfg = make_cfg_with_cache(cache_size);
    let metrics = Arc::clone(&cfg.metrics);
    let records_emitted = Arc::clone(&cfg.records_emitted);

    let inputs: Vec<BamTemplateBatch> = (0..n_batches).map(mixed_batch).collect();
    let kept: Arc<Mutex<Vec<BamTemplateBatch>>> = Arc::new(Mutex::new(Vec::new()));
    let rejects: Arc<Mutex<Vec<DecompressedBlock>>> = Arc::new(Mutex::new(Vec::new()));
    let (kept_sink, rejects_sink) = (Arc::clone(&kept), Arc::clone(&rejects));

    let builder = Pipeline::builder();
    let multi =
        builder.chain(ReplaySource::new(inputs)).chain(correct_step_with_rejects(cfg)).into_multi();
    multi.b0.chain(CollectSink { collected: kept_sink }).into_sink_marker();
    multi.b1.chain(CollectSink { collected: rejects_sink }).into_sink_marker();
    let pipeline = builder.build().expect("chain builds");
    pipeline.run(PipelineConfig { threads, ..Default::default() }).expect("chain runs");

    // Two templates kept per batch: the on-whitelist `exact` and `correctable`
    // (sorted: all `correctable` then all `exact`).
    let expected_kept: Vec<String> = std::iter::repeat_n("correctable".to_string(), n_batches)
        .chain(std::iter::repeat_n("exact".to_string(), n_batches))
        .collect();
    assert_eq!(sorted_kept_qnames(&kept), expected_kept, "kept branch keeps exact + correctable");

    // Rejects: the off-whitelist `offlist` and the missing-UMI `no_umi`, framed
    // once per batch. Concatenating the framed blocks yields a self-delimiting
    // sequence the decoder walks record by record.
    let reject_bytes: Vec<u8> = rejects
        .lock()
        .expect("rejects mutex not poisoned")
        .iter()
        .flat_map(|block| block.bytes.clone())
        .collect();
    let mut reject_names = reject_block_qnames(&reject_bytes);
    reject_names.sort();
    // Sorted: all `no_umi` then all `offlist`.
    let expected_rejects: Vec<String> = std::iter::repeat_n("no_umi".to_string(), n_batches)
        .chain(std::iter::repeat_n("offlist".to_string(), n_batches))
        .collect();
    assert_eq!(reject_names, expected_rejects, "rejects branch frames offlist + missing per batch");

    // Every input batch has identical records, so QNAME counts alone cannot
    // distinguish a dropped batch that a duplicate silently replaced. Pin the
    // per-batch identity by asserting both branches carry serials `0..n_batches`
    // in output order.
    let expected_serials: Vec<u64> = (0..n_batches as u64).collect();
    assert_eq!(
        kept_batch_serials(&kept),
        expected_serials,
        "kept branch carries a dense serial run"
    );
    assert_eq!(
        reject_block_serials(&rejects),
        expected_serials,
        "rejects branch emits one block per input serial, densely",
    );

    let n_batches = n_batches as u64;
    assert_eq!(
        records_emitted.load(Ordering::Relaxed),
        2 * n_batches,
        "records_emitted counts the two kept single-record templates per batch",
    );
    assert_eq!(sum_metrics(&metrics, |m| m.mismatched), n_batches, "one mismatched per batch");
    assert_eq!(sum_metrics(&metrics, |m| m.missing_umis), n_batches, "one missing-UMI per batch");
    assert_eq!(
        sum_metrics(&metrics, |m| m.templates_processed),
        4 * n_batches,
        "every template in every batch is counted once",
    );
}

/// A fully clean input (every template kept) driven end-to-end through
/// `correct_step_with_rejects`. The `run_batch_with_rejects_no_rejects_when_all_kept`
/// unit test asserts the same X5-001 invariant against `run_batch` directly, but
/// only the real pipeline exercises the factory body and the rejects branch's
/// `ByItemOrdinal` reorder stage: a regression that omitted the second output
/// (pushing `only_a`) would wedge that stage rather than fail an assert. This
/// case proves the pipeline completes and emits exactly one *empty* rejects
/// block per input serial.
#[rstest]
fn correct_step_with_rejects_pipeline_emits_empty_rejects_for_clean_input(
    #[values(1, 4)] threads: usize,
) {
    let n_batches: usize = 3;
    let cfg = make_default_cfg();
    let records_emitted = Arc::clone(&cfg.records_emitted);

    let inputs: Vec<BamTemplateBatch> = (0..n_batches).map(clean_batch).collect();
    let kept: Arc<Mutex<Vec<BamTemplateBatch>>> = Arc::new(Mutex::new(Vec::new()));
    let rejects: Arc<Mutex<Vec<DecompressedBlock>>> = Arc::new(Mutex::new(Vec::new()));
    let (kept_sink, rejects_sink) = (Arc::clone(&kept), Arc::clone(&rejects));

    let builder = Pipeline::builder();
    let multi =
        builder.chain(ReplaySource::new(inputs)).chain(correct_step_with_rejects(cfg)).into_multi();
    multi.b0.chain(CollectSink { collected: kept_sink }).into_sink_marker();
    multi.b1.chain(CollectSink { collected: rejects_sink }).into_sink_marker();
    let pipeline = builder.build().expect("chain builds");
    pipeline.run(PipelineConfig { threads, ..Default::default() }).expect("chain runs");

    // Both templates in every batch survive, so all `2 * n_batches` records are emitted.
    let expected_serials: Vec<u64> = (0..n_batches as u64).collect();
    assert_eq!(
        kept_batch_serials(&kept),
        expected_serials,
        "kept branch carries a dense serial run"
    );
    let kept_count: usize =
        kept.lock().expect("kept mutex not poisoned").iter().map(|b| b.templates().len()).sum();
    assert_eq!(kept_count, 2 * n_batches, "every clean template is kept");
    assert_eq!(records_emitted.load(Ordering::Relaxed), 2 * n_batches as u64);

    // The rejects branch emits exactly one block per input serial, each empty —
    // the X5-001 invariant that keeps `ByItemOrdinal` dense on a clean input.
    assert_eq!(
        reject_block_serials(&rejects),
        expected_serials,
        "one rejects block per input serial, densely",
    );
    assert!(
        rejects.lock().expect("rejects mutex not poisoned").iter().all(|b| b.bytes.is_empty()),
        "a clean batch emits a zero-byte rejects block, not a gap",
    );
}

/// The kept-only factory, driven through a real `Pipeline`, keeps the same
/// templates as the 2-output path but drops the rejected ones silently — there
/// is no rejects branch. Running it end-to-end exercises the kept-only
/// factory's `init` and `body` closures.
#[rstest]
fn correct_step_kept_only_pipeline_keeps_survivors(
    #[values(1, 4)] threads: usize,
    #[values(0, 100_000)] cache_size: usize,
) {
    let n_batches: usize = 3;
    let cfg = make_cfg_with_cache(cache_size);
    let metrics = Arc::clone(&cfg.metrics);
    let records_emitted = Arc::clone(&cfg.records_emitted);

    let inputs: Vec<BamTemplateBatch> = (0..n_batches).map(mixed_batch).collect();
    let kept: Arc<Mutex<Vec<BamTemplateBatch>>> = Arc::new(Mutex::new(Vec::new()));
    let kept_sink = Arc::clone(&kept);

    let builder = Pipeline::builder();
    builder
        .chain(ReplaySource::new(inputs))
        .chain(correct_step_kept_only(cfg))
        .chain(CollectSink { collected: kept_sink })
        .into_sink_marker();
    let pipeline = builder.build().expect("chain builds");
    pipeline.run(PipelineConfig { threads, ..Default::default() }).expect("chain runs");

    let expected_kept: Vec<String> = std::iter::repeat_n("correctable".to_string(), n_batches)
        .chain(std::iter::repeat_n("exact".to_string(), n_batches))
        .collect();
    assert_eq!(sorted_kept_qnames(&kept), expected_kept, "kept-only keeps exact + correctable");

    let n_batches = n_batches as u64;
    assert_eq!(
        records_emitted.load(Ordering::Relaxed),
        2 * n_batches,
        "records_emitted counts the two kept templates per batch",
    );
    // Rejected templates are dropped but still counted in metrics.
    assert_eq!(sum_metrics(&metrics, |m| m.mismatched), n_batches, "one mismatched per batch");
    assert_eq!(sum_metrics(&metrics, |m| m.missing_umis), n_batches, "one missing-UMI per batch");
}

/// `CorrectWorkerState` deliberately reports `heap_size = 0`: the `LruCache` it
/// wraps has no `HeapSize` impl, so its footprint is intentionally excluded from
/// the pipeline's memory bound (documented on the `impl`). Pin that for both the
/// cache-enabled and cache-disabled workers so a future switch to real
/// accounting has to update this consciously.
#[test]
fn correct_worker_state_reports_zero_heap_size() {
    assert_eq!(make_state().heap_size(), 0, "cache-enabled worker reports zero");
    assert_eq!(
        CorrectWorkerState { cache: None }.heap_size(),
        0,
        "cache-disabled worker reports zero",
    );
}
