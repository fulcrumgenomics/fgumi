use super::*;
use crate::commands::correct::CorrectOptions;
use crate::pipeline::core::queues::QueueSpec;
use crate::pipeline::core::step::StepKind;
use crate::sam::SamTag;
use crate::template::Template;
use fgumi_raw_bam::SamBuilder;

fn make_default_cfg() -> CorrectStepConfig {
    CorrectStepConfig {
        encoded_umi_set: Arc::new(EncodedUmiSet::new(&["AAAA".to_string()])),
        umi_length: 4,
        umi_tag: SamTag::RX.into(),
        opts: CorrectOptions::default(),
        metrics: PerThreadAccumulator::new(1),
        records_emitted: Arc::new(AtomicU64::new(0)),
        output_byte_limit: 1024,
        unmatched_umi: "NNNN".to_string(),
    }
}

#[test]
fn correct_step_with_rejects_profile() {
    let cfg = make_default_cfg();
    let step = correct_step_with_rejects(cfg);
    let profile = step.profile();
    assert_eq!(profile.name, "correct");
    assert!(matches!(profile.kind, StepKind::Parallel));
    assert_eq!(profile.output_queues.len(), 2);
    for q in &profile.output_queues {
        assert!(matches!(q, QueueSpec::ByteBounded { .. }));
    }
}

#[test]
fn correct_step_kept_only_profile() {
    let cfg = make_default_cfg();
    let step = correct_step_kept_only(cfg);
    let profile = step.profile();
    assert_eq!(profile.name, "correct");
    assert!(matches!(profile.kind, StepKind::Parallel));
    assert_eq!(profile.output_queues.len(), 1);
}

// --- Batch-routing tests for run_batch_with_rejects / run_batch_kept_only ---

/// Build a fresh per-worker state matching how the step's `init` closure
/// constructs it for the default (cache-enabled) `CorrectOptions`.
fn make_state() -> CorrectWorkerState {
    CorrectWorkerState {
        cache: Some(LruCache::new(
            NonZero::new(CorrectOptions::default().cache_size).expect("default cache_size > 0"),
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

    let (kept, rejects) =
        run_batch_with_rejects(&mut state, &cfg, batch).expect("run_batch_with_rejects");

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
    let rejects = rejects.expect("rejected template routed to rejects branch");
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

/// An all-on-whitelist batch produces no rejects (`None`) from the
/// 2-output path, so the framework emits a zero-byte rejects block instead.
#[test]
fn run_batch_with_rejects_no_rejects_when_all_kept() {
    let cfg = Arc::new(make_default_cfg());
    let mut state = make_state();

    let batch = BamTemplateBatch::new(3, vec![make_template_with_rx(b"exact", "AAAA")]);
    let (kept, rejects) =
        run_batch_with_rejects(&mut state, &cfg, batch).expect("run_batch_with_rejects");

    assert_eq!(kept.templates().len(), 1);
    assert!(rejects.is_none(), "no rejected records -> no rejects block");
    assert_eq!(cfg.records_emitted.load(Ordering::Relaxed), 1);
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

    let (kept, ()) = run_batch_kept_only(&mut state, &cfg, batch).expect("run_batch_kept_only");

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
