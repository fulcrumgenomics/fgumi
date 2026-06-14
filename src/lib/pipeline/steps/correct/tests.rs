use super::*;
use crate::commands::correct::CorrectOptions;
use crate::pipeline::core::queues::QueueSpec;
use crate::pipeline::core::step::StepKind;
use crate::sam::SamTag;

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
