//! Three independent validators for [`super::ChainSpec`]:
//!
//! 1. [`validate_stage_progression`] — ordering rules and mutual exclusions.
//! 2. [`validate_stage_opts_present`] — each stage has its options in the bag.
//! 3. [`validate_cross_stage_constraints`] — placeholder for Phase 3 cross-stage rules.

use anyhow::{Result, bail};

use crate::pipeline::chains::{ChainSpec, Stage};

/// Linear position of a stage in the canonical pipeline order. Used
/// by [`validate_stage_progression`] to enforce that stages appear
/// in a buildable order.
///
/// Mirrors `runall.rs::RunAllStage::ord` and extends it for the
/// standalone-only stages. Consensus stages share ord = 7 because
/// they're mutually exclusive in a chain.
///
/// Adjustments vs runall's ord:
/// - Extract = 0 (input-producing, must precede everything)
/// - Correct = 1 (was 0)
/// - Align = 2 (was 1, was `AlignAndMerge`)
/// - Zipper = 3 (was 2)
/// - Sort = 4 (was 3)
/// - Group = 5 (was 4)
/// - Clip, Dedup, Downsample = 6 (post-group standalone-only)
/// - Simplex, Duplex, Codec = 7 (terminal consensus, was 5)
/// - Filter = 8 (post-consensus filtering; follows the consensus caller
///   in a runall `consensus → filter` chain, so it sorts after the
///   consensus bucket)
fn stage_ord(stage: Stage) -> usize {
    match stage {
        Stage::Extract => 0,
        Stage::Correct => 1,
        Stage::Align => 2,
        Stage::Zipper => 3,
        Stage::Sort => 4,
        Stage::Group => 5,
        Stage::Clip | Stage::Dedup | Stage::Downsample => 6,
        Stage::Simplex | Stage::Duplex | Stage::Codec => 7,
        Stage::Filter => 8,
    }
}

/// Reject illegal stage orderings.
///
/// Rules:
/// - At least one stage (empty chain is meaningless).
/// - Each adjacent pair must be non-decreasing in `stage_ord`.
/// - At most one `Align` per chain (mutually exclusive with `Zipper`).
/// - At most one `Zipper` per chain (mutually exclusive with `Align`).
/// - At most one consensus stage (`Simplex`/`Duplex`/`Codec`).
/// - If a consensus stage is present, it must be either the last stage
///   or immediately followed by a single trailing `Filter` (the
///   `consensus → filter` runall chain) — nothing else may follow
///   consensus.
/// - `Extract` must be first when present.
///
/// # Errors
///
/// Returns the first violation with a human-readable message.
pub fn validate_stage_progression(spec: &ChainSpec) -> Result<()> {
    let stages = &spec.stages;
    if stages.is_empty() {
        bail!("ChainSpec.stages is empty — a chain needs at least one stage");
    }

    // Extract must be first when present.
    if let Some(extract_pos) = stages.iter().position(|s| *s == Stage::Extract) {
        if extract_pos != 0 {
            bail!(
                "Stage::Extract must be the first stage in the chain (got position {extract_pos})"
            );
        }
    }

    // Non-decreasing ord (allowing equality for stages that share an ord bucket,
    // e.g. multiple post-group standalone-only stages).
    for window in stages.windows(2) {
        let a = window[0];
        let b = window[1];
        if stage_ord(a) > stage_ord(b) {
            bail!(
                "Stage {a:?} (ord {}) cannot precede stage {b:?} (ord {}); \
                 stages must appear in canonical pipeline order",
                stage_ord(a),
                stage_ord(b),
            );
        }
    }

    // Mutually-exclusive: at most one Align, at most one Zipper.
    let align_count = stages.iter().filter(|s| **s == Stage::Align).count();
    let zipper_count = stages.iter().filter(|s| **s == Stage::Zipper).count();
    if align_count > 1 {
        bail!("At most one Stage::Align per chain; got {align_count}");
    }
    if zipper_count > 1 {
        bail!("At most one Stage::Zipper per chain; got {zipper_count}");
    }
    if align_count > 0 && zipper_count > 0 {
        bail!(
            "Stage::Align and Stage::Zipper are mutually exclusive in a chain; \
             got both"
        );
    }

    // At most one consensus stage, and it must be terminal.
    let consensus_positions: Vec<usize> =
        stages.iter().enumerate().filter(|(_, s)| s.is_consensus()).map(|(i, _)| i).collect();
    if consensus_positions.len() > 1 {
        bail!(
            "At most one consensus stage per chain; got {} ({:?})",
            consensus_positions.len(),
            consensus_positions.iter().map(|i| stages[*i]).collect::<Vec<_>>()
        );
    }
    if let Some(&pos) = consensus_positions.first() {
        let last = stages.len() - 1;
        // Consensus may be last, or second-to-last followed by a single
        // trailing Filter (the runall `consensus → filter` chain). Any
        // other stage after consensus is illegal.
        let consensus_terminal = pos == last;
        let consensus_then_filter = pos + 1 == last && stages[last] == Stage::Filter;
        if !consensus_terminal && !consensus_then_filter {
            bail!(
                "Consensus stage {:?} at position {pos} must be terminal \
                 (last in the chain) or immediately followed by a single \
                 trailing Filter; chain has {} stages total",
                stages[pos],
                stages.len()
            );
        }
    }

    Ok(())
}

/// Reject specs where a referenced stage has no options in the bag.
///
/// As of Phase 5 T5.3, eleven stages have their options in the bag:
/// Correct, Sort, Group, Zipper, Duplex, Codec, Dedup, Filter, Clip, Simplex,
/// and Extract. The remaining stage (Downsample) adds its slot incrementally
/// during its migration task (T2.19–T2.22). For now that stage skips the
/// options-presence check — once its slot lands in the bag, add the check here.
///
/// # Errors
///
/// Returns the first stage with a referenced-but-missing entry.
pub fn validate_stage_opts_present(spec: &ChainSpec) -> Result<()> {
    let bag = &spec.stage_opts;
    for stage in &spec.stages {
        // Downsample has no bag slot yet — fail early with a clear message.
        if *stage == Stage::Downsample {
            bail!(
                "Stage::Downsample has no bag slot yet — \
                 downsample isn't migrated to chains::build_for"
            );
        }

        // Check wired stages have their options present.
        let present = match stage {
            Stage::Correct => bag.correct.is_some(),
            Stage::Zipper => bag.zipper.is_some(),
            Stage::Sort => bag.sort.is_some(),
            Stage::Group => bag.group.is_some(),
            #[cfg(feature = "consensus")]
            Stage::Duplex => bag.duplex.is_some(),
            #[cfg(feature = "consensus")]
            Stage::Codec => bag.codec.is_some(),
            Stage::Dedup => bag.dedup.is_some(),
            Stage::Filter => bag.filter.is_some(),
            Stage::Clip => bag.clip.is_some(),
            #[cfg(feature = "consensus")]
            Stage::Simplex => bag.simplex.is_some(),
            // Without the `consensus` feature the consensus option-bag slots do
            // not exist, so this stage can never be satisfied. Bail with an
            // explicit feature-disabled error rather than returning `false`,
            // which would surface the generic "options not provided" message
            // below and tell reduced-feature callers to populate bag slots that
            // do not exist. Unreachable in practice — the consensus CLI commands
            // are compiled out too — but the message is correct if reached.
            #[cfg(not(feature = "consensus"))]
            Stage::Duplex | Stage::Codec | Stage::Simplex => {
                bail!("Stage {stage:?} requires building fgumi with the `consensus` feature");
            }
            Stage::Align => bag.aligner.is_some(),
            Stage::Extract => bag.extract.is_some(),
            Stage::Downsample => unreachable!(),
        };
        if !present {
            bail!(
                "Stage {stage:?} requested but its options are not provided \
                 in StageOptionsBag"
            );
        }
    }
    Ok(())
}

/// Reject specs that violate cross-stage constraints (e.g. duplex
/// requires group's strategy to be Paired).
///
/// Currently enforces two rules:
///
/// 1. `Stage::Zipper` requires `SourceSpec::PairedBams`.
/// 2. `Stage::Duplex` following `Stage::Group` in the same chain
///    requires `GroupOptions::strategy` to be `Strategy::Paired`.
///    Other strategies don't tag `/A` vs `/B` endpoints in the way
///    `DuplexConsensusCaller` expects.
///
/// **Rule 2 scope.** The rule fires only when both `Stage::Group` and
/// `Stage::Duplex` appear in the same chain. Standalone `fgumi duplex`
/// (chain = `[Stage::Duplex]`) skips the check because the input BAM is
/// already MI-tagged with `/A`/`/B` annotations from a prior group run
/// and the strategy is not re-checkable from the spec alone.
///
/// **Deferred rule.** Dedup's "input must be template-coordinate sorted"
/// rule is checked at chain-build time (`build_dedup_chain` reads the
/// input BAM header) — it cannot be determined from `spec.stages` alone,
/// so it intentionally stays out of this validator.
///
/// # Errors
///
/// Returns an error on the first violated constraint.
pub fn validate_cross_stage_constraints(spec: &ChainSpec) -> Result<()> {
    use crate::assigner::Strategy;
    use crate::pipeline::chains::{SinkSpec, SourceSpec};

    // Rule 1: Zipper requires a PairedBams source.
    if spec.stages.contains(&Stage::Zipper) {
        match &spec.source {
            SourceSpec::PairedBams { .. } => {}
            _ => {
                bail!(
                    "Stage::Zipper requires SourceSpec::PairedBams (unmapped + mapped + reference); \
                     got a different source variant"
                );
            }
        }
    }

    // Rule 2: Duplex (terminal) following Group (assigning MI ids) requires
    // Group's strategy to be Paired. Other strategies don't tag /A vs /B
    // endpoints in the way DuplexConsensusCaller expects.
    //
    // Standalone duplex (chain = [Stage::Duplex]) skips Group; the input is
    // already MI-tagged with /A/B annotations from a prior group run. We can't
    // validate that at spec-construct time, so the rule only fires when both
    // stages are present in the same chain.
    //
    // Dedup's "input must be template-coordinate sorted" rule is checked at
    // chain-build time (build_dedup_chain reads the input BAM header) — it
    // can't be checked from spec.stages alone, so it stays out of this
    // validator by design.
    if spec.stages.contains(&Stage::Duplex) && spec.stages.contains(&Stage::Group) {
        if let Some(group_opts) = spec.stage_opts.group.as_ref() {
            if !matches!(group_opts.strategy, Strategy::Paired) {
                bail!(
                    "Stage::Duplex requires Stage::Group to use Strategy::Paired \
                     (got Strategy::{:?})",
                    group_opts.strategy
                );
            }
        }
    }

    // Rule 3: `SinkSpec::BamWithIndex` requires the terminal stage to be
    // `Stage::Sort`. The BAI file format indexes BGZF virtual offsets of a
    // coordinate-sorted BAM; emitting one for a non-sort terminal stage (or
    // for a chain where sort is intermediate, leaving some other stage as
    // terminal) would either reference a file that isn't coordinate-sorted
    // or attempt to index an output type the hook isn't equipped to handle.
    // The per-command CLI check (`Sort::execute` rejects
    // `--write-index --order queryname`) handles the "coordinate
    // specifically, not template-coordinate or queryname" constraint.
    if matches!(spec.sink, SinkSpec::BamWithIndex(_)) {
        let terminal_is_sort = spec.stages.last().is_some_and(|s| *s == Stage::Sort);
        if !terminal_is_sort {
            // `validate_stage_progression` rejects empty `stages` before
            // this validator runs, so `last()` is guaranteed `Some` here.
            // The `unreachable!` documents the invariant inline.
            let terminal = spec.stages.last().ok_or_else(|| {
                anyhow::anyhow!(
                    "SinkSpec::BamWithIndex requires Stage::Sort as the terminal chain stage; \
                     got an empty chain (should have been caught by validate_stage_progression)",
                )
            })?;
            bail!(
                "SinkSpec::BamWithIndex requires Stage::Sort as the terminal chain stage; \
                 got terminal stage {terminal:?}",
            );
        }
    }

    // Rule 4: `Stage::Extract` requires a `SourceSpec::Fastqs` source.
    // Extract reads directly from FASTQ files; feeding it a BAM, SAM, or
    // paired-BAM source is a programmer error that would be caught late
    // (at chain-build time) without this guard. Catching it here at
    // spec-validation time gives a clearer error message.
    if spec.stages.contains(&Stage::Extract) && !matches!(spec.source, SourceSpec::Fastqs { .. }) {
        bail!("Stage::Extract requires SourceSpec::Fastqs; got {:?}", spec.source);
    }

    Ok(())
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::pipeline::chains::{SinkSpec, SourceSpec};
    use std::path::PathBuf;

    fn empty_spec(stages: Vec<Stage>) -> ChainSpec {
        use crate::commands::common::{
            CompressionOptions, QueueMemoryOptions, SchedulerOptions, ThreadingOptions,
        };
        use crate::pipeline::chains::StageOptionsBag;
        ChainSpec {
            stages,
            source: SourceSpec::Bam(PathBuf::from("in.bam")),
            sink: SinkSpec::Bam(PathBuf::from("out.bam")),
            stage_opts: StageOptionsBag::default(),
            threading: ThreadingOptions { threads: None },
            compression: CompressionOptions::default(),
            scheduler: SchedulerOptions::default(),
            queue_memory: QueueMemoryOptions::default(),
            async_reader: false,
            command_line: String::new(),
        }
    }

    #[test]
    fn empty_chain_rejected() {
        let err = validate_stage_progression(&empty_spec(vec![])).unwrap_err();
        assert!(err.to_string().contains("empty"));
    }

    #[test]
    fn good_chain_accepted() {
        let spec = empty_spec(vec![Stage::Sort, Stage::Group, Stage::Simplex]);
        validate_stage_progression(&spec).expect("good chain should validate");
    }

    #[test]
    fn out_of_order_rejected() {
        let spec = empty_spec(vec![Stage::Group, Stage::Sort]);
        let err = validate_stage_progression(&spec).unwrap_err();
        assert!(err.to_string().contains("canonical pipeline order"));
    }

    #[test]
    fn align_and_zipper_mutually_exclusive() {
        let spec = empty_spec(vec![Stage::Align, Stage::Zipper]);
        let err = validate_stage_progression(&spec).unwrap_err();
        assert!(err.to_string().contains("mutually exclusive"));
    }

    #[test]
    fn consensus_must_be_terminal() {
        // Consensus followed by a non-Filter stage (Clip) is illegal:
        // Clip (ord 6) cannot follow Simplex (ord 7), and Clip is not the
        // permitted trailing Filter. Either the ordering check or the
        // terminal check fires first — both are valid rejections.
        let spec = empty_spec(vec![Stage::Simplex, Stage::Clip]);
        let err = validate_stage_progression(&spec).unwrap_err();
        let msg = err.to_string();
        assert!(msg.contains("terminal") || msg.contains("canonical pipeline order"), "got: {msg}");
    }

    #[test]
    fn consensus_then_filter_is_valid() {
        // The runall `consensus → filter` chain: Simplex (ord 7) followed
        // by a single trailing Filter (ord 8) is permitted.
        let spec = empty_spec(vec![Stage::Group, Stage::Simplex, Stage::Filter]);
        assert!(
            validate_stage_progression(&spec).is_ok(),
            "group → simplex → filter must be a valid progression"
        );
    }

    #[test]
    fn multiple_consensus_rejected() {
        let spec = empty_spec(vec![Stage::Simplex, Stage::Duplex]);
        let err = validate_stage_progression(&spec).unwrap_err();
        // Caught either by "mutually exclusive" or "at most one consensus" or "terminal"
        // depending on which check fires first; either is acceptable.
        let msg = err.to_string();
        assert!(msg.contains("consensus") || msg.contains("terminal"), "got: {msg}");
    }

    #[test]
    fn extract_must_be_first() {
        let spec = empty_spec(vec![Stage::Correct, Stage::Extract]);
        let err = validate_stage_progression(&spec).unwrap_err();
        assert!(err.to_string().contains("first"));
    }

    #[test]
    fn missing_options_caught() {
        let spec = empty_spec(vec![Stage::Correct]);
        let err = validate_stage_opts_present(&spec).unwrap_err();
        assert!(err.to_string().contains("Correct"));
    }

    #[test]
    fn cross_stage_correct_accepted() {
        // Non-Zipper stages are unaffected by the PairedBams constraint.
        let spec = empty_spec(vec![Stage::Correct]);
        validate_cross_stage_constraints(&spec).expect("Correct with Bam source is valid");
    }

    #[test]
    fn cross_stage_zipper_requires_paired_bams() {
        use crate::pipeline::chains::SourceSpec;

        // Zipper with a Bam source should be rejected.
        let mut spec = empty_spec(vec![Stage::Zipper]);
        spec.source = SourceSpec::Bam(std::path::PathBuf::from("in.bam"));
        let err = validate_cross_stage_constraints(&spec).unwrap_err();
        let msg = err.to_string();
        assert!(msg.contains("PairedBams"), "expected PairedBams message, got: {msg}");
    }

    #[test]
    fn cross_stage_zipper_with_paired_bams_accepted() {
        use crate::pipeline::chains::SourceSpec;

        let mut spec = empty_spec(vec![Stage::Zipper]);
        spec.source = SourceSpec::PairedBams {
            unmapped: std::path::PathBuf::from("u.bam"),
            mapped: std::path::PathBuf::from("m.bam"),
            reference: std::path::PathBuf::from("ref.fa"),
        };
        validate_cross_stage_constraints(&spec).expect("Zipper with PairedBams is valid");
    }

    #[test]
    fn cross_stage_bam_with_index_accepted_on_terminal_sort() {
        let mut spec = empty_spec(vec![Stage::Sort]);
        spec.sink = SinkSpec::BamWithIndex(std::path::PathBuf::from("out.bam"));
        validate_cross_stage_constraints(&spec)
            .expect("BamWithIndex on a Stage::Sort-terminal chain is valid");
    }

    #[test]
    fn cross_stage_bam_with_index_rejected_on_non_sort_terminal() {
        // Sort → Group (Sort is intermediate, Group is terminal) — BAI makes
        // no sense because Group's output doesn't go through the indexer's
        // BGZF-virtual-offset assumptions.
        let mut spec = empty_spec(vec![Stage::Sort, Stage::Group]);
        spec.sink = SinkSpec::BamWithIndex(std::path::PathBuf::from("out.bam"));
        let err = validate_cross_stage_constraints(&spec).unwrap_err();
        let msg = err.to_string();
        assert!(msg.contains("BamWithIndex"), "expected BamWithIndex in error, got: {msg}");
        assert!(msg.contains("Stage::Sort"), "expected Stage::Sort in error, got: {msg}");
    }

    #[test]
    fn cross_stage_bam_with_index_rejected_on_no_sort_chain() {
        // Group only — no sort at all.
        let mut spec = empty_spec(vec![Stage::Group]);
        spec.sink = SinkSpec::BamWithIndex(std::path::PathBuf::from("out.bam"));
        let err = validate_cross_stage_constraints(&spec).unwrap_err();
        assert!(err.to_string().contains("BamWithIndex"));
    }

    #[test]
    fn cross_stage_bam_with_index_accepted_on_multi_stage_terminal_sort() {
        // Multi-stage chain with Sort as the terminal stage — e.g. a
        // hypothetical `[Stage::Correct, Stage::Sort]` runall invocation
        // that asks for BAI. Validator accepts; `add_sort`'s streaming
        // branch (where this chain lands because `current_tail.is_some()`
        // after the upstream Correct step) registers the hook in mirror of
        // the standalone branch.
        let mut spec = empty_spec(vec![Stage::Correct, Stage::Sort]);
        spec.sink = SinkSpec::BamWithIndex(std::path::PathBuf::from("out.bam"));
        validate_cross_stage_constraints(&spec)
            .expect("BamWithIndex on a multi-stage chain with Sort terminal is valid");
    }

    #[test]
    fn cross_stage_bam_with_index_error_message_unwraps_terminal_stage() {
        // Regression guard against formatting `spec.stages.last()` as
        // `Option<&Stage>` (which would emit "Some(Group)") instead of
        // the unwrapped stage name. The user-facing error should read
        // "got terminal stage Group", not "got terminal stage Some(Group)".
        let mut spec = empty_spec(vec![Stage::Group]);
        spec.sink = SinkSpec::BamWithIndex(std::path::PathBuf::from("out.bam"));
        let msg = validate_cross_stage_constraints(&spec).unwrap_err().to_string();
        assert!(!msg.contains("Some("), "error message wrapped terminal stage in Some(): {msg}");
        assert!(msg.contains("Group"), "expected stage name in error, got: {msg}");
    }

    #[test]
    fn opts_present_downsample_not_wired() {
        // Stage::Downsample has no bag slot yet; should fail at options-presence time.
        let spec = empty_spec(vec![Stage::Downsample]);
        let err = validate_stage_opts_present(&spec).unwrap_err();
        let msg = err.to_string();
        assert!(msg.contains("Downsample"), "expected Downsample in error, got: {msg}");
        assert!(msg.contains("no bag slot"), "expected 'no bag slot' in error, got: {msg}");
    }

    #[test]
    fn opts_present_extract_missing_fails() {
        // Stage::Extract is now wired (T5.3); missing options should fail with
        // "options not provided" — not "no bag slot".
        let spec = empty_spec(vec![Stage::Extract]);
        let err = validate_stage_opts_present(&spec).unwrap_err();
        let msg = err.to_string();
        assert!(msg.contains("Extract"), "expected Extract in error, got: {msg}");
        assert!(
            !msg.contains("no bag slot"),
            "expected options-missing error (not 'no bag slot') for Extract, got: {msg}"
        );
    }

    #[test]
    fn opts_present_align_missing_fails() {
        // Stage::Align is now wired; missing options should fail with a
        // "options not provided" error (not "no bag slot").
        let spec = empty_spec(vec![Stage::Align]);
        let err = validate_stage_opts_present(&spec).unwrap_err();
        let msg = err.to_string();
        assert!(msg.contains("Align"), "expected Align in error, got: {msg}");
        assert!(
            !msg.contains("no bag slot"),
            "expected options-missing error (not 'no bag slot') for Align, got: {msg}"
        );
    }

    // These three exercise the consensus option-bag slots, which only exist
    // with the `consensus` feature.
    #[cfg(feature = "consensus")]
    #[test]
    fn duplex_after_group_requires_paired_strategy() {
        use crate::assigner::Strategy;
        use crate::commands::duplex::DuplexOptions;
        use crate::commands::group::GroupOptions;

        let mut spec = empty_spec(vec![Stage::Group, Stage::Duplex]);
        spec.stage_opts.group =
            Some(GroupOptions { strategy: Strategy::Identity, ..Default::default() });
        spec.stage_opts.duplex = Some(DuplexOptions::default());

        let err = validate_cross_stage_constraints(&spec).unwrap_err();
        assert!(err.to_string().contains("Paired"), "got: {err}");
    }

    #[cfg(feature = "consensus")]
    #[test]
    fn duplex_after_group_with_paired_strategy_accepted() {
        use crate::assigner::Strategy;
        use crate::commands::duplex::DuplexOptions;
        use crate::commands::group::GroupOptions;

        let mut spec = empty_spec(vec![Stage::Group, Stage::Duplex]);
        spec.stage_opts.group =
            Some(GroupOptions { strategy: Strategy::Paired, ..Default::default() });
        spec.stage_opts.duplex = Some(DuplexOptions::default());

        validate_cross_stage_constraints(&spec).expect("paired strategy should pass");
    }

    #[cfg(feature = "consensus")]
    #[test]
    fn standalone_duplex_skips_paired_check() {
        // Chain has only Stage::Duplex (no Stage::Group). Rule 2 does not fire
        // because the check requires both Stage::Group AND Stage::Duplex to be
        // present in the same chain.
        use crate::commands::duplex::DuplexOptions;

        let mut spec = empty_spec(vec![Stage::Duplex]);
        spec.stage_opts.duplex = Some(DuplexOptions::default());

        validate_cross_stage_constraints(&spec).expect("standalone duplex skips paired check");
    }

    // ── Rule 4: Stage::Extract requires SourceSpec::Fastqs ──────────────────

    /// Helper: a minimal `ExtractOptions` suitable for tests.
    fn minimal_extract_opts() -> crate::commands::extract::ExtractOptions {
        use crate::commands::extract::{ExtractOptions, QualityEncoding};
        ExtractOptions {
            sample: "sample".to_string(),
            library: "library".to_string(),
            platform: None,
            platform_unit: None,
            read_group_id: "A".to_string(),
            comments: vec![],
            barcode: None,
            platform_model: None,
            sequencing_center: None,
            predicted_insert_size: None,
            description: None,
            run_date: None,
            quality_encoding: QualityEncoding::Standard,
            store_umi_quals: false,
            store_cell_quals: false,
            single_tag: None,
            annotate_read_names: false,
            extract_umis_from_read_names: false,
            store_sample_barcode_qualities: false,
            async_reader: false,
        }
    }

    #[test]
    fn cross_stage_extract_with_bam_source_rejected() {
        // Rule 4: Stage::Extract + SourceSpec::Bam must be rejected.
        let mut spec = empty_spec(vec![Stage::Extract]);
        spec.source = SourceSpec::Bam(std::path::PathBuf::from("in.bam"));
        spec.stage_opts.extract = Some(minimal_extract_opts());
        let err = validate_cross_stage_constraints(&spec).unwrap_err();
        let msg = err.to_string();
        assert!(msg.contains("requires SourceSpec::Fastqs"), "got: {msg}");
    }

    #[test]
    fn cross_stage_extract_with_fastqs_source_accepted() {
        // Rule 4: Stage::Extract + SourceSpec::Fastqs must be accepted.
        use read_structure::ReadStructure;
        use std::str::FromStr;

        let mut spec = empty_spec(vec![Stage::Extract]);
        spec.source = SourceSpec::Fastqs {
            paths: vec![std::path::PathBuf::from("r1.fq.gz"), std::path::PathBuf::from("r2.fq.gz")],
            read_structures: vec![
                ReadStructure::from_str("+T").unwrap(),
                ReadStructure::from_str("+T").unwrap(),
            ],
        };
        spec.stage_opts.extract = Some(minimal_extract_opts());
        validate_cross_stage_constraints(&spec)
            .expect("Stage::Extract + SourceSpec::Fastqs should be valid");
    }

    #[test]
    fn opts_present_extract_with_options_passes() {
        // Stage::Extract with its bag slot populated must pass the opts-present check.
        use read_structure::ReadStructure;
        use std::str::FromStr;

        let mut spec = empty_spec(vec![Stage::Extract]);
        spec.source = SourceSpec::Fastqs {
            paths: vec![std::path::PathBuf::from("r1.fq.gz")],
            read_structures: vec![ReadStructure::from_str("+T").unwrap()],
        };
        spec.stage_opts.extract = Some(minimal_extract_opts());
        validate_stage_opts_present(&spec)
            .expect("Stage::Extract with options populated should pass");
    }
}
