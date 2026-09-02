//! Three independent validators for [`super::ChainSpec`]:
//!
//! 1. [`validate_stage_progression`] — ordering rules and mutual exclusions.
//! 2. [`validate_stage_opts_present`] — each stage has its options in the bag.
//! 3. [`validate_cross_stage_constraints`] — Zipper source, Duplex/Group
//!    strategy, `BamWithIndex` terminal-sort, and Extract source rules.

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
        // Terminal BAM→FASTQ export. Only ever appears as the sole stage of a
        // `fgumi fastq` chain, so its rank is nominal; place it last.
        Stage::Fastq => 9,
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
    if let Some(extract_pos) = stages.iter().position(|s| *s == Stage::Extract)
        && extract_pos != 0
    {
        bail!("Stage::Extract must be the first stage in the chain (got position {extract_pos})");
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
/// As of Phase 5 T5.3, thirteen stages have their options in the bag:
/// Correct, Sort, Group, Zipper, Duplex, Codec, Dedup, Filter, Clip, Simplex,
/// Align, Extract, and Fastq. The remaining stage (Downsample) adds its slot
/// incrementally during its migration task (T2.19–T2.22). For now that stage
/// skips the options-presence check — once its slot lands in the bag, add the
/// check here.
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
            Stage::Fastq => bag.fastq.is_some(),
            // Downsample is rejected by the early `bail!` above (no bag slot
            // yet). Report "not provided" rather than `unreachable!()` so that
            // if the early guard is ever removed without adding a slot here, the
            // validator returns the clean "options not provided" error instead
            // of panicking on a valid `[Downsample]` spec (S5b2-003).
            Stage::Downsample => false,
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
/// rule is checked at chain-build time (the chain builder reads the
/// input BAM header) — it cannot be determined from `spec.stages` alone,
/// so it intentionally stays out of this validator.
///
/// # Errors
///
/// Returns an error on the first violated constraint.
pub fn validate_cross_stage_constraints(spec: &ChainSpec) -> Result<()> {
    use crate::assigner::Strategy;
    use crate::commands::sort::SortOrderArg;
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
    // chain-build time (the chain builder reads the input BAM header) — it
    // can't be checked from spec.stages alone, so it stays out of this
    // validator by design.
    if spec.stages.contains(&Stage::Duplex)
        && spec.stages.contains(&Stage::Group)
        && let Some(group_opts) = spec.stage_opts.group.as_ref()
        && !matches!(group_opts.strategy, Strategy::Paired)
    {
        bail!(
            "Stage::Duplex requires Stage::Group to use Strategy::Paired \
             (got Strategy::{:?})",
            group_opts.strategy
        );
    }

    // Rule 3: `SinkSpec::BamWithIndex` requires the terminal stage to be
    // `Stage::Sort`, sorting by *coordinate*. The BAI file format indexes BGZF
    // virtual offsets of a coordinate-sorted BAM; emitting one for a non-sort
    // terminal stage (or for a chain where sort is intermediate, leaving some
    // other stage as terminal) would either reference a file that isn't
    // coordinate-sorted or attempt to index an output type the hook isn't
    // equipped to handle. A queryname- or template-coordinate-sorted terminal
    // is just as wrong: the BAI would bin records that are not in coordinate
    // order, producing a structurally-invalid index.
    //
    // The per-command CLI check (`Sort::execute` rejects `--write-index` for a
    // non-coordinate `--order`) covers the sort command, and `add_sink` carries
    // a `debug_assert!` on the produced @HD SO tag — but that assert is skipped
    // in release. Enforcing the order here, unconditionally and before any sink
    // is constructed, is what protects a *direct* chain (built programmatically
    // rather than through the sort command) from writing an invalid BAI in a
    // release build.
    if matches!(spec.sink, SinkSpec::BamWithIndex(_)) {
        let terminal_is_sort = spec.stages.last().is_some_and(|s| *s == Stage::Sort);
        if !terminal_is_sort {
            // `validate_stage_progression` rejects empty `stages` before
            // this validator runs, so `last()` is guaranteed `Some` here.
            // Rather than `unreachable!()`, the invariant is documented with a
            // descriptive error via `ok_or_else` so an out-of-order validator
            // change degrades to a clean error instead of a panic.
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

        // Terminal is `Stage::Sort`; require it to sort by coordinate. Only
        // check when the sort options are present — `validate_stage_opts_present`
        // (which `build_for` runs before this validator) guarantees they are on
        // the build path, so the `None` arm is reachable only when this
        // validator is called in isolation, and skipping it there mirrors how
        // Rule 2 guards on `if let Some(group_opts)`.
        if let Some(sort_opts) = spec.stage_opts.sort.as_ref()
            && !matches!(sort_opts.order, SortOrderArg::Coordinate)
        {
            bail!(
                "SinkSpec::BamWithIndex requires the terminal Stage::Sort to sort by coordinate \
                 (the BAI format indexes a coordinate-sorted BAM); got sort order {:?}",
                sort_opts.order,
            );
        }
    }

    // Rule 3b: structural preconditions on `SinkSpec::BamWithIndex` that
    // `add_sink` also enforces (as `ensure!` guards) but only *after* earlier
    // stages have run. `build_for` adds the sink last, so a rejects-bearing
    // stage such as `add_correct` has already created and truncated its rejects
    // `WriteBgzfFile` by the time `add_sink` would return one of these errors —
    // leaving a clobbered rejects file behind on an input that was invalid from
    // the start. Validating here, before any stage or sink is constructed (so
    // before any `File::create` truncates a file), preserves the
    // check-then-truncate discipline; the `add_sink` guards remain as defense in
    // depth. This affects `Correct → Sort` and `Correct → Align → Sort`.
    if let SinkSpec::BamWithIndex(output) = &spec.sink {
        // An inline BAI is a sidecar file joined against the primary BAM's
        // compressed offsets; it cannot be produced for a streamed stdout
        // output (there is no seekable file to index). Mirrors `add_sink`'s
        // `--write-index cannot target stdout` guard.
        if fgumi_bam_io::is_stdout_path(output) {
            bail!(
                "SinkSpec::BamWithIndex cannot target stdout; an inline BAI index requires a \
                 seekable on-disk output"
            );
        }

        // `Stage::Align` defers the output header (it is rewritten once the
        // aligner has emitted its @SQ lines), and the inline BAI indexer
        // requires a resolved (non-deferred) header. Mirrors `add_sink`'s
        // `inline BAI indexing requires a resolved (non-deferred) header`
        // guard, which fires on `self.pending_header_handle` — set only by
        // `add_align`.
        if spec.stages.contains(&Stage::Align) {
            bail!(
                "SinkSpec::BamWithIndex is incompatible with Stage::Align: alignment defers the \
                 output header, but inline BAI indexing requires a resolved (non-deferred) header"
            );
        }

        // The derived `.bai` sidecar must not collide with the primary BAM or
        // any rejects branch. The sidecar path is *derived* (`bai_sidecar_path`),
        // so the command-layer output-collision check never saw it: two
        // `WriteBgzfFile` sinks aimed at one path would otherwise clobber each
        // other byte for byte.
        reject_index_sidecar_collisions(
            output,
            &fgumi_bam_io::bai_sidecar_path(output),
            &spec_rejects_paths(spec),
        )?;
    }

    // Rule 4: `Stage::Extract` requires a `SourceSpec::Fastqs` source.
    // Extract reads directly from FASTQ files; feeding it a BAM, SAM, or
    // paired-BAM source is a programmer error that would be caught late
    // (at chain-build time) without this guard. Catching it here at
    // spec-validation time gives a clearer error message.
    if spec.stages.contains(&Stage::Extract) && !matches!(spec.source, SourceSpec::Fastqs { .. }) {
        bail!("Stage::Extract requires SourceSpec::Fastqs; got {:?}", spec.source);
    }

    // Rule 5: the FASTQ sink and `Stage::Fastq` must go together. A `Fastq`
    // sink pushes raw bytes through `add_fastq_sink` with no BAM header, and
    // `Stage::Fastq` emits FASTQ text, not serialized BAM records — pairing
    // either with the other format's counterpart would silently write a
    // corrupt file (FASTQ text wrapped in a BAM header, or BAM bytes with no
    // header). Require the biconditional.
    let sink_is_fastq = matches!(spec.sink, SinkSpec::Fastq(_) | SinkSpec::FastqPaired { .. });
    let terminal_is_fastq = spec.stages.last() == Some(&Stage::Fastq);
    if sink_is_fastq != terminal_is_fastq {
        bail!(
            "SinkSpec::Fastq requires Stage::Fastq as the terminal chain stage and vice versa; \
             got sink={:?}, terminal stage={:?}",
            spec.sink,
            spec.stages.last()
        );
    }

    // Rule 6: a `Stage::Sort` that immediately precedes a consensus stage
    // (Simplex/Duplex/Codec) must sort by template-coordinate. Consensus callers
    // group a molecule's reads by *adjacency* on the MI tag (`GroupByMi`), which
    // only merges runs of consecutive identical MI values. In a fused
    // `Sort → <consensus>` chain the consensus stage's own
    // `check_consensus_sort_order` guard is skipped — it fires only when the
    // consensus stage is first and consumes the raw source — so a coordinate- or
    // queryname-ordered intermediate sort would scatter a molecule's reads across
    // the stream and silently split it into multiple consensus groups (a
    // success-exit correctness bug). Template-coordinate is the only order that
    // guarantees same-MI reads are contiguous. No CLI builds this chain today;
    // like Rules 3/3b this protects a directly-constructed programmatic chain.
    for window in spec.stages.windows(2) {
        if window[0] == Stage::Sort
            && window[1].is_consensus()
            && let Some(sort_opts) = spec.stage_opts.sort.as_ref()
            && !matches!(sort_opts.order, SortOrderArg::TemplateCoordinate)
        {
            bail!(
                "Stage::Sort immediately preceding consensus stage {:?} must sort by \
                 template-coordinate (consensus groups reads by adjacency on the MI tag, so \
                 same-MI reads must be contiguous); got sort order {:?}",
                window[1],
                sort_opts.order,
            );
        }
    }

    Ok(())
}

/// The rejects/secondary output paths a chain's stage options declare.
///
/// These are the writable files a chain produces besides its primary sink;
/// [`validate_cross_stage_constraints`] checks the derived `.bai` sidecar
/// against them (see [`reject_index_sidecar_collisions`]). Only stages that
/// actually emit a rejects branch carry one — the consensus stages are
/// `consensus`-gated, matching their `StageOptionsBag` slots.
fn spec_rejects_paths(spec: &ChainSpec) -> Vec<&std::path::Path> {
    let mut paths: Vec<&std::path::Path> = Vec::new();
    let opts = &spec.stage_opts;
    if let Some(correct) = &opts.correct
        && let Some(rejects) = &correct.rejects_path
    {
        paths.push(rejects.as_path());
    }
    if let Some(filter) = &opts.filter
        && let Some(rejects) = &filter.rejects
    {
        paths.push(rejects.as_path());
    }
    #[cfg(feature = "consensus")]
    {
        if let Some(duplex) = &opts.duplex
            && let Some(rejects) = &duplex.rejects_opts.rejects
        {
            paths.push(rejects.as_path());
        }
        if let Some(codec) = &opts.codec
            && let Some(rejects) = &codec.rejects_opts.rejects
        {
            paths.push(rejects.as_path());
        }
        if let Some(simplex) = &opts.simplex
            && let Some(rejects) = &simplex.rejects_opts.rejects
        {
            paths.push(rejects.as_path());
        }
    }
    paths
}

/// Reject a collision between the inline-index `.bai` sidecar and any other file
/// the chain writes (the primary BAM or a rejects branch).
///
/// The sidecar path is *derived* ([`fgumi_bam_io::bai_sidecar_path`]), so the
/// command-layer
/// [`reject_output_collisions`](crate::commands::common::reject_output_collisions)
/// check — which only sees the paths the CLI passed — never inspects it. Two
/// `WriteBgzfFile` sinks aimed at one path would truncate and clobber each other
/// byte for byte, exactly the corruption `reject_output_collisions` prevents.
fn reject_index_sidecar_collisions(
    output_path: &std::path::Path,
    bai_path: &std::path::Path,
    rejects_paths: &[&std::path::Path],
) -> Result<()> {
    let mut targets: Vec<(&std::path::Path, &str)> =
        vec![(output_path, "--output"), (bai_path, "--write-index (.bai sidecar)")];
    targets.extend(rejects_paths.iter().map(|p| (*p, "--rejects")));
    crate::commands::common::reject_output_collisions(&targets)
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::commands::sort::{SortOptions, SortOrderArg};
    use crate::pipeline::chains::{SinkSpec, SourceSpec};
    use rstest::rstest;
    use std::path::PathBuf;

    /// A `SortOptions` for a given order, covering every field the type
    /// declares (it has no `Default`). Only `order` varies across callers —
    /// the rest are the standalone-sort defaults, matching the integration
    /// helper in `tests/integration/test_chain_bam_with_index.rs`.
    fn sort_opts_with_order(order: SortOrderArg) -> SortOptions {
        use crate::commands::common::{MaxTempFiles, MemoryLimit, MemoryReserve};
        SortOptions {
            order,
            key_types: None,
            max_memory: MemoryLimit::Auto,
            memory_reserve: MemoryReserve::Auto,
            memory_per_thread: true,
            tmp_dirs: Vec::new(),
            sort_threads: None,
            merge_threads: None,
            temp_compression: 1,
            temp_codec: fgumi_sort::SpillCodec::default(),
            max_temp_files: MaxTempFiles::Auto,
            block_batch: 4,
            file_granularity: false,
        }
    }

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
            read_streams: fgumi_bam_io::ReadStreams::Fixed(1),
            verify_crc: true,
            command_line: String::new(),
        }
    }

    #[test]
    fn empty_chain_rejected() {
        let err = validate_stage_progression(&empty_spec(vec![])).unwrap_err();
        assert!(err.to_string().contains("empty"));
    }

    #[test]
    fn is_sort_terminal_only_for_lone_sort() {
        // The single-sourced predicate that both `build_for` (source/sink skip)
        // and `ChainBuilder::new` (header-open skip) rely on.
        assert!(empty_spec(vec![Stage::Sort]).is_sort_terminal());
        assert!(!empty_spec(vec![]).is_sort_terminal());
        assert!(!empty_spec(vec![Stage::Group, Stage::Sort]).is_sort_terminal());
        assert!(!empty_spec(vec![Stage::Sort, Stage::Group]).is_sort_terminal());
        assert!(!empty_spec(vec![Stage::Correct, Stage::Sort]).is_sort_terminal());
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

    #[rstest]
    #[case::queryname_lexicographic(SortOrderArg::Queryname)]
    #[case::queryname_natural(SortOrderArg::QuerynameNatural)]
    #[case::template_coordinate(SortOrderArg::TemplateCoordinate)]
    fn cross_stage_bam_with_index_rejects_non_coordinate_sort_order(#[case] order: SortOrderArg) {
        // A terminal `Stage::Sort` under `SinkSpec::BamWithIndex` must sort by
        // coordinate — the BAI format indexes a coordinate-sorted BAM. This is
        // a plain `#[test]` (not `#[cfg(debug_assertions)]`), so it runs in
        // release too, covering the case the `add_sink` `debug_assert!` skips:
        // a direct chain built programmatically with a non-coordinate order.
        let mut spec = empty_spec(vec![Stage::Sort]);
        spec.sink = SinkSpec::BamWithIndex(PathBuf::from("out.bam"));
        spec.stage_opts.sort = Some(sort_opts_with_order(order));
        let msg = validate_cross_stage_constraints(&spec).unwrap_err().to_string();
        assert!(msg.contains("BamWithIndex"), "expected BamWithIndex in error, got: {msg}");
        assert!(msg.contains("coordinate"), "expected 'coordinate' in error, got: {msg}");
        assert!(
            msg.contains(&format!("{order:?}")),
            "expected the rejected order {order:?} in error, got: {msg}"
        );
    }

    #[test]
    fn cross_stage_bam_with_index_accepts_coordinate_sort_order() {
        // The one accepted order, with the sort options present (the `None`
        // arm is exercised by `cross_stage_bam_with_index_accepted_on_terminal_sort`).
        let mut spec = empty_spec(vec![Stage::Sort]);
        spec.sink = SinkSpec::BamWithIndex(PathBuf::from("out.bam"));
        spec.stage_opts.sort = Some(sort_opts_with_order(SortOrderArg::Coordinate));
        validate_cross_stage_constraints(&spec)
            .expect("BamWithIndex with a coordinate terminal sort is valid");
    }

    #[rstest]
    // A `Sort → <consensus>` fused chain skips the consensus stage's own
    // `check_consensus_sort_order` guard, so the validator must reject any
    // intermediate sort order that does not leave same-MI reads contiguous.
    // Cover every non-template-coordinate order against every consensus stage.
    fn cross_stage_sort_before_consensus_rejects_non_template_coordinate(
        #[values(Stage::Simplex, Stage::Duplex, Stage::Codec)] consensus: Stage,
        #[values(
            SortOrderArg::Coordinate,
            SortOrderArg::Queryname,
            SortOrderArg::QuerynameNatural
        )]
        order: SortOrderArg,
    ) {
        let mut spec = empty_spec(vec![Stage::Sort, consensus]);
        spec.stage_opts.sort = Some(sort_opts_with_order(order));
        let msg = validate_cross_stage_constraints(&spec).unwrap_err().to_string();
        assert!(
            msg.contains("template-coordinate"),
            "expected 'template-coordinate' in error, got: {msg}"
        );
        assert!(
            msg.contains(&format!("{consensus:?}")),
            "expected consensus stage {consensus:?} in error, got: {msg}"
        );
        assert!(
            msg.contains(&format!("{order:?}")),
            "expected the rejected order {order:?} in error, got: {msg}"
        );
    }

    #[rstest]
    fn cross_stage_sort_before_consensus_accepts_template_coordinate(
        #[values(Stage::Simplex, Stage::Duplex, Stage::Codec)] consensus: Stage,
    ) {
        // Template-coordinate is the one intermediate order that keeps same-MI
        // reads contiguous, so a `Sort(template-coordinate) → <consensus>` chain
        // passes the cross-stage validator.
        let mut spec = empty_spec(vec![Stage::Sort, consensus]);
        spec.stage_opts.sort = Some(sort_opts_with_order(SortOrderArg::TemplateCoordinate));
        validate_cross_stage_constraints(&spec).unwrap_or_else(|e| {
            panic!("Sort(template-coordinate) → {consensus:?} should be valid, got: {e}")
        });
    }

    #[test]
    fn cross_stage_non_adjacent_sort_before_consensus_is_unconstrained() {
        // `Sort → Group → Simplex`: Group re-groups by MI, so the intermediate
        // sort order does not feed the consensus stage directly and Rule 6 must
        // not fire even for a coordinate sort.
        let mut spec = empty_spec(vec![Stage::Sort, Stage::Group, Stage::Simplex]);
        spec.stage_opts.sort = Some(sort_opts_with_order(SortOrderArg::Coordinate));
        validate_cross_stage_constraints(&spec)
            .expect("Sort → Group → Simplex is not a direct Sort → consensus feed");
    }

    #[rstest]
    #[case::dash("-")]
    #[case::dev_stdout("/dev/stdout")]
    fn cross_stage_bam_with_index_rejects_stdout_output(#[case] output: &str) {
        // `Correct → Sort` with an indexed sink aimed at stdout. `add_correct`
        // would create and truncate its rejects file during stage construction,
        // before `add_sink` returns the stdout error — so validation must reject
        // this up front, before `build_for` constructs anything.
        let mut spec = empty_spec(vec![Stage::Correct, Stage::Sort]);
        spec.sink = SinkSpec::BamWithIndex(PathBuf::from(output));
        let msg = validate_cross_stage_constraints(&spec).unwrap_err().to_string();
        assert!(msg.contains("BamWithIndex"), "expected BamWithIndex in error, got: {msg}");
        assert!(msg.contains("stdout"), "expected 'stdout' in error, got: {msg}");
    }

    #[test]
    fn cross_stage_bam_with_index_rejects_align_stage() {
        // `Correct → Align → Sort` with an indexed sink. `Stage::Align` defers
        // the output header, which the inline BAI indexer cannot use;
        // `add_sink` rejects the deferred header only after `add_correct` has
        // truncated its rejects file, so validation must reject up front.
        let mut spec = empty_spec(vec![Stage::Correct, Stage::Align, Stage::Sort]);
        spec.sink = SinkSpec::BamWithIndex(PathBuf::from("out.bam"));
        let msg = validate_cross_stage_constraints(&spec).unwrap_err().to_string();
        assert!(msg.contains("BamWithIndex"), "expected BamWithIndex in error, got: {msg}");
        assert!(msg.contains("Align"), "expected 'Align' in error, got: {msg}");
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
        use crate::read_structure::ReadStructure;
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
        use crate::read_structure::ReadStructure;
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

    // ── Rule 5: SinkSpec::{Fastq,FastqPaired} <-> Stage::Fastq biconditional ─

    #[test]
    fn cross_stage_fastq_paired_accepted_on_terminal_fastq() {
        let mut spec = empty_spec(vec![Stage::Fastq]);
        spec.sink = SinkSpec::FastqPaired {
            out1: std::path::PathBuf::from("r1.fq.gz"),
            out2: std::path::PathBuf::from("r2.fq.gz"),
            out0: Some(std::path::PathBuf::from("other.fq.gz")),
        };
        validate_cross_stage_constraints(&spec)
            .expect("FastqPaired + terminal Stage::Fastq must be valid");
    }

    #[test]
    fn cross_stage_fastq_paired_rejected_on_non_fastq_terminal() {
        // Sort -> Group (Group is terminal, not Fastq).
        let mut spec = empty_spec(vec![Stage::Sort, Stage::Group]);
        spec.sink = SinkSpec::FastqPaired {
            out1: std::path::PathBuf::from("r1.fq"),
            out2: std::path::PathBuf::from("r2.fq"),
            out0: None,
        };
        let err = validate_cross_stage_constraints(&spec).unwrap_err();
        assert!(err.to_string().contains("Stage::Fastq"), "got: {err}");
    }

    #[test]
    fn cross_stage_fastq_terminal_rejected_on_non_fastq_sink() {
        // The other direction of the Rule 5 biconditional: a Stage::Fastq
        // terminal with a non-FASTQ sink (the default Bam sink) must be
        // rejected, mirroring the FASTQ-sink/non-Fastq-terminal case above.
        let spec = empty_spec(vec![Stage::Fastq]); // empty_spec defaults sink = Bam
        assert!(matches!(spec.sink, SinkSpec::Bam(_)), "precondition: default sink is Bam");
        let err = validate_cross_stage_constraints(&spec).unwrap_err();
        assert!(err.to_string().contains("Stage::Fastq"), "got: {err}");
    }

    // --- Rule 3b: inline-index `.bai` sidecar collision ---

    /// The `reject_index_sidecar_collisions` helper: a sidecar colliding with a
    /// rejects branch is rejected, and the error names both flags.
    #[test]
    fn sidecar_helper_rejects_collision_with_a_rejects_branch() {
        use std::path::Path;
        let output = Path::new("out.bam");
        let bai = Path::new("out.bam.bai");
        let rejects = [Path::new("out.bam.bai")];
        let err = super::reject_index_sidecar_collisions(output, bai, &rejects)
            .expect_err("sidecar/rejects collision must be rejected");
        let msg = err.to_string();
        assert!(
            msg.contains("--write-index") && msg.contains("--rejects"),
            "error must name both colliding flags, got: {msg}"
        );
    }

    /// The helper also rejects a sidecar colliding with the primary BAM output.
    #[test]
    fn sidecar_helper_rejects_collision_with_the_output() {
        use std::path::Path;
        let output = Path::new("out.bam");
        let err = super::reject_index_sidecar_collisions(output, output, &[])
            .expect_err("sidecar/output collision must be rejected");
        assert!(err.to_string().contains("--output"), "error must name --output");
    }

    /// The normal case — sidecar derived from the output, distinct or absent
    /// rejects — has no collision.
    #[test]
    fn sidecar_helper_accepts_distinct_paths() {
        use std::path::Path;
        let output = Path::new("out.bam");
        let bai = Path::new("out.bam.bai");
        assert!(super::reject_index_sidecar_collisions(output, bai, &[]).is_ok());
        assert!(
            super::reject_index_sidecar_collisions(output, bai, &[Path::new("rejects.bam")])
                .is_ok(),
            "a distinct rejects path must not collide with the sidecar"
        );
    }

    /// A `[Correct, Sort]` chain whose `--rejects` names the derived `.bai`
    /// sidecar is rejected at validation — before `build_for` opens any sink, so
    /// no file is truncated. Pins that Rule 3b is wired into the validator, not
    /// just tested as a free helper.
    #[test]
    fn cross_stage_bam_with_index_rejects_sidecar_rejects_collision() {
        let output = PathBuf::from("out.bam");
        let sidecar = fgumi_bam_io::bai_sidecar_path(&output);
        let mut spec = empty_spec(vec![Stage::Correct, Stage::Sort]);
        spec.sink = SinkSpec::BamWithIndex(output);
        spec.stage_opts.correct = Some(correct_opts_with_rejects(Some(sidecar)));
        let err = validate_cross_stage_constraints(&spec)
            .expect_err("a rejects path colliding with the derived .bai must be rejected");
        let msg = err.to_string();
        assert!(
            msg.contains("--write-index") && msg.contains("--rejects"),
            "error must name both colliding flags, got: {msg}"
        );
    }

    /// The same chain with a distinct `--rejects` path validates cleanly — the
    /// rule must not false-positive on a normal rejects-bearing indexed sort.
    #[test]
    fn cross_stage_bam_with_index_distinct_rejects_is_valid() {
        let mut spec = empty_spec(vec![Stage::Correct, Stage::Sort]);
        spec.sink = SinkSpec::BamWithIndex(PathBuf::from("out.bam"));
        spec.stage_opts.correct =
            Some(correct_opts_with_rejects(Some(PathBuf::from("rejects.bam"))));
        validate_cross_stage_constraints(&spec)
            .expect("a distinct rejects path must not collide with the derived .bai");
    }

    /// `CorrectOptions` carrying only the `rejects_path` under test; all other
    /// fields are inert defaults (mirrors `correct::tests::make_default_opts`).
    fn correct_opts_with_rejects(
        rejects_path: Option<PathBuf>,
    ) -> crate::commands::correct::CorrectOptions {
        use crate::commands::correct::{CorrectOptions, Target};
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
            rejects_path,
        }
    }
}
