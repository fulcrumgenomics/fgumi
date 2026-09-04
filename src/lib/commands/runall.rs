//! Fused multi-stage `runall` pipeline command.
//!
//! `runall` fuses a contiguous slice of the fgumi pipeline (extract, UMI
//! correction, alignment + zipper-merge, sort, group, and consensus calling,
//! optionally followed by filter) into a single in-memory chain, so every
//! intermediate BAM a sequential run of the standalone commands would write
//! is elided.
//!
//! This module currently defines the stage-selection vocabulary shared by the
//! rest of the command: [`RunAllMode`] (which consensus algorithm to run),
//! [`RunAllStage`] (the linear pipeline stage order used by `--start-from`),
//! and [`StopAfter`] (the `--stop-after`-only subset of [`RunAllStage`], which
//! excludes `align` because raw aligner output before the zipper-merge has
//! lost every original tag).

use anyhow::{Result, bail};

use crate::assigner::Strategy;

/// Consensus mode selector for `runall`. Controls which consensus
/// caller runs after the fused group + MI-grouping stages.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash, clap::ValueEnum)]
pub enum RunAllMode {
    /// Single-strand consensus via `VanillaUmiConsensusCaller`. Pair
    /// with `--strategy adjacency` (or `identity`/`edit`).
    Simplex,
    /// CODEC duplex consensus via `CodecConsensusCaller`. Pair with
    /// `--strategy adjacency` (or `identity`/`edit`).
    Codec,
    /// Two-strand duplex consensus via `DuplexConsensusCaller`. Requires
    /// `--strategy paired` so MIs carry `/A`/`/B` suffixes.
    Duplex,
}

impl std::fmt::Display for RunAllMode {
    /// Lower-case the variant name so error messages match the CLI
    /// flag values the user typed (`--consensus simplex`, not
    /// `RunAllMode::Simplex`).
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        let name = match self {
            Self::Simplex => "simplex",
            Self::Codec => "codec",
            Self::Duplex => "duplex",
        };
        f.write_str(name)
    }
}

/// Pipeline stage marker for the `--start-from` / `--stop-after`
/// runall stage flags. Stages are linearly ordered:
/// `Extract < Correct < AlignAndMerge < Zipper < Sort < Group <
/// Consensus < Filter` (see [`RunAllStage::ord`]).
///
/// **Stage semantics — "the input is positioned AT this stage's
/// entry":**
///
/// * `Sort` — raw unsorted BAM, pre-sort. `--start-from sort`
///   includes the sort step. `--stop-after sort` writes the
///   sorted BAM and stops (apples-to-apples replacement for
///   standalone `fgumi sort`).
/// * `Group` — BAM is sorted by template-coordinate (output of
///   `fgumi sort`). `--start-from group` skips the sort step.
///   `--stop-after group` writes the grouped BAM (output of
///   `fgumi group`) and stops.
/// * `Consensus` — BAM is sorted and grouped (MI tags assigned,
///   output of `fgumi group`). `--start-from consensus` skips sort
///   and group. `--stop-after consensus` runs the consensus caller
///   chosen by `--consensus {simplex,duplex,codec}` and writes the
///   consensus BAM (output of `fgumi simplex`/`duplex`/`codec`). The
///   algorithm is selected by `--consensus`, not by a stage variant.
///
/// `AlignAndMerge` is a valid `--start-from` but not a valid
/// `--stop-after` (see [`StopAfter`]): raw aligner output before the
/// zipper-merge has lost every original tag.
///
/// **Case-insensitive parsing**: `sort`, `Sort`, `SORT` all match
/// because the `--start-from` / `--stop-after` flags use
/// `#[arg(value_enum, ignore_case = true)]`.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash, clap::ValueEnum)]
#[clap(rename_all = "kebab-case")]
pub enum RunAllStage {
    /// Extract UMIs from FASTQ files and produce an unmapped BAM.
    /// `--start-from extract` accepts FASTQ files via
    /// `--extract::inputs` and read structures via
    /// `--extract::read-structures`. `--start-from extract --stop-after
    /// extract` runs the extract-only chain and writes the unmapped BAM
    /// (record-identical to standalone `fgumi extract`, with matching `@HD`
    /// `SO`/`GO`/`SS`, `@SQ`, and `@RG`; runall's `@PG` provenance differs, and
    /// `@HD` `VN` / `@CO` are not compared). Otherwise Extract feeds into
    /// Correct (when `--correct::umi-files` is set) or Align.
    Extract,
    /// Correct UMIs in an unmapped BAM against a fixed whitelist of
    /// expected sequences. `--start-from correct` accepts an
    /// unmapped BAM via `--input` and requires at least one of
    /// `--correct::umis` / `--correct::umi-files`. `--stop-after
    /// correct` writes the corrected unmapped BAM and exits.
    Correct,
    /// Align an unmapped BAM with UMIs through an external aligner
    /// subprocess (bwa-mem3 / bwa / user-supplied command), then merge
    /// the resulting alignments with the original unmapped reads via
    /// the zipper-merge logic — all streamed end-to-end. `--start-from
    /// align` accepts an unmapped BAM via `--input` and requires
    /// `--ref` (the aligner reference) plus one of
    /// `--aligner::preset` or `--aligner::command`. Earliest stop
    /// point reachable is `--stop-after zipper` (raw aligner output
    /// without merge is not a supported stop point because it loses
    /// every original tag).
    ///
    /// CLI value is `align` (short form for "align-and-merge" — the
    /// stage does both, but the user-facing name is the verb).
    #[clap(name = "align")]
    AlignAndMerge,
    /// Zip an unmapped BAM with an aligned BAM, transferring tags.
    /// `--start-from zipper` accepts a queryname-sorted unmapped BAM
    /// (via `--unmapped`) and a queryname-sorted mapped BAM (via
    /// `--input`); `--stop-after zipper` writes the merged BAM.
    /// Requires `--reference` (the same FASTA + `.dict` the standalone
    /// `fgumi zipper` requires).
    Zipper,
    /// Raw unsorted BAM. Includes the sort step when used as
    /// `--start-from sort`; writes sorted BAM and stops when used
    /// as `--stop-after sort`.
    Sort,
    /// Sorted BAM (template-coordinate). Includes group + MI tag
    /// assignment when used as `--start-from group`; writes
    /// grouped BAM and stops when used as `--stop-after group`.
    Group,
    /// Consensus calling. The algorithm (simplex / duplex / CODEC) is
    /// selected by `--consensus`, not by this stage. `--start-from
    /// consensus` accepts a grouped (MI-tagged) BAM; `--stop-after
    /// consensus` runs the chosen consensus caller and writes the
    /// consensus BAM (record-identical to `fgumi simplex` / `duplex` / `codec`,
    /// with matching `@HD` `SO`/`GO`/`SS`, `@SQ`, and `@RG`; runall's `@PG`
    /// provenance differs, and `@HD` `VN` / `@CO` are not compared).
    /// `--stop-after filter` fuses consensus → filter into one
    /// in-memory pipeline (no intermediate consensus BAM).
    Consensus,
    /// Post-consensus filtering. `--start-from filter` accepts a
    /// consensus BAM (output of `fgumi simplex`/`duplex`/`codec`) via
    /// `--input`; `--stop-after filter` runs the filter stage and writes
    /// the filtered BAM (record-identical to standalone `fgumi filter`, with
    /// matching `@HD` `SO`/`GO`/`SS`, `@SQ`, and `@RG`; runall's `@PG`
    /// provenance differs, and `@HD` `VN` / `@CO` are not compared).
    /// Tunable via `--filter::*`. Filter is the last stage in the linear
    /// order — no stage follows it. It runs either as the filter self-pair
    /// (`--start-from filter --stop-after filter`, input = a consensus
    /// BAM) or fused after an upstream consensus caller in the same
    /// pipeline (`… → consensus → filter`), where the consensus output is
    /// decoded in-memory and streamed straight into filter.
    Filter,
}

/// The subset of [`RunAllStage`] values that are valid as `--stop-after`.
///
/// `align` (`AlignAndMerge`) is a valid `--start-from` but **not** a valid
/// `--stop-after`: raw aligner output before the zipper-merge has lost
/// every original tag (RX, QX, RG, …), so the earliest post-align stop
/// point is `zipper`. Giving `--stop-after` its own value-enum makes
/// `--help` list only the genuinely-valid stops and rejects
/// `--stop-after align` at the clap parse layer with a clear
/// "invalid value 'align'" error — instead of parsing it and failing
/// later at runtime validation.
///
/// Convert to [`RunAllStage`] with `RunAllStage::from(stop)` /
/// `stop.into()`; the two enums are kept in sync by
/// `stop_after_values_are_runallstage_minus_align`.
#[derive(Debug, Clone, Copy, PartialEq, Eq, clap::ValueEnum)]
#[clap(rename_all = "kebab-case")]
pub enum StopAfter {
    /// `--stop-after extract` — extract-only (FASTQ → unmapped BAM).
    Extract,
    /// `--stop-after correct` — write the corrected unmapped BAM.
    Correct,
    /// `--stop-after zipper` — merged-but-unsorted BAM (earliest
    /// post-align stop).
    Zipper,
    /// `--stop-after sort` — template-coordinate sorted BAM.
    Sort,
    /// `--stop-after group` — grouped (MI-tagged) BAM.
    Group,
    /// `--stop-after consensus` — consensus BAM (algorithm chosen by
    /// `--consensus`); same output as `fgumi simplex`/`duplex`/`codec`.
    Consensus,
    /// `--stop-after filter` — post-consensus filtered BAM (same output
    /// as `fgumi filter`). Tunable via `--filter::*`.
    Filter,
}

impl From<StopAfter> for RunAllStage {
    fn from(stop: StopAfter) -> Self {
        match stop {
            StopAfter::Extract => Self::Extract,
            StopAfter::Correct => Self::Correct,
            StopAfter::Zipper => Self::Zipper,
            StopAfter::Sort => Self::Sort,
            StopAfter::Group => Self::Group,
            StopAfter::Consensus => Self::Consensus,
            StopAfter::Filter => Self::Filter,
        }
    }
}

impl RunAllStage {
    /// Position in the linear stage order. `Extract = 0`,
    /// `Correct = 1`, `AlignAndMerge = 2`, `Zipper = 3`, `Sort = 4`,
    /// `Group = 5`, `Consensus = 6`, `Filter = 7`. Used to validate that
    /// `start_from <= stop_after` and that the chain is buildable.
    #[must_use]
    pub fn ord(self) -> usize {
        match self {
            Self::Extract => 0,
            Self::Correct => 1,
            Self::AlignAndMerge => 2,
            Self::Zipper => 3,
            Self::Sort => 4,
            Self::Group => 5,
            Self::Consensus => 6,
            Self::Filter => 7,
        }
    }

    /// `true` if this stage is the consensus stage. Used by the Group
    /// step-derivation exception (a consensus-start input is already
    /// MI-tagged, so no explicit Group stage is added).
    #[must_use]
    pub fn is_consensus(self) -> bool {
        matches!(self, Self::Consensus)
    }

    /// `true` if this stage is the extract stage.
    #[must_use]
    pub fn is_extract(self) -> bool {
        matches!(self, Self::Extract)
    }

    /// Validate that `self` (as `--start-from`) and `stop_after` form
    /// a buildable pipeline range.
    ///
    /// This is the sole ordinal-ordering check for the `--start-from` /
    /// `--stop-after` pair: it enforces `self.ord() <= stop_after.ord()`,
    /// i.e. that the pipeline does not start later than it stops.
    ///
    /// # Errors
    ///
    /// Returns `Err` if `self.ord() > stop_after.ord()`.
    pub fn validate_with(self, stop_after: Self) -> Result<()> {
        if self.ord() > stop_after.ord() {
            bail!(
                "--start-from {self} comes after --stop-after {stop_after}; \
                 stop-after must be at or after start-from in the pipeline order \
                 (extract < correct < align < zipper < sort < group < consensus < filter)"
            );
        }
        Ok(())
    }
}

impl std::fmt::Display for RunAllStage {
    /// Lower-case the variant name so error messages match the CLI
    /// flag values the user typed (`--start-from sort`, not
    /// `--start-from Sort`).
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        let name = match self {
            Self::Extract => "extract",
            Self::Correct => "correct",
            Self::AlignAndMerge => "align",
            Self::Zipper => "zipper",
            Self::Sort => "sort",
            Self::Group => "group",
            Self::Consensus => "consensus",
            Self::Filter => "filter",
        };
        f.write_str(name)
    }
}

/// Validate that the `--strategy` chosen by the user is compatible
/// with the terminal consensus mode derived from `--stop-after`.
///
/// `Duplex` requires `Strategy::Paired` (MIs need `/A`/`/B` suffixes);
/// `Simplex` and `Codec` require any non-paired strategy. Error
/// messages name `--stop-after` (Display, lowercase) so the user
/// sees the flag they typed.
///
/// Free function so the error-message Display formatting is unit-testable
/// without constructing a full `RunAll`.
///
/// Not yet called outside tests: the `RunAll` struct and its
/// `validate_strategy` wrapper land in a later task of this same PR-B
/// campaign, at which point this has a production call site.
#[allow(dead_code)]
fn validate_strategy_for_mode(mode: RunAllMode, strategy: Strategy) -> Result<()> {
    match (mode, strategy) {
        (RunAllMode::Duplex, Strategy::Paired) => Ok(()),
        (RunAllMode::Duplex, _) => {
            bail!("--consensus duplex requires --strategy paired (so MIs carry /A and /B suffixes)")
        }
        (RunAllMode::Simplex | RunAllMode::Codec, Strategy::Paired) => {
            bail!("--consensus {mode} requires a non-paired strategy (identity/edit/adjacency)")
        }
        _ => Ok(()),
    }
}

/// Validate the `--start-from` / `--stop-after` pair against the
/// structural (linear-ordering) rules. Delegates to
/// [`RunAllStage::validate_with`]; consensus is no longer terminal-only
/// (it fuses into `consensus → filter`), so only the ordinal order is
/// enforced.
///
/// Historically this function also enforced "not yet implemented"
/// gates while #33 tasks 4-9 incrementally wired each combination;
/// those gates are all lifted now. The function is kept as the
/// runall-side entry point so future validation rules (e.g.
/// strategy/stage compatibility beyond what `validate_with`
/// expresses) have a clear home.
///
/// Not yet called outside tests: `RunAll::validate_stages` (a later task of
/// this same PR-B campaign) is the production call site.
#[allow(dead_code)]
fn validate_stages_for(start_from: RunAllStage, stop_after: RunAllStage) -> Result<()> {
    start_from.validate_with(stop_after)?;
    // `--stop-after align-and-merge` would mean "stop after the
    // alignment but before zipper-merge". The design dropped that
    // stop point because raw aligner output without the zipper-merge
    // loses every original tag (RX, QX, RG, etc.) and isn't useful.
    // Earliest stop point reachable from `--start-from align-and-merge`
    // is `--stop-after zipper`. `validate_with`'s linear-order check
    // would have let this combo through (ord 0 → ord 0 is structurally
    // valid as a self-pair), so we reject it explicitly here.
    if stop_after == RunAllStage::AlignAndMerge {
        bail!(
            "--stop-after align is not supported; raw aligner output \
             without zipper-merge loses every original tag (RX, QX, RG, ...). \
             Use `--stop-after zipper` for merged-but-unsorted output."
        );
    }
    Ok(())
}

/// Derive the ordered [`Stage`](crate::pipeline::chains::Stage) list for a
/// given `(start_from, stop_after)` pair.
///
/// Free function so the stage-derivation logic is unit-testable without
/// constructing a full `RunAll` instance — same pattern as
/// [`validate_stages_for`]. `RunAll::derive_stages` is the thin wrapper
/// that reads `self` and calls this.
///
/// # Stage derivation rules
///
/// 1. **Correct**: included when `start_from == Correct`.
///    When `stop_after > Correct`, the chain always chains through AAM
///    (`Correct → Align` is the only supported downstream path; a corrected
///    unmapped BAM must be aligned before it can be sorted or grouped).
///
/// 2. **Align**: included when `start_from <= AlignAndMerge` AND
///    `stop_after >= Zipper`. `Stage::Align` encapsulates both the aligner
///    subprocess AND the zipper-merge; there is no separate
///    `--stop-after align` stop point (the validator rejects it).
///    Included when:
///    - `start_from == AlignAndMerge` (explicit AAM start), OR
///    - `start_from == Correct` and `stop_after > Correct` (correct feeds
///      into AAM as the mandatory next step).
///
/// 3. **Zipper**: included when `start_from == Zipper` (standalone zipper
///    start). Mutually exclusive with `Stage::Align` — both produce a merged
///    BAM from different source shapes.
///
/// 4. **Sort**: included when `start_from <= Sort` OR when `Align` /
///    `Zipper` is in the chain AND `stop_after >= Sort`.
///    **Sort-forcing rule**: `Stage::Align` and `Stage::Zipper` both emit
///    queryname-sorted output. Downstream `Stage::Group` requires
///    template-coordinate order — Sort is therefore mandatory between any
///    Align/Zipper stage and Group, even when `--start-from` is not `sort`.
///
/// 5. **Group**: included when `stop_after >= Group` AND `start_from` is NOT
///    a consensus stage. **Consensus-self-pair exception**: when
///    `start_from.is_consensus()`, the input BAM is already MI-tagged. The
///    standalone consensus chain builders handle MI-tag grouping internally
///    via `GroupByMi` — adding an explicit `Stage::Group` here would
///    double-group an already-grouped input. For all other start stages
///    (correct, align, zipper, sort, group), the Group step is required to
///    produce MI-tagged output before the consensus caller sees it.
///
/// 6. **Simplex / Duplex / Codec**: one of the three is appended when
///    `stop_after` names that consensus stage (terminal, mutually exclusive).
///
/// # Errors
///
/// Returns `Err` for any `(start, stop)` pair that `validate_stages` would
/// reject (should not be reached if `validate_stages` was called first), or
/// if the Stage derivation reaches a branch that is a programming error (e.g.
/// a non-consensus `stop` that makes it through all prior early-returns).
///
/// Not yet called outside tests: `RunAll::derive_stages` (a later task of
/// this same PR-B campaign) is the production call site.
#[allow(dead_code)]
fn derive_stages_for(
    start_from: RunAllStage,
    stop_after: RunAllStage,
    extract_wants_correct: bool,
    consensus_mode: Option<RunAllMode>,
) -> Result<Vec<crate::pipeline::chains::Stage>> {
    use crate::pipeline::chains::Stage;

    // `validate_stages` must have been called before this; the assert
    // catches misuse in tests or future refactors.
    debug_assert!(
        start_from.ord() <= stop_after.ord() && stop_after != RunAllStage::AlignAndMerge,
        "derive_stages_for called with invalid (start={start_from}, stop={stop_after}); \
         validate_stages must run first"
    );

    let mut stages: Vec<Stage> = Vec::with_capacity(7);

    // ── Step 0: Extract ──────────────────────────────────────────────────
    // Included when the chain starts at Extract.
    if start_from == RunAllStage::Extract {
        stages.push(Stage::Extract);
        // extract → extract: extract-only chain (FASTQ → unmapped BAM).
        // The terminal Extract stage serialises and writes the unmapped
        // BAM directly; there are no downstream stages.
        if stop_after == RunAllStage::Extract {
            return Ok(stages);
        }
    }

    // Spec §6: an extract-start that stops at correct without a UMI source has
    // no buildable chain (Correct would be skipped, and the ordinal fall-through
    // would misreport "--consensus required"). Return the actionable error here
    // so the pure fn is total over every (start, stop, wants_correct, mode).
    if start_from == RunAllStage::Extract
        && stop_after == RunAllStage::Correct
        && !extract_wants_correct
    {
        bail!(
            "--start-from extract --stop-after correct requires \
             --correct::umi-files or --correct::umis"
        );
    }

    // ── Filter self-pair ─────────────────────────────────────────────────
    // `--start-from filter --stop-after filter`: input is a consensus BAM;
    // run only the (terminal) filter stage. No group/consensus steps apply
    // (the standalone filter chain handles its own queryname grouping when
    // `--filter::filter-by-template` is set), so short-circuit here. The
    // validator guarantees `start_from == Filter` implies
    // `stop_after == Filter` (filter is the last stage in the linear order).
    if start_from == RunAllStage::Filter {
        stages.push(Stage::Filter);
        return Ok(stages);
    }

    // ── Step 1: Correct ──────────────────────────────────────────────────
    // Included when the chain starts at Correct. Also included when the
    // chain starts at Extract AND the user supplied --correct::umi-files
    // or --correct::umis (the extract_wants_correct flag).
    // For a Correct self-pair there are no downstream stages; return early.
    if start_from == RunAllStage::Correct
        || (start_from == RunAllStage::Extract && extract_wants_correct)
    {
        stages.push(Stage::Correct);
        if stop_after == RunAllStage::Correct {
            return Ok(stages);
        }
        // `Correct → stop > Correct` always chains through AAM.
        // Any stop beyond Correct (Zipper, Sort, Group, consensus)
        // requires align-and-merge — a corrected unmapped BAM cannot
        // be fed directly to Sort (no query-coordinate BAM).
    }

    // ── Step 2: Align (AAM = align + zipper-merge, fused) ────────────────
    // `Stage::Align` covers both the aligner subprocess and the
    // zipper-merge; there is no separate `--stop-after align` stop point.
    // Included when:
    //   (a) start_from == AlignAndMerge (explicit AAM start), OR
    //   (b) start_from == Correct and stop_after > Correct (correct feeds
    //       into AAM — see Step 1 above).
    //   (c) start_from == Extract and stop_after > Correct (extract feeds
    //       into Align, possibly through Correct first).
    let includes_align = start_from == RunAllStage::AlignAndMerge
        || (start_from == RunAllStage::Correct && stop_after != RunAllStage::Correct)
        || (start_from == RunAllStage::Extract && stop_after.ord() > RunAllStage::Correct.ord());
    if includes_align {
        stages.push(Stage::Align);
        // `--stop-after zipper` stops after the AAM step (which includes
        // the zipper-merge internally). Return early.
        if stop_after == RunAllStage::Zipper {
            return Ok(stages);
        }
    }

    // ── Step 3: Zipper (standalone zipper start) ─────────────────────────
    // Mutually exclusive with `Stage::Align`. Included only when
    // `start_from == Zipper` (the user supplies separate mapped + unmapped
    // BAMs and the zipper-merge runs as the first step).
    let includes_zipper = start_from == RunAllStage::Zipper;
    if includes_zipper {
        stages.push(Stage::Zipper);
        if stop_after == RunAllStage::Zipper {
            return Ok(stages);
        }
    }

    // ── Step 4: Sort ──────────────────────────────────────────────────────
    // Sort-forcing rule: Sort is ALWAYS included when Align or Zipper
    // precede anything past Sort, because both produce queryname-sorted
    // output and downstream Group requires template-coordinate order.
    // Sort is also included when start_from == Sort (explicit sort start).
    let sort_forced =
        (includes_align || includes_zipper) && stop_after.ord() >= RunAllStage::Sort.ord();
    let sort_explicit = start_from == RunAllStage::Sort;
    if sort_forced || sort_explicit {
        stages.push(Stage::Sort);
    }
    if stop_after == RunAllStage::Sort {
        return Ok(stages);
    }

    // ── Step 5: Group ─────────────────────────────────────────────────────
    // Included when stop_after >= Group AND start_from is NOT a consensus
    // stage. The consensus-self-pair exception: when start_from is already
    // a consensus stage (Simplex/Duplex/Codec), the input BAM is already
    // MI-tagged. The standalone consensus chain builder (build_simplex_chain
    // / build_duplex_chain / build_codec_chain) handles MI-tag grouping
    // internally via GroupByMi — it does NOT run GroupByPosition /
    // ProcessGroups / MiAssign. Adding an explicit Group stage here would
    // double-group an already-grouped input.
    //
    // For `--start-from group --stop-after {consensus}` or any upstream
    // start (sort, align, zipper, correct), Group IS required because the
    // input is not yet MI-tagged.
    let includes_group = !start_from.is_consensus()
        && (start_from == RunAllStage::Group || stop_after.ord() >= RunAllStage::Group.ord());
    if includes_group {
        stages.push(Stage::Group);
    }
    if stop_after == RunAllStage::Group {
        return Ok(stages);
    }

    // ── Step 6: Consensus ────────────────────────────────────────────────
    // Reached when stop_after is Consensus or Filter (all earlier stops
    // have returned, and the filter self-pair short-circuited above). The
    // chain stage is chosen by the `--consensus` algorithm. When
    // stop_after == Filter, consensus runs first, then Filter is appended
    // below (consensus → filter chain).
    debug_assert!(
        matches!(stop_after, RunAllStage::Consensus | RunAllStage::Filter),
        "derive_stages_for reached the consensus branch with unexpected \
         stop stage {stop_after}; validate_stages must run first"
    );
    let mode = consensus_mode.ok_or_else(|| {
        anyhow::anyhow!(
            "--consensus <simplex|duplex|codec> is required when --stop-after \
             reaches the consensus stage"
        )
    })?;
    match mode {
        RunAllMode::Simplex => stages.push(Stage::Simplex),
        RunAllMode::Duplex => stages.push(Stage::Duplex),
        RunAllMode::Codec => stages.push(Stage::Codec),
    }

    // ── Step 7: Filter (terminal) ────────────────────────────────────────
    // Appended after the consensus caller when `--stop-after filter` chains
    // consensus → filter. Filter is the last stage in the linear order, so
    // no stage follows it.
    if stop_after == RunAllStage::Filter {
        stages.push(Stage::Filter);
    }

    Ok(stages)
}

#[cfg(test)]
mod derive_tests {
    use super::*;
    use crate::assigner::Strategy;
    use crate::pipeline::chains::Stage;
    use RunAllStage::*;
    use rstest::rstest;

    #[rstest]
    // ── self-pairs ──
    #[case::extract_only(Extract, Extract, false, None, vec![Stage::Extract])]
    #[case::correct_only(Correct, Correct, false, None, vec![Stage::Correct])]
    #[case::sort_only(Sort, Sort, false, None, vec![Stage::Sort])]
    #[case::group_only(Group, Group, false, None, vec![Stage::Group])]
    #[case::filter_only(Filter, Filter, false, None, vec![Stage::Filter])]
    // ── extract-fed fan-out ──
    #[case::extract_to_group_no_umi(
        Extract, Group, false, None,
        vec![Stage::Extract, Stage::Align, Stage::Sort, Stage::Group])]
    #[case::extract_to_group_with_umi(
        Extract, Group, true, None,
        vec![Stage::Extract, Stage::Correct, Stage::Align, Stage::Sort, Stage::Group])]
    #[case::extract_to_zipper_allowed(
        Extract, Zipper, false, None,
        vec![Stage::Extract, Stage::Align])] // §6: ALLOW; align satisfies a zipper stop
    #[case::extract_to_simplex(
        Extract, Consensus, false, Some(RunAllMode::Simplex),
        vec![Stage::Extract, Stage::Align, Stage::Sort, Stage::Group, Stage::Simplex])]
    // ── correct-fed ──
    #[case::correct_to_sort(
        Correct, Sort, false, None,
        vec![Stage::Correct, Stage::Align, Stage::Sort])]
    // ── sort/group/consensus fed ──
    #[case::sort_to_group(Sort, Group, false, None, vec![Stage::Sort, Stage::Group])]
    #[case::sort_to_duplex(
        Sort, Consensus, false, Some(RunAllMode::Duplex),
        vec![Stage::Sort, Stage::Group, Stage::Duplex])]
    #[case::group_to_codec(
        Group, Consensus, false, Some(RunAllMode::Codec),
        vec![Stage::Group, Stage::Codec])]
    #[case::group_to_consensus_to_filter(
        Group, Filter, false, Some(RunAllMode::Simplex),
        vec![Stage::Group, Stage::Simplex, Stage::Filter])]
    // ── consensus self-pair: NO Stage::Group (already MI-tagged) ──
    #[case::consensus_self_pair(
        Consensus, Consensus, false, Some(RunAllMode::Simplex),
        vec![Stage::Simplex])]
    // ── zipper start ──
    #[case::zipper_to_sort(Zipper, Sort, false, None, vec![Stage::Zipper, Stage::Sort])]
    fn derive_stages_matches_spec(
        #[case] start: RunAllStage,
        #[case] stop: RunAllStage,
        #[case] wants_correct: bool,
        #[case] mode: Option<RunAllMode>,
        #[case] expected: Vec<Stage>,
    ) {
        let got = derive_stages_for(start, stop, wants_correct, mode).expect("derives");
        assert_eq!(got, expected);
    }

    #[rstest]
    #[case::extract_correct_no_umi(Extract, Correct, false, None)]
    fn derive_stages_errors(
        #[case] start: RunAllStage,
        #[case] stop: RunAllStage,
        #[case] wants_correct: bool,
        #[case] mode: Option<RunAllMode>,
    ) {
        let err = derive_stages_for(start, stop, wants_correct, mode).unwrap_err().to_string();
        assert!(err.contains("--correct::umi-files or --correct::umis"), "got: {err}");
    }

    #[rstest]
    #[case::backwards(Group, Sort)]
    fn validate_stages_rejects_backwards(#[case] start: RunAllStage, #[case] stop: RunAllStage) {
        assert!(validate_stages_for(start, stop).is_err());
    }

    #[rstest]
    #[case::duplex_needs_paired(RunAllMode::Duplex, Strategy::Adjacency, true)]
    #[case::duplex_paired_ok(RunAllMode::Duplex, Strategy::Paired, false)]
    #[case::simplex_paired_rejected(RunAllMode::Simplex, Strategy::Paired, true)]
    #[case::simplex_adjacency_ok(RunAllMode::Simplex, Strategy::Adjacency, false)]
    fn strategy_for_mode(#[case] mode: RunAllMode, #[case] strat: Strategy, #[case] is_err: bool) {
        assert_eq!(validate_strategy_for_mode(mode, strat).is_err(), is_err);
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn runallstage_ord_is_linear() {
        assert!(RunAllStage::Extract.ord() < RunAllStage::Correct.ord());
        assert!(RunAllStage::Correct.ord() < RunAllStage::AlignAndMerge.ord());
        assert!(RunAllStage::AlignAndMerge.ord() < RunAllStage::Zipper.ord());
        assert!(RunAllStage::Zipper.ord() < RunAllStage::Sort.ord());
        assert!(RunAllStage::Sort.ord() < RunAllStage::Group.ord());
        assert!(RunAllStage::Group.ord() < RunAllStage::Consensus.ord());
        assert!(RunAllStage::Consensus.ord() < RunAllStage::Filter.ord());
    }

    #[test]
    fn stopafter_maps_to_runallstage_minus_align() {
        // StopAfter has exactly RunAllStage minus AlignAndMerge; every StopAfter
        // round-trips to a RunAllStage that is NOT AlignAndMerge.
        for stop in [
            StopAfter::Extract,
            StopAfter::Correct,
            StopAfter::Zipper,
            StopAfter::Sort,
            StopAfter::Group,
            StopAfter::Consensus,
            StopAfter::Filter,
        ] {
            let s: RunAllStage = stop.into();
            assert_ne!(s, RunAllStage::AlignAndMerge);
        }
    }

    #[test]
    fn display_matches_cli_values() {
        assert_eq!(RunAllStage::AlignAndMerge.to_string(), "align");
        assert_eq!(RunAllStage::Consensus.to_string(), "consensus");
        assert_eq!(RunAllMode::Simplex.to_string(), "simplex");
    }

    #[test]
    fn validate_with_rejects_start_after_stop() {
        assert!(RunAllStage::Group.validate_with(RunAllStage::Sort).is_err());
        assert!(RunAllStage::Sort.validate_with(RunAllStage::Group).is_ok());
    }
}
