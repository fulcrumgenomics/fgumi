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

use std::path::PathBuf;

use anyhow::{Result, bail};
use clap::Parser;

use crate::assigner::Strategy;
use crate::commands::common::{
    CompressionOptions, QueueMemoryOptions, ReadGroupOptions, RejectsOptions, SchedulerOptions,
    StatsOptions, ThreadingOptions,
};

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
/// `RunAll::validate_stages` is the thin wrapper that calls this.
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
/// `RunAll::derive_stages` is the thin wrapper that calls this.
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

/// Fused runall command. Stages are selected by `--start-from` and
/// `--stop-after`; the accepted starts are
/// `{Extract, Correct, AlignAndMerge, Zipper, Sort, Group, Consensus,
/// Filter}` (the consensus algorithm is chosen by `--consensus`, not by a
/// per-algorithm start value — `simplex`/`duplex`/`codec` are NOT accepted
/// start values) and the accepted stops are
/// `{Extract, Correct, Zipper, Sort, Group, Consensus, Filter}`
/// (no `--stop-after align` — raw aligner output without zipper-merge
/// would lose every original tag).
#[derive(Debug, Parser)]
#[command(
    name = "runall",
    about = "\x1b[38;5;166m[UTILITIES]\x1b[0m      \x1b[36mFused multi-stage pipeline (correct + align + zipper + sort + group + consensus, no intermediate BAM)\x1b[0m",
    long_about = "\
Fuses extract, UMI correction, alignment + zipper-merge, sort, group, and \
consensus calling (optionally followed by filter) into a single in-memory \
pipeline. Select the slice to run with `--start-from` and `--stop-after` (both \
required); the stages run in the fixed order extract < correct < align < zipper \
< sort < group < consensus < filter. Records flow directly between stages in \
memory, so every intermediate BAM a sequential run would write is elided. The \
record stream matches running the equivalent standalone commands in turn, as do \
the `@HD` `SO`/`GO`/`SS` fields, the `@SQ` dictionary, and the `@RG` records. The \
output is NOT byte-for-byte identical: the fused chain writes a single `@PG` \
record while the staged run writes one `@PG` per command, and `@HD` `VN` and \
`@CO` are not compared.

`--start-from` declares the state of the input and selects the upstream work:

* `extract` reads FASTQ via `--extract::inputs` / `--extract::read-structures` \
  (no `--input`).
* `correct` reads an unmapped BAM via `--input` and corrects UMIs against \
  `--correct::umis` / `--correct::umi-files`.
* `align` reads an unmapped BAM via `--input`, runs the aligner subprocess \
  (`--aligner::preset` or `--aligner::command`) against `--ref`, and \
  zipper-merges the result.
* `zipper` reads a queryname-sorted mapped BAM via `--input` plus the matching \
  unmapped BAM via `--unmapped` and the reference (`--ref` + `.dict`).
* `sort` reads a raw unsorted aligned BAM via `--input`.
* `group` reads a template-coordinate-sorted BAM via `--input`.
* `consensus` reads a grouped (MI-tagged) BAM via `--input`.
* `filter` reads a consensus BAM via `--input`.

The consensus algorithm is chosen by `--consensus {simplex,duplex,codec}`, not \
by a stage name; `--stop-after consensus` produces the same BAM as the matching \
`fgumi {simplex,duplex,codec}` command. `--stop-after filter` fuses a \
post-consensus `fgumi filter` pass (tunable via `--filter::*`) into the same \
pipeline, so no intermediate consensus BAM is written. Tune any stage with its \
prefixed flags: `--sort::max-memory`, `--group::strategy`, `--duplex::min-reads`, \
and so on.

`--stop-after align` is not allowed: raw aligner output before the zipper-merge \
has lost every original tag, so the earliest stop after `--start-from align` is \
`zipper`. The validator checks the requested stage range but cannot verify that \
the on-disk input actually matches `--start-from` — make sure it does."
)]
pub struct RunAll {
    /// Input BAM/SAM file. Optional because runall's input contract is
    /// stage-dependent, unlike the standard single-required-BAM commands
    /// that flatten `BamIoOptions`: `--start-from extract` reads FASTQ
    /// from `--extract::inputs` (no `--input`); `--start-from zipper`
    /// pairs `--input` (mapped) with `--unmapped`; the BAM-source start
    /// stages (`correct`/`align`/`sort`/`group`/consensus) require
    /// `--input`. Presence is validated per start-stage in
    /// `RunAll::derive_source_spec`.
    #[arg(short = 'i', long = "input")]
    pub input: Option<PathBuf>,

    /// Output BAM file.
    #[arg(short = 'o', long = "output")]
    pub output: PathBuf,

    /// Wrap the input in a userspace async prefetch reader: a background
    /// thread reads ahead so disk I/O overlaps decompression/compute.
    /// Defaults to off.
    #[arg(long = "async-reader", default_value_t = false, hide = true)]
    pub async_reader: bool,

    /// Optional output for rejected reads.
    #[command(flatten)]
    pub rejects_opts: RejectsOptions,

    /// Optional output for statistics.
    #[command(flatten)]
    pub stats_opts: StatsOptions,

    /// Read group and read name prefix options.
    #[command(flatten)]
    pub read_group: ReadGroupOptions,

    /// Threading options.
    #[command(flatten)]
    pub threading: ThreadingOptions,

    /// Compression options for output BAM.
    #[command(flatten)]
    pub compression: CompressionOptions,

    /// Scheduler and pipeline statistics options.
    #[command(flatten)]
    pub scheduler_opts: SchedulerOptions,

    /// Queue memory options.
    #[command(flatten)]
    pub queue_memory: QueueMemoryOptions,

    /// Methylation-aware consensus calling mode (requires --ref).
    #[arg(long = "methylation-mode", value_enum)]
    pub methylation_mode: Option<crate::commands::common::MethylationModeArg>,

    /// Path to the reference FASTA file. Required when `--methylation-mode`
    /// is set (consensus stages), when `--start-from zipper` is used (zipper
    /// requires the FASTA + `.dict` to build the output BAM header from the
    /// reference dictionary), and whenever the derived stage chain includes
    /// `Stage::Align` (the aligner reference) — not only `--start-from
    /// align`, but also fused runs that reach align from an upstream start
    /// (e.g. `--start-from extract`/`--start-from correct` with a
    /// `--stop-after` past zipper). For those align-bearing upstream starts
    /// the reference is consumed by the align stage and is accepted without
    /// `--methylation-mode`.
    #[arg(long = "ref")]
    pub reference: Option<PathBuf>,

    // ───────── extract-side options (used only when --start-from=extract) ─────────
    /// Per-stage extract tuning, exposed as `--extract::inputs`,
    /// `--extract::read-structures`, `--extract::sample`,
    /// `--extract::library`, etc. via the `MultiExtractRunallOptions`
    /// companion struct (generated by `#[multi_options]` on
    /// `ExtractRunallOptions` in `commands::extract`).
    /// Ignored when `--start-from` is not `extract`.
    #[command(flatten)]
    pub extract_opts: crate::commands::extract::MultiExtractRunallOptions,

    // ───────── correct-side options (used only when --start-from <= correct) ─────────
    /// Per-stage UMI-correction tuning, exposed as `--correct::umis`,
    /// `--correct::umi-files`, `--correct::max-mismatches`,
    /// `--correct::min-distance`, `--correct::revcomp`, etc. via the
    /// `MultiCorrectOptions` companion struct (generated by
    /// `#[multi_options]` on `CorrectOptions` in `commands::correct`).
    /// `--correct::min-distance` is required when `--start-from
    /// correct` is selected (no default); at least one of
    /// `--correct::umis` / `--correct::umi-files` is required.
    /// Ignored when `--start-from` is not `correct`.
    #[command(flatten)]
    pub correct_opts: crate::commands::correct::MultiCorrectOptions,

    // ───────── aligner-side options (used whenever the chain includes Stage::Align) ─────────
    /// Per-stage aligner tuning, exposed as `--aligner::preset`,
    /// `--aligner::command`, `--aligner::threads`, `--aligner::chunk-size`
    /// via the `MultiAlignerOptions` companion struct (generated by
    /// `#[multi_options]` on `AlignerOptions` in `crate::aligner`).
    /// Used whenever the derived stage chain includes `Stage::Align` — not
    /// only on `--start-from=align`, but also on fused runs that reach align
    /// from an upstream start (e.g. `--start-from extract` or
    /// `--start-from correct` with a `--stop-after` past zipper). Ignored
    /// only when the chain has no align stage.
    #[command(flatten)]
    pub aligner_opts: crate::aligner::MultiAlignerOptions,

    /// Override path for the aligner binary (preset mode only). When
    /// unset, the preset's binary (`bwa-mem3` / `bwa`) is found via
    /// `which::which()` on `PATH`. Rejected with a clear error if
    /// `--aligner::command` is used (command mode owns its own
    /// binary). Used whenever the chain includes `Stage::Align` (see
    /// `--aligner::*` above); ignored only when the chain has no align stage.
    ///
    /// Note: this flag is at the runall top level (`--aligner-bin`),
    /// NOT inside the `--aligner::*` family. The `--aligner::*` flags
    /// are mode-orthogonal tuning knobs (chunk-size, threads); the
    /// binary identity is mode-shaped — command mode owns it inside
    /// `--aligner::command "..."`, preset mode owns it here. Putting
    /// the override in the `--aligner::*` family would suggest a
    /// non-existent `--aligner::bin` symmetry.
    #[arg(long = "aligner-bin")]
    pub aligner_bin: Option<PathBuf>,

    // ───────── zipper-side options (used only when --start-from=zipper) ─────────
    /// Path to the unmapped BAM (required when `--start-from zipper`).
    /// The standalone `fgumi zipper` reads this via `-u`/`--unmapped`;
    /// runall surfaces it as a top-level flag because the
    /// `BamIoOptions`-shared `--input` already names the mapped BAM.
    /// Ignored when `--start-from` is not `zipper`.
    #[arg(long = "unmapped")]
    pub unmapped: Option<PathBuf>,

    /// Per-stage zipper tuning, exposed as `--zipper::tags-to-remove`,
    /// `--zipper::buffer`, `--zipper::skip-tc-tags`, etc. via the
    /// `MultiZipperOptions` companion struct (generated by
    /// `#[multi_options]` on `ZipperOptions` in `commands::zipper`).
    /// Ignored when `--start-from` is not `zipper`.
    #[command(flatten)]
    pub zipper_opts: crate::commands::zipper::MultiZipperOptions,

    // ───────── group-side options ─────────
    /// Per-stage group tuning, exposed as `--group::strategy`,
    /// `--group::edits`, `--group::min-map-q`, etc. via the
    /// `MultiGroupOptions` companion struct (generated by
    /// `#[multi_options]` on `GroupOptions` in `commands::group`).
    /// `--group::strategy` is required when `--start-from {sort,group}`
    /// is selected (no default).
    #[command(flatten)]
    pub group_opts: crate::commands::group::MultiGroupOptions,

    // ───────── codec-specific options (used only when --consensus=codec) ─────────
    /// Per-stage codec tuning, exposed as `--codec::min-duplex-length`,
    /// `--codec::single-strand-qual`, `--codec::outer-bases-qual`,
    /// `--codec::outer-bases-length`, `--codec::max-duplex-disagreements`,
    /// `--codec::max-duplex-disagreement-rate`. Ignored when `--consensus`
    /// is not `codec`.
    ///
    /// `--codec::allow-unmapped` is NOT exposed: PR A implemented spec
    /// §12.2 option (b) (`CodecOptions` carries a skipped, always-disabled
    /// `AllowUnmappedOptions`), so every runall codec stage runs with
    /// `allow_unmapped: false` — a known, deferred capability gap versus
    /// the standalone `fgumi codec` command.
    #[cfg(feature = "consensus")]
    #[command(flatten)]
    pub codec_opts: crate::commands::codec::MultiCodecOptions,

    // ───────── simplex-specific options (used only with --consensus simplex) ─────────
    /// Per-stage simplex tuning, exposed as `--simplex::error-rate-pre-umi`,
    /// `--simplex::min-reads`, `--simplex::max-reads`,
    /// `--simplex::min-input-base-quality`, etc. via the
    /// `MultiSimplexOptions` companion struct (generated by
    /// `#[multi_options]` on `SimplexOptions` in `commands::simplex`).
    /// `--simplex::min-reads` is required when `--consensus simplex` is
    /// selected (no default). Ignored unless the chain reaches the
    /// consensus stage with `--consensus simplex`.
    ///
    /// `--simplex::allow-unmapped` is NOT exposed: PR A implemented spec
    /// §12.2 option (b) (`SimplexOptions` carries a skipped, always-disabled
    /// `AllowUnmappedOptions`), so every runall simplex stage runs with
    /// `allow_unmapped: false` — a known, deferred capability gap versus
    /// the standalone `fgumi simplex` command.
    #[cfg(feature = "consensus")]
    #[command(flatten)]
    pub simplex_opts: crate::commands::simplex::MultiSimplexOptions,

    // ───────── duplex-specific options (used only with --consensus duplex) ─────────
    /// Per-stage duplex tuning, exposed as `--duplex::error-rate-pre-umi`,
    /// `--duplex::min-reads`, `--duplex::max-reads-per-strand`,
    /// `--duplex::consensus-call-overlapping-bases`, etc. via the
    /// `MultiDuplexOptions` companion struct. Ignored unless the chain
    /// reaches the consensus stage with `--consensus duplex`.
    ///
    /// `--duplex::allow-unmapped` is NOT exposed: PR A implemented spec
    /// §12.2 option (b) (`DuplexOptions` carries a skipped, always-disabled
    /// `AllowUnmappedOptions`), so every runall duplex stage runs with
    /// `allow_unmapped: false` — a known, deferred capability gap versus
    /// the standalone `fgumi duplex` command.
    #[cfg(feature = "consensus")]
    #[command(flatten)]
    pub duplex_opts: crate::commands::duplex::MultiDuplexOptions,

    // ───────── pipeline stage control ─────────
    /// Pipeline stage to start from: extract, correct, align, zipper, sort,
    /// group, consensus, or filter (required, case-insensitive).
    ///
    /// Declares the state of the input and selects the upstream work
    /// runall performs before the terminal `--stop-after` stage. See the
    /// command description for what each start stage reads (FASTQ for
    /// `extract`; `--input` plus, for `zipper`, `--unmapped`; a
    /// template-coordinate-sorted BAM for `group`; a grouped MI-tagged
    /// BAM for `consensus`; and so on).
    ///
    /// **Consensus-start chain:** `--start-from consensus` (with
    /// `--consensus {simplex,duplex,codec}`) runs the standalone consensus
    /// chain, which groups by the *existing* `MI` tag. It does NOT add a
    /// position-grouping stage — doing so would double-group an
    /// already-grouped input (see `derive_stages_for`'s
    /// consensus-self-pair exception).
    ///
    /// Case-insensitive. Must satisfy `start_from.ord() <=
    /// stop_after.ord()` (see [`RunAllStage::validate_with`]).
    #[arg(long = "start-from", value_enum, ignore_case = true)]
    pub start_from: RunAllStage,

    /// Pipeline stage to stop after: extract, correct, zipper, sort, group,
    /// consensus, or filter (required, case-insensitive).
    ///
    /// Determines the terminal stage of the chain. `--stop-after
    /// consensus` writes the consensus BAM produced by the caller chosen
    /// with `--consensus`; `--stop-after filter` fuses a post-consensus
    /// filter pass into the same pipeline.
    ///
    /// `align` is intentionally not a valid stop: raw aligner output
    /// before the zipper-merge has lost every original tag, so the
    /// earliest stop after `--start-from align` is `zipper`. Must satisfy
    /// `start_from.ord() <= stop_after.ord()` (`RunAllStage::validate_with`);
    /// use `RunAll::stop_after_stage` to read this as a [`RunAllStage`].
    #[arg(long = "stop-after", value_enum, ignore_case = true)]
    pub stop_after: StopAfter,

    /// Consensus algorithm to run at the consensus stage: `simplex`,
    /// `duplex`, or `codec`. Required whenever the chain reaches the
    /// consensus stage (`--stop-after consensus`, or `--start-from
    /// consensus`); ignored otherwise. `duplex` requires
    /// `--strategy paired`. Selects which standalone caller's output the
    /// run reproduces (`fgumi simplex` / `duplex` / `codec`).
    #[arg(long = "consensus", value_enum, ignore_case = true)]
    pub consensus_mode: Option<RunAllMode>,

    /// Per-stage sort tuning, exposed as `--sort::max-memory`,
    /// `--sort::tmp-dir`, etc. Consumed only when `--start-from sort`
    /// is selected; otherwise ignored. The full standalone option
    /// set of `fgumi sort` is mirrored here via the
    /// `#[multi_options]` macro — see `SortOptions` in
    /// `commands::sort` for the field list.
    #[command(flatten)]
    pub sort_opts: crate::commands::sort::MultiSortOptions,

    /// Per-stage filter tuning, exposed as `--filter::min-reads`,
    /// `--filter::max-read-error-rate`, `--filter::ref`, … Used when the
    /// chain reaches the filter stage (`--stop-after filter` or
    /// `--start-from filter`); otherwise ignored. The full standalone
    /// option set of `fgumi filter` is mirrored here via the
    /// `#[multi_options]` macro — see `FilterOptions` in
    /// `commands::filter` for the field list.
    #[command(flatten)]
    pub filter_opts: crate::commands::filter::MultiFilterOptions,
}

impl RunAll {
    /// Returns `true` if the user passed enough `--correct::*` flags to
    /// make Correct a load-bearing stage in an extract-fed chain.
    /// Concretely: at least one of `--correct::umi-files` or
    /// `--correct::umis` is non-empty. Used by `derive_stages` to
    /// decide whether to splice Correct between Extract and Align.
    #[must_use]
    pub(crate) fn extract_chain_wants_correct(&self) -> bool {
        !self.correct_opts.correct_umis.is_empty()
            || !self.correct_opts.correct_umi_files.is_empty()
    }

    /// Returns the consensus algorithm for this run, from `--consensus`.
    /// Errors (does not panic) if `--consensus` was omitted — it is required
    /// whenever the chain reaches the consensus stage. Only called on chains
    /// that actually reach consensus; `--stop-after sort`/`group`/`filter` are
    /// accepted stops that may not require it.
    ///
    /// Not yet called outside tests: the production call site (the
    /// stage-options-bag builder) lands in a later task of this same PR-B
    /// campaign.
    #[allow(dead_code)]
    pub(crate) fn require_consensus_mode(&self) -> Result<RunAllMode> {
        self.consensus_mode.ok_or_else(|| {
            anyhow::anyhow!(
                "--consensus <simplex|duplex|codec> is required when the chain \
                 reaches the consensus stage"
            )
        })
    }

    /// Reads `--stop-after` as a [`RunAllStage`] (`StopAfter` converts via
    /// [`From<StopAfter> for RunAllStage`]).
    #[must_use]
    pub(crate) fn stop_after_stage(&self) -> RunAllStage {
        self.stop_after.into()
    }

    /// Returns `self.input`, or an actionable error naming `--start-from`
    /// when it is absent. Used by every BAM-source `--start-from` stage
    /// (`correct`/`align`/`zipper`(mapped)/`sort`/`group`/`consensus`/`filter`).
    ///
    /// # Errors
    ///
    /// Returns `Err` if `--input` was not supplied.
    ///
    /// Not yet called outside tests: the production call sites (`execute`
    /// and the stage-options-bag builder) land in a later task of this same
    /// PR-B campaign.
    #[allow(dead_code)]
    pub(crate) fn require_input(&self) -> Result<PathBuf> {
        self.input.clone().ok_or_else(|| {
            anyhow::anyhow!(
                "--input is required with --start-from {} (the input BAM)",
                self.start_from
            )
        })
    }

    /// Returns the reference FASTA to thread into the consensus caller's
    /// `#[arg(skip)]` `reference` slot, gated on `--methylation-mode` being
    /// set: `Some(self.reference.clone())` when methylation mode is
    /// requested, `None` otherwise (even if `--ref` was supplied — e.g. to
    /// feed an upstream align stage in the same chain).
    ///
    /// Not yet called outside tests: the production call site (the
    /// stage-options-bag builder) lands in a later task of this same PR-B
    /// campaign.
    #[must_use]
    #[allow(dead_code)]
    pub(crate) fn consensus_reference(&self) -> Option<PathBuf> {
        if self.methylation_mode.is_some() { self.reference.clone() } else { None }
    }

    /// Validate the `--start-from` / `--stop-after` pair. Thin wrapper over
    /// the free function [`validate_stages_for`] so the ordering rules stay
    /// unit-testable without constructing a full `RunAll`.
    ///
    /// # Errors
    ///
    /// Returns `Err` under the same conditions as [`validate_stages_for`].
    ///
    /// Not yet called outside tests: the production call site (`execute`)
    /// lands in a later task of this same PR-B campaign.
    #[allow(dead_code)]
    pub(crate) fn validate_stages(&self) -> Result<()> {
        validate_stages_for(self.start_from, self.stop_after_stage())
    }

    /// Derive the ordered stage list for this run. Thin wrapper over the
    /// free function [`derive_stages_for`] so the derivation rules stay
    /// unit-testable without constructing a full `RunAll`.
    ///
    /// # Errors
    ///
    /// Returns `Err` under the same conditions as [`derive_stages_for`].
    ///
    /// Not yet called outside tests: the production call site (`execute`)
    /// lands in a later task of this same PR-B campaign.
    #[allow(dead_code)]
    pub(crate) fn derive_stages(&self) -> Result<Vec<crate::pipeline::chains::Stage>> {
        derive_stages_for(
            self.start_from,
            self.stop_after_stage(),
            self.extract_chain_wants_correct(),
            self.consensus_mode,
        )
    }

    /// Build the [`SourceSpec`](crate::pipeline::chains::SourceSpec) for the
    /// `ChainSpec` derived from `--start-from`.
    ///
    /// | `--start-from`     | `SourceSpec` variant                                              |
    /// |--------------------|-------------------------------------------------------------------|
    /// | `extract`          | `Fastqs`/`InterleavedFastq` from `--extract::inputs`              |
    /// | `correct`          | `Bam(self.require_input()?)` — unmapped BAM feeding correct       |
    /// | `align`            | `Bam(self.require_input()?)` — unmapped BAM feeding AAM           |
    /// | `zipper`           | `PairedBams { unmapped, mapped: input, reference }`               |
    /// | `sort` / `group` / consensus / `filter` | `Bam(self.require_input()?)` — aligned BAM |
    ///
    /// # Errors
    ///
    /// - `Extract` start with `--input` set is rejected (FASTQ paths come
    ///   from `--extract::inputs`, not `-i`/`--input`).
    /// - `Extract` start requires a non-empty `--extract::inputs` and
    ///   `--extract::read-structures` of matching length.
    /// - `Zipper` start requires `--unmapped` and `--ref`; returns `Err` if
    ///   either is absent.
    ///
    /// Exercised directly by `source_tests` below; the production call site
    /// (`execute`) lands in a later task of this same PR-B campaign.
    #[allow(dead_code)]
    pub(crate) fn derive_source_spec(&self) -> Result<crate::pipeline::chains::SourceSpec> {
        use crate::pipeline::chains::SourceSpec;
        match self.start_from {
            RunAllStage::Extract => {
                // When --start-from extract, FASTQ paths come from
                // --extract::inputs, not from -i/--input.
                if self.input.is_some() {
                    bail!(
                        "when --start-from extract, pass --extract::inputs instead of \
                         -i/--input (--input is for BAM-source start-from values)"
                    );
                }
                let opts = self.extract_opts.clone().validate()?;
                anyhow::ensure!(
                    !opts.inputs.is_empty(),
                    "--extract::inputs is required when --start-from extract"
                );
                anyhow::ensure!(
                    !opts.read_structures.is_empty(),
                    "--extract::read-structures is required when --start-from extract"
                );
                anyhow::ensure!(
                    opts.inputs.len() == opts.read_structures.len(),
                    "--extract::inputs and --extract::read-structures must have the same count"
                );
                if opts.interleaved {
                    let path = opts.inputs.first().cloned().ok_or_else(|| {
                        anyhow::anyhow!("--extract::inputs is required when --start-from extract")
                    })?;
                    // Enforces exactly 2 read structures (the R1/R2 pair).
                    SourceSpec::interleaved_fastq(path, opts.read_structures)
                } else {
                    // Enforces paths.len() == read_structures.len().
                    SourceSpec::fastqs(opts.inputs, opts.read_structures)
                }
            }

            RunAllStage::Correct
            | RunAllStage::AlignAndMerge
            | RunAllStage::Sort
            | RunAllStage::Group
            | RunAllStage::Consensus
            | RunAllStage::Filter => Ok(SourceSpec::Bam(self.require_input()?)),

            RunAllStage::Zipper => {
                let unmapped = self.unmapped.clone().ok_or_else(|| {
                    anyhow::anyhow!(
                        "--unmapped is required with --start-from zipper (the second \
                         input BAM, queryname-sorted, that carries the source tags)"
                    )
                })?;
                let reference = self.reference.clone().ok_or_else(|| {
                    anyhow::anyhow!(
                        "--ref is required with --start-from zipper (the reference FASTA \
                         with an accompanying `.dict` for the output BAM header)"
                    )
                })?;
                Ok(SourceSpec::PairedBams { unmapped, mapped: self.require_input()?, reference })
            }
        }
    }

    /// Build the [`SinkSpec`](crate::pipeline::chains::SinkSpec) for the
    /// `ChainSpec`. Runall always writes a plain BAM via `SinkSpec::Bam` —
    /// `BamWithIndex` is never constructed here because runall's sort output
    /// is either intermediate or a final template-coordinate BAM, neither of
    /// which admits a valid BAI. Runall does not expose a
    /// `--sort::write-index` flag.
    ///
    /// # Errors
    ///
    /// Currently infallible; returns `Result` for API symmetry with the
    /// other `derive_*` helpers.
    ///
    /// Not yet called outside tests: the production call site (`execute`)
    /// lands in a later task of this same PR-B campaign.
    #[allow(dead_code)]
    pub(crate) fn derive_sink_spec(&self) -> Result<crate::pipeline::chains::SinkSpec> {
        Ok(crate::pipeline::chains::SinkSpec::Bam(self.output.clone()))
    }
}

#[cfg(test)]
mod source_tests {
    use super::*;
    use crate::pipeline::chains::SourceSpec;
    use clap::Parser;

    fn parse(args: &[&str]) -> RunAll {
        RunAll::try_parse_from(std::iter::once("runall").chain(args.iter().copied())).unwrap()
    }

    #[test]
    fn bam_start_uses_input() {
        let r = parse(&[
            "--start-from",
            "group",
            "--stop-after",
            "group",
            "-i",
            "in.bam",
            "-o",
            "out.bam",
            "--group::strategy",
            "adjacency",
        ]);
        assert!(matches!(r.derive_source_spec().unwrap(), SourceSpec::Bam(_)));
    }

    #[test]
    fn extract_start_rejects_input() {
        let r = parse(&[
            "--start-from",
            "extract",
            "--stop-after",
            "extract",
            "-i",
            "in.bam",
            "-o",
            "out.bam",
            "--extract::inputs",
            "r1.fq",
            "--extract::read-structures",
            "+T",
            "--extract::sample",
            "s",
            "--extract::library",
            "l",
        ]);
        let err = r.derive_source_spec().unwrap_err().to_string();
        assert!(err.contains("--extract::inputs instead of"), "got: {err}");
    }

    #[test]
    fn zipper_start_requires_unmapped_and_ref() {
        let r = parse(&[
            "--start-from",
            "zipper",
            "--stop-after",
            "zipper",
            "-i",
            "mapped.bam",
            "-o",
            "out.bam",
        ]);
        assert!(r.derive_source_spec().unwrap_err().to_string().contains("--unmapped is required"));
    }
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
