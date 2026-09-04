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
use crate::commands::command::Command;
use crate::commands::common::{
    BamIoOptions, CompressionOptions, QueueMemoryOptions, ReadGroupOptions, RejectsOptions,
    SchedulerOptions, StatsOptions, ThreadingOptions,
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
/// `stopafter_maps_to_runallstage_minus_align`.
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
    /// Called by [`Self::build_stage_options_bag`] to populate the
    /// `#[arg(skip)]` `reference` slot on the Simplex/Duplex consensus
    /// options structs.
    #[must_use]
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
                // The equal-count check is only correct for the non-interleaved
                // (one FASTQ per read structure) shape. Interleaved packs both
                // reads of a pair into a single FASTQ, so it is 1 input + 2 read
                // structures. Mirror `build_stage_options_bag`'s Extract arm and
                // extract.rs so the two agree.
                if opts.interleaved {
                    anyhow::ensure!(
                        opts.inputs.len() == 1,
                        "--extract::interleaved requires exactly one --extract::inputs; got {}",
                        opts.inputs.len()
                    );
                    anyhow::ensure!(
                        opts.read_structures.len() == 2,
                        "--extract::interleaved requires exactly two \
                         --extract::read-structures; got {}",
                        opts.read_structures.len()
                    );
                    let path = opts.inputs.first().cloned().ok_or_else(|| {
                        anyhow::anyhow!("--extract::inputs is required when --start-from extract")
                    })?;
                    SourceSpec::interleaved_fastq(path, opts.read_structures)
                } else {
                    anyhow::ensure!(
                        opts.inputs.len() == opts.read_structures.len(),
                        "--extract::inputs and --extract::read-structures must have the same count"
                    );
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
    pub(crate) fn derive_sink_spec(&self) -> Result<crate::pipeline::chains::SinkSpec> {
        Ok(crate::pipeline::chains::SinkSpec::Bam(self.output.clone()))
    }

    /// Validate the align-stage CLI surface without executing the chain.
    /// Runs the same [`crate::aligner::AlignerOptions::resolve`] path the
    /// eventual align step uses, so a misconfigured invocation fails before
    /// any input is read AND produces the same error message users will see
    /// at execute time.
    ///
    /// Called from [`Command::execute`] whenever the derived stage list
    /// contains `Stage::Align`.
    ///
    /// # Errors
    ///
    /// - `--ref` not set (the aligner needs a reference path).
    /// - `--methylation-mode` combined with an align-bearing chain (not yet
    ///   supported; see the error message for the EM-seq workaround).
    /// - No sequence-dictionary file (`.dict`) found next to `--ref`.
    /// - `--aligner::preset` and `--aligner::command` both unset, or both
    ///   set, or a preset-only flag paired with command mode, or preset-mode
    ///   index files missing, or a command-mode template missing `{ref}`
    ///   (delegated to [`crate::aligner::AlignerOptions::resolve`]).
    pub(crate) fn validate_align_and_merge(&self) -> Result<()> {
        let reference = self.reference.as_ref().ok_or_else(|| {
            anyhow::anyhow!(
                "a runall chain that includes align requires --ref (the aligner \
                 reference FASTA with its index files alongside)"
            )
        })?;
        // `--methylation-mode` drives the *consensus* stage (which AAM never
        // reaches alone) and would also conflict with the aligner's
        // reference handling. Reject explicitly rather than silently ignore;
        // methylation-aware AAM presets are a follow-up PR per the design doc.
        if self.methylation_mode.is_some() {
            bail!(
                "--methylation-mode is not yet supported for runall chains that include \
                 align. For EM-seq today, use `--aligner::command \"bwameth.py ...\"` (or \
                 `bwa-mem3 --methylation-mode em-seq ...`) in command mode and apply \
                 methylation downstream as a separate step."
            );
        }
        // The downstream zipper-merge step needs a `.dict` file alongside
        // the reference FASTA (used to populate the output BAM header).
        // Validate up front so a missing `.dict` fails before the aligner
        // runs — saves potentially hours of CPU.
        if crate::reference::find_dict_path(reference).is_none() {
            bail!(
                "no sequence-dictionary file found next to --ref {} \
                 (expected `<ref>.dict` or `<ref-without-extension>.dict`). \
                 Generate one with `samtools dict {} -o <ref>.dict` before running \
                 a runall chain that includes align.",
                reference.display(),
                reference.display(),
            );
        }
        let top_threads = self.threading.num_threads();
        // `resolve` consumes the options struct; we clone so the validator
        // can be called multiple times if needed.
        let _resolved = self.aligner_opts.clone().validate()?.resolve(
            reference,
            top_threads,
            self.aligner_bin.as_deref(),
        )?;
        Ok(())
    }

    /// Build the [`StageOptionsBag`](crate::pipeline::chains::StageOptionsBag)
    /// for the `ChainSpec` derived from the active stages in
    /// [`Self::derive_stages`].
    ///
    /// Each active stage's options slot is populated; inactive slots remain
    /// `None`. Fields annotated `#[arg(skip)]` on the per-stage options
    /// structs are populated here using the same logic the existing
    /// single-stage `execute_*_only` paths apply.
    ///
    /// # Option-population rules by stage
    ///
    /// * **Correct** — validates `MultiCorrectOptions`, then sets
    ///   `rejects_path`:
    ///   - Self-pair (`--stop-after correct`): honors the top-level
    ///     `--rejects`, falling back to `--correct::rejects`.
    ///   - Cross-stage (correct feeds AAM): `rejects_path` = `None` (the
    ///     fused chain uses the kept-only correct step; UMI rejects are
    ///     discarded — `RunAll::execute` emits a warning when `--rejects` is
    ///     set on a chained correct run).
    ///
    /// * **Align** — constructs
    ///   [`AlignOptions`](crate::pipeline::chains::options_bag::AlignOptions)
    ///   from `--aligner::*` + `--aligner-bin` + `--ref`. Requires `--ref` to
    ///   be set.
    ///
    /// * **Zipper** — validates `MultiZipperOptions`. No `#[arg(skip)]` fields.
    ///
    /// * **Sort** — validates `MultiSortOptions`, forces
    ///   `order = TemplateCoordinate` (runall's sort step always feeds
    ///   downstream group; coordinate order is not applicable), resolves
    ///   `tmp_dirs` against `FGUMI_TMP_DIRS`, and — when sort is fused with
    ///   another stage — overrides an explicit `--sort::max-memory auto` to
    ///   the standalone 768 MiB default (auto-detection is not supported in a
    ///   fused chain). Finally runs the inherent
    ///   [`crate::commands::sort::SortOptions::validate`] cross-field check
    ///   (rejects `--sort::temp-codec zstd` with `--sort::temp-compression 0`)
    ///   before the chain builder wires the spill compressor.
    ///
    /// * **Group** — validates `MultiGroupOptions`, then computes
    ///   `effective_strategy` / `effective_edits` via
    ///   [`crate::commands::group::GroupOptions::resolve_strategy_and_edits`]
    ///   and clears the three histogram/metrics paths (anti-goal documented
    ///   in the module-level doc).
    ///
    /// * **Simplex** / **Duplex** — validate the per-mode `Multi<X>`, then
    ///   populate the cross-cutting `#[arg(skip)]` fields
    ///   (`rejects_opts`/`stats_opts`/`read_group`/`methylation_mode`/
    ///   `reference`) from `self`.
    ///
    /// * **Codec** — rejects `--methylation-mode` up front (codec has no
    ///   methylation support), validates `MultiCodecOptions`, runs the
    ///   inherent numeric-bounds `CodecOptions::validate`, then populates
    ///   `rejects_opts`/`stats_opts`/`read_group`.
    ///
    /// * **Extract** — validates `MultiExtractRunallOptions`, runs the
    ///   inherent `ExtractRunallOptions::validate` (the 2 macro-dropped
    ///   conflicts), then the interleaved-vs-count check, template-count
    ///   validation, non-empty-read-structure check, and
    ///   `ExtractOptions::validate` (reserved-tag + store-umi-quals).
    ///
    /// * **Filter** — validates `MultiFilterOptions`. No cross-stage
    ///   rewiring is needed: filter is always terminal in a runall chain.
    ///
    /// # Errors
    ///
    /// Returns `Err` if any per-stage `validate()` call fails (e.g. missing
    /// required option like `--group::strategy` or `--correct::min-distance`),
    /// if `Stage::Align` is present but `--ref` is absent, if
    /// `--methylation-mode` is combined with `Stage::Codec`, or if `stages`
    /// contains a stage `derive_stages` never emits for a runall chain (an
    /// internal-error bail naming the unexpected stage).
    pub(crate) fn build_stage_options_bag(
        &self,
        stages: &[crate::pipeline::chains::Stage],
    ) -> Result<crate::pipeline::chains::StageOptionsBag> {
        use crate::commands::common::MemoryLimit;
        use crate::commands::sort::{SortOrderArg, TMP_DIRS_ENV, resolve_tmp_dirs};
        use crate::pipeline::chains::options_bag::AlignOptions;
        use crate::pipeline::chains::{Stage, StageOptionsBag};

        let mut bag = StageOptionsBag::default();

        for &stage in stages {
            match stage {
                Stage::Correct => {
                    let mut opts = self.correct_opts.clone().validate()?;
                    opts.rejects_path = if self.stop_after_stage() == RunAllStage::Correct {
                        // self-pair: honor top-level --rejects, falling back to --correct::rejects.
                        self.rejects_opts.rejects.clone().or(opts.rejects_path)
                    } else {
                        // fused kept-only correct discards UMI rejects (see execute() warn)
                        None
                    };
                    bag.correct = Some(opts);
                }

                Stage::Align => {
                    let reference = self.reference.clone().ok_or_else(|| {
                        anyhow::anyhow!(
                            "a runall chain that includes align requires --ref (the aligner \
                             reference FASTA with its index files alongside)"
                        )
                    })?;
                    let aligner = self.aligner_opts.clone().validate()?;
                    bag.aligner = Some(AlignOptions {
                        aligner,
                        reference,
                        aligner_bin: self.aligner_bin.clone(),
                    });
                }

                Stage::Zipper => {
                    let opts = self.zipper_opts.clone().validate()?;
                    bag.zipper = Some(opts);
                }

                Stage::Sort => {
                    let mut sort_opts = self.sort_opts.clone().validate()?;
                    // Runall's sort step always produces template-coordinate
                    // output — the only order compatible with downstream Group.
                    // Warn before overriding a user-supplied order so the
                    // override is not silent (mirrors the max-memory warn below).
                    if sort_opts.order != SortOrderArg::TemplateCoordinate {
                        log::warn!(
                            "--sort::order {:?} is ignored; runall forces template-coordinate \
                             sort order for its sort step",
                            sort_opts.order
                        );
                    }
                    sort_opts.order = SortOrderArg::TemplateCoordinate;
                    sort_opts.tmp_dirs = resolve_tmp_dirs(
                        &sort_opts.tmp_dirs,
                        std::env::var(TMP_DIRS_ENV).ok().as_deref(),
                    );
                    // --sort::max-memory=auto bails in any non-sole-[Sort] chain
                    // (builder.rs:2560), so override auto -> the standalone
                    // 768 MiB default when sort is fused with other stages. Use
                    // the exact expression SortOptions::default() uses so the
                    // value cannot diverge from the standalone default.
                    if stages != [Stage::Sort] && sort_opts.max_memory == MemoryLimit::Auto {
                        log::warn!(
                            "--sort::max-memory auto is not supported in a fused runall \
                             chain; using 768M"
                        );
                        sort_opts.max_memory =
                            crate::commands::common::parse_memory("768M").expect("valid default");
                    }
                    // Inherent cross-field check (mirrors the Codec/Extract
                    // branches below): reject `--sort::temp-codec zstd` paired
                    // with `--sort::temp-compression 0` here, before the chain
                    // builder wires `SpillBlockCompress`, rather than failing
                    // lazily on the first spill. The macro-generated
                    // `MultiSortOptions::validate` only converts staged-required
                    // fields; it does not run this check.
                    sort_opts.validate()?;
                    bag.sort = Some(sort_opts);
                }

                Stage::Group => {
                    let mut group_opts = self.group_opts.clone().validate()?;
                    // Mirror GroupReadsByUmi::execute: compute
                    // effective_strategy / effective_edits via the shared
                    // helper, then null the per-position metrics outputs
                    // (anti-goal documented in the module-level doc comment).
                    let (effective_strategy, effective_edits) =
                        group_opts.resolve_strategy_and_edits();
                    group_opts.effective_strategy = effective_strategy;
                    group_opts.effective_edits = effective_edits;
                    // `--index-threshold always` asserts indexing will happen;
                    // reject it when the resolved strategy/edits can never index,
                    // matching standalone `fgumi group` (group.rs
                    // GroupReadsByUmi::validate) rather than silently ignoring it.
                    crate::commands::common::validate_index_threshold(
                        group_opts.index_threshold,
                        effective_strategy,
                        effective_edits,
                    )?;
                    // Forcing these to None is spec-mandated (runall does not
                    // emit per-position group metrics in a fused run), but warn
                    // so the dropped flags are not silent.
                    if group_opts.family_size_histogram.is_some()
                        || group_opts.grouping_metrics.is_some()
                        || group_opts.metrics_prefix.is_some()
                    {
                        log::warn!(
                            "--group::family-size-histogram / --group::grouping-metrics / \
                             --group::metrics are ignored: runall does not emit per-position \
                             group metrics in a fused run"
                        );
                    }
                    group_opts.family_size_histogram = None;
                    group_opts.grouping_metrics = None;
                    group_opts.metrics_prefix = None;
                    bag.group = Some(group_opts);
                }

                #[cfg(feature = "consensus")]
                Stage::Simplex => {
                    let mut opts = self.simplex_opts.clone().validate()?;
                    opts.rejects_opts.clone_from(&self.rejects_opts);
                    opts.stats_opts.clone_from(&self.stats_opts);
                    opts.read_group.clone_from(&self.read_group);
                    opts.methylation_mode =
                        crate::commands::common::resolve_methylation_mode(self.methylation_mode);
                    opts.reference = self.consensus_reference();
                    bag.simplex = Some(opts);
                }

                #[cfg(feature = "consensus")]
                Stage::Duplex => {
                    let mut opts = self.duplex_opts.clone().validate()?;
                    // `add_duplex` (builder.rs) does not validate numeric bounds
                    // the way `add_simplex`/codec do, so run the same guards the
                    // standalone `fgumi duplex` applies (DuplexOptions::validate_numeric,
                    // extracted from duplex.rs Duplex::validate) — otherwise runall
                    // accepts degenerate configs that yield a silent empty BAM.
                    opts.validate_numeric()?;
                    opts.rejects_opts.clone_from(&self.rejects_opts);
                    opts.stats_opts.clone_from(&self.stats_opts);
                    opts.read_group.clone_from(&self.read_group);
                    opts.methylation_mode =
                        crate::commands::common::resolve_methylation_mode(self.methylation_mode);
                    opts.reference = self.consensus_reference();
                    bag.duplex = Some(opts);
                }

                #[cfg(feature = "consensus")]
                Stage::Codec => {
                    // Codec does not support methylation calling — fail loud
                    // rather than silently dropping the flag.
                    if self.methylation_mode.is_some() {
                        bail!(
                            "--methylation-mode is not supported with codec consensus; \
                             it is valid only for simplex/duplex"
                        );
                    }
                    // Validate the `--codec::*` flags into CodecOptions, then
                    // run the numeric/semantic checks (min-reads, qual
                    // ceilings, disagreement-rate, …) that standalone `fgumi
                    // codec` runs — otherwise runall would accept degenerate
                    // configs the standalone command rejects.
                    let mut opts = self.codec_opts.clone().validate()?;
                    opts.validate()?;
                    opts.rejects_opts.clone_from(&self.rejects_opts);
                    opts.stats_opts.clone_from(&self.stats_opts);
                    opts.read_group.clone_from(&self.read_group);
                    bag.codec = Some(opts);
                }

                Stage::Extract => {
                    use crate::commands::extract::validate_template_count;

                    let opts = self.extract_opts.clone().validate()?; // -> ExtractRunallOptions
                    opts.validate()?; // the 2 macro-dropped conflicts
                    let rs = &opts.read_structures;
                    if opts.interleaved {
                        anyhow::ensure!(
                            opts.inputs.len() == 1,
                            "--extract::interleaved requires exactly one --extract::inputs; got {}",
                            opts.inputs.len()
                        );
                        anyhow::ensure!(
                            rs.len() == 2,
                            "--extract::interleaved requires exactly two \
                             --extract::read-structures; got {}",
                            rs.len()
                        );
                    } else {
                        anyhow::ensure!(
                            opts.inputs.len() == rs.len(),
                            "--extract::inputs and --extract::read-structures must have the \
                             same count"
                        );
                    }
                    validate_template_count(rs)?;
                    for (i, r) in rs.iter().enumerate() {
                        anyhow::ensure!(
                            !r.segments().is_empty(),
                            "Read structure {} is empty",
                            i + 1
                        );
                    }
                    let extract_options = opts.to_extract_options();
                    extract_options.validate()?; // reserved-tag + store-umi-quals
                    bag.extract = Some(extract_options);
                }

                Stage::Filter => {
                    // The standalone filter chain builder reads `rejects` and
                    // `stats` straight off the bag's FilterOptions, so
                    // `--filter::rejects` / `--filter::stats` flow through
                    // unchanged. No cross-stage rewiring is needed: filter is
                    // always terminal in a runall chain.
                    let filter_opts = self.filter_opts.clone().validate()?;
                    bag.filter = Some(filter_opts);
                }

                // Stages runall never derives — a programming error if they appear.
                Stage::Clip
                | Stage::Dedup
                | Stage::Downsample
                | Stage::Fastq
                | Stage::CopyUmi
                | Stage::Retag => bail!(
                    "internal error: build_stage_options_bag encountered unexpected \
                     stage {stage:?} in a runall chain; this is a bug in derive_stages"
                ),
                // When built without the consensus feature, the
                // Simplex/Duplex/Codec arms above are compiled out; keep the
                // match exhaustive (mirrors validate.rs).
                #[cfg(not(feature = "consensus"))]
                Stage::Simplex | Stage::Duplex | Stage::Codec => {
                    bail!("Stage {stage:?} requires building fgumi with the `consensus` feature")
                }
            }
        }

        Ok(bag)
    }
}

impl Command for RunAll {
    fn execute(&self, command_line: &str) -> Result<()> {
        use crate::pipeline::chains::{ChainSpec, Stage};

        let timer = crate::logging::OperationTimer::new(&format!(
            "runall (--start-from {} --stop-after {})",
            self.start_from,
            self.stop_after_stage()
        ));

        // stdin is supported: the input flows into `SourceSpec::Bam`/`Fastqs`
        // and is opened once via `InputSource::open` (BGZF/SAM auto-detected,
        // stdin-aware). runall reads the source exactly once per run, so a
        // streamed `-` works for every start stage.

        self.validate_stages()?;

        log::info!(
            "Starting runall (--start-from {} --stop-after {})",
            self.start_from,
            self.stop_after_stage()
        );
        if let Some(input) = &self.input {
            log::info!("Input: {}", input.display());
        }
        log::info!("Output: {}", self.output.display());

        // Derive the ordered stage list for this (start_from, stop_after) pair.
        let stages = self.derive_stages()?;

        // `--ref` pairs with `--methylation-mode` for the consensus stages, but
        // it is also legitimately consumed *without* methylation when the
        // derived chain includes `Stage::Align` (the aligner reference) or
        // starts at `zipper` (the FASTA + `.dict` for the output BAM header).
        // For any other chain `--ref` is dead, so reject it without
        // `--methylation-mode` rather than silently ignore it. Keyed on the
        // derived chain — not the start stage — so a fused `extract`/`correct →
        // … → align` run is exempted while a non-aligning chain (e.g.
        // `correct → correct`) is not. This must run after `derive_stages()`.
        let ref_is_consumed_without_methylation =
            self.start_from == RunAllStage::Zipper || stages.contains(&Stage::Align);
        if !ref_is_consumed_without_methylation
            && self.reference.is_some()
            && self.methylation_mode.is_none()
        {
            bail!("--ref requires --methylation-mode to be set");
        }

        // Symmetric guard for the *other* flag: in runall, `--methylation-mode`
        // is wired only into the simplex/duplex consensus stages (see the
        // `simplex_opts`/`duplex_opts` bag population in
        // `build_stage_options_bag`; codec and align reject it with their own
        // messages). On a chain that reaches no consensus stage (e.g.
        // `group → group`, `correct → sort`) it is dead — silently ignored —
        // so reject it rather than mislead. Align chains are exempt here:
        // `validate_align_and_merge` (run just below for any chain containing
        // `Stage::Align`) owns the align-specific EM-seq message, which is
        // more actionable than this generic one.
        let chain_reaches_consensus = stages.iter().any(|s| s.is_consensus());
        let chain_includes_align = stages.contains(&Stage::Align);
        // `--methylation-mode` on a simplex/duplex chain requires `--ref`: the
        // consensus stage consumes the FASTA (via `consensus_reference()`), and
        // the standalone contract (`--methylation-mode` doc, "requires --ref")
        // demands it. Without this guard a non-align `group → simplex` run with
        // `--methylation-mode` but no `--ref` slips through (the dead-consensus
        // guard below does not fire because the chain *does* reach consensus),
        // and `build_stage_options_bag()` then threads `methylation_mode =
        // Some(...)` with `reference = None` into the stage — deferring the
        // failure to a later, less actionable path. Codec is excluded (it
        // rejects `--methylation-mode` outright below), and align chains are
        // exempt (`validate_align_and_merge` owns the more specific message).
        let chain_uses_methylation_capable_consensus =
            stages.iter().any(|s| matches!(s, Stage::Simplex | Stage::Duplex));
        if self.methylation_mode.is_some()
            && chain_uses_methylation_capable_consensus
            && !chain_includes_align
            && self.reference.is_none()
        {
            bail!("--methylation-mode requires --ref to be set");
        }
        if self.methylation_mode.is_some() && !chain_reaches_consensus && !chain_includes_align {
            bail!(
                "--methylation-mode is only consumed by the consensus stage; \
                 it is dead on a runall chain that stops before consensus"
            );
        }

        // D12: log auto-inserted stages for extract-fed chains so users
        // understand the chain expansion rules.
        if self.start_from == RunAllStage::Extract {
            if stages.contains(&Stage::Correct) {
                log::info!(
                    "Including Stage::Correct between Extract and downstream \
                     (because --correct::umi-files or --correct::umis is set)"
                );
            }
            if stages.contains(&Stage::Align) {
                log::info!(
                    "Including Stage::Align between Extract and downstream \
                     (FASTQ source requires alignment)"
                );
            }
            if stages.contains(&Stage::Sort) && self.start_from != RunAllStage::Sort {
                log::info!(
                    "Including Stage::Sort between Align and Group \
                     (template-coordinate sort precedes Group)"
                );
            }
        }

        // Align-stage pre-flight: run `validate_align_and_merge` so that
        // aligner option conflicts (preset vs command, missing --ref,
        // --methylation-mode with AAM) are caught before the pipeline starts.
        if stages.contains(&Stage::Align) {
            self.validate_align_and_merge()?;
        }

        // Cross-stage validation: when Group feeds a consensus stage, the
        // group strategy must be compatible with the consensus mode (e.g.
        // duplex requires `--strategy paired`).
        let has_group = stages.contains(&Stage::Group);
        let has_consensus = stages.iter().any(|s| s.is_consensus());
        if has_group && has_consensus {
            let consensus_mode = self.require_consensus_mode()?;
            let group_opts = self.group_opts.clone().validate()?;
            validate_strategy_for_mode(consensus_mode, group_opts.strategy)?;
            // Log group/consensus summary banner (mirrors old fused dispatcher).
            log::info!("Strategy: {:?}, edits: {}", group_opts.strategy, group_opts.edits);
        }

        // Surface a warning if --rejects is set with --start-from
        // extract/correct and a chained --stop-after. The fused chain uses
        // the kept-only correct step (no UMI-rejects branch), so the rejects
        // file collects only downstream rejects (post-consensus / etc.), not
        // correct's UMI rejects. Users who need UMI rejects should stage a
        // `fgumi correct --rejects` step separately.
        if matches!(self.start_from, RunAllStage::Extract | RunAllStage::Correct)
            && self.stop_after_stage() != RunAllStage::Correct
            && stages.contains(&Stage::Correct)
            && self.rejects_opts.rejects.is_some()
        {
            log::warn!(
                "--rejects with --start-from {} --stop-after {} discards correct's UMI rejects \
                 (the fused chain has no UMI-rejects branch). The rejects file will collect only \
                 downstream rejects. Use a separate `fgumi correct --rejects` step if you need \
                 UMI rejects captured.",
                self.start_from,
                self.stop_after_stage()
            );
        }

        // A7: warn on top-level --stats / --rejects that the derived chain
        // never consumes, so a silently-dead flag does not mislead. Top-level
        // --stats is cloned only into consensus stages; top-level --rejects is
        // wired only into a correct self-pair or a consensus stage. Filter reads
        // its own --filter::stats / --filter::rejects, so point the user there.
        // Use warn (not bail) — the run is still valid, just missing the output
        // the user expected on the wrong flag.
        if self.stats_opts.stats.is_some() && !chain_reaches_consensus {
            log::warn!(
                "--stats is consumed only by the consensus stage; it is dead on a runall chain \
                 that reaches no consensus stage. Use --filter::stats to capture filter statistics."
            );
        }
        if self.rejects_opts.rejects.is_some()
            && !chain_reaches_consensus
            && !stages.contains(&Stage::Correct)
        {
            log::warn!(
                "--rejects is wired nowhere on this runall chain (it is consumed only by a correct \
                 self-pair or a consensus stage). Use --filter::rejects to capture filter rejects."
            );
        }

        // Validate the input BAM exists, once every earlier guard (stage
        // ordering, --ref/--methylation-mode, align pre-flight) has already
        // passed — those checks are pure CLI-argument reasoning and should
        // fail before we require any file to be present on disk.
        // `--start-from extract` has no `--input` (FASTQ comes from
        // `--extract::inputs`, validated in `derive_source_spec`), so this is
        // conditional.
        if let Some(input) = &self.input {
            BamIoOptions::new(input, &self.output).validate()?;
        }

        // Construct the ChainSpec from the derived stages + spec-construction
        // helpers. Build `source`/`sink`/`stage_opts` into locals FIRST — the
        // bag borrows `stages`, so `stages` cannot be moved into the struct
        // literal before that borrow is done with.
        let verify_crc = match &self.input {
            Some(p) => BamIoOptions::new(p, &self.output).effective_check_crc(),
            // FASTQ source ignores verify_crc; value is inert.
            None => true,
        };
        let source = self.derive_source_spec()?;
        let sink = self.derive_sink_spec()?;
        let stage_opts = self.build_stage_options_bag(&stages)?;

        // A1: reject an --output that collides with any rejects/stats path the
        // derived chain will actually open as a *second* writer. Two writers on
        // one file silently corrupt it (see `reject_output_collisions` in
        // common.rs). Enumerate only the writer paths the bag actually wired —
        // reading them off the populated `stage_opts` so the guard fires on the
        // exact paths the chain will open. `--input`/`--unmapped` are read-only
        // and exempt; stdin/`-`/`/dev/null` are handled inside the guard.
        let mut write_targets: Vec<(&std::path::Path, &str)> =
            vec![(self.output.as_path(), "--output")];
        if let Some(rejects) = stage_opts.correct.as_ref().and_then(|c| c.rejects_path.as_ref()) {
            write_targets.push((rejects.as_path(), "--rejects"));
        }
        #[cfg(feature = "consensus")]
        {
            // Only one consensus stage is ever wired, but enumerate all three
            // uniformly; the absent bags contribute nothing.
            for rejects in [
                stage_opts.simplex.as_ref().and_then(|o| o.rejects_opts.rejects.as_ref()),
                stage_opts.duplex.as_ref().and_then(|o| o.rejects_opts.rejects.as_ref()),
                stage_opts.codec.as_ref().and_then(|o| o.rejects_opts.rejects.as_ref()),
            ]
            .into_iter()
            .flatten()
            {
                write_targets.push((rejects.as_path(), "--rejects"));
            }
            for stats in [
                stage_opts.simplex.as_ref().and_then(|o| o.stats_opts.stats.as_ref()),
                stage_opts.duplex.as_ref().and_then(|o| o.stats_opts.stats.as_ref()),
                stage_opts.codec.as_ref().and_then(|o| o.stats_opts.stats.as_ref()),
            ]
            .into_iter()
            .flatten()
            {
                write_targets.push((stats.as_path(), "--stats"));
            }
        }
        if let Some(filter) = stage_opts.filter.as_ref() {
            if let Some(rejects) = filter.rejects.as_ref() {
                write_targets.push((rejects.as_path(), "--filter::rejects"));
            }
            if let Some(stats) = filter.stats.as_ref() {
                write_targets.push((stats.as_path(), "--filter::stats"));
            }
        }
        crate::commands::common::reject_output_collisions(&write_targets)?;

        let spec = ChainSpec {
            stages,
            source,
            sink,
            stage_opts,
            threading: self.threading.clone(),
            compression: self.compression.clone(),
            scheduler: self.scheduler_opts.clone(),
            queue_memory: self.queue_memory.clone(),
            async_reader: self.async_reader,
            read_streams: fgumi_bam_io::ReadStreams::Fixed(1),
            verify_crc,
            command_line: command_line.to_string(),
        };

        crate::pipeline::chains::build_for(spec)?.run()?;

        timer.log_completion(0);
        log::info!("runall completed successfully");
        Ok(())
    }
}

#[cfg(test)]
mod execute_tests {
    use super::*;
    use clap::Parser;

    fn run(args: &[&str]) -> anyhow::Result<()> {
        RunAll::try_parse_from(std::iter::once("runall").chain(args.iter().copied()))
            .unwrap()
            .execute("runall test")
    }

    #[test]
    fn rejects_ref_without_methylation_on_nonalign_chain() {
        let e = run(&[
            "--start-from",
            "group",
            "--stop-after",
            "group",
            "-i",
            "in.bam",
            "-o",
            "o.bam",
            "--group::strategy",
            "adjacency",
            "--ref",
            "ref.fa",
        ])
        .unwrap_err()
        .to_string();
        assert!(e.contains("--ref requires --methylation-mode"), "got: {e}");
    }

    #[test]
    fn rejects_stop_after_align_at_clap_layer() {
        // StopAfter has no `align` variant → clap rejects it before execute().
        assert!(
            RunAll::try_parse_from([
                "runall",
                "--start-from",
                "align",
                "--stop-after",
                "align",
                "-i",
                "i.bam",
                "-o",
                "o.bam"
            ])
            .is_err()
        );
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
mod bag_tests {
    use super::*;
    use crate::assigner::Strategy;
    use crate::commands::sort::SortOrderArg;
    use crate::pipeline::chains::Stage;
    use clap::Parser;

    fn parse(args: &[&str]) -> RunAll {
        RunAll::try_parse_from(std::iter::once("runall").chain(args.iter().copied())).unwrap()
    }

    #[test]
    fn sort_order_is_forced_template_coordinate() {
        let r = parse(&[
            "--start-from",
            "sort",
            "--stop-after",
            "group",
            "-i",
            "in.bam",
            "-o",
            "out.bam",
            "--group::strategy",
            "adjacency",
            "--sort::order",
            "coordinate",
        ]);
        let bag = r.build_stage_options_bag(&[Stage::Sort, Stage::Group]).unwrap();
        assert_eq!(bag.sort.unwrap().order, SortOrderArg::TemplateCoordinate);
    }

    #[test]
    fn group_effective_strategy_resolved() {
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
            "--group::edits",
            "1",
        ]);
        let bag = r.build_stage_options_bag(&[Stage::Group]).unwrap();
        let g = bag.group.unwrap();
        assert_eq!(g.effective_strategy, Strategy::Adjacency);
        assert_eq!(g.effective_edits, 1);
    }

    #[test]
    fn fused_sort_max_memory_auto_is_overridden_to_fixed() {
        let r = parse(&[
            "--start-from",
            "sort",
            "--stop-after",
            "group",
            "-i",
            "in.bam",
            "-o",
            "out.bam",
            "--group::strategy",
            "adjacency",
            "--sort::max-memory",
            "auto",
        ]);
        let bag = r.build_stage_options_bag(&[Stage::Sort, Stage::Group]).unwrap();
        assert!(matches!(
            bag.sort.unwrap().max_memory,
            crate::commands::common::MemoryLimit::Fixed(_)
        ));
    }

    #[test]
    fn sort_rejects_zstd_with_zero_temp_compression() {
        // The inherent `SortOptions::validate` cross-field check must run in the
        // runall sort branch, rejecting the invalid zstd/level-0 pairing before
        // the chain builder wires the spill compressor (rather than failing
        // lazily on the first spill).
        let r = parse(&[
            "--start-from",
            "sort",
            "--stop-after",
            "group",
            "-i",
            "in.bam",
            "-o",
            "out.bam",
            "--group::strategy",
            "adjacency",
            "--sort::temp-codec",
            "zstd",
            "--sort::temp-compression",
            "0",
        ]);
        let Err(err) = r.build_stage_options_bag(&[Stage::Sort, Stage::Group]) else {
            panic!("expected --sort::temp-codec zstd + --sort::temp-compression 0 to be rejected");
        };
        let err = err.to_string();
        assert!(
            err.contains("--temp-compression 0 is only supported with --temp-codec bgzf"),
            "got: {err}"
        );
    }

    #[cfg(feature = "consensus")]
    #[test]
    fn codec_rejects_methylation_mode() {
        let r = parse(&[
            "--start-from",
            "group",
            "--stop-after",
            "consensus",
            "-i",
            "in.bam",
            "-o",
            "out.bam",
            "--group::strategy",
            "adjacency",
            "--consensus",
            "codec",
            "--codec::min-reads",
            "1",
            "--methylation-mode",
            "em-seq",
        ]);
        let Err(err) = r.build_stage_options_bag(&[Stage::Codec]) else {
            panic!("expected --methylation-mode + codec to be rejected");
        };
        let err = err.to_string();
        assert!(err.contains("not supported with codec consensus"), "got: {err}");
    }

    #[cfg(not(feature = "consensus"))]
    #[test]
    fn consensus_stages_rejected_without_consensus_feature() {
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
        for stage in [Stage::Simplex, Stage::Duplex, Stage::Codec] {
            let Err(err) = r.build_stage_options_bag(&[stage]) else {
                panic!("expected {stage:?} to be rejected without the consensus feature");
            };
            let err = err.to_string();
            assert!(err.contains("requires building fgumi with the `consensus` feature"));
        }
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

    /// The reverse of `stopafter_maps_to_runallstage_minus_align`: every
    /// `RunAllStage` variant EXCEPT `AlignAndMerge` has a corresponding
    /// `StopAfter` variant that round-trips back to it via `RunAllStage::from`.
    ///
    /// Enumerates all `RunAllStage` variants in a fixed array and maps each to
    /// its `StopAfter` counterpart through an EXHAUSTIVE match (no wildcard
    /// arm) — so if a new `RunAllStage` variant is ever added without also
    /// adding a matching `StopAfter` variant (and updating this match), this
    /// test fails to compile rather than silently passing.
    #[test]
    fn every_runallstage_except_align_has_a_stopafter() {
        let variants = [
            RunAllStage::Extract,
            RunAllStage::Correct,
            RunAllStage::AlignAndMerge,
            RunAllStage::Zipper,
            RunAllStage::Sort,
            RunAllStage::Group,
            RunAllStage::Consensus,
            RunAllStage::Filter,
        ];
        for stage in variants {
            let stop_after: Option<StopAfter> = match stage {
                RunAllStage::Extract => Some(StopAfter::Extract),
                RunAllStage::Correct => Some(StopAfter::Correct),
                RunAllStage::AlignAndMerge => None,
                RunAllStage::Zipper => Some(StopAfter::Zipper),
                RunAllStage::Sort => Some(StopAfter::Sort),
                RunAllStage::Group => Some(StopAfter::Group),
                RunAllStage::Consensus => Some(StopAfter::Consensus),
                RunAllStage::Filter => Some(StopAfter::Filter),
            };
            match stop_after {
                Some(stop) => {
                    assert_eq!(
                        RunAllStage::from(stop),
                        stage,
                        "{stage} should round-trip through its StopAfter counterpart"
                    );
                }
                None => {
                    assert_eq!(
                        stage,
                        RunAllStage::AlignAndMerge,
                        "only AlignAndMerge should have no StopAfter counterpart"
                    );
                }
            }
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
