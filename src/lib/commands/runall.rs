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
