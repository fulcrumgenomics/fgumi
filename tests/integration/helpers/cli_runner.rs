//! CLI-invocation helper for the runall parity tests.
//!
//! Centralises how each runall combination — and each equivalent
//! standalone / staged sequence — is invoked, so the parity tests can
//! treat both sides of the comparison uniformly:
//!
//!   * Class A (single-stage parity): one runall invocation vs one
//!     standalone command (`fgumi runall --start-from X --stop-after X`
//!     vs `fgumi X`).
//!   * Class B (composition / multi-stage parity): one runall
//!     invocation vs a sequence of standalone commands chained through
//!     intermediate BAMs (`fgumi runall --start-from S --stop-after T`
//!     vs `fgumi S → tmp.bam → ... → fgumi T`).
//!
//! Today's runall only accepts `{Sort, Group} → {Simplex, Duplex, Codec}`
//! (6 combos). The remaining 6 combos (5 Class A diagonal + Sort→Group)
//! are gated behind #33 tasks 5-9 and the corresponding tests are
//! `#[ignore = "TDD: blocked on task N of #33"]` — they compile, build
//! the right CLI invocation via this module, and turn green
//! automatically when the gate is lifted.

#![allow(dead_code)]

use std::ffi::OsString;
use std::path::{Path, PathBuf};
use std::process::{Command, Output};

/// Pipeline stages exposed by `--start-from` / `--stop-after`,
/// intentionally a subset of `RunAllStage`: the `run_runall`
/// helper below only supports BAM-in/BAM-out stages with no
/// extra-input requirements. `RunAllStage::Correct` (needs a
/// UMI whitelist), `AlignAndMerge` (needs a reference + aligner),
/// and `Zipper` (needs an unmapped second BAM) are excluded here
/// and tested via dedicated harnesses.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum Stage {
    AlignAndMerge,
    Zipper,
    Sort,
    Group,
    Simplex,
    Duplex,
    Codec,
}

impl Stage {
    pub fn cli_value(self) -> &'static str {
        match self {
            Self::AlignAndMerge => "align-and-merge",
            Self::Zipper => "zipper",
            Self::Sort => "sort",
            Self::Group => "group",
            Self::Simplex => "simplex",
            Self::Duplex => "duplex",
            Self::Codec => "codec",
        }
    }

    /// `true` for `Simplex`/`Duplex`/`Codec` (the three consensus
    /// algorithms, which `fgumi runall` selects via `--consensus`).
    pub fn is_consensus(self) -> bool {
        matches!(self, Self::Simplex | Self::Duplex | Self::Codec)
    }

    /// The `--start-from`/`--stop-after` value for `fgumi runall`. The
    /// three consensus algorithms collapse to the single `consensus`
    /// stage; the algorithm is passed separately via `--consensus`.
    pub fn runall_stage_value(self) -> &'static str {
        if self.is_consensus() { "consensus" } else { self.cli_value() }
    }

    /// The `--consensus` mode value when this stage is a consensus
    /// algorithm, else `None`.
    pub fn consensus_mode_value(self) -> Option<&'static str> {
        match self {
            Self::Simplex => Some("simplex"),
            Self::Duplex => Some("duplex"),
            Self::Codec => Some("codec"),
            _ => None,
        }
    }
}

/// Zipper-specific input triple. The standalone `fgumi zipper` and
/// `fgumi runall --start-from zipper` both consume a mapped BAM,
/// an unmapped BAM, and a reference FASTA (with `.dict`). The Class A
/// / Class B parity helpers below use this struct to bind the three
/// paths together; non-zipper stages still use the single-input
/// `run_runall` / `run_standalone` helpers.
#[derive(Debug, Clone)]
pub struct ZipperInputs {
    pub mapped: PathBuf,
    pub unmapped: PathBuf,
    pub reference: PathBuf,
}

/// Strategy / consensus knobs shared by every invocation in a parity
/// pair. Both sides of the parity comparison use the same values so
/// the test is apples-to-apples.
#[derive(Debug, Clone)]
pub struct ParityArgs {
    pub strategy: &'static str,
    pub edits: u32,
    pub min_reads: &'static str,
    pub threads: u32,
}

impl ParityArgs {
    /// Default args for simplex / codec single-strand flows.
    pub fn for_simplex_like() -> Self {
        Self { strategy: "adjacency", edits: 1, min_reads: "1", threads: 1 }
    }

    /// Default args for the duplex two-strand flow.
    pub fn for_duplex() -> Self {
        Self { strategy: "paired", edits: 1, min_reads: "1,1,0", threads: 1 }
    }
}

pub fn fgumi_binary() -> PathBuf {
    PathBuf::from(env!("CARGO_BIN_EXE_fgumi"))
}

pub fn fgumi(args: &[OsString]) -> Output {
    Command::new(fgumi_binary()).args(args).output().expect("failed to execute fgumi")
}

/// Run `fgumi runall --start-from <start> --stop-after <stop>` over
/// a single-input pipeline. `Stage::Zipper` is rejected because it
/// needs the (mapped, unmapped, reference) triple — use
/// `run_runall_zipper` for that.
pub fn run_runall(
    start: Stage,
    stop: Stage,
    input: &Path,
    output: &Path,
    args: &ParityArgs,
) -> Output {
    assert!(
        !matches!(start, Stage::Zipper | Stage::AlignAndMerge)
            && !matches!(stop, Stage::Zipper | Stage::AlignAndMerge),
        "run_runall does not support Stage::Zipper / Stage::AlignAndMerge — \
         use the dedicated runner (`run_runall_zipper` for zipper; AAM has no \
         runner yet — lands in #33 C4)"
    );
    // After the riker-style port (#33), per-stage tuning flags live
    // under `--<stage>::<flag>` on `fgumi runall`. Strategy and edits
    // are group-stage knobs, exposed as `--group::strategy` and
    // `--group::edits`. Per-mode consensus tuning (min-reads, error
    // rates, …) lives under `--simplex::*` / `--duplex::*` / `--codec::*`
    // — there is no shared top-level `--min-reads`. `--threads` stays
    // top-level.
    let mut cli_args: Vec<OsString> = vec![
        "runall".into(),
        "--start-from".into(),
        start.runall_stage_value().into(),
        "--stop-after".into(),
        stop.runall_stage_value().into(),
        "--input".into(),
        input.as_os_str().to_owned(),
        "--output".into(),
        output.as_os_str().to_owned(),
        "--group::strategy".into(),
        args.strategy.into(),
        "--group::edits".into(),
        args.edits.to_string().into(),
        "--threads".into(),
        args.threads.to_string().into(),
    ];
    // The consensus algorithm (simplex/duplex/codec) is selected with
    // `--consensus` when the chain reaches the consensus stage; min-reads
    // is then passed via that mode's `--<mode>::min-reads` flag.
    if let Some(mode) = stop.consensus_mode_value().or_else(|| start.consensus_mode_value()) {
        cli_args.push("--consensus".into());
        cli_args.push(mode.into());
        cli_args.push(format!("--{mode}::min-reads").into());
        cli_args.push(args.min_reads.into());
    }
    fgumi(&cli_args)
}

/// Run the fused `fgumi runall --start-from <start> --consensus <mode>
/// --stop-after filter` pipeline (consensus → filter in one pass, no
/// intermediate consensus BAM). `consensus` must be a consensus stage
/// (`Simplex`/`Duplex`/`Codec`); `start` is the pre-consensus start stage
/// (e.g. `Group`). Filter uses its default `--filter::*` parameters.
pub fn run_runall_consensus_to_filter(
    start: Stage,
    consensus: Stage,
    input: &Path,
    output: &Path,
    args: &ParityArgs,
) -> Output {
    let mode = consensus
        .consensus_mode_value()
        .expect("run_runall_consensus_to_filter: `consensus` must be a consensus stage");
    let cli_args: Vec<OsString> = vec![
        "runall".into(),
        "--start-from".into(),
        start.runall_stage_value().into(),
        "--stop-after".into(),
        "filter".into(),
        "--consensus".into(),
        mode.into(),
        "--input".into(),
        input.as_os_str().to_owned(),
        "--output".into(),
        output.as_os_str().to_owned(),
        "--group::strategy".into(),
        args.strategy.into(),
        "--group::edits".into(),
        args.edits.to_string().into(),
        format!("--{mode}::min-reads").into(),
        args.min_reads.into(),
        // `fgumi filter` has no default for `--min-reads` (it is required);
        // pass the same value the staged `run_standalone_filter` uses so the
        // two filter passes are apples-to-apples.
        "--filter::min-reads".into(),
        FILTER_MIN_READS.into(),
        "--threads".into(),
        args.threads.to_string().into(),
    ];
    fgumi(&cli_args)
}

/// Filter `--min-reads` value used by both sides of the consensus → filter
/// parity comparison. `fgumi filter` requires this flag (no clap default), so
/// the fused and staged invocations must pass the same value.
const FILTER_MIN_READS: &str = "1";

/// Run standalone `fgumi filter --input <input> --output <output>` with the
/// shared thread count and the required `--min-reads`. The staged counterpart
/// to the filter stage in [`run_runall_consensus_to_filter`].
pub fn run_standalone_filter(input: &Path, output: &Path, args: &ParityArgs) -> Output {
    let cli_args: Vec<OsString> = vec![
        "filter".into(),
        "--input".into(),
        input.as_os_str().to_owned(),
        "--output".into(),
        output.as_os_str().to_owned(),
        "--min-reads".into(),
        FILTER_MIN_READS.into(),
        "--threads".into(),
        args.threads.to_string().into(),
    ];
    fgumi(&cli_args)
}

/// Run a single standalone command for Class A parity. Each arm maps
/// `stage` to the equivalent `fgumi <stage>` invocation with the
/// command-specific flags every `ParityArgs` field controls.
///
/// Argument mapping rationale (verified against
/// `src/commands/<stage>.rs` clap-derive signatures):
///   * `sort` — input/output/threads only; the parity fixture is
///     template-coordinate, so use `--order template-coordinate`.
///   * `group` — strategy/edits/threads.
///   * `simplex`/`duplex`/`codec` — min-reads/threads.
pub fn run_standalone(stage: Stage, input: &Path, output: &Path, args: &ParityArgs) -> Output {
    let mut cli_args: Vec<OsString> = vec![
        stage.cli_value().into(),
        "--input".into(),
        input.as_os_str().to_owned(),
        "--output".into(),
        output.as_os_str().to_owned(),
        "--threads".into(),
        args.threads.to_string().into(),
    ];
    match stage {
        Stage::Sort => {
            cli_args.push("--order".into());
            cli_args.push("template-coordinate".into());
        }
        Stage::Group => {
            cli_args.push("--strategy".into());
            cli_args.push(args.strategy.into());
            cli_args.push("--edits".into());
            cli_args.push(args.edits.to_string().into());
        }
        Stage::Simplex | Stage::Duplex | Stage::Codec => {
            cli_args.push("--min-reads".into());
            cli_args.push(args.min_reads.into());
        }
        Stage::Zipper => panic!(
            "run_standalone does not support Stage::Zipper — zipper has \
             two BAM inputs + a reference. Use `run_standalone_zipper` instead."
        ),
        Stage::AlignAndMerge => panic!(
            "run_standalone does not support Stage::AlignAndMerge — there is no \
             standalone `fgumi align-and-merge` command. AAM is only reachable via \
             `fgumi runall --start-from align-and-merge`."
        ),
    }
    fgumi(&cli_args)
}

/// Run a multi-stage chain through standalone commands via intermediate
/// temp BAMs. Produces the same final output that a fused
/// `runall --start-from S --stop-after T` would produce.
///
/// `chain` lists the standalone commands in execution order. The first
/// command reads from `input`; intermediate commands chain through
/// temp BAMs in `scratch_dir`; the last command writes to `final_output`.
///
/// Returns the `Output` of the LAST command in the chain. If any
/// earlier command fails, this function panics with that command's
/// stderr (chain failures should surface immediately, not silently
/// poison `final_output`).
pub fn run_staged_chain(
    chain: &[Stage],
    input: &Path,
    final_output: &Path,
    scratch_dir: &Path,
    args: &ParityArgs,
) -> Output {
    assert!(!chain.is_empty(), "run_staged_chain: chain must be non-empty");
    let mut current_input: PathBuf = input.to_path_buf();
    for (i, stage) in chain.iter().enumerate() {
        let is_last = i + 1 == chain.len();
        let stage_output = if is_last {
            final_output.to_path_buf()
        } else {
            scratch_dir.join(format!("staged_{}_{}.bam", i, stage.cli_value()))
        };
        let output = run_standalone(*stage, &current_input, &stage_output, args);
        if is_last {
            return output;
        }
        assert!(
            output.status.success(),
            "staged_chain step {i} ({}) failed: stderr={}",
            stage.cli_value(),
            String::from_utf8_lossy(&output.stderr)
        );
        current_input = stage_output;
    }
    unreachable!("loop returns on the last iteration")
}

/// Run `fgumi runall --start-from zipper --stop-after <stop>` over a
/// zipper triple (mapped + unmapped + reference). When
/// `stop != Stage::Zipper`, runall chains zipper through the
/// downstream stages internally (via tempfiles, per
/// `execute_zipper_then_downstream`).
pub fn run_runall_zipper(
    stop: Stage,
    inputs: &ZipperInputs,
    output: &Path,
    args: &ParityArgs,
) -> Output {
    let mut cli_args: Vec<OsString> = vec![
        "runall".into(),
        "--start-from".into(),
        Stage::Zipper.cli_value().into(),
        "--stop-after".into(),
        stop.runall_stage_value().into(),
        "--input".into(),
        inputs.mapped.as_os_str().to_owned(),
        "--unmapped".into(),
        inputs.unmapped.as_os_str().to_owned(),
        "--ref".into(),
        inputs.reference.as_os_str().to_owned(),
        "--output".into(),
        output.as_os_str().to_owned(),
        "--threads".into(),
        args.threads.to_string().into(),
    ];
    // Downstream group stage needs the group knobs. They're ignored when
    // stop = Zipper but harmless.
    if !matches!(stop, Stage::Zipper) {
        cli_args.extend([
            "--group::strategy".into(),
            args.strategy.into(),
            "--group::edits".into(),
            args.edits.to_string().into(),
        ]);
    }
    // Select the consensus algorithm when the chain reaches consensus;
    // min-reads is passed via that mode's `--<mode>::min-reads` flag.
    if let Some(mode) = stop.consensus_mode_value() {
        cli_args.push("--consensus".into());
        cli_args.push(mode.into());
        cli_args.push(format!("--{mode}::min-reads").into());
        cli_args.push(args.min_reads.into());
    }
    fgumi(&cli_args)
}

/// Run the standalone `fgumi zipper` command against a zipper triple.
pub fn run_standalone_zipper(inputs: &ZipperInputs, output: &Path, args: &ParityArgs) -> Output {
    let cli_args: Vec<OsString> = vec![
        Stage::Zipper.cli_value().into(),
        "--input".into(),
        inputs.mapped.as_os_str().to_owned(),
        "--unmapped".into(),
        inputs.unmapped.as_os_str().to_owned(),
        "--reference".into(),
        inputs.reference.as_os_str().to_owned(),
        "--output".into(),
        output.as_os_str().to_owned(),
        "--threads".into(),
        args.threads.to_string().into(),
    ];
    fgumi(&cli_args)
}
