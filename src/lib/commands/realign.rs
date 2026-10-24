//! `realign` — re-align consensus reads, preferring the alignment inside each
//! molecule's raw-read window.
//!
//! Consensus callers emit unmapped reads. `fgumi group --realign-window` records
//! each molecule's raw-read window (`rw`), the consensus callers copy it onto the
//! consensus reads, and this command aligns them with the same backends as
//! `fgumi runall`'s align stage, zipper-merges the alignments onto the consensus
//! reads, then chooses each end's primary with the window: bwa's primary if it
//! lies in the window, else the first in-window `XA` alternative, else an
//! in-window supplementary (split) piece of the same read, else bwa's primary
//! unchanged. `rw` is removed from the output.

use std::path::{Path, PathBuf};
use std::str::FromStr;

use anyhow::{Result, bail};
use clap::Parser;

use crate::aligner::MultiAlignerOptions;
use crate::commands::command::Command;
use crate::commands::common::{
    BamIoOptions, CompressionOptions, MethylationModeArg, QueueMemoryOptions, SchedulerOptions,
    ThreadingOptions, reject_output_collisions, resolve_methylation_mode,
};
use crate::pipeline::chains::options_bag::AlignOptions;

/// bwa's `-h <primary>,<alt>`: the most alternative hits reported in `XA` (`alt`
/// applies to hits on ALT contigs). Only `XA` hits can be promoted to an
/// in-window primary, so widening this widens the candidate set.
///
/// Parse errors do not name the flag: clap prefixes them with the flag the user
/// passed (`--xa-max-hits` or `runall`'s `--realign::xa-max-hits`).
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct XaMaxHits {
    /// Maximum `XA` hits for a read with primary hits (`-h` first value).
    pub primary: u32,
    /// Maximum `XA` hits on ALT contigs (`-h` second value).
    pub alt: u32,
}

impl Default for XaMaxHits {
    /// bwa-mem3's default, `-h 5,200`.
    fn default() -> Self {
        Self { primary: 5, alt: 200 }
    }
}

impl FromStr for XaMaxHits {
    type Err = anyhow::Error;

    /// Parses `<primary>,<alt>`, or a single `<n>` (bwa sets both to `n`).
    fn from_str(s: &str) -> Result<Self> {
        // The flags take hyphenated words as their value, so that `-1,5` given
        // as its own word reaches the checks below; a word that looks like an
        // option instead means the value was left out.
        if looks_like_an_option(s) {
            bail!("a value is required, but '{s}' is an option; expected <primary>,<alt> or <n>");
        }
        let parse = |v: &str| -> Result<u32> {
            let n: u32 = v.trim().parse().map_err(|e| {
                anyhow::anyhow!(
                    "'{v}' is not a non-negative integer ({e}); expected <primary>,<alt> or <n>"
                )
            })?;
            if i32::try_from(n).is_err() {
                bail!("{n} exceeds {}", i32::MAX);
            }
            Ok(n)
        };
        if let Some((primary, alt)) = s.split_once(',') {
            Ok(Self { primary: parse(primary)?, alt: parse(alt)? })
        } else {
            let n = parse(s)?;
            Ok(Self { primary: n, alt: n })
        }
    }
}

/// Whether `word` reads as an option rather than an [`XaMaxHits`] value: it
/// starts with `-` and no digit follows (so `-1,5` is a value, and is rejected
/// as negative).
fn looks_like_an_option(word: &str) -> bool {
    word.strip_prefix('-').is_some_and(|rest| !rest.starts_with(|c: char| c.is_ascii_digit()))
}

/// The flags that take an [`XaMaxHits`], with the subcommand each belongs to.
const XA_MAX_HITS_FLAGS: [(&str, &str); 2] =
    [("--xa-max-hits", "realign"), ("--realign::xa-max-hits", "runall")];

/// The error for an `XaMaxHits` flag left without a value, as `(subcommand,
/// message)`, or `None` when every one in `args` has a value.
///
/// The flags take hyphenated words as their value (so `-1,5` reaches the
/// non-negative check), so clap takes the next option as the value. When that
/// option takes a value of its own, clap reports the stray value as an
/// unexpected argument before the hits value is ever checked, and that error
/// does not name the flag. Called on a failed parse, this names it instead: a
/// flag that is the last word, or is followed by a word that reads as an
/// option, is missing its value. Only a flag after its own subcommand's name
/// counts, so another subcommand's error is left as clap reports it.
#[must_use]
pub fn missing_xa_max_hits_value<S: AsRef<std::ffi::OsStr>>(
    args: &[S],
) -> Option<(&'static str, String)> {
    args.iter().enumerate().find_map(|(i, arg)| {
        let (flag, subcommand) = XA_MAX_HITS_FLAGS
            .iter()
            .find(|(flag, _)| arg.as_ref() == std::ffi::OsStr::new(flag))?;
        if !args[..i].iter().any(|word| word.as_ref() == std::ffi::OsStr::new(subcommand)) {
            return None;
        }
        let next = args.get(i + 1).map(|next| next.as_ref().to_string_lossy());
        let message = match next {
            None => format!(
                "a value is required for '{flag} <PRIMARY[,ALT]>' but none was supplied; \
                 expected <primary>,<alt> or <n>"
            ),
            Some(next) if looks_like_an_option(&next) => format!(
                "a value is required for '{flag} <PRIMARY[,ALT]>' but none was supplied \
                 ('{next}' is an option); expected <primary>,<alt> or <n>"
            ),
            Some(_) => return None,
        };
        Some((*subcommand, message))
    })
}

impl std::fmt::Display for XaMaxHits {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        write!(f, "{},{}", self.primary, self.alt)
    }
}

/// Re-align consensus reads and prefer the alignment inside each molecule's raw-read window.
#[derive(Debug, Parser)]
#[command(
    name = "realign",
    about = "\x1b[38;5;173m[POST-CONSENSUS]\x1b[0m \x1b[36mRe-align consensus reads, preferring each molecule's raw-read window\x1b[0m",
    long_about = r#"
Re-align consensus reads and choose each end's primary alignment inside the molecule's
raw-read window.

The input is the (optionally filtered) output of `fgumi simplex`, `fgumi duplex` or
`fgumi codec`, called from a BAM grouped with `fgumi group --realign-window`. That flag
records each molecule's raw-read window in the `rw` tag, which the consensus callers copy
onto the consensus reads. Only mapped raw reads give a window, so a consensus end whose
raw reads never mapped has no `rw`: it is aligned, left as bwa placed it, and counted as
without a window. An input that was never window-grouped is an error naming its first read:
realign fails if none of the first 1000 consensus templates (or, in a shorter input, none at
all) carries `rw`, and it decides this before any read is aligned or written, so a rejected
input leaves no records in the output.

Reads are aligned with the same backends as `fgumi runall`'s align stage
(`--aligner::preset bwa-mem3|bwa|bwa-mem3-inproc` or `--aligner::command`) and merged
onto the consensus reads as `fgumi zipper --tags-to-reverse Consensus
--tags-to-revcomp Consensus` would: the per-base consensus tags of a reverse-strand
read are reversed (and its per-base bases reverse-complemented) here, so do not pass
`fgumi filter --reverse-per-base-tags` on realign's output. For each end:

- bwa's primary overlaps the window (same contig and strand): kept as is.
- else the first `XA` hit that does is promoted: it becomes the primary with MAPQ 0,
  NM from the `XA` hit, MD regenerated from the reference, AS from the hit (below),
  the old primary moved into `XA`, `XS` set to the old primary's score, and `SA` and
  that end's split records dropped: its supplementary records, and the secondary
  records carrying `SA` that an aligner run with `-M` (as `--meth` is) writes for
  the pieces of a split. Secondary records without `SA` (`-a`) are kept.
- else, if a split piece of the read (a supplementary record, or under `-M` a
  secondary one carrying `SA`) overlaps the window, it swaps roles with the primary;
  of several, the one with the longest reference span, then the first. The piece
  becomes the primary with the full SEQ/QUAL in its strand's orientation (hard clips
  turned into soft clips) and keeps its own MAPQ, NM, MD, AS and XS as bwa computed
  them, without `XA`. The old primary becomes the piece, flagged and clipped (hard or
  soft) as the piece was. `SA` is rewritten on every piece of the read, and its other
  pieces stay.
- else bwa's primary is kept and counted as outside the window.

The proper-pair guard: neither swap is made when bwa flagged the pair proper and the
mate's primary lies in the mate's own window. The mate then agrees with bwa's
placement, so this end's window is the suspect (for example, raw reads placed on a
poly-G artefact that the consensus read no longer prefers); the end keeps bwa's
primary and is counted separately ("Primary kept outside window (proper pair with an
in-window mate)" in the summary). Both
ends are judged as bwa left them, before either is swapped. Fragments, and pairs whose
mate is unmapped, outside its window or without one, are not guarded.

After a swap the mate fields are recomputed (secondary and supplementary records
included, with their `tc` tags). The proper-pair flag is never set by a swap: it is
kept only if bwa had set it and the pair is still FR on one contig with an insert no
larger than before, and cleared otherwise. The windows say where each end belongs,
not that the molecule's raw reads were a concordant pair.

A consensus end whose mate is absent (removed by a per-record `fgumi filter
--filter-by-template false`) is sent to the aligner with a one-base placeholder mate,
so the aligner's batches (`-K`, cut only after an even number of reads) still hold
whole pairs; the placeholder is removed again. The end keeps its consensus tags and
its R1/R2 flag, with its mate marked unmapped and placed at the read itself.

`--xa-max-hits` sets bwa's `-h`: only hits listed in `XA` can be promoted, so widening
it widens the candidate set. bwa omits `XA` entirely when more than `-h` hits score
near the best, so the most repetitive reads have no candidates at all. Unset, the
aligner keeps its own default (bwa's `5,200`) and an `--aligner::command` is passed
through unchanged.

A malformed `XA` never fails the run: an `XA` tag that does not parse is treated as
absent and a hit whose CIGAR is invalid is skipped; both are counted in the summary.
A malformed `rw` does fail it, naming the read: `rw` is fgumi's own tag, so a value
that does not parse means the input was damaged or written by another tool.

A promoted alignment's AS is the score bwa lists for the hit in a six-field `XA`
(`bwa-mem3 mem -u`). With the usual four-field `XA` it is rescored from its CIGAR and
the `XA` hit's edit distance
(so bisulfite conversions bwa did not count are not penalised), with the aligner's
scoring, and never below 0 (bwa reports no negative AS): the in-process backend's
options, or for a subprocess preset or `--aligner::command` what the command passes
after its `mem` word: `-A`, `-B`, `-O`, `-E`, an `-x` preset, and `--meth` (whose
collapsed scoring, the default except for `--meth=taps`, sets `-B 2`; `--meth-scoring
genomic|neutral` keeps bwa's `-B`). Options it does not pass keep bwa-mem3's defaults.
An `--aligner::command` is read for these options (and for `-h`, above) by splitting
it on whitespace up to the first `|`, `;`, `&&` or `||`: shell quoting is not
interpreted, so an option-like word inside a quoted argument (say, in `-R '...'`) is
read as an option.

The `rw` tag is removed. Output is in input order (unsorted); pipe it to
`fgumi sort --order coordinate` for a coordinate-sorted BAM.
"#
)]
pub struct Realign {
    /// Input/output BAM paths and reader options.
    #[command(flatten)]
    pub io: BamIoOptions,

    /// Reference FASTA the aligner index was built from (needs a `.dict` alongside).
    #[arg(short = 'r', long = "ref")]
    pub reference: PathBuf,

    /// Aligner selection and tuning, exposed as `--aligner::*` exactly as on `fgumi runall`.
    #[command(flatten)]
    pub aligner_opts: MultiAlignerOptions,

    /// Override path for the aligner binary (preset mode only).
    #[arg(long = "aligner-bin")]
    pub aligner_bin: Option<PathBuf>,

    /// Methylation chemistry. With `em-seq`, the bwa-mem3 presets align
    /// bisulfite-aware (`bwa-mem3 mem --meth`, against a `bwa-mem3 index --meth` index),
    /// which suits simplex consensus reads. Do not pass it for duplex consensus
    /// reads: their sequence is the molecule's own (methylation is in `MM`/`ML`),
    /// so they align with a standard aligner.
    #[arg(long = "methylation-mode", value_enum)]
    pub methylation_mode: Option<MethylationModeArg>,

    /// bwa's `-h <primary>,<alt>`: the most alternative hits reported in `XA`
    /// (a single `<n>` sets both; `0` reports none, so nothing can be promoted).
    /// bwa omits `XA` entirely when more than this many hits score near the best,
    /// so the most repetitive reads get no candidates. Unset, the aligner keeps
    /// its default (bwa's `5,200`) and `--aligner::command` is left untouched.
    /// Set, `-h` goes right after the command's first `mem` word, or at its end
    /// (with a warning) when it has none; a command that already passes `-h`
    /// keeps its own (with a warning).
    #[arg(long = "xa-max-hits", value_name = "PRIMARY[,ALT]", allow_hyphen_values = true)]
    pub xa_max_hits: Option<XaMaxHits>,

    /// Hold the whole reference FASTA in memory instead of reading it on demand.
    /// Realign reads the reference only to rescore the alignments it promotes,
    /// which are few, so by default it reads just those regions through the FASTA's `.fai` index instead of holding the
    /// whole FASTA (on hg38 the default measured about 3.3 GB lower peak memory,
    /// with identical output). Each region is a small disk read, slower when the
    /// FASTA is not in the page cache. A compressed (gzip or BGZF) FASTA, one
    /// without a `.fai` (`samtools faidx`), or one that is not a regular file (a
    /// pipe, such as `<(cat ref.fa)`), cannot be read on demand: realign then
    /// loads it into memory (decompressing a compressed one) with a warning. Set
    /// this to load it into memory regardless. Output is identical either way.
    #[arg(long = "in-memory-reference", value_name = "true|false", default_value = "false", num_args = 0..=1, default_missing_value = "true", action = clap::ArgAction::Set, value_parser = clap::builder::BoolishValueParser::new(), hide_possible_values = true)]
    pub in_memory_reference: bool,

    /// Threading options for the parallel pipeline.
    #[command(flatten)]
    pub threading: ThreadingOptions,

    /// Compression options for the output BAM.
    #[command(flatten)]
    pub compression: CompressionOptions,

    /// Scheduler and pipeline statistics options.
    #[command(flatten)]
    pub scheduler_opts: SchedulerOptions,

    /// Queue memory options.
    #[command(flatten)]
    pub queue_memory: QueueMemoryOptions,
}

/// Options for [`Stage::Realign`](crate::pipeline::chains::Stage::Realign): the
/// aligner (resolved as the align stage resolves it) plus bwa's `-h`.
pub struct RealignOptions {
    /// The aligner, reference and binary override.
    pub align: AlignOptions,
    /// bwa's `-h`, applied only to the realign stage's aligner, and only when the
    /// user set it (`None` leaves the aligner's default).
    pub xa_max_hits: Option<XaMaxHits>,
    /// The flag `xa_max_hits` came from, for messages (`--xa-max-hits` or
    /// `--realign::xa-max-hits`).
    pub xa_max_hits_option: &'static str,
    /// Hold the whole reference in memory rather than reading regions on demand
    /// through its `.fai` index (`--in-memory-reference` /
    /// `--realign::in-memory-reference`).
    pub in_memory_reference: bool,
}

impl Realign {
    /// Project the parsed CLI flags into [`RealignOptions`].
    ///
    /// # Errors
    ///
    /// Returns an error if the `--aligner::*` options fail validation.
    pub fn to_realign_options(&self) -> Result<RealignOptions> {
        let mut aligner = self.aligner_opts.clone().validate()?;
        aligner.methylation_mode = resolve_methylation_mode(self.methylation_mode);
        Ok(RealignOptions {
            align: AlignOptions {
                aligner,
                reference: self.reference.clone(),
                aligner_bin: self.aligner_bin.clone(),
            },
            xa_max_hits: self.xa_max_hits,
            xa_max_hits_option: "--xa-max-hits",
            in_memory_reference: self.in_memory_reference,
        })
    }

    /// Run the realign stage on the declarative chain builder, which opens the
    /// source, logs the banner and wires the aligner itself.
    fn execute_chain(&self, command_line: &str) -> Result<()> {
        use crate::pipeline::chains::{
            ChainSpec, SingleStageContext, Stage, StageOptionsBag, build_for,
        };

        self.io.log_effective_check_crc();
        let stage_opts =
            StageOptionsBag { realign: Some(self.to_realign_options()?), ..Default::default() };
        let ctx = SingleStageContext {
            io: &self.io,
            threading: &self.threading,
            compression: &self.compression,
            scheduler: &self.scheduler_opts,
            queue_memory: &self.queue_memory,
            command_line,
        };
        let spec = ChainSpec::single_stage(Stage::Realign, stage_opts, &ctx);
        build_for(spec)?.run()
    }
}

impl Command for Realign {
    fn execute(&self, command_line: &str) -> Result<()> {
        self.io.validate()?;
        let outputs: [(&Path, &str); 1] = [(self.io.output.as_path(), "--output")];
        reject_output_collisions(&outputs)?;
        crate::commands::common::reject_writes_aliasing_inputs(
            &[(self.io.input.as_path(), "--input"), (self.reference.as_path(), "--ref")],
            &outputs,
        )?;
        self.execute_chain(command_line)
    }
}

#[cfg(test)]
mod tests {
    use rstest::rstest;

    use super::*;
    use crate::aligner::{AlignerPreset, DedupReads};
    use fgumi_consensus::MethylationMode;

    fn parse(extra: &[&str]) -> Realign {
        let mut args = vec!["realign", "-i", "in.bam", "-o", "out.bam", "-r", "ref.fa"];
        args.extend_from_slice(extra);
        Realign::try_parse_from(args).expect("parse")
    }

    #[test]
    fn to_realign_options_carries_every_tuning_flag() {
        let cmd = parse(&[
            "--aligner::preset",
            "bwa-mem3-inproc",
            "--aligner::chunk-size",
            "1000",
            "--aligner::sub-batch-templates",
            "64",
            "--aligner::dedup-reads",
            "off",
            "--aligner-bin",
            "/opt/bwa-mem3",
            "--methylation-mode",
            "em-seq",
            "--xa-max-hits",
            "20,300",
        ]);
        let opts = cmd.to_realign_options().expect("options");
        assert_eq!(opts.align.reference, PathBuf::from("ref.fa"));
        assert_eq!(opts.align.aligner_bin, Some(PathBuf::from("/opt/bwa-mem3")));
        assert_eq!(opts.align.aligner.preset, Some(AlignerPreset::BwaMem3InProc));
        assert_eq!(opts.align.aligner.command, None);
        assert_eq!(opts.align.aligner.chunk_size, 1000);
        assert_eq!(opts.align.aligner.sub_batch_templates, Some(64));
        assert_eq!(opts.align.aligner.dedup_reads, Some(DedupReads::Off));
        assert_eq!(opts.align.aligner.methylation_mode, MethylationMode::EmSeq);
        assert_eq!(opts.xa_max_hits, Some(XaMaxHits { primary: 20, alt: 300 }));
        assert_eq!(opts.xa_max_hits_option, "--xa-max-hits");
    }

    #[test]
    fn to_realign_options_carries_defaults() {
        let cmd = parse(&["--aligner::command", "bwa-mem3 mem -p {ref} /dev/stdin"]);
        let opts = cmd.to_realign_options().expect("options");
        let defaults = crate::aligner::AlignerOptions::default();
        assert_eq!(opts.align.reference, PathBuf::from("ref.fa"));
        assert_eq!(opts.align.aligner_bin, None);
        assert_eq!(opts.align.aligner.preset, None);
        assert_eq!(opts.align.aligner.command.as_deref(), Some("bwa-mem3 mem -p {ref} /dev/stdin"));
        assert_eq!(opts.align.aligner.threads, None);
        assert_eq!(opts.align.aligner.chunk_size, defaults.chunk_size);
        assert_eq!(opts.align.aligner.sub_batch_templates, None);
        assert_eq!(opts.align.aligner.dedup_reads, None);
        assert_eq!(opts.align.aligner.methylation_mode, MethylationMode::Disabled);
        assert_eq!(opts.xa_max_hits, None);
        assert!(!opts.in_memory_reference, "the reference is read on demand by default");
    }

    /// `--in-memory-reference` alone opts out of on-demand reads; an explicit
    /// value sets it either way.
    #[rstest]
    #[case::unset(&[], false)]
    #[case::bare(&["--in-memory-reference"], true)]
    #[case::explicit_true(&["--in-memory-reference", "true"], true)]
    #[case::explicit_false(&["--in-memory-reference", "false"], false)]
    fn in_memory_reference_parses(#[case] extra: &[&str], #[case] expected: bool) {
        let opts = parse(extra).to_realign_options().expect("options");
        assert_eq!(opts.in_memory_reference, expected);
    }

    #[rstest]
    #[case::pair("5,200", Some(XaMaxHits { primary: 5, alt: 200 }))]
    #[case::single("12", Some(XaMaxHits { primary: 12, alt: 12 }))]
    #[case::zero("0,0", Some(XaMaxHits { primary: 0, alt: 0 }))]
    #[case::spaces(" 3 , 9 ", Some(XaMaxHits { primary: 3, alt: 9 }))]
    fn xa_max_hits_parses(#[case] input: &str, #[case] expected: Option<XaMaxHits>) {
        assert_eq!(input.parse::<XaMaxHits>().ok(), expected);
    }

    /// Errors describe the value, not a flag name: clap prefixes them with the
    /// flag the user passed, which differs between `realign` and `runall`.
    #[rstest]
    #[case::negative("-1,5", "'-1' is not a non-negative integer")]
    #[case::too_large("2147483648,5", "2147483648 exceeds 2147483647")]
    #[case::not_a_number("five", "'five' is not a non-negative integer")]
    #[case::three_values("1,2,3", "'2,3' is not a non-negative integer")]
    #[case::an_option("--threads", "a value is required, but '--threads' is an option")]
    #[case::a_short_option("-t", "a value is required, but '-t' is an option")]
    fn xa_max_hits_rejects(#[case] input: &str, #[case] message: &str) {
        let err = input.parse::<XaMaxHits>().expect_err("invalid").to_string();
        assert!(err.contains(message), "{err}");
        assert!(!err.contains("--xa-max-hits"), "{err}");
    }

    /// clap names the flag the user passed; the value's own message follows,
    /// also for a negative value given as its own word (clap would otherwise
    /// read `-1,5` as an unknown option), and for a left-out value, where the
    /// next word is an option that takes no value or is the last word.
    #[rstest]
    #[case::not_a_number(&["five"], "'five' is not a non-negative integer")]
    #[case::negative(&["-1,5"], "'-1' is not a non-negative integer")]
    #[case::missing_before_a_flag(
        &["--in-memory-reference"],
        "a value is required, but '--in-memory-reference' is an option"
    )]
    #[case::missing_at_the_end(&["--threads"], "a value is required, but '--threads' is an option")]
    fn xa_max_hits_parse_error_names_the_flag(#[case] value: &[&str], #[case] message: &str) {
        let mut args =
            vec!["realign", "-i", "in.bam", "-o", "out.bam", "-r", "ref.fa", "--xa-max-hits"];
        args.extend_from_slice(value);
        let err = Realign::try_parse_from(args).expect_err("invalid value");
        let rendered = err.to_string();
        assert!(rendered.contains("--xa-max-hits"), "{rendered}");
        assert!(rendered.contains(message), "{rendered}");
    }

    /// A left-out value before an option that takes one never parses: clap
    /// reads the option as the value and then rejects the option's own value
    /// as an unexpected argument (it reports that before the value's own
    /// message), so `--threads` is never silently taken as the hits count.
    /// `fgumi` replaces that error with [`missing_xa_max_hits_value`]'s.
    #[test]
    fn xa_max_hits_missing_before_an_option_with_a_value_fails() {
        let args = [
            "realign",
            "-i",
            "in.bam",
            "-o",
            "out.bam",
            "-r",
            "ref.fa",
            "--xa-max-hits",
            "--threads",
            "4",
        ];
        let err = Realign::try_parse_from(args).expect_err("missing value");
        assert!(err.to_string().contains("unexpected argument '4'"), "{err}");
        let (subcommand, message) = missing_xa_max_hits_value(&args).expect("missing");
        assert_eq!(subcommand, "realign");
        assert!(message.contains("'--xa-max-hits <PRIMARY[,ALT]>'"), "{message}");
        assert!(message.contains("('--threads' is an option)"), "{message}");
    }

    /// An `XaMaxHits` flag is missing its value when it is the last word or the
    /// next word reads as an option; a value, negative or not, is left to the
    /// value's own checks. Only the flag of the subcommand being run, after it,
    /// is considered: another subcommand's error is left to clap.
    #[rstest]
    #[case::realign_before_option(&["realign", "--xa-max-hits", "--threads", "4"], Some("realign"))]
    #[case::realign_last(&["realign", "--xa-max-hits"], Some("realign"))]
    #[case::realign_short_option(&["realign", "--xa-max-hits", "-t", "4"], Some("realign"))]
    #[case::realign_value(&["realign", "--xa-max-hits", "5,200"], None)]
    #[case::realign_negative(&["realign", "--xa-max-hits", "-1,5"], None)]
    #[case::runall_before_option(
        &["runall", "--realign::xa-max-hits", "--threads", "4"],
        Some("runall")
    )]
    #[case::runall_negative(&["runall", "--realign::xa-max-hits", "-1,5"], None)]
    #[case::attached(&["realign", "--xa-max-hits=--threads", "4"], None)]
    #[case::absent(&["realign", "--threads", "4"], None)]
    #[case::other_subcommand(&["group", "--xa-max-hits", "--threads", "4"], None)]
    #[case::other_subcommands_flag(&["runall", "--xa-max-hits", "--threads", "4"], None)]
    #[case::before_the_subcommand(&["--xa-max-hits", "realign", "--threads", "4"], None)]
    fn missing_xa_max_hits_value_finds_a_left_out_value(
        #[case] args: &[&str],
        #[case] expected: Option<&str>,
    ) {
        let got = missing_xa_max_hits_value(args);
        assert_eq!(got.as_ref().map(|(subcommand, _)| *subcommand), expected, "{got:?}");
        if let Some((_, message)) = got {
            assert!(message.contains(&format!("'{} <PRIMARY[,ALT]>'", args[1])), "{message}");
        }
    }

    #[test]
    fn xa_max_hits_display_round_trips() {
        let hits = XaMaxHits { primary: 7, alt: 99 };
        assert_eq!(hits.to_string().parse::<XaMaxHits>().expect("parse"), hits);
    }
}
