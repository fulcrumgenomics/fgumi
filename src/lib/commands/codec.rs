//! # CODEC Consensus Calling Command
//!
//! Implementation of the `codec` command for calling consensus from CODEC sequencing data.
//!
//! CODEC (Bae et al 2023) is a sequencing protocol where each read-pair sequences both strands
//! of the original duplex molecule. R1 comes from one strand, R2 from the opposite strand,
//! allowing even a single read-pair to generate duplex consensus.
//!
//! When `--rejects` is set, the declarative chain (see [`crate::pipeline::chains`])
//! routes rejected records through its rejects fan-out branch. Both reject
//! paths (success and duplex-disagreement recovery) land in batch-input order.
//! The rejects BAM advertises the input header so raw-input RG/PG/contig
//! metadata is preserved. The pattern matches `commands::filter`,
//! `commands::correct`, `commands::simplex`, and `commands::duplex`.

use anyhow::{Result, bail};
use clap::Parser;
use std::path::Path;

use super::common::{
    AllowUnmappedOptions, BamIoOptions, CompressionOptions, ConsensusCallingOptions,
    QueueMemoryOptions, ReadGroupOptions, RejectsOptions, SchedulerOptions, StatsOptions,
    ThreadingOptions, reject_output_collisions,
};
use crate::commands::command::Command;
// Used only by unit tests that build stats structs directly, read back the
// `--stats` TSV, or set string tags on synthetic test records.
#[cfg(test)]
use crate::commands::consensus_runner::ConsensusStatsOps;
#[cfg(test)]
use crate::consensus::codec_caller::CodecConsensusStats;
#[cfg(test)]
use crate::sam::SamTag;
#[cfg(test)]
use fgoxide::io::DelimFile;

/// Call CODEC consensus reads from template-coordinate sorted BAM
///
/// CODEC is a sequencing protocol where a single read-pair sequences both strands of the
/// original duplex molecule. R1 sequences one strand, R2 sequences the opposite strand,
/// enabling duplex consensus calling from individual read-pairs.
///
/// The algorithm:
/// 1. Groups reads by molecule ID (MI tag)
/// 2. For each molecule:
///    - Clips read pairs where they extend past mate ends
///    - Filters R1s and R2s for compatible alignments
///    - Checks for sufficient overlap
///    - Calls single-strand consensus from R1s and R2s separately
///    - Combines into duplex consensus
///    - Applies quality masking
/// 3. Outputs unmapped consensus fragments (not paired-end)
///
/// ## Input Requirements
///
/// - BAM must be template-coordinate sorted (or queryname sorted)
/// - Reads must be aligned (mapped)
/// - Required tags:
///   - `RX`: Raw UMI bases
///   - `MI`: Molecule ID (from `group`)
///   - `CB`: Cell barcode (optional, for single-cell data)
///
/// ## Grouping Strategy
///
/// Must use `adjacency` or `identity` grouping (NOT `paired`).
/// MI tags should not have /A or /B suffixes.
///
/// ## Output
///
/// - Output BAM contains unmapped consensus fragments (not paired-end)
/// - Each consensus represents a full duplex molecule
/// - Includes rich per-base and per-read tags
///
/// Consensus reads have a number of additional optional tags set in the resulting BAM file.
/// The tag names follow a pattern where the first letter (a, b or c) denotes that the tag
/// applies to the first single strand consensus (a), second single-strand consensus (b) or
/// the final duplex consensus (c). The second letter captures the meaning of the tag
/// (e.g. d=depth, m=min depth, e=errors/error-rate) and is upper case for values that are
/// one per read and lower case for values that are one per base.
///
/// The tags break down into those that are single-valued per read:
///
///   consensus depth      \[aD,bD,cD\] (int)  : the maximum depth of raw reads
///   consensus min depth  \[aM,bM,cM\] (int)  : the minimum depth of raw reads
///   consensus error rate \[aE,bE,cE\] (float): the fraction of bases disagreeing with consensus
///
/// And those that have a value per base:
///
///   consensus depth  \[ad,bd\] (short[]): the count of bases contributing to consensus
///   consensus errors \[ae,be\] (short[]): the count of disagreeing bases
///   consensus bases  \[ac,bc\] (string) : the single-strand consensus bases
///   consensus quals  \[aq,bq\] (string) : the single-strand consensus qualities
///
/// ## Examples
///
/// Basic usage:
/// ```bash
/// fgumi codec -i grouped.bam -o codec_consensus.bam
/// ```
///
/// With quality thresholds and statistics:
/// ```bash
/// fgumi codec \
///   -i grouped.bam \
///   -o codec_consensus.bam \
///   -r rejects.bam \
///   -s stats.txt \
///   -m 15 \
///   -M 3 \
///   -d 10
/// ```
#[derive(Parser, Debug)]
#[command(
    name = "codec",
    about = "\x1b[38;5;180m[CONSENSUS]\x1b[0m      \x1b[36mCall CODEC consensus reads from grouped BAM\x1b[0m",
    long_about = r#"
Calls consensus reads from CODEC (Concatenating Original Duplex for Error Correction) data. CODEC
libraries place both strands of a source duplex molecule into a single read pair, so unlike `duplex`
-- which combines two separately sequenced single-strand consensus reads -- the two strands arrive
already paired: R1 carries one strand and R2 carries the other.

Consequently each input read pair yields at most one consensus fragment, and the duplex evidence
comes from comparing R1 against R2 rather than from grouping reads across the file. Prior to running
this tool, reads must have been grouped with `group`.

The consensus reads produced are unaligned fragments, so they should be aligned afterwards.

Quality masking
---------------

Several options overwrite base qualities in regions where the duplex evidence is weaker, rather
than discarding the bases outright. Each assigns its value unconditionally, so a position whose
quality is already lower is raised to it:

  --single-strand-qual         sets quality to this value in regions covered by only one strand,
                               where there is no duplex confirmation
  --outer-bases-qual /         sets quality to this value for the first and last
  --outer-bases-length         --outer-bases-length bases (5 by default) of the fragment, which
                               are the most error-prone

Duplex agreement filters
------------------------

  --min-duplex-length          minimum overlap, in bases, between the two strands for a fragment to
                               be emitted (default 1)
  --max-duplex-disagreement-rate  maximum fraction of overlapping positions where the strands may
                               disagree (default 1.0, i.e. no limit)
  --max-duplex-disagreements   maximum absolute count of such disagreements
  --legacy-overlap-window      use fgbio's legacy `[neg.start, pos.end]` overlap window, which
                               rejects dovetailed short-insert pairs as dovetails and matches
                               fgbio's current-release consensus output (off by default; the
                               corrected intersection window admits those pairs)

Note that `--methylation-mode` is not supported for CODEC data.
"#
)]
pub struct Codec {
    /// Input/output BAM options
    #[command(flatten)]
    pub io: BamIoOptions,

    /// Options for writing rejected reads
    #[command(flatten)]
    pub rejects_opts: RejectsOptions,

    /// Options for writing statistics
    #[command(flatten)]
    pub stats_opts: StatsOptions,

    /// Read group and name prefix options
    #[command(flatten)]
    pub read_group: ReadGroupOptions,

    /// Consensus calling options
    #[command(flatten)]
    pub consensus: ConsensusCallingOptions,

    /// Threading options for parallel processing
    #[command(flatten)]
    pub threading: ThreadingOptions,

    /// Compression options for output
    #[command(flatten)]
    pub compression: CompressionOptions,

    // --- CODEC-specific options below ---
    /// Minimum read pairs per strand to form consensus (fgbio's `--min-read-pairs`)
    #[arg(short = 'M', long = "min-reads", default_value = "1")]
    pub min_reads: usize,

    /// Maximum read pairs per strand, downsample if exceeded (fgbio's `--max-read-pairs`).
    ///
    /// If more than this many read pairs are present in a tag family, the family is downsampled
    /// to exactly this many read pairs.
    ///
    /// Which pairs are retained is determined by a hash of the read names, so the selection is
    /// reproducible across runs, thread counts, and execution modes.
    ///
    /// The cap is applied after filtering each strand to its most common alignment, and to each
    /// strand independently — matching fgbio. Mates share a read name and therefore a rank, so
    /// both ends of a template are normally retained or discarded together; they can diverge
    /// when one end survives the alignment filter and the other does not.
    #[arg(long = "max-reads")]
    pub max_reads: Option<usize>,

    /// Minimum duplex overlap length in bases
    #[arg(short = 'd', long = "min-duplex-length", default_value = "1")]
    pub min_duplex_length: usize,

    /// Use fgbio's legacy `[negative.start, positive.end]` duplex overlap window instead of the
    /// corrected intersection window (fgumi#761).
    ///
    /// Reproduces fgbio's CODEC consensus *output* for the current fgbio release: dovetailed FR
    /// pairs — common on short inserts, where each strand's alignment runs past the far end of
    /// the other — are rejected rather than consensed, so far fewer consensus reads are emitted.
    /// Those rejections are reported under `raw_reads_rejected_for_dovetail`, whereas fgbio
    /// buckets them under `indel_error_between_strands`, so the `--stats` breakdown is not
    /// byte-identical to fgbio's. Off by default (fulcrumgenomics/fgbio#1173).
    #[arg(long = "legacy-overlap-window")]
    pub legacy_overlap_window: bool,

    /// Set single-strand region quality to this value (0-93). Assigned
    /// unconditionally, so a lower quality is raised to it.
    /// Note: This uses a different short flag than duplex's -q for min-base-quality.
    #[arg(long = "single-strand-qual")]
    pub single_strand_qual: Option<u8>,

    /// Set outer bases quality to this value (0-93). Assigned unconditionally,
    /// so a lower quality is raised to it.
    #[arg(short = 'Q', long = "outer-bases-qual")]
    pub outer_bases_qual: Option<u8>,

    /// Number of outer bases to set the quality of
    #[arg(short = 'O', long = "outer-bases-length", default_value = "5")]
    pub outer_bases_length: usize,

    /// Maximum duplex disagreement rate (0.0-1.0)
    #[arg(short = 'x', long = "max-duplex-disagreement-rate", default_value = "1.0")]
    pub max_duplex_disagreement_rate: f64,

    /// Maximum number of duplex disagreements
    #[arg(short = 'X', long = "max-duplex-disagreements")]
    pub max_duplex_disagreements: Option<usize>,

    /// Whether to process unmapped reads (the shared `--allow-unmapped` flag).
    #[command(flatten)]
    pub allow_unmapped: AllowUnmappedOptions,

    /// Scheduler and pipeline statistics options.
    #[command(flatten)]
    pub scheduler_opts: SchedulerOptions,

    /// Queue memory options.
    #[command(flatten)]
    pub queue_memory: QueueMemoryOptions,
}

// ─────────────────────────────────────────────────────────────────────────────
// CodecOptions — the stage's tuning knobs, projected out of the CLI struct
// ─────────────────────────────────────────────────────────────────────────────

/// Codec-stage tuning, independent of how the values were supplied.
///
/// See [`crate::commands::zipper::ZipperOptions`] for why this is a plain
/// struct rather than a flattened `clap::Args`. Note that the consensus-calling
/// knobs are held **flat** here even though [`Codec`] nests them behind
/// `#[command(flatten)]` sub-structs: the chain builder wants one bag per stage,
/// not a re-run of the CLI's grouping.
#[derive(Debug, Clone)]
pub struct CodecOptions {
    /// Pre-UMI error rate (phred).
    pub error_rate_pre_umi: u8,
    /// Post-UMI error rate (phred).
    pub error_rate_post_umi: u8,
    /// Minimum input base quality.
    pub min_input_base_quality: u8,
    /// Emit per-base consensus tags.
    pub output_per_base_tags: bool,
    /// Trim consensus reads.
    pub trim: bool,
    /// Minimum consensus base quality.
    pub min_consensus_base_quality: u8,
    /// How to resolve a near-tie between the two most likely consensus bases.
    pub tie_rule: fgumi_consensus::TieRule,
    /// Minimum reads per consensus.
    pub min_reads: usize,
    /// Cap on reads per consensus.
    pub max_reads: Option<usize>,
    /// Minimum duplex overlap length.
    pub min_duplex_length: usize,
    /// Reproduce fgbio's legacy (pre-fgumi#761) overlap window for dovetailed
    /// FR pairs. Off by default; see [`Codec::legacy_overlap_window`].
    pub legacy_overlap_window: bool,
    /// Quality cap for single-strand positions.
    pub single_strand_qual: Option<u8>,
    /// Quality cap for outer bases.
    pub outer_bases_qual: Option<u8>,
    /// How many bases count as outer.
    pub outer_bases_length: usize,
    /// Maximum duplex disagreement rate.
    pub max_duplex_disagreement_rate: f64,
    /// Maximum duplex disagreements.
    pub max_duplex_disagreements: Option<usize>,
    /// Let fully-unmapped primary templates through the pre-group filter.
    ///
    /// Carried as the whole flattened sub-struct, like `io` / `rejects_opts` /
    /// `read_group`, rather than as a bare `bool`.
    pub allow_unmapped: AllowUnmappedOptions,
    /// Input/output paths and reader mode.
    pub io: BamIoOptions,
    /// Optional rejects output.
    pub rejects_opts: RejectsOptions,
    /// Optional stats output.
    pub stats_opts: StatsOptions,
    /// Read-group identity for emitted reads.
    pub read_group: ReadGroupOptions,
}

impl Codec {
    /// Project the parsed CLI flags into [`CodecOptions`].
    #[must_use]
    pub fn to_codec_options(&self) -> CodecOptions {
        CodecOptions {
            error_rate_pre_umi: self.consensus.error_rate_pre_umi,
            error_rate_post_umi: self.consensus.error_rate_post_umi,
            min_input_base_quality: self.consensus.min_input_base_quality,
            output_per_base_tags: self.consensus.output_per_base_tags,
            trim: self.consensus.trim,
            min_consensus_base_quality: self.consensus.min_consensus_base_quality,
            tie_rule: self.consensus.tie_rule.into(),
            min_reads: self.min_reads,
            max_reads: self.max_reads,
            min_duplex_length: self.min_duplex_length,
            legacy_overlap_window: self.legacy_overlap_window,
            single_strand_qual: self.single_strand_qual,
            outer_bases_qual: self.outer_bases_qual,
            outer_bases_length: self.outer_bases_length,
            max_duplex_disagreement_rate: self.max_duplex_disagreement_rate,
            max_duplex_disagreements: self.max_duplex_disagreements,
            allow_unmapped: self.allow_unmapped.clone(),
            io: self.io.clone(),
            rejects_opts: self.rejects_opts.clone(),
            stats_opts: self.stats_opts.clone(),
            read_group: self.read_group.clone(),
        }
    }
}

impl CodecOptions {
    /// Validate the codec consensus parameter bounds. Single source of truth,
    /// shared by `Codec::validate` (the CLI pre-flight) and `add_codec` (the
    /// chain builder), so the two paths cannot drift and every direct
    /// `StageOptionsBag { codec: .. }` caller is guarded identically to the CLI.
    /// Pure function of the options — no I/O, no side effects.
    pub fn validate(&self) -> Result<()> {
        // Validate error rates.
        if self.error_rate_pre_umi == 0 {
            bail!("error-rate-pre-umi must be > 0");
        }
        if self.error_rate_post_umi == 0 {
            bail!("error-rate-post-umi must be > 0");
        }

        // Validate min/max reads.
        if self.min_reads == 0 {
            bail!("min-reads must be >= 1");
        }
        if let Some(max) = self.max_reads
            && max < self.min_reads
        {
            bail!("max-reads ({}) must be >= min-reads ({})", max, self.min_reads);
        }

        // Validate duplex length.
        if self.min_duplex_length == 0 {
            bail!("min-duplex-length must be >= 1");
        }

        // Validate disagreement rate. The non-finite check must come first: `NaN`
        // compares `false` against both bounds, so a bare range test admits it, and
        // every downstream `rate > threshold` comparison is then false — silently
        // disabling duplex-disagreement rejection so molecules that should be
        // rejected are emitted as consensus reads.
        if !self.max_duplex_disagreement_rate.is_finite() {
            bail!(
                "max-duplex-disagreement-rate must be a finite number, got {}",
                self.max_duplex_disagreement_rate
            );
        }
        if self.max_duplex_disagreement_rate < 0.0 || self.max_duplex_disagreement_rate > 1.0 {
            bail!("max-duplex-disagreement-rate must be between 0.0 and 1.0");
        }

        // Validate the fixed quality overrides. These are written straight into the
        // consensus record's QUAL array, so a value above the maximum Phred score
        // would emit BAM quality bytes outside the legal SAM `!`..`~` range.
        let max_phred = ConsensusCallingOptions::MAX_PHRED;
        for (name, value) in [
            ("single-strand-qual", self.single_strand_qual),
            ("outer-bases-qual", self.outer_bases_qual),
        ] {
            if let Some(qual) = value
                && qual > max_phred
            {
                bail!("{name} ({qual}) exceeds maximum Phred score ({max_phred})");
            }
        }

        Ok(())
    }

    /// Reconstruct the shared [`ConsensusCallingOptions`] from the inlined flat
    /// fields, so the chain builder can read `.consensus().error_rate_pre_umi`
    /// etc. `tie_rule` is stored resolved on `CodecOptions`; convert it back to
    /// the CLI-facing [`crate::commands::common::TieRuleArg`] via the 1:1 mapping.
    #[must_use]
    pub fn consensus(&self) -> ConsensusCallingOptions {
        ConsensusCallingOptions {
            error_rate_pre_umi: self.error_rate_pre_umi,
            error_rate_post_umi: self.error_rate_post_umi,
            min_input_base_quality: self.min_input_base_quality,
            output_per_base_tags: self.output_per_base_tags,
            trim: self.trim,
            min_consensus_base_quality: self.min_consensus_base_quality,
            tie_rule: self.tie_rule.into(),
        }
    }
}

impl Command for Codec {
    fn execute(&self, command_line: &str) -> Result<()> {
        // Validate inputs
        self.validate()?;
        // Validate the input exists (stdin paths are exempt — the reader
        // streams them in a single pass).
        self.io.validate()?;
        let mut outputs: Vec<(&Path, &str)> = vec![(self.io.output.as_path(), "--output")];
        if let Some(path) = &self.rejects_opts.rejects {
            outputs.push((path.as_path(), "--rejects"));
        }
        if let Some(path) = &self.stats_opts.stats {
            outputs.push((path.as_path(), "--stats"));
        }
        reject_output_collisions(&outputs)?;

        // The declarative chain is the only execution path: `execute` runs the
        // reader-free pre-flight above and then always dispatches to
        // `execute_chain`, with or without `--threads` (absent `--threads` runs
        // the chain at a single worker). `add_codec` logs the `Calling CODEC
        // consensus` timer, the `Starting CODEC consensus calling` banner +
        // option lines, and the codec finalize hook logs the summary and the
        // `--stats` TSV. Running any of those here would double-log and
        // pre-consume stdin, so `execute` only does the pre-flight above.
        self.execute_chain(command_line)
    }
}
impl Codec {
    /// Run the codec stage on the declarative chain builder — the only
    /// execution path, with or without `--threads`.
    ///
    /// The chain opens its own source, validates the template-coordinate sort
    /// order, calls consensus, writes the output BAM, and writes the
    /// rejects/stats via the codec finalize hook. `--threads 1` (or absent
    /// `--threads`) runs the chain at a single worker, which is the in-process
    /// parity oracle for the multi-worker case (see
    /// `test_codec_chain_matches_single_threaded`).
    fn execute_chain(&self, command_line: &str) -> Result<()> {
        use crate::pipeline::chains::{
            ChainSpec, SingleStageContext, Stage, StageOptionsBag, build_for,
        };
        self.io.log_effective_check_crc();
        let stage_opts =
            StageOptionsBag { codec: Some(self.to_codec_options()), ..Default::default() };
        let ctx = SingleStageContext {
            io: &self.io,
            threading: &self.threading,
            compression: &self.compression,
            scheduler: &self.scheduler_opts,
            queue_memory: &self.queue_memory,
            command_line,
        };
        let spec = ChainSpec::single_stage(Stage::Codec, stage_opts, &ctx);
        build_for(spec)?.run()
    }

    /// Validates command-line arguments.
    ///
    /// Delegates to [`CodecOptions::validate`], the single source of truth for
    /// the codec parameter bounds, so this CLI pre-flight and the chain builder
    /// (`add_codec`, via the same `CodecOptions::validate`) cannot drift.
    fn validate(&self) -> Result<()> {
        self.to_codec_options().validate()
    }
}

// Helper methods moved to consensus_runner module

#[cfg(test)]
mod tests {
    use super::*;

    /// Every tuning flag must survive the projection into [`CodecOptions`].
    /// See the simplex counterpart for why this parses rather than constructs.
    ///
    /// Every value here is deliberately **non-default**, and every field of
    /// `CodecOptions` is asserted. Both halves matter: an assertion that
    /// compares a default against a default passes even when the projection
    /// ignores the parsed field entirely, and a field left unasserted can be
    /// dropped from the projection without any test noticing.
    #[test]
    fn to_codec_options_carries_every_tuning_flag() {
        let cmd = Codec::try_parse_from([
            "codec",
            "-i",
            "in.bam",
            "-o",
            "out.bam",
            "--error-rate-pre-umi",
            "42",
            "--error-rate-post-umi",
            "37",
            "--min-input-base-quality",
            "16",
            "--output-per-base-tags=false",
            "--trim=true",
            "--min-consensus-base-quality",
            "23",
            "--tie-rule",
            "ulp-relative",
            "--min-reads",
            "4",
            "--max-reads",
            "88",
            "--min-duplex-length",
            "9",
            "--legacy-overlap-window",
            "--single-strand-qual",
            "12",
            "--outer-bases-qual",
            "15",
            "--outer-bases-length",
            "7",
            "--max-duplex-disagreement-rate",
            "0.25",
            "--max-duplex-disagreements",
            "6",
            "--rejects",
            "rej.bam",
            "--stats",
            "stats.txt",
            "--read-group-id",
            "Z",
            "--read-name-prefix",
            "pfx",
            "--allow-unmapped=true",
        ])
        .expect("parses");

        let opts = cmd.to_codec_options();

        assert_eq!(opts.error_rate_pre_umi, 42);
        assert_eq!(opts.error_rate_post_umi, 37);
        assert_eq!(opts.min_input_base_quality, 16);
        assert!(!opts.output_per_base_tags, "an explicit false must not be lost");
        assert!(opts.trim);
        assert_eq!(opts.min_consensus_base_quality, 23);
        assert_eq!(
            opts.tie_rule,
            fgumi_consensus::TieRule::UlpRelative,
            "--tie-rule must reach the projection"
        );
        assert_eq!(opts.min_reads, 4);
        assert_eq!(opts.max_reads, Some(88));
        assert_eq!(opts.min_duplex_length, 9);
        assert!(opts.legacy_overlap_window, "--legacy-overlap-window must reach the projection");
        assert_eq!(opts.single_strand_qual, Some(12));
        assert_eq!(opts.outer_bases_qual, Some(15));
        assert_eq!(opts.outer_bases_length, 7);
        assert!((opts.max_duplex_disagreement_rate - 0.25).abs() < f64::EPSILON);
        assert_eq!(opts.max_duplex_disagreements, Some(6));
        assert!(opts.allow_unmapped.enabled, "--allow-unmapped must reach the projection");
        // The flattened sub-structs must come across whole, not field by field.
        assert_eq!(opts.io.input, std::path::PathBuf::from("in.bam"));
        assert_eq!(opts.io.output, std::path::PathBuf::from("out.bam"));
        assert_eq!(opts.rejects_opts.rejects, Some(std::path::PathBuf::from("rej.bam")));
        assert_eq!(opts.stats_opts.stats, Some(std::path::PathBuf::from("stats.txt")));
        assert_eq!(opts.read_group.read_group_id, "Z");
        assert_eq!(opts.read_group.read_name_prefix, Some("pfx".to_string()));
    }

    /// The projection must carry defaults faithfully too — a field hard-coded to
    /// the value the non-default test happens to pass would slip through it.
    #[test]
    fn to_codec_options_carries_defaults() {
        let cmd =
            Codec::try_parse_from(["codec", "-i", "in.bam", "-o", "out.bam"]).expect("parses");

        let opts = cmd.to_codec_options();

        assert_eq!(opts.error_rate_pre_umi, 45);
        assert_eq!(opts.error_rate_post_umi, 40);
        assert_eq!(opts.min_input_base_quality, 10);
        assert!(opts.output_per_base_tags);
        assert!(!opts.trim);
        assert_eq!(opts.min_consensus_base_quality, 2);
        assert_eq!(opts.tie_rule, fgumi_consensus::TieRule::FgbioCompat);
        assert_eq!(opts.min_reads, 1);
        assert_eq!(opts.max_reads, None);
        assert_eq!(opts.min_duplex_length, 1);
        assert!(!opts.legacy_overlap_window, "the default must be preserved, not hard-coded true");
        assert_eq!(opts.single_strand_qual, None);
        assert_eq!(opts.outer_bases_qual, None);
        assert_eq!(opts.outer_bases_length, 5);
        assert!((opts.max_duplex_disagreement_rate - 1.0).abs() < f64::EPSILON);
        assert_eq!(opts.max_duplex_disagreements, None);
        assert!(!opts.allow_unmapped.enabled);
        assert_eq!(opts.rejects_opts.rejects, None);
        assert_eq!(opts.stats_opts.stats, None);
        assert_eq!(opts.read_group.read_group_id, "A");
        assert_eq!(opts.read_group.read_name_prefix, None);
    }

    use noodles::sam::alignment::io::Write as AlignmentWrite;
    use rstest::rstest;
    use std::path::PathBuf;

    /// Helper to create a Codec with specified input/output paths
    fn create_codec_with_paths(input: PathBuf, output: PathBuf) -> Codec {
        Codec {
            io: BamIoOptions {
                input,
                output,
                async_reader: false,
                check_crc: false,
                no_check_crc: false,
            },
            rejects_opts: RejectsOptions::default(),
            stats_opts: StatsOptions::default(),
            read_group: ReadGroupOptions::default(),
            consensus: ConsensusCallingOptions {
                output_per_base_tags: false,
                min_consensus_base_quality: 0,
                ..ConsensusCallingOptions::default()
            },
            threading: ThreadingOptions::none(),
            compression: CompressionOptions { compression_level: 1 },
            min_reads: 1,
            max_reads: None,
            min_duplex_length: 1,
            legacy_overlap_window: false,
            single_strand_qual: None,
            outer_bases_qual: None,
            outer_bases_length: 5,
            max_duplex_disagreement_rate: 1.0,
            max_duplex_disagreements: None,
            allow_unmapped: AllowUnmappedOptions { enabled: false },
            scheduler_opts: SchedulerOptions::default(),
            queue_memory: QueueMemoryOptions::default(),
        }
    }

    fn create_test_codec() -> Codec {
        create_codec_with_paths(PathBuf::from("input.bam"), PathBuf::from("output.bam"))
    }

    #[rstest]
    #[case::typical(0.1, true)]
    #[case::lower_bound(0.0, true)]
    #[case::upper_bound(1.0, true)]
    #[case::negative(-0.1, false)]
    #[case::above_one(1.1, false)]
    #[case::nan(f64::NAN, false)]
    #[case::infinity(f64::INFINITY, false)]
    fn test_validate_max_duplex_disagreement_rate(#[case] rate: f64, #[case] expected_ok: bool) {
        let mut codec = create_test_codec();
        codec.max_duplex_disagreement_rate = rate;
        assert_eq!(
            codec.validate().is_ok(),
            expected_ok,
            "unexpected validation result for rate {rate}"
        );
    }

    #[rstest]
    #[case::ss_qual_valid(Some(40), None, true)]
    #[case::ss_qual_at_max(Some(93), None, true)]
    #[case::ss_qual_too_high(Some(94), None, false)]
    #[case::outer_qual_valid(None, Some(2), true)]
    #[case::outer_qual_at_max(None, Some(93), true)]
    #[case::outer_qual_too_high(None, Some(94), false)]
    #[case::both_unset(None, None, true)]
    fn test_validate_fixed_quality_overrides(
        #[case] single_strand_qual: Option<u8>,
        #[case] outer_bases_qual: Option<u8>,
        #[case] expected_ok: bool,
    ) {
        let mut codec = create_test_codec();
        codec.single_strand_qual = single_strand_qual;
        codec.outer_bases_qual = outer_bases_qual;
        assert_eq!(codec.validate().is_ok(), expected_ok);
    }

    /// The `--allow-unmapped` CLI flag defaults to `false` (fgbio parity: the
    /// pre-group filter is applied before consensus calling).
    #[test]
    fn test_codec_allow_unmapped_defaults_to_false() {
        let default_cmd =
            <Codec as clap::Parser>::try_parse_from(["codec", "-i", "in.bam", "-o", "out.bam"])
                .expect("codec should parse with only the required args");
        assert!(!default_cmd.allow_unmapped.enabled, "--allow-unmapped must default to false");

        let enabled_cmd = <Codec as clap::Parser>::try_parse_from([
            "codec",
            "-i",
            "in.bam",
            "-o",
            "out.bam",
            "--allow-unmapped",
        ])
        .expect("codec should parse with --allow-unmapped");
        assert!(enabled_cmd.allow_unmapped.enabled, "bare --allow-unmapped must enable the flag");
    }

    /// The codec pre-group filter drops secondary/supplementary alignments before
    /// grouping (fgbio parity), and `--allow-unmapped` must **not** relax that
    /// exclusion — it only relaxes the mapped-record rule.
    ///
    /// This is the codec counterpart to simplex's
    /// `test_allow_unmapped_gates_pregroup_filter` /
    /// duplex's `test_duplex_allow_unmapped_gates_pregroup_filter`. Those tests
    /// gate on unmapped, unpaired reads, but the codec caller itself rejects
    /// unmapped pairs (`codec_caller.rs`: "Reads must be mapped and overlap on the
    /// genome"), so an unmapped fixture yields zero consensus regardless of the
    /// flag and cannot exercise the gate. Instead we add a supplementary
    /// alignment to a valid mapped molecule and assert it never reaches the
    /// consensus caller (observed via `total_input_reads`, which counts
    /// post-filter reads) — with the flag both off and on. This pins the leak
    /// where a bypassed filter would let non-primary alignments into grouping.
    /// Runs the single-worker chain (absent `--threads`) and the multi-worker
    /// chain, which install the filter at the same `add_codec` site — worker
    /// count is the only difference between the two cases.
    #[rstest]
    #[case::no_threads_default(ThreadingOptions::none(), false)]
    #[case::no_threads_allow_unmapped(ThreadingOptions::none(), true)]
    #[case::threaded_default(ThreadingOptions::new(2), false)]
    #[case::threaded_allow_unmapped(ThreadingOptions::new(2), true)]
    fn test_codec_allow_unmapped_gates_pregroup_filter(
        #[case] threading: ThreadingOptions,
        #[case] allow_unmapped: bool,
    ) -> Result<()> {
        let dir = TempDir::new()?;
        let input_path = dir.path().join("input.bam");
        let output_path = dir.path().join("output.bam");
        let stats_path = dir.path().join("stats.txt");

        // A valid mapped, overlapping FR pair (one molecule) plus a supplementary
        // alignment of R1, all sharing one MI. The pre-group filter must drop the
        // supplementary before grouping, so only the two primaries reach the caller.
        let (r1, r2) = create_codec_fr_pair_overlapping("UMI001", 100, 105, 20, &[30; 20]);
        let supplementary = create_codec_supplementary_read("UMI001", 100, 20, &[30; 20]);
        write_codec_bam(&input_path, vec![r1, r2, supplementary])?;

        let mut cmd = create_codec_with_paths(input_path, output_path);
        cmd.threading = threading;
        cmd.allow_unmapped.enabled = allow_unmapped;
        cmd.stats_opts.stats = Some(stats_path.clone());
        cmd.read_group.read_name_prefix = Some("codec".to_string());
        cmd.outer_bases_length = 0;
        cmd.execute("test")?;

        // Codec emits fgbio's vertical key-value metrics format, so read the rows
        // back as `ConsensusKvMetric` and look up the two counts by key: the
        // supplementary read must be filtered before grouping, so exactly the two
        // primaries are considered (`raw_reads_considered`) and the single FR
        // molecule is emitted (`consensus_reads_emitted`).
        let kv_metrics: Vec<ConsensusKvMetric> = DelimFile::default().read_tsv(&stats_path)?;
        let value_for = |key: &str| -> u64 {
            kv_metrics
                .iter()
                .find(|m| m.key == key)
                .unwrap_or_else(|| panic!("stats file should contain a `{key}` row"))
                .value
                .parse()
                .unwrap_or_else(|_| panic!("`{key}` value should be an integer"))
        };
        assert_eq!(
            value_for("raw_reads_considered"),
            2,
            "the supplementary alignment must be dropped before grouping (allow_unmapped={allow_unmapped}): \
             only the two primary reads reach the codec caller"
        );
        assert_eq!(
            value_for("consensus_reads_emitted"),
            1,
            "the valid mapped FR pair is consensus-called (allow_unmapped={allow_unmapped})"
        );

        Ok(())
    }

    #[test]
    fn test_validation() {
        let mut cmd = create_test_codec();

        assert!(cmd.validate().is_ok());

        // Test invalid error rate
        cmd.consensus.error_rate_pre_umi = 0;
        assert!(cmd.validate().is_err());
        cmd.consensus.error_rate_pre_umi = 45;

        // Test invalid min reads
        cmd.min_reads = 0;
        assert!(cmd.validate().is_err());
        cmd.min_reads = 1;

        // Test invalid max < min
        cmd.max_reads = Some(0);
        assert!(cmd.validate().is_err());
        cmd.max_reads = None;

        // Test invalid duplex length
        cmd.min_duplex_length = 0;
        assert!(cmd.validate().is_err());
        cmd.min_duplex_length = 1;

        // Test invalid disagreement rate
        cmd.max_duplex_disagreement_rate = 1.5;
        assert!(cmd.validate().is_err());
        cmd.max_duplex_disagreement_rate = 1.0;

        assert!(cmd.validate().is_ok());
    }

    #[test]
    fn test_stats_to_metrics() {
        let stats = CodecConsensusStats {
            total_input_reads: 1000,
            consensus_reads_generated: 500,
            reads_filtered: 500,
            ..Default::default()
        };

        let metrics = stats.to_metrics();
        assert_eq!(metrics.total_input_reads, 1000);
        assert_eq!(metrics.consensus_reads, 500);
        assert_eq!(metrics.filtered_reads, 500);
    }

    // Integration tests
    use crate::metrics::consensus::ConsensusKvMetric;
    use fgumi_raw_bam::{
        SamBuilder as RawSamBuilder, flags, raw_record_to_record_buf, testutil::encode_op,
    };
    use noodles::sam::Header;
    use noodles::sam::alignment::record_buf::RecordBuf;
    use tempfile::TempDir;

    fn to_record_buf(raw: fgumi_raw_bam::RawRecord) -> RecordBuf {
        raw_record_to_record_buf(&raw, &Header::default())
            .expect("raw_record_to_record_buf failed in test")
    }

    /// Helper to create a CODEC-style FR read pair with proper overlap
    /// R1 forward at start1, R2 reverse at start2, with `read_len` bases each
    #[allow(clippy::cast_possible_truncation)]
    fn create_codec_fr_pair_overlapping(
        mi_value: &str,
        start1: usize,
        start2: usize,
        read_len: usize,
        quals: &[u8],
    ) -> (RecordBuf, RecordBuf) {
        // Use a simple reference-matching sequence
        let seq_forward = b"ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGT";
        let seq = &seq_forward[..read_len];

        // For R2 reverse, we need the reverse complement
        let seq_rc: Vec<u8> = seq
            .iter()
            .rev()
            .map(|&c| match c {
                b'A' => b'T',
                b'T' => b'A',
                b'C' => b'G',
                b'G' => b'C',
                other => other,
            })
            .collect();

        let insert_size: i32 = (start2 + read_len) as i32 - start1 as i32;
        let cigar = encode_op(0, read_len);
        // R1: PAIRED | PROPERLY_SEGMENTED | FIRST_SEGMENT | MATE_REVERSE = 0x1|0x2|0x40|0x20 = 99
        const PROPERLY_PAIRED: u16 = 0x2;
        let r1_flags = flags::PAIRED | PROPERLY_PAIRED | flags::FIRST_SEGMENT | flags::MATE_REVERSE;
        // R2: PAIRED | PROPERLY_SEGMENTED | LAST_SEGMENT | REVERSE = 0x1|0x2|0x80|0x10 = 147
        let r2_flags = flags::PAIRED | PROPERLY_PAIRED | flags::LAST_SEGMENT | flags::REVERSE;

        let mut b1 = RawSamBuilder::new();
        b1.read_name(format!("read_{mi_value}").as_bytes())
            .flags(r1_flags)
            .ref_id(0)
            .pos(start1 as i32 - 1)
            .mapq(60)
            .cigar_ops(&[cigar])
            .sequence(seq)
            .qualities(quals)
            .mate_ref_id(0)
            .mate_pos(start2 as i32 - 1)
            .template_length(insert_size);
        b1.add_string_tag(SamTag::MI, mi_value.as_bytes());
        b1.add_string_tag(SamTag::RG, b"A");
        let r1 = to_record_buf(b1.build());

        let mut b2 = RawSamBuilder::new();
        b2.read_name(format!("read_{mi_value}").as_bytes())
            .flags(r2_flags)
            .ref_id(0)
            .pos(start2 as i32 - 1)
            .mapq(60)
            .cigar_ops(&[cigar])
            .sequence(&seq_rc)
            .qualities(quals)
            .mate_ref_id(0)
            .mate_pos(start1 as i32 - 1)
            .template_length(-insert_size);
        b2.add_string_tag(SamTag::MI, mi_value.as_bytes());
        b2.add_string_tag(SamTag::RG, b"A");
        let r2 = to_record_buf(b2.build());

        (r1, r2)
    }

    /// A fully-unmapped primary pair (both ends `UNMAPPED`, both mates
    /// `MATE_UNMAPPED`) sharing `mi_value`. Used to verify that `--allow-unmapped`
    /// admits such a pair past the pre-group filter but the codec caller still
    /// declines to call it.
    fn create_codec_unmapped_pair(
        mi_value: &str,
        read_len: usize,
        quals: &[u8],
    ) -> (RecordBuf, RecordBuf) {
        let seq_forward = b"ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGT";
        let seq = &seq_forward[..read_len];

        let build = |flags: u16| {
            let mut b = RawSamBuilder::new();
            // ref_id/pos/mate_ref_id/mate_pos of -1 are BAM's "no alignment", and
            // an unmapped record carries no cigar.
            b.read_name(format!("read_{mi_value}").as_bytes())
                .flags(flags)
                .ref_id(-1)
                .pos(-1)
                .mapq(0)
                .sequence(seq)
                .qualities(quals)
                .mate_ref_id(-1)
                .mate_pos(-1)
                .template_length(0);
            b.add_string_tag(SamTag::MI, mi_value.as_bytes());
            b.add_string_tag(SamTag::RG, b"A");
            to_record_buf(b.build())
        };

        let r1 =
            build(flags::PAIRED | flags::FIRST_SEGMENT | flags::UNMAPPED | flags::MATE_UNMAPPED);
        let r2 =
            build(flags::PAIRED | flags::LAST_SEGMENT | flags::UNMAPPED | flags::MATE_UNMAPPED);
        (r1, r2)
    }

    /// `--allow-unmapped` admits unmapped primary records into grouping, but the
    /// codec caller requires a mapped primary FR pair (`is_primary_fr_pair_raw`
    /// rejects an unmapped read or an unmapped mate), so `codec` emits no
    /// consensus for a fully-unmapped pair however the flag is set.
    ///
    /// This pins the split documented on
    /// [`consensus_pregroup_keep_flags`](crate::commands::common::consensus_pregroup_keep_flags):
    /// the flag governs the *filter*, the caller governs whether a consensus is
    /// produced. `raw_reads_considered` counts post-filter reads, so the two
    /// assertions separate the two stages — with the flag on, both reads reach
    /// the caller and are still not called; with it off, the filter drops them
    /// before grouping.
    #[rstest]
    #[case::default_filters_before_grouping(false, 0)]
    #[case::allow_unmapped_admits_then_caller_rejects(true, 2)]
    fn test_codec_emits_no_consensus_for_unmapped_pair(
        #[case] allow_unmapped: bool,
        #[case] expected_reads_considered: u64,
    ) -> Result<()> {
        let dir = TempDir::new()?;
        let input_path = dir.path().join("input.bam");
        let output_path = dir.path().join("output.bam");
        let stats_path = dir.path().join("stats.txt");

        let (r1, r2) = create_codec_unmapped_pair("UMI001", 20, &[30; 20]);
        write_codec_bam(&input_path, vec![r1, r2])?;

        let mut cmd = create_codec_with_paths(input_path, output_path);
        cmd.allow_unmapped.enabled = allow_unmapped;
        cmd.stats_opts.stats = Some(stats_path.clone());
        cmd.read_group.read_name_prefix = Some("codec".to_string());
        cmd.outer_bases_length = 0;
        cmd.execute("test")?;

        let kv_metrics: Vec<ConsensusKvMetric> = DelimFile::default().read_tsv(&stats_path)?;
        let value_for = |key: &str| -> u64 {
            kv_metrics
                .iter()
                .find(|m| m.key == key)
                .unwrap_or_else(|| panic!("stats file should contain a `{key}` row"))
                .value
                .parse()
                .unwrap_or_else(|_| panic!("`{key}` value should be an integer"))
        };
        assert_eq!(
            value_for("raw_reads_considered"),
            expected_reads_considered,
            "--allow-unmapped governs only the pre-group filter (allow_unmapped={allow_unmapped})"
        );
        assert_eq!(
            value_for("consensus_reads_emitted"),
            0,
            "the codec caller never calls a fully-unmapped pair (allow_unmapped={allow_unmapped})"
        );

        Ok(())
    }

    /// A mapped supplementary alignment of R1 (`SUPPLEMENTARY` flag set) sharing
    /// `mi_value`. Used to verify the consensus pre-group filter drops
    /// non-primary alignments before grouping.
    #[allow(clippy::cast_possible_truncation)]
    fn create_codec_supplementary_read(
        mi_value: &str,
        start: usize,
        read_len: usize,
        quals: &[u8],
    ) -> RecordBuf {
        let seq_forward = b"ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGT";
        let seq = &seq_forward[..read_len];
        let cigar = encode_op(0, read_len);
        // PAIRED | FIRST_SEGMENT | SUPPLEMENTARY: a supplementary alignment of R1.
        let sup_flags = flags::PAIRED | flags::FIRST_SEGMENT | flags::SUPPLEMENTARY;
        let mut b = RawSamBuilder::new();
        b.read_name(format!("read_{mi_value}").as_bytes())
            .flags(sup_flags)
            .ref_id(0)
            .pos(start as i32 - 1)
            .mapq(60)
            .cigar_ops(&[cigar])
            .sequence(seq)
            .qualities(quals)
            .mate_ref_id(0)
            .mate_pos(start as i32 - 1);
        b.add_string_tag(SamTag::MI, mi_value.as_bytes());
        b.add_string_tag(SamTag::RG, b"A");
        to_record_buf(b.build())
    }

    fn write_codec_bam(path: &std::path::Path, records: Vec<RecordBuf>) -> Result<()> {
        use noodles::sam::header::record::value::Map;
        use noodles::sam::header::record::value::map::{ReferenceSequence, header::Version};
        use std::num::NonZeroUsize;

        let mut header = Header::builder()
            .set_header(Map::<noodles::sam::header::record::value::map::Header>::new(Version::new(
                1, 6,
            )))
            .build();

        // Add reference sequence
        let rs = Map::<ReferenceSequence>::new(
            NonZeroUsize::new(1000).expect("non-zero chromosome length"),
        );
        header.reference_sequences_mut().insert(bstr::BString::from("chr1"), rs);

        // Add read group
        let rg = Map::<noodles::sam::header::record::value::map::ReadGroup>::default();
        header.read_groups_mut().insert(bstr::BString::from("A"), rg);

        // Consensus callers require template-coordinate-sorted input (CONS-01).
        let header = fgumi_sam::header_as_template_coordinate(&header);

        let mut writer = noodles::bam::io::writer::Builder.build_from_path(path)?;
        writer.write_header(&header)?;

        for record in &records {
            writer.write_alignment_record(&header, record)?;
        }

        Ok(())
    }

    fn read_bam_records(path: &std::path::Path) -> Result<Vec<RecordBuf>> {
        let mut reader = noodles::bam::io::reader::Builder.build_from_path(path)?;
        let header = reader.read_header()?;
        let records: Vec<_> = reader.record_bufs(&header).collect::<std::io::Result<Vec<_>>>()?;
        Ok(records)
    }

    #[test]
    fn test_codec_execute_basic() -> Result<()> {
        let dir = TempDir::new()?;
        let input_path = dir.path().join("input.bam");
        let output_path = dir.path().join("output.bam");

        // Create FR pairs with same MI (molecule ID)
        // For proper overlap: R1 at pos 100 forward, R2 at pos 105 reverse, both 20bp
        // This gives overlap from pos 105-119 (15 bases overlap)
        let mut records = Vec::new();
        let (r1, r2) = create_codec_fr_pair_overlapping("UMI001", 100, 105, 20, &[30; 20]);
        records.push(r1);
        records.push(r2);

        write_codec_bam(&input_path, records)?;

        let mut cmd = create_codec_with_paths(input_path, output_path.clone());
        cmd.read_group.read_name_prefix = Some("codec".to_string());
        cmd.outer_bases_length = 0;

        // Execute should complete without errors (coverage is the goal)
        cmd.execute("test")?;

        // Output file should be created
        assert!(output_path.exists());

        Ok(())
    }

    #[rstest]
    #[case::fast_path(ThreadingOptions::none())]
    #[case::pipeline_1(ThreadingOptions::new(1))]
    #[case::pipeline_2(ThreadingOptions::new(2))]
    fn test_codec_execute_rejects_non_template_coordinate_input(
        #[case] threading: ThreadingOptions,
    ) -> Result<()> {
        // CONS-01: codec requires template-coordinate-sorted input. An ungrouped header must be
        // rejected by execute() (via check_consensus_sort_order) rather than silently
        // mis-grouping molecules; the accept branch is covered by the other execute tests.
        // `check_consensus_sort_order` is guarded on codec being the chain's source stage, so
        // parameterize over both the absent-`--threads` (single-worker) and `--threads` cases to
        // guard that the check still fires regardless of worker count.
        use noodles::sam::header::record::value::Map;
        use noodles::sam::header::record::value::map::{ReferenceSequence, header::Version};
        use std::num::NonZeroUsize;

        let dir = TempDir::new()?;
        let input_path = dir.path().join("input.bam");
        let output_path = dir.path().join("output.bam");

        // Plain header — deliberately NOT stamped template-coordinate (unlike write_codec_bam).
        let mut header = Header::builder()
            .set_header(Map::<noodles::sam::header::record::value::map::Header>::new(Version::new(
                1, 6,
            )))
            .build();
        header.reference_sequences_mut().insert(
            bstr::BString::from("chr1"),
            Map::<ReferenceSequence>::new(NonZeroUsize::new(1000).expect("non-zero length")),
        );
        header.read_groups_mut().insert(
            bstr::BString::from("A"),
            Map::<noodles::sam::header::record::value::map::ReadGroup>::default(),
        );
        let (r1, r2) = create_codec_fr_pair_overlapping("UMI001", 100, 105, 20, &[30; 20]);
        let mut writer = noodles::bam::io::writer::Builder.build_from_path(&input_path)?;
        writer.write_header(&header)?;
        writer.write_alignment_record(&header, &r1)?;
        writer.write_alignment_record(&header, &r2)?;
        drop(writer);

        let mut cmd = create_codec_with_paths(input_path, output_path);
        cmd.threading = threading;
        let err = cmd.execute("test").expect_err("codec must reject non-template-coordinate input");
        assert!(
            err.to_string().contains("template-coordinate"),
            "error must suggest template-coordinate sorting: {err}"
        );

        Ok(())
    }

    #[test]
    fn test_codec_execute_with_rejects() -> Result<()> {
        let dir = TempDir::new()?;
        let input_path = dir.path().join("input.bam");
        let output_path = dir.path().join("output.bam");
        let rejects_path = dir.path().join("rejects.bam");

        let mut records = Vec::new();
        // Create multiple FR pairs
        for i in 0..3 {
            let (r1, r2) = create_codec_fr_pair_overlapping(
                &format!("UMI{i:03}"),
                100 + i * 50,
                105 + i * 50,
                20,
                &[30; 20],
            );
            records.push(r1);
            records.push(r2);
        }

        write_codec_bam(&input_path, records)?;

        let mut cmd = create_codec_with_paths(input_path, output_path.clone());
        cmd.rejects_opts.rejects = Some(rejects_path.clone());
        cmd.read_group.read_name_prefix = Some("codec".to_string());
        cmd.outer_bases_length = 0;

        cmd.execute("test")?;

        // Rejects file should exist
        assert!(rejects_path.exists());

        Ok(())
    }

    #[test]
    fn test_codec_execute_with_stats() -> Result<()> {
        let dir = TempDir::new()?;
        let input_path = dir.path().join("input.bam");
        let output_path = dir.path().join("output.bam");
        let stats_path = dir.path().join("stats.txt");

        let mut records = Vec::new();
        let (r1, r2) = create_codec_fr_pair_overlapping("UMI001", 100, 105, 20, &[30; 20]);
        records.push(r1);
        records.push(r2);

        write_codec_bam(&input_path, records)?;

        let mut cmd = create_codec_with_paths(input_path, output_path.clone());
        cmd.stats_opts.stats = Some(stats_path.clone());
        cmd.read_group.read_name_prefix = Some("codec".to_string());
        cmd.outer_bases_length = 0;

        cmd.execute("test")?;

        // Stats file should exist and contain data
        assert!(stats_path.exists());
        let stats_content = std::fs::read_to_string(&stats_path)?;
        assert!(!stats_content.is_empty());

        Ok(())
    }

    #[test]
    fn test_codec_execute_multithreaded() -> Result<()> {
        let dir = TempDir::new()?;
        let input_path = dir.path().join("input.bam");
        let output_single = dir.path().join("output_single.bam");
        let output_multi = dir.path().join("output_multi.bam");
        let rejects_path = dir.path().join("rejects.bam");

        // Create many UMI groups to test parallel processing (25+ to trigger multiple batches with BATCH_SIZE=10)
        let mut records = Vec::new();
        for i in 0..25 {
            let (r1, r2) = create_codec_fr_pair_overlapping(
                &format!("UMI{i:03}"),
                100 + i * 50,
                105 + i * 50,
                20,
                &[30; 20],
            );
            records.push(r1);
            records.push(r2);
        }

        write_codec_bam(&input_path, records)?;

        // Single-threaded
        let mut cmd_single = create_codec_with_paths(input_path.clone(), output_single.clone());
        cmd_single.read_group.read_name_prefix = Some("codec".to_string());
        cmd_single.outer_bases_length = 0;

        cmd_single.execute("test")?;

        // Multi-threaded with rejects
        let mut cmd_multi = create_codec_with_paths(input_path, output_multi.clone());
        cmd_multi.rejects_opts.rejects = Some(rejects_path.clone());
        cmd_multi.read_group.read_name_prefix = Some("codec".to_string());
        cmd_multi.outer_bases_length = 0;
        cmd_multi.threading = ThreadingOptions::new(4);

        cmd_multi.execute("test")?;

        let records_single = read_bam_records(&output_single)?;
        let records_multi = read_bam_records(&output_multi)?;

        // Both should produce the same number of consensus reads
        assert_eq!(
            records_single.len(),
            records_multi.len(),
            "Single and multi-threaded should produce same number of reads"
        );

        Ok(())
    }

    #[test]
    fn test_codec_execute_with_per_base_tags() -> Result<()> {
        let dir = TempDir::new()?;
        let input_path = dir.path().join("input.bam");
        let output_path = dir.path().join("output.bam");

        let mut records = Vec::new();
        let (r1, r2) = create_codec_fr_pair_overlapping("UMI001", 100, 105, 20, &[30; 20]);
        records.push(r1);
        records.push(r2);

        write_codec_bam(&input_path, records)?;

        let mut cmd = create_codec_with_paths(input_path, output_path.clone());
        cmd.read_group.read_name_prefix = Some("codec".to_string());
        cmd.outer_bases_length = 0;
        cmd.consensus.output_per_base_tags = true; // Enable per-base tags

        cmd.execute("test")?;

        let output_records = read_bam_records(&output_path)?;
        assert!(!output_records.is_empty());

        Ok(())
    }

    #[test]
    fn test_codec_execute_with_trim() -> Result<()> {
        let dir = TempDir::new()?;
        let input_path = dir.path().join("input.bam");
        let output_path = dir.path().join("output.bam");

        let mut records = Vec::new();
        // Create reads with low quality at ends
        let quals = [5, 5, 30, 30, 30, 30, 30, 30, 30, 30, 30, 30, 30, 30, 30, 30, 5, 5, 5, 5]; // Low quality at ends
        let (r1, r2) = create_codec_fr_pair_overlapping("UMI001", 100, 105, 20, &quals);
        records.push(r1);
        records.push(r2);

        write_codec_bam(&input_path, records)?;

        let mut cmd = create_codec_with_paths(input_path, output_path.clone());
        cmd.read_group.read_name_prefix = Some("codec".to_string());
        cmd.outer_bases_length = 0;
        cmd.consensus.trim = true; // Enable trimming

        cmd.execute("test")?;

        // Should complete successfully (may or may not produce output depending on overlap after trim)
        assert!(output_path.exists());

        Ok(())
    }

    // Note: Overlapping consensus tests removed - CODEC does not support overlapping consensus
    // (matching fgbio's CallCodecConsensusReads which has no such option).

    #[test]
    fn test_codec_execute_with_max_reads() -> Result<()> {
        let dir = TempDir::new()?;
        let input_path = dir.path().join("input.bam");
        let output_path = dir.path().join("output.bam");

        // Create multiple read pairs for same UMI (to test max_reads downsampling)
        let mut records = Vec::new();
        for i in 0..5 {
            let (r1, r2) =
                create_codec_fr_pair_overlapping("UMI001", 100 + i, 105 + i, 20, &[30; 20]);
            records.push(r1);
            records.push(r2);
        }

        write_codec_bam(&input_path, records)?;

        let mut cmd = create_codec_with_paths(input_path, output_path.clone());
        cmd.read_group.read_name_prefix = Some("codec".to_string());
        cmd.outer_bases_length = 0;
        cmd.max_reads = Some(2); // Limit to 2 reads per strand

        cmd.execute("test")?;

        assert!(output_path.exists());

        Ok(())
    }

    #[test]
    fn test_codec_execute_validation_errors() {
        // Test various validation errors
        let mut cmd = create_test_codec();

        // Test post-UMI error rate = 0
        cmd.consensus.error_rate_post_umi = 0;
        assert!(cmd.validate().is_err());
        cmd.consensus.error_rate_post_umi = 40;

        // Test max < min reads
        cmd.min_reads = 5;
        cmd.max_reads = Some(2);
        assert!(cmd.validate().is_err());
        cmd.min_reads = 1;
        cmd.max_reads = None;

        // Test invalid disagreement rate
        cmd.max_duplex_disagreement_rate = -0.1;
        assert!(cmd.validate().is_err());
    }

    /// Parameterized test for all threading modes.
    ///
    /// Tests:
    /// - `None`: Single-threaded fast path, no pipeline
    /// - `Some(1)`: Pipeline with 1 thread
    /// - `Some(2)`: Pipeline with 2 threads
    #[rstest]
    #[case::fast_path(ThreadingOptions::none())]
    #[case::pipeline_1(ThreadingOptions::new(1))]
    #[case::pipeline_2(ThreadingOptions::new(2))]
    fn test_threading_modes(#[case] threading: ThreadingOptions) -> Result<()> {
        let dir = TempDir::new()?;
        let input_path = dir.path().join("input.bam");
        let output_path = dir.path().join("output.bam");

        // Create FR pairs with same MI (molecule ID)
        let mut records = Vec::new();
        let (r1, r2) = create_codec_fr_pair_overlapping("UMI001", 100, 105, 20, &[30; 20]);
        records.push(r1);
        records.push(r2);

        write_codec_bam(&input_path, records)?;

        let mut cmd = create_codec_with_paths(input_path, output_path.clone());
        cmd.read_group.read_name_prefix = Some("codec".to_string());
        cmd.outer_bases_length = 0;
        cmd.threading = threading;
        cmd.execute("test")?;

        assert!(output_path.exists());

        Ok(())
    }

    /// Asserts that single-threaded and multi-threaded codec produce the same number of
    /// consensus records and identical CB tag presence when some groups have CB and some do not.
    #[test]
    fn test_threading_parity_mixed_cb() -> Result<()> {
        use noodles::sam::alignment::record::data::field::Tag;
        use noodles::sam::alignment::record_buf::data::field::value::Value;

        let dir = TempDir::new()?;
        let input_path = dir.path().join("input.bam");

        // Group 1: MI001 with CB=CELL1
        let (mut r1, mut r2) = create_codec_fr_pair_overlapping("MI001", 100, 105, 20, &[30; 20]);
        let cb_tag = Tag::new(b'C', b'B');
        r1.data_mut().insert(cb_tag, Value::String("CELL1".into()));
        r2.data_mut().insert(cb_tag, Value::String("CELL1".into()));

        // Group 2: MI002 without CB
        let (r3, r4) = create_codec_fr_pair_overlapping("MI002", 200, 205, 20, &[30; 20]);

        write_codec_bam(&input_path, vec![r1, r2, r3, r4])?;

        // Run single-threaded
        let out_st = dir.path().join("out_st.bam");
        let stats_st = dir.path().join("stats_st.txt");
        let mut cmd_st = create_codec_with_paths(input_path.clone(), out_st.clone());
        cmd_st.outer_bases_length = 0;
        cmd_st.threading = ThreadingOptions::none();
        cmd_st.stats_opts.stats = Some(stats_st.clone());
        cmd_st.execute("test")?;
        let records_st = read_bam_records(&out_st)?;

        // Run multi-threaded
        let out_mt = dir.path().join("out_mt.bam");
        let stats_mt = dir.path().join("stats_mt.txt");
        let mut cmd_mt = create_codec_with_paths(input_path, out_mt.clone());
        cmd_mt.outer_bases_length = 0;
        cmd_mt.threading = ThreadingOptions::new(2);
        cmd_mt.stats_opts.stats = Some(stats_mt.clone());
        cmd_mt.execute("test")?;
        let records_mt = read_bam_records(&out_mt)?;

        // The metrics output must not depend on thread count, and must be the
        // fgbio KV format (seeded rejection rows), matching simplex/duplex.
        let stats_content = std::fs::read_to_string(&stats_st)?;
        assert_eq!(
            stats_content,
            std::fs::read_to_string(&stats_mt)?,
            "single-threaded and multi-threaded codec must write identical stats"
        );
        assert!(
            stats_content.contains("raw_reads_rejected_for_non_paired_reads"),
            "codec stats must be the fgbio KV format with the seeded non_paired_reads row"
        );

        assert_eq!(
            records_st.len(),
            records_mt.len(),
            "single-threaded and multi-threaded should produce the same number of consensus records"
        );

        // Assert record *identity*, not just count/CB presence: single- vs multi-threaded must
        // emit byte-identical consensus records, so content divergence (name/sequence/quality
        // reordering or corruption) is caught rather than masked by matching counts.
        // Tuple: (name, sequence, base qualities, CB-tag presence) per record.
        type RecordIdentity = (Vec<u8>, Vec<u8>, Vec<u8>, bool);
        let record_identity = |records: &[RecordBuf]| -> Vec<RecordIdentity> {
            records
                .iter()
                .map(|r| {
                    (
                        r.name().map(|n| n.to_vec()).unwrap_or_default(),
                        r.sequence().as_ref().to_vec(),
                        r.quality_scores().as_ref().to_vec(),
                        r.data().get(&cb_tag).is_some(),
                    )
                })
                .collect()
        };
        assert_eq!(
            record_identity(&records_st),
            record_identity(&records_mt),
            "single- and multi-threaded codec must emit byte-identical consensus records \
             (name, sequence, base qualities, and CB-tag presence)"
        );

        Ok(())
    }
}
