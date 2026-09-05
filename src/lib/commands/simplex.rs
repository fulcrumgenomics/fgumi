//! Call molecular consensus reads from UMI-grouped reads.
//!
//! This tool takes reads that have been grouped by UMI (via `group`) and generates
//! consensus reads using a likelihood-based model that accounts for sequencing errors and
//! errors introduced during sample preparation.
//!
//! When `--rejects` is set, the declarative chain (see [`crate::pipeline::chains`])
//! routes rejected records through its rejects fan-out branch. Both reject
//! sources (input-record skips below `--min-reads` and caller-emitted rejects)
//! land in batch-input order. The rejects BAM advertises the input header so
//! raw-input RG/PG/contig metadata is preserved. The pattern matches
//! `commands::filter` and `commands::correct`.

use anyhow::{Result, bail};
use clap::Parser;
use std::path::Path;

use crate::commands::command::Command;
use crate::commands::common::{
    AllowUnmappedOptions, BamIoOptions, CompressionOptions, ConsensusCallingOptions,
    OverlappingConsensusOptions, QueueMemoryOptions, ReadGroupOptions, RejectsOptions,
    SchedulerOptions, StatsOptions, ThreadingOptions, reject_output_collisions,
};
// Used only by unit tests that read back the `--stats` TSV.
#[cfg(test)]
use fgoxide::io::DelimFile;

/// Calls simplex consensus sequences from reads with the same unique molecular tag.
#[derive(Debug, Parser)]
#[command(
    name = "simplex",
    about = "\x1b[38;5;180m[CONSENSUS]\x1b[0m      \x1b[36mCall simplex consensus sequences from UMI-grouped reads\x1b[0m",
    long_about = r#"
Calls consensus sequences from reads with the same unique molecular tag.

Reads with the same unique molecular tag are examined base-by-base to assess the likelihood of each base in the
source molecule. The likelihood model is as follows:

1. First, the base qualities are adjusted. The base qualities are assumed to represent the probability of a
   sequencing error (i.e. the sequencer observed the wrong base present on the cluster/flowcell/well). The base
   quality scores are converted to probabilities incorporating a probability representing the chance of an error
   from the time the unique molecular tags were integrated to just prior to sequencing. The resulting probability
   is the error rate of all processes from right after integrating the molecular tag through to the end of
   sequencing.
2. Next, a consensus sequence is called for all reads with the same unique molecular tag base-by-base. For a
   given base position in the reads, the likelihoods that an A, C, G, or T is the base for the underlying
   source molecule respectively are computed by multiplying the likelihood of each read observing the base
   position being considered. The probability of error (from 1.) is used when the observed base does not match
   the hypothesized base for the underlying source molecule, while one minus that probability is used otherwise.
   The computed likelihoods are normalized by dividing them by the sum of all four likelihoods to produce a
   posterior probability, namely the probability that the source molecule was an A, C, G, or T from just after
   integrating molecular tag through to sequencing, given the observations. The base with the maximum posterior
   probability as the consensus call, and the posterior probability is used as its raw base quality.
3. Finally, the consensus raw base quality is modified by incorporating the probability of an error prior to
   integrating the unique molecular tags. Therefore, the probability used for the final consensus base
   quality is the posterior probability of the source molecule having the consensus base given the observed
   reads with the same molecular tag, all the way from sample extraction and through sample and library
   preparation, through preparing the library for sequencing (e.g. amplification, target selection), and finally,
   through sequencing.

This tool assumes that reads with the same tag are grouped together (consecutive in the file). Bases that
overlap within a read pair are consensus called jointly before UMI consensus calling; pass
--consensus-call-overlapping-bases false to call each end of a pair independently, as fgbio does. Insertion or
deletion errors in the reads are not considered in the consensus model.

The consensus reads produced are unaligned, due to the difficulty and error-prone nature of inferring the consensus
alignment. Consensus reads should therefore be aligned after, which should not be too expensive as likely there
are far fewer consensus reads than input raw reads.

Particular attention should be paid to setting the --min-reads parameter as this can have a dramatic effect on
both results and runtime. For libraries with low duplication rates (e.g. 100-300X exomes libraries) in which it
is desirable to retain singleton reads while making consensus reads from sets of duplicates, --min-reads=1 is
appropriate. For libraries with high duplication rates where it is desirable to only produce consensus reads
supported by 2+ reads to allow error correction, --min-reads=2 or higher is appropriate. After generation,
consensus reads can be further filtered using the filter tool. As such it is always safe to run
with --min-reads=1 and filter later, but filtering at this step can improve performance significantly.

Consensus reads have a number of additional optional tags set in the resulting BAM file. The tags break down into
those that are single-valued per read:

  consensus depth      [cD] (int)  : the maximum depth of raw reads at any point in the consensus read
  consensus min depth  [cM] (int)  : the minimum depth of raw reads at any point in the consensus read
  consensus error rate [cE] (float): the fraction of bases in raw reads disagreeing with the final consensus calls

And those that have a value per base:

  consensus depth  [cd] (short[]): the count of bases contributing to the consensus read at each position
  consensus errors [ce] (short[]): the number of bases from raw reads disagreeing with the final consensus base

The per base depths and errors are both capped at 32,767. In all cases no-calls (Ns) and bases below the
--min-input-base-quality are not counted in tag value calculations.
"#
)]
pub struct Simplex {
    /// Input/output BAM file paths
    #[command(flatten)]
    pub io: BamIoOptions,

    /// Optional output for rejected reads
    #[command(flatten)]
    pub rejects_opts: RejectsOptions,

    /// Optional output for statistics
    #[command(flatten)]
    pub stats_opts: StatsOptions,

    /// Read group and read name prefix options
    #[command(flatten)]
    pub read_group: ReadGroupOptions,

    /// Consensus calling options (error rates, quality thresholds)
    #[command(flatten)]
    pub consensus: ConsensusCallingOptions,

    /// Overlapping bases consensus options
    #[command(flatten)]
    pub overlapping: OverlappingConsensusOptions,

    /// Threading options for parallel processing
    #[command(flatten)]
    pub threading: ThreadingOptions,

    /// Compression options for output BAM.
    #[command(flatten)]
    pub compression: CompressionOptions,

    /// Minimum number of reads to produce a consensus (required, no default)
    /// Matches fgbio's `CallMolecularConsensusReads` which requires this argument
    #[arg(short = 'M', long = "min-reads")]
    pub min_reads: usize,

    /// Maximum reads to use per end (Fragment, R1, or R2) when building the consensus
    /// (downsample if exceeded).
    ///
    /// The cap is applied independently to each end rather than to the whole tag family, matching
    /// fgbio's per-end semantics (fgumi#723).
    ///
    /// Which reads are retained is determined by a hash of the read names, so the selection is
    /// reproducible across runs, thread counts, and execution modes. Mates share a read name and
    /// therefore a rank, so both ends of a template are normally retained or discarded together;
    /// because each end is capped independently, a template whose R1 survives upstream alignment
    /// filtering while its R2 does not can shift the cap boundary on one end alone, and then the
    /// two ends' retained subsets diverge.
    #[arg(long = "max-reads")]
    pub max_reads: Option<usize>,

    /// Whether to process unmapped reads (the shared `--allow-unmapped` flag).
    #[command(flatten)]
    pub allow_unmapped: AllowUnmappedOptions,

    /// Scheduler and pipeline statistics options.
    #[command(flatten)]
    pub scheduler_opts: SchedulerOptions,

    /// Queue memory options.
    #[command(flatten)]
    pub queue_memory: QueueMemoryOptions,

    /// Methylation-aware consensus calling mode.
    /// EM-Seq: C→T at ref-C = unmethylated (enzymatic conversion); TAPs: C→T at ref-C = methylated.
    /// Emits MM/ML methylation tags and cu/ct per-base count tags on consensus reads.
    /// Requires --ref.
    #[arg(long = "methylation-mode", value_enum)]
    pub methylation_mode: Option<crate::commands::common::MethylationModeArg>,

    /// Path to the reference FASTA file (required when --methylation-mode is set).
    #[arg(long = "ref")]
    pub reference: Option<std::path::PathBuf>,
}

// ─────────────────────────────────────────────────────────────────────────────
// SimplexOptions — the stage's tuning knobs, projected out of the CLI struct
// ─────────────────────────────────────────────────────────────────────────────

/// Simplex-stage tuning, independent of how the values were supplied.
///
/// See [`crate::commands::zipper::ZipperOptions`] for why this is a plain
/// struct rather than a flattened `clap::Args`. Note that the consensus-calling
/// knobs are held **flat** here even though [`Simplex`] nests them behind
/// `#[command(flatten)]` sub-structs: the chain builder wants one bag per stage,
/// not a re-run of the CLI's grouping.
#[derive(Debug, Clone)]
pub struct SimplexOptions {
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
    /// Call overlapping bases jointly.
    pub consensus_call_overlapping_bases: bool,
    /// Minimum reads per consensus.
    pub min_reads: usize,
    /// Cap on reads per consensus.
    pub max_reads: Option<usize>,
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
    /// Resolved methylation calling mode (`Disabled` when the flag is unset).
    pub methylation_mode: fgumi_consensus::MethylationMode,
    /// Reference FASTA for methylation-aware modes.
    pub reference: Option<std::path::PathBuf>,
}

impl Simplex {
    /// Project the parsed CLI flags into [`SimplexOptions`].
    #[must_use]
    pub fn to_simplex_options(&self) -> SimplexOptions {
        SimplexOptions {
            error_rate_pre_umi: self.consensus.error_rate_pre_umi,
            error_rate_post_umi: self.consensus.error_rate_post_umi,
            min_input_base_quality: self.consensus.min_input_base_quality,
            output_per_base_tags: self.consensus.output_per_base_tags,
            trim: self.consensus.trim,
            min_consensus_base_quality: self.consensus.min_consensus_base_quality,
            tie_rule: self.consensus.tie_rule.into(),
            consensus_call_overlapping_bases: self.overlapping.consensus_call_overlapping_bases,
            min_reads: self.min_reads,
            max_reads: self.max_reads,
            allow_unmapped: self.allow_unmapped.clone(),
            io: self.io.clone(),
            rejects_opts: self.rejects_opts.clone(),
            stats_opts: self.stats_opts.clone(),
            read_group: self.read_group.clone(),
            methylation_mode: crate::commands::common::resolve_methylation_mode(
                self.methylation_mode,
            ),
            reference: self.reference.clone(),
        }
    }
}

impl SimplexOptions {
    /// Reconstruct the shared [`ConsensusCallingOptions`] from the inlined flat
    /// fields, so the chain builder can read `.consensus().error_rate_pre_umi`
    /// etc. `tie_rule` is stored resolved on `SimplexOptions`; convert it back to
    /// the CLI-facing `TieRuleArg` via the 1:1 mapping.
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

    /// Reconstruct the shared [`OverlappingConsensusOptions`] from the inlined
    /// flat field.
    #[must_use]
    pub fn overlapping(&self) -> OverlappingConsensusOptions {
        OverlappingConsensusOptions {
            consensus_call_overlapping_bases: self.consensus_call_overlapping_bases,
        }
    }

    /// Validate the `--min-reads` / `--max-reads` range.
    ///
    /// `min_reads` must be `>= 1` (a value of 0 admits empty groups) and, when
    /// `max_reads` is set, it must be `>= min_reads`. Shared by the standalone
    /// `Simplex::execute` and the chain builder's `add_simplex` so `runall`
    /// rejects the same degenerate configurations the standalone command does.
    pub(crate) fn validate_read_bounds(&self) -> Result<()> {
        if self.min_reads == 0 {
            bail!("--min-reads must be >= 1 (a value of 0 admits empty groups)");
        }
        if let Some(max) = self.max_reads
            && max < self.min_reads
        {
            bail!("--max-reads ({}) must be >= --min-reads ({})", max, self.min_reads);
        }
        Ok(())
    }
}

impl Command for Simplex {
    fn execute(&self, command_line: &str) -> Result<()> {
        // ---- reader-free pre-flight (runs on BOTH paths) ----
        self.io.validate()?;
        let mut outputs: Vec<(&Path, &str)> = vec![(self.io.output.as_path(), "--output")];
        if let Some(path) = &self.rejects_opts.rejects {
            outputs.push((path.as_path(), "--rejects"));
        }
        if let Some(path) = &self.stats_opts.stats {
            outputs.push((path.as_path(), "--stats"));
        }
        reject_output_collisions(&outputs)?;

        self.validate_read_bounds()?;

        if self.reference.is_some() && self.methylation_mode.is_none() {
            bail!("--ref requires --methylation-mode to be set");
        }

        // The declarative chain is the only execution path: `execute` runs the
        // reader-free pre-flight above and then always dispatches to
        // `execute_chain`, with or without `--threads` (absent `--threads` runs
        // the chain at a single worker). `add_simplex` logs the `Calling simplex
        // consensus` timer, the `Starting Simplex` banner + Input/Output/Min
        // reads/Error rate lines, and `Processing reads and calling consensus
        // (streaming)…`; `SimplexFinalizeHook` logs the overlapping stats, the
        // summary, and the `--stats` TSV. Running any of those here first would
        // double-log and pre-consume stdin, so `execute` only does the
        // pre-flight above.
        self.execute_chain(command_line)
    }
}
impl Simplex {
    /// Validates the `--min-reads` / `--max-reads` family-size bounds.
    ///
    /// `--min-reads 0` is rejected because it admits empty groups, silently turning
    /// the minimum-family-size filter into a no-op (the `raw_records.len() < min_reads`
    /// test can never fire). `codec` already enforces the same lower bound.
    ///
    /// # Errors
    ///
    /// Returns an error if `min_reads` is 0, or if `max_reads` is below `min_reads`.
    fn validate_read_bounds(&self) -> Result<()> {
        self.to_simplex_options().validate_read_bounds()
    }

    /// Run the simplex stage on the declarative chain builder — the only
    /// execution path, with or without `--threads`.
    ///
    /// The chain opens its own source, validates the template-coordinate sort
    /// order, calls consensus, writes the output BAM, and writes the
    /// rejects/stats via `SimplexFinalizeHook`. `--threads 1` (or absent
    /// `--threads`) runs the chain at a single worker, which is the in-process
    /// parity oracle for the multi-worker case (see
    /// `test_simplex_chain_matches_single_threaded`).
    fn execute_chain(&self, command_line: &str) -> Result<()> {
        use crate::pipeline::chains::{
            ChainSpec, SingleStageContext, Stage, StageOptionsBag, build_for,
        };
        self.io.log_effective_check_crc();
        let stage_opts =
            StageOptionsBag { simplex: Some(self.to_simplex_options()), ..Default::default() };
        let ctx = SingleStageContext {
            io: &self.io,
            threading: &self.threading,
            compression: &self.compression,
            scheduler: &self.scheduler_opts,
            queue_memory: &self.queue_memory,
            command_line,
        };
        let spec = ChainSpec::single_stage(Stage::Simplex, stage_opts, &ctx);
        build_for(spec)?.run()
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    /// Every tuning flag must survive the projection into [`SimplexOptions`].
    ///
    /// Driven through `try_parse_from` rather than a struct literal: a literal
    /// would still compile if a flag were renamed or unwired, whereas parsing
    /// pins the whole path from command line to option struct. Non-default
    /// values throughout, so a field read from the wrong source fails rather
    /// than coincidentally matching its default.
    #[test]
    fn to_simplex_options_carries_every_tuning_flag() {
        let cmd = Simplex::try_parse_from([
            "simplex",
            "-i",
            "in.bam",
            "-o",
            "out.bam",
            "--error-rate-pre-umi",
            "40",
            "--error-rate-post-umi",
            "35",
            "--min-input-base-quality",
            "17",
            "--output-per-base-tags=false",
            "--trim=true",
            "--min-consensus-base-quality",
            "19",
            "--tie-rule",
            "ulp-relative",
            "--consensus-call-overlapping-bases=false",
            "--min-reads",
            "3",
            "--max-reads",
            "77",
            "--rejects",
            "rej.bam",
            "--stats",
            "stats.txt",
            "--read-group-id",
            "Z",
            "--read-name-prefix",
            "pfx",
            "--methylation-mode",
            "em-seq",
            "--ref",
            "ref.fa",
            "--allow-unmapped=true",
        ])
        .expect("parses");

        let opts = cmd.to_simplex_options();

        assert_eq!(opts.error_rate_pre_umi, 40);
        assert_eq!(opts.error_rate_post_umi, 35);
        assert_eq!(opts.min_input_base_quality, 17);
        assert!(!opts.output_per_base_tags, "an explicit false must not be lost");
        assert!(opts.trim);
        assert_eq!(opts.min_consensus_base_quality, 19);
        assert_eq!(
            opts.tie_rule,
            fgumi_consensus::TieRule::UlpRelative,
            "--tie-rule must reach the projection"
        );
        assert!(!opts.consensus_call_overlapping_bases);
        assert_eq!(opts.min_reads, 3);
        assert_eq!(opts.max_reads, Some(77));
        assert!(opts.allow_unmapped.enabled, "--allow-unmapped must reach the projection");
        assert_eq!(opts.reference, Some(std::path::PathBuf::from("ref.fa")));
        assert_eq!(
            opts.methylation_mode,
            fgumi_consensus::MethylationMode::EmSeq,
            "--methylation-mode must reach the projection",
        );
        // The flattened sub-structs must come across whole, not field by field.
        assert_eq!(opts.io.input, std::path::PathBuf::from("in.bam"));
        assert_eq!(opts.io.output, std::path::PathBuf::from("out.bam"));
        assert_eq!(opts.read_group.read_group_id, "Z");
        assert_eq!(opts.read_group.read_name_prefix, Some("pfx".to_string()));
        assert_eq!(opts.rejects_opts.rejects, Some(std::path::PathBuf::from("rej.bam")));
        assert_eq!(opts.stats_opts.stats, Some(std::path::PathBuf::from("stats.txt")));
    }

    /// The projection must carry defaults faithfully too — a field hard-coded to
    /// the value the non-default test happens to pass would slip through it.
    #[test]
    fn to_simplex_options_carries_defaults() {
        let cmd = Simplex::try_parse_from([
            "simplex",
            "-i",
            "in.bam",
            "-o",
            "out.bam",
            "--min-reads",
            "1",
        ])
        .expect("parses");

        let opts = cmd.to_simplex_options();

        assert_eq!(opts.error_rate_pre_umi, 45);
        assert_eq!(opts.error_rate_post_umi, 40);
        assert_eq!(opts.min_input_base_quality, 10);
        assert!(opts.output_per_base_tags);
        assert!(!opts.trim);
        assert_eq!(opts.min_consensus_base_quality, 2);
        assert_eq!(opts.tie_rule, fgumi_consensus::TieRule::FgbioCompat);
        assert!(opts.consensus_call_overlapping_bases);
        assert_eq!(opts.max_reads, None);
        assert!(!opts.allow_unmapped.enabled);
        assert_eq!(opts.methylation_mode, fgumi_consensus::MethylationMode::Disabled);
        assert_eq!(opts.reference, None);
        assert_eq!(opts.rejects_opts.rejects, None);
        assert_eq!(opts.stats_opts.stats, None);
        assert_eq!(opts.read_group.read_group_id, "A");
        assert_eq!(opts.read_group.read_name_prefix, None);
    }

    use crate::metrics::consensus::ConsensusKvMetric;
    use noodles::sam::alignment::record::data::field::Tag;
    use noodles::sam::alignment::record_buf::RecordBuf;
    use rstest::rstest;
    use std::path::PathBuf;

    /// Creates a test Simplex command instance with default parameters.
    ///
    /// Generates a Simplex command configured with standard test values,
    /// including file paths, UMI tag, error rates, and threading settings.
    ///
    /// # Returns
    ///
    /// A `Simplex` instance configured for testing
    fn create_test_simplex() -> Simplex {
        create_simplex_with_paths(PathBuf::from("test.bam"), PathBuf::from("out.bam"))
    }

    #[rstest]
    #[case::min_reads_zero_rejected(0, None, false)]
    #[case::min_reads_one(1, None, true)]
    #[case::min_reads_three(3, None, true)]
    #[case::max_below_min(3, Some(2), false)]
    #[case::max_equals_min(3, Some(3), true)]
    #[case::max_above_min(3, Some(9), true)]
    fn test_validate_read_bounds(
        #[case] min_reads: usize,
        #[case] max_reads: Option<usize>,
        #[case] expected_ok: bool,
    ) {
        let mut simplex = create_test_simplex();
        simplex.min_reads = min_reads;
        simplex.max_reads = max_reads;
        assert_eq!(
            simplex.validate_read_bounds().is_ok(),
            expected_ok,
            "unexpected result for min_reads={min_reads} max_reads={max_reads:?}"
        );
    }

    /// Creates a Simplex command with the given input/output paths and default parameters.
    fn create_simplex_with_paths(input: PathBuf, output: PathBuf) -> Simplex {
        Simplex {
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
            overlapping: OverlappingConsensusOptions { consensus_call_overlapping_bases: false },
            threading: ThreadingOptions::none(),
            compression: CompressionOptions { compression_level: 1 },
            queue_memory: QueueMemoryOptions::default(),
            min_reads: 1,
            max_reads: None,
            // These fixtures use `SamBuilder::new_unmapped()` synthetic reads to
            // exercise consensus behavior; opt back into consensus on unmapped
            // input so they survive the default fgbio pre-group filter. The
            // real CLI default (false) is pinned by
            // `test_allow_unmapped_defaults_to_false`.
            allow_unmapped: AllowUnmappedOptions { enabled: true },
            scheduler_opts: SchedulerOptions::default(),
            methylation_mode: None,
            reference: None,
        }
    }

    #[test]
    fn test_default_parameters() {
        let cmd = create_test_simplex();
        assert_eq!(cmd.consensus.error_rate_pre_umi, 45);
        assert_eq!(cmd.consensus.error_rate_post_umi, 40);
        assert_eq!(cmd.consensus.min_input_base_quality, 10);
        assert_eq!(cmd.min_reads, 1);
        assert_eq!(cmd.max_reads, None);
        assert!(!cmd.consensus.output_per_base_tags);
        assert!(cmd.threading.is_single_threaded());
        // Note: create_simplex_with_paths disables overlapping by default for tests
        assert!(!cmd.overlapping.consensus_call_overlapping_bases);
    }

    #[test]
    fn test_custom_parameters() {
        let mut cmd = create_test_simplex();
        cmd.min_reads = 3;
        cmd.max_reads = Some(10);
        cmd.consensus.output_per_base_tags = true;
        cmd.threading = ThreadingOptions::new(4);

        assert_eq!(cmd.min_reads, 3);
        assert_eq!(cmd.max_reads, Some(10));
        assert!(cmd.consensus.output_per_base_tags);
        assert_eq!(cmd.threading.threads, Some(4));
    }

    #[test]
    fn test_missing_input_file_fails() {
        let cmd = create_test_simplex();
        let result = cmd.execute("test");
        assert!(result.is_err());
        assert!(result.unwrap_err().to_string().contains("does not exist"));
    }

    #[test]
    fn test_rejects_non_template_coordinate_input() -> Result<()> {
        // CONS-01: consensus calling groups reads by molecule as they stream, so it
        // requires template-coordinate-sorted input. A bare (ungrouped) header must be
        // rejected — proceeding would silently split every molecule and corrupt output.
        let mut builder = SamBuilder::new_unmapped();
        let mut attrs = HashMap::new();
        attrs.insert("MI", BufValue::from("GATTACA:0"));
        attrs.insert("RX", BufValue::from("ACGT-TGCA"));
        builder.add_pair_with_attrs("READ:0", None, None, true, true, &attrs);

        let paths = TestPaths::new()?;
        // Deliberately NOT stamped with template-coordinate order.
        builder.write(&paths.input)?;

        let cmd = create_simplex_with_paths(paths.input.clone(), paths.output.clone());
        let result = cmd.execute("test");
        assert!(result.is_err(), "simplex must reject non-template-coordinate input");
        assert!(result.unwrap_err().to_string().contains("template-coordinate"));
        Ok(())
    }

    // ========================================================================
    // Integration Tests
    // ========================================================================

    use crate::sam::builder::SamBuilder;
    use noodles::sam::alignment::record_buf::data::field::Value as BufValue;
    use std::collections::HashMap;
    use tempfile::TempDir;

    /// Helper struct for managing temporary test file paths.
    struct TestPaths {
        #[allow(dead_code)]
        dir: TempDir,
        pub input: PathBuf,
        pub output: PathBuf,
        pub rejects: PathBuf,
        pub stats: PathBuf,
    }

    impl TestPaths {
        fn new() -> Result<Self> {
            let dir = TempDir::new()?;
            Ok(Self {
                input: dir.path().join("input.bam"),
                output: dir.path().join("output.bam"),
                rejects: dir.path().join("rejects.bam"),
                stats: dir.path().join("stats.txt"),
                dir,
            })
        }
    }

    /// Helper to read all records from a BAM file
    fn read_bam_records(path: &std::path::Path) -> Result<Vec<RecordBuf>> {
        let mut reader = noodles::bam::io::reader::Builder.build_from_path(path)?;
        let header = reader.read_header()?;
        let mut records = Vec::new();

        for result in reader.records() {
            let record = result?;
            let record_buf = RecordBuf::try_from_alignment_record(&header, &record)?;
            records.push(record_buf);
        }

        Ok(records)
    }

    /// Helper to get a string tag from a record
    fn get_string_tag(record: &RecordBuf, tag_name: &str) -> Option<String> {
        let tag_bytes = tag_name.as_bytes();
        let tag = Tag::from([tag_bytes[0], tag_bytes[1]]);

        record.data().get(&tag).and_then(|v| {
            if let noodles::sam::alignment::record_buf::data::field::Value::String(s) = v {
                Some(String::from_utf8_lossy(s).to_string())
            } else {
                None
            }
        })
    }

    /// Helper to get an integer tag from a record
    fn get_int_tag(record: &RecordBuf, tag_name: &str) -> Option<i64> {
        let tag_bytes = tag_name.as_bytes();
        let tag = Tag::from([tag_bytes[0], tag_bytes[1]]);

        record.data().get(&tag).and_then(|v| match v {
            noodles::sam::alignment::record_buf::data::field::Value::Int8(i) => Some(*i as i64),
            noodles::sam::alignment::record_buf::data::field::Value::UInt8(i) => Some(*i as i64),
            noodles::sam::alignment::record_buf::data::field::Value::Int16(i) => Some(*i as i64),
            noodles::sam::alignment::record_buf::data::field::Value::UInt16(i) => Some(*i as i64),
            noodles::sam::alignment::record_buf::data::field::Value::Int32(i) => Some(*i as i64),
            noodles::sam::alignment::record_buf::data::field::Value::UInt32(i) => Some(*i as i64),
            _ => None,
        })
    }

    #[test]
    fn test_end_to_end_paired_end_workflow() -> Result<()> {
        // Create test data similar to Scala test: 1000 UMI groups with 2 pairs each
        let mut builder = SamBuilder::new_unmapped();

        for idx in 0..1000 {
            let umi = format!("GATTACA:{idx}");
            let mut attrs = HashMap::new();
            attrs.insert("MI", BufValue::from(umi.clone()));
            attrs.insert("RX", BufValue::from("ACGT-TGCA"));

            // Add 2 pairs per UMI group
            builder.add_pair_with_attrs(
                &format!("READ:{}", 2 * idx),
                None,
                None,
                true,
                true,
                &attrs,
            );
            builder.add_pair_with_attrs(
                &format!("READ:{}", 2 * idx + 1),
                None,
                None,
                true,
                true,
                &attrs,
            );
        }

        let paths = TestPaths::new()?;
        builder.set_template_coordinate_sort_order().write(&paths.input)?;

        let mut cmd = create_simplex_with_paths(paths.input.clone(), paths.output.clone());
        cmd.read_group.read_group_id = "ABC".to_string();

        cmd.execute("test")?;

        // Verify output
        let records = read_bam_records(&paths.output)?;

        // Should have 2000 consensus reads (1000 UMI groups × 2 reads per pair)
        assert_eq!(records.len(), 2000, "Should have 2000 consensus reads");

        // Verify first and second of pair counts
        let first_count = records.iter().filter(|r| r.flags().is_first_segment()).count();
        let second_count = records.iter().filter(|r| r.flags().is_last_segment()).count();
        assert_eq!(first_count, 1000, "Should have 1000 first-of-pair reads");
        assert_eq!(second_count, 1000, "Should have 1000 second-of-pair reads");

        // Verify all reads have expected attributes
        for record in &records {
            // Check sequence
            assert_eq!(record.sequence().len(), 100, "Sequence length should be 100");

            // Check read group
            let rg = get_string_tag(record, "RG");
            assert_eq!(rg.as_deref(), Some("ABC"), "Read group should be ABC");

            // Check consensus tags are present
            let cd_tag = get_int_tag(record, "cD");
            assert!(cd_tag.is_some(), "cD tag should be present");
            assert_eq!(
                cd_tag.expect("cD tag should have a value"),
                2,
                "Depth should be 2 (2 reads per UMI)"
            );

            // MI and RX tags should be preserved by the consensus caller
            let mi_tag = get_string_tag(record, "MI");
            assert!(mi_tag.is_some(), "MI tag should be preserved");
        }

        Ok(())
    }

    #[rstest]
    #[case::fast_path(ThreadingOptions::none())]
    #[case::threaded(ThreadingOptions::new(2))]
    fn test_end_to_end_single_end_workflow(#[case] threading: ThreadingOptions) -> Result<()> {
        // Verify CB partitioning: two groups share the same MI but differ in CB.
        // They must produce two separate consensus reads, one per cell barcode.
        let mut builder = SamBuilder::new_unmapped();

        let mut attrs_a = HashMap::new();
        attrs_a.insert("RX", BufValue::from("ACGT"));
        attrs_a.insert("MI", BufValue::from("shared"));
        attrs_a.insert("CB", BufValue::from("AB"));

        // 3 fragments from cell AB
        builder.add_frag_with_attrs("a1", None, true, &attrs_a);
        builder.add_frag_with_attrs("a2", None, true, &attrs_a);
        builder.add_frag_with_attrs("a3", None, true, &attrs_a);

        let mut attrs_b = HashMap::new();
        attrs_b.insert("RX", BufValue::from("ACGT"));
        attrs_b.insert("MI", BufValue::from("shared"));
        attrs_b.insert("CB", BufValue::from("CD"));

        // 2 fragments from cell CD — same MI, different cell
        builder.add_frag_with_attrs("b1", None, true, &attrs_b);
        builder.add_frag_with_attrs("b2", None, true, &attrs_b);

        let paths = TestPaths::new()?;
        builder.set_template_coordinate_sort_order().write(&paths.input)?;

        let mut cmd = create_simplex_with_paths(paths.input.clone(), paths.output.clone());
        cmd.read_group.read_group_id = "ABC".to_string();
        cmd.threading = threading;

        cmd.execute("test")?;

        // Verify output
        let records = read_bam_records(&paths.output)?;

        // Should have 2 consensus reads — one per cell barcode
        assert_eq!(records.len(), 2, "Should have 2 consensus reads (one per cell barcode)");

        // All should be unpaired
        assert!(records.iter().all(|r| !r.flags().is_segmented()), "All reads should be unpaired");

        // Both cell barcodes must survive independently
        let observed_cbs: std::collections::HashSet<_> =
            records.iter().filter_map(|r| get_string_tag(r, "CB")).collect();
        assert_eq!(
            observed_cbs,
            ["AB".to_string(), "CD".to_string()].into_iter().collect(),
            "Both cell barcodes should produce independent consensus reads"
        );

        // Verify per-record attributes
        for record in &records {
            let rg = get_string_tag(record, "RG");
            assert_eq!(rg.as_deref(), Some("ABC"), "Read group should be ABC");

            let cd_tag = get_int_tag(record, "cD");
            assert!(cd_tag.is_some(), "cD tag should be present");
            assert!(cd_tag.expect("cD tag should have a value") >= 2, "Depth should be at least 2");

            assert_eq!(record.sequence().len(), 100, "Sequence length should be 100");
        }

        Ok(())
    }

    /// `--allow-unmapped` gates the fgbio `ConsensusCallingIterator` pre-group
    /// filter. By default (fgbio parity) unmapped, unpaired reads are dropped
    /// before consensus calling, so no consensus is produced; with
    /// `--allow-unmapped` they are consensus-called. Runs both the
    /// single-worker chain (absent `--threads`) and the multi-worker chain,
    /// which install the filter at the same `add_simplex` site — worker count
    /// is the only difference between the two cases.
    #[rstest]
    #[case::no_threads(ThreadingOptions::none())]
    #[case::threaded(ThreadingOptions::new(2))]
    fn test_allow_unmapped_gates_pregroup_filter(
        #[case] threading: ThreadingOptions,
    ) -> Result<()> {
        let mut builder = SamBuilder::new_unmapped();
        let mut attrs = HashMap::new();
        attrs.insert("RX", BufValue::from("ACGT"));
        attrs.insert("MI", BufValue::from("m1"));
        attrs.insert("CB", BufValue::from("AB"));
        // Two unmapped, unpaired fragments in one MI group.
        builder.add_frag_with_attrs("r1", None, true, &attrs);
        builder.add_frag_with_attrs("r2", None, true, &attrs);
        // Advertise template-coordinate order so the input clears the consensus
        // sort-order guard (`check_consensus_sort_order`); this test exercises the
        // pre-group filter, not the sort guard.
        builder.set_template_coordinate_sort_order();

        // Default (allow_unmapped = false): the pre-group filter drops
        // unmapped-unpaired reads -> no consensus records. (The shared test
        // helper enables the flag for its unmapped fixtures, so set it back to
        // the CLI default here explicitly.)
        let default_paths = TestPaths::new()?;
        builder.write(&default_paths.input)?;
        let mut default_cmd =
            create_simplex_with_paths(default_paths.input.clone(), default_paths.output.clone());
        default_cmd.threading = threading.clone();
        default_cmd.allow_unmapped.enabled = false;
        default_cmd.execute("test")?;
        let default_records = read_bam_records(&default_paths.output)?;
        assert_eq!(
            default_records.len(),
            0,
            "default (fgbio parity) drops unmapped, unpaired reads before consensus"
        );

        // --allow-unmapped: the unmapped reads are consensus-called (one
        // consensus read for the single MI group).
        let allow_paths = TestPaths::new()?;
        builder.write(&allow_paths.input)?;
        let mut allow_cmd =
            create_simplex_with_paths(allow_paths.input.clone(), allow_paths.output.clone());
        allow_cmd.threading = threading;
        allow_cmd.allow_unmapped.enabled = true;
        allow_cmd.execute("test")?;
        let allow_records = read_bam_records(&allow_paths.output)?;
        assert_eq!(
            allow_records.len(),
            1,
            "--allow-unmapped consensus-calls unmapped reads (one consensus per MI group)"
        );

        Ok(())
    }

    /// The `--allow-unmapped` CLI flag defaults to `false` (fgbio parity: the
    /// pre-group filter is on). Pinned separately from the behavioral tests
    /// because the shared test helper enables the flag for its unmapped
    /// fixtures.
    #[test]
    fn test_allow_unmapped_defaults_to_false() {
        let default_cmd = <Simplex as clap::Parser>::try_parse_from([
            "simplex", "-i", "in.bam", "-o", "out.bam", "-M", "1",
        ])
        .expect("simplex should parse with only the required args");
        assert!(!default_cmd.allow_unmapped.enabled, "--allow-unmapped must default to false");

        let enabled_cmd = <Simplex as clap::Parser>::try_parse_from([
            "simplex",
            "-i",
            "in.bam",
            "-o",
            "out.bam",
            "-M",
            "1",
            "--allow-unmapped",
        ])
        .expect("simplex should parse with --allow-unmapped");
        assert!(enabled_cmd.allow_unmapped.enabled, "bare --allow-unmapped must enable the flag");
    }

    #[test]
    fn test_min_reads_filtering() -> Result<()> {
        // Test that groups with fewer than min_reads are rejected
        let mut builder = SamBuilder::new_unmapped();

        // Group 1: 3 reads (should pass with min_reads=3)
        let mut attrs1 = HashMap::new();
        attrs1.insert("MI", BufValue::from("group1"));
        builder.add_frag_with_attrs("a1", None, true, &attrs1);
        builder.add_frag_with_attrs("a2", None, true, &attrs1);
        builder.add_frag_with_attrs("a3", None, true, &attrs1);

        // Group 2: 2 reads (should be rejected with min_reads=3)
        let mut attrs2 = HashMap::new();
        attrs2.insert("MI", BufValue::from("group2"));
        builder.add_frag_with_attrs("b1", None, true, &attrs2);
        builder.add_frag_with_attrs("b2", None, true, &attrs2);

        // Group 3: 1 read (should be rejected with min_reads=3)
        let mut attrs3 = HashMap::new();
        attrs3.insert("MI", BufValue::from("group3"));
        builder.add_frag_with_attrs("c1", None, true, &attrs3);

        let paths = TestPaths::new()?;
        builder.set_template_coordinate_sort_order().write(&paths.input)?;

        let mut cmd = create_simplex_with_paths(paths.input.clone(), paths.output.clone());
        cmd.min_reads = 3;

        cmd.execute("test")?;

        // Verify output
        let records = read_bam_records(&paths.output)?;

        // Should have only 1 consensus read (from group1)
        assert_eq!(records.len(), 1, "Should have 1 consensus read");

        // Check consensus depth tag is present (validates we got the group with 3 reads)
        let cd_tag = get_int_tag(&records[0], "cD");
        assert_eq!(
            cd_tag.expect("cD tag should have a value"),
            3,
            "Depth should be 3 for the passing group"
        );

        Ok(())
    }

    #[test]
    fn test_per_base_tags_generation() -> Result<()> {
        // Test that per-base tags (cd, ce) are generated when requested
        let mut builder = SamBuilder::new_unmapped();

        let mut attrs = HashMap::new();
        attrs.insert("MI", BufValue::from("test_umi"));

        // Add 5 reads to get sufficient depth
        for i in 0..5 {
            builder.add_frag_with_attrs(&format!("read{i}"), None, true, &attrs);
        }

        let paths = TestPaths::new()?;
        builder.set_template_coordinate_sort_order().write(&paths.input)?;

        let mut cmd = create_simplex_with_paths(paths.input.clone(), paths.output.clone());
        cmd.consensus.output_per_base_tags = true;

        cmd.execute("test")?;

        // Verify output
        let records = read_bam_records(&paths.output)?;
        assert_eq!(records.len(), 1, "Should have 1 consensus read");

        let record = &records[0];

        // Check per-read tags (cD, cM, cE) are present
        let cd_tag = get_int_tag(record, "cD");
        assert!(cd_tag.is_some(), "cD tag should be present");
        assert_eq!(cd_tag.expect("cD tag should have a value"), 5, "Max depth should be 5");

        let cm_tag = get_int_tag(record, "cM");
        assert!(cm_tag.is_some(), "cM tag should be present");
        assert_eq!(cm_tag.expect("cM tag should have a value"), 5, "Min depth should be 5");

        // Check per-base tags (cd, ce) are present
        let tag_bytes = "cd".as_bytes();
        let cd_array_tag = Tag::from([tag_bytes[0], tag_bytes[1]]);
        assert!(record.data().get(&cd_array_tag).is_some(), "cd per-base tag should be present");

        let tag_bytes = "ce".as_bytes();
        let ce_array_tag = Tag::from([tag_bytes[0], tag_bytes[1]]);
        assert!(record.data().get(&ce_array_tag).is_some(), "ce per-base tag should be present");

        Ok(())
    }

    #[test]
    fn test_multithreading() -> Result<()> {
        // Test that multithreading produces same results as single-threaded
        let mut builder = SamBuilder::new_unmapped();

        // Create 100 UMI groups with 2 reads each
        for idx in 0..100 {
            let mut attrs = HashMap::new();
            attrs.insert("MI", BufValue::from(format!("umi_{idx}")));
            builder.add_frag_with_attrs(&format!("read_{idx}a"), None, true, &attrs);
            builder.add_frag_with_attrs(&format!("read_{idx}b"), None, true, &attrs);
        }

        let paths = TestPaths::new()?;
        let output_multi_path = paths.dir.path().join("output_multi.bam");
        builder.set_template_coordinate_sort_order().write(&paths.input)?;

        // Run with single thread
        let cmd_single = create_simplex_with_paths(paths.input.clone(), paths.output.clone());
        cmd_single.execute("test")?;

        // Run with multiple threads
        let mut cmd_multi =
            create_simplex_with_paths(paths.input.clone(), output_multi_path.clone());
        cmd_multi.threading = ThreadingOptions::new(4);
        cmd_multi.execute("test")?;

        // Verify both outputs have same number of records
        let records_single = read_bam_records(&paths.output)?;
        let records_multi = read_bam_records(&output_multi_path)?;

        assert_eq!(
            records_single.len(),
            records_multi.len(),
            "Single and multi-threaded should produce same number of records"
        );
        assert_eq!(records_single.len(), 100, "Should have 100 consensus reads");

        Ok(())
    }

    #[test]
    fn test_max_reads_downsampling() -> Result<()> {
        // Test that max_reads properly downsamples large UMI groups
        let mut builder = SamBuilder::new_unmapped();

        let mut attrs = HashMap::new();
        attrs.insert("MI", BufValue::from("large_group"));

        // Add 20 reads to this group
        for i in 0..20 {
            builder.add_frag_with_attrs(&format!("read{i}"), None, true, &attrs);
        }

        let paths = TestPaths::new()?;
        builder.set_template_coordinate_sort_order().write(&paths.input)?;

        let mut cmd = create_simplex_with_paths(paths.input.clone(), paths.output.clone());
        cmd.max_reads = Some(10);

        cmd.execute("test")?;

        // Verify output
        let records = read_bam_records(&paths.output)?;
        assert_eq!(records.len(), 1, "Should have 1 consensus read");

        // Check depth tag shows downsampling was applied (max depth should be 10)
        let cd_tag = get_int_tag(&records[0], "cD");
        assert!(cd_tag.is_some(), "cD tag should be present");
        assert!(
            cd_tag.expect("cD tag should have a value") <= 10,
            "Max depth should be <= 10 due to downsampling"
        );

        Ok(())
    }

    #[test]
    fn test_max_reads_less_than_min_reads_fails() {
        // Create a dummy BAM file so validation gets to tag check
        let mut builder = SamBuilder::new_unmapped();
        let paths = TestPaths::new().expect("failed to create test paths");
        builder
            .set_template_coordinate_sort_order()
            .write(&paths.input)
            .expect("failed to write test BAM");

        let mut cmd = create_simplex_with_paths(paths.input.clone(), PathBuf::from("out.bam"));
        cmd.min_reads = 5;
        cmd.max_reads = Some(3); // Invalid: less than min_reads

        let result = cmd.execute("test");
        assert!(result.is_err());
        let error_msg = result.unwrap_err().to_string();
        assert!(error_msg.contains("--max-reads"));
        assert!(error_msg.contains("--min-reads"));
    }

    #[test]
    fn test_statistics_file_generation() -> Result<()> {
        let mut builder = SamBuilder::new_unmapped();

        // Group 1: passes (3 reads)
        let mut attrs1 = HashMap::new();
        attrs1.insert("MI", BufValue::from("pass"));
        builder.add_frag_with_attrs("p1", None, true, &attrs1);
        builder.add_frag_with_attrs("p2", None, true, &attrs1);
        builder.add_frag_with_attrs("p3", None, true, &attrs1);

        // Group 2: filtered (1 read, min_reads=2)
        let mut attrs2 = HashMap::new();
        attrs2.insert("MI", BufValue::from("fail"));
        builder.add_frag_with_attrs("f1", None, true, &attrs2);

        let paths = TestPaths::new()?;
        builder.set_template_coordinate_sort_order().write(&paths.input)?;

        let mut cmd = create_simplex_with_paths(paths.input.clone(), paths.output.clone());
        cmd.stats_opts.stats = Some(paths.stats.clone());
        cmd.min_reads = 2;

        cmd.execute("test")?;

        // Verify stats file was created and contains data
        assert!(&paths.stats.exists(), "Stats file should exist");

        // Read and verify TSV format (now vertical key-value-description format)
        let kv_metrics: Vec<ConsensusKvMetric> =
            DelimFile::default().read_tsv(&paths.stats).expect("Failed to read metrics file");

        // Should have multiple rows (one per metric)
        assert!(!kv_metrics.is_empty(), "Should have metrics");

        // Verify expected keys are present
        let keys: Vec<&str> = kv_metrics.iter().map(|m| m.key.as_str()).collect();
        assert!(keys.contains(&"raw_reads_considered"), "Should have raw_reads_considered");
        assert!(keys.contains(&"raw_reads_rejected"), "Should have raw_reads_rejected");
        assert!(keys.contains(&"consensus_reads_emitted"), "Should have consensus_reads_emitted");

        // Verify raw_reads_considered has a value
        let raw_reads = kv_metrics
            .iter()
            .find(|m| m.key == "raw_reads_considered")
            .expect("Should have raw_reads_considered");
        let count: u64 = raw_reads.value.parse().expect("Should be a number");
        assert!(count > 0, "Stats should have total reads");

        Ok(())
    }

    #[test]
    fn test_custom_read_name_prefix() -> Result<()> {
        let mut builder = SamBuilder::new_unmapped();

        let mut attrs = HashMap::new();
        attrs.insert("MI", BufValue::from("test"));
        builder.add_frag_with_attrs("read1", None, true, &attrs);
        builder.add_frag_with_attrs("read2", None, true, &attrs);

        let paths = TestPaths::new()?;
        builder.set_template_coordinate_sort_order().write(&paths.input)?;

        let mut cmd = create_simplex_with_paths(paths.input.clone(), paths.output.clone());
        cmd.read_group.read_name_prefix = Some("MYCONSENSUS".to_string());

        cmd.execute("test")?;

        // Verify output read name has custom prefix
        let records = read_bam_records(&paths.output)?;
        assert_eq!(records.len(), 1, "Should have 1 consensus read");

        let read_name = records[0].name().map(std::string::ToString::to_string).unwrap_or_default();
        assert!(read_name.starts_with("MYCONSENSUS"), "Read name should start with custom prefix");

        Ok(())
    }

    #[test]
    fn test_rejects_file_generation() -> Result<()> {
        // Test that rejected reads are written to rejects file
        let mut builder = SamBuilder::new_unmapped();

        // Group 1: passes (3 reads)
        let mut attrs1 = HashMap::new();
        attrs1.insert("MI", BufValue::from("pass"));
        builder.add_frag_with_attrs("p1", None, true, &attrs1);
        builder.add_frag_with_attrs("p2", None, true, &attrs1);
        builder.add_frag_with_attrs("p3", None, true, &attrs1);

        // Group 2: filtered (1 read, min_reads=2)
        let mut attrs2 = HashMap::new();
        attrs2.insert("MI", BufValue::from("fail"));
        builder.add_frag_with_attrs("f1", None, true, &attrs2);

        let paths = TestPaths::new()?;
        builder.set_template_coordinate_sort_order().write(&paths.input)?;

        let mut cmd = create_simplex_with_paths(paths.input.clone(), paths.output.clone());
        cmd.rejects_opts.rejects = Some(paths.rejects.clone());
        cmd.min_reads = 2;

        cmd.execute("test")?;

        // Verify rejects file exists and contains rejected reads
        assert!(&paths.rejects.exists(), "Rejects file should exist");

        let reject_records = read_bam_records(&paths.rejects)?;
        assert_eq!(reject_records.len(), 1, "Should have 1 rejected read");

        // Verify the rejected read has the correct UMI
        let umi = get_string_tag(&reject_records[0], "MI");
        assert_eq!(umi.as_deref(), Some("fail"), "Rejected read should have 'fail' UMI");

        Ok(())
    }

    #[test]
    fn test_multithreaded_with_rejects() -> Result<()> {
        // Test multi-threaded execution with rejects tracking
        let mut builder = SamBuilder::new_unmapped();

        // Create 50 UMI groups
        for idx in 0..50 {
            let mut attrs = HashMap::new();
            if idx % 3 == 0 {
                // Every 3rd group has only 1 read (will be rejected with min_reads=2)
                attrs.insert("MI", BufValue::from(format!("umi_{idx}")));
                builder.add_frag_with_attrs(&format!("read_{idx}"), None, true, &attrs);
            } else {
                // Other groups have 2 reads (will pass)
                attrs.insert("MI", BufValue::from(format!("umi_{idx}")));
                builder.add_frag_with_attrs(&format!("read_{idx}a"), None, true, &attrs);
                builder.add_frag_with_attrs(&format!("read_{idx}b"), None, true, &attrs);
            }
        }

        let paths = TestPaths::new()?;
        builder.set_template_coordinate_sort_order().write(&paths.input)?;

        let mut cmd = create_simplex_with_paths(paths.input.clone(), paths.output.clone());
        cmd.rejects_opts.rejects = Some(paths.rejects.clone());
        cmd.min_reads = 2;
        cmd.threading = ThreadingOptions::new(4);

        cmd.execute("test")?;

        // Verify outputs
        let output_records = read_bam_records(&paths.output)?;
        let reject_records = read_bam_records(&paths.rejects)?;

        // About 33 groups should pass (those with 2 reads)
        assert!(output_records.len() >= 30, "Should have at least 30 consensus reads");

        // About 17 groups should be rejected (those with 1 read)
        assert!(reject_records.len() >= 15, "Should have at least 15 rejected reads");

        Ok(())
    }

    #[test]
    fn test_multithreaded_with_stats() -> Result<()> {
        // Test multi-threaded execution with statistics file generation
        let mut builder = SamBuilder::new_unmapped();

        for idx in 0..100 {
            let mut attrs = HashMap::new();
            attrs.insert("MI", BufValue::from(format!("umi_{idx}")));
            builder.add_frag_with_attrs(&format!("read_{idx}a"), None, true, &attrs);
            builder.add_frag_with_attrs(&format!("read_{idx}b"), None, true, &attrs);
        }

        let paths = TestPaths::new()?;
        builder.set_template_coordinate_sort_order().write(&paths.input)?;

        let mut cmd = create_simplex_with_paths(paths.input.clone(), paths.output.clone());
        cmd.stats_opts.stats = Some(paths.stats.clone());
        cmd.threading = ThreadingOptions::new(4);

        cmd.execute("test")?;

        // Verify stats file exists and has expected content
        assert!(&paths.stats.exists(), "Stats file should exist");

        // Read and verify TSV format (now vertical key-value-description format)
        let kv_metrics: Vec<ConsensusKvMetric> =
            DelimFile::default().read_tsv(&paths.stats).expect("Failed to read metrics file");

        // Should have multiple rows (one per metric)
        assert!(!kv_metrics.is_empty(), "Should have metrics");

        // Verify raw_reads_considered has a value
        let raw_reads = kv_metrics
            .iter()
            .find(|m| m.key == "raw_reads_considered")
            .expect("Should have raw_reads_considered");
        let count: u64 = raw_reads.value.parse().expect("Should be a number");
        assert!(count > 0, "Stats should have total reads");

        // Verify consensus_reads_emitted has a value
        let consensus = kv_metrics
            .iter()
            .find(|m| m.key == "consensus_reads_emitted")
            .expect("Should have consensus_reads_emitted");
        let count: u64 = consensus.value.parse().expect("Should be a number");
        assert!(count > 0, "Stats should have consensus count");

        Ok(())
    }

    #[test]
    fn test_overlapping_consensus_calling_paired() -> Result<()> {
        // Test overlapping consensus calling on paired-end reads
        let mut builder = SamBuilder::new_unmapped();

        let mut attrs = HashMap::new();
        attrs.insert("MI", BufValue::from("test_umi"));

        // Add 2 pairs with overlapping bases (both are unmapped in this simple test)
        builder.add_pair_with_attrs("pair1", None, None, true, true, &attrs);
        builder.add_pair_with_attrs("pair2", None, None, true, true, &attrs);

        let paths = TestPaths::new()?;
        builder.set_template_coordinate_sort_order().write(&paths.input)?;

        let mut cmd = create_simplex_with_paths(paths.input.clone(), paths.output.clone());
        cmd.overlapping.consensus_call_overlapping_bases = true;

        cmd.execute("test")?;

        // Verify output exists (overlapping consensus was called)
        let records = read_bam_records(&paths.output)?;
        assert_eq!(records.len(), 2, "Should have 2 consensus reads (R1 and R2)");

        Ok(())
    }

    #[test]
    fn test_overlapping_consensus_unpaired_read() -> Result<()> {
        // Test that unpaired reads in overlapping mode are handled correctly
        let mut builder = SamBuilder::new_unmapped();

        let mut attrs = HashMap::new();
        attrs.insert("MI", BufValue::from("test_umi"));

        // Add 3 unpaired fragments
        builder.add_frag_with_attrs("frag1", None, true, &attrs);
        builder.add_frag_with_attrs("frag2", None, true, &attrs);
        builder.add_frag_with_attrs("frag3", None, true, &attrs);

        let paths = TestPaths::new()?;
        builder.set_template_coordinate_sort_order().write(&paths.input)?;

        let mut cmd = create_simplex_with_paths(paths.input.clone(), paths.output.clone());
        cmd.overlapping.consensus_call_overlapping_bases = true;

        cmd.execute("test")?;

        // Verify unpaired reads are processed correctly
        let records = read_bam_records(&paths.output)?;
        assert_eq!(records.len(), 1, "Should have 1 consensus read from unpaired fragments");

        Ok(())
    }

    #[test]
    fn test_trim_enabled() -> Result<()> {
        // Test quality trimming option
        let mut builder = SamBuilder::new_unmapped();

        let mut attrs = HashMap::new();
        attrs.insert("MI", BufValue::from("test_umi"));
        builder.add_frag_with_attrs("read1", None, true, &attrs);
        builder.add_frag_with_attrs("read2", None, true, &attrs);

        let paths = TestPaths::new()?;
        builder.set_template_coordinate_sort_order().write(&paths.input)?;

        let mut cmd = create_simplex_with_paths(paths.input.clone(), paths.output.clone());
        cmd.consensus.trim = true;

        cmd.execute("test")?;

        // Verify output exists with trimming enabled
        let records = read_bam_records(&paths.output)?;
        assert_eq!(records.len(), 1, "Should have 1 consensus read");

        Ok(())
    }

    #[test]
    fn test_min_consensus_base_quality() -> Result<()> {
        // Test minimum consensus base quality masking
        let mut builder = SamBuilder::new_unmapped();

        let mut attrs = HashMap::new();
        attrs.insert("MI", BufValue::from("test_umi"));
        builder.add_frag_with_attrs("read1", None, true, &attrs);
        builder.add_frag_with_attrs("read2", None, true, &attrs);

        let paths = TestPaths::new()?;
        builder.set_template_coordinate_sort_order().write(&paths.input)?;

        let mut cmd = create_simplex_with_paths(paths.input.clone(), paths.output.clone());
        cmd.consensus.min_consensus_base_quality = 30;

        cmd.execute("test")?;

        // Verify output exists with quality filtering
        let records = read_bam_records(&paths.output)?;
        assert_eq!(records.len(), 1, "Should have 1 consensus read");

        Ok(())
    }

    #[test]
    fn test_reads_without_umi_tag_skipped() -> Result<()> {
        // Test that reads without the UMI tag are skipped
        let mut builder = SamBuilder::new_unmapped();

        let mut attrs_with_umi = HashMap::new();
        attrs_with_umi.insert("MI", BufValue::from("has_umi"));
        builder.add_frag_with_attrs("with_umi1", None, true, &attrs_with_umi);
        builder.add_frag_with_attrs("with_umi2", None, true, &attrs_with_umi);

        // Add reads without UMI tag
        let attrs_no_umi = HashMap::new();
        builder.add_frag_with_attrs("no_umi1", None, true, &attrs_no_umi);
        builder.add_frag_with_attrs("no_umi2", None, true, &attrs_no_umi);

        let paths = TestPaths::new()?;
        builder.set_template_coordinate_sort_order().write(&paths.input)?;

        let cmd = create_simplex_with_paths(paths.input.clone(), paths.output.clone());

        cmd.execute("test")?;

        // Verify only reads with UMI tag generated consensus
        let records = read_bam_records(&paths.output)?;
        assert_eq!(records.len(), 1, "Should have 1 consensus read from reads with UMI only");

        Ok(())
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
        let paths = TestPaths::new()?;

        let mut builder = SamBuilder::with_single_ref("chr1", 100);
        let mut attrs = HashMap::new();
        attrs.insert("MI", BufValue::from("1"));
        // Add a few reads with the same MI tag
        builder.add_frag_with_attrs("read1", None, true, &attrs);
        builder.add_frag_with_attrs("read2", None, true, &attrs);
        builder.add_frag_with_attrs("read3", None, true, &attrs);
        builder.set_template_coordinate_sort_order().write(&paths.input)?;

        let mut cmd = create_simplex_with_paths(paths.input.clone(), paths.output.clone());
        cmd.threading = threading;
        cmd.execute("test")?;

        // Should produce consensus output
        let records = read_bam_records(&paths.output)?;
        assert_eq!(records.len(), 1, "Should have 1 consensus read");

        Ok(())
    }

    #[rstest]
    #[case::single_threaded(ThreadingOptions::none())]
    #[case::multi_threaded(ThreadingOptions::new(2))]
    fn test_simplex_em_seq_command(#[case] threading: ThreadingOptions) -> Result<()> {
        use crate::sam::builder::{Strand, create_test_fasta};

        // Create a FASTA with chr1 containing C bases at positions 0..20
        let ref_seq = "C".repeat(200);
        let ref_fasta = create_test_fasta(&[("chr1", &ref_seq)])?;

        // Create mapped reads with MI tag showing C bases (methylated) at ref-C positions
        let mut builder = SamBuilder::with_single_ref("chr1", 200);
        let mut attrs = HashMap::new();
        attrs.insert("MI", BufValue::from("1"));

        // Add 3 fragment reads at position 1 (1-based), all showing C (methylated)
        for i in 0..3 {
            let _ = builder
                .add_frag()
                .name(&format!("r{i}"))
                .start(1)
                .strand(Strand::Plus)
                .bases("CCCCCCCCCC")
                .attr("MI", BufValue::from("1"))
                .build();
        }

        let paths = TestPaths::new()?;
        builder.set_template_coordinate_sort_order().write(&paths.input)?;

        let mut cmd = create_simplex_with_paths(paths.input.clone(), paths.output.clone());
        cmd.methylation_mode = Some(crate::commands::common::MethylationModeArg::EmSeq);
        cmd.reference = Some(ref_fasta.path().to_path_buf());
        cmd.threading = threading;
        cmd.execute("test")?;

        let records = read_bam_records(&paths.output)?;
        assert_eq!(records.len(), 1, "Should have 1 consensus read");

        // Verify methylation tags are present and correct
        use noodles::sam::alignment::record_buf::data::field::value::{
            Array as BufArray, Value as BufValue,
        };
        let record = &records[0];
        let cu_tag = Tag::from([b'c', b'u']);
        let cu_value = record.data().get(&cu_tag).expect("cu tag should be present with EM-Seq");
        // All 3 reads show C at ref-C positions → unconverted counts should be non-zero
        if let BufValue::Array(BufArray::Int16(cu_vals)) = cu_value {
            assert!(
                cu_vals.iter().any(|&v| v > 0),
                "cu should have non-zero values for methylated reads"
            );
        } else {
            panic!("cu tag should be an i16 array");
        }

        let ct_tag = Tag::from([b'c', b't']);
        let ct_value = record.data().get(&ct_tag).expect("ct tag should be present with EM-Seq");
        // All reads show C (not T) → converted counts should be 0
        if let BufValue::Array(BufArray::Int16(ct_vals)) = ct_value {
            assert!(
                ct_vals.iter().all(|&v| v == 0),
                "ct should be all zeros for fully methylated reads"
            );
        } else {
            panic!("ct tag should be an i16 array");
        }

        Ok(())
    }
}
