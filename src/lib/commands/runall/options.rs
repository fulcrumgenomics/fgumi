//! Enums and option structs for pipeline CLI arguments.

use std::path::PathBuf;

use crate::aligner::AlignerPreset;
use crate::commands::common::parse_bool;
use crate::defaults;
use clap::{Args, ValueEnum};
use fgumi_derive::multi_options;
use fgumi_umi::Strategy;

/// Pipeline entry point: which stage to start from.
#[derive(Debug, Clone, Copy, ValueEnum)]
pub enum StartFrom {
    /// Start from grouped BAM (MI tags present). Runs: consensus -> filter -> write.
    Group,
    /// Start from coordinate-sorted BAM. Runs: sort -> group -> consensus -> filter -> write.
    Sort,
    /// Start from unmapped BAM with UMIs. Runs: correct -> fastq -> aligner -> zipper -> sort -> ...
    Correct,
    /// Start from unmapped BAM. Runs: fastq -> aligner -> zipper -> sort -> ...
    Fastq,
    /// Start from FASTQ. Runs: extract -> fastq -> aligner -> zipper -> sort -> ... (full pipeline).
    Extract,
}

/// Consensus algorithm to use in the pipeline.
#[derive(Debug, Clone, Copy, ValueEnum)]
pub enum ConsensusMode {
    /// Simplex (vanilla UMI) consensus calling.
    Simplex,
    /// Duplex consensus calling (requires paired UMI grouping with /A and /B suffixes).
    Duplex,
    /// CODEC consensus calling (R1/R2 from same molecule sequence opposite strands).
    Codec,
}

/// Aligner preset for `--aligner::preset` CLI option.
///
/// Each variant maps to an [`AlignerPreset`] that builds the aligner command string
/// and validates that the binary and index files are present.
#[derive(clap::ValueEnum, Clone, Debug)]
pub enum AlignerPresetArg {
    /// BWA-MEM aligner (`bwa mem`).
    BwaMem,
    /// BWA-MEM2 aligner (`bwa-mem2 mem`).
    BwaMem2,
}

impl AlignerPresetArg {
    /// Convert to the corresponding [`AlignerPreset`].
    pub(crate) fn to_preset(&self) -> AlignerPreset {
        match self {
            AlignerPresetArg::BwaMem => AlignerPreset::BwaMem,
            AlignerPresetArg::BwaMem2 => AlignerPreset::BwaMem2,
        }
    }
}

/// Pipeline stage to stop after (inclusive).
///
/// Controls how far the pipeline runs before writing output and exiting.
/// The value must be at or after the `--start-from` stage.
#[derive(clap::ValueEnum, Clone, Copy, Debug, PartialEq, Eq, PartialOrd, Ord)]
pub enum StopAfter {
    /// Stop after extraction; output is an unmapped BAM with UMI tags.
    Extract,
    /// Stop after UMI correction; output is a corrected unmapped BAM.
    Correct,
    /// Stop after FASTQ conversion; writes interleaved FASTQ to the output path.
    /// Valid with --start-from extract, correct, or fastq.
    Fastq,
    /// Stop after alignment; output is the raw aligner SAM output.
    Align,
    /// Stop after zipper merge (--start-from zipper or extract).
    Zipper,
    /// Stop after template-coordinate sort; output is a sorted BAM.
    Sort,
    /// Stop after UMI grouping; output is a sorted BAM with MI tags.
    Group,
    /// Stop after consensus calling; output is consensus reads (no min-reads filter).
    Consensus,
    /// Stop after filtering (default); output is filtered consensus reads.
    Filter,
}

// ============================================================================
// Option structs with multi_options for scoped CLI naming
// ============================================================================

/// Options for the external sort stage.
#[multi_options("sort", "Sort Options")]
#[derive(Args, Debug, Clone)]
pub struct SortOptions {
    /// Memory limit in bytes for external sort.
    #[arg(long, default_value_t = defaults::SORT_MEMORY_LIMIT)]
    pub memory_limit: usize,
    /// Temporary directory for sort spill files (default: system temp).
    #[arg(long)]
    pub temp_dir: Option<PathBuf>,
}

impl Default for SortOptions {
    fn default() -> Self {
        Self { memory_limit: defaults::SORT_MEMORY_LIMIT, temp_dir: None }
    }
}

/// Options for UMI grouping.
#[multi_options("group", "Group Options")]
#[derive(Args, Debug, Clone)]
pub struct GroupOptions {
    /// UMI assignment strategy.
    #[arg(long, default_value_t = Strategy::Adjacency, value_enum)]
    pub strategy: Strategy,
    /// Maximum allowed UMI mismatches for grouping.
    #[arg(long, default_value_t = defaults::GROUP_MAX_EDITS)]
    pub max_edits: u32,
    /// Minimum mapping quality for grouping. Templates where any mapped read has
    /// MAPQ below this threshold are discarded before UMI assignment.
    #[arg(long, default_value_t = defaults::GROUP_MIN_MAPQ)]
    pub min_mapq: u8,
    /// Allow fully unmapped templates (both reads unmapped) through grouping.
    #[arg(long, default_value_t = defaults::GROUP_ALLOW_UNMAPPED, num_args = 0..=1,
          default_missing_value = "true", action = clap::ArgAction::Set,
          value_parser = parse_bool)]
    pub allow_unmapped: bool,
}

impl Default for GroupOptions {
    fn default() -> Self {
        Self {
            strategy: Strategy::Adjacency,
            max_edits: defaults::GROUP_MAX_EDITS,
            min_mapq: defaults::GROUP_MIN_MAPQ,
            allow_unmapped: defaults::GROUP_ALLOW_UNMAPPED,
        }
    }
}

/// Options for consensus read filtering.
///
/// Defaults mirror the standalone `fgumi filter` command by pulling from the
/// same `crate::defaults::FILTER_*` constants, so running `fgumi filter`
/// standalone and `fgumi runall` with otherwise-identical CLI args produces
/// the same output.
#[multi_options("filter", "Filter Options")]
#[derive(Args, Debug, Clone)]
pub struct FilterOptions {
    /// Minimum number of reads to produce a consensus read. Accepts 1-3
    /// comma-separated values for duplex: `[total, AB, BA]`.
    ///
    /// No default: matches standalone `fgumi filter` which also requires an
    /// explicit `--min-reads`. Required when the plan actually runs the filter
    /// stage (see `Runall::validate`); optional otherwise (e.g. plans that
    /// stop before filter).
    #[arg(long, value_delimiter = ',')]
    pub min_reads: Vec<usize>,
    /// Minimum input base quality for consensus calling.
    #[arg(long, default_value_t = defaults::FILTER_MIN_BASE_QUALITY)]
    pub min_base_quality: u8,
    /// Maximum raw read error rate for a single-strand consensus base/read.
    /// Accepts 1-3 comma-separated values for duplex: `[total, AB, BA]`.
    #[arg(long, value_delimiter = ',', default_value = "0.025")]
    pub max_read_error_rate: Vec<f64>,
    /// Maximum base error rate across raw reads. Accepts 1-3 comma-separated
    /// values for duplex: `[total, AB, BA]`.
    #[arg(long, value_delimiter = ',', default_value = "0.1")]
    pub max_base_error_rate: Vec<f64>,
    /// Maximum fraction of no-call (N) bases per read.
    #[arg(long, default_value_t = defaults::FILTER_MAX_NO_CALL_FRACTION)]
    pub max_no_call_fraction: f64,
    /// Minimum mean base quality across the read (after masking). If unset,
    /// no mean-quality filter is applied.
    #[arg(long)]
    pub min_mean_base_quality: Option<f64>,
    /// Drop all records of a template when any primary record fails the filter.
    /// Matches the `--filter-by-template` flag on standalone `fgumi filter`
    /// (default `true`).
    #[arg(long, default_value = "true", num_args = 0..=1,
          default_missing_value = "true", action = clap::ArgAction::Set,
          value_parser = parse_bool)]
    pub filter_by_template: bool,
    /// Require single-strand agreement when masking duplex consensus bases.
    /// Only applies to duplex and codec consensus modes.
    #[arg(long, default_value = "false", num_args = 0..=1,
          default_missing_value = "true", action = clap::ArgAction::Set,
          value_parser = parse_bool)]
    pub require_ss_agreement: bool,
    /// Optional output path for filtering statistics TSV. Matches standalone
    /// `fgumi filter --stats`: a four-row TSV with `total_reads`,
    /// `passed_reads`, `failed_reads`, `pass_rate`.
    #[arg(long)]
    pub stats: Option<PathBuf>,
}

impl Default for FilterOptions {
    fn default() -> Self {
        Self {
            min_reads: Vec::new(),
            min_base_quality: defaults::FILTER_MIN_BASE_QUALITY,
            max_read_error_rate: vec![defaults::FILTER_MAX_READ_ERROR_RATE],
            max_base_error_rate: vec![defaults::FILTER_MAX_BASE_ERROR_RATE],
            max_no_call_fraction: defaults::FILTER_MAX_NO_CALL_FRACTION,
            min_mean_base_quality: None,
            filter_by_template: true,
            require_ss_agreement: false,
            stats: None,
        }
    }
}

/// Options for consensus calling.
///
/// Defaults mirror the standalone consensus commands' [`ConsensusCallingOptions`]
/// by pulling from the same `crate::defaults::*` constants, so running
/// `fgumi simplex` back-to-back with the standalone commands produces the same
/// output as `fgumi runall`.
///
/// [`ConsensusCallingOptions`]: crate::commands::common::ConsensusCallingOptions
#[multi_options("consensus", "Consensus Options")]
#[derive(Args, Debug, Clone)]
pub struct ConsensusOptions {
    /// Pre-UMI error rate in Phred scale.
    #[arg(long, default_value_t = defaults::CONSENSUS_ERROR_RATE_PRE_UMI)]
    pub error_rate_pre_umi: u8,
    /// Post-UMI error rate in Phred scale.
    #[arg(long, default_value_t = defaults::CONSENSUS_ERROR_RATE_POST_UMI)]
    pub error_rate_post_umi: u8,
    /// Minimum base quality in raw reads to use for consensus.
    #[arg(long, default_value_t = defaults::CONSENSUS_MIN_INPUT_BASE_QUALITY)]
    pub min_input_base_quality: u8,
    /// Minimum consensus base quality (output consensus bases below this are masked to N).
    #[arg(long, default_value_t = defaults::CONSENSUS_MIN_BASE_QUALITY)]
    pub min_consensus_base_quality: u8,
    /// Produce per-base tags (cd, ce) in addition to per-read tags.
    #[arg(long, default_value_t = defaults::CONSENSUS_OUTPUT_PER_BASE_TAGS,
          num_args = 0..=1, default_missing_value = "true",
          action = clap::ArgAction::Set, value_parser = parse_bool)]
    pub output_per_base_tags: bool,
    /// Quality-trim reads before consensus calling.
    #[arg(long, default_value_t = defaults::CONSENSUS_TRIM,
          num_args = 0..=1, default_missing_value = "true",
          action = clap::ArgAction::Set, value_parser = parse_bool)]
    pub trim: bool,
    /// Apply overlapping-bases consensus correction to paired reads before
    /// consensus calling (matches standalone `fgumi simplex`/`duplex`
    /// `--consensus-call-overlapping-bases`; default `true`).
    #[arg(long, default_value_t = defaults::CONSENSUS_CALL_OVERLAPPING_BASES,
          num_args = 0..=1, default_missing_value = "true",
          action = clap::ArgAction::Set, value_parser = parse_bool)]
    pub call_overlapping_bases: bool,
    /// Read group ID for consensus reads.
    #[arg(long, default_value_t = String::from(defaults::READ_GROUP_ID))]
    pub read_group_id: String,
    /// Prefix for consensus read names (derived from header if not set).
    #[arg(long)]
    pub read_name_prefix: Option<String>,
}

impl Default for ConsensusOptions {
    fn default() -> Self {
        Self {
            error_rate_pre_umi: defaults::CONSENSUS_ERROR_RATE_PRE_UMI,
            error_rate_post_umi: defaults::CONSENSUS_ERROR_RATE_POST_UMI,
            min_input_base_quality: defaults::CONSENSUS_MIN_INPUT_BASE_QUALITY,
            min_consensus_base_quality: defaults::CONSENSUS_MIN_BASE_QUALITY,
            output_per_base_tags: defaults::CONSENSUS_OUTPUT_PER_BASE_TAGS,
            trim: defaults::CONSENSUS_TRIM,
            call_overlapping_bases: defaults::CONSENSUS_CALL_OVERLAPPING_BASES,
            read_group_id: String::from(defaults::READ_GROUP_ID),
            read_name_prefix: None,
        }
    }
}

/// Options for FASTQ extraction from input reads.
#[multi_options("extract", "Extract Options")]
#[derive(Args, Debug, Clone)]
pub struct ExtractOptions {
    /// Sample name (required for --start-from extract).
    #[arg(long)]
    pub sample: Option<String>,
    /// Library name (required for --start-from extract).
    #[arg(long)]
    pub library: Option<String>,
    /// Read structures for FASTQ extraction (one per input FASTQ, --start-from extract).
    #[arg(long, short = 'R', num_args = 1..)]
    pub read_structures: Option<Vec<String>>,
    /// Extract UMIs from field 8 of the Illumina read name (colon-delimited).
    #[arg(long, default_value = "false", num_args = 0..=1,
          default_missing_value = "true", action = clap::ArgAction::Set,
          value_parser = parse_bool)]
    pub extract_umis_from_read_names: bool,

    // ---- Read group metadata ----
    /// Sequencing platform.
    #[arg(long, default_value_t = String::from(defaults::PLATFORM))]
    pub platform: String,
    /// Platform unit (e.g. 'flowcell-barcode.lane.sample-barcode').
    #[arg(long)]
    pub platform_unit: Option<String>,
    /// Platform model (e.g. miseq, hiseq2500, hiseqX).
    #[arg(long)]
    pub platform_model: Option<String>,
    /// Library or sample barcode sequence.
    #[arg(long)]
    pub barcode: Option<String>,
    /// Sequencing center from which the data originated.
    #[arg(long)]
    pub sequencing_center: Option<String>,
    /// Predicted median insert size for the read group header.
    #[arg(long)]
    pub predicted_insert_size: Option<u32>,
    /// Description of the read group.
    #[arg(long)]
    pub description: Option<String>,
    /// Date the run was produced.
    #[arg(long)]
    pub run_date: Option<String>,
    /// Comment(s) to include in the output file's header.
    #[arg(long, num_args = 1..)]
    pub comment: Option<Vec<String>>,

    // ---- Extract tag/behavior options ----
    /// Single tag to store all concatenated UMIs (in addition to per-segment tags).
    #[arg(long)]
    pub single_tag: Option<String>,
    /// Annotate read names with UMIs (appends "+UMIs" to read names).
    #[arg(long, default_value = "false", num_args = 0..=1,
          default_missing_value = "true", action = clap::ArgAction::Set,
          value_parser = parse_bool)]
    pub annotate_read_names: bool,
    /// Store sample barcode quality scores in the QT tag.
    #[arg(long, default_value = "false", num_args = 0..=1,
          default_missing_value = "true", action = clap::ArgAction::Set,
          value_parser = parse_bool)]
    pub store_sample_barcode_qualities: bool,
}

impl Default for ExtractOptions {
    fn default() -> Self {
        Self {
            sample: None,
            library: None,
            read_structures: None,
            extract_umis_from_read_names: false,
            platform: String::from(defaults::PLATFORM),
            platform_unit: None,
            platform_model: None,
            barcode: None,
            sequencing_center: None,
            predicted_insert_size: None,
            description: None,
            run_date: None,
            comment: None,
            single_tag: None,
            annotate_read_names: false,
            store_sample_barcode_qualities: false,
        }
    }
}

/// Default aligner chunk size in bases for BWA `-K` flag.
const DEFAULT_ALIGNER_CHUNK_SIZE: u64 = defaults::BWA_CHUNK_SIZE;

/// Options for the aligner subprocess.
#[multi_options("aligner", "Aligner Options")]
#[derive(Args, Debug, Clone)]
pub struct AlignerOptions {
    /// Aligner preset (conflicts with `--aligner::command`).
    #[arg(long, conflicts_with = "aligner_command")]
    pub preset: Option<AlignerPresetArg>,
    /// Aligner command template. Use {ref} for reference path and {threads} for thread count.
    #[arg(long, conflicts_with = "aligner_preset")]
    pub command: Option<String>,
    /// Threads for the aligner subprocess (default: 4).
    #[arg(long, default_value_t = 4)]
    pub threads: usize,
    /// Processing chunk size in bases for the aligner (BWA `-K` flag).
    #[arg(long, default_value_t = DEFAULT_ALIGNER_CHUNK_SIZE)]
    pub chunk_size: u64,
}

impl Default for AlignerOptions {
    fn default() -> Self {
        Self { preset: None, command: None, threads: 4, chunk_size: DEFAULT_ALIGNER_CHUNK_SIZE }
    }
}

/// Options for the zipper merge stage.
#[multi_options("zipper", "Zipper Options")]
#[derive(Args, Debug, Clone, Default)]
pub struct ZipperOptions {
    /// Tags to remove from mapped reads during zipper merge.
    #[arg(long, num_args = 1..)]
    pub tags_to_remove: Option<Vec<String>>,
    /// Tags to reverse for negative strand reads during zipper merge.
    #[arg(long, num_args = 1..)]
    pub tags_to_reverse: Option<Vec<String>>,
    /// Tags to reverse complement for negative strand reads during zipper merge.
    #[arg(long, num_args = 1..)]
    pub tags_to_revcomp: Option<Vec<String>>,
    /// Skip adding primary alignment (PA) tags during zipper merge.
    #[arg(long, default_value = "false", num_args = 0..=1,
          default_missing_value = "true", action = clap::ArgAction::Set,
          value_parser = parse_bool)]
    pub skip_pa_tags: bool,
}

/// Options for UMI correction.
#[multi_options("correct", "Correct Options")]
#[derive(Args, Debug, Clone)]
pub struct CorrectOptions {
    /// Known UMI sequences for correction (one per argument).
    #[arg(long, num_args = 1..)]
    pub umis: Option<Vec<String>>,
    /// Files containing known UMI sequences (one per line).
    #[arg(long, num_args = 1..)]
    pub umi_files: Option<Vec<PathBuf>>,
    /// Maximum mismatches for UMI correction.
    #[arg(long, default_value_t = defaults::CORRECT_MAX_MISMATCHES)]
    pub max_mismatches: usize,
    /// Minimum difference between best and second-best match.
    #[arg(long, default_value_t = 2)]
    pub min_distance: usize,
    /// LRU cache size per worker for UMI segment matching (0 = no cache).
    #[arg(long, default_value_t = defaults::CORRECT_CACHE_SIZE)]
    pub cache_size: usize,
    /// Output path for records that could not be corrected (rejects BAM).
    #[arg(long)]
    pub rejects: Option<PathBuf>,
    /// Output path for per-UMI correction metrics TSV.
    #[arg(long)]
    pub metrics: Option<PathBuf>,
    /// Reverse complement UMIs before matching.
    #[arg(long, default_value = "false", num_args = 0..=1,
          default_missing_value = "true", action = clap::ArgAction::Set,
          value_parser = parse_bool)]
    pub revcomp: bool,
    /// Do not store original UMIs in the OX tag when correcting.
    #[arg(long, default_value = "false", num_args = 0..=1,
          default_missing_value = "true", action = clap::ArgAction::Set,
          value_parser = parse_bool)]
    pub dont_store_original_umis: bool,
    /// Minimum fraction of reads that must be correctable (0.0-1.0).
    #[arg(long)]
    pub min_corrected: Option<f64>,
}

impl Default for CorrectOptions {
    fn default() -> Self {
        Self {
            umis: None,
            umi_files: None,
            max_mismatches: defaults::CORRECT_MAX_MISMATCHES,
            min_distance: 2,
            cache_size: defaults::CORRECT_CACHE_SIZE,
            rejects: None,
            metrics: None,
            revcomp: false,
            dont_store_original_umis: false,
            min_corrected: None,
        }
    }
}

/// Options for duplex consensus calling.
#[multi_options("duplex", "Duplex Options")]
#[derive(Args, Debug, Clone, Default)]
pub struct DuplexOptions {
    /// Minimum reads for duplex consensus calling (1-3 comma-separated values:
    /// [duplex] or [duplex, AB/BA] or [duplex, AB, BA]).
    #[arg(long, value_delimiter = ',', default_value = "1")]
    pub min_reads: Option<Vec<usize>>,
    /// Maximum reads per strand for duplex (downsample if exceeded).
    #[arg(long)]
    pub max_reads_per_strand: Option<usize>,
}

/// Options for CODEC consensus calling.
#[multi_options("codec", "Codec Options")]
#[derive(Args, Debug, Clone)]
pub struct CodecOptions {
    /// Minimum read pairs per strand for CODEC consensus.
    #[arg(long, default_value_t = defaults::CODEC_MIN_READS)]
    pub min_reads: usize,
    /// Maximum read pairs per strand for CODEC consensus.
    #[arg(long)]
    pub max_reads: Option<usize>,
    /// Minimum duplex overlap length in bases for CODEC consensus.
    #[arg(long, default_value_t = defaults::CODEC_MIN_DUPLEX_LENGTH)]
    pub min_duplex_length: usize,
}

impl Default for CodecOptions {
    fn default() -> Self {
        Self {
            min_reads: defaults::CODEC_MIN_READS,
            max_reads: None,
            min_duplex_length: defaults::CODEC_MIN_DUPLEX_LENGTH,
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn stop_after_ordering_matches_pipeline_stage_order() {
        assert!(StopAfter::Extract < StopAfter::Correct);
        assert!(StopAfter::Correct < StopAfter::Fastq);
        assert!(StopAfter::Fastq < StopAfter::Align);
        assert!(StopAfter::Align < StopAfter::Zipper);
        assert!(StopAfter::Zipper < StopAfter::Sort);
        assert!(StopAfter::Sort < StopAfter::Group);
        assert!(StopAfter::Group < StopAfter::Consensus);
        assert!(StopAfter::Consensus < StopAfter::Filter);
    }
}
