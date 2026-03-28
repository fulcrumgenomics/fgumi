//! Enums and option structs for pipeline CLI arguments.

use std::path::PathBuf;

use clap::{Args, ValueEnum};
use fgumi_derive::multi_options;
use fgumi_lib::defaults;
use fgumi_pipeline::stages::aligner::AlignerPreset;
use fgumi_umi::Strategy;

/// Pipeline entry point: which stage to start from.
#[derive(Debug, Clone, Copy, ValueEnum)]
pub enum StartFrom {
    /// Start from grouped BAM (MI tags present). Runs: consensus -> filter -> write.
    Group,
    /// Start from coordinate-sorted BAM. Runs: sort -> group -> consensus -> filter -> write.
    Sort,
    /// Start from unmapped + mapped BAMs. Runs: zipper -> sort -> group -> consensus -> filter.
    Zipper,
    /// Start from unmapped BAM with UMIs. Runs: correct -> fastq -> aligner -> zipper -> sort -> ...
    Correct,
    /// Start from unmapped BAM. Runs: fastq -> aligner -> zipper -> sort -> ...
    Fastq,
    /// Start from FASTQ file. Runs: aligner -> zipper -> sort -> ...
    /// Requires `--unmapped` for the unmapped BAM used in zipper merge.
    Align,
    /// Start from FASTQ. Runs: extract -> fastq -> aligner -> zipper -> sort -> ... (full pipeline).
    Extract,
}

/// Consensus algorithm to use in the pipeline.
#[derive(Debug, Clone, Copy, ValueEnum)]
pub enum ConsensusMode {
    /// Simplex (vanilla UMI) consensus calling.
    Simplex,
    /// Duplex consensus calling (requires paired UMI grouping with /A and /B suffixes).
    #[cfg(feature = "duplex")]
    Duplex,
    /// CODEC consensus calling (R1/R2 from same molecule sequence opposite strands).
    #[cfg(feature = "codec")]
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
    /// Stop after FASTQ conversion (--start-from extract only).
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
    /// Two-character BAM tag containing the raw UMI.
    #[arg(long, default_value_t = String::from(defaults::UMI_TAG))]
    pub umi_tag: String,
    /// Minimum mapping quality for grouping. Templates where any mapped read has
    /// MAPQ below this threshold are discarded before UMI assignment.
    #[arg(long, default_value_t = 1)]
    pub min_mapq: u8,
    /// Allow fully unmapped templates (both reads unmapped) through grouping.
    #[arg(long, default_value_t = false)]
    pub allow_unmapped: bool,
}

impl Default for GroupOptions {
    fn default() -> Self {
        Self {
            strategy: Strategy::Adjacency,
            max_edits: defaults::GROUP_MAX_EDITS,
            umi_tag: String::from(defaults::UMI_TAG),
            min_mapq: 1,
            allow_unmapped: false,
        }
    }
}

/// Options for consensus read filtering.
#[multi_options("filter", "Filter Options")]
#[derive(Args, Debug, Clone)]
pub struct FilterOptions {
    /// Minimum number of reads to produce a consensus read.
    #[arg(long, default_value_t = 1)]
    pub min_reads: usize,
    /// Minimum input base quality for consensus calling.
    #[arg(long, default_value_t = 2)]
    pub min_base_quality: u8,
}

impl Default for FilterOptions {
    fn default() -> Self {
        Self { min_reads: 1, min_base_quality: 2 }
    }
}

/// Options for consensus calling.
#[multi_options("consensus", "Consensus Options")]
#[derive(Args, Debug, Clone)]
pub struct ConsensusOptions {
    /// Pre-UMI error rate in Phred scale.
    #[arg(long, default_value_t = defaults::CONSENSUS_ERROR_RATE_PRE_UMI)]
    pub error_rate_pre_umi: u8,
    /// Post-UMI error rate in Phred scale.
    #[arg(long, default_value_t = defaults::CONSENSUS_ERROR_RATE_POST_UMI)]
    pub error_rate_post_umi: u8,
    /// Read group ID for consensus reads.
    #[arg(long, default_value_t = String::from(defaults::READ_GROUP_ID))]
    pub read_group_id: String,
    /// Prefix for consensus read names (derived from header if not set).
    #[arg(long)]
    pub read_name_prefix: Option<String>,
    /// SAM tag containing the molecular identifier.
    #[arg(long, default_value_t = String::from(defaults::MI_TAG))]
    pub mi_tag: String,
}

impl Default for ConsensusOptions {
    fn default() -> Self {
        Self {
            error_rate_pre_umi: defaults::CONSENSUS_ERROR_RATE_PRE_UMI,
            error_rate_post_umi: defaults::CONSENSUS_ERROR_RATE_POST_UMI,
            read_group_id: String::from(defaults::READ_GROUP_ID),
            read_name_prefix: None,
            mi_tag: String::from(defaults::MI_TAG),
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
    #[arg(long, short = 'R', num_args = 0..)]
    pub read_structures: Option<Vec<String>>,
    /// Extract UMIs from field 8 of the Illumina read name (colon-delimited).
    #[arg(long, default_value_t = false)]
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
    #[arg(long, num_args = 0..)]
    pub comment: Option<Vec<String>>,

    // ---- Extract tag/behavior options ----
    /// SAM tag for storing cell barcode sequences.
    #[arg(long, default_value_t = String::from(defaults::CELL_TAG))]
    pub cell_tag: String,
    /// SAM tag for storing UMI quality scores (e.g. QX).
    #[arg(long)]
    pub umi_qual_tag: Option<String>,
    /// SAM tag for storing cell barcode quality scores.
    #[arg(long)]
    pub cell_qual_tag: Option<String>,
    /// Single tag to store all concatenated UMIs (in addition to per-segment tags).
    #[arg(long)]
    pub single_tag: Option<String>,
    /// Annotate read names with UMIs (appends "+UMIs" to read names).
    #[arg(long, default_value_t = false)]
    pub annotate_read_names: bool,
    /// Store sample barcode quality scores in the QT tag.
    #[arg(long, default_value_t = false)]
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
            cell_tag: String::from(defaults::CELL_TAG),
            umi_qual_tag: None,
            cell_qual_tag: None,
            single_tag: None,
            annotate_read_names: false,
            store_sample_barcode_qualities: false,
        }
    }
}

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
}

impl Default for AlignerOptions {
    fn default() -> Self {
        Self { preset: None, command: None, threads: 4 }
    }
}

/// Options for the zipper merge stage.
#[multi_options("zipper", "Zipper Options")]
#[derive(Args, Debug, Clone, Default)]
pub struct ZipperOptions {
    /// Tags to remove from mapped reads during zipper merge.
    #[arg(long, num_args = 0..)]
    pub tags_to_remove: Option<Vec<String>>,
    /// Tags to reverse for negative strand reads during zipper merge.
    #[arg(long, num_args = 0..)]
    pub tags_to_reverse: Option<Vec<String>>,
    /// Tags to reverse complement for negative strand reads during zipper merge.
    #[arg(long, num_args = 0..)]
    pub tags_to_revcomp: Option<Vec<String>>,
    /// Skip adding primary alignment (PA) tags during zipper merge.
    #[arg(long, default_value_t = false)]
    pub skip_pa_tags: bool,
}

/// Options for UMI correction.
#[multi_options("correct", "Correct Options")]
#[derive(Args, Debug, Clone)]
pub struct CorrectOptions {
    /// Known UMI sequences for correction (one per argument).
    #[arg(long, num_args = 0..)]
    pub umis: Option<Vec<String>>,
    /// Files containing known UMI sequences (one per line).
    #[arg(long, num_args = 0..)]
    pub umi_files: Option<Vec<PathBuf>>,
    /// Maximum mismatches for UMI correction.
    #[arg(long, default_value_t = defaults::CORRECT_MAX_MISMATCHES)]
    pub max_mismatches: usize,
    /// LRU cache size per worker for UMI segment matching (0 = no cache).
    #[arg(long, default_value_t = defaults::CORRECT_CACHE_SIZE)]
    pub cache_size: usize,
}

impl Default for CorrectOptions {
    fn default() -> Self {
        Self {
            umis: None,
            umi_files: None,
            max_mismatches: defaults::CORRECT_MAX_MISMATCHES,
            cache_size: defaults::CORRECT_CACHE_SIZE,
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
