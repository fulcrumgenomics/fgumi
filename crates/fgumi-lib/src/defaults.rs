//! Shared default values for CLI options.
//!
//! These constants are the single source of truth for default values used by both
//! standalone commands and the runall pipeline. Both places import from here.

// -- SAM Tags --

/// Default UMI tag.
pub const UMI_TAG: &str = "RX";
/// Default molecule identifier tag.
pub const MI_TAG: &str = "MI";
/// Default cell barcode tag.
pub const CELL_TAG: &str = "CB";

// -- Extract --

/// Default read group ID (used by both extract and consensus stages).
pub const READ_GROUP_ID: &str = "A";
/// Default sequencing platform.
pub const PLATFORM: &str = "illumina";

// -- Correct --

/// Default maximum mismatches for UMI correction.
pub const CORRECT_MAX_MISMATCHES: usize = 2;
/// Default LRU cache size for UMI matching.
pub const CORRECT_CACHE_SIZE: usize = 100_000;

// -- Fastq / Aligner --

/// Default SAM flags to exclude when converting to FASTQ.
pub const FASTQ_EXCLUDE_FLAGS: u16 = 0x900;
/// Default SAM flags to require when converting to FASTQ.
pub const FASTQ_REQUIRE_FLAGS: u16 = 0;
/// Default BWA chunk size in bases (BWA `-K` parameter).
pub const BWA_CHUNK_SIZE: u64 = 150_000_000;

// -- Sort --

/// Default memory limit for external sort in bytes (768 MiB).
pub const SORT_MEMORY_LIMIT: usize = 768 * 1024 * 1024;
/// Display string for the sort memory limit (for `--help` text in standalone).
pub const SORT_MEMORY_LIMIT_DISPLAY: &str = "768M";
/// Default: scale memory limit per thread.
pub const SORT_MEMORY_PER_THREAD: bool = true;
/// Default compression level for sort temporary files.
pub const SORT_TEMP_COMPRESSION: u32 = 1;

// -- Group --

/// Default maximum edit distance between UMIs for grouping.
pub const GROUP_MAX_EDITS: u32 = 1;
/// Default minimum UMIs per position to use index-based grouping.
pub const GROUP_INDEX_THRESHOLD: usize = 100;

// -- Consensus --

/// Default Phred-scaled error rate prior to UMI integration.
pub const CONSENSUS_ERROR_RATE_PRE_UMI: u8 = 45;
/// Default Phred-scaled error rate after UMI integration.
pub const CONSENSUS_ERROR_RATE_POST_UMI: u8 = 40;
/// Default minimum input base quality for consensus calling.
pub const CONSENSUS_MIN_INPUT_BASE_QUALITY: u8 = 10;
/// Default: produce per-base tags (cd, ce) in consensus output.
pub const CONSENSUS_OUTPUT_PER_BASE_TAGS: bool = true;
/// Default: do not trim reads before consensus calling.
pub const CONSENSUS_TRIM: bool = false;
/// Default: consensus-call overlapping bases in read pairs.
pub const CONSENSUS_CALL_OVERLAPPING_BASES: bool = true;
/// Default minimum consensus base quality (bases below are masked to N).
pub const CONSENSUS_MIN_BASE_QUALITY: u8 = 2;

// -- Codec --

/// Default minimum reads per UMI group for codec consensus.
pub const CODEC_MIN_READS: usize = 1;
/// Default minimum duplex length for codec consensus.
pub const CODEC_MIN_DUPLEX_LENGTH: usize = 1;
/// Default outer bases length for codec consensus.
pub const CODEC_OUTER_BASES_LENGTH: usize = 5;
/// Default maximum duplex disagreement rate for codec consensus.
pub const CODEC_MAX_DUPLEX_DISAGREEMENT_RATE: f64 = 1.0;

// -- Filter --

/// Default maximum per-read error rate (string form for `Vec<f64>` defaults).
pub const FILTER_MAX_READ_ERROR_RATE: &str = "0.025";
/// Default maximum per-base error rate (string form for `Vec<f64>` defaults).
pub const FILTER_MAX_BASE_ERROR_RATE: &str = "0.1";
/// Default maximum fraction of no-call (N) bases.
pub const FILTER_MAX_NO_CALL_FRACTION: f64 = 0.2;

// -- Compression --

/// Default BAM compression level.
pub const COMPRESSION_LEVEL: u32 = 1;
