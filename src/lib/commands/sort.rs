//! Sort BAM files by various orderings.
//!
//! Uses high-performance raw-bytes sorting with radix sort for in-memory
//! chunks and O(1) merge comparisons via pre-computed sort keys.
//!
//! # Sort Orders
//!
//! - **Template-coordinate**: Groups paired-end reads by template position (for `fgumi group`)
//! - **Queryname**: Groups reads by read name (for `fgumi zipper`)
//! - **Coordinate**: Standard genomic coordinate order (for IGV, `fgumi review`)
//!
//! # Performance
//!
//! - 1.9x faster than samtools on template-coordinate sort
//! - Handles BAM files larger than available RAM via spill-to-disk
//! - Uses parallel sorting for in-memory chunks
//! - Configurable temp file compression (--temp-compression)
//!
//! # Verification
//!
//! Use `--verify` to check if a BAM file is correctly sorted without writing output.

use crate::logging::OperationTimer;
use crate::sam::SamTag;
use crate::validation::validate_input_exists;
use anyhow::{Result, bail};
use bytesize::ByteSize;
use clap::Parser;
use fgumi_bam_io::is_stdout_path;
use fgumi_sort::{KeyTypesSpec, QuerynameComparator, SortOrder, verify_sort_order};

use log::{debug, info, warn};
use std::path::{Path, PathBuf};

use crate::commands::command::Command;
use crate::commands::common::{
    CompressionOptions, MaxTempFiles, MemoryLimit, MemoryReserve, QueueMemoryOptions,
    SchedulerOptions, ThreadingOptions, parse_max_temp_files, parse_memory, parse_memory_reserve,
    resolve_memory_budget,
};
use crate::pipeline::chains::{ChainSpec, SinkSpec, SourceSpec, Stage, StageOptionsBag, build_for};

/// Sort order for BAM files.
///
/// Queryname sort supports sub-sort specification via `::` syntax:
/// - `queryname` — lexicographic ordering (default, fast)
/// - `queryname::lexicographic` — explicit lexicographic ordering
/// - `queryname::natural` — natural numeric ordering (samtools-compatible)
///
/// The `@HD SS` sub-sort for lexicographic queryname is emitted with the
/// SAM-spec spelling `lexicographical` (accepted back as `queryname::lexicographical`).
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum SortOrderArg {
    /// Coordinate sort (tid → pos → strand)
    Coordinate,
    /// Queryname sort with lexicographic ordering (default)
    Queryname,
    /// Queryname sort with natural numeric ordering
    QuerynameNatural,
    /// Template-coordinate sort (for UMI grouping)
    TemplateCoordinate,
}

impl SortOrderArg {
    /// Parse a sort order string, supporting `::` sub-sort syntax for queryname.
    ///
    /// Valid values:
    /// - `coordinate`
    /// - `queryname` (default: lexicographic)
    /// - `queryname::lex`
    /// - `queryname::lexicographic`
    /// - `queryname::lexicographical` (alias; fgumi emits `queryname:lexicographical` in `@HD` SS)
    /// - `queryname::natural`
    /// - `template-coordinate`
    ///
    /// # Errors
    ///
    /// Returns an error if the string is not a valid sort order or has an
    /// unrecognized sub-sort specifier.
    pub fn parse(s: &str) -> Result<Self, String> {
        match s {
            "coordinate" => Ok(Self::Coordinate),
            "queryname"
            | "queryname::lex"
            | "queryname::lexicographic"
            | "queryname::lexicographical" => Ok(Self::Queryname),
            "queryname::natural" => Ok(Self::QuerynameNatural),
            "template-coordinate" => Ok(Self::TemplateCoordinate),
            other => {
                if other.starts_with("queryname::") {
                    let sub =
                        other.strip_prefix("queryname::").expect("guarded by starts_with check");
                    Err(format!(
                        "unknown queryname sub-sort '{sub}', expected 'lex', 'lexicographic', \
                         'lexicographical', or 'natural'"
                    ))
                } else {
                    Err(format!(
                        "unknown sort order '{other}', expected 'coordinate', 'queryname', \
                         'queryname::lex', 'queryname::lexicographic', \
                         'queryname::lexicographical', 'queryname::natural', \
                         or 'template-coordinate'"
                    ))
                }
            }
        }
    }
}

impl From<SortOrderArg> for SortOrder {
    fn from(arg: SortOrderArg) -> Self {
        match arg {
            SortOrderArg::Coordinate => SortOrder::Coordinate,
            SortOrderArg::Queryname => SortOrder::Queryname(QuerynameComparator::Lexicographic),
            SortOrderArg::QuerynameNatural => SortOrder::Queryname(QuerynameComparator::Natural),
            SortOrderArg::TemplateCoordinate => SortOrder::TemplateCoordinate,
        }
    }
}

/// Sort a BAM file.
///
/// Sorts BAM files using high-performance external merge-sort, supporting
/// multiple sort orders required by the fgumi pipeline.
#[derive(Debug, Parser)]
#[command(
    name = "sort",
    about = "\x1b[38;5;72m[ALIGNMENT]\x1b[0m      \x1b[36mSort BAM file by coordinate, queryname, or template-coordinate\x1b[0m",
    long_about = r#"
Sort a BAM file using high-performance external merge-sort.

This tool provides efficient BAM sorting with support for multiple sort orders:

SORT ORDERS:

  coordinate              Standard genomic coordinate sort (tid → pos → strand).
                          Use for IGV visualization, variant calling, `fgumi review`.

  queryname                Lexicographic read name sort (fast, default sub-sort).
  queryname::lex           Short alias for lexicographic ordering (same as above).
  queryname::lexicographic Explicit lexicographic ordering (same as above).
  queryname::lexicographical Alias for the above; fgumi emits `queryname:lexicographical` in @HD SS.
  queryname::natural       Natural numeric ordering (samtools-compatible).
                          Use for `fgumi zipper`, template-level operations.

  template-coordinate      Template-level position sort for UMI grouping.
                          Use for `fgumi group`, `fgumi dedup`, and `fgumi downsample` input.

PERFORMANCE:

  - 1.9x faster than samtools on template-coordinate sort
  - Handles BAM files larger than available RAM via spill-to-disk
  - Uses parallel sorting (--threads) for in-memory chunks
  - Configurable temp file compression (--temp-compression)
  - Default 768M per-thread memory limit (samtools-compatible); pass
    `--max-memory auto` to detect system memory (opt-in)
  - Spilled runs are consolidated once they reach `--max-temp-files`
    (default "auto": sized to the process's open-file limit, `ulimit -n`);
    consolidation rewrites already-sorted data, so raising that limit
    avoids it on very large inputs

EXAMPLES:

  # Sort for fgumi group input
  fgumi sort -i aligned.bam -o sorted.bam --order template-coordinate

  # Sort by coordinate for IGV
  fgumi sort -i input.bam -o sorted.bam --order coordinate

  # Sort by queryname for zipper
  fgumi sort -i input.bam -o sorted.bam --order queryname

  # Multi-threaded sort (default 768M per thread)
  fgumi sort -i input.bam -o sorted.bam --order template-coordinate --threads 8

  # Override the per-thread memory limit
  fgumi sort -i input.bam -o sorted.bam -m 2GiB --threads 8

  # Opt in to auto-detected system memory (subtracts --memory-reserve)
  fgumi sort -i input.bam -o sorted.bam -m auto --threads 8

  # Reserve extra memory for bwa mem running in a pipeline
  fgumi sort -i input.bam -o sorted.bam --memory-reserve 12GiB --threads 4

  # Cede cores to the aligner during ingest, but keep the merge wide
  bwa mem -t 32 ref.fa r1.fq r2.fq | fgumi sort -i - -o sorted.bam -@ 8 --sort-threads 4

  # Allow more spilled runs before consolidating (fewer consolidation passes)
  fgumi sort -i input.bam -o sorted.bam --order coordinate --max-temp-files 512

  # Verify a BAM file is correctly sorted
  fgumi sort -i sorted.bam --verify --order template-coordinate

  # Spread spill chunks across multiple temp dirs (round-robin, free-space aware)
  fgumi sort -i in.bam -o out.bam -T /mnt/ssd1 -T /mnt/ssd2

  # Same via FGUMI_TMP_DIRS env var (PATH-style list)
  FGUMI_TMP_DIRS=/mnt/ssd1:/mnt/ssd2 fgumi sort -i in.bam -o out.bam
"#
)]
#[allow(clippy::struct_excessive_bools)]
pub struct Sort {
    /// Input BAM file.
    #[arg(short = 'i', long = "input")]
    pub input: PathBuf,

    /// Output BAM file (required unless --verify is used).
    #[arg(short = 'o', long = "output")]
    pub output: Option<PathBuf>,

    /// Verify the input file is correctly sorted (no output written).
    ///
    /// Reads records sequentially and checks that each record's sort key
    /// is >= the previous record's key. Exits 0 if sorted correctly,
    /// non-zero if any records are out of order.
    #[arg(long = "verify", value_name = "true|false", default_value = "false", num_args = 0..=1, default_missing_value = "true", action = clap::ArgAction::Set, value_parser = clap::builder::BoolishValueParser::new(), hide_possible_values = true)]
    pub verify: bool,

    /// Sort order.
    ///
    /// Queryname sort supports sub-sort specifiers:
    ///   `queryname`                  Lexicographic byte ordering (default, fast)
    ///   `queryname::lexicographic`   Explicit lexicographic ordering (alias: `queryname::lex`)
    ///   `queryname::lexicographical` Alias; written as `queryname:lexicographical` in `@HD` SS
    ///   `queryname::natural`         Natural numeric ordering (samtools-compatible)
    #[arg(long = "order", default_value = "template-coordinate", value_parser = SortOrderArg::parse)]
    pub order: SortOrderArg,

    /// Which optional lanes to keep in the template-coordinate sort key.
    ///
    /// Smaller keys use less memory and spill less. Only meaningful for
    /// `--order template-coordinate`; ignored for other orders.
    ///
    ///   (omitted)            Auto-detect from the first record + verify (default).
    ///   full                 Keep all lanes (CB + library/MI). Largest key.
    ///   none                 Drop all optional lanes (smallest, bulk pre-group).
    ///   cb,library,mi        Comma/space list; keep the named lanes.
    ///
    /// A record carrying a value in a dropped lane aborts the sort with a message
    /// naming the field and the token to re-include it.
    #[arg(long = "key-types", value_parser = parse_key_types)]
    pub key_types: Option<KeyTypesSpec>,

    /// Maximum memory for in-memory sorting.
    ///
    /// Default is "768M" per thread (matching samtools behavior). Pass "auto"
    /// to detect system memory and subtract --memory-reserve, leaving room
    /// for the OS and co-running processes (e.g. an aligner). Explicit values
    /// like "512M", "1G", "4GiB" are per-thread when --memory-per-thread is
    /// enabled (default).
    ///
    /// When the limit is reached, sorted chunks spill to temporary files.
    ///
    /// This bounds two budgets: the in-memory record buffer (the dominant
    /// consumer, reported on the "Max memory" log line) and the bytes held in the
    /// queues between pipeline stages (applied, but not separately logged). It is
    /// not a hard cap on total process RSS: ingest decompression/parse buffers add
    /// a few percent, and the final k-way merge's per-file read-ahead is separate
    /// again, so peak RSS runs somewhat above this value. When sorting in a pipe
    /// alongside an aligner, size the budget to leave headroom for the aligner's
    /// resident index on top of this.
    ///
    /// The two budgets scale by different thread counts: the record buffer by
    /// max(--threads, --sort-threads) (the sort phase fills it), the queue budget
    /// by --threads (the width of the ingest/output plumbing), so their resolved
    /// totals can differ for the same value when --sort-threads > --threads.
    /// Raising this value does not raise the queues' per-stage backpressure marks,
    /// so a large or `auto` value bounds the queues at those marks, not the full
    /// total; only a value below them tightens the queues.
    #[arg(short = 'm', long = "max-memory", default_value = "768M", value_parser = parse_memory)]
    pub max_memory: MemoryLimit,

    /// Memory to reserve for other processes when --max-memory=auto.
    ///
    /// "auto" (default) reserves min(10 GiB, 50% of system memory). Explicit
    /// values like "10G", "8GiB" set a fixed reservation. Set higher when
    /// running alongside a memory-intensive aligner (e.g. `bwa mem` with a
    /// human genome index uses ~8 GiB).
    ///
    /// Ignored when --max-memory is set to an explicit value.
    #[arg(long = "memory-reserve", default_value = "auto", value_parser = parse_memory_reserve)]
    pub memory_reserve: MemoryReserve,

    /// Scale memory limit by thread count (samtools behavior).
    ///
    /// When enabled (default), --max-memory specifies memory per thread. Total
    /// memory = `max_memory` × the larger of --threads and --sort-threads, since
    /// the sort phase is what fills the in-memory buffer. Disable for fixed total
    /// memory.
    ///
    /// This formula is for the in-memory sort buffer. --max-memory also bounds
    /// the inter-stage queue budget (see its docs), which scales by --threads
    /// alone, so the two totals differ when --sort-threads > --threads.
    #[arg(long = "memory-per-thread", value_name = "true|false", default_value = "true", num_args = 0..=1, default_missing_value = "true", action = clap::ArgAction::Set, value_parser = clap::builder::BoolishValueParser::new(), hide_possible_values = true)]
    pub memory_per_thread: bool,

    /// Temporary directory for intermediate files. Repeatable.
    ///
    /// Pass `-T <path>` one or more times to spread spill chunks across multiple
    /// directories in free-space-aware round-robin order. Useful when one
    /// filesystem is too small or slower than the aggregate of several.
    ///
    /// If no flags are given and the `FGUMI_TMP_DIRS` environment variable is
    /// set, its value is parsed as a `PATH`-style list (colon-separated on
    /// Unix, semicolon-separated on Windows) and used instead.
    ///
    /// If neither is provided, the system default temp directory is used.
    /// For best performance, use fast SSDs.
    #[arg(short = 'T', long = "tmp-dir", action = clap::ArgAction::Append)]
    pub tmp_dirs: Vec<PathBuf>,

    /// Number of threads for parallel operations.
    ///
    /// Used for parallel sorting of in-memory chunks and parallel temp-chunk
    /// (BGZF or zstd) compression.
    ///
    /// This is also the floor for the `--max-memory --memory-per-thread`
    /// multiplier: the sort runs on one global memory budget, sized by the larger
    /// of this and `--sort-threads` (the phase that fills the buffer).
    /// `--merge-threads` changes scheduling only and never resizes it.
    #[arg(short = '@', short_alias = 't', long = "threads", default_value = "1")]
    pub threads: usize,

    /// Number of threads for the sort phase (accumulate, sort, spill).
    ///
    /// Defaults to `--threads`. Lower this to cede cores to an upstream
    /// producer while keeping the merge wide -- with `-@ 8 --sort-threads 4`,
    /// ingest contends with the producer over only 4 threads, while the merge
    /// still uses 8 because it cannot start until the input is exhausted, by
    /// which point the producer has finished writing.
    ///
    /// The output is byte-identical, but this is not purely a scheduling knob:
    /// with --memory-per-thread enabled (default) the budget scales by the larger
    /// of --threads and --sort-threads, so raising this above --threads raises
    /// total memory by the same factor.
    #[arg(long = "sort-threads")]
    pub sort_threads: Option<usize>,

    /// Number of threads for the merge phase (k-way merge and output write).
    ///
    /// Defaults to `--threads`. This only changes scheduling; the output is
    /// byte-identical.
    #[arg(long = "merge-threads")]
    pub merge_threads: Option<usize>,

    /// Compression options for output BAM.
    #[command(flatten)]
    pub compression: CompressionOptions,

    /// Compression level for temporary chunk files (0-9).
    ///
    /// Applies to the codec selected by `--temp-codec`:
    ///   * For `bgzf`, level 0 produces uncompressed (stored) BGZF blocks
    ///     (fastest, uses most disk space); 1..=9 are libdeflate levels.
    ///   * For `zstd`, only 1..=9 are valid; level 0 is rejected because zstd
    ///     has no equivalent "stored" mode and silently remapping it to 1
    ///     would surprise users counting on uncompressed spill.
    ///
    /// Level 1 (default) provides fast compression with reasonable space savings.
    /// Higher levels (up to 9) provide better compression but are slower.
    #[arg(long = "temp-compression", default_value = "1", value_parser = clap::value_parser!(u32).range(0..=9))]
    pub temp_compression: u32,

    /// Codec used for temporary spill chunks: `zstd` (default) or `bgzf`.
    ///
    /// zstd is significantly faster than bgzf at comparable compression
    /// ratios for BAM-record data; we default to zstd because spill files
    /// are internal to the sort and never read by other tools. Pass `bgzf`
    /// to fall back to the legacy on-disk format.
    #[arg(long = "temp-codec", default_value = "zstd")]
    pub temp_codec: fgumi_sort::SpillCodec,

    /// Maximum number of temporary spill files kept before the oldest are
    /// consolidated into a single run.
    ///
    /// Large inputs spill many sorted runs to disk. When the number of runs
    /// reaches this limit, the oldest are merged together in a single pass so
    /// the final k-way merge opens fewer files at once. That merge is the only
    /// reason the limit exists: it opens every remaining run at once, so the
    /// limit bounds how many file descriptors the sort needs.
    ///
    /// Consolidation rewrites data that is already sorted, so it is pure
    /// overhead whenever the descriptor budget could have carried the runs.
    /// Raising this avoids it on very large inputs (at the cost of more open
    /// file descriptors during the final merge); lowering it keeps fewer files
    /// open. Must be at least 2; to effectively disable consolidation, pass a
    /// value larger than the number of runs you expect to spill.
    ///
    /// "auto" (default) sizes the limit to the process's soft open-file limit
    /// (`ulimit -n`), less a reserve for the input, output and index handles,
    /// and capped at a tested maximum. Explicit values like "64", "256" pin it;
    /// a pinned value larger than the open-file budget is reported at startup.
    /// Must be at least 2.
    #[arg(long = "max-temp-files", default_value = "auto", value_parser = parse_max_temp_files)]
    pub max_temp_files: MaxTempFiles,

    /// Write BAM index (.bai) alongside output.
    ///
    /// Only valid for coordinate sort, and not with output to stdout (`-` /
    /// `/dev/stdout`): an index needs a seekable file and a sidecar path, so
    /// that combination is rejected rather than silently skipped.
    ///
    /// The index file will be written to
    /// `<output>.bai`. Output BGZF compression stays multi-threaded (scales
    /// with `--threads`); the BAI virtual offsets are recovered from each BGZF
    /// block as it finalizes, so indexing does not serialize compression.
    #[arg(long = "write-index", value_name = "true|false", default_value = "false", num_args = 0..=1, default_missing_value = "true", action = clap::ArgAction::Set, value_parser = clap::builder::BoolishValueParser::new(), hide_possible_values = true)]
    pub write_index: bool,

    /// Concurrent positional reads to fetch the input BAM with: `auto`, or a
    /// fixed count (`1` = the plain sequential reader).
    ///
    /// One blocking read keeps only about one read-ahead window outstanding at
    /// the device, which on network-attached storage is far below what it can
    /// serve -- and no buffer size fixes it, because a bigger single request is
    /// not more concurrency. How much that matters depends entirely on the
    /// storage: on EBS gp3 a 38 GB coordinate sort ran ~14% faster at four
    /// streams, while on a direct-attached instance-store SSD one stream is
    /// already faster than four are on gp3 and adding more costs 1.8%.
    ///
    /// `auto` therefore measures instead of guessing. It reads the first several
    /// fills at a single stream -- exactly the pre-existing behaviour -- then
    /// commits once to a count of `ceil(target-throughput / measured)`, capped at
    /// 8. It lands on four for gp3 and stays at one for `NVMe` without being told
    /// which is which, and does not revisit the choice afterwards.
    ///
    /// Each active stream is served by a scoped OS thread spawned per fill window,
    /// bounded by the chosen count (at most 8) -- this is independent of
    /// `--threads`. Applies to the input BAM only, not the merge's spill files.
    /// Ignored for stdin and other non-seekable inputs, where positional reads do
    /// not exist.
    #[arg(long = "read-streams", default_value_t = fgumi_sort::ReadStreams::Auto)]
    pub read_streams: fgumi_sort::ReadStreams,

    /// Emit the sort's performance diagnostics.
    ///
    /// Prints spill geometry, per-phase timing, the merge's floor (which of
    /// consumer serial CPU, worker capacity, or coordination limits it, and how
    /// much is recoverable), worker utilization, and park attribution -- roughly
    /// a hundred lines. Off by default: it is instrumentation for performance
    /// work, read from a log with a grep, and it buried the handful of lines a
    /// normal run should show.
    #[arg(long = "sort-stats", default_value_t = false, hide = true)]
    pub sort_stats: bool,
}

/// Sort-stage tuning, projected out of the [`Sort`] CLI struct for the chain
/// builder's `add_sort` (which reads these values from the
/// [`crate::pipeline::chains::StageOptionsBag`] rather than from `Sort`
/// directly).
///
/// `block_batch` and `file_granularity` are chain-engine knobs that `Sort`
/// does not yet expose as CLI flags; until the sort command is rewired onto the
/// chain they take the engine defaults (`block_batch = 4`, the original
/// `MAX_BATCH_PER_CALL`; `file_granularity = false`, block-parallel).
#[derive(Debug, Clone)]
#[allow(clippy::struct_excessive_bools)]
pub struct SortOptions {
    /// Requested output sort order.
    pub order: SortOrderArg,
    /// Expert override for the provisioned template-coordinate key lanes.
    pub key_types: Option<KeyTypesSpec>,
    /// In-memory sort budget before spilling.
    pub max_memory: MemoryLimit,
    /// Memory reserved for other processes under `--max-memory=auto`.
    pub memory_reserve: MemoryReserve,
    /// Whether `max_memory` is per-thread (samtools behavior).
    pub memory_per_thread: bool,
    /// Temp directories for spill chunks (free-space-aware round-robin).
    pub tmp_dirs: Vec<PathBuf>,
    /// Worker threads for the accumulate/sort/spill phase (Phase 1).
    pub sort_threads: Option<usize>,
    /// Worker threads for the k-way merge phase (Phase 2).
    pub merge_threads: Option<usize>,
    /// Compression level for temporary spill files (0-9).
    pub temp_compression: u32,
    /// Codec for temporary spill files.
    pub temp_codec: fgumi_sort::SpillCodec,
    /// Spill-file consolidation limit (`--max-temp-files`). Resolved against the
    /// host `RLIMIT_NOFILE` when `Auto`; without this the chain sorter fell back
    /// to the engine's portable default and ignored the CLI value.
    pub max_temp_files: MaxTempFiles,
    /// Records batched per parallel sort call (chain engine; not a CLI flag).
    pub block_batch: usize,
    /// Spill at file rather than block granularity (chain engine; not a CLI flag).
    pub file_granularity: bool,
}

impl Sort {
    /// Projects the parsed CLI flags into [`SortOptions`] for the chain builder.
    #[must_use]
    pub fn to_sort_options(&self) -> SortOptions {
        SortOptions {
            order: self.order,
            key_types: self.key_types,
            max_memory: self.max_memory,
            memory_reserve: self.memory_reserve,
            memory_per_thread: self.memory_per_thread,
            tmp_dirs: self.tmp_dirs.clone(),
            sort_threads: self.sort_threads,
            merge_threads: self.merge_threads,
            temp_compression: self.temp_compression,
            temp_codec: self.temp_codec,
            max_temp_files: self.max_temp_files,
            // Chain-engine defaults (see `SortOptions` docs) — not yet CLI flags.
            block_batch: 4,
            file_granularity: false,
        }
    }

    /// Projects the sort command's single memory budget onto the pipeline's
    /// inter-stage queue budget ([`QueueMemoryOptions`]).
    ///
    /// The sort command exposes exactly one memory knob (`--max-memory`, plus
    /// `--memory-reserve` and `--memory-per-thread`), so that one knob bounds
    /// *both* the sorter's in-memory record buffer (via [`Self::to_sort_options`])
    /// and the bytes held in the queues between pipeline stages. Without this the
    /// chain hardcoded `QueueMemoryOptions::default()` (768 MiB/thread), so a small
    /// `--max-memory` shrank the record buffer while the queues stayed large.
    ///
    /// The two budgets do not resolve to the same *number*: the sorter scales by
    /// `max(--threads, --sort-threads)` (`memory_budget_threads`, since the
    /// sort phase is what fills the buffer), whereas the queue budget scales by
    /// plain `--threads` (the width of the ingest/output plumbing the queues
    /// track). So a run with `--sort-threads` > `--threads` resolves a larger
    /// per-command sorter budget than queue budget for the same `--max-memory`;
    /// that asymmetry is intentional, not a drift.
    ///
    /// Two consequences of sharing the flag are deliberate and benign:
    /// - The default queue budget now follows `--max-memory`'s default (`"768M"`
    ///   = 768 MB decimal via the size parser), marginally below the former
    ///   hardcoded `QueueMemoryOptions::default()` of 768 MiB. This aligns the
    ///   queue default with the sorter's own default rather than leaving the two
    ///   at different units.
    /// - A large or `auto` `--max-memory` does not balloon queue memory: raising
    ///   the queue *total* does not raise the per-stage backpressure high-water
    ///   marks (512 MiB / 256 MiB — see [`QueueMemoryOptions`] and issue #765),
    ///   which is what actually bounds in-flight queue bytes. A budget *below*
    ///   those marks tightens them (the point of a small `--max-memory`); a budget
    ///   above them leaves them alone.
    #[must_use]
    pub fn queue_memory_options(&self) -> QueueMemoryOptions {
        QueueMemoryOptions {
            max_memory: self.max_memory,
            memory_reserve: self.memory_reserve,
            memory_per_thread: self.memory_per_thread,
        }
    }

    /// Builds the declarative [`ChainSpec`] that `fgumi sort` runs.
    ///
    /// Extracted from [`Self::execute_sort`] so the wiring is unit-testable: the
    /// spec's `queue_memory` (see [`Self::queue_memory_options`]), `stages`,
    /// `threading`, and `write_index`-selected sink are all otherwise only
    /// exercised through a full run, where a dropped knob is invisible (it does
    /// not change output bytes).
    ///
    /// `resolved_tmp_dirs` is the already-resolved temp-dir list (CLI flags with
    /// the `FGUMI_TMP_DIRS` fallback applied); it is baked into the `SortOptions`
    /// verbatim because `add_sort` reads `SortOptions.tmp_dirs` and never consults
    /// the environment.
    fn build_sort_chain_spec(
        &self,
        output: &Path,
        resolved_tmp_dirs: Vec<PathBuf>,
        command_line: &str,
    ) -> ChainSpec {
        let mut sort_options = self.to_sort_options();
        sort_options.tmp_dirs = resolved_tmp_dirs;

        let stage_opts = StageOptionsBag { sort: Some(sort_options), ..Default::default() };
        let sink = if self.write_index {
            SinkSpec::BamWithIndex(output.to_path_buf())
        } else {
            SinkSpec::Bam(output.to_path_buf())
        };
        ChainSpec {
            stages: vec![Stage::Sort],
            source: SourceSpec::Bam(self.input.clone()),
            sink,
            stage_opts,
            threading: ThreadingOptions { threads: Some(self.threads) },
            compression: self.compression.clone(),
            // Match sibling chain commands (10s), not the derived Default (0 = monitor off).
            scheduler: SchedulerOptions { deadlock_timeout: 10, ..Default::default() },
            // The single `--max-memory` knob bounds the inter-stage queue budget too,
            // not just the sorter buffer (see `queue_memory_options`). The queue
            // budget scales by `--threads`; the sorter's by max(threads, sort_threads).
            queue_memory: self.queue_memory_options(),
            async_reader: false,
            // Thread the sort command's --read-streams into the chain's BAM
            // source: Auto probes the device and picks a concurrent-read count,
            // Fixed(n) pins it, Fixed(1) is the plain sequential reader.
            read_streams: self.read_streams,
            // The owned engine always verified CRC (incl. stdin); keep parity -- a future
            // PR can add --no-check-crc if opt-out is wanted. `effective_check_crc()` would
            // skip verification for stdin input (the file-vs-stdin default other commands
            // use), which is a silent regression from the owned sorter's behavior.
            verify_crc: true,
            command_line: command_line.to_string(),
        }
    }
}

/// Environment variable name for the fallback temp-dir list, parsed as a
/// `PATH`-style list when no `-T/--tmp-dir` flags are passed.
pub(crate) const TMP_DIRS_ENV: &str = "FGUMI_TMP_DIRS";

/// Resolve the final list of temp directories for a sort run.
///
/// Precedence: CLI flags (if non-empty) > `FGUMI_TMP_DIRS` env var > empty.
/// Empty strings and whitespace-only entries are filtered out of the env-var
/// value so that `FGUMI_TMP_DIRS=:` or trailing separators don't produce bogus
/// paths.
pub(crate) fn resolve_tmp_dirs(cli: &[PathBuf], env_value: Option<&str>) -> Vec<PathBuf> {
    if !cli.is_empty() {
        return cli.to_vec();
    }

    let Some(value) = env_value else { return Vec::new() };
    if value.is_empty() {
        return Vec::new();
    }

    std::env::split_paths(value)
        .filter(|p| !p.as_os_str().is_empty())
        .filter(|p| !p.to_string_lossy().trim().is_empty())
        .collect()
}

/// Parse the cell tag for template-coordinate sort/verify, returning `None`
/// for other sort orders.
pub(crate) fn parse_cell_tag(order: SortOrderArg) -> Result<Option<SamTag>> {
    if matches!(order, SortOrderArg::TemplateCoordinate) { Ok(Some(SamTag::CB)) } else { Ok(None) }
}

/// Parse `--key-types` into a [`KeyTypesSpec`].
///
/// `full` | `none`/`""` | comma/space list of `cb`,`library`,`mi`. Unknown
/// tokens are rejected. `library`, `mi`, and `library,mi` all map to the shared
/// `tertiary` lane.
pub(crate) fn parse_key_types(s: &str) -> Result<KeyTypesSpec, String> {
    let s = s.trim();
    if s.eq_ignore_ascii_case("full") {
        return Ok(KeyTypesSpec::Full);
    }
    if s.is_empty() || s.eq_ignore_ascii_case("none") {
        return Ok(KeyTypesSpec::None);
    }
    let mut cb = false;
    let mut tertiary = false;
    for tok in s.split([',', ' ']).filter(|t| !t.is_empty()) {
        match tok.to_ascii_lowercase().as_str() {
            "cb" => cb = true,
            "library" | "mi" => tertiary = true,
            other => {
                return Err(format!(
                    "unknown --key-types token '{other}', expected 'full', 'none', \
                     or a list of 'cb','library','mi'"
                ));
            }
        }
    }
    Ok(KeyTypesSpec::Explicit { cb, tertiary })
}

impl Command for Sort {
    fn execute(&self, command_line: &str) -> Result<()> {
        if self.verify && self.output.is_some() {
            bail!("--verify cannot be used with --output");
        }
        if self.verify && self.write_index {
            bail!("--write-index cannot be used with --verify");
        }
        // A BAI is written to a sidecar path beside a seekable file, and a pipe
        // is neither. Rejected up front rather than at the writer, so the run
        // fails before streaming a BAM it could never index.
        if self.write_index && self.output.as_deref().is_some_and(is_stdout_path) {
            bail!(
                "--write-index cannot be used with output to stdout: an index requires a \
                 seekable file, so give --output a path"
            );
        }

        // Validate inputs. Exempt stdin paths (`-` / `/dev/stdin`): the sort
        // engine reads stdin in a single pass (a `TeeReader` reads the header,
        // then a `ChainedReader` replays it and streams the rest), so a
        // file-existence check would spuriously reject it — matching the stdin
        // exemption in group/dedup/clip. `--verify` takes the same single-pass
        // path, so it is not an exception.
        validate_input_exists(&self.input, "Input BAM")?;

        // Either --output or --verify must be specified
        if !self.verify && self.output.is_none() {
            bail!("Either --output or --verify must be specified");
        }

        if self.verify {
            return self.execute_verify();
        }

        self.execute_sort(command_line)
    }
}

impl Sort {
    /// Thread count that `--memory-per-thread` multiplies `--max-memory` by.
    ///
    /// The budget sizes the in-memory accumulation buffer, which the sort phase
    /// fills, so `--sort-threads` is the count that should drive it. It is
    /// combined with `--threads` rather than replacing it: lowering only the sort
    /// phase (`-@ 32 --sort-threads 4`) is the documented way to cede cores to an
    /// upstream producer, and letting that shrink the buffer 8x would turn a
    /// scheduling hint into a throughput cliff. Taking the larger of the two
    /// raises the budget when the sort phase is the wider one -- the case where
    /// `--threads` alone under-counts -- and never lowers it.
    fn memory_budget_threads(&self) -> usize {
        // Single source of truth shared with the chain builder's `add_sort`
        // (`sort_budget_threads`) — see `common::sort_memory_budget_threads` — so
        // the banner path and the chain path cannot drift. `--threads 0` still
        // resolves to 0 there, so `resolve_memory_budget` rejects a zero-thread run.
        crate::commands::common::sort_memory_budget_threads(self.threads, self.sort_threads)
    }

    /// The flag [`memory_budget_threads`](Self::memory_budget_threads) took its
    /// count from, for the `Max memory:` log line.
    ///
    /// `--threads` is the floor, so it is the source unless `--sort-threads` was
    /// set strictly above it. Reported so that line and the per-phase `Threads:`
    /// line below it cannot be read as disagreeing.
    fn memory_budget_threads_flag(&self) -> &'static str {
        if self.sort_threads.is_some_and(|n| n > self.threads) {
            "--sort-threads"
        } else {
            "--threads"
        }
    }

    /// Effective Phase-1 (accumulate/sort/spill) worker count.
    ///
    /// Mirrors `RawExternalSorter::phase1_threads` in
    /// `fgumi-sort/src/external.rs` exactly. Duplicated rather than shared:
    /// the engine's own version is used only internally by
    /// `SortBuffer::from_sorter` after this method replaces the CLI's only
    /// other caller, so there is no third call site to justify a cross-crate
    /// export. Kept honest by `test_phase_threads_match_sorters_formula`
    /// (values) and `test_phase_threads_match_engine` (direct parity with
    /// the engine's formula) here.
    fn phase1_threads(&self) -> usize {
        self.sort_threads.unwrap_or(self.threads).max(1)
    }

    /// Effective Phase-2 (merge/write) worker count. See `phase1_threads`.
    fn phase2_threads(&self) -> usize {
        self.merge_threads.unwrap_or(self.threads).max(1)
    }

    /// The spill-file consolidation limit this command will use.
    ///
    /// The engine's own default is a fixed, portable value; sizing it to the
    /// host's descriptor budget is a decision the command line makes, so it is
    /// resolved here and passed down. Resolving in one place keeps the number
    /// that gets logged identical to the number the engine is handed.
    fn resolved_max_temp_files(&self, soft_nofile: Option<u64>) -> usize {
        match self.max_temp_files {
            MaxTempFiles::Auto => fgumi_sort::temp_file_limit_from_nofile(soft_nofile),
            MaxTempFiles::Fixed(limit) => limit,
        }
    }

    /// Parse the cell tag for template-coordinate sort/verify, returning `None`
    /// for other sort orders.
    fn parse_cell_tag(&self) -> Result<Option<SamTag>> {
        parse_cell_tag(self.order)
    }
}

/// The line naming the spill-file limit in force, for the run configuration log.
///
/// Returns the string rather than logging it so the wording — in particular
/// whether the number was asked for, derived from a budget this line then names,
/// or fell back to the portable default — is assertable without a log capture.
fn max_temp_files_log_line(
    setting: MaxTempFiles,
    resolved: usize,
    soft_nofile: Option<u64>,
) -> String {
    match (setting, soft_nofile) {
        (MaxTempFiles::Fixed(_), _) => format!("Max temp files: {resolved}"),
        (MaxTempFiles::Auto, Some(soft)) => {
            format!("Max temp files: {resolved} (derived from RLIMIT_NOFILE soft limit {soft})")
        }
        (MaxTempFiles::Auto, None) => format!("Max temp files: {resolved} (default)"),
    }
}

/// The warning for a spill-file limit that overruns the process's descriptor
/// budget, or `None` when `resolved` fits.
///
/// Checked for explicit limits too, not just derived ones. A derived limit is
/// clamped to the budget and can only fall short of it when the budget is below
/// the floor; an explicit limit can overrun it outright, which makes it the case
/// more likely to fail — so the advice names `--max-temp-files` only when the
/// user actually set it, and otherwise points at `ulimit -n`, the only lever
/// left.
fn fd_budget_warning(
    setting: MaxTempFiles,
    resolved: usize,
    soft_nofile: Option<u64>,
) -> Option<String> {
    if fgumi_sort::fits_nofile_budget(soft_nofile, resolved) {
        return None;
    }
    let advice = if matches!(setting, MaxTempFiles::Fixed(_)) {
        "lower --max-temp-files or raise the open file limit with `ulimit -n`"
    } else {
        "raise the open file limit with `ulimit -n`"
    };
    Some(format!(
        "A spill-file limit of {resolved} exceeds this process's open file budget; \
         {advice} if the sort fails with \"Too many open files\"."
    ))
}

impl Sort {
    /// Execute sort mode: read, sort, and write output.
    fn execute_sort(&self, command_line: &str) -> Result<()> {
        let output = self.output.as_ref().expect("output required for sort mode");

        // --write-index only valid for coordinate sort
        if self.write_index && !matches!(self.order, SortOrderArg::Coordinate) {
            bail!("--write-index is only valid for coordinate sort");
        }

        // zstd has no level-0 "stored" mode; silently remapping to 1 would
        // surprise users who pass --temp-compression 0 to disable temp
        // compression (which works for BGZF). Reject the combination
        // explicitly.
        if self.temp_compression == 0 && matches!(self.temp_codec, fgumi_sort::SpillCodec::Zstd) {
            bail!(
                "--temp-compression 0 is only supported with --temp-codec bgzf; \
                 zstd does not have an uncompressed mode. Pass --temp-codec bgzf \
                 to keep level-0 spill, or pick a zstd level >= 1."
            );
        }

        // The "Sorting BAM ..." start line and its completion line are both owned
        // by the chain's `SortSummaryFinalizeHook` timer (`add_sort` constructs an
        // `OperationTimer::new("Sorting BAM")` there). execute_sort no longer
        // builds its own timer — doing so logged the identical start line twice.

        // Resolve memory limit (auto-detect or fixed)
        let budget_threads = self.memory_budget_threads();
        let effective_memory = resolve_memory_budget(
            self.max_memory,
            self.memory_reserve,
            budget_threads,
            self.memory_per_thread,
        )?;

        let cell_tag = self.parse_cell_tag()?;

        // One reading of the descriptor budget serves the whole run: the limit
        // handed to the sorter, the budget named in the log line, and the check
        // the warning makes all derive from this value, so they cannot describe
        // different limits if `RLIMIT_NOFILE` changes mid-run.
        let soft_nofile = fgumi_sort::soft_nofile();

        debug!("Starting Sort");
        info!("Input: {}", self.input.display());
        info!("Output: {}", output.display());
        info!("Sort order: {:?}", self.order);
        if let Some(ct) = cell_tag {
            let ct_bytes = *ct;
            info!("Cell tag: {}{}", ct_bytes[0] as char, ct_bytes[1] as char);
        }
        if let MemoryLimit::Fixed(per_thread) = self.max_memory {
            if self.memory_per_thread {
                // `resolve_memory_budget` sizes one global budget for the whole
                // run, and the multiplier is the larger of `--threads` and
                // `--sort-threads` -- never `--merge-threads`. Name the flag it
                // came from so this line and the per-phase counts logged below
                // cannot be read as disagreeing.
                info!(
                    "Max memory: {} ({}/thread x {} threads, from {})",
                    ByteSize(effective_memory as u64),
                    ByteSize(per_thread as u64),
                    budget_threads,
                    self.memory_budget_threads_flag()
                );
            } else {
                info!("Max memory: {} (fixed)", ByteSize(effective_memory as u64));
            }
        }
        // Report the effective per-phase thread counts (not the raw --threads,
        // which is only where --sort-threads/--merge-threads default from):
        // logging the flag alone reports 1 thread for a run that asked for more.
        info!(
            "Threads: {}",
            fgumi_sort::format_thread_counts(self.phase1_threads(), self.phase2_threads())
        );
        info!("Temp compression level: {}", self.temp_compression);
        // The resolved max-temp-files the sort will actually consolidate at,
        // sourced the same way as the thread counts above so the banner cannot
        // drift from the value handed to the sort.
        let max_temp_files = self.resolved_max_temp_files(soft_nofile);
        info!("{}", max_temp_files_log_line(self.max_temp_files, max_temp_files, soft_nofile));
        if let Some(warning) = fd_budget_warning(self.max_temp_files, max_temp_files, soft_nofile) {
            warn!("{warning}");
        }
        if self.write_index {
            info!("Write index: enabled");
        }
        let env_value = std::env::var(TMP_DIRS_ENV).ok();
        let resolved_tmp_dirs = resolve_tmp_dirs(&self.tmp_dirs, env_value.as_deref());
        if !resolved_tmp_dirs.is_empty() {
            let joined = resolved_tmp_dirs
                .iter()
                .map(|p| p.display().to_string())
                .collect::<Vec<_>>()
                .join(", ");
            info!("Temp directories: {joined}");
        }

        // `--read-streams` is now threaded end-to-end (into the chain's BAM
        // source below); no warn shim. `--sort-stats` is still unthreaded.
        if self.sort_stats {
            warn!("--sort-stats is ignored by the sort chain (not yet threaded)");
        }

        // --- cutover: run via the chain instead of sorter.sort() ---
        // (`phase1_threads`/`phase2_threads`/`resolved_max_temp_files` above are retained
        //  only to source the banner's thread/temp-file numbers; the owned sorter is not
        //  constructed or executed here. PR B removes the owned engine.)
        //
        // `output` is already bound at the top of this method (validated in `execute`),
        // and `resolved_tmp_dirs` was resolved above for the banner; hand both to the
        // (unit-tested) spec builder. `resolved_tmp_dirs` is moved here — its last
        // borrow was the banner log above.
        let spec = self.build_sort_chain_spec(output, resolved_tmp_dirs, command_line);
        build_for(spec)?.run()
    }

    /// Execute verify mode: read records and check sort order.
    fn execute_verify(&self) -> Result<()> {
        use fgumi_sort::open_raw_bam_record_reader_with_header;
        use fgumi_sort::{
            LibraryLookup, RawQuerynameKey, RawQuerynameLexKey, RawSortKey, SortContext, cb_hasher,
            extract_coordinate_key_inline, extract_template_key_inline,
        };
        use std::cmp::Ordering;

        let cell_tag = self.parse_cell_tag()?;

        let timer = OperationTimer::new("Verifying BAM sort order");

        debug!("Starting Sort Verification");
        info!("Input: {}", self.input.display());
        info!("Expected order: {:?}", self.order);
        if let Some(ct) = cell_tag {
            let ct_bytes = *ct;
            info!("Cell tag: {}{}", ct_bytes[0] as char, ct_bytes[1] as char);
        }

        // One open yields both the header and the records: the header is parsed
        // through a tee and the consumed bytes replayed. Reading the path twice
        // is what used to make `--verify` reject a pipe.
        let (raw_reader, header) = open_raw_bam_record_reader_with_header(&self.input)?;

        let (total_records, violations, first_violation) = match self.order {
            SortOrderArg::Coordinate => {
                let nref = header.reference_sequences().len() as u32;
                verify_sort_order(
                    raw_reader,
                    |bam| extract_coordinate_key_inline(bam, nref),
                    |key, prev| key < prev,
                )?
            }
            SortOrderArg::Queryname => {
                let ctx = SortContext::from_header(&header);
                verify_sort_order(
                    raw_reader,
                    |bam| RawQuerynameLexKey::extract(bam, &ctx),
                    |key, prev| key < prev,
                )?
            }
            SortOrderArg::QuerynameNatural => {
                let ctx = SortContext::from_header(&header);
                verify_sort_order(
                    raw_reader,
                    |bam| RawQuerynameKey::extract(bam, &ctx),
                    |key, prev| key < prev,
                )?
            }
            SortOrderArg::TemplateCoordinate => {
                let lib_lookup = LibraryLookup::from_header(&header);
                let hasher = cb_hasher();
                verify_sort_order(
                    raw_reader,
                    |bam| extract_template_key_inline(bam, &lib_lookup, cell_tag, &hasher),
                    // Use core_cmp to ignore name_hash tie-breaker differences
                    // This allows both fgumi and samtools sorted files to pass
                    |key, prev| key.core_cmp(prev) == Ordering::Less,
                )?
            }
        };

        // Summary
        info!("=== Verification Summary ===");
        info!("Records checked: {total_records}");
        info!("Sort order violations: {violations}");

        if violations > 0 {
            if let Some((record_num, name)) = first_violation {
                info!("First violation at record {record_num}: {name}");
            }
            timer.log_completion(total_records);
            bail!(
                "BAM file is NOT correctly sorted by {:?}: {violations} violations found",
                self.order
            );
        }

        info!("Result: PASS - file is correctly sorted by {:?}", self.order);
        timer.log_completion(total_records);
        Ok(())
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    // Memory-budget helpers moved to `commands::common`; import the `pub(crate)`
    // items these tests exercise that are not re-exported through `super::*`.
    use crate::commands::common::{MIN_MEMORY_PER_THREAD, detect_total_memory, resolve_reserve};
    use clap::{CommandFactory, Parser};
    use rstest::rstest;

    // ========================================================================
    // Help-text tests
    // ========================================================================

    /// Returns true if `help` names the `fgumi` binary as a standalone token.
    ///
    /// Splits on every character that cannot appear inside a Rust identifier, so
    /// a bare `fgumi`, a trailing `fgumi`, and punctuated forms (`` `fgumi` ``,
    /// `fgumi.`, `fgumi-sort`) all match, while `fgumi` embedded in a longer word
    /// (`fgumidocs`) and identifier-shaped mentions (`fgumi_sort`) do not.
    fn names_the_binary(help: &str) -> bool {
        help.split(|c: char| !c.is_ascii_alphanumeric() && c != '_').any(|token| token == "fgumi")
    }

    #[rstest]
    #[case::bare_invocation("run fgumi sort to order records", true)]
    #[case::trailing_token("this flag mirrors fgumi", true)]
    #[case::followed_by_period("see fgumi.", true)]
    #[case::backticked("see `fgumi` for details", true)]
    #[case::hyphenated("see fgumi-sort", true)]
    #[case::embedded_in_word("see the fgumidocs site", false)]
    #[case::suffix_of_word("see myfgumi", false)]
    #[case::rust_identifier("handled by fgumi_sort::run", false)]
    #[case::no_mention("Sort records by template coordinate", false)]
    fn test_names_the_binary(#[case] help: &str, #[case] expected: bool) {
        assert_eq!(names_the_binary(help), expected, "help was: {help}");
    }

    /// `Sort` is `#[command(flatten)]`-ed into other binaries, so its
    /// per-argument help renders under a program name that is not `fgumi` — a
    /// concrete `fgumi ...` invocation there tells those users to run a command
    /// they do not have. Per-argument help must therefore describe the flag
    /// without naming the binary.
    ///
    /// The command-level `long_about` is deliberately exempt, and is not walked
    /// here: a wrapper replaces it wholesale, so its EXAMPLES block is the right
    /// home for worked `fgumi sort` invocations.
    #[test]
    fn test_arg_help_does_not_name_the_binary() {
        let command = Sort::command();

        // Guard against a vacuous pass: an empty (or trivially short) arg walk
        // would satisfy the loop below without checking anything.
        let arg_count = command.get_arguments().count();
        assert!(arg_count > 10, "expected Sort to expose its flags, walked only {arg_count}");

        for arg in command.get_arguments() {
            let help = format!(
                "{} {}",
                arg.get_help().map(ToString::to_string).unwrap_or_default(),
                arg.get_long_help().map(ToString::to_string).unwrap_or_default(),
            );
            assert!(
                !names_the_binary(&help),
                "help for `--{}` names the `fgumi` binary; describe the flag instead and put \
                 worked invocations in the command-level EXAMPLES block. Help was: {help}",
                arg.get_id()
            );
        }
    }

    /// The advanced performance-instrumentation flags are hidden from generated
    /// help. The performance-tuning guide documents them as available but hidden,
    /// and surfacing them in `fgumi sort --help` would bury the handful of
    /// options a normal run cares about. They stay parseable; only their help
    /// rendering is suppressed.
    #[rstest]
    #[case::sort_stats("sort-stats")]
    fn test_advanced_flags_are_hidden_from_help(#[case] long: &str) {
        let command = Sort::command();
        let arg = command
            .get_arguments()
            .find(|arg| arg.get_long() == Some(long))
            .unwrap_or_else(|| panic!("`--{long}` not found on `fgumi sort`"));
        assert!(arg.is_hide_set(), "`--{long}` must be hidden from generated help");
    }

    // Boolean-flag help is guarded repo-wide in `src/main.rs`'s `mod tests`,
    // which walks every command rather than `Sort` alone and identifies boolean
    // flags by their parser's accepted set rather than by `num_args`.

    // ========================================================================
    // Temp-dir resolution tests
    // ========================================================================

    #[test]
    fn test_resolve_tmp_dirs_empty() {
        assert!(resolve_tmp_dirs(&[], None).is_empty());
        assert!(resolve_tmp_dirs(&[], Some("")).is_empty());
    }

    #[test]
    fn test_resolve_tmp_dirs_cli_only() {
        let cli = vec![PathBuf::from("/tmp/a"), PathBuf::from("/tmp/b")];
        let got = resolve_tmp_dirs(&cli, None);
        assert_eq!(got, cli);
    }

    #[test]
    fn test_resolve_tmp_dirs_env_only() {
        #[cfg(unix)]
        let env = "/tmp/x:/tmp/y";
        #[cfg(windows)]
        let env = "C:/tmp/x;C:/tmp/y";

        let got = resolve_tmp_dirs(&[], Some(env));
        assert_eq!(got.len(), 2);
        assert!(got[0].to_string_lossy().ends_with('x'));
        assert!(got[1].to_string_lossy().ends_with('y'));
    }

    #[test]
    fn test_resolve_tmp_dirs_cli_overrides_env() {
        let cli = vec![PathBuf::from("/tmp/cli")];
        #[cfg(unix)]
        let env = "/tmp/env1:/tmp/env2";
        #[cfg(windows)]
        let env = "C:/tmp/env1;C:/tmp/env2";

        let got = resolve_tmp_dirs(&cli, Some(env));
        assert_eq!(got, cli, "CLI flags must take precedence over env var");
    }

    #[test]
    fn test_resolve_tmp_dirs_skips_empty_segments() {
        // Trailing separator / empty segment must not produce a bogus empty PathBuf.
        #[cfg(unix)]
        let env = "/tmp/a::/tmp/b:";
        #[cfg(windows)]
        let env = "C:/tmp/a;;C:/tmp/b;";

        let got = resolve_tmp_dirs(&[], Some(env));
        assert_eq!(got.len(), 2, "empty path segments must be filtered: {got:?}");
    }

    // ========================================================================
    // Clap parsing tests for repeatable -T flag
    // ========================================================================

    #[rstest]
    #[case::zero(&[], vec![])]
    #[case::single_short(&["-T", "/tmp/a"], vec![PathBuf::from("/tmp/a")])]
    #[case::multiple_short(
        &["-T", "/tmp/a", "-T", "/tmp/b", "-T", "/tmp/c"],
        vec![PathBuf::from("/tmp/a"), PathBuf::from("/tmp/b"), PathBuf::from("/tmp/c")],
    )]
    #[case::multiple_long(
        &["--tmp-dir", "/tmp/a", "--tmp-dir", "/tmp/b"],
        vec![PathBuf::from("/tmp/a"), PathBuf::from("/tmp/b")],
    )]
    fn test_clap_tmp_dir_repeatable(#[case] extra: &[&str], #[case] expected: Vec<PathBuf>) {
        let base = ["sort", "-i", "in.bam", "-o", "out.bam", "--order", "coordinate"];
        let args: Vec<&str> = base.iter().copied().chain(extra.iter().copied()).collect();
        let sort = Sort::try_parse_from(args).expect("parse should succeed");
        assert_eq!(sort.tmp_dirs, expected);
    }

    /// Pin the `--temp-codec` parser default and explicit override so a flipped
    /// default (or a removed `default_value`) fails at unit-test time.
    #[rstest]
    #[case::default_omitted(&[], fgumi_sort::SpillCodec::Zstd)]
    #[case::explicit_zstd(&["--temp-codec", "zstd"], fgumi_sort::SpillCodec::Zstd)]
    #[case::explicit_bgzf(&["--temp-codec", "bgzf"], fgumi_sort::SpillCodec::Bgzf)]
    fn test_clap_temp_codec_default(
        #[case] extra: &[&str],
        #[case] expected: fgumi_sort::SpillCodec,
    ) {
        let base = ["sort", "-i", "in.bam", "-o", "out.bam", "--order", "coordinate"];
        let args: Vec<&str> = base.iter().copied().chain(extra.iter().copied()).collect();
        let sort = Sort::try_parse_from(args).expect("parse should succeed");
        assert_eq!(sort.temp_codec, expected);
    }

    /// The per-phase flags must parse and reach the command struct, defaulting
    /// to `None` (which the builder resolves to `--threads`) when not passed.
    #[rstest]
    #[case::neither(&[], None, None)]
    #[case::sort_only(&["--sort-threads", "4"], Some(4), None)]
    #[case::merge_only(&["--merge-threads", "3"], None, Some(3))]
    #[case::both(&["--sort-threads", "2", "--merge-threads", "8"], Some(2), Some(8))]
    fn test_parse_per_phase_threads(
        #[case] extra: &[&str],
        #[case] expected_sort: Option<usize>,
        #[case] expected_merge: Option<usize>,
    ) {
        let base = ["sort", "-i", "in.bam", "-o", "out.bam", "--order", "coordinate"];
        let args: Vec<&str> = base.iter().copied().chain(extra.iter().copied()).collect();
        let sort = Sort::try_parse_from(args).expect("parse should succeed");
        assert_eq!(sort.sort_threads, expected_sort);
        assert_eq!(sort.merge_threads, expected_merge);
    }

    /// The single `--max-memory` knob must drive the pipeline's queue budget too,
    /// not just the sorter buffer. The helper projects the command's three memory
    /// flags 1:1 onto `QueueMemoryOptions`, so a run's queue budget agrees with its
    /// sorter budget on those fields (one flag, both budgets). Assert against the
    /// parsed command fields rather than hardcoded byte counts — the byte math is
    /// the memory-size parser's contract, not this projection's.
    #[rstest]
    #[case::defaults(&[])]
    #[case::small_fixed(&["-m", "64KiB"])]
    #[case::auto_with_reserve(&["-m", "auto", "--memory-reserve", "8GiB"])]
    #[case::fixed_total(&["-m", "2GiB", "--memory-per-thread", "false"])]
    fn test_queue_memory_options_projects_command_flags(#[case] extra: &[&str]) {
        let base = ["sort", "-i", "in.bam", "-o", "out.bam", "--order", "coordinate"];
        let args: Vec<&str> = base.iter().copied().chain(extra.iter().copied()).collect();
        let sort = Sort::try_parse_from(args).expect("parse should succeed");

        let queue = sort.queue_memory_options();
        assert_eq!(queue.max_memory, sort.max_memory);
        assert_eq!(queue.memory_reserve, sort.memory_reserve);
        assert_eq!(queue.memory_per_thread, sort.memory_per_thread);

        // One flag drives both budgets: the queue options must agree with the
        // sorter options on the three shared memory fields.
        let sorter = sort.to_sort_options();
        assert_eq!(queue.max_memory, sorter.max_memory);
        assert_eq!(queue.memory_reserve, sorter.memory_reserve);
        assert_eq!(queue.memory_per_thread, sorter.memory_per_thread);
    }

    /// Regression guard for the exact bug this follow-up fixes: before threading
    /// `--max-memory` into `queue_memory`, the chain hardcoded
    /// `QueueMemoryOptions::default()` (768 MiB/thread), so a small `--max-memory`
    /// shrank the sorter buffer while the queues stayed large. A small explicit
    /// value must now move the queue budget off that default.
    #[test]
    fn small_max_memory_moves_queue_budget_off_default() {
        let sort = Sort::try_parse_from([
            "sort",
            "-i",
            "in.bam",
            "-o",
            "out.bam",
            "--order",
            "coordinate",
            "-m",
            "64KiB",
        ])
        .expect("parse should succeed");

        let queue = sort.queue_memory_options();
        assert_ne!(
            queue.max_memory,
            QueueMemoryOptions::default().max_memory,
            "a small --max-memory must move the queue budget off the 768 MiB default",
        );
        assert_eq!(queue.max_memory, sort.max_memory);
    }

    /// The projection helper being correct is not enough: the whole point of this
    /// change is line-of-code wiring — that `execute_sort`'s `ChainSpec` actually
    /// carries `queue_memory_options()` (and the other flag-derived fields). That
    /// wiring is invisible to a full run (the queue budget does not change output
    /// bytes) and to the isolated helper tests, so reverting it to
    /// `QueueMemoryOptions::default()` would keep every other test green. Assert on
    /// the built spec directly. `QueueMemoryOptions` has no `PartialEq`, so compare
    /// its fields rather than the whole struct.
    #[rstest]
    #[case::plain(false)]
    #[case::write_index(true)]
    fn build_sort_chain_spec_wires_flags_into_the_spec(#[case] write_index: bool) {
        let mut args = vec![
            "sort",
            "-i",
            "in.bam",
            "-o",
            "out.bam",
            "--order",
            "coordinate",
            "-m",
            "64KiB",
            "-@",
            "3",
        ];
        if write_index {
            args.push("--write-index");
        }
        let sort = Sort::try_parse_from(args).expect("parse should succeed");

        let spec =
            sort.build_sort_chain_spec(Path::new("out.bam"), Vec::new(), "fgumi sort (test)");

        // The queue budget must carry the parsed --max-memory, not the default.
        assert_eq!(spec.queue_memory.max_memory, sort.max_memory);
        assert_eq!(spec.queue_memory.memory_reserve, sort.memory_reserve);
        assert_eq!(spec.queue_memory.memory_per_thread, sort.memory_per_thread);
        assert_ne!(
            spec.queue_memory.max_memory,
            QueueMemoryOptions::default().max_memory,
            "the spec must carry the parsed --max-memory, proving the wiring is live",
        );

        // The other flag-derived spec fields the run relies on.
        assert_eq!(spec.stages, vec![Stage::Sort]);
        assert_eq!(spec.threading.threads, Some(sort.threads));
        assert!(matches!(spec.source, SourceSpec::Bam(_)));
        assert!(spec.verify_crc);
        // --read-streams reaches the chain's BAM source (defaulted to Auto here),
        // proving the field is wired even without an explicit flag.
        assert_eq!(spec.read_streams, sort.read_streams);
        match (write_index, &spec.sink) {
            (true, SinkSpec::BamWithIndex(_)) | (false, SinkSpec::Bam(_)) => {}
            (wi, other) => panic!("write_index={wi} produced the wrong sink: {other:?}"),
        }
    }

    /// The sorter-free phase thread helpers `Sort::phase1_threads` /
    /// `Sort::phase2_threads` must match the engine's own formula, since
    /// they replace the CLI's only other caller of `RawExternalSorter`'s
    /// methods.
    #[rstest]
    #[case::defaults(None, None, 2, 2)]
    #[case::narrow_sort(Some(1), Some(4), 1, 4)]
    #[case::narrow_merge(Some(4), Some(1), 4, 1)]
    #[case::sort_only(Some(3), None, 3, 2)]
    #[case::merge_only(None, Some(3), 2, 3)]
    fn test_phase_threads_match_sorters_formula(
        #[case] sort_threads: Option<usize>,
        #[case] merge_threads: Option<usize>,
        #[case] expected_phase1: usize,
        #[case] expected_phase2: usize,
    ) {
        let mut sort = make_sort(SortOrderArg::Coordinate);
        sort.threads = 2;
        sort.sort_threads = sort_threads;
        sort.merge_threads = merge_threads;

        assert_eq!(sort.phase1_threads(), expected_phase1);
        assert_eq!(sort.phase2_threads(), expected_phase2);
    }

    /// The CLI's `phase1_threads`/`phase2_threads` are hand-copied from
    /// `RawExternalSorter`'s formula. Pin them against the engine directly so
    /// the two copies cannot drift silently (both would otherwise stay green
    /// against their own tables).
    #[rstest]
    #[case::defaults(2, None, None)]
    #[case::narrow_sort(2, Some(1), Some(4))]
    #[case::narrow_merge(2, Some(4), Some(1))]
    #[case::sort_only(2, Some(3), None)]
    #[case::merge_only(2, None, Some(3))]
    #[case::zero_sort_override_clamps(4, Some(0), None)]
    #[case::zero_merge_override_clamps(4, None, Some(0))]
    fn test_phase_threads_match_engine(
        #[case] threads: usize,
        #[case] sort_threads: Option<usize>,
        #[case] merge_threads: Option<usize>,
    ) {
        let mut sort = make_sort(SortOrderArg::Coordinate);
        sort.threads = threads;
        sort.sort_threads = sort_threads;
        sort.merge_threads = merge_threads;

        let mut engine = fgumi_sort::RawExternalSorter::new(sort.order.into()).threads(threads);
        if let Some(n) = sort_threads {
            engine = engine.sort_threads(n);
        }
        if let Some(n) = merge_threads {
            engine = engine.merge_threads(n);
        }
        assert_eq!(sort.phase1_threads(), engine.phase1_threads());
        assert_eq!(sort.phase2_threads(), engine.phase2_threads());
    }

    /// The per-thread budget sizes the buffer the sort phase fills, so it follows
    /// `--sort-threads` upward. `--sort-threads 8` with `--threads` unset used to
    /// resolve a one-thread budget for an eight-thread sort.
    ///
    /// It never follows `--sort-threads` downward: `-@ 32 --sort-threads 4` is the
    /// documented way to cede cores to an upstream producer, and must keep the
    /// 32-thread budget it resolves today.
    #[rstest]
    #[case::defaults(1, None, 1)]
    #[case::threads_only(4, None, 4)]
    #[case::sort_above_threads(1, Some(8), 8)]
    #[case::sort_below_threads(32, Some(4), 32)]
    #[case::sort_equals_threads(4, Some(4), 4)]
    #[case::zero_sort_override_keeps_threads(4, Some(0), 4)]
    #[case::zero_threads_is_not_clamped(0, Some(8), 0)]
    fn test_memory_budget_threads(
        #[case] threads: usize,
        #[case] sort_threads: Option<usize>,
        #[case] expected: usize,
    ) {
        let mut sort = make_sort(SortOrderArg::Coordinate);
        sort.threads = threads;
        sort.sort_threads = sort_threads;

        assert_eq!(sort.memory_budget_threads(), expected);
    }

    /// The `Max memory:` line names whichever flag supplied the multiplier.
    ///
    /// `--threads` is the floor, so it stays the attributed source whenever it
    /// ties or wins — including the `--sort-threads 0` case the engine clamps.
    #[rstest]
    #[case::defaults(1, None, "--threads")]
    #[case::threads_only(4, None, "--threads")]
    #[case::sort_above_threads(1, Some(8), "--sort-threads")]
    #[case::sort_below_threads(32, Some(4), "--threads")]
    #[case::sort_equals_threads(4, Some(4), "--threads")]
    #[case::zero_sort_override(4, Some(0), "--threads")]
    fn test_memory_budget_threads_flag(
        #[case] threads: usize,
        #[case] sort_threads: Option<usize>,
        #[case] expected: &str,
    ) {
        let mut sort = make_sort(SortOrderArg::Coordinate);
        sort.threads = threads;
        sort.sort_threads = sort_threads;

        assert_eq!(sort.memory_budget_threads_flag(), expected);
    }

    /// The budget the sort-phase count actually resolves to, end to end through
    /// `resolve_memory_budget`.
    ///
    /// The `Max memory:` log line is printed from the same local that feeds
    /// `resolve_memory_budget`, so the integration tests that assert on it cannot
    /// tell a correctly wired budget from one that logs `budget_threads` and
    /// resolves `self.threads`. This pins the resolved byte count instead.
    #[rstest]
    #[case::sort_above_threads(1, Some(8), 800_000_000)]
    #[case::sort_below_threads(4, Some(2), 400_000_000)]
    #[case::threads_only(4, None, 400_000_000)]
    #[case::zero_sort_override_keeps_threads(4, Some(0), 400_000_000)]
    fn test_memory_budget_threads_resolves_the_scaled_budget(
        #[case] threads: usize,
        #[case] sort_threads: Option<usize>,
        #[case] expected_bytes: usize,
    ) {
        let mut sort = make_sort(SortOrderArg::Coordinate);
        sort.threads = threads;
        sort.sort_threads = sort_threads;

        let resolved = resolve_memory_budget(
            MemoryLimit::Fixed(100_000_000),
            MemoryReserve::Auto,
            sort.memory_budget_threads(),
            true,
        )
        .expect("budget should resolve");

        assert_eq!(resolved, expected_bytes);
    }

    /// `--max-temp-files` parses onto the struct and defaults to `auto`, which
    /// is also spellable explicitly -- the point of the literal is that the
    /// host-derived behaviour has a name rather than being the absence of a flag.
    #[rstest]
    #[case::unset(&[], MaxTempFiles::Auto)]
    #[case::auto(&["--max-temp-files", "auto"], MaxTempFiles::Auto)]
    #[case::auto_mixed_case(&["--max-temp-files", "AUTO"], MaxTempFiles::Auto)]
    #[case::minimum(&["--max-temp-files", "2"], MaxTempFiles::Fixed(2))]
    #[case::explicit(&["--max-temp-files", "256"], MaxTempFiles::Fixed(256))]
    #[case::large(&["--max-temp-files", "100000"], MaxTempFiles::Fixed(100_000))]
    fn test_parse_max_temp_files(#[case] extra: &[&str], #[case] expected: MaxTempFiles) {
        let base = ["sort", "-i", "in.bam", "-o", "out.bam", "--order", "coordinate"];
        let args: Vec<&str> = base.iter().copied().chain(extra.iter().copied()).collect();
        let sort = Sort::try_parse_from(args).expect("parse should succeed");
        assert_eq!(sort.max_temp_files, expected);
    }

    /// The config log must say where the number came from, because the same
    /// value means different things depending on its provenance: an explicit
    /// limit is what the user asked for, a derived one is a function of a budget
    /// the line then names, and the fallback is neither.
    #[rstest]
    #[case::fixed_omits_provenance(
        MaxTempFiles::Fixed(256),
        256,
        Some(1024),
        "Max temp files: 256"
    )]
    #[case::fixed_ignores_the_soft_limit(
        MaxTempFiles::Fixed(256),
        256,
        None,
        "Max temp files: 256"
    )]
    #[case::auto_names_the_budget_it_came_from(
        MaxTempFiles::Auto,
        992,
        Some(1024),
        "Max temp files: 992 (derived from RLIMIT_NOFILE soft limit 1024)"
    )]
    #[case::auto_without_a_readable_budget_is_the_default(
        MaxTempFiles::Auto,
        64,
        None,
        "Max temp files: 64 (default)"
    )]
    fn test_max_temp_files_log_line(
        #[case] setting: MaxTempFiles,
        #[case] resolved: usize,
        #[case] soft_nofile: Option<u64>,
        #[case] expected: &str,
    ) {
        assert_eq!(max_temp_files_log_line(setting, resolved, soft_nofile), expected);
    }

    /// A limit that fits the descriptor budget must not warn, and one that
    /// cannot must name the lever the user actually has: `--max-temp-files` only
    /// when they set it, and `ulimit -n` either way.
    ///
    /// The budget is supplied rather than read from the host, so every branch
    /// is exercised on every target: a 1024-descriptor budget admits a
    /// merge-width limit of 8 and rejects `usize::MAX`, and a `None` budget
    /// (the non-Unix case) admits everything by design.
    #[rstest]
    #[case::fits_when_explicit(MaxTempFiles::Fixed(8), 8, Some(1024), None)]
    #[case::fits_when_derived(MaxTempFiles::Auto, 8, Some(1024), None)]
    #[case::explicit_overrun_points_at_the_flag(
        MaxTempFiles::Fixed(usize::MAX),
        usize::MAX,
        Some(1024),
        Some("lower --max-temp-files or raise the open file limit with `ulimit -n`")
    )]
    #[case::derived_overrun_points_only_at_ulimit(
        MaxTempFiles::Auto,
        usize::MAX,
        Some(1024),
        Some("raise the open file limit with `ulimit -n`")
    )]
    #[case::unreadable_budget_never_warns(MaxTempFiles::Fixed(usize::MAX), usize::MAX, None, None)]
    fn test_fd_budget_warning(
        #[case] setting: MaxTempFiles,
        #[case] resolved: usize,
        #[case] soft_nofile: Option<u64>,
        #[case] expected_advice: Option<&str>,
    ) {
        let warning = fd_budget_warning(setting, resolved, soft_nofile);
        match expected_advice {
            None => assert_eq!(warning, None, "a limit within the budget must not warn"),
            Some(advice) => {
                let warning = warning.expect("an oversized limit must warn");
                assert!(warning.contains(advice), "advice missing from: {warning}");
                assert!(
                    warning.contains(&format!("limit of {resolved}")),
                    "warning must name the limit: {warning}"
                );
            }
        }
    }

    /// The derived-overrun advice must not mention `--max-temp-files`: the user
    /// never set it, so telling them to lower it is advice they cannot act on.
    #[test]
    fn test_derived_overrun_advice_omits_the_flag_the_user_did_not_set() {
        let warning = fd_budget_warning(MaxTempFiles::Auto, usize::MAX, Some(1024))
            .expect("an oversized limit must warn");
        assert!(!warning.contains("lower --max-temp-files"), "unexpected advice: {warning}");
    }

    /// Invalid values are rejected at parse time rather than silently clamped or
    /// truncated: values below the floor of 2 (a merge needs at least two
    /// inputs), and non-numeric, negative, or overflowing input.
    #[rstest]
    #[case::zero("0")]
    #[case::one("1")]
    #[case::non_numeric("abc")]
    #[case::negative("-1")]
    #[case::overflow("99999999999999999999999999")]
    fn test_parse_max_temp_files_rejects_invalid(#[case] value: &str) {
        let args = [
            "sort",
            "-i",
            "in.bam",
            "-o",
            "out.bam",
            "--order",
            "coordinate",
            "--max-temp-files",
            value,
        ];
        assert!(Sort::try_parse_from(args).is_err(), "--max-temp-files {value} must be rejected");
    }

    /// An explicit limit passes through untouched; `auto` derives one from the
    /// descriptor budget it is handed. The budget is supplied rather than read
    /// from the host, so the expected derived value is a fixed number instead of
    /// one that depends on the machine the test runs on — and it is deliberately
    /// *not* the engine's own default, which stays a fixed portable value now
    /// that the derivation lives in the command.
    #[rstest]
    #[case::set(MaxTempFiles::Fixed(256), Some(1024), 256)]
    #[case::auto_derives_from_the_budget(MaxTempFiles::Auto, Some(1024), 992)]
    #[case::auto_without_a_budget_falls_back(
        MaxTempFiles::Auto,
        None,
        fgumi_sort::FALLBACK_MAX_TEMP_FILES
    )]
    fn test_resolved_max_temp_files(
        #[case] max_temp_files: MaxTempFiles,
        #[case] soft_nofile: Option<u64>,
        #[case] expected: usize,
    ) {
        let mut sort = make_sort(SortOrderArg::Coordinate);
        sort.max_temp_files = max_temp_files;
        assert_eq!(sort.resolved_max_temp_files(soft_nofile), expected);
    }

    /// The accessor test above only proves the flag reaches the builder. This
    /// one proves `--max-temp-files` is purely a *how-we-consolidate* knob and
    /// never changes *what* is written: with a 1 KiB memory budget (no floor in
    /// `Fixed` mode) and enough records to spill many runs, `Some(2)` forces the
    /// consolidation merge to fire repeatedly (oldest runs are merged once their
    /// count reaches the limit), while a limit far above the spill count never
    /// consolidates and merges every run in the final k-way pass. Both paths
    /// must emit the same records in the same order — this guards the
    /// consolidation path against a stable-order or record-loss regression the
    /// accessor test cannot catch.
    ///
    /// The control arm pins an explicit large limit rather than passing `None`.
    /// The engine default is derived from the host's `ulimit -n`, so a `None`
    /// control arm could itself consolidate on a low-limit machine, quietly
    /// turning this into a consolidation-vs-consolidation comparison that still
    /// passes while no longer testing what it claims.
    ///
    /// We compare decoded records rather than the raw BGZF bytes on purpose:
    /// the two write paths flush BGZF blocks at different points, so the
    /// compressed framing legitimately differs even when the record stream is
    /// identical. The contract consolidation must uphold is record identity,
    /// not byte identity of the container.
    #[test]
    fn test_max_temp_files_consolidation_is_record_identical() {
        use fgumi_sam::SamBuilder;

        // Reverse-ordered, all-distinct starts so the coordinate sort is a total
        // order (no ties whose order could depend on merge grouping); 60 pairs
        // at a 1 KiB budget spill well over a dozen runs.
        let mut builder = SamBuilder::new();
        for i in 0..60 {
            let _ = builder
                .add_pair()
                .name(&format!("read{i:04}"))
                .start1((60 - i) * 100 + 1)
                .start2((60 - i) * 100 + 51)
                .build();
        }
        let dir = tempfile::tempdir().expect("tempdir");
        let input = dir.path().join("in.bam");
        builder.write_bam(&input).expect("write bam");

        // 1 KiB memory forces spilling; this mirrors the order/memory/threads/
        // max-temp-files fields the command's owned-engine wiring used to set,
        // so it exercises the same consolidation behavior.
        let run = |max_temp_files: MaxTempFiles, out: PathBuf| {
            let mut sort = make_sort(SortOrderArg::Coordinate);
            sort.input = input.clone();
            sort.output = Some(out.clone());
            sort.threads = 1;
            sort.max_temp_files = max_temp_files;

            let sorter = fgumi_sort::RawExternalSorter::new(sort.order.into())
                .memory_limit(1024)
                .threads(1)
                .max_temp_files(sort.resolved_max_temp_files(fgumi_sort::soft_nofile()));
            let stats = sorter.sort(&input, &out).expect("sort should succeed");

            let mut reader =
                noodles::bam::io::reader::Builder.build_from_path(&out).expect("open output");
            let header = reader.read_header().expect("read header");
            let records = reader
                .record_bufs(&header)
                .collect::<std::io::Result<Vec<_>>>()
                .expect("read records");
            (stats, records)
        };

        let (default_stats, default_records) =
            run(MaxTempFiles::Fixed(usize::MAX), dir.path().join("default.bam"));
        let (limited_stats, limited_records) =
            run(MaxTempFiles::Fixed(2), dir.path().join("limited.bam"));

        // Precondition: the tiny budget really did spill multiple runs, so the
        // `Some(2)` run exercised consolidation rather than a trivial in-memory
        // sort. Without this the identity assertion below could pass vacuously.
        assert!(
            default_stats.runs_written >= 2,
            "test must spill multiple runs to exercise consolidation, got {} run(s)",
            default_stats.runs_written,
        );
        assert_eq!(
            limited_stats.output_records, default_stats.output_records,
            "consolidation must not drop or duplicate records",
        );
        assert_eq!(
            limited_records, default_records,
            "--max-temp-files is a consolidation knob and must not change the emitted records",
        );
    }

    /// Separately, the knob must not change what is written. Covers the same
    /// asymmetric splits end-to-end through `execute`.
    #[rstest]
    #[case::defaults(None, None)]
    #[case::narrow_sort(Some(1), Some(4))]
    #[case::narrow_merge(Some(4), Some(1))]
    fn test_execute_with_per_phase_threads_is_output_identical(
        #[case] sort_threads: Option<usize>,
        #[case] merge_threads: Option<usize>,
    ) {
        use fgumi_sam::SamBuilder;

        let mut builder = SamBuilder::new();
        for i in 0..60 {
            let _ = builder
                .add_pair()
                .name(&format!("read{i:04}"))
                .start1((60 - i) * 100 + 1)
                .start2((60 - i) * 100 + 51)
                .build();
        }
        let dir = tempfile::tempdir().expect("tempdir");
        let input = dir.path().join("in.bam");
        builder.write_bam(&input).expect("write bam");

        let run = |st: Option<usize>, mt: Option<usize>, out: PathBuf| {
            let mut sort = make_sort(SortOrderArg::Coordinate);
            sort.input = input.clone();
            sort.output = Some(out.clone());
            sort.threads = 2;
            sort.sort_threads = st;
            sort.merge_threads = mt;
            sort.execute("fgumi sort (test)").expect("sort should succeed");
            std::fs::read(&out).expect("read output")
        };

        let baseline = run(None, None, dir.path().join("baseline.bam"));
        let actual = run(sort_threads, merge_threads, dir.path().join("actual.bam"));
        assert_eq!(
            actual, baseline,
            "per-phase threads must not change output bytes (sort={sort_threads:?}, \
             merge={merge_threads:?})",
        );
    }

    /// Helper to construct a `Sort` struct with a given order.
    fn make_sort(order: SortOrderArg) -> Sort {
        Sort {
            input: PathBuf::from("test.bam"),
            output: None,
            verify: false,
            order,
            key_types: None,
            max_memory: MemoryLimit::Fixed(512 * 1024 * 1024),
            memory_reserve: MemoryReserve::Auto,
            memory_per_thread: true,
            tmp_dirs: Vec::new(),
            threads: 1,
            sort_threads: None,
            merge_threads: None,
            compression: CompressionOptions::default(),
            temp_compression: 1,
            temp_codec: fgumi_sort::SpillCodec::default(),
            max_temp_files: MaxTempFiles::Auto,
            write_index: false,
            read_streams: fgumi_sort::ReadStreams::Auto,
            sort_stats: false,
        }
    }

    #[rstest]
    #[case(SortOrderArg::TemplateCoordinate, Some(SamTag::CB))]
    #[case(SortOrderArg::Coordinate, None)]
    #[case(SortOrderArg::Queryname, None)]
    fn test_parse_cell_tag(#[case] order: SortOrderArg, #[case] expected: Option<SamTag>) {
        let sort = make_sort(order);
        assert_eq!(sort.parse_cell_tag().expect("parse_cell_tag should succeed"), expected);
    }

    #[test]
    fn test_parse_memory_auto() {
        assert_eq!(
            parse_memory("auto").expect("parse_memory should succeed for 'auto'"),
            MemoryLimit::Auto
        );
        assert_eq!(
            parse_memory("AUTO").expect("parse_memory should succeed for 'AUTO'"),
            MemoryLimit::Auto
        );
        assert_eq!(
            parse_memory("Auto").expect("parse_memory should succeed for 'Auto'"),
            MemoryLimit::Auto
        );
    }

    #[test]
    fn test_parse_memory_plain_numbers_as_mb() {
        // Plain numbers are interpreted as MB (via parse_memory_size)
        assert_eq!(
            parse_memory("768").expect("parse_memory should succeed for 768"),
            MemoryLimit::Fixed(768 * 1024 * 1024)
        );
        assert_eq!(
            parse_memory("1").expect("parse_memory should succeed for 1"),
            MemoryLimit::Fixed(1024 * 1024)
        );
    }

    #[test]
    fn test_parse_memory_human_readable() {
        // Suffixed values use ByteSize (decimal: G=1000^3, M=1000^2)
        assert_eq!(
            parse_memory("512MB").expect("parse_memory should succeed for 512MB"),
            MemoryLimit::Fixed(512 * 1000 * 1000)
        );
        assert_eq!(
            parse_memory("1G").expect("parse_memory should succeed for 1G"),
            MemoryLimit::Fixed(1_000_000_000)
        );
        assert_eq!(
            parse_memory("2GB").expect("parse_memory should succeed for 2GB"),
            MemoryLimit::Fixed(2_000_000_000)
        );
        // Binary suffixes (GiB, MiB)
        assert_eq!(
            parse_memory("1GiB").expect("parse_memory should succeed for 1GiB"),
            MemoryLimit::Fixed(1024 * 1024 * 1024)
        );
        assert_eq!(
            parse_memory("512MiB").expect("parse_memory should succeed for 512MiB"),
            MemoryLimit::Fixed(512 * 1024 * 1024)
        );
    }

    #[test]
    fn test_parse_memory_case_insensitive() {
        assert_eq!(
            parse_memory("512mb").expect("parse_memory should succeed for lowercase 512mb"),
            MemoryLimit::Fixed(512 * 1000 * 1000)
        );
        assert_eq!(
            parse_memory("1gb").expect("parse_memory should succeed for lowercase 1gb"),
            MemoryLimit::Fixed(1_000_000_000)
        );
    }

    #[test]
    fn test_parse_memory_decimal_with_suffix() {
        assert_eq!(
            parse_memory("1.5GB").expect("parse_memory should succeed for 1.5GB"),
            MemoryLimit::Fixed(1_500_000_000)
        );
    }

    #[test]
    fn test_parse_memory_invalid() {
        assert!(parse_memory("").is_err());
        assert!(parse_memory("abc").is_err());
        assert!(parse_memory("-1G").is_err());
    }

    #[test]
    fn test_resolve_memory_limit_fixed() {
        let fixed = MemoryLimit::Fixed(1024 * 1024 * 1024); // 1 GiB
        // Reserve is ignored for fixed limits
        let resolved =
            resolve_memory_budget(fixed, MemoryReserve::Auto, 4, true).expect("should succeed");
        // Fixed + memory_per_thread: total = 1 GiB * 4 = 4 GiB
        assert_eq!(resolved, 4 * 1024 * 1024 * 1024);
    }

    #[test]
    fn test_resolve_memory_limit_fixed_no_per_thread() {
        let fixed = MemoryLimit::Fixed(4 * 1024 * 1024 * 1024); // 4 GiB
        let resolved =
            resolve_memory_budget(fixed, MemoryReserve::Auto, 4, false).expect("should succeed");
        // Fixed + no per-thread: total = 4 GiB
        assert_eq!(resolved, 4 * 1024 * 1024 * 1024);
    }

    #[test]
    fn test_resolve_memory_limit_auto() {
        let total = detect_total_memory();

        let resolved = resolve_memory_budget(MemoryLimit::Auto, MemoryReserve::Auto, 4, true)
            .expect("should succeed");
        // Must be at least the per-thread minimum floor (256 MiB * 4 threads)
        let min_expected = MIN_MEMORY_PER_THREAD.saturating_mul(4).min(total);
        assert!(
            resolved >= min_expected,
            "auto resolved to {resolved} bytes, expected at least {min_expected}"
        );
        // And not more than effective total memory (cgroup-aware), unless total
        // memory is so small that the per-thread floor already exceeds it.
        if total >= MIN_MEMORY_PER_THREAD.saturating_mul(4) {
            assert!(resolved <= total);
        }
    }

    #[test]
    fn test_resolve_memory_limit_auto_no_per_thread() {
        let resolved = resolve_memory_budget(MemoryLimit::Auto, MemoryReserve::Auto, 8, false)
            .expect("should succeed");
        // Auto + no per-thread: should be total budget, not divided by threads
        // Must be at least the minimum floor (256 MB)
        assert!(resolved >= 256 * 1024 * 1024);
    }

    #[test]
    fn test_resolve_reserve_auto() {
        let gib = 1024 * 1024 * 1024;
        // 32 GiB system: min(10 GiB, 16 GiB) = 10 GiB
        assert_eq!(resolve_reserve(MemoryReserve::Auto, 32 * gib), 10 * gib);
        // 16 GiB system: min(10 GiB, 8 GiB) = 8 GiB
        assert_eq!(resolve_reserve(MemoryReserve::Auto, 16 * gib), 8 * gib);
        // 8 GiB system: min(10 GiB, 4 GiB) = 4 GiB
        assert_eq!(resolve_reserve(MemoryReserve::Auto, 8 * gib), 4 * gib);
        // 128 GiB system: min(10 GiB, 64 GiB) = 10 GiB
        assert_eq!(resolve_reserve(MemoryReserve::Auto, 128 * gib), 10 * gib);
    }

    #[test]
    fn test_resolve_reserve_fixed() {
        let gib = 1024 * 1024 * 1024;
        assert_eq!(resolve_reserve(MemoryReserve::Fixed(12 * gib), 64 * gib), 12 * gib);
    }

    #[test]
    fn test_parse_memory_reserve() {
        assert_eq!(parse_memory_reserve("auto").expect("should parse 'auto'"), MemoryReserve::Auto,);
        assert_eq!(
            parse_memory_reserve("10GiB").expect("should parse '10GiB'"),
            MemoryReserve::Fixed(10 * 1024 * 1024 * 1024),
        );
        assert_eq!(
            parse_memory_reserve("8G").expect("should parse '8G'"),
            MemoryReserve::Fixed(8_000_000_000),
        );
    }

    #[test]
    fn test_resolve_memory_limit_auto_with_fixed_reserve() {
        // With a larger fixed reserve, auto should return less memory.
        // Use modest reserve sizes (128 MiB vs 512 MiB) to stay well within
        // CI runner RAM and avoid the per-thread floor clamping both results.
        let large_reserve = resolve_memory_budget(
            MemoryLimit::Auto,
            MemoryReserve::Fixed(512 * 1024 * 1024),
            4,
            true,
        )
        .expect("should succeed");
        let small_reserve = resolve_memory_budget(
            MemoryLimit::Auto,
            MemoryReserve::Fixed(128 * 1024 * 1024),
            4,
            true,
        )
        .expect("should succeed");
        assert!(large_reserve < small_reserve);
    }

    #[test]
    fn test_sort_order_conversion() {
        assert_eq!(SortOrder::from(SortOrderArg::Coordinate), SortOrder::Coordinate);
        assert_eq!(
            SortOrder::from(SortOrderArg::Queryname),
            SortOrder::Queryname(QuerynameComparator::Lexicographic)
        );
        assert_eq!(
            SortOrder::from(SortOrderArg::QuerynameNatural),
            SortOrder::Queryname(QuerynameComparator::Natural)
        );
        assert_eq!(
            SortOrder::from(SortOrderArg::TemplateCoordinate),
            SortOrder::TemplateCoordinate
        );
    }

    // ========================================================================
    // SortOrderArg::parse tests
    // ========================================================================

    #[rstest]
    #[case("coordinate", Ok(SortOrderArg::Coordinate))]
    #[case("queryname", Ok(SortOrderArg::Queryname))]
    #[case("queryname::lexicographic", Ok(SortOrderArg::Queryname))]
    // SORT3-10: the spec-spelling alias is accepted so a user can pass back the
    // `@HD SS:queryname:lexicographical` value fgumi now emits.
    #[case("queryname::lexicographical", Ok(SortOrderArg::Queryname))]
    #[case("queryname::lex", Ok(SortOrderArg::Queryname))]
    #[case("queryname::natural", Ok(SortOrderArg::QuerynameNatural))]
    #[case("template-coordinate", Ok(SortOrderArg::TemplateCoordinate))]
    #[case("queryname::fast", Err("unknown queryname sub-sort 'fast'"))]
    #[case("random", Err("unknown sort order 'random'"))]
    #[case("queryname::", Err("unknown queryname sub-sort ''"))]
    fn test_parse_sort_order(#[case] input: &str, #[case] expected: Result<SortOrderArg, &str>) {
        match expected {
            Ok(order) => assert_eq!(
                SortOrderArg::parse(input).expect("parse should succeed for valid sort order"),
                order
            ),
            Err(msg) => {
                let err = SortOrderArg::parse(input)
                    .expect_err("parse should fail for invalid sort order");
                assert!(err.contains(msg), "expected error containing {msg:?}, got: {err}");
            }
        }
    }

    // ========================================================================
    // Header sub-sort tag tests
    // ========================================================================

    #[test]
    fn test_queryname_lex_header_has_subsort() {
        let order = SortOrder::from(SortOrderArg::Queryname);
        assert_eq!(order.header_so_tag(), "queryname");
        // SORT3-10: the emitted SS sub-sort value uses the spec spelling "lexicographical"
        // (with the `<sort-order>:` prefix retained).
        assert_eq!(order.header_ss_tag(), Some("queryname:lexicographical"));
    }

    #[test]
    fn test_queryname_natural_header_has_subsort() {
        let order = SortOrder::from(SortOrderArg::QuerynameNatural);
        assert_eq!(order.header_so_tag(), "queryname");
        assert_eq!(order.header_ss_tag(), Some("queryname:natural"));
    }

    #[test]
    fn test_coordinate_header_no_subsort() {
        let order = SortOrder::from(SortOrderArg::Coordinate);
        assert_eq!(order.header_so_tag(), "coordinate");
        assert_eq!(order.header_ss_tag(), None);
    }

    #[test]
    fn test_template_coordinate_header_subsort() {
        let order = SortOrder::from(SortOrderArg::TemplateCoordinate);
        assert_eq!(order.header_so_tag(), "unsorted");
        assert_eq!(order.header_ss_tag(), Some("unsorted:template-coordinate"));
    }

    #[test]
    fn test_verify_conflicts_with_output() {
        let sort = Sort {
            verify: true,
            output: Some(PathBuf::from("out.bam")),
            ..make_sort(SortOrderArg::Coordinate)
        };
        let err = sort.execute("test").unwrap_err();
        assert!(err.to_string().contains("--verify cannot be used with --output"));
    }

    #[test]
    fn test_verify_conflicts_with_write_index() {
        let sort = Sort { verify: true, write_index: true, ..make_sort(SortOrderArg::Coordinate) };
        let err = sort.execute("test").unwrap_err();
        assert!(err.to_string().contains("--write-index cannot be used with --verify"));
    }

    #[test]
    fn test_temp_compression_zero_with_zstd_rejected() {
        let sort = Sort {
            output: Some(PathBuf::from("out.bam")),
            temp_compression: 0,
            temp_codec: fgumi_sort::SpillCodec::Zstd,
            ..make_sort(SortOrderArg::Coordinate)
        };
        // Call execute_sort directly to bypass the input-file existence check
        // in execute(); the codec validation lives inside execute_sort.
        let err = sort.execute_sort("test").unwrap_err();
        let msg = err.to_string();
        assert!(
            msg.contains("--temp-compression 0 is only supported with --temp-codec bgzf"),
            "unexpected error: {msg}"
        );
    }

    #[test]
    fn test_verify_coordinate_fails_on_unsorted() -> Result<()> {
        use fgumi_sort::RawBamRecordReader;
        use fgumi_sort::extract_coordinate_key_inline;

        // Build BAM with records deliberately out of coordinate order
        let mut builder = crate::sam::builder::SamBuilder::new();
        let _ = builder.add_pair().name("a").contig(1).start1(100).build();
        let _ = builder.add_pair().name("b").contig(0).start1(200).build();

        let dir = tempfile::tempdir()?;
        let bam_path = dir.path().join("unsorted.bam");
        builder.write_bam(&bam_path)?;

        let file = std::fs::File::open(&bam_path)?;
        let (_, header) = fgumi_bam_io::create_bam_reader(&bam_path, 1)?;
        let mut reader = RawBamRecordReader::new(file)?;
        reader.skip_header()?;

        let nref = header.reference_sequences().len() as u32;
        let (total, violations, _) = verify_sort_order(
            reader,
            |bam| extract_coordinate_key_inline(bam, nref),
            |key, prev| key < prev,
        )?;

        assert!(total > 0);
        assert!(violations > 0, "unsorted file should fail coordinate verify");
        Ok(())
    }

    /// Verifies that template-coordinate sort groups reads by CB when pairs share
    /// the same outer genomic coordinates. CB is used as a hash-based tiebreaker,
    /// so reads with the same CB value must appear contiguously in the output.
    #[test]
    fn test_template_coordinate_sorts_by_cell_barcode() -> Result<()> {
        use crate::commands::command::Command;
        use crate::sam::builder::SamBuilder;
        use bstr::ByteSlice;

        let dir = tempfile::tempdir()?;
        let input = dir.path().join("input.bam");
        let output = dir.path().join("output.bam");

        let mut builder = SamBuilder::new();
        // Three pairs at the same position: two with CB=A and one with CB=B interleaved.
        // After sorting, the two CB=A pairs must be adjacent (not split by CB=B).
        let _ = builder.add_pair().name("pair_a1").contig(0).start1(100).attr("CB", "A").build();
        let _ = builder.add_pair().name("pair_b").contig(0).start1(100).attr("CB", "B").build();
        let _ = builder.add_pair().name("pair_a2").contig(0).start1(100).attr("CB", "A").build();
        builder.write_bam(&input)?;

        let mut sort = make_sort(SortOrderArg::TemplateCoordinate);
        sort.input = input;
        sort.output = Some(output.clone());
        sort.execute("test")?;

        let mut reader = noodles::bam::io::reader::Builder.build_from_path(&output)?;
        let header = reader.read_header()?;
        let records: Vec<_> = reader.record_bufs(&header).collect::<std::io::Result<Vec<_>>>()?;

        assert_eq!(records.len(), 6, "should have 6 records (3 pairs × 2 reads)");

        // Collect output record names and find positions of the two CB=A pairs.
        // They must be contiguous — i.e. not interleaved with the CB=B pair.
        let names: Vec<String> = records
            .iter()
            .map(|r| {
                r.name()
                    .map(|n| String::from_utf8_lossy(n.as_bytes()).into_owned())
                    .unwrap_or_default()
            })
            .collect();
        let a_positions: Vec<usize> = names
            .iter()
            .enumerate()
            .filter(|(_, n)| n.starts_with("pair_a"))
            .map(|(i, _)| i)
            .collect();
        assert_eq!(a_positions.len(), 4, "expected 4 reads for the two CB=A pairs");
        let min = a_positions[0];
        let max = *a_positions.last().unwrap();
        assert_eq!(
            max - min,
            3,
            "CB=A reads must be grouped together; got positions {a_positions:?}"
        );

        Ok(())
    }

    // ========================================================================
    // parse_key_types tests
    // ========================================================================

    #[rstest]
    #[case("full", KeyTypesSpec::Full)]
    #[case("none", KeyTypesSpec::None)]
    #[case("", KeyTypesSpec::None)]
    #[case("cb", KeyTypesSpec::Explicit { cb: true, tertiary: false })]
    #[case("library", KeyTypesSpec::Explicit { cb: false, tertiary: true })]
    #[case("mi", KeyTypesSpec::Explicit { cb: false, tertiary: true })]
    #[case("library,mi", KeyTypesSpec::Explicit { cb: false, tertiary: true })]
    #[case("cb,mi", KeyTypesSpec::Explicit { cb: true, tertiary: true })]
    #[case("cb library", KeyTypesSpec::Explicit { cb: true, tertiary: true })]
    #[case("FULL", KeyTypesSpec::Full)]
    #[case("None", KeyTypesSpec::None)]
    #[case("CB", KeyTypesSpec::Explicit { cb: true, tertiary: false })]
    #[case("Cb,MI", KeyTypesSpec::Explicit { cb: true, tertiary: true })]
    fn test_parse_key_types_ok(#[case] input: &str, #[case] expected: KeyTypesSpec) {
        assert_eq!(parse_key_types(input).expect("valid"), expected);
    }

    #[rstest]
    #[case("bogus")]
    #[case("cb,bogus")]
    fn test_parse_key_types_err(#[case] input: &str) {
        assert!(parse_key_types(input).is_err());
    }

    #[test]
    fn test_key_types_clap_default_is_none_option() {
        let sort = Sort::try_parse_from(["sort", "-i", "in.bam", "-o", "out.bam"]).expect("parse");
        assert!(sort.key_types.is_none());
    }
}
