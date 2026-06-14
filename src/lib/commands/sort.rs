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
use crate::validation::validate_file_exists;
use anyhow::{Result, bail};
use bytesize::ByteSize;
use clap::Parser;
use fgumi_bam_io::create_raw_bam_reader;
use fgumi_sort::{QuerynameComparator, SortOrder};

use log::info;
use std::path::PathBuf;

use crate::commands::command::Command;
use crate::commands::common::{CompressionOptions, detect_total_memory, parse_bool};
use crate::validation::parse_memory_size;

/// Sort order for BAM files.
///
/// Queryname sort supports sub-sort specification via `::` syntax:
/// - `queryname` — lexicographic ordering (default, fast)
/// - `queryname::lexicographic` — explicit lexicographic ordering
/// - `queryname::natural` — natural numeric ordering (samtools-compatible)
#[derive(Debug, Clone, Copy, PartialEq, Eq, Default)]
pub enum SortOrderArg {
    /// Coordinate sort (tid → pos → strand)
    Coordinate,
    /// Queryname sort with lexicographic ordering (default)
    Queryname,
    /// Queryname sort with natural numeric ordering
    QuerynameNatural,
    /// Template-coordinate sort (for UMI grouping)
    #[default]
    TemplateCoordinate,
}

impl SortOrderArg {
    /// Parse a sort order string, supporting `::` sub-sort syntax for queryname.
    ///
    /// Valid values:
    /// - `coordinate`
    /// - `queryname` (default: lexicographic)
    /// - `queryname::lexicographic`
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
            "queryname" | "queryname::lex" | "queryname::lexicographic" => Ok(Self::Queryname),
            "queryname::natural" => Ok(Self::QuerynameNatural),
            "template-coordinate" => Ok(Self::TemplateCoordinate),
            other => {
                if other.starts_with("queryname::") {
                    let sub =
                        other.strip_prefix("queryname::").expect("guarded by starts_with check");
                    Err(format!(
                        "unknown queryname sub-sort '{sub}', expected 'lex', 'lexicographic', or 'natural'"
                    ))
                } else {
                    Err(format!(
                        "unknown sort order '{other}', expected 'coordinate', 'queryname', \
                         'queryname::lex', 'queryname::lexicographic', 'queryname::natural', \
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

    /// Wrap the input in a userspace async prefetch reader: a background
    /// thread reads ahead so disk I/O overlaps sorting/compression. Helps
    /// on slow or networked storage. Defaults to off.
    #[arg(long = "async-reader", default_value_t = false, hide = true)]
    pub async_reader: bool,

    /// Verify the input file is correctly sorted (no output written).
    ///
    /// Reads records sequentially and checks that each record's sort key
    /// is >= the previous record's key. Exits 0 if sorted correctly,
    /// non-zero if any records are out of order.
    #[arg(long = "verify", default_value = "false", num_args = 0..=1, default_missing_value = "true", action = clap::ArgAction::Set, value_parser = parse_bool)]
    pub verify: bool,

    /// Sort order.
    ///
    /// Queryname sort supports sub-sort specifiers:
    ///   `queryname`                  Lexicographic byte ordering (default, fast)
    ///   `queryname::lexicographic`   Explicit lexicographic ordering (alias: `queryname::lex`)
    ///   `queryname::natural`         Natural numeric ordering (samtools-compatible)
    #[arg(long = "order", default_value = "template-coordinate", value_parser = SortOrderArg::parse)]
    pub order: SortOrderArg,

    /// Per-stage sort tuning knobs (max memory, tmp dirs, etc.).
    /// Flattened here so `fgumi sort` exposes them as unprefixed
    /// flags (`--max-memory`, `-T`, …) and `fgumi runall` exposes
    /// them as prefixed `--sort::*` via `MultiSortOptions`.
    #[command(flatten)]
    pub options: SortOptions,

    /// Number of threads for parallel operations.
    ///
    /// Used for parallel sorting of in-memory chunks and
    /// multi-threaded BGZF compression.
    #[arg(short = '@', short_alias = 't', long = "threads", default_value = "1")]
    pub threads: usize,

    /// Compression options for output BAM.
    #[command(flatten)]
    pub compression: CompressionOptions,

    /// Write BAM index (.bai) alongside output.
    ///
    /// Only valid for coordinate sort. The index file is written to
    /// `<output>.bam.bai` after sort completes, by re-reading the
    /// finished BAM (a post-pipeline pass). The BAM itself uses the
    /// same multi-threaded BGZF compression as a non-`--write-index`
    /// run, so output bytes are identical between the two; only the
    /// extra BAI sidecar distinguishes them.
    #[arg(long = "write-index", default_value = "false", num_args = 0..=1, default_missing_value = "true", action = clap::ArgAction::Set, value_parser = parse_bool)]
    pub write_index: bool,
}

/// Represents the memory limit configuration for sorting.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum MemoryLimit {
    /// Automatically detect system memory and compute an optimal limit.
    Auto,
    /// Use a fixed memory limit in bytes.
    Fixed(usize),
}

impl Default for MemoryLimit {
    /// Matches the clap `default_value = "768M"` on `Sort::max_memory`
    /// (also used by `SortOptions::max_memory`).
    fn default() -> Self {
        Self::Fixed(768 * 1024 * 1024)
    }
}

impl std::fmt::Display for MemoryLimit {
    /// Round-trips through `parse_memory`. `Fixed(N)` is rendered in
    /// the largest unit that divides cleanly (`G` → `M` → `K` → `B`)
    /// so `--help` shows e.g. `"768M"` rather than `"805306368B"`.
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        match self {
            Self::Auto => f.write_str("auto"),
            Self::Fixed(bytes) => format_binary_bytes(*bytes, f),
        }
    }
}

/// Represents how much memory to reserve for other processes (OS, aligners, etc.)
/// when `--max-memory=auto`.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum MemoryReserve {
    /// Automatic: `min(10 GiB, 50% of system memory)`.
    Auto,
    /// Reserve a fixed number of bytes.
    Fixed(usize),
}

impl Default for MemoryReserve {
    /// Matches the clap `default_value = "auto"` on `Sort::memory_reserve`.
    fn default() -> Self {
        Self::Auto
    }
}

impl std::fmt::Display for MemoryReserve {
    /// Round-trips through `parse_memory_reserve`. See [`MemoryLimit`]'s
    /// `Display` for the formatting rule.
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        match self {
            Self::Auto => f.write_str("auto"),
            Self::Fixed(bytes) => format_binary_bytes(*bytes, f),
        }
    }
}

/// Format a byte count in the largest binary unit that divides
/// cleanly (`G` → `M` → `K` → `B`). Used by `MemoryLimit::Display`
/// and `MemoryReserve::Display` so that `parse_memory` /
/// `parse_memory_reserve` round-trip the result.
fn format_binary_bytes(bytes: usize, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
    const K: usize = 1024;
    const M: usize = K * 1024;
    const G: usize = M * 1024;
    if bytes >= G && bytes.is_multiple_of(G) {
        write!(f, "{}G", bytes / G)
    } else if bytes >= M && bytes.is_multiple_of(M) {
        write!(f, "{}M", bytes / M)
    } else if bytes >= K && bytes.is_multiple_of(K) {
        write!(f, "{}K", bytes / K)
    } else {
        write!(f, "{bytes}B")
    }
}

/// Per-stage sort tuning options, flattened into both the standalone
/// `Sort` command (as bare `--max-memory`, `-T`, … flags) and the
/// `RunAll` command (as prefixed `--sort::max-memory`, `--sort::tmp-dir`
/// flags via the `MultiSortOptions` companion struct generated by
/// `#[multi_options]`).
///
/// `Default` must match each field's clap `default_value` exactly —
/// `MultiSortOptions`' generated `default_value_t` reads from
/// `SortOptions::default().field`, and clap requires the default
/// shown in `--help` to round-trip through the value-parser.
#[fgumi_cli_macros::multi_options("sort", "Sort Options")]
#[derive(clap::Args, Debug, Clone)]
pub struct SortOptions {
    /// Maximum memory for in-memory sorting.
    ///
    /// Default is "768M" per thread (matching samtools behavior). Pass "auto"
    /// to detect system memory and subtract --memory-reserve, leaving room
    /// for the OS and co-running processes (e.g. an aligner). Explicit values
    /// like "512M", "1G", "4GiB" are per-thread when --memory-per-thread is
    /// enabled (default).
    ///
    /// When the limit is reached, sorted chunks spill to temporary files.
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
    /// When enabled (default), --max-memory specifies memory per thread.
    /// Total memory = `max_memory` × threads. Disable for fixed total memory.
    #[arg(long = "memory-per-thread", default_value = "true", num_args = 0..=1, default_missing_value = "true", action = clap::ArgAction::Set, value_parser = parse_bool)]
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

    /// Compression level for temporary chunk files (0-9).
    ///
    /// Level 0 disables compression (fastest, uses most disk space).
    /// Level 1 (default) provides fast compression with reasonable space savings.
    /// Higher levels (up to 9) provide better compression but are slower.
    #[arg(long = "temp-compression", default_value = "1", value_parser = clap::value_parser!(u32).range(0..=9))]
    pub temp_compression: u32,

    /// Sort order (chain-builder slot).
    ///
    /// Carried here so the chain builder can read it from the bag without
    /// needing a separate out-of-band parameter. Populated by `Sort::execute`
    /// from `Sort::order` before constructing the `ChainSpec`.
    #[arg(skip)]
    pub order: SortOrderArg,
}

impl Default for SortOptions {
    fn default() -> Self {
        Self {
            max_memory: MemoryLimit::default(),
            memory_reserve: MemoryReserve::default(),
            memory_per_thread: true,
            tmp_dirs: Vec::new(),
            temp_compression: 1,
            order: SortOrderArg::TemplateCoordinate,
        }
    }
}

/// Parse a memory size string into `usize` bytes, suitable for use in clap
/// value parsers.
///
/// Delegates to [`parse_memory_size`] for numeric parsing. Plain numbers are
/// interpreted as MiB (e.g. "768" = 768 MiB). Supports human-readable formats
/// like "2GB", "1GiB", "512MiB". See [`parse_memory_size`] for full details.
fn parse_memory_bytes(s: &str, label: &str) -> Result<usize, String> {
    let bytes = parse_memory_size(s).map_err(|e| e.to_string())?;
    usize::try_from(bytes).map_err(|_| format!("{label} too large: {bytes}"))
}

/// Parse memory size string (e.g., "512M", "1G", "2G", "auto").
pub(crate) fn parse_memory(s: &str) -> Result<MemoryLimit, String> {
    let s = s.trim();

    if s.eq_ignore_ascii_case("auto") {
        return Ok(MemoryLimit::Auto);
    }

    Ok(MemoryLimit::Fixed(parse_memory_bytes(s, "Memory size")?))
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

/// Parse memory reserve string (e.g., "10G", "auto").
pub(crate) fn parse_memory_reserve(s: &str) -> Result<MemoryReserve, String> {
    let s = s.trim();

    if s.eq_ignore_ascii_case("auto") {
        return Ok(MemoryReserve::Auto);
    }

    let bytes = parse_memory_bytes(s, "Memory reserve")?;
    Ok(MemoryReserve::Fixed(bytes))
}

/// Parse the cell tag for template-coordinate sort/verify, returning `None`
/// for other sort orders.
pub(crate) fn parse_cell_tag(order: SortOrderArg) -> Result<Option<SamTag>> {
    if matches!(order, SortOrderArg::TemplateCoordinate) { Ok(Some(SamTag::CB)) } else { Ok(None) }
}

use fgumi_sort::verify_sort_order;

impl Command for Sort {
    fn execute(&self, command_line: &str) -> Result<()> {
        if self.verify && self.output.is_some() {
            bail!("--verify cannot be used with --output");
        }
        if self.verify && self.write_index {
            bail!("--write-index cannot be used with --verify");
        }

        // Validate inputs. Exempt stdin paths (`-` / `/dev/stdin`): the sort
        // engine reads stdin in a single pass (a `TeeReader` reads the header,
        // then a `ChainedReader` replays it and streams the rest), so a
        // file-existence check would spuriously reject it — matching the stdin
        // exemption in group/dedup/clip. `--verify` is the exception: it
        // re-scans the input (header probe + a fresh `File::open` record pass),
        // which a non-seekable stdin can't satisfy, so reject stdin there up
        // front with a clear message.
        if fgumi_bam_io::is_stdin_path(&self.input) {
            if self.verify {
                bail!(
                    "fgumi sort --verify cannot read from stdin (it re-scans the input); \
                     provide a file path instead"
                );
            }
        } else {
            validate_file_exists(&self.input, "Input BAM")?;
        }

        // Either --output or --verify must be specified
        if !self.verify && self.output.is_none() {
            bail!("Either --output or --verify must be specified");
        }

        if self.verify {
            // --verify is a standalone read-and-check path: no
            // `ChainSpec`, no `Pipeline`, just a raw record-stream walk
            // that compares each key to the previous (see
            // `execute_verify`). It is intentionally NOT modeled as a
            // `SinkSpec` variant because there is no BAM/BAI output —
            // the path emits only log lines + a nonzero exit on
            // violations.
            return self.execute_verify();
        }

        // Route the non-verify sort path through chains::build_for.
        let output = self.output.as_ref().expect("output required for sort mode");

        // BAI is only defined for coordinate sort. The cross-stage validator
        // (Rule 3 in `chains::validate`) catches the spec-level constraint
        // "BamWithIndex requires Stage::Sort terminal"; this CLI check
        // enforces the file-format-physics constraint "BAI indexes
        // coordinate-sorted BGZF offsets specifically, not template-coordinate
        // or queryname". Together they cover both axes (chain shape and sort
        // order) of the indexability invariant.
        if self.write_index && !matches!(self.order, SortOrderArg::Coordinate) {
            bail!("--write-index is only valid for coordinate sort");
        }

        // Copy order into SortOptions so the chain builder can read it from
        // the bag. `--write-index` no longer rides on SortOptions: it routes
        // through the sink variant below so the chain-builder's
        // `IndexBamFinalizeHook` (registered by `add_sort`) owns BAI
        // generation as a post-pipeline step, decoupled from the BGZF
        // compression path.
        let mut sort_opts = self.options.clone();
        sort_opts.order = self.order;

        let sink = if self.write_index {
            crate::pipeline::chains::SinkSpec::BamWithIndex(output.clone())
        } else {
            crate::pipeline::chains::SinkSpec::Bam(output.clone())
        };

        let spec = crate::pipeline::chains::ChainSpec {
            stages: vec![crate::pipeline::chains::Stage::Sort],
            source: crate::pipeline::chains::SourceSpec::Bam(self.input.clone()),
            sink,
            stage_opts: crate::pipeline::chains::StageOptionsBag {
                sort: Some(sort_opts),
                ..Default::default()
            },
            threading: crate::commands::common::ThreadingOptions::new(self.threads),
            compression: self.compression.clone(),
            scheduler: crate::commands::common::SchedulerOptions::default(),
            queue_memory: crate::commands::common::QueueMemoryOptions::default(),
            async_reader: self.async_reader,
            command_line: command_line.to_string(),
        };
        crate::pipeline::chains::build_for(spec)?.run()
    }
}

/// The minimum per-thread memory budget (256 MiB).
const MIN_MEMORY_PER_THREAD: usize = 256 * 1024 * 1024;

/// Default auto-reserve cap: 10 GiB.
const AUTO_RESERVE_CAP: usize = 10 * 1024 * 1024 * 1024;

/// Resolve a [`MemoryReserve`] to a concrete byte count given total system memory.
fn resolve_reserve(reserve: MemoryReserve, total_memory: usize) -> usize {
    match reserve {
        MemoryReserve::Fixed(bytes) => bytes,
        MemoryReserve::Auto => {
            // min(10 GiB, 50% of system memory)
            AUTO_RESERVE_CAP.min(total_memory / 2)
        }
    }
}

/// Resolves a [`MemoryLimit`] to a concrete byte count.
///
/// For [`MemoryLimit::Auto`]: detects total system memory (cgroup-aware via
/// [`detect_total_memory`]), subtracts the reserve (via [`MemoryReserve`]),
/// and divides by thread count (when `memory_per_thread` is true). The result
/// is clamped to a minimum of 256 MiB per thread.
///
/// For [`MemoryLimit::Fixed`]: applies the same `memory_per_thread`
/// multiplication as before. The reserve is ignored.
///
/// # Note
///
/// Calls [`detect_total_memory`] exactly once regardless of `limit` variant.
/// That function invokes `sysinfo` (creates a `System` instance and refreshes
/// memory counters), so repeated calls add unnecessary overhead.
pub(crate) fn resolve_memory_limit(
    limit: MemoryLimit,
    reserve: MemoryReserve,
    threads: usize,
    memory_per_thread: bool,
) -> Result<usize> {
    if threads == 0 {
        bail!("--threads must be at least 1");
    }

    // Call once — detect_total_memory() invokes sysinfo, which is not free.
    let total = detect_total_memory();

    let total_budget = match limit {
        MemoryLimit::Fixed(bytes) => {
            if memory_per_thread {
                bytes
                    .checked_mul(threads)
                    .ok_or_else(|| anyhow::anyhow!("--max-memory × --threads overflowed"))?
            } else {
                bytes
            }
        }
        MemoryLimit::Auto => {
            let margin = resolve_reserve(reserve, total);
            let available = total.saturating_sub(margin);

            if memory_per_thread {
                let per_thread = (available / threads).max(MIN_MEMORY_PER_THREAD);
                let budget = per_thread
                    .checked_mul(threads)
                    .ok_or_else(|| anyhow::anyhow!("auto memory budget overflowed"))?;
                if budget > available {
                    log::warn!(
                        "Auto memory: total budget {} exceeds available {} \
                         ({}/thread x {} threads, reserve {}); may spill to disk earlier than expected",
                        ByteSize(budget as u64),
                        ByteSize(available as u64),
                        ByteSize(per_thread as u64),
                        threads,
                        ByteSize(margin as u64),
                    );
                }
                info!(
                    "Auto memory: using {} of {} ({}/thread x {} threads, reserve {})",
                    ByteSize(budget as u64),
                    ByteSize(total as u64),
                    ByteSize(per_thread as u64),
                    threads,
                    ByteSize(margin as u64),
                );
                budget
            } else {
                let budget = available.max(MIN_MEMORY_PER_THREAD);
                info!(
                    "Auto memory: using {} of {} (fixed total, reserve {})",
                    ByteSize(budget as u64),
                    ByteSize(total as u64),
                    ByteSize(margin as u64),
                );
                budget
            }
        }
    };

    // Post-resolution sanity check: warn if budget exceeds the effective memory limit.
    if total_budget > total {
        log::warn!(
            "Memory budget {} exceeds total system memory {}; spill-to-disk is likely",
            ByteSize(total_budget as u64),
            ByteSize(total as u64),
        );
    }

    Ok(total_budget)
}

impl Sort {
    /// Parse the cell tag for template-coordinate sort/verify, returning `None`
    /// for other sort orders.
    fn parse_cell_tag(&self) -> Result<Option<SamTag>> {
        parse_cell_tag(self.order)
    }

    /// Execute verify mode: read records and check sort order.
    fn execute_verify(&self) -> Result<()> {
        use fgumi_sort::RawBamRecordReader;
        use fgumi_sort::{
            LibraryLookup, RawQuerynameKey, RawQuerynameLexKey, RawSortKey, SortContext, cb_hasher,
            extract_coordinate_key_inline, extract_template_key_inline,
        };
        use std::cmp::Ordering;
        use std::fs::File;

        let cell_tag = self.parse_cell_tag()?;

        let timer = OperationTimer::new("Verifying BAM sort order");

        info!("Starting Sort Verification");
        info!("Input: {}", self.input.display());
        info!("Expected order: {:?}", self.order);
        if let Some(ct) = cell_tag {
            let ct_bytes = *ct;
            info!("Cell tag: {}{}", ct_bytes[0] as char, ct_bytes[1] as char);
        }

        // Get header via the raw-byte reader, then re-open for raw record iteration.
        let (_, header) = create_raw_bam_reader(&self.input, 1)?;

        let file = File::open(&self.input)?;
        let mut raw_reader = RawBamRecordReader::new(file)?;
        raw_reader.skip_header()?;

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
    use clap::Parser;
    use rstest::rstest;

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
        assert_eq!(sort.options.tmp_dirs, expected);
    }

    /// Helper to construct a `Sort` struct with a given order.
    fn make_sort(order: SortOrderArg) -> Sort {
        Sort {
            input: PathBuf::from("test.bam"),
            output: None,
            async_reader: false,
            verify: false,
            order,
            options: SortOptions {
                max_memory: MemoryLimit::Fixed(512 * 1024 * 1024),
                memory_reserve: MemoryReserve::Auto,
                memory_per_thread: true,
                tmp_dirs: Vec::new(),
                temp_compression: 1,
                order,
            },
            threads: 1,
            compression: CompressionOptions::default(),
            write_index: false,
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
            resolve_memory_limit(fixed, MemoryReserve::Auto, 4, true).expect("should succeed");
        // Fixed + memory_per_thread: total = 1 GiB * 4 = 4 GiB
        assert_eq!(resolved, 4 * 1024 * 1024 * 1024);
    }

    #[test]
    fn test_resolve_memory_limit_fixed_no_per_thread() {
        let fixed = MemoryLimit::Fixed(4 * 1024 * 1024 * 1024); // 4 GiB
        let resolved =
            resolve_memory_limit(fixed, MemoryReserve::Auto, 4, false).expect("should succeed");
        // Fixed + no per-thread: total = 4 GiB
        assert_eq!(resolved, 4 * 1024 * 1024 * 1024);
    }

    #[test]
    fn test_resolve_memory_limit_auto() {
        let total = detect_total_memory();

        let resolved = resolve_memory_limit(MemoryLimit::Auto, MemoryReserve::Auto, 4, true)
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
        let resolved = resolve_memory_limit(MemoryLimit::Auto, MemoryReserve::Auto, 8, false)
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
        let large_reserve = resolve_memory_limit(
            MemoryLimit::Auto,
            MemoryReserve::Fixed(512 * 1024 * 1024),
            4,
            true,
        )
        .expect("should succeed");
        let small_reserve = resolve_memory_limit(
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
        assert_eq!(order.header_ss_tag(), Some("lexicographic"));
    }

    #[test]
    fn test_queryname_natural_header_has_subsort() {
        let order = SortOrder::from(SortOrderArg::QuerynameNatural);
        assert_eq!(order.header_so_tag(), "queryname");
        assert_eq!(order.header_ss_tag(), Some("natural"));
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
        assert_eq!(order.header_ss_tag(), Some("template-coordinate"));
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
}
