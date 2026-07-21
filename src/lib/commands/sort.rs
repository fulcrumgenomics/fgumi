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
use fgumi_sort::{
    KeyTypesSpec, QuerynameComparator, RawExternalSorter, SortOrder, verify_sort_order,
};

use log::{debug, info};
use std::path::PathBuf;

use crate::commands::command::Command;
use crate::commands::common::{
    CompressionOptions, MemoryLimit, MemoryReserve, parse_bool, parse_memory, parse_memory_reserve,
    resolve_memory_budget,
};

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
    ///   `queryname::lexicographical` Alias; fgumi emits `queryname:lexicographical` in `@HD` SS
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

    /// Number of threads for parallel operations.
    ///
    /// Used for parallel sorting of in-memory chunks and parallel temp-chunk
    /// (BGZF or zstd) compression.
    #[arg(short = '@', short_alias = 't', long = "threads", default_value = "1")]
    pub threads: usize,

    /// Number of threads for the sort phase (accumulate, sort, spill).
    ///
    /// Defaults to `--threads`. Lower this to cede cores to an upstream
    /// producer while keeping the merge wide -- e.g. in
    /// `bwa mem -t 32 ... | fgumi sort -@ 8 --sort-threads 4`, ingest competes
    /// with the aligner but the merge does not.
    ///
    /// This only changes scheduling; the output is byte-identical.
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

    /// Write BAM index (.bai) alongside output.
    ///
    /// Only valid for coordinate sort. The index file will be written to
    /// `<output>.bai`. Output BGZF compression stays multi-threaded (scales
    /// with `--threads`); the BAI virtual offsets are recovered from each BGZF
    /// block as it finalizes, so indexing does not serialize compression.
    #[arg(long = "write-index", default_value = "false", num_args = 0..=1, default_missing_value = "true", action = clap::ArgAction::Set, value_parser = parse_bool)]
    pub write_index: bool,

    /// Enable async userspace prefetch on the input BAM.
    ///
    /// Spawns a dedicated I/O thread that reads raw bytes into a bounded
    /// queue ahead of the decompression step, so pool workers do not block
    /// on disk. Prototype flag; defaults to off.
    #[arg(long = "async-reader", default_value_t = false, hide = true)]
    pub async_reader: bool,
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
            return self.execute_verify();
        }

        self.execute_sort(command_line)
    }
}

impl Sort {
    /// Construct the sorter from this command's options.
    ///
    /// Extracted from `execute` so the option-to-builder wiring is testable on
    /// its own. That matters most for `--sort-threads` / `--merge-threads`:
    /// they only affect scheduling, so an end-to-end run cannot distinguish a
    /// correctly wired flag from one that was parsed and then dropped.
    fn build_sorter(&self, effective_memory: usize, command_line: &str) -> RawExternalSorter {
        let mut sorter = RawExternalSorter::new(self.order.into())
            .memory_limit(effective_memory)
            .threads(self.threads)
            .output_compression(self.compression.compression_level)
            .temp_compression(self.temp_compression)
            .spill_codec(self.temp_codec)
            .write_index(self.write_index)
            .async_reader(self.async_reader)
            .pg_info(crate::version::VERSION.to_string(), command_line.to_string());

        // Each per-phase override is optional and falls back to `--threads`.
        if let Some(n) = self.sort_threads {
            sorter = sorter.sort_threads(n);
        }
        if let Some(n) = self.merge_threads {
            sorter = sorter.merge_threads(n);
        }
        sorter
    }

    /// Parse the cell tag for template-coordinate sort/verify, returning `None`
    /// for other sort orders.
    fn parse_cell_tag(&self) -> Result<Option<SamTag>> {
        parse_cell_tag(self.order)
    }

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

        let timer = OperationTimer::new("Sorting BAM");

        // Resolve memory limit (auto-detect or fixed)
        let effective_memory = resolve_memory_budget(
            self.max_memory,
            self.memory_reserve,
            self.threads,
            self.memory_per_thread,
        )?;

        let cell_tag = self.parse_cell_tag()?;

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
                info!(
                    "Max memory: {} ({}/thread x {} threads)",
                    ByteSize(effective_memory as u64),
                    ByteSize(per_thread as u64),
                    self.threads
                );
            } else {
                info!("Max memory: {} (fixed)", ByteSize(effective_memory as u64));
            }
        }
        info!("Threads: {}", self.threads);
        info!("Temp compression level: {}", self.temp_compression);
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

        // Sort using raw-bytes sorter for optimal memory efficiency and speed
        let mut sorter = self.build_sorter(effective_memory, command_line);

        // For auto mode, cap initial buffer pre-allocation at 768 MiB/thread
        // (matching samtools default) to avoid huge upfront allocations.
        // The buffer will grow on demand up to memory_limit.
        if matches!(self.max_memory, MemoryLimit::Auto) {
            let init = 768_usize
                .checked_mul(1024 * 1024)
                .and_then(|b| b.checked_mul(self.threads))
                .ok_or_else(|| anyhow::anyhow!("initial auto buffer size overflowed"))?;
            sorter = sorter.initial_capacity(effective_memory.min(init));
        }

        if let Some(ct) = cell_tag {
            sorter = sorter.cell_tag(ct);
        }

        let key_types = self.key_types.unwrap_or_default(); // Auto
        if !matches!(self.order, SortOrderArg::TemplateCoordinate) && self.key_types.is_some() {
            info!("--key-types is ignored for --order {:?}", self.order);
        }
        sorter = sorter.key_types(key_types);

        if !resolved_tmp_dirs.is_empty() {
            sorter = sorter.temp_dirs(resolved_tmp_dirs);
        }

        let stats = sorter.sort(&self.input, output)?;
        let (total_records, output_records, chunks_written) =
            (stats.total_records, stats.output_records, stats.chunks_written);

        // Summary
        info!("=== Summary ===");
        info!("Records processed: {total_records}");
        info!("Records written: {output_records}");
        if chunks_written > 0 {
            info!("Temporary chunks: {chunks_written}");
        }
        info!("Output: {}", output.display());

        timer.log_completion(total_records);
        Ok(())
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

        debug!("Starting Sort Verification");
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
    // Memory-budget helpers moved to `commands::common`; import the `pub(crate)`
    // items these tests exercise that are not re-exported through `super::*`.
    use crate::commands::common::{MIN_MEMORY_PER_THREAD, detect_total_memory, resolve_reserve};
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

    /// The flags must actually reach the sorter. This cannot be checked through
    /// an end-to-end run: the per-phase counts only affect scheduling, so a flag
    /// that parsed and was then silently dropped produces byte-identical output
    /// to one that was wired correctly. Assert the resolved per-phase counts on
    /// the built sorter instead.
    #[rstest]
    #[case::defaults(None, None, 2, 2)]
    #[case::narrow_sort(Some(1), Some(4), 1, 4)]
    #[case::narrow_merge(Some(4), Some(1), 4, 1)]
    #[case::sort_only(Some(3), None, 3, 2)]
    #[case::merge_only(None, Some(3), 2, 3)]
    fn test_build_sorter_wires_per_phase_threads(
        #[case] sort_threads: Option<usize>,
        #[case] merge_threads: Option<usize>,
        #[case] expected_phase1: usize,
        #[case] expected_phase2: usize,
    ) {
        let mut sort = make_sort(SortOrderArg::Coordinate);
        sort.threads = 2;
        sort.sort_threads = sort_threads;
        sort.merge_threads = merge_threads;

        let sorter = sort.build_sorter(512 * 1024 * 1024, "fgumi sort (test)");
        assert_eq!(sorter.phase1_threads(), expected_phase1);
        assert_eq!(sorter.phase2_threads(), expected_phase2);
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
            write_index: false,
            async_reader: false,
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
