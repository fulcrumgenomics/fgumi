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

use std::path::PathBuf;

use anyhow::{Result, bail};
use clap::Parser;
use fgumi_bam_io::create_raw_bam_reader;
use fgumi_cli_common::{
    CompressionOptions, MemoryLimit, MemoryReserve, OperationTimer, parse_bool, parse_memory,
    parse_memory_reserve,
};
use fgumi_sam::SamTag;
use fgumi_sort::{KeyTypesSpec, QuerynameComparator, SortOrder, SpillCodec, verify_sort_order};
use log::info;

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
                if let Some(sub) = other.strip_prefix("queryname::") {
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
  - Default 768MiB per-thread memory limit (samtools-compatible); pass
    `--max-memory auto` to detect host memory (opt-in)

EXAMPLES:

  # Sort for fgumi group input
  fgumi sort -i aligned.bam -o sorted.bam --order template-coordinate

  # Sort by coordinate for IGV
  fgumi sort -i input.bam -o sorted.bam --order coordinate

  # Sort by queryname for zipper
  fgumi sort -i input.bam -o sorted.bam --order queryname

  # Multi-threaded sort (default 768MiB per thread)
  fgumi sort -i input.bam -o sorted.bam --order template-coordinate --threads 8

  # Override the per-thread memory limit
  fgumi sort -i input.bam -o sorted.bam -m 2GiB --threads 8

  # Opt in to auto-detected host memory (subtracts --memory-reserve)
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

    /// Print detailed pipeline occupancy statistics at completion (per-step /
    /// -worker busy-vs-idle, pool utilisation, and the off-pool "detached"
    /// coordination/I-O lines). Diagnostic; sort now runs through the typed-step
    /// pipeline, so this matches its siblings (correct/filter/clip/duplex). Also
    /// settable via `FGUMI_PIPELINE_STATS=1`.
    #[arg(long = "pipeline-stats", default_value = "false", num_args = 0..=1, default_missing_value = "true", action = clap::ArgAction::Set, value_parser = parse_bool, hide = true)]
    pub pipeline_stats: bool,

    /// Per-edge instrumentation level: off | summary | timeline | deep.
    ///
    /// `summary` adds per-edge throughput/occupancy/latency + a bottleneck
    /// verdict to the end-of-run report; `timeline` also writes a per-tick TSV
    /// (see `--pipeline-trace-out`); `deep` adds direct dwell/park-time latency.
    /// Diagnostic only — throughput is slightly depressed under tracing. Also
    /// settable via `FGUMI_PIPELINE_TRACE`.
    #[arg(long = "pipeline-trace", default_value = "off", value_parser = ["off", "summary", "timeline", "deep"], hide = true)]
    pub pipeline_trace: String,

    /// Path for the `--pipeline-trace timeline` per-tick TSV
    /// (default `pipeline-trace.tsv` in the working directory).
    #[arg(long = "pipeline-trace-out", hide = true)]
    pub pipeline_trace_out: Option<std::path::PathBuf>,
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
    /// Default is "768MiB" per thread (matching samtools' 768 MiB). Pass "auto"
    /// to detect host memory and subtract --memory-reserve, leaving room
    /// for the OS and co-running processes (e.g. an aligner). Explicit values
    /// like "512MiB", "1GiB", "4GiB" are per-thread when --memory-per-thread is
    /// enabled (default). Note bare "M"/"G" are decimal (1000ⁿ); "MiB"/"GiB" are
    /// binary (1024ⁿ).
    ///
    /// When the limit is reached, sorted chunks spill to temporary files.
    #[arg(short = 'm', long = "max-memory", default_value = "768MiB", value_parser = parse_memory)]
    pub max_memory: MemoryLimit,

    /// Memory to reserve for other processes when --max-memory=auto.
    ///
    /// "auto" (default) reserves min(10 GiB, 50% of host memory). Explicit
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

    /// Codec for temporary spill files: `zstd` (default) or `bgzf`.
    ///
    /// Spill files are internal to the sort and never read by other tools, so
    /// `zstd` is the default — it is faster than `bgzf` on every axis. Pass
    /// `bgzf` for the legacy on-disk format, or to use `--temp-compression 0`
    /// (uncompressed spill), which zstd does not support.
    #[arg(long = "temp-codec", default_value = "zstd")]
    pub temp_codec: SpillCodec,

    /// Expert override for which sort-key lanes to provision (template-coordinate
    /// order only).
    ///
    /// Omitted (default) auto-detects the narrowest key from the first record.
    /// Accepts `full`, `none` (or empty), or a comma/space list of `cb`,
    /// `library`, `mi` (`library`/`mi` both map to the shared tertiary lane).
    /// Ignored for non-template-coordinate orders.
    #[arg(long = "key-types", value_parser = parse_key_types)]
    pub key_types: Option<KeyTypesSpec>,

    /// Worker threads for the accumulation/sort/spill phase (Phase 1).
    ///
    /// Defaults to `--threads`. Lower this to cede CPU to an upstream producer
    /// (e.g. `bwa mem | fgumi sort`, or runall's align→sort chain); use
    /// `--merge-threads` for the merge/write phase. Output is byte-identical
    /// regardless of this value.
    #[arg(long = "sort-threads")]
    pub sort_threads: Option<usize>,

    /// Worker threads for the Phase-2 spill-decompression stage.
    ///
    /// Defaults to `--threads`. Every sort — standalone `fgumi sort` and every
    /// runall sort (sole-stage or fused) — runs the same streaming pipeline:
    /// this caps how many spill files decompress concurrently; the merge step
    /// itself is always serial. Output is byte-identical regardless of this
    /// value.
    #[arg(long = "merge-threads")]
    pub merge_threads: Option<usize>,

    /// Phase-2 spill decompression granularity (expert tuning).
    ///
    /// `false` (default) — block-level parallel: workers hold the per-file reader
    /// lock only for the READ (sequence-tagged), release it, then decompress
    /// OUTSIDE the lock, so multiple workers decompress one file's blocks
    /// concurrently, reassembled by a per-slot reorder buffer. Hardened by the
    /// loom model (`fgumi-sort` `tests/loom_merge_slots.rs`, driving the real
    /// `SortMergeSlot`) and the external-watchdog soak matrix.
    ///
    /// `true` — one worker owns a spill file's decompression at a time: it reads
    /// AND decompresses blocks inline under the reader lock, in order, into a
    /// plain FIFO. The older, single-worker-per-file path; pass `true` to fall
    /// back to it.
    #[arg(long = "file-granularity", default_value = "false", num_args = 0..=1, default_missing_value = "true", action = clap::ArgAction::Set, value_parser = parse_bool, hide = true)]
    pub file_granularity: bool,

    /// Raw blocks claimed per reader-lock acquisition during Phase-2
    /// decompression (expert tuning). Default `4` (the original
    /// `MAX_BATCH_PER_CALL`; a fleet bench will tune it). Must be >= 1.
    #[arg(long = "block-batch", default_value = "4", value_parser = clap::value_parser!(usize), hide = true)]
    pub block_batch: usize,

    /// Sort order (chain-builder slot).
    ///
    /// Carried here so the chain builder can read it from the bag without
    /// needing a separate out-of-band parameter. Populated by `Sort::execute`
    /// from `Sort::order` before constructing the step.
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
            temp_codec: SpillCodec::default(),
            key_types: None,
            sort_threads: None,
            merge_threads: None,
            file_granularity: false,
            block_batch: 4,
            order: SortOrderArg::TemplateCoordinate,
        }
    }
}

impl SortOptions {
    /// Validate the order-independent semantic constraints on sort options that
    /// clap's per-field value parser cannot express. This is the **single**
    /// validator for sort options — both the standalone `fgumi sort` path
    /// (`execute_sort_command`) and the `fgumi runall` path (after
    /// `MultiSortOptions::validate`) call it, so the two CLI surfaces share one
    /// set of rules and cannot drift (order-*dependent* rules like
    /// `--write-index` requires coordinate order stay with the standalone
    /// command, which owns the sink/order surface).
    ///
    /// Rules:
    /// - `--block-batch` must be `>= 1`. A value of `0` would read zero raw
    ///   blocks per reader-lock acquisition — the inline (file-granularity)
    ///   decompress path would declare a phantom EOF after reading nothing
    ///   (silent record loss), and the block-parallel path would livelock
    ///   (`queue_eof` never finalizes). `SortSpillDecompress::new` clamps to
    ///   `>= 1` as defense-in-depth, but the CLI rejects `0` up front rather
    ///   than letting a silent clamp be the *only* guard.
    /// - `--temp-compression 0` requires `--temp-codec bgzf`: zstd has no
    ///   uncompressed mode, so level 0 + zstd is rejected.
    ///
    /// # Errors
    ///
    /// Returns an error if any option is out of range or mutually inconsistent.
    pub fn validate(&self) -> anyhow::Result<()> {
        // Name both surfaces in messages: `fgumi sort` exposes the bare flags,
        // `fgumi runall` the prefixed `--sort::*` flags (this struct is shared).
        anyhow::ensure!(
            self.block_batch >= 1,
            "--block-batch (--sort::block-batch under runall) must be >= 1 (got {})",
            self.block_batch
        );
        anyhow::ensure!(
            !(self.temp_compression == 0 && matches!(self.temp_codec, SpillCodec::Zstd)),
            "--temp-compression 0 (--sort::temp-compression 0 under runall) is only supported \
             with --temp-codec bgzf; zstd has no uncompressed mode. Pass --temp-codec bgzf to \
             keep uncompressed spill, or use a zstd level >= 1."
        );
        Ok(())
    }
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

/// Environment variable name for the fallback temp-dir list, parsed as a
/// `PATH`-style list when no `-T/--tmp-dir` flags are passed.
pub const TMP_DIRS_ENV: &str = "FGUMI_TMP_DIRS";

/// Resolve the final list of temp directories for a sort run.
///
/// Precedence: CLI flags (if non-empty) > `FGUMI_TMP_DIRS` env var > empty.
/// Empty strings and whitespace-only entries are filtered out of the env-var
/// value so that `FGUMI_TMP_DIRS=:` or trailing separators don't produce bogus
/// paths.
#[must_use]
pub fn resolve_tmp_dirs(cli: &[PathBuf], env_value: Option<&str>) -> Vec<PathBuf> {
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
///
/// # Errors
///
/// Currently infallible, but returns `Result` for forward compatibility with
/// configurable cell tags.
pub fn parse_cell_tag(order: SortOrderArg) -> Result<Option<SamTag>> {
    if matches!(order, SortOrderArg::TemplateCoordinate) { Ok(Some(SamTag::CB)) } else { Ok(None) }
}

impl Sort {
    /// Parse the cell tag for template-coordinate sort/verify, returning `None`
    /// for other sort orders.
    fn parse_cell_tag(&self) -> Result<Option<SamTag>> {
        parse_cell_tag(self.order)
    }

    /// Execute verify mode: read records and check sort order.
    ///
    /// # Errors
    ///
    /// Returns an error if the input cannot be opened or read, or if any
    /// record is found out of order relative to the requested sort order.
    pub fn execute_verify(&self) -> Result<()> {
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
                let nref = u32::try_from(header.reference_sequences().len()).unwrap_or(u32::MAX);
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
                    // so both fgumi- and samtools-sorted files pass.
                    |key, prev| key.core_cmp(prev) == Ordering::Less,
                )?
            }
        };

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
    fn default_decompress_tuning_is_block_parallel_batch_of_four() {
        // block_batch default is 4 (the original `MAX_BATCH_PER_CALL`, pending a
        // fleet bench); file_granularity default is false (block-parallel is the
        // hardened production default, P5c). Both clap `default_value`s must
        // round-trip through `SortOptions::default`.
        let defaults = SortOptions::default();
        assert_eq!(defaults.block_batch, 4);
        assert!(
            !defaults.file_granularity,
            "block-parallel is the default (file_granularity=false)"
        );

        let sort = Sort::try_parse_from(["fgumi-sort", "-i", "in.bam", "-o", "out.bam"])
            .expect("parse with defaults");
        assert_eq!(sort.options.block_batch, 4, "CLI default block_batch must be 4");
        assert!(
            !sort.options.file_granularity,
            "CLI default must be block-parallel (file_granularity=false)",
        );
    }

    #[rstest]
    #[case(1)]
    #[case(2)]
    #[case(4)]
    fn validate_accepts_positive_block_batch(#[case] block_batch: usize) {
        let opts = SortOptions { block_batch, ..SortOptions::default() };
        assert!(opts.validate().is_ok(), "block_batch={block_batch} must be accepted");
    }

    #[test]
    fn validate_rejects_zero_block_batch() {
        let opts = SortOptions { block_batch: 0, ..SortOptions::default() };
        let err = opts.validate().expect_err("block_batch 0 must be rejected");
        let msg = err.to_string();
        assert!(msg.contains("--block-batch"), "error should name the standalone flag: {msg}");
        assert!(
            msg.contains("--sort::block-batch"),
            "error should also name the runall flag: {msg}"
        );
        assert!(msg.contains(">= 1"), "error should state the constraint: {msg}");
    }

    #[test]
    fn validate_rejects_uncompressed_zstd_spill() {
        // temp-compression 0 + zstd has no uncompressed mode. This rule is now in
        // the shared validator, so BOTH `fgumi sort` and `fgumi runall` reject it
        // (previously runall slipped it through — the guard-parity gap this closes).
        let opts = SortOptions {
            temp_compression: 0,
            temp_codec: SpillCodec::Zstd,
            ..SortOptions::default()
        };
        let err = opts.validate().expect_err("temp-compression 0 + zstd must be rejected");
        let msg = err.to_string();
        assert!(msg.contains("temp-compression 0"), "error should name the flag: {msg}");
        assert!(msg.contains("bgzf"), "error should point at the bgzf fix: {msg}");
    }

    #[test]
    fn validate_accepts_uncompressed_bgzf_and_compressed_zstd() {
        // The two valid temp-spill configs: uncompressed bgzf, and zstd >= 1.
        let bgzf0 = SortOptions {
            temp_compression: 0,
            temp_codec: SpillCodec::Bgzf,
            ..SortOptions::default()
        };
        assert!(bgzf0.validate().is_ok(), "temp-compression 0 + bgzf is valid");
        let zstd1 = SortOptions {
            temp_compression: 1,
            temp_codec: SpillCodec::Zstd,
            ..SortOptions::default()
        };
        assert!(zstd1.validate().is_ok(), "zstd level 1 is valid");
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

    #[rstest]
    #[case("coordinate", SortOrderArg::Coordinate)]
    #[case("queryname", SortOrderArg::Queryname)]
    #[case("queryname::lex", SortOrderArg::Queryname)]
    #[case("queryname::lexicographic", SortOrderArg::Queryname)]
    #[case("queryname::natural", SortOrderArg::QuerynameNatural)]
    #[case("template-coordinate", SortOrderArg::TemplateCoordinate)]
    fn test_sort_order_parse_valid(#[case] input: &str, #[case] expected: SortOrderArg) {
        assert_eq!(SortOrderArg::parse(input).unwrap(), expected);
    }

    #[test]
    fn test_sort_order_parse_invalid() {
        assert!(SortOrderArg::parse("bogus").is_err());
        assert!(SortOrderArg::parse("queryname::bogus").is_err());
    }

    #[test]
    fn test_sort_order_to_sort_order() {
        assert_eq!(SortOrder::from(SortOrderArg::Coordinate), SortOrder::Coordinate);
        assert_eq!(
            SortOrder::from(SortOrderArg::TemplateCoordinate),
            SortOrder::TemplateCoordinate
        );
    }

    #[test]
    fn test_parse_cell_tag() {
        assert_eq!(parse_cell_tag(SortOrderArg::TemplateCoordinate).unwrap(), Some(SamTag::CB));
        assert_eq!(parse_cell_tag(SortOrderArg::Coordinate).unwrap(), None);
        assert_eq!(parse_cell_tag(SortOrderArg::Queryname).unwrap(), None);
    }

    #[test]
    fn test_sort_options_default_matches_clap() {
        // The Sort command must parse with only required args, taking all
        // SortOptions defaults.
        let sort = Sort::try_parse_from(["sort", "-i", "in.bam", "-o", "out.bam"]).unwrap();
        let defaults = SortOptions::default();
        assert_eq!(sort.options.max_memory, defaults.max_memory);
        assert_eq!(sort.options.memory_reserve, defaults.memory_reserve);
        assert_eq!(sort.options.memory_per_thread, defaults.memory_per_thread);
        assert_eq!(sort.options.temp_compression, defaults.temp_compression);
        assert_eq!(sort.options.temp_codec, defaults.temp_codec);
        assert_eq!(sort.options.temp_codec, SpillCodec::Zstd, "default spill codec is zstd");
        assert_eq!(sort.options.key_types, defaults.key_types);
        assert_eq!(sort.options.key_types, None, "--key-types defaults to Auto (None)");
    }

    #[rstest]
    #[case::default_omitted(&["sort", "-i", "in.bam", "-o", "out.bam"], SpillCodec::Zstd)]
    #[case::explicit_zstd(
        &["sort", "-i", "in.bam", "-o", "out.bam", "--temp-codec", "zstd"],
        SpillCodec::Zstd
    )]
    #[case::explicit_bgzf(
        &["sort", "-i", "in.bam", "-o", "out.bam", "--temp-codec", "bgzf"],
        SpillCodec::Bgzf
    )]
    fn test_temp_codec_parses(#[case] args: &[&str], #[case] expected: SpillCodec) {
        let sort = Sort::try_parse_from(args).expect("parse");
        assert_eq!(sort.options.temp_codec, expected);
    }

    #[test]
    fn test_parse_key_types() {
        assert_eq!(parse_key_types("full").unwrap(), KeyTypesSpec::Full);
        assert_eq!(parse_key_types("none").unwrap(), KeyTypesSpec::None);
        assert_eq!(parse_key_types("").unwrap(), KeyTypesSpec::None);
        assert_eq!(
            parse_key_types("cb").unwrap(),
            KeyTypesSpec::Explicit { cb: true, tertiary: false }
        );
        // `library` and `mi` both map to the shared tertiary lane.
        assert_eq!(
            parse_key_types("library").unwrap(),
            KeyTypesSpec::Explicit { cb: false, tertiary: true }
        );
        assert_eq!(
            parse_key_types("mi").unwrap(),
            KeyTypesSpec::Explicit { cb: false, tertiary: true }
        );
        assert_eq!(
            parse_key_types("cb, mi").unwrap(),
            KeyTypesSpec::Explicit { cb: true, tertiary: true }
        );
        assert!(parse_key_types("bogus").is_err(), "unknown token must be rejected");
    }

    #[test]
    fn test_phase_thread_flags_parse() {
        use clap::Parser;
        // Standalone: flags now live on SortOptions, still spelled --sort-threads / --merge-threads.
        let cmd = Sort::parse_from([
            "sort",
            "-i",
            "in.bam",
            "-o",
            "out.bam",
            "--sort-threads",
            "2",
            "--merge-threads",
            "5",
        ]);
        assert_eq!(cmd.options.sort_threads, Some(2));
        assert_eq!(cmd.options.merge_threads, Some(5));

        // Defaults: unset => None (engine falls back to --threads).
        let cmd = Sort::parse_from(["sort", "-i", "in.bam", "-o", "out.bam"]);
        assert_eq!(cmd.options.sort_threads, None);
        assert_eq!(cmd.options.merge_threads, None);
    }

    #[test]
    fn test_runall_prefixed_phase_thread_flags_parse() {
        use clap::Parser;
        #[derive(clap::Parser)]
        struct Probe {
            #[command(flatten)]
            sort: MultiSortOptions,
        }
        let p = Probe::parse_from(["x", "--sort::sort-threads", "3", "--sort::merge-threads", "7"]);
        let opts = p.sort.validate().expect("validate");
        assert_eq!(opts.sort_threads, Some(3));
        assert_eq!(opts.merge_threads, Some(7));
    }
}
