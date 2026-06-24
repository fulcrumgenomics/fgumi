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
use std::sync::Arc;

use anyhow::{Result, bail};
use clap::Parser;
use fgumi_bam_io::create_raw_bam_reader;
use fgumi_cli_common::{
    Command, CompressionOptions, MemoryLimit, MemoryReserve, OperationTimer, parse_bool,
    parse_memory, parse_memory_reserve, resolve_memory_budget, validate_file_exists,
};
use fgumi_pipeline_core::{FinalizeHook, Pipeline, PipelineConfig};
use fgumi_sam::SamTag;
use fgumi_sort::{KeyTypesSpec, QuerynameComparator, SortOrder, SpillCodec, verify_sort_order};
use log::info;
use parking_lot::Mutex;

use crate::chains::{
    IndexBamFinalizeHook, SortFinalizeHook, SortStepCaptures, build_sort_step, log_sort_start,
};

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

    /// Worker threads for the merge/write phase (Phase 2).
    ///
    /// Defaults to `--threads`. In standalone `fgumi sort` (and runall's
    /// sole-stage sort) this caps the merge/write workers of the sort engine's
    /// pool. In a multi-stage runall chain it caps the number of concurrent
    /// spill-decompression workers (the merge step itself is serial). Output is
    /// byte-identical regardless of this value.
    #[arg(long = "merge-threads")]
    pub merge_threads: Option<usize>,

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
            order: SortOrderArg::TemplateCoordinate,
        }
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

impl Command for Sort {
    fn execute(&self, command_line: &str) -> Result<()> {
        if self.verify && self.output.is_some() {
            bail!("--verify cannot be used with --output");
        }
        if self.verify && self.write_index {
            bail!("--write-index cannot be used with --verify");
        }

        // Validate inputs. Exempt stdin paths (`-` / `/dev/stdin`): the
        // streaming sort path reads stdin once (the sort engine's reader
        // handles stdin directly), so a file-existence check would spuriously
        // reject it. `--verify` is the exception: it re-scans the input, which
        // a non-seekable stdin can't satisfy, so reject stdin there up front.
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
            // --verify is a standalone read-and-check path: no pipeline, just a
            // raw record-stream walk that compares each key to the previous.
            return self.execute_verify();
        }

        let output = self.output.as_ref().expect("output required for sort mode");

        // BAI is only defined for coordinate sort.
        if self.write_index && !matches!(self.order, SortOrderArg::Coordinate) {
            bail!("--write-index is only valid for coordinate sort");
        }

        // Up-front guard (S3-002): zstd has no uncompressed mode, so reject the
        // `--temp-compression 0` + `--temp-codec zstd` combination before the run
        // starts rather than failing late. (The engine also enforces this via
        // `ensure_spill_codec_compat`, but catching it here gives a CLI-level
        // message and avoids any setup work.)
        if self.options.temp_compression == 0 && matches!(self.options.temp_codec, SpillCodec::Zstd)
        {
            bail!(
                "--temp-compression 0 is only supported with --temp-codec bgzf; \
                 zstd does not have an uncompressed mode. Pass --temp-codec bgzf \
                 to keep uncompressed spill, or use a zstd level >= 1."
            );
        }

        // `--key-types` only applies to template-coordinate sort; warn if set for
        // any other order so the override is not silently ignored.
        if self.options.key_types.is_some()
            && !matches!(self.order, SortOrderArg::TemplateCoordinate)
        {
            info!("--key-types is ignored for --order {:?}", self.order);
        }

        // Copy order into SortOptions so the step factory can read it.
        let mut sort_opts = self.options.clone();
        sort_opts.order = self.order;

        self.execute_sort(sort_opts, output, command_line)
    }
}

impl Sort {
    /// Execute the non-verify sort path by building the typed-step pipeline
    /// directly (no chains monolith): a single `Exclusive` `SortBamFile` step
    /// registered as a source, run on one framework driver thread while the
    /// sort engine's own worker pool provides internal concurrency. Finalize
    /// actions (summary log, optional BAI index) run after `Pipeline::run`.
    fn execute_sort(
        &self,
        sort_opts: SortOptions,
        output: &std::path::Path,
        command_line: &str,
    ) -> Result<()> {
        let input_path = self.input.clone();
        let output_path = output.to_path_buf();
        let num_sorter_threads = self.threads.max(1);
        // Per-phase worker counts; each defaults to --threads. Memory budget
        // intentionally stays tied to --threads (these knobs are CPU-only). Read
        // from `sort_opts` (the resolved options used by the rest of the method)
        // so the phase counts cannot diverge from the logged / built config.
        let num_phase1_threads = sort_opts.sort_threads.unwrap_or(num_sorter_threads).max(1);
        let num_phase2_threads = sort_opts.merge_threads.unwrap_or(num_sorter_threads).max(1);

        let effective_memory = resolve_memory_budget(
            sort_opts.max_memory,
            sort_opts.memory_reserve,
            num_sorter_threads,
            sort_opts.memory_per_thread,
        )?;

        let timer = OperationTimer::new("Sorting BAM");
        log_sort_start(
            &sort_opts,
            &input_path,
            &output_path,
            num_sorter_threads,
            num_phase1_threads,
            num_phase2_threads,
            effective_memory,
        );

        let stats_slot = Arc::new(Mutex::new(None::<fgumi_sort::SortStats>));

        let sort_step = build_sort_step(SortStepCaptures {
            sort: sort_opts,
            input_path,
            output_path: output_path.clone(),
            num_sorter_threads,
            num_phase1_threads,
            num_phase2_threads,
            effective_memory,
            output_compression: self.compression.compression_level,
            command_line: command_line.to_string(),
            stats_slot: Arc::clone(&stats_slot),
            async_reader: self.async_reader,
        })?;

        // SortBamFile has `Input = ()` and `Outputs = ()`, so it registers via
        // `append_source` and the builder won't flag unwired branches. A single
        // Exclusive step needs only one framework driver thread.
        let builder = Pipeline::builder();
        let _ = builder.append_source(sort_step);
        let pipeline =
            builder.build().map_err(|e| anyhow::anyhow!("Pipeline build failed: {e}"))?;
        let config = PipelineConfig { threads: 1, ..Default::default() };

        // Always drain the summary finalize action, even on the error path, so
        // the summary logs. The pipeline-run error takes precedence.
        let run_result = pipeline.run(config).map_err(|e| anyhow::anyhow!("Pipeline::run: {e:?}"));

        // Drain the summary hook even on the error path (it is infallible —
        // always `Ok(())`), then let the pipeline-run error take precedence.
        let summary_result =
            Box::new(SortFinalizeHook { stats_slot, output_path: output_path.clone(), timer })
                .finalize();

        // Only write the BAM index after a successful sort: a failed run leaves
        // the output BAM incomplete, so emitting (or overwriting) `.bam.bai`
        // would publish a stale index for a partial file.
        run_result?;
        summary_result?;
        if self.write_index {
            Box::new(IndexBamFinalizeHook { output_path }).finalize()?;
        }
        Ok(())
    }

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
    fn test_temp_compression_zero_with_zstd_rejected_up_front() {
        // zstd has no uncompressed mode — the combination must be rejected
        // before the sort runs (S3-002). Use a valid input so the up-front
        // config guard, not input validation, is what fires.
        let dir = tempfile::tempdir().expect("tempdir");
        let input = dir.path().join("in.bam");
        let output = dir.path().join("out.bam");
        write_coordinate_bam(&input);

        let sort = Sort::try_parse_from([
            "sort",
            "-i",
            input.to_str().unwrap(),
            "-o",
            output.to_str().unwrap(),
            "--temp-compression",
            "0",
            "--temp-codec",
            "zstd",
        ])
        .expect("parse");
        let err = sort.execute("fgumi sort").expect_err("zstd + level 0 must be rejected");
        assert!(
            err.to_string().contains("only supported with --temp-codec bgzf"),
            "unexpected error: {err}"
        );
        assert!(!output.exists(), "no output should be written when the config is rejected");
    }

    #[test]
    fn test_temp_compression_zero_with_bgzf_runs() {
        // bgzf supports uncompressed (level-0) spill, so the combination is
        // accepted and the sort completes.
        let dir = tempfile::tempdir().expect("tempdir");
        let input = dir.path().join("in.bam");
        let output = dir.path().join("out.bam");
        write_coordinate_bam(&input);

        let sort = Sort::try_parse_from([
            "sort",
            "-i",
            input.to_str().unwrap(),
            "-o",
            output.to_str().unwrap(),
            "--order",
            "coordinate",
            // Bound the memory reservation: the default is 768 MiB per thread,
            // which a single-record test should not reserve (especially under
            // nextest's concurrent execution).
            "--max-memory",
            "16MiB",
            "--temp-compression",
            "0",
            "--temp-codec",
            "bgzf",
        ])
        .expect("parse");
        sort.execute("fgumi sort").expect("bgzf level-0 spill must sort successfully");
        assert!(output.exists(), "sorted output not written");
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

    /// Write a minimal valid coordinate-sorted BAM (one mapped record) to
    /// `path`, so a pre-placed output file looks indexable to the BAI hook.
    #[allow(clippy::cast_possible_truncation)]
    fn write_coordinate_bam(path: &std::path::Path) {
        use noodles::sam::Header;
        use noodles::sam::header::record::value::{Map, map::ReferenceSequence};
        use std::num::NonZeroUsize;

        let mut builder = Header::builder();
        builder = builder.add_reference_sequence(
            b"chr1",
            Map::<ReferenceSequence>::new(NonZeroUsize::new(1000).expect("non-zero")),
        );
        let header =
            fgumi_sort::create_output_header(fgumi_sort::SortOrder::Coordinate, &builder.build());

        // One 10-base mapped record on chr1:100 (BAM record body bytes).
        let name = b"r1";
        let name_with_null = name.len() + 1;
        let padding = (4 - (name_with_null % 4)) % 4;
        let mut record = Vec::with_capacity(64);
        record.extend_from_slice(&0_i32.to_le_bytes()); // ref_id
        record.extend_from_slice(&100_i32.to_le_bytes()); // pos
        record.push((name_with_null + padding) as u8); // l_read_name
        record.push(60_u8); // mapq
        record.extend_from_slice(&4681_u16.to_le_bytes()); // bin
        record.extend_from_slice(&1_u16.to_le_bytes()); // n_cigar_op
        record.extend_from_slice(&0_u16.to_le_bytes()); // flag
        record.extend_from_slice(&10_u32.to_le_bytes()); // l_seq
        record.extend_from_slice(&(-1_i32).to_le_bytes()); // next_ref_id
        record.extend_from_slice(&(-1_i32).to_le_bytes()); // next_pos
        record.extend_from_slice(&0_i32.to_le_bytes()); // tlen
        record.extend_from_slice(name);
        record.push(0);
        record.extend(std::iter::repeat_n(0_u8, padding));
        record.extend_from_slice(&(10_u32 << 4).to_le_bytes()); // 10M cigar
        record.extend_from_slice(&[0x11_u8; 5]); // packed seq
        record.extend_from_slice(&[30_u8; 10]); // qualities

        let mut writer = fgumi_bam_io::create_raw_bam_writer(path, &header, 1, 1)
            .expect("create_raw_bam_writer");
        writer.write_raw_record(&record).expect("write record");
        writer.finish().expect("finish");
    }

    /// A failed sort must not write (or overwrite) the `.bam.bai`: the index
    /// finalize hook only runs once the pipeline run succeeds.
    #[test]
    fn execute_sort_skips_index_when_sort_fails() {
        let dir = tempfile::tempdir().expect("tempdir");
        let missing_input = dir.path().join("does-not-exist.bam");
        let output = dir.path().join("out.bam");

        // Pre-place a valid coordinate BAM at the output path. The sort opens its
        // (missing) input before touching the output, so this file survives —
        // making the BAI's presence a clean signal of whether the index hook ran.
        write_coordinate_bam(&output);
        let bai = fgumi_bam_io::bai_sidecar_path(&output);
        assert!(!bai.exists(), "precondition: no index before the failed sort");

        // A Sort configured to write the index on success.
        let sort = Sort::try_parse_from([
            "sort",
            "--input",
            missing_input.to_str().expect("utf8"),
            "--output",
            output.to_str().expect("utf8"),
            "--order",
            "coordinate",
            "--write-index",
            "true",
        ])
        .expect("parse Sort");

        let mut sort_opts = sort.options.clone();
        sort_opts.order = sort.order;
        let result = sort.execute_sort(sort_opts, &output, "fgumi sort test");

        assert!(result.is_err(), "sort should fail when the input is missing");
        assert!(!bai.exists(), "index hook must not run after a failed sort");
    }
}
