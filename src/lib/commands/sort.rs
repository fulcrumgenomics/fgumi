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

use crate::bam_io::create_bam_reader;
use crate::logging::OperationTimer;
use crate::sam::SamTag;
use crate::sort::{QuerynameComparator, RawExternalSorter, SortOrder};
use crate::validation::validate_file_exists;
use anyhow::{Result, bail};
use bytesize::ByteSize;
use clap::Parser;

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
  - Auto memory detection with configurable reserve for co-running processes

EXAMPLES:

  # Sort for fgumi group input
  fgumi sort -i aligned.bam -o sorted.bam --order template-coordinate

  # Sort by coordinate for IGV
  fgumi sort -i input.bam -o sorted.bam --order coordinate

  # Sort by queryname for zipper
  fgumi sort -i input.bam -o sorted.bam --order queryname

  # Multi-threaded sort (auto-detects available memory)
  fgumi sort -i input.bam -o sorted.bam --order template-coordinate --threads 8

  # Override auto memory with an explicit per-thread limit
  fgumi sort -i input.bam -o sorted.bam -m 2GiB --threads 8

  # Reserve extra memory for bwa mem running in a pipeline
  fgumi sort -i input.bam -o sorted.bam --memory-reserve 12GiB --threads 4

  # Verify a BAM file is correctly sorted
  fgumi sort -i sorted.bam --verify --order template-coordinate
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
    ///   `queryname::natural`         Natural numeric ordering (samtools-compatible)
    #[arg(long = "order", default_value = "template-coordinate", value_parser = SortOrderArg::parse)]
    pub order: SortOrderArg,

    /// Maximum memory for in-memory sorting.
    ///
    /// "auto" (default) detects system memory and subtracts --memory-reserve
    /// to leave room for the OS and co-running processes (e.g. an aligner).
    /// Explicit values like "512M", "1G", "4GiB" are per-thread when
    /// --memory-per-thread is enabled (default).
    ///
    /// When the limit is reached, sorted chunks spill to temporary files.
    #[arg(short = 'm', long = "max-memory", default_value = "auto", value_parser = parse_memory)]
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

    /// Temporary directory for intermediate files.
    ///
    /// If not specified, uses the system default temp directory.
    /// For best performance, use a fast SSD.
    #[arg(short = 'T', long = "tmp-dir")]
    pub tmp_dir: Option<PathBuf>,

    /// Number of threads for parallel operations.
    ///
    /// Used for parallel sorting of in-memory chunks and
    /// multi-threaded BGZF compression.
    #[arg(short = '@', short_alias = 't', long = "threads", default_value = "1")]
    pub threads: usize,

    /// Compression options for output BAM.
    #[command(flatten)]
    pub compression: CompressionOptions,

    /// Compression level for temporary chunk files (0-9).
    ///
    /// Level 0 disables compression (fastest, uses most disk space).
    /// Level 1 (default) provides fast compression with reasonable space savings.
    /// Higher levels (up to 9) provide better compression but are slower.
    #[arg(long = "temp-compression", default_value = "1", value_parser = clap::value_parser!(u32).range(0..=9))]
    pub temp_compression: u32,

    /// Write BAM index (.bai) alongside output.
    ///
    /// Only valid for coordinate sort. The index file will be written to
    /// `<output>.bai`. Uses single-threaded compression for accurate virtual
    /// position tracking.
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

/// Represents how much memory to reserve for other processes (OS, aligners, etc.)
/// when `--max-memory=auto`.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum MemoryReserve {
    /// Automatic: `min(10 GiB, 50% of system memory)`.
    Auto,
    /// Reserve a fixed number of bytes.
    Fixed(usize),
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
pub(crate) fn parse_cell_tag(order: SortOrderArg) -> Result<Option<[u8; 2]>> {
    if matches!(order, SortOrderArg::TemplateCoordinate) { Ok(Some(*SamTag::CB)) } else { Ok(None) }
}

/// Summary of sort-order verification: `(total_records, violations, first_violation)`.
type VerifySummary = (u64, u64, Option<(u64, String)>);

/// Verify that records from a raw BAM reader are in sorted order.
///
/// Iterates all records, extracting a sort key from each and checking that
/// consecutive keys satisfy the ordering invariant (no violations).
fn verify_sort_order<K>(
    raw_reader: crate::sort::raw_bam_reader::RawBamRecordReader<std::fs::File>,
    extract_key: impl Fn(&[u8]) -> K,
    is_violation: impl Fn(&K, &K) -> bool,
) -> Result<VerifySummary> {
    let mut total_records: u64 = 0;
    let mut violations: u64 = 0;
    let mut first_violation: Option<(u64, String)> = None;
    let mut prev_key: Option<K> = None;

    for result in raw_reader {
        let record_bytes = result?;
        total_records += 1;
        let bam = record_bytes.as_slice();

        let key = extract_key(bam);

        if let Some(ref prev) = prev_key {
            if is_violation(&key, prev) {
                violations += 1;
                if first_violation.is_none() {
                    let name = String::from_utf8_lossy(fgumi_raw_bam::read_name(bam)).to_string();
                    first_violation = Some((total_records, name));
                }
            }
        }
        prev_key = Some(key);
    }

    Ok((total_records, violations, first_violation))
}

impl Command for Sort {
    fn execute(&self, command_line: &str) -> Result<()> {
        if self.verify && self.output.is_some() {
            bail!("--verify cannot be used with --output");
        }
        if self.verify && self.write_index {
            bail!("--write-index cannot be used with --verify");
        }

        // Validate inputs
        validate_file_exists(&self.input, "Input BAM")?;

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
fn resolve_memory_limit(
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
    fn parse_cell_tag(&self) -> Result<Option<[u8; 2]>> {
        parse_cell_tag(self.order)
    }

    /// Execute sort mode: read, sort, and write output.
    fn execute_sort(&self, command_line: &str) -> Result<()> {
        let output = self.output.as_ref().expect("output required for sort mode");

        // --write-index only valid for coordinate sort
        if self.write_index && !matches!(self.order, SortOrderArg::Coordinate) {
            bail!("--write-index is only valid for coordinate sort");
        }

        let timer = OperationTimer::new("Sorting BAM");

        // Resolve memory limit (auto-detect or fixed)
        let effective_memory = resolve_memory_limit(
            self.max_memory,
            self.memory_reserve,
            self.threads,
            self.memory_per_thread,
        )?;

        let cell_tag = self.parse_cell_tag()?;

        info!("Starting Sort");
        info!("Input: {}", self.input.display());
        info!("Output: {}", output.display());
        info!("Sort order: {:?}", self.order);
        if let Some(ct) = cell_tag {
            info!("Cell tag: {}{}", ct[0] as char, ct[1] as char);
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
        if let Some(ref tmp) = self.tmp_dir {
            info!("Temp directory: {}", tmp.display());
        }

        // Sort using raw-bytes sorter for optimal memory efficiency and speed
        let mut sorter = RawExternalSorter::new(self.order.into())
            .memory_limit(effective_memory)
            .threads(self.threads)
            .output_compression(self.compression.compression_level)
            .temp_compression(self.temp_compression)
            .write_index(self.write_index)
            .pg_info(crate::version::VERSION.to_string(), command_line.to_string());

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

        if let Some(ref tmp) = self.tmp_dir {
            sorter = sorter.temp_dir(tmp.clone());
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
        use crate::sort::raw_bam_reader::RawBamRecordReader;
        use crate::sort::{
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
            info!("Cell tag: {}{}", ct[0] as char, ct[1] as char);
        }

        // Get header using noodles reader, then use raw reader for records
        let (_, header) = create_bam_reader(&self.input, 1)?;

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
                    |bam| extract_template_key_inline(bam, &lib_lookup, cell_tag.as_ref(), &hasher),
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
    use rstest::rstest;

    /// Helper to construct a `Sort` struct with a given order.
    fn make_sort(order: SortOrderArg) -> Sort {
        Sort {
            input: PathBuf::from("test.bam"),
            output: None,
            verify: false,
            order,
            max_memory: MemoryLimit::Fixed(512 * 1024 * 1024),
            memory_reserve: MemoryReserve::Auto,
            memory_per_thread: true,
            tmp_dir: None,
            threads: 1,
            compression: CompressionOptions::default(),
            temp_compression: 1,
            write_index: false,
        }
    }

    #[rstest]
    #[case(SortOrderArg::TemplateCoordinate, Some(*SamTag::CB))]
    #[case(SortOrderArg::Coordinate, None)]
    #[case(SortOrderArg::Queryname, None)]
    fn test_parse_cell_tag(#[case] order: SortOrderArg, #[case] expected: Option<[u8; 2]>) {
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
    fn test_verify_sort_order_sorted() -> Result<()> {
        use crate::sam::builder::SamBuilder;
        use crate::sort::raw_bam_reader::RawBamRecordReader;

        let mut builder = SamBuilder::new();
        // Add records with names in sorted order
        let _ = builder.add_pair().name("aaa").build();
        let _ = builder.add_pair().name("bbb").build();
        let _ = builder.add_pair().name("ccc").build();

        let dir = tempfile::tempdir()?;
        let bam_path = dir.path().join("sorted.bam");
        builder.write_bam(&bam_path)?;

        let file = std::fs::File::open(&bam_path)?;
        let mut reader = RawBamRecordReader::new(file)?;
        reader.skip_header()?;

        let (total, violations, first_violation) = verify_sort_order(
            reader,
            |bam| fgumi_raw_bam::read_name(bam).to_vec(),
            |key, prev| key < prev,
        )?;

        assert_eq!(total, 6); // 3 pairs = 6 records
        assert_eq!(violations, 0);
        assert!(first_violation.is_none());
        Ok(())
    }

    #[test]
    fn test_verify_sort_order_unsorted() -> Result<()> {
        use crate::sam::builder::SamBuilder;
        use crate::sort::raw_bam_reader::RawBamRecordReader;

        let mut builder = SamBuilder::new();
        // Add records with names out of order (pairs are interleaved)
        let _ = builder.add_pair().name("ccc").build();
        let _ = builder.add_pair().name("aaa").build();
        let _ = builder.add_pair().name("bbb").build();

        let dir = tempfile::tempdir()?;
        let bam_path = dir.path().join("unsorted.bam");
        builder.write_bam(&bam_path)?;

        let file = std::fs::File::open(&bam_path)?;
        let mut reader = RawBamRecordReader::new(file)?;
        reader.skip_header()?;

        let (total, violations, first_violation) = verify_sort_order(
            reader,
            |bam| fgumi_raw_bam::read_name(bam).to_vec(),
            |key, prev| key < prev,
        )?;

        assert_eq!(total, 6);
        assert!(violations > 0);
        assert!(first_violation.is_some());
        let (record_num, name) =
            first_violation.expect("first violation should be present for unsorted file");
        assert!(record_num > 1); // violation can't be on first record
        assert!(!name.is_empty());
        Ok(())
    }

    #[test]
    fn test_verify_sort_order_empty() -> Result<()> {
        use crate::sam::builder::SamBuilder;
        use crate::sort::raw_bam_reader::RawBamRecordReader;

        let builder = SamBuilder::new();

        let dir = tempfile::tempdir()?;
        let bam_path = dir.path().join("empty.bam");
        builder.write_bam(&bam_path)?;

        let file = std::fs::File::open(&bam_path)?;
        let mut reader = RawBamRecordReader::new(file)?;
        reader.skip_header()?;

        let (total, violations, first_violation) = verify_sort_order(
            reader,
            |bam| fgumi_raw_bam::read_name(bam).to_vec(),
            |key, prev| key < prev,
        )?;

        assert_eq!(total, 0);
        assert_eq!(violations, 0);
        assert!(first_violation.is_none());
        Ok(())
    }

    // ========================================================================
    // Verify sort order tests (passing and failing) for all orders + sub-sorts
    // ========================================================================

    /// Helper: build a BAM, sort it with the given order, then verify it passes.
    fn sort_and_verify_pass(order_str: &str) -> Result<()> {
        use crate::sort::raw_bam_reader::RawBamRecordReader;
        use crate::sort::{
            RawExternalSorter, RawQuerynameKey, RawQuerynameLexKey, RawSortKey, SortContext,
            extract_coordinate_key_inline,
        };

        let mut builder = crate::sam::builder::SamBuilder::new();
        // Add pairs with names that sort differently under natural vs lexicographic
        let _ = builder.add_pair().name("read2").contig(0).start1(200).build();
        let _ = builder.add_pair().name("read10").contig(0).start1(100).build();
        let _ = builder.add_pair().name("read1").contig(1).start1(50).build();

        let dir = tempfile::tempdir()?;
        let input_bam = dir.path().join("input.bam");
        let sorted_bam = dir.path().join("sorted.bam");
        builder.write_bam(&input_bam)?;

        let order_arg =
            SortOrderArg::parse(order_str).expect("parse should succeed for valid sort order");
        let sort_order: SortOrder = order_arg.into();

        let sorter = RawExternalSorter::new(sort_order).threads(1).output_compression(6);
        sorter.sort(&input_bam, &sorted_bam)?;

        // Now verify
        let file = std::fs::File::open(&sorted_bam)?;
        let (_, header) = crate::bam_io::create_bam_reader(&sorted_bam, 1)?;
        let mut reader = RawBamRecordReader::new(file)?;
        reader.skip_header()?;

        let (total, violations, _) = match order_arg {
            SortOrderArg::Coordinate => {
                let nref = header.reference_sequences().len() as u32;
                verify_sort_order(
                    reader,
                    |bam| extract_coordinate_key_inline(bam, nref),
                    |key, prev| key < prev,
                )?
            }
            SortOrderArg::Queryname => {
                let ctx = SortContext::from_header(&header);
                verify_sort_order(
                    reader,
                    |bam| RawQuerynameLexKey::extract(bam, &ctx),
                    |key, prev| key < prev,
                )?
            }
            SortOrderArg::QuerynameNatural => {
                let ctx = SortContext::from_header(&header);
                verify_sort_order(
                    reader,
                    |bam| RawQuerynameKey::extract(bam, &ctx),
                    |key, prev| key < prev,
                )?
            }
            SortOrderArg::TemplateCoordinate => {
                // Skip template-coordinate verify — requires more complex setup
                return Ok(());
            }
        };

        assert!(total > 0, "should have records for {order_str}");
        assert_eq!(violations, 0, "should be sorted for {order_str}");
        Ok(())
    }

    #[test]
    fn test_verify_coordinate_sorted_pass() -> Result<()> {
        sort_and_verify_pass("coordinate")
    }

    #[test]
    fn test_verify_queryname_default_sorted_pass() -> Result<()> {
        sort_and_verify_pass("queryname")
    }

    #[test]
    fn test_verify_queryname_lexicographic_sorted_pass() -> Result<()> {
        sort_and_verify_pass("queryname::lexicographic")
    }

    #[test]
    fn test_verify_queryname_natural_sorted_pass() -> Result<()> {
        sort_and_verify_pass("queryname::natural")
    }

    #[test]
    fn test_verify_queryname_lex_fails_with_natural_verifier() -> Result<()> {
        use crate::sort::raw_bam_reader::RawBamRecordReader;
        use crate::sort::{RawExternalSorter, RawQuerynameKey, RawSortKey, SortContext};

        // Build BAM with names that differ between lex and natural ordering
        let mut builder = crate::sam::builder::SamBuilder::new();
        let _ = builder.add_pair().name("read2").contig(0).start1(100).build();
        let _ = builder.add_pair().name("read10").contig(0).start1(200).build();

        let dir = tempfile::tempdir()?;
        let input_bam = dir.path().join("input.bam");
        let sorted_bam = dir.path().join("sorted.bam");
        builder.write_bam(&input_bam)?;

        // Sort with lexicographic order (read10 < read2 in lex)
        let sorter =
            RawExternalSorter::new(SortOrder::Queryname(QuerynameComparator::Lexicographic))
                .threads(1)
                .output_compression(6);
        sorter.sort(&input_bam, &sorted_bam)?;

        // Verify with natural comparator — should detect violations
        // because lex order puts read10 before read2, but natural expects read2 before read10
        let file = std::fs::File::open(&sorted_bam)?;
        let (_, header) = crate::bam_io::create_bam_reader(&sorted_bam, 1)?;
        let mut reader = RawBamRecordReader::new(file)?;
        reader.skip_header()?;

        let ctx = SortContext::from_header(&header);
        let (total, violations, _) = verify_sort_order(
            reader,
            |bam| RawQuerynameKey::extract(bam, &ctx),
            |key, prev| key < prev,
        )?;

        assert!(total > 0);
        assert!(violations > 0, "lex-sorted file should fail natural verify");
        Ok(())
    }

    #[test]
    fn test_verify_queryname_natural_fails_with_lex_verifier() -> Result<()> {
        use crate::sort::raw_bam_reader::RawBamRecordReader;
        use crate::sort::{RawExternalSorter, RawQuerynameLexKey, RawSortKey, SortContext};

        // Build BAM with names that differ between lex and natural ordering
        let mut builder = crate::sam::builder::SamBuilder::new();
        let _ = builder.add_pair().name("read2").contig(0).start1(100).build();
        let _ = builder.add_pair().name("read10").contig(0).start1(200).build();

        let dir = tempfile::tempdir()?;
        let input_bam = dir.path().join("input.bam");
        let sorted_bam = dir.path().join("sorted.bam");
        builder.write_bam(&input_bam)?;

        // Sort with natural order (read2 < read10 in natural)
        let sorter = RawExternalSorter::new(SortOrder::Queryname(QuerynameComparator::Natural))
            .threads(1)
            .output_compression(6);
        sorter.sort(&input_bam, &sorted_bam)?;

        // Verify with lexicographic comparator — should detect violations
        // because natural puts read2 before read10, but lex expects read10 before read2
        let file = std::fs::File::open(&sorted_bam)?;
        let (_, header) = crate::bam_io::create_bam_reader(&sorted_bam, 1)?;
        let mut reader = RawBamRecordReader::new(file)?;
        reader.skip_header()?;

        let ctx = SortContext::from_header(&header);
        let (total, violations, _) = verify_sort_order(
            reader,
            |bam| RawQuerynameLexKey::extract(bam, &ctx),
            |key, prev| key < prev,
        )?;

        assert!(total > 0);
        assert!(violations > 0, "natural-sorted file should fail lex verify");
        Ok(())
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
        use crate::sort::extract_coordinate_key_inline;
        use crate::sort::raw_bam_reader::RawBamRecordReader;

        // Build BAM with records deliberately out of coordinate order
        let mut builder = crate::sam::builder::SamBuilder::new();
        let _ = builder.add_pair().name("a").contig(1).start1(100).build();
        let _ = builder.add_pair().name("b").contig(0).start1(200).build();

        let dir = tempfile::tempdir()?;
        let bam_path = dir.path().join("unsorted.bam");
        builder.write_bam(&bam_path)?;

        let file = std::fs::File::open(&bam_path)?;
        let (_, header) = crate::bam_io::create_bam_reader(&bam_path, 1)?;
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
