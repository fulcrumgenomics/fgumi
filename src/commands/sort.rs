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

use anyhow::{Result, bail};
use clap::Parser;
use fgumi_lib::bam_io::create_bam_reader;
use fgumi_lib::logging::OperationTimer;
use fgumi_lib::sort::{QuerynameComparator, RawExternalSorter, SortOrder};
use fgumi_lib::validation::{string_to_tag, validate_file_exists};
use log::info;
use std::path::PathBuf;

use crate::commands::command::Command;
use crate::commands::common::{CompressionOptions, parse_bool};

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
  - Memory scales with threads by default (768MB/thread, like samtools)

EXAMPLES:

  # Sort for fgumi group input
  fgumi sort -i aligned.bam -o sorted.bam --order template-coordinate

  # Sort by coordinate for IGV
  fgumi sort -i input.bam -o sorted.bam --order coordinate

  # Sort by queryname for zipper
  fgumi sort -i input.bam -o sorted.bam --order queryname

  # High-performance sort with more memory and threads
  fgumi sort -i input.bam -o sorted.bam --order template-coordinate \
    --max-memory 8G --threads 8

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

    /// Maximum memory to use for in-memory sorting (per thread when --memory-per-thread).
    ///
    /// Accepts values like "512M", "1G", "2G". When the memory limit
    /// is reached, sorted chunks are written to temporary files and
    /// merged at the end.
    ///
    /// With --memory-per-thread (default), this is multiplied by thread count
    /// (like samtools' -m option). With 8 threads and 768M, total = 6GB.
    #[arg(short = 'm', long = "max-memory", default_value = "768M", value_parser = parse_memory)]
    pub max_memory: usize,

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

    /// Cell barcode tag for template-coordinate sort.
    ///
    /// When sorting in template-coordinate order, this tag is included in the
    /// sort key so that templates from different cells at the same locus are
    /// not interleaved. Use the same value as `--cell-tag` in `fgumi group`
    /// and `fgumi dedup`. Only used for template-coordinate sort.
    #[arg(short = 'c', long = "cell-tag", default_value = "CB")]
    pub cell_tag: String,
}

/// Parse memory size string (e.g., "512M", "1G", "2G").
pub(crate) fn parse_memory(s: &str) -> Result<usize, String> {
    let s = s.trim().to_uppercase();

    if s.is_empty() {
        return Err("Empty memory specification".to_string());
    }

    let (num_str, multiplier) = if s.ends_with('G') {
        (&s[..s.len() - 1], 1024 * 1024 * 1024)
    } else if s.ends_with('M') {
        (&s[..s.len() - 1], 1024 * 1024)
    } else if s.ends_with('K') {
        (&s[..s.len() - 1], 1024)
    } else {
        // Assume bytes
        (s.as_str(), 1)
    };

    let num: f64 = num_str.parse().map_err(|_| format!("Invalid number: {num_str}"))?;

    if num < 0.0 {
        return Err("Memory size must be positive".to_string());
    }

    Ok((num * multiplier as f64) as usize)
}

/// Parse the cell tag for template-coordinate sort/verify, returning `None`
/// for other sort orders.
pub(crate) fn parse_cell_tag(order: SortOrderArg, cell_tag: &str) -> Result<Option<[u8; 2]>> {
    if matches!(order, SortOrderArg::TemplateCoordinate) {
        let tag = string_to_tag(cell_tag, "cell-tag")?;
        Ok(Some([tag.as_ref()[0], tag.as_ref()[1]]))
    } else {
        Ok(None)
    }
}

/// Summary of sort-order verification: `(total_records, violations, first_violation)`.
type VerifySummary = (u64, u64, Option<(u64, String)>);

/// Verify that records from a raw BAM reader are in sorted order.
///
/// Iterates all records, extracting a sort key from each and checking that
/// consecutive keys satisfy the ordering invariant (no violations).
fn verify_sort_order<K>(
    raw_reader: fgumi_lib::sort::raw_bam_reader::RawBamRecordReader<std::fs::File>,
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

impl Sort {
    /// Parse the cell tag for template-coordinate sort/verify, returning `None`
    /// for other sort orders.
    fn parse_cell_tag(&self) -> Result<Option<[u8; 2]>> {
        parse_cell_tag(self.order, &self.cell_tag)
    }

    /// Execute sort mode: read, sort, and write output.
    fn execute_sort(&self, command_line: &str) -> Result<()> {
        let output = self.output.as_ref().expect("output required for sort mode");

        if self.max_memory == 0 {
            bail!("--max-memory must be greater than 0");
        }

        // --write-index only valid for coordinate sort
        if self.write_index && !matches!(self.order, SortOrderArg::Coordinate) {
            bail!("--write-index is only valid for coordinate sort");
        }

        let timer = OperationTimer::new("Sorting BAM");

        // Calculate effective memory: per-thread or fixed total
        let effective_memory = if self.memory_per_thread {
            self.max_memory.saturating_mul(self.threads.max(1))
        } else {
            self.max_memory
        };

        let cell_tag = self.parse_cell_tag()?;

        info!("Starting Sort");
        info!("Input: {}", self.input.display());
        info!("Output: {}", output.display());
        info!("Sort order: {:?}", self.order);
        if let Some(ct) = cell_tag {
            info!("Cell tag: {}{}", ct[0] as char, ct[1] as char);
        }
        if self.memory_per_thread {
            info!(
                "Max memory: {} MB ({} MB/thread × {} threads)",
                effective_memory / (1024 * 1024),
                self.max_memory / (1024 * 1024),
                self.threads
            );
        } else {
            info!("Max memory: {} MB (fixed)", effective_memory / (1024 * 1024));
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
        use fgumi_lib::sort::raw_bam_reader::RawBamRecordReader;
        use fgumi_lib::sort::{
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

    /// Helper to construct a `Sort` struct with a given order and cell tag.
    fn make_sort(order: SortOrderArg, cell_tag: &str) -> Sort {
        Sort {
            input: PathBuf::from("test.bam"),
            output: None,
            verify: false,
            order,
            max_memory: 512 * 1024 * 1024,
            memory_per_thread: true,
            tmp_dir: None,
            threads: 1,
            compression: CompressionOptions::default(),
            temp_compression: 1,
            write_index: false,
            cell_tag: cell_tag.to_string(),
        }
    }

    #[rstest]
    #[case(SortOrderArg::TemplateCoordinate, "CB", Some(*b"CB"))]
    #[case(SortOrderArg::TemplateCoordinate, "CR", Some(*b"CR"))]
    #[case(SortOrderArg::Coordinate, "CB", None)]
    #[case(SortOrderArg::Queryname, "CB", None)]
    fn test_parse_cell_tag(
        #[case] order: SortOrderArg,
        #[case] cell_tag: &str,
        #[case] expected: Option<[u8; 2]>,
    ) {
        let sort = make_sort(order, cell_tag);
        assert_eq!(sort.parse_cell_tag().expect("parse_cell_tag should succeed"), expected);
    }

    #[test]
    fn test_parse_cell_tag_invalid_single_char() {
        let sort = make_sort(SortOrderArg::TemplateCoordinate, "X");
        assert!(sort.parse_cell_tag().is_err());
    }

    #[test]
    fn test_parse_cell_tag_invalid_too_long() {
        let sort = make_sort(SortOrderArg::TemplateCoordinate, "ABC");
        assert!(sort.parse_cell_tag().is_err());
    }

    #[test]
    fn test_parse_memory_megabytes() {
        assert_eq!(
            parse_memory("512M").expect("parse_memory should succeed for 512M"),
            512 * 1024 * 1024
        );
        assert_eq!(
            parse_memory("1024M").expect("parse_memory should succeed for 1024M"),
            1024 * 1024 * 1024
        );
    }

    #[test]
    fn test_parse_memory_gigabytes() {
        assert_eq!(
            parse_memory("1G").expect("parse_memory should succeed for 1G"),
            1024 * 1024 * 1024
        );
        assert_eq!(
            parse_memory("2G").expect("parse_memory should succeed for 2G"),
            2 * 1024 * 1024 * 1024
        );
    }

    #[test]
    fn test_parse_memory_kilobytes() {
        assert_eq!(
            parse_memory("1024K").expect("parse_memory should succeed for 1024K"),
            1024 * 1024
        );
    }

    #[test]
    fn test_parse_memory_bytes() {
        assert_eq!(
            parse_memory("1048576").expect("parse_memory should succeed for bare bytes"),
            1_048_576
        );
    }

    #[test]
    fn test_parse_memory_lowercase() {
        assert_eq!(
            parse_memory("512m").expect("parse_memory should succeed for lowercase 512m"),
            512 * 1024 * 1024
        );
        assert_eq!(
            parse_memory("1g").expect("parse_memory should succeed for lowercase 1g"),
            1024 * 1024 * 1024
        );
    }

    #[test]
    fn test_parse_memory_decimal() {
        assert_eq!(
            parse_memory("1.5G").expect("parse_memory should succeed for 1.5G"),
            (1.5 * 1024.0 * 1024.0 * 1024.0) as usize
        );
    }

    #[test]
    fn test_parse_memory_invalid() {
        assert!(parse_memory("").is_err());
        assert!(parse_memory("abc").is_err());
        assert!(parse_memory("-1G").is_err());
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
        use fgumi_lib::sam::builder::SamBuilder;
        use fgumi_lib::sort::raw_bam_reader::RawBamRecordReader;

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
        use fgumi_lib::sam::builder::SamBuilder;
        use fgumi_lib::sort::raw_bam_reader::RawBamRecordReader;

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
        use fgumi_lib::sam::builder::SamBuilder;
        use fgumi_lib::sort::raw_bam_reader::RawBamRecordReader;

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
        use fgumi_lib::sort::raw_bam_reader::RawBamRecordReader;
        use fgumi_lib::sort::{
            RawExternalSorter, RawQuerynameKey, RawQuerynameLexKey, RawSortKey, SortContext,
            extract_coordinate_key_inline,
        };

        let mut builder = fgumi_lib::sam::builder::SamBuilder::new();
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
        let (_, header) = fgumi_lib::bam_io::create_bam_reader(&sorted_bam, 1)?;
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
        use fgumi_lib::sort::raw_bam_reader::RawBamRecordReader;
        use fgumi_lib::sort::{RawExternalSorter, RawQuerynameKey, RawSortKey, SortContext};

        // Build BAM with names that differ between lex and natural ordering
        let mut builder = fgumi_lib::sam::builder::SamBuilder::new();
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
        let (_, header) = fgumi_lib::bam_io::create_bam_reader(&sorted_bam, 1)?;
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
        use fgumi_lib::sort::raw_bam_reader::RawBamRecordReader;
        use fgumi_lib::sort::{RawExternalSorter, RawQuerynameLexKey, RawSortKey, SortContext};

        // Build BAM with names that differ between lex and natural ordering
        let mut builder = fgumi_lib::sam::builder::SamBuilder::new();
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
        let (_, header) = fgumi_lib::bam_io::create_bam_reader(&sorted_bam, 1)?;
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
            ..make_sort(SortOrderArg::Coordinate, "CB")
        };
        let err = sort.execute("test").unwrap_err();
        assert!(err.to_string().contains("--verify cannot be used with --output"));
    }

    #[test]
    fn test_verify_conflicts_with_write_index() {
        let sort =
            Sort { verify: true, write_index: true, ..make_sort(SortOrderArg::Coordinate, "CB") };
        let err = sort.execute("test").unwrap_err();
        assert!(err.to_string().contains("--write-index cannot be used with --verify"));
    }

    #[test]
    fn test_verify_coordinate_fails_on_unsorted() -> Result<()> {
        use fgumi_lib::sort::extract_coordinate_key_inline;
        use fgumi_lib::sort::raw_bam_reader::RawBamRecordReader;

        // Build BAM with records deliberately out of coordinate order
        let mut builder = fgumi_lib::sam::builder::SamBuilder::new();
        let _ = builder.add_pair().name("a").contig(1).start1(100).build();
        let _ = builder.add_pair().name("b").contig(0).start1(200).build();

        let dir = tempfile::tempdir()?;
        let bam_path = dir.path().join("unsorted.bam");
        builder.write_bam(&bam_path)?;

        let file = std::fs::File::open(&bam_path)?;
        let (_, header) = fgumi_lib::bam_io::create_bam_reader(&bam_path, 1)?;
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
}
