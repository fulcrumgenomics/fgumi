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
use clap::{Parser, ValueEnum};
use fgumi_lib::bam_io::create_bam_reader;
use fgumi_lib::logging::OperationTimer;
use fgumi_lib::sort::{RawExternalSorter, SortOrder};
use fgumi_lib::validation::validate_file_exists;
use log::info;
use std::path::PathBuf;

use crate::commands::command::Command;
use crate::commands::common::CompressionOptions;

/// Sort order for BAM files.
#[derive(Debug, Clone, Copy, ValueEnum)]
pub enum SortOrderArg {
    /// Coordinate sort (tid → pos → strand)
    Coordinate,
    /// Queryname sort (read name with natural ordering)
    Queryname,
    /// Template-coordinate sort (for UMI grouping)
    TemplateCoordinate,
}

impl From<SortOrderArg> for SortOrder {
    fn from(arg: SortOrderArg) -> Self {
        match arg {
            SortOrderArg::Coordinate => SortOrder::Coordinate,
            SortOrderArg::Queryname => SortOrder::Queryname,
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

  coordinate          Standard genomic coordinate sort (tid → pos → strand).
                      Use for IGV visualization, variant calling, `fgumi review`.

  queryname           Read name sort with natural numeric ordering.
                      Use for `fgumi zipper`, template-level operations.

  template-coordinate Template-level position sort for UMI grouping.
                      Use for `fgumi group` input. Equivalent to:
                        samtools sort --template-coordinate

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
pub struct Sort {
    /// Input BAM file.
    #[arg(short = 'i', long = "input")]
    pub input: PathBuf,

    /// Output BAM file (required unless --verify is used).
    #[arg(short = 'o', long = "output", conflicts_with = "verify")]
    pub output: Option<PathBuf>,

    /// Verify the input file is correctly sorted (no output written).
    ///
    /// Reads records sequentially and checks that each record's sort key
    /// is >= the previous record's key. Exits 0 if sorted correctly,
    /// non-zero if any records are out of order.
    #[arg(long = "verify", conflicts_with = "output")]
    pub verify: bool,

    /// Sort order.
    #[arg(long = "order", value_enum, default_value = "template-coordinate")]
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
    #[arg(long = "memory-per-thread", default_value = "true", action = clap::ArgAction::Set)]
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
    /// <output>.bai. Uses single-threaded compression for accurate virtual
    /// position tracking.
    #[arg(long = "write-index", default_value = "false")]
    pub write_index: bool,
}

/// Parse memory size string (e.g., "512M", "1G", "2G").
fn parse_memory(s: &str) -> Result<usize, String> {
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

impl Command for Sort {
    fn execute(&self, command_line: &str) -> Result<()> {
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

        info!("Starting Sort");
        info!("Input: {}", self.input.display());
        info!("Output: {}", output.display());
        info!("Sort order: {:?}", self.order);
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
            LibraryLookup, RawQuerynameKey, RawSortKey, SortContext, TemplateKey,
            extract_coordinate_key_inline, extract_template_key_inline,
        };
        use std::cmp::Ordering;
        use std::fs::File;

        let timer = OperationTimer::new("Verifying BAM sort order");

        info!("Starting Sort Verification");
        info!("Input: {}", self.input.display());
        info!("Expected order: {:?}", self.order);

        // Get header using noodles reader, then use raw reader for records
        let (_, header) = create_bam_reader(&self.input, 1)?;

        let mut total_records: u64 = 0;
        let mut violations: u64 = 0;
        let mut first_violation: Option<(u64, String)> = None;

        // Helper to extract name from raw BAM bytes
        let extract_name = |bam: &[u8]| -> String {
            let l_read_name = bam.get(8).copied().unwrap_or(0) as usize;
            let name_len = l_read_name.saturating_sub(1);
            if name_len > 0 && 32 + name_len <= bam.len() {
                String::from_utf8_lossy(&bam[32..32 + name_len]).to_string()
            } else {
                "<unnamed>".to_string()
            }
        };

        match self.order {
            SortOrderArg::Coordinate => {
                // Use raw BAM reading with same key extraction as sort
                let file = File::open(&self.input)?;
                let mut raw_reader = RawBamRecordReader::new(file)?;
                raw_reader.skip_header()?;
                let nref = header.reference_sequences().len() as u32;

                let mut prev_key: Option<u64> = None;

                for result in raw_reader {
                    let record_bytes = result?;
                    total_records += 1;
                    let bam = record_bytes.as_slice();

                    let key = extract_coordinate_key_inline(bam, nref);

                    if let Some(prev) = prev_key {
                        if key < prev {
                            violations += 1;
                            if first_violation.is_none() {
                                first_violation = Some((total_records, extract_name(bam)));
                            }
                        }
                    }
                    prev_key = Some(key);
                }
            }
            SortOrderArg::Queryname => {
                // Use raw BAM reading with same key extraction as sort
                let file = File::open(&self.input)?;
                let mut raw_reader = RawBamRecordReader::new(file)?;
                raw_reader.skip_header()?;
                let ctx = SortContext::from_header(&header);

                let mut prev_key: Option<RawQuerynameKey> = None;

                for result in raw_reader {
                    let record_bytes = result?;
                    total_records += 1;
                    let bam = record_bytes.as_slice();

                    let key = RawQuerynameKey::extract(bam, &ctx);

                    if let Some(ref prev) = prev_key {
                        if key < *prev {
                            violations += 1;
                            if first_violation.is_none() {
                                first_violation = Some((total_records, extract_name(bam)));
                            }
                        }
                    }
                    prev_key = Some(key);
                }
            }
            SortOrderArg::TemplateCoordinate => {
                // Use raw BAM reading with same key extraction as sort
                let file = File::open(&self.input)?;
                let mut raw_reader = RawBamRecordReader::new(file)?;
                raw_reader.skip_header()?;
                let lib_lookup = LibraryLookup::from_header(&header);

                let mut prev_key: Option<TemplateKey> = None;

                for result in raw_reader {
                    let record_bytes = result?;
                    total_records += 1;
                    let bam = record_bytes.as_slice();

                    let key = extract_template_key_inline(bam, &lib_lookup);

                    if let Some(ref prev) = prev_key {
                        // Use core_cmp to ignore name_hash tie-breaker differences
                        // This allows both fgumi and samtools sorted files to pass
                        if key.core_cmp(prev) == Ordering::Less {
                            violations += 1;
                            if first_violation.is_none() {
                                first_violation = Some((total_records, extract_name(bam)));
                            }
                        }
                    }
                    prev_key = Some(key);
                }
            }
        }

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

    #[test]
    fn test_parse_memory_megabytes() {
        assert_eq!(parse_memory("512M").unwrap(), 512 * 1024 * 1024);
        assert_eq!(parse_memory("1024M").unwrap(), 1024 * 1024 * 1024);
    }

    #[test]
    fn test_parse_memory_gigabytes() {
        assert_eq!(parse_memory("1G").unwrap(), 1024 * 1024 * 1024);
        assert_eq!(parse_memory("2G").unwrap(), 2 * 1024 * 1024 * 1024);
    }

    #[test]
    fn test_parse_memory_kilobytes() {
        assert_eq!(parse_memory("1024K").unwrap(), 1024 * 1024);
    }

    #[test]
    fn test_parse_memory_bytes() {
        assert_eq!(parse_memory("1048576").unwrap(), 1_048_576);
    }

    #[test]
    fn test_parse_memory_lowercase() {
        assert_eq!(parse_memory("512m").unwrap(), 512 * 1024 * 1024);
        assert_eq!(parse_memory("1g").unwrap(), 1024 * 1024 * 1024);
    }

    #[test]
    fn test_parse_memory_decimal() {
        assert_eq!(parse_memory("1.5G").unwrap(), (1.5 * 1024.0 * 1024.0 * 1024.0) as usize);
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
        assert_eq!(SortOrder::from(SortOrderArg::Queryname), SortOrder::Queryname);
        assert_eq!(
            SortOrder::from(SortOrderArg::TemplateCoordinate),
            SortOrder::TemplateCoordinate
        );
    }
}
