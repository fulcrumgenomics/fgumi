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

use anyhow::{Result, bail};
use clap::{Parser, ValueEnum};
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
  - Configurable memory limit (--max-memory)

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
"#
)]
pub struct Sort {
    /// Input BAM file.
    #[arg(short = 'i', long = "input")]
    pub input: PathBuf,

    /// Output BAM file.
    #[arg(short = 'o', long = "output")]
    pub output: PathBuf,

    /// Sort order.
    #[arg(long = "order", value_enum, default_value = "template-coordinate")]
    pub order: SortOrderArg,

    /// Maximum memory to use for in-memory sorting.
    ///
    /// Accepts values like "512M", "1G", "2G". When the memory limit
    /// is reached, sorted chunks are written to temporary files and
    /// merged at the end.
    #[arg(short = 'm', long = "max-memory", default_value = "512M", value_parser = parse_memory)]
    pub max_memory: usize,

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

    /// Compression level for temporary chunk files (0-12).
    ///
    /// Level 0 disables compression (fastest, uses most disk space).
    /// Level 1 (default) provides fast compression with reasonable space savings.
    /// Higher levels provide better compression but are slower.
    #[arg(long = "temp-compression", default_value = "1")]
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

        if self.max_memory == 0 {
            bail!("--max-memory must be greater than 0");
        }

        // --write-index only valid for coordinate sort
        if self.write_index && !matches!(self.order, SortOrderArg::Coordinate) {
            bail!("--write-index is only valid for coordinate sort");
        }

        let timer = OperationTimer::new("Sorting BAM");

        info!("Starting Sort");
        info!("Input: {}", self.input.display());
        info!("Output: {}", self.output.display());
        info!("Sort order: {:?}", self.order);
        info!("Max memory: {} MB", self.max_memory / (1024 * 1024));
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
            .memory_limit(self.max_memory)
            .threads(self.threads)
            .output_compression(self.compression.compression_level)
            .temp_compression(self.temp_compression)
            .write_index(self.write_index)
            .pg_info(crate::version::VERSION.to_string(), command_line.to_string());

        if let Some(ref tmp) = self.tmp_dir {
            sorter = sorter.temp_dir(tmp.clone());
        }

        let stats = sorter.sort(&self.input, &self.output)?;
        let (total_records, output_records, chunks_written) =
            (stats.total_records, stats.output_records, stats.chunks_written);

        // Summary
        info!("=== Summary ===");
        info!("Records processed: {total_records}");
        info!("Records written: {output_records}");
        if chunks_written > 0 {
            info!("Temporary chunks: {chunks_written}");
        }
        info!("Output: {}", self.output.display());

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
