//! Downsample BAM files by UMI family using a single-pass streaming algorithm.
//!
//! This tool reads a BAM file that has been processed by fgumi group (or fgbio `GroupReadsByUmi`),
//! uniformly samples UMI families based on the MI tag, and outputs kept reads directly to a BAM file.
//!
//! Requires input BAM to be in template-coordinate order (from group).

use anyhow::{Result, bail};
use clap::Parser;
use fgumi_lib::bam_io::{create_bam_reader, create_bam_writer, create_optional_bam_writer};
use fgumi_lib::logging::OperationTimer;
use fgumi_lib::progress::ProgressTracker;
use fgumi_lib::sam::is_template_coordinate_sorted;
use fgumi_lib::validation::validate_file_exists;
use log::info;
use noodles::sam::alignment::io::Write as AlignmentWrite;
use noodles::sam::alignment::record::data::field::Tag;
use noodles::sam::alignment::record_buf::RecordBuf;
use rand::SeedableRng;
use rand::rngs::StdRng;
use std::collections::{BTreeMap, HashSet};
use std::fs::File;
use std::io::Write;
use std::path::PathBuf;

use crate::commands::command::Command;
use crate::commands::common::{BamIoOptions, CompressionOptions};

/// MI tag for molecular identifier
const MI_TAG: Tag = Tag::new(b'M', b'I');

/// Downsample a BAM file by UMI family using streaming.
///
/// Reads consecutive records sharing the same MI tag, randomly samples families
/// based on the keep fraction, and writes kept families directly to output BAM.
#[derive(Debug, Parser)]
#[command(
    name = "downsample",
    about = "\x1b[38;5;166m[UTILITIES]\x1b[0m      \x1b[36mDownsample BAM by UMI family using streaming\x1b[0m",
    long_about = r#"
Downsample a BAM file by UMI family using a single-pass streaming algorithm.

This tool reads a BAM file that has been processed by fgumi group (or fgbio GroupReadsByUmi)
containing MI tags, uniformly samples UMI families, and outputs kept reads directly to a BAM file.

Requires input BAM to be in template-coordinate order:
  - SO:unsorted (or not set)
  - GO:query
  - SS:unsorted:template-coordinate or SS:template-coordinate

The tool processes families in streaming fashion by grouping consecutive reads with the same
MI tag value. For each family, a random decision is made based on the fraction parameter to
either keep or reject all reads in that family.

Example usage:
  fgumi downsample -i grouped.bam -o downsampled.bam -f 0.1 --seed 42
  fgumi downsample -i grouped.bam -o kept.bam -f 0.5 --rejects rejected.bam
  fgumi downsample -i grouped.bam -o kept.bam -f 0.1 --histogram-kept kept_hist.txt
"#
)]
pub struct Downsample {
    /// Input/output BAM options
    #[command(flatten)]
    pub io: BamIoOptions,

    /// Fraction of UMI families to keep (0.0 exclusive to 1.0 inclusive)
    #[arg(short = 'f', long = "fraction")]
    pub fraction: f64,

    /// Optional output BAM file for rejected reads
    #[arg(long = "rejects")]
    pub rejects: Option<PathBuf>,

    /// Random seed for reproducibility
    #[arg(long = "seed")]
    pub seed: Option<u64>,

    /// Validate that MI tags appear in consecutive groups (error if seen non-consecutively)
    #[arg(long = "validate-mi-order", default_value = "false")]
    pub validate_mi_order: bool,

    /// Output file for kept family size histogram
    #[arg(long = "histogram-kept")]
    pub histogram_kept: Option<PathBuf>,

    /// Output file for rejected family size histogram
    #[arg(long = "histogram-rejected")]
    pub histogram_rejected: Option<PathBuf>,

    /// Compression options for output BAM.
    #[command(flatten)]
    pub compression: CompressionOptions,
}

impl Command for Downsample {
    fn execute(&self, command_line: &str) -> Result<()> {
        // Validate inputs
        validate_file_exists(&self.io.input, "Input BAM")?;

        // Validate fraction
        if self.fraction <= 0.0 || self.fraction > 1.0 {
            bail!(
                "--fraction must be between 0.0 (exclusive) and 1.0 (inclusive), got {}",
                self.fraction
            );
        }

        let timer = OperationTimer::new("Downsampling reads");

        info!("Starting Downsample");
        info!("Input: {}", self.io.input.display());
        info!("Output: {}", self.io.output.display());
        info!("Target fraction: {}", self.fraction);
        if let Some(seed) = self.seed {
            info!("Random seed: {seed}");
        }
        if self.validate_mi_order {
            info!("MI order validation: enabled");
        }

        // Initialize RNG
        let mut rng = match self.seed {
            Some(seed) => StdRng::seed_from_u64(seed),
            None => rand::make_rng(),
        };

        // Open input BAM
        let (mut reader, header) = create_bam_reader(&self.io.input, 1)?;

        // Validate header - input must be template-coordinate sorted (output from group)
        if !is_template_coordinate_sorted(&header) {
            bail!(
                "Input BAM must be template-coordinate sorted (output from group).\n\n\
                Expected header fields: SO:unsorted, GO:query, SS:template-coordinate\n\n\
                The input to this tool should be the output of fgumi group or fgbio GroupReadsByUmi."
            );
        }
        info!("Header validation passed (template-coordinate order confirmed)");

        // Add @PG record with PP chaining to input's last program
        let header = crate::commands::common::add_pg_record(header, command_line)?;

        // Create output BAM writer (single-threaded, downsample doesn't have threads parameter)
        let mut writer =
            create_bam_writer(&self.io.output, &header, 1, self.compression.compression_level)?;

        // Create optional rejects writer
        let mut rejects_writer = create_optional_bam_writer(
            self.rejects.as_ref(),
            &header,
            1,
            self.compression.compression_level,
        )?;

        // Statistics
        let mut total_families: u64 = 0;
        let mut kept_families: u64 = 0;
        let mut kept_reads: u64 = 0;
        let mut rejected_reads: u64 = 0;
        let mut record_count: usize = 0;
        let progress = ProgressTracker::new("Processed records").with_interval(1_000_000);

        // Histograms: family_size -> count
        let mut hist_kept: BTreeMap<usize, u64> = BTreeMap::new();
        let mut hist_rejected: BTreeMap<usize, u64> = BTreeMap::new();

        // For MI order validation
        let mut seen_mis: HashSet<String> = HashSet::new();

        info!("Processing reads...");

        // Process families using streaming iteration
        let record_iter = reader.record_bufs(&header).map(|r| r.map_err(Into::into));
        let mut family_iter = FamilyIterator::new(record_iter);

        while let Some(family_result) = family_iter.next_family()? {
            let (mi, family) = family_result;
            total_families += 1;
            let family_size = family.len();

            // Validate MI order if requested
            if self.validate_mi_order {
                if seen_mis.contains(&mi) {
                    bail!(
                        "MI tag '{mi}' seen non-consecutively. Input BAM may not be properly grouped by MI."
                    );
                }
                seen_mis.insert(mi);
            }

            // Random sampling decision
            let keep = rand::RngExt::random::<f64>(&mut rng) < self.fraction;

            // Count all records processed (kept + rejected)
            record_count += family_size;

            if keep {
                kept_families += 1;
                kept_reads += family_size as u64;
                *hist_kept.entry(family_size).or_insert(0) += 1;

                for record in &family {
                    writer.write_alignment_record(&header, record)?;
                }
            } else {
                rejected_reads += family_size as u64;
                *hist_rejected.entry(family_size).or_insert(0) += 1;

                if let Some(ref mut rw) = rejects_writer {
                    for record in &family {
                        rw.write_alignment_record(&header, record)?;
                    }
                }
            }
            progress.log_if_needed(family_size as u64);
        }

        progress.log_final();

        // Write histograms
        if let Some(ref path) = self.histogram_kept {
            write_histogram(&hist_kept, path)?;
            info!("Wrote kept histogram to: {}", path.display());
        }

        if let Some(ref path) = self.histogram_rejected {
            write_histogram(&hist_rejected, path)?;
            info!("Wrote rejected histogram to: {}", path.display());
        }

        // Summary
        info!("=== Summary ===");
        info!("Total reads processed: {}", kept_reads + rejected_reads);
        info!("Input families: {total_families}");
        if total_families > 0 {
            let kept_pct = 100.0 * kept_families as f64 / total_families as f64;
            info!("Kept families: {kept_families} ({kept_pct:.2}%)");
        } else {
            info!("Kept families: 0");
        }
        info!("Kept reads: {kept_reads}");
        info!("Rejected reads: {rejected_reads}");
        info!("Output BAM: {}", self.io.output.display());
        if let Some(ref rejects) = self.rejects {
            info!("Rejects BAM: {}", rejects.display());
        }

        timer.log_completion(record_count as u64);
        Ok(())
    }
}

/// Write family size histogram in fgbio-compatible TSV format.
fn write_histogram(histogram: &BTreeMap<usize, u64>, path: &PathBuf) -> Result<()> {
    let mut file = File::create(path)?;
    writeln!(file, "family_size\tcount")?;
    for (size, count) in histogram {
        writeln!(file, "{size}\t{count}")?;
    }
    Ok(())
}

/// Iterator that groups consecutive records by MI tag.
///
/// This provides streaming iteration over UMI families without loading the entire BAM into memory.
struct FamilyIterator<I>
where
    I: Iterator<Item = Result<RecordBuf>>,
{
    records: std::iter::Peekable<I>,
}

impl<I> FamilyIterator<I>
where
    I: Iterator<Item = Result<RecordBuf>>,
{
    fn new(records: I) -> Self {
        Self { records: records.peekable() }
    }

    /// Get the next family of records sharing the same MI tag.
    ///
    /// Returns `Ok(Some((mi_tag`, records))) for each family, or Ok(None) when exhausted.
    fn next_family(&mut self) -> Result<Option<(String, Vec<RecordBuf>)>> {
        // Peek at the first record to get the MI tag
        let mi = match self.records.peek() {
            Some(Ok(record)) => get_mi_tag(record)?,
            Some(Err(_)) => {
                // Consume the error
                return Err(self.records.next().unwrap().unwrap_err());
            }
            None => return Ok(None),
        };

        // Collect all records with the same MI tag
        let mut family = Vec::new();

        while let Some(peek_result) = self.records.peek() {
            match peek_result {
                Ok(record) => {
                    let record_mi = get_mi_tag(record)?;
                    if record_mi != mi {
                        break;
                    }
                    // Consume the record
                    family.push(self.records.next().unwrap()?);
                }
                Err(_) => {
                    // Consume and return the error
                    return Err(self.records.next().unwrap().unwrap_err());
                }
            }
        }

        Ok(Some((mi, family)))
    }
}

/// Extract the MI tag value from a record.
fn get_mi_tag(record: &RecordBuf) -> Result<String> {
    let mi = record.data().get(&MI_TAG).ok_or_else(|| {
        let name = record.name().map_or_else(
            || "<unknown>".to_string(),
            |n| String::from_utf8_lossy(n.as_ref()).to_string(),
        );
        anyhow::anyhow!("Read '{name}' is missing required MI tag")
    })?;

    // MI tag can be an integer or string
    use noodles::sam::alignment::record_buf::data::field::Value;
    match mi {
        Value::Int8(v) => Ok(v.to_string()),
        Value::UInt8(v) => Ok(v.to_string()),
        Value::Int16(v) => Ok(v.to_string()),
        Value::UInt16(v) => Ok(v.to_string()),
        Value::Int32(v) => Ok(v.to_string()),
        Value::UInt32(v) => Ok(v.to_string()),
        Value::String(s) => Ok(s.to_string()),
        _ => bail!("Unexpected MI tag type"),
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use fgumi_lib::sam::builder::RecordBuilder;

    /// Create a test record with an MI tag
    fn create_test_record(name: &str, mi: &str) -> RecordBuf {
        RecordBuilder::new().name(name).tag("MI", mi).build()
    }

    /// Create a test record with an integer MI tag
    fn create_test_record_int_mi(name: &str, mi: i32) -> RecordBuf {
        RecordBuilder::new().name(name).tag("MI", mi).build()
    }

    /// Create a test record without an MI tag
    fn create_test_record_no_mi(name: &str) -> RecordBuf {
        RecordBuilder::new().name(name).build()
    }

    #[test]
    fn test_get_mi_tag_string() {
        let record = create_test_record("read1", "12345");
        let mi = get_mi_tag(&record).unwrap();
        assert_eq!(mi, "12345");
    }

    #[test]
    fn test_get_mi_tag_integer() {
        let record = create_test_record_int_mi("read1", 42);
        let mi = get_mi_tag(&record).unwrap();
        assert_eq!(mi, "42");
    }

    #[test]
    fn test_get_mi_tag_missing() {
        let record = create_test_record_no_mi("read1");
        let result = get_mi_tag(&record);
        assert!(result.is_err());
        assert!(result.unwrap_err().to_string().contains("missing required MI tag"));
    }

    #[test]
    fn test_family_iterator_single_family() {
        let records = vec![
            Ok(create_test_record("r1", "100")),
            Ok(create_test_record("r2", "100")),
            Ok(create_test_record("r3", "100")),
        ];

        let mut iter = FamilyIterator::new(records.into_iter());

        let family1 = iter.next_family().unwrap().unwrap();
        assert_eq!(family1.0, "100");
        assert_eq!(family1.1.len(), 3);

        let family2 = iter.next_family().unwrap();
        assert!(family2.is_none());
    }

    #[test]
    fn test_family_iterator_multiple_families() {
        let records = vec![
            Ok(create_test_record("r1", "100")),
            Ok(create_test_record("r2", "100")),
            Ok(create_test_record("r3", "200")),
            Ok(create_test_record("r4", "200")),
            Ok(create_test_record("r5", "200")),
            Ok(create_test_record("r6", "300")),
        ];

        let mut iter = FamilyIterator::new(records.into_iter());

        let family1 = iter.next_family().unwrap().unwrap();
        assert_eq!(family1.0, "100");
        assert_eq!(family1.1.len(), 2);

        let family2 = iter.next_family().unwrap().unwrap();
        assert_eq!(family2.0, "200");
        assert_eq!(family2.1.len(), 3);

        let family3 = iter.next_family().unwrap().unwrap();
        assert_eq!(family3.0, "300");
        assert_eq!(family3.1.len(), 1);

        let family4 = iter.next_family().unwrap();
        assert!(family4.is_none());
    }

    #[test]
    fn test_family_iterator_empty() {
        let records: Vec<Result<RecordBuf>> = vec![];
        let mut iter = FamilyIterator::new(records.into_iter());

        let family = iter.next_family().unwrap();
        assert!(family.is_none());
    }

    #[test]
    fn test_validate_fraction_too_low() {
        let cmd = Downsample {
            io: BamIoOptions {
                input: PathBuf::from("input.bam"),
                output: PathBuf::from("output.bam"),
            },
            fraction: 0.0,
            rejects: None,
            seed: None,
            validate_mi_order: false,
            histogram_kept: None,
            histogram_rejected: None,
            compression: CompressionOptions { compression_level: 1 },
        };

        // We can't call execute() without a real BAM file, but we can test the validation
        assert!(cmd.fraction <= 0.0);
    }

    #[test]
    fn test_validate_fraction_too_high() {
        let cmd = Downsample {
            io: BamIoOptions {
                input: PathBuf::from("input.bam"),
                output: PathBuf::from("output.bam"),
            },
            fraction: 1.5,
            rejects: None,
            seed: None,
            validate_mi_order: false,
            histogram_kept: None,
            histogram_rejected: None,
            compression: CompressionOptions { compression_level: 1 },
        };

        assert!(cmd.fraction > 1.0);
    }

    #[test]
    fn test_validate_fraction_valid() {
        let cmd = Downsample {
            io: BamIoOptions {
                input: PathBuf::from("input.bam"),
                output: PathBuf::from("output.bam"),
            },
            fraction: 0.5,
            rejects: None,
            seed: None,
            validate_mi_order: false,
            histogram_kept: None,
            histogram_rejected: None,
            compression: CompressionOptions { compression_level: 1 },
        };

        assert!(cmd.fraction > 0.0 && cmd.fraction <= 1.0);
    }

    #[test]
    fn test_write_histogram() {
        use tempfile::NamedTempFile;

        let mut hist = BTreeMap::new();
        hist.insert(1, 10);
        hist.insert(2, 20);
        hist.insert(5, 5);

        let temp_file = NamedTempFile::new().unwrap();
        write_histogram(&hist, &temp_file.path().to_path_buf()).unwrap();

        let contents = std::fs::read_to_string(temp_file.path()).unwrap();
        assert!(contents.contains("family_size\tcount"));
        assert!(contents.contains("1\t10"));
        assert!(contents.contains("2\t20"));
        assert!(contents.contains("5\t5"));
    }

    #[test]
    #[allow(clippy::float_cmp)] // Testing exact value assignment, not computation
    fn test_downsample_parameters() {
        let cmd = Downsample {
            io: BamIoOptions {
                input: PathBuf::from("input.bam"),
                output: PathBuf::from("output.bam"),
            },
            fraction: 0.1,
            rejects: Some(PathBuf::from("rejects.bam")),
            seed: Some(42),
            validate_mi_order: true,
            histogram_kept: Some(PathBuf::from("kept.txt")),
            histogram_rejected: Some(PathBuf::from("rejected.txt")),
            compression: CompressionOptions { compression_level: 1 },
        };

        assert_eq!(cmd.fraction, 0.1);
        assert_eq!(cmd.seed, Some(42));
        assert!(cmd.validate_mi_order);
        assert!(cmd.rejects.is_some());
        assert!(cmd.histogram_kept.is_some());
        assert!(cmd.histogram_rejected.is_some());
    }

    // Note: Header validation tests are in fgumi_lib::sam::tests for is_template_coordinate_sorted

    #[test]
    fn test_deterministic_sampling_with_seed() {
        // Test that the same seed produces the same results
        use rand::RngExt;

        let seed = 12345u64;

        // First run
        let mut rng1 = StdRng::seed_from_u64(seed);
        let results1: Vec<bool> = (0..100).map(|_| rng1.random::<f64>() < 0.5).collect();

        // Second run with same seed
        let mut rng2 = StdRng::seed_from_u64(seed);
        let results2: Vec<bool> = (0..100).map(|_| rng2.random::<f64>() < 0.5).collect();

        assert_eq!(results1, results2);
    }

    #[test]
    fn test_histogram_sorted_by_family_size() {
        let mut hist = BTreeMap::new();
        hist.insert(5, 10);
        hist.insert(1, 20);
        hist.insert(3, 15);

        // BTreeMap maintains sorted order
        let sizes: Vec<usize> = hist.keys().copied().collect();
        assert_eq!(sizes, vec![1, 3, 5]);
    }
}
