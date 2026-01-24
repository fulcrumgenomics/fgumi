//! Gzipped FASTQ file writer for simulation output.
//!
//! Provides a simple interface for writing FASTQ records to gzip-compressed files.
//! Supports both single-threaded and multi-threaded parallel compression.

use super::parallel_gzip_writer::{ParallelGzipConfig, ParallelGzipWriter};
use anyhow::{Context, Result};
use flate2::Compression;
use flate2::write::GzEncoder;
use std::fs::File;
use std::io::{BufWriter, Write};
use std::path::Path;

/// A writer for gzip-compressed FASTQ files.
///
/// Writes FASTQ records in standard format with Phred+33 quality encoding.
/// Supports both single-threaded and multi-threaded parallel compression.
///
/// # Examples
///
/// ```no_run
/// use fgumi_lib::simulate::FastqWriter;
///
/// # fn main() -> anyhow::Result<()> {
/// // Single-threaded writer
/// let mut writer = FastqWriter::new("output.fastq.gz")?;
///
/// // Multi-threaded writer with 4 compression threads
/// let mut mt_writer = FastqWriter::with_threads("output_mt.fastq.gz", 4)?;
///
/// // Write a record with numeric quality scores (will be converted to Phred+33)
/// writer.write_record("read1", b"ACGTACGT", &[30, 30, 30, 30, 30, 30, 30, 30])?;
///
/// // Important: call finish() to ensure all data is written
/// writer.finish()?;
/// # Ok(())
/// # }
/// ```
pub struct FastqWriter {
    inner: FastqWriterInner,
}

enum FastqWriterInner {
    SingleThreaded(GzEncoder<BufWriter<File>>),
    MultiThreaded(ParallelGzipWriter),
}

impl FastqWriter {
    /// Create a new single-threaded FASTQ writer for the given path.
    ///
    /// The output will be gzip-compressed with default compression level.
    ///
    /// # Arguments
    ///
    /// * `path` - Path to the output FASTQ.gz file
    ///
    /// # Errors
    ///
    /// Returns an error if the file cannot be created.
    pub fn new<P: AsRef<Path>>(path: P) -> Result<Self> {
        let path = path.as_ref();
        let file =
            File::create(path).with_context(|| format!("Failed to create {}", path.display()))?;
        let buf = BufWriter::new(file);
        let gz = GzEncoder::new(buf, Compression::default());
        Ok(Self { inner: FastqWriterInner::SingleThreaded(gz) })
    }

    /// Create a new FASTQ writer with specified thread count.
    ///
    /// Uses parallel gzip compression when threads > 1.
    ///
    /// # Arguments
    ///
    /// * `path` - Path to the output FASTQ.gz file
    /// * `threads` - Number of compression threads (1 = single-threaded)
    ///
    /// # Errors
    ///
    /// Returns an error if the file cannot be created.
    pub fn with_threads<P: AsRef<Path>>(path: P, threads: usize) -> Result<Self> {
        let path = path.as_ref();
        let file =
            File::create(path).with_context(|| format!("Failed to create {}", path.display()))?;

        if threads <= 1 {
            let buf = BufWriter::new(file);
            let gz = GzEncoder::new(buf, Compression::default());
            Ok(Self { inner: FastqWriterInner::SingleThreaded(gz) })
        } else {
            let config = ParallelGzipConfig::with_threads(threads);
            let writer = ParallelGzipWriter::new(file, config)
                .with_context(|| "Failed to create parallel gzip writer")?;
            Ok(Self { inner: FastqWriterInner::MultiThreaded(writer) })
        }
    }

    /// Write a FASTQ record.
    ///
    /// # Arguments
    ///
    /// * `name` - Read name (without the leading '@')
    /// * `seq` - DNA sequence as bytes (A, C, G, T, N)
    /// * `qual` - Quality scores as numeric Phred values (0-41)
    ///
    /// Quality scores are automatically converted to Phred+33 ASCII encoding.
    ///
    /// # Errors
    ///
    /// Returns an error if writing fails.
    pub fn write_record(&mut self, name: &str, seq: &[u8], qual: &[u8]) -> Result<()> {
        match &mut self.inner {
            FastqWriterInner::SingleThreaded(writer) => write_record_to(writer, name, seq, qual),
            FastqWriterInner::MultiThreaded(writer) => write_record_to(writer, name, seq, qual),
        }
    }

    /// Finish writing and flush all data.
    ///
    /// This must be called to ensure all data is written and the gzip stream
    /// is properly terminated.
    ///
    /// # Errors
    ///
    /// Returns an error if flushing fails.
    pub fn finish(self) -> Result<()> {
        match self.inner {
            FastqWriterInner::SingleThreaded(writer) => {
                writer.finish().context("Failed to finish gzip stream")?;
            }
            FastqWriterInner::MultiThreaded(writer) => {
                writer.finish().context("Failed to finish parallel gzip stream")?;
            }
        }
        Ok(())
    }
}

/// Write a FASTQ record to any writer.
fn write_record_to<W: Write>(writer: &mut W, name: &str, seq: &[u8], qual: &[u8]) -> Result<()> {
    // Write header line
    writeln!(writer, "@{name}")?;

    // Write sequence
    writer.write_all(seq)?;
    writeln!(writer)?;

    // Write separator
    writeln!(writer, "+")?;

    // Write quality scores (convert numeric to Phred+33 ASCII)
    for &q in qual {
        let ascii_q = q.saturating_add(33).min(126);
        writer.write_all(&[ascii_q])?;
    }
    writeln!(writer)?;

    Ok(())
}

#[cfg(test)]
mod tests {
    use super::*;
    use std::io::Read;
    use tempfile::NamedTempFile;

    #[test]
    fn test_write_single_record() -> Result<()> {
        let temp = NamedTempFile::new()?;
        let path = temp.path();

        {
            let mut writer = FastqWriter::new(path)?;
            writer.write_record("read1", b"ACGT", &[30, 30, 30, 30])?;
            writer.finish()?;
        }

        // Read and decompress
        let file = File::open(path)?;
        let mut decoder = flate2::read::GzDecoder::new(file);
        let mut content = String::new();
        decoder.read_to_string(&mut content)?;

        assert!(content.contains("@read1"));
        assert!(content.contains("ACGT"));
        assert!(content.contains("????")); // Q30 + 33 = 63 = '?'

        Ok(())
    }

    #[test]
    fn test_write_multiple_records() -> Result<()> {
        let temp = NamedTempFile::new()?;
        let path = temp.path();

        {
            let mut writer = FastqWriter::new(path)?;
            writer.write_record("read1", b"AAAA", &[10, 20, 30, 40])?;
            writer.write_record("read2", b"CCCC", &[40, 30, 20, 10])?;
            writer.finish()?;
        }

        let file = File::open(path)?;
        let mut decoder = flate2::read::GzDecoder::new(file);
        let mut content = String::new();
        decoder.read_to_string(&mut content)?;

        assert!(content.contains("@read1"));
        assert!(content.contains("@read2"));
        assert!(content.contains("AAAA"));
        assert!(content.contains("CCCC"));

        Ok(())
    }

    #[test]
    fn test_quality_encoding() -> Result<()> {
        let temp = NamedTempFile::new()?;
        let path = temp.path();

        {
            let mut writer = FastqWriter::new(path)?;
            writer.write_record("test", b"ACGT", &[0, 10, 30, 41])?;
            writer.finish()?;
        }

        let file = File::open(path)?;
        let mut decoder = flate2::read::GzDecoder::new(file);
        let mut content = String::new();
        decoder.read_to_string(&mut content)?;

        let lines: Vec<&str> = content.lines().collect();
        let qual_line = lines[3];

        // Q0 + 33 = '!', Q10 + 33 = '+', Q30 + 33 = '?', Q41 + 33 = 'J'
        assert_eq!(qual_line, "!+?J");

        Ok(())
    }

    #[test]
    fn test_multi_threaded_writer() -> Result<()> {
        let temp = NamedTempFile::new()?;
        let path = temp.path();

        {
            let mut writer = FastqWriter::with_threads(path, 4)?;
            writer.write_record("read1", b"ACGT", &[30, 30, 30, 30])?;
            writer.write_record("read2", b"TGCA", &[35, 35, 35, 35])?;
            writer.finish()?;
        }

        let file = File::open(path)?;
        let mut decoder = flate2::read::GzDecoder::new(file);
        let mut content = String::new();
        decoder.read_to_string(&mut content)?;

        assert!(content.contains("@read1"));
        assert!(content.contains("@read2"));
        assert!(content.contains("ACGT"));
        assert!(content.contains("TGCA"));

        Ok(())
    }

    #[test]
    fn test_multi_threaded_large_output() -> Result<()> {
        let temp = NamedTempFile::new()?;
        let path = temp.path();

        {
            let mut writer = FastqWriter::with_threads(path, 4)?;
            // Write enough records to trigger multiple compression blocks
            for i in 0..10000 {
                let name = format!("read{i:05}");
                let seq = b"ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGT";
                let qual = vec![30u8; seq.len()];
                writer.write_record(&name, seq, &qual)?;
            }
            writer.finish()?;
        }

        let file = File::open(path)?;
        // Use MultiGzDecoder since parallel gzip produces concatenated streams
        let mut decoder = flate2::read::MultiGzDecoder::new(file);
        let mut content = String::new();
        decoder.read_to_string(&mut content)?;

        // Verify first and last records
        assert!(content.contains("@read00000"));
        assert!(content.contains("@read09999"));

        Ok(())
    }
}
