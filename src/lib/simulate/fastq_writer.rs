//! Gzipped FASTQ file writer for simulation output.
//!
//! Provides a simple interface for writing FASTQ records to gzip-compressed files.
//! Supports both single-threaded and multi-threaded parallel compression.

use super::parallel_gzip_writer::{ParallelGzipConfig, ParallelGzipWriter};
use anyhow::{Context, Result};
use fgumi_bam_io::{OutputFile, OutputSink, close_buffered};
use flate2::Compression;
use flate2::write::GzEncoder;
use std::io::{BufWriter, Write};
use std::path::{Path, PathBuf};

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
    /// The output path, named in close errors.
    path: PathBuf,
}

enum FastqWriterInner {
    SingleThreaded(GzEncoder<BufWriter<Box<dyn OutputSink>>>),
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
        Self::with_threads(path, 1)
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
        let file = OutputFile::create(path)
            .with_context(|| format!("Failed to create {}", path.display()))?;
        Self::from_sink(Box::new(file), path, threads)
    }

    /// Build the writer over an already-open sink; `path` names it in errors.
    fn from_sink(sink: Box<dyn OutputSink>, path: &Path, threads: usize) -> Result<Self> {
        let inner = if threads <= 1 {
            let buf = BufWriter::new(sink);
            FastqWriterInner::SingleThreaded(GzEncoder::new(buf, Compression::default()))
        } else {
            let config = ParallelGzipConfig::with_threads(threads);
            let writer = ParallelGzipWriter::new(sink, &config)
                .with_context(|| "Failed to create parallel gzip writer")?;
            FastqWriterInner::MultiThreaded(writer)
        };
        Ok(Self { inner, path: path.to_path_buf() })
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

    /// Finish writing, flush all data, and sync and close the file.
    ///
    /// This must be called to ensure all data is written and the gzip stream
    /// is properly terminated.
    ///
    /// # Errors
    ///
    /// Returns an error if flushing, syncing, or closing the file fails.
    pub fn finish(self) -> Result<()> {
        let path = self.path.display();
        match self.inner {
            FastqWriterInner::SingleThreaded(writer) => {
                let buf = writer
                    .finish()
                    .with_context(|| format!("Failed to finish gzip stream {path}"))?;
                close_buffered(buf).with_context(|| format!("Failed to sync/close {path}"))?;
            }
            FastqWriterInner::MultiThreaded(writer) => {
                writer
                    .finish()
                    .with_context(|| format!("Failed to finish parallel gzip stream {path}"))?;
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
    use std::fs::File;
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

    /// Both arms close the file once, and a close failure fails `finish`
    /// naming the output; the bytes written are a complete gzip stream.
    #[rstest::rstest]
    #[case::single_threaded(1)]
    #[case::multi_threaded(2)]
    fn test_finish_closes_once_and_reports_close_error(#[case] threads: usize) -> Result<()> {
        use crate::simulate::test_support::CloseProbe;

        let path = Path::new("reads.fq.gz");
        let probe = CloseProbe::default();
        let mut writer = FastqWriter::from_sink(Box::new(probe.clone()), path, threads)?;
        writer.write_record("r1", b"ACGT", &[30; 4])?;
        writer.finish()?;
        assert_eq!(probe.closes(), 1);
        assert_eq!(probe.decoded(), b"@r1\nACGT\n+\n????\n");

        let probe = CloseProbe::failing_close();
        let mut writer = FastqWriter::from_sink(Box::new(probe.clone()), path, threads)?;
        writer.write_record("r1", b"ACGT", &[30; 4])?;
        let err = writer.finish().expect_err("a failed close must surface");
        let msg = format!("{err:#}");
        assert!(msg.contains(crate::simulate::test_support::CLOSE_ERROR), "{msg}");
        assert!(msg.contains("reads.fq.gz"), "the error must name the output: {msg}");
        Ok(())
    }
}
