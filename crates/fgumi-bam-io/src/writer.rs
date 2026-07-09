//! BAM writer factories and related types.
//!
//! This module provides writer factories for BAM output, including standard BAM writers,
//! raw-bytes BAM writers, and BAI sidecar helpers (`write_bai_index` for a prebuilt index,
//! `write_bai_sidecar` to index a finished BAM and write `<bam>.bai` in one call, and
//! `bai_sidecar_path` for the path convention). Incremental indexing during write was
//! removed in Phase 4 of issue #330; BAI generation now lives in the typed-step
//! `IndexBamFinalizeHook` as a post-pipeline pass that delegates to `write_bai_sidecar`.

use anyhow::{Context, Result};
use bgzf::CompressionLevel;
use noodles::bam::bai;
use noodles::sam::Header;
use noodles_bgzf::io::Writer as BgzfWriter;
use std::fs::File;
use std::io::{self, Write};
use std::num::NonZero;
use std::path::{Path, PathBuf};

use crate::paths::is_stdout_path;
use crate::vendored::{MultithreadedWriter, MultithreadedWriterBuilder};

/// Opaque wrapper over either a single- or multi-threaded BGZF writer.
///
/// Uses `Box<dyn Write + Send>` to support both file and stdout output. The
/// internal representation is intentionally private so that the choice of
/// underlying multi-threaded BGZF writer is not part of the public API.
pub struct BgzfWriterEnum {
    inner: BgzfWriterImpl,
}

enum BgzfWriterImpl {
    SingleThreaded(BgzfWriter<Box<dyn Write + Send>>),
    MultiThreaded(MultithreadedWriter<Box<dyn Write + Send>>),
}

impl BgzfWriterEnum {
    pub(crate) fn single_threaded(writer: BgzfWriter<Box<dyn Write + Send>>) -> Self {
        Self { inner: BgzfWriterImpl::SingleThreaded(writer) }
    }

    pub(crate) fn multi_threaded(writer: MultithreadedWriter<Box<dyn Write + Send>>) -> Self {
        Self { inner: BgzfWriterImpl::MultiThreaded(writer) }
    }
}

impl Write for BgzfWriterEnum {
    fn write(&mut self, buf: &[u8]) -> io::Result<usize> {
        match &mut self.inner {
            BgzfWriterImpl::SingleThreaded(w) => w.write(buf),
            BgzfWriterImpl::MultiThreaded(w) => w.write(buf),
        }
    }

    fn flush(&mut self) -> io::Result<()> {
        match &mut self.inner {
            BgzfWriterImpl::SingleThreaded(w) => w.flush(),
            BgzfWriterImpl::MultiThreaded(w) => w.flush(),
        }
    }
}

impl BgzfWriterEnum {
    /// Finish writing and close the writer properly.
    /// This is especially important for multi-threaded writers to ensure
    /// all data is flushed and the EOF marker is written.
    ///
    /// # Errors
    /// Returns an error if flushing or finalizing the writer fails.
    pub fn finish(self) -> io::Result<()> {
        match self.inner {
            BgzfWriterImpl::SingleThreaded(mut w) => {
                w.flush()?;
                // Single-threaded writer writes EOF on drop
                Ok(())
            }
            BgzfWriterImpl::MultiThreaded(mut w) => {
                w.finish().map_err(|e| io::Error::other(e.to_string()))?;
                Ok(())
            }
        }
    }
}

/// Type alias for a BAM writer that supports both single and multi-threaded BGZF.
pub type BamWriter = noodles::bam::io::Writer<BgzfWriterEnum>;

/// Create a [`BgzfWriterEnum`] from a writer, selecting single- or multi-threaded based on
/// `threads`.
#[allow(clippy::cast_possible_truncation)]
pub(crate) fn make_bgzf_writer(
    output: Box<dyn Write + Send>,
    threads: usize,
    compression_level: u32,
) -> BgzfWriterEnum {
    if threads > 1 {
        let worker_count = NonZero::new(threads).expect("threads > 1 checked above");
        let mut builder = MultithreadedWriterBuilder::default().set_worker_count(worker_count);

        if let Ok(cl) = CompressionLevel::new(compression_level as u8) {
            builder = builder.set_compression_level(cl);
        }

        BgzfWriterEnum::multi_threaded(builder.build_from_writer(output))
    } else {
        let level = noodles_bgzf::io::writer::CompressionLevel::new(compression_level as u8)
            .unwrap_or_default();
        let writer = noodles_bgzf::io::writer::Builder::default()
            .set_compression_level(level)
            .build_from_writer(output);
        BgzfWriterEnum::single_threaded(writer)
    }
}

/// Create a BGZF-compressed writer over an arbitrary `output`, selecting a
/// multi-threaded encoder when `threads > 1`.
///
/// This exposes the same BGZF backend fgumi uses for BAM output so other output
/// formats (e.g. gzipped FASTQ) get identical, block-gzip (BGZF) framing and
/// multi-threaded compression rather than a separate implementation.
///
/// `compression_level` is the BGZF/DEFLATE level (0–12, where 0 means
/// uncompressed); values outside the valid range fall back to the backend
/// default.
#[must_use]
pub fn create_bgzf_writer(
    output: Box<dyn Write + Send>,
    threads: usize,
    compression_level: u32,
) -> BgzfWriterEnum {
    make_bgzf_writer(output, threads, compression_level)
}

/// Write a BAM header (magic, SAM text, reference sequences) to any writer.
///
/// Shared implementation used by the BAM writer factories (e.g. [`RawBamWriter`]).
///
/// Note: we use a manual implementation rather than delegating to
/// `noodles::bam::io::Writer::write_header` because the noodles encoder
/// produces subtly different output that causes read-back failures with
/// some writer backends (e.g. `MultithreadedWriter`).
///
/// # Errors
///
/// Returns an [`io::Error`] if writing to `writer` fails or if noodles
/// cannot serialize the SAM header text.
#[allow(clippy::cast_possible_truncation, clippy::cast_possible_wrap)]
pub fn write_bam_header(writer: &mut impl Write, header: &Header) -> io::Result<()> {
    // BAM magic
    writer.write_all(fgumi_raw_bam::BAM_MAGIC)?;

    // Header text (SAM header serialized using noodles)
    let mut sam_writer = noodles::sam::io::Writer::new(Vec::new());
    sam_writer.write_header(header)?;
    let header_bytes = sam_writer.into_inner();
    let l_text = header_bytes.len() as i32;
    writer.write_all(&l_text.to_le_bytes())?;
    writer.write_all(&header_bytes)?;

    // Reference sequences
    let n_ref = header.reference_sequences().len() as i32;
    writer.write_all(&n_ref.to_le_bytes())?;

    for (name, map) in header.reference_sequences() {
        // l_name: length of name + null terminator
        let l_name = (name.len() + 1) as u32;
        writer.write_all(&l_name.to_le_bytes())?;
        writer.write_all(name)?;
        writer.write_all(&[0u8])?; // null terminator

        // l_ref: reference length
        let l_ref = map.length().get() as i32;
        writer.write_all(&l_ref.to_le_bytes())?;
    }

    Ok(())
}

/// Raw BAM writer for writing raw record bytes directly.
///
/// This bypasses the noodles Record encoding path for maximum performance
/// when you already have raw BAM bytes. Writes records as:
/// - 4-byte `block_size` (little-endian)
/// - raw BAM record bytes
pub struct RawBamWriter {
    inner: BgzfWriterEnum,
}

impl RawBamWriter {
    /// Create a new raw BAM writer from a BGZF writer.
    #[must_use]
    pub fn new(inner: BgzfWriterEnum) -> Self {
        Self { inner }
    }

    /// Write the BAM header.
    ///
    /// # Errors
    /// Returns an error if writing to the underlying writer fails.
    pub fn write_header(&mut self, header: &Header) -> io::Result<()> {
        write_bam_header(&mut self.inner, header)
    }

    /// Write a raw BAM record.
    ///
    /// The bytes should be the raw BAM record data (without the 4-byte `block_size` prefix).
    ///
    /// # Errors
    /// Returns an error if writing to the underlying writer fails.
    #[inline]
    #[allow(clippy::cast_possible_truncation)]
    pub fn write_raw_record(&mut self, record_bytes: &[u8]) -> io::Result<()> {
        use std::io::Write;

        // Write block_size (4 bytes, little-endian)
        let block_size = record_bytes.len() as u32;
        self.inner.write_all(&block_size.to_le_bytes())?;

        // Write raw record bytes
        self.inner.write_all(record_bytes)
    }

    /// Write pre-formatted BAM record bytes directly to the BGZF stream.
    ///
    /// Unlike `write_raw_record`, this does NOT add a length prefix — the data
    /// must already contain properly formatted BAM records (4-byte LE size + record bytes).
    /// Used by the pipeline for bulk secondary output.
    ///
    /// # Errors
    /// Returns an error if writing to the underlying writer fails.
    pub fn write_raw_bytes(&mut self, data: &[u8]) -> io::Result<()> {
        use std::io::Write;
        self.inner.write_all(data)
    }

    /// Finish writing and close the writer.
    ///
    /// # Errors
    /// Returns an error if finalizing the writer fails.
    pub fn finish(self) -> io::Result<()> {
        self.inner.finish()
    }
}

/// Create a raw BAM writer for writing raw record bytes directly.
///
/// This is more efficient than `create_bam_writer` when you already have raw BAM bytes
/// because it bypasses the noodles Record encoding path.
///
/// # Errors
/// Returns an error if the file cannot be created or the header cannot be written.
///
/// # Panics
/// Panics if `threads > 1` but `NonZero::new` fails (should not happen).
pub fn create_raw_bam_writer<P: AsRef<Path>>(
    path: P,
    header: &Header,
    threads: usize,
    compression_level: u32,
) -> Result<RawBamWriter> {
    let path_ref = path.as_ref();
    let output = open_output_writer(path_ref)?;

    let bgzf_writer = make_bgzf_writer(output, threads, compression_level);

    let mut writer = RawBamWriter::new(bgzf_writer);
    writer
        .write_header(header)
        .with_context(|| format!("Failed to write header to: {}", path_ref.display()))?;
    Ok(writer)
}

/// Write a BAI index to a file.
///
/// # Arguments
/// * `path` - Path for the output BAI file (typically `output.bam.bai`)
/// * `index` - The BAI index to write
///
/// # Errors
/// Returns an error if the file cannot be created or writing the index fails.
pub fn write_bai_index<P: AsRef<Path>>(path: P, index: &bai::Index) -> Result<()> {
    let path_ref = path.as_ref();
    let file = File::create(path_ref)
        .with_context(|| format!("Failed to create index file: {}", path_ref.display()))?;
    let mut writer = bai::io::Writer::new(file);
    writer
        .write_index(index)
        .with_context(|| format!("Failed to write index to: {}", path_ref.display()))?;
    Ok(())
}

/// Append `.bai` to a BAM path to form its conventional sidecar index path.
///
/// This appends to the **full** path (samtools convention) rather than
/// replacing the final extension:
///
/// - `foo.bam` → `foo.bam.bai`
/// - `foo` (no extension) → `foo.bai`
/// - `foo.sorted` → `foo.sorted.bai`
///
/// Using `Path::with_extension("bam.bai")` instead would replace whatever
/// follows the last `.`, producing the correct result only for paths that end
/// in exactly `.bam` and mis-naming every other path (`foo` → `foo.bam.bai`,
/// `foo.sorted` → `foo.bam.bai`).
#[must_use]
pub fn bai_sidecar_path<P: AsRef<Path>>(bam_path: P) -> PathBuf {
    let mut index_os = bam_path.as_ref().as_os_str().to_owned();
    index_os.push(".bai");
    PathBuf::from(index_os)
}

/// Build a BAI index for a finished coordinate-sorted BAM and write it to the
/// conventional `<bam_path>.bai` sidecar, returning the path written.
///
/// The sidecar path is derived by [`bai_sidecar_path`] (append `.bai` to the
/// full path), so BAMs whose path does not end in `.bam` still get a correctly
/// named sidecar.
///
/// **Invariant:** the BAM at `bam_path` must be writer-closed before this is
/// called — [`noodles::bam::fs::index`] re-reads the file from disk to compute
/// virtual offsets, so any buffered writes must already be flushed.
///
/// # Errors
/// Returns an error if indexing the BAM or writing the BAI sidecar fails.
pub fn write_bai_sidecar<P: AsRef<Path>>(bam_path: P) -> Result<PathBuf> {
    let bam_path = bam_path.as_ref();

    let index = noodles::bam::fs::index(bam_path)
        .with_context(|| format!("Failed to index {}", bam_path.display()))?;

    let index_path = bai_sidecar_path(bam_path);
    write_bai_index(&index_path, &index)
        .with_context(|| format!("Failed to write BAI to {}", index_path.display()))?;

    Ok(index_path)
}

/// Create a BAM writer and write the header in one operation.
///
/// # Arguments
/// * `path` - Path for the output BAM file (use `-` or `/dev/stdout` for stdout)
/// * `header` - SAM header to write
/// * `threads` - Number of threads for BGZF compression (1 = single-threaded)
/// * `compression_level` - Compression level (0-12; 0 = uncompressed)
///
/// # Returns
/// A BAM writer with the header already written.
///
/// # Errors
/// Returns an error if the file cannot be created or the header cannot be written.
///
/// # Example
/// ```no_run
/// use fgumi_bam_io::writer::create_bam_writer;
/// use noodles::sam::Header;
/// use std::path::Path;
///
/// let header = Header::default();
/// // Single-threaded with fast compression
/// let mut writer = create_bam_writer(Path::new("output.bam"), &header, 1, 1).unwrap();
/// // Multi-threaded with 4 compression threads
/// let mut writer = create_bam_writer(Path::new("output.bam"), &header, 4, 1).unwrap();
/// ```
pub fn create_bam_writer<P: AsRef<Path>>(
    path: P,
    header: &Header,
    threads: usize,
    compression_level: u32,
) -> Result<BamWriter> {
    let path_ref = path.as_ref();
    let output = open_output_writer(path_ref)?;

    let mut bgzf_writer = make_bgzf_writer(output, threads, compression_level);

    // Write the header with the libdeflate-backed helper rather than the
    // noodles encoder: noodles' `Writer::write_header` corrupts the leading
    // BAM bytes when the underlying writer is `MultithreadedWriter`.
    write_bam_header(&mut bgzf_writer, header)
        .with_context(|| format!("Failed to write header to: {}", path_ref.display()))?;
    Ok(noodles::bam::io::Writer::from(bgzf_writer))
}

/// Create an optional BAM writer for rejects or filtered reads.
///
/// This is a convenience function for commands that optionally output rejected/filtered reads.
/// If the path is `None`, returns `None`. If the path is `Some`, creates a writer and writes
/// the header.
///
/// # Arguments
/// * `path` - Optional path for the output BAM file
/// * `header` - SAM header to write
/// * `threads` - Number of threads for BGZF compression (1 = single-threaded)
/// * `compression_level` - Compression level (0-12; 0 = uncompressed)
///
/// # Returns
/// An optional BAM writer, or None if path is None.
///
/// # Errors
/// Returns an error if the file cannot be created or the header cannot be written.
///
/// # Example
/// ```no_run
/// use fgumi_bam_io::writer::create_optional_bam_writer;
/// use noodles::sam::Header;
/// use std::path::PathBuf;
///
/// let header = Header::default();
/// let rejects = Some(PathBuf::from("rejects.bam"));
/// let mut writer = create_optional_bam_writer(rejects, &header, 1, 1).unwrap();
/// ```
pub fn create_optional_bam_writer<P: AsRef<Path>>(
    path: Option<P>,
    header: &Header,
    threads: usize,
    compression_level: u32,
) -> Result<Option<BamWriter>> {
    match path {
        Some(p) => Ok(Some(create_bam_writer(p, header, threads, compression_level)?)),
        None => Ok(None),
    }
}

/// Open an output writer for a path, supporting stdout via `-` or `/dev/stdout`.
pub(crate) fn open_output_writer<P: AsRef<Path>>(path: P) -> Result<Box<dyn Write + Send>> {
    let path_ref = path.as_ref();
    if is_stdout_path(path_ref) {
        Ok(Box::new(std::io::stdout()))
    } else {
        let file = File::create(path_ref)
            .with_context(|| format!("Failed to create output BAM: {}", path_ref.display()))?;
        Ok(Box::new(file))
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use noodles::sam::header::record::value::{Map, map::ReferenceSequence};
    use std::num::NonZeroUsize;
    use tempfile::NamedTempFile;

    fn create_test_header() -> Header {
        let mut builder = Header::builder();
        let ref_seq = Map::<ReferenceSequence>::new(
            NonZeroUsize::new(100).expect("100 is non-zero constant"),
        );
        builder = builder.add_reference_sequence(b"chr1", ref_seq);
        builder.build()
    }

    #[test]
    fn bai_sidecar_path_appends_dot_bai_to_full_path() {
        // The samtools convention is to append `.bai` to the whole BAM path,
        // not to replace the final extension. The previous
        // `with_extension("bam.bai")` only produced the right answer for paths
        // ending in exactly `.bam`; these cases pin the fix (S3-006 / S5b2-008).
        assert_eq!(bai_sidecar_path("foo.bam"), PathBuf::from("foo.bam.bai"));
        // Extensionless path: must become `foo.bai`, not `foo.bam.bai`.
        assert_eq!(bai_sidecar_path("foo"), PathBuf::from("foo.bai"));
        // Non-`.bam` extension: must keep the original name and append `.bai`.
        assert_eq!(bai_sidecar_path("foo.sorted"), PathBuf::from("foo.sorted.bai"));
        // Path with directory components is preserved.
        assert_eq!(bai_sidecar_path("/tmp/run1/out.bam"), PathBuf::from("/tmp/run1/out.bam.bai"));
    }

    #[test]
    fn test_create_bam_writer() -> Result<()> {
        let temp_file = NamedTempFile::new()?;
        let header = create_test_header();

        let writer = create_bam_writer(temp_file.path(), &header, 1, 6);
        assert!(writer.is_ok());

        Ok(())
    }

    #[test]
    fn test_create_bam_writer_multithreaded() -> Result<()> {
        let temp_file = NamedTempFile::new()?;
        let header = create_test_header();

        let writer = create_bam_writer(temp_file.path(), &header, 4, 6);
        assert!(writer.is_ok());

        Ok(())
    }

    #[test]
    fn test_create_bam_writer_invalid_path() {
        let header = create_test_header();
        let result = create_bam_writer("/invalid/path/output.bam", &header, 1, 6);
        assert!(result.is_err());
        if let Err(e) = result {
            let err_msg = e.to_string();
            assert!(err_msg.contains("Failed to create output BAM"));
        }
    }

    #[test]
    fn test_create_optional_bam_writer_some() -> Result<()> {
        let temp_file = NamedTempFile::new()?;
        let header = create_test_header();

        let writer = create_optional_bam_writer(Some(temp_file.path()), &header, 1, 6)?;
        assert!(writer.is_some());

        Ok(())
    }

    #[test]
    fn test_create_optional_bam_writer_none() -> Result<()> {
        let header = create_test_header();

        let writer = create_optional_bam_writer(None::<&str>, &header, 1, 6)?;
        assert!(writer.is_none());

        Ok(())
    }

    #[test]
    fn test_create_optional_bam_writer_invalid_path() {
        let header = create_test_header();
        let result = create_optional_bam_writer(Some("/invalid/path/output.bam"), &header, 1, 6);
        assert!(result.is_err());
    }

    #[test]
    fn test_bgzf_writer_flush_single_threaded() -> Result<()> {
        let temp_file = NamedTempFile::new()?;
        let output_file: Box<dyn Write + Send> = Box::new(File::create(temp_file.path())?);
        let mut writer = BgzfWriterEnum::single_threaded(BgzfWriter::new(output_file));

        // Write some data and flush
        writer.write_all(b"test data")?;
        let result = writer.flush();
        assert!(result.is_ok());

        Ok(())
    }

    #[test]
    fn test_bgzf_writer_flush_multithreaded() -> Result<()> {
        let temp_file = NamedTempFile::new()?;
        let output_file: Box<dyn Write + Send> = Box::new(File::create(temp_file.path())?);
        let worker_count = NonZero::new(2).expect("2 is non-zero");
        let compression_level = CompressionLevel::new(6).expect("valid compression level");
        let mut writer = BgzfWriterEnum::multi_threaded(MultithreadedWriter::with_worker_count(
            worker_count,
            output_file,
            compression_level,
        ));

        // Write some data and flush
        writer.write_all(b"test data")?;
        let result = writer.flush();
        assert!(result.is_ok());

        Ok(())
    }

    #[test]
    fn test_bgzf_writer_finish_single_threaded() -> Result<()> {
        let temp_file = NamedTempFile::new()?;
        let output_file: Box<dyn Write + Send> = Box::new(File::create(temp_file.path())?);
        let mut writer = BgzfWriterEnum::single_threaded(BgzfWriter::new(output_file));

        // Write some data
        writer.write_all(b"test data")?;

        // Finish the writer - this should flush and close properly
        let result = writer.finish();
        assert!(result.is_ok());

        Ok(())
    }

    #[test]
    fn test_bgzf_writer_finish_multithreaded() -> Result<()> {
        let temp_file = NamedTempFile::new()?;
        let output_file: Box<dyn Write + Send> = Box::new(File::create(temp_file.path())?);
        let worker_count = NonZero::new(2).expect("2 is non-zero");
        let compression_level = CompressionLevel::new(6).expect("valid compression level");
        let mut writer = BgzfWriterEnum::multi_threaded(MultithreadedWriter::with_worker_count(
            worker_count,
            output_file,
            compression_level,
        ));

        // Write some data
        writer.write_all(b"test data")?;

        // Finish the writer - this ensures proper finalization
        let result = writer.finish();
        assert!(result.is_ok());

        Ok(())
    }

    #[test]
    fn test_bgzf_writer_write_single_threaded() -> Result<()> {
        let temp_file = NamedTempFile::new()?;
        let output_file: Box<dyn Write + Send> = Box::new(File::create(temp_file.path())?);
        let mut writer = BgzfWriterEnum::single_threaded(BgzfWriter::new(output_file));

        // Test writing via the Write trait
        let bytes_written = writer.write(b"test")?;
        assert_eq!(bytes_written, 4);

        Ok(())
    }

    #[test]
    fn test_bgzf_writer_write_multithreaded() -> Result<()> {
        let temp_file = NamedTempFile::new()?;
        let output_file: Box<dyn Write + Send> = Box::new(File::create(temp_file.path())?);
        let worker_count = NonZero::new(2).expect("2 is non-zero");
        let compression_level = CompressionLevel::new(6).expect("valid compression level");
        let mut writer = BgzfWriterEnum::multi_threaded(MultithreadedWriter::with_worker_count(
            worker_count,
            output_file,
            compression_level,
        ));

        // Test writing via the Write trait
        let bytes_written = writer.write(b"test")?;
        assert_eq!(bytes_written, 4);

        Ok(())
    }

    #[test]
    fn test_is_stdout_path() {
        use crate::paths::is_stdout_path;
        use std::path::Path;
        assert!(is_stdout_path("-"));
        assert!(is_stdout_path("/dev/stdout"));
        assert!(is_stdout_path(Path::new("-")));
        assert!(is_stdout_path(Path::new("/dev/stdout")));

        assert!(!is_stdout_path("output.bam"));
        assert!(!is_stdout_path("/path/to/file.bam"));
        assert!(!is_stdout_path(""));
        assert!(!is_stdout_path("/dev/null"));
    }

    #[test]
    fn test_write_raw_bytes() -> Result<()> {
        let temp = tempfile::NamedTempFile::new()?;
        let header = noodles::sam::Header::default();
        let mut writer = create_raw_bam_writer(temp.path(), &header, 1, 1)?;

        let record_bytes = vec![1, 2, 3, 4, 5];
        let mut formatted = Vec::new();
        #[allow(clippy::cast_possible_truncation)]
        let len = record_bytes.len() as u32;
        formatted.extend_from_slice(&len.to_le_bytes());
        formatted.extend_from_slice(&record_bytes);

        writer.write_raw_bytes(&formatted)?;
        writer.finish()?;

        assert!(temp.path().metadata()?.len() > 0);
        Ok(())
    }
}
