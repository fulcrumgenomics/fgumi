//! Shared raw alignment output, selected by the user-visible destination name.

use std::io;
use std::path::Path;

use anyhow::{Context, Result};
use noodles::sam::Header;

use crate::writer::RawBamWriter;

/// A raw alignment writer with its header already written.
///
/// BAM uses the requested BGZF worker count. CRAM deliberately packs, encodes,
/// and writes containers inline, without creating an htslib thread pool. This
/// serial fallback does not build an index. Always call [`Self::finish`] to
/// check finalization errors; dropping CRAM output safely abandons the stream.
pub struct RawAlignmentWriter {
    inner: AlignmentWriterBackend,
}

enum AlignmentWriterBackend {
    Bam(RawBamWriter),
    #[cfg(feature = "cram")]
    Cram(fgumi_htslib::CramWriter),
}

impl RawAlignmentWriter {
    /// Write a BAM record body, without its four-byte length prefix.
    ///
    /// # Errors
    /// Returns an error if writing fails, or CRAM rejects the record or fails
    /// to encode a completed container.
    pub fn write_raw_record(&mut self, record: &[u8]) -> io::Result<()> {
        match &mut self.inner {
            AlignmentWriterBackend::Bam(writer) => writer.write_raw_record(record),
            #[cfg(feature = "cram")]
            AlignmentWriterBackend::Cram(writer) => {
                writer.write_raw(record)?;
                for job in writer.drain_jobs() {
                    job.encode()?.write()?;
                }
                Ok(())
            }
        }
    }

    /// Write complete framed BAM records (four-byte little-endian size + body).
    ///
    /// No partial frame is retained between calls. Records preceding an invalid
    /// frame may already have been written when an error is returned.
    ///
    /// # Errors
    /// Returns an error for a truncated prefix or body, a body shorter than the
    /// BAM core, or any record-write or container-encode failure.
    pub fn write_raw_bytes(&mut self, mut data: &[u8]) -> io::Result<()> {
        while !data.is_empty() {
            let (prefix, remaining) = data.split_first_chunk::<4>().ok_or_else(|| {
                io::Error::new(io::ErrorKind::UnexpectedEof, "truncated BAM record length prefix")
            })?;
            let size = usize::try_from(u32::from_le_bytes(*prefix)).map_err(|_| {
                io::Error::new(
                    io::ErrorKind::InvalidData,
                    "BAM record length exceeds address space",
                )
            })?;
            if size < 32 {
                return Err(io::Error::new(
                    io::ErrorKind::InvalidData,
                    "BAM record length is shorter than the 32-byte core",
                ));
            }
            let (record, rest) = remaining.split_at_checked(size).ok_or_else(|| {
                io::Error::new(io::ErrorKind::UnexpectedEof, "truncated BAM record body")
            })?;
            self.write_raw_record(record)?;
            data = rest;
        }
        Ok(())
    }

    /// Finish encoding, write the EOF marker, and close the output.
    ///
    /// # Errors
    /// Returns an error if final encoding, writing, or closing fails.
    pub fn finish(self) -> io::Result<()> {
        match self.inner {
            AlignmentWriterBackend::Bam(writer) => writer.finish(),
            #[cfg(feature = "cram")]
            AlignmentWriterBackend::Cram(writer) => {
                let (jobs, finisher) = writer.seal()?;
                for job in jobs {
                    job.encode()?.write()?;
                }
                finisher.finish()
            }
        }
    }
}

/// Create an alignment writer and write its header.
///
/// `path` is the actual write path (possibly a staged temporary), while `named`
/// is the user-visible destination: its `.cram` extension selects CRAM, its
/// basename becomes the CRAM file ID, and errors identify this name. Every other
/// name selects BAM, including stdout (`-` or `/dev/stdout`). `threads` and
/// `compression_level` apply only to BAM; CRAM uses serial encoding and default
/// codecs. CRAM references follow the `cram_input` module's configuration when the
/// `cram` feature is enabled; otherwise CRAM names fail before opening `path`.
/// Neither format produces an index.
///
/// # Errors
/// Returns an error if CRAM is unavailable, the output cannot be opened, the
/// reference cannot be used, or the header cannot be written.
pub fn create_raw_alignment_writer(
    path: &Path,
    named: &Path,
    header: &Header,
    threads: usize,
    compression_level: u32,
) -> Result<RawAlignmentWriter> {
    validate_alignment_output(named)?;
    if crate::paths::is_cram_path(named) {
        #[cfg(feature = "cram")]
        {
            return create_serial_cram_writer(path, named, header, compression_level)
                .map(|writer| RawAlignmentWriter { inner: AlignmentWriterBackend::Cram(writer) })
                .with_context(|| format!("Failed to create CRAM output: {}", named.display()));
        }
    }
    let output = crate::output::open_output_sink(path)
        .with_context(|| format!("Failed to create output BAM: {}", named.display()))?;
    let mut writer =
        RawBamWriter::new(crate::writer::make_bgzf_writer(output, threads, compression_level));
    writer
        .write_header(header)
        .with_context(|| format!("Failed to write header to: {}", named.display()))?;
    Ok(RawAlignmentWriter { inner: AlignmentWriterBackend::Bam(writer) })
}

/// Validate an alignment destination without opening or changing it.
///
/// CRAM names (case-insensitive `.cram`) require the `cram` feature. Other
/// destinations, including stdout, are accepted. This validates format support,
/// not filesystem permissions or path existence.
///
/// # Errors
/// Returns an error for a CRAM destination when CRAM support is disabled.
pub fn validate_alignment_output(path: &Path) -> Result<()> {
    #[cfg(not(feature = "cram"))]
    crate::reject_cram_output(path)?;
    #[cfg(feature = "cram")]
    let _ = path;
    Ok(())
}

/// Resolved CRAM output policy shared by serial and pipeline writers.
///
/// The actual write path is deliberately separate: only the user-visible name
/// supplies the file ID and warning destination.
#[cfg(feature = "cram")]
pub struct CramOutputPolicy {
    /// Explicit reference, or `None` to leave lookup to htslib.
    pub reference: Option<std::path::PathBuf>,
    /// Store bases verbatim when neither an explicit reference nor lookup is configured.
    pub no_ref: bool,
    /// User-visible basename embedded in the CRAM file definition.
    pub file_id: Option<String>,
}

#[cfg(feature = "cram")]
impl CramOutputPolicy {
    /// Resolve reference selection and emit applicable output warnings.
    ///
    /// `may_be_mapped` includes headers deferred until an aligner runs. The
    /// default BAM compression level (1) is silently ignored for CRAM; a
    /// non-default request warns because CRAM uses its own default codecs.
    #[must_use]
    pub fn new(named: &Path, may_be_mapped: bool, compression_level: u32) -> Self {
        use crate::cram_input::{NO_REFERENCE_HINT, cram_reference, cram_reference_configured};

        let reference = cram_reference();
        let no_ref = !cram_reference_configured();
        if no_ref && may_be_mapped {
            log::warn!(
                "Writing {} as a reference-free CRAM ({}): bases are stored verbatim, so \
                 mapped reads compress far less than with a reference",
                named.display(),
                NO_REFERENCE_HINT
            );
        }
        if compression_level != 1 {
            log::warn!(
                "--compression-level applies to BGZF (BAM) output only; CRAM output {} uses \
                 htslib's CRAM codecs at their defaults",
                named.display()
            );
        }
        Self {
            reference,
            no_ref,
            file_id: named.file_name().map(|name| name.to_string_lossy().into_owned()),
        }
    }
}

/// Set up the serial executor with the shared CRAM reference and naming policy.
#[cfg(feature = "cram")]
fn create_serial_cram_writer(
    path: &Path,
    named: &Path,
    header: &Header,
    compression_level: u32,
) -> Result<fgumi_htslib::CramWriter> {
    let policy =
        CramOutputPolicy::new(named, !header.reference_sequences().is_empty(), compression_level);
    let options = fgumi_htslib::CramWriteOptions {
        reference: policy.reference.as_deref(),
        no_ref: policy.no_ref,
        file_id: policy.file_id.as_deref(),
        ..fgumi_htslib::CramWriteOptions::default()
    };
    let text = crate::writer::sam_header_text(header)?;
    Ok(fgumi_htslib::CramWriter::create(path, std::str::from_utf8(&text)?, &options)?)
}

#[cfg(test)]
mod tests {
    use super::*;
    use fgumi_raw_bam::{RawRecord, SamBuilder};
    use fgumi_tag::SamTag;
    use noodles::sam::header::record::value::{Map, map::ReferenceSequence};
    use rstest::rstest;
    use std::num::NonZeroUsize;

    fn header() -> Header {
        Header::builder()
            .add_reference_sequence(
                b"chr1",
                Map::<ReferenceSequence>::new(NonZeroUsize::new(1000).unwrap()),
            )
            .build()
    }

    fn record(unmapped: bool) -> RawRecord {
        indexed_record(unmapped, 0)
    }

    fn indexed_record(unmapped: bool, index: usize) -> RawRecord {
        let name = format!("read{index:05}");
        let rx = format!("ACGT{index:05}");
        let quality = u8::try_from(index % 40).unwrap();
        SamBuilder::new()
            .read_name(name.as_bytes())
            .ref_id(if unmapped { -1 } else { 0 })
            .pos(if unmapped { -1 } else { 10 })
            .flags(if unmapped { fgumi_raw_bam::flags::UNMAPPED } else { 0 })
            // Deliberately nonzero MAPQ/CIGAR on unmapped input: CRAM drops them.
            .mapq(42)
            .cigar_ops(&[4 << 4])
            .sequence(b"ACGT")
            .qualities(&[quality, quality + 1, quality + 2, quality + 3])
            .add_string_tag(SamTag::RX, rx.as_bytes())
            .build()
    }

    fn framed(record: &RawRecord) -> Vec<u8> {
        let mut bytes = u32::try_from(record.len()).unwrap().to_le_bytes().to_vec();
        bytes.extend_from_slice(record);
        bytes
    }

    fn read_records(path: &Path, cram: bool) -> Result<(Header, Vec<RawRecord>)> {
        #[cfg(feature = "cram")]
        if cram {
            let mut reader = fgumi_htslib::CramReader::open(path, None)?;
            let header = reader.header_text().parse()?;
            let mut records = Vec::new();
            while let Some(job) = reader.next_slice()? {
                let decoded = job.decode()?;
                let mut bytes = std::io::Cursor::new(decoded.bytes);
                let mut record = RawRecord::new();
                while fgumi_raw_bam::read_raw_record(&mut bytes, &mut record)? != 0 {
                    records.push(record.clone());
                }
            }
            return Ok((header, records));
        }
        assert!(!cram);
        let (mut reader, header) = crate::create_raw_bam_reader(path, 1)?;
        let mut records = Vec::new();
        let mut record = RawRecord::new();
        while reader.read_record(&mut record)? != 0 {
            records.push(record.clone());
        }
        Ok((header, records))
    }

    fn assert_records(decoded: &[RawRecord], expected: &[RawRecord], cram: bool, unmapped: bool) {
        assert_eq!(decoded.len(), expected.len());
        for (actual, record) in decoded.iter().zip(expected) {
            if !cram {
                assert_eq!(actual.as_ref(), record.as_ref());
                continue;
            }
            assert_eq!(actual.read_name(), record.read_name());
            assert_eq!(actual.flags(), record.flags());
            assert_eq!(actual.ref_id(), record.ref_id());
            assert_eq!(actual.pos(), record.pos());
            assert_eq!(actual.mate_ref_id(), record.mate_ref_id());
            assert_eq!(actual.mate_pos(), record.mate_pos());
            assert_eq!(actual.template_length(), record.template_length());
            assert_eq!(actual.sequence_vec(), record.sequence_vec());
            assert_eq!(actual.quality_scores(), record.quality_scores());
            assert_eq!(actual.tags().get(SamTag::RX), record.tags().get(SamTag::RX));
            assert_eq!(actual.mapq(), if unmapped { 0 } else { record.mapq() });
            assert_eq!(actual.cigar_to_string(), if unmapped { "" } else { "4M" });
            // Reference-based CRAM regenerates MD/NM; reference-free CRAM
            // need not. Compare the other tags independently of aux byte order.
            let mut preserved_tags: Vec<_> = actual
                .tags()
                .iter_typed()
                .filter(|(tag, _)| {
                    *tag != <[u8; 2]>::from(SamTag::MD) && *tag != <[u8; 2]>::from(SamTag::NM)
                })
                .collect();
            preserved_tags.sort_by_key(|(tag, _)| *tag);
            let mut expected_tags: Vec<_> = record.tags().iter_typed().collect();
            expected_tags.sort_by_key(|(tag, _)| *tag);
            assert_eq!(preserved_tags, expected_tags);
        }
    }

    #[rstest]
    #[case::bam_mapped("out.bam", false)]
    #[case::bam_unmapped("out.bam", true)]
    fn staged_bam_round_trip(
        #[case] name: &str,
        #[case] unmapped: bool,
        #[values(false, true)] bulk: bool,
        #[values(1, 2)] threads: usize,
    ) -> Result<()> {
        staged_round_trip(name, unmapped, bulk, threads)
    }

    #[cfg(feature = "cram")]
    #[rstest]
    #[case::cram_mapped("out.cram", false)]
    #[case::cram_unmapped("out.cram", true)]
    fn staged_cram_round_trip(
        #[case] name: &str,
        #[case] unmapped: bool,
        #[values(false, true)] bulk: bool,
        #[values(1, 2)] threads: usize,
    ) -> Result<()> {
        staged_round_trip(name, unmapped, bulk, threads)
    }

    fn staged_round_trip(name: &str, unmapped: bool, bulk: bool, threads: usize) -> Result<()> {
        let dir = tempfile::tempdir()?;
        let staged = dir.path().join(".staged.tmp");
        let named = dir.path().join(name);
        let header = header();
        let expected: Vec<_> = (0..3).map(|i| indexed_record(unmapped, i)).collect();
        let mut writer = create_raw_alignment_writer(&staged, &named, &header, threads, 1)?;
        if bulk {
            let bytes: Vec<_> = expected.iter().flat_map(framed).collect();
            writer.write_raw_bytes(&[])?;
            writer.write_raw_bytes(&bytes)?;
        } else {
            for record in &expected {
                writer.write_raw_record(record)?;
            }
        }
        writer.finish()?;
        let cram = crate::is_cram_path(Path::new(name));
        if cram {
            let bytes = std::fs::read(&staged)?;
            assert_eq!(&bytes[..4], b"CRAM");
            assert_eq!(&bytes[6..14], b"out.cram");
        }
        let (decoded_header, decoded) = read_records(&staged, cram)?;
        assert_eq!(decoded_header.reference_sequences(), header.reference_sequences());
        assert_records(&decoded, &expected, cram, unmapped);
        Ok(())
    }

    #[rstest]
    #[case::bam("out.bam")]
    #[cfg_attr(feature = "cram", case::cram("out.cram"))]
    fn empty_output(#[case] name: &str) -> Result<()> {
        let dir = tempfile::tempdir()?;
        let path = dir.path().join(".staged.tmp");
        create_raw_alignment_writer(&path, Path::new(name), &Header::default(), 1, 1)?.finish()?;
        let (header, records) = read_records(&path, crate::is_cram_path(Path::new(name)))?;
        assert_eq!(header, Header::default());
        assert!(records.is_empty());
        Ok(())
    }

    #[cfg_attr(feature = "cram", rstest(name => ["out.bam", "out.cram"]))]
    #[cfg_attr(not(feature = "cram"), rstest(name => ["out.bam"]))]
    #[case::short_prefix(vec![1, 2, 3], io::ErrorKind::UnexpectedEof)]
    #[case::short_body(40_u32.to_le_bytes().to_vec(), io::ErrorKind::UnexpectedEof)]
    #[case::huge_body(u32::MAX.to_le_bytes().to_vec(), io::ErrorKind::UnexpectedEof)]
    #[case::short_core(31_u32.to_le_bytes().to_vec(), io::ErrorKind::InvalidData)]
    #[case::zero_length(0_u32.to_le_bytes().to_vec(), io::ErrorKind::InvalidData)]
    fn malformed_frame(
        #[case] bytes: Vec<u8>,
        #[case] expected: io::ErrorKind,
        name: &str,
    ) -> Result<()> {
        let dir = tempfile::tempdir()?;
        let path = dir.path().join(".staged.tmp");
        let mut writer = create_raw_alignment_writer(&path, Path::new(name), &header(), 1, 1)?;
        let error = writer.write_raw_bytes(&bytes).unwrap_err();
        assert_eq!(error.kind(), expected);
        Ok(())
    }

    #[cfg(feature = "cram")]
    #[rstest]
    #[case::record_at_a_time(false)]
    #[case::framed_bytes(true)]
    fn cram_drains_completed_containers_and_finishes_partial_container(
        #[case] bulk: bool,
    ) -> Result<()> {
        let dir = tempfile::tempdir()?;
        let path = dir.path().join(".staged.tmp");
        let mut writer =
            create_raw_alignment_writer(&path, Path::new("out.cram"), &header(), 4, 1)?;
        let expected: Vec<_> = (0..25_001).map(|i| indexed_record(false, i)).collect();
        if bulk {
            writer.write_raw_bytes(&expected.iter().flat_map(framed).collect::<Vec<_>>())?;
        } else {
            for record in &expected {
                writer.write_raw_record(record)?;
            }
        }
        writer.finish()?;
        let (_, records) = read_records(&path, true)?;
        assert_records(&records, &expected, true, false);
        Ok(())
    }

    #[rstest]
    #[case::single_threaded(1)]
    #[case::multi_threaded(2)]
    fn bam_matches_existing_factory(#[case] threads: usize) -> Result<()> {
        let dir = tempfile::tempdir()?;
        // The user-visible name, not the staging suffix, selects the format.
        let staged = dir.path().join("staged.cram");
        let legacy = dir.path().join("legacy.bam");
        let header = header();
        let record = record(false);
        let mut writer =
            create_raw_alignment_writer(&staged, Path::new("out.bam"), &header, threads, 1)?;
        let mut original = crate::create_raw_bam_writer(&legacy, &header, threads, 1)?;
        writer.write_raw_record(&record)?;
        original.write_raw_record(&record)?;
        writer.finish()?;
        original.finish()?;
        assert_eq!(std::fs::read(staged)?, std::fs::read(legacy)?);
        Ok(())
    }

    #[rstest]
    #[case::bam("out.bam")]
    #[cfg_attr(feature = "cram", case::cram("out.cram"))]
    fn creation_error_names_destination(#[case] name: &str) -> Result<()> {
        let dir = tempfile::tempdir()?;
        let missing = dir.path().join("missing").join(".staged.tmp");
        let error = create_raw_alignment_writer(&missing, Path::new(name), &header(), 1, 1)
            .err()
            .expect("missing parent must fail");
        assert!(error.to_string().contains(name), "{error:#}");
        if !crate::is_cram_path(Path::new(name)) {
            assert!(error.to_string().starts_with("Failed to create output BAM:"));
        }
        Ok(())
    }

    #[cfg(feature = "cram")]
    #[test]
    fn dropping_cram_abandons_without_finalizing() -> Result<()> {
        let dir = tempfile::tempdir()?;
        let path = dir.path().join(".staged.tmp");
        let mut writer =
            create_raw_alignment_writer(&path, Path::new("out.cram"), &header(), 1, 1)?;
        writer.write_raw_record(&record(false))?;
        drop(writer);
        // Abandonment is safe but does not make incomplete output look complete.
        let mut reader = fgumi_htslib::CramReader::open(&path, None)?;
        assert!(reader.next_slice().is_err());
        Ok(())
    }

    #[cfg(feature = "cram")]
    #[rstest]
    #[case::short_prefix(vec![1, 2, 3], io::ErrorKind::UnexpectedEof)]
    #[case::short_body(40_u32.to_le_bytes().to_vec(), io::ErrorKind::UnexpectedEof)]
    #[case::short_core(31_u32.to_le_bytes().to_vec(), io::ErrorKind::InvalidData)]
    fn malformed_suffix_abandons_completed_containers_without_eof(
        #[case] suffix: Vec<u8>,
        #[case] error_kind: io::ErrorKind,
    ) -> Result<()> {
        let dir = tempfile::tempdir()?;
        let path = dir.path().join(".staged.tmp");
        let mut writer =
            create_raw_alignment_writer(&path, Path::new("out.cram"), &header(), 1, 1)?;
        let expected: Vec<_> = (0..25_001).map(|i| indexed_record(false, i)).collect();
        let mut bytes: Vec<_> = expected.iter().flat_map(framed).collect();
        bytes.extend_from_slice(&suffix);
        assert_eq!(writer.write_raw_bytes(&bytes).unwrap_err().kind(), error_kind);
        drop(writer);

        let mut reader = fgumi_htslib::CramReader::open(&path, None)?;
        let mut decoded = Vec::new();
        let error = loop {
            match reader.next_slice() {
                Ok(Some(job)) => {
                    let mut bytes = io::Cursor::new(job.decode()?.bytes);
                    let mut record = RawRecord::new();
                    while fgumi_raw_bam::read_raw_record(&mut bytes, &mut record)? != 0 {
                        decoded.push(record.clone());
                    }
                }
                Ok(None) => panic!("abandoned CRAM must not report clean EOF"),
                Err(error) => break error,
            }
        };
        assert!(error.to_string().contains("truncated"), "{error}");
        // Completed containers survive; the pending final container is abandoned.
        assert!(!decoded.is_empty());
        assert!(decoded.len() < expected.len());
        assert_records(&decoded, &expected[..decoded.len()], true, false);
        Ok(())
    }

    #[rstest]
    #[case::bam("out.bam", true)]
    #[case::stdout("-", true)]
    #[case::cram("out.cram", cfg!(feature = "cram"))]
    #[case::uppercase_cram("out.CRAM", cfg!(feature = "cram"))]
    fn output_validation_does_not_open_destination(
        #[case] name: &str,
        #[case] supported: bool,
    ) -> Result<()> {
        let dir = tempfile::tempdir()?;
        let path = dir.path().join(name);
        assert_eq!(crate::validate_alignment_output(&path).is_ok(), supported);
        assert!(!path.exists());
        std::fs::write(&path, b"preserve existing output")?;
        assert_eq!(crate::validate_alignment_output(&path).is_ok(), supported);
        assert_eq!(std::fs::read(path)?, b"preserve existing output");
        Ok(())
    }

    #[cfg(not(feature = "cram"))]
    #[rstest]
    #[case::staged(".staged.tmp")]
    #[case::direct("out.cram")]
    fn cram_refused_before_opening(#[case] basename: &str) -> Result<()> {
        let dir = tempfile::tempdir()?;
        let path = dir.path().join(basename);
        std::fs::write(&path, b"keep existing output")?;
        let error = create_raw_alignment_writer(&path, Path::new("out.cram"), &header(), 1, 1)
            .err()
            .expect("CRAM must require the feature");
        let message = error.to_string();
        assert!(message.contains("out.cram"));
        assert!(message.contains("--features cram"));
        assert_eq!(std::fs::read(path)?, b"keep existing output");
        Ok(())
    }
}
