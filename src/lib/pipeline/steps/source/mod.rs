//! BAM, SAM, and (future) FASTQ source steps.

pub mod pair_fastq;
pub mod parse_fastq;
pub mod parse_zip_fastq;
pub mod read_bam;
pub mod read_fastq;
pub mod read_sam_chunks;
pub mod zip_fastq;

use std::fs::File;
use std::io::{self, BufRead, BufReader, Read};
use std::path::Path;

use fgumi_bam_io::{ChainedReader, TeeReader};
use noodles::bgzf::io::Reader as BgzfReader;
use noodles::sam;

/// Opaque, opened input source ready to feed into the typed-step pipeline.
///
/// The two variants encode the BAM vs. SAM choice. The BAM variant carries
/// a `Box<dyn Read + Send>` whose first bytes are the BAM header (file
/// seeked to byte 0 for on-disk inputs, or a `TeeReader`/`ChainedReader`
/// replay for stdin) and the parsed header. The SAM variant carries the
/// SAM reader already positioned past its header — record reads go
/// straight into [`read_sam_chunks::ReadSamChunks`].
///
/// Each command's new-pipeline path opens its input via [`InputSource::open`]
/// once, threads the [`sam::Header`] into upstream validation /
/// `@PG`-mutation, then `match`es on the variant to build the appropriate
/// chain preamble (3-step BAM vs. 2-step SAM).
pub enum InputSource {
    /// BAM input (BGZF-compressed). `reader` is positioned at byte 0 of
    /// the BAM stream — the `read_bam::ReadBgzfBlocks` source emits the
    /// header bytes as part of the first block(s); `FindBamBoundaries::new`
    /// strips them downstream.
    Bam { reader: Box<dyn Read + Send>, header: sam::Header },
    /// SAM (text) input. `reader` has consumed the `@`-prefixed header
    /// lines; subsequent `read_record` calls yield the first record.
    /// Pair with `FindBamBoundaries::new_no_header` (the SAM source
    /// emits record bodies directly, no header bytes to strip).
    Sam { reader: sam::io::Reader<Box<dyn BufRead + Send>>, header: sam::Header },
}

impl InputSource {
    /// Open `path` and detect the input shape.
    ///
    /// Detection: `is_stdin_path` chooses stdin vs file; then the file
    /// extension (`.sam`/`.bam`/none) plus a peek at the first byte
    /// (BGZF magic `\x1f`) selects SAM vs BAM. For stdin without a
    /// suffix hint, the BGZF-magic peek is authoritative — SAM text
    /// always starts with `@HD` (or another `@`-prefixed line, or a
    /// non-`@` record), neither of which collides with BGZF.
    ///
    /// # Errors
    ///
    /// Returns I/O errors from file open, stdin read, or header parse.
    pub fn open<P: AsRef<Path>>(path: P) -> io::Result<Self> {
        Self::open_with_opts(path, fgumi_bam_io::PipelineReaderOpts::default())
    }

    /// [`InputSource::open`] with a [`fgumi_bam_io::PipelineReaderOpts`].
    ///
    /// # Errors
    ///
    /// Returns I/O errors from file open, stdin read, or header parse.
    pub fn open_with_opts<P: AsRef<Path>>(
        path: P,
        opts: fgumi_bam_io::PipelineReaderOpts,
    ) -> io::Result<Self> {
        let path = path.as_ref();
        let is_stdin = fgumi_bam_io::is_stdin_path(path);

        // Suffix hint: explicit `.sam` → SAM; `.bam` → BAM. Otherwise we
        // peek at the first byte (BGZF blocks start with `\x1f`; SAM
        // never does).
        let extension = path.extension();
        let suffix_says_sam = extension.is_some_and(|e| e.eq_ignore_ascii_case("sam"));
        let suffix_is_bam = extension.is_some_and(|e| e.eq_ignore_ascii_case("bam"));

        if is_stdin {
            // For stdin we need to peek at the first byte to disambiguate
            // (the path is just `-`, no extension to look at). Use a
            // `BufReader<Stdin>` so the peek is non-destructive.
            let mut buf = BufReader::with_capacity(2 * 1024 * 1024, io::stdin());
            let bgzf_magic = peek_is_bgzf(&mut buf)?;
            if bgzf_magic {
                let _ = opts; // PipelineReaderOpts not threaded through the stdin path (no async-reader wiring)
                open_bam_from_stdin_boxed_buf(Box::new(buf))
            } else {
                open_sam_from_boxed_buf(Box::new(buf))
            }
        } else if suffix_says_sam {
            let file = File::open(path).map_err(io::Error::other)?;
            let buf: Box<dyn BufRead + Send> =
                Box::new(BufReader::with_capacity(2 * 1024 * 1024, file));
            open_sam_from_boxed_buf(buf)
        } else if suffix_is_bam {
            // BAM file: parse header via the existing helper, which
            // handles the seek-to-0 + opt.async_reader wiring.
            let (reader, header) =
                fgumi_bam_io::create_bam_reader_for_pipeline_with_opts(path, opts)
                    .map_err(|e| io::Error::other(format!("open BAM: {e}")))?;
            Ok(InputSource::Bam { reader, header })
        } else {
            // No suffix hint: peek to disambiguate.
            let file = File::open(path).map_err(io::Error::other)?;
            let mut buf = BufReader::with_capacity(2 * 1024 * 1024, file);
            let bgzf_magic = peek_is_bgzf(&mut buf)?;
            if bgzf_magic {
                // Drop the BufReader and let the helper re-open via path
                // so it preserves the async-reader wiring.
                drop(buf);
                let (reader, header) =
                    fgumi_bam_io::create_bam_reader_for_pipeline_with_opts(path, opts)
                        .map_err(|e| io::Error::other(format!("open BAM: {e}")))?;
                Ok(InputSource::Bam { reader, header })
            } else {
                open_sam_from_boxed_buf(Box::new(buf))
            }
        }
    }

    /// Borrow the parsed header. Common for both variants.
    #[must_use]
    pub fn header(&self) -> &sam::Header {
        match self {
            InputSource::Bam { header, .. } | InputSource::Sam { header, .. } => header,
        }
    }
}

/// Peek the first byte of `reader` (via `fill_buf`, so the byte stays in
/// the buffer for the next read) and return `true` if it matches the BGZF
/// magic `\x1f` (the first byte of the gzip header that wraps every BGZF
/// block).
fn peek_is_bgzf<R: BufRead>(reader: &mut R) -> io::Result<bool> {
    let head = reader.fill_buf()?;
    Ok(head.first().copied() == Some(0x1f))
}

/// Wrap an already-buffered stdin-shape BAM reader into a BAM
/// `InputSource`. The BGZF reader's header-probe goes through a
/// `TeeReader` so we can replay the consumed header bytes ahead of the
/// remaining stdin stream via `ChainedReader`. File inputs take the
/// re-open-at-byte-0 path via `create_bam_reader_for_pipeline_with_opts`
/// and never reach this function.
fn open_bam_from_stdin_boxed_buf(buf: Box<dyn BufRead + Send>) -> io::Result<InputSource> {
    let tee = TeeReader::new(buf);
    let bgzf = BgzfReader::new(tee);
    let mut bam_reader = noodles::bam::io::Reader::from(bgzf);
    let header = bam_reader
        .read_header()
        .map_err(|e| io::Error::other(format!("read BAM header from stdin: {e}")))?;
    let bgzf = bam_reader.into_inner();
    let tee = bgzf.into_inner();
    let (buffered, remaining) = tee.into_parts();
    let chained: Box<dyn Read + Send> = Box::new(ChainedReader::new(buffered, remaining));
    Ok(InputSource::Bam { reader: chained, header })
}

/// Wrap an already-buffered SAM-shape reader (file or stdin) into a SAM
/// `InputSource`. Parses the header up front so subsequent record reads
/// start at the first record line.
fn open_sam_from_boxed_buf(buf: Box<dyn BufRead + Send>) -> io::Result<InputSource> {
    let mut sam_reader = sam::io::Reader::new(buf);
    let header =
        sam_reader.read_header().map_err(|e| io::Error::other(format!("read SAM header: {e}")))?;
    Ok(InputSource::Sam { reader: sam_reader, header })
}

#[cfg(test)]
mod tests {
    use super::*;
    use std::io::Write;

    const SAM_TEXT: &str = "\
@HD\tVN:1.6\tSO:unsorted\n\
@SQ\tSN:chr1\tLN:1000\n\
read1\t0\tchr1\t10\t60\t5M\t*\t0\t0\tACGTA\tIIIII\n";

    #[test]
    fn opens_bam_file_by_extension() {
        let path = tempfile::NamedTempFile::with_suffix(".bam").unwrap().into_temp_path();
        let header = noodles::sam::Header::default();
        let writer = fgumi_bam_io::create_raw_bam_writer(&path, &header, 1, 1).unwrap();
        writer.finish().unwrap();

        let source = InputSource::open(&path).unwrap();
        assert!(matches!(source, InputSource::Bam { .. }));
    }

    #[test]
    fn opens_sam_file_by_extension() {
        let path = tempfile::NamedTempFile::with_suffix(".sam").unwrap().into_temp_path();
        std::fs::File::create(&path).unwrap().write_all(SAM_TEXT.as_bytes()).unwrap();

        let source = InputSource::open(&path).unwrap();
        let InputSource::Sam { header, .. } = source else {
            panic!("expected SAM source");
        };
        assert_eq!(header.reference_sequences().len(), 1);
    }

    #[test]
    fn falls_back_to_magic_byte_when_extension_missing() {
        // SAM content, no `.sam` suffix.
        let path = tempfile::NamedTempFile::new().unwrap().into_temp_path();
        std::fs::File::create(&path).unwrap().write_all(SAM_TEXT.as_bytes()).unwrap();

        let source = InputSource::open(&path).unwrap();
        assert!(matches!(source, InputSource::Sam { .. }));
    }
}
