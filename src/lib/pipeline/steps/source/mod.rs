//! BAM, SAM, and FASTQ source steps.
//!
//! - `read_bam` / `read_sam_chunks` — open a BAM/SAM input and emit blocks
//!   or line-aligned chunks.
//! - `read_fastq` — read raw FASTQ bytes, cut on whole-record boundaries.
//! - `parse_fastq` / `parse_zip_fastq` — parse those bytes into records.
//! - `pair_fastq` / `zip_fastq` — join per-stream chunks into templates.

#[cfg(test)]
mod chain_tests;
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

/// Wrap a `File::open` failure with the path while **preserving the original
/// [`io::ErrorKind`]**.
///
/// `io::Error::other` would flatten every kind to [`io::ErrorKind::Other`], so a
/// missing input would stop reporting as [`io::ErrorKind::NotFound`] and callers
/// (and users) would lose the distinction between "no such file" and
/// "permission denied".
fn open_error(path: &std::path::Path, source: &io::Error) -> io::Error {
    io::Error::new(source.kind(), format!("open {}: {source}", path.display()))
}

/// The BAM-opener counterpart to [`open_error`], preserving the underlying
/// [`io::ErrorKind`] across an `anyhow::Error`.
///
/// `create_bam_reader_for_pipeline_with_opts` returns `anyhow::Error`, so
/// `io::Error::other` on it would flatten `NotFound` to `Other` — exactly the
/// loss `open_error` exists to prevent on the plain `File::open` paths, just one
/// error type further away.
///
/// The `io::Error` is usually **not** the top of the chain: the opener adds
/// `.with_context(..)`, so a top-level `downcast` would miss it. Walking
/// `chain()` finds the first I/O error wherever it sits, and a non-I/O failure
/// (a malformed header, say) still falls back to `Other`, which is accurate for
/// it.
fn open_bam_error(path: &std::path::Path, source: &anyhow::Error) -> io::Error {
    let message = format!("open BAM {}: {source}", path.display());
    match source.chain().find_map(|e| e.downcast_ref::<io::Error>()) {
        Some(io_err) => io::Error::new(io_err.kind(), message),
        None => io::Error::other(message),
    }
}

/// Parse a whole-record-aligned FASTQ chunk into records, returning them with
/// the summed `HeapSize` of the parsed records.
///
/// `data` must contain a whole number of 4-line FASTQ records, each terminated
/// by `\n` (as produced by `read_fastq_raw_bytes_from_bufread`). Record
/// boundaries are located with the SIMD lexer (`find_record_offsets`) — a single
/// vectorized pass over the chunk rather than a per-byte scalar scan.
///
/// `step_name` appears in the truncation error so the two callers
/// (`ParseFastqChunks` and `ParseAndZipFastq`) stay distinguishable in logs.
///
/// # Errors
///
/// Returns `InvalidData` when the chunk ends mid-record, i.e. the reader's
/// whole-record alignment guarantee was violated upstream.
pub(crate) fn parse_fastq_chunk(
    data: &[u8],
    step_name: &str,
) -> io::Result<(Vec<crate::fastq_parse::FastqRecord>, usize)> {
    use crate::pipeline::core::item::HeapSize as _;

    // `find_record_offsets` returns `[0, end_1, .., end_n]` where `end_n` is one
    // byte past the final *complete* record; trailing incomplete bytes are
    // excluded. The reader guarantees whole-record-aligned chunks, so the last
    // offset must equal `data.len()` — anything else is a truncated record.
    let offsets = fgumi_simd_fastq::find_record_offsets(data);
    let last = offsets.last().copied().unwrap_or(0);
    if last != data.len() {
        return Err(io::Error::new(
            io::ErrorKind::InvalidData,
            format!(
                "{step_name}: truncated FASTQ record ({} trailing bytes after the last \
                 complete record)",
                data.len() - last
            ),
        ));
    }

    let num_records = offsets.len().saturating_sub(1);
    let mut records = Vec::with_capacity(num_records);
    let mut total_bytes = 0usize;
    for window in offsets.windows(2) {
        let record = crate::fastq_parse::FastqRecord::from_slice(&data[window[0]..window[1]])?;
        total_bytes += record.heap_size();
        records.push(record);
    }
    Ok((records, total_bytes))
}

use fgumi_bam_io::{ChainedReader, InputFormat, TeeReader};
use flate2::read::MultiGzDecoder;
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
    /// Detection: `is_stdin_path` chooses stdin vs file; a `.bam` extension
    /// selects the BAM reader (itself content-aware — see below); anything
    /// else (no extension, or `.sam`) is classified by content through
    /// [`fgumi_bam_io::classify_input`].
    ///
    /// Content classification never reduces to a magic-byte test. A leading
    /// `\x1f` says "gzip", not "BGZF", and plain gzip is a supported input on
    /// every other fgumi path — so treating it as BAM would hand a single
    /// deflate stream to a block reader that rejects it. Reopenable inputs get
    /// this for free from `create_bam_reader_for_pipeline_with_opts`, which
    /// normalizes; pipes are classified (and decompressed if needed) in
    /// `open_stream_by_content`.
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

        // Suffix hint: `.bam` → BAM, handed to `create_bam_reader_for_pipeline_with_opts`,
        // which is itself content-aware (its `normalize_to_bgzf` transcodes SAM text or
        // plain gzip found under a `.bam` name, so a mislabeled BAM self-corrects there).
        // A `.sam` suffix carries no equivalent self-correction below it — `sam::io::Reader`
        // has no fallback for BGZF/BAM bytes — so it is deliberately NOT special-cased here;
        // it falls through to the same content classification as no extension at all. Content
        // always decides in the end, including for `.gz`, whose extension says only that it
        // is compressed, not by what.
        let extension = path.extension();
        let suffix_is_bam = extension.is_some_and(|e| e.eq_ignore_ascii_case("bam"));

        if is_stdin {
            // The path is just `-`, so there is no suffix to consult and the
            // content decides. `PipelineReaderOpts` is not threaded through here:
            // there is no path to reopen, so there is no async-reader wiring and
            // no positional (scatter) reads. Warn if `--read-streams` explicitly
            // asked for concurrency stdin can't provide (the `Auto` default falls
            // back silently) so this path gives the same feedback as the others.
            fgumi_bam_io::scatter_reader::warn_read_streams_unavailable(
                opts.read_streams,
                "stdin",
                "is not a seekable regular file",
            );
            open_stream_by_content(Box::new(io::stdin()))
        } else if suffix_is_bam {
            // BAM file: parse header via the existing helper, which
            // handles the seek-to-0 + opt.async_reader wiring.
            let (reader, header) =
                fgumi_bam_io::create_bam_reader_for_pipeline_with_opts(path, opts)
                    .map_err(|e| open_bam_error(path, &e))?;
            Ok(InputSource::Bam { reader, header })
        } else {
            // No suffix hint, or a `.sam` suffix that content may contradict:
            // the content decides.
            let file = File::open(path).map_err(|e| open_error(path, &e))?;
            // Whether re-opening `path` can replay the bytes classification
            // consumes. Only a regular file can: for a FIFO or a `/dev/fd/N`
            // handed over by process substitution, reading those bytes drains
            // them out of the pipe for good, and a second open yields the *rest*
            // of the stream, not the start of it — silently losing the header.
            //
            // An `fstat` on an already-open descriptor essentially cannot fail, but
            // `unwrap_or(false)` picks the safe direction if it somehow does: the
            // in-stream path always yields correct bytes and merely forgoes the
            // async-reader wiring, whereas a wrong "regular" answer corrupts input.
            if !file.metadata().map(|m| m.is_file()).unwrap_or(false) {
                // Non-regular: classify in-stream. Same trade-off the stdin branch
                // makes — no async-reader wiring and no positional reads, because
                // there is no reopenable path. Warn on an explicit `--read-streams`
                // the input can't honor (the `Auto` default stays silent).
                fgumi_bam_io::scatter_reader::warn_read_streams_unavailable(
                    opts.read_streams,
                    &path.display().to_string(),
                    "is not a seekable regular file",
                );
                return open_stream_by_content(Box::new(file));
            }

            // Regular file: classify here, then choose between keeping the text
            // in the SAM reader and handing the path to the normalizing helper.
            let mut file = file;
            let (prefix, format) = take_format_prefix(&mut file)?;
            match format {
                // SAM text stays SAM: `ReadSamChunks` parses it directly, and
                // routing it through the helper would transcode it to BGZF and
                // throw that away.
                InputFormat::Text => {
                    open_sam_from_boxed_buf(buffered(ChainedReader::new(prefix, file)))
                }
                // Compressed: re-open through the helper, which seeks to byte 0,
                // preserves the async-reader wiring, and normalizes — so a plain
                // gzip member is decompressed there rather than handed to a block
                // reader that would reject it. The bytes consumed above are simply
                // re-read; a regular file can always be replayed.
                InputFormat::Bgzf | InputFormat::Gzip => {
                    drop(file);
                    let (reader, header) =
                        fgumi_bam_io::create_bam_reader_for_pipeline_with_opts(path, opts)
                            .map_err(|e| open_bam_error(path, &e))?;
                    Ok(InputSource::Bam { reader, header })
                }
                InputFormat::Empty => Err(io::Error::new(
                    io::ErrorKind::UnexpectedEof,
                    format!("input is empty: {} (expected BAM or SAM records)", path.display()),
                )),
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

/// Consume up to [`fgumi_bam_io::FORMAT_PREFIX_LEN`] bytes and classify them.
///
/// The bytes are returned rather than pushed back so the caller can replay them
/// through a [`ChainedReader`]: these inputs are pipes, which cannot be rewound
/// or reopened. Short reads are looped over, because one `read` on a pipe may
/// return fewer bytes than are available and a short prefix classifies as
/// [`fgumi_bam_io::InputFormat::Gzip`] regardless of what the stream really is
/// — BGZF cannot be recognized from fewer than a full header.
fn take_format_prefix<R: Read + ?Sized>(reader: &mut R) -> io::Result<(Vec<u8>, InputFormat)> {
    let mut prefix = vec![0u8; fgumi_bam_io::FORMAT_PREFIX_LEN];
    let mut filled = 0;
    while filled < prefix.len() {
        match reader.read(&mut prefix[filled..]) {
            Ok(0) => break,
            Ok(n) => filled += n,
            // A signal can interrupt a `read` on a stream that is still
            // perfectly readable. `read_exact` retries this kind rather than
            // propagating it, and a hand-rolled fill loop has to do the same —
            // otherwise a signal during startup fails the open outright.
            Err(e) if e.kind() == io::ErrorKind::Interrupted => {}
            Err(e) => return Err(e),
        }
    }
    prefix.truncate(filled);
    let format = fgumi_bam_io::classify_input(&prefix);
    Ok((prefix, format))
}

/// Open a stream that cannot be reopened — stdin, a FIFO, a `/dev/fd/N` from
/// process substitution — by classifying its content.
///
/// Regular files do not come through here: they go to
/// `create_bam_reader_for_pipeline_with_opts`, which normalizes via
/// `normalize_to_bgzf` and so already decompresses plain gzip and transcodes
/// SAM. A pipe has no such second pass, so the classification and the gzip
/// decompression have to happen here.
///
/// Classification is by content through the shared
/// [`fgumi_bam_io::classify_input`], never by a hand-rolled magic-byte test: a
/// bare `\x1f` check cannot tell BGZF from plain gzip, and answering "BGZF" for
/// a plain-gzip member hands it to a block reader that rejects it — the
/// misleading diagnosis that classifier exists to prevent.
fn open_stream_by_content(reader: Box<dyn Read + Send>) -> io::Result<InputSource> {
    let mut reader = reader;
    let (prefix, format) = take_format_prefix(&mut reader)?;
    match format {
        InputFormat::Empty => Err(io::Error::new(
            io::ErrorKind::UnexpectedEof,
            "input stream is empty (expected BAM or SAM records)",
        )),
        InputFormat::Bgzf => {
            open_bam_from_stdin_boxed_buf(buffered(ChainedReader::new(prefix, reader)))
        }
        InputFormat::Text => open_sam_from_boxed_buf(buffered(ChainedReader::new(prefix, reader))),
        // Plain gzip: one deflate stream, so no block-parallel decode. Decompress
        // at this boundary and classify what comes out, exactly as
        // `normalize_to_bgzf` does for reopenable inputs.
        InputFormat::Gzip => {
            log::warn!(
                "input stream is gzip-compressed rather than BGZF — decompressing it \
                 single-threaded, so --threads will not speed up reading. Recompress \
                 with `bgzip` for block-parallel decode."
            );
            let mut decoded: Box<dyn Read + Send> =
                Box::new(MultiGzDecoder::new(ChainedReader::new(prefix, reader)));
            let (inner_prefix, inner_format) = take_format_prefix(&mut *decoded)?;
            match inner_format {
                InputFormat::Bgzf => open_bam_from_stdin_boxed_buf(buffered(ChainedReader::new(
                    inner_prefix,
                    decoded,
                ))),
                InputFormat::Text => {
                    open_sam_from_boxed_buf(buffered(ChainedReader::new(inner_prefix, decoded)))
                }
                InputFormat::Gzip => Err(io::Error::new(
                    io::ErrorKind::InvalidData,
                    "input stream is gzip-compressed more than once; decompress it before \
                     piping it in",
                )),
                InputFormat::Empty => Err(io::Error::new(
                    io::ErrorKind::UnexpectedEof,
                    "input stream decompresses to an empty gzip member",
                )),
            }
        }
    }
}

/// Box a reader as the buffered stream both openers take.
fn buffered<R: Read + Send + 'static>(reader: R) -> Box<dyn BufRead + Send> {
    Box::new(BufReader::with_capacity(2 * 1024 * 1024, reader))
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
        // Preserve the source kind: a truncated header arrives as
        // `UnexpectedEof`, and flattening it to `Other` loses the one signal a
        // caller can branch on. Same policy as `open_error`/`open_bam_error`.
        .map_err(|e| io::Error::new(e.kind(), format!("read BAM header from stdin: {e}")))?;
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
    let header = sam_reader.read_header().map_err(|e| {
        // See the BAM header read above: keep the source `ErrorKind`.
        io::Error::new(e.kind(), format!("read SAM header: {e}"))
    })?;
    // A SAM header parse consumes only lines beginning with `@`, so arbitrary
    // text yields an *empty* header rather than an error — mirroring
    // `fgumi_bam_io::sam_input`'s `SamToBamStream`, which has the identical
    // check for the identical reason. Left unchecked, this would silently
    // treat a truncated or wrong-format file as a valid, headerless SAM
    // stream and only fail later — confusingly, on the first "record" line
    // (which is really just more non-SAM text) — instead of naming the input
    // itself as the problem. Before SAM was accepted this input was rejected
    // outright for not being BGZF, and it stays rejected.
    if header.is_empty() {
        return Err(io::Error::new(
            io::ErrorKind::InvalidData,
            "Input is neither BAM nor SAM (no SAM header found; SAM input must begin with \
             an `@` header line)",
        ));
    }
    Ok(InputSource::Sam { reader: sam_reader, header })
}

#[cfg(test)]
mod tests {
    use rstest::rstest;

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

    /// Write a minimal single-reference BAM and return its path.
    fn write_bam(suffix: &str) -> tempfile::TempPath {
        let path = tempfile::NamedTempFile::with_suffix(suffix).unwrap().into_temp_path();
        let header = noodles::sam::Header::builder()
            .add_reference_sequence(
                bstr::BString::from("chr1"),
                noodles::sam::header::record::value::Map::<
                    noodles::sam::header::record::value::map::ReferenceSequence,
                >::new(std::num::NonZeroUsize::new(1000).unwrap()),
            )
            .build();
        let writer = fgumi_bam_io::create_raw_bam_writer(&path, &header, 1, 1).unwrap();
        writer.finish().unwrap();
        path
    }

    /// The magic-byte fallback must resolve BGZF content to BAM, not just
    /// resolve text content to SAM. Without this the suffix-less BAM branch —
    /// which re-opens through the async-reader-aware helper — is never taken.
    #[test]
    fn falls_back_to_magic_byte_for_bam_without_extension() {
        let path = write_bam("");
        let source = InputSource::open(&path).unwrap();
        assert!(matches!(source, InputSource::Bam { .. }));
        assert_eq!(source.header().reference_sequences().len(), 1);
    }

    /// Serve `bytes` over a FIFO at a suffix-less path and open it.
    ///
    /// The FIFO is what makes the input non-reopenable, which is the branch
    /// under test: a regular file would be handed to the normalizing helper
    /// instead. Opening a FIFO for writing blocks until a reader opens it, so
    /// the writer runs concurrently with the open.
    #[cfg(unix)]
    fn open_via_fifo(bytes: Vec<u8>, dir: &tempfile::TempDir) -> io::Result<InputSource> {
        let fifo = dir.path().join("stream"); // deliberately no extension
        let status = std::process::Command::new("mkfifo").arg(&fifo).status().expect("mkfifo runs");
        assert!(status.success(), "mkfifo failed for {}", fifo.display());

        let writer_path = fifo.clone();
        let writer = std::thread::spawn(move || {
            let mut f = std::fs::File::create(&writer_path).expect("open fifo for write");
            let _ = f.write_all(&bytes); // the reader may reject before draining
        });

        let opened = InputSource::open(&fifo);
        let _ = writer.join();
        opened
    }

    /// Plain-gzip SAM over a pipe must be decompressed and read as SAM.
    ///
    /// A bare `\x1f` check calls this BGZF and hands it to the BAM block reader,
    /// which fails with a header error on a stream that is perfectly readable —
    /// the misleading diagnosis `classify_input` exists to prevent. Plain gzip is
    /// a supported verdict on every other fgumi input path, so it must be one
    /// here too.
    #[cfg(unix)]
    #[test]
    fn suffix_less_plain_gzip_sam_over_a_pipe_is_read_as_sam() {
        use flate2::write::GzEncoder;

        let mut encoder = GzEncoder::new(Vec::new(), flate2::Compression::fast());
        encoder.write_all(SAM_TEXT.as_bytes()).expect("gzip the SAM text");
        let gzipped = encoder.finish().expect("finish gzip member");
        assert_eq!(gzipped[0], 0x1f, "the fixture must start with the gzip magic byte");

        let dir = tempfile::tempdir().unwrap();
        let source = open_via_fifo(gzipped, &dir).expect("gzipped SAM opens");

        let InputSource::Sam { header, .. } = source else {
            panic!("expected a SAM source; plain gzip must not be routed to the BAM reader");
        };
        assert_eq!(header.reference_sequences().len(), 1);
    }

    /// The error arms of `open_stream_by_content` and of the suffix-less regular
    /// -file branch, which the happy-path tests above never reach.
    ///
    /// Each is a named diagnosis rather than a downstream parse failure, so the
    /// message is what is asserted — a caller piping the wrong thing in should be
    /// told what is wrong with it, not handed "invalid BGZF header".
    #[cfg(unix)]
    #[test]
    fn an_empty_stream_is_named_as_empty() {
        let dir = tempfile::tempdir().unwrap();
        let Err(err) = open_via_fifo(Vec::new(), &dir) else {
            panic!("an empty stream must be an error");
        };

        assert_eq!(err.kind(), io::ErrorKind::UnexpectedEof);
        assert!(err.to_string().contains("empty"), "got: {err}");
    }

    /// gzip wrapped around gzip: the inner classification sees gzip again, which
    /// no fgumi reader decodes. Named explicitly rather than left to fail as a
    /// corrupt BAM header.
    #[cfg(unix)]
    #[test]
    fn a_doubly_gzipped_stream_is_rejected_by_name() {
        use flate2::write::GzEncoder;

        let gzip_once = |bytes: &[u8]| {
            let mut e = GzEncoder::new(Vec::new(), flate2::Compression::fast());
            e.write_all(bytes).expect("gzip");
            e.finish().expect("finish gzip member")
        };
        let doubled = gzip_once(&gzip_once(SAM_TEXT.as_bytes()));

        let dir = tempfile::tempdir().unwrap();
        let Err(err) = open_via_fifo(doubled, &dir) else {
            panic!("a doubly-gzipped stream must be an error");
        };

        assert_eq!(err.kind(), io::ErrorKind::InvalidData);
        assert!(err.to_string().contains("more than once"), "got: {err}");
    }

    /// A well-formed gzip member carrying nothing: the outer classification
    /// succeeds and the *inner* one is `Empty`, a different branch from the
    /// empty-stream case above.
    #[cfg(unix)]
    #[test]
    fn a_gzip_member_that_decompresses_to_nothing_is_named() {
        use flate2::write::GzEncoder;

        let encoder = GzEncoder::new(Vec::new(), flate2::Compression::fast());
        let empty_member = encoder.finish().expect("finish an empty gzip member");
        assert_eq!(empty_member[0], 0x1f, "still a well-formed gzip member");

        let dir = tempfile::tempdir().unwrap();
        let Err(err) = open_via_fifo(empty_member, &dir) else {
            panic!("an empty gzip member must be an error");
        };

        assert_eq!(err.kind(), io::ErrorKind::UnexpectedEof);
        assert!(err.to_string().contains("empty gzip member"), "got: {err}");
    }

    /// The regular-file branch has its own `Empty` arm, reached without any pipe.
    #[test]
    fn an_empty_suffix_less_file_is_named_as_empty() {
        let path = tempfile::NamedTempFile::new().unwrap().into_temp_path();
        let Err(err) = InputSource::open(&path) else {
            panic!("an empty file must be an error");
        };

        assert_eq!(err.kind(), io::ErrorKind::UnexpectedEof);
        assert!(err.to_string().contains("input is empty"), "got: {err}");
    }

    /// The same path must still resolve gzipped *BAM* to BAM: the verdict comes
    /// from what the member decompresses to, not from the outer wrapper.
    #[cfg(unix)]
    #[test]
    fn suffix_less_plain_gzip_bam_over_a_pipe_is_read_as_bam() {
        use flate2::write::GzEncoder;

        let bam_bytes = std::fs::read(write_bam("")).unwrap();
        let mut encoder = GzEncoder::new(Vec::new(), flate2::Compression::fast());
        encoder.write_all(&bam_bytes).expect("gzip the BAM bytes");
        let gzipped = encoder.finish().expect("finish gzip member");

        let dir = tempfile::tempdir().unwrap();
        let source = open_via_fifo(gzipped, &dir).expect("gzipped BAM opens");

        assert!(matches!(source, InputSource::Bam { .. }));
        assert_eq!(source.header().reference_sequences().len(), 1);
    }

    /// A suffix-less BGZF input that is **not** a regular file must be read
    /// through the buffer that already holds the peeked bytes. Re-opening the
    /// path — which is what a regular file gets, to keep the async-reader
    /// wiring — would hand back the stream *minus* the bytes `fill_buf` already
    /// drained out of the pipe, so the BAM header would be gone. This is the
    /// `<(...)` process-substitution and named-pipe case.
    ///
    /// A regression here manifests as a hang rather than an assertion failure:
    /// the re-open blocks forever waiting for a second writer that never comes.
    #[cfg(unix)]
    #[test]
    fn suffix_less_bam_from_a_fifo_keeps_the_peeked_bytes() {
        let bam_bytes = std::fs::read(write_bam("")).unwrap();

        let dir = tempfile::tempdir().unwrap();
        let fifo = dir.path().join("stream"); // deliberately no extension
        let status = std::process::Command::new("mkfifo").arg(&fifo).status().expect("mkfifo runs");
        assert!(status.success(), "mkfifo failed for {}", fifo.display());

        // Opening a FIFO for writing blocks until a reader opens it, so the
        // writer has to run concurrently with the `open` below. The payload is
        // far under the pipe buffer, so the writer never blocks on a full pipe.
        let writer_path = fifo.clone();
        let writer = std::thread::spawn(move || {
            let mut f = std::fs::File::create(&writer_path).expect("open fifo for write");
            f.write_all(&bam_bytes).expect("write bam to fifo");
        });

        let source = InputSource::open(&fifo).expect("fifo opens as BAM");
        assert!(matches!(source, InputSource::Bam { .. }));
        assert_eq!(source.header().reference_sequences().len(), 1);

        writer.join().expect("writer thread joins");
    }

    /// `header()` is the variant-agnostic accessor commands use before they
    /// branch on the input shape, so it must resolve for both variants.
    #[test]
    fn header_accessor_resolves_for_both_variants() {
        let sam_path = tempfile::NamedTempFile::with_suffix(".sam").unwrap().into_temp_path();
        std::fs::File::create(&sam_path).unwrap().write_all(SAM_TEXT.as_bytes()).unwrap();
        let sam = InputSource::open(&sam_path).unwrap();
        assert!(matches!(sam, InputSource::Sam { .. }));
        assert_eq!(sam.header().reference_sequences().len(), 1);

        let bam = InputSource::open(write_bam(".bam")).unwrap();
        assert!(matches!(bam, InputSource::Bam { .. }));
        assert_eq!(bam.header().reference_sequences().len(), 1);
    }

    /// The stdin BAM path cannot re-open by path, so it tees the header-probe
    /// bytes and replays them ahead of the remaining stream. Driving the helper
    /// directly over an in-memory BAM covers that replay without needing a real
    /// stdin: the parsed header must come back, and the returned reader must
    /// still begin at byte 0 of the BAM stream (magic included) so the
    /// downstream boundary step can strip the header itself.
    #[test]
    fn stdin_bam_path_replays_the_consumed_header_bytes() {
        let bytes = std::fs::read(write_bam(".bam")).unwrap();
        let buf: Box<dyn BufRead + Send> = Box::new(BufReader::new(io::Cursor::new(bytes.clone())));

        let source = open_bam_from_stdin_boxed_buf(buf).unwrap();
        let InputSource::Bam { mut reader, header } = source else {
            panic!("expected BAM source");
        };
        assert_eq!(header.reference_sequences().len(), 1);

        // The replayed stream must be byte-identical to the original file: the
        // header probe consumed bytes that `ChainedReader` has to hand back.
        let mut replayed = Vec::new();
        reader.read_to_end(&mut replayed).unwrap();
        assert_eq!(replayed, bytes, "the teed header bytes must be replayed intact");
    }

    /// A `.sam`-suffixed file whose header does not parse must surface an error
    /// rather than yielding a source with an empty header.
    #[test]
    fn a_malformed_sam_header_is_an_error() {
        let path = tempfile::NamedTempFile::with_suffix(".sam").unwrap().into_temp_path();
        std::fs::File::create(&path).unwrap().write_all(b"@HD\tthis is not a valid tag\n").unwrap();
        assert!(InputSource::open(&path).is_err());
    }

    /// Opening a path that does not exist must be an error, not a panic — and
    /// must keep the original [`io::ErrorKind`] and name the path.
    ///
    /// `open_error` exists precisely so `io::Error::other` does not flatten
    /// `NotFound` into `Other`, which would cost callers the distinction between
    /// "no such file" and "permission denied". Asserting only `is_err()` would
    /// let a regression back to `io::Error::other` pass unnoticed, so both
    /// properties are pinned here. One case per branch that can fail to open: the
    /// `.sam` and suffix-less branches go through `open_error`, and the `.bam`
    /// branch through `open_bam_error`, which has to recover the kind from an
    /// `anyhow::Error` chain rather than from an `io::Error` directly.
    #[rstest]
    #[case::sam_suffix("/nonexistent/fgumi-test-input.sam")]
    #[case::bam_suffix("/nonexistent/fgumi-test-input.bam")]
    #[case::suffix_less_path("/nonexistent/fgumi-test-input")]
    fn a_missing_file_reports_not_found_with_its_path(#[case] path: &str) {
        let Err(err) = InputSource::open(path) else {
            panic!("opening {path} must be an error");
        };

        assert_eq!(err.kind(), io::ErrorKind::NotFound, "kind must survive the wrap: {err}");
        assert!(err.to_string().contains(path), "the message must name the path; got: {err}");
    }
}
