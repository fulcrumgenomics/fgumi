//! BAM, SAM, and FASTQ source steps.
//!
//! - `read_bam` / `read_sam_chunks` — open a BAM/SAM input and emit blocks
//!   or line-aligned chunks.
//! - `read_fastq` — read raw FASTQ bytes, cut on whole-record boundaries.
//! - `fastq_zip` / `parse_zip_fastq` — shared K-way zip primitives and the
//!   parallel parse-and-zip step.
//! - `zip_raw_fastq_k` — the K-input lockstep aligner (and the 1-input adapter)
//!   that joins per-stream raw chunks for the parallel parse-and-zip step.

#[cfg(test)]
mod chain_tests;
pub mod fastq_bgzf;
pub mod fastq_zip;
pub mod find_fastq_boundaries;
pub mod parse_zip_fastq;
pub mod read_bam;
pub mod read_fastq;
pub mod read_sam_chunks;
pub mod zip_raw_fastq_k;

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

/// The BAM/CRAM-opener counterpart to [`open_error`], preserving the underlying
/// [`io::ErrorKind`] across an `anyhow::Error`; `format` (`"BAM"`, `"CRAM"`)
/// names what was being opened.
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
fn open_input_error(format: &str, path: &std::path::Path, source: &anyhow::Error) -> io::Error {
    let message = format!("open {format} {}: {source:#}", path.display());
    match source.chain().find_map(|e| e.downcast_ref::<io::Error>()) {
        Some(io_err) => io::Error::new(io_err.kind(), message),
        None => io::Error::other(message),
    }
}

/// The error for CRAM input `path` in a build without the `cram` feature.
#[cfg(not(feature = "cram"))]
fn cram_unsupported(path: &Path) -> io::Error {
    io::Error::new(
        io::ErrorKind::InvalidData,
        fgumi_bam_io::cram_input_unsupported(path).to_string(),
    )
}

/// Parse a whole-record-aligned FASTQ chunk into records, returning them with
/// the summed `HeapSize` of the parsed records.
///
/// `data` must contain a whole number of 4-line FASTQ records, each terminated
/// by `\n` (as produced by `read_fastq_raw_bytes_from_bufread`). Record
/// boundaries are located with the SIMD lexer (`find_record_offsets`) — a single
/// vectorized pass over the chunk rather than a per-byte scalar scan.
///
/// `step_name` appears in the truncation error so the caller
/// (`ParseAndZipFastqN`) is identifiable in logs.
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
    /// CRAM input (feature `cram`), positioned after its header. Read slice by
    /// slice by `ReadCramSlices` and decoded on the pool by `DecodeCramSlices`
    /// into header-less framed BAM records (pair with
    /// `FindBamBoundaries::new_no_header`).
    #[cfg(feature = "cram")]
    Cram { reader: fgumi_htslib::CramReader, header: sam::Header, label: String },
}

impl InputSource {
    /// Open `path` and detect the input shape.
    ///
    /// Detection: `is_stdin_path` chooses stdin vs file; a `.bam` extension
    /// selects the BAM reader (itself content-aware — see below) unless, with
    /// the `cram` feature, the file's bytes are CRAM; anything else (no
    /// extension, or `.sam`) is classified by content through
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
        Self::open_with_read_buffer(path, opts, 0)
    }

    /// [`InputSource::open_with_opts`] that reads a streamed input (stdin, a
    /// FIFO, a `/dev/fd/N` from process substitution — see
    /// [`fgumi_bam_io::is_stream_input`]) through a
    /// [`fgumi_bam_io::read_buffer::ReadBuffer`] of `read_buffer` bytes, so the
    /// process writing it is not blocked while the consumer pauses. `0` reads
    /// it directly. A regular file, including stdin redirected from one,
    /// ignores `read_buffer`: nothing upstream waits on it.
    ///
    /// # Errors
    ///
    /// Returns I/O errors from file open, stdin read, header parse, or starting
    /// the read buffer's thread.
    pub(crate) fn open_with_read_buffer<P: AsRef<Path>>(
        path: P,
        opts: fgumi_bam_io::PipelineReaderOpts,
        read_buffer: usize,
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
            // content decides. Stdin redirected from a regular file has no
            // writer to wait on, so it gets no read buffer.
            let read_buffer = if read_buffer > 0 && fgumi_bam_io::is_stream_input(path) {
                read_buffer
            } else {
                0
            };
            open_stream_by_content(
                stream_reader(Box::new(io::stdin()), read_buffer, opts.async_reader)?,
                "stdin",
            )
        } else if bam_by_suffix(path, suffix_is_bam) {
            // BAM file: parse header via the existing helper, which
            // handles the seek-to-0 + opt.async_reader wiring. A FIFO named
            // `*.bam` goes through content classification below instead, retaining
            // its original format and prefix without reopening the stream.
            let (reader, header) =
                fgumi_bam_io::create_bam_reader_for_pipeline_with_opts(path, opts)
                    .map_err(|e| open_input_error("BAM", path, &e))?;
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
            // `is_ok_and` picks the safe direction if it somehow does: the
            // in-stream path always yields correct bytes and merely forgoes
            // reopening, whereas a wrong "regular" answer corrupts input.
            if !file.metadata().is_ok_and(|m| m.is_file()) {
                // Non-regular: classify in-stream, from the descriptor already
                // open, as the stdin branch does.
                return open_stream_by_content(
                    stream_reader(Box::new(file), read_buffer, opts.async_reader)?,
                    &path.display().to_string(),
                );
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
                // CRAM is read slice by slice on the pipeline (feature `cram`).
                #[cfg(feature = "cram")]
                InputFormat::Cram => {
                    let replay = Box::new(ChainedReader::new(prefix, file));
                    open_cram_source(replay, path)
                }
                // Without `cram`, refuse CRAM with a clear error.
                #[cfg(not(feature = "cram"))]
                InputFormat::Cram => Err(cram_unsupported(path)),
                // Compressed: re-open through the helper, which seeks to byte 0,
                // preserves the async-reader wiring, and normalizes — so a plain
                // gzip member is decompressed there rather than handed to a block
                // reader that would reject it. The bytes consumed above are simply
                // re-read; a regular file can always be replayed.
                InputFormat::Bgzf | InputFormat::Gzip => {
                    drop(file);
                    let (reader, header) =
                        fgumi_bam_io::create_bam_reader_for_pipeline_with_opts(path, opts)
                            .map_err(|e| open_input_error("BAM", path, &e))?;
                    Ok(InputSource::Bam { reader, header })
                }
                InputFormat::Empty => Err(io::Error::new(
                    io::ErrorKind::UnexpectedEof,
                    format!("input is empty: {} (expected BAM or SAM records)", path.display()),
                )),
            }
        }
    }

    /// Parse the BAM header from an already-open BGZF file, read from its
    /// current offset (the sort's native input path, which reads the body
    /// positionally from the same open file, so header and body cannot come
    /// from two different files). `path` names the file in errors.
    ///
    /// # Errors
    ///
    /// Returns I/O errors from the read or the BAM-header parse.
    pub fn open_bgzf_file(file: File, path: &Path) -> io::Result<Self> {
        open_bam_from_boxed_buf(buffered(file), &path.display().to_string())
    }

    /// Borrow the parsed header. Common for both variants.
    #[must_use]
    pub fn header(&self) -> &sam::Header {
        match self {
            InputSource::Bam { header, .. } | InputSource::Sam { header, .. } => header,
            #[cfg(feature = "cram")]
            InputSource::Cram { header, .. } => header,
        }
    }
}

/// Whether `path` takes the `.bam` suffix fast path: a regular file named
/// `*.bam` (a stream named so is classified in-stream instead). With the `cram`
/// feature, a regular `*.bam` file whose bytes are CRAM does not: content
/// decides, so it opens as `InputSource::Cram` and the source reports CRAM.
fn bam_by_suffix(path: &Path, suffix_is_bam: bool) -> bool {
    let by_suffix = suffix_is_bam && !fgumi_bam_io::is_stream_input(path);
    #[cfg(feature = "cram")]
    let by_suffix = by_suffix && !regular_file_is_cram(path);
    by_suffix
}

/// Whether `path` is a regular file whose content (not its name) is CRAM.
///
/// Streams (stdin, FIFOs, `/dev/fd/N`) answer `false` without being read:
/// consuming their prefix here would lose it before the source reads them.
/// An unreadable file also answers `false`; opening the source reports why.
#[cfg(feature = "cram")]
pub(crate) fn regular_file_is_cram(path: &Path) -> bool {
    if fgumi_bam_io::is_stdin_path(path) || !std::fs::metadata(path).is_ok_and(|m| m.is_file()) {
        return false;
    }
    let Ok(mut file) = File::open(path) else { return false };
    take_format_prefix(&mut file).is_ok_and(|(_, format)| format == InputFormat::Cram)
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

/// Open a CRAM stream (its sniffed leading bytes replayed) as a pipeline
/// [`InputSource::Cram`], with the configured reference.
#[cfg(feature = "cram")]
fn open_cram_source(reader: Box<dyn Read + Send>, path: &Path) -> io::Result<InputSource> {
    let (reader, header) = fgumi_bam_io::cram_input::open_cram_reader(reader, path)
        .map_err(|e| open_input_error("CRAM", path, &e))?;
    Ok(InputSource::Cram { reader, header, label: path.display().to_string() })
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
///
/// `subject` names the stream in a BAM-header error (`stdin`, or the FIFO's
/// path).
fn open_stream_by_content(reader: Box<dyn Read + Send>, subject: &str) -> io::Result<InputSource> {
    let mut reader = reader;
    let (prefix, format) = take_format_prefix(&mut reader)?;
    match format {
        InputFormat::Empty => Err(io::Error::new(
            io::ErrorKind::UnexpectedEof,
            "input stream is empty (expected BAM or SAM records)",
        )),
        InputFormat::Bgzf => {
            open_bam_from_boxed_buf(buffered(ChainedReader::new(prefix, reader)), subject)
        }
        InputFormat::Text => open_sam_from_boxed_buf(buffered(ChainedReader::new(prefix, reader))),
        // CRAM: with the `cram` feature, read slice by slice on the pipeline
        // (`InputSource::Cram`); without it, the shared normalizer reports an
        // error naming the feature.
        #[cfg(feature = "cram")]
        InputFormat::Cram => {
            open_cram_source(Box::new(ChainedReader::new(prefix, reader)), Path::new(subject))
        }
        #[cfg(not(feature = "cram"))]
        InputFormat::Cram => Err(cram_unsupported(Path::new(subject))),
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
                InputFormat::Bgzf => open_bam_from_boxed_buf(
                    buffered(ChainedReader::new(inner_prefix, decoded)),
                    subject,
                ),
                // A BAM whose BGZF framing a gzip decompressor stripped: re-frame
                // it, as `normalize_to_bgzf` does, rather than parse it as SAM.
                InputFormat::Text if fgumi_bam_io::is_unframed_bam(&inner_prefix) => {
                    open_bam_from_boxed_buf(
                        buffered(fgumi_bam_io::frame_unframed_bam(ChainedReader::new(
                            inner_prefix,
                            decoded,
                        ))),
                        subject,
                    )
                }
                InputFormat::Text => {
                    open_sam_from_boxed_buf(buffered(ChainedReader::new(inner_prefix, decoded)))
                }
                InputFormat::Gzip => Err(io::Error::new(
                    io::ErrorKind::InvalidData,
                    "input stream is gzip-compressed more than once; decompress it before \
                     piping it in",
                )),
                InputFormat::Cram => Err(io::Error::new(
                    io::ErrorKind::InvalidData,
                    "input stream is a gzip-compressed CRAM; decompress it before piping it in",
                )),
                InputFormat::Empty => Err(io::Error::new(
                    io::ErrorKind::UnexpectedEof,
                    "input stream decompresses to an empty gzip member",
                )),
            }
        }
    }
}

/// Read a streamed input through a [`fgumi_bam_io::read_buffer::ReadBuffer`]
/// of `budget` bytes. With no budget (`0`), `async_reader` (`--async-reader`)
/// prefetches it on a [`fgumi_bam_io::prefetch_reader::PrefetchReader`]
/// thread, as the path opener does for a regular file; otherwise it is read
/// directly.
///
/// # Errors
///
/// Returns the error from starting the read buffer's thread.
fn stream_reader(
    reader: Box<dyn Read + Send>,
    budget: usize,
    async_reader: bool,
) -> io::Result<Box<dyn Read + Send>> {
    if budget > 0 {
        Ok(Box::new(fgumi_bam_io::read_buffer::ReadBuffer::try_new(reader, budget)?))
    } else if async_reader {
        Ok(Box::new(fgumi_bam_io::prefetch_reader::PrefetchReader::new(reader)))
    } else {
        Ok(reader)
    }
}

/// Box a reader as the buffered stream both openers take.
fn buffered<R: Read + Send + 'static>(reader: R) -> Box<dyn BufRead + Send> {
    Box::new(BufReader::with_capacity(2 * 1024 * 1024, reader))
}

/// Wrap an already-buffered BAM stream (stdin, a FIFO, or an open file) into a
/// BAM `InputSource`; `subject` names it in the header-parse error. The BGZF
/// reader's header-probe goes through a `TeeReader` so we can replay the
/// consumed header bytes ahead of the remaining stream via `ChainedReader`.
/// Reopenable files otherwise take the re-open-at-byte-0 path via
/// `create_bam_reader_for_pipeline_with_opts`.
fn open_bam_from_boxed_buf(buf: Box<dyn BufRead + Send>, subject: &str) -> io::Result<InputSource> {
    let tee = TeeReader::new(buf);
    let bgzf = BgzfReader::new(tee);
    let mut bam_reader = noodles::bam::io::Reader::from(bgzf);
    let header = bam_reader
        .read_header()
        // Preserve the source kind: a truncated header arrives as
        // `UnexpectedEof`, and flattening it to `Other` loses the one signal a
        // caller can branch on. Same policy as `open_error`/`open_input_error`.
        .map_err(|e| io::Error::new(e.kind(), format!("read BAM header from {subject}: {e}")))?;
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
        open_via_named_fifo(bytes, dir, "stream") // deliberately no extension
    }

    /// [`open_via_fifo`] with the FIFO named `name`.
    #[cfg(unix)]
    fn open_via_named_fifo(
        bytes: Vec<u8>,
        dir: &tempfile::TempDir,
        name: &str,
    ) -> io::Result<InputSource> {
        let fifo = dir.path().join(name);
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

    /// Counts the bytes taken from a stream, to tell read-ahead from on-demand reads.
    struct CountingStream {
        inner: io::Cursor<Vec<u8>>,
        taken: std::sync::Arc<std::sync::atomic::AtomicUsize>,
    }

    impl Read for CountingStream {
        fn read(&mut self, buf: &mut [u8]) -> io::Result<usize> {
            let n = self.inner.read(buf)?;
            self.taken.fetch_add(n, std::sync::atomic::Ordering::SeqCst);
            Ok(n)
        }
    }

    /// `--async-reader` prefetches a stream (a FIFO or stdin) on its own thread
    /// when no read buffer is set, as it does for a regular file; without it,
    /// the stream is read only as the consumer asks.
    #[rstest]
    #[case::async_reader(true, true)]
    #[case::sync_reader(false, false)]
    fn stream_reader_prefetches_only_with_async_reader(
        #[case] async_reader: bool,
        #[case] expect_read_ahead: bool,
    ) {
        let taken = std::sync::Arc::new(std::sync::atomic::AtomicUsize::new(0));
        let stream = CountingStream {
            inner: io::Cursor::new(vec![7u8; 1 << 20]),
            taken: std::sync::Arc::clone(&taken),
        };
        let mut reader = stream_reader(Box::new(stream), 0, async_reader).unwrap();
        let mut one = [0u8; 1];
        reader.read_exact(&mut one).unwrap();
        let taken = taken.load(std::sync::atomic::Ordering::SeqCst);
        assert_eq!(taken > 1, expect_read_ahead, "{taken} bytes taken for a 1-byte read");
    }

    /// Plain gzip around an *unframed* BAM byte stream (what `gzip -dc` makes of
    /// a BAM, recompressed) is re-framed and read as BAM over a pipe, whatever
    /// the pipe is named: a `.bam` FIFO is classified in-stream too, and must
    /// not lose the re-framing a regular `.bam` file gets from the normalizer.
    #[cfg(unix)]
    #[rstest]
    #[case::suffix_less("stream")]
    #[case::bam_suffix("stream.bam")]
    fn plain_gzip_of_unframed_bam_over_a_pipe_is_read_as_bam(#[case] name: &str) {
        use flate2::read::MultiGzDecoder;
        use flate2::write::GzEncoder;
        use std::io::Read as _;

        let mut unframed = Vec::new();
        MultiGzDecoder::new(std::fs::File::open(write_bam("")).unwrap())
            .read_to_end(&mut unframed)
            .expect("strip the BGZF framing");
        assert!(unframed.starts_with(b"BAM\x01"), "the fixture is an unframed BAM stream");
        let mut encoder = GzEncoder::new(Vec::new(), flate2::Compression::fast());
        encoder.write_all(&unframed).expect("gzip the unframed BAM");
        let gzipped = encoder.finish().expect("finish gzip member");

        let dir = tempfile::tempdir().unwrap();
        let source = open_via_named_fifo(gzipped, &dir, name).expect("gzipped unframed BAM opens");

        assert!(matches!(source, InputSource::Bam { .. }));
        // The fixture is header-only, so its one reference is the whole content.
        let refs = source.header().reference_sequences();
        assert_eq!(refs.len(), 1);
        let (name, reference) = refs.first().expect("one reference");
        assert_eq!(name.as_slice(), b"chr1");
        assert_eq!(reference.length().get(), 1000);
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

        let source = open_bam_from_boxed_buf(buf, "stdin").unwrap();
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
    /// branch through `open_input_error`, which has to recover the kind from an
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

    /// Write a one-record, reference-free CRAM to `path`.
    #[cfg(feature = "cram")]
    fn write_tiny_cram(path: &Path) {
        use fgumi_htslib::{CramWriteOptions, CramWriter};
        let options = CramWriteOptions { no_ref: true, ..CramWriteOptions::default() };
        let mut writer = CramWriter::create(path, "@HD\tVN:1.6\n", &options).unwrap();
        let record = fgumi_raw_bam::SamBuilder::new()
            .read_name(b"r0")
            .flags(4)
            .sequence(b"ACGT")
            .qualities(&[30; 4])
            .build();
        writer.write_raw(record.as_ref()).unwrap();
        let (jobs, finisher) = writer.seal().unwrap();
        for job in jobs {
            job.encode().unwrap().write().unwrap();
        }
        finisher.finish().unwrap();
    }

    /// A regular file's content, not its name, decides whether it is CRAM:
    /// a CRAM named `*.bam` opens as CRAM (and reports so), like one named
    /// `*.cram` or with no suffix.
    #[cfg(feature = "cram")]
    #[rstest]
    #[case::mislabeled_bam("input.bam")]
    #[case::cram_suffix("input.cram")]
    #[case::no_suffix("input")]
    fn cram_content_opens_as_cram_whatever_the_name(#[case] name: &str) {
        let dir = tempfile::tempdir().unwrap();
        let path = dir.path().join(name);
        write_tiny_cram(&path);
        assert!(regular_file_is_cram(&path));
        let source = InputSource::open(&path).unwrap();
        assert!(matches!(source, InputSource::Cram { .. }), "{name} must open as CRAM");
    }

    /// Non-CRAM regular files, and paths that are not regular files, are not
    /// CRAM (and a missing path is not opened to find out).
    #[cfg(feature = "cram")]
    #[rstest]
    #[case::bam_file(true)]
    #[case::missing(false)]
    fn non_cram_paths_are_not_cram(#[case] write_bam: bool) {
        let dir = tempfile::tempdir().unwrap();
        let path = dir.path().join("input.bam");
        if write_bam {
            let header = noodles::sam::Header::default();
            let writer = fgumi_bam_io::create_raw_bam_writer(&path, &header, 1, 1).unwrap();
            writer.finish().unwrap();
        }
        assert!(!regular_file_is_cram(&path));
        assert!(!regular_file_is_cram(dir.path()), "a directory is not a regular file");
    }
}
