//! Transparent support for uncompressed SAM input.
//!
//! Every fgumi reader factory is a BGZF reader: the single-threaded fast
//! paths hand out raw BAM record bytes, and the multi-threaded pipelines hand
//! out a compressed byte stream that downstream stages BGZF-decode
//! themselves. Teaching each of those consumers to accept a second, text
//! record source would mean a second orchestration per command — the exact
//! kind of split that lets one input format work under one flag and not
//! another.
//!
//! Instead this module normalizes at the boundary: SAM text is parsed and
//! re-encoded on the fly into the BGZF byte stream the consumers already
//! expect, so `-i in.sam` and `-i in.bam` are the same stream by the time any
//! consumer sees them. The transcode uses [`CompressionLevel::NONE`] — the
//! bytes are decoded again a few microseconds later in the same process, so
//! there is nothing to gain from deflating them.
//!
//! Detection is by content, never by file extension: BGZF begins with the
//! gzip magic bytes and SAM text cannot, so a misnamed `reads.bam` holding
//! SAM text and a `reads.sam` holding BGZF both read correctly.

use anyhow::{Context, Result, bail};
use noodles::sam;
use noodles::sam::alignment::io::Write as _;
use noodles_bgzf::io::writer::CompressionLevel;
use std::io::{self, BufReader, Read, Write};
use std::path::Path;
use std::sync::{Arc, Mutex};

use crate::reader::ChainedReader;

/// gzip magic, which every BGZF block starts with.
///
/// SAM text cannot begin with these bytes: a SAM file starts either with `@`
/// (a header line) or with the first character of a read name, and the SAM
/// specification restricts read names to printable ASCII.
const GZIP_MAGIC: [u8; 2] = [0x1f, 0x8b];

/// Read buffer wrapped around SAM text input. SAM is line-oriented and much
/// larger than the equivalent BAM, so it is read in big gulps.
const SAM_READ_BUFFER_SIZE: usize = 2 * 1024 * 1024;

/// Upper bound on how many SAM records are transcoded per refill of the output
/// buffer.
///
/// The BGZF encoder emits nothing until a block fills, so a refill stops at the
/// first emitted bytes and only falls back on this bound while the encoder is
/// still filling its first block. That keeps the number of round trips through
/// [`Read::read`] small for short reads without letting long reads — where a
/// single record can be hundreds of kilobytes — stage a whole batch before any
/// byte is handed back.
const RECORDS_PER_BATCH: usize = 1024;

/// A byte sink shared between the BGZF encoder and the reader draining it.
///
/// [`Arc`]/[`Mutex`] rather than [`std::rc::Rc`]/[`std::cell::RefCell`]
/// because the resulting stream is handed out as `Box<dyn Read + Send>`. The
/// lock is uncontended — the encoder and the drain both run on the caller's
/// thread, inside [`Read::read`].
#[derive(Clone, Default)]
struct SharedBuffer(Arc<Mutex<Vec<u8>>>);

impl SharedBuffer {
    /// Whether the encoder has yet to emit any bytes.
    fn is_empty(&self) -> bool {
        self.0.lock().expect("SAM transcode buffer mutex poisoned").is_empty()
    }

    /// Take everything written so far, leaving the buffer empty.
    fn take(&self) -> Vec<u8> {
        let mut guard = self.0.lock().expect("SAM transcode buffer mutex poisoned");
        std::mem::take(&mut *guard)
    }
}

impl Write for SharedBuffer {
    fn write(&mut self, buf: &[u8]) -> io::Result<usize> {
        let mut guard = self.0.lock().expect("SAM transcode buffer mutex poisoned");
        guard.extend_from_slice(buf);
        Ok(buf.len())
    }

    fn flush(&mut self) -> io::Result<()> {
        Ok(())
    }
}

/// Where the transcode has got to.
enum Stage {
    /// The BAM header has not been emitted yet.
    Header,
    /// Records are being transcoded.
    Records,
    /// The input is exhausted and the BGZF EOF block has been written.
    Done,
}

/// Presents an uncompressed SAM stream as the BGZF-compressed BAM byte stream
/// that fgumi's readers consume.
///
/// Records are pulled from the SAM parser only as the consumer reads, so a
/// SAM file is streamed rather than materialized.
pub struct SamToBamStream<R> {
    reader: sam::io::Reader<BufReader<R>>,
    header: sam::Header,
    writer: noodles::bam::io::Writer<noodles_bgzf::io::Writer<SharedBuffer>>,
    sink: SharedBuffer,
    /// Encoded bytes not yet handed to the consumer.
    staged: Vec<u8>,
    /// Read cursor into `staged`.
    offset: usize,
    stage: Stage,
}

impl<R: Read> SamToBamStream<R> {
    /// Parse the SAM header from `inner` and prepare the transcode.
    ///
    /// # Errors
    ///
    /// Returns an error if the SAM header cannot be parsed.
    pub fn new(inner: R) -> Result<Self> {
        let mut reader =
            sam::io::Reader::new(BufReader::with_capacity(SAM_READ_BUFFER_SIZE, inner));
        let header = reader.read_header().context("Failed to parse SAM header")?;

        let sink = SharedBuffer::default();
        let bgzf = noodles_bgzf::io::writer::Builder::default()
            .set_compression_level(CompressionLevel::NONE)
            .build_from_writer(sink.clone());

        Ok(Self {
            reader,
            header,
            writer: noodles::bam::io::Writer::from(bgzf),
            sink,
            staged: Vec::new(),
            offset: 0,
            stage: Stage::Header,
        })
    }

    /// Advance the transcode by one unit of work and stage whatever the BGZF
    /// encoder emitted.
    fn produce(&mut self) -> io::Result<()> {
        match self.stage {
            Stage::Header => {
                self.writer.write_header(&self.header)?;
                self.stage = Stage::Records;
            }
            Stage::Records => {
                let mut record = sam::Record::default();
                for _ in 0..RECORDS_PER_BATCH {
                    if self.reader.read_record(&mut record)? == 0 {
                        // try_finish flushes the trailing partial block and
                        // writes the BGZF EOF block, so the transcoded stream
                        // is a complete BAM rather than one that trips
                        // truncated-file checks.
                        self.writer.try_finish()?;
                        self.stage = Stage::Done;
                        break;
                    }
                    self.writer.write_alignment_record(&self.header, &record)?;
                    // Once a block has been emitted there are bytes to hand
                    // back, and transcoding further records would only stage
                    // output the consumer has not asked for.
                    if !self.sink.is_empty() {
                        break;
                    }
                }
            }
            Stage::Done => {}
        }

        let mut encoded = self.sink.take();
        if self.offset >= self.staged.len() {
            self.staged.clear();
            self.offset = 0;
        }
        self.staged.append(&mut encoded);

        Ok(())
    }
}

impl<R: Read> Read for SamToBamStream<R> {
    fn read(&mut self, buf: &mut [u8]) -> io::Result<usize> {
        while self.offset >= self.staged.len() {
            if matches!(self.stage, Stage::Done) {
                return Ok(0);
            }
            self.produce()?;
        }

        let n = (&self.staged[self.offset..]).read(buf)?;
        self.offset += n;
        Ok(n)
    }
}

/// Consume up to [`GZIP_MAGIC`]`.len()` bytes from `reader`.
///
/// Returns fewer bytes only at end of input. The bytes are returned rather
/// than pushed back so the caller can replay them via [`ChainedReader`] —
/// `reader` may be a pipe, which cannot be rewound or reopened.
fn read_magic<R: Read>(reader: &mut R) -> io::Result<Vec<u8>> {
    let mut magic = [0u8; GZIP_MAGIC.len()];
    let mut filled = 0;

    while filled < magic.len() {
        match reader.read(&mut magic[filled..])? {
            0 => break,
            n => filled += n,
        }
    }

    Ok(magic[..filled].to_vec())
}

/// Normalize `reader` to a BGZF byte stream, transcoding it if it is SAM text.
///
/// `path` is used only to describe the input in errors.
///
/// The magic bytes are consumed from `reader` and replayed, so this is safe
/// for non-seekable inputs (stdin, FIFOs) and never reopens `path`.
///
/// # Errors
///
/// Returns an error if the input is empty or its SAM header cannot be parsed.
pub fn normalize_to_bgzf(
    mut reader: Box<dyn Read + Send>,
    path: &Path,
) -> Result<Box<dyn Read + Send>> {
    let magic = read_magic(&mut reader)
        .with_context(|| format!("Failed to read from: {}", path.display()))?;

    if magic.is_empty() {
        bail!("Input is empty: {} (expected BAM or SAM records)", path.display());
    }

    if magic.starts_with(&GZIP_MAGIC) {
        return Ok(Box::new(ChainedReader::new(magic, reader)));
    }

    let sam = SamToBamStream::new(ChainedReader::new(magic, reader))
        .with_context(|| format!("Failed to read SAM input: {}", path.display()))?;
    Ok(Box::new(sam))
}

#[cfg(test)]
mod tests {
    use super::*;
    use std::fmt::Write as _;

    const SAM_TEXT: &str = "@HD\tVN:1.6\tSO:unsorted\n\
                            @SQ\tSN:chr1\tLN:100\n\
                            r1\t0\tchr1\t1\t60\t4M\t*\t0\t0\tACGT\tIIII\n";

    /// The transcoded bytes must be a complete BAM stream — magic, header and
    /// records — so that consumers that only know how to read BAM are
    /// unaffected by the input having been SAM.
    #[test]
    fn transcodes_sam_text_into_a_readable_bam_stream() -> Result<()> {
        let stream = SamToBamStream::new(io::Cursor::new(SAM_TEXT))?;

        let mut reader = noodles::bam::io::Reader::new(BufReader::new(stream));
        let header = reader.read_header()?;
        let records = reader.records().collect::<io::Result<Vec<_>>>()?;

        let refs = header.reference_sequences();
        assert_eq!(refs.len(), 1, "header @SQ lost");
        let (name, seq) = refs.get_index(0).expect("one @SQ");
        assert_eq!(AsRef::<[u8]>::as_ref(name), b"chr1", "@SQ SN");
        assert_eq!(usize::from(seq.length()), 100, "@SQ LN");

        // Field-for-field identity, not just count: prove the SAM record survived
        // the transcode intact rather than that *some* record was emitted.
        assert_eq!(records.len(), 1, "record count");
        let r = &records[0];
        assert_eq!(r.name().map(std::convert::AsRef::as_ref), Some(b"r1" as &[u8]), "QNAME");
        assert_eq!(r.flags().bits(), 0, "FLAG");
        assert_eq!(r.alignment_start().expect("mapped").expect("pos").get(), 1, "POS");
        let cigar: Vec<_> = r
            .cigar()
            .iter()
            .map(|op| op.map(|o| (o.kind(), o.len())))
            .collect::<io::Result<_>>()
            .expect("cigar");
        assert_eq!(cigar, [(noodles::sam::alignment::record::cigar::op::Kind::Match, 4)], "CIGAR");
        assert_eq!(r.sequence().iter().collect::<Vec<u8>>(), b"ACGT", "SEQ");
        assert_eq!(r.quality_scores().as_ref(), &[40, 40, 40, 40], "QUAL");
        Ok(())
    }

    /// A one-byte-at-a-time consumer must see the same bytes as a bulk one;
    /// the staging buffer is refilled across `read` calls.
    #[test]
    fn byte_at_a_time_reads_match_bulk_reads() -> Result<()> {
        let mut bulk = Vec::new();
        SamToBamStream::new(io::Cursor::new(SAM_TEXT))?.read_to_end(&mut bulk)?;

        let mut stream = SamToBamStream::new(io::Cursor::new(SAM_TEXT))?;
        let mut byte_at_a_time = Vec::new();
        let mut byte = [0u8; 1];
        while stream.read(&mut byte)? == 1 {
            byte_at_a_time.push(byte[0]);
        }

        assert_eq!(byte_at_a_time, bulk);
        Ok(())
    }

    /// Long-read SAM, where one record can be hundreds of kilobytes: a batch
    /// of them is far larger than a BGZF block.
    fn long_read_sam_text(records: usize, read_length: usize) -> String {
        let bases = "ACGT".repeat(read_length / 4);
        let quals = "I".repeat(bases.len());
        let mut text = format!("@HD\tVN:1.6\tSO:unsorted\n@SQ\tSN:chr1\tLN:{}\n", read_length * 2);
        for i in 0..records {
            writeln!(text, "r{i}\t4\t*\t0\t0\t*\t*\t0\t0\t{bases}\t{quals}")
                .expect("writing to a String cannot fail");
        }
        text
    }

    /// A consumer that asks for one byte must not pay for a whole batch: the
    /// transcode stops at the first bytes the encoder emits, so what is staged
    /// is bounded by the BGZF block size rather than by the record size times
    /// [`RECORDS_PER_BATCH`].
    #[test]
    fn a_small_read_stages_only_the_first_block() -> Result<()> {
        const READ_LENGTH: usize = 20_000;
        let text = long_read_sam_text(RECORDS_PER_BATCH, READ_LENGTH);

        let mut stream = SamToBamStream::new(io::Cursor::new(text))?;
        let mut byte = [0u8; 1];
        assert_eq!(stream.read(&mut byte)?, 1);

        // A single record already exceeds the 64 KiB block, so one block's
        // worth is the floor; a whole batch would be ~20 MB.
        assert!(
            stream.staged.len() < 4 * READ_LENGTH,
            "one read staged {} bytes, expected roughly one BGZF block",
            stream.staged.len()
        );
        Ok(())
    }

    /// A SAM file with a header but no records is legal and must transcode to
    /// an empty-but-valid BAM rather than erroring.
    #[test]
    fn header_only_sam_transcodes_to_an_empty_bam() -> Result<()> {
        let text = "@HD\tVN:1.6\tSO:unsorted\n@SQ\tSN:chr1\tLN:100\n";
        let stream = SamToBamStream::new(io::Cursor::new(text))?;

        let mut reader = noodles::bam::io::Reader::new(BufReader::new(stream));
        let header = reader.read_header()?;
        let records = reader.records().collect::<io::Result<Vec<_>>>()?;

        assert_eq!(header.reference_sequences().len(), 1);
        assert!(records.is_empty());
        Ok(())
    }

    /// Detection is by content, so BGZF passes through untouched.
    #[test]
    fn bgzf_input_is_passed_through_unchanged() -> Result<()> {
        let mut bgzf = Vec::new();
        {
            let mut writer = noodles::bam::io::Writer::new(&mut bgzf);
            writer.write_header(&sam::Header::default())?;
            writer.try_finish()?;
        }

        let mut normalized =
            normalize_to_bgzf(Box::new(io::Cursor::new(bgzf.clone())), Path::new("in.bam"))?;
        let mut round_tripped = Vec::new();
        normalized.read_to_end(&mut round_tripped)?;

        assert_eq!(round_tripped, bgzf);
        Ok(())
    }

    /// An empty input is neither BAM nor SAM; it must be reported rather than
    /// silently read as a headerless, record-less SAM file.
    #[test]
    fn empty_input_is_rejected() {
        let Err(err) =
            normalize_to_bgzf(Box::new(io::Cursor::new(Vec::new())), Path::new("in.bam"))
        else {
            panic!("empty input should be rejected");
        };

        assert!(
            format!("{err:#}").contains("empty"),
            "error should say the input is empty, got: {err:#}"
        );
    }
}
