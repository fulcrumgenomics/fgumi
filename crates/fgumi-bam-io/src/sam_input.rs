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
//! Detection is by content, never by file extension (see [`crate::format`]), so
//! a misnamed `reads.bam` holding SAM text and a `reads.sam` holding BGZF both
//! read correctly. gzip that is *not* BGZF is named as such rather than handed
//! to the block decoder, which would fail as though the file were corrupt.

use anyhow::{Context, Result, bail};
use noodles::sam;
use noodles::sam::alignment::io::Write as _;
use noodles_bgzf::io::writer::CompressionLevel;
use std::io::{self, BufReader, Read, Write};
use std::path::Path;
use std::sync::atomic::{AtomicUsize, Ordering};
use std::sync::{Arc, Mutex};

use crate::format::{FORMAT_PREFIX_LEN, InputFormat, classify_input};
use crate::reader::ChainedReader;

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
/// lock is uncontended: the encoder and the drain both run inside
/// [`Read::read`], so whichever thread is driving the stream holds both sides.
/// Which thread that is varies — the caller's under a single-threaded reader,
/// a `MultithreadedReader` worker under `--threads > 1`, the `fgumi-prefetch`
/// thread under `--async-reader` — but there is only ever one of them.
#[derive(Clone, Default)]
struct SharedBuffer {
    bytes: Arc<Mutex<Vec<u8>>>,
    /// Mirrors `bytes.len()`, so [`SharedBuffer::is_empty`] — which
    /// `transcode_batch` calls once per record — is an atomic load rather than a
    /// mutex acquisition. An uncontended lock is only ~10-20ns, but at a billion
    /// records that is ~10-20 seconds spent asking a question the encoder
    /// already knows the answer to.
    ///
    /// Every mutation of `bytes` must update this in the same critical section;
    /// the `debug_assert!`s below fail loudly if the two ever drift.
    len: Arc<AtomicUsize>,
}

impl SharedBuffer {
    /// Whether the encoder has yet to emit any bytes.
    ///
    /// Read with `Acquire` against the `Release` stores in [`Self::write`] and
    /// [`Self::swap_into`], so a reader that observes a non-zero length also
    /// observes the bytes behind it. The answer is only a hint either way — a
    /// stale `true` transcodes one more record, a stale `false` produces an
    /// empty batch the caller's refill loop retries — but the ordering keeps it
    /// from being a data race in the handoff between reader threads.
    fn is_empty(&self) -> bool {
        self.len.load(Ordering::Acquire) == 0
    }

    /// Exchange the encoder's buffer with `other`, leaving the encoder holding
    /// `other`'s (drained) allocation.
    ///
    /// A swap rather than a take so neither buffer is reallocated: `take` would
    /// leave the sink at zero capacity, forcing it to regrow across the ~14
    /// writes the BGZF encoder issues per block.
    fn swap_into(&self, other: &mut Vec<u8>) {
        let mut guard = self.bytes.lock().expect("SAM transcode buffer mutex poisoned");
        debug_assert_eq!(guard.len(), self.len.load(Ordering::Relaxed), "sink length drifted");
        std::mem::swap(&mut *guard, other);
        self.len.store(guard.len(), Ordering::Release);
    }
}

impl Write for SharedBuffer {
    fn write(&mut self, buf: &[u8]) -> io::Result<usize> {
        let mut guard = self.bytes.lock().expect("SAM transcode buffer mutex poisoned");
        debug_assert_eq!(guard.len(), self.len.load(Ordering::Relaxed), "sink length drifted");
        guard.extend_from_slice(buf);
        self.len.store(guard.len(), Ordering::Release);
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
    /// Reused across records. `read_record` clears and refills it, so keeping it
    /// here preserves its capacity; a fresh `Record` per call would regrow the
    /// buffer every time. On long reads `produce` runs once per record, so that
    /// regrowth would be per record rather than per block.
    scratch: sam::Record,
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
            scratch: sam::Record::default(),
            stage: Stage::Header,
        })
    }

    /// Read one SAM record, labelling a parse failure as a record failure.
    ///
    /// The consumer sees this stream as a BAM byte stream and wraps whatever it
    /// gets in its own context, so without this the raw noodles message
    /// (`unexpected EOL`) would be the only clue that the fault was in a
    /// record rather than in the framing around it.
    fn read_sam_record(&mut self, record: &mut sam::Record) -> io::Result<usize> {
        self.reader
            .read_record(record)
            .map_err(|e| io::Error::new(e.kind(), format!("Failed to parse SAM record: {e}")))
    }

    /// Advance the transcode by one unit of work and stage whatever the BGZF
    /// encoder emitted.
    fn produce(&mut self) -> io::Result<()> {
        match self.stage {
            Stage::Header => {
                self.writer.write_header(&self.header)?;
                // Seal the header into its own BGZF block. Without this the
                // encoder emits nothing until a block fills, so the consumer's
                // `read_header()` keeps pulling and ends up transcoding
                // records — and a malformed record then surfaces as a *header*
                // failure, pointing whoever reads the error at the wrong end
                // of the file. Costs one short block per stream.
                self.writer.get_mut().flush()?;
                self.stage = Stage::Records;
            }
            Stage::Records => {
                // Moved out so `self` stays borrowable inside the loop; put back
                // below so the allocation survives to the next call.
                let mut record = std::mem::take(&mut self.scratch);
                let result = self.transcode_batch(&mut record);
                self.scratch = record;
                result?;
            }
            Stage::Done => {}
        }

        self.stage_encoded();

        Ok(())
    }

    /// Transcode up to [`RECORDS_PER_BATCH`] records into the BGZF encoder,
    /// stopping early once the encoder has emitted a block.
    fn transcode_batch(&mut self, record: &mut sam::Record) -> io::Result<()> {
        for _ in 0..RECORDS_PER_BATCH {
            if self.read_sam_record(record)? == 0 {
                // try_finish flushes the trailing partial block and writes the
                // BGZF EOF block, so the transcoded stream is a complete BAM
                // rather than one that trips truncated-file checks.
                self.writer.try_finish()?;
                self.stage = Stage::Done;
                break;
            }
            self.writer.write_alignment_record(&self.header, record)?;
            // Once a block has been emitted there are bytes to hand back, and
            // transcoding further records would only stage output the consumer
            // has not asked for.
            if !self.sink.is_empty() {
                break;
            }
        }
        Ok(())
    }

    /// Move whatever the encoder emitted into the staging buffer.
    ///
    /// `produce` is only ever called from `read` with the staging buffer already
    /// drained, so this swaps rather than appends: the sink hands its allocation
    /// to `staged` and inherits the drained one, which keeps both buffers warm
    /// and copies nothing.
    fn stage_encoded(&mut self) {
        debug_assert!(
            self.offset >= self.staged.len(),
            "produce() ran with {} unread staged bytes; clearing here would drop them",
            self.staged.len() - self.offset
        );
        self.staged.clear();
        self.offset = 0;
        self.sink.swap_into(&mut self.staged);
    }
}

impl<R: Read> Read for SamToBamStream<R> {
    fn read(&mut self, buf: &mut [u8]) -> io::Result<usize> {
        // A zero-length request leaves the staging buffer drained, which the
        // refill loop below would read as "transcode more" — encoding a whole
        // BGZF block to satisfy a request for nothing. Nothing is lost without
        // the guard (`offset` does not advance), but the work is wasted, and
        // `ChainedReader::read` guards the same way at the same point in the
        // same chain.
        if buf.is_empty() {
            return Ok(0);
        }

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

/// Consume up to [`FORMAT_PREFIX_LEN`] bytes from `reader`.
///
/// Returns fewer bytes only at end of input. The bytes are returned rather
/// than pushed back so the caller can replay them via [`ChainedReader`] —
/// `reader` may be a pipe, which cannot be rewound or reopened.
fn read_format_prefix<R: Read>(reader: &mut R) -> io::Result<Vec<u8>> {
    let mut prefix = [0u8; FORMAT_PREFIX_LEN];
    let filled = crate::reader::read_prefix(reader, &mut prefix)?;
    Ok(prefix[..filled].to_vec())
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
/// Returns an error if the input is empty, is gzip but not BGZF, or is text
/// whose SAM header cannot be parsed (or is absent).
pub fn normalize_to_bgzf(
    mut reader: Box<dyn Read + Send>,
    path: &Path,
) -> Result<Box<dyn Read + Send>> {
    let magic = read_format_prefix(&mut reader)
        .with_context(|| format!("Failed to read from: {}", path.display()))?;

    match classify_input(&magic) {
        InputFormat::Empty => {
            bail!("Input is empty: {} (expected BAM or SAM records)", path.display())
        }
        InputFormat::Bgzf => return Ok(Box::new(ChainedReader::new(magic, reader))),
        // Distinguished from BGZF so the error names the real problem. Left
        // undetected, a plain-gzipped file reaches the BGZF decoder and fails
        // with "Failed to read header", which reads as a corrupt BAM.
        //
        // FUTURE: support plain-gzip SAM rather than rejecting it. Now that the
        // transcode exists this is a small arm — decode with
        // `flate2::read::MultiGzDecoder` and hand the result to
        // `SamToBamStream`, landing on the same normalized stream. Two things to
        // settle first:
        //
        //   1. It costs block-parallel decode. gzip is one deflate stream, so
        //      `--threads` would silently stop helping on the read side; that
        //      wants a one-time `log::warn!` rather than a silent cliff.
        //   2. Gzipped *SAM* is a legitimate text file and worth accepting;
        //      gzipped *BAM* is a malformed artifact that samtools cannot read
        //      either, so accepting it would just move the failure downstream.
        //      The arm should therefore transcode as SAM and still reject if the
        //      decompressed bytes turn out to be BGZF/BAM.
        //
        // `extract` already accepts plain-gzip FASTQ, so the asymmetry users see
        // today is real; it is defensible only while it stays explained.
        // Phrased as "not readable as BGZF" rather than "not BGZF": the verdict also
        // covers a gzip member carrying a `BC` subfield in a layout fgumi's decoder
        // rejects, which is BGZF by the letter of RFC 1952 (see
        // `classify_input`). The remedy is the same for both, so the message does not
        // split them.
        InputFormat::Gzip => bail!(
            "Input is gzip-compressed but not readable as BGZF: {} — fgumi reads BGZF (bgzip) \
             so it can decompress blocks in parallel. Recompress with `bgzip`, or pipe it in: \
             `gzip -dc {} | fgumi ... -i -`",
            path.display(),
            path.display()
        ),
        InputFormat::Text => {}
    }

    let sam = SamToBamStream::new(ChainedReader::new(magic, reader))
        .with_context(|| format!("Failed to read SAM input: {}", path.display()))?;

    // A SAM header parse consumes only lines beginning with `@`, so arbitrary
    // text yields an empty header instead of an error. Left alone that would
    // transcode a truncated or wrong-format file into a valid-looking BAM with
    // no reference sequences at all. Before SAM was accepted this input was
    // rejected outright for not being BGZF, and it stays rejected.
    if sam.header.is_empty() {
        bail!(
            "Input is neither BAM nor SAM: {} (no SAM header found; SAM input must begin \
             with an `@` header line)",
            path.display()
        );
    }

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

    /// A zero-length read must not drive a transcode. `read` refills whenever
    /// the staging buffer is drained, and a zero-length request leaves it
    /// drained, so without a guard a consumer asking for nothing pays for a
    /// whole BGZF block. [`crate::reader::ChainedReader::read`] guards the same
    /// way at the same point in the same chain.
    #[test]
    fn a_zero_length_read_transcodes_nothing() -> Result<()> {
        let mut stream = SamToBamStream::new(io::Cursor::new(SAM_TEXT))?;

        assert_eq!(stream.read(&mut [])?, 0, "zero-length read");
        assert!(
            stream.staged.is_empty(),
            "a zero-length read staged {} bytes",
            stream.staged.len()
        );
        assert!(matches!(stream.stage, Stage::Header), "a zero-length read advanced the transcode");

        // Returning early must not cost the stream anything: the same reader
        // still delivers the whole BAM afterwards.
        let mut after_guard = Vec::new();
        stream.read_to_end(&mut after_guard)?;
        let mut whole = Vec::new();
        SamToBamStream::new(io::Cursor::new(SAM_TEXT))?.read_to_end(&mut whole)?;
        assert_eq!(after_guard, whole, "bytes read after a zero-length read");
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

    /// The header must be readable without pulling a single record through the
    /// transcoder. Otherwise a malformed *record* surfaces inside the
    /// consumer's `read_header()` and gets reported as a header failure —
    /// sending whoever reads the error to the wrong end of the file.
    #[test]
    fn a_malformed_record_does_not_fail_the_header_read() -> Result<()> {
        let text = "@HD\tVN:1.6\tSO:unsorted\n\
                    @SQ\tSN:chr1\tLN:100\n\
                    r1\t0\tchr1\t1\t60\t4M\t*\t0\t0\tACGT\tIIII\n\
                    NOT_A_RECORD\n";

        let stream = SamToBamStream::new(io::Cursor::new(text))?;
        let mut reader = noodles::bam::io::Reader::new(BufReader::new(stream));

        let header = reader.read_header().expect("header must parse despite the bad record");
        assert_eq!(header.reference_sequences().len(), 1);

        let err = reader
            .records()
            .collect::<io::Result<Vec<_>>>()
            .expect_err("the malformed record must still be an error");
        assert!(
            format!("{err:#}").contains("SAM record"),
            "error should name the record, got: {err:#}"
        );
        Ok(())
    }

    /// Text that is neither BGZF nor SAM parses to an *empty* SAM header, which
    /// would otherwise transcode into a BAM with no reference sequences at all.
    /// Before SAM input was accepted this was a hard "not BGZF" rejection, and
    /// it must stay a rejection.
    #[test]
    fn non_sam_text_is_rejected() {
        let Err(err) = normalize_to_bgzf(
            Box::new(io::Cursor::new(b"this is not a SAM file\n".to_vec())),
            Path::new("in.txt"),
        ) else {
            panic!("non-SAM text should be rejected");
        };

        let msg = format!("{err:#}");
        assert!(
            msg.contains("neither BAM nor SAM"),
            "error should say the input is not BAM or SAM, got: {msg}"
        );
    }

    /// The rejection above keys on the *header*, not on the record parse, so a
    /// SAM file whose header is present still transcodes even if it declares no
    /// reference sequences (a legal unmapped-only SAM).
    #[test]
    fn sam_with_a_header_but_no_reference_sequences_is_accepted() -> Result<()> {
        let text = "@HD\tVN:1.6\tSO:unsorted\n\
                    r1\t4\t*\t0\t0\t*\t*\t0\t0\tACGT\tIIII\n";

        let normalized = normalize_to_bgzf(
            Box::new(io::Cursor::new(text.as_bytes().to_vec())),
            Path::new("in.sam"),
        )?;
        let mut reader = noodles::bam::io::Reader::new(BufReader::new(normalized));
        let header = reader.read_header()?;
        let records = reader.records().collect::<io::Result<Vec<_>>>()?;

        assert!(header.reference_sequences().is_empty(), "no @SQ expected");
        assert_eq!(records.len(), 1, "the unmapped record must survive");
        Ok(())
    }
}
