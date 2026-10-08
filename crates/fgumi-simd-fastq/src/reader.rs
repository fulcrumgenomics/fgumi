//! Buffered FASTQ reader using SIMD-accelerated parsing.
//!
//! [`SimdFastqReader`] wraps a `BufRead` and yields owned FASTQ records,
//! serving as a drop-in replacement for `seq_io::fastq::Reader`.

use std::io::{self, BufRead};

use crate::parser::{self, try_parse_single_record};

/// Default internal buffer size (1 MiB), matching `seq_io`'s default.
const DEFAULT_BUFFER_SIZE: usize = 1 << 20;

/// An owned FASTQ record with heap-allocated name, sequence, and quality.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct OwnedFastqRecord {
    /// Read name (without leading `@`).
    pub name: Vec<u8>,
    /// Sequence bases.
    pub sequence: Vec<u8>,
    /// Quality scores (Phred-encoded ASCII).
    pub quality: Vec<u8>,
}

/// Buffered FASTQ reader that uses SIMD-accelerated record boundary detection.
///
/// Reads chunks from the underlying `BufRead`, finds record boundaries with
/// [`find_record_offsets`](crate::find_record_offsets), and yields owned records.
///
/// Errors do not stall iteration. A malformed record, including a truncated or
/// unparseable tail at EOF, is consumed when its `InvalidData` error is returned,
/// so the next call moves on (to `None` after a bad tail). An I/O error from the
/// source keeps every byte already buffered, so a caller that retries resumes the
/// stream where it stopped.
///
/// # Example
///
/// ```
/// use fgumi_simd_fastq::SimdFastqReader;
/// use std::io::Cursor;
///
/// let data = b"@r1\nACGT\n+\nIIII\n@r2\nTTTT\n+\nJJJJ\n";
/// let mut reader = SimdFastqReader::new(Cursor::new(&data[..]));
///
/// let rec = reader.next().unwrap().unwrap();
/// assert_eq!(rec.name, b"r1");
/// assert_eq!(rec.sequence, b"ACGT");
/// ```
pub struct SimdFastqReader<R: BufRead> {
    inner: R,
    /// Internal buffer holding data read from the source.
    buffer: Vec<u8>,
    /// Pre-computed record boundary offsets within `buffer[..valid]`.
    offsets: Vec<usize>,
    /// Index into `offsets` for the next record to yield.
    next_record_idx: usize,
    /// Number of valid bytes in `buffer`.
    valid: usize,
    /// True when the underlying reader has returned 0 bytes.
    at_eof: bool,
}

impl<R: BufRead> SimdFastqReader<R> {
    /// Create a new reader with the default buffer size (1 MiB).
    pub fn new(inner: R) -> Self {
        Self::with_capacity(inner, DEFAULT_BUFFER_SIZE)
    }

    /// Create a new reader with a custom buffer capacity.
    pub fn with_capacity(inner: R, capacity: usize) -> Self {
        Self {
            inner,
            buffer: Vec::with_capacity(capacity),
            offsets: Vec::new(),
            next_record_idx: 0,
            valid: 0,
            at_eof: false,
        }
    }

    /// Consume the reader and return the underlying source.
    ///
    /// Any bytes this reader pulled from the source but has not yielded as
    /// records are **discarded** — they live in the internal buffer, not in the
    /// source. A caller that needs the whole stream afterwards must therefore
    /// have been capturing the source's bytes independently (e.g. through a tee)
    /// rather than relying on the source's own position.
    pub fn into_inner(self) -> R {
        self.inner
    }

    /// Fill the internal buffer, preserving any leftover bytes from incomplete records.
    ///
    /// Returns `true` if the buffer holds any bytes for `next()` to consume: complete
    /// records, or a final leftover at EOF that is yielded or reported as an error.
    fn fill_buffer(&mut self) -> io::Result<bool> {
        // Determine leftover: bytes from the last complete record offset to end of valid data.
        let leftover_start = if self.offsets.is_empty() {
            0
        } else {
            // The last offset in `offsets` is the start of the first incomplete record
            // (or the end of the last complete record, which is the same thing).
            self.offsets.last().copied().unwrap_or(0)
        };

        // Move leftover bytes to the front of the buffer
        if leftover_start > 0 && leftover_start < self.valid {
            self.buffer.copy_within(leftover_start..self.valid, 0);
            self.valid -= leftover_start;
        } else if leftover_start >= self.valid {
            self.valid = 0;
        }
        // The old offsets index the pre-compaction buffer. Drop them now, so that if the
        // read below fails and the caller retries, the next fill treats all of
        // `buffer[..valid]` as leftover instead of re-applying stale offsets to it.
        self.offsets.clear();
        self.next_record_idx = 0;

        // If the buffer is full of leftover (no complete records found), grow it
        // so we can read more data and find the end of the current record.
        if self.valid >= self.buffer.capacity() {
            self.buffer.reserve(self.buffer.capacity().max(4096));
        }
        self.buffer.resize(self.buffer.capacity(), 0);

        // Read new data into the buffer after the leftover
        let mut total_read = 0;
        while self.valid + total_read < self.buffer.len() {
            let buf = &mut self.buffer[self.valid + total_read..];
            if buf.is_empty() {
                break;
            }
            match self.inner.read(buf) {
                Ok(0) => {
                    self.at_eof = true;
                    break;
                }
                Ok(n) => total_read += n,
                Err(e) if e.kind() == io::ErrorKind::Interrupted => {}
                Err(e) => {
                    // Keep the bytes this fill already read, so a retry continues the
                    // stream rather than skipping them.
                    self.valid += total_read;
                    self.buffer.truncate(self.valid);
                    return Err(e);
                }
            }
        }

        self.valid += total_read;
        self.buffer.truncate(self.valid);

        // Find record boundaries in the buffer
        self.offsets = parser::find_record_offsets(&self.buffer[..self.valid]);
        self.next_record_idx = 0;

        // Report data whenever there are complete records OR any buffered bytes. At EOF a
        // non-empty buffer with no complete record is a final leftover, which `next()`
        // must still see: it is either a complete-but-unterminated record (yielded) or a
        // truncated one (an error). Returning `false` there would end iteration and drop
        // that leftover silently, e.g. a lone record with no trailing newline, or a final
        // record whose bytes all arrive in the last fill.
        Ok(self.offsets.len() > 1 || self.valid > 0)
    }
}

impl<R: BufRead> Iterator for SimdFastqReader<R> {
    type Item = io::Result<OwnedFastqRecord>;

    fn next(&mut self) -> Option<Self::Item> {
        loop {
            if self.next_record_idx + 1 < self.offsets.len() {
                let start = self.offsets[self.next_record_idx];
                let end = self.offsets[self.next_record_idx + 1];
                self.next_record_idx += 1;

                let borrowed = match try_parse_single_record(&self.buffer[start..end]) {
                    Ok(rec) => rec,
                    Err(e) => {
                        // Preserve the typed `FastqParseError` as the io::Error's inner error
                        // so downstream callers can inspect or downcast it via `get_ref()`.
                        return Some(Err(io::Error::new(io::ErrorKind::InvalidData, e)));
                    }
                };
                return Some(Ok(OwnedFastqRecord {
                    name: borrowed.name.to_vec(),
                    sequence: borrowed.sequence.to_vec(),
                    quality: borrowed.quality.to_vec(),
                }));
            }

            if self.at_eof {
                // Leftover bytes after the last complete-record boundary. A final
                // record with no trailing newline is still *complete* — it has the
                // three internal newlines after name/seq/`+`, only the fourth
                // (terminating) newline is absent — so yield it rather than
                // erroring (fgbio's `lines.take(4)` accepts an unterminated final
                // record). Only a leftover with fewer than three newlines is
                // genuinely truncated.
                let leftover_start = self.offsets.last().copied().unwrap_or(0);
                if leftover_start < self.valid {
                    let leftover = &self.buffer[leftover_start..self.valid];
                    // A complete-but-unterminated final record has exactly the three
                    // internal newlines after name/seq/`+`; cap the scan at four so a
                    // longer malformed run still fails the `== 3` check. It must also
                    // NOT end in a newline: a genuine unterminated record ends with
                    // its quality bytes, whereas a leftover like `@n\nseq\n+\n` has a
                    // terminated `+` line and an *absent* quality line — three
                    // newlines but no quality — which must error, not be parsed (its
                    // empty quality span is not a valid record).
                    let internal_newlines =
                        leftover.iter().filter(|&&b| b == b'\n').take(4).count();
                    if leftover.first() == Some(&b'@')
                        && leftover.last() != Some(&b'\n')
                        && internal_newlines == 3
                    {
                        // Route through the fallible parser like the main record path
                        // above: the coarse guard admits the leftover, but a typed
                        // `FastqParseError` (e.g. length mismatch) still surfaces as a
                        // clean io::Error rather than a panic.
                        let borrowed = match try_parse_single_record(leftover) {
                            Ok(rec) => rec,
                            Err(e) => {
                                // Consume the leftover so a caller that continues past
                                // the error terminates rather than seeing it repeatedly.
                                self.valid = leftover_start;
                                return Some(Err(io::Error::new(io::ErrorKind::InvalidData, e)));
                            }
                        };
                        let record = OwnedFastqRecord {
                            name: borrowed.name.to_vec(),
                            sequence: borrowed.sequence.to_vec(),
                            quality: borrowed.quality.to_vec(),
                        };
                        // Mark the leftover consumed so the next call terminates.
                        self.valid = leftover_start;
                        return Some(Ok(record));
                    }
                    let leftover_len = self.valid - leftover_start;
                    // Consume the leftover so a caller that continues past the error
                    // terminates rather than seeing it repeatedly.
                    self.valid = leftover_start;
                    return Some(Err(io::Error::new(
                        io::ErrorKind::InvalidData,
                        format!("Truncated FASTQ record at EOF ({leftover_len} leftover bytes)"),
                    )));
                }
                return None;
            }

            match self.fill_buffer() {
                Ok(true) => {}
                Ok(false) => return None,
                Err(e) => return Some(Err(e)),
            }
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use rstest::rstest;
    use std::io::Cursor;

    #[test]
    fn test_reader_single_record() {
        let data = b"@r1\nACGT\n+\nIIII\n";
        let mut reader = SimdFastqReader::new(Cursor::new(&data[..]));

        let rec = reader
            .next()
            .expect("reader should yield a record")
            .expect("record should parse successfully");
        assert_eq!(rec.name, b"r1");
        assert_eq!(rec.sequence, b"ACGT");
        assert_eq!(rec.quality, b"IIII");

        assert!(reader.next().is_none());
    }

    #[test]
    fn test_reader_multiple_records() {
        let data = b"@r1\nACGT\n+\nIIII\n@r2\nTTTT\n+\nJJJJ\n";
        let mut reader = SimdFastqReader::new(Cursor::new(&data[..]));

        let rec1 = reader
            .next()
            .expect("reader should yield a record")
            .expect("record should parse successfully");
        assert_eq!(rec1.name, b"r1");

        let rec2 = reader
            .next()
            .expect("reader should yield a record")
            .expect("record should parse successfully");
        assert_eq!(rec2.name, b"r2");

        assert!(reader.next().is_none());
    }

    /// EXT3-07: a final record with no trailing newline is complete (it has the
    /// three internal newlines after name/seq/`+`); the reader must yield it
    /// rather than erroring "Truncated FASTQ record at EOF". fgbio's
    /// `lines.take(4)` likewise accepts an unterminated final record.
    #[test]
    fn reader_accepts_final_record_without_trailing_newline() {
        let data = b"@r1\nACGT\n+\nIIII\n@r2\nTTTT\n+\nJJJJ"; // note: no trailing \n
        let mut reader = SimdFastqReader::new(Cursor::new(&data[..]));

        let rec1 = reader.next().expect("first record").expect("first record parses");
        assert_eq!(rec1.name, b"r1");
        let rec2 = reader
            .next()
            .expect("unterminated final record must still be yielded")
            .expect("unterminated final record must parse, not error");
        assert_eq!(rec2.name, b"r2");
        assert_eq!(rec2.sequence, b"TTTT");
        assert_eq!(rec2.quality, b"JJJJ");
        assert!(reader.next().is_none(), "reader must terminate cleanly after the last record");
    }

    /// A malformed final record must surface an `InvalidData` error rather than
    /// being silently dropped or mis-parsed, across every shape the coarse
    /// end-of-stream gate can admit:
    /// * `missing_quality_line` — genuinely truncated (fewer than the three
    ///   internal newlines); the `+`/quality lines are absent entirely.
    /// * `newline_terminated_plus_absent_quality` — the `+` line IS
    ///   newline-terminated but the quality line is then entirely absent
    ///   (`@n\nseq\n+\n`): exactly three internal newlines and a leading `@`, yet
    ///   NOT a complete record (parsing it would compute `qual_start > qual_end`
    ///   and panic). A genuine unterminated final record ends with its quality
    ///   bytes, never a newline.
    /// * `length_mismatch_via_eof` — passes the coarse gate (leading `@`, three
    ///   internal newlines, no trailing `\n`) but the final record's quality is
    ///   shorter than its sequence, so the coarse-gate→`try_parse_single_record`
    ///   propagation path must reject it.
    #[rstest]
    #[case::missing_quality_line(b"@r1\nACGT\n+\nIIII\n@r2\nTTTT\n+")]
    #[case::newline_terminated_plus_absent_quality(b"@r1\nACGT\n+\nIIII\n@r2\nTTTT\n+\n")]
    #[case::length_mismatch_via_eof(b"@r1\nACGT\n+\nIIII\n@r2\nTTTT\n+\nJJ")]
    fn reader_errors_on_malformed_final_record(#[case] data: &[u8]) {
        let mut reader = SimdFastqReader::new(Cursor::new(data));
        let _ = reader.next().expect("first record").expect("first record parses");
        assert_next_is_invalid_data(&mut reader);
        // The erroring leftover must be consumed, so a caller that continues past the
        // error terminates instead of receiving the same error forever.
        assert!(reader.next().is_none(), "reader must terminate after the EOF-leftover error");
    }

    #[test]
    fn test_reader_tiny_buffer() {
        // Use a very small buffer to force multiple refills
        let data = b"@r1\nACGT\n+\nIIII\n@r2\nTTTT\n+\nJJJJ\n";
        let mut reader = SimdFastqReader::with_capacity(Cursor::new(&data[..]), 20);

        let rec1 = reader
            .next()
            .expect("reader should yield a record")
            .expect("record should parse successfully");
        assert_eq!(rec1.name, b"r1");
        assert_eq!(rec1.sequence, b"ACGT");

        let rec2 = reader
            .next()
            .expect("reader should yield a record")
            .expect("record should parse successfully");
        assert_eq!(rec2.name, b"r2");
        assert_eq!(rec2.sequence, b"TTTT");

        assert!(reader.next().is_none());
    }

    #[test]
    fn test_reader_empty_input() {
        let data = b"";
        let mut reader = SimdFastqReader::new(Cursor::new(&data[..]));
        assert!(reader.next().is_none());
    }

    /// A final record with no trailing newline must be yielded even when no complete
    /// (newline-terminated) record shares its buffer fill: either it is the only record
    /// in the input (`lone`), or a buffer exactly one record wide leaves it alone in the
    /// last fill (`after_buffer_boundary`, 16 bytes = `@r1\nACGT\n+\nIIII\n`). Both
    /// were previously dropped silently, because `fill_buffer` reported no data at EOF.
    #[rstest]
    #[case::lone(b"@r2\nTTTT\n+\nJJJJ", DEFAULT_BUFFER_SIZE, &[])]
    #[case::after_buffer_boundary(b"@r1\nACGT\n+\nIIII\n@r2\nTTTT\n+\nJJJJ", 16, &["r1"])]
    fn reader_yields_unterminated_final_record_alone_in_last_fill(
        #[case] data: &[u8],
        #[case] capacity: usize,
        #[case] leading_names: &[&str],
    ) {
        let mut reader = SimdFastqReader::with_capacity(Cursor::new(data), capacity);
        for &name in leading_names {
            let rec = reader.next().expect("leading record").expect("leading record parses");
            assert_eq!(rec.name, name.as_bytes());
        }
        let rec = reader
            .next()
            .expect("unterminated final record must be yielded, not dropped")
            .expect("unterminated final record must parse");
        assert_eq!(
            rec,
            OwnedFastqRecord {
                name: b"r2".to_vec(),
                sequence: b"TTTT".to_vec(),
                quality: b"JJJJ".to_vec(),
            }
        );
        assert!(reader.next().is_none(), "reader must terminate after the last record");
    }

    /// A truncated leftover that is alone in the last buffer fill must be an error, not
    /// silently dropped: the `lone` input has no complete record at all, and the
    /// `after_buffer_boundary` input's truncated tail arrives on its own in the final fill.
    #[rstest]
    #[case::lone(b"@r2\nTTTT\n+", DEFAULT_BUFFER_SIZE, &[], 10)]
    #[case::after_buffer_boundary(b"@r1\nACGT\n+\nIIII\n@r2\nTTTT\n+", 16, &["r1"], 10)]
    fn reader_errors_on_truncated_leftover_alone_in_last_fill(
        #[case] data: &[u8],
        #[case] capacity: usize,
        #[case] leading_names: &[&str],
        #[case] leftover_bytes: usize,
    ) {
        let mut reader = SimdFastqReader::with_capacity(Cursor::new(data), capacity);
        for &name in leading_names {
            let rec = reader.next().expect("leading record").expect("leading record parses");
            assert_eq!(rec.name, name.as_bytes());
        }
        let err = reader
            .next()
            .expect("a truncated leftover must surface, not end iteration")
            .expect_err("a truncated leftover must be an error");
        assert_eq!(err.kind(), io::ErrorKind::InvalidData);
        assert_eq!(
            err.to_string(),
            format!("Truncated FASTQ record at EOF ({leftover_bytes} leftover bytes)")
        );
        assert!(
            reader.next().is_none(),
            "reader must terminate after the truncated-leftover error"
        );
    }

    // ── Malformed-record validation (mirrors fgbio's FastqSource throwing
    // behavior in FastqIoTest, adapted to fgumi's fallible `io::Result` item
    // contract: the streaming reader surfaces an `InvalidData` error rather than
    // panicking or silently yielding a corrupt record). ──

    fn assert_next_is_invalid_data<R: BufRead>(reader: &mut SimdFastqReader<R>) {
        let err = reader
            .next()
            .expect("reader should yield an item")
            .expect_err("malformed record should be an error, not Ok");
        assert_eq!(err.kind(), io::ErrorKind::InvalidData, "unexpected error: {err}");
    }

    #[rstest]
    // Header line does not start with '@'.
    #[case::missing_at_prefix(b"r1\nACGT\n+\nIIII\n")]
    // Third line (quality header) does not start with '+'.
    #[case::missing_plus_line(b"@r1\nACGT\n?\nIIII\n")]
    // Sequence has 4 bases but quality has 2 scores.
    #[case::seq_qual_length_mismatch(b"@r1\nACGT\n+\nII\n")]
    fn test_reader_rejects_malformed_record(#[case] data: &[u8]) {
        let mut reader = SimdFastqReader::new(Cursor::new(data));
        assert_next_is_invalid_data(&mut reader);
        assert!(reader.next().is_none(), "reader must terminate after the malformed record");
    }

    #[test]
    fn test_reader_error_preserves_parse_error_source() {
        // Sequence has 4 bases but quality has 2 scores.
        let data = b"@r1\nACGT\n+\nII\n";
        let mut reader = SimdFastqReader::new(Cursor::new(&data[..]));
        let err = reader
            .next()
            .expect("reader should yield an item")
            .expect_err("malformed record should be an error, not Ok");

        // The typed FastqParseError must survive as the io::Error's inner error so callers
        // can downcast it, rather than being flattened into a string.
        let inner = err.get_ref().expect("io::Error should carry an inner error");
        let parse_err = inner
            .downcast_ref::<parser::FastqParseError>()
            .expect("inner error should downcast to FastqParseError");
        assert_eq!(*parse_err, parser::FastqParseError::LengthMismatch { seq_len: 4, qual_len: 2 });
    }

    #[test]
    fn test_reader_long_records() {
        // A record four times the reader's capacity, so it can only be produced by
        // refilling the buffer mid-record. The bases and qualities vary with position
        // rather than repeating one byte: a misplaced boundary that dropped, duplicated,
        // or reordered a chunk would still yield 500 bytes, so only comparing the exact
        // bytes back can catch it.
        let seq: String = (0..500).map(|i| ['A', 'C', 'G', 'T'][i % 4]).collect();
        let qual: String =
            (0..500).map(|i| char::from(b'!' + u8::try_from(i % 60).unwrap())).collect();
        let data = format!("@longread\n{seq}\n+\n{qual}\n");
        let mut reader = SimdFastqReader::with_capacity(Cursor::new(data.as_bytes()), 256);

        let rec = reader
            .next()
            .expect("reader should yield a record")
            .expect("record should parse successfully");
        assert_eq!(rec.name, b"longread");
        assert_eq!(rec.sequence, seq.as_bytes(), "sequence mismatch across buffer refills");
        assert_eq!(rec.quality, qual.as_bytes(), "quality mismatch across buffer refills");

        assert!(reader.next().is_none());
    }
    /// A source that replays a fixed script of reads: each `Ok` chunk is returned whole
    /// by one `read` call, each `Err` is returned once, and the source is at EOF after.
    struct ScriptedReader(std::collections::VecDeque<io::Result<Vec<u8>>>);

    impl io::Read for ScriptedReader {
        fn read(&mut self, buf: &mut [u8]) -> io::Result<usize> {
            match self.0.pop_front() {
                None => Ok(0),
                Some(Err(e)) => Err(e),
                Some(Ok(chunk)) => {
                    assert!(chunk.len() <= buf.len(), "test chunk larger than the read buffer");
                    buf[..chunk.len()].copy_from_slice(&chunk);
                    Ok(chunk.len())
                }
            }
        }
    }

    /// An I/O error from the source must not lose or corrupt buffered bytes: a caller that
    /// retries after the error must still receive every record, in order. `partial_read`
    /// fails after a short read inside one fill (those bytes must be kept);
    /// `after_leftover_compaction` fails on the fill that follows a buffer ending
    /// mid-record (the compacted leftover `@r2\n` must survive the retry).
    #[rstest]
    #[case::partial_read(DEFAULT_BUFFER_SIZE)]
    #[case::after_leftover_compaction(20)]
    fn reader_recovers_after_source_io_error(#[case] capacity: usize) {
        let data: &[u8] = b"@r1\nACGT\n+\nIIII\n@r2\nTTTT\n+\nJJJJ\n";
        let script = vec![
            Ok(data[..20].to_vec()),
            Err(io::Error::other("transient read failure")),
            Ok(data[20..].to_vec()),
        ];
        let source = io::BufReader::with_capacity(1, ScriptedReader(script.into()));
        let reader = SimdFastqReader::with_capacity(source, capacity);

        let mut names = Vec::new();
        let mut errors = 0;
        for item in reader {
            match item {
                Ok(rec) => names.push(String::from_utf8(rec.name).expect("ASCII name")),
                Err(e) => {
                    assert_eq!(e.kind(), io::ErrorKind::Other, "unexpected error: {e}");
                    errors += 1;
                    assert!(errors <= 1, "the transient error must surface exactly once");
                }
            }
        }
        assert_eq!(errors, 1, "the transient error must surface");
        assert_eq!(names, ["r1", "r2"]);
    }
}
