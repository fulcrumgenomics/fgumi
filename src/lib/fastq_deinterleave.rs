//! De-interleave a single interleaved paired-end FASTQ stream into two logical
//! byte streams (R1 and R2).
//!
//! An interleaved FASTQ stream carries the two reads of each template back to
//! back — `R1, R2, R1, R2, …`. Some streaming trimmers and aligners emit this
//! shape so a single pipe can carry a pair (e.g. `bwa mem -p`). [`deinterleave`]
//! splits such a stream into two [`BufRead`] halves so downstream code that
//! expects two separate FASTQ inputs (like `fgumi extract`'s paired path) can
//! consume it unchanged.
//!
//! The split is done on **raw record bytes**, not parsed records: each 4-line
//! FASTQ record is routed whole to its side, so the read name — including any
//! trailing comment that carries a UMI — is preserved verbatim. Only 4-line
//! FASTQ is supported (one line each for header, sequence, `+`, and quality),
//! which every modern sequencer and every tool that would pipe into extract
//! emits.
//!
//! The de-interleaver itself does **not** validate record structure — it only
//! counts lines in groups of four. Its safety against malformed input rests on
//! matching the framing of the downstream FASTQ reader, which also delimits
//! records every fourth newline: a multi-line (wrapped) record is mis-framed
//! here into pieces the downstream parser then rejects (a piece whose third line
//! is not `+`, or a truncated final piece), rather than silently completing.
//! Pairing this splitter with a laxer, `@`-scanning consumer would reintroduce
//! silent mis-splitting.
//!
//! The two halves share the source behind a mutex, so a blocking read on one
//! side briefly serializes the other — an interleaved single stream cannot be
//! read in parallel the way two separate files can. Buffering of the not-yet-
//! consumed side is bounded by how far the two consumers drift apart, not by the
//! input size (see `MAX_PENDING_RECORDS`).

use std::collections::VecDeque;
use std::io::{self, BufRead, BufReader, Read};
use std::sync::{Arc, Mutex};

/// Lines per FASTQ record: `@header`, sequence, `+`, quality.
const FASTQ_LINES_PER_RECORD: usize = 4;

/// Safety valve on the per-side buffer. The two halves are meant to be read
/// roughly in lockstep (the single-threaded consumer alternates strictly; the
/// threaded pipeline keeps its streams within one batch via boundary sync), so
/// the buffer normally holds at most a handful of records. If one half is read
/// this far ahead of the other, something is wrong upstream and unbounded
/// buffering would exhaust memory — fail loudly instead.
const MAX_PENDING_RECORDS: usize = 1_000_000;

/// Shared source state behind the two reader halves.
struct DeinterleaverCore<R: BufRead> {
    /// The interleaved source (already decompressed).
    src: R,
    /// Records pulled from `src` that belong to a side which has not yet asked
    /// for them. `pending[0]` buffers R1 records, `pending[1]` buffers R2.
    pending: [VecDeque<Vec<u8>>; 2],
    /// Index of the next record to be pulled from `src`; even → R1, odd → R2.
    next_index: usize,
    /// Set once `src` reaches a clean EOF at a record boundary.
    eof: bool,
}

impl<R: BufRead> DeinterleaverCore<R> {
    /// Read the next whole 4-line record's raw bytes from `src`.
    ///
    /// Returns `Ok(None)` on a clean EOF at a record boundary, and an error if
    /// the stream ends partway through a record.
    fn read_one_record(&mut self) -> io::Result<Option<Vec<u8>>> {
        let mut record = Vec::new();
        for line_no in 0..FASTQ_LINES_PER_RECORD {
            let n = self.src.read_until(b'\n', &mut record)?;
            if n == 0 {
                if line_no == 0 {
                    return Ok(None); // clean EOF at a record boundary
                }
                return Err(io::Error::new(
                    io::ErrorKind::UnexpectedEof,
                    format!(
                        "interleaved FASTQ ended mid-record: read {line_no} of \
                         {FASTQ_LINES_PER_RECORD} lines"
                    ),
                ));
            }
        }
        Ok(Some(record))
    }

    /// Return the next record belonging to `side`, pulling (and buffering the
    /// other side's) records from `src` as needed. `Ok(None)` at EOF for `side`.
    fn next_for(&mut self, side: usize) -> io::Result<Option<Vec<u8>>> {
        if let Some(rec) = self.pending[side].pop_front() {
            return Ok(Some(rec));
        }
        loop {
            if self.eof {
                return Ok(None);
            }
            match self.read_one_record()? {
                None => {
                    self.eof = true;
                    return Ok(None);
                }
                Some(rec) => {
                    let rec_side = self.next_index % 2;
                    // Advance `next_index` only after the record is committed
                    // (returned or buffered), so an (allocator) panic in
                    // `push_back` cannot leave the index advanced past a dropped
                    // record and flip R1/R2 parity for everything after it.
                    if rec_side == side {
                        self.next_index += 1;
                        return Ok(Some(rec));
                    }
                    self.pending[rec_side].push_back(rec);
                    self.next_index += 1;
                    if self.pending[rec_side].len() > MAX_PENDING_RECORDS {
                        return Err(io::Error::other(format!(
                            "interleaved de-interleave buffer exceeded {MAX_PENDING_RECORDS} \
                             records: the two reads are being consumed too far out of sync"
                        )));
                    }
                }
            }
        }
    }
}

/// One half (R1 or R2) of a de-interleaved stream, yielding the raw bytes of its
/// records in order. Implements [`Read`]; the [`deinterleave`] constructor wraps
/// it in a [`BufReader`] for [`BufRead`].
struct DeinterleavedReader<R: BufRead> {
    core: Arc<Mutex<DeinterleaverCore<R>>>,
    side: usize,
    /// The current record's bytes, being drained into caller buffers.
    current: io::Cursor<Vec<u8>>,
}

impl<R: BufRead> Read for DeinterleavedReader<R> {
    fn read(&mut self, out: &mut [u8]) -> io::Result<usize> {
        loop {
            let n = self.current.read(out)?;
            if n > 0 {
                return Ok(n);
            }
            // Current record drained; fetch this side's next record under lock.
            // The read of `src` happens under this lock, so the two halves briefly
            // serialize here — an interleaved single stream cannot be read in
            // parallel the way two independent files can. Fail closed on a
            // poisoned core: a panic while the lock was held may have left the
            // framing state (`next_index`, `pending`) inconsistent, and resuming
            // could silently mis-pair every later record, violating the
            // byte-exact split contract.
            let next = {
                let mut core = self.core.lock().map_err(|_| {
                    io::Error::other(
                        "interleaved de-interleaver poisoned: a reader thread panicked, so \
                         record framing may be inconsistent",
                    )
                })?;
                core.next_for(self.side)?
            };
            match next {
                Some(rec) => self.current = io::Cursor::new(rec),
                None => return Ok(0), // EOF for this side
            }
        }
    }
}

/// Split one interleaved paired-end FASTQ stream into `(r1, r2)` readers.
///
/// Records alternate `R1, R2, R1, R2, …`; the first reader yields the R1 records
/// (indices 0, 2, 4, …) and the second the R2 records (1, 3, 5, …). Each record's
/// raw bytes are preserved exactly, so read-name comments (e.g. an inline UMI)
/// survive. An odd number of records leaves the final R1 without an R2 mate,
/// which downstream paired handling surfaces as an out-of-sync error.
///
/// The two readers share the source behind a mutex and may be read concurrently
/// from different threads.
pub fn deinterleave<R: BufRead + Send + 'static>(
    src: R,
) -> (Box<dyn BufRead + Send>, Box<dyn BufRead + Send>) {
    let core = Arc::new(Mutex::new(DeinterleaverCore {
        src,
        pending: [VecDeque::new(), VecDeque::new()],
        next_index: 0,
        eof: false,
    }));
    let r1 = DeinterleavedReader {
        core: Arc::clone(&core),
        side: 0,
        current: io::Cursor::new(Vec::new()),
    };
    let r2 = DeinterleavedReader { core, side: 1, current: io::Cursor::new(Vec::new()) };
    (Box::new(BufReader::new(r1)), Box::new(BufReader::new(r2)))
}

#[cfg(test)]
mod tests {
    use super::*;

    /// Build an interleaved FASTQ blob from `(name, seq, qual)` triples.
    fn fastq(records: &[(&str, &str, &str)]) -> Vec<u8> {
        let mut out = Vec::new();
        for (name, seq, qual) in records {
            out.extend_from_slice(b"@");
            out.extend_from_slice(name.as_bytes());
            out.extend_from_slice(b"\n");
            out.extend_from_slice(seq.as_bytes());
            out.extend_from_slice(b"\n+\n");
            out.extend_from_slice(qual.as_bytes());
            out.extend_from_slice(b"\n");
        }
        out
    }

    fn read_to_string(mut r: Box<dyn BufRead + Send>) -> String {
        let mut s = String::new();
        r.read_to_string(&mut s).expect("read");
        s
    }

    #[test]
    fn splits_alternating_records_into_two_streams() {
        let blob = fastq(&[
            ("r0/1", "ACGT", "IIII"),
            ("r0/2", "TGCA", "JJJJ"),
            ("r1/1", "AAAA", "IIII"),
            ("r1/2", "TTTT", "JJJJ"),
        ]);
        let (r1, r2) = deinterleave(io::Cursor::new(blob));
        assert_eq!(read_to_string(r1), "@r0/1\nACGT\n+\nIIII\n@r1/1\nAAAA\n+\nIIII\n");
        assert_eq!(read_to_string(r2), "@r0/2\nTGCA\n+\nJJJJ\n@r1/2\nTTTT\n+\nJJJJ\n");
    }

    #[test]
    fn preserves_read_name_comments_verbatim() {
        // The UMI-in-name use case: the comment after the space must survive so
        // `--extract-umis-from-read-names` still sees it.
        let blob = fastq(&[
            ("read1 1:N:0:ACGT+TGCA", "ACGT", "IIII"),
            ("read1 2:N:0:ACGT+TGCA", "TGCA", "JJJJ"),
        ]);
        let (r1, r2) = deinterleave(io::Cursor::new(blob));
        assert!(read_to_string(r1).starts_with("@read1 1:N:0:ACGT+TGCA\n"));
        let (r1b, r2b) = deinterleave(io::Cursor::new(fastq(&[
            ("read1 1:N:0:ACGT+TGCA", "ACGT", "IIII"),
            ("read1 2:N:0:ACGT+TGCA", "TGCA", "JJJJ"),
        ])));
        let _ = r1b;
        assert!(read_to_string(r2).starts_with("@read1 2:N:0:ACGT+TGCA\n"));
        let _ = r2b;
    }

    #[test]
    fn qual_line_starting_with_at_is_not_a_false_boundary() {
        // '@' is a valid Phred quality character; line-counting (not '@'-scanning)
        // must route by 4-line blocks regardless of content.
        let blob = fastq(&[("a/1", "ACGT", "@@@@"), ("a/2", "TGCA", "@@@@")]);
        let (r1, r2) = deinterleave(io::Cursor::new(blob));
        assert_eq!(read_to_string(r1), "@a/1\nACGT\n+\n@@@@\n");
        assert_eq!(read_to_string(r2), "@a/2\nTGCA\n+\n@@@@\n");
    }

    #[test]
    fn empty_input_yields_two_empty_streams() {
        let (r1, r2) = deinterleave(io::Cursor::new(Vec::new()));
        assert_eq!(read_to_string(r1), "");
        assert_eq!(read_to_string(r2), "");
    }

    #[test]
    fn odd_record_count_gives_r1_the_unpaired_tail() {
        let blob = fastq(&[("a/1", "AA", "II"), ("a/2", "TT", "JJ"), ("b/1", "CC", "II")]);
        let (r1, r2) = deinterleave(io::Cursor::new(blob));
        assert_eq!(read_to_string(r1), "@a/1\nAA\n+\nII\n@b/1\nCC\n+\nII\n");
        assert_eq!(read_to_string(r2), "@a/2\nTT\n+\nJJ\n");
    }

    #[test]
    fn missing_final_newline_still_routes_last_record() {
        let blob = b"@a/1\nAA\n+\nII\n@a/2\nTT\n+\nJJ".to_vec();
        let (r1, r2) = deinterleave(io::Cursor::new(blob));
        assert_eq!(read_to_string(r1), "@a/1\nAA\n+\nII\n");
        assert_eq!(read_to_string(r2), "@a/2\nTT\n+\nJJ");
    }

    #[test]
    fn truncated_record_is_an_error() {
        // A record with only 2 of 4 lines before EOF must not be silently split.
        // Whichever side reads far enough to reach the truncated record surfaces
        // the error (R1 gets there first here, reading ahead for its next record),
        // so assert order-independently that draining the stream errors.
        let blob = b"@a/1\nAA\n+\nII\n@a/2\nTT\n".to_vec();
        let (mut r1, mut r2) = deinterleave(io::Cursor::new(blob));
        let mut s1 = String::new();
        let mut s2 = String::new();
        let e1 = r1.read_to_string(&mut s1);
        let e2 = r2.read_to_string(&mut s2);
        let errs: Vec<_> = [e1, e2].into_iter().filter_map(Result::err).collect();
        assert!(!errs.is_empty(), "a truncated record must surface an error");
        assert!(
            errs.iter().all(|e| e.kind() == io::ErrorKind::UnexpectedEof),
            "truncation errors must be UnexpectedEof, got {errs:?}"
        );
    }

    #[test]
    fn concurrent_reads_from_two_threads_are_consistent() {
        let n = 2000;
        let mut recs: Vec<(String, String, String)> = Vec::new();
        for i in 0..n {
            recs.push((format!("r{i}/1"), "ACGTACGT".into(), "IIIIIIII".into()));
            recs.push((format!("r{i}/2"), "TGCATGCA".into(), "JJJJJJJJ".into()));
        }
        let record_refs: Vec<(&str, &str, &str)> =
            recs.iter().map(|(a, b, c)| (a.as_str(), b.as_str(), c.as_str())).collect();
        let (r1, r2) = deinterleave(io::Cursor::new(fastq(&record_refs)));
        let h1 = std::thread::spawn(move || read_to_string(r1));
        let h2 = std::thread::spawn(move || read_to_string(r2));
        let out1 = h1.join().unwrap();
        let out2 = h2.join().unwrap();
        assert_eq!(out1.matches("@r").count(), n, "R1 record count");
        assert_eq!(out2.matches("@r").count(), n, "R2 record count");
        assert!(out1.contains(&format!("@r{}/1\n", n - 1)));
        assert!(out2.contains(&format!("@r{}/2\n", n - 1)));
    }
}
