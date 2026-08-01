//! Comparison engines for `fgumi compare bams`.
//!
//! An "engine" pairs records from the two input streams (see
//! [`positional`](self::positional), added alongside the positional-alignment work) and
//! decides whether each pair is content-equal via a [`content::ContentPredicate`].
//! Splitting pairing (positional) from equality (content) keeps the pairing logic honest:
//! it can never quietly resync past a mismatch just because two *unrelated* records
//! happen to compare equal only because a lenient predicate ignores an actual divergence.
//!
//! [`header`] is the one engine-adjacent module that isn't about record pairing at all: it
//! compares the two inputs' `noodles::sam::Header` (`@HD`/`@SQ`/`@RG`, normalizing away
//! `@PG`/`@CO`), since a record-level match says nothing about whether the two files declare
//! the same reference dictionary, read groups, or sort order. Every engine below wires its
//! result in as an additional divergence source alongside its own pairing/content checks.

pub mod content;
pub mod header;
pub(crate) mod molecule_join;
pub mod positional;
pub mod sort_verify;

use std::path::{Path, PathBuf};

use fgumi_raw_bam::RawRecord;
use fgumi_sort::RawReadAheadReader;
use noodles::sam::Header;

/// A record stream that reports read failures as `Err`, not as end-of-stream.
///
/// [`RawReadAheadReader`] decodes on a background thread and therefore cannot
/// return an `io::Error` from `Iterator::next`: it ends iteration and parks the
/// error in a slot for the consumer to collect via
/// [`RawReadAheadReader::take_error`] *after* the loop. That convention is a
/// hazard in a comparison oracle. A truncated or CRC-corrupt BAM would read as a
/// **shorter** file rather than an error — reported as a record-count `DIFFER`,
/// or, when both inputs are equally corrupt, as a false `IDENTICAL`. Both are
/// silently wrong answers from the tool whose entire job is deciding whether two
/// files agree.
///
/// Rather than requiring every engine loop to remember a trailing `take_error()`
/// check, this adapter converts the error back into a final `Err` item at the
/// point the stream ends. Engines keep using ordinary `rec?` propagation and are
/// correct by construction; forgetting the check is not an available mistake.
///
/// A reader thread that dies mid-input parks an error of its own, so the same
/// `take_error()` call also catches a stream that stopped without ever reaching
/// its end — the failure `positional_compare` names "reader disconnected before
/// EOF". Every way this stream can end short of the input is therefore an `Err`
/// here, never a `None`.
pub(crate) struct CheckedRecords {
    inner: RawReadAheadReader,
    /// Which input this stream reads, so a failure names the offending file.
    ///
    /// The background reader reports a truncated stream as a bare
    /// `failed to fill whole buffer`, which in a two-input comparison does not say
    /// *which* input failed. Carrying the path lets the error identify it, matching
    /// the diagnostics the inline reader used to produce.
    label: String,
    /// Set once the inner iterator is exhausted, so a collected error is yielded
    /// exactly once and the stream stays fused afterwards.
    finished: bool,
}

impl CheckedRecords {
    fn new(inner: RawReadAheadReader, label: String) -> Self {
        Self { inner, label, finished: false }
    }
}

impl Iterator for CheckedRecords {
    type Item = std::io::Result<RawRecord>;

    fn next(&mut self) -> Option<Self::Item> {
        if self.finished {
            return None;
        }
        if let Some(record) = self.inner.next() {
            return Some(Ok(record));
        }
        self.finished = true;
        // End of stream: a parked error means the stream did not end cleanly, so
        // surface it instead of reporting EOF, naming the input it came from.
        self.inner.take_error().map(|e| {
            Err(std::io::Error::new(e.kind(), format!("reading records from {}: {e}", self.label)))
        })
    }
}

/// One input, opened exactly once: its header and its records from the same stream.
///
/// Passed to the single-pass engines instead of a path so the caller's open is the
/// only one. `CompareBams::execute` must read both headers before dispatching — an
/// incompatible pair has to hard-exit ahead of any engine — and if the engine then
/// opened the path again for its records, the command as a whole would require a
/// re-openable input even though each half of it reads the stream once. A FIFO, a
/// process substitution and stdin are all opened once and never again, so that
/// second open is the difference between accepting them and not.
///
/// `path` is carried only to name the input in diagnostics; nothing reopens it.
pub(crate) struct OpenedInput {
    /// Records, positioned just past the header.
    pub(crate) reader: CheckedRecords,
    /// The header parsed off the front of `reader`'s stream.
    pub(crate) header: Header,
    /// What to call this input in messages.
    pub(crate) path: PathBuf,
}

/// Split a `--threads` budget between the two inputs a comparison reads.
///
/// `--threads` is one budget for the command, so giving the full count to *each*
/// of two inputs would request twice the parallelism asked for: on a 12-core host
/// `--threads 8` would spawn 16 BGZF workers plus two read-ahead threads and the
/// main thread, paying context switches for no extra decode throughput.
///
/// `zipper` splits asymmetrically — 1 thread for its unmapped input, the full
/// count for its mapped input — because its two inputs differ in cost. A
/// comparison must fully decode both inputs and they are the same size, so an even
/// split is the right shape here.
///
/// The per-input read-ahead thread is deliberately *not* charged against this
/// budget: it overlaps I/O rather than decompressing, and `content` mode has
/// always spawned one reader thread per input regardless of `--threads` (see
/// `start_raw_batch_reader`).
fn decompression_threads_per_input(threads: usize) -> usize {
    (threads / 2).max(1)
}

impl OpenedInput {
    /// Open `path` once with single-threaded BGZF decode.
    ///
    /// For callers that have no thread budget to spend. Records still arrive via a
    /// read-ahead thread, so the decode overlaps the consumer's work.
    ///
    /// # Errors
    ///
    /// Returns an error if the input cannot be opened, is neither BAM nor SAM, or
    /// its header cannot be parsed.
    pub(crate) fn open(path: &Path) -> anyhow::Result<Self> {
        Self::open_one(path, 1)
    }

    /// Open both inputs of a comparison, splitting `threads` between them.
    ///
    /// The single entry point for the two-input engines, so the budget split lives
    /// in one place instead of being re-derived (or forgotten) per call site.
    ///
    /// # Errors
    ///
    /// Returns an error if either input cannot be opened, is neither BAM nor SAM,
    /// or its header cannot be parsed.
    pub(crate) fn open_pair(
        bam1: &Path,
        bam2: &Path,
        threads: usize,
    ) -> anyhow::Result<(Self, Self)> {
        let per_input = decompression_threads_per_input(threads);
        Ok((Self::open_one(bam1, per_input)?, Self::open_one(bam2, per_input)?))
    }

    /// Open `path` once with `decompression_threads` BGZF workers, reading ahead on
    /// a background thread.
    ///
    /// `decompression_threads > 1` gives the multithreaded BGZF decoder `content`
    /// mode has always used, and the read-ahead thread overlaps this input's decode
    /// with the other input's and with the comparison itself. Before this, the
    /// `sort` and `grouping` engines decoded both inputs inline on the single main
    /// thread — `--threads` created no threads at all (issue #686), leaving BGZF
    /// inflate plus its per-block CRC32 (~61% of the run) serialized on one core.
    ///
    /// Takes an already-divided per-input count, not the user's `--threads`; use
    /// [`Self::open_pair`] for the two-input case so the split happens once.
    ///
    /// Still exactly one open per input, so a FIFO, a process substitution, or
    /// stdin all still work: [`fgumi_bam_io::create_raw_bam_reader`] parses the
    /// header off the front of the same stream it then yields records from, which
    /// is the property this type exists to preserve (see the type docs).
    ///
    /// # Errors
    ///
    /// Returns an error if the input cannot be opened, is neither BAM nor SAM, or
    /// its header cannot be parsed.
    fn open_one(path: &Path, decompression_threads: usize) -> anyhow::Result<Self> {
        use anyhow::Context as _;

        let (reader, header) = fgumi_bam_io::create_raw_bam_reader(path, decompression_threads)
            .with_context(|| format!("opening raw BAM reader for {}", path.display()))?;
        Ok(Self {
            reader: CheckedRecords::new(
                RawReadAheadReader::new(reader),
                path.display().to_string(),
            ),
            header,
            path: path.to_path_buf(),
        })
    }
}

/// Append `msg()` to `details` unless it is already at the `max_diffs` cap — lazily, so
/// callers only pay for building the message string when it will actually be kept.
///
/// Shared by every engine's diff-collection loop (`molecule_join`, `metrics`, and, via
/// [`header::fold_header_diffs`], `positional`/`sort_verify`) — the capping/lazy-build
/// discipline is identical everywhere `diff_details` is populated.
pub(crate) fn push_diff(details: &mut Vec<String>, max_diffs: usize, msg: impl FnOnce() -> String) {
    if details.len() < max_diffs {
        details.push(msg());
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use rstest::rstest;

    /// Write a small valid BAM so a test can then damage its bytes.
    fn write_temp_bam(records: &[RawRecord]) -> tempfile::NamedTempFile {
        use noodles::sam::alignment::io::Write as _;
        let tmp = tempfile::NamedTempFile::new().expect("create temp BAM");
        let header = Header::default();
        let mut writer =
            noodles::bam::io::Writer::new(std::fs::File::create(tmp.path()).expect("create BAM"));
        writer.write_header(&header).expect("write header");
        for record in records {
            let buf = fgumi_raw_bam::raw_record_to_record_buf(record, &header)
                .expect("raw_record_to_record_buf");
            writer.write_alignment_record(&header, &buf).expect("write record");
        }
        writer.try_finish().expect("finish BAM");
        tmp
    }

    /// A damaged input must surface as an `Err` from the record stream, never as a
    /// clean end-of-stream.
    ///
    /// This is the invariant that makes the read-ahead reader safe to use in a
    /// comparison oracle. [`RawReadAheadReader`] decodes on a background thread and
    /// signals failure by ending iteration and parking the error for
    /// `take_error()`, so a consumer that only watches for `None` sees a damaged
    /// file as a *shorter* file. In `compare bams` that silently becomes a
    /// record-count `DIFFER` — or a false `IDENTICAL` when both inputs are damaged
    /// the same way — which is a wrong answer from the tool whose only job is
    /// deciding whether two files agree. [`CheckedRecords`] converts the parked
    /// error back into a final `Err`; this test is what keeps that true.
    ///
    /// Both damage modes are covered because they fail in different layers:
    /// truncation starves the block reader mid-record, while flipped payload bytes
    /// pass length checks and fail the per-block CRC32.
    #[rstest]
    #[case::truncated_mid_block(true)]
    #[case::corrupt_payload_crc(false)]
    fn damaged_input_yields_err_not_clean_eof(#[case] truncate: bool) {
        let records: Vec<RawRecord> = (0..2000)
            .map(|i| {
                fgumi_raw_bam::SamBuilder::new()
                    .read_name(format!("read{i:06}").as_bytes())
                    .flags(fgumi_raw_bam::flags::FIRST_SEGMENT)
                    .build()
            })
            .collect();
        let tmp = write_temp_bam(&records);
        let bytes = std::fs::read(tmp.path()).expect("read temp BAM");

        let damaged = if truncate {
            // Cut mid-BGZF-block: no EOF marker and an incomplete final block.
            bytes[..bytes.len() * 3 / 5].to_vec()
        } else {
            // Flip bytes deep inside a compressed block: still parses as a block,
            // cannot match its recorded CRC32.
            let mut v = bytes.clone();
            let off = v.len() / 2;
            for b in &mut v[off..off + 64] {
                *b ^= 0xFF;
            }
            v
        };
        let dmg = tempfile::NamedTempFile::new().expect("create damaged BAM");
        std::fs::write(dmg.path(), &damaged).expect("write damaged BAM");

        // Opening may already fail (the damage can land in the header region);
        // that is an acceptable hard failure. What must never happen is a
        // successful open whose stream then ends cleanly, hiding the damage.
        let Ok(opened) = OpenedInput::open(dmg.path()) else {
            return;
        };
        let mut saw_err = false;
        for item in opened.reader {
            if item.is_err() {
                saw_err = true;
                break;
            }
        }
        assert!(
            saw_err,
            "damaged input must report an error, not a clean EOF (truncate={truncate}): a \
             silently-short stream turns a corrupt file into a record-count DIFFER, or a false \
             IDENTICAL when both sides are damaged alike"
        );
    }

    /// `--threads` is one budget for the command, so two symmetric inputs split it
    /// rather than each taking the full count (which would double the requested
    /// parallelism). Always at least one worker per input.
    #[rstest]
    #[case::single(1, 1)]
    #[case::two(2, 1)]
    #[case::four(4, 2)]
    #[case::eight(8, 4)]
    #[case::odd_rounds_down(7, 3)]
    #[case::zero_still_yields_one(0, 1)]
    fn threads_are_split_between_the_two_inputs(#[case] threads: usize, #[case] expected: usize) {
        assert_eq!(decompression_threads_per_input(threads), expected);
    }
}
