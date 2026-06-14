//! `ReadSamChunks` source step + supporting helpers.
//!
//! 2-step parallel-parse SAM ingest. See
//! `docs/design/refactor/unified-chain-builder-architecture.md` (the
//! "Source resolver (delivered in Phase 1)" diagram) for how this step
//! fits into the per-command chain:
//!
//! ```text
//! ReadSamChunks (Serial + Reader)
//!     ↓ SamChunk { batch_serial, bytes, line_offsets }
//! ParseSamChunk (Parallel)
//!     ↓ DecodedRecordBatch
//! ```
//!
//! The Read step does NOT parse records — it reads bytes, finds line
//! boundaries with `memchr`, and emits a chunk plus its line-offset table.
//! Per-line parsing happens in parallel downstream.

use std::collections::VecDeque;
use std::io::{self, BufRead};
use std::sync::Arc;
use std::sync::atomic::{AtomicBool, Ordering};

use parking_lot::Mutex;

use crate::pipeline::core::Unpushed;
use crate::pipeline::core::held::HeldSlot;
use crate::pipeline::core::outputs::OrderedBytesSingle;
use crate::pipeline::core::queues::QueueSpec;
use crate::pipeline::core::reorder::BranchOrdering;
use crate::pipeline::core::step::{Affinity, Step, StepCtx, StepKind, StepOutcome, StepProfile};
use crate::pipeline::steps::types::SamChunk;

/// Target bytes per emitted [`SamChunk`]. Picked to match a typical L1 cache
/// working set (256 KB) — each parallel parser worker chews through one
/// chunk while the next one is in flight.
pub const DEFAULT_SAM_CHUNK_BYTES: usize = 256 * 1024;

/// Upper bound on the carryover buffer while still searching for the first
/// newline. A single SAM alignment record — even a long-read record with a
/// multi-megabase CIGAR/sequence — is realistically a few MB of text, so a
/// 1 GiB ceiling is enormously generous. Exceeding it without ever seeing a
/// `\n` means the input is not newline-delimited SAM (e.g. a truncated/binary
/// stream); we surface an error rather than letting `leftover` grow until the
/// process is OOM-killed.
const MAX_RECORD_BYTES: usize = 1024 * 1024 * 1024;

/// Mutable per-source state: the SAM byte reader, a small carryover buffer
/// for the trailing partial line, and an EOF flag.
///
/// Read at the byte level — header parsing is the caller's job (typically
/// via [`crate::pipeline::steps::source::InputSource::open`], which
/// drives `noodles::sam::io::Reader::read_header` then hands the post-header
/// reader's underlying buffered stream into [`ReadSamState::new`]).
pub struct ReadSamState {
    reader: Box<dyn BufRead + Send>,
    /// Trailing partial line from the previous `read_next_chunk` call. The
    /// next call prepends these bytes before scanning so records never get
    /// split across chunks. Private (encapsulation); use `leftover_len`
    /// and `leftover_bytes` for in-crate test assertions.
    leftover: Vec<u8>,
    eof: bool,
}

impl ReadSamState {
    #[must_use]
    pub fn new(reader: Box<dyn BufRead + Send>) -> Self {
        Self { reader, leftover: Vec::new(), eof: false }
    }

    /// Read the next chunk: append fresh bytes to `leftover`, split at the
    /// last newline, return the complete-records prefix plus its
    /// sentinel-form offset table. The trailing partial line stays in
    /// `self.leftover` for the next call.
    ///
    /// `target_bytes` is a **soft lower bound**: the loop reads at least
    /// that many bytes AND waits for at least one newline before
    /// splitting. This guarantees forward progress even when a single
    /// SAM line exceeds `target_bytes` (otherwise an upper-bound loop
    /// would stall — leftover would fill with bytes that contain no
    /// newline and the splitter would emit empty chunks forever).
    /// Practical implication: chunks may be larger than `target_bytes`
    /// if a record spans several reads, but never smaller than one
    /// complete record once data is available.
    ///
    /// Returns `Ok(None)` only when the reader has signalled EOF AND no
    /// bytes remain in `leftover` — i.e. the stream is fully consumed.
    /// If EOF is reached with leftover bytes present, a synthetic `\n`
    /// is appended so the final partial line becomes a complete record
    /// (SAM tools traditionally accept a missing trailing newline).
    ///
    /// # Panics (debug)
    ///
    /// Debug-asserts that `target_bytes <= u32::MAX`. Line offsets emitted
    /// by `split_complete_lines` are stored as `u32`; the production
    /// entry point [`ReadSamChunks::new`] also asserts this at
    /// construction. Direct callers (e.g. unit tests) should respect
    /// the same cap.
    ///
    /// # Errors
    ///
    /// Returns the underlying reader's I/O error.
    pub fn read_next_chunk(
        &mut self,
        target_bytes: usize,
    ) -> io::Result<Option<(Vec<u8>, Vec<u32>)>> {
        debug_assert!(
            u32::try_from(target_bytes).is_ok(),
            "target_bytes ({target_bytes}) exceeds u32::MAX; \
             line offsets in the emitted chunk are stored as u32"
        );
        // Soft-lower-bound read loop: keep pulling from the reader until
        // we have at least `target_bytes` AND at least one newline (so
        // the splitter can emit a complete record). `read_until(b'\n', ...)`
        // would force per-record I/O sync — we want larger reads to
        // amortize syscalls, so we pull raw bytes via `fill_buf`/`consume`.
        while !self.eof {
            // Stop when we have enough bytes for the target AND at least
            // one complete record sits in the buffer.
            if self.leftover.len() >= target_bytes
                && memchr::memchr(b'\n', &self.leftover).is_some()
            {
                break;
            }
            let chunk = self.reader.fill_buf()?;
            if chunk.is_empty() {
                self.eof = true;
                break;
            }
            let take = chunk.len();
            self.leftover.extend_from_slice(&chunk[..take]);
            self.reader.consume(take);

            // Guard against unbounded growth: if `leftover` exceeds the
            // per-record ceiling while still containing no newline, the input
            // is not newline-delimited SAM. Bail out instead of buffering the
            // whole stream into memory. (Once a newline is present the loop
            // breaks above, so this only fires on a genuinely newline-less span.)
            if self.leftover.len() > MAX_RECORD_BYTES
                && memchr::memchr(b'\n', &self.leftover).is_none()
            {
                return Err(io::Error::new(
                    io::ErrorKind::InvalidData,
                    format!(
                        "SAM record exceeds {MAX_RECORD_BYTES} bytes with no newline; \
                         input does not appear to be newline-delimited SAM text"
                    ),
                ));
            }
        }

        if self.leftover.is_empty() {
            // Reader is fully drained and nothing was held over.
            return Ok(None);
        }

        // At true EOF with leftover bytes: append a synthetic `\n` so the
        // trailing partial line becomes a complete record. SAM tools
        // traditionally accept files without a trailing newline.
        if self.eof && !self.leftover.ends_with(b"\n") {
            self.leftover.push(b'\n');
        }

        let (offsets, leftover_slice) = split_complete_lines(&self.leftover);
        let split_at = self.leftover.len() - leftover_slice.len();
        let new_leftover = self.leftover.split_off(split_at);
        let bytes = std::mem::replace(&mut self.leftover, new_leftover);
        Ok(Some((bytes, offsets)))
    }

    /// Test-only accessor: length of the carryover buffer holding the
    /// trailing partial line. Mostly useful for asserting that a
    /// mid-stream read held the right partial bytes.
    #[cfg(test)]
    fn leftover_len(&self) -> usize {
        self.leftover.len()
    }

    /// Test-only accessor: bytes of the carryover buffer.
    #[cfg(test)]
    fn leftover_bytes(&self) -> &[u8] {
        &self.leftover
    }
}

/// `Serial + sticky + Affinity::Reader` SAM-chunk source. Reads SAM text
/// bytes, splits at newline boundaries, and emits `SamChunk` items with
/// inline line-offset tables. Per-line parsing is delegated to the
/// downstream `Parallel` parser step.
pub struct ReadSamChunks {
    state: Arc<Mutex<Option<ReadSamState>>>,
    next_serial: u64,
    pending: VecDeque<SamChunk>,
    held: HeldSlot<Unpushed<SamChunk>>,
    target_chunk_bytes: usize,
    output_byte_limit: u64,
    finished: Arc<AtomicBool>,
}

impl ReadSamChunks {
    /// Build a SAM chunk source from a buffered reader positioned past the
    /// SAM header (typically obtained via `InputSource::open`).
    ///
    /// # Panics
    ///
    /// Panics if `target_chunk_bytes` exceeds `u32::MAX` — line offsets
    /// inside an emitted `SamChunk` are held as `u32`. In production all
    /// callsites pass `DEFAULT_SAM_CHUNK_BYTES` (256 KB) so this is
    /// unreachable; the assert exists to catch a configuration mistake
    /// at construction rather than panicking deep in the read loop.
    #[must_use]
    pub fn new(
        reader: Box<dyn BufRead + Send>,
        target_chunk_bytes: usize,
        output_byte_limit: u64,
    ) -> Self {
        assert!(
            u32::try_from(target_chunk_bytes).is_ok(),
            "target_chunk_bytes ({target_chunk_bytes}) exceeds u32::MAX; \
             line offsets in SamChunk are stored as u32"
        );
        Self {
            state: Arc::new(Mutex::new(Some(ReadSamState::new(reader)))),
            next_serial: 0,
            pending: VecDeque::new(),
            held: HeldSlot::new(),
            target_chunk_bytes: target_chunk_bytes.max(1),
            output_byte_limit,
            finished: Arc::new(AtomicBool::new(false)),
        }
    }
}

impl Step for ReadSamChunks {
    type Input = ();
    type Outputs = OrderedBytesSingle<SamChunk>;

    fn profile(&self) -> StepProfile {
        StepProfile {
            name: "ReadSamChunks",
            kind: StepKind::Serial,
            sticky: true,
            output_queues: vec![QueueSpec::ByteBounded { limit_bytes: self.output_byte_limit }],
            branch_ordering: vec![BranchOrdering::ByItemOrdinal],
        }
    }

    fn affinity(&self) -> Affinity {
        Affinity::Reader
    }

    fn try_run(&mut self, ctx: &mut StepCtx<'_, Self>) -> io::Result<StepOutcome> {
        // 1. Drain the held slot first.
        if let Some(unpushed) = self.held.take() {
            match ctx.outputs.retry(unpushed) {
                Ok(()) => {}
                Err(again) => {
                    self.held.put(again);
                    return Ok(StepOutcome::Contention);
                }
            }
        }

        // 2. Drain pending chunks.
        if let Some(chunk) = self.pending.pop_front() {
            match ctx.outputs.push(chunk) {
                Ok(()) => return Ok(StepOutcome::Progress),
                Err(unpushed) => {
                    self.held.put(unpushed);
                    return Ok(StepOutcome::Progress);
                }
            }
        }

        if self.finished.load(Ordering::Acquire) {
            return Ok(StepOutcome::Finished);
        }

        // 3. Pull the next chunk from the reader.
        let next = {
            let mut guard = self.state.lock();
            let state = guard.as_mut().expect("ReadSamChunks: state missing — was clone() called?");
            state.read_next_chunk(self.target_chunk_bytes)?
        };

        let Some((bytes, line_offsets)) = next else {
            self.finished.store(true, Ordering::Release);
            return Ok(StepOutcome::Finished);
        };

        let serial = self.next_serial;
        self.next_serial += 1;
        let chunk = SamChunk { batch_serial: serial, bytes, line_offsets };
        match ctx.outputs.push(chunk) {
            Ok(()) => Ok(StepOutcome::Progress),
            Err(unpushed) => {
                self.held.put(unpushed);
                Ok(StepOutcome::Progress)
            }
        }
    }
}

/// Split `data` into complete `\n`-terminated records.
///
/// Returns `(line_offsets, leftover)` where:
///
/// - `line_offsets` is the sentinel-form table — `N+1` entries describing
///   `N` complete lines. Line `i` is `data[line_offsets[i]..line_offsets[i+1]]`.
///   Each line's bytes include its trailing `\n`. Empty if `data` contains
///   no complete line.
/// - `leftover` is the trailing partial line (bytes after the last `\n`).
///   Empty if `data` ends on a `\n` boundary. The caller carries this
///   forward into the next read so records aren't split mid-line.
///
/// # Panics
///
/// Panics if `data.len()` exceeds `u32::MAX` (we hold offsets as `u32`
/// because individual chunks are sized in the hundreds-of-KB range).
#[must_use]
pub fn split_complete_lines(data: &[u8]) -> (Vec<u32>, &[u8]) {
    // Find the position just past the last newline. Everything up to that
    // position is "complete lines"; the remainder is the partial-line
    // leftover for the next read.
    let Some(last_nl) = memchr::memrchr(b'\n', data) else {
        return (Vec::new(), data);
    };
    let split = last_nl + 1; // include the trailing \n in the complete region
    let (complete, leftover) = data.split_at(split);

    // Walk newline positions to build the sentinel-form offset table.
    let mut offsets: Vec<u32> = Vec::new();
    offsets.push(0);
    for nl in memchr::memchr_iter(b'\n', complete) {
        offsets.push(u32::try_from(nl + 1).expect("chunk size fits in u32"));
    }
    (offsets, leftover)
}

#[cfg(test)]
mod tests {
    use super::*;
    use std::io::{BufReader, Cursor};

    fn reader_from(bytes: &[u8]) -> Box<dyn std::io::BufRead + Send> {
        Box::new(BufReader::new(Cursor::new(bytes.to_vec())))
    }

    #[test]
    fn split_complete_lines_empty_input_returns_empty_offsets_and_empty_leftover() {
        let (offsets, leftover) = split_complete_lines(b"");
        assert!(offsets.is_empty(), "expected no line offsets for empty input");
        assert!(leftover.is_empty(), "expected no leftover for empty input");
    }

    #[test]
    fn split_complete_lines_single_line_terminated_by_newline_no_leftover() {
        // "abc\n" — one complete record, ends at newline. Sentinel-form
        // offsets are [0, 4]: line 0 is bytes[0..4] = "abc\n".
        let (offsets, leftover) = split_complete_lines(b"abc\n");
        assert_eq!(offsets, vec![0, 4]);
        assert!(leftover.is_empty(), "expected no leftover when input ends on \\n");
    }

    #[test]
    fn split_complete_lines_multiple_lines_no_trailing_newline_leaves_partial_leftover() {
        // "abc\nde\nf" — two complete lines ("abc\n", "de\n") and a partial
        // line "f" carried over. Offsets [0, 4, 7] index into bytes[..7].
        let input = b"abc\nde\nf";
        let (offsets, leftover) = split_complete_lines(input);
        assert_eq!(offsets, vec![0, 4, 7]);
        assert_eq!(leftover, b"f");
        // Verify the offsets actually slice out the lines.
        assert_eq!(&input[offsets[0] as usize..offsets[1] as usize], b"abc\n");
        assert_eq!(&input[offsets[1] as usize..offsets[2] as usize], b"de\n");
    }

    #[test]
    fn split_complete_lines_only_partial_line_returns_empty_offsets() {
        // "abc" — no newline at all. Nothing complete; everything is leftover.
        let (offsets, leftover) = split_complete_lines(b"abc");
        assert!(offsets.is_empty());
        assert_eq!(leftover, b"abc");
    }

    #[test]
    fn read_next_chunk_emits_complete_lines_and_holds_partial_as_leftover() {
        // Soft-lower-bound semantics: with target=10 against a 16-byte
        // cursor that delivers everything in one `fill_buf`, the loop
        // reads it all (target satisfied + at least one `\n` present),
        // splits at the last newline, and holds the trailing partial
        // line "linX" as the carryover for the next call.
        let mut state = ReadSamState::new(reader_from(b"line1\nline2\nlinX"));
        let (bytes, offsets) =
            state.read_next_chunk(10).expect("read ok").expect("at least one chunk");
        assert_eq!(offsets, vec![0, 6, 12]);
        assert_eq!(&bytes[0..6], b"line1\n");
        assert_eq!(&bytes[6..12], b"line2\n");
        assert_eq!(state.leftover_bytes(), b"linX");
    }

    #[test]
    fn read_next_chunk_flushes_partial_line_at_eof_with_synthetic_newline() {
        // Input without a trailing newline. A large target reads everything
        // in one go, then EOF triggers the synthetic-`\n` flush so the
        // final partial line becomes a complete record.
        let mut state = ReadSamState::new(reader_from(b"only_one"));
        let (bytes, offsets) =
            state.read_next_chunk(1024).expect("read ok").expect("at least one chunk");
        assert_eq!(offsets, vec![0, 9]);
        assert_eq!(bytes, b"only_one\n");
        assert_eq!(state.leftover_len(), 0, "EOF flush should drain leftover");
    }

    #[test]
    fn read_next_chunk_returns_none_after_full_drain() {
        // After EOF + flushed final record, the next call returns None.
        let mut state = ReadSamState::new(reader_from(b"a\nb\n"));
        let _first = state.read_next_chunk(1024).expect("read ok").expect("first chunk");
        let second = state.read_next_chunk(1024).expect("read ok");
        assert!(second.is_none(), "expected None on second call after drain");
    }

    #[test]
    fn read_next_chunk_target_smaller_than_line_does_not_stall() {
        // Regression test for the target-bytes-as-upper-bound stall: if
        // `target_bytes` was smaller than the smallest line, the old
        // upper-bound loop terminated with leftover holding partial
        // bytes and no newline. Subsequent calls returned empty chunks
        // forever (the exit condition was already satisfied).
        //
        // Fix treats `target_bytes` as a soft LOWER bound: read at least
        // that many bytes AND wait for a newline before splitting.
        //
        // Input: two 5-byte records against a 4-byte target. First call
        // must emit BOTH records (the Cursor delivers everything in one
        // `fill_buf`, the second `\n` satisfies the lower bound). A
        // follow-up call must return None — drained.
        let mut state = ReadSamState::new(reader_from(b"hello\nworld\n"));
        let (bytes, offsets) =
            state.read_next_chunk(4).expect("read ok").expect("at least one chunk");
        assert_eq!(bytes, b"hello\nworld\n", "must drain both complete records, not stall");
        assert_eq!(offsets, vec![0, 6, 12], "two records: hello\\n at [0..6], world\\n at [6..12]");
        assert_eq!(state.leftover_len(), 0, "leftover should be empty after draining");

        // Second call: reader is drained, leftover empty → returns None.
        let second = state.read_next_chunk(4).expect("read ok");
        assert!(second.is_none(), "second call must return None — stream fully consumed");
    }

    #[test]
    fn read_next_chunk_returns_none_for_empty_input() {
        let mut state = ReadSamState::new(reader_from(b""));
        let result = state.read_next_chunk(1024).expect("read ok");
        assert!(result.is_none());
    }

    #[test]
    fn split_complete_lines_consecutive_newlines_yield_empty_lines() {
        // "\n\n" — two empty records, both terminated. The parser should
        // accept zero-length lines; SAM doesn't naturally produce them but
        // the splitter shouldn't drop them either (drives the invariant
        // that offsets.windows(2) cover every byte of `complete`).
        let (offsets, leftover) = split_complete_lines(b"\n\n");
        assert_eq!(offsets, vec![0, 1, 2]);
        assert!(leftover.is_empty());
    }

    #[test]
    fn read_sam_chunks_profile_advertises_serial_reader_byordinal() {
        let step = ReadSamChunks::new(reader_from(b""), 64 * 1024, 1024 * 1024);
        let profile = step.profile();
        assert_eq!(profile.name, "ReadSamChunks");
        assert_eq!(profile.kind, crate::pipeline::core::step::StepKind::Serial);
        assert!(profile.sticky);
        assert_eq!(step.affinity(), crate::pipeline::core::step::Affinity::Reader);
        assert_eq!(
            profile.branch_ordering,
            vec![crate::pipeline::core::reorder::BranchOrdering::ByItemOrdinal]
        );
        assert!(matches!(
            profile.output_queues[0],
            crate::pipeline::core::queues::QueueSpec::ByteBounded { .. }
        ));
    }
}
