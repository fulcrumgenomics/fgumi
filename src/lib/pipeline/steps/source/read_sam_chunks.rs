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
    /// Ceiling on a single newline-free span before the input is rejected as
    /// not newline-delimited. A field rather than a bare constant so tests can
    /// drive the rejection path with a small cap instead of buffering a
    /// gigabyte; production callers get [`MAX_RECORD_BYTES`] via `new`.
    max_record_bytes: usize,
    /// Batch cut mode. [`BatchCut::Record`] (default) splits at any line
    /// boundary; [`BatchCut::Queryname`] carries the trailing partial queryname
    /// run so every emitted chunk is closed under queryname (mirrors the BAM
    /// `FindBamBoundaries` cut).
    cut: crate::pipeline::steps::boundaries::state::BatchCut,
    /// Coalescing target for `BatchCut::Queryname`: keep the closed prefix in
    /// `leftover` until it reaches this many bytes before emitting. `0`
    /// disables coalescing. Ignored for `Record`.
    min_emit_bytes: usize,
    /// QNAME of the last line in the most recently emitted chunk, under
    /// `BatchCut::Queryname`. The next emitted chunk's first line must have a
    /// different QNAME or the closed-under-queryname invariant is broken; the
    /// self-check errors if not. Empty when nothing has been emitted yet.
    last_emitted_qname: Vec<u8>,
}

impl ReadSamState {
    #[must_use]
    pub fn new(reader: Box<dyn BufRead + Send>) -> Self {
        Self {
            reader,
            leftover: Vec::new(),
            eof: false,
            max_record_bytes: MAX_RECORD_BYTES,
            cut: crate::pipeline::steps::boundaries::state::BatchCut::Record,
            min_emit_bytes: 0,
            last_emitted_qname: Vec::new(),
        }
    }

    /// Select the batch cut mode + coalescing target (see [`BatchCut`]).
    #[must_use]
    fn with_cut(
        mut self,
        cut: crate::pipeline::steps::boundaries::state::BatchCut,
        min_emit_bytes: usize,
    ) -> Self {
        self.cut = cut;
        self.min_emit_bytes = min_emit_bytes;
        self
    }

    /// Construct with a custom newline-free span ceiling. Test-facing: it makes
    /// the `InvalidData` rejection reachable without a gigabyte of input.
    #[cfg(test)]
    fn with_max_record_bytes(reader: Box<dyn BufRead + Send>, max_record_bytes: usize) -> Self {
        Self {
            reader,
            leftover: Vec::new(),
            eof: false,
            max_record_bytes,
            cut: crate::pipeline::steps::boundaries::state::BatchCut::Record,
            min_emit_bytes: 0,
            last_emitted_qname: Vec::new(),
        }
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
        // Scan for the first newline incrementally. Re-running `memchr` over the
        // whole of `leftover` on every iteration re-reads bytes already known to
        // be newline-free, which is quadratic in the length of a newline-less
        // span — exactly the input the ceiling below exists to reject, so the
        // rejection path was the slowest one. `scanned` marks how far the search
        // has already reached; everything before it is known newline-free.
        // In queryname mode, a chunk whose complete lines are all one open run
        // closes nothing, so the cut carries everything and we must read MORE
        // before we can emit. The read-loop break gate is `leftover.len() >=
        // effective_target`, so once the cut carries everything we raise the
        // target by `target_bytes` and loop, forcing another read (otherwise the
        // gate stays satisfied by the already-buffered leftover and the reader
        // never advances to EOF — an infinite loop). `Record` mode never carries
        // a whole chunk, so `effective_target` stays at `target_bytes` there.
        let mut effective_target = target_bytes;
        loop {
            let mut newline_at = memchr::memchr(b'\n', &self.leftover);
            let mut scanned = if newline_at.is_some() { 0 } else { self.leftover.len() };

            while !self.eof {
                // Stop when we have enough bytes for the target AND at least
                // one complete record sits in the buffer.
                if self.leftover.len() >= effective_target && newline_at.is_some() {
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

                if newline_at.is_none() {
                    newline_at =
                        memchr::memchr(b'\n', &self.leftover[scanned..]).map(|off| scanned + off);
                    if newline_at.is_none() {
                        scanned = self.leftover.len();
                    }
                }

                // Guard against unbounded growth: if `leftover` exceeds the
                // per-record ceiling while still containing no newline, the input
                // is not newline-delimited SAM. Bail out instead of buffering the
                // whole stream into memory. (Once a newline is present the loop
                // breaks above, so this only fires on a genuinely newline-less span.)
                if newline_at.is_none() && self.leftover.len() > self.max_record_bytes {
                    let max = self.max_record_bytes;
                    return Err(io::Error::new(
                        io::ErrorKind::InvalidData,
                        format!(
                            "SAM record exceeds {max} bytes with no newline; \
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

            let (mut offsets, leftover_slice) = split_complete_lines(&self.leftover);
            let mut split_at = self.leftover.len() - leftover_slice.len();

            // Queryname cut: pull `split_at`/`offsets` back so the emitted chunk is
            // closed under queryname — the trailing partial queryname run stays in
            // `leftover` with the partial line, to be emitted once a differing
            // QNAME proves the run ended (or at EOF). Mirrors the BAM cut.
            if self.cut == crate::pipeline::steps::boundaries::state::BatchCut::Queryname {
                self.apply_queryname_cut(&mut offsets, &mut split_at)?;
                // Cut carried everything (nothing closed) and we are not at EOF:
                // raise the target and read more so the reader advances toward a
                // queryname boundary or EOF. Guaranteed to terminate — each pass
                // either reads (draining toward EOF) or the growth reveals a
                // boundary.
                if split_at == 0 && !self.eof {
                    effective_target = effective_target.saturating_add(target_bytes);
                    continue;
                }
            }

            let new_leftover = self.leftover.split_off(split_at);
            let bytes = std::mem::replace(&mut self.leftover, new_leftover);
            return Ok(Some((bytes, offsets)));
        } // end outer loop
    }

    /// Pull `offsets`/`split_at` back so the emitted lines are closed under
    /// queryname. `offsets` is the sentinel-form table for the complete lines
    /// in `self.leftover[..split_at]`; on return it describes only the emitted
    /// prefix, and `split_at` marks where `leftover` is cut (carry = the rest).
    ///
    /// Walks back from the last complete line while adjacent QNAMEs match, then
    /// (if a coalescing target is set and the closed prefix is short, and we are
    /// not at EOF) carries everything to accumulate more. Emits nothing (offsets
    /// truncated to the leading sentinel, `split_at = 0`) when the whole chunk
    /// is one open run — `read_next_chunk`'s caller loop keeps reading.
    ///
    /// # Errors
    ///
    /// `InvalidData` if the closed-under-queryname self-check trips.
    fn apply_queryname_cut(
        &mut self,
        offsets: &mut Vec<u32>,
        split_at: &mut usize,
    ) -> io::Result<()> {
        let n_lines = offsets.len().saturating_sub(1);
        if n_lines == 0 {
            return Ok(());
        }
        // QNAME of line `i` (bytes before the first tab; the line includes its
        // trailing `\n`). A tab-less line is treated as its own whole name.
        let qname = |i: usize| -> &[u8] {
            let mut line = &self.leftover[offsets[i] as usize..offsets[i + 1] as usize];
            if line.last() == Some(&b'\n') {
                line = &line[..line.len() - 1];
            }
            let end = memchr::memchr(b'\t', line).unwrap_or(line.len());
            &line[..end]
        };

        // Choose how many leading lines to emit (`emit`): the closed prefix.
        let emit = if self.eof {
            // At EOF the trailing run is complete (nothing can extend it), so
            // emit every complete line. This also prevents an infinite carry
            // when the whole final stream is one QNAME (the backward walk would
            // otherwise pick 0 and the caller would re-read the same leftover
            // forever).
            n_lines
        } else {
            // `j` = first line of the trailing (open) queryname run to carry.
            let mut j = n_lines - 1;
            while j > 0 && qname(j - 1) == qname(j) {
                j -= 1;
            }
            // Coalescing: carry everything until the closed prefix is big enough.
            if self.min_emit_bytes > 0 && (offsets[j] as usize) < self.min_emit_bytes {
                0
            } else {
                j
            }
        };

        if emit == 0 {
            // Nothing closed to emit: carry all complete lines + the partial
            // tail. `split_at = 0` leaves everything in `leftover`.
            offsets.clear();
            *split_at = 0;
            return Ok(());
        }

        // Self-check: the first emitted line's QNAME must differ from the
        // previous emitted chunk's last QNAME (else a run straddled a chunk).
        let first_q = qname(0).to_vec();
        if !self.last_emitted_qname.is_empty() && self.last_emitted_qname == first_q {
            return Err(io::Error::new(
                io::ErrorKind::InvalidData,
                format!(
                    "ReadSamChunks: queryname cut invariant violated \
                     (a run of QNAME {:?} straddled a chunk boundary)",
                    String::from_utf8_lossy(&first_q),
                ),
            ));
        }
        self.last_emitted_qname = qname(emit - 1).to_vec();

        // Emit lines 0..emit; carry the rest. When `emit == n_lines` (EOF or a
        // run boundary exactly at the last complete line) this keeps everything
        // — `split_at`/`offsets` are already correct, and the truncate is a
        // no-op.
        *split_at = offsets[emit] as usize;
        offsets.truncate(emit + 1);
        Ok(())
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
            held: HeldSlot::new(),
            target_chunk_bytes: target_chunk_bytes.max(1),
            output_byte_limit,
            finished: Arc::new(AtomicBool::new(false)),
        }
    }

    /// Like [`Self::new`] but emits chunks closed under queryname (every
    /// queryname run lies entirely within one chunk), coalescing the emitted
    /// prefix to ~`min_emit_bytes`. Used when the first stage groups by
    /// queryname so the downstream grouper can run as a parallel map. `Record`
    /// callers use [`Self::new`] unchanged.
    ///
    /// # Panics
    ///
    /// Panics if `target_chunk_bytes` exceeds `u32::MAX` (see [`Self::new`]).
    #[must_use]
    pub fn new_queryname_cut(
        reader: Box<dyn BufRead + Send>,
        target_chunk_bytes: usize,
        output_byte_limit: u64,
        min_emit_bytes: usize,
    ) -> Self {
        assert!(
            u32::try_from(target_chunk_bytes).is_ok(),
            "target_chunk_bytes ({target_chunk_bytes}) exceeds u32::MAX"
        );
        let state = ReadSamState::new(reader).with_cut(
            crate::pipeline::steps::boundaries::state::BatchCut::Queryname,
            min_emit_bytes,
        );
        Self {
            state: Arc::new(Mutex::new(Some(state))),
            next_serial: 0,
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

        // Under the queryname cut, a chunk whose whole content was one open run
        // is carried entirely to `leftover` and comes back with no complete
        // lines. Don't mint a serial for it (that would leave a gap in the dense
        // ordinal the downstream reorder needs) — report Progress; the reader
        // advanced, so the next dispatch continues from the carried bytes.
        if line_offsets.len() <= 1 {
            return Ok(StepOutcome::Progress);
        }

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

    /// A `BufRead` that releases at most `piece` bytes per `fill_buf`.
    ///
    /// `BufReader<Cursor<_>>` hands the whole input over on the first
    /// `fill_buf`, so every test built on `reader_from` leaves the read loop
    /// after one iteration — with `scanned` still `0`. That makes the
    /// incremental-scan arithmetic in `read_next_chunk` (which converts a
    /// slice-relative newline index back to an absolute one) unreachable. This
    /// reader forces the loop to append repeatedly and resume the search from
    /// `scanned`.
    struct Dribble {
        data: Vec<u8>,
        pos: usize,
        piece: usize,
    }

    impl Dribble {
        fn boxed(data: Vec<u8>, piece: usize) -> Box<dyn std::io::BufRead + Send> {
            Box::new(Self { data, pos: 0, piece })
        }
    }

    impl io::Read for Dribble {
        fn read(&mut self, out: &mut [u8]) -> io::Result<usize> {
            let n = self.piece.min(out.len()).min(self.data.len() - self.pos);
            out[..n].copy_from_slice(&self.data[self.pos..self.pos + n]);
            self.pos += n;
            Ok(n)
        }
    }

    impl io::BufRead for Dribble {
        fn fill_buf(&mut self) -> io::Result<&[u8]> {
            let end = (self.pos + self.piece).min(self.data.len());
            Ok(&self.data[self.pos..end])
        }

        fn consume(&mut self, amt: usize) {
            self.pos = (self.pos + amt).min(self.data.len());
        }
    }

    /// The newline search must resume at `scanned` and still report an
    /// **absolute** offset. The first `\n` sits at byte 20, well past the 8-byte
    /// piece size, so it is found on the third `fill_buf` in a slice that starts
    /// at 16 — a sign error or a missing `scanned +` would place it at 4.
    #[test]
    fn read_next_chunk_resumes_the_newline_search_across_fill_buf_calls() {
        let mut data = vec![b'A'; 20];
        data.push(b'\n');
        data.extend_from_slice(b"sec");

        let mut state = ReadSamState::new(Dribble::boxed(data, 8));

        let (bytes, offsets) =
            state.read_next_chunk(4).expect("read succeeds").expect("a chunk is emitted");

        // Three 8-byte pieces are pulled before the newline is seen, so the
        // buffer holds 24 bytes; only the 21 through the newline are complete.
        let mut expected = vec![b'A'; 20];
        expected.push(b'\n');
        assert_eq!(bytes, expected);
        assert_eq!(offsets, vec![0, 21]);
        assert_eq!(state.leftover_bytes(), b"sec");
    }

    /// The newline-free ceiling must also fire when bytes trickle in, since the
    /// check runs inside the same loop the incremental scan advances.
    #[test]
    fn the_newline_free_cap_still_fires_on_small_pieces() {
        const CAP: usize = 64;

        let mut state =
            ReadSamState::with_max_record_bytes(Dribble::boxed(vec![b'A'; CAP * 4], 7), CAP);

        let err = state.read_next_chunk(16).expect_err("newline-free input is rejected");

        assert_eq!(err.kind(), io::ErrorKind::InvalidData);
        assert!(err.to_string().contains("no newline"), "got: {err}");
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

    /// A newline-free span past the ceiling must be rejected as `InvalidData`
    /// rather than buffered. This is the path the cap exists for — a
    /// non-newline-delimited input (a binary file handed to the SAM reader)
    /// would otherwise pull the whole stream into `leftover`.
    ///
    /// Driven through a small cap so the rejection is reachable without
    /// allocating the production-sized ceiling.
    #[test]
    fn a_newline_free_span_past_the_cap_is_rejected() {
        const CAP: usize = 4 * 1024;
        let input = vec![b'A'; CAP * 4]; // no newline anywhere
        let mut state = ReadSamState::with_max_record_bytes(
            Box::new(std::io::BufReader::new(std::io::Cursor::new(input))),
            CAP,
        );

        let err = state
            .read_next_chunk(1024)
            .expect_err("a newline-free input must be rejected, not buffered");
        assert_eq!(err.kind(), io::ErrorKind::InvalidData);
        assert!(
            err.to_string().contains("newline-delimited"),
            "error should name the actual problem, got: {err}"
        );
    }

    /// The cap must not fire on legitimate input whose records simply exceed the
    /// soft target: a single long line under the ceiling still parses.
    #[test]
    fn a_long_but_terminated_line_under_the_cap_is_accepted() {
        const CAP: usize = 4 * 1024;
        let mut input = vec![b'A'; CAP / 2];
        input.push(b'\n');
        let expected = input.clone();
        let mut state = ReadSamState::with_max_record_bytes(
            Box::new(std::io::BufReader::new(std::io::Cursor::new(input))),
            CAP,
        );

        let (bytes, offsets) = state
            .read_next_chunk(64)
            .expect("a terminated line under the cap must parse")
            .expect("one complete record");
        assert_eq!(bytes, expected);
        assert_eq!(offsets, vec![0, u32::try_from(expected.len()).unwrap()]);
    }

    // -- Queryname cut (BatchCut::Queryname) --

    /// One SAM line: `qname\t<rest>\n`. The rest is a fixed opaque field; only
    /// the QNAME (bytes before the first tab) matters to the cut.
    fn sam_line(qname: &str) -> Vec<u8> {
        format!("{qname}\t0\t*\t0\t0\t*\t*\t0\t0\t*\t*\n").into_bytes()
    }

    /// Drive a queryname-cut `ReadSamState` over `input` (via a `Dribble` so the
    /// read loop iterates) with the given chunk target, returning the QNAMEs of
    /// every emitted line, grouped per emitted chunk.
    fn qn_sam_chunks(input: Vec<u8>, target: usize, min_emit: usize) -> Vec<Vec<String>> {
        let mut state = ReadSamState::new(Dribble::boxed(input, 17))
            .with_cut(crate::pipeline::steps::boundaries::state::BatchCut::Queryname, min_emit);
        let mut out: Vec<Vec<String>> = Vec::new();
        loop {
            match state.read_next_chunk(target).expect("read ok") {
                None => break,
                Some((bytes, offsets)) => {
                    if offsets.len() <= 1 {
                        continue; // absorbed (whole chunk one open run)
                    }
                    let mut names = Vec::new();
                    for w in offsets.windows(2) {
                        let mut line = &bytes[w[0] as usize..w[1] as usize];
                        if line.last() == Some(&b'\n') {
                            line = &line[..line.len() - 1];
                        }
                        let end = line.iter().position(|&b| b == b'\t').unwrap_or(line.len());
                        names.push(String::from_utf8_lossy(&line[..end]).into_owned());
                    }
                    out.push(names);
                }
            }
        }
        out
    }

    fn qn_sam_flat(chunks: &[Vec<String>]) -> Vec<String> {
        chunks.iter().flatten().cloned().collect()
    }

    fn qn_sam_no_straddle(chunks: &[Vec<String>]) {
        for w in chunks.windows(2) {
            assert_ne!(w[0].last(), w[1].first(), "a QNAME run straddled a chunk boundary");
        }
    }

    #[test]
    fn sam_queryname_cut_carries_a_run_across_chunks() {
        let mut input = Vec::new();
        for q in ["a", "a", "b", "b", "b", "c"] {
            input.extend(sam_line(q));
        }
        // Small target so the reader emits several chunks mid-run.
        let chunks = qn_sam_chunks(input, 40, 0);
        qn_sam_no_straddle(&chunks);
        assert_eq!(qn_sam_flat(&chunks), vec!["a", "a", "b", "b", "b", "c"],);
    }

    #[test]
    fn sam_queryname_cut_single_run_whole_stream() {
        let mut input = Vec::new();
        for _ in 0..6 {
            input.extend(sam_line("solo"));
        }
        // Every line is QNAME "solo": nothing closes until EOF, then all emit.
        let chunks = qn_sam_chunks(input, 20, 0);
        assert_eq!(qn_sam_flat(&chunks), vec!["solo"; 6]);
    }

    #[test]
    fn sam_queryname_cut_coalesces_to_min_emit() {
        let mut input = Vec::new();
        for i in 0..10 {
            input.extend(sam_line(&format!("q{i}")));
        }
        // Large min_emit holds everything until EOF: one emitted chunk.
        let chunks = qn_sam_chunks(input, 20, 1 << 20);
        assert_eq!(chunks.len(), 1);
        assert_eq!(chunks[0].len(), 10);
    }

    #[test]
    fn sam_queryname_cut_missing_trailing_newline_flushes_at_eof() {
        let mut input = Vec::new();
        input.extend(sam_line("a"));
        input.extend(sam_line("a"));
        input.extend(b"b\t0\t*\t0\t0\t*\t*\t0\t0\t*\t*"); // no trailing \n
        let chunks = qn_sam_chunks(input, 30, 0);
        qn_sam_no_straddle(&chunks);
        assert_eq!(qn_sam_flat(&chunks), vec!["a", "a", "b"]);
    }

    #[test]
    fn sam_queryname_cut_tabless_line_is_its_own_name() {
        // A line with no tab is treated as a whole-line QNAME; adjacent
        // identical tab-less lines group, distinct ones split.
        let input = b"xx\nxx\nyy\n".to_vec();
        let chunks = qn_sam_chunks(input, 4, 0);
        qn_sam_no_straddle(&chunks);
        assert_eq!(qn_sam_flat(&chunks), vec!["xx", "xx", "yy"]);
    }
}
