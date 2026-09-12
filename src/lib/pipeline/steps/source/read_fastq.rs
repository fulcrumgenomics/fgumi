//! `ReadFastqInputs` source step: reads (and gzip-decompresses) raw FASTQ
//! bytes from one or more `BufRead` streams in round-robin order and emits
//! `FastqRawChunk` items. Parsing into `FastqRecord`s is deferred to the
//! downstream `Parallel` `ParseAndZipFastqN` step so the decode/parse work
//! can fan out across worker threads rather than bottlenecking on the
//! single reader thread.
//!
//! ## Per-stream vs. all-streams instantiation
//!
//! For the common paired-end case (R1 + R2) the gzip decompression is the
//! bottleneck, and a single reader serializing both streams cannot keep up.
//! The framework runs distinct `Serial` steps on distinct workers
//! concurrently (a `Serial` step holds a per-step mutex, not a global one),
//! so the fix is to instantiate **one reader per stream** and join their
//! outputs K-way (see [`ZipRawFastqK`]). Each per-stream reader is built with
//! [`ReadFastqInputs::new_single`], is pinned to its **own** worker
//! (`Affinity::Worker(i)` for stream `i`), and draws its chunks' `ordinal` from
//! its OWN [`FastqOrdinalSequence`]. Each reader is a distinct `ByItemOrdinal`
//! producer edge into the K-way join, so its per-edge ordinal must be dense on
//! its own; cross-stream row alignment is by `chunk_serial`, not by a shared
//! ordinal. Handing all readers one shared sequence would make each edge sparse
//! and stall its reorder stage.
//!
//! The distinct-worker pinning is load-bearing, not a tuning choice: a `Serial`
//! source dispatched by more than one worker hits a `try_lock`-after-`Finished`
//! race in the runtime's source-drain path, so concurrent per-stream readers
//! must never use [`Affinity::None`]. See [`ReadFastqInputs::new_single`].
//!
//! [`ReadFastqInputs::new`] ([`Affinity::Reader`]) retains the original
//! single-reader, all-streams, round-robin model — a single worker-0-sticky
//! thread driving all I/O reads, same as `ReadBgzfBlocks`. The unified K-way
//! decode front no longer wires it (each stream gets its own `new_single`
//! reader), but it is kept for callers that want an all-streams reader.
//!
//! [`ZipRawFastqK`]: super::zip_raw_fastq_k::ZipRawFastqK

use std::collections::VecDeque;
use std::io::{self, BufRead};
use std::sync::Arc;
use std::sync::atomic::{AtomicU64, Ordering};

use parking_lot::Mutex;

use crate::pipeline::core::Unpushed;
use crate::pipeline::core::held::HeldSlot;
use crate::pipeline::core::item::{HeapSize, Ordered};
use crate::pipeline::core::outputs::OrderedBytesSingle;
use crate::pipeline::core::queues::QueueSpec;
use crate::pipeline::core::reorder::BranchOrdering;
use crate::pipeline::core::step::{Affinity, Step, StepCtx, StepKind, StepOutcome, StepProfile};

// ─────────────────────────────────────────────────────────────────────────────
// FastqRawChunk — one chunk of raw (decompressed, unparsed) FASTQ bytes from
// a single stream.
// ─────────────────────────────────────────────────────────────────────────────

/// A monotonic ordinal sequence for one FASTQ producer edge's chunks.
///
/// Every chunk emitted on one producer edge draws its ordinal from a single
/// counter, so the sequence on that edge is both **unique** and **gap-free**
/// (dense). In the per-stream (per-edge) K-way decode front each stream owns
/// its OWN sequence — its edge into `ZipRawFastqK`/`WrapRawFastq1` is
/// `ByItemOrdinal`, so that per-edge ordinal must be dense on its own;
/// cross-stream alignment is by `chunk_serial`, not by a shared ordinal.
///
/// Density is the load-bearing half. `BranchOrdering::ByItemOrdinal` releases
/// only in contiguous ordinal order, so a hole on an edge stalls its reorder
/// stage forever — a silent hang, not an error. Drawing an edge's ordinals from
/// its own single counter (rather than an arithmetic partition like
/// `chunk_serial * n_streams + stream_idx`, which is only dense when every
/// stream yields the same chunk count) makes each edge dense by construction.
/// An unequal-length desync no longer strands an ordinal: each edge stays
/// dense, and `ZipRawFastqK` reports "FASTQ sources out of sync" when the short
/// stream drains with the front row still incomplete.
///
/// The cost is one uncontended `fetch_add` per *chunk* (hundreds of KB of
/// FASTQ), which is far below the noise floor of decompressing that chunk.
#[derive(Debug, Clone, Default)]
pub struct FastqOrdinalSequence(Arc<AtomicU64>);

impl FastqOrdinalSequence {
    /// A fresh sequence starting at 0. Share one clone per reader instance
    /// feeding the same downstream branch.
    #[must_use]
    pub fn new() -> Self {
        Self::default()
    }

    /// Take the next ordinal. `Relaxed` is sufficient: the value only has to be
    /// unique and monotonic, and the reorder stage does its own ordering.
    pub(crate) fn next(&self) -> u64 {
        self.0.fetch_add(1, Ordering::Relaxed)
    }
}

/// A chunk of decompressed, whole-record-aligned raw FASTQ bytes emitted by
/// [`ReadFastqInputs`] (or the BGZF `FindFastqBoundaries`) and consumed by the
/// K-way aligner `ZipRawFastqK` / the 1-way adapter `WrapRawFastq1`.
///
/// `ordinal` is a per-edge, monotonically-increasing counter assigned once per
/// emitted chunk on this stream's edge. It is the key the framework uses to
/// reorder that edge — `stream_idx` and `chunk_serial` are NOT unique across
/// streams (two streams emit chunks with the same `chunk_serial`), so they must
/// not be used as a reorder key. The `(stream_idx, chunk_serial)` pair is
/// preserved so `ZipRawFastqK` can re-join the per-stream chunks by row.
pub struct FastqRawChunk {
    /// Globally-unique monotonic ordinal, for the framework's reorder buffer.
    pub ordinal: u64,
    /// Index of the originating FASTQ stream (0-based).
    pub stream_idx: usize,
    /// Per-stream round-robin cycle serial, for the `ZipRawFastqK` join.
    pub chunk_serial: u64,
    /// Decompressed raw FASTQ bytes, aligned to whole records (4 lines each).
    ///
    /// Rebuild this chunk rather than mutating it in place: an in-place edit
    /// that breaks whole-record (4-line) alignment desyncs downstream parsing.
    pub data: Vec<u8>,
}

impl HeapSize for FastqRawChunk {
    fn heap_size(&self) -> usize {
        // Allocated capacity, not logical length: byte-bounded queues budget the
        // memory actually held, and these buffers are read into with spare
        // capacity. `PairFastqStreams` accounts for the same chunks and must use
        // the same measure, or `pending_total_bytes` drifts against the queues.
        self.data.capacity()
    }
}

impl Ordered for FastqRawChunk {
    fn ordinal(&self) -> u64 {
        self.ordinal
    }
}

// ─────────────────────────────────────────────────────────────────────────────
// FASTQ reading helper
// ─────────────────────────────────────────────────────────────────────────────

/// Read up to `count` complete FASTQ records' worth of raw bytes from
/// `reader` into a single `Vec<u8>`, returning the bytes plus the number of
/// complete records read.
///
/// Each record is 4 lines (`@name`, sequence, `+`, quality); the bytes are
/// always whole-record-aligned (the function never splits mid-record). A
/// returned count of 0 with empty bytes signals EOF for the stream.
///
/// Decompression (when the underlying `BufRead` is a gzip decoder) happens
/// here; FASTQ structure validation and `FastqRecord` parsing are deferred to
/// `ParseAndZipFastqN`. The only structural check performed here is that each
/// record begins with `@`, so the round-robin reader fails fast on a malformed
/// stream rather than emitting garbage to the parse workers.
fn read_fastq_raw_bytes_from_bufread(
    reader: &mut dyn BufRead,
    count: usize,
) -> io::Result<(Vec<u8>, usize)> {
    let mut data: Vec<u8> = Vec::with_capacity(count * 300);
    let mut records_read = 0;

    for _ in 0..count {
        // Line 1: @name
        let line_start = data.len();
        let bytes_read = reader.read_until(b'\n', &mut data)?;
        if bytes_read == 0 {
            break;
        }
        if data[line_start] != b'@' {
            // A lone newline (blank line between records) is a common
            // hand-edit/corruption artifact; give it a clearer message than
            // the raw "got '\n'" the generic branch would produce.
            let message = if data[line_start] == b'\n' {
                "Unexpected blank line in FASTQ stream (expected a record \
                 starting with '@')"
                    .to_string()
            } else {
                format!(
                    "Expected FASTQ record to start with '@', got '{}'",
                    data[line_start] as char
                )
            };
            return Err(io::Error::new(io::ErrorKind::InvalidData, message));
        }
        if data.last() != Some(&b'\n') {
            data.push(b'\n');
        }

        // Line 2: sequence
        let bytes_read = reader.read_until(b'\n', &mut data)?;
        if bytes_read == 0 {
            return Err(io::Error::new(
                io::ErrorKind::UnexpectedEof,
                "Unexpected EOF reading FASTQ sequence line",
            ));
        }
        if data.last() != Some(&b'\n') {
            data.push(b'\n');
        }

        // Line 3: +
        let bytes_read = reader.read_until(b'\n', &mut data)?;
        if bytes_read == 0 {
            return Err(io::Error::new(
                io::ErrorKind::UnexpectedEof,
                "Unexpected EOF reading FASTQ separator line",
            ));
        }
        if data.last() != Some(&b'\n') {
            data.push(b'\n');
        }

        // Line 4: quality
        let bytes_read = reader.read_until(b'\n', &mut data)?;
        if bytes_read == 0 {
            return Err(io::Error::new(
                io::ErrorKind::UnexpectedEof,
                "Unexpected EOF reading FASTQ quality line",
            ));
        }
        if data.last() != Some(&b'\n') {
            data.push(b'\n');
        }

        records_read += 1;
    }

    Ok((data, records_read))
}

// ─────────────────────────────────────────────────────────────────────────────
// ReadFastqInputs — Serial+Reader source step
// ─────────────────────────────────────────────────────────────────────────────

type FastqReaders = Vec<Box<dyn BufRead + Send>>;

pub struct ReadFastqInputs {
    readers: Arc<Mutex<Option<FastqReaders>>>,
    /// Number of streams this instance reads (the length of `readers`). For
    /// `new_single` this is always 1; for `new` it is the total stream count.
    n_streams: usize,
    /// Global index of this instance's first (local) stream. The emitted
    /// `FastqRawChunk::stream_idx` is `stream_idx_base + local_idx`, so a
    /// per-stream reader for R2 reports `stream_idx = 1` even though its only
    /// reader is at local index 0.
    stream_idx_base: usize,
    /// Shared ordinal source. Every reader feeding the same downstream branch
    /// holds a clone, so the emitted ordinals are unique *and* gap-free even
    /// when the streams differ in length. See [`FastqOrdinalSequence`].
    ordinals: FastqOrdinalSequence,
    /// Scheduling hint reported by `affinity()`. `Affinity::Worker(i)` for the
    /// per-stream concurrent readers (each pinned to a distinct worker);
    /// `Affinity::Reader` for the single all-streams fallback.
    affinity: Affinity,
    /// Whether the step advertises `sticky` in its profile. The all-streams
    /// fallback is `sticky` (the legacy worker-0 burst-read model). The
    /// per-stream concurrent readers are NOT sticky: a sticky reader would
    /// monopolize its worker via the scheduler's sticky re-entry loop
    /// (spinning on the reader while it makes Progress), starving the
    /// `Parallel` `BgzfCompress` of that worker between read bursts and
    /// regressing throughput at low `--threads`. Non-sticky lets the
    /// affinity-pinned worker interleave compression work between reads.
    sticky: bool,
    batch_record_count: usize,
    next_chunk_serial: u64,
    pending: VecDeque<FastqRawChunk>,
    held: HeldSlot<Unpushed<FastqRawChunk>>,
    output_byte_limit: u64,
    exhausted: Vec<bool>,
    all_done: bool,
}

impl ReadFastqInputs {
    /// All-streams reader: reads every stream in `readers` round-robin on a
    /// single worker (`Affinity::Reader`). Used for the N≥3 fallback where
    /// the per-stream concurrent model's coordination overhead is not worth
    /// it. `stream_idx_base` is 0, so emitted chunks carry their natural
    /// stream index. This is the only reader on its branch, so it owns a
    /// private [`FastqOrdinalSequence`].
    #[must_use]
    pub fn new(
        readers: Vec<Box<dyn BufRead + Send>>,
        batch_record_count: usize,
        output_byte_limit: u64,
    ) -> Self {
        Self::build(
            readers,
            0,
            FastqOrdinalSequence::new(),
            Affinity::Reader,
            /* sticky */ true,
            batch_record_count,
            output_byte_limit,
        )
    }

    /// Single-stream reader: reads exactly one stream so multiple instances
    /// (e.g. one for R1, one for R2) can run concurrently on different
    /// workers. `global_stream_idx` is this stream's index in the full
    /// pipeline, carried on each chunk so `ZipRawFastqK` can re-join the
    /// streams.
    ///
    /// `ordinals` is this reader's OWN [`FastqOrdinalSequence`]. In the K-way
    /// decode front each single-stream reader is its own producer edge into
    /// `ZipRawFastqK`/`WrapRawFastq1`, and each such edge is `ByItemOrdinal`, so
    /// its per-edge ordinal must be dense (`0, 1, 2, …`) on its own. Cross-stream
    /// row alignment is by `chunk_serial`, not by a shared ordinal — so each
    /// reader gets a fresh sequence, NOT one shared across streams (a shared
    /// counter would make each edge's ordinals sparse and stall its reorder).
    ///
    /// `affinity` pins this reader to a specific worker. Each per-stream
    /// reader is given a **distinct** worker (e.g. `Affinity::Worker(0)` for
    /// R1, `Affinity::Worker(1)` for R2) so the two gzip decoders run
    /// concurrently AND no two workers ever contend the same source mutex.
    /// The latter is load-bearing: a `Serial` source dispatched by more than
    /// one worker hits a `try_lock`-after-`Finished` race in the runtime's
    /// source-drain path, so concurrent readers must use disjoint
    /// single-worker affinities, never `Affinity::None`.
    ///
    #[must_use]
    pub fn new_single(
        reader: Box<dyn BufRead + Send>,
        global_stream_idx: usize,
        ordinals: FastqOrdinalSequence,
        affinity: Affinity,
        batch_record_count: usize,
        output_byte_limit: u64,
    ) -> Self {
        Self::build(
            vec![reader],
            global_stream_idx,
            ordinals,
            affinity,
            /* sticky */ false,
            batch_record_count,
            output_byte_limit,
        )
    }

    fn build(
        readers: Vec<Box<dyn BufRead + Send>>,
        stream_idx_base: usize,
        ordinals: FastqOrdinalSequence,
        affinity: Affinity,
        sticky: bool,
        batch_record_count: usize,
        output_byte_limit: u64,
    ) -> Self {
        // Each reader draws its own edge's ordinals from the `ordinals`
        // sequence it was handed. In the K-way front that is a fresh per-stream
        // sequence (one edge per stream), so each edge's ordinals are dense on
        // their own; in the legacy all-streams `new` reader it is a single
        // private sequence over the round-robin output. Either way the sequence
        // is per-edge, so there is no cross-reader collision to guard. Density
        // (not uniqueness) is the load-bearing half — a sparse edge stalls its
        // reorder into a silent hang. See [`FastqOrdinalSequence`].
        let n_streams = readers.len();

        Self {
            readers: Arc::new(Mutex::new(Some(readers))),
            n_streams,
            stream_idx_base,
            ordinals,
            affinity,
            sticky,
            batch_record_count: batch_record_count.max(1),
            next_chunk_serial: 0,
            pending: VecDeque::new(),
            held: HeldSlot::new(),
            output_byte_limit,
            exhausted: vec![false; n_streams],
            all_done: false,
        }
    }
}

impl Step for ReadFastqInputs {
    type Input = ();
    type Outputs = OrderedBytesSingle<FastqRawChunk>;

    fn profile(&self) -> StepProfile {
        StepProfile {
            name: "ReadFastqInputs",
            kind: StepKind::Serial,
            sticky: self.sticky,
            output_queues: vec![QueueSpec::ByteBounded { limit_bytes: self.output_byte_limit }],
            branch_ordering: vec![BranchOrdering::ByItemOrdinal],
        }
    }

    fn affinity(&self) -> Affinity {
        self.affinity
    }

    fn try_run(&mut self, ctx: &mut StepCtx<'_, Self>) -> io::Result<StepOutcome> {
        // 1. Drain held slot.
        if let Some(unpushed) = self.held.take() {
            match ctx.outputs.retry(unpushed) {
                Ok(()) => {}
                Err(again) => {
                    self.held.put(again);
                    return Ok(StepOutcome::Contention);
                }
            }
        }

        // 2. Drain pending batches.
        if let Some(batch) = self.pending.pop_front() {
            match ctx.outputs.push(batch) {
                Ok(()) => return Ok(StepOutcome::Progress),
                Err(unpushed) => {
                    self.held.put(unpushed);
                    return Ok(StepOutcome::Progress);
                }
            }
        }

        if self.all_done {
            return Ok(StepOutcome::Finished);
        }

        // 3. Read one round-robin cycle: one chunk per stream at the current serial.
        let chunk_serial = self.next_chunk_serial;
        let mut any_read = false;

        {
            let mut guard = self.readers.lock();
            let readers = guard.as_mut().expect("ReadFastqInputs: readers missing");
            debug_assert_eq!(
                readers.len(),
                self.n_streams,
                "reader count must match the declared stream count"
            );

            // Read one chunk from every non-exhausted stream this cycle. Each
            // cycle reads all streams, so there is no rotation offset.
            for (idx, reader) in readers.iter_mut().enumerate() {
                if self.exhausted[idx] {
                    continue;
                }

                let (data, records_read) =
                    read_fastq_raw_bytes_from_bufread(reader.as_mut(), self.batch_record_count)?;

                if records_read == 0 {
                    self.exhausted[idx] = true;
                    continue;
                }

                any_read = true;
                // Global stream index (== local idx for the all-streams
                // reader, which has stream_idx_base == 0).
                let stream_idx = self.stream_idx_base + idx;
                // Drawn per emitted chunk, so a stream that exhausts early
                // leaves no hole in the sequence.
                let ordinal = self.ordinals.next();
                self.pending.push_back(FastqRawChunk { ordinal, stream_idx, chunk_serial, data });
            }
        }

        // Advance to the next round-robin cycle.
        self.next_chunk_serial += 1;

        if !any_read {
            self.all_done = true;
            return Ok(StepOutcome::Finished);
        }

        // Emit one batch from the freshly-filled pending queue.
        if let Some(batch) = self.pending.pop_front() {
            match ctx.outputs.push(batch) {
                Ok(()) => Ok(StepOutcome::Progress),
                Err(unpushed) => {
                    self.held.put(unpushed);
                    Ok(StepOutcome::Progress)
                }
            }
        } else {
            Ok(StepOutcome::NoProgress)
        }
    }
}

#[cfg(test)]
mod tests {
    use rstest::rstest;

    use super::*;
    use crate::pipeline::core::item::Ordered;

    /// Every way `read_fastq_raw_bytes_from_bufread` can reject a malformed
    /// record, one case per rejection. The reader exists to fail fast on
    /// degenerate input, and without these the only covered paths were the happy
    /// path, the count cap, and empty input — so a regression that reordered the
    /// `read_until` / `bytes_read == 0` checks would still have passed.
    ///
    /// Note that the *first* line's `bytes_read == 0` is a clean `break` (end of
    /// stream between records), not an error, so it is covered by the empty-input
    /// test rather than here.
    #[rstest]
    #[case::eof_before_sequence(b"@read1\n".as_slice(), io::ErrorKind::UnexpectedEof, "sequence line")]
    #[case::eof_before_separator(
        b"@read1\nACGT\n".as_slice(),
        io::ErrorKind::UnexpectedEof,
        "separator line"
    )]
    #[case::eof_before_quality(
        b"@read1\nACGT\n+\n".as_slice(),
        io::ErrorKind::UnexpectedEof,
        "quality line"
    )]
    #[case::blank_line_between_records(
        b"\n@read1\nACGT\n+\nIIII\n".as_slice(),
        io::ErrorKind::InvalidData,
        "blank line"
    )]
    #[case::name_line_missing_at_sign(
        b"read1\nACGT\n+\nIIII\n".as_slice(),
        io::ErrorKind::InvalidData,
        "start with '@'"
    )]
    fn read_fastq_raw_bytes_rejects_malformed_records(
        #[case] input: &[u8],
        #[case] expected_kind: io::ErrorKind,
        #[case] expected_message_fragment: &str,
    ) {
        let mut reader = io::Cursor::new(input.to_vec());
        let err = read_fastq_raw_bytes_from_bufread(&mut reader, 2)
            .expect_err("malformed input must be rejected");

        assert_eq!(err.kind(), expected_kind, "wrong error kind for: {err}");
        assert!(
            err.to_string().contains(expected_message_fragment),
            "message {err:?} should mention {expected_message_fragment:?}",
        );
    }

    /// A final quality line with no trailing newline is normalized rather than
    /// rejected, so the emitted chunk stays 4-line aligned for the parser
    /// downstream. This is the boundary case the `data.last() != Some(&b'\n')`
    /// guards exist for, and it sits right next to the EOF rejections above.
    #[test]
    fn read_fastq_raw_bytes_appends_a_missing_final_newline() {
        let mut reader = io::Cursor::new(b"@read1\nACGT\n+\nIIII".to_vec());

        let (bytes, records_read) =
            read_fastq_raw_bytes_from_bufread(&mut reader, 1).expect("record is complete");

        assert_eq!(records_read, 1);
        assert_eq!(bytes, b"@read1\nACGT\n+\nIIII\n");
    }

    #[test]
    fn fastq_raw_chunk_heap_size_and_ordinal() {
        let chunk =
            FastqRawChunk { ordinal: 3, stream_idx: 1, chunk_serial: 0, data: vec![0u8; 42] };
        assert_eq!(chunk.heap_size(), 42);
        assert_eq!(chunk.ordinal(), 3);
    }

    /// `heap_size` must report allocation capacity, not live length: the reader
    /// builds each chunk with `Vec::with_capacity(count * 300)`, so a short
    /// chunk — the final partial batch, or reads under the per-record estimate —
    /// holds materially more memory than it spans.
    ///
    /// The test above cannot pin that: `vec![0u8; 42]` has `capacity() == len()`,
    /// so it passes for either measure and a regression to `len()` would stay
    /// green while the byte-bounded queue under-counted what it holds. Mirrors
    /// `pair_fastq::tests::paired_batch_heap_size_counts_spare_capacity`.
    #[test]
    fn fastq_raw_chunk_heap_size_counts_spare_capacity() {
        let mut data = Vec::with_capacity(300);
        data.extend_from_slice(&[0u8; 42]);
        let capacity = data.capacity();

        let chunk = FastqRawChunk { ordinal: 3, stream_idx: 1, chunk_serial: 0, data };

        assert_eq!(chunk.heap_size(), capacity);
        assert!(
            chunk.heap_size() > 42,
            "heap_size must exceed the 42 live bytes; got {}",
            chunk.heap_size()
        );
    }

    #[test]
    fn profile_advertises_serial_reader_byordinal() {
        let readers: Vec<Box<dyn BufRead + Send>> =
            vec![Box::new(io::Cursor::new(Vec::<u8>::new()))];
        let step = ReadFastqInputs::new(readers, 400, 1024 * 1024);
        let profile = step.profile();
        assert_eq!(profile.name, "ReadFastqInputs");
        assert_eq!(profile.kind, StepKind::Serial);
        assert!(profile.sticky);
        assert_eq!(step.affinity(), Affinity::Reader);
        assert_eq!(profile.branch_ordering, vec![BranchOrdering::ByItemOrdinal]);
        assert!(matches!(profile.output_queues[0], QueueSpec::ByteBounded { .. }));
    }

    /// The per-stream constructor keeps `global_stream_idx` purely as chunk
    /// metadata for `ZipRawFastqK` to re-join on.
    ///
    /// Each per-stream reader is its own producer edge with its own dense
    /// ordinal sequence, so there is no cross-reader collision to guard.
    #[test]
    fn new_single_records_its_global_stream_index() {
        let reader: Box<dyn BufRead + Send> = Box::new(io::Cursor::new(Vec::<u8>::new()));
        let step = ReadFastqInputs::new_single(
            reader,
            1,
            FastqOrdinalSequence::new(),
            Affinity::Worker(1),
            400,
            1024,
        );
        assert_eq!(step.stream_idx_base, 1);
    }

    #[test]
    fn new_single_pins_to_requested_worker() {
        // Per-stream readers pin to a distinct worker so the two gzip
        // decoders run concurrently (R1 on worker 0, R2 on worker 1) without
        // contending the same source mutex.
        let reader: Box<dyn BufRead + Send> = Box::new(io::Cursor::new(Vec::<u8>::new()));
        let step = ReadFastqInputs::new_single(
            reader,
            1,
            FastqOrdinalSequence::new(),
            Affinity::Worker(1),
            400,
            1024,
        );
        assert_eq!(step.affinity(), Affinity::Worker(1));
        assert_eq!(step.stream_idx_base, 1);
        assert_eq!(step.n_streams, 1);
        // Per-stream readers are NOT sticky — a sticky reader would monopolize
        // its worker and starve the parallel compressor.
        assert!(!step.profile().sticky);
    }

    /// The shared sequence must be unique *and* gap-free across every reader
    /// drawing from it — including when one stream exhausts long before
    /// another, which is precisely the case the old strided formula could not
    /// express.
    ///
    /// Density is the half that matters: `BranchOrdering::ByItemOrdinal`
    /// releases only in contiguous order, so a hole stalls the reorder stage
    /// forever rather than erroring.
    #[test]
    fn shared_ordinal_sequence_is_dense_across_readers() {
        let seq = FastqOrdinalSequence::new();
        // Two "readers" sharing one sequence, drawing at wildly different
        // rates — stream 0 keeps going long after stream 1 stops.
        let r1 = seq.clone();
        let r2 = seq.clone();
        let mut minted = Vec::new();
        for cycle in 0..10u64 {
            minted.push(r1.next());
            if cycle < 3 {
                minted.push(r2.next());
            }
        }

        let expected: Vec<u64> = (0..minted.len() as u64).collect();
        let mut sorted = minted.clone();
        sorted.sort_unstable();
        assert_eq!(
            sorted, expected,
            "ordinals must be unique and gap-free however unevenly the readers draw",
        );
    }

    #[test]
    fn read_fastq_raw_bytes_from_bufread_reads_whole_records() {
        let data = b"@read1\nACGT\n+\nIIII\n@read2\nTGCA\n+\nJJJJ\n";
        let mut reader = io::Cursor::new(data.to_vec());
        let (bytes, records_read) = read_fastq_raw_bytes_from_bufread(&mut reader, 10).unwrap();
        assert_eq!(records_read, 2);
        // Raw bytes are returned verbatim (whole-record-aligned).
        assert_eq!(bytes, data);
    }

    #[test]
    fn read_fastq_raw_bytes_from_bufread_respects_count() {
        let data = b"@read1\nACGT\n+\nIIII\n@read2\nTGCA\n+\nJJJJ\n";
        let mut reader = io::Cursor::new(data.to_vec());
        let (bytes, records_read) = read_fastq_raw_bytes_from_bufread(&mut reader, 1).unwrap();
        assert_eq!(records_read, 1);
        assert_eq!(bytes, b"@read1\nACGT\n+\nIIII\n");
    }

    #[test]
    fn read_fastq_raw_bytes_from_bufread_empty_input() {
        let mut reader = io::Cursor::new(Vec::<u8>::new());
        let (bytes, records_read) = read_fastq_raw_bytes_from_bufread(&mut reader, 10).unwrap();
        assert_eq!(records_read, 0);
        assert!(bytes.is_empty());
    }
}
