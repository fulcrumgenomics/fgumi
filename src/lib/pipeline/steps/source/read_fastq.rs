//! `ReadFastqInputs` source step: reads (and gzip-decompresses) raw FASTQ
//! bytes from one or more `BufRead` streams in round-robin order and emits
//! `FastqRawChunk` items. Parsing into `FastqRecord`s is deferred to the
//! downstream `Parallel` `ParseFastqChunks` step so the decode/parse work
//! can fan out across worker threads rather than bottlenecking on the
//! single reader thread.
//!
//! ## Per-stream vs. all-streams instantiation
//!
//! For the common paired-end case (R1 + R2) the gzip decompression is the
//! bottleneck, and a single reader serializing both streams cannot keep up.
//! The framework runs distinct `Serial` steps on distinct workers
//! concurrently (a `Serial` step holds a per-step mutex, not a global one),
//! so the fix is to instantiate **one reader per stream** and pair their
//! outputs 2-way (see `PairRawFastq`). Each per-stream reader is built
//! with [`ReadFastqInputs::new_single`], runs with [`Affinity::None`] so any
//! free worker can drive it, and computes its chunks' globally-unique
//! `ordinal` independently from its own `global_stream_idx` /
//! `n_streams_total` — no shared counter, no cross-reader collision.
//!
//! The N≥3 fallback keeps the original single-reader, all-streams,
//! round-robin model ([`ReadFastqInputs::new`], [`Affinity::Reader`]): a
//! single worker-0-sticky thread drives all I/O reads, same as
//! `ReadBgzfBlocks`.

use std::collections::VecDeque;
use std::io::{self, BufRead};
use std::sync::Arc;

use parking_lot::Mutex;

use crate::fastq_parse::FastqRecord;
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

/// A chunk of decompressed, whole-record-aligned raw FASTQ bytes emitted by
/// [`ReadFastqInputs`] and consumed by `ParseFastqChunks`.
///
/// `ordinal` is a globally-unique, monotonically-increasing counter assigned
/// once per emitted chunk (regardless of stream). It is the key the framework
/// uses to reorder the `Parallel` `ParseFastqChunks` output — `stream_idx` and
/// `chunk_serial` are NOT globally unique (two streams emit chunks with the
/// same `chunk_serial`), so they must not be used as a reorder key. The
/// `(stream_idx, chunk_serial)` pair is preserved purely so `ZipFastqRecords`
/// can re-join the per-stream chunks.
pub struct FastqRawChunk {
    /// Globally-unique monotonic ordinal, for the framework's reorder buffer.
    pub ordinal: u64,
    /// Index of the originating FASTQ stream (0-based).
    pub stream_idx: usize,
    /// Per-stream round-robin cycle serial, for the `ZipFastqRecords` join.
    pub chunk_serial: u64,
    /// Decompressed raw FASTQ bytes, aligned to whole records (4 lines each).
    pub data: Vec<u8>,
}

impl HeapSize for FastqRawChunk {
    fn heap_size(&self) -> usize {
        self.data.len()
    }
}

impl Ordered for FastqRawChunk {
    fn ordinal(&self) -> u64 {
        self.ordinal
    }
}

// ─────────────────────────────────────────────────────────────────────────────
// FastqChunkBatch — one chunk of parsed records from a single FASTQ stream.
// ─────────────────────────────────────────────────────────────────────────────

pub struct FastqChunkBatch {
    /// Globally-unique monotonic ordinal, propagated from the originating
    /// [`FastqRawChunk`]. Used by the framework to reorder the output of the
    /// `Parallel` `ParseFastqChunks` step. Distinct from `chunk_serial`,
    /// which is only per-stream unique and is reused across streams.
    pub ordinal: u64,
    pub stream_idx: usize,
    pub chunk_serial: u64,
    pub records: Vec<FastqRecord>,
    total_bytes: usize,
}

impl FastqChunkBatch {
    #[must_use]
    pub fn new(
        ordinal: u64,
        stream_idx: usize,
        chunk_serial: u64,
        records: Vec<FastqRecord>,
        total_bytes: usize,
    ) -> Self {
        Self { ordinal, stream_idx, chunk_serial, records, total_bytes }
    }
}

impl HeapSize for FastqChunkBatch {
    fn heap_size(&self) -> usize {
        self.total_bytes
    }
}

impl Ordered for FastqChunkBatch {
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
/// `ParseFastqChunks`. The only structural check performed here is that each
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
            return Err(io::Error::new(
                io::ErrorKind::InvalidData,
                format!(
                    "Expected FASTQ record to start with '@', got '{}'",
                    data[line_start] as char
                ),
            ));
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
    /// Total number of FASTQ streams across the whole pipeline (all
    /// concurrent readers combined). Used in the ordinal formula so every
    /// chunk has a globally-unique ordinal regardless of which reader
    /// produced it.
    n_streams_total: usize,
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
    current_stream: usize,
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
    /// it. `stream_idx_base` is 0 and `n_streams_total` equals the number of
    /// readers, so emitted chunks carry their natural stream index and a
    /// `chunk_serial`-derived globally-unique ordinal.
    #[must_use]
    pub fn new(
        readers: Vec<Box<dyn BufRead + Send>>,
        batch_record_count: usize,
        output_byte_limit: u64,
    ) -> Self {
        let n_streams = readers.len();
        Self::build(
            readers,
            0,
            n_streams,
            Affinity::Reader,
            /* sticky */ true,
            batch_record_count,
            output_byte_limit,
        )
    }

    /// Single-stream reader: reads exactly one stream so multiple instances
    /// (e.g. one for R1, one for R2) can run concurrently on different
    /// workers. `global_stream_idx` is this stream's index in the full
    /// pipeline and `n_streams_total` is the total stream count; together they
    /// make this reader's chunk ordinals globally unique without any shared
    /// counter:
    /// `ordinal = chunk_serial * n_streams_total + global_stream_idx`.
    ///
    /// `affinity` pins this reader to a specific worker. Each per-stream
    /// reader is given a **distinct** worker (e.g. `Affinity::Worker(0)` for
    /// R1, `Affinity::Worker(1)` for R2) so the two gzip decoders run
    /// concurrently AND no two workers ever contend the same source mutex.
    /// The latter is load-bearing: a `Serial` source dispatched by more than
    /// one worker hits a `try_lock`-after-`Finished` race in the runtime's
    /// source-drain path, so concurrent readers must use disjoint
    /// single-worker affinities, never `Affinity::None`.
    #[must_use]
    pub fn new_single(
        reader: Box<dyn BufRead + Send>,
        global_stream_idx: usize,
        n_streams_total: usize,
        affinity: Affinity,
        batch_record_count: usize,
        output_byte_limit: u64,
    ) -> Self {
        Self::build(
            vec![reader],
            global_stream_idx,
            n_streams_total,
            affinity,
            /* sticky */ false,
            batch_record_count,
            output_byte_limit,
        )
    }

    fn build(
        readers: Vec<Box<dyn BufRead + Send>>,
        stream_idx_base: usize,
        n_streams_total: usize,
        affinity: Affinity,
        sticky: bool,
        batch_record_count: usize,
        output_byte_limit: u64,
    ) -> Self {
        let n_streams = readers.len();
        Self {
            readers: Arc::new(Mutex::new(Some(readers))),
            n_streams,
            stream_idx_base,
            n_streams_total: n_streams_total.max(1),
            affinity,
            sticky,
            batch_record_count: batch_record_count.max(1),
            next_chunk_serial: 0,
            current_stream: 0,
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

            // Read from each stream starting at current_stream, wrapping around.
            for offset in 0..self.n_streams {
                let idx = (self.current_stream + offset) % self.n_streams;
                if self.exhausted[idx] {
                    continue;
                }

                let (data, records_read) = read_fastq_raw_bytes_from_bufread(
                    readers[idx].as_mut(),
                    self.batch_record_count,
                )?;

                if records_read == 0 {
                    self.exhausted[idx] = true;
                    continue;
                }

                any_read = true;
                // Global stream index (== local idx for the all-streams
                // reader, which has stream_idx_base == 0).
                let stream_idx = self.stream_idx_base + idx;
                // Globally-unique across every (chunk_serial, stream) pair:
                // for a fixed chunk_serial the streams occupy
                // [serial*N, serial*N + N), disjoint from the next serial's
                // window. For N==1 this reduces to ordinal == chunk_serial.
                let ordinal = chunk_serial * self.n_streams_total as u64 + stream_idx as u64;
                self.pending.push_back(FastqRawChunk { ordinal, stream_idx, chunk_serial, data });
            }
        }

        // Advance round-robin and serial.
        self.current_stream = 0;
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
    use super::*;
    use crate::pipeline::core::item::Ordered;

    #[test]
    fn fastq_chunk_batch_heap_size_returns_total_bytes() {
        // ordinal (7) and chunk_serial (5) are intentionally different to
        // confirm `ordinal()` reports the unique ordinal, not chunk_serial.
        let batch = FastqChunkBatch::new(7, 0, 5, vec![], 1234);
        assert_eq!(batch.heap_size(), 1234);
        assert_eq!(batch.ordinal(), 7);
    }

    #[test]
    fn fastq_raw_chunk_heap_size_and_ordinal() {
        let chunk =
            FastqRawChunk { ordinal: 3, stream_idx: 1, chunk_serial: 0, data: vec![0u8; 42] };
        assert_eq!(chunk.heap_size(), 42);
        assert_eq!(chunk.ordinal(), 3);
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

    #[test]
    fn new_single_pins_to_requested_worker() {
        // Per-stream readers pin to a distinct worker so the two gzip
        // decoders run concurrently (R1 on worker 0, R2 on worker 1) without
        // contending the same source mutex.
        let reader: Box<dyn BufRead + Send> = Box::new(io::Cursor::new(Vec::<u8>::new()));
        let step = ReadFastqInputs::new_single(reader, 1, 2, Affinity::Worker(1), 400, 1024);
        assert_eq!(step.affinity(), Affinity::Worker(1));
        assert_eq!(step.stream_idx_base, 1);
        assert_eq!(step.n_streams_total, 2);
        assert_eq!(step.n_streams, 1);
        // Per-stream readers are NOT sticky — a sticky reader would monopolize
        // its worker and starve the parallel compressor.
        assert!(!step.profile().sticky);
    }

    /// The ordinal formula must produce a globally-unique value for every
    /// `(chunk_serial, global_stream_idx)` pair so the `Parallel`
    /// `ParseFastqChunks` reorder buffer (keyed on `ordinal`) never collides
    /// across the concurrent per-stream readers. This pins the formula and
    /// its key invariants: uniqueness across streams within a serial, and
    /// reduction to `ordinal == chunk_serial` for the single-stream (N==1)
    /// case.
    #[test]
    fn ordinal_formula_is_globally_unique() {
        let ordinal = |chunk_serial: u64, stream_idx: usize, n_total: usize| -> u64 {
            chunk_serial * n_total as u64 + stream_idx as u64
        };

        // N == 1: ordinal collapses to chunk_serial.
        for s in 0..5 {
            assert_eq!(ordinal(s, 0, 1), s);
        }

        // N == 2 (R1 + R2): R1 and R2 at the same chunk_serial differ, and
        // no ordinal repeats across the full (serial, stream) grid.
        let mut seen = std::collections::HashSet::new();
        for serial in 0..100u64 {
            let r1 = ordinal(serial, 0, 2);
            let r2 = ordinal(serial, 1, 2);
            assert_ne!(r1, r2, "R1/R2 ordinals collided at serial {serial}");
            assert!(seen.insert(r1), "duplicate ordinal {r1}");
            assert!(seen.insert(r2), "duplicate ordinal {r2}");
        }

        // N == 3 fallback: all three streams across many serials stay unique.
        let mut seen3 = std::collections::HashSet::new();
        for serial in 0..100u64 {
            for stream_idx in 0..3 {
                assert!(
                    seen3.insert(ordinal(serial, stream_idx, 3)),
                    "duplicate ordinal at serial {serial} stream {stream_idx}"
                );
            }
        }
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
