//! `FindFastqBoundaries` — `Serial` per-stream record-seam resolver for the
//! parallel BGZF FASTQ decode split.
//!
//! The `Parallel` [`FastqDecompress`](super::fastq_bgzf::FastqDecompress) step
//! emits [`FastqDecompressedBlock`](super::fastq_bgzf::FastqDecompressedBlock)s
//! whose bytes end wherever the compressor cut the block — almost never on a
//! FASTQ record (4-line) boundary. This step is the FASTQ analogue of
//! `FindBamBoundaries`: it re-frames the concatenated decompressed stream into
//! whole-record chunks and emits the existing
//! [`FastqRawChunk`](super::read_fastq::FastqRawChunk) type, so everything
//! downstream (`ZipRawFastqK` / `WrapRawFastq1` / `ParseAndZipFastqN`) is
//! unchanged.
//!
//! ## Fixed-record-count chunks (block-layout independence)
//!
//! The chunks are cut at a **fixed record count** (`batch_record_count`), NOT
//! one-per-block. This is load-bearing for the K-way join: `ZipRawFastqK`
//! matches every stream's chunk *N* by `chunk_serial` and assumes they all hold
//! the same record range. Block boundaries follow
//! *compression*, which differs between R1 and R2 (e.g. R1 in one BGZF block,
//! R2 one-record-per-block), so emitting one chunk per block would desync the
//! join. Accumulating decompressed bytes and cutting every `batch_record_count`
//! records reproduces the fused
//! [`ReadFastqInputs`](super::read_fastq::ReadFastqInputs) reader's chunk cadence
//! (which reads a fixed record count per chunk), keeping chunk *N* = records
//! `[N·batch, N·batch + batch)` on both streams regardless of block layout.
//!
//! ## Ordinals
//!
//! The upstream per-stream block ordinal (on `FastqDecompressedBlock`) is only
//! the reorder key for the parallel decompress; it is **not** the dense chunk
//! ordinal the downstream reorder stage needs. So this step mints the emitted
//! `FastqRawChunk.ordinal` from its own per-stream
//! [`FastqOrdinalSequence`](super::read_fastq::FastqOrdinalSequence). Each
//! stream is its own producer edge into `ZipRawFastqK`/`WrapRawFastq1`, and that
//! edge is `ByItemOrdinal`, so its ordinals must be dense (`0, 1, 2, …`) on
//! their own — hence a fresh sequence per stream, not one shared across
//! streams. `chunk_serial` is the per-stream monotonic cycle counter that
//! aligns rows ACROSS streams in the `ZipRawFastqK` join.
//!
//! Record seams are found with [`fgumi_simd_fastq::find_record_offsets`], whose
//! offsets mark whole-record (4-line) boundaries; the emitted chunk bytes
//! therefore always hold a whole number of records, the invariant
//! `parse_fastq_chunk` relies on.

use std::io;

use fgumi_simd_fastq::find_record_offsets;

use super::fastq_bgzf::FastqDecompressedBlock;
use super::read_fastq::{FastqOrdinalSequence, FastqRawChunk};
use crate::pipeline::core::Unpushed;
use crate::pipeline::core::held::HeldSlot;
use crate::pipeline::core::outputs::OrderedBytesSingle;
use crate::pipeline::core::queues::QueueSpec;
use crate::pipeline::core::reorder::BranchOrdering;
use crate::pipeline::core::step::{Step, StepCtx, StepKind, StepOutcome, StepProfile};

/// Max decompressed blocks consumed per `try_run`, amortizing the `Serial`
/// mutex acquisition. Mirrors `FindBamBoundaries::MAX_BATCHES_PER_LOCK`.
const MAX_BLOCKS_PER_LOCK: usize = 8;

/// `Serial` per-stream FASTQ record-seam resolver. Accumulates decompressed
/// bytes and emits fixed-record-count chunks (see the module docs for why the
/// count — not the block boundary — sets the cut).
pub struct FindFastqBoundaries {
    /// Global stream index (0 = R1, 1 = R2, …), stamped on every emitted chunk
    /// so `ZipRawFastqK` can re-join the streams.
    stream_idx: usize,
    /// Records per emitted chunk. Matches the fused reader's `batch_record_count`
    /// so both paths — and both streams — cut chunks at the same record cadence.
    batch_record_count: usize,
    /// This stream's OWN dense ordinal source. Each stream is a distinct
    /// producer edge into the K-way join, and that edge is `ByItemOrdinal`, so
    /// its ordinals must be dense on their own (a hole hangs the reorder stage);
    /// the sequence is therefore per-stream, not shared across streams.
    ordinals: FastqOrdinalSequence,
    /// Per-stream monotonic cycle serial for the `ZipRawFastqK` join.
    next_chunk_serial: u64,
    /// Accumulated decompressed bytes not yet emitted: a run of complete records
    /// (at the front) plus a possibly-incomplete trailing record.
    acc: Vec<u8>,
    /// Cached record-boundary table over `acc` (`find_record_offsets(&acc)`):
    /// `[0, .., e_last]`, `e_last` one byte past the last complete record.
    /// Refreshed on every mutation of `acc` (`ingest` / `take_full_chunk`) so the
    /// hot path scans `acc` once per change rather than 2–3× per emitted chunk.
    offsets: Vec<usize>,
    /// Pending emit held after a rejected push.
    held: HeldSlot<Unpushed<FastqRawChunk>>,
    output_byte_limit: u64,
    finished: bool,
}

impl FindFastqBoundaries {
    /// Construct a per-stream boundary resolver emitting `batch_record_count`
    /// records per chunk.
    ///
    /// `ordinals` is this stream's OWN [`FastqOrdinalSequence`] (a fresh
    /// `FastqOrdinalSequence::new()` per stream). Each stream is a distinct
    /// `ByItemOrdinal` producer edge into the K-way join, so its ordinals must
    /// be dense on their own; a sequence shared across streams would make each
    /// edge sparse and hang its reorder stage.
    #[must_use]
    pub fn new(
        stream_idx: usize,
        batch_record_count: usize,
        ordinals: FastqOrdinalSequence,
        output_byte_limit: u64,
    ) -> Self {
        Self {
            stream_idx,
            batch_record_count: batch_record_count.max(1),
            ordinals,
            next_chunk_serial: 0,
            acc: Vec::new(),
            offsets: vec![0],
            held: HeldSlot::new(),
            output_byte_limit,
            finished: false,
        }
    }

    /// Ingest one decompressed block: append its bytes to `acc` and rescan the
    /// complete-record boundaries over the accumulator.
    ///
    /// The caller (`try_run` / `flush_full_chunks`) emits every full chunk before
    /// the next ingest, so `acc` holds fewer than `batch_record_count` complete
    /// records (the residue below one chunk) plus at most one `MAX_BLOCKS_PER_LOCK`
    /// burst of freshly-appended blocks. It is therefore bounded by config
    /// (`batch_record_count`) and the block/burst size, never by the input length
    /// — the per-stage-memory-is-a-function-of-config rule. The scan is cached in
    /// `offsets` so `take_full_chunk`/`records_ready` reuse it without re-scanning.
    fn ingest(&mut self, block_data: &[u8]) {
        self.acc.extend_from_slice(block_data);
        self.rescan();
    }

    /// Recompute `offsets` (the cached record-boundary table) over the current
    /// `acc`. `offsets = [0, e1, .., e_last]` where `e_last` is one byte past the
    /// final COMPLETE record; `acc[e_last..]` is the trailing partial record (or
    /// empty). `offsets.len() - 1` is the complete-record count.
    fn rescan(&mut self) {
        self.offsets = find_record_offsets(&self.acc);
    }

    /// Number of complete records currently buffered, from the cached scan.
    fn records_ready(&self) -> usize {
        self.offsets.len().saturating_sub(1)
    }

    /// If at least `batch_record_count` complete records are buffered, cut the
    /// first `batch_record_count` off the front and return their bytes.
    ///
    /// Reuses the cached `offsets` for the cut, then rescans the (much smaller)
    /// remainder exactly once. Returns `None` when fewer than a full chunk's
    /// records are ready (the caller waits for more blocks, or flushes at drain).
    fn take_full_chunk(&mut self) -> Option<Vec<u8>> {
        if self.records_ready() < self.batch_record_count {
            return None;
        }
        // offsets[batch_record_count] is one byte past the Nth complete record.
        let cut = self.offsets[self.batch_record_count];
        let chunk: Vec<u8> = self.acc[..cut].to_vec();
        // Retain the remainder (complete records beyond the cut + trailing
        // partial) and rescan just that tail.
        self.acc.drain(..cut);
        self.rescan();
        Some(chunk)
    }

    /// Consume the final `acc` remainder (< `batch_record_count` records) at
    /// drain, returning whole-record-aligned bytes or an error if the tail is a
    /// genuinely truncated record.
    ///
    /// Only ONE normalization is legitimate here: a FASTQ file whose last record
    /// omits its trailing `\n`. That record has all four lines but only three
    /// newlines, so the 4-newline scan leaves it in `acc`; appending the missing
    /// `\n` completes it. Any OTHER residue — a record cut mid-line, or fewer
    /// than four lines — is a truncated stream, and the split must fail closed
    /// like the fused reader's `read_fastq_raw_bytes_from_bufread` (which errors
    /// `UnexpectedEof` on a short record) rather than silently synthesize a
    /// malformed record by tacking on a newline. Papering over truncation is
    /// worse than no finding: it emits corrupt records downstream.
    fn finalize_remainder(&mut self) -> io::Result<Vec<u8>> {
        let mut data = std::mem::take(&mut self.acc);
        self.offsets = vec![0];
        if data.last() != Some(&b'\n') {
            data.push(b'\n');
        }
        // After normalization the bytes must form a WHOLE number of records: the
        // last boundary equals the byte length. If not, the tail was truncated
        // mid-record (not merely missing a final newline) — surface it.
        let offs = find_record_offsets(&data);
        let last = offs.last().copied().unwrap_or(0);
        if last != data.len() {
            return Err(io::Error::new(
                io::ErrorKind::UnexpectedEof,
                format!(
                    "FindFastqBoundaries: truncated FASTQ record at end of stream \
                     ({} trailing byte(s) after the last complete record)",
                    data.len() - last
                ),
            ));
        }
        Ok(data)
    }

    /// Package `data` (whole-record-aligned bytes) into a `FastqRawChunk` with a
    /// freshly-minted per-stream (per-edge) ordinal and per-stream chunk serial.
    fn make_chunk(&mut self, data: Vec<u8>) -> FastqRawChunk {
        let ordinal = self.ordinals.next();
        let chunk_serial = self.next_chunk_serial;
        self.next_chunk_serial += 1;
        FastqRawChunk { ordinal, stream_idx: self.stream_idx, chunk_serial, data }
    }

    /// Push a chunk, holding it on rejection. Returns `true` if a chunk was
    /// produced (emitted or held) this call.
    fn emit(&mut self, data: Vec<u8>, ctx: &mut StepCtx<'_, Self>) {
        let chunk = self.make_chunk(data);
        if let Err(unpushed) = ctx.outputs.push(chunk) {
            self.held.put(unpushed);
        }
    }
}

impl Step for FindFastqBoundaries {
    type Input = FastqDecompressedBlock;
    type Outputs = OrderedBytesSingle<FastqRawChunk>;

    fn profile(&self) -> StepProfile {
        StepProfile {
            name: "FindFastqBoundaries",
            kind: StepKind::Serial,
            sticky: false,
            output_queues: vec![QueueSpec::ByteBounded { limit_bytes: self.output_byte_limit }],
            branch_ordering: vec![BranchOrdering::ByItemOrdinal],
        }
    }

    fn try_run(&mut self, ctx: &mut StepCtx<'_, Self>) -> io::Result<StepOutcome> {
        // 1. Drain a held emit first.
        if let Some(unpushed) = self.held.take() {
            match ctx.outputs.retry(unpushed) {
                Ok(()) => {}
                Err(again) => {
                    self.held.put(again);
                    return Ok(StepOutcome::Contention);
                }
            }
        }

        // 2. Emit a full chunk if one is already buffered, before taking more
        // input (bounds `acc` to ~one chunk + one block).
        if let Some(chunk) = self.take_full_chunk() {
            self.emit(chunk, ctx);
            return Ok(StepOutcome::Progress);
        }

        // 3. Ingest up to MAX_BLOCKS_PER_LOCK decompressed blocks, emitting a
        // chunk as soon as a full one accumulates.
        let mut did_work = false;
        for _ in 0..MAX_BLOCKS_PER_LOCK {
            let Some(block) = ctx.input.pop() else { break };
            did_work = true;
            self.ingest(&block.data);
            if let Some(chunk) = self.take_full_chunk() {
                self.emit(chunk, ctx);
                return Ok(StepOutcome::Progress);
            }
        }

        if did_work {
            return Ok(StepOutcome::Progress);
        }

        // 4. No input this call. On drained input, first drain any remaining
        // FULL chunks (step 2 already emitted the head, but `acc` may still hold
        // ≥ batch_record_count records if the last ingest overshot), one per
        // call — the held slot holds a single item, so multiple full chunks must
        // span several `try_run` calls. Emitting them full (not as one giant
        // chunk) keeps R1 and R2 cutting at the SAME record cadence, which the
        // chunk_serial join requires: a stream that dumped its residue as one
        // oversized final chunk would desync against a stream that cut it into
        // full chunks + a short remainder.
        if ctx.input.is_drained() {
            if let Some(chunk) = self.take_full_chunk() {
                self.emit(chunk, ctx);
                return Ok(StepOutcome::Progress);
            }
            // Then flush the short final remainder (< batch_record_count records)
            // exactly once. `finished` guards against a second flush across the
            // multi-pass completion drain.
            if !self.finished {
                self.finished = true;
                if !self.acc.is_empty() {
                    let data = self.finalize_remainder()?;
                    self.emit(data, ctx);
                    return Ok(StepOutcome::Progress);
                }
            }
            return Ok(StepOutcome::Finished);
        }
        Ok(StepOutcome::NoProgress)
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::pipeline::core::item::Ordered;

    fn rec(name: &str) -> Vec<u8> {
        format!("@{name}\nACGT\n+\nIIII\n").into_bytes()
    }

    #[test]
    fn profile_is_serial_byordinal() {
        let s = FindFastqBoundaries::new(0, 400, FastqOrdinalSequence::new(), 1 << 20);
        let p = s.profile();
        assert_eq!(p.name, "FindFastqBoundaries");
        assert_eq!(p.kind, StepKind::Serial);
        assert_eq!(p.branch_ordering, vec![BranchOrdering::ByItemOrdinal]);
    }

    /// Chunks cut at a fixed record count regardless of how blocks are fed:
    /// two records per chunk even when records arrive one block at a time.
    #[test]
    fn cuts_at_fixed_record_count_across_blocks() {
        let mut s = FindFastqBoundaries::new(0, 2, FastqOrdinalSequence::new(), 1 << 20);

        // Feed 5 records one "block" at a time; a full chunk (2 records) should
        // become available after the 2nd and 4th records.
        s.ingest(&rec("r0"));
        assert!(s.take_full_chunk().is_none(), "1 record < batch of 2");
        s.ingest(&rec("r1"));
        let c0 = s.take_full_chunk().expect("2 records ready");
        assert_eq!(c0, [rec("r0"), rec("r1")].concat());

        s.ingest(&rec("r2"));
        s.ingest(&rec("r3"));
        let c1 = s.take_full_chunk().expect("2 more ready");
        assert_eq!(c1, [rec("r2"), rec("r3")].concat());

        // One record remains; not a full chunk.
        s.ingest(&rec("r4"));
        assert!(s.take_full_chunk().is_none());
        // The remainder is exactly r4, whole-record-aligned, ready for the
        // drain-time flush.
        assert_eq!(s.acc, rec("r4"));
    }

    /// A block ending mid-record carries the partial forward; the next block
    /// completes it with no record split.
    #[test]
    fn carries_partial_record_across_blocks() {
        let mut s = FindFastqBoundaries::new(0, 1, FastqOrdinalSequence::new(), 1 << 20);

        // First block: one whole record + the first half of a second.
        let mut b1 = rec("r0");
        b1.extend_from_slice(b"@r1\nAC");
        s.ingest(&b1);
        let c0 = s.take_full_chunk().expect("one complete record");
        assert_eq!(c0, rec("r0"));
        assert!(s.take_full_chunk().is_none(), "r1 still partial");

        // Second block completes r1.
        s.ingest(b"GT\n+\nII\n");
        let c1 = s.take_full_chunk().expect("r1 now complete");
        assert_eq!(c1, b"@r1\nACGT\n+\nII\n".to_vec());
    }

    /// `finalize_remainder` normalizes a final record missing only its trailing
    /// `\n` (all four lines present, three newlines) into a whole-record chunk —
    /// the one legitimate normalization, mirroring the fused reader.
    #[test]
    fn finalize_remainder_normalizes_missing_final_newline() {
        let mut s = FindFastqBoundaries::new(0, 400, FastqOrdinalSequence::new(), 1 << 20);
        s.ingest(b"@r0\nACGT\n+\nIIII"); // no trailing newline
        assert!(s.take_full_chunk().is_none(), "not a complete record by the 4-newline scan");
        let out = s.finalize_remainder().expect("a record missing only its final newline is fine");
        assert_eq!(out, b"@r0\nACGT\n+\nIIII\n");
        // Whole-record aligned after normalization.
        let offs = find_record_offsets(&out);
        assert_eq!(*offs.last().unwrap(), out.len());
    }

    /// `finalize_remainder` FAILS CLOSED on a genuinely truncated tail (a record
    /// cut mid-line, or with fewer than four lines): appending one newline does
    /// not make it a whole record, so it must surface as `UnexpectedEof` rather
    /// than being silently emitted as a malformed record. This is the fail-fast
    /// parity with the fused reader's `read_fastq_raw_bytes_from_bufread`.
    #[test]
    fn finalize_remainder_errors_on_truncated_tail() {
        // Two whole records + a third truncated after only its name+seq lines.
        for truncated in [
            &b"@r0\nACGT\n+\nIIII\n@r1\nAC"[..], // cut mid-record (2 lines only)
            &b"@r0\nACGT\n+"[..],                // cut mid-record, no quality
        ] {
            let mut s = FindFastqBoundaries::new(0, 400, FastqOrdinalSequence::new(), 1 << 20);
            s.ingest(truncated);
            let err = s
                .finalize_remainder()
                .expect_err("a mid-record truncation must fail closed, not normalize");
            assert_eq!(err.kind(), io::ErrorKind::UnexpectedEof, "got: {err}");
            assert!(err.to_string().contains("truncated"), "got: {err}");
        }
    }

    /// The residue below one full chunk stays whole-record-aligned in `acc`
    /// (the trailing partial, if any, is carried; complete records sit at the
    /// front). Pins that the cached `offsets` and `acc` stay consistent after a
    /// cut so the residue is exactly the un-emitted records.
    #[test]
    fn residue_after_cut_is_the_unemitted_records() {
        let mut s = FindFastqBoundaries::new(0, 2, FastqOrdinalSequence::new(), 1 << 20);
        // 5 complete records in one ingest.
        let all: Vec<u8> = (0..5).flat_map(|i| rec(&format!("r{i}"))).collect();
        s.ingest(&all);
        let c0 = s.take_full_chunk().expect("first 2");
        assert_eq!(c0, [rec("r0"), rec("r1")].concat());
        let c1 = s.take_full_chunk().expect("next 2");
        assert_eq!(c1, [rec("r2"), rec("r3")].concat());
        assert!(s.take_full_chunk().is_none(), "only 1 record left, < batch of 2");
        assert_eq!(s.acc, rec("r4"), "residue is exactly the un-emitted record");
        assert_eq!(s.records_ready(), 1, "cached offsets agree with acc");
    }

    /// Each stream owns its OWN ordinal sequence (its edge into the K-way join
    /// is `ByItemOrdinal`, so its per-edge ordinals must be dense on their own).
    /// Both streams therefore mint `0, 1, 2, …` independently; cross-stream row
    /// alignment is by `chunk_serial`, not by the ordinal.
    #[test]
    fn each_stream_mints_its_own_dense_ordinal() {
        let mut r1 = FindFastqBoundaries::new(0, 400, FastqOrdinalSequence::new(), 1 << 20);
        let mut r2 = FindFastqBoundaries::new(1, 400, FastqOrdinalSequence::new(), 1 << 20);

        // Stream 0: rows 0, 1. Stream 1: row 0.
        let a0 = r1.make_chunk(rec("a"));
        let b0 = r2.make_chunk(rec("b"));
        let a1 = r1.make_chunk(rec("c"));

        // Per-edge ordinals are dense per stream (not shared/interleaved).
        assert_eq!(a0.ordinal(), 0);
        assert_eq!(a1.ordinal(), 1);
        assert_eq!(b0.ordinal(), 0, "stream 1 has its own sequence starting at 0");

        // chunk_serial is the per-stream row index that aligns across streams.
        assert_eq!(a0.chunk_serial, 0);
        assert_eq!(a1.chunk_serial, 1);
        assert_eq!(b0.chunk_serial, 0);
        assert_eq!(a0.stream_idx, 0);
        assert_eq!(b0.stream_idx, 1);
    }
}
