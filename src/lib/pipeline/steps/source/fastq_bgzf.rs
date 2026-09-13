//! Parallel BGZF FASTQ decode front — the three-step split that mirrors the
//! BAM path (`ReadBgzfBlocks → BgzfDecompress → FindBamBoundaries`) for
//! **bgzip-compressed** FASTQ input.
//!
//! The fused [`ReadFastqInputs`](super::read_fastq::ReadFastqInputs) step reads
//! *and* gzip-decodes in one `Serial` step, so decode caps at one thread per
//! stream (two for a paired R1/R2 run). For a plain-gzip stream that is the only
//! parallelism available — a single continuous DEFLATE stream cannot be
//! block-parallel-decoded. But **BGZF** is block-structured: each block is an
//! independent DEFLATE stream, so decode can fan across the whole pool exactly
//! as the BAM path does.
//!
//! This module provides the pieces for that split, one per stream:
//!
//! ```text
//! ReadFastqBlocks (Serial, per stream)   → FastqRawBlock          (raw compressed block)
//!   → FastqDecompress (Parallel, pooled)  → FastqDecompressedBlock (inflated bytes)
//!   → FindFastqBoundaries (Serial)        → FastqRawChunk          (4-line-aligned records)
//! ```
//!
//! [`FindFastqBoundaries`](super::find_fastq_boundaries::FindFastqBoundaries)
//! lives in its own module (it carries the cross-block record-seam state, the
//! analogue of `FindBamBoundaries`); the raw-block reader and the parallel
//! decompressor live here.
//!
//! ## Why this only applies to BGZF
//!
//! [`ReadFastqBlocks`] frames the input with
//! [`fgumi_bgzf::reader::read_raw_blocks`], which requires a valid BGZF block
//! header on every block. Plain gzip has no such framing, so the chain builder
//! only routes bgzip-detected inputs here; plain-gzip / uncompressed inputs keep
//! the fused [`ReadFastqInputs`](super::read_fastq::ReadFastqInputs) path.

use std::collections::VecDeque;
use std::io::{self, Read};

use fgumi_bgzf::reader::{decompress_block_slice_into_opts, read_raw_blocks};
use libdeflater::Decompressor;

use crate::pipeline::core::Unpushed;
use crate::pipeline::core::held::HeldSlot;
use crate::pipeline::core::item::{HeapSize, Ordered};
use crate::pipeline::core::outputs::OrderedBytesSingle;
use crate::pipeline::core::queues::QueueSpec;
use crate::pipeline::core::reorder::BranchOrdering;
use crate::pipeline::core::step::{Affinity, Step, StepCtx, StepKind, StepOutcome, StepProfile};

/// Max raw BGZF blocks read per `ReadFastqBlocks::try_run`. Amortizes the
/// `Serial` mutex acquisition, mirroring `read_bam::DEFAULT_BLOCKS_PER_BATCH`.
const BLOCKS_PER_BATCH: usize = 16;

/// Per-worker decompression scratch capacity. Matches `BgzfDecompress`
/// (256 KiB ≈ 4 blocks × 64 KiB) so freed buffers land on one mimalloc size
/// class and the thread-local cache can recycle them.
const DECOMPRESS_SCRATCH_CAPACITY: usize = 256 * 1024;

// ─────────────────────────────────────────────────────────────────────────────
// FastqRawBlock — one raw (still-compressed) BGZF block from one FASTQ stream.
// ─────────────────────────────────────────────────────────────────────────────

/// A raw BGZF block read from a single FASTQ stream, not yet decompressed.
///
/// `ordinal` is a per-stream monotonic block counter used purely as the
/// framework's reorder key for the `Parallel` [`FastqDecompress`] step (so the
/// pool can inflate blocks out of order and the reorder stage restores order).
/// It is **not** the gap-free chunk ordinal the downstream
/// [`FastqRawChunk`](super::read_fastq::FastqRawChunk) needs — that is minted
/// later by [`FindFastqBoundaries`](super::find_fastq_boundaries::FindFastqBoundaries)
/// from the shared
/// [`FastqOrdinalSequence`](super::read_fastq::FastqOrdinalSequence), because
/// only the serial boundary step sees the final per-stream chunk stream.
///
/// `stream_idx` rides along so the downstream join (`ZipRawFastqK`) can
/// re-associate the streams by row, exactly as it does for the fused reader's
/// chunks.
pub struct FastqRawBlock {
    /// Per-stream monotonic ordinal, the reorder key for `FastqDecompress`.
    pub ordinal: u64,
    /// Index of the originating FASTQ stream (0-based; 0 = R1, 1 = R2).
    pub stream_idx: usize,
    /// Complete raw BGZF block bytes (header + compressed + footer).
    pub data: Vec<u8>,
}

impl HeapSize for FastqRawBlock {
    fn heap_size(&self) -> usize {
        self.data.capacity()
    }
}

impl Ordered for FastqRawBlock {
    fn ordinal(&self) -> u64 {
        self.ordinal
    }
}

// ─────────────────────────────────────────────────────────────────────────────
// FastqDecompressedBlock — one inflated BGZF block from one FASTQ stream.
// ─────────────────────────────────────────────────────────────────────────────

/// A decompressed BGZF block from one FASTQ stream. The bytes are *not*
/// record-aligned — a block ends wherever the compressor cut it, mid-record —
/// so [`FindFastqBoundaries`](super::find_fastq_boundaries::FindFastqBoundaries)
/// re-frames the stream on whole-record (4-line) seams downstream.
///
/// Carries the same per-stream `ordinal` and `stream_idx` as the
/// [`FastqRawBlock`] it came from, so the reorder key and stream identity
/// survive the parallel decompress.
pub struct FastqDecompressedBlock {
    /// Per-stream monotonic ordinal, propagated from the source
    /// [`FastqRawBlock`]; the reorder key for `FastqDecompress`'s output.
    pub ordinal: u64,
    /// Index of the originating FASTQ stream (0-based).
    pub stream_idx: usize,
    /// Decompressed block bytes (arbitrary record alignment).
    pub data: Vec<u8>,
}

impl HeapSize for FastqDecompressedBlock {
    fn heap_size(&self) -> usize {
        self.data.capacity()
    }
}

impl Ordered for FastqDecompressedBlock {
    fn ordinal(&self) -> u64 {
        self.ordinal
    }
}

// ─────────────────────────────────────────────────────────────────────────────
// ReadFastqBlocks — Serial, per-stream raw BGZF block reader.
// ─────────────────────────────────────────────────────────────────────────────

/// `Serial` per-stream reader that emits raw (still-compressed) BGZF blocks.
///
/// This is the thin I/O half of the FASTQ decode split — the analogue of the
/// BAM path's `ReadBgzfBlocks` source. It does no decoding: it reads whole BGZF
/// blocks off the file and hands them to the `Parallel` [`FastqDecompress`]
/// step, which fans the inflate work across the pool.
///
/// One instance per stream, pinned to a distinct worker via [`Affinity`], for
/// the same reason the fused per-stream `ReadFastqInputs` readers are: a
/// `Serial` source dispatched by more than one worker races on the runtime's
/// `try_lock`-after-`Finished` source-drain path.
pub struct ReadFastqBlocks {
    reader: Option<Box<dyn Read + Send>>,
    stream_idx: usize,
    affinity: Affinity,
    next_ordinal: u64,
    pending: VecDeque<FastqRawBlock>,
    held: HeldSlot<Unpushed<FastqRawBlock>>,
    output_byte_limit: u64,
    finished: bool,
}

impl ReadFastqBlocks {
    /// Construct a per-stream BGZF block reader.
    ///
    /// `reader` is the **raw compressed** file stream (not a gzip decoder) —
    /// the chain builder opens the file directly for the BGZF split path.
    /// `stream_idx` is the global stream index (0 = R1, 1 = R2), stamped on
    /// every emitted block. `affinity` pins this reader to a distinct worker.
    #[must_use]
    pub fn new(
        reader: Box<dyn Read + Send>,
        stream_idx: usize,
        affinity: Affinity,
        output_byte_limit: u64,
    ) -> Self {
        Self {
            reader: Some(reader),
            stream_idx,
            affinity,
            next_ordinal: 0,
            pending: VecDeque::new(),
            held: HeldSlot::new(),
            output_byte_limit,
            finished: false,
        }
    }
}

impl Step for ReadFastqBlocks {
    type Input = ();
    type Outputs = OrderedBytesSingle<FastqRawBlock>;

    fn profile(&self) -> StepProfile {
        StepProfile {
            name: "ReadFastqBlocks",
            kind: StepKind::Serial,
            // Sticky like the other per-stream readers would starve the pool's
            // Parallel decompress of this worker; keep it non-sticky so the
            // affinity-pinned worker interleaves decode work between reads.
            sticky: false,
            output_queues: vec![QueueSpec::ByteBounded { limit_bytes: self.output_byte_limit }],
            branch_ordering: vec![BranchOrdering::ByItemOrdinal],
        }
    }

    fn affinity(&self) -> Affinity {
        self.affinity
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

        // 2. Drain one pending block.
        if let Some(block) = self.pending.pop_front() {
            match ctx.outputs.push(block) {
                Ok(()) => return Ok(StepOutcome::Progress),
                Err(unpushed) => {
                    self.held.put(unpushed);
                    return Ok(StepOutcome::Progress);
                }
            }
        }

        if self.finished {
            return Ok(StepOutcome::Finished);
        }

        // 3. Read up to BLOCKS_PER_BATCH raw BGZF blocks.
        let raw_blocks = {
            let reader = self
                .reader
                .as_mut()
                .expect("ReadFastqBlocks: try_run called after the source reported Finished");
            read_raw_blocks(reader.as_mut(), BLOCKS_PER_BATCH)?
        };

        if raw_blocks.is_empty() {
            self.finished = true;
            // Release the reader (and its buffer) at end of stream.
            self.reader = None;
            return Ok(StepOutcome::Finished);
        }

        for raw in raw_blocks {
            let ordinal = self.next_ordinal;
            self.next_ordinal += 1;
            self.pending.push_back(FastqRawBlock {
                ordinal,
                stream_idx: self.stream_idx,
                data: raw.data,
            });
        }

        if let Some(block) = self.pending.pop_front() {
            match ctx.outputs.push(block) {
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

// ─────────────────────────────────────────────────────────────────────────────
// FastqDecompress — Parallel BGZF block decompressor.
// ─────────────────────────────────────────────────────────────────────────────

/// `Parallel` decompressor: each pool worker holds its own
/// `libdeflater::Decompressor` (`Clone` builds a fresh one) and inflates one
/// BGZF block per `try_run`, emitting a [`FastqDecompressedBlock`]. This is the
/// step that lets FASTQ decode fan across the whole pool — the direct analogue
/// of [`BgzfDecompress`](crate::pipeline::steps::bgzf::decompress::BgzfDecompress).
///
/// A single `FastqDecompress` declaration is appended to *each* stream's chain;
/// because it is `Parallel`, all its per-worker clones draw from the one shared
/// pool, so per-stream declarations do not multiply thread usage — they only
/// keep each stream on its own edge (independent backpressure, independent
/// reorder domain).
pub struct FastqDecompress {
    decompressor: Decompressor,
    output_scratch: Vec<u8>,
    held: HeldSlot<Unpushed<FastqDecompressedBlock>>,
    output_byte_limit: u64,
    /// Whether to verify each block's CRC32. The chain threads the command's
    /// `--check-crc` / `--no-check-crc` policy here, matching the BAM path.
    verify_crc: bool,
}

impl FastqDecompress {
    /// Construct a decompressor with an explicit CRC-verification policy.
    #[must_use]
    pub fn new(output_byte_limit: u64, verify_crc: bool) -> Self {
        Self {
            decompressor: Decompressor::new(),
            output_scratch: Vec::with_capacity(DECOMPRESS_SCRATCH_CAPACITY),
            held: HeldSlot::new(),
            output_byte_limit,
            verify_crc,
        }
    }
}

impl Clone for FastqDecompress {
    fn clone(&self) -> Self {
        Self {
            decompressor: Decompressor::new(),
            output_scratch: Vec::with_capacity(DECOMPRESS_SCRATCH_CAPACITY),
            held: HeldSlot::new(),
            output_byte_limit: self.output_byte_limit,
            verify_crc: self.verify_crc,
        }
    }
}

impl Step for FastqDecompress {
    type Input = FastqRawBlock;
    type Outputs = OrderedBytesSingle<FastqDecompressedBlock>;

    fn profile(&self) -> StepProfile {
        StepProfile {
            name: "FastqDecompress",
            kind: StepKind::Parallel,
            sticky: false,
            output_queues: vec![QueueSpec::ByteBounded { limit_bytes: self.output_byte_limit }],
            branch_ordering: vec![BranchOrdering::ByItemOrdinal],
        }
    }

    fn try_run(&mut self, ctx: &mut StepCtx<'_, Self>) -> io::Result<StepOutcome> {
        if let Some(unpushed) = self.held.take() {
            match ctx.outputs.retry(unpushed) {
                Ok(()) => {}
                Err(again) => {
                    self.held.put(again);
                    // Contention (not NoProgress) keeps this worker alive to
                    // retry — NoProgress could let the framework Skip it with a
                    // held item still buffered.
                    return Ok(StepOutcome::Contention);
                }
            }
        }

        let Some(block) = ctx.input.pop() else {
            if ctx.input.is_drained() {
                return Ok(StepOutcome::Finished);
            }
            return Ok(StepOutcome::NoProgress);
        };

        // Decompress into the per-worker scratch, then swap it out with a fresh
        // fixed-capacity buffer (single mimalloc size class, hot TLS cache).
        decompress_block_slice_into_opts(
            &block.data,
            &mut self.decompressor,
            &mut self.output_scratch,
            self.verify_crc,
        )?;
        let bytes = std::mem::replace(
            &mut self.output_scratch,
            Vec::with_capacity(DECOMPRESS_SCRATCH_CAPACITY),
        );

        let out = FastqDecompressedBlock {
            ordinal: block.ordinal,
            stream_idx: block.stream_idx,
            data: bytes,
        };

        match ctx.outputs.push(out) {
            Ok(()) => Ok(StepOutcome::Progress),
            Err(unpushed) => {
                self.held.put(unpushed);
                Ok(StepOutcome::Progress)
            }
        }
    }

    fn new_worker_copy(&self) -> Self {
        self.clone()
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn read_fastq_blocks_profile_is_serial_byordinal() {
        let reader: Box<dyn Read + Send> = Box::new(io::Cursor::new(Vec::<u8>::new()));
        let s = ReadFastqBlocks::new(reader, 0, Affinity::Worker(0), 1 << 20);
        let p = s.profile();
        assert_eq!(p.name, "ReadFastqBlocks");
        assert_eq!(p.kind, StepKind::Serial);
        assert_eq!(s.affinity(), Affinity::Worker(0));
        assert_eq!(p.branch_ordering, vec![BranchOrdering::ByItemOrdinal]);
    }

    #[test]
    fn fastq_decompress_profile_is_parallel_byordinal() {
        let s = FastqDecompress::new(1 << 20, true);
        let p = s.profile();
        assert_eq!(p.name, "FastqDecompress");
        assert_eq!(p.kind, StepKind::Parallel);
        assert!(!p.sticky);
        assert_eq!(p.branch_ordering, vec![BranchOrdering::ByItemOrdinal]);
    }

    #[test]
    fn fastq_decompress_clone_constructs_fresh_decompressor() {
        let s = FastqDecompress::new(1 << 20, false);
        let _cloned = s.clone();
    }

    #[test]
    fn fastq_raw_block_heap_size_and_ordinal() {
        let block = FastqRawBlock { ordinal: 4, stream_idx: 1, data: vec![0u8; 30] };
        assert_eq!(block.ordinal(), 4);
        assert!(block.heap_size() >= 30);
    }

    #[test]
    fn fastq_decompressed_block_heap_size_and_ordinal() {
        let block = FastqDecompressedBlock { ordinal: 9, stream_idx: 0, data: vec![0u8; 12] };
        assert_eq!(block.ordinal(), 9);
        assert!(block.heap_size() >= 12);
    }
}
