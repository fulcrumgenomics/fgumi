//! BGZF compression of one output block, wherever it runs.
//!
//! [`OutputCompressor`] turns a [`DecompressedBlock`] of framed BAM records
//! into one logical [`BgzfBlock`]: the concatenation of the physical 64 KiB
//! BGZF blocks the compressor produced, and, when indexing, the
//! [`BamIndexManifest`] an inline BAI indexer needs (built from the
//! uncompressed bytes). The pool's `BgzfCompress` step runs one per worker,
//! and a range-mode sort merge runs one of its own to compress the blocks its
//! output edge refuses ([`OutputBlock::Compressed`]); libdeflate at a fixed
//! level is deterministic on its input, so both produce the same bytes and
//! manifests for a block.

use std::io;
use std::sync::atomic::{AtomicU64, Ordering};

use fgumi_bgzf::writer::InlineBgzfCompressor;
use fgumi_pipeline_core::{HeapSize, Ordered};

use crate::types::{BamIndexManifest, BgzfBlock, DecompressedBlock, RecordIndexEntry};

/// Per-compressor output scratch capacity. The compressed output of one
/// block is the concatenation of physical 64 KiB BGZF blocks; a batch
/// typically produces 1-4 of them. Sized to the typical case so the same
/// mimalloc size class is hit every call.
const COMPRESS_SCRATCH_CAPACITY: usize = 256 * 1024;

/// Compresses output blocks into [`BgzfBlock`]s at one level, with or
/// without the inline BAI manifest.
pub struct OutputCompressor {
    compressor: InlineBgzfCompressor,
    level: u32,
    /// Output scratch for the multi-block case; `mem::replace`d on each use
    /// so the replacement is always the same size class.
    scratch: Vec<u8>,
    /// Build a [`BamIndexManifest`] for each block.
    index_bam: bool,
}

impl OutputCompressor {
    /// A compressor at `level`, building manifests when `index_bam`.
    #[must_use]
    pub fn new(level: u32, index_bam: bool) -> Self {
        Self {
            compressor: InlineBgzfCompressor::new(level),
            level,
            scratch: Vec::with_capacity(COMPRESS_SCRATCH_CAPACITY),
            index_bam,
        }
    }

    /// Compress `block` into one [`BgzfBlock`] with the same `batch_serial`.
    /// The manifest (when indexing) is built from the uncompressed framed
    /// records, never from the compressed output.
    ///
    /// # Errors
    ///
    /// The compressor's own error.
    pub fn compress(&mut self, block: DecompressedBlock) -> io::Result<BgzfBlock> {
        let DecompressedBlock { batch_serial, bytes: input } = block;
        let uncompressed_size = u32::try_from(input.len()).unwrap_or(u32::MAX);
        self.compressor.write_all(&input)?;
        self.compressor.flush()?;
        let (bytes, index) = if self.index_bam {
            let records = build_record_entries(&input);
            let (bytes, phys_comp_len) = self.assemble_output_bytes_with_lens();
            debug_assert_eq!(
                phys_comp_len.iter().map(|&n| u64::from(n)).sum::<u64>(),
                bytes.len() as u64,
                "physical block lengths must sum to the total compressed byte count"
            );
            debug_assert_eq!(
                phys_comp_len.len(),
                (uncompressed_size as usize).div_ceil(fgumi_bgzf::BGZF_MAX_BLOCK_SIZE),
                "one physical block per BGZF_MAX_BLOCK_SIZE-sized chunk of uncompressed input"
            );
            let index =
                (!bytes.is_empty()).then(|| Box::new(BamIndexManifest { phys_comp_len, records }));
            (bytes, index)
        } else {
            (self.assemble_output_bytes(), None)
        };
        Ok(BgzfBlock { batch_serial, bytes: bytes.into(), uncompressed_size, index })
    }

    /// Drain the compressor's physical BGZF blocks into one byte buffer.
    ///
    /// INVARIANT: each physical block's `serial` is `InlineBgzfCompressor`'s
    /// per-compressor counter, unrelated to the pipeline's `batch_serial`; it
    /// is dropped here and must never enter the framework's ordinal stream.
    ///
    /// One physical block (the common case) is moved out as is; several are
    /// concatenated into the scratch (then `mem::replace`d out) and their
    /// buffers handed back to the compressor's pool. No physical block (an
    /// empty input) yields an unallocated `Vec`, not the scratch.
    fn assemble_output_bytes(&mut self) -> Vec<u8> {
        let children = self.compressor.take_blocks();
        self.concat(children)
    }

    /// [`Self::assemble_output_bytes`] plus each physical block's compressed
    /// length, in emission order (a manifest's `phys_comp_len`).
    fn assemble_output_bytes_with_lens(&mut self) -> (Vec<u8>, Vec<u32>) {
        let children = self.compressor.take_blocks();
        let lens = children
            .iter()
            .map(|c| u32::try_from(c.data.len()).expect("block len fits u32"))
            .collect();
        (self.concat(children), lens)
    }

    fn concat(&mut self, mut children: Vec<fgumi_bgzf::writer::CompressedBlock>) -> Vec<u8> {
        match children.len() {
            0 => Vec::new(),
            1 => children.pop().expect("len == 1").data,
            _ => {
                for child in children {
                    self.scratch.extend_from_slice(&child.data);
                    self.compressor.recycle_buffer(child.data);
                }
                std::mem::replace(&mut self.scratch, Vec::with_capacity(COMPRESS_SCRATCH_CAPACITY))
            }
        }
    }
}

impl Clone for OutputCompressor {
    /// A fresh compressor with the same settings (no shared state).
    fn clone(&self) -> Self {
        Self::new(self.level, self.index_bam)
    }
}

/// Walk a batch of `[u32 le len][body]*` frames into per-record index
/// entries. `uoffset` is the frame start (offset of the length prefix from
/// the start of the batch's uncompressed bytes); `len` is `4 +` the body
/// length.
#[must_use]
pub(crate) fn build_record_entries(batch: &[u8]) -> Vec<RecordIndexEntry> {
    // Records in a block are similar in size, so the first frame's length
    // gives a cheap estimate of the record count that avoids most of the
    // doubling reallocations; a short estimate just grows as before.
    let estimate = batch.get(..4).map_or(0, |prefix| {
        let first = 4 + u32::from_le_bytes(prefix.try_into().expect("4 bytes")) as usize;
        batch.len() / first
    });
    let mut out = Vec::with_capacity(estimate);
    let mut u = 0usize;
    while u + 4 <= batch.len() {
        let block_size = u32::from_le_bytes(batch[u..u + 4].try_into().expect("4 bytes")) as usize;
        let body = &batch[u + 4..u + 4 + block_size];
        out.push(RecordIndexEntry {
            uoffset: u32::try_from(u).expect("batch offset fits u32"),
            len: u32::try_from(4 + block_size).expect("record len fits u32"),
            ctx: fgumi_bam_io::extract_alignment_context(body),
        });
        u += 4 + block_size;
    }
    debug_assert_eq!(u, batch.len(), "frames must exactly cover the batch");
    out
}

/// The depth of a sort merge's output queue: blocks the merge pushed onto
/// the edge to `BgzfCompress` that no compressor has popped yet. The merge
/// counts its accepted pushes, the compressor its pops; the difference is
/// compression's backlog whatever the edge's capacity (the queue budget may
/// grow or shrink it while blocks wait). Range mode gates the pool's
/// merge-side claims on it.
#[derive(Debug, Default)]
pub struct OutputDepth {
    pushed: AtomicU64,
    popped: AtomicU64,
}

impl OutputDepth {
    /// A counter pair at zero.
    #[must_use]
    pub fn new() -> Self {
        Self::default()
    }

    /// Count `n` blocks the producer pushed onto the edge.
    pub fn note_pushed(&self, n: u64) {
        self.pushed.fetch_add(n, Ordering::AcqRel);
    }

    /// Count one block a compressor popped from the edge.
    pub fn note_popped(&self) {
        self.popped.fetch_add(1, Ordering::AcqRel);
    }

    /// Blocks on the edge not yet popped, saturating at 0.
    ///
    /// The count can transiently read low: the producer counts its pushes
    /// after it pushes them (`RangeModeState::flush_output`), so a compressor
    /// can pop, and count, a block before its push is counted. A low read
    /// only opens the workers' claim gate early for a moment, which is
    /// harmless. It also reads 0 when pops ran ahead of the pushes counted
    /// for good (blocks a serial merge pushed without counting).
    #[must_use]
    pub fn depth(&self) -> u64 {
        let popped = self.popped.load(Ordering::Acquire);
        self.pushed.load(Ordering::Acquire).saturating_sub(popped)
    }

    /// Blocks popped so far.
    #[must_use]
    pub fn popped(&self) -> u64 {
        self.popped.load(Ordering::Acquire)
    }

    /// Pushes counted so far (tests): unlike [`Self::depth`], not saturated,
    /// so a push the producer never counted shows.
    #[cfg(test)]
    #[must_use]
    pub(crate) fn pushed_for_test(&self) -> u64 {
        self.pushed.load(Ordering::Acquire)
    }
}

/// An output block on its way to the BGZF writer: framed records still to
/// compress, or already compressed by the producer (a range-mode sort
/// merge's thread, while its output edge refused it).
#[derive(Debug)]
pub enum OutputBlock {
    /// Framed BAM records, for `BgzfCompress` to compress.
    Raw(DecompressedBlock),
    /// Compressed by the producer with an [`OutputCompressor`] of the
    /// writer's settings; `BgzfCompress` passes it through.
    Compressed(BgzfBlock),
}

impl OutputBlock {
    /// The block's serial (either variant's).
    #[must_use]
    pub fn batch_serial(&self) -> u64 {
        match self {
            Self::Raw(b) => b.batch_serial,
            Self::Compressed(b) => b.batch_serial,
        }
    }

    /// Compress a `Raw` block in place with `c` (same serial); `false` when
    /// it was already compressed.
    ///
    /// # Errors
    ///
    /// The compressor's own error.
    pub fn compress_in_place(&mut self, c: &mut OutputCompressor) -> io::Result<bool> {
        let Self::Raw(raw) = self else { return Ok(false) };
        let raw = std::mem::replace(
            raw,
            DecompressedBlock { batch_serial: raw.batch_serial, bytes: Vec::new() },
        );
        *self = Self::Compressed(c.compress(raw)?);
        Ok(true)
    }
}

impl HeapSize for OutputBlock {
    fn heap_size(&self) -> usize {
        match self {
            Self::Raw(b) => b.heap_size(),
            Self::Compressed(b) => b.heap_size(),
        }
    }
}

impl Ordered for OutputBlock {
    fn ordinal(&self) -> u64 {
        self.batch_serial()
    }
}

/// The input of a BGZF compression step: a block to compress, or one its
/// producer already compressed.
pub trait CompressInput: Send + HeapSize + Ordered + 'static {
    /// The block as an [`OutputBlock`].
    fn into_output_block(self) -> OutputBlock;
}

impl CompressInput for DecompressedBlock {
    fn into_output_block(self) -> OutputBlock {
        OutputBlock::Raw(self)
    }
}

impl CompressInput for OutputBlock {
    fn into_output_block(self) -> OutputBlock {
        self
    }
}

#[cfg(test)]
mod tests {
    use fgumi_sort::range_merge_test_support::frame_bodies as framed;
    use rstest::rstest;

    use super::*;

    /// The plain concatenation of a fresh compressor's physical blocks.
    fn reference_concat(level: u32, input: &[u8]) -> Vec<u8> {
        let mut c = InlineBgzfCompressor::new(level);
        c.write_all(input).unwrap();
        c.flush().unwrap();
        c.take_blocks().into_iter().flat_map(|b| b.data).collect()
    }

    /// Bytes of 200 KiB: more than one physical block.
    fn big() -> Vec<u8> {
        (0..(200 * 1024u32)).map(|i| u8::try_from(i % 251).unwrap()).collect()
    }

    /// A mapped record body at `pos` with a `5M` CIGAR.
    fn record_body(pos: i32) -> Vec<u8> {
        fgumi_raw_bam::SamBuilder::new()
            .ref_id(0)
            .pos(pos)
            .flags(0)
            .cigar_ops(&[5 << 4])
            .build()
            .to_vec()
    }

    /// An `OutputCompressor` emits exactly a fresh compressor's
    /// concatenation, at every level, for one physical block (moved out) and
    /// for several (concatenated through the scratch), and keeps doing so
    /// after its pool is primed by recycled buffers.
    #[rstest]
    fn compress_equals_a_fresh_compressor(
        #[values(0, 1, 6)] level: u32,
        #[values(false, true)] multi: bool,
        #[values(false, true)] index_bam: bool,
    ) {
        // Framed records (the manifest walks frames): 4 fit one physical
        // block, 4000 span several.
        let records = if multi { 4000 } else { 4 };
        let input = framed(&(0..records).map(|i| record_body(i * 10)).collect::<Vec<_>>());
        assert_eq!(input.len() > fgumi_bgzf::BGZF_MAX_BLOCK_SIZE, multi);
        let mut c = OutputCompressor::new(level, index_bam);
        for serial in 0..2 {
            let block = c
                .compress(DecompressedBlock { batch_serial: serial, bytes: input.clone() })
                .unwrap();
            assert_eq!(block.batch_serial, serial);
            assert_eq!(block.uncompressed_size as usize, input.len());
            assert_eq!(&block.bytes[..], &reference_concat(level, &input)[..], "call {serial}");
            assert_eq!(&block.bytes[..2], &[0x1f, 0x8b], "BGZF magic");
        }
    }

    #[test]
    fn empty_block_returns_unallocated_bytes() {
        let mut c = OutputCompressor::new(1, false);
        // Prime the scratch with a multi-block batch first.
        c.compress(DecompressedBlock { batch_serial: 0, bytes: big() }).unwrap();
        let empty = c.compress(DecompressedBlock { batch_serial: 1, bytes: Vec::new() }).unwrap();
        assert!(empty.bytes.is_empty());
        assert_eq!(empty.heap_size(), 0, "an empty output carries no scratch allocation");
        assert!(empty.index.is_none());
    }

    #[test]
    fn build_record_entries_matches_frames() {
        let bodies = vec![record_body(0), record_body(100)];
        let entries = build_record_entries(&framed(&bodies));
        assert_eq!(entries.len(), 2);
        assert_eq!(entries[0].uoffset, 0);
        assert_eq!(entries[1].uoffset as usize, 4 + bodies[0].len());
        for (e, body) in entries.iter().zip(&bodies) {
            assert_eq!(e.len as usize, 4 + body.len());
            assert!(e.ctx.is_some());
        }
    }

    /// Two compressors build the same manifest for the same block; its
    /// physical lengths sum to the bytes and number one per 64 KiB chunk.
    #[rstest]
    #[case::one_block(4)]
    #[case::several_blocks(4000)]
    fn manifests_agree_across_compressors(#[case] records: i32) {
        let bodies: Vec<Vec<u8>> = (0..records).map(|i| record_body(i * 10)).collect();
        let bytes = framed(&bodies);
        let block = || DecompressedBlock { batch_serial: 3, bytes: bytes.clone() };
        let a = OutputCompressor::new(1, true).compress(block()).unwrap();
        let mut primed = OutputCompressor::new(1, true);
        let primer = framed(&(0..4000).map(|i| record_body(i * 7)).collect::<Vec<_>>());
        primed.compress(DecompressedBlock { batch_serial: 0, bytes: primer }).unwrap();
        let b = primed.compress(block()).unwrap();
        assert_eq!(&a.bytes[..], &b.bytes[..]);
        let (ma, mb) = (a.index.expect("a manifest"), b.index.expect("a manifest"));
        assert_eq!(ma.phys_comp_len, mb.phys_comp_len);
        assert_eq!(
            ma.phys_comp_len.iter().map(|&n| n as usize).sum::<usize>(),
            a.bytes.len(),
            "physical lengths sum to the bytes"
        );
        assert_eq!(ma.phys_comp_len.len(), bytes.len().div_ceil(fgumi_bgzf::BGZF_MAX_BLOCK_SIZE));
        let key = |m: &BamIndexManifest| -> Vec<(u32, u32, String)> {
            m.records.iter().map(|r| (r.uoffset, r.len, format!("{:?}", r.ctx))).collect()
        };
        assert_eq!(key(&ma), key(&mb));
    }

    /// Compressing a `Raw` block in place keeps its serial and size, equals
    /// the compressor's output, and is done once.
    #[test]
    fn compress_in_place_once() {
        let bodies: Vec<Vec<u8>> = (0..200).map(|i| record_body(i * 10)).collect();
        let raw = DecompressedBlock { batch_serial: 7, bytes: framed(&bodies) };
        let len = raw.bytes.len();
        let want = OutputCompressor::new(6, true)
            .compress(DecompressedBlock { batch_serial: 7, bytes: raw.bytes.clone() })
            .unwrap();
        let mut block = OutputBlock::Raw(raw);
        let mut c = OutputCompressor::new(6, true);
        assert!(block.compress_in_place(&mut c).unwrap());
        assert!(!block.compress_in_place(&mut c).unwrap(), "already compressed");
        assert_eq!(block.batch_serial(), 7);
        assert_eq!(block.ordinal(), 7);
        let OutputBlock::Compressed(got) = block else { panic!("compressed") };
        assert_eq!(got.uncompressed_size as usize, len);
        assert_eq!(&got.bytes[..], &want.bytes[..]);
        assert_eq!(got.index.map(|m| m.phys_comp_len), want.index.map(|m| m.phys_comp_len));
    }

    /// A clone is a fresh compressor with the original's level and indexing:
    /// it emits a fresh compressor's bytes at that level and the manifest
    /// exactly when the original builds one.
    #[rstest]
    #[case::level_1_indexed(1, true)]
    #[case::level_6_unindexed(6, false)]
    fn clone_keeps_settings(#[case] level: u32, #[case] index_bam: bool) {
        let input = framed(&(0..200).map(|i| record_body(i * 10)).collect::<Vec<_>>());
        let original = OutputCompressor::new(level, index_bam);
        let got = original
            .clone()
            .compress(DecompressedBlock { batch_serial: 0, bytes: input.clone() })
            .unwrap();
        assert_eq!(&got.bytes[..], &reference_concat(level, &input)[..]);
        assert_eq!(got.index.is_some(), index_bam);
    }

    /// The depth is pushes minus pops, never below zero.
    /// Columns: pushes, pops, expected depth.
    #[rstest]
    #[case::empty(0, 0, 0)]
    #[case::waiting(5, 2, 3)]
    #[case::drained(4, 4, 0)]
    #[case::uncounted_pushes(1, 3, 0)]
    fn output_depth_counts_waiting_blocks(
        #[case] pushes: u64,
        #[case] pops: u64,
        #[case] depth: u64,
    ) {
        let d = OutputDepth::new();
        d.note_pushed(pushes);
        for _ in 0..pops {
            d.note_popped();
        }
        assert_eq!(d.depth(), depth);
        assert_eq!(d.popped(), pops);
    }

    #[test]
    fn inputs_become_output_blocks() {
        let raw = DecompressedBlock { batch_serial: 2, bytes: vec![1, 2] };
        assert!(matches!(raw.into_output_block(), OutputBlock::Raw(b) if b.batch_serial == 2));
        let done = OutputCompressor::new(1, false)
            .compress(DecompressedBlock { batch_serial: 5, bytes: vec![3] })
            .unwrap();
        let block = OutputBlock::Compressed(done).into_output_block();
        assert!(matches!(block, OutputBlock::Compressed(b) if b.batch_serial == 5));
    }
}
