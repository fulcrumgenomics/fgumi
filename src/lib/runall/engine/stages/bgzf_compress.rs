//! BGZF compression stage.
//!
//! Parallel pool stage: one [`BgzfCompress`] instance per worker,
//! constructed via a factory at pool start. Each worker owns a
//! persistent [`InlineBgzfCompressor`] reused across batches — the
//! compressor's internal buffer pool recycles 64 KB buffers instead of
//! reallocating them per call, matching `unified_pipeline`'s per-worker
//! compressor model.
//!
//! ## Input
//!
//! [`SerializedBatch`] — uncompressed, length-prefixed BAM record bytes
//! with an optional secondary stream (typically rejects).
//!
//! ## Output
//!
//! [`CompressedBatch`] — concatenated complete BGZF blocks, ready to
//! append to an output file after the BAM header. Primary and
//! secondary streams are compressed independently; the `ordinal` field
//! is passed through unchanged from input to output.
//!
//! ## Ordering guarantees
//!
//! None imposed by this stage. `ordinal` is preserved verbatim and is
//! consumed downstream by
//! [`crate::runall::engine::sink::bam_file_write::BamFileWrite`]'s
//! [`crate::runall::engine::reorder::ReorderBuffer`] to restore the order
//! established upstream.
//!
//! ## Memory model
//!
//! 1:1 with input: one output batch per input batch. The output buffer
//! is freshly allocated at `total_len` capacity per call; the
//! compressor's internal buffer pool recycles its own 64 KB scratch.
//!
//! ## Determinism
//!
//! Byte-identical per batch: given identical input bytes and
//! compression level, `libdeflate` produces the same BGZF output.
//! Across a multi-worker run, per-batch compressed bytes are
//! deterministic; only the ordering of batches in the output stream
//! depends on the reorder buffer at the sink (which is deterministic).

use anyhow::{Context, Result};
use fgumi_bgzf::InlineBgzfCompressor;

use crate::runall::engine::output_types::{
    CompressedBatch, CompressedBytes, RawBytes, SerializedBatch,
};
use crate::runall::engine::stage::{Parallelism, Stage};

/// Parallel stage: compresses [`SerializedBatch`] BAM-record bytes into
/// concatenated BGZF blocks, preserving ordinal and record counts.
///
/// Holds a persistent [`InlineBgzfCompressor`] that is reused across all
/// batches processed by this worker instance. Not `Clone` because the
/// compressor holds internal state (buffer pool, zlib deflate context).
pub struct BgzfCompress {
    compressor: InlineBgzfCompressor,
}

impl BgzfCompress {
    /// Construct a new stage with the given BGZF compression level (1..=12).
    #[must_use]
    pub fn new(compression_level: u32) -> Self {
        Self { compressor: InlineBgzfCompressor::new(compression_level) }
    }
}

impl Stage for BgzfCompress {
    type Input = SerializedBatch;
    type Output = CompressedBatch;

    #[tracing::instrument(name = "bgzf_compress", skip_all)]
    fn process(&mut self, input: Self::Input, out: &mut dyn FnMut(Self::Output)) -> Result<()> {
        let primary = compress_stream(&mut self.compressor, &input.primary)
            .context("BgzfCompress: primary stream")?;
        let secondary = match input.secondary.as_ref() {
            Some(raw) => Some(
                compress_stream(&mut self.compressor, raw)
                    .context("BgzfCompress: secondary stream")?,
            ),
            None => None,
        };
        out(CompressedBatch { primary, secondary, ordinal: input.ordinal });
        Ok(())
    }

    fn parallelism(&self) -> Parallelism {
        Parallelism::Parallel
    }

    fn output_memory_estimate(&self, out: &Self::Output) -> usize {
        out.primary.data.capacity() + out.secondary.as_ref().map_or(0, |s| s.data.capacity())
    }

    fn name(&self) -> &'static str {
        "BgzfCompress"
    }
}

/// Compress one [`RawBytes`] stream into a single contiguous `Vec<u8>` of
/// BGZF blocks, preserving the record count. Reuses the given compressor's
/// internal state — callers must drain `take_blocks()` between streams so
/// primary and secondary outputs don't mix.
fn compress_stream(
    compressor: &mut InlineBgzfCompressor,
    raw: &RawBytes,
) -> Result<CompressedBytes> {
    compressor.write_all(&raw.data).context("InlineBgzfCompressor::write_all failed")?;
    compressor.flush().context("InlineBgzfCompressor::flush failed")?;

    let blocks = compressor.take_blocks();
    let total_len: usize = blocks.iter().map(|b| b.data.len()).sum();
    let mut data = Vec::with_capacity(total_len);
    for block in blocks {
        data.extend_from_slice(&block.data);
    }

    Ok(CompressedBytes { data, record_count: raw.record_count })
}

#[cfg(test)]
mod tests {
    use super::*;
    use fgumi_bgzf::{decompress_block, read_raw_blocks};
    use libdeflater::Decompressor;
    use std::io::Cursor;

    /// Decompress a contiguous stream of BGZF blocks into a single Vec<u8>.
    fn decompress_bgzf_stream(blocks_data: &[u8]) -> Vec<u8> {
        let mut out = Vec::new();
        let mut reader = Cursor::new(blocks_data);
        let mut decompressor = Decompressor::new();
        while let Some(block) =
            read_raw_blocks(&mut reader, 1).expect("reading BGZF blocks").into_iter().next()
        {
            let bytes = decompress_block(&block, &mut decompressor).expect("decompress block");
            out.extend_from_slice(&bytes);
        }
        out
    }

    fn raw(bytes: Vec<u8>, count: u64) -> RawBytes {
        RawBytes { data: bytes, record_count: count }
    }

    fn process_emit(stage: &mut BgzfCompress, input: SerializedBatch) -> CompressedBatch {
        let mut captured = None;
        stage
            .process(input, &mut |v| {
                captured = Some(v);
            })
            .expect("compression must succeed");
        captured.expect("stage must emit")
    }

    #[test]
    fn test_bgzf_compress_round_trip_primary_only() {
        let mut stage = BgzfCompress::new(1);
        let payload: Vec<u8> = (0..5000u32).flat_map(u32::to_le_bytes).collect();
        let input =
            SerializedBatch { primary: raw(payload.clone(), 42), secondary: None, ordinal: 3 };
        let out = process_emit(&mut stage, input);

        assert_eq!(out.ordinal, 3);
        assert!(out.secondary.is_none());
        assert_eq!(out.primary.record_count, 42);
        let decoded = decompress_bgzf_stream(&out.primary.data);
        assert_eq!(decoded, payload, "round-trip byte mismatch");
    }

    #[test]
    fn test_bgzf_compress_empty_input_produces_no_blocks() {
        let mut stage = BgzfCompress::new(1);
        let input = SerializedBatch { primary: raw(vec![], 0), secondary: None, ordinal: 0 };
        let out = process_emit(&mut stage, input);
        assert_eq!(out.primary.record_count, 0);
        assert!(out.primary.data.is_empty(), "empty input must produce zero-byte output");
    }

    #[test]
    fn test_bgzf_compress_secondary_present_round_trips_both() {
        let mut stage = BgzfCompress::new(1);
        let primary_payload = b"PRIMARY-STREAM-DATA".repeat(100);
        let secondary_payload = b"SECONDARY-STREAM".repeat(50);
        let input = SerializedBatch {
            primary: raw(primary_payload.clone(), 10),
            secondary: Some(raw(secondary_payload.clone(), 5)),
            ordinal: 11,
        };
        let out = process_emit(&mut stage, input);
        assert_eq!(out.ordinal, 11);
        assert_eq!(out.primary.record_count, 10);
        assert_eq!(out.secondary.as_ref().unwrap().record_count, 5);
        assert_eq!(decompress_bgzf_stream(&out.primary.data), primary_payload);
        assert_eq!(decompress_bgzf_stream(&out.secondary.unwrap().data), secondary_payload);
    }

    #[test]
    fn test_bgzf_compress_preserves_ordinal() {
        let mut stage = BgzfCompress::new(1);
        for ord in [0u64, 1, 42, u64::MAX - 1] {
            let out = process_emit(
                &mut stage,
                SerializedBatch { primary: raw(vec![1, 2, 3], 1), secondary: None, ordinal: ord },
            );
            assert_eq!(out.ordinal, ord);
        }
    }

    #[test]
    fn test_bgzf_compress_declares_parallel() {
        let stage = BgzfCompress::new(6);
        assert_eq!(stage.parallelism(), Parallelism::Parallel);
        assert_eq!(stage.name(), "BgzfCompress");
    }

    /// The persistent compressor must round-trip consecutive batches correctly:
    /// `take_blocks()` drains previous output so the next batch starts clean.
    #[test]
    fn test_bgzf_compress_reuses_compressor_across_batches() {
        let mut stage = BgzfCompress::new(1);

        let payload_a = b"PAYLOAD-A-0123456789".repeat(100);
        let out_a = process_emit(
            &mut stage,
            SerializedBatch { primary: raw(payload_a.clone(), 10), secondary: None, ordinal: 0 },
        );
        assert_eq!(decompress_bgzf_stream(&out_a.primary.data), payload_a);

        let payload_b = b"PAYLOAD-B-abcdefghij".repeat(200);
        let out_b = process_emit(
            &mut stage,
            SerializedBatch { primary: raw(payload_b.clone(), 20), secondary: None, ordinal: 1 },
        );
        assert_eq!(decompress_bgzf_stream(&out_b.primary.data), payload_b);

        // Batch B's output must not contain any of batch A's bytes.
        assert_ne!(out_a.primary.data, out_b.primary.data);
    }
}
