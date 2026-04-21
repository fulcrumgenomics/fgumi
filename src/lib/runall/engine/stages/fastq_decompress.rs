//! Parallel BGZF decompression stage.
//!
//! Parallel pool stage. Each instance holds its own
//! `libdeflater::Decompressor` (expensive to allocate, cheap to reuse).
//! Cloned per worker — `Clone` creates a fresh decompressor.
//!
//! ## Input
//!
//! [`PerStreamChunk`] with raw BGZF-compressed bytes
//! (`offsets == None`), tagged by `stream_idx` and `batch_num`.
//!
//! ## Output
//!
//! [`PerStreamChunk`] with decompressed bytes (`offsets == None`;
//! record boundaries are found in the downstream parse step). The
//! `stream_idx`, `batch_num`, and `is_last` fields are preserved.
//!
//! ## Ordering guarantees
//!
//! Per-stream order is preserved via the passed-through `batch_num`.
//! No cross-stream ordering is imposed; stream coordination happens in
//! the downstream parse/merge stages.
//!
//! ## Memory model
//!
//! 1:1 with input. Decompressed output is freshly allocated at
//! `raw.len() * 4` initial capacity (BGZF's worst-case inflation is
//! bounded by `BGZF_MAX_BLOCK_SIZE`).
//!
//! ## Determinism
//!
//! Byte-identical: zlib decompression is a deterministic function of
//! the compressed bytes.

use anyhow::Result;

use crate::runall::engine::fastq_types::PerStreamChunk;
use crate::runall::engine::stage::{Parallelism, Stage};

pub struct BgzfDecompress {
    decompressor: libdeflater::Decompressor,
}

impl BgzfDecompress {
    #[must_use]
    pub fn new() -> Self {
        Self { decompressor: libdeflater::Decompressor::new() }
    }
}

impl Default for BgzfDecompress {
    fn default() -> Self {
        Self::new()
    }
}

impl Clone for BgzfDecompress {
    fn clone(&self) -> Self {
        Self::new()
    }
}

impl Stage for BgzfDecompress {
    type Input = PerStreamChunk;
    type Output = PerStreamChunk;

    #[tracing::instrument(name = "bgzf_decompress", skip_all)]
    fn process(&mut self, input: Self::Input, out: &mut dyn FnMut(Self::Output)) -> Result<()> {
        if input.is_last && input.data.is_empty() {
            out(input);
            return Ok(());
        }
        let decompressed = decompress_bgzf_chunk(&input.data, &mut self.decompressor)?;
        out(PerStreamChunk {
            stream_idx: input.stream_idx,
            batch_num: input.batch_num,
            data: decompressed,
            offsets: None,
            is_last: input.is_last,
        });
        Ok(())
    }

    fn parallelism(&self) -> Parallelism {
        Parallelism::Parallel
    }
    fn output_memory_estimate(&self, output: &Self::Output) -> usize {
        output.data.len()
    }
    fn name(&self) -> &'static str {
        "BgzfDecompress"
    }
}

/// Decompress a buffer containing one or more concatenated BGZF blocks.
///
/// Mirrors `src/lib/unified_pipeline/fastq.rs::decompress_bgzf_chunk`.
fn decompress_bgzf_chunk(
    raw: &[u8],
    decompressor: &mut libdeflater::Decompressor,
) -> Result<Vec<u8>> {
    let mut out = Vec::with_capacity(raw.len() * 4);
    let mut pos = 0;
    while pos < raw.len() {
        // Each BGZF block has: 18-byte fixed header + CDATA + 8-byte footer.
        // The header's BSIZE (bytes 16-17, little-endian) + 1 is the total block size.
        // The footer's ISIZE (bytes block_end-4..block_end, little-endian) is the uncompressed size.
        anyhow::ensure!(pos + 18 <= raw.len(), "truncated BGZF at {pos} (no room for header)");
        let bsize = u16::from_le_bytes([raw[pos + 16], raw[pos + 17]]);
        let block_total = bsize as usize + 1;
        let block_end = pos + block_total;
        anyhow::ensure!(
            block_end <= raw.len(),
            "truncated BGZF block at {pos}: need {block_total} bytes, have {}",
            raw.len() - pos
        );
        fgumi_bgzf::decompress_block_slice_into(&raw[pos..block_end], decompressor, &mut out)?;
        pos = block_end;
    }
    Ok(out)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_passes_through_eof_sentinel() {
        let mut stage = BgzfDecompress::new();
        let input = PerStreamChunk {
            stream_idx: 0,
            batch_num: 42,
            data: Vec::new(),
            offsets: None,
            is_last: true,
        };
        let mut captured = None;
        stage
            .process(input, &mut |v| {
                captured = Some(v);
            })
            .unwrap();
        let output = captured.expect("stage must emit");
        assert!(output.is_last);
        assert!(output.data.is_empty());
        assert_eq!(output.batch_num, 42);
    }

    #[test]
    fn test_decompresses_real_bgzf_block() {
        // Create a BGZF block using the inline compressor
        use fgumi_bgzf::InlineBgzfCompressor;

        let original_data = b"hello world";
        let mut compressor = InlineBgzfCompressor::new(6);
        compressor.write_all(original_data).unwrap();
        compressor.flush().unwrap();
        let blocks = compressor.take_blocks();
        assert_eq!(blocks.len(), 1);

        let raw = blocks[0].data.clone();

        let mut stage = BgzfDecompress::new();
        let input = PerStreamChunk {
            stream_idx: 0,
            batch_num: 0,
            data: raw,
            offsets: None,
            is_last: false,
        };
        let mut captured = None;
        stage
            .process(input, &mut |v| {
                captured = Some(v);
            })
            .unwrap();
        let output = captured.expect("stage must emit");
        assert_eq!(&output.data, original_data);
    }
}
