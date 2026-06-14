//! End-to-end no-op BAM round-trip helper. Wires the full block-level
//! BAM chain (`ReadBgzfBlocks → BgzfDecompress → FindBamBoundaries →
//! DecodeRecords → GroupBam → SerializeBamRecords → BgzfCompress →
//! WriteBgzfFile`) and runs it against an input BAM, producing an output
//! BAM. Convenience for the `fgumi compare bam-roundtrip <bam>`
//! gate (Phase 3 bench-prep) and for any caller that wants to drive the
//! framework end-to-end on a real BAM.

use std::io;
use std::path::Path;

use fgumi_bam_io::PipelineReaderOpts;

use crate::pipeline::core::builder::{Pipeline, PipelineConfig};
use crate::pipeline::steps::bgzf::compress::BgzfCompress;
use crate::pipeline::steps::bgzf::decompress::BgzfDecompress;
use crate::pipeline::steps::boundaries::bam::FindBamBoundaries;
use crate::pipeline::steps::group::bam::GroupBam;
use crate::pipeline::steps::parse::decode::DecodeRecords;
use crate::pipeline::steps::serialize::SerializeBamRecords;
use crate::pipeline::steps::sink::write_bgzf::WriteBgzfFile;
use crate::pipeline::steps::source::read_bam::read_bam;
use crate::pipeline::steps::tuning::BamPipelineTuning;
use fgumi_bam_io::GroupKeyConfig;

/// Configuration for [`run_bam_roundtrip`]. Backed by
/// [`BamPipelineTuning`] which mirrors legacy `auto_tuned` defaults
/// (`blocks_per_batch` 16-64 by thread count, `template_batch_size=500`,
/// etc.).
#[derive(Debug, Clone, Copy)]
pub struct RoundtripConfig {
    pub tuning: BamPipelineTuning,
}

impl RoundtripConfig {
    /// Auto-tune for `threads` workers. Equivalent to
    /// `RoundtripConfig { tuning: BamPipelineTuning::auto_tuned(threads) }`.
    #[must_use]
    pub fn auto_tuned(threads: usize) -> Self {
        Self { tuning: BamPipelineTuning::auto_tuned(threads) }
    }

    /// Override the BGZF compression level on the underlying tuning.
    #[must_use]
    pub fn with_compression_level(mut self, level: u32) -> Self {
        self.tuning = self.tuning.with_compression_level(level);
        self
    }
}

impl Default for RoundtripConfig {
    fn default() -> Self {
        Self::auto_tuned(4)
    }
}

/// Run the full block-level BAM chain on `input` and write the output to
/// `output`. The output BAM is record-equivalent to the input modulo
/// header `@PG` (Phase 3 has no `@PG` mutation).
///
/// # Errors
///
/// Returns I/O errors from file open / read / write, BAM parsing errors,
/// or `BuildError` propagated as `io::Error::other` if the pipeline graph
/// is malformed.
pub fn run_bam_roundtrip(input: &Path, output: &Path, cfg: RoundtripConfig) -> io::Result<()> {
    let t = cfg.tuning;
    let opts = PipelineReaderOpts::default();
    let (read_step, header) = read_bam(input, opts, t.blocks_per_batch, t.per_step_byte_limit)?;

    let builder = Pipeline::builder();
    builder
        .chain(read_step)
        .chain(BgzfDecompress::new(t.per_step_byte_limit))
        .chain(FindBamBoundaries::new(t.per_step_byte_limit))
        // DecodeRecords does parse + per-record GroupKey computation
        // (CIGAR walk + tag scan + library/cell-barcode hashing) in one
        // parallel pass. Matches the legacy pipeline's combined Decode
        // step; drops the intermediate Parse → Decode queue.
        .chain(DecodeRecords::new(GroupKeyConfig::default(), t.per_step_byte_limit))
        .chain(GroupBam::new(t.template_batch_size, t.per_step_byte_limit))
        .chain(SerializeBamRecords::new(t.per_step_byte_limit))
        .chain(BgzfCompress::new(t.compression_level, t.per_step_byte_limit))
        .chain(
            WriteBgzfFile::new(output, &header, t.compression_level)
                .map_err(|e| io::Error::other(format!("WriteBgzfFile::new: {e}")))?,
        )
        .into_sink_marker();
    let pipeline =
        builder.build().map_err(|e| io::Error::other(format!("Pipeline::build: {e:?}")))?;
    pipeline
        .run(PipelineConfig { threads: t.threads, ..Default::default() })
        .map_err(|e| io::Error::other(format!("Pipeline::run: {e:?}")))?;
    Ok(())
}
