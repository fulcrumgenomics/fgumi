#![deny(unsafe_code)]

//! BAM-pipeline I/O layer for the fgumi typed-step pipeline.
//!
//! Provides the source, sink, and sort typed-step building blocks plus the
//! record-batch and BGZF-block buffer types shared across the fgumi pipeline:
//!
//!   * [`source`] — BAM ingest steps ([`ReadBgzfBlocks`] and the
//!     `read_bam*` helpers) that turn a reader into decompressed blocks.
//!   * [`sink`] — BAM output steps ([`WriteBgzfFile`]).
//!   * [`sort`] — in-pipeline sort steps ([`SortAndSpill`], [`SortMerge`],
//!     [`SortBamFile`], [`SortSpillDecompress`]).
//!   * [`types`] — the record-batch / decompressed-block buffers
//!     ([`RecordBatch`], [`DecompressedBlock`], [`BgzfBlock`], …) exchanged
//!     between steps.

pub mod sink;
pub mod sort;
pub mod source;
pub mod types;

pub use fgumi_pipeline_core::HeaderHandle;
pub use sink::write_bgzf::WriteBgzfFile;
pub use sort::{SortAndSpill, SortBamFile, SortMerge, SortSpillDecompress};
pub use source::read_bam::{
    DEFAULT_BLOCKS_PER_BATCH, ReadBgzfBlocks, read_bam, read_bam_auto, read_bam_from_reader,
    read_bam_stdin,
};
pub use types::{BgzfBlock, DecompressedBlock, RecordBatch, RecordBatchBuilder};
