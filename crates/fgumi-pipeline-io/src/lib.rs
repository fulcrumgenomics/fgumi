#![deny(unsafe_code)]

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
