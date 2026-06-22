//! Re-export shim. The `ReadBgzfBlocks` source step and `read_bam*` helpers
//! now live in the `fgumi-pipeline-io` crate.
pub use fgumi_pipeline_io::source::read_bam::{
    DEFAULT_BLOCKS_PER_BATCH, ReadBgzfBlocks, read_bam, read_bam_auto, read_bam_from_reader,
    read_bam_stdin,
};
