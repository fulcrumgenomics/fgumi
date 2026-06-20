//! Re-export shim. The sort typed-steps now live in the `fgumi-pipeline-io`
//! crate (`fgumi_pipeline_io::sort`).
pub use fgumi_pipeline_io::sort::{
    SortAndSpill, SortBamFile, SortMerge, SortSpillDecompress, protocol,
};

pub mod and_spill {
    pub use fgumi_pipeline_io::sort::and_spill::*;
}
pub mod merge {
    pub use fgumi_pipeline_io::sort::merge::*;
}
pub mod spill_decompress {
    pub use fgumi_pipeline_io::sort::spill_decompress::*;
}
