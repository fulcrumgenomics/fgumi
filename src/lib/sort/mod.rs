//! High-performance BAM sorting module.
//!
//! This module provides efficient BAM file sorting with support for multiple sort orders:
//! - **Template-coordinate**: Groups paired-end reads by template position (for `fgumi group`)
//! - **Queryname**: Groups reads by read name (for `fgumi zipper`)
//! - **Coordinate**: Standard genomic coordinate order (for IGV, `fgumi review`)
//!
//! # Performance Features
//!
//! - **External merge-sort**: Handles BAM files larger than available RAM via spill-to-disk
//! - **Lazy decoding**: Only parses fields needed for sort key extraction
//! - **Parallel sorting**: Uses rayon for in-memory parallel sort
//! - **Buffer recycling**: Reuses buffers via channel-based allocation patterns
//! - **Fast compression**: Uses libdeflate or ISA-L for temporary file compression
//!
//! # Architecture
//!
//! The sorting process follows this pipeline:
//!
//! 1. **Read phase**: Stream BAM records, extract sort keys lazily
//! 2. **Accumulate phase**: Buffer records until memory limit reached
//! 3. **Sort phase**: Parallel sort in-memory records using rayon
//! 4. **Spill phase**: Compress and write sorted chunk to temp file
//! 5. **Merge phase**: K-way merge of sorted temp files using min-heap

pub mod bam_fields;
pub mod external;
pub mod inline_buffer;
pub mod keys;
pub mod pipeline;
pub mod radix;
pub mod raw;
pub mod read_ahead;

pub use external::ExternalSorter;
pub use inline_buffer::TemplateKey;
pub use keys::{
    CoordinateKey, QuerynameKey, RawCoordinateKey, RawQuerynameKey, RawSortKey, SortContext,
    SortKey, SortOrder, TemplateCoordinateKey,
};
pub use pipeline::{ParallelMergeConfig, parallel_merge, parallel_merge_buffered};
pub use raw::RawExternalSorter;
