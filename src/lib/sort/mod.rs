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
//! - **Fast compression**: Uses libdeflate for temporary file compression
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

use std::path::Path;

use anyhow::{Context, Result};
use tempfile::TempDir;

pub use fgumi_raw_bam as bam_fields;
pub mod external;
pub mod inline_buffer;
pub mod keys;
pub mod pipeline;
pub mod radix;
pub mod raw;
pub mod raw_bam_reader;
pub mod read_ahead;

/// Buffer size for `BufReader` during merge phase.
const MERGE_BUFFER_SIZE: usize = 64 * 1024;

/// Statistics from a sort operation.
#[derive(Default, Debug)]
pub struct SortStats {
    /// Total records read from input.
    pub total_records: u64,
    /// Records written to output.
    pub output_records: u64,
    /// Number of temporary chunk files written.
    pub chunks_written: usize,
}

/// Create a temporary directory for sort spill files.
fn create_temp_dir(base: Option<&Path>) -> Result<TempDir> {
    match base {
        Some(base) => {
            std::fs::create_dir_all(base)?;
            TempDir::new_in(base).context("Failed to create temp directory")
        }
        None => TempDir::new().context("Failed to create temp directory"),
    }
}

pub use external::ExternalSorter;
pub use inline_buffer::{TemplateKey, extract_coordinate_key_inline};
pub use keys::{
    CoordinateKey, PA_TAG, PrimaryAlignmentInfo, QuerynameKey, RawCoordinateKey, RawQuerynameKey,
    RawSortKey, SortContext, SortKey, SortOrder, TemplateCoordinateKey,
};
pub use pipeline::{ParallelMergeConfig, parallel_merge, parallel_merge_buffered};
pub use raw::{LibraryLookup, RawExternalSorter, extract_template_key_inline};

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_create_temp_dir_default() {
        let dir = create_temp_dir(None).unwrap();
        assert!(dir.path().exists());
    }

    #[test]
    fn test_create_temp_dir_with_base() {
        let base = tempfile::tempdir().unwrap();
        let subdir = base.path().join("sort_spill");
        let dir = create_temp_dir(Some(&subdir)).unwrap();
        assert!(dir.path().exists());
        assert!(dir.path().starts_with(&subdir));
    }

    #[test]
    fn test_sort_stats_default() {
        let stats = SortStats::default();
        assert_eq!(stats.total_records, 0);
        assert_eq!(stats.output_records, 0);
        assert_eq!(stats.chunks_written, 0);
    }
}
