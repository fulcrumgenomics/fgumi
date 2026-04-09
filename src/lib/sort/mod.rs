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
//! - **N+2 worker pool**: Dedicated reader, writer, and N sort workers stream
//!   chunks through the pipeline; in-memory sort uses parallel radix sort
//! - **Buffer recycling**: Reuses buffers via channel-based allocation patterns
//! - **Fast compression**: Uses libdeflate for temporary file compression
//!
//! # Architecture
//!
//! The sorting process follows this pipeline:
//!
//! 1. **Read phase**: Stream BAM records, extract sort keys lazily
//! 2. **Accumulate phase**: Buffer records until memory limit reached
//! 3. **Sort phase**: Parallel radix sort over the in-memory record buffer
//! 4. **Spill phase**: Compress and write sorted chunk to temp file
//! 5. **Merge phase**: K-way merge of sorted temp files using a loser tree

use std::path::Path;

use anyhow::{Context, Result};
use bstr::BString;
use noodles::sam::Header;
use noodles::sam::header::record::value::Map;
use noodles::sam::header::record::value::map::header::tag as header_tag;
use tempfile::TempDir;

pub use fgumi_raw_bam as bam_fields;
pub mod bgzf_io;
pub mod inline_buffer;
pub mod keys;
pub mod loser_tree;
pub(crate) mod memory_probe;
pub mod pipeline;
pub mod pooled_bam_writer;
pub mod pooled_chunk_writer;
pub mod radix;
pub mod raw;
pub mod raw_bam_reader;
pub mod read_ahead;
pub(crate) mod segmented_buf;
pub mod worker_pool;

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

/// Create an output header with appropriate sort order tags (SO, GO, SS).
///
/// Preserves all existing header content (reference sequences, read groups, programs,
/// comments, and `@HD` fields like `VN`), then overwrites only the sort-related tags
/// (`SO`, `GO`, `SS`) based on the requested sort order.
pub(crate) fn create_output_header(sort_order: keys::SortOrder, header: &Header) -> Header {
    let mut builder = Header::builder();

    // Copy reference sequences
    for (name, seq) in header.reference_sequences() {
        builder = builder.add_reference_sequence(name.as_slice(), seq.clone());
    }

    // Copy read groups
    for (id, rg) in header.read_groups() {
        builder = builder.add_read_group(id.as_slice(), rg.clone());
    }

    // Copy programs
    for (id, pg) in header.programs().as_ref() {
        builder = builder.add_program(id.as_slice(), pg.clone());
    }

    // Copy comments
    for comment in header.comments() {
        builder = builder.add_comment(comment.clone());
    }

    // Start from the existing @HD record (preserving VN and any other fields),
    // or create a fresh one if the input has no @HD line.
    let mut hd = header.header().cloned().unwrap_or_else(|| {
        Map::<noodles::sam::header::record::value::map::Header>::builder()
            .build()
            .expect("valid default header")
    });

    // Clear sort-related tags before setting new ones, so stale values from a
    // previous sort order don't leak through.
    hd.other_fields_mut().swap_remove(&header_tag::SORT_ORDER);
    hd.other_fields_mut().swap_remove(&header_tag::GROUP_ORDER);
    hd.other_fields_mut().swap_remove(&header_tag::SUBSORT_ORDER);

    match sort_order {
        keys::SortOrder::Coordinate => {
            hd.other_fields_mut().insert(header_tag::SORT_ORDER, BString::from("coordinate"));
        }
        keys::SortOrder::Queryname(_) => {
            hd.other_fields_mut().insert(header_tag::SORT_ORDER, BString::from("queryname"));
            if let Some(ss) = sort_order.header_ss_tag() {
                hd.other_fields_mut().insert(header_tag::SUBSORT_ORDER, BString::from(ss));
            }
        }
        keys::SortOrder::TemplateCoordinate => {
            hd.other_fields_mut().insert(header_tag::SORT_ORDER, BString::from("unsorted"));
            hd.other_fields_mut().insert(header_tag::GROUP_ORDER, BString::from("query"));
            hd.other_fields_mut()
                .insert(header_tag::SUBSORT_ORDER, BString::from("template-coordinate"));
        }
    }

    builder = builder.set_header(hd);
    builder.build()
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

pub use inline_buffer::{TemplateKey, extract_coordinate_key_inline};
pub use keys::{
    CoordinateKey, PA_TAG, PrimaryAlignmentInfo, QuerynameComparator, QuerynameKey,
    RawCoordinateKey, RawQuerynameKey, RawQuerynameLexKey, RawSortKey, SortContext, SortKey,
    SortOrder, natural_compare, normalize_natural_key,
};
pub use pipeline::{ParallelMergeConfig, parallel_merge, parallel_merge_buffered};
pub use raw::{LibraryLookup, RawExternalSorter, cb_hasher, extract_template_key_inline};

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_create_temp_dir_default() {
        let dir = create_temp_dir(None).expect("creating temp dir should succeed");
        assert!(dir.path().exists());
    }

    #[test]
    fn test_create_temp_dir_with_base() {
        let base = tempfile::tempdir().expect("creating temp file/dir should succeed");
        let subdir = base.path().join("sort_spill");
        let dir = create_temp_dir(Some(&subdir)).expect("creating temp dir should succeed");
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

    #[test]
    fn test_create_output_header_preserves_vn() {
        // Build a header with VN:1.6 in the @HD line.
        let hd = Map::<noodles::sam::header::record::value::map::Header>::new(
            noodles::sam::header::record::value::map::header::Version::new(1, 6),
        );
        let header = Header::builder().set_header(hd).build();

        let output = create_output_header(keys::SortOrder::Coordinate, &header);

        let hd = output.header().expect("should have @HD");
        assert_eq!(
            hd.version(),
            noodles::sam::header::record::value::map::header::Version::new(1, 6)
        );
        let so = hd.other_fields().get(b"SO").expect("should have SO");
        assert_eq!(<_ as AsRef<[u8]>>::as_ref(so), b"coordinate");
    }

    #[test]
    fn test_create_output_header_clears_stale_sort_tags() {
        // Start with template-coordinate tags (SO:unsorted, GO:query, SS:template-coordinate).
        let hd = Map::<noodles::sam::header::record::value::map::Header>::builder()
            .insert(header_tag::SORT_ORDER, BString::from("unsorted"))
            .insert(header_tag::GROUP_ORDER, BString::from("query"))
            .insert(header_tag::SUBSORT_ORDER, BString::from("template-coordinate"))
            .build()
            .expect("valid header");
        let header = Header::builder().set_header(hd).build();

        // Re-sort as coordinate — GO and SS should be removed.
        let output = create_output_header(keys::SortOrder::Coordinate, &header);

        let hd = output.header().expect("should have @HD");
        let so = hd.other_fields().get(b"SO").expect("should have SO");
        assert_eq!(<_ as AsRef<[u8]>>::as_ref(so), b"coordinate");
        assert!(hd.other_fields().get(b"GO").is_none(), "GO should be cleared");
        assert!(hd.other_fields().get(b"SS").is_none(), "SS should be cleared");
    }

    #[test]
    fn test_create_output_header_no_existing_hd() {
        // Header with no @HD line at all.
        let header = Header::builder().build();

        let output = create_output_header(keys::SortOrder::Coordinate, &header);

        let hd = output.header().expect("should have @HD");
        let so = hd.other_fields().get(b"SO").expect("should have SO");
        assert_eq!(<_ as AsRef<[u8]>>::as_ref(so), b"coordinate");
    }
}
