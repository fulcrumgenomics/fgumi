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
//! - **Fast compression**: Uses zstd (default) or BGZF/libdeflate for temporary file compression
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

// Unsafe code is prohibited at the crate level. The sole approved exception is
// `memory_probe::platform_ffi`, which is guarded by `#[allow(unsafe_code)]` on
// that inner module. See CLAUDE.md "Approved non-stdlib FFI exceptions".
#![deny(unsafe_code)]
#![deny(missing_docs)]

use std::path::Path;

use anyhow::{Context, Result};
use bstr::BString;
use noodles::sam::Header;
use noodles::sam::header::record::value::Map;
use noodles::sam::header::record::value::map::header::tag as header_tag;
use tempfile::TempDir;

// All sub-modules are crate-private. Items intended for external consumers are
// re-exported at the crate root below.
pub mod arena_pool;
pub(crate) mod bgzf_io;
pub(crate) mod block_offsets;
pub(crate) mod chunk_sorter;
pub mod codec;
pub(crate) mod external;
pub(crate) mod inline;
pub(crate) mod keys;
pub(crate) mod loser_tree;
pub(crate) mod memory_probe;
pub(crate) mod merge_slots;
pub(crate) mod pipeline;
pub(crate) mod pooled_bam_writer;
pub(crate) mod pooled_chunk_writer;
pub(crate) mod radix;
pub(crate) mod read_ahead;
pub(crate) mod reader;
pub mod ref_sort;
pub mod segmented_buf;
pub mod template_arena;
// Block-granular spill compression kernel for the block-parallel spill steps
// (`SpillGather`/`SpillBlockCompress`/`SpillWrite`). Surfaced via the re-exports
// below.
pub(crate) mod spill_block;
pub(crate) mod spill_block_reader;
// Synchronous inline spill compress for the P6 `CompressSpill` step (retires the
// `SortWorkerPool` compress path). Surfaced via `write_sorted_chunk` below.
pub(crate) mod sync_spill_writer;
pub(crate) mod tmp_dir_alloc;
pub(crate) mod verify;
pub(crate) mod worker_pool;
pub(crate) mod zspill_stream;

/// Print mimalloc allocator statistics to stderr.
///
/// This is a thin re-export of the internal FFI call, exposed for the main
/// `fgumi` crate's memory-debug monitor loop. Only available when the
/// `memory-debug` feature is enabled.
#[cfg(feature = "memory-debug")]
pub use memory_probe::print_mi_stats;

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
///
/// # Panics
///
/// Panics only if constructing a default `@HD` map fails for an input header that
/// has no `@HD` line — this cannot happen for the empty-map default used here.
#[must_use]
pub fn create_output_header(sort_order: keys::SortOrder, header: &Header) -> Header {
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

pub use arena_pool::{ArenaPool, PooledSegmentedBuf};
pub use chunk_sorter::{CoordinateChunkSorter, TemplateChunkSorter};
pub use codec::SpillCodec;
pub use external::MemorySources;
pub use external::{
    KeyTypesSpec, LibraryLookup, MergeDriver, MergeDriverDyn, MergeStep, RawExternalSorter,
    cb_hasher, extract_template_key_inline, open_spill_slot,
};
pub use inline::{
    CbKey32, InMemoryChunk, PackedCoordinateKey, RecordBuffer, RecordRef, SORT_SEGMENT_SIZE,
    TemplateKey, TemplateKey24, TertKey32, extract_coordinate_key_inline, radix_sort_record_refs,
    radix_sort_record_refs_with_max,
};
pub use keys::{
    QuerynameComparator, RawCoordinateKey, RawQuerynameKey, RawQuerynameLexKey, RawSortKey,
    SortContext, SortOrder, natural_compare, natural_compare_nul, normalize_natural_key,
};
pub use merge_slots::{PHASE2_DECOMP_CAP, SortMergeSlot};
pub use reader::RawBamRecordReader;
pub use ref_sort::{
    coordinate_chunk_from_arena_refs, coordinate_chunk_from_refs, queryname_chunk_from_arena_refs,
};
pub use segmented_buf::SegmentedBuf;
pub use spill_block::{SpillBlockCompressor, frame_keyed_record_into, spill_magic, spill_trailer};
pub use spill_block_reader::SpillBlockDecompressor;
pub use sync_spill_writer::{write_sorted_chunk, write_sorted_chunk_inmem};
pub use template_arena::{
    TemplateArenaAccumulator, TemplateMemChunk, template_chunk_from_arena_refs,
};
pub use tmp_dir_alloc::TmpDirAllocator;
pub use verify::{VerifySummary, verify_sort_order};

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
