//! Engine components for sorting and merging BAM records.
//!
//! The pieces `fgumi sort`, `fgumi merge` and `fgumi simulate` are built from:
//! sort keys and their inline extractors (coordinate, queryname lexicographic
//! and natural, template-coordinate with its `--key-types` lanes), the
//! in-memory record buffers, arenas and parallel radix sorts, the spill-file
//! block format and its decompressor, the merge slots the spill read planner
//! feeds, the k-way merge driver and its loser tree, and a sort-order verifier.
//! The sort itself is a chain-native pipeline in the `fgumi` crate
//! (`fgumi_lib::sort::sort_bam`, or the `fgumi sort` command); this crate does
//! not ship a standalone sorter.

// Unsafe code is prohibited at the crate level. Approved exceptions carry
// #[allow(unsafe_code)] on the item: memory_probe::platform_ffi (CLAUDE.md
// "Approved non-stdlib FFI exceptions") and the sort hot paths listed in CLAUDE.md
// "Approved hot-path unsafe (sort engine)".
#![deny(unsafe_code)]
#![deny(missing_docs)]

use bstr::BString;
use noodles::sam::Header;
use noodles::sam::header::record::value::Map;
use noodles::sam::header::record::value::map::header::tag as header_tag;

// Sub-modules are crate-private except where a consumer needs to name the type
// itself rather than just receive it (`arena_pool`, `segmented_buf`, `codec`,
// `prefetch`, `ref_sort`, `sort_pool`, `template_arena`, `template_key`).
// Everything else is re-exported at the crate root below.
pub mod arena_pool;
pub(crate) mod chunk_sorter;
pub mod codec;
pub(crate) mod fd_limit;
pub(crate) mod inline;
pub(crate) mod keys;
pub(crate) mod loser_tree;
pub(crate) mod memory_probe;
// The merge's demand facts (awaited file, consumer wake, demand counters),
// shared by the spill supply and the merge. Surfaced via the re-export below.
pub(crate) mod merge_demand;
// The k-way merge driver and slot parser. Surfaced via the re-export below.
pub(crate) mod merge_driver;
pub(crate) mod merge_slots;
pub mod prefetch;
pub(crate) mod queryname_ref;
pub(crate) mod radix;
pub(crate) mod read_ahead;
pub(crate) mod reader;
pub mod ref_sort;
// A merge or sort output staged beside its destination and persisted on
// success. Surfaced via the re-export below.
pub(crate) mod output_target;
// Per-order run boundary key for the arena spill run-former (extend-vs-new-run
// decision). Surfaced via the re-export below.
pub(crate) mod run_bound;
pub mod segmented_buf;
pub mod sort_pool;
pub mod template_arena;
pub mod template_key;
// Block-granular spill compression kernel for the block-parallel spill steps
// (`SpillGather`/`SpillCompress`/`SpillWrite`). Surfaced via the re-exports
// below.
pub(crate) mod spill_block;
pub(crate) mod spill_block_reader;
pub(crate) mod tmp_dir_alloc;
pub(crate) mod verify;

/// Print mimalloc allocator statistics to stderr.
///
/// This is a thin re-export of the internal FFI call, exposed for the main
/// `fgumi` crate's memory-debug monitor loop. Only available when the
/// `memory-debug` feature is enabled.
#[cfg(feature = "memory-debug")]
pub use memory_probe::print_mi_stats;
/// Sample the current process's resident set size, in bytes.
///
/// Thin re-export of the internal RSS probe (reads `VmRSS` from
/// `/proc/self/status` on Linux — no FFI; `mach2::task_info(TASK_VM_INFO)` FFI on
/// macOS; `None` elsewhere), exposed unconditionally — unlike `print_mi_stats`,
/// this is not gated behind `memory-debug` — for the main `fgumi` crate's
/// pipeline thread-telemetry sampler, which needs an RSS reading on every tick
/// regardless of whether that feature is enabled.
pub use memory_probe::process_rss_bytes;
/// Tune mimalloc to keep freed pages instead of returning them to the OS, and
/// read the setting back. Re-exported for the in-process aligner, which frees
/// and reallocates large per-batch buffers on every pool thread.
pub use memory_probe::{RETAINED_PURGE_DELAY_MS, mi_purge_delay_ms, retain_freed_memory};
/// Background read-ahead record reader, re-exported for `fgumi compare bams`,
/// which reads two inputs concurrently and needs each decode off the main thread.
pub use read_ahead::RawReadAheadReader;

/// Statistics from a sort operation.
#[derive(Default, Debug, Clone, Copy, PartialEq, Eq)]
pub struct SortStats {
    /// Total records read from input.
    pub total_records: u64,
    /// Records written to output.
    pub output_records: u64,
    /// Number of spill *runs* written, before any `--max-temp-files`
    /// consolidation merged some of them together.
    ///
    /// A run is one temp file, built from one or more sorted chunks: consecutive
    /// chunks that are already in order extend the open run rather than each
    /// starting a file. So an input already in the requested order spills a single
    /// run, and this is below the number of chunks that were spilled.
    pub runs_written: usize,
}

/// Create an output header with appropriate sort order tags (SO, GO, SS).
///
/// Preserves all existing header content (reference sequences, read groups, programs,
/// comments, and `@HD` fields like `VN`), then overwrites only the sort-related tags
/// (`SO`, `GO`, `SS`) based on the requested sort order.
///
/// # Panics
///
/// Panics only if constructing the default `@HD` map fails, which cannot happen for
/// the hard-coded default built here.
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
            if let Some(ss) = sort_order.header_ss_tag() {
                hd.other_fields_mut().insert(header_tag::SUBSORT_ORDER, BString::from(ss));
            }
        }
    }

    builder = builder.set_header(hd);
    builder.build()
}

pub use arena_pool::{ArenaPool, PooledSegmentedBuf};
pub use chunk_sorter::{CoordinateChunkSorter, TemplateChunkSorter};
pub use codec::SpillCodec;
pub use fd_limit::{
    FALLBACK_MAX_TEMP_FILES, FD_RESERVE, fits_nofile_budget, soft_nofile,
    temp_file_limit_from_nofile,
};
pub use inline::{
    CbKey32, ChunkRef, InMemoryChunk, PackedCoordinateKey, RecordBuffer, RecordRef,
    SORT_SEGMENT_SIZE, TemplateKey, TemplateKey24, TemplateLaneKey, TemplateRecordRef, TertKey32,
    extract_coordinate_key_inline, radix_sort_record_refs, radix_sort_record_refs_with_max,
};
/// The template run sort's kernels, for the `seal_allocations` test.
#[cfg(feature = "test-utils")]
pub use inline::{radix_sort_template_refs, sort_template_refs_into_records};
pub use keys::{
    QuerynameComparator, RawCoordinateKey, RawQuerynameKey, RawQuerynameLexKey, RawSortKey,
    SortContext, SortOrder, natural_compare, natural_compare_nul, normalize_natural_key,
};
pub use merge_demand::{
    AwaitedSlotState, HotIds, MAX_HOT_IDS, MergeDemand, MergeDemandSnapshot, MergeDemandStats,
    MergePoolRequest,
};
pub use merge_driver::{
    MergeCtxErased, MergeDriver, MergeDriverDyn, MergeStep, MergeStepFramed, RecordRewrite,
    input_order_violation, merge_driver_for_kind, open_bam_slot, open_spill_slot,
    open_spill_slot_as,
};
pub use merge_slots::{
    CONSOLIDATION_SLOT_IDS_START, ClaimedBlock, IngestOutcome, PHASE2_DECOMP_CAP, SortMergeSlot,
};
pub use output_target::MergeOutputTarget;
pub use prefetch::KEY_PREFETCH_DISTANCE;
pub use queryname_ref::{
    QN_CLASS_BATCH_RECORDS, QN_MAX_CLASSES, QnBatchLcp, QnClassTable, QnCmpStats, QnComparator,
    QnKey, QnLex, QnMergeKey, QnNatural, QnRef, QnRefBuilder, QnSortState, cmp_qn_refs,
    queryname_chunk_from_qn_refs, queryname_chunks_from_qn_refs,
};
pub use reader::{
    OwnedRawBamRecordReader, RawBamRecordReader, open_raw_bam_record_reader,
    open_raw_bam_record_reader_with_header,
};
pub use ref_sort::{
    coordinate_chunk_from_arena_refs, coordinate_chunk_from_refs, coordinate_chunks_from_refs,
    radix_sorts_in_parallel,
};
pub use run_bound::{RunBound, format_run_formation};
pub use segmented_buf::SegmentedBuf;
pub use sort_pool::BoundedSortPool;
pub use spill_block::{
    SpillBlockCompressor, SpillKeyKind, frame_embedded_record_into, frame_keyed_record_into,
    spill_magic, spill_trailer,
};
pub use spill_block_reader::{RawBlock, SpillBlockDecompressor};
pub use template_arena::{
    TemplateArenaAccumulator, TemplateArenaRefs, TemplateMemChunk, seal_template_refs,
    seal_template_refs_pieces, template_chunk_from_arena_refs,
};
pub use template_key::{
    DroppedLaneViolation, KeyTypesSpec, LibraryLookup, TemplateKeyVariant, cb_hasher,
    extract_template_key_inline,
};
pub use template_key::{ExtractedRefs, KeyViolation, TemplateKeyContext, TemplateKeySeed};
pub use tmp_dir_alloc::TmpDirAllocator;
pub use tmp_dir_alloc::create_sort_temp_dirs;
pub use verify::{VerifySummary, verify_sort_order};

/// Render the effective per-phase worker counts for a log line.
///
/// Collapses to a single number when the phases agree (the common case, where
/// neither `--sort-threads` nor `--merge-threads` was overridden) and spells
/// both out otherwise. Shared by the engine's end-of-run summary and the `sort`
/// command's config dump so the two emitters cannot drift.
#[must_use]
pub fn format_thread_counts(sort_threads: usize, merge_threads: usize) -> String {
    if sort_threads == merge_threads {
        format!("{sort_threads}")
    } else {
        format!("sort {sort_threads}, merge {merge_threads}")
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_sort_stats_default() {
        let stats = SortStats::default();
        assert_eq!(stats.total_records, 0);
        assert_eq!(stats.output_records, 0);
        assert_eq!(stats.runs_written, 0);
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

    #[test]
    fn test_create_output_header_ss_has_sort_order_prefix() {
        // fgbio/samtools write SS as `<sort-order>:<sub-sort>`; pin that fgumi
        // does too, so fgbio's `SamOrder.apply` re-recognizes the header.
        let empty = Header::builder().build();
        let ss_of = |hd: &Map<noodles::sam::header::record::value::map::Header>| {
            <_ as AsRef<[u8]>>::as_ref(hd.other_fields().get(b"SS").expect("SS")).to_vec()
        };

        // queryname (natural) -> SO:queryname, SS:queryname:natural
        let out = create_output_header(
            keys::SortOrder::Queryname(keys::QuerynameComparator::Natural),
            &empty,
        );
        let hd = out.header().expect("@HD");
        assert_eq!(
            <_ as AsRef<[u8]>>::as_ref(hd.other_fields().get(b"SO").expect("SO")),
            b"queryname"
        );
        assert_eq!(ss_of(hd), b"queryname:natural");

        // template-coordinate -> SO:unsorted, GO:query, SS:unsorted:template-coordinate
        let out = create_output_header(keys::SortOrder::TemplateCoordinate, &empty);
        let hd = out.header().expect("@HD");
        assert_eq!(
            <_ as AsRef<[u8]>>::as_ref(hd.other_fields().get(b"SO").expect("SO")),
            b"unsorted"
        );
        assert_eq!(<_ as AsRef<[u8]>>::as_ref(hd.other_fields().get(b"GO").expect("GO")), b"query");
        assert_eq!(ss_of(hd), b"unsorted:template-coordinate");
    }
}
