//! Integration test: `output_compression` must reach the BAM the sorter writes,
//! on both the spilling and the fully-in-memory merge paths.
//!
//! Output BGZF blocks are compressed by the worker pool, and which compressor a
//! worker reaches for is chosen by the pool's phase. Phase 2 used to be entered
//! only when there were spill files to merge, so an in-memory sort left the pool
//! in Phase 1 — where the only compress step is the *spill* one, using
//! `temp_compression`. The user's `--compression-level` was silently replaced by
//! the temp level.
//!
//! The spilling merge had the mirror-image defect at the other end: it returned
//! the pool to `LEGACY` *before* finishing the writer, so the output blocks
//! still queued at teardown — plus the writer's final flush — went through the
//! spill compressor and the tail of the BAM came out at `temp_compression`.
//! Both defects are the same coupling seen from opposite sides: which
//! compressor a block gets is decided by the pool's phase when a worker pops
//! it, not by where the block came from.
//!
//! Hence the two tests below: the first asserts `output_compression` reaches
//! the file, the second that `temp_compression` never does.

use fgumi_raw_bam::{IndexedRawBamReader, RawRecord};
use fgumi_sam::SamBuilder;
use fgumi_sort::{
    LibraryLookup, QuerynameComparator, RawBamRecordReader, RawExternalSorter, RawQuerynameKey,
    RawQuerynameLexKey, RawSortKey, SortContext, SortOrder, SpillCodec, cb_hasher,
    extract_coordinate_key_inline, extract_template_key_inline, verify_sort_order,
};
use noodles::core::Region;
use rstest::rstest;

/// Number of pairs to sort. Large enough that the output spans many BGZF blocks,
/// so the compression level has a clearly measurable effect on file size.
const N_PAIRS: usize = 5_000;

/// Stride used to emit the pairs in a scrambled but deterministic order. Coprime
/// with [`N_PAIRS`] (5000 = 2^3·5^4, 2001 = 3·23·29), so `i * STRIDE % N_PAIRS`
/// walks every index exactly once.
const EMIT_STRIDE: usize = 2001;

/// Memory limit small enough that the sort has to spill chunk files.
const SPILL_MEMORY_LIMIT: usize = 1024 * 1024;

/// Builds a BAM whose records compress well, so level 0 and level 9 differ a lot.
///
/// Pairs are emitted in scrambled index order so the input is *not* already in
/// any of the sort orders under test. The sorter therefore has to do real work,
/// and the order check in [`sort_at_level`] cannot pass just by echoing the
/// input.
fn build_input(path: &std::path::Path) {
    let mut builder = SamBuilder::new();
    for i in 0..N_PAIRS {
        let pair = i * EMIT_STRIDE % N_PAIRS;
        let _ = builder
            .add_pair()
            .name(&format!("read{pair:06}"))
            .start1(pair * 20 + 1)
            .start2(pair * 20 + 101)
            .build();
    }
    builder.write_bam(path).expect("write_bam");
}

/// Decodes every record from a BAM into an owned, in-order `Vec<RawRecord>`,
/// re-reading through `RawBamRecordReader` rather than the writer that produced
/// the file.
fn decode_records(path: &std::path::Path) -> Vec<RawRecord> {
    let file = std::fs::File::open(path).expect("open bam");
    let mut reader = RawBamRecordReader::new(file).expect("bam reader");
    reader.skip_header().expect("skip header");
    let mut records = Vec::new();
    while let Some(record) = reader.next_record().expect("read record") {
        records.push(record);
    }
    records
}

/// Canonicalises a decoded record stream into a comparable multiset: the record
/// payloads, ordered by their bytes so file order does not matter.
fn record_multiset(records: &[RawRecord]) -> Vec<Vec<u8>> {
    let mut payloads: Vec<Vec<u8>> = records.iter().map(|record| record.to_vec()).collect();
    payloads.sort_unstable();
    payloads
}

/// Reads the BAI sidecar written next to `path` and uses it to query the file,
/// returning `(index_bytes, records_resolved)`.
///
/// The query covers all of `chr1` — the only reference [`build_input`] places
/// records on — so `records_resolved` is the whole file, record for record,
/// when the index is correct. Going through the index (rather than just loading
/// it) is what makes this a check on the *virtual offsets*: a `.bai` whose
/// offsets point into the wrong BGZF blocks still parses, but seeking to them
/// yields short reads or garbage, so records go missing or the read errors out.
/// Callers compare the returned records against a linear scan of the same file,
/// so a query that merely returns the right *number* of records is not enough.
///
/// The raw index bytes are returned alongside so two runs can be compared
/// directly: identical output blocks must produce an identical index.
fn query_via_index(path: &std::path::Path) -> (Vec<u8>, Vec<RawRecord>) {
    let bai_path = fgumi_bam_io::bai_sidecar_path(path);
    // Parsed from the bytes already in hand rather than re-reading the sidecar.
    let index_bytes = std::fs::read(&bai_path)
        .unwrap_or_else(|e| panic!("expected a BAI sidecar at {}: {e}", bai_path.display()));
    let index = noodles::bam::bai::io::Reader::new(&index_bytes[..])
        .read_index()
        .expect("BAI must be loadable");

    let mut reader =
        IndexedRawBamReader::from_path(path, index).expect("open indexed reader over sorted BAM");
    let header = reader.read_header().expect("read header through the indexed reader");
    let region: Region = "chr1".parse().expect("chr1 is a valid region");
    let resolved: Vec<RawRecord> = reader
        .query(&header, &region)
        .expect("query chr1")
        .map(|record| record.expect("indexed query must resolve to a readable record"))
        .collect();

    (index_bytes, resolved)
}

/// Counts sort-order violations in `path` under `sort_order`.
///
/// This runs the exported verifier (`verify_sort_order` plus the public key
/// extractors) over the finished file — the same check `fgumi sort --verify`
/// performs, and a different code path from the sort engine that produced it.
fn count_sort_violations(sort_order: SortOrder, path: &std::path::Path) -> u64 {
    let (_, header) = fgumi_bam_io::create_bam_reader(path, 1).expect("open bam for header");
    let new_reader = || {
        let file = std::fs::File::open(path).expect("open bam");
        let mut reader = RawBamRecordReader::new(file).expect("bam reader");
        reader.skip_header().expect("skip header");
        reader
    };

    let (_, violations, _) = match sort_order {
        SortOrder::Coordinate => {
            let n_ref = u32::try_from(header.reference_sequences().len())
                .expect("reference sequence count fits in u32");
            verify_sort_order(
                new_reader(),
                |bam| extract_coordinate_key_inline(bam, n_ref),
                |key, prev| key < prev,
            )
        }
        SortOrder::Queryname(QuerynameComparator::Natural) => {
            let ctx = SortContext::from_header(&header);
            verify_sort_order(
                new_reader(),
                |bam| RawQuerynameKey::extract(bam, &ctx),
                |key, prev| key < prev,
            )
        }
        SortOrder::Queryname(QuerynameComparator::Lexicographic) => {
            let ctx = SortContext::from_header(&header);
            verify_sort_order(
                new_reader(),
                |bam| RawQuerynameLexKey::extract(bam, &ctx),
                |key, prev| key < prev,
            )
        }
        SortOrder::TemplateCoordinate => {
            let lib_lookup = LibraryLookup::from_header(&header);
            let hasher = cb_hasher();
            verify_sort_order(
                new_reader(),
                |bam| extract_template_key_inline(bam, &lib_lookup, None, &hasher),
                |key, prev| key < prev,
            )
        }
    }
    .expect("verify sort order");
    violations
}

/// The knobs a single [`sort_at_level`] run varies.
///
/// `write_index` picks the finalization path under test: `false` routes the
/// merge through `merge_chunks_generic` (`writer.finish()`), `true` through
/// `merge_chunks_with_index` (`writer.finish_index()`), which additionally pins
/// BAI virtual offsets to the block layout the tail compression produces.
#[derive(Clone, Copy)]
struct SortRun {
    memory_limit: usize,
    output_compression: u32,
    temp_compression: u32,
    write_index: bool,
}

/// Sorts `input` in `sort_order` under `run` and returns
/// `(output_size_bytes, chunks_written, decoded_records)`.
///
/// `expected` is the independent baseline the output is checked against — the
/// multiset of the *input* records, decoded straight from the input BAM without
/// going through the sorter. Combined with the sort-order verification below,
/// it pins the output exactly: a correct sort emits every input record, no
/// others, in non-decreasing key order. Neither check involves any other sort
/// run, so two runs cannot drift the same way.
fn sort_at_level(
    sort_order: SortOrder,
    input: &std::path::Path,
    output: &std::path::Path,
    run: SortRun,
    expected: &[Vec<u8>],
) -> (u64, usize, Vec<RawRecord>) {
    let SortRun { memory_limit, output_compression, temp_compression, write_index } = run;

    let stats = RawExternalSorter::new(sort_order)
        .threads(2)
        .memory_limit(memory_limit)
        .spill_codec(SpillCodec::Bgzf)
        .temp_compression(temp_compression)
        .output_compression(output_compression)
        .write_index(write_index)
        .sort(input, output)
        .expect("sort should succeed");

    assert_eq!(
        stats.output_records,
        (N_PAIRS * 2) as u64,
        "all records must survive the sort at level {output_compression}"
    );
    let records = decode_records(output);
    // Independent of the sorter's own `output_records` bookkeeping: count the
    // records that actually decode out of the written BAM.
    assert_eq!(
        records.len(),
        N_PAIRS * 2,
        "decoded record count must match the input at level {output_compression}"
    );
    // Content: exactly the input records, byte for byte, with nothing dropped,
    // duplicated or corrupted by the output compressor.
    assert!(
        record_multiset(&records) == expected,
        "level {output_compression} output is not the input record set ({sort_order:?})"
    );
    // Order: the output really is sorted, so the content check above cannot be
    // satisfied by an arbitrary permutation.
    assert_eq!(
        count_sort_violations(sort_order, output),
        0,
        "level {output_compression} output is not in {sort_order:?} order"
    );
    let size = std::fs::metadata(output).expect("output metadata").len();
    (size, stats.chunks_written, records)
}

/// `output_compression` must change the size of the written BAM.
///
/// Both merge paths are covered: a memory limit large enough to hold everything
/// (no spill, `chunks_written == 0`) and one small enough to force spill files.
/// Only the in-memory case was broken — it never entered Phase 2, so the pool
/// compressed output blocks with the Phase 1 spill compressor.
///
/// Every sort order that routes its in-memory output through `PooledBamWriter`
/// is exercised, because the `begin_phase2` transition is duplicated per order
/// and a regression could reintroduce the bug in just one of them. That is:
/// `Coordinate` (`sort_coordinate_optimized`), `Queryname`
/// (`sort_queryname_keyed`), and `TemplateCoordinate`
/// (`sort_template_coordinate_impl`). The `--write-index` coordinate path
/// (`sort_coordinate_with_index`) is deliberately excluded: its in-memory
/// branch bypasses the pool and writes at `output_compression` directly, so it
/// is not at risk here.
///
/// The payload is checked against an oracle built independently of both
/// compression runs — the input record multiset, plus `verify_sort_order` over
/// the finished file — so a routing change that silently mangled records could
/// not pass by mangling them the same way at both levels. See
/// [`sort_at_level`].
#[rstest]
fn test_output_compression_level_reaches_the_output_bam(
    #[values(
        SortOrder::Coordinate,
        SortOrder::Queryname(QuerynameComparator::Natural),
        SortOrder::TemplateCoordinate
    )]
    sort_order: SortOrder,
    #[values((256 * 1024 * 1024, false), (SPILL_MEMORY_LIMIT, true))] memory: (usize, bool),
) {
    let (memory_limit, expect_spill) = memory;
    let dir = tempfile::tempdir().expect("tempdir");
    let input = dir.path().join("input.bam");
    build_input(&input);

    // The baseline both outputs are checked against: the input records, decoded
    // from the input BAM before any sort runs.
    let expected = record_multiset(&decode_records(&input));

    // The temp level is held constant across both runs, and is deliberately
    // different from both output levels under test: if the output were
    // compressed at the temp level, both runs would produce the same size.
    let run_at = |output_compression| SortRun {
        memory_limit,
        output_compression,
        temp_compression: 1,
        write_index: false,
    };

    let uncompressed = dir.path().join("level0.bam");
    let (size_level_0, chunks_0, records_level_0) =
        sort_at_level(sort_order, &input, &uncompressed, run_at(0), &expected);

    let compressed = dir.path().join("level9.bam");
    let (size_level_9, chunks_9, records_level_9) =
        sort_at_level(sort_order, &input, &compressed, run_at(9), &expected);

    assert_eq!(
        chunks_0 > 0,
        expect_spill,
        "expected spill={expect_spill}, got chunks_written={chunks_0} at level 0 ({sort_order:?})"
    );
    assert_eq!(
        chunks_9 > 0,
        expect_spill,
        "expected spill={expect_spill}, got chunks_written={chunks_9} at level 9 ({sort_order:?})"
    );

    // Both runs were already checked against the independent baseline above;
    // this additionally pins tie ordering, which "same records, sorted" leaves
    // free: the two runs differ only in the output compressor, so the compressor
    // must not perturb the emitted sequence at all.
    assert_eq!(
        records_level_0, records_level_9,
        "level 0 and level 9 output must decode to identical records in order ({sort_order:?}, \
         spill={expect_spill})"
    );

    assert!(
        size_level_9 < size_level_0,
        "output_compression was ignored for {sort_order:?}: level 9 wrote {size_level_9} bytes and \
         level 0 wrote {size_level_0} bytes (spill={expect_spill}). Equal sizes mean the output was \
         compressed at temp_compression instead."
    );
}

/// `temp_compression` must not reach the output BAM.
///
/// The mirror image of the check above, at the other end of the sort: see this
/// file's module header for the defect. The output must be finalized while the
/// pool is still in Phase 2, or the trailing blocks go through the *temp*
/// compressor — so raising `temp_compression` alone would change the output.
///
/// So: sort twice at `output_compression = 0`, changing only the temp level.
/// The written BAM must be exactly the same size both times. Level 0 is used
/// for the output because it makes any block that slipped through the temp
/// compressor strictly smaller, and therefore visible in the total.
///
/// Only the spilling path is exercised — the in-memory path never builds a
/// `Phase2Guard`, so it has no teardown to get wrong.
///
/// Both spilling finalizations are covered: `merge_chunks_generic`
/// (`writer.finish()`) and, via the `coordinate_indexed` case,
/// `merge_chunks_with_index` (`writer.finish_index()`). Each sort order is
/// covered because each routes to the merge through its own entry point. The
/// indexed path carries the extra risk that the tail block's compression level
/// moves the BGZF boundaries the BAI records, so that case additionally pins the
/// index itself — see [`query_via_index`]. `write_index` is only meaningful for
/// `SortOrder::Coordinate`; the other orders never reach
/// `sort_coordinate_with_index`, so pairing them with it would just re-run the
/// unindexed case.
#[rstest]
#[case::coordinate(SortOrder::Coordinate, false)]
#[case::coordinate_indexed(SortOrder::Coordinate, true)]
#[case::queryname_natural(SortOrder::Queryname(QuerynameComparator::Natural), false)]
#[case::template_coordinate(SortOrder::TemplateCoordinate, false)]
fn test_temp_compression_does_not_reach_the_output_bam(
    #[case] sort_order: SortOrder,
    #[case] write_index: bool,
) {
    let dir = tempfile::tempdir().expect("tempdir");
    let input = dir.path().join("input.bam");
    build_input(&input);
    let expected = record_multiset(&decode_records(&input));

    // Returns the output size and, for an indexed run, the BAI bytes. The BAI is
    // checked here rather than by the caller because this is where `scanned` —
    // the linear-scan oracle `sort_at_level` already produced — is in hand.
    let sort_with_temp = |tag: &str, temp_compression: u32| -> (u64, Option<Vec<u8>>) {
        let output = dir.path().join(format!("{tag}.bam"));
        let run = SortRun {
            memory_limit: SPILL_MEMORY_LIMIT,
            output_compression: 0,
            temp_compression,
            write_index,
        };
        let (size, chunks, scanned) = sort_at_level(sort_order, &input, &output, run, &expected);
        assert!(
            chunks > 0,
            "test must exercise the spilling path, but no chunks were written \
             ({sort_order:?}, temp_compression={temp_compression})"
        );

        // The index is built from the virtual positions of the blocks the pool
        // compressed, so a tail block that slipped through the temp compressor
        // would move every offset after it.
        //
        // The oracle for "resolved correctly" is `scanned`, a linear scan of the
        // same BAM that does not consult the index at all: an indexed query over
        // the only reference in the file must yield exactly those records, in
        // that order. Comparing counts alone would pass on an index that pointed
        // at the right number of wrong records.
        let bai = write_index.then(|| {
            let (bai, queried) = query_via_index(&output);
            assert!(
                queried == scanned,
                "indexed query at temp_compression={temp_compression} resolved {} records against \
                 {} in the file: the BAI's virtual offsets do not match the blocks in the BAM",
                queried.len(),
                scanned.len()
            );
            bai
        });

        (size, bai)
    };

    let (size_temp_0, bai_temp_0) = sort_with_temp("temp0", 0);
    let (size_temp_9, bai_temp_9) = sort_with_temp("temp9", 9);

    assert_eq!(
        size_temp_0, size_temp_9,
        "temp_compression leaked into the output for {sort_order:?}: the same sort wrote \
         {size_temp_0} bytes at temp_compression=0 and {size_temp_9} bytes at \
         temp_compression=9, but only output_compression (0 for both) may affect the output"
    );

    // Identical output blocks must produce an identical index. Compared as raw
    // bytes rather than with `assert_eq!`, which would dump both whole indexes
    // into the failure message.
    assert!(
        bai_temp_0 == bai_temp_9,
        "temp_compression changed the BAI for {sort_order:?} ({:?} and {:?} bytes): the indexed \
         merge compressed part of the output at the temp level, shifting the virtual offsets the \
         index records",
        bai_temp_0.as_ref().map(Vec::len),
        bai_temp_9.as_ref().map(Vec::len),
    );
}
