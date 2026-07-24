//! Integration test: `output_compression` must reach the BAM the sorter writes,
//! on both the spilling and the fully-in-memory merge paths.
//!
//! Output BGZF blocks are compressed by the worker pool, and which compressor a
//! worker reaches for is chosen by the pool's phase. Phase 2 used to be entered
//! only when there were spill files to merge, so an in-memory sort left the pool
//! in Phase 1 — where the only compress step is the *spill* one, using
//! `temp_compression`. The user's `--compression-level` was silently replaced by
//! the temp level.

use fgumi_raw_bam::RawRecord;
use fgumi_sam::SamBuilder;
use fgumi_sort::{
    LibraryLookup, QuerynameComparator, RawBamRecordReader, RawExternalSorter, RawQuerynameKey,
    RawQuerynameLexKey, RawSortKey, SortContext, SortOrder, SpillCodec, cb_hasher,
    extract_coordinate_key_inline, extract_template_key_inline, verify_sort_order,
};
use rstest::rstest;

/// Number of pairs to sort. Large enough that the output spans many BGZF blocks,
/// so the compression level has a clearly measurable effect on file size.
const N_PAIRS: usize = 5_000;

/// Stride used to emit the pairs in a scrambled but deterministic order. Coprime
/// with [`N_PAIRS`] (5000 = 2^3·5^4, 2001 = 3·23·29), so `i * STRIDE % N_PAIRS`
/// walks every index exactly once.
const EMIT_STRIDE: usize = 2001;

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

/// Sorts `input` in `sort_order` at the given output compression level and
/// returns `(output_size_bytes, chunks_written, decoded_records)`.
///
/// `expected` is the independent baseline the output is checked against — the
/// multiset of the *input* records, decoded straight from the input BAM without
/// going through the sorter. Combined with the sort-order verification below,
/// it pins the output exactly: a correct sort emits every input record, no
/// others, in non-decreasing key order. Neither check involves the other
/// compression run, so they cannot both drift the same way.
fn sort_at_level(
    sort_order: SortOrder,
    input: &std::path::Path,
    output: &std::path::Path,
    memory_limit: usize,
    output_compression: u32,
    expected: &[Vec<u8>],
) -> (u64, usize, Vec<RawRecord>) {
    let stats = RawExternalSorter::new(sort_order)
        .threads(2)
        .memory_limit(memory_limit)
        .spill_codec(SpillCodec::Bgzf)
        // Held constant across levels, and deliberately different from both
        // output levels under test: if the output is compressed at the temp
        // level, both runs produce the same size.
        .temp_compression(1)
        .output_compression(output_compression)
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
    #[values((256 * 1024 * 1024, false), (1024 * 1024, true))] memory: (usize, bool),
) {
    let (memory_limit, expect_spill) = memory;
    let dir = tempfile::tempdir().expect("tempdir");
    let input = dir.path().join("input.bam");
    build_input(&input);

    // The baseline both outputs are checked against: the input records, decoded
    // from the input BAM before any sort runs.
    let expected = record_multiset(&decode_records(&input));

    let uncompressed = dir.path().join("level0.bam");
    let (size_level_0, chunks_0, records_level_0) =
        sort_at_level(sort_order, &input, &uncompressed, memory_limit, 0, &expected);

    let compressed = dir.path().join("level9.bam");
    let (size_level_9, chunks_9, records_level_9) =
        sort_at_level(sort_order, &input, &compressed, memory_limit, 9, &expected);

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
