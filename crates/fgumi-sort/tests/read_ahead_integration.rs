//! Integration test: exercises the pool-integrated BAM reader path via
//! `RawExternalSorter`, which internally uses `worker_pool` and
//! `read_ahead::PooledInputStream`.
//!
//! Originally, the logic for `create_raw_bam_reader_pool_integrated` and
//! `skip_bam_header` was excised from `src/lib/bam_io.rs` in the main
//! `fgumi` crate (Phase B, Task B.7) and placed in the sort module because
//! it depends on `SortWorkerPool` and `PooledInputStream`. Those helpers now
//! live in `fgumi_sort::external` (previously `sort/raw.rs`). This test
//! exercises the integrated reader path through a complete sort round-trip.

use fgumi_sam::SamBuilder;
use fgumi_sort::{RawExternalSorter, SortOrder};

/// Builds a small BAM file with `n_pairs` paired-end reads, sorts it
/// coordinate-order via the pool-integrated reader path, and verifies all
/// records survive.
#[test]
fn test_pool_integrated_reader_round_trip_coordinate() {
    let n_pairs: usize = 10;
    let mut builder = SamBuilder::new();
    for i in 0..n_pairs {
        let _ = builder
            .add_pair()
            .name(&format!("read{i}"))
            .start1(i * 200 + 1)
            .start2(i * 200 + 101)
            .build();
    }

    let dir = tempfile::tempdir().expect("tempdir");
    let input = dir.path().join("input.bam");
    let output = dir.path().join("output.bam");
    builder.write_bam(&input).expect("write_bam");

    let stats = RawExternalSorter::new(SortOrder::Coordinate)
        .threads(2)
        .memory_limit(64 * 1024 * 1024) // 64 MiB — fits in memory, no spill
        .temp_compression(0)
        .output_compression(0)
        .sort(&input, &output)
        .expect("sort should succeed");

    let expected = (n_pairs * 2) as u64;
    assert_eq!(
        stats.total_records, expected,
        "expected {} records, got {}",
        expected, stats.total_records
    );
    assert_eq!(stats.output_records, expected);
    assert_eq!(stats.chunks_written, 0, "tiny BAM should fit in memory");
}

/// Same round-trip but with template-coordinate sort order to exercise the
/// `TemplateKey` extraction path.
#[test]
fn test_pool_integrated_reader_round_trip_template_coordinate() {
    let n_pairs: usize = 10;
    let mut builder = SamBuilder::new();
    for i in 0..n_pairs {
        let _ = builder
            .add_pair()
            .name(&format!("pair{i}"))
            .start1(i * 300 + 1)
            .start2(i * 300 + 151)
            .build();
    }

    let dir = tempfile::tempdir().expect("tempdir");
    let input = dir.path().join("input.bam");
    let output = dir.path().join("output_tc.bam");
    builder.write_bam(&input).expect("write_bam");

    let stats = RawExternalSorter::new(SortOrder::TemplateCoordinate)
        .threads(2)
        .memory_limit(64 * 1024 * 1024)
        .temp_compression(0)
        .output_compression(0)
        .sort(&input, &output)
        .expect("template-coordinate sort should succeed");

    let expected = (n_pairs * 2) as u64;
    assert_eq!(stats.total_records, expected);
    assert_eq!(stats.output_records, expected);
}
