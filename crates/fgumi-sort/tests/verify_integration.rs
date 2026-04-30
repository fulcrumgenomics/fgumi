//! Integration tests for `fgumi_sort::verify_sort_order`.
//!
//! Originally in `src/lib/commands/sort.rs` of the main `fgumi` crate;
//! relocated here when `verify_sort_order` was lifted into `fgumi-sort`.

use anyhow::Result;
use fgumi_sam::SamBuilder;
use fgumi_sort::{
    RawBamRecordReader, RawExternalSorter, RawQuerynameKey, RawQuerynameLexKey, RawSortKey,
    SortContext, SortOrder, extract_coordinate_key_inline, verify_sort_order,
};
use fgumi_sort::keys::QuerynameComparator;

// ============================================================================
// Basic verify_sort_order tests
// ============================================================================

#[test]
fn test_verify_sort_order_sorted() -> Result<()> {
    let mut builder = SamBuilder::new();
    let _ = builder.add_pair().name("aaa").build();
    let _ = builder.add_pair().name("bbb").build();
    let _ = builder.add_pair().name("ccc").build();

    let dir = tempfile::tempdir()?;
    let bam_path = dir.path().join("sorted.bam");
    builder.write_bam(&bam_path)?;

    let file = std::fs::File::open(&bam_path)?;
    let mut reader = RawBamRecordReader::new(file)?;
    reader.skip_header()?;

    let (total, violations, first_violation) = verify_sort_order(
        reader,
        |bam| fgumi_raw_bam::RawRecordView::new(bam).read_name().to_vec(),
        |key, prev| key < prev,
    )?;

    assert_eq!(total, 6); // 3 pairs = 6 records
    assert_eq!(violations, 0);
    assert!(first_violation.is_none());
    Ok(())
}

#[test]
fn test_verify_sort_order_unsorted() -> Result<()> {
    let mut builder = SamBuilder::new();
    let _ = builder.add_pair().name("ccc").build();
    let _ = builder.add_pair().name("aaa").build();
    let _ = builder.add_pair().name("bbb").build();

    let dir = tempfile::tempdir()?;
    let bam_path = dir.path().join("unsorted.bam");
    builder.write_bam(&bam_path)?;

    let file = std::fs::File::open(&bam_path)?;
    let mut reader = RawBamRecordReader::new(file)?;
    reader.skip_header()?;

    let (total, violations, first_violation) = verify_sort_order(
        reader,
        |bam| fgumi_raw_bam::RawRecordView::new(bam).read_name().to_vec(),
        |key, prev| key < prev,
    )?;

    assert_eq!(total, 6);
    assert!(violations > 0);
    assert!(first_violation.is_some());
    let (record_num, name) =
        first_violation.expect("first violation should be present for unsorted file");
    assert!(record_num > 1);
    assert!(!name.is_empty());
    Ok(())
}

#[test]
fn test_verify_sort_order_empty() -> Result<()> {
    let builder = SamBuilder::new();

    let dir = tempfile::tempdir()?;
    let bam_path = dir.path().join("empty.bam");
    builder.write_bam(&bam_path)?;

    let file = std::fs::File::open(&bam_path)?;
    let mut reader = RawBamRecordReader::new(file)?;
    reader.skip_header()?;

    let (total, violations, first_violation) = verify_sort_order(
        reader,
        |bam| fgumi_raw_bam::RawRecordView::new(bam).read_name().to_vec(),
        |key, prev| key < prev,
    )?;

    assert_eq!(total, 0);
    assert_eq!(violations, 0);
    assert!(first_violation.is_none());
    Ok(())
}

// ============================================================================
// Sort-and-verify round-trip helpers and tests
// ============================================================================

/// Build a BAM, sort it with the given order, then verify it passes.
fn sort_and_verify_pass(order: SortOrder) -> Result<()> {
    let mut builder = SamBuilder::new();
    let _ = builder.add_pair().name("read2").contig(0).start1(200).build();
    let _ = builder.add_pair().name("read10").contig(0).start1(100).build();
    let _ = builder.add_pair().name("read1").contig(1).start1(50).build();

    let dir = tempfile::tempdir()?;
    let input_bam = dir.path().join("input.bam");
    let sorted_bam = dir.path().join("sorted.bam");
    builder.write_bam(&input_bam)?;

    RawExternalSorter::new(order).threads(1).output_compression(6).sort(&input_bam, &sorted_bam)?;

    let (_, header) = fgumi_bam_io::create_bam_reader(&sorted_bam, 1)?;
    let file = std::fs::File::open(&sorted_bam)?;
    let mut reader = RawBamRecordReader::new(file)?;
    reader.skip_header()?;

    let (total, violations, _) = match order {
        SortOrder::Coordinate => {
            let nref = header.reference_sequences().len() as u32;
            verify_sort_order(
                reader,
                |bam| extract_coordinate_key_inline(bam, nref),
                |key, prev| key < prev,
            )?
        }
        SortOrder::Queryname(QuerynameComparator::Lexicographic) => {
            let ctx = SortContext::from_header(&header);
            verify_sort_order(
                reader,
                |bam| RawQuerynameLexKey::extract(bam, &ctx),
                |key, prev| key < prev,
            )?
        }
        SortOrder::Queryname(QuerynameComparator::Natural) => {
            let ctx = SortContext::from_header(&header);
            verify_sort_order(
                reader,
                |bam| RawQuerynameKey::extract(bam, &ctx),
                |key, prev| key < prev,
            )?
        }
        SortOrder::TemplateCoordinate => return Ok(()),
    };

    assert!(total > 0, "should have records for {:?}", order);
    assert_eq!(violations, 0, "should be sorted for {:?}", order);
    Ok(())
}

#[test]
fn test_verify_coordinate_sorted_pass() -> Result<()> {
    sort_and_verify_pass(SortOrder::Coordinate)
}

#[test]
fn test_verify_queryname_default_sorted_pass() -> Result<()> {
    sort_and_verify_pass(SortOrder::Queryname(QuerynameComparator::Lexicographic))
}

#[test]
fn test_verify_queryname_natural_sorted_pass() -> Result<()> {
    sort_and_verify_pass(SortOrder::Queryname(QuerynameComparator::Natural))
}

#[test]
fn test_verify_queryname_lex_fails_with_natural_verifier() -> Result<()> {
    let mut builder = SamBuilder::new();
    let _ = builder.add_pair().name("read2").contig(0).start1(100).build();
    let _ = builder.add_pair().name("read10").contig(0).start1(200).build();

    let dir = tempfile::tempdir()?;
    let input_bam = dir.path().join("input.bam");
    let sorted_bam = dir.path().join("sorted.bam");
    builder.write_bam(&input_bam)?;

    RawExternalSorter::new(SortOrder::Queryname(QuerynameComparator::Lexicographic))
        .threads(1)
        .output_compression(6)
        .sort(&input_bam, &sorted_bam)?;

    let (_, header) = fgumi_bam_io::create_bam_reader(&sorted_bam, 1)?;
    let file = std::fs::File::open(&sorted_bam)?;
    let mut reader = RawBamRecordReader::new(file)?;
    reader.skip_header()?;

    let ctx = SortContext::from_header(&header);
    let (total, violations, _) = verify_sort_order(
        reader,
        |bam| RawQuerynameKey::extract(bam, &ctx),
        |key, prev| key < prev,
    )?;

    assert!(total > 0);
    assert!(violations > 0, "lex-sorted file should fail natural verify");
    Ok(())
}

#[test]
fn test_verify_queryname_natural_fails_with_lex_verifier() -> Result<()> {
    let mut builder = SamBuilder::new();
    let _ = builder.add_pair().name("read2").contig(0).start1(100).build();
    let _ = builder.add_pair().name("read10").contig(0).start1(200).build();

    let dir = tempfile::tempdir()?;
    let input_bam = dir.path().join("input.bam");
    let sorted_bam = dir.path().join("sorted.bam");
    builder.write_bam(&input_bam)?;

    RawExternalSorter::new(SortOrder::Queryname(QuerynameComparator::Natural))
        .threads(1)
        .output_compression(6)
        .sort(&input_bam, &sorted_bam)?;

    let (_, header) = fgumi_bam_io::create_bam_reader(&sorted_bam, 1)?;
    let file = std::fs::File::open(&sorted_bam)?;
    let mut reader = RawBamRecordReader::new(file)?;
    reader.skip_header()?;

    let ctx = SortContext::from_header(&header);
    let (total, violations, _) = verify_sort_order(
        reader,
        |bam| RawQuerynameLexKey::extract(bam, &ctx),
        |key, prev| key < prev,
    )?;

    assert!(total > 0);
    assert!(violations > 0, "natural-sorted file should fail lex verify");
    Ok(())
}

#[test]
fn test_verify_coordinate_fails_on_unsorted() -> Result<()> {
    let mut builder = SamBuilder::new();
    let _ = builder.add_pair().name("a").contig(1).start1(100).build();
    let _ = builder.add_pair().name("b").contig(0).start1(200).build();

    let dir = tempfile::tempdir()?;
    let bam_path = dir.path().join("unsorted.bam");
    builder.write_bam(&bam_path)?;

    let (_, header) = fgumi_bam_io::create_bam_reader(&bam_path, 1)?;
    let file = std::fs::File::open(&bam_path)?;
    let mut reader = RawBamRecordReader::new(file)?;
    reader.skip_header()?;

    let nref = header.reference_sequences().len() as u32;
    let (total, violations, _) = verify_sort_order(
        reader,
        |bam| extract_coordinate_key_inline(bam, nref),
        |key, prev| key < prev,
    )?;

    assert!(total > 0);
    assert!(violations > 0, "unsorted file should fail coordinate verify");
    Ok(())
}
