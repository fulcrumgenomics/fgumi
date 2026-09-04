//! End-to-end tests for the dedup command.
//!
//! These tests invoke `MarkDuplicates::execute()` in-process and validate:
//! 1. Basic duplicate marking
//! 2. Metrics output
//! 3. Remove duplicates mode

use clap::Parser;
use rstest::rstest;

use fgoxide::io::DelimFile;
use fgumi_lib::commands::command::Command;
use fgumi_lib::commands::dedup::MarkDuplicates;
use fgumi_lib::metrics::DeduplicationMetrics;
use fgumi_lib::sam::SamTag;
use fgumi_raw_bam::{RawRecord, SamBuilder, flags};
use noodles::bam;
use noodles::sam::alignment::io::Write as AlignmentWrite;
use noodles::sam::alignment::record::data::field::Tag;
use noodles::sam::alignment::record_buf::data::field::Value;
use std::fs;
use std::path::Path;
use tempfile::TempDir;

use crate::helpers::assertions::assert_bam_sorted;
use crate::helpers::bam_generator::{create_minimal_header, to_record_buf};

/// Create a template-coordinate sorted BAM with UMI-tagged reads.
///
/// Template-coordinate sort groups reads by position, then by name within each position.
/// The header must have SO:unsorted GO:query SS:template-coordinate tags.
fn create_sorted_bam(path: &Path, records: Vec<RawRecord>) {
    // Delegate to the header-aware variant with the default minimal header so the
    // two helpers share a single implementation and can never drift apart.
    let header = create_minimal_header("chr1", 10000);
    create_sorted_bam_with_header(path, &header, records);
}

/// Create a group of paired-end reads at the same position with the same UMI
/// (simulating PCR duplicates).
fn create_duplicate_group(base_name: &str, umi: &str, count: usize, start: i32) -> Vec<RawRecord> {
    create_duplicate_group_inner(base_name, umi, count, start, None, 60, 30)
}

/// A single paired-end template (R1 + R2) at `start` with UMI `umi`, every base
/// at quality `base_qual`. `dedup` scores a template by the sum of its primary
/// reads' base qualities, so building several same-coordinate/same-UMI templates
/// with *distinct* `base_qual`s makes which one wins the representative slot (and
/// which are marked duplicate) deterministic — the unique maximum wins, with no
/// score tie whose resolution would depend on stream order.
fn create_template_qual(name: &str, umi: &str, start: i32, base_qual: u8) -> Vec<RawRecord> {
    create_duplicate_group_inner(name, umi, 1, start, None, 60, base_qual)
}

/// Like [`create_duplicate_group`] but with an explicit `mapq`, so a fixture can
/// carry sub-threshold reads that `--min-map-q` filters (exercising the filter
/// funnel on the chain path).
fn create_duplicate_group_mapq(
    base_name: &str,
    umi: &str,
    count: usize,
    start: i32,
    mapq: u8,
) -> Vec<RawRecord> {
    create_duplicate_group_inner(base_name, umi, count, start, None, mapq, 30)
}

/// Shared implementation for [`create_duplicate_group`] and
/// [`create_duplicate_group_with_rg`]: builds `count` paired-end duplicate
/// templates, tagging each record with `RG:Z:{rg_id}` only when `rg_id` is
/// `Some`. Keeping one implementation means the exact record shape the per-library
/// and ladder tests derive their template counts from can never diverge between
/// the RG and non-RG variants.
fn create_duplicate_group_inner(
    base_name: &str,
    umi: &str,
    count: usize,
    start: i32,
    rg_id: Option<&str>,
    mapq: u8,
    base_qual: u8,
) -> Vec<RawRecord> {
    let mut records = Vec::new();
    for i in 0..count {
        let name = format!("{base_name}_{i}");

        let r1 = {
            let mut b = SamBuilder::new();
            b.read_name(name.as_bytes())
                .sequence(b"ACGTACGT")
                .qualities(&[base_qual; 8])
                .flags(flags::PAIRED | flags::FIRST_SEGMENT)
                .ref_id(0)
                .pos(start - 1)
                .mapq(mapq)
                .cigar_ops(&[8 << 4]) // 8M
                .mate_ref_id(0)
                .mate_pos(start + 99)
                .template_length(108)
                .add_string_tag(SamTag::RX, umi.as_bytes())
                .add_string_tag(SamTag::MC, b"8M");
            if let Some(rg_id) = rg_id {
                b.add_string_tag(SamTag::RG, rg_id.as_bytes());
            }
            b.build()
        };

        let r2 = {
            let mut b = SamBuilder::new();
            b.read_name(name.as_bytes())
                .sequence(b"ACGTACGT")
                .qualities(&[base_qual; 8])
                .flags(flags::PAIRED | flags::LAST_SEGMENT | flags::REVERSE)
                .ref_id(0)
                .pos(start + 99)
                .mapq(mapq)
                .cigar_ops(&[8 << 4]) // 8M
                .mate_ref_id(0)
                .mate_pos(start - 1)
                .template_length(-108)
                .add_string_tag(SamTag::RX, umi.as_bytes())
                .add_string_tag(SamTag::MC, b"8M");
            if let Some(rg_id) = rg_id {
                b.add_string_tag(SamTag::RG, rg_id.as_bytes());
            }
            b.build()
        };

        records.push(r1);
        records.push(r2);
    }
    records
}

/// Test basic dedup command (mark duplicates).
#[test]
fn test_dedup_command_basic() {
    let temp_dir = TempDir::new().unwrap();
    let input_bam = temp_dir.path().join("input.bam");
    let output_bam = temp_dir.path().join("output.bam");

    // Create 3 duplicate pairs at position 100 with same UMI, and 2 at position 500
    let mut records = create_duplicate_group("dup1", "ACGTACGT", 3, 100);
    records.extend(create_duplicate_group("dup2", "TGCATGCA", 2, 500));
    create_sorted_bam(&input_bam, records);

    let cmd = MarkDuplicates::try_parse_from([
        "dedup",
        "--input",
        input_bam.to_str().unwrap(),
        "--output",
        output_bam.to_str().unwrap(),
        "--strategy",
        "identity",
        "--compression-level",
        "1",
    ])
    .expect("failed to parse dedup args");
    cmd.execute("fgumi dedup").expect("Dedup command failed");
    assert!(output_bam.exists(), "Output BAM not created");

    // dedup consumes and re-emits template-coordinate order; gate that its output
    // still verifies as template-coordinate sorted (SS:template-coordinate).
    assert_bam_sorted(&output_bam, "template-coordinate", Some("mi"));

    // All reads should be present (duplicates are marked, not removed)
    let mut reader = bam::io::Reader::new(fs::File::open(&output_bam).unwrap());
    let _header = reader.read_header().unwrap();
    let count = reader.records().count();
    assert_eq!(count, 10, "All reads should be in output (marked, not removed)");
}

/// End-to-end guard for fgumi#739. With a mapping-quality threshold no read can
/// meet, every template is dropped -- and the metrics file must say so, rather
/// than reporting a row of zeros with no indication anything was discarded.
///
/// This is the regression test for the reported bug: before the fix, the counts
/// were collected and merged correctly, then dropped at the serialization
/// boundary because the output struct had no field to hold them.
///
/// The same input is run twice -- once unfiltered, once with the threshold -- and
/// the templates that reach the output in the first run must equal the templates
/// reported as filtered in the second. Asserting the two against each other,
/// rather than against a hard-coded count, keeps the test about the accounting
/// invariant instead of about how this fixture happens to group into templates.
#[test]
fn test_dedup_metrics_report_filtered_templates() {
    fn run_dedup(temp_dir: &TempDir, label: &str, min_map_q: &str) -> DeduplicationMetrics {
        let input_bam = temp_dir.path().join(format!("{label}.in.bam"));
        let output_bam = temp_dir.path().join(format!("{label}.out.bam"));
        let metrics_path = temp_dir.path().join(format!("{label}.metrics.txt"));

        create_sorted_bam(&input_bam, create_duplicate_group("dup1", "ACGTACGT", 3, 100));

        let cmd = MarkDuplicates::try_parse_from([
            "dedup",
            "--input",
            input_bam.to_str().unwrap(),
            "--output",
            output_bam.to_str().unwrap(),
            "--strategy",
            "identity",
            "--min-map-q",
            min_map_q,
            "--metrics",
            metrics_path.to_str().unwrap(),
            "--compression-level",
            "1",
        ])
        .expect("failed to parse dedup args");
        cmd.execute("fgumi dedup").expect("Dedup command failed");

        let rows: Vec<DeduplicationMetrics> =
            DelimFile::default().read_tsv(&metrics_path).expect("failed to read dedup metrics");
        rows.into_iter().next().expect("one metrics row")
    }

    let temp_dir = TempDir::new().unwrap();

    // Baseline: nothing is filtered, so every template reaches the output.
    let kept = run_dedup(&temp_dir, "kept", "0");
    assert_eq!(kept.filtered_templates, 0, "baseline should filter nothing");
    assert!(kept.total_templates > 0, "baseline should emit templates");

    // Same input, mapping-quality threshold no read can meet.
    let dropped = run_dedup(&temp_dir, "dropped", "61");

    assert_eq!(
        dropped.filtered_templates, kept.total_templates,
        "every template the baseline emitted should be reported as filtered"
    );
    assert_eq!(
        dropped.filtered_low_mapping_quality, kept.total_templates,
        "and all of them attributed to low mapping quality"
    );
    assert_eq!(dropped.total_templates, 0, "nothing should reach the output");
    assert_eq!(dropped.filtered_ns_in_umi, 0, "no other reason should be credited");
    assert_eq!(dropped.filtered_missing_umi, 0);
    assert_eq!(dropped.filtered_unmapped, 0);
}

/// Test dedup command with metrics output.
#[test]
fn test_dedup_command_with_metrics() {
    let temp_dir = TempDir::new().unwrap();
    let input_bam = temp_dir.path().join("input.bam");
    let output_bam = temp_dir.path().join("output.bam");
    let metrics_path = temp_dir.path().join("metrics.txt");

    let records = create_duplicate_group("dup1", "ACGTACGT", 3, 100);
    create_sorted_bam(&input_bam, records);

    let cmd = MarkDuplicates::try_parse_from([
        "dedup",
        "--input",
        input_bam.to_str().unwrap(),
        "--output",
        output_bam.to_str().unwrap(),
        "--strategy",
        "identity",
        "--metrics",
        metrics_path.to_str().unwrap(),
        "--compression-level",
        "1",
    ])
    .expect("failed to parse dedup args");
    cmd.execute("fgumi dedup").expect("Dedup command with metrics failed");
    assert!(metrics_path.exists(), "Metrics file not created");
}

/// End-to-end guard for the #804 pair/orphan breakdown. The fixture is entirely
/// mapped paired-end duplicates. `dedup` dedups each mate in its own position
/// group, so the 3 read-pairs become 6 single-mate templates (4 marked duplicate:
/// 2 per group); classified by flags they are 6 pair halves -> `mapped_pairs = 3`,
/// `duplicate_pairs = 2`, with nothing in the orphan or unmapped categories. Also
/// asserts the documented read-unit reconciliation invariants.
#[test]
fn test_dedup_metrics_pair_orphan_breakdown_all_mapped_pairs() {
    let temp_dir = TempDir::new().unwrap();
    let input_bam = temp_dir.path().join("input.bam");
    let output_bam = temp_dir.path().join("output.bam");
    let metrics_path = temp_dir.path().join("metrics.txt");

    // 3 mapped paired-end duplicates: 1 kept + 2 marked duplicate.
    create_sorted_bam(&input_bam, create_duplicate_group("dup1", "ACGTACGT", 3, 100));

    let cmd = MarkDuplicates::try_parse_from([
        "dedup",
        "--input",
        input_bam.to_str().unwrap(),
        "--output",
        output_bam.to_str().unwrap(),
        "--strategy",
        "identity",
        "--metrics",
        metrics_path.to_str().unwrap(),
        "--compression-level",
        "1",
    ])
    .expect("failed to parse dedup args");
    cmd.execute("fgumi dedup").expect("Dedup command failed");

    let rows: Vec<DeduplicationMetrics> =
        DelimFile::default().read_tsv(&metrics_path).expect("failed to read dedup metrics");
    // One per-library row plus the aggregate "All Reads" row, so a regression in
    // `DedupCounts::merge` (which builds the aggregate) cannot hide behind the first row.
    assert_eq!(rows.len(), 2, "one library row plus the All Reads aggregate");

    // The pair/orphan breakdown and both read-unit reconciliation invariants must hold on
    // every row, including the aggregate. This fixture has a single library, so the aggregate
    // repeats the per-library counts — but asserting it exercises the merge path all the same.
    let assert_breakdown = |row: &DeduplicationMetrics, label: &str| {
        // 3 read-pairs -> 6 single-mate templates, 4 marked duplicate (2 per group).
        assert_eq!(row.total_templates, 6, "{label}: 3 pairs -> 6 single-mate templates");
        assert_eq!(row.duplicate_templates, 4, "{label}: duplicate templates");
        // Every read is a pair half, halved to whole pairs.
        assert_eq!(row.mapped_pairs, 3, "{label}: 6 pair halves -> 3 mapped pairs");
        assert_eq!(row.duplicate_pairs, 2, "{label}: 4 duplicate pair halves -> 2 duplicate pairs");
        assert_eq!(row.mapped_orphans, 0, "{label}: mapped orphans");
        assert_eq!(row.duplicate_orphans, 0, "{label}: duplicate orphans");
        assert_eq!(row.unmapped_pairs, 0, "{label}: unmapped pairs");
        assert_eq!(row.unmapped_orphans, 0, "{label}: unmapped orphans");
        assert_eq!(row.unmated_templates, 0, "{label}: no single-end reads in this fixture");

        // Documented read-unit reconciliation invariants (each pair == two templates).
        assert_eq!(
            2 * row.mapped_pairs + row.mapped_orphans + row.unmapped_pairs + row.unmapped_orphans,
            row.total_templates,
            "{label}: mapping buckets must reconcile with total_templates in read units"
        );
        assert_eq!(
            2 * row.duplicate_pairs + row.duplicate_orphans,
            row.duplicate_templates,
            "{label}: duplicate buckets must reconcile with duplicate_templates in read units"
        );
    };

    let all_reads = rows
        .iter()
        .find(|r| r.library == "All Reads")
        .expect("metrics file must contain an All Reads aggregate row");
    let per_library = rows
        .iter()
        .find(|r| r.library != "All Reads")
        .expect("metrics file must contain a per-library row");

    assert_breakdown(per_library, "per-library row");
    assert_breakdown(all_reads, "All Reads row");
}

/// Build a template-coordinate-sorted SAM header with two `@RG` lines that have
/// distinct `LB` values ("libA"/"libB"), for the per-library metrics test below.
fn create_multi_library_header(ref_name: &str, ref_len: usize) -> noodles::sam::Header {
    // Both read groups share one sample, so the resolved `sample` column is a
    // single value; the distinct-sample ordering is exercised separately.
    create_multi_library_header_with_samples(ref_name, ref_len, "sampleX", "sampleX")
}

/// Like [`create_multi_library_header`], but lets each read group declare its
/// own `@RG SM:` value so a test can pin how distinct samples resolve.
fn create_multi_library_header_with_samples(
    ref_name: &str,
    ref_len: usize,
    sample_rg1: &str,
    sample_rg2: &str,
) -> noodles::sam::Header {
    use bstr::BString;
    use noodles::sam::header::record::value::Map;
    use noodles::sam::header::record::value::map::Map as HeaderRecordMap;
    use noodles::sam::header::record::value::map::header::tag::Tag as HeaderTag;
    use noodles::sam::header::record::value::map::read_group::tag as rg_tag;
    use noodles::sam::header::record::value::map::{
        Header as HeaderRecord, ReadGroup, ReferenceSequence,
    };
    use std::num::NonZeroUsize;

    let mut header_builder = HeaderRecordMap::<HeaderRecord>::builder();
    for &(tag_bytes, value) in
        &[(*b"SO", "unsorted"), (*b"GO", "query"), (*b"SS", "template-coordinate")]
    {
        let HeaderTag::Other(tag) = HeaderTag::from(tag_bytes) else { unreachable!() };
        header_builder = header_builder.insert(tag, value);
    }
    let header_map = header_builder.build().expect("valid header map");

    let reference_sequence = Map::<ReferenceSequence>::new(
        NonZeroUsize::new(ref_len).expect("reference length must be non-zero"),
    );

    let rg_a = Map::<ReadGroup>::builder()
        .insert(rg_tag::LIBRARY, String::from("libA"))
        .insert(rg_tag::SAMPLE, String::from(sample_rg1))
        .build()
        .expect("building read group RG1 should succeed");
    let rg_b = Map::<ReadGroup>::builder()
        .insert(rg_tag::LIBRARY, String::from("libB"))
        .insert(rg_tag::SAMPLE, String::from(sample_rg2))
        .build()
        .expect("building read group RG2 should succeed");

    noodles::sam::Header::builder()
        .set_header(header_map)
        .add_reference_sequence(BString::from(ref_name), reference_sequence)
        .add_read_group(BString::from("RG1"), rg_a)
        .add_read_group(BString::from("RG2"), rg_b)
        .build()
}

/// Like [`create_duplicate_group`], but tags every record with `RG:Z:{rg_id}`
/// so it is attributed to a specific library at dedup time.
fn create_duplicate_group_with_rg(
    base_name: &str,
    umi: &str,
    count: usize,
    start: i32,
    rg_id: &str,
) -> Vec<RawRecord> {
    create_duplicate_group_inner(base_name, umi, count, start, Some(rg_id), 60, 30)
}

/// Shared implementation for [`create_sorted_bam`]: writes `records` against the
/// caller-supplied `header`, then template-coordinate sorts the result into
/// `path` with fgumi's own sorter so the dedup input is *genuinely* in
/// template-coordinate order rather than merely labelled as such. (dedup trusts
/// the header's SS tag; feeding it mislabelled-but-unsorted input produces
/// mislabelled-but-unsorted output.) [`create_sorted_bam`] passes the default
/// [`create_minimal_header`]; the per-library tests pass a multi-`@RG` header so
/// each read is attributed to a specific library at dedup time.
fn create_sorted_bam_with_header(
    path: &Path,
    header: &noodles::sam::Header,
    records: Vec<RawRecord>,
) {
    let unsorted = path.with_file_name("dedup_input_unsorted.bam");
    let mut writer =
        bam::io::Writer::new(fs::File::create(&unsorted).expect("Failed to create BAM file"));
    writer.write_header(header).expect("Failed to write header");
    for record in records {
        writer
            .write_alignment_record(header, &to_record_buf(&record))
            .expect("Failed to write record");
    }
    writer.try_finish().expect("Failed to finish BAM");

    let status = std::process::Command::new(env!("CARGO_BIN_EXE_fgumi"))
        .args([
            "sort",
            "-i",
            unsorted.to_str().unwrap(),
            "-o",
            path.to_str().unwrap(),
            "--order",
            "template-coordinate",
            "--key-types",
            "mi",
        ])
        .status()
        .expect("failed to spawn `fgumi sort` for dedup test input");
    assert!(status.success(), "failed to template-coordinate sort dedup test input");
}

/// Reads the dedup metrics TSV back as one `column name -> value` map per data
/// row, in file order.
fn read_dedup_metrics_rows(path: &Path) -> Vec<std::collections::HashMap<String, String>> {
    let content = fs::read_to_string(path).expect("failed to read metrics file");
    let mut lines = content.lines();
    let header: Vec<String> =
        lines.next().expect("metrics header row").split('\t').map(String::from).collect();
    lines
        .map(|line| {
            let values: Vec<&str> = line.split('\t').collect();
            header.iter().cloned().zip(values.iter().map(ToString::to_string)).collect()
        })
        .collect()
}

/// Validate one library's `--duplication-ladder` rows and return their count.
///
/// Checks: `templates_seen` matches the derived interval-crossing sequence
/// exactly (which also pins the snapshot-row count) and is strictly increasing;
/// cumulative and per-window `duplicate_fraction` both in `[0, 1]`; each
/// `window_templates` equals the increment over the previous snapshot and the
/// windows sum to the true per-library total; the final snapshot lands at that
/// total with a cumulative fraction matching the metrics file's
/// `duplicate_rate`.
/// Assert `row[field]` parses to the exact fraction `num/den` (within a small
/// float tolerance), so a pinned expected table catches wrong saturation values,
/// not merely out-of-range ones.
fn assert_fraction_eq(
    row: &std::collections::HashMap<String, String>,
    field: &str,
    (num, den): (u64, u64),
    label: &str,
) {
    let actual: f64 = row[field].parse().unwrap_or_else(|_| panic!("{field} is f64"));
    #[expect(
        clippy::cast_precision_loss,
        reason = "numerator and denominator are small test constants (< 100)"
    )]
    let expected = num as f64 / den as f64;
    assert!(
        (actual - expected).abs() < 1e-9,
        "{label} {field} ({actual}) must equal {num}/{den} ({expected}): {row:?}"
    );
}

#[expect(clippy::type_complexity, reason = "an inline pinned table of expected rows")]
fn assert_ladder_library_ok(
    ladder_rows: &[std::collections::HashMap<String, String>],
    library: &str,
    expected_rows: &[(u64, (u64, u64), (u64, u64))],
    library_metrics: &std::collections::HashMap<String, String>,
) -> usize {
    let rows: Vec<_> = ladder_rows.iter().filter(|r| r["library"] == library).collect();

    // Pin every row's snapshot exactly, in file order (not merely by count), so a
    // regression that emits the right number of rows at the wrong interval
    // crossings — or with the wrong cumulative/marginal saturation values — is
    // caught.
    assert_eq!(
        rows.len(),
        expected_rows.len(),
        "{library}: expected {} ladder rows: {rows:?}",
        expected_rows.len()
    );

    let mut previous_seen = 0u64;
    let mut window_sum = 0u64;
    for (row, &(expected_seen, cumulative, window)) in rows.iter().zip(expected_rows) {
        let templates_seen: u64 = row["templates_seen"].parse().expect("templates_seen is u64");
        assert_eq!(
            templates_seen, expected_seen,
            "{library}: snapshot must land on the derived interval crossing: {row:?}"
        );
        assert!(
            templates_seen > previous_seen,
            "{library}: templates_seen must be strictly increasing: {rows:?}"
        );

        // Cumulative duplicate_fraction, pinned exactly (not merely in [0, 1]).
        assert_fraction_eq(
            row,
            "duplicate_fraction",
            cumulative,
            &format!("{library}: cumulative"),
        );

        // Marginal (window) columns: this window's depth is the increment over
        // the previous snapshot, and its fraction is pinned exactly.
        let window_templates: u64 =
            row["window_templates"].parse().expect("window_templates is u64");
        assert_eq!(
            window_templates,
            templates_seen - previous_seen,
            "{library}: window_templates must equal the increment since the last snapshot: {row:?}"
        );
        assert_fraction_eq(row, "window_duplicate_fraction", window, &format!("{library}: window"));
        window_sum += window_templates;
        previous_seen = templates_seen;
    }

    let last_row = rows.last().unwrap_or_else(|| panic!("{library}: no ladder rows emitted"));
    assert_eq!(
        last_row["templates_seen"], library_metrics["total_templates"],
        "{library}: final ladder snapshot must land at the true per-library total: {last_row:?}"
    );
    // The per-window depths partition the cumulative total exactly.
    assert_eq!(
        window_sum.to_string(),
        library_metrics["total_templates"],
        "{library}: window_templates must sum to the per-library total"
    );
    let expected_fraction: f64 =
        library_metrics["duplicate_rate"].parse().expect("duplicate_rate is f64");
    let actual_fraction: f64 =
        last_row["duplicate_fraction"].parse().expect("duplicate_fraction is f64");
    assert!(
        (actual_fraction - expected_fraction).abs() < 1e-9,
        "{library}: final duplicate_fraction ({actual_fraction}) must match the metrics \
         file's cumulative duplicate_rate ({expected_fraction})"
    );
    rows.len()
}

/// Dedup over a two-library input (library A: 3 duplicate pairs at one
/// position; library B: 2 duplicate pairs at a different position) and verify
/// the metrics file has one row per library, sorted by name, plus a final
/// "All Reads" total row summing across libraries — with correct per-library
/// `duplicate_templates`/`percent_duplication`, and `estimated_library_size`
/// populated for both (both libraries have duplicates).
///
/// Each mate's position group is dedup'd independently (R1's group and R2's
/// group are template-coordinate-adjacent but distinct `RawPositionGroup`s,
/// each holding single-mate templates), so a library's `count` duplicate
/// pairs contribute `2 * count` templates/reads to its row: `count` from R1's
/// group (1 unique + `count - 1` duplicate) and `count` from R2's group
/// (same split) — e.g. library A's 3 pairs give `total_templates = 6`,
/// `unique_templates = 2`, `duplicate_templates = 4`.
#[test]
fn test_dedup_command_per_library_metrics() {
    let temp_dir = TempDir::new().unwrap();
    let input_bam = temp_dir.path().join("input.bam");
    let output_bam = temp_dir.path().join("output.bam");
    let metrics_path = temp_dir.path().join("metrics.txt");

    let header = create_multi_library_header("chr1", 10000);
    let mut records = create_duplicate_group_with_rg("libA_dup", "ACGTACGT", 3, 100, "RG1");
    records.extend(create_duplicate_group_with_rg("libB_dup", "TGCATGCA", 2, 500, "RG2"));
    create_sorted_bam_with_header(&input_bam, &header, records);

    let cmd = MarkDuplicates::try_parse_from([
        "dedup",
        "--input",
        input_bam.to_str().unwrap(),
        "--output",
        output_bam.to_str().unwrap(),
        "--strategy",
        "identity",
        "--metrics",
        metrics_path.to_str().unwrap(),
        "--compression-level",
        "1",
    ])
    .expect("failed to parse dedup args");
    cmd.execute("fgumi dedup").expect("Dedup command with per-library metrics failed");

    let rows = read_dedup_metrics_rows(&metrics_path);
    assert_eq!(rows.len(), 3, "expected libA row + libB row + All Reads total row: {rows:?}");

    let libraries: Vec<&str> = rows.iter().map(|r| r["library"].as_str()).collect();
    assert_eq!(
        libraries,
        vec!["libA", "libB", "All Reads"],
        "rows must be library-name-sorted, with the total row last: {libraries:?}"
    );

    // The `sample` column is resolved from the header's @RG SM tags (both read
    // groups declare SM:sampleX) and is the same on every row.
    for row in &rows {
        assert_eq!(row["sample"], "sampleX", "sample resolved from @RG SM: {row:?}");
    }

    // Each position group in `fgumi dedup` holds a single mate (R1's group and
    // R2's group are template-coordinate-adjacent but distinct groups), so a
    // library's `count` duplicate pairs contribute `2 * count` templates/reads
    // (one per mate's own position group), each with exactly one record.
    let lib_a = &rows[0];
    assert_eq!(lib_a["total_templates"], "6");
    assert_eq!(lib_a["unique_templates"], "2");
    assert_eq!(lib_a["duplicate_templates"], "4");
    assert_eq!(lib_a["total_reads"], "6");
    assert_eq!(lib_a["duplicate_reads"], "4");
    assert!(
        (lib_a["percent_duplication"].parse::<f64>().unwrap() - 4.0 / 6.0).abs() < 1e-6,
        "lib_a percent_duplication (read-level 4/6): {lib_a:?}"
    );
    assert!(
        !lib_a["estimated_library_size"].is_empty(),
        "lib_a has duplicates, so the library-size estimate must be populated: {lib_a:?}"
    );

    let lib_b = &rows[1];
    assert_eq!(lib_b["total_templates"], "4");
    assert_eq!(lib_b["unique_templates"], "2");
    assert_eq!(lib_b["duplicate_templates"], "2");
    assert_eq!(lib_b["total_reads"], "4");
    assert_eq!(lib_b["duplicate_reads"], "2");
    assert!(
        (lib_b["percent_duplication"].parse::<f64>().unwrap() - 2.0 / 4.0).abs() < 1e-6,
        "lib_b percent_duplication (read-level 2/4): {lib_b:?}"
    );
    assert!(
        !lib_b["estimated_library_size"].is_empty(),
        "lib_b has duplicates, so the library-size estimate must be populated: {lib_b:?}"
    );

    let total = &rows[2];
    assert_eq!(total["total_templates"], "10");
    assert_eq!(total["unique_templates"], "4");
    assert_eq!(total["duplicate_templates"], "6");
    assert_eq!(total["total_reads"], "10");
    assert_eq!(total["duplicate_reads"], "6");
    // The aggregate row spans two distinct libraries, so a pooled Lander-Waterman
    // estimate would be meaningless: it must be left empty (the per-library rows
    // above carry the meaningful estimates). Matches dupblaster.
    assert!(
        total["estimated_library_size"].is_empty(),
        "the aggregate row spans >1 library, so its estimate must be empty: {total:?}"
    );
}

/// Distinct `@RG SM:` values resolve deterministically: the `sample` column is
/// the sorted, comma-joined set of unique samples, identical on every row,
/// regardless of the order the read groups declare them in the header. The
/// header here lists `sampleB` before `sampleA` so a non-sorting resolver would
/// emit `sampleB,sampleA`.
#[test]
fn test_dedup_command_metrics_sample_column_sorts_distinct_samples() {
    let temp_dir = TempDir::new().unwrap();
    let input_bam = temp_dir.path().join("input.bam");
    let output_bam = temp_dir.path().join("output.bam");
    let metrics_path = temp_dir.path().join("metrics.txt");

    // RG1 declares the lexically-later sample, RG2 the earlier one, so the
    // header order is the reverse of the expected sorted output.
    let header = create_multi_library_header_with_samples("chr1", 10000, "sampleB", "sampleA");
    let mut records = create_duplicate_group_with_rg("libA_dup", "ACGTACGT", 3, 100, "RG1");
    records.extend(create_duplicate_group_with_rg("libB_dup", "TGCATGCA", 2, 500, "RG2"));
    create_sorted_bam_with_header(&input_bam, &header, records);

    let cmd = MarkDuplicates::try_parse_from([
        "dedup",
        "--input",
        input_bam.to_str().unwrap(),
        "--output",
        output_bam.to_str().unwrap(),
        "--strategy",
        "identity",
        "--metrics",
        metrics_path.to_str().unwrap(),
        "--compression-level",
        "1",
    ])
    .expect("failed to parse dedup args");
    cmd.execute("fgumi dedup").expect("Dedup command with distinct samples failed");

    let rows = read_dedup_metrics_rows(&metrics_path);
    assert_eq!(rows.len(), 3, "expected libA row + libB row + All Reads total row: {rows:?}");

    let libraries: Vec<&str> = rows.iter().map(|r| r["library"].as_str()).collect();
    assert_eq!(
        libraries,
        vec!["libA", "libB", "All Reads"],
        "rows must be library-name-sorted, with the total row last: {libraries:?}"
    );

    for row in &rows {
        assert_eq!(
            row["sample"], "sampleA,sampleB",
            "the sample column must be the sorted, comma-joined unique @RG SM values, \
             independent of header order: {row:?}"
        );
    }
}

/// `--sample` overrides the header-derived `@RG SM:` value on every metrics row.
#[test]
fn test_dedup_command_sample_override_wins_over_header() {
    let temp_dir = TempDir::new().unwrap();
    let input_bam = temp_dir.path().join("input.bam");
    let output_bam = temp_dir.path().join("output.bam");
    let metrics_path = temp_dir.path().join("metrics.txt");

    // The header declares @RG SM:sampleX; --sample must win.
    let header = create_multi_library_header("chr1", 10000);
    let records = create_duplicate_group_with_rg("libA_dup", "ACGTACGT", 3, 100, "RG1");
    create_sorted_bam_with_header(&input_bam, &header, records);

    let cmd = MarkDuplicates::try_parse_from([
        "dedup",
        "--input",
        input_bam.to_str().unwrap(),
        "--output",
        output_bam.to_str().unwrap(),
        "--strategy",
        "identity",
        "--metrics",
        metrics_path.to_str().unwrap(),
        "--sample",
        "OVERRIDE",
        "--compression-level",
        "1",
    ])
    .expect("failed to parse dedup args");
    cmd.execute("fgumi dedup").expect("Dedup command with --sample override failed");

    let rows = read_dedup_metrics_rows(&metrics_path);
    assert_eq!(rows.len(), 2, "expected libA row + All Reads total row: {rows:?}");

    let libraries: Vec<&str> = rows.iter().map(|r| r["library"].as_str()).collect();
    assert_eq!(
        libraries,
        vec!["libA", "All Reads"],
        "rows must be library-name-sorted, with the total row last: {libraries:?}"
    );

    for row in &rows {
        assert_eq!(row["sample"], "OVERRIDE", "--sample must override @RG SM: {row:?}");
    }
}

/// Builds a two-`@RG` header whose read groups declare `LB` but no `SM:`, for
/// the empty-`sample` resolution test below.
fn create_multi_library_header_without_samples(
    ref_name: &str,
    ref_len: usize,
) -> noodles::sam::Header {
    use bstr::BString;
    use noodles::sam::header::record::value::Map;
    use noodles::sam::header::record::value::map::Map as HeaderRecordMap;
    use noodles::sam::header::record::value::map::header::tag::Tag as HeaderTag;
    use noodles::sam::header::record::value::map::read_group::tag as rg_tag;
    use noodles::sam::header::record::value::map::{
        Header as HeaderRecord, ReadGroup, ReferenceSequence,
    };
    use std::num::NonZeroUsize;

    let mut header_builder = HeaderRecordMap::<HeaderRecord>::builder();
    for &(tag_bytes, value) in
        &[(*b"SO", "unsorted"), (*b"GO", "query"), (*b"SS", "template-coordinate")]
    {
        let HeaderTag::Other(tag) = HeaderTag::from(tag_bytes) else { unreachable!() };
        header_builder = header_builder.insert(tag, value);
    }
    let header_map = header_builder.build().expect("valid header map");

    let reference_sequence = Map::<ReferenceSequence>::new(
        NonZeroUsize::new(ref_len).expect("reference length must be non-zero"),
    );

    // Read groups carry LB but deliberately omit SM.
    let rg_a = Map::<ReadGroup>::builder()
        .insert(rg_tag::LIBRARY, String::from("libA"))
        .build()
        .expect("building read group RG1 should succeed");
    let rg_b = Map::<ReadGroup>::builder()
        .insert(rg_tag::LIBRARY, String::from("libB"))
        .build()
        .expect("building read group RG2 should succeed");

    noodles::sam::Header::builder()
        .set_header(header_map)
        .add_reference_sequence(BString::from(ref_name), reference_sequence)
        .add_read_group(BString::from("RG1"), rg_a)
        .add_read_group(BString::from("RG2"), rg_b)
        .build()
}

/// When neither `--sample` nor any `@RG SM:` value is present, the `sample`
/// column resolves to the empty string on every row.
#[test]
fn test_dedup_command_sample_empty_when_no_rg_sm() {
    let temp_dir = TempDir::new().unwrap();
    let input_bam = temp_dir.path().join("input.bam");
    let output_bam = temp_dir.path().join("output.bam");
    let metrics_path = temp_dir.path().join("metrics.txt");

    // Read groups declare LB but no @RG SM:, and --sample is not passed.
    let header = create_multi_library_header_without_samples("chr1", 10000);
    let mut records = create_duplicate_group_with_rg("libA_dup", "ACGTACGT", 3, 100, "RG1");
    records.extend(create_duplicate_group_with_rg("libB_dup", "TGCATGCA", 2, 500, "RG2"));
    create_sorted_bam_with_header(&input_bam, &header, records);

    let cmd = MarkDuplicates::try_parse_from([
        "dedup",
        "--input",
        input_bam.to_str().unwrap(),
        "--output",
        output_bam.to_str().unwrap(),
        "--strategy",
        "identity",
        "--metrics",
        metrics_path.to_str().unwrap(),
        "--compression-level",
        "1",
    ])
    .expect("failed to parse dedup args");
    cmd.execute("fgumi dedup").expect("Dedup command with no @RG SM: failed");

    let rows = read_dedup_metrics_rows(&metrics_path);
    assert_eq!(rows.len(), 3, "expected libA row + libB row + All Reads total row: {rows:?}");

    let libraries: Vec<&str> = rows.iter().map(|r| r["library"].as_str()).collect();
    assert_eq!(
        libraries,
        vec!["libA", "libB", "All Reads"],
        "rows must be library-name-sorted, with the total row last: {libraries:?}"
    );

    for row in &rows {
        assert_eq!(
            row["sample"], "",
            "sample must be empty when neither --sample nor any @RG SM: is present: {row:?}"
        );
    }
}

/// Row order with named libraries AND reads lacking `@RG`/`LB`: named
/// libraries first (alphabetical), then the "Unknown Library" catch-all, then
/// the "All Reads" total last. The catch-all is a sentinel, not a real
/// library, so it is grouped after the named libraries rather than sorted in
/// among them by name.
#[test]
fn test_dedup_command_metrics_orders_unknown_after_named_before_total() {
    let temp_dir = TempDir::new().unwrap();
    let input_bam = temp_dir.path().join("input.bam");
    let output_bam = temp_dir.path().join("output.bam");
    let metrics_path = temp_dir.path().join("metrics.txt");

    let header = create_multi_library_header("chr1", 10000);
    // Named libraries (RG1 -> libA, RG2 -> libB) plus reads with no @RG, which
    // fall into the idx-0 "Unknown Library" bucket. Distinct positions keep them
    // in separate position groups.
    let mut records = create_duplicate_group_with_rg("libA_dup", "ACGTACGT", 3, 100, "RG1");
    records.extend(create_duplicate_group_with_rg("libB_dup", "TGCATGCA", 2, 500, "RG2"));
    records.extend(create_duplicate_group("noRG_dup", "GGCCGGCC", 2, 900));
    create_sorted_bam_with_header(&input_bam, &header, records);

    let cmd = MarkDuplicates::try_parse_from([
        "dedup",
        "--input",
        input_bam.to_str().unwrap(),
        "--output",
        output_bam.to_str().unwrap(),
        "--strategy",
        "identity",
        "--metrics",
        metrics_path.to_str().unwrap(),
        "--compression-level",
        "1",
    ])
    .expect("failed to parse dedup args");
    cmd.execute("fgumi dedup").expect("Dedup command with mixed-library metrics failed");

    let rows = read_dedup_metrics_rows(&metrics_path);
    let libraries: Vec<&str> = rows.iter().map(|r| r["library"].as_str()).collect();
    assert_eq!(
        libraries,
        vec!["libA", "libB", "Unknown Library", "All Reads"],
        "named libraries first (alphabetical), then Unknown Library, then the total: {libraries:?}"
    );
}

/// Single-library input (no `@RG`/`LB` at all) must still produce two rows:
/// the "Unknown Library" bucket, then the "All Reads" total — the metrics
/// writer never collapses to a single row, regardless of library count, so
/// the output shape is uniform across single- and multi-library inputs.
#[test]
fn test_dedup_command_single_library_metrics_has_unknown_and_total_rows() {
    let temp_dir = TempDir::new().unwrap();
    let input_bam = temp_dir.path().join("input.bam");
    let output_bam = temp_dir.path().join("output.bam");
    let metrics_path = temp_dir.path().join("metrics.txt");

    let records = create_duplicate_group("dup1", "ACGTACGT", 3, 100);
    create_sorted_bam(&input_bam, records);

    let cmd = MarkDuplicates::try_parse_from([
        "dedup",
        "--input",
        input_bam.to_str().unwrap(),
        "--output",
        output_bam.to_str().unwrap(),
        "--strategy",
        "identity",
        "--metrics",
        metrics_path.to_str().unwrap(),
        "--compression-level",
        "1",
    ])
    .expect("failed to parse dedup args");
    cmd.execute("fgumi dedup").expect("Dedup command with metrics failed");

    let rows = read_dedup_metrics_rows(&metrics_path);
    let libraries: Vec<&str> = rows.iter().map(|r| r["library"].as_str()).collect();
    assert_eq!(
        libraries,
        vec!["Unknown Library", "All Reads"],
        "no @RG/LB in the input: expect the Unknown Library bucket, then the total: {libraries:?}"
    );
    assert_eq!(rows[0]["total_templates"], rows[1]["total_templates"]);
}

/// Dedup over a two-library input with `--duplication-ladder` and a small
/// `--ladder-interval`, verifying per-library snapshot rows are correct.
///
/// libA gets three duplicate-group positions with `count` = 3, 2, 2 pairs
/// (spaced far enough apart to fall in distinct position groups). Per
/// [`test_dedup_command_per_library_metrics`]'s documented shape, each
/// position's mate groups (R1 then R2, in coordinate order) each contribute
/// `count` templates, so libA's cumulative per-group `templates_seen`
/// sequence is 3, 6, 8, 10, 12, 14 (true total 14). With `--ladder-interval
/// 4`, thresholds are crossed at cumulative values 6, 8, 12 — three interval
/// rows — and the true total (14) does not land on a crossing, so a fourth,
/// final row is appended at 14: four rows total.
///
/// libB gets two duplicate-group positions with `count` = 2 each, giving the
/// cumulative sequence 2, 4, 6, 8 (true total 8). Thresholds are crossed at 4
/// and 8 — two interval rows — and the second crossing lands exactly on the
/// true total, so no extra final row is appended: two rows total.
#[test]
fn test_dedup_duplication_ladder_multi_library() {
    let temp_dir = TempDir::new().unwrap();
    let input_bam = temp_dir.path().join("input.bam");
    let output_bam = temp_dir.path().join("output.bam");
    let metrics_path = temp_dir.path().join("metrics.txt");
    let ladder_path = temp_dir.path().join("ladder.txt");

    let header = create_multi_library_header("chr1", 10000);
    let mut records = create_duplicate_group_with_rg("libA_p1", "AAAAAAAA", 3, 100, "RG1");
    records.extend(create_duplicate_group_with_rg("libA_p2", "CCCCCCCC", 2, 500, "RG1"));
    records.extend(create_duplicate_group_with_rg("libA_p3", "GGGGGGGG", 2, 900, "RG1"));
    records.extend(create_duplicate_group_with_rg("libB_p1", "TTTTTTTT", 2, 1300, "RG2"));
    records.extend(create_duplicate_group_with_rg("libB_p2", "ACACACAC", 2, 1700, "RG2"));
    create_sorted_bam_with_header(&input_bam, &header, records);

    let cmd = MarkDuplicates::try_parse_from([
        "dedup",
        "--input",
        input_bam.to_str().unwrap(),
        "--output",
        output_bam.to_str().unwrap(),
        "--strategy",
        "identity",
        "--metrics",
        metrics_path.to_str().unwrap(),
        "--duplication-ladder",
        ladder_path.to_str().unwrap(),
        "--ladder-interval",
        "4",
        "--compression-level",
        "1",
    ])
    .expect("failed to parse dedup args");
    cmd.execute("fgumi dedup").expect("Dedup command with duplication ladder failed");

    assert!(ladder_path.exists(), "duplication ladder file must be created when requested");

    let metrics_rows = read_dedup_metrics_rows(&metrics_path);
    let metrics_by_library: std::collections::HashMap<
        &str,
        &std::collections::HashMap<String, String>,
    > = metrics_rows.iter().map(|r| (r["library"].as_str(), r)).collect();

    let ladder_rows = read_dedup_metrics_rows(&ladder_path);

    // Expected snapshot sequence per library, pinned exactly (not merely by row
    // count) so a regression that emitted the right number of rows at the wrong
    // interval crossings is still caught: libA crosses the interval at 6, 8, 12,
    // then gets a final row at its true total (14), which doesn't land on a
    // crossing — 4 rows. libB crosses at 4 and 8, with 8 landing exactly on its
    // true total, so no extra final row — 2 rows.
    // Every snapshot row is pinned exactly — both `templates_seen` and the
    // cumulative `duplicate_fraction` (as `(numerator, denominator)`) — so a
    // regression that emits the right row count at the wrong interval crossings,
    // reorders the library rows, or writes the wrong intermediate saturation
    // values is caught, not just one with a wrong final row. libA crosses the
    // interval at 6, 8, 12, then gets a final row at its true total (14), which
    // doesn't land on a crossing — 4 rows. libB crosses at 4 and 8, with 8
    // landing exactly on its true total, so no extra final row — 2 rows.
    // Per snapshot row, keyed by library: (templates_seen, cumulative (num,
    // den), window (num, den)), pinned so the table reads as a spec. The window
    // duplicate count is the cumulative duplicates minus the previous snapshot's,
    // over the window's templates: libA cumulative dups run 4, 5, 7, 8 over seen
    // 6, 8, 12, 14, so the windows are 4/6, 1/2, 2/4, 1/2; libB runs 2, 4 over 4,
    // 8, so 2/4, 2/4.
    #[expect(clippy::type_complexity, reason = "an inline pinned table of expected rows")]
    let expected_sequences: [(&str, &[(u64, (u64, u64), (u64, u64))]); 2] = [
        (
            "libA",
            &[
                (6, (4, 6), (4, 6)),
                (8, (5, 8), (1, 2)),
                (12, (7, 12), (2, 4)),
                (14, (8, 14), (1, 2)),
            ],
        ),
        ("libB", &[(4, (2, 4), (2, 4)), (8, (4, 8), (2, 4))]),
    ];
    let mut matched_row_count = 0usize;

    // For each library: every row's templates_seen, cumulative duplicate_fraction,
    // window_templates and window_duplicate_fraction match the pinned sequence in
    // file order; the windows sum to the true per-library total; and the last
    // row's cumulative fraction matches the metrics file's template-level
    // duplicate_rate exactly.
    for (library, expected_rows) in expected_sequences {
        matched_row_count += assert_ladder_library_ok(
            &ladder_rows,
            library,
            expected_rows,
            metrics_by_library[library],
        );
    }

    assert_eq!(
        matched_row_count,
        ladder_rows.len(),
        "the ladder is inherently per-library: no unexpected library names: {ladder_rows:?}"
    );
}

/// Without `--duplication-ladder`, no ladder file is written and dedup's
/// ordinary behavior is unaffected (mirrors [`test_dedup_command_basic`]).
#[test]
fn test_dedup_duplication_ladder_off_by_default() {
    let temp_dir = TempDir::new().unwrap();
    let input_bam = temp_dir.path().join("input.bam");
    let output_bam = temp_dir.path().join("output.bam");
    let ladder_path = temp_dir.path().join("ladder_never_written.txt");

    let mut records = create_duplicate_group("dup1", "ACGTACGT", 3, 100);
    records.extend(create_duplicate_group("dup2", "TGCATGCA", 2, 500));
    create_sorted_bam(&input_bam, records);

    let cmd = MarkDuplicates::try_parse_from([
        "dedup",
        "--input",
        input_bam.to_str().unwrap(),
        "--output",
        output_bam.to_str().unwrap(),
        "--strategy",
        "identity",
        "--compression-level",
        "1",
    ])
    .expect("failed to parse dedup args");
    cmd.execute("fgumi dedup").expect("Dedup command without duplication ladder failed");

    assert!(output_bam.exists(), "Output BAM not created");
    assert!(
        !ladder_path.exists(),
        "duplication ladder file must not be created when --duplication-ladder is not passed"
    );

    // Record-identity oracle: dedup marks duplicates in place, so the output must
    // hold exactly the input's 10 records (both segments of all five templates),
    // none dropped, added, or renamed — and carrying the duplicate marks the
    // identity strategy assigns. A bare count of 10 would still pass if a record
    // were dropped and another duplicated, or the wrong templates were marked.
    let mut reader = bam::io::Reader::new(fs::File::open(&output_bam).unwrap());
    let header = reader.read_header().unwrap();
    let mut names: Vec<String> = Vec::new();
    // Per template name, the duplicate flag of each of its segments.
    let mut duplicate_flags_by_template: std::collections::BTreeMap<String, Vec<bool>> =
        std::collections::BTreeMap::new();
    for record in reader.record_bufs(&header) {
        let record = record.expect("read output record");
        let name =
            String::from_utf8(record.name().expect("output record must have a name").to_vec())
                .expect("read name is UTF-8");
        duplicate_flags_by_template
            .entry(name.clone())
            .or_default()
            .push(record.flags().is_duplicate());
        names.push(name);
    }

    names.sort_unstable();
    assert_eq!(
        names,
        vec![
            "dup1_0", "dup1_0", "dup1_1", "dup1_1", "dup1_2", "dup1_2", "dup2_0", "dup2_0",
            "dup2_1", "dup2_1",
        ],
        "output must hold exactly the input records (both segments of every template): {names:?}"
    );

    // Both segments of a template must carry the same duplicate mark, and the
    // identity strategy keeps exactly one primary per position group: dup1's three
    // templates yield two duplicates, dup2's two yield one.
    for (name, flags) in &duplicate_flags_by_template {
        assert_eq!(flags.len(), 2, "template {name} must have both segments in output");
        assert_eq!(
            flags[0], flags[1],
            "template {name}: both segments must share the duplicate mark"
        );
    }
    let marked_duplicate_templates = |prefix: &str| {
        duplicate_flags_by_template
            .iter()
            .filter(|(name, _)| name.starts_with(prefix))
            .filter(|(_, flags)| flags[0])
            .count()
    };
    assert_eq!(
        marked_duplicate_templates("dup1"),
        2,
        "dup1: 3 templates at one position -> 1 primary kept, 2 marked duplicate"
    );
    assert_eq!(
        marked_duplicate_templates("dup2"),
        1,
        "dup2: 2 templates at one position -> 1 primary kept, 1 marked duplicate"
    );
}

/// Regression test: a single position group whose template increment is
/// large enough to leap past more than one `--ladder-interval` multiple in a
/// single `DuplicationLadderRecorder::record` call must still emit only one
/// row per call, not one row per crossed multiple.
///
/// A single 10-pair duplicate group (no `@RG`, so library = "Unknown
/// Library") with `--ladder-interval 3` produces two position groups (R1
/// then R2, in coordinate order), each contributing 10 templates in one
/// `record()` call — each call alone crosses three interval multiples (the
/// first call crosses 3, 6, 9; the second crosses 12, 15, 18). Before the fix
/// this produced duplicate rows `[10, 10, 10, 20, 20, 20]`; after the fix it
/// must produce exactly `[10, 20]`.
#[test]
fn test_dedup_duplication_ladder_single_group_exceeds_interval() {
    let temp_dir = TempDir::new().unwrap();
    let input_bam = temp_dir.path().join("input.bam");
    let output_bam = temp_dir.path().join("output.bam");
    let metrics_path = temp_dir.path().join("metrics.txt");
    let ladder_path = temp_dir.path().join("ladder.txt");

    let records = create_duplicate_group("dup1", "ACGTACGT", 10, 100);
    create_sorted_bam(&input_bam, records);

    let cmd = MarkDuplicates::try_parse_from([
        "dedup",
        "--input",
        input_bam.to_str().unwrap(),
        "--output",
        output_bam.to_str().unwrap(),
        "--strategy",
        "identity",
        "--metrics",
        metrics_path.to_str().unwrap(),
        "--duplication-ladder",
        ladder_path.to_str().unwrap(),
        "--ladder-interval",
        "3",
        "--compression-level",
        "1",
    ])
    .expect("failed to parse dedup args");
    cmd.execute("fgumi dedup").expect("Dedup command with duplication ladder failed");

    let metrics_rows = read_dedup_metrics_rows(&metrics_path);
    let unknown_library = metrics_rows
        .iter()
        .find(|r| r["library"] == "Unknown Library")
        .expect("Unknown Library row present");

    let ladder_rows = read_dedup_metrics_rows(&ladder_path);
    let templates_seen: Vec<u64> = ladder_rows
        .iter()
        .map(|r| r["templates_seen"].parse().expect("templates_seen is u64"))
        .collect();

    assert_eq!(
        templates_seen,
        vec![10, 20],
        "a single group whose increment leaps past more than one interval \
         multiple must still emit exactly one row per record() call, not one \
         row per crossed multiple (regression: previously produced \
         [10, 10, 10, 20, 20, 20]): {ladder_rows:?}"
    );

    // Redundant with the exact-sequence check above, but states the general
    // invariant explicitly: no duplicate/non-increasing templates_seen.
    for window in templates_seen.windows(2) {
        assert!(
            window[1] > window[0],
            "templates_seen must be strictly increasing: {templates_seen:?}"
        );
    }

    // Pin the marginal (window) columns on the multi-threshold path too, not just
    // the cumulative fraction: a defect in `window_templates` or
    // `window_duplicate_fraction` — e.g. one that mis-attributes depth when a
    // single `record()` call crosses several interval multiples — can pass a
    // cumulative-only check. Each position group contributes exactly 10 templates
    // with 9 duplicates, so every window is 10 templates at 9/10.
    for row in &ladder_rows {
        assert_eq!(
            row["window_templates"], "10",
            "each window covers exactly one 10-template position group: {row:?}"
        );
        assert_fraction_eq(row, "window_duplicate_fraction", (9, 10), "multi-threshold window");
    }

    // Pin the intermediate row's fraction too, not just the final one: the
    // first snapshot lands at 10 templates_seen with 9 cumulative duplicates.
    let first_row = ladder_rows.first().expect("at least one ladder row");
    let first_fraction: f64 =
        first_row["duplicate_fraction"].parse().expect("duplicate_fraction is f64");
    assert!(
        (first_fraction - 9.0 / 10.0).abs() < 1e-9,
        "the first snapshot's duplicate_fraction must be 9/10: {first_row:?}"
    );

    let last_row = ladder_rows.last().expect("at least one ladder row");
    assert_eq!(
        last_row["templates_seen"], unknown_library["total_templates"],
        "final ladder snapshot must land at the true total: {last_row:?}"
    );
    let expected_fraction: f64 =
        unknown_library["duplicate_rate"].parse().expect("duplicate_rate is f64");
    let actual_fraction: f64 =
        last_row["duplicate_fraction"].parse().expect("duplicate_fraction is f64");
    assert!(
        (actual_fraction - expected_fraction).abs() < 1e-9,
        "final duplicate_fraction ({actual_fraction}) must match the metrics \
         file's cumulative duplicate_rate ({expected_fraction})"
    );
}

/// Test dedup command with remove-duplicates flag.
#[test]
fn test_dedup_command_remove_duplicates() {
    let temp_dir = TempDir::new().unwrap();
    let input_bam = temp_dir.path().join("input.bam");
    let output_bam = temp_dir.path().join("output.bam");

    // 3 duplicate pairs → should keep 1 pair (2 records), remove 2 pairs (4 records)
    let records = create_duplicate_group("dup1", "ACGTACGT", 3, 100);
    create_sorted_bam(&input_bam, records);

    let cmd = MarkDuplicates::try_parse_from([
        "dedup",
        "--input",
        input_bam.to_str().unwrap(),
        "--output",
        output_bam.to_str().unwrap(),
        "--strategy",
        "identity",
        "--remove-duplicates",
        "--compression-level",
        "1",
    ])
    .expect("failed to parse dedup args");
    cmd.execute("fgumi dedup").expect("Dedup command with --remove-duplicates failed");
    assert!(output_bam.exists(), "Output BAM not created");

    // With remove-duplicates, only the best pair should remain
    let mut reader = bam::io::Reader::new(fs::File::open(&output_bam).unwrap());
    let _header = reader.read_header().unwrap();
    let count = reader.records().count();
    assert!(count < 6, "Remove-duplicates should produce fewer reads than input");
    assert!(count >= 2, "Should keep at least one pair");
}

/// Create a paired-end template whose *both* mates are unmapped (`ref_id=-1`, `pos=-1`),
/// simulating a read pair that failed to align.
fn create_unmapped_pair(name: &str, umi: &str) -> Vec<RawRecord> {
    let make = |segment_flag: u16| {
        let mut b = SamBuilder::new();
        b.read_name(name.as_bytes())
            .sequence(b"ACGTACGT")
            .qualities(&[30; 8])
            .flags(flags::PAIRED | flags::UNMAPPED | flags::MATE_UNMAPPED | segment_flag)
            .add_string_tag(SamTag::RX, umi.as_bytes());
        b.build()
    };
    vec![make(flags::FIRST_SEGMENT), make(flags::LAST_SEGMENT)]
}

/// `--include-unmapped` emits templates with no mapped read untouched, while the default
/// drops them. Input is 3 mapped duplicate pairs (6 reads) + 1 fully-unmapped pair (2 reads);
/// the default keeps only the 6 mapped reads, and `--include-unmapped` keeps all 8.
/// `expected_passthrough` / `expected_filtered_unmapped` pin where the no-mapped-read
/// template is *accounted*, not just whether its reads survive: the two modes must route
/// it to exactly one of the two columns, never both and never neither.
#[rstest]
#[case::default_drops_unmapped(false, 6, 0, 0, 1, 0)]
#[case::include_unmapped_keeps(true, 8, 2, 1, 0, 2)]
fn test_dedup_include_unmapped(
    #[case] include_unmapped: bool,
    #[case] expected_total: usize,
    #[case] expected_unmapped: usize,
    #[case] expected_passthrough: u64,
    #[case] expected_filtered_unmapped: u64,
    #[case] expected_unmapped_pairs: u64,
) {
    use noodles::sam::alignment::record::data::field::Tag;
    use noodles::sam::alignment::record_buf::RecordBuf;
    use noodles::sam::alignment::record_buf::data::field::value::Value as DataValue;

    let temp_dir = TempDir::new().unwrap();
    let input_bam = temp_dir.path().join("input.bam");
    let output_bam = temp_dir.path().join("output.bam");
    let metrics_path = temp_dir.path().join("metrics.txt");

    let mut records = create_duplicate_group("dup1", "ACGTACGT", 3, 100);
    records.extend(create_unmapped_pair("unmapped1", "TGCATGCA"));
    create_sorted_bam(&input_bam, records);

    let mut args: Vec<String> = vec![
        "dedup".into(),
        "--input".into(),
        input_bam.to_str().unwrap().into(),
        "--output".into(),
        output_bam.to_str().unwrap().into(),
        "--strategy".into(),
        "identity".into(),
        "--compression-level".into(),
        "1".into(),
        "--metrics".into(),
        metrics_path.to_str().unwrap().into(),
    ];
    if include_unmapped {
        args.push("--include-unmapped".into());
    }
    let cmd = MarkDuplicates::try_parse_from(&args).expect("failed to parse dedup args");
    cmd.execute("fgumi dedup").expect("Dedup command failed");

    let mut reader = bam::io::Reader::new(fs::File::open(&output_bam).unwrap());
    let out_header = reader.read_header().unwrap();

    let mi_tag = Tag::from(SamTag::MI);
    let mut total = 0usize;
    let mut unmapped_seen = 0usize;
    for result in reader.record_bufs(&out_header) {
        let record: RecordBuf = result.expect("read record");
        total += 1;
        let record_flags = record.flags();
        if record_flags.is_unmapped() {
            unmapped_seen += 1;
            // Pass-through reads are emitted untouched: never duplicate-marked, never MI-tagged.
            assert!(
                !record_flags.is_duplicate(),
                "unmapped pass-through read must not be duplicate-marked"
            );
            assert!(
                !matches!(record.data().get(&mi_tag), Some(DataValue::String(_))),
                "unmapped pass-through read must not carry an MI tag"
            );
        }
    }

    assert_eq!(total, expected_total, "unexpected output record count");
    assert_eq!(unmapped_seen, expected_unmapped, "unexpected unmapped-read count");

    // The no-mapped-read template must be accounted for in exactly one column: it is
    // either passed through untouched or filtered as unmapped, never both, never neither.
    let rows: Vec<DeduplicationMetrics> =
        DelimFile::default().read_tsv(&metrics_path).expect("failed to read dedup metrics");
    let metrics = rows.first().expect("one metrics row");
    assert_eq!(
        metrics.passthrough_templates, expected_passthrough,
        "unexpected passthrough_templates"
    );
    assert_eq!(
        metrics.filtered_unmapped, expected_filtered_unmapped,
        "unexpected filtered_unmapped"
    );

    // The fully-unmapped pass-through pair is the multi-primary case the
    // pair/orphan reconciliation invariant explicitly excludes: it is emitted as
    // ONE template holding both mates, so it adds 1 to `total_templates` but 2 to
    // `unmapped_pairs` (one per primary read). Under `--include-unmapped` that is
    // two unmapped pairs from a single template; by default the pair is filtered
    // out before pair/orphan counting, so the counter stays zero.
    assert_eq!(metrics.unmapped_pairs, expected_unmapped_pairs, "unexpected unmapped_pairs");

    // The single-primary-read reconciliation `2*mapped_pairs + mapped_orphans +
    // unmapped_pairs + unmapped_orphans == total_templates` holds only for
    // deduplicated (single-primary) templates. This fixture's sole pass-through is
    // a fully-unmapped *pair* (two primary reads in one template), so the left
    // side runs ahead of `total_templates` by exactly `passthrough_templates`.
    let reconciliation_lhs = 2 * metrics.mapped_pairs
        + metrics.mapped_orphans
        + metrics.unmapped_pairs
        + metrics.unmapped_orphans;
    assert_eq!(
        reconciliation_lhs,
        metrics.total_templates + metrics.passthrough_templates,
        "pair/orphan reconciliation must exceed total_templates by the fully-unmapped \
         pass-through pairs: {metrics:?}"
    );
}

/// Create a secondary alignment for `name` with no `tc` tag.
///
/// `ref_id`/`pos` of `None` leave it unplaced, so it travels in the same
/// template-coordinate group as an unmapped primary pair of the same name.
fn create_secondary_without_tc(name: &str, umi: &str, placed: bool) -> RawRecord {
    let mut b = SamBuilder::new();
    b.read_name(name.as_bytes())
        .sequence(b"ACGTACGT")
        .qualities(&[30; 8])
        .flags(flags::PAIRED | flags::SECONDARY | flags::FIRST_SEGMENT)
        .add_string_tag(SamTag::RX, umi.as_bytes());
    if placed {
        // Same coordinate as `create_duplicate_group`'s R1, so it travels in that
        // template's position group.
        b.ref_id(0).pos(99).mapq(60).cigar_ops(&[fgumi_raw_bam::testutil::encode_op(0, 8)]);
    } else {
        b.flags(flags::PAIRED | flags::SECONDARY | flags::FIRST_SEGMENT | flags::UNMAPPED);
    }
    b.build()
}

/// A secondary/supplementary read with no `tc` tag is a hard failure, and
/// `--include-unmapped` must not be a way around it.
///
/// `template_is_unmapped_passthrough` decides on the primary `r1`/`r2` alone, so a
/// template whose primaries are both unmapped is routed to the pass-through partition
/// carrying whatever non-primary records came with it. Those records are still counted
/// as `secondary_reads`/`supplementary_reads`, so if the pass-through loop skipped the
/// `tc` check, `--include-unmapped` would report the read and waive the failure that
/// the same read triggers on any other template.
#[rstest]
#[case::passthrough_template("unmapped1", false)]
#[case::filtered_template("dup1_0", true)]
fn test_dedup_missing_tc_fails_on_every_template(#[case] name: &str, #[case] placed: bool) {
    let temp_dir = TempDir::new().unwrap();
    let input_bam = temp_dir.path().join("input.bam");
    let output_bam = temp_dir.path().join("output.bam");

    let mut records = create_duplicate_group("dup1", "ACGTACGT", 3, 100);
    records.extend(create_unmapped_pair("unmapped1", "TGCATGCA"));
    records.push(create_secondary_without_tc(name, "TGCATGCA", placed));
    create_sorted_bam(&input_bam, records);

    let cmd = MarkDuplicates::try_parse_from([
        "dedup",
        "--input",
        input_bam.to_str().unwrap(),
        "--output",
        output_bam.to_str().unwrap(),
        "--strategy",
        "identity",
        "--include-unmapped",
        "--compression-level",
        "1",
    ])
    .expect("failed to parse dedup args");

    let err = cmd.execute("fgumi dedup").expect_err("a tc-less secondary read must fail dedup");
    assert!(
        err.to_string().contains("missing the `tc` tag"),
        "expected the tc-tag failure, got: {err}"
    );
}

/// Regression test for OOM with large position groups in `--no-umi` mode.
///
/// WES data can have extreme depth pileups at capture targets, creating positions
/// with thousands of reads. With `--no-umi` (identity strategy), ALL reads at the
/// same position form ONE group. This test exercises a 5,000-template position
/// group to verify the pipeline completes without unbounded memory growth.
#[test]
#[allow(clippy::too_many_lines)]
fn test_dedup_no_umi_large_position_group() {
    use bstr::BString;
    use fgumi_raw_bam::raw_record_to_record_buf;
    use fgumi_raw_bam::{SamBuilder as RawSamBuilder, flags, testutil::encode_op};
    use noodles::sam::Header;
    use noodles::sam::alignment::record::data::field::Tag;
    use noodles::sam::alignment::record_buf::RecordBuf;
    use noodles::sam::alignment::record_buf::data::field::value::Value as DataValue;
    use noodles::sam::header::record::value::map::Map as HeaderRecordMap;
    use noodles::sam::header::record::value::map::header::tag::Tag as HeaderTag;
    use noodles::sam::header::record::value::{
        Map, map::Header as HeaderRecord, map::ReferenceSequence,
    };
    use std::collections::HashSet;
    use std::num::NonZeroUsize;

    const NUM_TEMPLATES: usize = 5_000;
    // PROPERLY_PAIRED flag = 0x2 (not exposed as a named constant in flags module).
    const PROPERLY_PAIRED: u16 = 0x2;
    // MATE_REVERSE flag = 0x20 (not exposed as a named constant in flags module).
    // R1 is forward / R2 is reverse; we set MATE_REVERSE on R1 so the pair is FR
    // self-consistent (mate reverse flag matches mate's REVERSE).
    const MATE_REVERSE: u16 = 0x20;

    let temp_dir = TempDir::new().unwrap();
    let input_bam = temp_dir.path().join("input.bam");
    let output_bam = temp_dir.path().join("output.bam");

    // Build a template-coordinate sorted header with one reference.
    let header = {
        let mut builder = HeaderRecordMap::<HeaderRecord>::builder();
        for (tag_bytes, value) in
            [(*b"SO", "unsorted"), (*b"GO", "query"), (*b"SS", "template-coordinate")]
        {
            let HeaderTag::Other(tag) = HeaderTag::from(tag_bytes) else { unreachable!() };
            builder = builder.insert(tag, value);
        }
        Header::builder()
            .set_header(builder.build().expect("valid header map"))
            .add_reference_sequence(
                BString::from("chr1"),
                Map::<ReferenceSequence>::new(NonZeroUsize::new(10_000).expect("non-zero length")),
            )
            .build()
    };

    // Write 5,000 paired-end templates all at position 100 with no UMI tag.
    // Template-coordinate sort groups by position then name, so we write
    // R1 and R2 together for each template in name-sorted order.
    {
        let mut writer = bam::io::Writer::new(fs::File::create(&input_bam).expect("create BAM"));
        writer.write_header(&header).expect("write header");

        for i in 0..NUM_TEMPLATES {
            let name = format!("read_{i:05}");

            let r1 = {
                let mut b = RawSamBuilder::new();
                b.read_name(name.as_bytes())
                    .sequence(b"ACGTACGT")
                    .qualities(&[30; 8])
                    .flags(flags::PAIRED | PROPERLY_PAIRED | flags::FIRST_SEGMENT | MATE_REVERSE)
                    .ref_id(0)
                    .pos(99)
                    .mapq(60)
                    .cigar_ops(&[encode_op(0, 8)])
                    .mate_ref_id(0)
                    .mate_pos(199)
                    .template_length(108);
                b.add_string_tag(SamTag::MC, b"8M");
                b.build()
            };

            let r2 = {
                let mut b = RawSamBuilder::new();
                b.read_name(name.as_bytes())
                    .sequence(b"ACGTACGT")
                    .qualities(&[30; 8])
                    .flags(flags::PAIRED | PROPERLY_PAIRED | flags::REVERSE | flags::LAST_SEGMENT)
                    .ref_id(0)
                    .pos(199)
                    .mapq(60)
                    .cigar_ops(&[encode_op(0, 8)])
                    .mate_ref_id(0)
                    .mate_pos(99)
                    .template_length(-108);
                b.add_string_tag(SamTag::MC, b"8M");
                b.build()
            };

            let r1_buf = raw_record_to_record_buf(&r1, &header).expect("convert R1");
            let r2_buf = raw_record_to_record_buf(&r2, &header).expect("convert R2");
            writer.write_alignment_record(&header, &r1_buf).expect("write R1");
            writer.write_alignment_record(&header, &r2_buf).expect("write R2");
        }
        writer.try_finish().expect("finish BAM");
    }

    // Run dedup with --no-umi via MarkDuplicates::execute().
    let cmd = MarkDuplicates::try_parse_from([
        "dedup",
        "--input",
        input_bam.to_str().unwrap(),
        "--output",
        output_bam.to_str().unwrap(),
        "--no-umi",
        "--compression-level",
        "1",
    ])
    .expect("failed to parse dedup args");
    cmd.execute("fgumi dedup").expect("dedup --no-umi should succeed with large position group");

    assert!(output_bam.exists(), "output BAM should be created");

    // Read back the output and verify duplicate marking and MI tags.
    let mut reader = bam::io::Reader::new(fs::File::open(&output_bam).expect("open output BAM"));
    let out_header = reader.read_header().expect("read output header");

    let mut total_records = 0usize;
    let mut duplicate_records = 0usize;
    let mut non_duplicate_records = 0usize;
    let mut mi_values = HashSet::new();
    let mut mi_count = 0usize;
    let mut non_dup_names = HashSet::new();

    let mi_tag = Tag::from(SamTag::MI);

    for result in reader.record_bufs(&out_header) {
        let record: RecordBuf = result.expect("read record");
        total_records += 1;

        let is_dup = record.flags().is_duplicate();
        if is_dup {
            duplicate_records += 1;
        } else {
            non_duplicate_records += 1;
            if let Some(name) = record.name() {
                non_dup_names.insert(name.to_owned());
            }
        }

        if let Some(DataValue::String(mi)) = record.data().get(&mi_tag) {
            mi_values.insert(mi.to_owned());
            mi_count += 1;
        }
    }

    // All 10,000 records (5,000 pairs) should be present.
    assert_eq!(total_records, NUM_TEMPLATES * 2, "all records should be present in output");

    // Exactly one template (2 records) should NOT be marked as duplicate.
    assert_eq!(
        non_duplicate_records, 2,
        "exactly one pair should be non-duplicate (the best-scoring template)"
    );
    assert_eq!(
        duplicate_records,
        (NUM_TEMPLATES - 1) * 2,
        "all other pairs should be marked as duplicates"
    );

    // The two non-duplicate records should share the same read name.
    assert_eq!(non_dup_names.len(), 1, "non-duplicate records should be from one template");

    // Every record should have an MI tag.
    assert_eq!(mi_count, NUM_TEMPLATES * 2, "all records should have an MI tag");

    // All templates should share the same MI value (one group in identity strategy).
    assert_eq!(mi_values.len(), 1, "all records should share a single MI tag value");
}

/// `dedup` carries its own copy of `--index-threshold`, so it needs the same
/// unsatisfiable-request check `group` has: `always` under a strategy with no index
/// code path must be rejected, not silently ignored.
#[test]
fn test_dedup_rejects_index_threshold_always_when_index_unreachable() {
    let temp_dir = TempDir::new().unwrap();
    let input_bam = temp_dir.path().join("input.bam");
    let output_bam = temp_dir.path().join("output.bam");
    create_sorted_bam(&input_bam, create_duplicate_group("dup1", "ACGTACGT", 3, 100));

    let cmd = MarkDuplicates::try_parse_from([
        "dedup",
        "--input",
        input_bam.to_str().unwrap(),
        "--output",
        output_bam.to_str().unwrap(),
        "--strategy",
        "identity",
        "--index-threshold",
        "always",
    ])
    .expect("failed to parse dedup args");

    let err = cmd.execute("fgumi dedup").expect_err("must reject an unsatisfiable request");
    let message = err.to_string();
    assert!(
        message.contains("--index-threshold always")
            && message.contains("never uses the UMI index"),
        "error should say why the request cannot be honoured, got: {message}"
    );
}

/// Read back (name, flags, MI) per record, in file order.
///
/// `never` only ever changes *how* neighbours are discovered, so this is the granularity
/// the contract is stated at: which molecule each read landed in, and whether it was
/// marked a duplicate.
fn dedup_output_summary(path: &Path) -> Vec<(String, u16, Option<String>)> {
    let mi_tag = Tag::from(SamTag::MI);
    let mut reader = bam::io::Reader::new(fs::File::open(path).expect("failed to open output BAM"));
    let header = reader.read_header().expect("failed to read output header");
    reader
        .record_bufs(&header)
        .map(|record| {
            let record = record.expect("failed to read output record");
            let mi = record.data().get(&mi_tag).map(|value| match value {
                Value::String(mi) => mi.to_string(),
                other => panic!("MI must be a string tag, got {other:?}"),
            });
            let name = record.name().expect("output record must be named").to_string();
            (name, record.flags().bits(), mi)
        })
        .collect()
}

/// Run dedup on `input_bam`, appending `extra_args`, and summarise the output.
fn run_dedup(
    input_bam: &Path,
    output_bam: &Path,
    extra_args: &[&str],
) -> Vec<(String, u16, Option<String>)> {
    let mut args = vec![
        "dedup",
        "--input",
        input_bam.to_str().unwrap(),
        "--output",
        output_bam.to_str().unwrap(),
        "--strategy",
        "identity",
        "--compression-level",
        "1",
    ];
    args.extend_from_slice(extra_args);

    MarkDuplicates::try_parse_from(args)
        .expect("failed to parse dedup args")
        .execute("fgumi dedup")
        .expect("dedup must succeed");
    dedup_output_summary(output_bam)
}

/// `never` is always satisfiable, so it must be accepted everywhere -- including under
/// a strategy that would not have indexed anyway.
///
/// The index is a pure optimisation: it changes which UMI pairs are *compared*, never
/// which reads group together. So the contract is output identity with the default
/// threshold, not merely "a file appeared" -- and the absolute grouping is pinned too,
/// since comparing the two runs alone would pass if both regressed the same way.
#[test]
fn test_dedup_accepts_index_threshold_never() {
    let temp_dir = TempDir::new().unwrap();
    let input_bam = temp_dir.path().join("input.bam");
    let never_bam = temp_dir.path().join("never.bam");
    let default_bam = temp_dir.path().join("default.bam");
    // Three pairs sharing one UMI at one position: one molecule, six records.
    create_sorted_bam(&input_bam, create_duplicate_group("dup1", "ACGTACGT", 3, 100));

    let never = run_dedup(&input_bam, &never_bam, &["--index-threshold", "never"]);
    let default = run_dedup(&input_bam, &default_bam, &[]);

    assert_eq!(
        never, default,
        "--index-threshold never must not change the output; it only suppresses the index"
    );

    // Absolute expectations, so this still fails if both runs regress together.
    assert_eq!(never.len(), 6, "all six records must be emitted (marked, not removed)");

    // All six records, pinned by name, flags and molecule. The two ends sit in
    // different template-coordinate position groups (99 and 199), so each end gets its
    // own molecule id -- the three copies of each end share one, which is the grouping
    // the index would have had to get wrong to matter here.
    let expected = vec![
        ("dup1_0".to_string(), flags::PAIRED | flags::FIRST_SEGMENT, Some("0".to_string())),
        (
            "dup1_1".to_string(),
            flags::PAIRED | flags::FIRST_SEGMENT | flags::DUPLICATE,
            Some("0".to_string()),
        ),
        (
            "dup1_2".to_string(),
            flags::PAIRED | flags::FIRST_SEGMENT | flags::DUPLICATE,
            Some("0".to_string()),
        ),
        (
            "dup1_0".to_string(),
            flags::PAIRED | flags::LAST_SEGMENT | flags::REVERSE,
            Some("1".to_string()),
        ),
        (
            "dup1_1".to_string(),
            flags::PAIRED | flags::LAST_SEGMENT | flags::REVERSE | flags::DUPLICATE,
            Some("1".to_string()),
        ),
        (
            "dup1_2".to_string(),
            flags::PAIRED | flags::LAST_SEGMENT | flags::REVERSE | flags::DUPLICATE,
            Some("1".to_string()),
        ),
    ];
    assert_eq!(never, expected, "--index-threshold never changed the grouping");

    // One template is the primary; the other two pairs are marked duplicates.
    let duplicates = never.iter().filter(|(_, flags, _)| flags & flags::DUPLICATE != 0).count();
    assert_eq!(duplicates, 4, "two of the three pairs must be marked duplicate: {never:?}");
}

/// Corrupt the CRC32 footer of the LAST non-EOF BGZF block in `path`, in place.
///
/// BAM writers flush the compressor after the header, so with enough alignment
/// records the file spans multiple BGZF blocks and the last one holds only
/// record data. Corrupting the *first* block instead would risk landing inside
/// the header's block: the single-threaded raw reader (`FgumiBgzfReader`) applies
/// `verify_crc` uniformly to every block, so a corrupted block 0 fails while the
/// reader is being constructed -- before the intended read/count assertion can
/// run, which would defeat the point of this test.
fn corrupt_last_block_crc(path: &Path) {
    let mut bytes = fs::read(path).expect("read bam for corruption");
    let mut cursor: &[u8] = &bytes;
    let blocks = fgumi_lib::bgzf_reader::read_raw_blocks(&mut cursor, 10_000)
        .expect("read bgzf blocks from test bam");
    assert!(
        blocks.len() >= 2,
        "test input must span at least 2 BGZF blocks so the corrupted block isn't also the \
         header's block; got {} -- generate more records",
        blocks.len()
    );
    let offset: usize =
        blocks[..blocks.len() - 1].iter().map(fgumi_lib::bgzf_reader::RawBgzfBlock::len).sum();
    let last = blocks.last().expect("checked len >= 2 above");
    // `read_raw_blocks` drops every BGZF EOF marker, so summing the returned
    // (real) block lengths yields the last block's on-disk offset only when no
    // marker sits *between* real blocks. Guard that: everything past the last
    // framed block must be whole trailing EOF markers (a writer may emit more
    // than one). An intermediate marker would leave real data here instead and
    // shift `crc_off` onto an unrelated byte, which the `>= 2` guard above cannot
    // detect.
    let eof = &fgumi_lib::bgzf_reader::BGZF_EOF;
    let tail = &bytes[offset + last.len()..];
    assert!(
        tail.len().is_multiple_of(eof.len())
            && tail.chunks_exact(eof.len()).all(|chunk| chunk == &eof[..]),
        "bytes after the last framed block must be only trailing BGZF EOF markers; \
         an intermediate marker would invalidate the CRC offset"
    );
    let crc_off = offset + last.len() - fgumi_lib::bgzf_reader::BGZF_FOOTER_SIZE;
    bytes[crc_off] ^= 0x01;
    fs::write(path, bytes).expect("write corrupted bam");
}

/// A file input with a corrupted BGZF CRC32 must fail dedup by default:
/// `--check-crc`/`--no-check-crc` are unset, and the dupblaster policy
/// verifies for file input. Also confirms the failure is the CRC32 check
/// (not, say, a coincidental header-parse failure from the corruption).
#[test]
fn test_dedup_rejects_corrupted_crc_on_file_input_by_default() {
    let temp_dir = TempDir::new().unwrap();
    let input_bam = temp_dir.path().join("input.bam");
    let output_bam = temp_dir.path().join("output.bam");

    // Enough duplicate pairs to span multiple BGZF blocks (see
    // `corrupt_last_block_crc`).
    create_sorted_bam(&input_bam, create_duplicate_group("dup", "ACGTACGT", 400, 100));
    corrupt_last_block_crc(&input_bam);

    let cmd = MarkDuplicates::try_parse_from([
        "dedup",
        "--input",
        input_bam.to_str().unwrap(),
        "--output",
        output_bam.to_str().unwrap(),
        "--strategy",
        "identity",
    ])
    .expect("failed to parse dedup args");

    let err = cmd
        .execute("fgumi dedup")
        .expect_err("default (verify-on for file input) must reject a corrupted BGZF CRC32");
    let message = format!("{err:#}");
    assert!(message.to_uppercase().contains("CRC32"), "error should mention CRC32: {message}");
}

/// `--no-check-crc` on a file input must accept the same corrupted BGZF CRC32
/// that the default (verify-on) run above rejects.
#[test]
fn test_dedup_no_check_crc_accepts_corrupted_crc_on_file_input() {
    let temp_dir = TempDir::new().unwrap();
    let input_bam = temp_dir.path().join("input.bam");
    let output_bam = temp_dir.path().join("output.bam");
    let clean_bam = temp_dir.path().join("clean.bam");
    let clean_output = temp_dir.path().join("clean_output.bam");

    create_sorted_bam(&input_bam, create_duplicate_group("dup", "ACGTACGT", 400, 100));
    // Keep a pristine copy before corrupting the CRC32, so we can dedup both and
    // compare the decoded output record-for-record.
    fs::copy(&input_bam, &clean_bam).expect("copy pristine input");
    corrupt_last_block_crc(&input_bam);

    let cmd = MarkDuplicates::try_parse_from([
        "dedup",
        "--input",
        input_bam.to_str().unwrap(),
        "--output",
        output_bam.to_str().unwrap(),
        "--strategy",
        "identity",
        "--no-check-crc",
    ])
    .expect("failed to parse dedup args");

    cmd.execute("fgumi dedup")
        .expect("--no-check-crc must accept a corrupted BGZF CRC32 and still complete");
    assert!(output_bam.exists(), "output BAM not created");

    // The CRC32 footer is metadata: skipping its check must not alter a single
    // decoded byte, so the corrupted run must reproduce the clean run exactly.
    // A record-identity oracle catches dropped, duplicated, or reordered records
    // that a bare count of 800 would miss.
    MarkDuplicates::try_parse_from([
        "dedup",
        "--input",
        clean_bam.to_str().unwrap(),
        "--output",
        clean_output.to_str().unwrap(),
        "--strategy",
        "identity",
    ])
    .expect("failed to parse dedup args")
    .execute("fgumi dedup")
    .expect("clean input must dedup");

    let read_records = |path: &Path| {
        let mut reader = bam::io::Reader::new(fs::File::open(path).unwrap());
        let header = reader.read_header().unwrap();
        reader
            .record_bufs(&header)
            .collect::<std::io::Result<Vec<_>>>()
            .expect("read deduped records")
    };
    let actual = read_records(&output_bam);
    let expected = read_records(&clean_output);
    assert_eq!(actual.len(), 800, "all records should be in output (marked, not removed)");
    assert_eq!(actual, expected, "--no-check-crc changed the decoded records");
}

//////////////////////////////////////////////////////////////////////////////
// Chain-path (`--threads`) parity tests
//
// `dedup --threads N` routes through the declarative chain builder; the
// no-`--threads` path keeps the hand-rolled unified pipeline and is the
// in-process oracle these tests diff against.
//////////////////////////////////////////////////////////////////////////////

/// Read a BAM's records back as decoded `RecordBuf`s, for record-for-record
/// comparison (order-sensitive, catches drops/dupes/reorders that a bare count
/// would miss).
fn read_deduped_records(path: &Path) -> Vec<noodles::sam::alignment::RecordBuf> {
    let mut reader = bam::io::Reader::new(fs::File::open(path).unwrap());
    let header = reader.read_header().unwrap();
    reader.record_bufs(&header).collect::<std::io::Result<Vec<_>>>().expect("read deduped records")
}

/// Run `dedup` on `input` writing `output`, with `extra` args appended
/// (e.g. `--threads`, `--strategy`, output flags). Asserts the run succeeds.
fn dedup_run(input: &Path, output: &Path, extra: &[&str]) {
    let mut args =
        vec!["dedup", "--input", input.to_str().unwrap(), "--output", output.to_str().unwrap()];
    args.extend_from_slice(extra);
    MarkDuplicates::try_parse_from(args)
        .expect("failed to parse dedup args")
        .execute("fgumi dedup")
        .expect("dedup run failed");
}

/// The chain (`--threads N`) path produces output records record-for-record
/// identical to the non-chain (no-`--threads`) path. Run at both `--threads 1`
/// (the minimal chain engine) and `--threads 4` (genuinely parallel) — dedup's
/// output is deterministic, so both must equal the single oracle.
#[rstest]
#[case::threads_1(&["--threads", "1"])]
#[case::threads_4(&["--threads", "4"])]
fn test_dedup_chain_matches_single_threaded(#[case] thread_args: &[&str]) {
    let temp_dir = TempDir::new().unwrap();
    let input_bam = temp_dir.path().join("input.bam");

    // Several distinct position groups so the chain sees multiple batches.
    let mut records = Vec::new();
    for i in 0..16 {
        records.extend(create_duplicate_group(&format!("g{i}"), "ACGTACGT", 3, 100 + i * 200));
    }
    create_sorted_bam(&input_bam, records);

    let oracle_out = temp_dir.path().join("oracle.bam");
    dedup_run(&input_bam, &oracle_out, &["--strategy", "identity"]);

    let chain_out = temp_dir.path().join("chain.bam");
    let mut chain_args = vec!["--strategy", "identity"];
    chain_args.extend_from_slice(thread_args);
    dedup_run(&input_bam, &chain_out, &chain_args);

    let expected = read_deduped_records(&oracle_out);
    let actual = read_deduped_records(&chain_out);
    assert!(!expected.is_empty(), "oracle output must be non-empty (guard against a vacuous pass)");
    assert_eq!(
        actual, expected,
        "chain {thread_args:?} output must match the non-chain path record-for-record"
    );
}

/// Build one position group at `start` mixing near-duplicate and distant UMIs.
/// The `ACGTACGT`/`ACGTACGA` pair is edit-distance 1, so under `--strategy
/// adjacency` (the CLI default) they cluster into one molecule while `TTTTTTTT`
/// stays separate — non-trivial clustering that `--strategy identity` (which
/// keeps all three distinct) cannot exercise. This is what makes the default-
/// strategy parity case below meaningful rather than a vacuous pass.
fn create_chain_parity_group(base: &str, start: i32) -> Vec<RawRecord> {
    let mut records = create_duplicate_group(&format!("{base}a"), "ACGTACGT", 3, start);
    records.extend(create_duplicate_group(&format!("{base}b"), "ACGTACGA", 2, start));
    records.extend(create_duplicate_group(&format!("{base}c"), "TTTTTTTT", 2, start));
    // A mapq-10 subfamily: passes at the default `--min-map-q 0` but is filtered
    // by the `--min-map-q 30` case, so the chain path's filter funnel (the ported
    // "Filtered out N templates" diagnostic + the filter-count metric columns) is
    // exercised with a non-zero count instead of only the all-zero case.
    records.extend(create_duplicate_group_mapq(&format!("{base}d"), "GGGGGGGG", 2, start, 10));
    records
}

/// Read a BAM's `@HD` record (declared sort order). The `@PG` command-line field
/// legitimately differs between the chain and non-chain invocations (different
/// `--threads` args), so parity checks compare `@HD` rather than the whole header.
fn read_bam_hd(path: &Path) -> Option<String> {
    let mut reader = bam::io::Reader::new(fs::File::open(path).unwrap());
    let header = reader.read_header().unwrap();
    header.header().map(|hd| format!("{hd:?}"))
}

/// The chain (`--threads`) path matches the non-chain path across the
/// output-changing knobs the identity-only parity tests leave uncovered:
/// - the CLI-default `adjacency` strategy (the common `dedup --threads N`
///   invocation) and `--strategy edit`, on mixed UMIs so both do real
///   clustering work rather than collapsing to identity (vacuous);
/// - `--remove-duplicates` (the serialize-step drop path);
/// - `--no-umi`, whose handling this diff specifically changed (it forces
///   identity/edits=0 and flips `filter_config.no_umi` on the chain path);
/// - `--min-map-q 30`, which filters the fixture's mapq-10 subfamily, exercising
///   the chain's ported filter funnel with a non-zero filtered count.
/// The `@HD` sort-order header is compared too (the `@PG` command-line field
/// legitimately differs between the two invocations, so it is excluded).
#[rstest]
#[case::identity(&["--strategy", "identity"])]
#[case::adjacency_default(&[])]
#[case::edit(&["--strategy", "edit"])]
#[case::adjacency_remove(&["--strategy", "adjacency", "--remove-duplicates"])]
#[case::identity_remove(&["--strategy", "identity", "--remove-duplicates"])]
#[case::no_umi(&["--no-umi"])]
#[case::min_map_q_filters(&["--strategy", "identity", "--min-map-q", "30"])]
fn test_dedup_chain_matches_non_chain_across_knobs(#[case] extra: &[&str]) {
    let temp_dir = TempDir::new().unwrap();
    let input_bam = temp_dir.path().join("input.bam");

    let mut records = Vec::new();
    for i in 0..12 {
        records.extend(create_chain_parity_group(&format!("g{i}"), 100 + i * 200));
    }
    create_sorted_bam(&input_bam, records);

    let oracle_out = temp_dir.path().join("oracle.bam");
    dedup_run(&input_bam, &oracle_out, extra);

    let chain_out = temp_dir.path().join("chain.bam");
    let mut chain_args = extra.to_vec();
    chain_args.extend_from_slice(&["--threads", "4"]);
    dedup_run(&input_bam, &chain_out, &chain_args);

    let expected = read_deduped_records(&oracle_out);
    let actual = read_deduped_records(&chain_out);
    assert!(!expected.is_empty(), "oracle output must be non-empty (guard against a vacuous pass)");
    assert_eq!(
        actual, expected,
        "chain --threads 4 output must match the non-chain path for knobs {extra:?}"
    );
    assert_eq!(
        read_bam_hd(&chain_out),
        read_bam_hd(&oracle_out),
        "chain and non-chain must declare the same @HD sort order for knobs {extra:?}"
    );
}

/// Non-vacuity guard for the adjacency parity cases above: on the mixed-UMI
/// fixture, `--strategy adjacency` (which merges the edit-distance-1 UMI pair)
/// must produce *different* output than `--strategy identity` (which keeps all
/// three UMIs distinct). If this ever fails, the fixture stopped exercising
/// adjacency's distinguishing logic and `case_2_adjacency_default` /
/// `case_3_adjacency_remove` would be passing vacuously.
#[test]
fn test_dedup_mixed_umi_fixture_distinguishes_adjacency_from_identity() {
    let temp_dir = TempDir::new().unwrap();
    let input_bam = temp_dir.path().join("input.bam");

    let mut records = Vec::new();
    for i in 0..12 {
        records.extend(create_chain_parity_group(&format!("g{i}"), 100 + i * 200));
    }
    create_sorted_bam(&input_bam, records);

    let identity_out = temp_dir.path().join("identity.bam");
    dedup_run(&input_bam, &identity_out, &["--strategy", "identity"]);
    let adjacency_out = temp_dir.path().join("adjacency.bam");
    dedup_run(&input_bam, &adjacency_out, &["--strategy", "adjacency"]);

    assert_ne!(
        read_deduped_records(&identity_out),
        read_deduped_records(&adjacency_out),
        "the mixed-UMI fixture must make adjacency cluster differently than identity, \
         else the adjacency parity cases pass vacuously"
    );
}

/// The `--check-crc` / `--no-check-crc` / default CRC policy is honored on the
/// chain (`--threads`) path. The accept case asserts record identity against an
/// intact-file baseline run *also on the chain path* — not merely a non-empty
/// output (the vacuous-pass trap).
#[rstest]
#[case::no_check_crc_accepts(&["--no-check-crc"], true)]
#[case::check_crc_rejects(&["--check-crc"], false)]
#[case::default_rejects(&[], false)]
fn test_dedup_threaded_crc_policy(#[case] crc_args: &[&str], #[case] expect_ok: bool) {
    let temp_dir = TempDir::new().unwrap();
    let input_bam = temp_dir.path().join("input.bam");
    let clean_bam = temp_dir.path().join("clean.bam");

    // Enough duplicate pairs to span multiple BGZF blocks (see
    // `corrupt_last_block_crc`).
    create_sorted_bam(&input_bam, create_duplicate_group("dup", "ACGTACGT", 400, 100));
    fs::copy(&input_bam, &clean_bam).expect("copy pristine input");

    // For the accept case, capture the intact-file chain output as the baseline
    // BEFORE corrupting, so the accept assertion is record-identity, not merely
    // "non-empty".
    let baseline = expect_ok.then(|| {
        let baseline_out = temp_dir.path().join("baseline.bam");
        dedup_run(&clean_bam, &baseline_out, &["--strategy", "identity", "--threads", "4"]);
        read_deduped_records(&baseline_out)
    });

    corrupt_last_block_crc(&input_bam);

    let output_bam = temp_dir.path().join("output.bam");
    let mut args = vec![
        "dedup",
        "--input",
        input_bam.to_str().unwrap(),
        "--output",
        output_bam.to_str().unwrap(),
        "--strategy",
        "identity",
        "--threads",
        "4",
    ];
    args.extend_from_slice(crc_args);
    let result = MarkDuplicates::try_parse_from(args)
        .expect("failed to parse dedup args")
        .execute("fgumi dedup");

    if expect_ok {
        result.expect("--no-check-crc on the chain path must accept a corrupted BGZF CRC32");
        let records = read_deduped_records(&output_bam);
        assert_eq!(
            baseline.expect("baseline is captured for the accept case"),
            records,
            "chain --no-check-crc output must be identical to the intact-file run"
        );
    } else {
        let err = result.expect_err("chain path must reject a corrupted BGZF CRC32");
        let message = format!("{err:#}");
        assert!(message.to_uppercase().contains("CRC32"), "error should mention CRC32: {message}");
    }
}

/// The `--duplication-ladder` (and `--metrics` / `--family-size-histogram`)
/// output from the chain path is byte-identical to the non-chain path.
///
/// This MUST run multi-threaded: the ladder is order-sensitive (it samples a
/// saturation curve at cumulative-template intervals), and its recording seam
/// is the chain's serial `MiAssign` step, which only actually reorders at
/// `--threads > 1`. At `--threads 1` a reorder regression would slip through.
/// The input carries several hundred position groups (so batches span many
/// in-flight batches at `--threads 4`) across two libraries (so per-library
/// ladder rows are exercised), with a small `--ladder-interval` so rows are
/// sampled mid-stream rather than only at `finish()`.
#[test]
fn test_dedup_threaded_duplication_ladder_parity() {
    let temp_dir = TempDir::new().unwrap();
    let input_bam = temp_dir.path().join("input.bam");

    let header = create_multi_library_header("chr1", 10000);
    let mut records = Vec::new();
    for i in 0..400i32 {
        let start = 100 + i * 20;
        let (rg, umi) = if i % 2 == 0 { ("RG1", "AAAAAAAA") } else { ("RG2", "CCCCCCCC") };
        // Non-uniform group sizes (2..6 templates) are essential: a ladder built
        // from uniform increments is insensitive to group order, so it could not
        // detect a reordering regression. Varying the per-group (templates,
        // duplicates) makes the sampled saturation curve depend on the exact
        // coordinate order the recorder observes.
        let count = 2 + usize::try_from(i % 5).expect("0..5 fits usize");
        records.extend(create_duplicate_group_with_rg(&format!("p{i}"), umi, count, start, rg));
    }
    create_sorted_bam_with_header(&input_bam, &header, records);

    // Run dedup writing the BAM plus all three metric outputs; `extra` carries
    // the thread flags (empty = non-chain oracle).
    let run = |tag: &str, extra: &[&str]| {
        let out = temp_dir.path().join(format!("{tag}.bam"));
        let ladder = temp_dir.path().join(format!("{tag}.ladder.txt"));
        let metrics = temp_dir.path().join(format!("{tag}.metrics.txt"));
        let hist = temp_dir.path().join(format!("{tag}.hist.txt"));
        let mut args = vec![
            "--strategy",
            "identity",
            "--duplication-ladder",
            ladder.to_str().unwrap(),
            "--ladder-interval",
            "7",
            "--metrics",
            metrics.to_str().unwrap(),
            "--family-size-histogram",
            hist.to_str().unwrap(),
            "--compression-level",
            "1",
        ];
        args.extend_from_slice(extra);
        dedup_run(&input_bam, &out, &args);
        (out, ladder, metrics, hist)
    };

    let (oracle_out, oracle_ladder, oracle_metrics, oracle_hist) = run("oracle", &[]);
    let (chain_out, chain_ladder, chain_metrics, chain_hist) = run("chain", &["--threads", "4"]);

    // Sanity: the ladder is non-vacuous (mid-stream rows actually landed).
    assert!(
        !read_dedup_metrics_rows(&oracle_ladder).is_empty(),
        "ladder must have rows, else parity is vacuous"
    );

    // The ladder is the order-sensitive file — its byte-identity is the point of
    // this test. Metrics and family-size histogram are order-insensitive
    // aggregates but checked too for completeness.
    assert_eq!(
        fs::read(&oracle_ladder).unwrap(),
        fs::read(&chain_ladder).unwrap(),
        "duplication ladder diverged between the chain and non-chain paths"
    );
    assert_eq!(
        fs::read(&oracle_metrics).unwrap(),
        fs::read(&chain_metrics).unwrap(),
        "metrics diverged between the chain and non-chain paths"
    );
    assert_eq!(
        fs::read(&oracle_hist).unwrap(),
        fs::read(&chain_hist).unwrap(),
        "family-size histogram diverged between the chain and non-chain paths"
    );
    assert_eq!(
        read_deduped_records(&oracle_out),
        read_deduped_records(&chain_out),
        "output records diverged between the chain and non-chain paths"
    );
}

/// Complement to the `group` #901 regression, on a fixture that exercises the
/// full PCR-duplicate marking contract rather than mere retention. `dedup` must:
///
/// 1. Mark the correct read in each duplicate family: the representative (highest
///    base-quality template) is kept unflagged and every other template in the
///    family gets the `PCR/optical duplicate` flag (`0x400`).
/// 2. *Retain* a `tc`-keyed secondary/supplementary read (the opposite of `group`,
///    which filters it) and propagate its primary's duplicate status onto it —
///    guarding against a regression that hoists the secondary/supplementary filter
///    above the `include_secondary_supplementary` guard, or that marks primaries
///    but skips the coalesced non-primary record.
///
/// The fixture is 10 paired templates across four independent duplicate families
/// (distinct coordinate + UMI, so they never cross-contaminate) of sizes 1, 4, 3,
/// and 2. Within each family every template gets a *distinct* base quality, so the
/// representative is the unique maximum — deterministic, with no score tie whose
/// resolution would depend on stream order. Every one of the 20 primary reads is
/// asserted to carry exactly the flag its family membership dictates. The flagged
/// alignment is attached to a marked duplicate (`famD_dup`), so its output flag
/// pins the propagation this test's ancestor could not: the single-template
/// fixture it replaced never marked any record duplicate.
#[rstest]
#[case::secondary(flags::SECONDARY)]
#[case::supplementary(flags::SUPPLEMENTARY)]
fn test_dedup_marks_families_and_flags_tc_keyed_secondary_supplementary(#[case] extra_flag: u16) {
    let temp_dir = TempDir::new().unwrap();
    let input_bam = temp_dir.path().join("input.bam");
    let output_bam = temp_dir.path().join("output.bam");

    // Four duplicate families at distinct coordinates and UMIs, of sizes 1/4/3/2.
    // Each entry is (UMI, start position, templates), where a template is a
    // (QNAME base, base quality) pair — `create_template_qual` builds one template
    // per pair, named `{base}_0`. Within a family the unique top quality (40) is the
    // representative and the rest (35/30/25) are marked duplicate.
    #[expect(clippy::type_complexity, reason = "an inline fixture table of duplicate families")]
    let families: [(&str, i32, &[(&str, u8)]); 4] = [
        // Family A: a lone template — a singleton is never a duplicate.
        ("AAAAAAAA", 100, &[("famA_rep", 40)]),
        // Family B: 4 templates -> 1 kept + 3 marked.
        (
            "CCCCCCCC",
            500,
            &[("famB_rep", 40), ("famB_dup0", 35), ("famB_dup1", 30), ("famB_dup2", 25)],
        ),
        // Family C: 3 templates -> 1 kept + 2 marked.
        ("GGGGGGGG", 900, &[("famC_rep", 40), ("famC_dup0", 35), ("famC_dup1", 30)]),
        // Family D: 2 templates -> 1 kept + 1 marked. This is the highest
        // coordinate, so family D's reads sort last; the flagged alignment below
        // (own position 5000, past every primary) sorts to the very end and
        // coalesces — by QNAME adjacency — onto whichever family-D primary sorts
        // last. `famD_dup` is named so it sorts after `famD_rep` (the
        // template-coordinate tie-break below the shared coordinate is `name_hash`,
        // not lexicographic), making the marked duplicate the trailing record. If
        // that ordering ever changes, the flagged read coalesces onto the wrong
        // primary (or is dropped) and the assertions below fail loudly.
        ("TTTTTTTT", 1300, &[("famD_rep", 40), ("famD_dup", 35)]),
    ];

    // Names that must be kept (representatives) vs. marked duplicate, derived from
    // the family table so the expectation and the fixture cannot drift.
    let mut expected_kept: Vec<String> = Vec::new();
    let mut expected_duplicate: Vec<String> = Vec::new();
    let mut records = Vec::new();
    for (umi, start, templates) in families {
        // Representative = the unique maximum base quality in the family.
        let rep_qual = templates.iter().map(|&(_, q)| q).max().expect("family is non-empty");
        for &(name, qual) in templates {
            records.extend(create_template_qual(name, umi, start, qual));
            let read_name = format!("{name}_0");
            if qual == rep_qual {
                expected_kept.push(read_name);
            } else {
                expected_duplicate.push(read_name);
            }
        }
    }

    // A `tc`-keyed secondary/supplementary alignment of `famD_dup` (a marked
    // duplicate). Its `tc` tag carries family D's template coordinate — R1 fwd
    // unclipped-5' = 1300, R2 rev unclipped-5' = 1407 -> [0,1300,0,0,1407,1] — the
    // resolved key dedup uses to place it in family D's molecule. dedup coalesces
    // it into `famD_dup`'s template, where it must inherit that template's
    // duplicate flag.
    let supp_primary = "famD_dup_0";
    let mut supp = SamBuilder::new();
    supp.read_name(supp_primary.as_bytes())
        .sequence(b"ACGTACGT")
        .qualities(&[30; 8])
        .flags(flags::PAIRED | flags::FIRST_SEGMENT | extra_flag)
        .ref_id(0)
        .pos(5000)
        .mapq(60)
        .cigar_ops(&[8u32 << 4])
        .mate_ref_id(0)
        .mate_pos(1399)
        .add_string_tag(SamTag::RX, b"TTTTTTTT")
        .add_string_tag(SamTag::MC, b"8M")
        .add_array_i32(SamTag::TC, &[0, 1300, 0, 0, 1407, 1]);
    records.push(supp.build());

    create_sorted_bam(&input_bam, records);
    let cmd = MarkDuplicates::try_parse_from([
        "dedup",
        "--input",
        input_bam.to_str().unwrap(),
        "--output",
        output_bam.to_str().unwrap(),
        "--strategy",
        "identity",
        "--compression-level",
        "1",
    ])
    .expect("failed to parse dedup args");
    cmd.execute("fgumi dedup").expect("dedup should succeed");

    // Collect the full output: primary reads keyed by QNAME -> the duplicate flag on
    // each of that template's records (R1 and R2 must agree), plus the retained
    // secondary/supplementary alignment.
    let mut reader = bam::io::Reader::new(fs::File::open(&output_bam).unwrap());
    let header = reader.read_header().unwrap();
    let mut primary_flags: std::collections::HashMap<String, Vec<bool>> =
        std::collections::HashMap::new();
    let mut secondary_supplementary: Vec<(String, bool)> = Vec::new();
    let mut total = 0usize;
    for record in reader.record_bufs(&header) {
        let record = record.expect("read dedup output record");
        let flag = record.flags();
        let name = record.name().map(ToString::to_string).expect("output record has a name");
        total += 1;
        if flag.is_secondary() || flag.is_supplementary() {
            secondary_supplementary.push((name, flag.is_duplicate()));
        } else {
            primary_flags.entry(name).or_default().push(flag.is_duplicate());
        }
    }

    // 10 templates * 2 mates + 1 flagged alignment, nothing removed (mark mode).
    assert_eq!(total, 21, "mark mode retains all 20 primary reads plus the flagged alignment");
    assert_eq!(primary_flags.len(), 10, "every one of the 10 templates must be present by QNAME");

    // Every primary template carries exactly the flag its family membership dictates,
    // on BOTH mates — checked read-by-read, not merely by count.
    for name in &expected_kept {
        assert_eq!(
            primary_flags.get(name).map(Vec::as_slice),
            Some([false, false].as_slice()),
            "representative {name} must keep both mates unflagged"
        );
    }
    for name in &expected_duplicate {
        assert_eq!(
            primary_flags.get(name).map(Vec::as_slice),
            Some([true, true].as_slice()),
            "duplicate {name} must have the PCR-duplicate flag on both mates"
        );
    }

    // The flagged alignment is retained, belongs to its `famD_dup` primary, and
    // inherits that primary's duplicate flag — the propagation this test pins.
    assert_eq!(
        secondary_supplementary,
        vec![(supp_primary.to_string(), true)],
        "the tc-keyed sec/supp read must be retained as {supp_primary}'s and marked duplicate"
    );
}

// ============================================================================
// --verify (strict template-coordinate sort-order gate)
// ============================================================================

/// Build `count` paired-end duplicate templates with INTERNALLY CONSISTENT
/// mate-strand flags: R1 forward with `MATE_REVERSE` set, R2 reverse with its
/// mate forward. Unlike [`create_duplicate_group`], whose R1 omits
/// `MATE_REVERSE`, this makes R1 and R2 resolve to the *same* template
/// coordinate — a prerequisite for the strict `--verify` order check, which keys
/// on the exact `fgumi sort --order template-coordinate` key. (An inconsistent
/// pair splits into two different coordinates, which the strict check then reads
/// as an out-of-order file.)
fn create_consistent_pair_group(
    base_name: &str,
    umi: &str,
    count: usize,
    start: i32,
) -> Vec<RawRecord> {
    let mut records = Vec::new();
    for i in 0..count {
        let name = format!("{base_name}_{i}");
        let r1 = {
            let mut b = SamBuilder::new();
            b.read_name(name.as_bytes())
                .sequence(b"ACGTACGT")
                .qualities(&[30; 8])
                .flags(flags::PAIRED | flags::FIRST_SEGMENT | flags::MATE_REVERSE)
                .ref_id(0)
                .pos(start - 1)
                .mapq(60)
                .cigar_ops(&[8 << 4]) // 8M
                .mate_ref_id(0)
                .mate_pos(start + 99)
                .template_length(108)
                .add_string_tag(SamTag::RX, umi.as_bytes())
                .add_string_tag(SamTag::MC, b"8M");
            b.build()
        };
        let r2 = {
            let mut b = SamBuilder::new();
            b.read_name(name.as_bytes())
                .sequence(b"ACGTACGT")
                .qualities(&[30; 8])
                .flags(flags::PAIRED | flags::LAST_SEGMENT | flags::REVERSE)
                .ref_id(0)
                .pos(start + 99)
                .mapq(60)
                .cigar_ops(&[8 << 4]) // 8M
                .mate_ref_id(0)
                .mate_pos(start - 1)
                .template_length(-108)
                .add_string_tag(SamTag::RX, umi.as_bytes())
                .add_string_tag(SamTag::MC, b"8M");
            b.build()
        };
        records.push(r1);
        records.push(r2);
    }
    records
}

/// Write `records` verbatim (NO sorting) under a header that advertises
/// `SS:template-coordinate`. Unlike [`create_sorted_bam`], which shells out to
/// `fgumi sort`, this preserves the caller's record order, so the header-level
/// ordering check passes while the records can be genuinely out of order — the
/// exact shape `--verify` exists to catch.
fn write_bam_tc_header_unsorted(path: &Path, records: &[RawRecord]) {
    let header = create_minimal_header("chr1", 10000);
    let mut writer =
        bam::io::Writer::new(fs::File::create(path).expect("Failed to create BAM file"));
    writer.write_header(&header).expect("Failed to write header");
    for record in records {
        writer
            .write_alignment_record(&header, &to_record_buf(record))
            .expect("Failed to write record");
    }
    writer.try_finish().expect("Failed to finish BAM");
}

/// Run `dedup --verify` (plus `threads`), returning the raw execute result so a
/// test can assert on either success or the ordering error.
fn dedup_verify_result(input: &Path, output: &Path, threads: &[&str]) -> anyhow::Result<()> {
    let mut args = vec![
        "dedup",
        "--input",
        input.to_str().unwrap(),
        "--output",
        output.to_str().unwrap(),
        "--verify",
    ];
    args.extend_from_slice(threads);
    MarkDuplicates::try_parse_from(args).expect("failed to parse dedup args").execute("fgumi dedup")
}

/// `--verify` on a correctly template-coordinate-sorted input must succeed and
/// produce output record-for-record identical to the same run without `--verify`
/// — the flag is a precondition gate, not a mode, so it must not perturb the
/// result (and `--output` is still written). Covers the non-chain path and the
/// `--threads N` chain path.
#[rstest]
#[case::single_threaded(&[])]
#[case::chain_threads2(&["--threads", "2"])]
fn test_dedup_verify_accepts_sorted_input(#[case] threads: &[&str]) {
    let temp_dir = TempDir::new().unwrap();
    let input_bam = temp_dir.path().join("input.bam");
    // Ascending start positions → correctly template-coordinate sorted.
    let mut records = create_consistent_pair_group("dup1", "ACGTACGT", 3, 100);
    records.extend(create_consistent_pair_group("dup2", "TGCATGCA", 2, 500));
    create_sorted_bam(&input_bam, records);

    let baseline_out = temp_dir.path().join("baseline.bam");
    dedup_run(&input_bam, &baseline_out, threads);
    let verified_out = temp_dir.path().join("verified.bam");
    let mut extra = threads.to_vec();
    extra.push("--verify");
    dedup_run(&input_bam, &verified_out, &extra);

    let baseline = read_deduped_records(&baseline_out);
    assert_eq!(baseline.len(), 10, "sanity: all 10 reads present (duplicates marked, not removed)");
    assert_eq!(
        baseline,
        read_deduped_records(&verified_out),
        "--verify must not change dedup output on a correctly-sorted input",
    );
}

/// `--verify` on an input that is out of template-coordinate order (despite a
/// header that advertises it) must abort with a non-zero exit and an error that
/// names the ordering requirement. The same input WITHOUT `--verify` is accepted
/// — proving the gate is what rejects it. Covers the non-chain and `--threads N`
/// chain paths.
#[rstest]
#[case::single_threaded(&[])]
#[case::chain_threads2(&["--threads", "2"])]
fn test_dedup_verify_rejects_out_of_order_input(#[case] threads: &[&str]) {
    let temp_dir = TempDir::new().unwrap();
    let input_bam = temp_dir.path().join("input.bam");
    // Descending start positions → genuinely out of order under a TC-advertising
    // header (the group at 500 is written before the group at 100). Written
    // verbatim — NOT through `fgumi sort`, which would reorder it into order.
    let mut records = create_consistent_pair_group("late", "ACGTACGT", 2, 500);
    records.extend(create_consistent_pair_group("early", "TGCATGCA", 2, 100));
    write_bam_tc_header_unsorted(&input_bam, &records);

    // Without --verify the out-of-order input is accepted (no ordering guard).
    let plain_out = temp_dir.path().join("plain.bam");
    dedup_run(&input_bam, &plain_out, threads);
    assert!(
        !read_deduped_records(&plain_out).is_empty(),
        "sanity: without --verify the input is deduplicated without error",
    );

    // With --verify it is rejected before completing.
    let verified_out = temp_dir.path().join("verified.bam");
    let err = dedup_verify_result(&input_bam, &verified_out, threads)
        .expect_err("out-of-order input must be rejected under --verify");
    let message = format!("{err:#}");
    assert!(
        message.contains("template-coordinate sort order"),
        "error should name the ordering requirement, got: {message}",
    );
}

/// dedup's chain verify path must accept a correctly-sorted input that spans many
/// pipeline batches under multi-worker reordering — the sibling of group's
/// `test_group_verify_accepts_sorted_multi_batch_chain`. Both commands wire the
/// same serial `GroupByPosition` step, so the reorder-across-batches risk is
/// identical; dedup additionally uses `with_secondary_supplementary`, so this also
/// exercises verify on that grouper variant across batch boundaries.
#[test]
fn test_dedup_verify_accepts_sorted_multi_batch_chain() {
    let temp_dir = TempDir::new().unwrap();
    let input_bam = temp_dir.path().join("input.bam");
    // 600 single-pair templates at distinct ascending positions → many position
    // groups spanning multiple pipeline batches (both the position-group
    // accumulator and the upstream decode batch).
    let mut records = Vec::new();
    for i in 0..600i32 {
        records.extend(create_consistent_pair_group(&format!("t{i}"), "ACGTACGT", 1, 100 + i * 10));
    }
    create_sorted_bam(&input_bam, records);

    let baseline_out = temp_dir.path().join("baseline.bam");
    dedup_run(&input_bam, &baseline_out, &["--threads", "4"]);
    let verified_out = temp_dir.path().join("verified.bam");
    dedup_run(&input_bam, &verified_out, &["--threads", "4", "--verify"]);

    let baseline = read_deduped_records(&baseline_out);
    assert_eq!(
        baseline.len(),
        1200,
        "600 pairs → 1200 reads present (duplicates marked, not removed)"
    );
    assert_eq!(
        baseline,
        read_deduped_records(&verified_out),
        "--verify must accept a correctly-sorted input spanning many pipeline batches under \
         multi-worker reordering, and must not change the output",
    );
}

/// Build a single-end mapped read carrying `RG` and `CB` tags at `start`
/// (1-based). Used to exercise the cell-barcode lane of the template-coordinate
/// sort key, which `--verify` must key on to match `fgumi sort`'s default.
fn single_end_rg_cb(name: &str, start: i32, rg_id: &str, cb: &str) -> RawRecord {
    let mut b = SamBuilder::new();
    b.read_name(name.as_bytes())
        .sequence(b"ACGTACGT")
        .qualities(&[30; 8])
        .flags(0) // unpaired, mapped, forward
        .ref_id(0)
        .pos(start - 1)
        .mapq(60)
        .cigar_ops(&[8 << 4]) // 8M
        .add_string_tag(SamTag::RX, b"ACGTACGT")
        .add_string_tag(SamTag::RG, rg_id.as_bytes())
        .add_string_tag(SamTag::CB, cb.as_bytes());
    b.build()
}

/// Write `records` under `header`, then sort with `fgumi sort --order
/// template-coordinate` using the DEFAULT `--key-types` (auto-detect, which
/// INCLUDES the CB lane for CB-tagged input). Unlike [`create_sorted_bam`], which
/// passes `--key-types mi` and drops the CB lane, this reproduces the exact order
/// the command `--verify`'s help/error tells users to run.
fn sort_default_keytypes(
    temp_dir: &Path,
    header: &noodles::sam::Header,
    records: &[RawRecord],
) -> std::path::PathBuf {
    let unsorted = temp_dir.join("cb_unsorted.bam");
    let mut writer =
        bam::io::Writer::new(fs::File::create(&unsorted).expect("Failed to create BAM file"));
    writer.write_header(header).expect("Failed to write header");
    for r in records {
        writer.write_alignment_record(header, &to_record_buf(r)).expect("Failed to write record");
    }
    writer.try_finish().expect("Failed to finish BAM");

    let sorted = temp_dir.join("cb_sorted.bam");
    let status = std::process::Command::new(env!("CARGO_BIN_EXE_fgumi"))
        .args([
            "sort",
            "-i",
            unsorted.to_str().unwrap(),
            "-o",
            sorted.to_str().unwrap(),
            "--order",
            "template-coordinate",
        ])
        .status()
        .expect("failed to spawn `fgumi sort` for CB test input");
    assert!(status.success(), "failed to template-coordinate sort CB test input");
    sorted
}

/// Regression for the CB-keying bug: `fgumi sort --order template-coordinate` keys
/// on the `CB` tag by default (`parse_cell_tag`), and `cb_hash` outranks the
/// library lane in `TemplateKey::core_cmp`. Two same-position reads whose `CB`
/// order is the OPPOSITE of their library order are emitted by that sort in `CB`
/// order; `--verify` must accept it because it keys on `CB` too. A CB-blind
/// verifier (`cell_tag = None`) would see the library lane out of order and
/// spuriously reject. `cb_hash(AAAAAAAA) < cb_hash(TTTTTTTT)` under fgumi's
/// fixed-seed hasher, so pairing `CB=AAAAAAAA` with the higher library ordinal
/// (`libB` = `RG2`) makes the CB-sorted order run library-descending — the exact
/// shape that trips a `None` verifier. (Sorted via `fgumi sort` itself, so the
/// input is genuinely what the recommended command produces, whatever the hash.)
#[test]
fn test_dedup_verify_accepts_cb_sorted_multi_library() {
    let temp_dir = TempDir::new().unwrap();
    let header = create_multi_library_header("chr1", 10000);
    let records = [
        single_end_rg_cb("readA", 501, "RG1", "TTTTTTTT"), // libA (ordinal 1), high cb_hash
        single_end_rg_cb("readB", 501, "RG2", "AAAAAAAA"), // libB (ordinal 2), low cb_hash
    ];
    let sorted = sort_default_keytypes(temp_dir.path(), &header, &records);

    let out = temp_dir.path().join("out.bam");
    dedup_verify_result(&sorted, &out, &[])
        .expect("--verify must accept the CB-sorted multi-library input `fgumi sort` produced");
}

/// A coordinate-ordered (mate-interleaved) paired file is a realistic mislabel:
/// records physically in `SO:coordinate` order but under a header advertising
/// `SS:template-coordinate`. That is NOT template-coordinate order — coordinate
/// sort interleaves mates by each read's own leftmost position instead of
/// grouping a template's two reads at their shared canonical coordinate — so
/// `--verify` must reject it. (An *honestly*-labeled non-TC file, i.e. an
/// `SO:coordinate`/`SO:queryname` header, is already rejected upstream by the
/// header-ordering check before `--verify` runs; this exercises the record-order
/// discontinuity that check cannot see.)
#[test]
fn test_dedup_verify_rejects_coordinate_ordered_input() {
    let temp_dir = TempDir::new().unwrap();
    let input_bam = temp_dir.path().join("input.bam");
    // Two templates. Template-coordinate order keeps each mate pair adjacent
    // (T1.R1,T1.R2,T2.R1,T2.R2); coordinate order interleaves them by leftmost
    // position: R1@100, R1@150, R2@200, R2@250.
    let t1 = create_consistent_pair_group("T1", "ACGTACGT", 1, 100); // [R1@100, R2@200]
    let t2 = create_consistent_pair_group("T2", "TGCATGCA", 1, 150); // [R1@150, R2@250]
    let coordinate_order = vec![t1[0].clone(), t2[0].clone(), t1[1].clone(), t2[1].clone()];
    write_bam_tc_header_unsorted(&input_bam, &coordinate_order);

    let out = temp_dir.path().join("out.bam");
    let err = dedup_verify_result(&input_bam, &out, &[])
        .expect_err("coordinate-ordered (mate-interleaved) input must be rejected under --verify");
    assert!(
        format!("{err:#}").contains("template-coordinate sort order"),
        "error should name the ordering requirement, got: {err:#}",
    );
}
