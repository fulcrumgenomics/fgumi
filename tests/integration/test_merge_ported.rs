//! The retired standalone engine's merge tests (`fgumi-sort`'s `external.rs`),
//! ported to the `fgumi merge` command: each fixture is sorted with `fgumi
//! sort` and merged with `fgumi merge`; every assertion is the original's, and
//! each merge is also checked record for record against the stable sort of its
//! inputs' concatenation (`assert_merge_is_a_stable_sort`), so a dropped or
//! duplicated record fails even where the original checked only counts or
//! order.

use std::path::{Path, PathBuf};

use fgumi_bam_io::create_raw_bam_reader;
use fgumi_sam::SamBuilder;
use fgumi_sort::{RawReadAheadReader, cb_hasher, extract_template_key_inline};

use crate::helpers::chain_spec::logged_count;
use crate::helpers::merge_fixtures::{merge_cli, read_raw, sort_for};

/// `num_pairs` pairs named `{prefix}_read{i:04}` starting at `start_offset`,
/// sorted into `order` (a `fgumi sort --order` value) by `fgumi sort`.
fn create_sorted_bam(
    dir: &Path,
    prefix: &str,
    num_pairs: usize,
    start_offset: usize,
    order: &str,
) -> (PathBuf, Vec<String>) {
    let mut builder = SamBuilder::new();
    let mut names = Vec::with_capacity(num_pairs);
    for i in 0..num_pairs {
        let name = format!("{prefix}_read{i:04}");
        names.push(name.clone());
        let _ = builder
            .add_pair()
            .name(&name)
            .start1((start_offset + i * 200) + 1)
            .start2((start_offset + i * 200) + 101)
            .build();
    }
    let unsorted = dir.join(format!("{prefix}_unsorted.bam"));
    let sorted = dir.join(format!("{prefix}_sorted.bam"));
    builder.write_bam(&unsorted).expect("failed to write BAM");
    let st = std::process::Command::new(env!("CARGO_BIN_EXE_fgumi"))
        .args(["sort", "-i"])
        .arg(&unsorted)
        .arg("-o")
        .arg(&sorted)
        .args(["--order", order, "--compression-level", "0"])
        .status()
        .expect("run fgumi sort");
    assert!(st.success(), "fgumi sort failed");
    (sorted, names)
}

fn count_bam_records(path: &Path) -> u64 {
    let (reader, _) = create_raw_bam_reader(path, 1).expect("failed to create raw BAM reader");
    RawReadAheadReader::new(reader).count() as u64
}

fn collect_read_names(path: &Path) -> Vec<String> {
    let (reader, _) = create_raw_bam_reader(path, 1).expect("failed to create raw BAM reader");
    RawReadAheadReader::new(reader)
        .map(|rec| {
            let name_bytes = fgumi_raw_bam::RawRecordView::new(rec.as_ref()).read_name();
            String::from_utf8(name_bytes.to_vec()).expect("read name should be valid UTF-8")
        })
        .collect()
}

/// Collect (`ref_id`, pos) tuples for every record in a BAM.
fn collect_positions(path: &Path) -> Vec<(i32, i32)> {
    let (reader, _) = create_raw_bam_reader(path, 1).expect("failed to create raw BAM reader");
    RawReadAheadReader::new(reader)
        .map(|rec| {
            let v = fgumi_raw_bam::RawRecordView::new(rec.as_ref());
            (v.ref_id(), v.pos())
        })
        .collect()
}

/// A record's identity in a merge: (name, flags, `ref_id`, pos).
type RecordId = (Vec<u8>, u16, i32, i32);

fn record_id(record: &[u8]) -> RecordId {
    let v = fgumi_raw_bam::RawRecordView::new(record);
    (v.read_name().to_vec(), v.flags(), v.ref_id(), v.pos())
}

/// The merge's record-identity contract: `merged` is, record for record, the
/// stable sort in `order` of `inputs` concatenated in input order (equal keys
/// keep input order), so a dropped, duplicated or substituted record fails.
fn assert_merge_is_a_stable_sort(inputs: &[PathBuf], merged: &Path, order: &str) {
    let sort_order = match order {
        "coordinate" => fgumi_sort::SortOrder::Coordinate,
        "queryname" => {
            fgumi_sort::SortOrder::Queryname(fgumi_sort::QuerynameComparator::Lexicographic)
        }
        "queryname::natural" => {
            fgumi_sort::SortOrder::Queryname(fgumi_sort::QuerynameComparator::Natural)
        }
        "template-coordinate" => fgumi_sort::SortOrder::TemplateCoordinate,
        other => panic!("unknown order {other}"),
    };
    let mut header = None;
    let mut concat: Vec<fgumi_raw_bam::RawRecord> = Vec::new();
    for input in inputs {
        let (h, records) = read_raw(input);
        header.get_or_insert(h);
        concat.extend(records.into_iter().map(fgumi_raw_bam::RawRecord::from));
    }
    let header = header.expect("at least one input");
    let want: Vec<RecordId> =
        sort_for(sort_order, &header, &concat).iter().map(|r| record_id(r)).collect();
    let got: Vec<RecordId> = read_raw(merged).1.iter().map(|r| record_id(r)).collect();
    assert_eq!(got.len(), want.len(), "record count");
    for (i, (g, w)) in got.iter().zip(&want).enumerate() {
        assert_eq!(g, w, "record {i} differs from the stable sort of the inputs");
    }
}

/// Merge `inputs` in `order`, check the merge's record identity, and return
/// the logged record count.
fn merge(inputs: &[PathBuf], merged: &Path, order: &str) -> u64 {
    let log = merge_cli(inputs, merged, order, 1);
    assert_merge_is_a_stable_sort(inputs, merged, order);
    logged_count(&log, "Records merged:").expect("merge logs its record count")
}

#[test]
fn test_merge_ported_coordinate_sort() {
    let dir = tempfile::tempdir().expect("failed to create temp directory");
    let (bam_a, _) = create_sorted_bam(dir.path(), "a", 10, 0, "coordinate");
    let (bam_b, _) = create_sorted_bam(dir.path(), "b", 10, 10_000, "coordinate");

    let merged = dir.path().join("merged.bam");
    let count = merge(&[bam_a, bam_b], &merged, "coordinate");

    // 10 pairs * 2 records * 2 inputs = 40
    assert_eq!(count, 40);
    assert_eq!(count_bam_records(&merged), 40);

    // Verify coordinate sort order: (ref_id, pos) is non-decreasing
    let positions = collect_positions(&merged);
    for w in positions.windows(2) {
        assert!(w[0] <= w[1], "coordinate sort violated: {:?} > {:?}", w[0], w[1]);
    }
}

#[test]
fn test_merge_ported_template_coordinate_sort() {
    let dir = tempfile::tempdir().expect("failed to create temp directory");
    let (bam_a, _) = create_sorted_bam(dir.path(), "a", 10, 0, "template-coordinate");
    let (bam_b, _) = create_sorted_bam(dir.path(), "b", 10, 10_000, "template-coordinate");

    let merged = dir.path().join("merged.bam");
    let count = merge(&[bam_a, bam_b], &merged, "template-coordinate");

    assert_eq!(count, 40);
    assert_eq!(count_bam_records(&merged), 40);

    // Verify template-coordinate order/identity via an independent oracle:
    // re-extract the template key from every merged record and assert
    // non-decreasing core order. A count-only check would pass even if the
    // merge reordered or duplicated records.
    let (_, hdr) = create_raw_bam_reader(&merged, 1).expect("header");
    let lib = fgumi_sort::LibraryLookup::from_header(&hdr);
    let hasher = cb_hasher();
    let file = std::fs::File::open(&merged).expect("open merged");
    let mut raw_reader = fgumi_sort::RawBamRecordReader::new(file).expect("reader");
    raw_reader.skip_header().expect("skip header");
    let (total, violations, first) = fgumi_sort::verify_sort_order(
        raw_reader,
        |bam| extract_template_key_inline(bam, &lib, None, &hasher),
        |cur: &fgumi_sort::TemplateKey, prev: &fgumi_sort::TemplateKey| {
            cur.core_cmp(prev) == std::cmp::Ordering::Less
        },
    )
    .expect("verify runs");
    assert_eq!(total, 40, "verify must see every merged record");
    assert_eq!(
        violations, 0,
        "template-coordinate order violated after merge (first={first:?}, total={total})"
    );
}

#[test]
fn test_merge_ported_queryname_sort() {
    let dir = tempfile::tempdir().expect("failed to create temp directory");
    let (bam_a, _) = create_sorted_bam(dir.path(), "a", 10, 0, "queryname");
    let (bam_b, _) = create_sorted_bam(dir.path(), "b", 10, 10_000, "queryname");

    let merged = dir.path().join("merged.bam");
    let count = merge(&[bam_a, bam_b], &merged, "queryname");

    assert_eq!(count, 40);
    assert_eq!(count_bam_records(&merged), 40);

    // Verify queryname sort order: read names are non-decreasing
    let names = collect_read_names(&merged);
    for w in names.windows(2) {
        assert!(w[0] <= w[1], "queryname sort violated: {:?} > {:?}", w[0], w[1]);
    }
}

#[test]
fn test_merge_ported_single_input() {
    let dir = tempfile::tempdir().expect("failed to create temp directory");
    let (bam_a, _) = create_sorted_bam(dir.path(), "a", 15, 0, "coordinate");

    let merged = dir.path().join("merged.bam");
    let count = merge(&[bam_a], &merged, "coordinate");

    // 15 pairs * 2 = 30
    assert_eq!(count, 30);
    assert_eq!(count_bam_records(&merged), 30);
}

#[test]
fn test_merge_ported_preserves_all_records() {
    let dir = tempfile::tempdir().expect("failed to create temp directory");
    let (bam_a, names_a) = create_sorted_bam(dir.path(), "a", 5, 0, "queryname");
    let (bam_b, names_b) = create_sorted_bam(dir.path(), "b", 5, 10_000, "queryname");

    let merged = dir.path().join("merged.bam");
    merge(&[bam_a, bam_b], &merged, "queryname");

    // Every expected name (from both inputs) appears exactly twice (once per
    // mate) in the merged output, and no other name does.
    let mut merged_names = collect_read_names(&merged);
    merged_names.sort();
    let mut want: Vec<String> =
        names_a.iter().chain(names_b.iter()).flat_map(|n| [n.clone(), n.clone()]).collect();
    want.sort();
    assert_eq!(merged_names, want, "every record of both inputs, once");
}

#[test]
fn test_merge_ported_many_inputs() {
    let dir = tempfile::tempdir().expect("failed to create temp directory");
    let k = 8;
    let pairs_per_input = 5;
    let mut inputs = Vec::with_capacity(k);

    for i in 0..k {
        let (bam, _) = create_sorted_bam(
            dir.path(),
            &format!("in{i}"),
            pairs_per_input,
            i * 50_000,
            "coordinate",
        );
        inputs.push(bam);
    }

    let merged = dir.path().join("merged.bam");
    let count = merge(&inputs, &merged, "coordinate");

    let expected = (k * pairs_per_input * 2) as u64; // 8 * 5 * 2 = 80
    assert_eq!(count, expected);
    assert_eq!(count_bam_records(&merged), expected);

    // Verify coordinate sort order
    let positions = collect_positions(&merged);
    for w in positions.windows(2) {
        assert!(w[0] <= w[1], "coordinate sort violated with k={k}: {:?} > {:?}", w[0], w[1]);
    }
}

#[test]
fn test_merge_ported_queryname_natural_sort() {
    let dir = tempfile::tempdir().expect("failed to create temp directory");
    let (bam_a, _) = create_sorted_bam(dir.path(), "a", 10, 0, "queryname::natural");
    let (bam_b, _) = create_sorted_bam(dir.path(), "b", 10, 10_000, "queryname::natural");

    let merged = dir.path().join("merged.bam");
    let count = merge(&[bam_a, bam_b], &merged, "queryname::natural");

    assert_eq!(count, 40);
    assert_eq!(count_bam_records(&merged), 40);

    // Verify natural queryname order: names are non-decreasing
    let names = collect_read_names(&merged);
    for w in names.windows(2) {
        assert!(w[0] <= w[1], "natural queryname sort violated: {:?} > {:?}", w[0], w[1]);
    }
}
