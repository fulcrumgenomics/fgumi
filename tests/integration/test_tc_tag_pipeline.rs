//! End-to-end tests for the `tc` tag across `extract → zipper → sort → dedup`.
//!
//! zipper stamps each secondary and supplementary record with its template's
//! primary coordinate in the `tc` tag. The template-coordinate sort must use it
//! to place the record with its primaries, because `dedup` only sees a template's
//! records when they are adjacent and discards a record it finds alone as a
//! template without primary reads.

use std::ffi::OsStr;
use std::fmt::Write as _;
use std::fs;
use std::path::Path;

use noodles::sam::alignment::RecordBuf;
use rstest::rstest;
use tempfile::TempDir;

use crate::helpers::cli::{run_args_and_capture_logs, utf8_path as p};
use crate::helpers::{create_plain_fastq, create_test_reference, read_bam_output};

/// Length of `chr1` in the reference `create_test_reference` writes.
const CONTIG_LEN: usize = 10_000;

/// Writes the inputs for templates `names` into `dir`: the shared test
/// reference, R1/R2 FASTQs (R1 = 4 bp UMI + 10 bp insert, R2 = 10 bp insert),
/// and the aligner's output `mapped.sam` with one record per `alignments` row of
/// `(name, flag, 1-based pos, mapq, cigar, 1-based mate pos)`.
fn write_inputs(dir: &TempDir, names: &[&str], alignments: &[(&str, u16, u32, u8, &str, u32)]) {
    create_test_reference(dir.path());
    let fastq = |seq: &'static str, qual: &'static str| -> Vec<(&str, &str, &str)> {
        names.iter().map(|&n| (n, seq, qual)).collect()
    };
    create_plain_fastq(dir, "r1.fq", &fastq("AAAAAAAAAAAAAA", "IIIIIIIIIIIIII"));
    create_plain_fastq(dir, "r2.fq", &fastq("AAAAAAAAAA", "IIIIIIIIII"));

    let mut sam = format!("@HD\tVN:1.6\tSO:unsorted\n@SQ\tSN:chr1\tLN:{CONTIG_LEN}\n");
    for (name, flag, pos, mapq, cigar, mate_pos) in alignments {
        writeln!(
            sam,
            "{name}\t{flag}\tchr1\t{pos}\t{mapq}\t{cigar}\t=\t{mate_pos}\t0\tAAAAAAAAAA\tIIIIIIIIII"
        )
        .unwrap();
    }
    fs::write(dir.path().join("mapped.sam"), sam).unwrap();
}

/// Runs `fgumi <args>...`, asserting it succeeds.
fn fgumi(args: &[&str]) {
    run_args_and_capture_logs(&args.iter().map(OsStr::new).collect::<Vec<_>>());
}

/// Runs `extract → zipper → sort --order template-coordinate → dedup` in `dir`
/// and returns the sorted records and the `dedup` output records.
fn run_pipeline(dir: &Path) -> (Vec<RecordBuf>, Vec<RecordBuf>) {
    let file = |name: &str| dir.join(name);
    let (fq1, fq2, unmapped) = (file("r1.fq"), file("r2.fq"), file("unmapped.bam"));
    let (mapped, reference, zipped) = (file("mapped.sam"), file("ref.fa"), file("zipped.bam"));
    let (sorted, deduped, metrics) = (file("sorted.bam"), file("dedup.bam"), file("dedup.txt"));

    fgumi(
        &["extract", "-i", p(&fq1), p(&fq2), "-r", "4M+T", "+T", "-o", p(&unmapped)]
            .into_iter()
            .chain(["--sample", "s", "--library", "l"])
            .collect::<Vec<_>>(),
    );
    fgumi(&["zipper", "-i", p(&mapped), "-u", p(&unmapped), "-r", p(&reference), "-o", p(&zipped)]);
    fgumi(&["sort", "-i", p(&zipped), "-o", p(&sorted), "--order", "template-coordinate"]);
    fgumi(&["dedup", "-i", p(&sorted), "-o", p(&deduped), "-m", p(&metrics)]);

    (read_bam_output(&sorted).1, read_bam_output(&deduped).1)
}

/// The read names of `records`, in order.
fn names(records: &[RecordBuf]) -> Vec<String> {
    records.iter().map(|r| r.name().expect("named record").to_string()).collect()
}

/// The flag bits that say which alignment of which segment a record is; zipper
/// rewrites mate flags and `dedup` sets the duplicate flag, but neither may
/// change these.
const IDENTITY_FLAGS: u16 = 0x4 | 0x40 | 0x80 | 0x100 | 0x800;

/// Each record as `(name, identity flags, 1-based position)`, sorted, so two
/// record sets compare equal exactly when they hold the same alignments.
fn identities(records: &[RecordBuf]) -> Vec<(String, u16, usize)> {
    let mut ids: Vec<_> = records
        .iter()
        .map(|r| {
            let name = r.name().expect("named record").to_string();
            let pos = r.alignment_start().map_or(0, usize::from);
            (name, r.flags().bits() & IDENTITY_FLAGS, pos)
        })
        .collect();
    ids.sort();
    ids
}

/// Asserts each template's records are contiguous in `names`.
fn assert_templates_contiguous(names: &[String]) {
    let mut seen: Vec<&str> = Vec::new();
    for name in names {
        if seen.last() != Some(&name.as_str()) {
            assert!(!seen.contains(&name.as_str()), "template {name} is split in {names:?}");
            seen.push(name);
        }
    }
}

/// A secondary or supplementary record must sort with its template's primaries,
/// and `dedup` must keep it.
#[rstest]
// #1088: A's supplementary was keyed 1 bp right of A, so B (starting 1 bp right
// of A) sorted between them and `dedup` dropped the supplementary.
#[case::template_one_bp_to_the_right(
    &["A", "B"],
    &[
        ("A", 99, 100, 60, "10M", 301),
        ("A", 147, 301, 60, "10M", 100),
        ("A", 2147, 2001, 60, "10M", 301),
        ("B", 99, 101, 60, "10M", 191),
        ("B", 147, 191, 60, "10M", 101),
    ],
)]
// #1089: C has one mapped primary (R2 unmapped). Its supplementary's `tc` copied
// that primary into both slots while the primary is keyed with "no mate", so the
// supplementary sorted before D and C's primaries after it.
#[case::mate_unmapped_template(
    &["C", "D"],
    &[
        ("C", 73, 1001, 60, "10M", 1001),
        ("C", 133, 1001, 0, "*", 1001),
        ("C", 2121, 2501, 60, "10M", 1001),
        ("D", 99, 1201, 60, "10M", 1401),
        ("D", 147, 1401, 60, "10M", 1201),
    ],
)]
fn test_tc_tagged_record_sorts_with_its_template_and_survives_dedup(
    #[case] template_names: &[&str],
    #[case] alignments: &[(&str, u16, u32, u8, &str, u32)],
) {
    let dir = tempfile::tempdir().unwrap();
    write_inputs(&dir, template_names, alignments);

    let (sorted, deduped) = run_pipeline(dir.path());

    let mut input: Vec<_> = alignments
        .iter()
        .map(|&(name, flag, pos, ..)| (name.to_string(), flag & IDENTITY_FLAGS, pos as usize))
        .collect();
    input.sort();
    assert_eq!(identities(&sorted), input, "sort must keep exactly the input alignments");
    assert_templates_contiguous(&names(&sorted));
    assert_eq!(identities(&deduped), input, "dedup must keep exactly the input alignments");
}
