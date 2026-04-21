//! Byte-comparison: runall correction matches v1 `fgumi correct`.

use std::io::Write;
use std::path::Path;

fn fgumi_bin() -> &'static str {
    env!("CARGO_BIN_EXE_fgumi")
}

const DNA_BASES: [u8; 4] = [b'A', b'C', b'G', b'T'];

fn write_paired_gzip_fastq(r1: &Path, r2: &Path, num_pairs: usize) {
    use flate2::{Compression, write::GzEncoder};
    let f1 = std::fs::File::create(r1).unwrap();
    let mut e1 = GzEncoder::new(f1, Compression::default());
    let f2 = std::fs::File::create(r2).unwrap();
    let mut e2 = GzEncoder::new(f2, Compression::default());
    for i in 0..num_pairs {
        let umi_bytes: [u8; 8] = std::array::from_fn(|k| DNA_BASES[(i + k) % 4]);
        let umi = std::str::from_utf8(&umi_bytes).unwrap();
        let r1_seq = format!("{umi}CGATCGATCGATCGAT");
        let r2_seq = "TAGCTAGCTAGCTAGC";
        writeln!(e1, "@read{i}/1").unwrap();
        writeln!(e1, "{r1_seq}").unwrap();
        writeln!(e1, "+").unwrap();
        writeln!(e1, "{}", "I".repeat(r1_seq.len())).unwrap();
        writeln!(e2, "@read{i}/2").unwrap();
        writeln!(e2, "{r2_seq}").unwrap();
        writeln!(e2, "+").unwrap();
        writeln!(e2, "{}", "I".repeat(r2_seq.len())).unwrap();
    }
    e1.finish().unwrap();
    e2.finish().unwrap();
}

/// Write a single-UMI reference file.
///
/// Only one UMI is listed so that reads whose UMI differs by > 1 mismatch
/// from `ACGTACGT` are rejected. With the rotating 4-pattern fixture this
/// yields ~25 kept and ~75 rejected per 100 pairs, guaranteeing that both
/// v1 and v2 write a non-empty rejects BAM — needed for a meaningful
/// byte-comparison of the rejects file.
fn write_umi_reference(path: &Path) {
    let mut f = std::fs::File::create(path).unwrap();
    writeln!(f, "ACGTACGT").unwrap();
}

fn run_v1_extract(r1: &Path, r2: &Path, out: &Path) {
    let output = std::process::Command::new(fgumi_bin())
        .args([
            "extract",
            "--inputs",
            r1.to_str().unwrap(),
            r2.to_str().unwrap(),
            "--output",
            out.to_str().unwrap(),
            "--sample",
            "S",
            "--library",
            "L",
            "--read-structures",
            "8M+T",
            "+T",
            "--compression-level",
            "1",
        ])
        .output()
        .expect("run v1 fgumi extract");
    assert!(
        output.status.success(),
        "v1 extract failed: {:?}",
        String::from_utf8_lossy(&output.stderr)
    );
}

fn run_v1_correct(input: &Path, umi_file: &Path, out: &Path, rejects: &Path, metrics: &Path) {
    let output = std::process::Command::new(fgumi_bin())
        .args([
            "correct",
            "--input",
            input.to_str().unwrap(),
            "--output",
            out.to_str().unwrap(),
            "--rejects",
            rejects.to_str().unwrap(),
            "--metrics",
            metrics.to_str().unwrap(),
            "--umi-files",
            umi_file.to_str().unwrap(),
            "--min-distance",
            "2",
            "--max-mismatches",
            "1",
            "--compression-level",
            "1",
            "--threads",
            "4",
        ])
        .output()
        .expect("run v1 fgumi correct");
    assert!(
        output.status.success(),
        "v1 correct failed: {:?}",
        String::from_utf8_lossy(&output.stderr)
    );
}

fn run_v2_inline(
    r1: &Path,
    r2: &Path,
    umi_file: &Path,
    out: &Path,
    rejects: &Path,
    metrics: &Path,
) {
    let output = std::process::Command::new(fgumi_bin())
        .args([
            "runall",
            "--start-from",
            "extract",
            "--stop-after",
            "correct",
            "--input",
            r1.to_str().unwrap(),
            r2.to_str().unwrap(),
            "--output",
            out.to_str().unwrap(),
            "--extract::sample",
            "S",
            "--extract::library",
            "L",
            "--extract::read-structures",
            "8M+T",
            "+T",
            "--correct::umi-files",
            umi_file.to_str().unwrap(),
            "--correct::max-mismatches",
            "1",
            "--correct::min-distance",
            "2",
            "--correct::rejects",
            rejects.to_str().unwrap(),
            "--correct::metrics",
            metrics.to_str().unwrap(),
            "--compression-level",
            "1",
            "--threads",
            "4",
        ])
        .output()
        .expect("run v2 inline runall");
    assert!(
        output.status.success(),
        "v2 inline runall failed: {:?}",
        String::from_utf8_lossy(&output.stderr)
    );
}

fn run_v2_standalone(input: &Path, umi_file: &Path, out: &Path, rejects: &Path, metrics: &Path) {
    let output = std::process::Command::new(fgumi_bin())
        .args([
            "runall",
            "--start-from",
            "correct",
            "--stop-after",
            "correct",
            "--input",
            input.to_str().unwrap(),
            "--output",
            out.to_str().unwrap(),
            "--correct::umi-files",
            umi_file.to_str().unwrap(),
            "--correct::max-mismatches",
            "1",
            "--correct::min-distance",
            "2",
            "--correct::rejects",
            rejects.to_str().unwrap(),
            "--correct::metrics",
            metrics.to_str().unwrap(),
            "--compression-level",
            "1",
            "--threads",
            "4",
        ])
        .output()
        .expect("run v2 standalone runall");
    assert!(
        output.status.success(),
        "v2 standalone runall failed: {:?}",
        String::from_utf8_lossy(&output.stderr)
    );
}

fn compare_bams(a: &Path, b: &Path) -> bool {
    let output = std::process::Command::new(fgumi_bin())
        .args(["compare", "bams", a.to_str().unwrap(), b.to_str().unwrap(), "--command", "correct"])
        .output()
        .expect("run fgumi compare bams");
    if !output.status.success() {
        eprintln!("compare stdout:\n{}", String::from_utf8_lossy(&output.stdout));
        eprintln!("compare stderr:\n{}", String::from_utf8_lossy(&output.stderr));
    }
    output.status.success()
}

#[test]
fn test_correct_v2_inline_matches_v1() {
    let tmp = tempfile::tempdir().unwrap();
    let r1 = tmp.path().join("R1.fastq.gz");
    let r2 = tmp.path().join("R2.fastq.gz");
    let umi_file = tmp.path().join("umis.txt");
    let v1_extracted = tmp.path().join("v1_extracted.bam");
    let v1_corrected = tmp.path().join("v1_corrected.bam");
    let v1_rejects = tmp.path().join("v1_rejects.bam");
    let v1_metrics = tmp.path().join("v1_metrics.tsv");
    let v2_corrected = tmp.path().join("v2_corrected.bam");
    let v2_rejects = tmp.path().join("v2_rejects.bam");
    let v2_metrics = tmp.path().join("v2_metrics.tsv");

    write_paired_gzip_fastq(&r1, &r2, 100);
    write_umi_reference(&umi_file);

    run_v1_extract(&r1, &r2, &v1_extracted);
    run_v1_correct(&v1_extracted, &umi_file, &v1_corrected, &v1_rejects, &v1_metrics);

    run_v2_inline(&r1, &r2, &umi_file, &v2_corrected, &v2_rejects, &v2_metrics);

    assert!(compare_bams(&v1_corrected, &v2_corrected), "primary BAMs differ");
    assert!(compare_bams(&v1_rejects, &v2_rejects), "rejects BAMs differ");
    assert_eq!(
        std::fs::read(&v1_metrics).unwrap(),
        std::fs::read(&v2_metrics).unwrap(),
        "metrics TSV differs"
    );
}

#[test]
fn test_correct_v2_standalone_matches_v1() {
    let tmp = tempfile::tempdir().unwrap();
    let r1 = tmp.path().join("R1.fastq.gz");
    let r2 = tmp.path().join("R2.fastq.gz");
    let umi_file = tmp.path().join("umis.txt");
    let extracted = tmp.path().join("extracted.bam");
    let v1_corrected = tmp.path().join("v1_corrected.bam");
    let v1_rejects = tmp.path().join("v1_rejects.bam");
    let v1_metrics = tmp.path().join("v1_metrics.tsv");
    let v2_corrected = tmp.path().join("v2_corrected.bam");
    let v2_rejects = tmp.path().join("v2_rejects.bam");
    let v2_metrics = tmp.path().join("v2_metrics.tsv");

    write_paired_gzip_fastq(&r1, &r2, 100);
    write_umi_reference(&umi_file);

    run_v1_extract(&r1, &r2, &extracted);
    run_v1_correct(&extracted, &umi_file, &v1_corrected, &v1_rejects, &v1_metrics);

    run_v2_standalone(&extracted, &umi_file, &v2_corrected, &v2_rejects, &v2_metrics);

    assert!(compare_bams(&v1_corrected, &v2_corrected), "primary BAMs differ");
    assert!(compare_bams(&v1_rejects, &v2_rejects), "rejects BAMs differ");
    assert_eq!(
        std::fs::read(&v1_metrics).unwrap(),
        std::fs::read(&v2_metrics).unwrap(),
        "metrics TSV differs"
    );
}
