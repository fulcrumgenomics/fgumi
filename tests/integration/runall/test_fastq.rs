//! Byte-comparison: runall fastq entry matches `fgumi fastq` output.

use std::io::Write;
use std::path::Path;

fn fgumi_bin() -> &'static str {
    env!("CARGO_BIN_EXE_fgumi")
}

/// Paired-end gzip FASTQ fixture. Same shape as `test_extract::write_paired_gzip_fastq`
/// / `test_correct` but kept local to avoid cross-file coupling. 100 pairs,
/// 8bp UMI + 16bp template.
fn write_paired_gzip_fastq(r1: &Path, r2: &Path, num_pairs: usize) {
    use flate2::{Compression, write::GzEncoder};
    const DNA_BASES: [u8; 4] = [b'A', b'C', b'G', b'T'];
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

/// Runs `fgumi fastq -i <bam>` and captures stdout to `out`.
fn run_v1_fastq(input: &Path, out: &Path) {
    let output = std::process::Command::new(fgumi_bin())
        .args(["fastq", "--input", input.to_str().unwrap()])
        .output()
        .expect("run v1 fgumi fastq");
    assert!(
        output.status.success(),
        "v1 fastq failed: {:?}",
        String::from_utf8_lossy(&output.stderr)
    );
    std::fs::write(out, &output.stdout).unwrap();
}

fn run_v2_fastq(input: &Path, out: &Path) {
    let output = std::process::Command::new(fgumi_bin())
        .args([
            "runall",
            "--start-from",
            "fastq",
            "--stop-after",
            "fastq",
            "--input",
            input.to_str().unwrap(),
            "--output",
            out.to_str().unwrap(),
            "--threads",
            "4",
        ])
        .output()
        .expect("run runall fastq");
    assert!(
        output.status.success(),
        "runall fastq failed: {:?}",
        String::from_utf8_lossy(&output.stderr)
    );
}

#[test]
fn test_fastq_v2_matches_v1() {
    let tmp = tempfile::tempdir().unwrap();
    let r1 = tmp.path().join("R1.fastq.gz");
    let r2 = tmp.path().join("R2.fastq.gz");
    let extracted = tmp.path().join("extracted.bam");
    let v1_fastq = tmp.path().join("v1.fastq");
    let v2_fastq = tmp.path().join("v2.fastq");

    write_paired_gzip_fastq(&r1, &r2, 100);
    run_v1_extract(&r1, &r2, &extracted);
    run_v1_fastq(&extracted, &v1_fastq);
    run_v2_fastq(&extracted, &v2_fastq);

    let v1_bytes = std::fs::read(&v1_fastq).unwrap();
    let v2_bytes = std::fs::read(&v2_fastq).unwrap();
    assert_eq!(
        v1_bytes.len(),
        v2_bytes.len(),
        "FASTQ lengths differ: v1={} v2={}",
        v1_bytes.len(),
        v2_bytes.len()
    );
    assert_eq!(v1_bytes, v2_bytes, "FASTQ bytes differ");
}
