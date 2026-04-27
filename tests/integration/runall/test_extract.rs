//! Byte-comparison: runall --start-from extract matches v1 `fgumi extract`.
//!
//! Tests subprocess `fgumi extract` (v1) vs `fgumi runall --start-from extract --stop-after extract`
//! (v2) on synthetic paired-end FASTQ fixtures. Uses `fgumi compare bams --command extract`
//! for record-by-record comparison (headers are ignored).

use std::io::Write;
use std::path::Path;

/// Binary path, resolved by cargo at test compile time.
fn fgumi_bin() -> &'static str {
    env!("CARGO_BIN_EXE_fgumi")
}

/// DNA alphabet used to build deterministic synthetic UMI sequences.
const DNA_BASES: [u8; 4] = [b'A', b'C', b'G', b'T'];

/// Write a paired-end gzip-compressed FASTQ fixture with a known UMI structure.
/// Record structure: R1 = 8bp UMI + 16bp template, R2 = 16bp template.
fn write_paired_gzip_fastq(r1: &Path, r2: &Path, num_pairs: usize) {
    use flate2::{Compression, write::GzEncoder};

    let f1 = std::fs::File::create(r1).unwrap();
    let mut e1 = GzEncoder::new(f1, Compression::default());
    let f2 = std::fs::File::create(r2).unwrap();
    let mut e2 = GzEncoder::new(f2, Compression::default());

    for i in 0..num_pairs {
        // 8bp UMI from a fixed DNA alphabet so every base is valid regardless of i.
        let umi_bytes: [u8; 8] = std::array::from_fn(|k| DNA_BASES[(i + k) % 4]);
        let umi = std::str::from_utf8(&umi_bytes).unwrap();
        let r1_seq = format!("{umi}CGATCGATCGATCGAT"); // 8 UMI + 16 template
        let r2_seq = "TAGCTAGCTAGCTAGC"; // 16 template

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

/// Run `fgumi extract` (v1 path) on paired FASTQ inputs.
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
        "v1 fgumi extract failed (exit {}):\nstdout:\n{}\nstderr:\n{}",
        output.status,
        String::from_utf8_lossy(&output.stdout),
        String::from_utf8_lossy(&output.stderr),
    );
}

/// Run `fgumi runall --start-from extract --stop-after extract` (v2 path).
fn run_v2_extract(r1: &Path, r2: &Path, out: &Path) {
    let output = std::process::Command::new(fgumi_bin())
        .args([
            "runall",
            "--start-from",
            "extract",
            "--stop-after",
            "extract",
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
            "--compression-level",
            "1",
        ])
        .output()
        .expect("run v2 fgumi runall");
    assert!(
        output.status.success(),
        "v2 fgumi runall failed (exit {}):\nstdout:\n{}\nstderr:\n{}",
        output.status,
        String::from_utf8_lossy(&output.stdout),
        String::from_utf8_lossy(&output.stderr),
    );
}

/// Invoke `fgumi compare bams --command extract` as a subprocess.
/// Returns true on exit status 0 (BAMs equivalent).
fn compare_bams(a: &Path, b: &Path) -> bool {
    let output = std::process::Command::new(fgumi_bin())
        .args(["compare", "bams", a.to_str().unwrap(), b.to_str().unwrap(), "--command", "extract"])
        .output()
        .expect("run fgumi compare bams");
    if !output.status.success() {
        eprintln!("fgumi compare bams stdout:\n{}", String::from_utf8_lossy(&output.stdout));
        eprintln!("fgumi compare bams stderr:\n{}", String::from_utf8_lossy(&output.stderr));
    }
    output.status.success()
}

/// Byte-comparison gate: `fgumi runall --start-from extract --stop-after extract`
/// must produce records identical to `fgumi extract` on the same paired-end gzip
/// FASTQ input.
///
/// `--mode content` compares all BAM record fields and tag values; BAM headers
/// (including the @PG line, which differs because the two commands have
/// different command lines) are ignored. Record order must match because
/// extract is deterministic.
#[test]
fn test_extract_v2_matches_v1_paired_gzip() {
    let tmp = tempfile::tempdir().unwrap();
    let r1 = tmp.path().join("R1.fastq.gz");
    let r2 = tmp.path().join("R2.fastq.gz");
    let v1_out = tmp.path().join("v1.bam");
    let v2_out = tmp.path().join("v2.bam");

    write_paired_gzip_fastq(&r1, &r2, 100);
    run_v1_extract(&r1, &r2, &v1_out);
    run_v2_extract(&r1, &r2, &v2_out);

    assert!(
        compare_bams(&v1_out, &v2_out),
        "runall --start-from extract output must match v1 fgumi extract record-for-record"
    );
}
