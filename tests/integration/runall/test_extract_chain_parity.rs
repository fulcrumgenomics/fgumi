//! Real-data parity gate for the `(Extract, Extract)` chain runner.
//!
//! Runs subprocess `fgumi extract` and subprocess
//! `fgumi runall --start-from extract --stop-after extract` against the
//! same paired-end gzip FASTQ fixture, then asserts
//! `fgumi compare bams --command extract` reports the two BAMs as
//! identical record-for-record (headers ignored by the compare preset).

use std::io::Write;
use std::path::Path;
use std::process::Command;

use flate2::Compression;
use flate2::write::GzEncoder;
use tempfile::TempDir;

/// Resolved path to the `fgumi` binary built by Cargo for this test.
fn fgumi_bin() -> &'static str {
    env!("CARGO_BIN_EXE_fgumi")
}

const DNA_BASES: [u8; 4] = [b'A', b'C', b'G', b'T'];

/// Write a paired-end gzip-compressed FASTQ fixture.
///
/// R1 = 8 bp UMI + 16 bp template, R2 = 16 bp template. The UMI is a
/// rotation of `ACGT` derived from the record index so successive reads
/// have varied UMIs, exercising both extraction and BAM record bytes.
fn write_paired_gzip_fastq(r1: &Path, r2: &Path, num_pairs: usize, r1_template_len: usize) {
    let f1 = std::fs::File::create(r1).unwrap();
    let mut e1 = GzEncoder::new(f1, Compression::default());
    let f2 = std::fs::File::create(r2).unwrap();
    let mut e2 = GzEncoder::new(f2, Compression::default());

    for i in 0..num_pairs {
        let umi_bytes: [u8; 8] = std::array::from_fn(|k| DNA_BASES[(i + k) % 4]);
        let umi = std::str::from_utf8(&umi_bytes).unwrap();
        let template_r1 = "CGATCGATCGATCGAT";
        let r1_seq = format!("{umi}{}", &template_r1[..r1_template_len]);
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

/// Subprocess `fgumi extract` (the standalone command).
///
/// `extra_args` are appended to the base command line so individual
/// parity cases can flip flags like `--store-umi-quals` while sharing
/// the same scaffolding.
fn run_standalone_extract(r1: &Path, r2: &Path, out: &Path, extra_args: &[&str]) {
    let mut args: Vec<&str> = vec![
        "extract",
        "--threads",
        "4",
        "--inputs",
        r1.to_str().unwrap(),
        r2.to_str().unwrap(),
        "--output",
        out.to_str().unwrap(),
        "--read-structures",
        "8M+T",
        "+T",
        "--sample",
        "S",
        "--library",
        "L",
        "--compression-level",
        "1",
    ];
    args.extend_from_slice(extra_args);
    let output = Command::new(fgumi_bin()).args(&args).output().expect("spawn fgumi extract");
    assert!(
        output.status.success(),
        "fgumi extract failed (exit {}):\nstdout:\n{}\nstderr:\n{}",
        output.status,
        String::from_utf8_lossy(&output.stdout),
        String::from_utf8_lossy(&output.stderr),
    );
}

/// Subprocess `fgumi runall --start-from extract --stop-after extract`.
///
/// `extra_args` are appended to the base command line so individual
/// parity cases can flip flags like `--extract::store-umi-quals` while
/// sharing the same scaffolding.
fn run_runall_extract_chain(r1: &Path, r2: &Path, out: &Path, extra_args: &[&str]) {
    let mut args: Vec<&str> = vec![
        "runall",
        "--threads",
        "4",
        "--start-from",
        "extract",
        "--stop-after",
        "extract",
        "--input",
        r1.to_str().unwrap(),
        r2.to_str().unwrap(),
        "--extract::sample",
        "S",
        "--extract::library",
        "L",
        "--extract::read-structures",
        "8M+T",
        "+T",
        "--output",
        out.to_str().unwrap(),
        "--compression-level",
        "1",
    ];
    args.extend_from_slice(extra_args);
    let output = Command::new(fgumi_bin()).args(&args).output().expect("spawn fgumi runall");
    assert!(
        output.status.success(),
        "fgumi runall failed (exit {}):\nstdout:\n{}\nstderr:\n{}",
        output.status,
        String::from_utf8_lossy(&output.stdout),
        String::from_utf8_lossy(&output.stderr),
    );
}

/// `fgumi compare bams --command extract` exit status (zero = identical
/// under the extract preset: content mode, headers ignored, in-order).
fn compare_bams_extract(a: &Path, b: &Path) -> std::process::Output {
    Command::new(fgumi_bin())
        .args(["compare", "bams", a.to_str().unwrap(), b.to_str().unwrap(), "--command", "extract"])
        .output()
        .expect("spawn fgumi compare bams")
}

#[test]
fn extract_chain_parity_matches_standalone_on_paired_gzip() {
    let dir = TempDir::new().unwrap();
    let r1 = dir.path().join("R1.fastq.gz");
    let r2 = dir.path().join("R2.fastq.gz");
    let standalone_out = dir.path().join("standalone.bam");
    let runall_out = dir.path().join("runall.bam");

    write_paired_gzip_fastq(&r1, &r2, 100, 16);
    run_standalone_extract(&r1, &r2, &standalone_out, &[]);
    run_runall_extract_chain(&r1, &r2, &runall_out, &[]);

    // Sanity: both runs must produce a non-empty BAM.
    let standalone_size = std::fs::metadata(&standalone_out).unwrap().len();
    let runall_size = std::fs::metadata(&runall_out).unwrap().len();
    assert!(standalone_size > 0, "standalone BAM is empty");
    assert!(runall_size > 0, "runall BAM is empty");

    let cmp = compare_bams_extract(&standalone_out, &runall_out);
    assert!(
        cmp.status.success(),
        "fgumi compare bams reported a mismatch (exit {}):\nstdout:\n{}\nstderr:\n{}",
        cmp.status,
        String::from_utf8_lossy(&cmp.stdout),
        String::from_utf8_lossy(&cmp.stderr),
    );
}

/// Parity gate exercising `--store-umi-quals` (standalone) vs
/// `--extract::store-umi-quals true` (runall). The two BAMs must be
/// byte-identical under the `extract` compare preset, including the
/// added QX aux tag.
#[test]
fn extract_chain_parity_with_store_umi_quals_matches_standalone() {
    let dir = TempDir::new().unwrap();
    let r1 = dir.path().join("R1.fastq.gz");
    let r2 = dir.path().join("R2.fastq.gz");
    let standalone_out = dir.path().join("standalone.bam");
    let runall_out = dir.path().join("runall.bam");

    write_paired_gzip_fastq(&r1, &r2, 100, 16);
    run_standalone_extract(&r1, &r2, &standalone_out, &["--store-umi-quals"]);
    run_runall_extract_chain(&r1, &r2, &runall_out, &["--extract::store-umi-quals", "true"]);

    let cmp = compare_bams_extract(&standalone_out, &runall_out);
    assert!(
        cmp.status.success(),
        "store-umi-quals parity mismatch (exit {}):\nstdout:\n{}\nstderr:\n{}",
        cmp.status,
        String::from_utf8_lossy(&cmp.stdout),
        String::from_utf8_lossy(&cmp.stderr),
    );
}
