//! Real-data parity gate for the `(Extract, Correct)` fused chain runner.
//!
//! Runs subprocess `fgumi extract` followed by `fgumi correct` (the
//! pipechain materialised through an intermediate BAM, run sequentially)
//! and subprocess
//! `fgumi runall --start-from extract --stop-after correct` on the same
//! synthetic paired-end gzip FASTQ fixture, then asserts
//! `fgumi compare bams --command correct` reports the two BAMs as
//! identical record-for-record.

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

/// Cycle of UMIs used by the synthetic fixture: four ACGT-homopolymer UMIs
/// (in the allowed list) plus an N8 UMI that lives outside the allowed
/// list. With `--max-mismatches 0 --min-distance 1`, the N8 templates are
/// dropped, exercising the rejection path; the others pass through as
/// perfect matches, exercising the kept path.
const UMI_CYCLE: &[&[u8]] = &[b"AAAAAAAA", b"CCCCCCCC", b"GGGGGGGG", b"TTTTTTTT", b"NNNNNNNN"];

/// Allowed UMI list used by both standalone and runall.
const ALLOWED_UMIS: &[&str] = &["AAAAAAAA", "CCCCCCCC", "GGGGGGGG", "TTTTTTTT"];

/// Write a paired-end gzip-compressed FASTQ fixture (same shape as the
/// correct-chain parity test). 80% allowed UMIs, 20% N8.
fn write_paired_gzip_fastq_with_cycling_umi(r1: &Path, r2: &Path, num_pairs: usize) {
    let f1 = std::fs::File::create(r1).unwrap();
    let mut e1 = GzEncoder::new(f1, Compression::default());
    let f2 = std::fs::File::create(r2).unwrap();
    let mut e2 = GzEncoder::new(f2, Compression::default());

    for i in 0..num_pairs {
        let umi = UMI_CYCLE[i % UMI_CYCLE.len()];
        let umi_str = std::str::from_utf8(umi).unwrap();
        let r1_seq = format!("{umi_str}AAAAAAAA");
        let r2_seq = "CCCCCCCC";

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

/// Subprocess `fgumi extract` to materialise the intermediate BAM.
fn run_extract(r1: &Path, r2: &Path, out: &Path) {
    let output = Command::new(fgumi_bin())
        .args([
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
        ])
        .output()
        .expect("spawn fgumi extract");
    assert!(
        output.status.success(),
        "fgumi extract failed (exit {}):\nstdout:\n{}\nstderr:\n{}",
        output.status,
        String::from_utf8_lossy(&output.stdout),
        String::from_utf8_lossy(&output.stderr),
    );
}

/// Subprocess `fgumi correct` (the standalone command) on `input`.
///
/// Standalone `fgumi correct --umis` takes a single value per occurrence, so
/// we repeat the flag per allowed UMI rather than passing multiple values to
/// one flag.
fn run_standalone_correct(input: &Path, out: &Path, max_mismatches: usize) {
    let mut args: Vec<String> = vec![
        "correct".into(),
        "--threads".into(),
        "4".into(),
        "--input".into(),
        input.to_str().unwrap().into(),
        "--output".into(),
        out.to_str().unwrap().into(),
        "--max-mismatches".into(),
        max_mismatches.to_string(),
        "--min-distance".into(),
        "1".into(),
        "--compression-level".into(),
        "1".into(),
    ];
    for umi in ALLOWED_UMIS {
        args.push("--umis".into());
        args.push((*umi).to_string());
    }
    let output = Command::new(fgumi_bin()).args(&args).output().expect("spawn fgumi correct");
    assert!(
        output.status.success(),
        "fgumi correct failed (exit {}):\nstdout:\n{}\nstderr:\n{}",
        output.status,
        String::from_utf8_lossy(&output.stdout),
        String::from_utf8_lossy(&output.stderr),
    );
}

/// Subprocess `fgumi runall --start-from extract --stop-after correct`.
fn run_runall_extract_correct_chain(r1: &Path, r2: &Path, out: &Path, max_mismatches: usize) {
    let mut args: Vec<String> = vec![
        "runall".into(),
        "--threads".into(),
        "4".into(),
        "--start-from".into(),
        "extract".into(),
        "--stop-after".into(),
        "correct".into(),
        "--input".into(),
        r1.to_str().unwrap().into(),
        r2.to_str().unwrap().into(),
        "--extract::sample".into(),
        "S".into(),
        "--extract::library".into(),
        "L".into(),
        "--extract::read-structures".into(),
        "8M+T".into(),
        "+T".into(),
        "--output".into(),
        out.to_str().unwrap().into(),
        "--compression-level".into(),
        "1".into(),
        "--correct::max-mismatches".into(),
        max_mismatches.to_string(),
        "--correct::min-distance".into(),
        "1".into(),
        "--correct::umis".into(),
    ];
    args.extend(ALLOWED_UMIS.iter().map(|s| (*s).to_string()));
    let output = Command::new(fgumi_bin()).args(&args).output().expect("spawn fgumi runall");
    assert!(
        output.status.success(),
        "fgumi runall failed (exit {}):\nstdout:\n{}\nstderr:\n{}",
        output.status,
        String::from_utf8_lossy(&output.stdout),
        String::from_utf8_lossy(&output.stderr),
    );
}

/// `fgumi compare bams --command correct` exit status (zero = identical).
fn compare_bams_correct(a: &Path, b: &Path) -> std::process::Output {
    Command::new(fgumi_bin())
        .args(["compare", "bams", a.to_str().unwrap(), b.to_str().unwrap(), "--command", "correct"])
        .output()
        .expect("spawn fgumi compare bams")
}

#[test]
fn extract_correct_chain_parity_matches_pipechain_on_paired_gzip() {
    let dir = TempDir::new().unwrap();
    let r1 = dir.path().join("R1.fastq.gz");
    let r2 = dir.path().join("R2.fastq.gz");
    let extracted = dir.path().join("extracted.bam");
    let pipechain_out = dir.path().join("pipechain.bam");
    let runall_out = dir.path().join("runall.bam");

    write_paired_gzip_fastq_with_cycling_umi(&r1, &r2, 100);

    // Pipechain: extract → correct, sequentially via intermediate BAM.
    run_extract(&r1, &r2, &extracted);
    run_standalone_correct(&extracted, &pipechain_out, 0);

    // In-process: runall fused chain.
    run_runall_extract_correct_chain(&r1, &r2, &runall_out, 0);

    // Sanity: both runs must produce a non-empty BAM.
    let pipechain_size = std::fs::metadata(&pipechain_out).unwrap().len();
    let runall_size = std::fs::metadata(&runall_out).unwrap().len();
    assert!(pipechain_size > 0, "pipechain BAM is empty");
    assert!(runall_size > 0, "runall BAM is empty");

    let cmp = compare_bams_correct(&pipechain_out, &runall_out);
    assert!(
        cmp.status.success(),
        "fgumi compare bams reported a mismatch (exit {}):\nstdout:\n{}\nstderr:\n{}",
        cmp.status,
        String::from_utf8_lossy(&cmp.stdout),
        String::from_utf8_lossy(&cmp.stderr),
    );
}
