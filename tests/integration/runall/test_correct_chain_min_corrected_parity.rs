//! Real-data parity gate for `--correct::min-corrected` on the
//! `(Correct, Correct)` and `(Extract, Correct)` chain runners.
//!
//! Two cases per chain:
//! 1. Low threshold (0.50): both standalone and runall must SUCCEED on a
//!    fixture where ~83% of templates pass correction; we additionally
//!    assert byte-identical output BAMs between standalone and the BAM-
//!    linear chain runner.
//! 2. High threshold (0.99): both must FAIL with an error message
//!    matching the standalone error text byte-for-byte (the runall
//!    runners replicate the exact `bail!()` string from
//!    `run_correct_pipeline`).

use std::io::Write;
use std::path::Path;
use std::process::{Command, Output};

use flate2::Compression;
use flate2::write::GzEncoder;
use tempfile::TempDir;

fn fgumi_bin() -> &'static str {
    env!("CARGO_BIN_EXE_fgumi")
}

const ALLOWED_UMIS: &[&str] = &["AAAAAAAA", "CCCCCCCC", "GGGGGGGG", "TTTTTTTT"];

/// 6-UMI cycle: 5 allowed (4 perfect + 1 single-mismatch off `AAAAAAAA`)
/// and 1 unmatchable `NNNNNNNN`. With `--max-mismatches 1
/// --min-distance 1`, 5 of every 6 templates are kept (~83%) and 1 is
/// dropped — the kept fraction sits well above 0.50 (low threshold) and
/// well below 0.99 (high threshold).
const UMI_CYCLE: &[&[u8]] =
    &[b"AAAAAAAA", b"CCCCCCCC", b"GGGGGGGG", b"TTTTTTTT", b"AAAAAAAT", b"NNNNNNNN"];

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

fn run_extract_to_bam(r1: &Path, r2: &Path, out: &Path) {
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

fn run_standalone_correct(input: &Path, out: &Path, min_corrected: f64) -> Output {
    let mut args: Vec<String> = vec![
        "correct".into(),
        "--threads".into(),
        "4".into(),
        "--input".into(),
        input.to_str().unwrap().into(),
        "--output".into(),
        out.to_str().unwrap().into(),
        "--min-corrected".into(),
        format!("{min_corrected}"),
        "--max-mismatches".into(),
        "1".into(),
        "--min-distance".into(),
        "1".into(),
        "--compression-level".into(),
        "1".into(),
    ];
    for umi in ALLOWED_UMIS {
        args.push("--umis".into());
        args.push((*umi).to_string());
    }
    Command::new(fgumi_bin()).args(&args).output().expect("spawn fgumi correct")
}

fn run_runall_correct_chain(input: &Path, out: &Path, min_corrected: f64) -> Output {
    let mut args: Vec<String> = vec![
        "runall".into(),
        "--threads".into(),
        "4".into(),
        "--start-from".into(),
        "correct".into(),
        "--stop-after".into(),
        "correct".into(),
        "--input".into(),
        input.to_str().unwrap().into(),
        "--output".into(),
        out.to_str().unwrap().into(),
        "--compression-level".into(),
        "1".into(),
        "--correct::min-corrected".into(),
        format!("{min_corrected}"),
        "--correct::max-mismatches".into(),
        "1".into(),
        "--correct::min-distance".into(),
        "1".into(),
        "--correct::umis".into(),
    ];
    args.extend(ALLOWED_UMIS.iter().map(|s| (*s).to_string()));
    Command::new(fgumi_bin()).args(&args).output().expect("spawn fgumi runall")
}

fn run_runall_extract_correct_chain(
    r1: &Path,
    r2: &Path,
    out: &Path,
    min_corrected: f64,
) -> Output {
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
        "--output".into(),
        out.to_str().unwrap().into(),
        "--compression-level".into(),
        "1".into(),
        "--extract::sample".into(),
        "S".into(),
        "--extract::library".into(),
        "L".into(),
        "--extract::read-structures".into(),
        "8M+T".into(),
        "+T".into(),
        "--correct::min-corrected".into(),
        format!("{min_corrected}"),
        "--correct::max-mismatches".into(),
        "1".into(),
        "--correct::min-distance".into(),
        "1".into(),
        "--correct::umis".into(),
    ];
    args.extend(ALLOWED_UMIS.iter().map(|s| (*s).to_string()));
    Command::new(fgumi_bin()).args(&args).output().expect("spawn fgumi runall")
}

fn compare_bams_correct(a: &Path, b: &Path) -> Output {
    Command::new(fgumi_bin())
        .args(["compare", "bams", a.to_str().unwrap(), b.to_str().unwrap(), "--command", "correct"])
        .output()
        .expect("spawn fgumi compare bams")
}

#[test]
fn correct_chain_min_corrected_low_threshold_succeeds_with_bam_parity() {
    let dir = TempDir::new().unwrap();
    let r1 = dir.path().join("R1.fastq.gz");
    let r2 = dir.path().join("R2.fastq.gz");
    let extracted = dir.path().join("extracted.bam");
    let standalone_out = dir.path().join("standalone.bam");
    let runall_out = dir.path().join("runall.bam");

    write_paired_gzip_fastq_with_cycling_umi(&r1, &r2, 120);
    run_extract_to_bam(&r1, &r2, &extracted);

    let standalone = run_standalone_correct(&extracted, &standalone_out, 0.50);
    assert!(
        standalone.status.success(),
        "standalone correct must succeed at min-corrected=0.50; stderr:\n{}",
        String::from_utf8_lossy(&standalone.stderr),
    );
    let runall = run_runall_correct_chain(&extracted, &runall_out, 0.50);
    assert!(
        runall.status.success(),
        "runall correct chain must succeed at --correct::min-corrected=0.50; stderr:\n{}",
        String::from_utf8_lossy(&runall.stderr),
    );

    let cmp = compare_bams_correct(&standalone_out, &runall_out);
    assert!(
        cmp.status.success(),
        "fgumi compare bams reported a mismatch (exit {}):\nstdout:\n{}\nstderr:\n{}",
        cmp.status,
        String::from_utf8_lossy(&cmp.stdout),
        String::from_utf8_lossy(&cmp.stderr),
    );
}

#[test]
fn correct_chain_min_corrected_high_threshold_fails_with_matching_error() {
    let dir = TempDir::new().unwrap();
    let r1 = dir.path().join("R1.fastq.gz");
    let r2 = dir.path().join("R2.fastq.gz");
    let extracted = dir.path().join("extracted.bam");
    let standalone_out = dir.path().join("standalone.bam");
    let runall_out = dir.path().join("runall.bam");

    write_paired_gzip_fastq_with_cycling_umi(&r1, &r2, 120);
    run_extract_to_bam(&r1, &r2, &extracted);

    let standalone = run_standalone_correct(&extracted, &standalone_out, 0.99);
    assert!(
        !standalone.status.success(),
        "standalone correct must fail at min-corrected=0.99; stderr:\n{}",
        String::from_utf8_lossy(&standalone.stderr),
    );
    let runall = run_runall_correct_chain(&extracted, &runall_out, 0.99);
    assert!(
        !runall.status.success(),
        "runall correct chain must fail at --correct::min-corrected=0.99; stderr:\n{}",
        String::from_utf8_lossy(&runall.stderr),
    );

    // Both runners route through `run_correct_pipeline`, so the
    // user-facing error text must be byte-identical. Look for the
    // signature substring in both stderr streams.
    let needle = "Final ratio of reads kept / total was";
    let s_err = String::from_utf8_lossy(&standalone.stderr);
    let r_err = String::from_utf8_lossy(&runall.stderr);
    assert!(s_err.contains(needle), "standalone stderr missing min-corrected error: {s_err}");
    assert!(r_err.contains(needle), "runall stderr missing min-corrected error: {r_err}");
}

#[test]
fn extract_correct_chain_min_corrected_low_threshold_succeeds() {
    let dir = TempDir::new().unwrap();
    let r1 = dir.path().join("R1.fastq.gz");
    let r2 = dir.path().join("R2.fastq.gz");
    let fused_out = dir.path().join("fused.bam");

    write_paired_gzip_fastq_with_cycling_umi(&r1, &r2, 120);

    let fused = run_runall_extract_correct_chain(&r1, &r2, &fused_out, 0.50);
    assert!(
        fused.status.success(),
        "fused extract→correct chain must succeed at --correct::min-corrected=0.50; stderr:\n{}",
        String::from_utf8_lossy(&fused.stderr),
    );
}

#[test]
fn extract_correct_chain_min_corrected_high_threshold_fails_with_matching_error() {
    let dir = TempDir::new().unwrap();
    let r1 = dir.path().join("R1.fastq.gz");
    let r2 = dir.path().join("R2.fastq.gz");
    let extracted = dir.path().join("extracted.bam");
    let standalone_out = dir.path().join("standalone.bam");
    let fused_out = dir.path().join("fused.bam");

    write_paired_gzip_fastq_with_cycling_umi(&r1, &r2, 120);
    run_extract_to_bam(&r1, &r2, &extracted);

    let standalone = run_standalone_correct(&extracted, &standalone_out, 0.99);
    assert!(
        !standalone.status.success(),
        "standalone correct must fail at min-corrected=0.99; stderr:\n{}",
        String::from_utf8_lossy(&standalone.stderr),
    );
    let fused = run_runall_extract_correct_chain(&r1, &r2, &fused_out, 0.99);
    assert!(
        !fused.status.success(),
        "fused extract→correct chain must fail at --correct::min-corrected=0.99; stderr:\n{}",
        String::from_utf8_lossy(&fused.stderr),
    );

    let needle = "Final ratio of reads kept / total was";
    let s_err = String::from_utf8_lossy(&standalone.stderr);
    let f_err = String::from_utf8_lossy(&fused.stderr);
    assert!(s_err.contains(needle), "standalone stderr missing min-corrected error: {s_err}");
    assert!(f_err.contains(needle), "fused chain stderr missing min-corrected error: {f_err}");
}
