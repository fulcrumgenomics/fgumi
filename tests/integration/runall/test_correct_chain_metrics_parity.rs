//! Real-data parity gate for `--correct::metrics` on the `(Correct, Correct)`
//! and `(Extract, Correct)` chain runners.
//!
//! Builds a tiny BAM (or paired FASTQ for the fused chain) with a realistic
//! UMI distribution: some UMIs are perfect matches, some are 8-mers within
//! the per-segment edit distance budget that need correction, and some are
//! N-only sentinels that fall outside the allowed list. We then run
//! standalone `fgumi correct --metrics standalone.tsv` and
//! `fgumi runall --correct::metrics runall.tsv` on the same input and assert
//! the two TSVs are byte-for-byte identical (after sorting lines, since the
//! standalone aggregation order isn't guaranteed across threads).

use std::collections::BTreeSet;
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

/// Allowed UMI list. The synthetic fixture cycles through these plus
/// `NNNNNNNN` and `AAAAAAAT` (a 1-mismatch off `AAAAAAAA`); with
/// `--max-mismatches 1 --min-distance 1` the latter exercises the
/// 1-mismatch correction bucket and the former exercises the rejection
/// bucket.
const ALLOWED_UMIS: &[&str] = &["AAAAAAAA", "CCCCCCCC", "GGGGGGGG", "TTTTTTTT"];

/// UMI cycle: 4 perfect matches, 1 single-mismatch, 1 unmatchable.
/// Layout chosen so the rounded counts are predictable per six-record block.
const UMI_CYCLE: &[&[u8]] =
    &[b"AAAAAAAA", b"CCCCCCCC", b"GGGGGGGG", b"TTTTTTTT", b"AAAAAAAT", b"NNNNNNNN"];

/// Write a paired-end gzip-compressed FASTQ fixture.
///
/// R1 = 8 bp UMI (cycled) + 8 bp template, R2 = 8 bp template.
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

/// Subprocess `fgumi extract` to produce an input BAM with RX tags.
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

/// Subprocess `fgumi correct --metrics ...` (standalone).
fn run_standalone_correct_metrics(input: &Path, out: &Path, metrics: &Path) {
    let mut args: Vec<String> = vec![
        "correct".into(),
        "--threads".into(),
        "4".into(),
        "--input".into(),
        input.to_str().unwrap().into(),
        "--output".into(),
        out.to_str().unwrap().into(),
        "--metrics".into(),
        metrics.to_str().unwrap().into(),
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
    let output = Command::new(fgumi_bin()).args(&args).output().expect("spawn fgumi correct");
    assert!(
        output.status.success(),
        "fgumi correct failed (exit {}):\nstdout:\n{}\nstderr:\n{}",
        output.status,
        String::from_utf8_lossy(&output.stdout),
        String::from_utf8_lossy(&output.stderr),
    );
}

/// Subprocess `fgumi runall --start-from correct --stop-after correct
/// --correct::metrics ...`.
fn run_runall_correct_chain_metrics(input: &Path, out: &Path, metrics: &Path) {
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
        "--correct::metrics".into(),
        metrics.to_str().unwrap().into(),
        "--correct::max-mismatches".into(),
        "1".into(),
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

/// Subprocess `fgumi runall --start-from extract --stop-after correct
/// --correct::metrics ...` (the fused chain).
fn run_runall_extract_correct_chain_metrics(r1: &Path, r2: &Path, out: &Path, metrics: &Path) {
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
        "--correct::metrics".into(),
        metrics.to_str().unwrap().into(),
        "--correct::max-mismatches".into(),
        "1".into(),
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

/// Read a TSV file and return its lines, preserving the header (first
/// line) and sorting the data rows lexicographically. Standalone
/// correct's per-UMI rows have a stable sort by UMI in
/// `finalize_correct_metrics`, but we sort defensively in case any
/// future change makes the order thread-dependent.
fn normalize_tsv(path: &Path) -> (String, BTreeSet<String>) {
    let text = std::fs::read_to_string(path)
        .unwrap_or_else(|e| panic!("read TSV {} failed: {e}", path.display()));
    let mut lines = text.lines();
    let header = lines.next().unwrap_or_default().to_string();
    let rows: BTreeSet<String> = lines.map(str::to_string).collect();
    (header, rows)
}

#[test]
fn correct_chain_metrics_parity_matches_standalone() {
    let dir = TempDir::new().unwrap();
    let r1 = dir.path().join("R1.fastq.gz");
    let r2 = dir.path().join("R2.fastq.gz");
    let extracted = dir.path().join("extracted.bam");

    let standalone_out = dir.path().join("standalone.bam");
    let standalone_metrics = dir.path().join("standalone-metrics.tsv");
    let runall_out = dir.path().join("runall.bam");
    let runall_metrics = dir.path().join("runall-metrics.tsv");

    // 120 templates → exactly 20 templates per UMI in the cycle.
    write_paired_gzip_fastq_with_cycling_umi(&r1, &r2, 120);
    run_extract_to_bam(&r1, &r2, &extracted);

    run_standalone_correct_metrics(&extracted, &standalone_out, &standalone_metrics);
    run_runall_correct_chain_metrics(&extracted, &runall_out, &runall_metrics);

    let (s_header, s_rows) = normalize_tsv(&standalone_metrics);
    let (r_header, r_rows) = normalize_tsv(&runall_metrics);

    assert_eq!(s_header, r_header, "TSV header rows must match");
    assert_eq!(
        s_rows, r_rows,
        "per-UMI metrics rows must match exactly between standalone and runall"
    );
    // Sanity: at minimum we expect one row per allowed UMI plus the
    // unmatched-UMI sentinel, so 5 rows.
    assert!(
        s_rows.len() > ALLOWED_UMIS.len(),
        "standalone TSV has unexpectedly few rows: {}",
        s_rows.len()
    );
}

#[test]
fn extract_correct_chain_metrics_parity_matches_standalone() {
    // Same input, but compare the fused extract→correct chain's TSV
    // against the standalone-correct TSV produced from the equivalent
    // pre-extracted BAM. The two pipelines do the same per-template
    // correction work, so the per-UMI counts must agree.
    let dir = TempDir::new().unwrap();
    let r1 = dir.path().join("R1.fastq.gz");
    let r2 = dir.path().join("R2.fastq.gz");
    let extracted = dir.path().join("extracted.bam");

    let standalone_out = dir.path().join("standalone.bam");
    let standalone_metrics = dir.path().join("standalone-metrics.tsv");
    let fused_out = dir.path().join("fused.bam");
    let fused_metrics = dir.path().join("fused-metrics.tsv");

    write_paired_gzip_fastq_with_cycling_umi(&r1, &r2, 120);
    run_extract_to_bam(&r1, &r2, &extracted);

    run_standalone_correct_metrics(&extracted, &standalone_out, &standalone_metrics);
    run_runall_extract_correct_chain_metrics(&r1, &r2, &fused_out, &fused_metrics);

    let (s_header, s_rows) = normalize_tsv(&standalone_metrics);
    let (f_header, f_rows) = normalize_tsv(&fused_metrics);

    assert_eq!(s_header, f_header, "TSV header rows must match");
    assert_eq!(
        s_rows, f_rows,
        "per-UMI metrics rows must match exactly between standalone correct and the \
         fused extract→correct chain"
    );
}
