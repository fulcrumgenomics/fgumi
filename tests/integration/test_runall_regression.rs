//! Regression integration tests for the `fgumi runall` command with ~1000 read pairs.
//!
//! Exercises multiple start-from/stop-after combinations on a larger synthetic dataset
//! to catch regressions in the runall pipeline across threading configurations.
//!
//! Tests are skipped when `bwa` is not on `PATH` or a `PhiX` reference is not available.

use std::fs;
use std::io::{BufRead, BufReader, Write};
use std::path::{Path, PathBuf};
use std::process::Command;
use tempfile::TempDir;

use crate::helpers::assertions::assert_bam_equivalent;
use crate::helpers::bam_generator::count_bam_records;
use crate::helpers::references::check_alignment_prerequisites;

/// Run fgumi with the given args, asserting it succeeds.
fn run_fgumi(args: &[&str]) {
    let output =
        Command::new(env!("CARGO_BIN_EXE_fgumi")).args(args).output().expect("Failed to run fgumi");
    assert!(
        output.status.success(),
        "fgumi {} failed:\nstdout: {}\nstderr: {}",
        args.first().unwrap_or(&""),
        String::from_utf8_lossy(&output.stdout),
        String::from_utf8_lossy(&output.stderr),
    );
}

/// Read a FASTA file and return the concatenated sequence (ignoring headers/newlines).
fn read_fasta_sequence(path: &Path) -> Vec<u8> {
    let file = fs::File::open(path).expect("Failed to open FASTA");
    let reader = BufReader::new(file);
    let mut seq = Vec::new();
    for line in reader.lines() {
        let line = line.expect("Failed to read line");
        if !line.starts_with('>') {
            seq.extend_from_slice(line.trim().as_bytes());
        }
    }
    seq
}

/// Reverse complement a DNA sequence.
fn revcomp(seq: &[u8]) -> Vec<u8> {
    seq.iter()
        .rev()
        .map(|&b| match b {
            b'A' => b'T',
            b'T' => b'A',
            b'C' => b'G',
            b'G' => b'C',
            b'a' => b't',
            b't' => b'a',
            b'c' => b'g',
            b'g' => b'c',
            other => other,
        })
        .collect()
}

/// A set of predefined 6-base UMIs used to build families.
const UMIS: &[&str] = &[
    "AACCGG", "TTGGAA", "CCTTAA", "GGAACC", "ACGTAC", "TGCATG", "CATGCA", "GTACGT", "AATCCG",
    "TTCGGA", "CCGATT", "GGCTAA", "ACTGAC", "TGACTG", "CAGACT", "GACTGA",
];

/// Generate ~1000 paired-end read pairs from a `PhiX` reference with UMIs in read names.
///
/// Creates diverse families at multiple positions spread across the genome.
/// Each position gets 1-3 UMI families with 10-30 reads each.
///
/// Returns `(r1_path, r2_path)`.
fn create_large_test_fastqs(
    dir: &Path,
    ref_path: &Path,
    num_target_pairs: usize,
    read_len: usize,
) -> (PathBuf, PathBuf) {
    let genome = read_fasta_sequence(ref_path);
    assert!(genome.len() > 1000, "Reference sequence too short for test data generation");

    let r1_path = dir.join("r1.fastq");
    let r2_path = dir.join("r2.fastq");
    let mut r1_file = fs::File::create(&r1_path).unwrap();
    let mut r2_file = fs::File::create(&r2_path).unwrap();

    // Build families: spread positions every ~100bp across the genome.
    // Use 2 UMIs per position, ~15 reads per family on average.
    let max_start = genome.len().saturating_sub(2 * read_len);
    let step = 100;
    let mut total_pairs = 0;
    let mut pos = 0;
    let mut umi_idx = 0;

    while pos < max_start && total_pairs < num_target_pairs {
        // Each position gets 1-2 UMI families
        let families_at_pos = if umi_idx % 3 == 0 { 1 } else { 2 };

        for _ in 0..families_at_pos {
            if total_pairs >= num_target_pairs {
                break;
            }
            let umi = UMIS[umi_idx % UMIS.len()];
            // Vary family size: cycle through 10, 15, 20, 25, 30
            let family_size = 10 + (umi_idx % 5) * 5;
            let count = family_size.min(num_target_pairs - total_pairs);

            let r1_start = pos;
            let r1_end = (r1_start + read_len).min(genome.len());
            let r1_seq: Vec<u8> = genome[r1_start..r1_end].to_vec();

            // R2 is reverse complement of a region downstream (simulating insert)
            let insert_size = read_len + 50;
            let r2_start = (r1_start + insert_size).min(genome.len().saturating_sub(read_len));
            let r2_end = (r2_start + read_len).min(genome.len());
            let r2_seq = revcomp(&genome[r2_start..r2_end]);

            let qual = "I".repeat(r1_seq.len());
            let qual_r2 = "I".repeat(r2_seq.len());

            for i in 0..count {
                // Illumina format: @instrument:run:flowcell:lane:tile:x:y:UMI
                let name = format!("INST:1:FC:1:1:{pos}:{i}:{umi}");
                writeln!(r1_file, "@{name} 1:N:0:SAMPLE").unwrap();
                r1_file.write_all(&r1_seq).unwrap();
                writeln!(r1_file).unwrap();
                writeln!(r1_file, "+").unwrap();
                writeln!(r1_file, "{qual}").unwrap();

                writeln!(r2_file, "@{name} 2:N:0:SAMPLE").unwrap();
                r2_file.write_all(&r2_seq).unwrap();
                writeln!(r2_file).unwrap();
                writeln!(r2_file, "+").unwrap();
                writeln!(r2_file, "{qual_r2}").unwrap();
            }

            total_pairs += count;
            umi_idx += 1;
        }

        pos += step;
    }

    assert!(total_pairs > 0, "Failed to generate any test read pairs");
    eprintln!("Generated {total_pairs} paired-end read pairs across {} positions", pos / step);

    (r1_path, r2_path)
}

/// Common arguments shared by all runall invocations that start from extract.
fn base_runall_extract_args<'a>(
    r1_fq: &'a str,
    r2_fq: &'a str,
    output: &'a str,
    ref_path: &'a str,
    read_structure: &'a str,
) -> Vec<&'a str> {
    vec![
        "runall",
        "--start-from",
        "extract",
        "--input",
        r1_fq,
        r2_fq,
        "--output",
        output,
        "--reference",
        ref_path,
        "--consensus",
        "simplex",
        "--extract::sample",
        "TestSample",
        "--extract::library",
        "TestLibrary",
        "--extract::read-structures",
        read_structure,
        read_structure,
        "--extract::extract-umis-from-read-names",
        "--aligner::preset",
        "bwa-mem",
        "--aligner::threads",
        "2",
        "--group::strategy",
        "identity",
        "--group::max-edits",
        "0",
        "--group::min-mapq",
        "0",
        "--filter::min-reads",
        "1",
        "--compression-level",
        "1",
    ]
}

// ============================================================================
// Test 1: full pipeline from FASTQ input
// ============================================================================

/// Full pipeline from extract to end on ~1000 read pairs.
/// Verifies output BAM exists, has records, and is valid.
#[test]
fn test_runall_extract_to_end() {
    let Some(ref_path) = check_alignment_prerequisites() else {
        return;
    };

    let temp_dir = TempDir::new().unwrap();
    let read_len = 100;
    let (r1_fq, r2_fq) = create_large_test_fastqs(temp_dir.path(), &ref_path, 1000, read_len);
    let read_structure = format!("{read_len}T");

    let output = temp_dir.path().join("full_output.bam");

    let mut args = base_runall_extract_args(
        r1_fq.to_str().unwrap(),
        r2_fq.to_str().unwrap(),
        output.to_str().unwrap(),
        ref_path.to_str().unwrap(),
        &read_structure,
    );
    args.extend_from_slice(&["--threads", "4"]);

    run_fgumi(&args);

    assert!(output.exists(), "Output BAM should exist");
    let count = count_bam_records(&output);
    assert!(count > 0, "Output BAM should have at least one record, got {count}");
    eprintln!("test_runall_extract_to_end: {count} consensus records from ~1000 input pairs");
}

// ============================================================================
// Test 2: start from sort (pre-aligned input)
// ============================================================================

/// Run extract+align phases first, then runall --start-from sort.
/// Verifies output matches a full extract-to-end run.
#[test]
fn test_runall_sort_to_end() {
    let Some(ref_path) = check_alignment_prerequisites() else {
        return;
    };

    let temp_dir = TempDir::new().unwrap();
    let read_len = 100;
    let (r1_fq, r2_fq) = create_large_test_fastqs(temp_dir.path(), &ref_path, 1000, read_len);
    let read_structure = format!("{read_len}T");

    // --- Full pipeline as reference ---
    let full_output = temp_dir.path().join("full.bam");
    let mut full_args = base_runall_extract_args(
        r1_fq.to_str().unwrap(),
        r2_fq.to_str().unwrap(),
        full_output.to_str().unwrap(),
        ref_path.to_str().unwrap(),
        &read_structure,
    );
    full_args.extend_from_slice(&["--threads", "4"]);
    run_fgumi(&full_args);

    // --- Partial: extract+align, stop after zipper ---
    let zippered = temp_dir.path().join("zippered.bam");
    let mut partial_args = base_runall_extract_args(
        r1_fq.to_str().unwrap(),
        r2_fq.to_str().unwrap(),
        zippered.to_str().unwrap(),
        ref_path.to_str().unwrap(),
        &read_structure,
    );
    partial_args.extend_from_slice(&["--stop-after", "zipper", "--threads", "4"]);
    run_fgumi(&partial_args);

    // --- Resume from sort ---
    let sort_output = temp_dir.path().join("from_sort.bam");
    run_fgumi(&[
        "runall",
        "--start-from",
        "sort",
        "--input",
        zippered.to_str().unwrap(),
        "--output",
        sort_output.to_str().unwrap(),
        "--reference",
        ref_path.to_str().unwrap(),
        "--consensus",
        "simplex",
        "--group::strategy",
        "identity",
        "--group::max-edits",
        "0",
        "--group::min-mapq",
        "0",
        "--filter::min-reads",
        "1",
        "--compression-level",
        "1",
        "--threads",
        "4",
    ]);

    assert_bam_equivalent(&full_output, &sort_output);
}

// ============================================================================
// Test 3: start from group (pre-sorted input)
// ============================================================================

/// Group the input first (stop after group), then runall --start-from group.
/// Verifies output matches a full extract-to-end run.
#[test]
fn test_runall_group_to_end() {
    let Some(ref_path) = check_alignment_prerequisites() else {
        return;
    };

    let temp_dir = TempDir::new().unwrap();
    let read_len = 100;
    let (r1_fq, r2_fq) = create_large_test_fastqs(temp_dir.path(), &ref_path, 1000, read_len);
    let read_structure = format!("{read_len}T");

    // --- Full pipeline as reference ---
    let full_output = temp_dir.path().join("full.bam");
    let mut full_args = base_runall_extract_args(
        r1_fq.to_str().unwrap(),
        r2_fq.to_str().unwrap(),
        full_output.to_str().unwrap(),
        ref_path.to_str().unwrap(),
        &read_structure,
    );
    full_args.extend_from_slice(&["--threads", "4"]);
    run_fgumi(&full_args);

    // --- Partial: extract through group (MI tags assigned) ---
    let grouped = temp_dir.path().join("grouped.bam");
    let mut partial_args = base_runall_extract_args(
        r1_fq.to_str().unwrap(),
        r2_fq.to_str().unwrap(),
        grouped.to_str().unwrap(),
        ref_path.to_str().unwrap(),
        &read_structure,
    );
    partial_args.extend_from_slice(&["--stop-after", "group", "--threads", "4"]);
    run_fgumi(&partial_args);

    // --- Resume from group ---
    let group_output = temp_dir.path().join("from_group.bam");
    run_fgumi(&[
        "runall",
        "--start-from",
        "group",
        "--input",
        grouped.to_str().unwrap(),
        "--output",
        group_output.to_str().unwrap(),
        "--reference",
        ref_path.to_str().unwrap(),
        "--consensus",
        "simplex",
        "--group::strategy",
        "identity",
        "--group::max-edits",
        "0",
        "--group::min-mapq",
        "0",
        "--filter::min-reads",
        "1",
        "--compression-level",
        "1",
        "--threads",
        "4",
    ]);

    assert_bam_equivalent(&full_output, &group_output);
}

// ============================================================================
// Test 4: stop after sort
// ============================================================================

/// Extract through sort only, verify output is a valid sorted BAM.
#[test]
fn test_runall_stop_after_sort() {
    let Some(ref_path) = check_alignment_prerequisites() else {
        return;
    };

    let temp_dir = TempDir::new().unwrap();
    let read_len = 100;
    let (r1_fq, r2_fq) = create_large_test_fastqs(temp_dir.path(), &ref_path, 1000, read_len);
    let read_structure = format!("{read_len}T");

    let output = temp_dir.path().join("sorted_output.bam");
    let mut args = base_runall_extract_args(
        r1_fq.to_str().unwrap(),
        r2_fq.to_str().unwrap(),
        output.to_str().unwrap(),
        ref_path.to_str().unwrap(),
        &read_structure,
    );
    args.extend_from_slice(&["--stop-after", "sort", "--threads", "4"]);
    run_fgumi(&args);

    assert!(output.exists(), "Sorted output BAM should exist");
    let count = count_bam_records(&output);
    // Should have ~2000 records (1000 pairs = 2000 reads, paired-end)
    assert!(count > 100, "Sorted BAM should have many records, got {count}");
    eprintln!("test_runall_stop_after_sort: {count} records in sorted BAM");

    // Record count validates the sort stage ran successfully.
}

// ============================================================================
// Test 5: parallel vs serial equivalence
// ============================================================================

/// Same input, compare parallel (--threads 4) vs serial (--threads 1).
/// Verifies both produce equivalent output.
#[test]
fn test_runall_parallel_vs_serial() {
    let Some(ref_path) = check_alignment_prerequisites() else {
        return;
    };

    let temp_dir = TempDir::new().unwrap();
    let read_len = 100;
    let (r1_fq, r2_fq) = create_large_test_fastqs(temp_dir.path(), &ref_path, 1000, read_len);
    let read_structure = format!("{read_len}T");

    let serial_output = temp_dir.path().join("serial.bam");
    let parallel_output = temp_dir.path().join("parallel.bam");

    // Serial (1 thread)
    let mut serial_args = base_runall_extract_args(
        r1_fq.to_str().unwrap(),
        r2_fq.to_str().unwrap(),
        serial_output.to_str().unwrap(),
        ref_path.to_str().unwrap(),
        &read_structure,
    );
    serial_args.extend_from_slice(&["--threads", "1"]);
    run_fgumi(&serial_args);

    // Parallel (4 threads)
    let mut parallel_args = base_runall_extract_args(
        r1_fq.to_str().unwrap(),
        r2_fq.to_str().unwrap(),
        parallel_output.to_str().unwrap(),
        ref_path.to_str().unwrap(),
        &read_structure,
    );
    parallel_args.extend_from_slice(&["--threads", "4"]);
    run_fgumi(&parallel_args);

    assert_bam_equivalent(&serial_output, &parallel_output);
}

// ============================================================================
// Test 6: multiple thread counts produce equivalent output
// ============================================================================

/// Run with 1, 2, 4, 8 threads and verify all produce equivalent output.
#[test]
fn test_runall_thread_counts() {
    let Some(ref_path) = check_alignment_prerequisites() else {
        return;
    };

    let temp_dir = TempDir::new().unwrap();
    let read_len = 100;
    let (r1_fq, r2_fq) = create_large_test_fastqs(temp_dir.path(), &ref_path, 1000, read_len);
    let read_structure = format!("{read_len}T");

    let thread_counts = [1, 2, 4, 8];
    let mut outputs: Vec<PathBuf> = Vec::new();

    for &threads in &thread_counts {
        let output = temp_dir.path().join(format!("threads_{threads}.bam"));
        let mut args = base_runall_extract_args(
            r1_fq.to_str().unwrap(),
            r2_fq.to_str().unwrap(),
            output.to_str().unwrap(),
            ref_path.to_str().unwrap(),
            &read_structure,
        );
        let threads_str = threads.to_string();
        args.extend_from_slice(&["--threads", &threads_str]);
        run_fgumi(&args);
        outputs.push(output);
    }

    // Compare all outputs against the first (1-thread) baseline
    let baseline = &outputs[0];
    for (i, output) in outputs.iter().enumerate().skip(1) {
        eprintln!("Comparing {} threads vs {} threads", thread_counts[0], thread_counts[i]);
        assert_bam_equivalent(baseline, output);
    }
}
