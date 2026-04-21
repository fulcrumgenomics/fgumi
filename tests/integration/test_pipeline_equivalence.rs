//! Equivalence integration tests for the `fgumi pipeline` command.
//!
//! Verifies that `fgumi pipeline --start-from <X>` produces identical (or record-level
//! equivalent) output to running the same stages as individual sequential commands.
//!
//! # Test matrix
//!
//! | Test                                          | Sequential path             | Pipeline path                         | External deps |
//! |-----------------------------------------------|-----------------------------|------------------------------------|---------------|
//! | `test_sort_group_simplex_equivalence`         | sort -> group -> simplex    | pipeline --start-from sort         | none          |
//! | `test_sort_group_simplex_parallel_vs_serial`  | pipeline threads=1          | pipeline threads=3                 | none          |
//! | `test_stop_after_sort_equivalence`            | sort (standalone)           | pipeline --stop-after sort         | none          |
//! | `test_stop_after_group_equivalence`           | sort -> group               | pipeline --stop-after group        | none          |
//! | `test_full_pipeline_with_bwa`                 | extract -> bwa -> zipper -> ... | pipeline --start-from extract   | bwa, PhiX ref |

use noodles::bam;
use noodles::sam::alignment::io::Write as AlignmentWrite;
use noodles::sam::alignment::record::data::field::Tag;
use std::fs;
use std::io::{BufRead, BufReader, Write};
use std::path::{Path, PathBuf};
use std::process::Command;
use tempfile::TempDir;

use crate::helpers::assertions::assert_bam_equivalent;
use crate::helpers::bam_generator::{
    count_bam_records, create_coordinate_sorted_header, create_mapped_bam_for_sort,
    create_test_reference, read_bam_records, read_bam_sequences,
};
use crate::helpers::references::check_alignment_prerequisites;
use fgumi_lib::sam::builder::RecordBuilder;

/// Build a mapped BAM with more position groups and varying family sizes for richer tests.
///
/// Creates reads at three positions with different UMIs and depths to exercise
/// the grouping and consensus logic more thoroughly.
fn create_mapped_bam_multi_group(path: &Path) {
    let header = create_coordinate_sorted_header("chr1", 10000);
    let mut writer =
        bam::io::Writer::new(fs::File::create(path).expect("Failed to create BAM file"));
    writer.write_header(&header).expect("Failed to write header");

    // Position 100: 4 reads, UMI "AAAA"
    for i in 0..4 {
        let rec = RecordBuilder::new()
            .name(&format!("pos100_{i}"))
            .sequence("ACGTACGTACGT")
            .qualities(&[35u8; 12])
            .reference_sequence_id(0)
            .alignment_start(100)
            .mapping_quality(60)
            .tag("RX", "AAAA")
            .build();
        writer.write_alignment_record(&header, &rec).expect("write");
    }

    // Position 300: 6 reads, UMI "TTTT"
    for i in 0..6 {
        let rec = RecordBuilder::new()
            .name(&format!("pos300_{i}"))
            .sequence("GGGGCCCCAAAA")
            .qualities(&[32u8; 12])
            .reference_sequence_id(0)
            .alignment_start(300)
            .mapping_quality(60)
            .tag("RX", "TTTT")
            .build();
        writer.write_alignment_record(&header, &rec).expect("write");
    }

    // Position 700: 3 reads, UMI "GGGG"
    for i in 0..3 {
        let rec = RecordBuilder::new()
            .name(&format!("pos700_{i}"))
            .sequence("AAAAAAAAAAAA")
            .qualities(&[28u8; 12])
            .reference_sequence_id(0)
            .alignment_start(700)
            .mapping_quality(60)
            .tag("RX", "GGGG")
            .build();
        writer.write_alignment_record(&header, &rec).expect("write");
    }

    writer.try_finish().expect("Failed to finish BAM");
}

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

/// Create synthetic paired-end FASTQ reads that align to a given reference.
///
/// Generates families of reads at specified positions with UMIs embedded in the
/// Illumina-style read name (field 8, colon-separated).
///
/// Returns `(r1_path, r2_path)`.
fn create_test_fastqs(dir: &Path, ref_path: &Path, read_len: usize) -> (PathBuf, PathBuf) {
    let genome = read_fasta_sequence(ref_path);
    assert!(genome.len() > 1000, "Reference sequence too short for test data generation");

    let r1_path = dir.join("r1.fastq");
    let r2_path = dir.join("r2.fastq");
    let mut r1_file = fs::File::create(&r1_path).unwrap();
    let mut r2_file = fs::File::create(&r2_path).unwrap();

    // Define families: (position, umi, count)
    let families: &[(usize, &str, usize)] = &[(100, "ACGT", 3), (400, "TTCC", 4), (700, "GGAA", 2)];

    for &(pos, umi, count) in families {
        let end = (pos + read_len).min(genome.len());
        let r1_seq: Vec<u8> = genome[pos..end].to_vec();
        let r2_seq = revcomp(&r1_seq);
        let qual = "I".repeat(r1_seq.len());

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
            writeln!(r2_file, "{qual}").unwrap();
        }
    }

    (r1_path, r2_path)
}

// ============================================================================
// Test 1: sort -> group -> simplex equivalence
// ============================================================================

/// Sequential sort -> group -> simplex should produce the same consensus reads as
/// `fgumi pipeline --start-from sort`.
///
/// Uses `--threads 3` for the runall path to exercise the parallel pipeline, and
/// `identity` grouping strategy with `--max-edits 0` for deterministic output.
#[test]
fn test_sort_group_simplex_equivalence() {
    let temp_dir = TempDir::new().unwrap();
    let input_bam = temp_dir.path().join("input.bam");
    let ref_path = create_test_reference(temp_dir.path());

    create_mapped_bam_multi_group(&input_bam);

    // --- Sequential path ---
    let sorted_bam = temp_dir.path().join("sorted.bam");
    let grouped_bam = temp_dir.path().join("grouped.bam");
    let sequential_output = temp_dir.path().join("sequential.bam");

    run_fgumi(&[
        "sort",
        "-i",
        input_bam.to_str().unwrap(),
        "-o",
        sorted_bam.to_str().unwrap(),
        "--order",
        "template-coordinate",
        "--max-memory",
        "64M",
    ]);

    run_fgumi(&[
        "group",
        "--input",
        sorted_bam.to_str().unwrap(),
        "--output",
        grouped_bam.to_str().unwrap(),
        "--strategy",
        "identity",
        "--edits",
        "0",
        "--compression-level",
        "1",
    ]);

    run_fgumi(&[
        "simplex",
        "--input",
        grouped_bam.to_str().unwrap(),
        "--output",
        sequential_output.to_str().unwrap(),
        "--min-reads",
        "1",
        "--compression-level",
        "1",
    ]);

    // --- Pipeline path (parallel, >= 3 threads) ---
    let pipeline_output = temp_dir.path().join("pipeline.bam");
    run_fgumi(&[
        "runall",
        "--start-from",
        "sort",
        "--input",
        input_bam.to_str().unwrap(),
        "--output",
        pipeline_output.to_str().unwrap(),
        "--reference",
        ref_path.to_str().unwrap(),
        "--consensus",
        "simplex",
        "--group::strategy",
        "identity",
        "--group::max-edits",
        "0",
        "--filter::min-reads",
        "1",
        "--compression-level",
        "1",
        "--threads",
        "3",
    ]);

    // --- Compare ---
    assert_bam_equivalent(&sequential_output, &pipeline_output);
}

// ============================================================================
// Test 2: parallel vs serial pipeline equivalence
// ============================================================================

/// `pipeline --start-from sort --threads 3` should produce the same consensus reads
/// as `pipeline --start-from sort --threads 6`.
///
/// Both use the parallel pipeline path (minimum 3 threads). Verifies that increasing
/// parallelism does not change the consensus output.
#[test]
fn test_sort_group_simplex_parallel_thread_counts() {
    let temp_dir = TempDir::new().unwrap();
    let input_bam = temp_dir.path().join("input.bam");
    let ref_path = create_test_reference(temp_dir.path());

    create_mapped_bam_multi_group(&input_bam);

    let output_3t = temp_dir.path().join("output_3t.bam");
    let output_6t = temp_dir.path().join("output_6t.bam");

    // 3 threads (minimum for parallel path)
    run_fgumi(&[
        "runall",
        "--start-from",
        "sort",
        "--input",
        input_bam.to_str().unwrap(),
        "--output",
        output_3t.to_str().unwrap(),
        "--reference",
        ref_path.to_str().unwrap(),
        "--consensus",
        "simplex",
        "--group::strategy",
        "identity",
        "--group::max-edits",
        "0",
        "--filter::min-reads",
        "1",
        "--compression-level",
        "1",
        "--threads",
        "3",
    ]);

    // 6 threads
    run_fgumi(&[
        "runall",
        "--start-from",
        "sort",
        "--input",
        input_bam.to_str().unwrap(),
        "--output",
        output_6t.to_str().unwrap(),
        "--reference",
        ref_path.to_str().unwrap(),
        "--consensus",
        "simplex",
        "--group::strategy",
        "identity",
        "--group::max-edits",
        "0",
        "--filter::min-reads",
        "1",
        "--compression-level",
        "1",
        "--threads",
        "6",
    ]);

    assert_bam_equivalent(&output_3t, &output_6t);
}

// ============================================================================
// Test 3: --stop-after sort equivalence
// ============================================================================

/// `fgumi sort --order template-coordinate` should produce the same records as
/// `fgumi pipeline --start-from sort --stop-after sort`.
#[test]
#[ignore = "plan 2c-2: pipeline does not yet honor --stop-after sort."]
fn test_stop_after_sort_equivalence() {
    let temp_dir = TempDir::new().unwrap();
    let input_bam = temp_dir.path().join("input.bam");
    let ref_path = create_test_reference(temp_dir.path());

    create_mapped_bam_for_sort(&input_bam);

    // Sequential: standalone sort
    let standalone_sorted = temp_dir.path().join("standalone_sorted.bam");
    run_fgumi(&[
        "sort",
        "-i",
        input_bam.to_str().unwrap(),
        "-o",
        standalone_sorted.to_str().unwrap(),
        "--order",
        "template-coordinate",
        "--max-memory",
        "64M",
    ]);

    // Runall: pipeline --stop-after sort
    let pipeline_sorted = temp_dir.path().join("pipeline_sorted.bam");
    run_fgumi(&[
        "runall",
        "--start-from",
        "sort",
        "--stop-after",
        "sort",
        "--input",
        input_bam.to_str().unwrap(),
        "--output",
        pipeline_sorted.to_str().unwrap(),
        "--reference",
        ref_path.to_str().unwrap(),
        "--consensus",
        "simplex",
        "--compression-level",
        "1",
    ]);

    // Both should have all 10 records
    let standalone_count = count_bam_records(&standalone_sorted);
    let pipeline_count = count_bam_records(&pipeline_sorted);
    assert_eq!(standalone_count, 10, "Standalone sort should preserve all 10 records");
    assert_eq!(pipeline_count, 10, "Pipeline sort should preserve all 10 records");

    // Compare sequences (sorted order should match since both use template-coordinate)
    let standalone_seqs = read_bam_sequences(&standalone_sorted);
    let pipeline_seqs = read_bam_sequences(&pipeline_sorted);
    assert_eq!(
        standalone_seqs, pipeline_seqs,
        "Standalone sort and pipeline sort should produce records in the same order"
    );
}

// ============================================================================
// Test 4: --stop-after group equivalence
// ============================================================================

/// `sort -> group` (sequential) should produce the same MI tags as
/// `pipeline --start-from sort --stop-after group`.
#[test]
#[ignore = "plan 2c-2: pipeline does not yet honor --stop-after group."]
fn test_stop_after_group_equivalence() {
    let temp_dir = TempDir::new().unwrap();
    let input_bam = temp_dir.path().join("input.bam");
    let ref_path = create_test_reference(temp_dir.path());

    create_mapped_bam_for_sort(&input_bam);

    // Sequential: sort -> group
    let sorted_bam = temp_dir.path().join("sorted.bam");
    let sequential_grouped = temp_dir.path().join("sequential_grouped.bam");

    run_fgumi(&[
        "sort",
        "-i",
        input_bam.to_str().unwrap(),
        "-o",
        sorted_bam.to_str().unwrap(),
        "--order",
        "template-coordinate",
        "--max-memory",
        "64M",
    ]);

    run_fgumi(&[
        "group",
        "--input",
        sorted_bam.to_str().unwrap(),
        "--output",
        sequential_grouped.to_str().unwrap(),
        "--strategy",
        "identity",
        "--edits",
        "0",
        "--compression-level",
        "1",
    ]);

    // Runall: pipeline --stop-after group
    let pipeline_grouped = temp_dir.path().join("pipeline_grouped.bam");
    run_fgumi(&[
        "runall",
        "--start-from",
        "sort",
        "--stop-after",
        "group",
        "--input",
        input_bam.to_str().unwrap(),
        "--output",
        pipeline_grouped.to_str().unwrap(),
        "--reference",
        ref_path.to_str().unwrap(),
        "--consensus",
        "simplex",
        "--group::strategy",
        "identity",
        "--group::max-edits",
        "0",
        "--compression-level",
        "1",
    ]);

    // Both should have all 10 records with MI tags
    let seq_records = read_bam_records(&sequential_grouped);
    let pipe_records = read_bam_records(&pipeline_grouped);
    assert_eq!(seq_records.len(), 10, "Sequential grouped should have 10 records");
    assert_eq!(pipe_records.len(), 10, "Pipeline grouped should have 10 records");

    let mi_tag = Tag::new(b'M', b'I');

    // Every record should have an MI tag in both outputs
    for rec in &seq_records {
        assert!(rec.data().get(&mi_tag).is_some(), "Sequential: every record needs MI tag");
    }
    for rec in &pipe_records {
        assert!(rec.data().get(&mi_tag).is_some(), "Pipeline: every record needs MI tag");
    }

    // Compare sequences to ensure same records are present
    let mut seq_seqs = read_bam_sequences(&sequential_grouped);
    let mut pipe_seqs = read_bam_sequences(&pipeline_grouped);
    seq_seqs.sort();
    pipe_seqs.sort();
    assert_eq!(seq_seqs, pipe_seqs, "Grouped outputs should contain the same sequences");
}

// ============================================================================
// Test 5: full pipeline with bwa
// ============================================================================

/// Run the sequential extract -> bwa -> zipper -> sort -> group -> simplex path.
///
/// Returns the path to the final consensus BAM.
fn run_sequential_extract_to_consensus(
    temp_dir: &Path,
    r1_fq: &Path,
    r2_fq: &Path,
    ref_path: &Path,
    read_len: usize,
) -> PathBuf {
    let read_structure = format!("{read_len}T");

    // Step 1: extract
    let extracted_bam = temp_dir.join("extracted.bam");
    run_fgumi(&[
        "extract",
        "--inputs",
        r1_fq.to_str().unwrap(),
        r2_fq.to_str().unwrap(),
        "--output",
        extracted_bam.to_str().unwrap(),
        "--read-structures",
        &read_structure,
        &read_structure,
        "--sample",
        "TestSample",
        "--library",
        "TestLibrary",
        "--extract-umis-from-read-names",
    ]);

    // Step 2: fastq (writes to stdout)
    let fastq_output = Command::new(env!("CARGO_BIN_EXE_fgumi"))
        .args(["fastq", "--input", extracted_bam.to_str().unwrap()])
        .output()
        .expect("Failed to run fastq");
    assert!(fastq_output.status.success(), "fgumi fastq failed");
    let fastq_for_align = temp_dir.join("for_align.fastq");
    fs::write(&fastq_for_align, &fastq_output.stdout).unwrap();

    // Step 3: bwa mem
    let aligned_sam = temp_dir.join("aligned.sam");
    let bwa_output = Command::new("bwa")
        .args([
            "mem",
            "-t",
            "2",
            "-p",
            "-K",
            "150000000",
            "-Y",
            ref_path.to_str().unwrap(),
            fastq_for_align.to_str().unwrap(),
        ])
        .output()
        .expect("Failed to run bwa mem");
    assert!(bwa_output.status.success(), "bwa mem failed");
    fs::write(&aligned_sam, &bwa_output.stdout).unwrap();

    // Step 4: zipper
    let zippered_bam = temp_dir.join("zippered.bam");
    run_fgumi(&[
        "zipper",
        "--input",
        aligned_sam.to_str().unwrap(),
        "--unmapped",
        extracted_bam.to_str().unwrap(),
        "--reference",
        ref_path.to_str().unwrap(),
        "--output",
        zippered_bam.to_str().unwrap(),
    ]);

    // Step 5: sort
    let sorted_bam = temp_dir.join("sorted.bam");
    run_fgumi(&[
        "sort",
        "-i",
        zippered_bam.to_str().unwrap(),
        "-o",
        sorted_bam.to_str().unwrap(),
        "--order",
        "template-coordinate",
        "--max-memory",
        "64M",
    ]);

    // Step 6: group
    let grouped_bam = temp_dir.join("grouped.bam");
    run_fgumi(&[
        "group",
        "--input",
        sorted_bam.to_str().unwrap(),
        "--output",
        grouped_bam.to_str().unwrap(),
        "--strategy",
        "identity",
        "--edits",
        "0",
        "--min-map-q",
        "0",
        "--compression-level",
        "1",
    ]);

    // Step 7: simplex
    let sequential_output = temp_dir.join("sequential_consensus.bam");
    run_fgumi(&[
        "simplex",
        "--input",
        grouped_bam.to_str().unwrap(),
        "--output",
        sequential_output.to_str().unwrap(),
        "--min-reads",
        "1",
        "--compression-level",
        "1",
    ]);

    sequential_output
}

/// End-to-end: sequential extract -> bwa mem -> zipper -> sort -> group -> simplex
/// vs. `fgumi pipeline --start-from extract --aligner::preset bwa-mem`.
///
/// Skipped when bwa is not on PATH or `PhiX` reference is not available.
#[test]
#[ignore = "plan 2c-2: pipeline does not yet implement --start-from extract."]
fn test_full_pipeline_with_bwa() {
    let Some(ref_path) = check_alignment_prerequisites() else {
        return;
    };

    let temp_dir = TempDir::new().unwrap();
    let read_len = 100;
    let (r1_fq, r2_fq) = create_test_fastqs(temp_dir.path(), &ref_path, read_len);
    let read_structure = format!("{read_len}T");

    // Sequential path
    let sequential_output =
        run_sequential_extract_to_consensus(temp_dir.path(), &r1_fq, &r2_fq, &ref_path, read_len);

    // Pipeline path
    let pipeline_output = temp_dir.path().join("pipeline_consensus.bam");
    run_fgumi(&[
        "runall",
        "--start-from",
        "extract",
        "--input",
        r1_fq.to_str().unwrap(),
        r2_fq.to_str().unwrap(),
        "--output",
        pipeline_output.to_str().unwrap(),
        "--reference",
        ref_path.to_str().unwrap(),
        "--consensus",
        "simplex",
        "--extract::sample",
        "TestSample",
        "--extract::library",
        "TestLibrary",
        "--extract::read-structures",
        &read_structure,
        &read_structure,
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
        "--threads",
        "3",
    ]);

    assert_bam_equivalent(&sequential_output, &pipeline_output);
}

// ============================================================================
// Test 6: adjacency strategy equivalence
// ============================================================================

/// Verify that adjacency grouping strategy produces the same output through both
/// sequential and runall paths. This is a more realistic scenario than identity
/// since most production runs use adjacency.
#[test]
fn test_sort_group_simplex_adjacency_equivalence() {
    let temp_dir = TempDir::new().unwrap();
    let input_bam = temp_dir.path().join("input.bam");
    let ref_path = create_test_reference(temp_dir.path());

    create_mapped_bam_multi_group(&input_bam);

    // --- Sequential path ---
    let sorted_bam = temp_dir.path().join("sorted.bam");
    let grouped_bam = temp_dir.path().join("grouped.bam");
    let sequential_output = temp_dir.path().join("sequential.bam");

    run_fgumi(&[
        "sort",
        "-i",
        input_bam.to_str().unwrap(),
        "-o",
        sorted_bam.to_str().unwrap(),
        "--order",
        "template-coordinate",
        "--max-memory",
        "64M",
    ]);

    run_fgumi(&[
        "group",
        "--input",
        sorted_bam.to_str().unwrap(),
        "--output",
        grouped_bam.to_str().unwrap(),
        "--strategy",
        "adjacency",
        "--edits",
        "1",
        "--compression-level",
        "1",
    ]);

    run_fgumi(&[
        "simplex",
        "--input",
        grouped_bam.to_str().unwrap(),
        "--output",
        sequential_output.to_str().unwrap(),
        "--min-reads",
        "2",
        "--compression-level",
        "1",
    ]);

    // --- Pipeline path ---
    let pipeline_output = temp_dir.path().join("pipeline.bam");
    run_fgumi(&[
        "runall",
        "--start-from",
        "sort",
        "--input",
        input_bam.to_str().unwrap(),
        "--output",
        pipeline_output.to_str().unwrap(),
        "--reference",
        ref_path.to_str().unwrap(),
        "--consensus",
        "simplex",
        "--group::strategy",
        "adjacency",
        "--group::max-edits",
        "1",
        "--filter::min-reads",
        "2",
        "--compression-level",
        "1",
        "--threads",
        "3",
    ]);

    assert_bam_equivalent(&sequential_output, &pipeline_output);
}

// ============================================================================
// Test 7: duplex consensus equivalence
// ============================================================================

/// Create a mapped BAM with paired-end reads and paired UMIs for duplex testing.
///
/// Produces reads at one position with A-B and B-A UMI pairing, suitable for
/// the `paired` grouping strategy. Each strand has reads with complementary
/// paired UMIs (e.g., "ACGT-TTGG" and "TTGG-ACGT") to form duplex families.
///
/// All reads carry CIGAR strings and MC (mate CIGAR) tags, which are required
/// by the standalone `group` command for paired-end template-coordinate grouping.
fn create_mapped_bam_duplex(path: &Path) {
    let header = create_coordinate_sorted_header("chr1", 10000);
    let mut writer =
        bam::io::Writer::new(fs::File::create(path).expect("Failed to create BAM file"));
    writer.write_header(&header).expect("Failed to write header");

    // Position 100: 3 read pairs with UMI "ACGT-TTGG" (A strand, FR orientation)
    for i in 0..3 {
        let name = format!("pos100_ab_{i}");
        // R1: forward at 100
        let r1 = RecordBuilder::new()
            .name(&name)
            .sequence("ACGTACGTACGT")
            .qualities(&[35u8; 12])
            .paired(true)
            .first_segment(true)
            .reverse_complement(false)
            .mate_reverse_complement(true)
            .reference_sequence_id(0)
            .alignment_start(100)
            .mate_reference_sequence_id(0)
            .mate_alignment_start(200)
            .template_length(112)
            .mapping_quality(60)
            .cigar("12M")
            .tag("RX", "ACGT-TTGG")
            .tag("MC", "12M")
            .build();
        writer.write_alignment_record(&header, &r1).expect("write");
        // R2: reverse at 200
        let r2 = RecordBuilder::new()
            .name(&name)
            .sequence("ACGTACGTACGT")
            .qualities(&[35u8; 12])
            .paired(true)
            .first_segment(false)
            .reverse_complement(true)
            .mate_reverse_complement(false)
            .reference_sequence_id(0)
            .alignment_start(200)
            .mate_reference_sequence_id(0)
            .mate_alignment_start(100)
            .template_length(-112)
            .mapping_quality(60)
            .cigar("12M")
            .tag("RX", "ACGT-TTGG")
            .tag("MC", "12M")
            .build();
        writer.write_alignment_record(&header, &r2).expect("write");
    }

    // Position 100: 3 read pairs with UMI "TTGG-ACGT" (B strand, RF orientation)
    for i in 0..3 {
        let name = format!("pos100_ba_{i}");
        // R1: reverse at 200
        let r1 = RecordBuilder::new()
            .name(&name)
            .sequence("ACGTACGTACGT")
            .qualities(&[35u8; 12])
            .paired(true)
            .first_segment(true)
            .reverse_complement(true)
            .mate_reverse_complement(false)
            .reference_sequence_id(0)
            .alignment_start(200)
            .mate_reference_sequence_id(0)
            .mate_alignment_start(100)
            .template_length(-112)
            .mapping_quality(60)
            .cigar("12M")
            .tag("RX", "TTGG-ACGT")
            .tag("MC", "12M")
            .build();
        writer.write_alignment_record(&header, &r1).expect("write");
        // R2: forward at 100
        let r2 = RecordBuilder::new()
            .name(&name)
            .sequence("ACGTACGTACGT")
            .qualities(&[35u8; 12])
            .paired(true)
            .first_segment(false)
            .reverse_complement(false)
            .mate_reverse_complement(true)
            .reference_sequence_id(0)
            .alignment_start(100)
            .mate_reference_sequence_id(0)
            .mate_alignment_start(200)
            .template_length(112)
            .mapping_quality(60)
            .cigar("12M")
            .tag("RX", "TTGG-ACGT")
            .tag("MC", "12M")
            .build();
        writer.write_alignment_record(&header, &r2).expect("write");
    }

    writer.try_finish().expect("Failed to finish BAM");
}

/// Sequential sort -> group(paired) -> duplex should produce the same consensus reads
/// as `fgumi pipeline --start-from sort --consensus duplex --group::strategy paired`.
#[test]
#[ignore = "plan 2c-2: pipeline only supports --consensus simplex; duplex is deferred."]
fn test_duplex_consensus_equivalence() {
    let temp_dir = TempDir::new().unwrap();
    let input_bam = temp_dir.path().join("input.bam");
    let ref_path = create_test_reference(temp_dir.path());

    create_mapped_bam_duplex(&input_bam);

    // --- Sequential path ---
    let sorted_bam = temp_dir.path().join("sorted.bam");
    let grouped_bam = temp_dir.path().join("grouped.bam");
    let sequential_output = temp_dir.path().join("sequential.bam");

    run_fgumi(&[
        "sort",
        "-i",
        input_bam.to_str().unwrap(),
        "-o",
        sorted_bam.to_str().unwrap(),
        "--order",
        "template-coordinate",
        "--max-memory",
        "64M",
    ]);

    run_fgumi(&[
        "group",
        "--input",
        sorted_bam.to_str().unwrap(),
        "--output",
        grouped_bam.to_str().unwrap(),
        "--strategy",
        "paired",
        "--edits",
        "1",
        "--compression-level",
        "1",
    ]);

    run_fgumi(&[
        "duplex",
        "--input",
        grouped_bam.to_str().unwrap(),
        "--output",
        sequential_output.to_str().unwrap(),
        "--min-reads",
        "1",
        "--compression-level",
        "1",
    ]);

    // --- Pipeline path ---
    let pipeline_output = temp_dir.path().join("pipeline.bam");
    run_fgumi(&[
        "runall",
        "--start-from",
        "sort",
        "--input",
        input_bam.to_str().unwrap(),
        "--output",
        pipeline_output.to_str().unwrap(),
        "--reference",
        ref_path.to_str().unwrap(),
        "--consensus",
        "duplex",
        "--group::strategy",
        "paired",
        "--group::max-edits",
        "1",
        "--duplex::min-reads",
        "1",
        "--filter::min-reads",
        "1",
        "--compression-level",
        "1",
        "--threads",
        "3",
    ]);

    // Both should produce consensus reads
    let seq_count = count_bam_records(&sequential_output);
    let pipeline_count = count_bam_records(&pipeline_output);
    assert!(seq_count > 0, "Sequential duplex should produce consensus reads");
    assert_eq!(
        seq_count, pipeline_count,
        "Duplex consensus: sequential ({seq_count}) vs pipeline ({pipeline_count}) record count mismatch"
    );

    // Compare sorted sequence multisets
    let mut seq_seqs = read_bam_sequences(&sequential_output);
    let mut pipeline_seqs = read_bam_sequences(&pipeline_output);
    seq_seqs.sort();
    pipeline_seqs.sort();
    assert_eq!(
        seq_seqs, pipeline_seqs,
        "Duplex consensus sequences should match between sequential and pipeline paths"
    );
}

// ============================================================================
// Test 8: codec consensus equivalence
// ============================================================================

/// Create a mapped BAM with paired-end reads suitable for CODEC consensus.
///
/// CODEC reads are FR-oriented pairs where R1 and R2 sequence opposite strands
/// of the same molecule. Even a single read pair can produce duplex consensus.
///
/// All reads carry CIGAR strings and MC (mate CIGAR) tags, required by the
/// standalone `group` command for paired-end template-coordinate grouping.
fn create_mapped_bam_codec(path: &Path) {
    let header = create_coordinate_sorted_header("chr1", 10000);
    let mut writer =
        bam::io::Writer::new(fs::File::create(path).expect("Failed to create BAM file"));
    writer.write_header(&header).expect("Failed to write header");

    // 4 read pairs with the same UMI at position 100, FR orientation
    for i in 0..4 {
        let name = format!("codec_{i}");
        // R1: forward at 100
        let r1 = RecordBuilder::new()
            .name(&name)
            .sequence("ACGTACGTACGT")
            .qualities(&[30u8; 12])
            .paired(true)
            .first_segment(true)
            .reverse_complement(false)
            .mate_reverse_complement(true)
            .reference_sequence_id(0)
            .alignment_start(100)
            .mate_reference_sequence_id(0)
            .mate_alignment_start(100)
            .template_length(12)
            .mapping_quality(60)
            .cigar("12M")
            .tag("RX", "AAAA")
            .tag("MC", "12M")
            .build();
        writer.write_alignment_record(&header, &r1).expect("write");
        // R2: reverse at 100 (same position, overlapping for CODEC)
        let r2 = RecordBuilder::new()
            .name(&name)
            .sequence("ACGTACGTACGT")
            .qualities(&[30u8; 12])
            .paired(true)
            .first_segment(false)
            .reverse_complement(true)
            .mate_reverse_complement(false)
            .reference_sequence_id(0)
            .alignment_start(100)
            .mate_reference_sequence_id(0)
            .mate_alignment_start(100)
            .template_length(-12)
            .mapping_quality(60)
            .cigar("12M")
            .tag("RX", "AAAA")
            .tag("MC", "12M")
            .build();
        writer.write_alignment_record(&header, &r2).expect("write");
    }

    writer.try_finish().expect("Failed to finish BAM");
}

/// Sequential sort -> group -> codec should produce the same consensus reads
/// as `fgumi pipeline --start-from sort --consensus codec`.
#[test]
#[ignore = "plan 2c-2: pipeline only supports --consensus simplex; codec is deferred."]
fn test_codec_consensus_equivalence() {
    let temp_dir = TempDir::new().unwrap();
    let input_bam = temp_dir.path().join("input.bam");
    let ref_path = create_test_reference(temp_dir.path());

    create_mapped_bam_codec(&input_bam);

    // --- Sequential path ---
    let sorted_bam = temp_dir.path().join("sorted.bam");
    let grouped_bam = temp_dir.path().join("grouped.bam");
    let sequential_output = temp_dir.path().join("sequential.bam");

    run_fgumi(&[
        "sort",
        "-i",
        input_bam.to_str().unwrap(),
        "-o",
        sorted_bam.to_str().unwrap(),
        "--order",
        "template-coordinate",
        "--max-memory",
        "64M",
    ]);

    run_fgumi(&[
        "group",
        "--input",
        sorted_bam.to_str().unwrap(),
        "--output",
        grouped_bam.to_str().unwrap(),
        "--strategy",
        "identity",
        "--edits",
        "0",
        "--compression-level",
        "1",
    ]);

    run_fgumi(&[
        "codec",
        "--input",
        grouped_bam.to_str().unwrap(),
        "--output",
        sequential_output.to_str().unwrap(),
        "--min-reads",
        "1",
        "--compression-level",
        "1",
    ]);

    // --- Pipeline path ---
    let pipeline_output = temp_dir.path().join("pipeline.bam");
    run_fgumi(&[
        "runall",
        "--start-from",
        "sort",
        "--input",
        input_bam.to_str().unwrap(),
        "--output",
        pipeline_output.to_str().unwrap(),
        "--reference",
        ref_path.to_str().unwrap(),
        "--consensus",
        "codec",
        "--group::strategy",
        "identity",
        "--group::max-edits",
        "0",
        "--codec::min-reads",
        "1",
        "--filter::min-reads",
        "1",
        "--compression-level",
        "1",
        "--threads",
        "3",
    ]);

    // Both should produce consensus reads
    let seq_count = count_bam_records(&sequential_output);
    let pipeline_count = count_bam_records(&pipeline_output);
    assert!(seq_count > 0, "Sequential codec should produce consensus reads");
    assert_eq!(
        seq_count, pipeline_count,
        "Codec consensus: sequential ({seq_count}) vs pipeline ({pipeline_count}) record count mismatch"
    );

    let mut seq_seqs = read_bam_sequences(&sequential_output);
    let mut pipeline_seqs = read_bam_sequences(&pipeline_output);
    seq_seqs.sort();
    pipeline_seqs.sort();
    assert_eq!(
        seq_seqs, pipeline_seqs,
        "Codec consensus sequences should match between sequential and pipeline paths"
    );
}

// ============================================================================
// Test 10: --stop-after consensus
// ============================================================================

/// Verify that `--stop-after consensus` produces consensus reads. Since the serial
/// path (threads < 4) does not skip the filter stage for `--stop-after consensus`,
/// we verify that output is produced and contains consensus tags.
#[test]
fn test_stop_after_consensus() {
    let temp_dir = TempDir::new().unwrap();
    let input_bam = temp_dir.path().join("input.bam");
    let ref_path = create_test_reference(temp_dir.path());

    create_mapped_bam_multi_group(&input_bam);

    let output_bam = temp_dir.path().join("consensus_only.bam");
    run_fgumi(&[
        "runall",
        "--start-from",
        "sort",
        "--stop-after",
        "consensus",
        "--input",
        input_bam.to_str().unwrap(),
        "--output",
        output_bam.to_str().unwrap(),
        "--reference",
        ref_path.to_str().unwrap(),
        "--consensus",
        "simplex",
        "--group::strategy",
        "identity",
        "--group::max-edits",
        "0",
        "--filter::min-reads",
        "1",
        "--compression-level",
        "1",
        "--threads",
        "3",
    ]);

    // Verify output has consensus reads with cD (max depth) tags
    let records = read_bam_records(&output_bam);
    assert!(!records.is_empty(), "Stop-after consensus should produce output records");

    let cd_tag = Tag::new(b'c', b'D');
    for record in &records {
        assert!(
            record.data().get(&cd_tag).is_some(),
            "Consensus records should have cD (max depth) tag"
        );
    }
}

// ============================================================================
// Test 11: empty input
// ============================================================================

/// An empty BAM (valid header, zero records) should produce zero-record output
/// without crashing.
#[test]
fn test_empty_input_produces_empty_output() {
    let temp_dir = TempDir::new().unwrap();
    let input_bam = temp_dir.path().join("empty.bam");
    let ref_path = create_test_reference(temp_dir.path());

    // Write a valid BAM with header but no records
    let header = create_coordinate_sorted_header("chr1", 10000);
    let mut writer =
        bam::io::Writer::new(fs::File::create(&input_bam).expect("Failed to create BAM file"));
    writer.write_header(&header).expect("Failed to write header");
    writer.try_finish().expect("Failed to finish BAM");

    let output_bam = temp_dir.path().join("output.bam");
    run_fgumi(&[
        "runall",
        "--start-from",
        "sort",
        "--input",
        input_bam.to_str().unwrap(),
        "--output",
        output_bam.to_str().unwrap(),
        "--reference",
        ref_path.to_str().unwrap(),
        "--consensus",
        "simplex",
        "--group::strategy",
        "identity",
        "--group::max-edits",
        "0",
        "--filter::min-reads",
        "1",
        "--compression-level",
        "1",
    ]);

    let count = count_bam_records(&output_bam);
    assert_eq!(count, 0, "Empty input should produce zero output records");
}

// ============================================================================
// Test 12: single-read families filtered by min-reads
// ============================================================================

/// With `--filter::min-reads 2`, position groups containing only single-read
/// families should produce no consensus reads, while groups with >= 2 reads
/// should produce consensus.
#[test]
fn test_single_read_families_filtered_by_min_reads() {
    let temp_dir = TempDir::new().unwrap();
    let input_bam = temp_dir.path().join("input.bam");
    let ref_path = create_test_reference(temp_dir.path());

    // Create a BAM with:
    //  - 1 read at pos 100 with UMI "AAAA" (single-read family, should be filtered)
    //  - 3 reads at pos 500 with UMI "CCCC" (multi-read family, should pass)
    let header = create_coordinate_sorted_header("chr1", 10000);
    let mut writer =
        bam::io::Writer::new(fs::File::create(&input_bam).expect("Failed to create BAM file"));
    writer.write_header(&header).expect("Failed to write header");

    // Single-read family
    let rec = RecordBuilder::new()
        .name("single_0")
        .sequence("ACGTACGTACGT")
        .qualities(&[35u8; 12])
        .reference_sequence_id(0)
        .alignment_start(100)
        .mapping_quality(60)
        .tag("RX", "AAAA")
        .build();
    writer.write_alignment_record(&header, &rec).expect("write");

    // Multi-read family
    for i in 0..3 {
        let rec = RecordBuilder::new()
            .name(&format!("multi_{i}"))
            .sequence("GGGGCCCCAAAA")
            .qualities(&[32u8; 12])
            .reference_sequence_id(0)
            .alignment_start(500)
            .mapping_quality(60)
            .tag("RX", "CCCC")
            .build();
        writer.write_alignment_record(&header, &rec).expect("write");
    }
    writer.try_finish().expect("Failed to finish BAM");

    let output_bam = temp_dir.path().join("output.bam");
    run_fgumi(&[
        "runall",
        "--start-from",
        "sort",
        "--input",
        input_bam.to_str().unwrap(),
        "--output",
        output_bam.to_str().unwrap(),
        "--reference",
        ref_path.to_str().unwrap(),
        "--consensus",
        "simplex",
        "--group::strategy",
        "identity",
        "--group::max-edits",
        "0",
        "--filter::min-reads",
        "2",
        "--compression-level",
        "1",
    ]);

    // Should produce exactly 1 consensus read (from the 3-read family)
    let count = count_bam_records(&output_bam);
    assert_eq!(
        count, 1,
        "With min-reads=2, only the 3-read family should produce consensus (got {count})"
    );

    // The consensus sequence should come from the multi-read family
    let seqs = read_bam_sequences(&output_bam);
    assert_eq!(seqs.len(), 1);
    assert_eq!(seqs[0], b"GGGGCCCCAAAA", "Consensus sequence should match the multi-read family");
}
