//! Integration tests for the `fgumi pipeline` command.
//!
//! Tests that the pipeline command produces correct output and that it can be
//! run successfully with various options.

use noodles::bam;
use noodles::sam::alignment::io::Write as AlignmentWrite;
use noodles::sam::alignment::record::data::field::Tag;
use noodles::sam::alignment::record_buf::RecordBuf;
use std::fs;
use std::path::Path;
use std::process::Command;
use tempfile::TempDir;

use crate::helpers::bam_generator::{
    count_bam_records, create_coordinate_sorted_header, create_mapped_bam_for_sort,
    create_minimal_header, create_test_reference, create_umi_family, read_bam_records,
    read_bam_sequences,
};

/// Write grouped BAM file (reads grouped by MI tag).
fn create_grouped_bam(path: &Path, families: Vec<(&str, Vec<RecordBuf>)>) {
    let header = create_minimal_header("chr1", 10000);
    let mut writer =
        bam::io::Writer::new(fs::File::create(path).expect("Failed to create BAM file"));
    writer.write_header(&header).expect("Failed to write header");

    for (mi, records) in families {
        for mut record in records {
            use noodles::sam::alignment::record::data::field::Tag;
            use noodles::sam::alignment::record_buf::data::field::Value;
            let mi_tag = Tag::new(b'M', b'I');
            record.data_mut().insert(mi_tag, Value::from(mi));
            writer.write_alignment_record(&header, &record).expect("Failed to write record");
        }
    }
    writer.try_finish().expect("Failed to finish BAM");
}

/// Test that `fgumi pipeline --start-from group` produces valid output.
#[test]
fn test_pipeline_from_group_basic() {
    let temp_dir = TempDir::new().unwrap();
    let input_bam = temp_dir.path().join("input.bam");
    let output_bam = temp_dir.path().join("output.bam");

    // Create two families of 5 reads each with distinct UMIs
    let family1 = create_umi_family("ACGT", 5, "fam1", "ACGTACGT", 30);
    let family2 = create_umi_family("TGCA", 5, "fam2", "TTTTAAAA", 30);
    create_grouped_bam(&input_bam, vec![("1", family1), ("2", family2)]);

    let status = Command::new(env!("CARGO_BIN_EXE_fgumi"))
        .args([
            "runall",
            "--start-from",
            "group",
            "--input",
            input_bam.to_str().unwrap(),
            "--output",
            output_bam.to_str().unwrap(),
            "--consensus",
            "simplex",
            "--filter::min-reads",
            "1",
            "--compression-level",
            "1",
        ])
        .status()
        .expect("Failed to run pipeline command");

    assert!(status.success(), "Pipeline command failed");
    assert!(output_bam.exists(), "Output BAM not created");

    // Verify output has consensus reads
    let count = count_bam_records(&output_bam);
    assert!(count > 0, "Pipeline should produce consensus reads, got 0");
}

/// Test that `fgumi pipeline --start-from group` produces the same output as
/// `fgumi simplex` when given the same input and parameters.
#[test]
fn test_pipeline_equivalence_with_simplex() {
    let temp_dir = TempDir::new().unwrap();
    let input_bam = temp_dir.path().join("input.bam");
    let pipeline_output = temp_dir.path().join("pipeline_output.bam");
    let simplex_output = temp_dir.path().join("simplex_output.bam");

    // Create three families with different depths
    let family1 = create_umi_family("AAAA", 5, "fam1", "ACGTACGT", 30);
    let family2 = create_umi_family("CCCC", 3, "fam2", "TTTTAAAA", 30);
    let family3 = create_umi_family("GGGG", 7, "fam3", "GGGGCCCC", 30);
    create_grouped_bam(&input_bam, vec![("1", family1), ("2", family2), ("3", family3)]);

    // Run pipeline command
    let pipeline_status = Command::new(env!("CARGO_BIN_EXE_fgumi"))
        .args([
            "runall",
            "--start-from",
            "group",
            "--input",
            input_bam.to_str().unwrap(),
            "--output",
            pipeline_output.to_str().unwrap(),
            "--consensus",
            "simplex",
            "--filter::min-reads",
            "2",
            "--compression-level",
            "1",
        ])
        .status()
        .expect("Failed to run pipeline command");
    assert!(pipeline_status.success(), "Pipeline command failed");

    // Run simplex command with same parameters
    let simplex_status = Command::new(env!("CARGO_BIN_EXE_fgumi"))
        .args([
            "simplex",
            "--input",
            input_bam.to_str().unwrap(),
            "--output",
            simplex_output.to_str().unwrap(),
            "--min-reads",
            "2",
            "--compression-level",
            "1",
        ])
        .status()
        .expect("Failed to run simplex command");
    assert!(simplex_status.success(), "Simplex command failed");

    // Compare: both should produce the same number of consensus reads
    let pipeline_count = count_bam_records(&pipeline_output);
    let simplex_count = count_bam_records(&simplex_output);
    assert_eq!(
        pipeline_count, simplex_count,
        "Pipeline and simplex should produce the same number of consensus reads"
    );

    // Compare sequences
    let pipeline_seqs = read_bam_sequences(&pipeline_output);
    let simplex_seqs = read_bam_sequences(&simplex_output);
    assert_eq!(
        pipeline_seqs, simplex_seqs,
        "Pipeline and simplex should produce identical consensus sequences"
    );
}

/// Test that the pipeline command writes metrics when --metrics is specified.
#[test]
fn test_pipeline_metrics_output() {
    let temp_dir = TempDir::new().unwrap();
    let input_bam = temp_dir.path().join("input.bam");
    let output_bam = temp_dir.path().join("output.bam");
    let metrics_path = temp_dir.path().join("metrics.tsv");

    let family = create_umi_family("ACGT", 5, "fam1", "ACGTACGT", 30);
    create_grouped_bam(&input_bam, vec![("1", family)]);

    let status = Command::new(env!("CARGO_BIN_EXE_fgumi"))
        .args([
            "runall",
            "--start-from",
            "group",
            "--input",
            input_bam.to_str().unwrap(),
            "--output",
            output_bam.to_str().unwrap(),
            "--consensus",
            "simplex",
            "--filter::min-reads",
            "1",
            "--compression-level",
            "1",
            "--metrics",
            metrics_path.to_str().unwrap(),
        ])
        .status()
        .expect("Failed to run pipeline command");

    assert!(status.success(), "Pipeline command failed");
    assert!(metrics_path.exists(), "Metrics file not created");

    let metrics_content = fs::read_to_string(&metrics_path).unwrap();
    assert!(
        metrics_content.contains("stage\twall_time_secs\trecords_in\trecords_out"),
        "Metrics file should have header row"
    );
    assert!(metrics_content.contains("total\t"), "Metrics file should have total row");
}

/// Test that the pipeline temp file is cleaned up on success (atomic output).
#[test]
fn test_pipeline_atomic_output_no_tmp_on_success() {
    let temp_dir = TempDir::new().unwrap();
    let input_bam = temp_dir.path().join("input.bam");
    let output_bam = temp_dir.path().join("output.bam");
    let tmp_output = temp_dir.path().join("output.bam.tmp");

    let family = create_umi_family("ACGT", 5, "fam1", "ACGTACGT", 30);
    create_grouped_bam(&input_bam, vec![("1", family)]);

    let status = Command::new(env!("CARGO_BIN_EXE_fgumi"))
        .args([
            "runall",
            "--start-from",
            "group",
            "--input",
            input_bam.to_str().unwrap(),
            "--output",
            output_bam.to_str().unwrap(),
            "--consensus",
            "simplex",
            "--filter::min-reads",
            "1",
            "--compression-level",
            "1",
        ])
        .status()
        .expect("Failed to run pipeline command");

    assert!(status.success(), "Pipeline command failed");
    assert!(output_bam.exists(), "Final output should exist");
    assert!(!tmp_output.exists(), "Temp file should not remain after success");
}

// ============================================================================
// Helpers for --start-from sort tests
// ============================================================================

// ============================================================================
// Test 1: --start-from sort equivalence with sequential commands
// ============================================================================

/// Compare runall pipeline (--start-from sort) against sequential sort -> group -> simplex.
///
/// Both paths should produce the same number of consensus reads with the same sequences.
/// Uses --threads 3 to trigger the parallel pipeline path (serial path requires >= 3 threads).
#[test]
fn test_pipeline_start_from_sort_equivalence() {
    let temp_dir = TempDir::new().unwrap();
    let input_bam = temp_dir.path().join("input.bam");
    let ref_path = create_test_reference(temp_dir.path());

    create_mapped_bam_for_sort(&input_bam);

    // --- Sequential path: sort -> group -> simplex ---
    let sorted_bam = temp_dir.path().join("sorted.bam");
    let grouped_bam = temp_dir.path().join("grouped.bam");
    let sequential_output = temp_dir.path().join("sequential.bam");

    // Step 1: sort
    let status = Command::new(env!("CARGO_BIN_EXE_fgumi"))
        .args([
            "sort",
            "-i",
            input_bam.to_str().unwrap(),
            "-o",
            sorted_bam.to_str().unwrap(),
            "--order",
            "template-coordinate",
            "--max-memory",
            "64M",
        ])
        .status()
        .expect("Failed to run sort");
    assert!(status.success(), "sort failed");

    // Step 2: group
    let status = Command::new(env!("CARGO_BIN_EXE_fgumi"))
        .args([
            "group",
            "--input",
            sorted_bam.to_str().unwrap(),
            "--output",
            grouped_bam.to_str().unwrap(),
            "--raw-tag",
            "RX",
            "--assign-tag",
            "MI",
            "--strategy",
            "adjacency",
            "--edits",
            "1",
            "--compression-level",
            "1",
        ])
        .status()
        .expect("Failed to run group");
    assert!(status.success(), "group failed");

    // Step 3: simplex (includes consensus + filter)
    let status = Command::new(env!("CARGO_BIN_EXE_fgumi"))
        .args([
            "simplex",
            "--input",
            grouped_bam.to_str().unwrap(),
            "--output",
            sequential_output.to_str().unwrap(),
            "--min-reads",
            "1",
            "--compression-level",
            "1",
        ])
        .status()
        .expect("Failed to run simplex");
    assert!(status.success(), "simplex failed");

    // --- Pipeline path: runall --start-from sort (parallel, threads >= 3) ---
    let pipeline_output = temp_dir.path().join("pipeline.bam");
    let status = Command::new(env!("CARGO_BIN_EXE_fgumi"))
        .args([
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
            "1",
            "--compression-level",
            "1",
            "--threads",
            "3",
        ])
        .status()
        .expect("Failed to run pipeline");
    assert!(status.success(), "pipeline --start-from sort failed");

    // --- Compare ---
    let seq_count = count_bam_records(&sequential_output);
    let pipeline_count = count_bam_records(&pipeline_output);
    assert!(seq_count > 0, "Sequential path should produce consensus reads");
    assert_eq!(
        seq_count, pipeline_count,
        "Pipeline and sequential paths should produce the same number of consensus reads"
    );

    let seq_seqs = read_bam_sequences(&sequential_output);
    let pipeline_seqs = read_bam_sequences(&pipeline_output);
    assert_eq!(
        seq_seqs, pipeline_seqs,
        "Pipeline and sequential paths should produce identical consensus sequences"
    );
}

// ============================================================================
// Test 2: --stop-after sort
// ============================================================================

/// Verify that --start-from sort --stop-after sort produces a valid
/// template-coordinate sorted BAM (more records out than zero, valid BAM).
#[test]
fn test_pipeline_stop_after_sort() {
    let temp_dir = TempDir::new().unwrap();
    let input_bam = temp_dir.path().join("input.bam");
    let output_bam = temp_dir.path().join("output.bam");
    let ref_path = create_test_reference(temp_dir.path());

    create_mapped_bam_for_sort(&input_bam);

    let status = Command::new(env!("CARGO_BIN_EXE_fgumi"))
        .args([
            "runall",
            "--start-from",
            "sort",
            "--stop-after",
            "sort",
            "--input",
            input_bam.to_str().unwrap(),
            "--output",
            output_bam.to_str().unwrap(),
            "--reference",
            ref_path.to_str().unwrap(),
            "--consensus",
            "simplex",
            "--compression-level",
            "1",
        ])
        .status()
        .expect("Failed to run pipeline");
    assert!(status.success(), "pipeline --stop-after sort failed");
    assert!(output_bam.exists(), "Output BAM not created");

    // All 10 input records should be present in the sorted output
    let count = count_bam_records(&output_bam);
    assert_eq!(count, 10, "Sorted output should contain all 10 input records");
}

// ============================================================================
// Test 3: --stop-after group
// ============================================================================

/// Verify that --start-from sort --stop-after group produces a BAM with MI tags.
#[test]
fn test_pipeline_stop_after_group() {
    let temp_dir = TempDir::new().unwrap();
    let input_bam = temp_dir.path().join("input.bam");
    let output_bam = temp_dir.path().join("output.bam");
    let ref_path = create_test_reference(temp_dir.path());

    create_mapped_bam_for_sort(&input_bam);

    let status = Command::new(env!("CARGO_BIN_EXE_fgumi"))
        .args([
            "runall",
            "--start-from",
            "sort",
            "--stop-after",
            "group",
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
            "--compression-level",
            "1",
        ])
        .status()
        .expect("Failed to run pipeline");
    assert!(status.success(), "pipeline --stop-after group failed");
    assert!(output_bam.exists(), "Output BAM not created");

    // Verify all records have MI tags
    let records = read_bam_records(&output_bam);
    assert_eq!(records.len(), 10, "Grouped output should contain all 10 input records");

    let mi_tag = Tag::new(b'M', b'I');
    for record in &records {
        assert!(record.data().get(&mi_tag).is_some(), "Every record should have an MI tag");
    }
}

// ============================================================================
// Test 4: validation errors
// ============================================================================

/// Test that the pipeline fails when the input file does not exist.
#[test]
fn test_pipeline_validation_missing_input() {
    let temp_dir = TempDir::new().unwrap();
    let output_bam = temp_dir.path().join("output.bam");

    let output = Command::new(env!("CARGO_BIN_EXE_fgumi"))
        .args([
            "runall",
            "--start-from",
            "group",
            "--input",
            "/nonexistent/path.bam",
            "--output",
            output_bam.to_str().unwrap(),
            "--consensus",
            "simplex",
        ])
        .output()
        .expect("Failed to run pipeline");
    assert!(!output.status.success(), "Pipeline should fail for missing input file");
}

/// Test that --start-from sort fails when --reference is not provided.
#[test]
fn test_pipeline_validation_missing_reference() {
    let temp_dir = TempDir::new().unwrap();
    let input_bam = temp_dir.path().join("dummy_input.bam");
    let output_bam = temp_dir.path().join("output.bam");

    // Create a minimal valid BAM so validation gets past the file-exists check
    let header = create_coordinate_sorted_header("chr1", 10000);
    let mut writer = bam::io::Writer::new(fs::File::create(&input_bam).expect("create BAM"));
    writer.write_header(&header).expect("write header");
    writer.try_finish().expect("finish BAM");

    let output = Command::new(env!("CARGO_BIN_EXE_fgumi"))
        .args([
            "runall",
            "--start-from",
            "sort",
            "--input",
            input_bam.to_str().unwrap(),
            "--output",
            output_bam.to_str().unwrap(),
            "--consensus",
            "simplex",
        ])
        .output()
        .expect("Failed to run pipeline");
    assert!(
        !output.status.success(),
        "Pipeline should fail when --reference is missing for --start-from sort"
    );
    let stderr = String::from_utf8_lossy(&output.stderr);
    assert!(
        stderr.contains("reference") || stderr.contains("Reference"),
        "Error should mention reference; got: {stderr}"
    );
}

/// Test that --stop-after before --start-from is rejected.
#[test]
fn test_pipeline_validation_stop_before_start() {
    let temp_dir = TempDir::new().unwrap();
    let input_bam = temp_dir.path().join("dummy_input.bam");
    let output_bam = temp_dir.path().join("output.bam");
    let ref_path = create_test_reference(temp_dir.path());

    let header = create_coordinate_sorted_header("chr1", 10000);
    let mut writer = bam::io::Writer::new(fs::File::create(&input_bam).expect("create BAM"));
    writer.write_header(&header).expect("write header");
    writer.try_finish().expect("finish BAM");

    let output = Command::new(env!("CARGO_BIN_EXE_fgumi"))
        .args([
            "runall",
            "--start-from",
            "sort",
            "--stop-after",
            "extract",
            "--input",
            input_bam.to_str().unwrap(),
            "--output",
            output_bam.to_str().unwrap(),
            "--reference",
            ref_path.to_str().unwrap(),
            "--consensus",
            "simplex",
        ])
        .output()
        .expect("Failed to run pipeline");
    assert!(
        !output.status.success(),
        "Pipeline should fail when --stop-after is before --start-from"
    );
    let stderr = String::from_utf8_lossy(&output.stderr);
    assert!(
        stderr.contains("before") || stderr.contains("stop"),
        "Error should mention stop-after ordering; got: {stderr}"
    );
}

/// Test that --start-from extract fails when aligner command is not provided.
#[test]
fn test_pipeline_validation_missing_aligner() {
    let temp_dir = TempDir::new().unwrap();
    let output_bam = temp_dir.path().join("output.bam");
    let dummy_fq = temp_dir.path().join("dummy.fastq");
    fs::write(&dummy_fq, "@read1\nACGT\n+\nIIII\n").unwrap();
    let ref_path = create_test_reference(temp_dir.path());

    let output = Command::new(env!("CARGO_BIN_EXE_fgumi"))
        .args([
            "runall",
            "--start-from",
            "extract",
            "--input",
            dummy_fq.to_str().unwrap(),
            "--output",
            output_bam.to_str().unwrap(),
            "--reference",
            ref_path.to_str().unwrap(),
            "--consensus",
            "simplex",
            "--extract::read-structures",
            "8M",
            "--extract::sample",
            "S1",
            "--extract::library",
            "L1",
        ])
        .output()
        .expect("Failed to run pipeline");
    assert!(
        !output.status.success(),
        "Pipeline should fail when aligner command is missing for --start-from extract"
    );
    let stderr = String::from_utf8_lossy(&output.stderr);
    assert!(
        stderr.contains("aligner") || stderr.contains("Aligner"),
        "Error should mention aligner; got: {stderr}"
    );
}

// ============================================================================
// Test 5: parallel vs serial equivalence
// ============================================================================

/// Run the same --start-from sort input with different thread counts (both in the
/// parallel path, which requires >= 3 threads). Both should produce the same consensus reads.
#[test]
fn test_pipeline_parallel_thread_count_equivalence() {
    let temp_dir = TempDir::new().unwrap();
    let input_bam = temp_dir.path().join("input.bam");
    let ref_path = create_test_reference(temp_dir.path());
    let output_3t = temp_dir.path().join("output_3t.bam");
    let output_6t = temp_dir.path().join("output_6t.bam");

    create_mapped_bam_for_sort(&input_bam);

    // 3 threads (minimum for parallel path)
    let status = Command::new(env!("CARGO_BIN_EXE_fgumi"))
        .args([
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
        ])
        .status()
        .expect("Failed to run 3-thread pipeline");
    assert!(status.success(), "3-thread pipeline failed");

    // 6 threads
    let status = Command::new(env!("CARGO_BIN_EXE_fgumi"))
        .args([
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
        ])
        .status()
        .expect("Failed to run 6-thread pipeline");
    assert!(status.success(), "6-thread pipeline failed");

    // Compare record counts
    let count_3t = count_bam_records(&output_3t);
    let count_6t = count_bam_records(&output_6t);
    assert!(count_3t > 0, "3-thread pipeline should produce consensus reads");
    assert_eq!(
        count_3t, count_6t,
        "Different thread counts should produce the same number of consensus reads"
    );

    // Compare sequences (sort by sequence since parallel output order may differ)
    let mut seqs_3t = read_bam_sequences(&output_3t);
    let mut seqs_6t = read_bam_sequences(&output_6t);
    seqs_3t.sort();
    seqs_6t.sort();
    assert_eq!(
        seqs_3t, seqs_6t,
        "Different thread counts should produce identical consensus sequences"
    );
}
