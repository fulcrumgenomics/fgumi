//! Error path integration tests.
//!
//! These tests verify that error conditions are handled correctly,
//! including validation failures, missing files, and invalid inputs.
//!
//! Two groups of tests live here:
//!
//! 1. Library-level error paths for `ConsensusSequence` and `FilterConfig`
//!    that cover invalid-argument panics and empty-input safety.
//! 2. CLI-level error paths for `fgumi runall` that spawn the compiled
//!    binary via `CARGO_BIN_EXE_fgumi`, feed deliberately broken inputs,
//!    and assert that the process fails with a clear, user-facing error
//!    message (not a panic or a generic "error"). These guard the four
//!    "obvious" failure modes surfaced by M1: corrupt BAM input, missing
//!    reference, invalid read structure, and a stage-ordering conflict
//!    between `--start-from` and `--stop-after`.

use std::process::Command;

use fgumi_lib::consensus::{ConsensusSequence, FilterConfig, FilterThresholds};

/// Absolute path to the compiled `fgumi` binary under test.
fn fgumi_bin() -> &'static str {
    env!("CARGO_BIN_EXE_fgumi")
}

// ==================== ConsensusSequence Error Paths ====================

#[test]
fn test_consensus_sequence_empty_operations() {
    let seq = ConsensusSequence::new();

    // Operations on empty sequence should be safe
    assert_eq!(seq.max_depth(), 0);
    assert_eq!(seq.min_depth(), 0);
    assert!((seq.error_rate() - 0.0).abs() < f32::EPSILON);
    assert!(seq.bases().is_empty());
    assert!(seq.quals().is_empty());
}

#[test]
#[should_panic(expected = "new_length")]
fn test_consensus_sequence_padded_shorter_panics() {
    let seq = ConsensusSequence::from_vecs(
        vec![b'A', b'C', b'G', b'T'],
        vec![30, 30, 30, 30],
        vec![10, 10, 10, 10],
        vec![0, 0, 0, 0],
    );

    // Trying to pad to shorter length should panic
    let _ = seq.padded(2, false, b'N', 2);
}

// ==================== FilterConfig Error Paths ====================

#[test]
#[should_panic(expected = "min-reads values must be specified high to low")]
fn test_filter_config_for_duplex_invalid_strand_greater_than_duplex() {
    let duplex =
        FilterThresholds { min_reads: 5, max_read_error_rate: 0.1, max_base_error_rate: 0.2 };
    let strand =
        FilterThresholds { min_reads: 10, max_read_error_rate: 0.1, max_base_error_rate: 0.2 };

    // strand.min_reads (10) > duplex.min_reads (5) should panic
    let _ = FilterConfig::for_duplex(duplex, strand, Some(13), None, 0.1);
}

#[test]
#[should_panic(expected = "max-read-error-rate for AB")]
fn test_filter_config_for_duplex_asymmetric_invalid_error_rate() {
    let duplex =
        FilterThresholds { min_reads: 10, max_read_error_rate: 0.05, max_base_error_rate: 0.1 };
    let ab = FilterThresholds { min_reads: 5, max_read_error_rate: 0.2, max_base_error_rate: 0.1 };
    let ba = FilterThresholds { min_reads: 3, max_read_error_rate: 0.1, max_base_error_rate: 0.2 };

    // ab.max_read_error_rate (0.2) > ba.max_read_error_rate (0.1) should panic
    let _ = FilterConfig::for_duplex_asymmetric(duplex, ab, ba, Some(13), None, 0.1);
}

// Note: ConsensusCallingOptions validation tests are in the binary crate's
// unit tests (src/commands/common.rs) since they require access to the
// command module which is not exposed from the library.

// ==================== runall CLI Error Paths ====================
//
// These tests spawn the real `fgumi` binary with deliberately broken inputs
// and assert on distinctive substrings in stderr to prove the *specific*
// error branch fired (not a panic, not a generic failure).

/// `runall --start-from group` must fail with a readable BAM/header error
/// when handed a file that is not a valid BAM.
#[test]
fn runall_rejects_corrupt_bam_input() {
    let tmp = tempfile::tempdir().unwrap();
    let bam = tmp.path().join("corrupt.bam");
    std::fs::write(&bam, b"not a bam file").unwrap();
    let out = tmp.path().join("out.bam");

    let output = Command::new(fgumi_bin())
        .args(["runall", "--start-from", "group"])
        // --filter::min-reads is required for the default stop-after (filter);
        // supply a value so validation doesn't short-circuit before the BAM
        // header-read error this test is asserting on.
        .args(["--filter::min-reads", "1"])
        .arg("--input")
        .arg(&bam)
        .arg("--output")
        .arg(&out)
        .output()
        .expect("spawn fgumi");

    assert!(!output.status.success(), "expected non-zero exit on corrupt BAM input");
    let stderr = String::from_utf8_lossy(&output.stderr);
    // fgumi currently reports "Failed to read header from: <path>" for a
    // non-BAM file; assert on that distinctive phrase so we catch
    // regressions where the error becomes a panic or a generic message.
    assert!(
        stderr.contains("Failed to read header"),
        "expected header-read error on corrupt BAM; got stderr:\n{stderr}"
    );
    assert!(
        stderr.contains("corrupt.bam"),
        "expected the offending input path in the error; got stderr:\n{stderr}"
    );
}

/// `runall --start-from sort --reference <missing>` must fail with a
/// reference-not-found error before any processing begins.
#[test]
fn runall_rejects_missing_reference() {
    let tmp = tempfile::tempdir().unwrap();
    // We need a valid input path so validation reaches the reference check;
    // the file can be empty because reference validation runs first after
    // `validate_file_exists(input, ...)`.
    let bam = tmp.path().join("input.bam");
    std::fs::write(&bam, b"").unwrap();
    let missing_ref = tmp.path().join("does-not-exist.fa");
    let out = tmp.path().join("out.bam");

    let output = Command::new(fgumi_bin())
        .args(["runall", "--start-from", "sort"])
        .arg("--input")
        .arg(&bam)
        .arg("--output")
        .arg(&out)
        .arg("--reference")
        .arg(&missing_ref)
        .output()
        .expect("spawn fgumi");

    assert!(!output.status.success(), "expected non-zero exit on missing reference");
    let stderr = String::from_utf8_lossy(&output.stderr);
    assert!(
        stderr.contains("Reference FASTA"),
        "expected 'Reference FASTA' in error; got stderr:\n{stderr}"
    );
    assert!(
        stderr.contains("does not exist"),
        "expected 'does not exist' in error; got stderr:\n{stderr}"
    );
}

/// `runall --start-from extract --extract::read-structures NOTAVALIDSTRUCTURE ...`
/// must fail with a parse error that identifies the offending read structure.
#[test]
fn runall_rejects_invalid_read_structure() {
    use std::io::Write;

    use flate2::{Compression, write::GzEncoder};

    let tmp = tempfile::tempdir().unwrap();
    // Create two empty (but valid) gzipped FASTQ files so validation passes
    // the existence check and reaches read-structure parsing.
    let r1 = tmp.path().join("r1.fastq.gz");
    let r2 = tmp.path().join("r2.fastq.gz");
    for path in [&r1, &r2] {
        let f = std::fs::File::create(path).unwrap();
        let mut enc = GzEncoder::new(f, Compression::default());
        enc.write_all(b"").unwrap();
        enc.finish().unwrap();
    }
    let out = tmp.path().join("out.bam");

    let output = Command::new(fgumi_bin())
        .args([
            "runall",
            "--start-from",
            "extract",
            "--stop-after",
            "extract",
            "--extract::sample",
            "S",
            "--extract::library",
            "L",
            "--extract::read-structures",
            "NOTAVALIDSTRUCTURE",
            "+T",
        ])
        .arg("--input")
        .arg(&r1)
        .arg(&r2)
        .arg("--output")
        .arg(&out)
        .output()
        .expect("spawn fgumi");

    assert!(!output.status.success(), "expected non-zero exit on invalid read structure");
    let stderr = String::from_utf8_lossy(&output.stderr);
    // The underlying `read-structure` crate reports
    // "Read structure missing length information: ..." for strings like
    // "NOTAVALIDSTRUCTURE" that have letters but no length spec.
    assert!(
        stderr.contains("Read structure"),
        "expected 'Read structure' in parse error; got stderr:\n{stderr}"
    );
    assert!(
        stderr.contains("OTAVALIDSTRUCTURE"),
        "expected the offending token (or its tail) in the error; got stderr:\n{stderr}"
    );
}

/// `runall --start-from group --stop-after sort` must fail with a
/// stage-ordering error because `sort` precedes `group` in the pipeline.
#[test]
fn runall_rejects_stop_after_before_start_from() {
    let tmp = tempfile::tempdir().unwrap();
    // Input content is irrelevant — validation of the stage ordering runs
    // before any BAM parsing, but `validate_file_exists` needs the path to
    // resolve. An empty file suffices.
    let bam = tmp.path().join("input.bam");
    std::fs::write(&bam, b"").unwrap();
    let out = tmp.path().join("out.bam");

    let output = Command::new(fgumi_bin())
        .args(["runall", "--start-from", "group", "--stop-after", "sort"])
        .arg("--input")
        .arg(&bam)
        .arg("--output")
        .arg(&out)
        .output()
        .expect("spawn fgumi");

    assert!(
        !output.status.success(),
        "expected non-zero exit when --stop-after precedes --start-from"
    );
    let stderr = String::from_utf8_lossy(&output.stderr);
    assert!(
        stderr.contains("--stop-after"),
        "expected '--stop-after' in stage-ordering error; got stderr:\n{stderr}"
    );
    assert!(
        stderr.contains("before --start-from") || stderr.contains("at or after the start stage"),
        "expected stage-ordering phrasing; got stderr:\n{stderr}"
    );
}
