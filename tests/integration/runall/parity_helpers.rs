//! Shared fixture helpers for the runall consensus parity tests
//! (`test_consensus_codec.rs`, `test_consensus_simplex.rs`,
//! `test_consensus_duplex.rs`).
//!
//! Each parity test verifies that runall-with-`--consensus <mode>` produces
//! the same output as the equivalent standalone subprocess chain. The earlier
//! versions of these tests used a synthetic 1000bp `ACGTACGT…` repeating
//! reference, which made bwa-mem assign MAPQ=0 to every alignment; the
//! default `--group::min-mapq=1` then dropped every record at `GroupAssign`
//! and both pipelines emitted empty BAMs. The byte comparison "passed"
//! trivially, so any consensus-stage divergence (e.g. the codec
//! overlapping-bases bug) was invisible to the gate.
//!
//! These helpers fix that by:
//!
//! 1. Pointing `--reference` at the tracked `PhiX` FASTA (5386bp, non-repetitive,
//!    pre-indexed with `bwa index` and `samtools dict` in the repo).
//! 2. Exposing `assert_records_gt_zero` so each parity test can guard its
//!    own output: zero-record BAMs are a fixture failure, not a pass.
//!
//! Both helpers are feature-gated identically to the parity test files
//! themselves (via `runall/mod.rs`).

#![allow(dead_code, reason = "helpers used selectively across the three parity test files")]

use std::path::{Path, PathBuf};
use std::process::Command;

/// Locate the tracked `PhiX` reference FASTA (with bwa+samtools indexes).
///
/// Returns `None` and either skips locally or panics in CI if `PhiX` is not
/// available — matches the behavior of the per-test bwa/samtools tool checks.
pub(crate) fn phix_reference_or_skip(test_label: &str) -> Option<PathBuf> {
    let phix = crate::helpers::references::phix_reference();
    if phix.is_none() {
        crate::helpers::references::skip_or_fail(&format!(
            "`PhiX` reference not available at tests/data/references/phix.fasta; \
             skipping {test_label}"
        ));
    }
    phix
}

/// Read a single-contig FASTA file into a raw byte vector (sequence only,
/// header line stripped, newlines removed). Used to slice short read templates
/// out of the `PhiX` reference.
///
/// # Panics
///
/// Panics if the file cannot be read or contains no sequence.
pub(crate) fn read_fasta_sequence(path: &Path) -> Vec<u8> {
    let raw = std::fs::read_to_string(path).expect("read FASTA");
    let mut seq = Vec::with_capacity(raw.len());
    for line in raw.lines() {
        if line.starts_with('>') {
            continue;
        }
        seq.extend_from_slice(line.as_bytes());
    }
    assert!(!seq.is_empty(), "FASTA at {} contained no sequence", path.display());
    seq
}

/// Count records in a BAM via `samtools view -c`. Returns 0 on tool failure
/// (so the caller's `assert_records_gt_zero` produces a useful diagnostic).
pub(crate) fn count_bam_records(path: &Path) -> usize {
    let output = Command::new("samtools")
        .args(["view", "-c", path.to_str().unwrap()])
        .output()
        .expect("spawn samtools view -c");
    if !output.status.success() {
        return 0;
    }
    String::from_utf8_lossy(&output.stdout).trim().parse::<usize>().unwrap_or(0)
}

/// Guard: assert that a BAM contains at least one record.
///
/// A consensus parity gate that compares two empty BAMs is a vacuous pass —
/// any divergence in the consensus stage is invisible because nothing reaches
/// it. Every parity test must run this on both v1 and v2 outputs to ensure
/// the fixture is actually exercising the pipeline through filter.
///
/// # Panics
///
/// Panics with a clear message naming the BAM and pipeline label when the
/// count is zero.
pub(crate) fn assert_records_gt_zero(path: &Path, label: &str) {
    let n = count_bam_records(path);
    assert!(
        n > 0,
        "{label} produced 0 records at {} — parity fixture is too weak \
         to gate the consensus stage; check upstream stages (alignment, \
         grouping) for silent record loss",
        path.display()
    );
}
