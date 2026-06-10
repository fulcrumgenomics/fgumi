//! BAM-equivalence assertion for runall parity tests.
//!
//! Two BAM files are "parity-equivalent" when they contain the same
//! sequence of alignment records (`RecordBuf` byte-equal). Header
//! differences are tolerated for fields that legitimately diverge
//! between invocations (`@PG` carries the literal command-line string;
//! `@HD` may carry slight SO/GO/SS variation depending on which
//! pipeline produced the file).
//!
//! For the parity tests' purposes, the RECORDS are what defines the
//! contract: same input + same parameters must yield the same record
//! stream regardless of CLI entry point.

#![allow(dead_code)]

use std::fs;
use std::path::Path;

use noodles::bam;
use noodles::sam::alignment::record_buf::RecordBuf;

/// Read all records from a BAM file into a `Vec<RecordBuf>`, in file
/// order. Used by `assert_bams_record_equivalent` and by tests that
/// need to inspect record streams directly.
pub fn read_bam_records(path: &Path) -> Vec<RecordBuf> {
    let mut reader = bam::io::Reader::new(
        fs::File::open(path).unwrap_or_else(|e| panic!("open {}: {e}", path.display())),
    );
    let header =
        reader.read_header().unwrap_or_else(|e| panic!("read header from {}: {e}", path.display()));
    reader
        .record_bufs(&header)
        .collect::<std::io::Result<Vec<_>>>()
        .unwrap_or_else(|e| panic!("read records from {}: {e}", path.display()))
}

/// Assert that two BAMs contain the same record stream (count + each
/// record byte-equal as `RecordBuf`).
///
/// Header content is intentionally NOT compared — `@PG` lines record
/// the literal command-line string, which differs between runall and
/// the equivalent standalone/staged chain. The parity contract is
/// about the record stream, not the metadata.
///
/// # Panics
///
/// Panics with a diagnostic on the first mismatch (record count or
/// first divergent record index).
pub fn assert_bams_record_equivalent(a: &Path, b: &Path) {
    let recs_a = read_bam_records(a);
    let recs_b = read_bam_records(b);
    assert_eq!(
        recs_a.len(),
        recs_b.len(),
        "record count differs: {} has {}, {} has {}",
        a.display(),
        recs_a.len(),
        b.display(),
        recs_b.len(),
    );
    for (i, (ra, rb)) in recs_a.iter().zip(recs_b.iter()).enumerate() {
        assert_eq!(ra, rb, "record {i} differs between {} and {}", a.display(), b.display());
    }
}
