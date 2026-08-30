//! Smoke test for the sort-command cutover (Task 3): `Sort::execute` runs its
//! coordinate sort through the declarative chain builder (`build_for(spec)?.run()`)
//! instead of the owned `RawExternalSorter` engine, and still produces
//! correctly sorted output.
//!
//! **Not a genuine RED/GREEN gate.** Pre-cutover, `execute_sort`'s owned-engine
//! path already sorts correctly, so this test passes before Task 3's change
//! too — there is no way to observe a real regression here that a broken
//! cutover wouldn't also need to break record count or ordering outright. It
//! exists as a basic end-to-end guard through the cutover, kept green on both
//! sides. The real parity gate — the chain's output is byte-identical to the
//! owned engine's, across sort orders and spill configurations — is Task 4's
//! job, not this test's.

use std::ffi::OsStr;

use clap::Parser;
use fgumi_lib::commands::command::Command as FgumiCommand;
use fgumi_lib::commands::sort::Sort;
use fgumi_raw_bam::{RawRecord, SamBuilder};
use tempfile::TempDir;

use crate::helpers::bam_generator::{create_minimal_header, write_bam};
use crate::helpers::read_bam_output;

/// `n` mapped, single-end records on one reference, placed at *descending*
/// positions so the input is deliberately unsorted — the sort stage has real
/// work to do before this test's assertions can distinguish a working
/// cutover from a no-op one.
fn unsorted_records(n: usize) -> Vec<RawRecord> {
    (0..n)
        .map(|i| {
            let pos = i32::try_from((n - i) * 100).expect("pos fits i32");
            let mut b = SamBuilder::new();
            b.read_name(format!("read{i}").as_bytes())
                .ref_id(0)
                .pos(pos)
                .mapq(60)
                .flags(0)
                .cigar_ops(&[4u32 << 4]) // 4M
                .sequence(b"ACGT")
                .qualities(&[30u8; 4]);
            b.build()
        })
        .collect()
}

#[test]
fn sort_command_produces_coordinate_sorted_output_via_chain() {
    let dir = TempDir::new().expect("create temp dir");
    let input_bam = dir.path().join("in.bam");
    let output_bam = dir.path().join("out.bam");

    let header = create_minimal_header("chr1", 1_000_000);
    let records = unsorted_records(500);
    let input_count = records.len();
    write_bam(&input_bam, &header, &records);

    // Pass `&OsStr` directly (not `&str`) so `try_parse_from` doesn't UTF-8
    // round-trip the temp-dir paths, matching the convention in
    // `test_sort_write_index.rs`.
    let cmd = Sort::try_parse_from([
        OsStr::new("sort"),
        OsStr::new("-i"),
        input_bam.as_os_str(),
        OsStr::new("-o"),
        output_bam.as_os_str(),
        OsStr::new("--order"),
        OsStr::new("coordinate"),
    ])
    .expect("failed to parse sort args");

    cmd.execute("fgumi sort").expect("sort command should succeed via the chain");

    let (_, out_records) = read_bam_output(&output_bam);
    assert_eq!(out_records.len(), input_count, "output record count must match input record count");

    let keys: Vec<(usize, usize)> = out_records
        .iter()
        .map(|r| {
            (
                r.reference_sequence_id().unwrap_or(usize::MAX),
                r.alignment_start().map_or(0, usize::from),
            )
        })
        .collect();
    assert!(
        keys.windows(2).all(|w| w[0] <= w[1]),
        "output is not coordinate-sorted (non-decreasing (ref_id, pos) expected): {keys:?}"
    );
}
