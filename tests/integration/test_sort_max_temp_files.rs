//! `--max-temp-files` on the sort command: the live spill runs are bounded by
//! consolidating them, and consolidation never changes the output.
//!
//! The arena sort path used to resolve and log the limit without enforcing it,
//! so every spilled run was held open for one merge as wide as the run count
//! (fgumi#991). These tests pin both halves of the fix: the sort consolidates
//! when it spills more runs than the limit allows, and its records are
//! byte-identical to the same sort run without a binding limit — including the
//! order of records whose sort keys tie, which is where a consolidation that
//! merged runs out of order would show.

use std::path::{Path, PathBuf};
use std::process::Command;

use fgumi_lib::sam::SamTag;
use fgumi_raw_bam::{RawRecord, SamBuilder, flags};
use rstest::rstest;
use tempfile::TempDir;

use crate::helpers::bam_generator::{create_minimal_header, write_bam};
use crate::helpers::cutover::decompressed_records_without_pg;

/// A limit far above any run count these fixtures produce, so it never binds.
const UNBOUNDED: &str = "100000";

/// Unsorted single-end records engineered to tie under every sort order: a small
/// pool of positions (coordinate ties), a small pool of read names sharing flags
/// (queryname ties), and a unique `RX` per record so that two tied records are
/// still distinguishable in the output bytes. Every 17th record is unmapped
/// (no reference, no position), so all of those tie on the coordinate sort's
/// unmapped sentinel key and only run order separates them.
fn tie_heavy_records(n: usize) -> Vec<RawRecord> {
    records(n, |i| format!("q{}", (i * 31) % (n / 4)))
}

/// Like [`tie_heavy_records`] but with a unique name per record, for commands
/// (such as `group`) that reject two primary reads sharing a name.
fn uniquely_named_records(n: usize) -> Vec<RawRecord> {
    records(n, |i| format!("q{i}"))
}

fn records<F: Fn(usize) -> String>(n: usize, name_of: F) -> Vec<RawRecord> {
    (0..n)
        .map(|i| {
            let name = name_of(i);
            let mut b = SamBuilder::new();
            b.read_name(name.as_bytes()).sequence(b"ACGTACGTAC").qualities(&[30u8; 10]);
            if i % 17 == 0 {
                b.ref_id(-1).pos(-1).mapq(0).flags(flags::UNMAPPED);
            } else {
                let pos = i32::try_from(1 + (i * 7919 % 500) * 10).expect("pos fits i32");
                let strand = if i % 2 == 0 { 0 } else { flags::REVERSE };
                b.ref_id(0).pos(pos).mapq(60).flags(strand).cigar_ops(&[10u32 << 4]); // 10M
            }
            b.add_string_tag(SamTag::RX, format!("u{i}").as_bytes());
            b.build()
        })
        .collect()
}

fn write_fixture(dir: &Path, n: usize) -> PathBuf {
    let input = dir.join("unsorted.bam");
    write_bam(&input, &create_minimal_header("chr1", 1_000_000), &tie_heavy_records(n));
    input
}

/// Run `fgumi sort` as a subprocess (so its logs can be read) and return the
/// output path and stderr.
fn sort(
    input: &Path,
    dir: &Path,
    name: &str,
    order: &str,
    threads: &str,
    limit: &str,
    extra_args: &[&str],
) -> (PathBuf, String) {
    let output = dir.join(format!("{name}.bam"));
    let result = Command::new(env!("CARGO_BIN_EXE_fgumi"))
        .env("RUST_LOG", "info")
        .args(["sort", "-i"])
        .arg(input)
        .arg("-o")
        .arg(&output)
        .args(["--order", order, "--threads", threads, "--max-temp-files", limit])
        // A tiny, fixed total budget forces many spilled runs from a small file.
        .args(["--max-memory", "64K", "--memory-per-thread", "false"])
        .args(extra_args)
        .output()
        .expect("run fgumi sort");
    let stderr = String::from_utf8_lossy(&result.stderr).into_owned();
    assert!(result.status.success(), "fgumi sort failed:\n{stderr}");
    (output, stderr)
}

/// The integer after `label` on the first stderr line containing it.
fn logged_count(stderr: &str, label: &str) -> Option<u64> {
    let line = stderr.lines().find(|l| l.contains(label))?;
    let rest = &line[line.find(label)? + label.len()..];
    rest.split_whitespace().next()?.parse().ok()
}

#[rstest]
fn consolidation_bounds_merge_sources_without_changing_the_output(
    #[values(
        "coordinate",
        "queryname::lexicographical",
        "queryname::natural",
        "template-coordinate"
    )]
    order: &str,
    // 2 and 4 give 2-wide merges; 8 gives merges of up to 4 runs.
    #[values("2", "4", "8")] limit: &str,
    #[values("1", "4")] threads: &str,
) {
    let dir = TempDir::new().expect("tempdir");
    let input = write_fixture(dir.path(), 20_000);

    assert_consolidation_preserves_output(&input, dir.path(), order, threads, limit, &[]);
}

/// Template-coordinate chooses its key width at runtime, so each lane is a
/// separate consolidation kernel. `--key-types` forces each one.
#[rstest]
fn consolidation_preserves_every_template_coordinate_key_lane(
    #[values("none", "cb", "mi", "full")] key_types: &str,
) {
    let dir = TempDir::new().expect("tempdir");
    let input = write_fixture(dir.path(), 20_000);
    assert_consolidation_preserves_output(
        &input,
        dir.path(),
        "template-coordinate",
        "4",
        "8",
        &["--key-types", key_types],
    );
}

/// Sort `input` without a binding limit and with `limit`, then assert the bounded
/// sort consolidated, kept its merge sources under the limit, and produced the
/// same records as the reference.
fn assert_consolidation_preserves_output(
    input: &Path,
    dir: &Path,
    order: &str,
    threads: &str,
    limit: &str,
    extra_args: &[&str],
) {
    let (reference, reference_log) =
        sort(input, dir, "reference", order, threads, UNBOUNDED, extra_args);
    let (bounded, bounded_log) = sort(input, dir, "bounded", order, threads, limit, extra_args);

    let limit_n: u64 = limit.parse().unwrap();
    let runs = logged_count(&bounded_log, "Spill runs:").expect("the fixture must spill");
    assert!(
        runs > limit_n,
        "{runs} spilled runs cannot exercise a limit of {limit}:\n{bounded_log}"
    );
    assert!(
        bounded_log.contains("Consolidating "),
        "{runs} runs under a limit of {limit} must consolidate:\n{bounded_log}"
    );
    let sources = logged_count(&bounded_log, "Merge sources:").expect("merge sources are logged");
    assert!(sources < limit_n, "{sources} merge sources exceed the limit of {limit}");
    assert!(
        !reference_log.contains("Consolidating "),
        "the unbounded reference must not consolidate:\n{reference_log}"
    );

    assert!(
        decompressed_records_without_pg(&bounded) == decompressed_records_without_pg(&reference),
        "consolidating under --max-temp-files {limit} changed the {order} output at {threads} \
         threads ({extra_args:?})"
    );
}

/// The limit exists to keep the sort inside the process's open-file budget. With
/// `ulimit -n` far below the number of spilled runs, the sort must still finish
/// when `--max-temp-files` is set below it — which it could not while every run
/// was held open until the merge.
#[cfg(unix)]
#[test]
fn a_sort_spilling_more_runs_than_the_descriptor_limit_completes() {
    // Low enough that the unbounded reference below (one descriptor per run)
    // still fits a common inherited soft limit of 256.
    const DESCRIPTOR_LIMIT: u64 = 32;
    let dir = TempDir::new().expect("tempdir");
    let input = write_fixture(dir.path(), 80_000);
    let output = dir.path().join("sorted.bam");

    // Count the runs without a descriptor limit first, so the assertion below
    // cannot pass vacuously on a fixture that spills too little.
    let (reference, log) = sort(&input, dir.path(), "count", "coordinate", "2", UNBOUNDED, &[]);
    let runs = logged_count(&log, "Spill runs:").expect("the fixture must spill");
    assert!(runs > 2 * DESCRIPTOR_LIMIT, "only {runs} runs; the descriptor limit would not bind");

    let script = format!(
        "ulimit -n {DESCRIPTOR_LIMIT} && exec \"$0\" sort -i \"$1\" -o \"$2\" --order coordinate \
         --threads 2 --max-temp-files 8 --max-memory 64K --memory-per-thread false"
    );
    let result = Command::new("sh")
        .args(["-c", &script, env!("CARGO_BIN_EXE_fgumi")])
        .arg(&input)
        .arg(&output)
        .env("RUST_LOG", "info")
        .output()
        .expect("run fgumi sort under a descriptor limit");
    let stderr = String::from_utf8_lossy(&result.stderr);
    assert!(result.status.success(), "sort failed under ulimit -n {DESCRIPTOR_LIMIT}:\n{stderr}");
    // Finishing is not enough: the bounded sort must produce the same records.
    assert!(
        decompressed_records_without_pg(&output) == decompressed_records_without_pg(&reference),
        "the sort under ulimit -n {DESCRIPTOR_LIMIT} changed the output"
    );
}

/// In a fused chain the spill writer is a pool-scheduled step rather than the
/// standalone sort's dedicated thread; consolidation must behave the same there.
#[test]
fn consolidation_in_a_fused_runall_chain_preserves_the_output() {
    let dir = TempDir::new().expect("tempdir");
    let input = dir.path().join("unsorted.bam");
    write_bam(&input, &create_minimal_header("chr1", 1_000_000), &uniquely_named_records(20_000));
    let run = |name: &str, limit: &str| -> (PathBuf, String) {
        let output = dir.path().join(format!("{name}.bam"));
        let result = Command::new(env!("CARGO_BIN_EXE_fgumi"))
            .env("RUST_LOG", "info")
            .args(["runall", "--start-from", "sort", "--stop-after", "group", "-i"])
            .arg(&input)
            .arg("-o")
            .arg(&output)
            .args(["--group::strategy", "identity", "--group::edits", "0", "--threads", "4"])
            .args(["--sort::max-temp-files", limit, "--sort::max-memory", "64K"])
            .args(["--memory-per-thread", "false"])
            .output()
            .expect("run fgumi runall");
        let stderr = String::from_utf8_lossy(&result.stderr).into_owned();
        assert!(result.status.success(), "fgumi runall failed:\n{stderr}");
        (output, stderr)
    };
    let (reference, reference_log) = run("reference", UNBOUNDED);
    let (bounded, bounded_log) = run("bounded", "2");
    assert!(
        bounded_log.contains("Consolidating "),
        "the bounded runall sort must consolidate:\n{bounded_log}"
    );
    assert!(!reference_log.contains("Consolidating "), "the reference must not consolidate");
    assert!(
        decompressed_records_without_pg(&bounded) == decompressed_records_without_pg(&reference),
        "consolidating inside runall changed the sort -> group output"
    );
}
