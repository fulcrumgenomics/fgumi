//! End-to-end checks of the range-parallel final merge in `fgumi sort`.
//!
//! The range merge must be byte-identical to the serial merge of the same
//! binary (BAM and `.bai`), must fall back to the serial merge, with the
//! documented reason, when it cannot run, and must be counted in the
//! `=== Sort Phase Timing ===` k-way merge row.

use std::ffi::OsStr;
use std::fmt::Write as _;
use std::path::Path;
use std::process::{Command, Output};

use fgumi_lib::sam::SamTag;
use fgumi_raw_bam::{RawRecord, SamBuilder, flags};
use rstest::rstest;
use tempfile::TempDir;

use crate::helpers::bam_generator::{create_minimal_header, write_bam};

const FGUMI: &str = env!("CARGO_BIN_EXE_fgumi");

/// `num_pairs` unsorted, paired, mapped read pairs on one reference (positions
/// scramble within the reference) plus a few unmapped pairs.
fn unsorted_pairs(num_pairs: usize) -> Vec<RawRecord> {
    let seq = b"ACGTACGTAC";
    let qual = [30u8; 10];
    let cigar_op = 10u32 << 4; // 10M
    let mut records = Vec::with_capacity(num_pairs * 2 + 8);
    for i in 0..num_pairs {
        // Multiplying by a prime co-prime with the modulus scrambles positions.
        let pos = i32::try_from((i * 7919) % 9973 * 50).expect("pos fits i32");
        let mate_pos = pos + 300;
        let name = format!("read{i}");
        let mut r1 = SamBuilder::new();
        r1.read_name(name.as_bytes())
            .ref_id(0)
            .pos(pos)
            .mapq(60)
            .flags(flags::PAIRED | flags::FIRST_SEGMENT | flags::MATE_REVERSE)
            .cigar_ops(&[cigar_op])
            .sequence(seq)
            .qualities(&qual)
            .mate_ref_id(0)
            .mate_pos(mate_pos)
            .template_length(310)
            .add_string_tag(SamTag::MC, b"10M");
        records.push(r1.build());
        let mut r2 = SamBuilder::new();
        r2.read_name(name.as_bytes())
            .ref_id(0)
            .pos(mate_pos)
            .mapq(60)
            .flags(flags::PAIRED | flags::LAST_SEGMENT | flags::REVERSE)
            .cigar_ops(&[cigar_op])
            .sequence(seq)
            .qualities(&qual)
            .mate_ref_id(0)
            .mate_pos(pos)
            .template_length(-310)
            .add_string_tag(SamTag::MC, b"10M");
        records.push(r2.build());
    }
    for i in 0..4 {
        for segment in [flags::FIRST_SEGMENT, flags::LAST_SEGMENT] {
            let mut b = SamBuilder::new();
            b.read_name(format!("unmapped{i}").as_bytes())
                .flags(flags::PAIRED | flags::UNMAPPED | flags::MATE_UNMAPPED | segment)
                .sequence(seq)
                .qualities(&qual);
            records.push(b.build());
        }
    }
    records
}

/// A generated input BAM of `num_pairs` pairs in its own temp dir.
fn input_bam(dir: &TempDir, num_pairs: usize) -> std::path::PathBuf {
    let path = dir.path().join("in.bam");
    write_bam(&path, &create_minimal_header("chr1", 1_000_000), &unsorted_pairs(num_pairs));
    path
}

/// Input pairs and `--max-memory` (per thread) for the range-mode runs: a
/// few spill runs (2 at `-@ 8`, 5 or 6 at `-@ 4`), so the window fits the
/// 64 MiB read-ahead floor, and about 10 MB of run data, so the planned
/// range target (a few MiB at that budget) still cuts several ranges.
const RANGE_PAIRS: usize = 60_000;
const RANGE_MEMORY: &str = "512K";

/// `--max-memory` per thread for a range-mode run at `threads`: the same
/// 2 MiB in total below `-@ 4` (so as few runs, and as many ranges, as at
/// `-@ 4`), else [`RANGE_MEMORY`].
fn range_memory(threads: usize) -> &'static str {
    match threads {
        1 => "2M",
        2 => "1M",
        _ => RANGE_MEMORY,
    }
}

/// `--max-memory` for the unmapped-tail runs: their input spills 8 runs at
/// `RANGE_MEMORY` and `-@ 8`, where the phase-1 front's worst-case estimate
/// (`S_c = S`, each range spanning its target plus one boundary frame per
/// run) does not fit the 64 MiB floor, so the final run would stay resident
/// and the merge serial. At 1 MiB per thread it spills fewer runs and fits.
const UNMAPPED_TAIL_MEMORY: &str = "1M";

/// Input pairs and `--max-memory` for the fallback runs.
const FALLBACK_PAIRS: usize = 20_000;
const MEMORY: &str = "64K";

/// `--max-memory` that spills too many runs (10 or more at `-@ 8`) for any
/// range window to fit the 64 MiB read-ahead floor.
const TOO_MANY_RUNS_MEMORY: &str = "16K";

/// Runs a spilled `fgumi sort --sort-stats --max-memory {memory}` writing to
/// `output`; asserts success.
fn run_sort(
    input: &Path,
    output: &Path,
    tmp: &Path,
    order: &str,
    memory: &str,
    extra: &[&str],
    envs: &[(&str, &str)],
) -> String {
    let mut cmd = Command::new(FGUMI);
    cmd.env("RUST_LOG", "info").envs(envs.iter().copied());
    cmd.args([OsStr::new("sort"), OsStr::new("-i"), input.as_os_str()]);
    cmd.args([OsStr::new("-o"), output.as_os_str()]);
    cmd.args(["--order", order, "--max-memory", memory, "--sort-stats"]);
    cmd.args([OsStr::new("--tmp-dir"), tmp.as_os_str()]);
    cmd.args(extra);
    let out: Output = cmd.output().expect("failed to spawn fgumi sort");
    let stderr = String::from_utf8_lossy(&out.stderr).into_owned();
    assert!(out.status.success(), "sort {order} {extra:?} failed:\n{stderr}");
    stderr
}

fn verify_sorted(path: &Path, order: &str) {
    let status = Command::new(FGUMI)
        .args([OsStr::new("sort"), OsStr::new("--verify"), OsStr::new("-i")])
        .arg(path)
        .args(["--order", order])
        .status()
        .expect("failed to spawn fgumi sort --verify");
    assert!(status.success(), "{} is not {order}-sorted", path.display());
}

/// The `Spill runs: N` summary count in `stderr` (the run-formation line
/// `Spill runs: R from C chunks` is not a bare count and is skipped).
fn spill_runs(stderr: &str) -> usize {
    stderr
        .lines()
        .find_map(|l| l.split_once("Spill runs: ").and_then(|(_, n)| n.trim().parse().ok()))
        .unwrap_or_else(|| panic!("no `Spill runs: N` line:\n{stderr}"))
}

/// The range count of the `Range merge: N ranges` line in `stderr`.
fn range_count(stderr: &str) -> usize {
    stderr
        .split("Range merge: ")
        .nth(1)
        .and_then(|rest| rest.split_whitespace().next())
        .and_then(|n| n.parse().ok())
        .unwrap_or_else(|| panic!("no `Range merge:` line:\n{stderr}"))
}

fn bai(path: &Path) -> std::path::PathBuf {
    let mut p = path.as_os_str().to_owned();
    p.push(".bai");
    p.into()
}

/// The range merge's BAM (and `.bai`) equal the serial merge's at every
/// pool width, with merge-thread compression enabled. How many blocks the
/// merge thread compresses depends on scheduling (none at one thread in local
/// runs, a varying share above), so this asserts byte identity with
/// compression on, not that any block was compressed; the unit test
/// `merge_thread_compression_is_byte_identical` (in `fgumi-pipeline-io`)
/// forces refused blocks and asserts the merge thread compressed some.
#[rstest]
#[case::coordinate("coordinate", true)]
#[case::template_coordinate("template-coordinate", false)]
fn range_merge_matches_serial_merge(
    #[case] order: &str,
    #[case] write_index: bool,
    #[values(1, 2, 4, 8)] threads: usize,
) {
    let dir = TempDir::new().expect("temp dir");
    let input = input_bam(&dir, RANGE_PAIRS);
    let tmp = TempDir::new().expect("tmp dir");
    let out = dir.path().join("out.bam");
    let memory = range_memory(threads);
    let threads = threads.to_string();
    let mut extra = vec!["-@", threads.as_str()];
    if write_index {
        extra.push("--write-index");
    }

    let range_stderr = run_sort(&input, &out, tmp.path(), order, memory, &extra, &[]);
    assert!(range_stderr.contains("Sort merge mode: range"), "{range_stderr}");
    assert!(range_stderr.contains("index stride 1,"), "{range_stderr}");
    assert!(range_stderr.contains("merge-thread compression on)"), "{range_stderr}");
    let ranges = range_count(&range_stderr);
    assert!(ranges > 1, "expected several key ranges, got {ranges}:\n{range_stderr}");
    verify_sorted(&out, order);
    let range_bam = std::fs::read(&out).expect("read range output");
    let range_bai = write_index.then(|| std::fs::read(bai(&out)).expect("read range bai"));

    let serial_stderr = run_sort(
        &input,
        &out,
        tmp.path(),
        order,
        memory,
        &extra,
        &[("FGUMI_SORT_SERIAL_MERGE", "1")],
    );
    assert!(
        serial_stderr.contains("Sort merge mode: serial (forced by FGUMI_SORT_SERIAL_MERGE)"),
        "{serial_stderr}"
    );
    // The range merge spills the final run the serial merge keeps in memory.
    assert_eq!(
        spill_runs(&range_stderr),
        spill_runs(&serial_stderr) + 1,
        "range:\n{range_stderr}\nserial:\n{serial_stderr}"
    );
    let serial_bam = std::fs::read(&out).expect("read serial output");
    assert!(range_bam == serial_bam, "range and serial BAMs differ byte for byte");
    if let Some(range_bai) = range_bai {
        let serial_bai = std::fs::read(bai(&out)).expect("read serial bai");
        assert!(range_bai == serial_bai, "range and serial .bai differ byte for byte");
    }
}

/// `mapped_pairs` pairs of [`unsorted_pairs`] interleaved with
/// `unmapped_pairs` fully unmapped pairs, so every spill run holds part of
/// the unmapped tail (one key in coordinate order). Unmapped reads are
/// 150 bp, so the tail spans many range targets.
fn unmapped_heavy_records(mapped_pairs: usize, unmapped_pairs: usize) -> Vec<RawRecord> {
    let seq = b"ACGTACGTAC".repeat(15);
    let qual = [30u8; 150];
    let mapped = unsorted_pairs(mapped_pairs);
    let mut unmapped = Vec::with_capacity(unmapped_pairs * 2);
    for i in 0..unmapped_pairs {
        for segment in [flags::FIRST_SEGMENT, flags::LAST_SEGMENT] {
            let mut b = SamBuilder::new();
            b.read_name(format!("tail{i}").as_bytes())
                .flags(flags::PAIRED | flags::UNMAPPED | flags::MATE_UNMAPPED | segment)
                .sequence(&seq)
                .qualities(&qual);
            unmapped.push(b.build());
        }
    }
    let per_mapped = unmapped.len().div_ceil(mapped.len().max(1));
    let mut unmapped = unmapped.into_iter();
    let mut records = Vec::with_capacity(mapped.len() + unmapped_pairs * 2);
    for r in mapped {
        records.push(r);
        records.extend(unmapped.by_ref().take(per_mapped));
    }
    records.extend(unmapped);
    records
}

/// The `largest range N% of records` figure of the `Range merge:` line.
fn largest_range_pct(stderr: &str) -> f64 {
    stderr
        .split("largest range ")
        .nth(1)
        .and_then(|rest| rest.split('%').next())
        .and_then(|n| n.parse().ok())
        .unwrap_or_else(|| panic!("no `largest range` figure:\n{stderr}"))
}

/// A coordinate sort whose unmapped tail (every record one key) is five
/// sixths of the records is cut into many ranges, bounds falling inside the
/// tail, so no range holds the whole tail (one range held all of it when
/// bounds fell only between distinct keys); the output is byte-identical to
/// the serial merge. Template order has no such tie and must stay identical.
#[rstest]
#[case::coordinate("coordinate")]
#[case::template_coordinate("template-coordinate")]
fn range_merge_splits_an_unmapped_tail(#[case] order: &str) {
    let dir = TempDir::new().expect("temp dir");
    let input = dir.path().join("in.bam");
    let records = unmapped_heavy_records(RANGE_PAIRS / 6, 5 * RANGE_PAIRS / 6);
    write_bam(&input, &create_minimal_header("chr1", 1_000_000), &records);
    let tmp = TempDir::new().expect("tmp dir");
    let out = dir.path().join("out.bam");
    let extra = ["-@", "8"];

    let range_stderr = run_sort(&input, &out, tmp.path(), order, UNMAPPED_TAIL_MEMORY, &extra, &[]);
    assert!(range_stderr.contains("Sort merge mode: range"), "{range_stderr}");
    assert!(range_count(&range_stderr) > 1, "{range_stderr}");
    let largest = largest_range_pct(&range_stderr);
    assert!(largest < 30.0, "largest range {largest}% of records:\n{range_stderr}");
    verify_sorted(&out, order);
    let range_bam = std::fs::read(&out).expect("read range output");

    run_sort(
        &input,
        &out,
        tmp.path(),
        order,
        UNMAPPED_TAIL_MEMORY,
        &extra,
        &[("FGUMI_SORT_SERIAL_MERGE", "1")],
    );
    assert!(
        range_bam == std::fs::read(&out).expect("read forced-serial output"),
        "range output differs from the forced-serial output"
    );
}

/// Each fallback merges serially, with output byte-identical to the
/// forced-serial merge, and its runs are the forced-serial (`main`) layout
/// exactly: the run indexes live in the runs' footers, never in the sort
/// buffer, so indexing (or turning it off part-way) does not move a run
/// boundary. None of these cases consolidates (the default temp-file limit is
/// not reached). Columns: order, memory, extra arguments, the serial reason
/// (the whole reason, or for a reason that names budgets its distinguishing
/// tail), whether the front indexed runs before turning indexing off.
#[rstest]
#[case::bgzf_spill(
    "coordinate",
    MEMORY,
    &["-@", "8", "--temp-codec", "bgzf"],
    "Sort merge mode: serial (spill run 0 has no run index)",
    false
)]
#[case::queryname(
    "queryname",
    MEMORY,
    &["-@", "8"],
    "Sort merge mode: serial (spill run 0 has no run index)",
    false
)]
// Ten or more spill runs at `-@ 8`: even one range per worker at decode
// amplification 1.25 needs 8 x k x 128 KiB of running frames plus
// 8 x 3 x (4 x k x 64 KiB), over the range budget (the larger of the 64 MiB
// read-ahead floor `R` and the merge memory). The front's estimate finds no
// fit once enough runs spilled and turns indexing off, so the final run stays
// resident, which the merge reports (the writer leaves the footers unread).
#[case::indexes_dropped_when_too_many_runs(
    "coordinate",
    TOO_MANY_RUNS_MEMORY,
    &["-@", "8"],
    "(range merge needs the final run spilled but",
    true
)]
// The gauntlet's reproduction (278 runs against 236 on `main` when the
// index was charged to the buffer): the layout is now `main`'s exactly.
#[case::template_tiny_memory_matches_main(
    "template-coordinate",
    "32K",
    &["-@", "8"],
    "(range merge needs the final run spilled but",
    true
)]
fn serial_fallbacks(
    #[case] order: &str,
    #[case] memory: &str,
    #[case] extra: &[&str],
    #[case] reason: &str,
    #[case] indexed_then_dropped: bool,
) {
    let dir = TempDir::new().expect("temp dir");
    let input = input_bam(&dir, FALLBACK_PAIRS);
    let tmp = TempDir::new().expect("tmp dir");
    let out = dir.path().join("out.bam");

    let stderr = run_sort(&input, &out, tmp.path(), order, memory, extra, &[]);
    assert!(
        stderr.contains("Sort merge mode: serial (") && stderr.contains(reason),
        "expected serial fallback with {reason:?}:\n{stderr}"
    );
    assert!(!stderr.contains("Consolidating"), "{stderr}");
    verify_sorted(&out, order);
    let fallback_bam = std::fs::read(&out).expect("read fallback output");

    let serial_stderr = run_sort(
        &input,
        &out,
        tmp.path(),
        order,
        memory,
        extra,
        &[("FGUMI_SORT_SERIAL_MERGE", "1")],
    );
    assert_eq!(spill_runs(&stderr), spill_runs(&serial_stderr), "{stderr}\n{serial_stderr}");
    assert_eq!(stderr.contains("Run indexes dropped after"), indexed_then_dropped, "{stderr}");
    assert!(
        fallback_bam == std::fs::read(&out).expect("read forced-serial output"),
        "fallback output differs from the forced-serial output"
    );
}

/// Consolidated runs are indexed, so a sort that consolidates (too many temp
/// files) still merges by range, with output identical to the serial merge.
#[test]
fn range_merge_after_consolidation() {
    let dir = TempDir::new().expect("temp dir");
    let input = input_bam(&dir, FALLBACK_PAIRS);
    let tmp = TempDir::new().expect("tmp dir");
    let out = dir.path().join("out.bam");
    let extra = ["-@", "8", "--max-temp-files", "2"];

    let stderr = run_sort(&input, &out, tmp.path(), "coordinate", MEMORY, &extra, &[]);
    assert!(stderr.contains("Consolidating"), "{stderr}");
    assert!(stderr.contains("Sort merge mode: range"), "{stderr}");
    assert!(range_count(&stderr) > 1, "{stderr}");
    verify_sorted(&out, "coordinate");
    let range_bam = std::fs::read(&out).expect("read range output");

    run_sort(
        &input,
        &out,
        tmp.path(),
        "coordinate",
        MEMORY,
        &extra,
        &[("FGUMI_SORT_SERIAL_MERGE", "1")],
    );
    assert!(
        range_bam == std::fs::read(&out).expect("read forced-serial output"),
        "range output differs from the forced-serial output"
    );
}

#[test]
fn range_merge_counts_in_kway_phase_timing() {
    let dir = TempDir::new().expect("temp dir");
    let input = input_bam(&dir, RANGE_PAIRS);
    let tmp = TempDir::new().expect("tmp dir");
    let out = dir.path().join("out.bam");

    let stderr = run_sort(
        &input,
        &out,
        tmp.path(),
        "coordinate",
        RANGE_MEMORY,
        &["-@", "8"],
        &[("FGUMI_PIPELINE_STATS", "1")],
    );
    assert!(stderr.contains("Sort merge mode: range"), "{stderr}");
    assert!(stderr.contains("=== Sort Phase Timing ==="), "{stderr}");
    let row = stderr.lines().find(|l| l.contains("k-way merge")).expect("k-way merge row");
    let secs: f64 = row
        .split_whitespace()
        .find_map(|t| t.strip_suffix('s').and_then(|n| n.parse().ok()))
        .unwrap_or_else(|| panic!("no seconds in row {row:?}"));
    assert!(secs > 0.0, "k-way merge row is zero in range mode: {row:?}");
}

/// A SAM rendering of [`unsorted_pairs`]`(num_pairs)`.
fn unsorted_pairs_sam(num_pairs: usize) -> String {
    let mut sam = String::from("@HD\tVN:1.6\tSO:unsorted\n@SQ\tSN:chr1\tLN:1000000\n");
    for i in 0..num_pairs {
        let pos = (i * 7919) % 9973 * 50 + 1;
        let mate_pos = pos + 300;
        writeln!(
            sam,
            "read{i}\t97\tchr1\t{pos}\t60\t10M\t=\t{mate_pos}\t310\tACGTACGTAC\tIIIIIIIIII\tMC:Z:10M\n\
             read{i}\t145\tchr1\t{mate_pos}\t60\t10M\t=\t{pos}\t-310\tACGTACGTAC\tIIIIIIIIII\tMC:Z:10M"
        )
        .expect("write SAM");
    }
    for i in 0..4 {
        for flag in [77, 141] {
            writeln!(sam, "unmapped{i}\t{flag}\t*\t0\t0\t*\t*\t0\t0\tACGTACGTAC\tIIIIIIIIII")
                .expect("write SAM");
        }
    }
    sam
}

/// A SAM input sorts through `SortBuffer`, not the BAM block front; it still
/// spills its final run and merges by range, identically to the serial merge.
#[test]
fn sam_input_range_merges() {
    let dir = TempDir::new().expect("temp dir");
    let input = dir.path().join("in.sam");
    std::fs::write(&input, unsorted_pairs_sam(RANGE_PAIRS)).expect("write SAM");
    let tmp = TempDir::new().expect("tmp dir");
    let out = dir.path().join("out.bam");
    let extra = ["-@", "8"];

    let stderr = run_sort(&input, &out, tmp.path(), "coordinate", RANGE_MEMORY, &extra, &[]);
    assert!(spill_runs(&stderr) >= 2, "{stderr}");
    assert!(stderr.contains("Sort merge mode: range"), "{stderr}");
    assert!(range_count(&stderr) > 1, "{stderr}");
    verify_sorted(&out, "coordinate");
    let range_bam = std::fs::read(&out).expect("read range output");

    let serial_stderr = run_sort(
        &input,
        &out,
        tmp.path(),
        "coordinate",
        RANGE_MEMORY,
        &extra,
        &[("FGUMI_SORT_SERIAL_MERGE", "1")],
    );
    // `SortBuffer` spilled the final run the serial merge keeps in memory.
    assert_eq!(spill_runs(&stderr), spill_runs(&serial_stderr) + 1, "{stderr}\n{serial_stderr}");
    assert!(
        range_bam == std::fs::read(&out).expect("read forced-serial output"),
        "range output differs from the forced-serial output"
    );
}

/// Sorts the generated BAM piped through `-i /dev/stdin` to `out`.
fn run_sort_stream(input: &Path, out: &Path, tmp: &Path, envs: &[(&str, &str)]) -> String {
    let mut child = Command::new(FGUMI)
        .env("RUST_LOG", "info")
        .envs(envs.iter().copied())
        .args(["sort", "-i", "/dev/stdin", "-o"])
        .arg(out)
        .args(["--order", "coordinate", "--max-memory", RANGE_MEMORY, "--sort-stats", "-@", "8"])
        .arg("--tmp-dir")
        .arg(tmp)
        .stdin(std::fs::File::open(input).expect("open input"))
        .stderr(std::process::Stdio::piped())
        .spawn()
        .expect("spawn fgumi sort");
    let mut stderr = String::new();
    std::io::Read::read_to_string(&mut child.stderr.take().expect("stderr"), &mut stderr)
        .expect("read stderr");
    assert!(child.wait().expect("wait").success(), "{stderr}");
    stderr
}

/// A piped (non-seekable) BAM input still merges by range, with output
/// identical to the forced-serial merge of the same stream. (The file input's
/// `@PG` command line names different paths, so its bytes are not compared.)
#[test]
fn stream_input_range_merges() {
    let dir = TempDir::new().expect("temp dir");
    let input = input_bam(&dir, RANGE_PAIRS);
    let tmp = TempDir::new().expect("tmp dir");
    let out = dir.path().join("out.bam");

    let stderr = run_sort_stream(&input, &out, tmp.path(), &[]);
    assert!(stderr.contains("Sort merge mode: range"), "{stderr}");
    let range_bam = std::fs::read(&out).expect("read range output");
    run_sort_stream(&input, &out, tmp.path(), &[("FGUMI_SORT_SERIAL_MERGE", "1")]);
    assert!(
        range_bam == std::fs::read(&out).expect("read forced-serial output"),
        "stream range output differs from the forced-serial output"
    );
}
