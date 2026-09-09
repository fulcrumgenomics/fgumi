//! Worker-count cutover parity for inline consensus metrics.
//!
//! Follows this repo's existing `test_*_cutover_parity.rs` precedent (e.g.
//! `crc_policy_parity_serial_and_chain` in `test_copy_umi_command.rs`): run
//! the *same* input through the *same* command at two different
//! `--threads` values and assert the metrics output is byte-for-byte
//! identical. Unlike `copy-umi`'s serial-oracle-vs-chain-builder split,
//! `simplex`/`duplex`/`codec` route every invocation through the same
//! chain/`Pipeline` builder regardless of thread count (Task 8's
//! `CoordinateGroupCollector` is wired in as a step, not a separate serial
//! engine) — so this is a pure worker-count determinism check: 1 worker vs.
//! 8 workers must derive the identical coordinate-group boundaries and
//! reducer totals from Task 8's `CoordinateGroupCollector` (T2) and Task 11's
//! T1 fused adapter alike.
//!
//! Uses a 60-distinct-MI single coordinate group so the run actually spans
//! `GroupByMi`'s 50-MI-group batch boundary (`DEFAULT_TARGET_BATCH_COUNT`) —
//! the shape most likely to expose a worker-count-dependent batching bug,
//! rather than a trivial single-family case where thread count cannot
//! matter.

use crate::helpers::bam_generator::{create_minimal_header, write_bam};
use fgumi_raw_bam::{RawRecord, SamBuilder, SamTag, flags};
use rstest::rstest;
use std::path::{Path, PathBuf};
use std::process::Command;

fn run_fgumi(args: &[&str]) {
    let status = Command::new(env!("CARGO_BIN_EXE_fgumi"))
        .args(args)
        .status()
        .unwrap_or_else(|e| panic!("failed to spawn fgumi {args:?}: {e}"));
    assert!(status.success(), "fgumi {args:?} failed with {status}");
}

fn suffixed(prefix: &Path, suffix: &str) -> PathBuf {
    PathBuf::from(format!("{}.{suffix}", prefix.display()))
}

fn assert_metrics_file_eq(a_prefix: &Path, b_prefix: &Path, suffix: &str, context: &str) {
    let a = std::fs::read_to_string(suffixed(a_prefix, suffix))
        .unwrap_or_else(|e| panic!("{context}: missing {suffix} at {}: {e}", a_prefix.display()));
    let b = std::fs::read_to_string(suffixed(b_prefix, suffix))
        .unwrap_or_else(|e| panic!("{context}: missing {suffix} at {}: {e}", b_prefix.display()));
    assert_eq!(a, b, "{context}: {suffix} differs across thread counts");
}

/// Base-4 (`ACGT`) encoding of `i`, unique for any `i < 4^8`. `--edits 0` is
/// used on every `group` call below, so uniqueness (not edit-distance
/// separation) is all that is required.
fn indexed_umi(i: usize) -> String {
    const BASES: [u8; 4] = [b'A', b'C', b'G', b'T'];
    let mut n = i;
    let mut out = [0u8; 8];
    for slot in out.iter_mut().rev() {
        *slot = BASES[n % 4];
        n /= 4;
    }
    String::from_utf8(out.to_vec()).expect("ACGT bytes are valid utf8")
}

const N_MIS: usize = 60;

/// `N_MIS` distinct single-index UMIs, one single-end read each, all at one
/// position — one coordinate group split across `N_MIS` MIs, straddling
/// `GroupByMi`'s 50-MI batch boundary.
fn simplex_many_mi_records() -> Vec<RawRecord> {
    (0..N_MIS)
        .map(|f| {
            let umi = indexed_umi(f);
            let mut b = SamBuilder::new();
            b.read_name(format!("f{f}").as_bytes())
                .sequence(b"ACGTACGTAC")
                .qualities(&[30; 10])
                .flags(0)
                .ref_id(0)
                .pos(100)
                .mapq(60)
                .cigar_ops(&[10 << 4]);
            b.add_string_tag(SamTag::RX, umi.as_bytes());
            b.build()
        })
        .collect()
}

/// One FR read pair sharing a dual UMI, no `MI` tag.
fn duplex_pair(name: &str, umi: &str, r1_pos: i32, r2_pos: i32) -> (RawRecord, RawRecord) {
    let seq = vec![b'A'; 10];
    let mut b1 = SamBuilder::new();
    b1.read_name(name.as_bytes())
        .sequence(&seq)
        .qualities(&[30u8; 10])
        .flags(flags::PAIRED | flags::FIRST_SEGMENT | flags::MATE_REVERSE)
        .ref_id(0)
        .pos(r1_pos)
        .mapq(60)
        .cigar_ops(&[10 << 4])
        .mate_ref_id(0)
        .mate_pos(r2_pos)
        .template_length(r2_pos + 10 - r1_pos);
    b1.add_string_tag(SamTag::RX, umi.as_bytes());
    b1.add_string_tag(SamTag::MC, b"10M");
    let r1 = b1.build();

    let mut b2 = SamBuilder::new();
    b2.read_name(name.as_bytes())
        .sequence(&seq)
        .qualities(&[30u8; 10])
        .flags(flags::PAIRED | flags::LAST_SEGMENT | flags::REVERSE)
        .ref_id(0)
        .pos(r2_pos)
        .mapq(60)
        .cigar_ops(&[10 << 4])
        .mate_ref_id(0)
        .mate_pos(r1_pos)
        .template_length(-(r2_pos + 10 - r1_pos));
    b2.add_string_tag(SamTag::RX, umi.as_bytes());
    b2.add_string_tag(SamTag::MC, b"10M");
    let r2 = b2.build();

    (r1, r2)
}

/// `N_MIS` distinct single-strand duplex molecules (dual UMI, one strand
/// each), all at one position — same shape as
/// [`simplex_many_mi_records`] but paired-end, for the duplex/codec cases.
fn duplex_many_mi_records() -> Vec<RawRecord> {
    let mut records = Vec::new();
    for f in 0..N_MIS {
        let umi = format!("{}-{}", indexed_umi(f), indexed_umi(f + 10_000));
        let (r1, r2) = duplex_pair(&format!("f{f}"), &umi, 100, 200);
        records.push(r1);
        records.push(r2);
    }
    records
}

// ---------------------------------------------------------------------------
// T1 (fused) worker-count cutover
// ---------------------------------------------------------------------------

fn run_simplex_fused(dir: &Path, input: &Path, tag: &str, threads: usize) -> PathBuf {
    let prefix = dir.join(format!("fused_{tag}"));
    let out = dir.join(format!("fused_{tag}_out.bam"));
    let threads_arg = threads.to_string();
    run_fgumi(&[
        "runall",
        "-i",
        input.to_str().unwrap(),
        "-o",
        out.to_str().unwrap(),
        "--start-from",
        "group",
        "--stop-after",
        "consensus",
        "--consensus",
        "simplex",
        "--group::strategy",
        "adjacency",
        "--group::edits",
        "0",
        "--simplex::min-reads",
        "1",
        "--simplex::metrics",
        prefix.to_str().unwrap(),
        "--threads",
        &threads_arg,
    ]);
    prefix
}

fn run_duplex_fused(dir: &Path, input: &Path, tag: &str, threads: usize) -> PathBuf {
    let prefix = dir.join(format!("fused_{tag}"));
    let out = dir.join(format!("fused_{tag}_out.bam"));
    let threads_arg = threads.to_string();
    run_fgumi(&[
        "runall",
        "-i",
        input.to_str().unwrap(),
        "-o",
        out.to_str().unwrap(),
        "--start-from",
        "group",
        "--stop-after",
        "consensus",
        "--consensus",
        "duplex",
        "--group::strategy",
        "paired",
        "--group::edits",
        "0",
        "--duplex::min-reads",
        "1",
        "--duplex::metrics",
        prefix.to_str().unwrap(),
        "--threads",
        &threads_arg,
    ]);
    prefix
}

fn run_codec_fused(dir: &Path, input: &Path, tag: &str, threads: usize) -> PathBuf {
    let prefix = dir.join(format!("fused_{tag}"));
    let out = dir.join(format!("fused_{tag}_out.bam"));
    let threads_arg = threads.to_string();
    run_fgumi(&[
        "runall",
        "-i",
        input.to_str().unwrap(),
        "-o",
        out.to_str().unwrap(),
        "--start-from",
        "group",
        "--stop-after",
        "consensus",
        "--consensus",
        "codec",
        "--group::strategy",
        "adjacency",
        "--group::edits",
        "0",
        "--codec::min-reads",
        "1",
        "--codec::metrics",
        prefix.to_str().unwrap(),
        "--threads",
        &threads_arg,
    ]);
    prefix
}

const SIMPLEX_SUFFIXES: [&str; 3] =
    ["family_sizes.txt", "umi_counts.txt", "simplex_yield_metrics.txt"];
const DUPLEX_SUFFIXES: [&str; 4] =
    ["family_sizes.txt", "duplex_family_sizes.txt", "umi_counts.txt", "duplex_yield_metrics.txt"];

#[test]
fn simplex_fused_metrics_match_across_thread_counts() {
    let dir = tempfile::tempdir().expect("tempdir");
    let header = create_minimal_header("chr1", 10_000);
    let bam_path = dir.path().join("in.bam");
    write_bam(&bam_path, &header, &simplex_many_mi_records());

    let one = run_simplex_fused(dir.path(), &bam_path, "t1", 1);
    let eight = run_simplex_fused(dir.path(), &bam_path, "t8", 8);

    for suffix in SIMPLEX_SUFFIXES {
        assert_metrics_file_eq(&one, &eight, suffix, "simplex fused: threads=1 vs threads=8");
    }
}

#[test]
fn duplex_fused_metrics_match_across_thread_counts() {
    let dir = tempfile::tempdir().expect("tempdir");
    let header = create_minimal_header("chr1", 10_000);
    let bam_path = dir.path().join("in.bam");
    write_bam(&bam_path, &header, &duplex_many_mi_records());

    let one = run_duplex_fused(dir.path(), &bam_path, "t1", 1);
    let eight = run_duplex_fused(dir.path(), &bam_path, "t8", 8);

    for suffix in DUPLEX_SUFFIXES {
        assert_metrics_file_eq(&one, &eight, suffix, "duplex fused: threads=1 vs threads=8");
    }
}

#[test]
fn codec_fused_metrics_match_across_thread_counts() {
    let dir = tempfile::tempdir().expect("tempdir");
    let header = create_minimal_header("chr1", 10_000);
    let bam_path = dir.path().join("in.bam");
    write_bam(&bam_path, &header, &duplex_many_mi_records());

    let one = run_codec_fused(dir.path(), &bam_path, "t1", 1);
    let eight = run_codec_fused(dir.path(), &bam_path, "t8", 8);

    for suffix in DUPLEX_SUFFIXES {
        assert_metrics_file_eq(&one, &eight, suffix, "codec fused: threads=1 vs threads=8");
    }
}

// ---------------------------------------------------------------------------
// T2 (standalone) worker-count cutover — the collector design's own
// correctness anchor, per the task brief.
// ---------------------------------------------------------------------------

fn run_simplex_standalone(dir: &Path, grouped: &Path, tag: &str, threads: usize) -> PathBuf {
    let prefix = dir.join(format!("standalone_{tag}"));
    let out = dir.join(format!("standalone_{tag}_out.bam"));
    let threads_arg = threads.to_string();
    run_fgumi(&[
        "simplex",
        "-i",
        grouped.to_str().unwrap(),
        "-o",
        out.to_str().unwrap(),
        "--min-reads",
        "1",
        "--metrics",
        prefix.to_str().unwrap(),
        "--threads",
        &threads_arg,
    ]);
    prefix
}

fn run_duplex_standalone(dir: &Path, grouped: &Path, tag: &str, threads: usize) -> PathBuf {
    let prefix = dir.join(format!("standalone_{tag}"));
    let out = dir.join(format!("standalone_{tag}_out.bam"));
    let threads_arg = threads.to_string();
    run_fgumi(&[
        "duplex",
        "-i",
        grouped.to_str().unwrap(),
        "-o",
        out.to_str().unwrap(),
        "--min-reads",
        "1",
        "--metrics",
        prefix.to_str().unwrap(),
        "--threads",
        &threads_arg,
    ]);
    prefix
}

fn run_codec_standalone(dir: &Path, grouped: &Path, tag: &str, threads: usize) -> PathBuf {
    let prefix = dir.join(format!("standalone_{tag}"));
    let out = dir.join(format!("standalone_{tag}_out.bam"));
    let threads_arg = threads.to_string();
    run_fgumi(&[
        "codec",
        "-i",
        grouped.to_str().unwrap(),
        "-o",
        out.to_str().unwrap(),
        "--min-reads",
        "1",
        "--metrics",
        prefix.to_str().unwrap(),
        "--threads",
        &threads_arg,
    ]);
    prefix
}

#[rstest]
#[case::simplex(2)]
#[case::simplex_wide(8)]
fn simplex_standalone_metrics_match_across_thread_counts(#[case] wide: usize) {
    let dir = tempfile::tempdir().expect("tempdir");
    let header = create_minimal_header("chr1", 10_000);
    let bam_path = dir.path().join("in.bam");
    write_bam(&bam_path, &header, &simplex_many_mi_records());

    let grouped = dir.path().join("grouped.bam");
    run_fgumi(&[
        "group",
        "-i",
        bam_path.to_str().unwrap(),
        "-o",
        grouped.to_str().unwrap(),
        "-s",
        "adjacency",
        "-e",
        "0",
    ]);

    let one = run_simplex_standalone(dir.path(), &grouped, "t1", 1);
    let other = run_simplex_standalone(dir.path(), &grouped, "twide", wide);

    for suffix in SIMPLEX_SUFFIXES {
        assert_metrics_file_eq(
            &one,
            &other,
            suffix,
            &format!("simplex standalone: threads=1 vs threads={wide}"),
        );
    }
}

#[test]
fn duplex_standalone_metrics_match_across_thread_counts() {
    let dir = tempfile::tempdir().expect("tempdir");
    let header = create_minimal_header("chr1", 10_000);
    let bam_path = dir.path().join("in.bam");
    write_bam(&bam_path, &header, &duplex_many_mi_records());

    let grouped = dir.path().join("grouped.bam");
    run_fgumi(&[
        "group",
        "-i",
        bam_path.to_str().unwrap(),
        "-o",
        grouped.to_str().unwrap(),
        "-s",
        "paired",
        "-e",
        "0",
    ]);

    let one = run_duplex_standalone(dir.path(), &grouped, "t1", 1);
    let eight = run_duplex_standalone(dir.path(), &grouped, "t8", 8);

    for suffix in DUPLEX_SUFFIXES {
        assert_metrics_file_eq(&one, &eight, suffix, "duplex standalone: threads=1 vs threads=8");
    }
}

#[test]
fn codec_standalone_metrics_match_across_thread_counts() {
    let dir = tempfile::tempdir().expect("tempdir");
    let header = create_minimal_header("chr1", 10_000);
    let bam_path = dir.path().join("in.bam");
    write_bam(&bam_path, &header, &duplex_many_mi_records());

    let grouped = dir.path().join("grouped.bam");
    run_fgumi(&[
        "group",
        "-i",
        bam_path.to_str().unwrap(),
        "-o",
        grouped.to_str().unwrap(),
        "-s",
        "adjacency",
        "-e",
        "0",
    ]);

    let one = run_codec_standalone(dir.path(), &grouped, "t1", 1);
    let eight = run_codec_standalone(dir.path(), &grouped, "t8", 8);

    for suffix in DUPLEX_SUFFIXES {
        assert_metrics_file_eq(&one, &eight, suffix, "codec standalone: threads=1 vs threads=8");
    }
}
