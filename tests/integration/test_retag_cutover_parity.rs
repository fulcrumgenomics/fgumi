//! Parity gate for the `retag` command's single-threaded-path retirement (C4):
//! `Retag::execute` no longer has a serial in-process loop reached when
//! `--threads` is absent — it *always* routes through the declarative chain
//! builder. This test proves that cutover lost nothing user-observable.
//!
//! Two independent things are checked here:
//!
//! 1. **The cutover actually happened** (`retag_no_threads_routes_through_chain`).
//!    A no-`--threads` run now emits the chain-only `"Using pipeline with N
//!    threads"` banner that `ChainBuilder::add_retag` logs and the retired serial
//!    tail never did. This is the genuine RED/GREEN discriminator: before the
//!    removal a no-`--threads` run took the serial path and printed no such line;
//!    after it, the chain does.
//!
//! 2. **Output parity with the pre-removal serial path**
//!    (`cutover_matches_baseline`). The current build's `retag` output — records
//!    (byte-identical, modulo the `@PG` line) and the `--metrics` TSV — must match
//!    the frozen owned-serial baseline binary. The baseline path comes from
//!    `FGUMI_BASELINE_BIN`; when it is unset (or names a missing file) the case
//!    degrades to a self-consistency oracle (expected tag rewrites, record
//!    survival, and metrics content asserted directly) rather than skipping — the
//!    exact fallback discipline of `test_sort_cutover_parity.rs`. Cases cover
//!    multiple operations, an operation that matches zero records (pinning the
//!    zero-match warning), and a run with `--metrics` (pinning the TSV).
//!
//! **Not a RED/GREEN gate for the parity half.** Because the removed serial loop
//! and the chain were already output-equivalent, the baseline byte-parity check
//! passes on both sides of the change — like the sort cutover it guards
//! equivalence, it does not observe a regression the cutover introduces. The
//! chain-banner check in (1) is the part that flips RED→GREEN across the removal.

use std::ffi::OsStr;
use std::path::Path;
use std::process::Command;

use fgumi_raw_bam::{RawRecord, SamBuilder, flags};
use rstest::rstest;
use tempfile::TempDir;

use crate::helpers::bam_generator::{create_minimal_header, write_bam};
use crate::helpers::cutover::{baseline_bin, decompressed_records_without_pg};
use crate::helpers::read_bam_output;
use fgumi_lib::sam::SamTag;

/// One mapped record: name, optional `RX`, optional pre-existing `BX` (to
/// exercise `copy`/`move` overwrite), optional `MI`. Every third record is a
/// mapped paired read (real CIGAR, mate info, `MC`/`RG`) so the chain's decode
/// runs over non-trivial records, not just unmapped primaries — mirroring the
/// varied corpus in the in-process retag parity unit tests.
fn build_varied_records() -> Vec<RawRecord> {
    (0..60u32)
        .map(|i| {
            let mut b = SamBuilder::new();
            let name = format!("r{i:04}");
            b.read_name(name.as_bytes());
            // Record index 41 deliberately carries NO RX, so `RX::*` skips it
            // (src_missing) — and no record ever carries `ZZ`, pinning the
            // zero-match op path.
            if i != 41 {
                let rx = format!("ACGT{i:04}");
                b.add_string_tag(SamTag::RX, rx.as_bytes());
            }
            if i % 2 == 0 {
                b.add_string_tag(SamTag::MI, i.to_string().as_bytes());
            }
            if i % 5 == 0 {
                b.add_string_tag("BX".parse::<SamTag>().expect("BX tag"), b"OLD");
            }
            let idx = i32::try_from(i).expect("record index fits i32");
            if i % 3 == 0 {
                let mut f = flags::PAIRED;
                if i % 6 == 0 {
                    f |= flags::SECONDARY;
                }
                b.ref_id(0)
                    .pos(idx * 7 + 1)
                    .mapq(60)
                    .flags(f)
                    .mate_ref_id(0)
                    .mate_pos(idx * 7 + 50)
                    .cigar_ops(&[8 << 4]) // 8M
                    .sequence(b"ACGTACGT")
                    .qualities(&[30u8; 8])
                    .add_string_tag(SamTag::MC, b"8M")
                    .add_string_tag(SamTag::RG, b"A");
            } else {
                b.flags(0).ref_id(0).pos(idx + 1).sequence(b"ACGT").qualities(&[30u8; 4]);
            }
            b.build()
        })
        .collect()
}

/// Writes the varied corpus to `dir/in.bam` and returns its path.
fn write_varied_input(dir: &Path) -> std::path::PathBuf {
    let input = dir.join("in.bam");
    write_bam(&input, &create_minimal_header("chr1", 10_000), &build_varied_records());
    input
}

/// The operations every parity case runs: fan `RX` out to `BX` (overwriting the
/// pre-existing `BX` on every 5th record) and to `CB`, move `MI` to `XM`, then a
/// `ZZ::delete` that matches zero records. Covers copy/move/delete, overwrite,
/// `src_missing`, and the zero-match warn in one shot.
const PARITY_OPS: &[&str] = &["RX::copy::BX", "RX::copy::CB", "MI::move::XM", "ZZ::delete"];

/// Runs `<bin> retag -i <input> -o <output> [-M <metrics>] <ops...>` (no
/// `--threads`, so the current build takes the post-cutover chain path and the
/// baseline takes its serial path) with `RUST_LOG=info`, and returns the process
/// output for stderr assertions.
fn run_retag(
    bin: &Path,
    input: &Path,
    output: &Path,
    metrics: Option<&Path>,
    ops: &[&str],
) -> std::process::Output {
    let mut cmd = Command::new(bin);
    cmd.env("RUST_LOG", "info").args([
        OsStr::new("retag"),
        OsStr::new("-i"),
        input.as_os_str(),
        OsStr::new("-o"),
        output.as_os_str(),
    ]);
    if let Some(m) = metrics {
        cmd.args([OsStr::new("-M"), m.as_os_str()]);
    }
    cmd.args(ops.iter().map(OsStr::new));
    cmd.output().unwrap_or_else(|e| panic!("failed to spawn `{}` retag: {e}", bin.display()))
}

/// A no-`--threads` run now routes through the declarative chain, which logs the
/// `"Using pipeline with N threads"` banner from `ChainBuilder::add_retag`. The
/// retired serial tail logged no such line, so this is the RED (pre-removal) →
/// GREEN (post-removal) discriminator for the cutover.
#[test]
fn retag_no_threads_routes_through_chain() {
    let dir = TempDir::new().expect("temp dir");
    let input = write_varied_input(dir.path());
    let output = dir.path().join("out.bam");

    let current_bin = Path::new(env!("CARGO_BIN_EXE_fgumi"));
    let out = run_retag(current_bin, &input, &output, None, PARITY_OPS);
    let stderr = String::from_utf8_lossy(&out.stderr);
    assert!(out.status.success(), "no-threads retag must succeed; stderr:\n{stderr}");
    assert!(
        stderr.contains("Using pipeline with 1 threads"),
        "a no-`--threads` retag must route through the chain (which logs the pipeline banner); \
         the serial path is retired. stderr:\n{stderr}"
    );
    assert!(
        stderr.contains("Starting Retag"),
        "the chain must still emit the `Starting Retag` banner; stderr:\n{stderr}"
    );
}

/// Output parity of the post-cutover chain against the pre-removal serial
/// baseline binary, across representative operation sets — plus the
/// always-available self-consistency oracle when no baseline is set.
///
/// `#[case]` args: a label, the ops to run, and whether `--metrics` is written.
#[rstest]
#[case::multi_op_no_metrics(PARITY_OPS, false)]
#[case::multi_op_with_metrics(PARITY_OPS, true)]
#[case::zero_match_only(&["ZZ::delete"], true)]
fn cutover_matches_baseline(#[case] ops: &[&str], #[case] with_metrics: bool) {
    let dir = TempDir::new().expect("temp dir");
    let input = write_varied_input(dir.path());

    let current_bin = Path::new(env!("CARGO_BIN_EXE_fgumi"));
    let current_out = dir.path().join("current.bam");
    let current_tsv = dir.path().join("current.tsv");
    let current_metrics = with_metrics.then(|| current_tsv.clone());
    let current = run_retag(current_bin, &input, &current_out, current_metrics.as_deref(), ops);
    let current_stderr = String::from_utf8_lossy(&current.stderr);
    assert!(current.status.success(), "current retag must succeed; stderr:\n{current_stderr}");

    // A zero-match op must still warn on the chain (verbatim message from the
    // retired serial tail). ZZ is present on no record in every case's ops.
    if ops.contains(&"ZZ::delete") {
        assert!(
            current_stderr.contains(
                "operation 'ZZ::delete' matched zero records: no record carried the source tag 'ZZ'"
            ),
            "the zero-match warning must survive the cutover; stderr:\n{current_stderr}"
        );
    }

    if let Some(baseline) = baseline_bin() {
        let baseline_out = dir.path().join("baseline.bam");
        let baseline_tsv = dir.path().join("baseline.tsv");
        let baseline_metrics = with_metrics.then(|| baseline_tsv.clone());
        let base = run_retag(&baseline, &input, &baseline_out, baseline_metrics.as_deref(), ops);
        assert!(
            base.status.success(),
            "baseline retag failed; stderr:\n{}",
            String::from_utf8_lossy(&base.stderr)
        );

        assert_eq!(
            decompressed_records_without_pg(&current_out),
            decompressed_records_without_pg(&baseline_out),
            "chain retag output diverges from the pre-removal serial baseline binary ({}) \
             after stripping @PG — a real cutover parity bug, not something to relax",
            baseline.display(),
        );
        if with_metrics {
            assert_eq!(
                std::fs::read_to_string(&current_tsv).expect("current tsv"),
                std::fs::read_to_string(&baseline_tsv).expect("baseline tsv"),
                "chain --metrics TSV diverges from the serial baseline binary"
            );
        }
    } else {
        eprintln!(
            "SKIP baseline half of cutover_matches_baseline[{ops:?}]: FGUMI_BASELINE_BIN is \
             unset or does not name an existing file — running self-consistency oracle instead"
        );
        assert_self_consistent(&current_out, current_metrics.as_deref(), ops);
    }
}

/// Always-available oracle used when no baseline binary is set — it is the only
/// retag header oracle that runs in plain CI (the pre-removal serial path was the
/// old one). It asserts the chain output: (1) synthesizes `@HD VN:1.6 SO:unsorted`
/// and stamps a retag `@PG` provenance record on the header; (2) keeps every input
/// record; (3) applies exactly the tag rewrites `ops` prescribe (spot-checked on a
/// known record for the standard op set); and (4) writes one `--metrics` row per
/// op, in order, with the zero-match op reporting `records_applied == 0`. It does
/// NOT recompute every count from the input — that full byte-parity check is the
/// baseline-binary half of `cutover_matches_baseline`.
fn assert_self_consistent(output: &Path, metrics: Option<&Path>, ops: &[&str]) {
    assert_output_header(output);

    let input_records = build_varied_records();
    let (_, out_records) = read_bam_output(output);
    assert_eq!(
        out_records.len(),
        input_records.len(),
        "retag must keep every record (retag never adds or drops records)"
    );

    // For the standard multi-op set, spot-check the rewrites on a known record.
    if ops == PARITY_OPS {
        // Record 0 carries RX="ACGT0000", MI="0", and a pre-existing BX="OLD".
        let r0 = &out_records[0];
        assert_eq!(
            tag_string(r0, "BX").as_deref(),
            Some("ACGT0000"),
            "RX::copy::BX must overwrite the pre-existing BX"
        );
        assert_eq!(tag_string(r0, "CB").as_deref(), Some("ACGT0000"), "RX::copy::CB must add CB");
        assert_eq!(tag_string(r0, "XM").as_deref(), Some("0"), "MI::move::XM must add XM");
        assert!(tag_string(r0, "MI").is_none(), "MI::move::XM must drop MI");
        assert!(tag_string(r0, "ZZ").is_none(), "ZZ never present");
    }

    if let Some(path) = metrics {
        let tsv = std::fs::read_to_string(path).expect("read metrics");
        let lines: Vec<&str> = tsv.lines().collect();
        assert_eq!(
            lines[0], "operation\tkind\trecords_applied\tdst_overwritten\tsrc_missing",
            "metrics header"
        );
        assert_eq!(lines.len(), ops.len() + 1, "one metrics row per operation");
        // Every op named in the input appears as a row in order.
        for (op, line) in ops.iter().zip(&lines[1..]) {
            assert!(
                line.starts_with(&format!("{op}\t")),
                "metrics row order/label mismatch: {line} (expected {op})"
            );
        }
        // The zero-match op reports records_applied == 0.
        if let Some(zz) = lines.iter().find(|l| l.starts_with("ZZ::delete\t")) {
            let fields: Vec<&str> = zz.split('\t').collect();
            assert_eq!(fields.get(2).copied(), Some("0"), "ZZ::delete records_applied must be 0");
        }
    }
}

/// Asserts the chain output header carries the synthesized `@HD VN:1.6
/// SO:unsorted` and a retag `@PG` provenance record.
///
/// Reads the header un-normalized (a fresh `noodles::bam` reader, NOT
/// `read_bam_output`, which rewrites the `@PG` `CL` to `<normalized>`) so the
/// `@PG` command line — which names `retag` — is still intact to assert on. When
/// the serial path was removed this became the only retag test that pins @HD
/// synthesis / @PG stamping in plain CI.
fn assert_output_header(output: &Path) {
    use noodles::sam::header::record::value::map::header::{Version, tag as header_tag};
    use noodles::sam::header::record::value::map::program::tag as program_tag;

    let mut reader = noodles::bam::io::reader::Builder
        .build_from_path(output)
        .unwrap_or_else(|e| panic!("open {}: {e}", output.display()));
    let header = reader.read_header().expect("read output header");

    let hd = header.header().expect("output BAM must carry a synthesized @HD line");
    assert_eq!(hd.version(), Version::new(1, 6), "synthesized @HD VN must be 1.6");
    assert_eq!(
        hd.other_fields().get(&header_tag::SORT_ORDER).map(|v| v.to_vec()),
        Some(b"unsorted".to_vec()),
        "synthesized @HD SO must be 'unsorted'"
    );

    // The input carries no @PG, so any @PG present is the one retag stamped. Its
    // CL names `retag`, which pins that this provenance record came from this run.
    let stamped_retag_pg = header.programs().as_ref().values().any(|program| {
        program
            .other_fields()
            .get(&program_tag::COMMAND_LINE)
            .is_some_and(|cl| cl.windows(5).any(|w| w == b"retag"))
    });
    assert!(stamped_retag_pg, "output header must carry a retag @PG provenance record");
}

/// Reads a string aux tag (named by its two-char SAM tag) off a parsed
/// `RecordBuf`, or `None` when the tag is absent or not a string value.
fn tag_string(record: &noodles::sam::alignment::RecordBuf, tag: &str) -> Option<String> {
    use noodles::sam::alignment::record::data::field::Tag;
    use noodles::sam::alignment::record_buf::data::field::Value;
    let sam_tag: SamTag = tag.parse().expect("valid two-char tag");
    match record.data().get(&Tag::from(sam_tag))? {
        Value::String(s) => Some(s.to_string()),
        _ => None,
    }
}
