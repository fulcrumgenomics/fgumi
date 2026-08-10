//! Determinism regression tests for consensus downsampling.
//!
//! Downsampling used to select which reads to keep by shuffling with a seeded RNG held on
//! the consensus caller. The RNG is stateful, so the surviving subset for a family depended
//! on how many families that caller instance had already processed.
//!
//! fgumi never had fgbio's thread-count bug — the unified pipeline builds a fresh caller per
//! batch at a fixed batch size, so `--threads 1` and `--threads 8` agreed. But the
//! no-`--threads` path holds a single caller for the entire run, putting its RNG at a
//! different position for every family, so the two execution *modes* disagreed under a cap.
//! On a ~77k-molecule library that diverged on 566 of 77,804 duplex records (0.73%).
//!
//! Ranking reads by a Murmur3 hash of the read name makes the retained subset a pure
//! function of the family, so every mode agrees. These tests pin that.
//!
//! Two properties of the fixtures are load-bearing, and a test that loses either would pass
//! vacuously — including against the old shuffle:
//!
//! 1. **Reads within a family differ in sequence.** If every read were identical, any
//!    selection would yield the same consensus and no divergence could be observed.
//! 2. **The family count exceeds the pipeline's per-batch group count** (50 for simplex, 100
//!    for duplex). The caller is constructed once per batch, so a single-batch input would
//!    put every family at the same RNG position and hide the divergence entirely.
//!
//! Uncapped runs always agreed across modes; the uncapped cases guard the default
//! configuration against a regression.

use clap::Parser;
use fgumi_lib::commands::command::Command;
use fgumi_lib::commands::duplex::Duplex;
use fgumi_lib::commands::simplex::Simplex;
use fgumi_lib::sam::SamTag;
use fgumi_raw_bam::{RawRecord, SamBuilder, flags};
use noodles::bam;
use noodles::sam::Header;
use noodles::sam::alignment::io::Write as AlignmentWrite;
use noodles::sam::alignment::record::data::field::Value;
use noodles::sam::alignment::record::data::field::value::Array;
use noodles::sam::alignment::record::data::field::value::array::Values;
use rstest::rstest;
use std::fs;
use std::path::{Path, PathBuf};
use tempfile::TempDir;

use crate::helpers::bam_generator::{create_minimal_header, to_record_buf};

/// The fields that must match across execution modes for two runs to be identical.
///
/// Compared field by field rather than by `Debug`-formatting the whole record, so a failure
/// names the field that diverged and the assertion cannot be weakened by a change to any
/// type's `Debug` implementation. Sequence and the consensus tags are what a change in
/// *which* reads were retained actually moves.
#[derive(Debug, PartialEq, Eq)]
struct RecordFingerprint {
    name: Vec<u8>,
    flags: u16,
    sequence: Vec<u8>,
    quality_scores: Vec<u8>,
    tags: Vec<(String, String)>,
}

/// Renders an aux value canonically, by matching the value's own type rather than formatting
/// it with `Debug`.
///
/// Keeps the fingerprint's promise that no type's `Debug` implementation can weaken the
/// comparison: every variant is rendered from its actual contents, floats through their bit
/// pattern so the comparison stays exact.
fn render_tag_value(value: &Value<'_>) -> String {
    match value {
        Value::Character(c) => format!("A:{}", char::from(*c)),
        Value::Int8(_)
        | Value::UInt8(_)
        | Value::Int16(_)
        | Value::UInt16(_)
        | Value::Int32(_)
        | Value::UInt32(_) => {
            format!("i:{}", value.as_int().expect("integer variant has an integer value"))
        }
        Value::Float(f) => format!("f:{:#010x}", f.to_bits()),
        Value::String(s) => format!("Z:{s}"),
        Value::Hex(s) => format!("H:{s}"),
        Value::Array(array) => format!("B:{}", render_array(array)),
    }
}

/// Renders an aux array as `<subtype>:<comma-separated values>`.
fn render_array(array: &Array<'_>) -> String {
    fn join<N: std::fmt::Display>(values: &(dyn Values<'_, N> + '_)) -> String {
        values
            .iter()
            .map(|v| v.expect("read array value").to_string())
            .collect::<Vec<_>>()
            .join(",")
    }
    match array {
        Array::Int8(v) => format!("c:{}", join(v.as_ref())),
        Array::UInt8(v) => format!("C:{}", join(v.as_ref())),
        Array::Int16(v) => format!("s:{}", join(v.as_ref())),
        Array::UInt16(v) => format!("S:{}", join(v.as_ref())),
        Array::Int32(v) => format!("i:{}", join(v.as_ref())),
        Array::UInt32(v) => format!("I:{}", join(v.as_ref())),
        Array::Float(v) => {
            let bits: Vec<String> =
                v.iter().map(|f| format!("{:#010x}", f.expect("read value").to_bits())).collect();
            format!("f:{}", bits.join(","))
        }
    }
}

/// The largest value of an integer aux tag across every record of a BAM.
///
/// `cD` is the number of raw reads that contributed to a consensus, so this is the absolute
/// oracle for "the cap actually capped": no consensus may be built from more reads than the cap
/// allows, whatever execution mode produced it.
fn max_int_tag(path: &Path, tag: SamTag) -> Option<i64> {
    let mut reader = bam::io::Reader::new(fs::File::open(path).expect("open output BAM"));
    let _header = reader.read_header().expect("read header");
    reader
        .records()
        .filter_map(|result| {
            let record = result.expect("read record");
            record.data().get(tag.as_ref()).map(|value| {
                value.expect("read aux field").as_int().expect("depth tag is an integer")
            })
        })
        .max()
}

/// Reads every record of a BAM into a comparable fingerprint, in file order.
fn fingerprints(path: &Path) -> Vec<RecordFingerprint> {
    let mut reader = bam::io::Reader::new(fs::File::open(path).expect("open output BAM"));
    let _header = reader.read_header().expect("read header");
    reader
        .records()
        .map(|result| {
            let record = result.expect("read record");
            let mut tags: Vec<(String, String)> = record
                .data()
                .iter()
                .map(|field| {
                    let (tag, value) = field.expect("read aux field");
                    (String::from_utf8_lossy(tag.as_ref()).into_owned(), render_tag_value(&value))
                })
                .collect();
            // Aux field ordering is not part of the contract; the set of values is.
            tags.sort();
            RecordFingerprint {
                name: record.name().map(|n| n.to_vec()).unwrap_or_default(),
                flags: record.flags().bits(),
                sequence: record.sequence().iter().collect::<Vec<u8>>(),
                quality_scores: record.quality_scores().as_ref().to_vec(),
                tags,
            }
        })
        .collect()
}

/// An 8bp sequence carrying a single `C` at `offset % 8`, the rest `A`.
///
/// Gives each read in a family distinct content so *which* reads survive a cap is visible in
/// the consensus, not merely in its depth.
fn distinguishing_sequence(offset: usize) -> String {
    let mut bases = vec![b'A'; 8];
    bases[offset % 8] = b'C';
    String::from_utf8(bases).expect("ASCII bases")
}

/// Writes `records` as a BAM with a minimal single-reference header.
fn write_bam(path: &Path, header: &Header, records: &[RawRecord]) {
    let mut writer = bam::io::Writer::new(fs::File::create(path).expect("create BAM"));
    writer.write_header(header).expect("write header");
    for raw in records {
        writer.write_alignment_record(header, &to_record_buf(raw)).expect("write record");
    }
    writer.try_finish().expect("finish BAM");
}

/// One `MI`-tagged FR read pair for the simplex fixture.
fn simplex_pair(name: &str, mi: &str, sequence: &str) -> (RawRecord, RawRecord) {
    let seq = sequence.as_bytes();
    let read_len = seq.len();
    let cigar_op = u32::try_from(read_len).expect("read_len fits u32") << 4;
    let span = i32::try_from(read_len + 100).expect("template span fits i32");

    let build = |first: bool| {
        let (pos, mate_pos, rev, mate_rev) =
            if first { (100, 200, false, true) } else { (200, 100, true, false) };
        let segment = if first { flags::FIRST_SEGMENT } else { flags::LAST_SEGMENT };
        let mut b = SamBuilder::new();
        b.read_name(name.as_bytes())
            .sequence(seq)
            .qualities(&vec![30u8; read_len])
            .flags(
                flags::PAIRED
                    | segment
                    | if rev { flags::REVERSE } else { 0 }
                    | if mate_rev { flags::MATE_REVERSE } else { 0 },
            )
            .ref_id(0)
            .pos(pos - 1)
            .mapq(60)
            .cigar_ops(&[cigar_op])
            .mate_ref_id(0)
            .mate_pos(mate_pos - 1)
            .template_length(if first { span } else { -span })
            .add_string_tag(SamTag::MI, mi.as_bytes());
        b.build()
    };

    (build(true), build(false))
}

/// One `MI`-tagged read pair on the `/A` (FR) or `/B` (RF) strand for the duplex fixture.
///
/// The orientation is the point: the duplex caller pairs an FR template with the RF template
/// at the same coordinates, so wrong flags or mate positions silently yield two single-strand
/// molecules instead of one duplex.
fn duplex_pair(name: &str, mi: &str, sequence: &str, is_b_strand: bool) -> (RawRecord, RawRecord) {
    let seq = sequence.as_bytes();
    let read_len = seq.len();
    let cigar_op = u32::try_from(read_len).expect("read_len fits u32") << 4;
    let span = i32::try_from(read_len + 100).expect("template span fits i32");

    let (r1_start, r2_start, r1_rev, r2_rev) =
        if is_b_strand { (200, 100, true, false) } else { (100, 200, false, true) };
    let tlen = if is_b_strand { -span } else { span };

    let build = |first: bool| {
        let (pos, mate_pos, rev, mate_rev) = if first {
            (r1_start, r2_start, r1_rev, r2_rev)
        } else {
            (r2_start, r1_start, r2_rev, r1_rev)
        };
        let segment = if first { flags::FIRST_SEGMENT } else { flags::LAST_SEGMENT };
        let mut b = SamBuilder::new();
        b.read_name(name.as_bytes())
            .sequence(seq)
            .qualities(&vec![30u8; read_len])
            .flags(
                flags::PAIRED
                    | segment
                    | if rev { flags::REVERSE } else { 0 }
                    | if mate_rev { flags::MATE_REVERSE } else { 0 },
            )
            .ref_id(0)
            .pos(pos - 1)
            .mapq(60)
            .cigar_ops(&[cigar_op])
            .mate_ref_id(0)
            .mate_pos(mate_pos - 1)
            .template_length(if first { tlen } else { -tlen })
            .add_string_tag(SamTag::RX, b"AAA-TTT")
            .add_string_tag(SamTag::MI, mi.as_bytes());
        b.build()
    };

    (build(true), build(false))
}

/// Builds a grouped BAM of `families` MI families, each with `depth` read pairs.
fn write_simplex_input(path: &Path, families: usize, depth: usize) {
    let header = create_minimal_header("chr1", 10_000);
    let mut records: Vec<RawRecord> = Vec::new();
    for family in 0..families {
        for read in 0..depth {
            let (r1, r2) = simplex_pair(
                &format!("f{family}:{read}"),
                &format!("{family}"),
                &distinguishing_sequence(read),
            );
            records.push(r1);
            records.push(r2);
        }
    }
    write_bam(path, &header, &records);
}

/// Builds a duplex-grouped BAM: `families` molecules, each with `depth_per_strand` read
/// pairs on `/A` and on `/B`.
fn write_duplex_input(path: &Path, families: usize, depth_per_strand: usize) {
    let header = create_minimal_header("chr1", 10_000);
    let mut records: Vec<RawRecord> = Vec::new();
    for family in 0..families {
        for (suffix, is_b_strand) in [("A", false), ("B", true)] {
            for read in 0..depth_per_strand {
                let (r1, r2) = duplex_pair(
                    &format!("m{family}{suffix}:{read}"),
                    &format!("{family}/{suffix}"),
                    &distinguishing_sequence(read),
                    is_b_strand,
                );
                records.push(r1);
                records.push(r2);
            }
        }
    }
    write_bam(path, &header, &records);
}

/// Where a labelled run of `command` writes its output.
fn output_path(dir: &Path, command: &str, label: &str) -> PathBuf {
    dir.join(format!("{command}.{label}.bam"))
}

/// Asserts that a cap actually constrained the output, in two independent ways.
///
/// Without this, the cross-mode assertions alone are satisfied by an implementation that ignores
/// the cap in *every* mode — dropping the option from both the sequential and the threaded
/// construction sites leaves all three runs equal and uncapped, and equally green.
///
/// The depth tags are the absolute check (no consensus may be built from more raw reads than
/// the cap allows); the comparison against the uncapped baseline is the relative one, and catches
/// a cap that is honoured but has no effect on this fixture.
///
/// `depth_tags` must name the depths the cap actually bounds. For `simplex --max-reads` that is
/// the combined `cD`; for `duplex --max-reads-per-strand` it is the per-strand `aD` and `bD` —
/// `cD` sums both strands and legitimately reaches twice the cap.
fn assert_cap_took_effect(
    capped: &[RecordFingerprint],
    uncapped: &[RecordFingerprint],
    capped_bam: &Path,
    cap: &str,
    depth_tags: &[SamTag],
) {
    let cap: i64 = cap.parse().expect("cap is numeric");
    for &tag in depth_tags {
        let name = String::from_utf8_lossy(tag.as_ref()).into_owned();
        let max_depth =
            max_int_tag(capped_bam, tag).unwrap_or_else(|| panic!("no record carries {name}"));
        assert!(
            max_depth <= cap,
            "a consensus used {max_depth} raw reads ({name}), exceeding the cap of {cap}"
        );
    }
    assert_ne!(
        capped, uncapped,
        "capped output is identical to the uncapped run, so the cap changed nothing"
    );
}

/// Runs `simplex` with the given cap and thread setting, returning the output fingerprints.
fn run_simplex(
    dir: &Path,
    input: &Path,
    label: &str,
    cap: Option<&str>,
    threads: Option<&str>,
) -> Vec<RecordFingerprint> {
    let output = output_path(dir, "simplex", label);
    let mut args: Vec<String> = vec![
        "simplex".into(),
        "--input".into(),
        input.to_str().expect("utf-8 path").into(),
        "--output".into(),
        output.to_str().expect("utf-8 path").into(),
        "--min-reads".into(),
        "1".into(),
        "--compression-level".into(),
        "1".into(),
    ];
    if let Some(cap) = cap {
        args.push("--max-reads".into());
        args.push(cap.into());
    }
    if let Some(threads) = threads {
        args.push("--threads".into());
        args.push(threads.into());
    }
    let cmd = Simplex::try_parse_from(&args).expect("parse simplex args");
    cmd.execute("fgumi simplex").expect("simplex should succeed");
    fingerprints(&output)
}

/// Runs `duplex` with the given cap and thread setting, returning the output fingerprints.
fn run_duplex(
    dir: &Path,
    input: &Path,
    label: &str,
    cap: Option<&str>,
    threads: Option<&str>,
) -> Vec<RecordFingerprint> {
    let output = output_path(dir, "duplex", label);
    let mut args: Vec<String> = vec![
        "duplex".into(),
        "--input".into(),
        input.to_str().expect("utf-8 path").into(),
        "--output".into(),
        output.to_str().expect("utf-8 path").into(),
        "--min-reads".into(),
        "1".into(),
        "--compression-level".into(),
        "1".into(),
    ];
    if let Some(cap) = cap {
        args.push("--max-reads-per-strand".into());
        args.push(cap.into());
    }
    if let Some(threads) = threads {
        args.push("--threads".into());
        args.push(threads.into());
    }
    let cmd = Duplex::try_parse_from(&args).expect("parse duplex args");
    cmd.execute("fgumi duplex").expect("duplex should succeed");
    fingerprints(&output)
}

/// `simplex` output must not depend on which execution mode produced it.
///
/// The capped case is the regression: before hash-based selection, `--threads N` and the
/// no-`--threads` path retained different reads. The uncapped case guards the default
/// configuration, which was always mode-independent and must stay so.
/// The odd cap is not redundant with the even one: reads within a family arrive interleaved
/// R1/R2 and mates share a rank, so a cap of 2 lands exactly on a tie boundary and never
/// consults the stable sort, while a cap of 3 lands *inside* one — making arrival order, the
/// last selection input that is not a pure function of the family, observable across modes.
#[rstest]
#[case::capped(Some("2"))]
#[case::capped_odd(Some("3"))]
#[case::uncapped(None)]
fn simplex_output_is_identical_across_execution_modes(#[case] cap: Option<&str>) {
    let temp = TempDir::new().expect("temp dir");
    let input = temp.path().join("grouped.bam");
    write_simplex_input(&input, 60, 6);

    let single = run_simplex(temp.path(), &input, "t1", cap, Some("1"));
    let many = run_simplex(temp.path(), &input, "t8", cap, Some("8"));
    let no_threads = run_simplex(temp.path(), &input, "none", cap, None);

    assert!(!single.is_empty(), "expected consensus reads to be produced");
    assert_eq!(single, many, "--threads 1 and --threads 8 must produce identical output");
    assert_eq!(single, no_threads, "--threads and no-threads must produce identical output");

    if let Some(cap) = cap {
        let uncapped = run_simplex(temp.path(), &input, "uncapped", None, Some("1"));
        assert_cap_took_effect(
            &single,
            &uncapped,
            &output_path(temp.path(), "simplex", "t1"),
            cap,
            &[SamTag::CD],
        );
    }
}

/// `duplex` output must not depend on which execution mode produced it.
///
/// This is the path where the divergence was measured on real data: 566 of 77,804 records
/// differed between `--threads N` and the no-threads path at `--max-reads-per-strand 2`.
#[rstest]
#[case::capped(Some("2"))]
#[case::capped_odd(Some("3"))]
#[case::uncapped(None)]
fn duplex_output_is_identical_across_execution_modes(#[case] cap: Option<&str>) {
    let temp = TempDir::new().expect("temp dir");
    let input = temp.path().join("duplex-grouped.bam");
    write_duplex_input(&input, 120, 4);

    let single = run_duplex(temp.path(), &input, "t1", cap, Some("1"));
    let many = run_duplex(temp.path(), &input, "t8", cap, Some("8"));
    let no_threads = run_duplex(temp.path(), &input, "none", cap, None);

    assert!(!single.is_empty(), "expected duplex consensus reads to be produced");
    assert_eq!(single, many, "--threads 1 and --threads 8 must produce identical output");
    assert_eq!(single, no_threads, "--threads and no-threads must produce identical output");

    if let Some(cap) = cap {
        let uncapped = run_duplex(temp.path(), &input, "uncapped", None, Some("1"));
        assert_cap_took_effect(
            &single,
            &uncapped,
            &output_path(temp.path(), "duplex", "t1"),
            cap,
            &[SamTag::AD, SamTag::BD],
        );
    }
}
