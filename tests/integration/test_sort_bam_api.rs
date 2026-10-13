//! `fgumi_lib::sort::sort_bam`, the library entry to the sort, against the
//! `fgumi sort` CLI; and the staged output both publish only after a
//! successful run.

use std::path::{Path, PathBuf};

use fgumi_lib::commands::sort::SortOrderArg;
use fgumi_lib::sam::SamTag;
use fgumi_lib::sort::{KeyTypesSpec, SortConfig};
use fgumi_raw_bam::{RawRecord, SamBuilder, flags};
use rstest::rstest;
use tempfile::TempDir;

use crate::helpers::bam_generator::{create_minimal_header, write_bam};
use crate::helpers::cutover::decompressed_records_without_pg;

/// `n` unsorted paired records whose `CB` and `MI` vary, so a template key
/// that dropped either lane orders them differently (or fails the sort).
fn lane_varying_records(n: usize) -> Vec<RawRecord> {
    (0..n)
        .map(|i| {
            let mut b = SamBuilder::new();
            let pos = i32::try_from(1 + (i * 7919 % 400) * 10).expect("pos fits i32");
            let mate = if i % 2 == 0 { flags::FIRST_SEGMENT } else { flags::LAST_SEGMENT };
            b.read_name(format!("t{}", i / 2).as_bytes())
                .sequence(b"ACGTACGTAC")
                .qualities(&[30u8; 10])
                .ref_id(0)
                .pos(pos)
                .mapq(60)
                .flags(flags::PAIRED | mate)
                .cigar_ops(&[10u32 << 4]);
            b.add_string_tag(SamTag::CB, ["AAAC", "CCAG", "GGTA"][i % 3].as_bytes());
            b.add_string_tag(SamTag::MI, format!("{}", i % 7).as_bytes());
            b.build()
        })
        .collect()
}

fn config(input: &Path, output: &Path, order: SortOrderArg, threads: usize) -> SortConfig {
    SortConfig {
        input: input.to_path_buf(),
        output: output.to_path_buf(),
        order,
        threads,
        max_memory: Some(64 * 1024),
        tmp_dirs: vec![input.parent().unwrap().to_path_buf()],
        compression: 1,
        write_index: false,
        key_types: KeyTypesSpec::default(),
    }
}

/// `fgumi sort` with `args` from `input` into `output`, returning whether it
/// succeeded.
fn cli_sort(input: &Path, output: &Path, args: &[&str]) -> bool {
    std::process::Command::new(env!("CARGO_BIN_EXE_fgumi"))
        .args(["sort", "-i"])
        .arg(input)
        .arg("-o")
        .arg(output)
        .args(["--max-memory", "64K", "--memory-per-thread", "false"])
        .args(args)
        .status()
        .unwrap()
        .success()
}

/// The library entry builds the CLI's chain, key types included: on records
/// whose `CB` and `MI` vary, every `--key-types` either writes the same records
/// through both or fails through both, as each case expects. A key that drops
/// a lane that varies (`none`, `cb`, and `auto`, which takes its lanes from
/// the first record) fails naming that lane. Catches: `sort_bam` ignoring
/// `key_types`, the order or the thread count, and a regression that fails
/// both paths the same way.
#[rstest]
fn sort_bam_matches_the_cli(
    #[values(1usize, 4)] threads: usize,
    #[values(
        ("coordinate", None, Ok(())),
        ("queryname", None, Ok(())),
        ("template-coordinate", None, Err("MI")),
        ("template-coordinate", Some(("full", KeyTypesSpec::Full)), Ok(())),
        ("template-coordinate", Some(("none", KeyTypesSpec::None)), Err("CB")),
        (
            "template-coordinate",
            Some(("cb", KeyTypesSpec::Explicit { cb: true, tertiary: false })),
            Err("MI")
        )
    )]
    case: (&str, Option<(&str, KeyTypesSpec)>, Result<(), &str>),
) {
    let (order, key_types, expected) = case;
    let dir = TempDir::new().unwrap();
    let input = dir.path().join("in.bam");
    write_bam(&input, &create_minimal_header("chr1", 1_000_000), &lane_varying_records(4_000));
    let via_cli = dir.path().join("cli.bam");
    let threads_arg = threads.to_string();
    let mut args = vec!["--order", order, "--threads", threads_arg.as_str()];
    if let Some((token, _)) = key_types {
        args.extend(["--key-types", token]);
    }
    let cli_ok = cli_sort(&input, &via_cli, &args);
    assert_eq!(cli_ok, expected.is_ok(), "{order} {key_types:?} t{threads}: the CLI");
    let via_api = dir.path().join("api.bam");
    let mut cfg = config(&input, &via_api, SortOrderArg::parse(order).unwrap(), threads);
    if let Some((_, spec)) = key_types {
        cfg.key_types = spec;
    }
    let api = fgumi_lib::sort::sort_bam(&cfg);
    assert_eq!(api.is_ok(), cli_ok, "{order} {key_types:?} t{threads}: {api:?}");
    if let (Err(lane), Err(e)) = (expected, &api) {
        let msg = format!("{e:#}");
        assert!(msg.contains(&format!("carries a {lane} value absent")), "{msg}");
    }
    if let Ok(stats) = api {
        assert_eq!(stats.total_records, 4_000);
        assert!(
            decompressed_records_without_pg(&via_api) == decompressed_records_without_pg(&via_cli),
            "{order} {key_types:?} t{threads}"
        );
    }
}

/// `write_index: true` writes a loadable BAI beside the output and the same
/// records as without it. Catches: `sort_bam` ignoring `write_index`.
#[test]
fn sort_bam_writes_a_loadable_index() {
    let dir = TempDir::new().unwrap();
    let input = dir.path().join("in.bam");
    write_bam(&input, &create_minimal_header("chr1", 1_000_000), &lane_varying_records(4_000));
    let (plain, indexed) = (dir.path().join("plain.bam"), dir.path().join("indexed.bam"));
    fgumi_lib::sort::sort_bam(&config(&input, &plain, SortOrderArg::Coordinate, 2)).unwrap();
    let mut cfg = config(&input, &indexed, SortOrderArg::Coordinate, 2);
    cfg.write_index = true;
    fgumi_lib::sort::sort_bam(&cfg).unwrap();
    let index = noodles::bam::bai::fs::read(fgumi_bam_io::bai_sidecar_path(&indexed)).unwrap();
    assert_eq!(index.reference_sequences().len(), 1);
    assert!(decompressed_records_without_pg(&plain) == decompressed_records_without_pg(&indexed));
}

/// `compression` reaches the written BAM. Catches: `sort_bam` ignoring it.
#[test]
fn sort_bam_honors_the_compression_level() {
    let dir = TempDir::new().unwrap();
    let input = dir.path().join("in.bam");
    write_bam(&input, &create_minimal_header("chr1", 1_000_000), &lane_varying_records(4_000));
    let size = |level: u32, name: &str| {
        let out = dir.path().join(name);
        let mut cfg = config(&input, &out, SortOrderArg::Coordinate, 2);
        cfg.compression = level;
        fgumi_lib::sort::sort_bam(&cfg).unwrap();
        std::fs::metadata(&out).unwrap().len()
    };
    let (stored, packed) = (size(0, "l0.bam"), size(9, "l9.bam"));
    assert!(packed < stored / 2, "level 9 ({packed}) did not shrink level 0 ({stored})");
}

/// Catches: `write_index` accepted for an order that cannot be indexed.
#[test]
fn sort_bam_rejects_an_index_for_a_non_coordinate_order() {
    let dir = TempDir::new().unwrap();
    let mut cfg =
        config(&dir.path().join("in.bam"), &dir.path().join("out.bam"), SortOrderArg::Queryname, 1);
    cfg.write_index = true;
    let err = fgumi_lib::sort::sort_bam(&cfg).unwrap_err();
    assert!(err.to_string().contains("only valid for coordinate"), "{err:#}");
}

/// A valid BAM cut off inside its records: the header parses (the chain
/// builds), and the run fails mid-way.
fn truncated_input(dir: &Path) -> PathBuf {
    let input = dir.join("truncated.bam");
    write_bam(&input, &create_minimal_header("chr1", 1_000_000), &lane_varying_records(20_000));
    let bytes = std::fs::read(&input).unwrap();
    std::fs::write(&input, &bytes[..bytes.len() * 3 / 5]).unwrap();
    input
}

/// A sort that fails mid-run publishes nothing, through the CLI and the API, with
/// and without an index: a new destination is not created and an existing one
/// (and its index) keeps its bytes. Catches: the chain writing its output in
/// place (a truncated BAM left at the destination, as the CLI did since the
/// cutover).
#[rstest]
fn a_failed_sort_leaves_the_destination_as_it_was(
    #[values(true, false)] via_cli: bool,
    #[values(true, false)] existing: bool,
    #[values(true, false)] write_index: bool,
) {
    let dir = TempDir::new().unwrap();
    let input = truncated_input(dir.path());
    let out = dir.path().join("out.bam");
    let bai = fgumi_bam_io::bai_sidecar_path(&out);
    if existing {
        std::fs::write(&out, b"previous bam").unwrap();
        std::fs::write(&bai, b"previous bai").unwrap();
    }
    let ok = if via_cli {
        let mut args = vec!["--order", "coordinate", "--threads", "2"];
        if write_index {
            args.push("--write-index");
        }
        cli_sort(&input, &out, &args)
    } else {
        let mut cfg = config(&input, &out, SortOrderArg::Coordinate, 2);
        cfg.write_index = write_index;
        fgumi_lib::sort::sort_bam(&cfg).is_ok()
    };
    assert!(!ok, "a truncated input must fail the sort");
    if existing {
        assert_eq!(std::fs::read(&out).unwrap(), b"previous bam", "the destination changed");
        assert_eq!(std::fs::read(&bai).unwrap(), b"previous bai", "the index changed");
    } else {
        assert!(!out.exists() && !bai.exists(), "a failed sort created its output");
    }
    let leftovers: Vec<_> = std::fs::read_dir(dir.path())
        .unwrap()
        .filter_map(Result::ok)
        .filter(|e| e.file_name().to_string_lossy().starts_with(".fgumi-"))
        .collect();
    assert!(leftovers.is_empty(), "staged files left behind: {leftovers:?}");
}
