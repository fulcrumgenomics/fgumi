//! Integration tests for streaming input support (stdin/pipes).
//!
//! These tests spawn cat processes whose stdout is piped to fgumi commands.
//! The child processes are properly cleaned up when their stdout is consumed.
#![allow(clippy::zombie_processes)]

use noodles::bam;
use noodles::sam::alignment::io::Write as AlignmentWrite;
use rstest::rstest;
use std::fs;
use std::path::PathBuf;
use std::process::{Command, Stdio};

use tempfile::TempDir;

use crate::helpers::bam_generator::{create_minimal_header, create_umi_family, to_record_buf};

/// Test that the group command works correctly with piped input.
#[test]
fn test_group_command_with_piped_input() {
    let temp_dir = TempDir::new().expect("Failed to create temp dir");
    let input_bam = temp_dir.path().join("input.bam");
    let output_from_file = temp_dir.path().join("output_file.bam");
    let output_from_pipe = temp_dir.path().join("output_pipe.bam");

    // Create test BAM file
    create_test_input_bam(&input_bam);

    // Run group command with file input (baseline)
    let status = Command::new(env!("CARGO_BIN_EXE_fgumi"))
        .args([
            "group",
            "--input",
            input_bam.to_str().unwrap(),
            "--output",
            output_from_file.to_str().unwrap(),
            "--strategy",
            "identity",
            "--edits",
            "0",
            "--compression-level",
            "1",
            "--threads",
            "2",
        ])
        .status()
        .expect("Failed to run group command with file input");
    assert!(status.success(), "Group command with file input failed");
    assert!(output_from_file.exists(), "Output BAM from file not created");

    // Run group command with piped input (cat file | fgumi group)
    let cat_child = Command::new("cat")
        .arg(input_bam.to_str().unwrap())
        .stdout(Stdio::piped())
        .spawn()
        .expect("Failed to spawn cat");

    let status = Command::new(env!("CARGO_BIN_EXE_fgumi"))
        .args([
            "group",
            "--input",
            "-", // stdin
            "--output",
            output_from_pipe.to_str().unwrap(),
            "--strategy",
            "identity",
            "--edits",
            "0",
            "--compression-level",
            "1",
            "--threads",
            "2",
        ])
        .stdin(cat_child.stdout.unwrap())
        .status()
        .expect("Failed to run group command with piped input");
    assert!(status.success(), "Group command with piped input failed");
    assert!(output_from_pipe.exists(), "Output BAM from pipe not created");

    // Compare outputs - records should be identical (headers may differ due to @PG command line)
    compare_bam_records(&output_from_file, &output_from_pipe);
}

/// Test that /dev/stdin path also works.
#[test]
fn test_group_command_with_dev_stdin_path() {
    let temp_dir = TempDir::new().expect("Failed to create temp dir");
    let input_bam = temp_dir.path().join("input.bam");
    let output_bam = temp_dir.path().join("output.bam");

    // Create test BAM file
    create_test_input_bam(&input_bam);

    // Run group command using /dev/stdin
    let cat_child = Command::new("cat")
        .arg(input_bam.to_str().unwrap())
        .stdout(Stdio::piped())
        .spawn()
        .expect("Failed to spawn cat");

    let status = Command::new(env!("CARGO_BIN_EXE_fgumi"))
        .args([
            "group",
            "--input",
            "/dev/stdin",
            "--output",
            output_bam.to_str().unwrap(),
            "--strategy",
            "identity",
            "--edits",
            "0",
            "--compression-level",
            "1",
            "--threads",
            "2",
        ])
        .stdin(cat_child.stdout.unwrap())
        .status()
        .expect("Failed to run group command with /dev/stdin");

    assert!(status.success(), "Group command with /dev/stdin failed");
    assert!(output_bam.exists(), "Output BAM not created");
}

/// S5c2-002: the downsample command must accept stdin input (`-`) instead of
/// failing the up-front existence check with a spurious "file not found".
/// Downsample requires template-coordinate-sorted input, so we first sort the
/// fixture, then pipe the sorted BAM into `fgumi downsample -i -`.
#[test]
fn test_downsample_command_with_piped_input() {
    let temp_dir = TempDir::new().expect("Failed to create temp dir");
    let input_bam = temp_dir.path().join("input.bam");
    let sorted_bam = temp_dir.path().join("sorted.bam");
    let output_bam = temp_dir.path().join("output.bam");

    // Downsample requires grouped (MI-tagged) reads sorted by template
    // coordinate.
    create_grouped_test_bam(&input_bam);

    // Produce a template-coordinate-sorted BAM (downsample's required input).
    let status = Command::new(env!("CARGO_BIN_EXE_fgumi"))
        .args([
            "sort",
            "--input",
            input_bam.to_str().unwrap(),
            "--output",
            sorted_bam.to_str().unwrap(),
            "--order",
            "template-coordinate",
            "--threads",
            "1",
        ])
        .status()
        .expect("Failed to run sort");
    assert!(status.success(), "sort failed");

    // Pipe the sorted BAM into downsample via stdin.
    let cat_child = Command::new("cat")
        .arg(sorted_bam.to_str().unwrap())
        .stdout(Stdio::piped())
        .spawn()
        .expect("Failed to spawn cat");

    let output = Command::new(env!("CARGO_BIN_EXE_fgumi"))
        .args([
            "downsample",
            "--input",
            "-", // stdin
            "--output",
            output_bam.to_str().unwrap(),
            "--fraction",
            "1.0",
            "--seed",
            "42",
            "--compression-level",
            "1",
        ])
        .stdin(cat_child.stdout.unwrap())
        .output()
        .expect("Failed to run downsample with piped input");
    let stderr = String::from_utf8_lossy(&output.stderr);
    assert!(
        !stderr.contains("not found") && !stderr.contains("Input BAM"),
        "downsample stdin failed the existence check; stderr: {stderr}"
    );
    assert!(output.status.success(), "downsample with piped input failed; stderr: {stderr}");
    assert!(output_bam.exists(), "Output BAM from pipe not created");

    // Oracle: run the SAME downsample directly against the sorted file (no
    // stdin). With `--fraction 1.0` and a fixed seed, the stdin path must
    // produce a byte-equivalent record stream to the direct-file path —
    // `output_bam.exists()` alone would pass even if the stdin path silently
    // dropped or reordered records.
    let direct_out = temp_dir.path().join("direct_out.bam");
    let direct = Command::new(env!("CARGO_BIN_EXE_fgumi"))
        .args([
            "downsample",
            "--input",
            sorted_bam.to_str().unwrap(),
            "--output",
            direct_out.to_str().unwrap(),
            "--fraction",
            "1.0",
            "--seed",
            "42",
            "--compression-level",
            "1",
        ])
        .output()
        .expect("Failed to run downsample with direct file input");
    assert!(
        direct.status.success(),
        "downsample with direct input failed; stderr: {}",
        String::from_utf8_lossy(&direct.stderr)
    );
    crate::helpers::parity::assert_bams_record_equivalent(&output_bam, &direct_out);
}

/// Test simplex command with piped input.
#[test]
fn test_simplex_command_with_piped_input() {
    let temp_dir = TempDir::new().expect("Failed to create temp dir");
    let input_bam = temp_dir.path().join("input.bam");
    let output_from_file = temp_dir.path().join("output_file.bam");
    let output_from_pipe = temp_dir.path().join("output_pipe.bam");

    // Create test BAM file with grouped reads (MI tag)
    create_grouped_test_bam(&input_bam);

    // Run simplex with file input
    let status = Command::new(env!("CARGO_BIN_EXE_fgumi"))
        .args([
            "simplex",
            "--input",
            input_bam.to_str().unwrap(),
            "--output",
            output_from_file.to_str().unwrap(),
            "--min-reads",
            "1",
            "--min-consensus-base-quality",
            "0",
            "--compression-level",
            "1",
            "--threads",
            "2",
        ])
        .status()
        .expect("Failed to run simplex command with file input");
    assert!(status.success(), "Simplex command with file input failed");
    assert!(output_from_file.exists(), "Output BAM from file not created");

    // Run simplex with piped input
    let cat_child = Command::new("cat")
        .arg(input_bam.to_str().unwrap())
        .stdout(Stdio::piped())
        .spawn()
        .expect("Failed to spawn cat");

    let status = Command::new(env!("CARGO_BIN_EXE_fgumi"))
        .args([
            "simplex",
            "--input",
            "-",
            "--output",
            output_from_pipe.to_str().unwrap(),
            "--min-reads",
            "1",
            "--min-consensus-base-quality",
            "0",
            "--compression-level",
            "1",
            "--threads",
            "2",
        ])
        .stdin(cat_child.stdout.unwrap())
        .status()
        .expect("Failed to run simplex command with piped input");
    assert!(status.success(), "Simplex command with piped input failed");
    assert!(output_from_pipe.exists(), "Output BAM from pipe not created");

    // Compare outputs - records should be identical (headers may differ due to @PG command line)
    compare_bam_records(&output_from_file, &output_from_pipe);
}

/// Single-threaded `simplex` (no `--threads`) must read stdin exactly once.
/// The single-threaded path previously pre-sniffed SAM-vs-BAM via `InputSource`
/// (consuming stdin) and then re-opened the input, so `-i -` failed outright and
/// `-i /dev/stdin` silently dropped the leading bytes. The multi-threaded test
/// above does not exercise this path — which is why the double-open shipped.
/// Reading the same grouped BAM from a file vs piped stdin must yield
/// byte-identical consensus records.
#[test]
fn test_simplex_single_threaded_reads_stdin_once() {
    let temp_dir = TempDir::new().expect("Failed to create temp dir");
    let input_bam = temp_dir.path().join("input.bam");
    let output_from_file = temp_dir.path().join("output_file.bam");
    let output_from_pipe = temp_dir.path().join("output_pipe.bam");

    create_grouped_test_bam(&input_bam);

    // Baseline: single-threaded (no `--threads`) from a file.
    let status = Command::new(env!("CARGO_BIN_EXE_fgumi"))
        .args([
            "simplex",
            "--input",
            input_bam.to_str().unwrap(),
            "--output",
            output_from_file.to_str().unwrap(),
            "--min-reads",
            "1",
            "--min-consensus-base-quality",
            "0",
            "--compression-level",
            "1",
        ])
        .status()
        .expect("Failed to run single-threaded simplex with file input");
    assert!(status.success(), "single-threaded simplex with file input failed");

    // Piped: `cat grouped.bam | fgumi simplex -i -` (single-threaded).
    let cat_child = Command::new("cat")
        .arg(input_bam.to_str().unwrap())
        .stdout(Stdio::piped())
        .spawn()
        .expect("Failed to spawn cat");
    let status = Command::new(env!("CARGO_BIN_EXE_fgumi"))
        .args([
            "simplex",
            "--input",
            "-",
            "--output",
            output_from_pipe.to_str().unwrap(),
            "--min-reads",
            "1",
            "--min-consensus-base-quality",
            "0",
            "--compression-level",
            "1",
        ])
        .stdin(cat_child.stdout.unwrap())
        .status()
        .expect("Failed to run single-threaded simplex with piped input");
    assert!(
        status.success(),
        "single-threaded simplex with piped input failed (stdin double-open regression?)"
    );

    compare_bam_records(&output_from_file, &output_from_pipe);
}

/// SAM (non-BGZF) input to the single-threaded codec/simplex/duplex paths must
/// fail at the header read with a `--threads` hint — those paths are BAM-only,
/// and the multi-threaded path is the SAM-capable one. Guards the `with_context`
/// message that replaced the removed SAM pre-sniff (codec/simplex) and the
/// matching hint added to duplex.
#[test]
fn test_single_threaded_consensus_rejects_sam_with_threads_hint() {
    let temp_dir = TempDir::new().expect("Failed to create temp dir");
    let bam = temp_dir.path().join("grouped.bam");
    let sam = temp_dir.path().join("grouped.sam");
    create_grouped_test_bam(&bam);
    convert_bam_to_sam(&bam, &sam);
    let input = sam.to_str().unwrap();

    // Per-command minimal args; each runs single-threaded (no `--threads`).
    let extra_args: [(&str, &[&str]); 3] = [
        ("codec", &["--min-reads", "1", "--min-duplex-length", "1"]),
        ("simplex", &["--min-reads", "1"]),
        ("duplex", &["--min-reads", "1"]),
    ];
    for (cmd, extra) in extra_args {
        let out = temp_dir.path().join(format!("{cmd}_out.bam"));
        let mut args: Vec<&str> = vec![
            cmd,
            "--input",
            input,
            "--output",
            out.to_str().unwrap(),
            "--compression-level",
            "1",
        ];
        args.extend_from_slice(extra);
        let output = Command::new(env!("CARGO_BIN_EXE_fgumi"))
            .args(&args)
            .output()
            .unwrap_or_else(|e| panic!("failed to spawn {cmd}: {e}"));
        assert!(!output.status.success(), "single-threaded {cmd} must reject SAM input");
        let stderr = String::from_utf8_lossy(&output.stderr);
        assert!(
            stderr.contains("--threads"),
            "single-threaded {cmd} SAM error should hint at --threads; stderr: {stderr}"
        );
    }
}

/// `--async-reader` over stdin must spawn the prefetch thread (the stdin branch
/// of `open_bgzf_reader` honors async like the file path) and still read the
/// input correctly — file vs piped output must match.
#[test]
fn test_simplex_single_threaded_async_reader_over_stdin() {
    let temp_dir = TempDir::new().expect("Failed to create temp dir");
    let input_bam = temp_dir.path().join("input.bam");
    let output_from_file = temp_dir.path().join("output_file.bam");
    let output_from_pipe = temp_dir.path().join("output_pipe.bam");
    create_grouped_test_bam(&input_bam);

    // Baseline: single-threaded + --async-reader from a file.
    let status = Command::new(env!("CARGO_BIN_EXE_fgumi"))
        .args([
            "simplex",
            "--input",
            input_bam.to_str().unwrap(),
            "--output",
            output_from_file.to_str().unwrap(),
            "--min-reads",
            "1",
            "--min-consensus-base-quality",
            "0",
            "--compression-level",
            "1",
            "--async-reader",
        ])
        .status()
        .expect("Failed to run async-reader simplex with file input");
    assert!(status.success(), "async-reader simplex from file failed");

    // Piped: `cat grouped.bam | fgumi simplex -i - --async-reader` (single-threaded).
    let cat_child = Command::new("cat")
        .arg(input_bam.to_str().unwrap())
        .stdout(Stdio::piped())
        .spawn()
        .expect("Failed to spawn cat");
    let status = Command::new(env!("CARGO_BIN_EXE_fgumi"))
        .args([
            "simplex",
            "--input",
            "-",
            "--output",
            output_from_pipe.to_str().unwrap(),
            "--min-reads",
            "1",
            "--min-consensus-base-quality",
            "0",
            "--compression-level",
            "1",
            "--async-reader",
        ])
        .stdin(cat_child.stdout.unwrap())
        .status()
        .expect("Failed to run async-reader simplex with piped input");
    assert!(
        status.success(),
        "async-reader simplex from stdin failed (prefetch-on-stdin regression?)"
    );

    compare_bam_records(&output_from_file, &output_from_pipe);
}

/// Convert a BAM file to a parallel SAM text file containing the same
/// records. Lets `fgumi group --input foo.sam` reuse the same record
/// set as the BAM-input baseline for parity tests.
fn convert_bam_to_sam(bam_path: &PathBuf, sam_path: &PathBuf) {
    use noodles::sam::alignment::io::Write as _;
    let mut bam_reader = bam::io::reader::Builder.build_from_path(bam_path).expect("open bam");
    let header = bam_reader.read_header().expect("read bam header");
    let mut sam_writer =
        noodles::sam::io::Writer::new(fs::File::create(sam_path).expect("create sam"));
    sam_writer.write_header(&header).expect("write sam header");
    for record in bam_reader.records() {
        let record = record.expect("bam record");
        sam_writer.write_alignment_record(&header, &record).expect("write sam record");
    }
}

/// SAM-input parity for the typed-step pipeline. Generates a SAM file
/// from the same record set as the BAM baseline, runs `fgumi group
/// --input foo.sam`, and compares the output BAMs.
/// Exercises the `ReadSamChunks` + `ParseSamChunk` parallel-parse path.
#[test]
fn test_group_command_with_sam_input_new_pipeline_matches_bam_baseline() {
    let temp_dir = TempDir::new().expect("Failed to create temp dir");
    let bam_path = temp_dir.path().join("input.bam");
    let sam_path = temp_dir.path().join("input.sam");
    let out_bam_baseline = temp_dir.path().join("output_bam.bam");
    let out_sam = temp_dir.path().join("output_sam.bam");

    create_test_input_bam(&bam_path);
    convert_bam_to_sam(&bam_path, &sam_path);

    let status = Command::new(env!("CARGO_BIN_EXE_fgumi"))
        .args([
            "group",
            "--input",
            bam_path.to_str().unwrap(),
            "--output",
            out_bam_baseline.to_str().unwrap(),
            "--strategy",
            "identity",
            "--edits",
            "0",
            "--compression-level",
            "1",
            "--threads",
            "2",
        ])
        .status()
        .expect("Failed to run group with BAM input");
    assert!(status.success(), "BAM-input baseline failed");

    let status = Command::new(env!("CARGO_BIN_EXE_fgumi"))
        .args([
            "group",
            "--input",
            sam_path.to_str().unwrap(),
            "--output",
            out_sam.to_str().unwrap(),
            "--strategy",
            "identity",
            "--edits",
            "0",
            "--compression-level",
            "1",
            "--threads",
            "2",
        ])
        .status()
        .expect("Failed to run group with SAM input");
    assert!(status.success(), "SAM-input run failed");

    compare_bam_records(&out_bam_baseline, &out_sam);
}

/// SAM-input parity for simplex. Generates a SAM from the same grouped
/// records as the BAM baseline (with MI tags), runs simplex on both,
/// compares the consensus outputs.
#[test]
fn test_simplex_command_with_sam_input_new_pipeline_matches_bam_baseline() {
    let temp_dir = TempDir::new().expect("Failed to create temp dir");
    let bam_path = temp_dir.path().join("input.bam");
    let sam_path = temp_dir.path().join("input.sam");
    let out_bam_baseline = temp_dir.path().join("output_bam.bam");
    let out_sam = temp_dir.path().join("output_sam.bam");

    create_grouped_test_bam(&bam_path);
    convert_bam_to_sam(&bam_path, &sam_path);

    let status = Command::new(env!("CARGO_BIN_EXE_fgumi"))
        .args([
            "simplex",
            "--input",
            bam_path.to_str().unwrap(),
            "--output",
            out_bam_baseline.to_str().unwrap(),
            "--min-reads",
            "1",
            "--min-consensus-base-quality",
            "0",
            "--compression-level",
            "1",
            "--threads",
            "2",
        ])
        .status()
        .expect("BAM baseline run");
    assert!(status.success(), "simplex BAM baseline failed");

    let status = Command::new(env!("CARGO_BIN_EXE_fgumi"))
        .args([
            "simplex",
            "--input",
            sam_path.to_str().unwrap(),
            "--output",
            out_sam.to_str().unwrap(),
            "--min-reads",
            "1",
            "--min-consensus-base-quality",
            "0",
            "--compression-level",
            "1",
            "--threads",
            "2",
        ])
        .status()
        .expect("SAM run");
    assert!(status.success(), "simplex SAM run failed");

    compare_bam_records(&out_bam_baseline, &out_sam);
}

/// SAM-input parity for `correct` — exercises the `(Sam, true)` arm of
/// correct.rs's 4-way match (the `into_multi()` shape that emits both
/// kept records and rejects on parallel branches).
#[test]
fn test_correct_command_with_sam_input_new_pipeline_with_rejects_matches_bam_baseline() {
    let temp_dir = TempDir::new().expect("Failed to create temp dir");
    let bam_path = temp_dir.path().join("input.bam");
    let sam_path = temp_dir.path().join("input.sam");
    let umi_file = temp_dir.path().join("umis.txt");
    let out_bam_baseline = temp_dir.path().join("output_bam.bam");
    let out_bam_rejects_baseline = temp_dir.path().join("rejects_bam.bam");
    let out_sam = temp_dir.path().join("output_sam.bam");
    let out_sam_rejects = temp_dir.path().join("rejects_sam.bam");

    create_test_input_bam(&bam_path);
    convert_bam_to_sam(&bam_path, &sam_path);
    // The UMIs used by `create_umi_family` are AAAAAAAA and CCCCCCCC.
    std::fs::write(&umi_file, b"AAAAAAAA\nCCCCCCCC\n").expect("write umi file");

    for (input, output, rejects, label) in [
        (&bam_path, &out_bam_baseline, &out_bam_rejects_baseline, "BAM"),
        (&sam_path, &out_sam, &out_sam_rejects, "SAM"),
    ] {
        let status = Command::new(env!("CARGO_BIN_EXE_fgumi"))
            .args([
                "correct",
                "--input",
                input.to_str().unwrap(),
                "--output",
                output.to_str().unwrap(),
                "--rejects",
                rejects.to_str().unwrap(),
                "--umi-files",
                umi_file.to_str().unwrap(),
                "--max-mismatches",
                "1",
                "--min-distance",
                "1",
                "--compression-level",
                "1",
                "--threads",
                "2",
            ])
            .status()
            .unwrap_or_else(|_| panic!("Failed to run correct {label}"));
        assert!(status.success(), "correct {label} failed");
    }
    // Non-vacuous guard (S9b-007): this is a SAM-vs-BAM parity check between two
    // fgumi outputs, so a both-empty pass is possible if correct silently dropped
    // every record on BOTH paths. Assert the BAM baseline kept output is
    // non-empty before the parity comparison so the test cannot pass vacuously.
    // (The rejects output is legitimately empty here — this fixture corrects all
    // UMIs, rejecting none — so only the kept output is guarded.)
    assert!(count_bam_records(&out_bam_baseline) > 0, "BAM baseline kept output is empty");
    compare_bam_records(&out_bam_baseline, &out_sam);
    compare_bam_records(&out_bam_rejects_baseline, &out_sam_rejects);
}

/// Build an unmapped-consensus BAM with depth tags so filter accepts it
/// without `--ref` (mapped reads would require ref for NM/UQ/MD tag
/// regeneration). Mirrors the fixture in `test_filter_command.rs`.
fn create_unmapped_consensus_bam(path: &PathBuf) {
    use fgumi_lib::sam::SamTag;
    use fgumi_raw_bam::{SamBuilder as RawSamBuilder, flags};

    let header = create_minimal_header("chr1", 10000);
    let records: Vec<_> = ["good1", "good2", "low_depth"]
        .iter()
        .map(|name| {
            let depth = if *name == "low_depth" { 1u16 } else { 10u16 };
            let mut b = RawSamBuilder::new();
            b.read_name(name.as_bytes())
                .flags(flags::UNMAPPED)
                .sequence(b"ACGTACGT")
                .qualities(&[35; 8]);
            b.add_array_u16(SamTag::CD_BASES, &[depth; 8]).add_array_u16(SamTag::CE_BASES, &[0; 8]);
            b.build()
        })
        .collect();

    let mut writer =
        bam::io::Writer::new(fs::File::create(path).expect("Failed to create BAM file"));
    writer.write_header(&header).expect("Failed to write header");
    for raw in records {
        let rec =
            fgumi_raw_bam::raw_record_to_record_buf(&raw, &header).expect("raw -> record_buf");
        writer.write_alignment_record(&header, &rec).expect("Failed to write record");
    }
    writer.try_finish().expect("Failed to finish BAM");
}

/// SAM-input parity for `filter --rejects` — exercises the `(Sam, ...)`
/// arm of `execute_single_read_with_rejects_new` (the `into_multi()`
/// fan-out shape that emits kept + rejected on parallel branches).
#[test]
fn test_filter_command_with_rejects_sam_input_matches_bam_baseline() {
    let temp_dir = TempDir::new().expect("Failed to create temp dir");
    let bam_path = temp_dir.path().join("input.bam");
    let sam_path = temp_dir.path().join("input.sam");
    let out_bam_baseline = temp_dir.path().join("output_bam.bam");
    let out_bam_rejects_baseline = temp_dir.path().join("rejects_bam.bam");
    let out_sam = temp_dir.path().join("output_sam.bam");
    let out_sam_rejects = temp_dir.path().join("rejects_sam.bam");

    create_unmapped_consensus_bam(&bam_path);
    convert_bam_to_sam(&bam_path, &sam_path);

    for (input, output, rejects, label) in [
        (&bam_path, &out_bam_baseline, &out_bam_rejects_baseline, "BAM"),
        (&sam_path, &out_sam, &out_sam_rejects, "SAM"),
    ] {
        let status = Command::new(env!("CARGO_BIN_EXE_fgumi"))
            .args([
                "filter",
                "--input",
                input.to_str().unwrap(),
                "--output",
                output.to_str().unwrap(),
                "--rejects",
                rejects.to_str().unwrap(),
                "--min-mean-base-quality",
                "0",
                "--min-reads",
                "3",
                "--compression-level",
                "1",
                "--threads",
                "2",
            ])
            .status()
            .unwrap_or_else(|_| panic!("Failed to run filter {label}"));
        assert!(status.success(), "filter {label} failed");
    }
    // Non-vacuous guard (S9b-007): this compares two fgumi outputs, so a
    // both-empty pass is possible if filter silently dropped every record on
    // BOTH paths. The fixture has two depth-10 reads (kept) and one depth-1 read
    // (rejected at --min-reads 3), so both baseline outputs must be non-empty.
    assert!(count_bam_records(&out_bam_baseline) > 0, "BAM baseline kept output is empty");
    assert!(
        count_bam_records(&out_bam_rejects_baseline) > 0,
        "BAM baseline rejects output is empty"
    );
    compare_bam_records(&out_bam_baseline, &out_sam);
    compare_bam_records(&out_bam_rejects_baseline, &out_sam_rejects);
}

/// SAM-on-stdin parity (`cat foo.sam | fgumi group --input -`). The
/// `InputSource::open` magic-byte peek picks the SAM-chunks source for
/// stdin (no `\x1f` BGZF magic, falls through to text parsing).
#[test]
fn test_group_command_with_sam_stdin_new_pipeline_matches_bam_baseline() {
    let temp_dir = TempDir::new().expect("Failed to create temp dir");
    let bam_path = temp_dir.path().join("input.bam");
    let sam_path = temp_dir.path().join("input.sam");
    let out_bam_baseline = temp_dir.path().join("output_bam.bam");
    let out_sam_stdin = temp_dir.path().join("output_sam_stdin.bam");

    create_test_input_bam(&bam_path);
    convert_bam_to_sam(&bam_path, &sam_path);

    // BAM baseline (file input).
    let status = Command::new(env!("CARGO_BIN_EXE_fgumi"))
        .args([
            "group",
            "--input",
            bam_path.to_str().unwrap(),
            "--output",
            out_bam_baseline.to_str().unwrap(),
            "--strategy",
            "identity",
            "--edits",
            "0",
            "--compression-level",
            "1",
            "--threads",
            "2",
        ])
        .status()
        .expect("BAM baseline run");
    assert!(status.success(), "BAM-input baseline failed");

    // SAM piped through stdin.
    let cat_child = Command::new("cat")
        .arg(sam_path.to_str().unwrap())
        .stdout(Stdio::piped())
        .spawn()
        .expect("spawn cat");
    let status = Command::new(env!("CARGO_BIN_EXE_fgumi"))
        .args([
            "group",
            "--input",
            "-",
            "--output",
            out_sam_stdin.to_str().unwrap(),
            "--strategy",
            "identity",
            "--edits",
            "0",
            "--compression-level",
            "1",
            "--threads",
            "2",
        ])
        .stdin(cat_child.stdout.unwrap())
        .status()
        .expect("SAM stdin run");
    assert!(status.success(), "SAM-stdin run failed");

    compare_bam_records(&out_bam_baseline, &out_sam_stdin);
}

/// Same shape as `test_group_command_with_piped_input` but routed through
/// the typed-step pipeline. Exercises the
/// `read_bam_auto` dispatcher introduced for issue #330 Phase 1 T1.3.4:
/// when the input path is `-`, the source step is `read_bam_stdin`
/// (which buffers the BAM header via `TeeReader` and replays it ahead
/// of the remaining stdin stream).
#[test]
fn test_group_command_with_piped_input_new_pipeline() {
    let temp_dir = TempDir::new().expect("Failed to create temp dir");
    let input_bam = temp_dir.path().join("input.bam");
    let output_from_file = temp_dir.path().join("output_file.bam");
    let output_from_pipe = temp_dir.path().join("output_pipe.bam");

    create_test_input_bam(&input_bam);

    let status = Command::new(env!("CARGO_BIN_EXE_fgumi"))
        .args([
            "group",
            "--input",
            input_bam.to_str().unwrap(),
            "--output",
            output_from_file.to_str().unwrap(),
            "--strategy",
            "identity",
            "--edits",
            "0",
            "--compression-level",
            "1",
            "--threads",
            "2",
        ])
        .status()
        .expect("Failed to run group with file input");
    assert!(status.success(), "group (file) failed");
    assert!(output_from_file.exists(), "file-input output not created");

    let cat_child = Command::new("cat")
        .arg(input_bam.to_str().unwrap())
        .stdout(Stdio::piped())
        .spawn()
        .expect("Failed to spawn cat");

    let status = Command::new(env!("CARGO_BIN_EXE_fgumi"))
        .args([
            "group",
            "--input",
            "-",
            "--output",
            output_from_pipe.to_str().unwrap(),
            "--strategy",
            "identity",
            "--edits",
            "0",
            "--compression-level",
            "1",
            "--threads",
            "2",
        ])
        .stdin(cat_child.stdout.unwrap())
        .status()
        .expect("Failed to run group with piped input");
    assert!(status.success(), "group (stdin) failed");
    assert!(output_from_pipe.exists(), "stdin output not created");

    compare_bam_records(&output_from_file, &output_from_pipe);
}

/// Helper to compare BAM records (ignoring header differences like @PG command line).
///
/// Decodes to eager `RecordBuf`s and compares them for full equality — every
/// field including bases, qualities, flags, CIGAR, positions, and aux tags — so
/// a record-level corruption from mishandled input (dropped or garbled bytes)
/// fails here, not just a count/name mismatch.
fn compare_bam_records(path1: &PathBuf, path2: &PathBuf) {
    use noodles::sam::alignment::RecordBuf;

    let mut reader1 =
        bam::io::reader::Builder.build_from_path(path1).expect("Failed to open BAM 1");
    let header1 = reader1.read_header().expect("Failed to read header 1");

    let mut reader2 =
        bam::io::reader::Builder.build_from_path(path2).expect("Failed to open BAM 2");
    let header2 = reader2.read_header().expect("Failed to read header 2");

    let records1: Vec<RecordBuf> =
        reader1.record_bufs(&header1).map(|r| r.expect("Failed to read record 1")).collect();
    let records2: Vec<RecordBuf> =
        reader2.record_bufs(&header2).map(|r| r.expect("Failed to read record 2")).collect();

    assert_eq!(records1.len(), records2.len(), "BAM files have different record counts");

    for (i, (r1, r2)) in records1.iter().zip(records2.iter()).enumerate() {
        assert_eq!(r1, r2, "Record {i} differs (full record identity)");
    }
}

/// Create a test BAM with UMI families for grouping.
fn create_test_input_bam(path: &PathBuf) {
    let header = create_minimal_header("chr1", 10000);

    let mut writer =
        bam::io::Writer::new(fs::File::create(path).expect("Failed to create BAM file"));
    writer.write_header(&header).expect("Failed to write header");

    // Create multiple UMI families
    let family1 = create_umi_family("AAAAAAAA", 5, "family1", "ACGTACGT", 30);
    let family2 = create_umi_family("CCCCCCCC", 3, "family2", "TGCATGCA", 30);

    for record in family1.iter().chain(family2.iter()) {
        writer
            .write_alignment_record(&header, &to_record_buf(record))
            .expect("Failed to write record");
    }

    writer.try_finish().expect("Failed to finish BAM");
}

/// Create a test BAM with already-grouped reads (MI tag set).
fn create_grouped_test_bam(path: &PathBuf) {
    use noodles::sam::alignment::record::data::field::Tag;
    use noodles::sam::alignment::record_buf::data::field::Value;

    let header = create_minimal_header("chr1", 10000);

    let mut writer =
        bam::io::Writer::new(fs::File::create(path).expect("Failed to create BAM file"));
    writer.write_header(&header).expect("Failed to write header");

    // Create grouped reads with MI tag
    let mi_tag = Tag::from(fgumi_lib::sam::SamTag::MI);
    let records = create_umi_family("AAAAAAAA", 5, "mol1", "ACGTACGT", 30);

    // Add MI tag to each record (convert to RecordBuf first to enable tag mutation)
    let mi_value = 1;
    for raw in &records {
        let mut rec = to_record_buf(raw);
        rec.data_mut().insert(mi_tag, Value::from(mi_value));
        writer.write_alignment_record(&header, &rec).expect("Failed to write record");
    }

    // Second family with different MI
    let mi_value2 = 2;
    let records2 = create_umi_family("CCCCCCCC", 3, "mol2", "TGCATGCA", 30);
    for raw in &records2 {
        let mut rec = to_record_buf(raw);
        rec.data_mut().insert(mi_tag, Value::from(mi_value2));
        writer.write_alignment_record(&header, &rec).expect("Failed to write record");
    }

    writer.try_finish().expect("Failed to finish BAM");
}

// ===========================================================================
// Comprehensive stdin-input coverage (issue #330 follow-up).
//
// Every command that accepts a BAM input must read stdin exactly once and
// produce output byte-identical to reading the same BAM from a file. These
// tests assert file-vs-stdin parity for both `-` and `/dev/stdin`, exercising
// the single-threaded and multi-threaded reader paths where a command has both
// (codec/duplex/clip have a single-threaded fast path; the rest route through
// the chain builder regardless). Parity (not correctness) is the contract: the
// command need only run to completion and emit the same records either way.
// ===========================================================================

/// Run `fgumi` with `args`; if `stdin_bam` is `Some`, pipe that BAM to its
/// stdin via `cat`. Asserts the process exits successfully.
fn run_fgumi_checked(args: &[String], stdin_bam: Option<&std::path::Path>) {
    let mut command = Command::new(env!("CARGO_BIN_EXE_fgumi"));
    command.args(args);
    let status = if let Some(bam) = stdin_bam {
        let cat = Command::new("cat")
            .arg(bam.to_str().unwrap())
            .stdout(Stdio::piped())
            .spawn()
            .expect("Failed to spawn cat");
        command.stdin(cat.stdout.unwrap()).status()
    } else {
        command.status()
    };
    let status = status.unwrap_or_else(|e| panic!("Failed to run fgumi {args:?}: {e}"));
    assert!(status.success(), "fgumi {args:?} failed (stdin_piped={})", stdin_bam.is_some());
}

/// Assert a command emits byte-identical records reading `input_bam` from a
/// file vs piped stdin (`-`) vs `/dev/stdin`. `build_args(input, output)`
/// returns the full argv for one invocation; `label` namespaces temp outputs.
fn assert_stdin_input_parity(
    input_bam: &std::path::Path,
    temp: &std::path::Path,
    label: &str,
    build_args: impl Fn(&str, &str) -> Vec<String>,
) {
    let out_file = temp.join(format!("{label}_file.bam"));
    let out_dash = temp.join(format!("{label}_dash.bam"));
    let out_dev = temp.join(format!("{label}_dev.bam"));

    run_fgumi_checked(&build_args(input_bam.to_str().unwrap(), out_file.to_str().unwrap()), None);
    run_fgumi_checked(&build_args("-", out_dash.to_str().unwrap()), Some(input_bam));
    run_fgumi_checked(&build_args("/dev/stdin", out_dev.to_str().unwrap()), Some(input_bam));

    // Guard against a vacuous pass: if the command emitted no records, a
    // file-vs-stdin parity check would be trivially true even if stdin were
    // dropped entirely. Require the file baseline to produce real output.
    assert!(
        count_bam_records(&out_file) > 0,
        "{label}: file baseline produced no records — parity check would be vacuous"
    );

    compare_bam_records(&out_file, &out_dash);
    compare_bam_records(&out_file, &out_dev);
}

/// Count records in a BAM (eager decode). Used to reject vacuous parity tests.
fn count_bam_records(path: &std::path::Path) -> usize {
    let mut reader = bam::io::reader::Builder.build_from_path(path).expect("open BAM for count");
    let header = reader.read_header().expect("read header for count");
    let mut count = 0;
    for record in reader.record_bufs(&header) {
        record.expect("read record for count");
        count += 1;
    }
    count
}

fn s(v: &str) -> String {
    v.to_string()
}

#[rstest]
#[case("1")]
#[case("2")]
fn test_sort_reads_stdin_once(#[case] threads: &str) {
    let dir = TempDir::new().unwrap();
    let input = dir.path().join("input.bam");
    create_test_input_bam(&input);
    assert_stdin_input_parity(&input, dir.path(), "sort", |i, o| {
        vec![
            s("sort"),
            s("--input"),
            s(i),
            s("--output"),
            s(o),
            s("--order"),
            s("coordinate"),
            s("--threads"),
            s(threads),
            s("--compression-level"),
            s("1"),
        ]
    });
}

/// `sort --verify` reads its input twice (header probe + a fresh record
/// re-scan), which a non-seekable stdin stream can't satisfy, so it must be
/// rejected up front with a clear message rather than failing deep in a
/// re-open. (Plain `sort` streams stdin once and IS supported — above.)
#[test]
fn test_sort_verify_rejects_stdin_with_clear_message() {
    let dir = TempDir::new().unwrap();
    let input = dir.path().join("input.bam");
    create_test_input_bam(&input);
    let cat = Command::new("cat")
        .arg(input.to_str().unwrap())
        .stdout(Stdio::piped())
        .spawn()
        .expect("spawn cat");
    let output = Command::new(env!("CARGO_BIN_EXE_fgumi"))
        .args(["sort", "--verify", "--input", "-", "--order", "coordinate"])
        .stdin(cat.stdout.unwrap())
        .output()
        .expect("run sort --verify -i -");
    assert!(!output.status.success(), "sort --verify from stdin must be rejected");
    let stderr = String::from_utf8_lossy(&output.stderr);
    assert!(
        stderr.contains("stdin"),
        "sort --verify stdin rejection should mention stdin; stderr: {stderr}"
    );
}

/// Push `--threads <n>` only when `threads` is `Some`; `None` exercises a
/// command's single-threaded fast path (where it has one).
fn push_threads(args: &mut Vec<String>, threads: Option<&str>) {
    if let Some(t) = threads {
        args.push(s("--threads"));
        args.push(s(t));
    }
}

#[rstest]
#[case("1")]
#[case("2")]
fn test_group_reads_stdin_once(#[case] threads: &str) {
    let dir = TempDir::new().unwrap();
    let input = dir.path().join("input.bam");
    create_test_input_bam(&input);
    assert_stdin_input_parity(&input, dir.path(), "group", |i, o| {
        vec![
            s("group"),
            s("--input"),
            s(i),
            s("--output"),
            s(o),
            s("--strategy"),
            s("identity"),
            s("--edits"),
            s("0"),
            s("--threads"),
            s(threads),
            s("--compression-level"),
            s("1"),
        ]
    });
}

#[rstest]
#[case("1")]
#[case("2")]
fn test_dedup_reads_stdin_once(#[case] threads: &str) {
    let dir = TempDir::new().unwrap();
    let input = dir.path().join("input.bam");
    create_test_input_bam(&input);
    assert_stdin_input_parity(&input, dir.path(), "dedup", |i, o| {
        vec![
            s("dedup"),
            s("--input"),
            s(i),
            s("--output"),
            s(o),
            s("--strategy"),
            s("identity"),
            s("--edits"),
            s("0"),
            s("--threads"),
            s(threads),
            s("--compression-level"),
            s("1"),
        ]
    });
}

#[rstest]
#[case("1")]
#[case("2")]
fn test_correct_reads_stdin_once(#[case] threads: &str) {
    let dir = TempDir::new().unwrap();
    let input = dir.path().join("input.bam");
    create_test_input_bam(&input);
    // create_test_input_bam tags reads with RX = AAAAAAAA / CCCCCCCC.
    assert_stdin_input_parity(&input, dir.path(), "correct", |i, o| {
        vec![
            s("correct"),
            s("--input"),
            s(i),
            s("--output"),
            s(o),
            s("--umis"),
            s("AAAAAAAA"),
            s("--umis"),
            s("CCCCCCCC"),
            s("--min-distance"),
            s("1"),
            s("--threads"),
            s(threads),
            s("--compression-level"),
            s("1"),
        ]
    });
}

/// codec has a single-threaded fast path (no `--threads`) plus the
/// multi-threaded chain path; cover both.
// Uses helpers from the consensus-gated `test_codec_command` module.
#[cfg(feature = "consensus")]
#[rstest]
#[case::single_threaded(None)]
#[case::multi_threaded(Some("2"))]
fn test_codec_reads_stdin_once(#[case] threads: Option<&str>) {
    let dir = TempDir::new().unwrap();
    let input = dir.path().join("codec.bam");
    // CODEC needs paired duplex reads (R1/R2 same molecule) to emit consensus;
    // reuse the proven builder so the parity check is non-vacuous.
    let pairs: Vec<_> = (0..3)
        .map(|i| {
            crate::test_codec_command::create_codec_read_pair(
                &format!("read{i}"),
                b"ACGTACGT",
                b"ACGTACGT",
                &[30; 8],
                &[30; 8],
                100,
                "UMI001",
                None,
            )
        })
        .collect();
    crate::test_codec_command::create_codec_test_bam(&input, pairs);
    assert_stdin_input_parity(&input, dir.path(), "codec", move |i, o| {
        let mut args = vec![
            s("codec"),
            s("--input"),
            s(i),
            s("--output"),
            s(o),
            s("--min-reads"),
            s("1"),
            s("--min-duplex-length"),
            s("1"),
            s("--compression-level"),
            s("1"),
        ];
        push_threads(&mut args, threads);
        args
    });
}

/// duplex has a single-threaded fast path (no `--threads`) plus the
/// multi-threaded chain path; cover both.
// Uses helpers from the consensus-gated `test_duplex_command` module.
#[cfg(feature = "consensus")]
#[rstest]
#[case::single_threaded(None)]
#[case::multi_threaded(Some("2"))]
fn test_duplex_reads_stdin_once(#[case] threads: Option<&str>) {
    let dir = TempDir::new().unwrap();
    let input = dir.path().join("duplex.bam");
    // Duplex needs A/B-strand grouped molecules to emit consensus; reuse the
    // proven builder so the parity check is non-vacuous.
    let molecule =
        crate::test_duplex_command::create_duplex_molecule("MI001", "ACGTACGT", 30, 100, 3);
    crate::test_duplex_command::create_duplex_bam(&input, vec![molecule]);
    assert_stdin_input_parity(&input, dir.path(), "duplex", move |i, o| {
        let mut args = vec![
            s("duplex"),
            s("--input"),
            s(i),
            s("--output"),
            s(o),
            s("--min-reads"),
            s("1"),
            s("--compression-level"),
            s("1"),
        ];
        push_threads(&mut args, threads);
        args
    });
}

/// Unlike codec/duplex/clip, filter has NO single-threaded fast path: `Filter::execute`
/// always routes through `chains::build_for(spec).run()` regardless of `--threads`
/// (see `src/lib/commands/filter.rs`). So there is no no-`--threads` code path to cover —
/// only the explicit `--threads 1` / `--threads 2` cases exist (no
/// `#[case::single_threaded(None)]`).
#[rstest]
#[case("1")]
#[case("2")]
fn test_filter_reads_stdin_once(#[case] threads: &str) {
    use crate::helpers::bam_generator::create_test_reference;
    let dir = TempDir::new().unwrap();
    let reference = create_test_reference(dir.path()).to_str().unwrap().to_string();
    let input = dir.path().join("consensus.bam");
    // Mapped consensus reads with CD/CE per-base tags that pass the filter;
    // reuse the proven builder so the parity check is non-vacuous.
    crate::test_filter_command::write_filter_consensus_bam(&input);

    assert_stdin_input_parity(&input, dir.path(), "filter", move |i, o| {
        // --ref is required when filtering mapped reads (NM/UQ/MD).
        vec![
            s("filter"),
            s("--input"),
            s(i),
            s("--output"),
            s(o),
            s("--ref"),
            reference.clone(),
            s("--min-reads"),
            s("1"),
            s("--max-no-call-fraction"),
            s("1.0"),
            s("--threads"),
            s(threads),
            s("--compression-level"),
            s("1"),
        ]
    });
}

/// clip has a single-threaded fast path (no `--threads`) plus the
/// multi-threaded chain path; cover both.
#[rstest]
#[case::single_threaded(None)]
#[case::multi_threaded(Some("2"))]
fn test_clip_reads_stdin_once(#[case] threads: Option<&str>) {
    use crate::helpers::bam_generator::{create_paired_umi_family, create_test_reference};
    let dir = TempDir::new().unwrap();
    let reference = create_test_reference(dir.path()).to_str().unwrap().to_string();
    let input = dir.path().join("paired.bam");
    let header = create_minimal_header("chr1", 10000);
    let mut writer = bam::io::Writer::new(fs::File::create(&input).expect("create paired BAM"));
    writer.write_header(&header).expect("write header");
    for raw in &create_paired_umi_family("ACGTACGT", 3, "clip", "ACGTACGTAC", "TGCATGCATG", 30) {
        writer.write_alignment_record(&header, &to_record_buf(raw)).expect("write paired record");
    }
    writer.try_finish().expect("finish BAM");

    assert_stdin_input_parity(&input, dir.path(), "clip", move |i, o| {
        let mut args = vec![
            s("clip"),
            s("--input"),
            s(i),
            s("--output"),
            s(o),
            s("--reference"),
            reference.clone(),
            s("--clip-overlapping-reads"),
            s("--compression-level"),
            s("1"),
        ];
        push_threads(&mut args, threads);
        args
    });
}
