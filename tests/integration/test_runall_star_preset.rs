//! `runall --aligner::preset star` against a stand-in `STAR` on `PATH`.
//!
//! The stand-in checks the streaming flags the preset must pass, reads the
//! unaligned SAM pairs it is fed on stdin, and reports them under the `@SQ`
//! lines of its genome directory's `chrNameLength.txt`: the first pair mapped
//! with a spliced R1, the second with only R1 mapped, and every other read
//! unmapped. No real STAR is needed.

use std::path::{Path, PathBuf};
use std::process::{Command, Output, Stdio};

use fgumi_raw_bam::{RawRecord, SamBuilder, flags};
use tempfile::TempDir;

use crate::helpers::bam_generator::{create_minimal_header, create_test_reference, write_bam};
use crate::helpers::cli::utf8_path as p;
use crate::helpers::read_bam_output;

/// The stand-in `STAR`. Exits 2 when a flag the preset owns is missing or
/// wrong, or when `--outTmpDir` already exists (STAR refuses that); exits 4 on
/// input that is not unaligned SAM pairs. Like STAR, which maps a chunk of
/// input before it writes that chunk's alignments, it reads all of its input
/// before writing any record (a stand-in that wrote each record as it read it
/// could fill its stdout pipe while the align stage is still writing a batch).
/// Appends its arguments to
/// `$FAKE_STAR_ARGS_LOG` when set, writes `Log.out` under
/// `--outFileNamePrefix`, and sleeps `$FAKE_STAR_LOAD_SECS` before reading its
/// input, standing in for a genome load. With `$FAKE_STAR_EARLY_HEADER` set it
/// writes its header before that sleep instead of after it.
const FAKE_STAR: &str = r#"#!/bin/sh
fail() { echo "fake STAR: $*" >&2; exit 2; }
all_args="$*"
while [ $# -gt 0 ]; do
  case "$1" in
    --runMode) [ "$2" = alignReads ] || fail "--runMode $2"; run_mode=1; shift 2 ;;
    --runThreadN) threads=$2; shift 2 ;;
    --genomeDir) genome=$2; shift 2 ;;
    --readFilesType) [ "$2 $3" = "SAM PE" ] || fail "--readFilesType $2 $3"; files_type=1; shift 3 ;;
    --readFilesIn) [ "$2" = /dev/stdin ] || fail "--readFilesIn $2"; files_in=1; shift 2 ;;
    --readNameSeparator) [ "$2" = none ] || fail "--readNameSeparator $2"; separator=1; shift 2 ;;
    --outStd) [ "$2" = SAM ] || fail "--outStd $2"; out_std=1; shift 2 ;;
    --outSAMtype) [ "$2" = SAM ] || fail "--outSAMtype $2"; sam_type=1; shift 2 ;;
    --outSAMorder) [ "$2" = PairedKeepInputOrder ] || fail "--outSAMorder $2"; order=1; shift 2 ;;
    --outSAMunmapped) [ "$2" = Within ] || fail "--outSAMunmapped $2"; unmapped=1; shift 2 ;;
    --outTmpDir) tmp_dir=$2; shift 2 ;;
    --outFileNamePrefix) prefix=$2; shift 2 ;;
    *) shift ;;
  esac
done
for flag in run_mode files_type files_in separator out_std sam_type order unmapped; do
  eval "[ -n \"\$$flag\" ]" || fail "missing flag: $flag"
done
[ -n "$threads" ] || fail "missing --runThreadN"
[ -f "$genome/chrNameLength.txt" ] || fail "no chrNameLength.txt in --genomeDir $genome"
[ -n "$tmp_dir" ] || fail "missing --outTmpDir"
[ -e "$tmp_dir" ] && fail "--outTmpDir $tmp_dir already exists"
[ -n "$prefix" ] || fail "missing --outFileNamePrefix"
mkdir "$tmp_dir" || fail "cannot create --outTmpDir $tmp_dir"
echo "$all_args" > "${prefix}Log.out" || fail "cannot write ${prefix}Log.out"
[ -n "$FAKE_STAR_ARGS_LOG" ] && echo "$all_args" >> "$FAKE_STAR_ARGS_LOG"
header() {
  awk -v sizes="$genome/chrNameLength.txt" 'BEGIN {
    OFS = "\t"; print "@HD", "VN:1.4"
    while ((getline line < sizes) > 0) { split(line, f, "\t"); print "@SQ", "SN:" f[1], "LN:" f[2] }
    print "@PG", "ID:STAR", "PN:STAR"
  }'
}
[ -n "$FAKE_STAR_EARLY_HEADER" ] && header
# Detached from the pipes, so killing this script closes them.
sleep "${FAKE_STAR_LOAD_SECS:-0}" </dev/null >/dev/null 2>&1
[ -z "$FAKE_STAR_EARLY_HEADER" ] && header
awk 'BEGIN { OFS = "\t" }
  /^@/ { next }
  NF != 11 || ($2 != 77 && $2 != 141) { exit 4 }
  {
    pair = int(n / 2); n++
    if (pair == 0 && $2 == 77) rec[n] = $1 OFS 65 OFS "chr1" OFS 100 OFS 255 OFS "10M100N10M" OFS "=" OFS 400 OFS 0 OFS $10 OFS $11 OFS "NH:i:1" OFS "HI:i:1"
    else if (pair == 0) rec[n] = $1 OFS 129 OFS "chr1" OFS 400 OFS 255 OFS "20M" OFS "=" OFS 100 OFS 0 OFS $10 OFS $11 OFS "NH:i:1" OFS "HI:i:1"
    else if (pair == 1 && $2 == 77) rec[n] = $1 OFS 73 OFS "chr1" OFS 200 OFS 255 OFS "20M" OFS "=" OFS 200 OFS 0 OFS $10 OFS $11 OFS "NH:i:1" OFS "HI:i:1"
    else if (pair == 1) rec[n] = $1 OFS 133 OFS "chr1" OFS 200 OFS 0 OFS "*" OFS "=" OFS 200 OFS 0 OFS $10 OFS $11 OFS "NH:i:0"
    else rec[n] = $1 OFS $2 OFS "*\t0\t0\t*\t*\t0\t0" OFS $10 OFS $11 OFS "NH:i:0"
  }
  END { for (i = 1; i <= n; i++) print rec[i] }'
status=$?
rm -rf "$tmp_dir"
exit $status
"#;

/// The flags `--aligner::preset star` passes, spelled out for the
/// `--aligner::command` route (without the per-run paths).
const STAR_STREAMING_FLAGS: &str = "--runMode alignReads --readFilesType SAM PE \
    --readFilesIn /dev/stdin --readNameSeparator none --outStd SAM --outSAMtype SAM \
    --outSAMorder PairedKeepInputOrder --outSAMunmapped Within";

/// A STAR-style genome directory for [`create_test_reference`]'s one contig, and
/// a directory holding the stand-in `STAR` executable.
struct Fixture {
    tmp: TempDir,
    reference: PathBuf,
    index: PathBuf,
    bin_dir: PathBuf,
}

impl Fixture {
    fn new() -> Self {
        use std::os::unix::fs::PermissionsExt;
        let tmp = TempDir::new().expect("tempdir");
        let reference = create_test_reference(tmp.path());
        let index = tmp.path().join("genome.star");
        std::fs::create_dir(&index).unwrap();
        for file in ["Genome", "SA", "SAindex", "genomeParameters.txt", "chrStart.txt"] {
            std::fs::write(index.join(file), b"").unwrap();
        }
        std::fs::write(index.join("chrNameLength.txt"), "chr1\t10000\n").unwrap();
        std::fs::write(index.join("chrName.txt"), "chr1\n").unwrap();
        std::fs::write(index.join("chrLength.txt"), "10000\n").unwrap();
        let bin_dir = tmp.path().join("bin");
        std::fs::create_dir(&bin_dir).unwrap();
        let star = bin_dir.join("STAR");
        std::fs::write(&star, FAKE_STAR).unwrap();
        std::fs::set_permissions(&star, std::fs::Permissions::from_mode(0o755)).unwrap();
        Self { tmp, reference, index, bin_dir }
    }

    fn path(&self) -> &Path {
        self.tmp.path()
    }

    /// `PATH` with the stand-in `STAR` first.
    fn path_env(&self) -> String {
        format!("{}:{}", self.bin_dir.display(), std::env::var("PATH").unwrap_or_default())
    }

    /// Run `runall --start-from align --stop-after zipper` on `input` with the
    /// given aligner arguments, from the working directory `cwd`.
    fn runall(&self, input: &Path, output: &Path, cwd: &Path, aligner: &[&str]) -> Command {
        let mut command = Command::new(env!("CARGO_BIN_EXE_fgumi"));
        command
            .args(["runall", "--start-from", "align", "--stop-after", "zipper"])
            .args(["-i", p(input), "--ref", p(&self.reference), "-o", p(output)])
            .args(aligner)
            .env("PATH", self.path_env())
            .current_dir(cwd)
            .stdin(Stdio::null());
        command
    }
}

/// Writes `pairs` unmapped, UMI-tagged 20 bp read pairs, named in queryname
/// order.
fn write_unmapped_pairs(path: &Path, pairs: usize) {
    let header = create_minimal_header("chr1", 10_000);
    let records: Vec<RawRecord> = (0..pairs)
        .flat_map(|i| {
            let name = format!("q{i:06}");
            let rx: Vec<u8> = (0..8).rev().map(|shift| b"ACGT"[(i >> (2 * shift)) & 3]).collect();
            [
                (flags::FIRST_SEGMENT, b"ACGTACGTACGTACGTACGT", 30),
                (flags::LAST_SEGMENT, b"TTGGCCAAGTTTGGCCAAGT", 20),
            ]
            .map(|(segment, seq, qual)| {
                let mut b = SamBuilder::new();
                b.read_name(name.as_bytes())
                    .flags(flags::PAIRED | flags::UNMAPPED | flags::MATE_UNMAPPED | segment)
                    .sequence(seq)
                    .qualities(&[qual; 20]);
                b.add_string_tag(fgumi_lib::sam::SamTag::RX, &rx);
                b.build()
            })
        })
        .collect();
    write_bam(path, &header, &records);
}

fn assert_success(output: &Output, context: &str) {
    assert!(
        output.status.success(),
        "{context} failed:\n{}",
        String::from_utf8_lossy(&output.stderr)
    );
}

#[rstest::rstest]
#[case::first(0, "AAAAAAAA")]
#[case::second(1, "AAAAAAAC")]
#[case::third(2, "AAAAAAAG")]
fn unmapped_fixture_distinguishes_template_tags(#[case] pair: usize, #[case] expected: &str) {
    use noodles::sam::alignment::record::data::field::Tag;
    use noodles::sam::alignment::record_buf::data::field::Value;
    let tmp = TempDir::new().unwrap();
    let input = tmp.path().join("input.bam");
    write_unmapped_pairs(&input, 3);
    let (_, records) = read_bam_output(&input);
    for record in &records[pair * 2..pair * 2 + 2] {
        assert_eq!(record.data().get(&Tag::new(b'R', b'X')), Some(&Value::String(expected.into())));
    }
}

#[rstest::rstest]
#[case::unpaired(0)]
#[case::r1_only(flags::PAIRED | flags::FIRST_SEGMENT)]
#[case::r2_only(flags::PAIRED | flags::LAST_SEGMENT)]
fn star_rejects_an_incomplete_pair_by_name(#[case] mate_flags: u16) {
    let fx = Fixture::new();
    let input = fx.path().join("orphan.bam");
    let mut builder = SamBuilder::new();
    builder
        .read_name(b"orphan_template")
        .flags(mate_flags | flags::UNMAPPED)
        .sequence(b"ACGT")
        .qualities(&[30; 4]);
    write_bam(&input, &create_minimal_header("chr1", 10_000), &[builder.build()]);
    let output = fx
        .runall(
            &input,
            &fx.path().join("out.bam"),
            fx.path(),
            &["--aligner::preset", "star", "--aligner::index", p(&fx.index)],
        )
        .output()
        .expect("spawn fgumi");
    let stderr = String::from_utf8_lossy(&output.stderr);
    assert!(!output.status.success(), "{stderr}");
    assert!(
        stderr.contains("complete primary pair") && stderr.contains("orphan_template"),
        "{stderr}"
    );
}

/// The preset runs STAR with the streaming flags, extra arguments and thread
/// count, keeps STAR's side files in a scratch directory under
/// `FGUMI_TMP_DIRS` that it removes afterwards (not in the working directory),
/// merges STAR's mapped, spliced and half-mapped records with the unmapped
/// input, and gives record-for-record the output of the same STAR run through
/// `--aligner::command` with SAM input.
#[test]
fn star_preset_matches_the_command_route() {
    let fx = Fixture::new();
    let input = fx.path().join("unmapped.bam");
    write_unmapped_pairs(&input, 200);
    let cwd = fx.path().join("cwd");
    std::fs::create_dir(&cwd).unwrap();
    let args_log = fx.path().join("star-args.txt");
    let tmp_root = fx.path().join("fgumi-tmp");
    std::fs::create_dir(&tmp_root).unwrap();

    let preset_out = fx.path().join("preset.bam");
    let output = fx
        .runall(
            &input,
            &preset_out,
            &cwd,
            &[
                "--aligner::preset",
                "star",
                "--aligner::index",
                p(&fx.index),
                "--aligner::threads",
                "3",
                "--aligner::extra-args",
                "--outSAMmapqUnique 60 --outSAMattributes NH HI",
            ],
        )
        .env("FAKE_STAR_ARGS_LOG", &args_log)
        .env("FGUMI_TMP_DIRS", &tmp_root)
        .output()
        .expect("spawn fgumi");
    assert_success(&output, "runall --aligner::preset star");

    let star_args = std::fs::read_to_string(&args_log).expect("STAR ran");
    for needle in [
        "--runThreadN 3".to_string(),
        "--outSAMmapqUnique 60 --outSAMattributes NH HI".to_string(),
        format!("--outFileNamePrefix {}/fgumi-star.", tmp_root.display()),
    ] {
        assert!(star_args.contains(&needle), "STAR args lack {needle:?}: {star_args}");
    }
    let leftovers: Vec<_> = std::fs::read_dir(&cwd).unwrap().collect();
    assert!(leftovers.is_empty(), "STAR's files must not land in the working directory");
    let scratch: Vec<_> = std::fs::read_dir(&tmp_root).unwrap().collect();
    assert!(scratch.is_empty(), "the STAR scratch directory must be removed: {scratch:?}");

    let command_out = fx.path().join("command.bam");
    let command_prefix = fx.path().join("command-star");
    std::fs::create_dir(&command_prefix).unwrap();
    let command = format!(
        "STAR {STAR_STREAMING_FLAGS} --runThreadN 3 --genomeDir {} --outTmpDir {} \
         --outFileNamePrefix {}/ # {{ref}}",
        fx.index.display(),
        fx.path().join("command-tmp").display(),
        command_prefix.display()
    );
    let output = fx
        .runall(
            &input,
            &command_out,
            &cwd,
            &["--aligner::command", &command, "--aligner::input-format", "sam"],
        )
        .output()
        .expect("spawn fgumi");
    assert_success(&output, "runall --aligner::command STAR");

    let (preset_header, preset_records) = read_bam_output(&preset_out);
    let (command_header, command_records) = read_bam_output(&command_out);
    let (_, input_records) = read_bam_output(&input);
    assert_star_records_merged(&preset_records, &input_records);
    assert_star_records_merged(&command_records, &input_records);
    assert_eq!(preset_records, command_records);
    assert_eq!(preset_header.reference_sequences(), command_header.reference_sequences());
}

/// The SAM character of the CIGAR operations the stand-in writes.
fn cigar_op_char(kind: noodles::sam::alignment::record::cigar::op::Kind) -> char {
    use noodles::sam::alignment::record::cigar::op::Kind;
    match kind {
        Kind::Match => 'M',
        Kind::Skip => 'N',
        other => panic!("unexpected CIGAR operation {other:?}"),
    }
}

/// Checks STAR's first two pairs' alignments and every input record's ordered
/// identity, mate, SEQ, QUAL and `RX`, including the unmapped tail.
fn assert_star_records_merged(
    records: &[noodles::sam::alignment::RecordBuf],
    input_records: &[noodles::sam::alignment::RecordBuf],
) {
    use noodles::sam::alignment::record::data::field::Tag;

    assert_eq!(records.len(), input_records.len(), "every input read must reach the output");
    for (i, (record, input)) in records.iter().zip(input_records).enumerate() {
        assert_eq!(record.name(), input.name(), "record {i} QNAME");
        assert_eq!(
            record.flags().is_first_segment(),
            input.flags().is_first_segment(),
            "record {i} first segment"
        );
        assert_eq!(
            record.flags().is_last_segment(),
            input.flags().is_last_segment(),
            "record {i} last segment"
        );
        assert_eq!(record.sequence(), input.sequence(), "record {i} SEQ");
        assert_eq!(record.quality_scores(), input.quality_scores(), "record {i} QUAL");
        assert_eq!(
            record.data().get(&Tag::new(b'R', b'X')),
            input.data().get(&Tag::new(b'R', b'X')),
            "record {i} RX"
        );
        if i >= 4 {
            assert_eq!(record.flags(), input.flags(), "record {i} unmapped flags");
        }
    }

    // (name, first segment?, flags, 1-based start, CIGAR, MAPQ) for the reads
    // the stand-in mapped, as it wrote them. STAR's unique-read MAPQ of 255 is
    // SAM's "unavailable", which noodles reads as `None`.
    let expected = [
        ("q000000", true, 65, Some(100), "10M100N10M", None),
        ("q000000", false, 129, Some(400), "20M", None),
        ("q000001", true, 73, Some(200), "20M", None),
        ("q000001", false, 133, Some(200), "", Some(0)),
    ];
    for (name, first, flags, start, cigar, mapq) in expected {
        let record = records
            .iter()
            .find(|r| {
                r.name().map(AsRef::as_ref) == Some(name.as_bytes())
                    && r.flags().is_first_segment() == first
                    && !r.flags().is_secondary()
                    && !r.flags().is_supplementary()
            })
            .unwrap_or_else(|| panic!("{name} (first={first}) missing from the output"));
        assert_eq!(record.flags().bits(), flags, "{name} first={first} flags");
        assert_eq!(record.alignment_start().map(usize::from), start, "{name} first={first} start");
        let rendered = record.cigar().as_ref().iter().fold(String::new(), |mut out, op| {
            use std::fmt::Write as _;
            write!(out, "{}{}", op.len(), cigar_op_char(op.kind())).expect("write to a String");
            out
        });
        assert_eq!(rendered, cigar, "{name} first={first} CIGAR");
        assert_eq!(record.mapping_quality().map(u8::from), mapq, "{name} first={first} MAPQ");
    }
}

/// The stand-in really checks its flags: a command without
/// `--outSAMorder PairedKeepInputOrder` fails rather than passing vacuously.
#[test]
fn fake_star_rejects_a_missing_streaming_flag() {
    let fx = Fixture::new();
    let input = fx.path().join("unmapped.bam");
    write_unmapped_pairs(&input, 10);
    let command = format!(
        "STAR --runMode alignReads --readFilesType SAM PE --readFilesIn /dev/stdin \
         --readNameSeparator none --outStd SAM \
         --outSAMtype SAM --outSAMunmapped Within --runThreadN 1 --genomeDir {} \
         --outTmpDir {} --outFileNamePrefix {}/ # {{ref}}",
        fx.index.display(),
        fx.path().join("tmp").display(),
        fx.path().display()
    );
    let output = fx
        .runall(
            &input,
            &fx.path().join("out.bam"),
            fx.path(),
            &["--aligner::command", &command, "--aligner::input-format", "sam"],
        )
        .output()
        .expect("spawn fgumi");
    assert!(!output.status.success(), "the stand-in must reject a missing flag");
    assert!(String::from_utf8_lossy(&output.stderr).contains("missing flag: order"));
}

/// An index naming other sequences than the `--ref` dictionary fails before
/// STAR is started.
#[test]
fn star_preset_rejects_a_mismatched_index_before_starting_star() {
    let fx = Fixture::new();
    let input = fx.path().join("unmapped.bam");
    write_unmapped_pairs(&input, 10);
    std::fs::write(fx.index.join("chrNameLength.txt"), "chr1\t10000\nchrM\t16569\n").unwrap();
    let args_log = fx.path().join("star-args.txt");
    let output = fx
        .runall(
            &input,
            &fx.path().join("out.bam"),
            fx.path(),
            &["--aligner::preset", "star", "--aligner::index", p(&fx.index)],
        )
        .env("FAKE_STAR_ARGS_LOG", &args_log)
        .output()
        .expect("spawn fgumi");
    assert!(!output.status.success());
    let stderr = String::from_utf8_lossy(&output.stderr);
    assert!(stderr.contains("does not match the sequence dictionary"), "{stderr}");
    assert!(!args_log.exists(), "STAR must not be started");
}

/// Seconds the stand-in spends "loading the genome" in the slow-start tests:
/// longer than the watchdog's fatal window at `--deadlock-timeout 1` (6 s).
const SLOW_LOAD_SECS: &str = "9";

/// Template pairs for the slow-start tests: enough that the input backs up
/// into the pipeline's queues while STAR loads, so the watchdog sees work in
/// flight.
const SLOW_START_PAIRS: usize = 50_000;

/// A STAR genome load longer than the deadlock watchdog's fatal window does
/// not fail a star-preset run: the wait for STAR's first output is an
/// external wait. Through `--aligner::command` the same silence is a wedge.
#[rstest::rstest]
#[case::preset(true)]
#[case::command(false)]
fn slow_star_start_is_not_a_deadlock_for_the_preset(#[case] preset: bool) {
    let fx = Fixture::new();
    let input = fx.path().join("unmapped.bam");
    write_unmapped_pairs(&input, SLOW_START_PAIRS);
    let command = format!(
        "STAR {STAR_STREAMING_FLAGS} --runThreadN 1 --genomeDir {} --outTmpDir {} \
         --outFileNamePrefix {}/ # {{ref}}",
        fx.index.display(),
        fx.path().join("command-tmp").display(),
        fx.path().display()
    );
    let aligner: Vec<&str> = if preset {
        vec!["--aligner::preset", "star", "--aligner::index", p(&fx.index)]
    } else {
        vec!["--aligner::command", &command, "--aligner::input-format", "sam"]
    };
    let output = fx
        .runall(&input, &fx.path().join("out.bam"), fx.path(), &aligner)
        .args(["--deadlock-timeout", "1"])
        .env("FAKE_STAR_LOAD_SECS", SLOW_LOAD_SECS)
        .output()
        .expect("spawn fgumi");
    let stderr = String::from_utf8_lossy(&output.stderr);
    if preset {
        assert_success(&output, "runall --aligner::preset star with a slow genome load");
        let (_, records) = read_bam_output(&fx.path().join("out.bam"));
        let (_, input_records) = read_bam_output(&input);
        assert_star_records_merged(&records, &input_records);
    } else {
        assert!(!output.status.success(), "a silent command aligner is a wedge:\n{stderr}");
        assert!(stderr.contains("deadlock"), "{stderr}");
    }
}

/// A configured startup deadline bounds even a STAR that neither reads stdin
/// nor writes stdout, while a sufficiently generous deadline allows the load.
#[rstest::rstest]
#[case::expires("1", false)]
#[case::allows_load("5", true)]
fn star_startup_timeout_is_configurable(#[case] timeout: &str, #[case] succeeds: bool) {
    let fx = Fixture::new();
    let input = fx.path().join("unmapped.bam");
    write_unmapped_pairs(&input, 10);
    let output_path = fx.path().join("out.bam");
    let tmp_root = fx.path().join("scratch");
    let output = fx
        .runall(
            &input,
            &output_path,
            fx.path(),
            &[
                "--aligner::preset",
                "star",
                "--aligner::index",
                p(&fx.index),
                "--aligner::startup-timeout",
                timeout,
            ],
        )
        .env("FAKE_STAR_LOAD_SECS", "2")
        .env("FGUMI_TMP_DIRS", &tmp_root)
        .output()
        .expect("spawn fgumi");
    let stderr = String::from_utf8_lossy(&output.stderr);
    assert_eq!(output.status.success(), succeeds, "{stderr}");
    if succeeds {
        let (_, records) = read_bam_output(&output_path);
        let (_, input_records) = read_bam_output(&input);
        assert_star_records_merged(&records, &input_records);
    } else {
        assert!(stderr.contains("aligner startup timed out"), "{stderr}");
        assert!(stderr.contains("--aligner::startup-timeout"), "{stderr}");
    }
    assert_eq!(std::fs::read_dir(tmp_root).unwrap().count(), 0, "scratch must be cleaned");
}

#[rstest::rstest]
#[case::zero("0")]
#[case::negative("-1")]
#[case::fraction("1.5")]
#[case::overflow("18446744073709551616")]
fn star_startup_timeout_rejects_invalid_seconds(#[case] seconds: &str) {
    let fx = Fixture::new();
    let input = fx.path().join("unmapped.bam");
    write_unmapped_pairs(&input, 10);
    let output = fx
        .runall(
            &input,
            &fx.path().join("out.bam"),
            fx.path(),
            &["--aligner::preset", "star", "--aligner::index", p(&fx.index)],
        )
        .arg(format!("--aligner::startup-timeout={seconds}"))
        .output()
        .expect("spawn fgumi");
    let stderr = String::from_utf8_lossy(&output.stderr);
    assert_eq!(output.status.code(), Some(2), "{stderr}");
    assert!(stderr.contains("invalid value"), "{stderr}");
    assert!(stderr.contains("--aligner::startup-timeout"), "{stderr}");
}

/// `--aligner::preset star` is rejected for realignment, before STAR is
/// started, by standalone `fgumi realign` and by `runall --realign`.
#[rstest::rstest]
#[case::standalone_realign(true, true)]
#[case::runall_realign(false, true)]
#[case::runall_without_index(false, false)]
fn star_preset_cannot_realign(#[case] standalone: bool, #[case] with_index: bool) {
    let fx = Fixture::new();
    let input = fx.path().join("unmapped.bam");
    write_unmapped_pairs(&input, 10);
    let output_bam = fx.path().join("out.bam");
    let args_log = fx.path().join("star-args.txt");
    let mut star = vec!["--aligner::preset", "star"];
    if with_index {
        star.extend(["--aligner::index", p(&fx.index)]);
    }
    let mut command = if standalone {
        let mut command = Command::new(env!("CARGO_BIN_EXE_fgumi"));
        command
            .args(["realign", "-i", p(&input), "-o", p(&output_bam), "-r", p(&fx.reference)])
            .args(star)
            .env("PATH", fx.path_env())
            .stdin(Stdio::null());
        command
    } else {
        let mut command = Command::new(env!("CARGO_BIN_EXE_fgumi"));
        command
            .args(["runall", "--start-from", "align", "--stop-after", "consensus"])
            .args(["--consensus", "simplex", "--realign"])
            .args(["-i", p(&input), "--ref", p(&fx.reference), "-o", p(&output_bam)])
            .args(star)
            .env("PATH", fx.path_env())
            .stdin(Stdio::null());
        command
    };
    let output = command.env("FAKE_STAR_ARGS_LOG", &args_log).output().expect("spawn fgumi");
    let stderr = String::from_utf8_lossy(&output.stderr);
    assert!(!output.status.success(), "{stderr}");
    assert!(stderr.contains("--aligner::preset star cannot realign"), "{stderr}");
    assert!(!args_log.exists(), "STAR must not be started");
}

/// The start-up wait the preset excuses ends at STAR's first output: a STAR
/// that writes its header and then goes silent past the watchdog's fatal window
/// with input pending is reported as a deadlock, as through
/// `--aligner::command`.
#[test]
fn silence_after_the_star_header_is_a_deadlock_for_the_preset() {
    let fx = Fixture::new();
    let input = fx.path().join("unmapped.bam");
    write_unmapped_pairs(&input, SLOW_START_PAIRS);
    let output = fx
        .runall(
            &input,
            &fx.path().join("out.bam"),
            fx.path(),
            &["--aligner::preset", "star", "--aligner::index", p(&fx.index)],
        )
        .args(["--deadlock-timeout", "1"])
        .env("FAKE_STAR_LOAD_SECS", SLOW_LOAD_SECS)
        .env("FAKE_STAR_EARLY_HEADER", "1")
        .output()
        .expect("spawn fgumi");
    let stderr = String::from_utf8_lossy(&output.stderr);
    assert!(!output.status.success(), "a silent STAR after its header is a wedge:\n{stderr}");
    assert!(stderr.contains("deadlock"), "{stderr}");
}
