//! Parity gate for the `extract` command's single-threaded-path retirement (C4):
//! `Extract::execute` no longer has a serial in-process loop reached when
//! `--threads` is absent — it *always* routes through the declarative chain
//! builder. This test proves that cutover lost nothing user-observable.
//!
//! Two independent things are checked here:
//!
//! 1. **The cutover actually happened** (`extract_no_threads_routes_through_chain`).
//!    A no-`--threads` run now emits the chain's `"Starting Extract"` banner (and
//!    the `Input:`/`Output:` lines) that `ChainBuilder::add_extract` logs and the
//!    retired serial tail never did. This is the genuine RED/GREEN discriminator:
//!    before the removal a no-`--threads` run took the serial path and logged only
//!    the `Extracting UMIs` timer + `Processed records` progress; after it, the
//!    chain logs the banner.
//!
//! 2. **Output parity with the pre-removal serial path**
//!    (`cutover_matches_baseline`). The current build's `extract` output — records
//!    (byte-identical, modulo the `@PG` line) — must match the frozen serial
//!    baseline binary. The baseline path comes from `FGUMI_BASELINE_BIN`; when it
//!    is unset (or names a missing file) the case degrades to a self-consistency
//!    oracle (every record's template sequence and combined `RX`/`RG` tags
//!    asserted directly) rather than skipping — the exact fallback discipline of
//!    `test_sort_cutover_parity.rs`. Cases
//!    cover the three extract input shapes the chain must support: two-file input,
//!    interleaved input (`-p/--interleaved`), and BGZF input with `--check-crc`,
//!    each over a read structure carrying a UMI segment (`5M+T`).
//!
//! **Not a RED/GREEN gate for the parity half.** Because the removed serial loop
//! and the chain were already output-equivalent, the baseline byte-parity check
//! passes on both sides of the change — like the sort/retag cutovers it guards
//! equivalence, it does not observe a regression the cutover introduces. The
//! banner check in (1) is the part that flips RED→GREEN across the removal.

use std::ffi::OsStr;
use std::fs::File;
use std::io::Write;
use std::path::{Path, PathBuf};
use std::process::Command;

use noodles::sam::alignment::RecordBuf;
use noodles::sam::alignment::record::data::field::Tag;
use noodles::sam::alignment::record_buf::data::field::Value;
use noodles_bgzf::io::Writer as BgzfWriter;
use rstest::rstest;
use tempfile::TempDir;

use crate::helpers::read_bam_output;
use fgumi_lib::sam::SamTag;

/// Resolves the saved pre-removal serial baseline binary to compare against.
///
/// The path comes solely from `FGUMI_BASELINE_BIN`; when unset (or naming a
/// missing file) this returns `None` — "no baseline oracle available", never a
/// silent pass. Callers layer the baseline byte-parity check on top of the
/// always-available self-consistency oracle; a missing baseline drops only the
/// byte-parity half, it never skips the case. No hardcoded fallback: a baseline
/// binary is host-specific and must never be a path committed into the repo.
fn baseline_bin() -> Option<PathBuf> {
    let path = PathBuf::from(std::env::var_os("FGUMI_BASELINE_BIN")?);
    if path.is_file() {
        return Some(path);
    }
    eprintln!(
        "FGUMI_BASELINE_BIN={} does not name an existing file; baseline oracle unavailable",
        path.display()
    );
    None
}

/// Removes every `@PG` line from a SAM header text blob.
///
/// Both the current build and the baseline binary stamp a single `@PG` line
/// whose `VN` (git-describe version) and `CL` (command line, naming argv[0] — a
/// different binary path for each side — and the per-run output path) necessarily
/// differ between two independent invocations. Stripping the whole line is
/// correct here because neither side emits more than one `@PG` for these inputs
/// (extract synthesizes a fresh header with exactly one `@PG`).
fn strip_pg_lines(text: &str) -> String {
    if text.is_empty() {
        return String::new();
    }
    let mut lines: Vec<&str> = text.split('\n').collect();
    let had_trailing_newline = lines.last() == Some(&"");
    if had_trailing_newline {
        lines.pop();
    }
    lines.retain(|line| !line.starts_with("@PG"));
    let mut out = lines.join("\n");
    if had_trailing_newline {
        out.push('\n');
    }
    out
}

/// Reads a BAM's raw BGZF stream, decompresses it, strips the `@PG` line(s) from
/// the embedded SAM header text, and returns the resulting bytes (new header +
/// unmodified `n_ref`/reference-list/record bytes).
///
/// Comparing the *decompressed* BAM binary — rather than `RecordBuf`-parsed
/// records or raw file bytes — keeps everything except the `@PG` line an exact,
/// uninterpreted byte comparison: BGZF block boundaries differ between the two
/// writers even for identical logical content (so raw file bytes never match),
/// while parsing into `RecordBuf` and re-encoding could mask a real tag-order or
/// binary-layout regression by normalizing it away. extract's contract is the
/// exact UMI/quality/tag bytes it writes, so that masking risk is what must be
/// avoided.
fn decompressed_records_without_pg(path: &Path) -> Vec<u8> {
    let file = File::open(path).unwrap_or_else(|e| panic!("open {}: {e}", path.display()));
    let mut reader = noodles::bgzf::io::Reader::new(std::io::BufReader::new(file));
    let mut raw = Vec::new();
    std::io::Read::read_to_end(&mut reader, &mut raw)
        .unwrap_or_else(|e| panic!("decompress BGZF stream for {}: {e}", path.display()));

    assert!(
        raw.len() >= 8 && &raw[0..4] == b"BAM\x01",
        "{} does not decompress to a BAM binary stream (missing magic)",
        path.display()
    );
    let l_text = i32::from_le_bytes(raw[4..8].try_into().expect("4 bytes"));
    let l_text = usize::try_from(l_text).expect("l_text is non-negative");
    let text_start = 8;
    let text_end = text_start + l_text;
    assert!(raw.len() >= text_end, "{} header text runs past end of stream", path.display());

    let text = String::from_utf8_lossy(&raw[text_start..text_end]);
    let stripped_text = strip_pg_lines(&text);

    let mut out = Vec::with_capacity(raw.len());
    out.extend_from_slice(b"BAM\x01");
    let new_l_text = i32::try_from(stripped_text.len()).expect("stripped header text fits i32");
    out.extend_from_slice(&new_l_text.to_le_bytes());
    out.extend_from_slice(stripped_text.as_bytes());
    out.extend_from_slice(&raw[text_end..]); // n_ref, reference list, and all records, untouched
    out
}

/// A pair of `(name, sequence, quality)` FASTQ records: R1 carries UMI `AAAAA`
/// (template `CCCCC`), R2 carries UMI `GGGGG` (template `TTTTT`), so a `5M+T`
/// read structure yields template `CCCCC`/`TTTTT` and a combined `RX` of
/// `AAAAA-GGGGG` on both mates.
const R1_RECORDS: &[(&str, &str, &str)] =
    &[("pair0", "AAAAACCCCC", "IIIIIIIIII"), ("pair1", "TTTTTGGGGG", "IIIIIIIIII")];
const R2_RECORDS: &[(&str, &str, &str)] =
    &[("pair0", "GGGGGTTTTT", "IIIIIIIIII"), ("pair1", "CCCCCAAAAA", "IIIIIIIIII")];

/// Write plain-text FASTQ records to `dir/name`, returning the path.
fn write_plain_fastq(dir: &Path, name: &str, records: &[(&str, &str, &str)]) -> PathBuf {
    let path = dir.join(name);
    let mut file = File::create(&path).unwrap_or_else(|e| panic!("create {name}: {e}"));
    for (rname, seq, qual) in records {
        writeln!(file, "@{rname}\n{seq}\n+\n{qual}").expect("write fastq record");
    }
    path
}

/// Write BGZF-compressed FASTQ records to `dir/name`, returning the path. Used by
/// the `--check-crc` case, where the input must be BGZF for the CRC policy to bite.
fn write_bgzf_fastq(dir: &Path, name: &str, records: &[(&str, &str, &str)]) -> PathBuf {
    let path = dir.join(name);
    let file = File::create(&path).unwrap_or_else(|e| panic!("create {name}: {e}"));
    let mut writer = BgzfWriter::new(file);
    for (rname, seq, qual) in records {
        writeln!(writer, "@{rname}\n{seq}\n+\n{qual}").expect("write bgzf fastq record");
    }
    writer.finish().expect("finish bgzf fastq");
    path
}

/// The extract input shape a case exercises. Each is a distinct chain-source
/// topology (`SourceSpec::Fastqs` vs `SourceSpec::InterleavedFastq`) and/or a
/// distinct reader-open policy (`--check-crc`), so all three must retain parity.
#[derive(Clone, Copy, Debug)]
enum InputShape {
    /// Two plain FASTQ files (R1 + R2).
    TwoFile,
    /// One interleaved plain FASTQ file (`-p/--interleaved`).
    Interleaved,
    /// Two BGZF FASTQ files run with `--check-crc`.
    CheckCrcBgzf,
}

/// Materializes the input FASTQ(s) for `shape` under `dir` and returns
/// `(input_paths, interleaved, check_crc)` ready to pass to `run_extract`.
fn build_inputs(shape: InputShape, dir: &Path) -> (Vec<PathBuf>, bool, bool) {
    match shape {
        InputShape::TwoFile => {
            let r1 = write_plain_fastq(dir, "r1.fq", R1_RECORDS);
            let r2 = write_plain_fastq(dir, "r2.fq", R2_RECORDS);
            (vec![r1, r2], false, false)
        }
        InputShape::Interleaved => {
            // R1/R2 records alternate in one physical stream.
            let interleaved: Vec<(&str, &str, &str)> =
                R1_RECORDS.iter().zip(R2_RECORDS.iter()).flat_map(|(a, b)| [*a, *b]).collect();
            let path = write_plain_fastq(dir, "interleaved.fq", &interleaved);
            (vec![path], true, false)
        }
        InputShape::CheckCrcBgzf => {
            let r1 = write_bgzf_fastq(dir, "r1.fq.gz", R1_RECORDS);
            let r2 = write_bgzf_fastq(dir, "r2.fq.gz", R2_RECORDS);
            (vec![r1, r2], false, true)
        }
    }
}

/// Runs `<bin> extract --inputs <inputs...> --output <output> --read-structures
/// 5M+T [5M+T] --sample s --library l --compression-level 1 [--interleaved]
/// [--check-crc]` (no `--threads`, so the current build takes the post-cutover
/// chain path and the baseline takes its serial path) with `RUST_LOG=info`, and
/// returns the process output for stderr assertions.
fn run_extract(
    bin: &Path,
    inputs: &[PathBuf],
    output: &Path,
    interleaved: bool,
    check_crc: bool,
) -> std::process::Output {
    let mut cmd = Command::new(bin);
    cmd.env("RUST_LOG", "info");
    cmd.arg("extract").arg("--inputs");
    cmd.args(inputs.iter().map(|p| p.as_os_str()));
    cmd.args([OsStr::new("--output"), output.as_os_str(), OsStr::new("--read-structures")]);
    // One `5M+T` per read: two structures for a two-read template (two-file and
    // interleaved both describe two reads), matching validation's requirement.
    if interleaved || inputs.len() == 2 {
        cmd.args([OsStr::new("5M+T"), OsStr::new("5M+T")]);
    } else {
        cmd.arg("5M+T");
    }
    cmd.args([
        OsStr::new("--sample"),
        OsStr::new("s"),
        OsStr::new("--library"),
        OsStr::new("l"),
        OsStr::new("--compression-level"),
        OsStr::new("1"),
    ]);
    if interleaved {
        cmd.arg("--interleaved");
    }
    if check_crc {
        cmd.arg("--check-crc");
    }
    cmd.output().unwrap_or_else(|e| panic!("failed to spawn `{}` extract: {e}", bin.display()))
}

/// A no-`--threads` run now routes through the declarative chain, which logs the
/// `"Starting Extract"` banner (and `Input:`/`Output:` lines) from
/// `ChainBuilder::add_extract`. The retired serial tail logged no such line — it
/// emitted only the `Extracting UMIs` timer + `Processed records` progress — so
/// this is the RED (pre-removal) → GREEN (post-removal) discriminator.
#[test]
fn extract_no_threads_routes_through_chain() {
    let dir = TempDir::new().expect("temp dir");
    let (inputs, interleaved, check_crc) = build_inputs(InputShape::TwoFile, dir.path());
    let output = dir.path().join("out.bam");

    let current_bin = Path::new(env!("CARGO_BIN_EXE_fgumi"));
    let out = run_extract(current_bin, &inputs, &output, interleaved, check_crc);
    let stderr = String::from_utf8_lossy(&out.stderr);
    assert!(out.status.success(), "no-threads extract must succeed; stderr:\n{stderr}");
    assert!(
        stderr.contains("Starting Extract"),
        "a no-`--threads` extract must route through the chain (which logs the \
         `Starting Extract` banner); the serial path is retired. stderr:\n{stderr}"
    );
}

/// Output parity of the post-cutover chain against the pre-removal serial
/// baseline binary, across the three extract input shapes — plus the
/// always-available self-consistency oracle when no baseline is set.
#[rstest]
#[case::two_file(InputShape::TwoFile)]
#[case::interleaved(InputShape::Interleaved)]
#[case::check_crc_bgzf(InputShape::CheckCrcBgzf)]
fn cutover_matches_baseline(#[case] shape: InputShape) {
    let dir = TempDir::new().expect("temp dir");
    let (inputs, interleaved, check_crc) = build_inputs(shape, dir.path());

    let current_bin = Path::new(env!("CARGO_BIN_EXE_fgumi"));
    let current_out = dir.path().join("current.bam");
    let current = run_extract(current_bin, &inputs, &current_out, interleaved, check_crc);
    let current_stderr = String::from_utf8_lossy(&current.stderr);
    assert!(
        current.status.success(),
        "current extract ({shape:?}) must succeed; stderr:\n{current_stderr}"
    );

    if let Some(baseline) = baseline_bin() {
        let baseline_out = dir.path().join("baseline.bam");
        let base = run_extract(&baseline, &inputs, &baseline_out, interleaved, check_crc);
        assert!(
            base.status.success(),
            "baseline extract ({shape:?}) failed; stderr:\n{}",
            String::from_utf8_lossy(&base.stderr)
        );

        assert_eq!(
            decompressed_records_without_pg(&current_out),
            decompressed_records_without_pg(&baseline_out),
            "chain extract output ({shape:?}) diverges from the pre-removal serial baseline \
             binary ({}) after stripping @PG — a real cutover parity bug, not something to relax",
            baseline.display(),
        );
    } else {
        eprintln!(
            "SKIP baseline half of cutover_matches_baseline[{shape:?}]: FGUMI_BASELINE_BIN is \
             unset or does not name an existing file — running self-consistency oracle instead"
        );
        assert_self_consistent(&current_out);
    }
}

/// Always-available oracle used when no baseline binary is set (the default in
/// CI): the chain's output must carry exactly the records the known input +
/// `5M+T` read structure prescribe — for *every* input template pair, two mates
/// carrying the template bases after the 5 bp UMI, both tagged with the combined
/// R1-R2 UMI in `RX` and the default `RG`. Every produced record is checked
/// (record identity, not just count), so a per-record divergence in any pair
/// fails the case.
fn assert_self_consistent(output: &Path) {
    let (_, recs) = read_bam_output(output);
    assert_eq!(recs.len(), R1_RECORDS.len() * 2, "two BAM records per input template pair");

    let string_tag = |rec: &RecordBuf, tag: SamTag| -> Vec<u8> {
        match rec.data().get(&Tag::from(tag)) {
            Some(Value::String(s)) => s.to_vec(),
            other => panic!("expected a string {tag:?} tag, got {other:?}"),
        }
    };

    // Derive the expected output from the known input `5M+T` structure so the
    // oracle reads as a spec: each read splits into a 5 bp UMI and the template
    // remainder; the two mates of a pair share `RX = <R1 UMI>-<R2 UMI>`.
    for (pair_idx, (r1, r2)) in R1_RECORDS.iter().zip(R2_RECORDS.iter()).enumerate() {
        let (r1_umi, r1_template) = r1.1.split_at(5);
        let (r2_umi, r2_template) = r2.1.split_at(5);
        let expected_rx = format!("{r1_umi}-{r2_umi}").into_bytes();

        let r1_out = &recs[pair_idx * 2];
        let r2_out = &recs[pair_idx * 2 + 1];
        assert_eq!(
            r1_out.sequence().as_ref(),
            r1_template.as_bytes(),
            "pair {pair_idx}: R1 template after 5M UMI"
        );
        assert_eq!(
            r2_out.sequence().as_ref(),
            r2_template.as_bytes(),
            "pair {pair_idx}: R2 template after 5M UMI"
        );
        for rec in [r1_out, r2_out] {
            assert_eq!(
                string_tag(rec, SamTag::RX),
                expected_rx,
                "pair {pair_idx}: combined R1-R2 UMI in RX"
            );
            assert_eq!(string_tag(rec, SamTag::RG), b"A", "pair {pair_idx}: default read group id");
        }
    }
}
