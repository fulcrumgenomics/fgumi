//! Parity gate for the consensus trio's single-threaded-path retirement (C4):
//! `Simplex::execute`, `Duplex::execute`, and `Codec::execute` no longer have a
//! serial in-process loop reached when `--threads` is absent — they *always*
//! route through the declarative chain builder, with or without `--threads`.
//! This test proves that cutover lost nothing user-observable.
//!
//! Two independent things are checked here:
//!
//! 1. **The cutover actually happened** (`no_threads_routes_through_chain`).
//!    A no-`--threads` run now emits the chain-only `"Total MI groups processed"`
//!    summary line that the chain finalize hooks log; the retired serial tail
//!    instead logged `"Total records processed"`. This is the genuine RED/GREEN
//!    discriminator: before the removal a no-`--threads` run took the serial fast
//!    path and printed the records-processed line; after it, the chain prints the
//!    MI-groups line. (The `"Using pipeline with N threads"` banner is emitted by
//!    the simplex/duplex chains but not the codec chain, so it is not a uniform
//!    discriminator across the trio; the MI-groups summary line is.)
//!
//! 2. **Output parity with the pre-removal serial path** (`cutover_matches_baseline`).
//!    The current build's consensus output — records (byte-identical, modulo the
//!    `@PG` line) and the `--stats` TSV — must match the frozen pre-removal
//!    serial baseline binary. The baseline path comes from `FGUMI_BASELINE_BIN`;
//!    when it is unset (or names a missing file) the case degrades to a
//!    self-consistency oracle (non-empty consensus output + a well-formed stats
//!    TSV asserted directly) rather than skipping — the exact fallback discipline
//!    of `test_sort_cutover_parity.rs`.
//!
//! **Not a RED/GREEN gate for the parity half.** Because the removed serial loop
//! and the chain were already output-equivalent (proven in-process by the
//! `test_*_chain_matches_single_threaded` suites), the baseline byte-parity check
//! passes on both sides of the change — it guards equivalence, it does not observe
//! a regression the cutover introduces. The chain-banner check in (1) is the part
//! that flips RED→GREEN across the removal.
//!
//! Only compiled with the `consensus` feature (gated at the `mod` declaration
//! in `main.rs`) — `simplex`/`duplex`/`codec` (and the chain machinery they
//! depend on) only exist in a `consensus` build, so the cutover this test
//! guards has no other configuration to run under.

use std::ffi::OsStr;
use std::path::{Path, PathBuf};
use std::process::Command;

use fgumi_raw_bam::RawRecord;
use rstest::rstest;
use tempfile::TempDir;

use crate::helpers::bam_generator::{create_minimal_header, create_umi_family};
use crate::helpers::read_bam_output;

/// The three consensus commands whose serial fast path this PR retires.
#[derive(Clone, Copy, Debug)]
enum ConsensusCmd {
    Simplex,
    Duplex,
    Codec,
}

impl ConsensusCmd {
    /// The subcommand name passed to the binary.
    fn name(self) -> &'static str {
        match self {
            ConsensusCmd::Simplex => "simplex",
            ConsensusCmd::Duplex => "duplex",
            ConsensusCmd::Codec => "codec",
        }
    }

    /// The `Starting …` banner the chain builder logs for this command, used to
    /// assert the chain path ran (not just that some pipeline banner appeared).
    fn starting_banner(self) -> &'static str {
        match self {
            ConsensusCmd::Simplex => "Starting Simplex",
            ConsensusCmd::Duplex => "Starting Duplex",
            ConsensusCmd::Codec => "Starting CODEC consensus calling",
        }
    }
}

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
/// (the test BAMs carry no pre-existing `@PG`).
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
/// binary-layout regression by normalizing it away. Consensus output carries
/// consensus quality/per-base tags whose exact byte layout is the contract, so
/// that masking risk is exactly what must be avoided.
fn decompressed_records_without_pg(path: &Path) -> Vec<u8> {
    let file = std::fs::File::open(path).unwrap_or_else(|e| panic!("open {}: {e}", path.display()));
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

/// Writes an MI-grouped input BAM for `simplex`: several deep single-strand UMI
/// families, each tagged with a distinct `MI`, deep enough to survive
/// `--min-reads 2`.
fn write_simplex_input(dir: &Path) -> PathBuf {
    use noodles::sam::alignment::io::Write as _;
    use noodles::sam::alignment::record::data::field::Tag;
    use noodles::sam::alignment::record_buf::data::field::Value;

    let input = dir.join("in.bam");
    let header = create_minimal_header("chr1", 10_000);
    let mut writer =
        noodles::bam::io::Writer::new(std::fs::File::create(&input).expect("create simplex input"));
    writer.write_header(&header).expect("write header");

    let families = [
        ("1", create_umi_family("ACGT", 5, "fam1", "ACGTACGTAC", 30)),
        ("2", create_umi_family("TGCA", 4, "fam2", "TTTTAAAAGG", 30)),
        ("3", create_umi_family("GGCC", 3, "fam3", "CCCCGGGGTT", 30)),
    ];
    let mi_tag = Tag::from(fgumi_lib::sam::SamTag::MI);
    for (mi, records) in families {
        for raw in &records {
            let mut record = crate::helpers::bam_generator::to_record_buf(raw);
            record.data_mut().insert(mi_tag, Value::from(mi));
            writer.write_alignment_record(&header, &record).expect("write record");
        }
    }
    writer.try_finish().expect("finish simplex input");
    input
}

/// Writes an AB/BA duplex input BAM for `duplex`: two double-stranded molecules,
/// each with matched top/bottom strands deep enough for a duplex call.
fn write_duplex_input(dir: &Path) -> PathBuf {
    let input = dir.join("in.bam");
    let molecules = vec![
        crate::test_duplex_command::create_duplex_molecule("1", "ACGTACGT", 30, 100, 4),
        crate::test_duplex_command::create_duplex_molecule("2", "TTGGCCAA", 30, 400, 3),
    ];
    crate::test_duplex_command::create_duplex_bam(&input, molecules);
    input
}

/// Writes a CODEC input BAM for `codec`: one molecule of several overlapping
/// read pairs sharing an `MI`, deep enough for a CODEC consensus call.
fn write_codec_input(dir: &Path) -> PathBuf {
    let input = dir.join("in.bam");
    let mut pairs: Vec<(RawRecord, RawRecord)> = Vec::new();
    for i in 0..4 {
        pairs.push(crate::test_codec_command::create_codec_read_pair(
            &format!("read{i}"),
            b"ACGTACGT",
            b"ACGTACGT",
            &[30; 8],
            &[30; 8],
            100,
            "UMI001",
            None,
        ));
    }
    crate::test_codec_command::create_codec_test_bam(&input, pairs);
    input
}

/// Writes the appropriate input BAM for `cmd` and returns its path.
fn write_input(cmd: ConsensusCmd, dir: &Path) -> PathBuf {
    match cmd {
        ConsensusCmd::Simplex => write_simplex_input(dir),
        ConsensusCmd::Duplex => write_duplex_input(dir),
        ConsensusCmd::Codec => write_codec_input(dir),
    }
}

/// Runs `<bin> <cmd> -i <input> -o <output> [--stats <stats>] --min-reads 2
/// --compression-level 1` (no `--threads`, so the current build takes the
/// post-cutover chain path and the baseline takes its serial path) with
/// `RUST_LOG=info`, and returns the process output for stderr assertions.
fn run_consensus(
    bin: &Path,
    cmd: ConsensusCmd,
    input: &Path,
    output: &Path,
    stats: Option<&Path>,
) -> std::process::Output {
    let mut command = Command::new(bin);
    command.env("RUST_LOG", "info").args([
        OsStr::new(cmd.name()),
        OsStr::new("-i"),
        input.as_os_str(),
        OsStr::new("-o"),
        output.as_os_str(),
        OsStr::new("--min-reads"),
        OsStr::new("2"),
        OsStr::new("--compression-level"),
        OsStr::new("1"),
    ]);
    if let Some(stats) = stats {
        command.args([OsStr::new("--stats"), stats.as_os_str()]);
    }
    command
        .output()
        .unwrap_or_else(|e| panic!("failed to spawn `{}` {}: {e}", bin.display(), cmd.name()))
}

/// A no-`--threads` run now routes through the declarative chain, whose finalize
/// hooks log the `"Total MI groups processed"` summary line (and the chain
/// builder logs the `Starting …` banner). The retired serial tail logged
/// `"Total records processed"` instead, so the MI-groups line is the RED
/// (pre-removal) → GREEN (post-removal) discriminator for the cutover.
#[rstest]
#[case::simplex(ConsensusCmd::Simplex)]
#[case::duplex(ConsensusCmd::Duplex)]
#[case::codec(ConsensusCmd::Codec)]
fn no_threads_routes_through_chain(#[case] cmd: ConsensusCmd) {
    let dir = TempDir::new().expect("temp dir");
    let input = write_input(cmd, dir.path());
    let output = dir.path().join("out.bam");

    let current_bin = Path::new(env!("CARGO_BIN_EXE_fgumi"));
    let out = run_consensus(current_bin, cmd, &input, &output, None);
    let stderr = String::from_utf8_lossy(&out.stderr);
    assert!(out.status.success(), "no-threads {} must succeed; stderr:\n{stderr}", cmd.name());
    assert!(
        stderr.contains("Total MI groups processed"),
        "a no-`--threads` {} must route through the chain (whose finalize hook logs \
         `Total MI groups processed`); the serial path (which logged `Total records processed`) \
         is retired. stderr:\n{stderr}",
        cmd.name()
    );
    assert!(
        !stderr.contains("Total records processed"),
        "a no-`--threads` {} must NOT emit the retired serial tail's `Total records processed` \
         line. stderr:\n{stderr}",
        cmd.name()
    );
    assert!(
        stderr.contains(cmd.starting_banner()),
        "the chain must still emit the `{}` banner; stderr:\n{stderr}",
        cmd.starting_banner()
    );
}

/// Output parity of the post-cutover chain against the pre-removal serial
/// baseline binary, across all three commands with and without `--stats` — plus
/// the always-available self-consistency oracle when no baseline is set.
#[rstest]
#[case::simplex_no_stats(ConsensusCmd::Simplex, false)]
#[case::simplex_with_stats(ConsensusCmd::Simplex, true)]
#[case::duplex_no_stats(ConsensusCmd::Duplex, false)]
#[case::duplex_with_stats(ConsensusCmd::Duplex, true)]
#[case::codec_no_stats(ConsensusCmd::Codec, false)]
#[case::codec_with_stats(ConsensusCmd::Codec, true)]
fn cutover_matches_baseline(#[case] cmd: ConsensusCmd, #[case] with_stats: bool) {
    let dir = TempDir::new().expect("temp dir");
    let input = write_input(cmd, dir.path());

    let current_bin = Path::new(env!("CARGO_BIN_EXE_fgumi"));
    let current_out = dir.path().join("current.bam");
    let current_tsv = dir.path().join("current.stats.tsv");
    let current_stats = with_stats.then(|| current_tsv.clone());
    let current = run_consensus(current_bin, cmd, &input, &current_out, current_stats.as_deref());
    let current_stderr = String::from_utf8_lossy(&current.stderr);
    assert!(
        current.status.success(),
        "current {} must succeed; stderr:\n{current_stderr}",
        cmd.name()
    );

    if let Some(baseline) = baseline_bin() {
        let baseline_out = dir.path().join("baseline.bam");
        let baseline_tsv = dir.path().join("baseline.stats.tsv");
        let baseline_stats = with_stats.then(|| baseline_tsv.clone());
        let base = run_consensus(&baseline, cmd, &input, &baseline_out, baseline_stats.as_deref());
        assert!(
            base.status.success(),
            "baseline {} failed; stderr:\n{}",
            cmd.name(),
            String::from_utf8_lossy(&base.stderr)
        );

        assert_eq!(
            decompressed_records_without_pg(&current_out),
            decompressed_records_without_pg(&baseline_out),
            "chain {} output diverges from the pre-removal serial baseline binary ({}) after \
             stripping @PG — a real cutover parity bug, not something to relax",
            cmd.name(),
            baseline.display(),
        );
        if with_stats {
            assert_eq!(
                std::fs::read_to_string(&current_tsv).expect("current stats tsv"),
                std::fs::read_to_string(&baseline_tsv).expect("baseline stats tsv"),
                "chain {} --stats TSV diverges from the serial baseline binary",
                cmd.name()
            );
        }
    } else {
        eprintln!(
            "SKIP baseline half of cutover_matches_baseline[{}]: FGUMI_BASELINE_BIN is unset or \
             does not name an existing file — running self-consistency oracle instead",
            cmd.name()
        );
        assert_self_consistent(cmd, &current_out, current_stats.as_deref());
    }
}

/// Always-available oracle used when no baseline binary is set: since
/// `FGUMI_BASELINE_BIN` is unset in CI, this fallback is the ONLY oracle CI
/// ever runs for this test. It asserts real content, not just "non-empty":
/// the exact consensus record count and each record's exact `SEQ`, both
/// fully determined by the fixture (every family / molecule agrees on every
/// base, so the consensus reproduces it exactly with no injected
/// disagreement to reason about) — plus, with `--stats`, a non-empty,
/// tab-delimited stats file with a header row.
///
/// Expected values were confirmed by running the chain against each fixture
/// and reading back the output (not derived by inspection alone):
/// - Simplex (`write_simplex_input`): 3 MI groups (`fam1`/`fam2`/`fam3`), one
///   consensus fragment each, `SEQ` equal to the family's uniform sequence.
/// - Duplex (`write_duplex_input`): 2 molecules, each yielding one consensus
///   template (R1 + R2), so 4 records; both ends reproduce the molecule's
///   uniform sequence since A/B strands agree exactly.
/// - Codec (`write_codec_input`): 1 molecule (4 overlapping pairs sharing one
///   `MI`), one consensus fragment, `SEQ` equal to the pairs' uniform R1/R2
///   sequence.
fn assert_self_consistent(cmd: ConsensusCmd, output: &Path, stats: Option<&Path>) {
    let (_, out_records) = read_bam_output(output);

    let expected_seqs: &[&str] = match cmd {
        ConsensusCmd::Simplex => &["ACGTACGTAC", "TTTTAAAAGG", "CCCCGGGGTT"],
        ConsensusCmd::Duplex => &["ACGTACGT", "ACGTACGT", "TTGGCCAA", "TTGGCCAA"],
        ConsensusCmd::Codec => &["ACGTACGT"],
    };

    assert_eq!(
        out_records.len(),
        expected_seqs.len(),
        "{} chain output must contain exactly {} consensus record(s) for this fixture; got {} \
         record(s)",
        cmd.name(),
        expected_seqs.len(),
        out_records.len()
    );

    for (i, (record, expected_seq)) in out_records.iter().zip(expected_seqs).enumerate() {
        let seq = String::from_utf8_lossy(record.sequence().as_ref()).into_owned();
        assert_eq!(
            &seq,
            expected_seq,
            "{} consensus record {i} must reproduce the fixture's known agreeing sequence \
             (every input read agrees at every base, so the consensus is fully determined); \
             got {seq}",
            cmd.name()
        );
    }

    if let Some(path) = stats {
        let tsv = std::fs::read_to_string(path).expect("read stats");
        let mut lines = tsv.lines();
        let header = lines.next().unwrap_or("");
        assert!(
            header.contains('\t'),
            "{} stats TSV must have a tab-delimited header row; got:\n{tsv}",
            cmd.name()
        );
        assert!(
            lines.next().is_some(),
            "{} stats TSV must have at least one data row; got:\n{tsv}",
            cmd.name()
        );
    }
}
