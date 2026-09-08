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

use fgumi_lib::sam::SamTag;
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

/// Runs `<bin> <cmd> -i <input> -o <output> [--stats <stats>] [--rejects <rejects>]
/// --min-reads 2 --compression-level 1` (no `--threads`, so the current build takes
/// the post-cutover chain path and the baseline takes its serial path) with
/// `RUST_LOG=info`, and returns the process output for stderr assertions.
fn run_consensus(
    bin: &Path,
    cmd: ConsensusCmd,
    input: &Path,
    output: &Path,
    stats: Option<&Path>,
    rejects: Option<&Path>,
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
    if let Some(rejects) = rejects {
        command.args([OsStr::new("--rejects"), rejects.as_os_str()]);
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
    let out = run_consensus(current_bin, cmd, &input, &output, None, None);
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
    // Always exercise the `--rejects` fan-out path (a distinct output the chain
    // writes) so a cutover regression in it is observable. These fixtures are
    // all-agreeing and above `--min-reads 2`, so the rejects BAM is expected to
    // be header-only (zero records) — asserted below.
    let current_rejects = dir.path().join("current.rejects.bam");
    let current_stats = with_stats.then(|| current_tsv.clone());
    let current = run_consensus(
        current_bin,
        cmd,
        &input,
        &current_out,
        current_stats.as_deref(),
        Some(&current_rejects),
    );
    let current_stderr = String::from_utf8_lossy(&current.stderr);
    assert!(
        current.status.success(),
        "current {} must succeed; stderr:\n{current_stderr}",
        cmd.name()
    );

    if let Some(baseline) = baseline_bin() {
        let baseline_out = dir.path().join("baseline.bam");
        let baseline_tsv = dir.path().join("baseline.stats.tsv");
        let baseline_rejects = dir.path().join("baseline.rejects.bam");
        let baseline_stats = with_stats.then(|| baseline_tsv.clone());
        let base = run_consensus(
            &baseline,
            cmd,
            &input,
            &baseline_out,
            baseline_stats.as_deref(),
            Some(&baseline_rejects),
        );
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
        // The rejects fan-out must match the serial baseline too (header + the
        // empty-rejects case here); reject *record* parity on reject-producing
        // input is pinned by the per-command suites (e.g. `test_duplex_chain_rejects_parity`).
        assert_eq!(
            decompressed_records_without_pg(&current_rejects),
            decompressed_records_without_pg(&baseline_rejects),
            "chain {} rejects output diverges from the serial baseline binary after stripping @PG",
            cmd.name(),
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
        assert_self_consistent(cmd, &current_out, current_stats.as_deref(), Some(&current_rejects));
    }
}

/// One consensus record's fully fixture-determined content, used as the
/// self-consistency oracle's expectation. Every field is pinned to an exact
/// value so a regression in any of them — not just `SEQ` — fails the test.
struct ExpectedRecord {
    /// Consensus read name (`<prefix>:<counter>`; the test prefix is empty, so
    /// e.g. `:1`).
    name: &'static str,
    /// Consensus `SEQ` — the family/molecule's uniform sequence.
    seq: &'static str,
    /// Uniform per-base Phred quality carried across the whole read (the inputs
    /// agree at every base, so the caller emits one capped value throughout).
    qual: u8,
    /// Exact SAM FLAG bits. Consensus reads are unmapped: a simplex/codec
    /// fragment is `UNMAPPED` alone; each duplex template emits an
    /// `UNMAPPED | MATE_UNMAPPED | SEGMENTED` pair with `FIRST_SEGMENT` (R1) or
    /// `LAST_SEGMENT` (R2). Pinning the exact bits catches a regression in
    /// pairing or the mapped/unmapped state that field-subset checks miss.
    flags: u16,
    /// Total consensus depth (`cD`): the number of raw reads supporting the
    /// molecule. Also every entry of the per-base depth array.
    total_depth: i64,
    /// For duplex/codec, the per-strand `(aD, bD)` depths; `None` for simplex
    /// (which instead carries the single-strand `cd` per-base depth array).
    strand_depths: Option<(i64, i64)>,
}

/// `UNMAPPED` (0x4): the flag a simplex/codec consensus fragment carries.
const F_UNMAPPED: u16 = 0x4;
/// A duplex R1 consensus: `SEGMENTED | UNMAPPED | MATE_UNMAPPED | FIRST_SEGMENT`.
const F_DUPLEX_R1: u16 = 0x1 | 0x4 | 0x8 | 0x40;
/// A duplex R2 consensus: `SEGMENTED | UNMAPPED | MATE_UNMAPPED | LAST_SEGMENT`.
const F_DUPLEX_R2: u16 = 0x1 | 0x4 | 0x8 | 0x80;

/// The expected records for each command's fixture, confirmed by running the
/// chain against each fixture and reading the output back (not by inspection):
/// - Simplex (`write_simplex_input`): 3 MI groups (`fam1`/`fam2`/`fam3`) of
///   depth 5/4/3, one consensus fragment each.
/// - Duplex (`write_duplex_input`): 2 molecules (depth 8/6, split 4+4 / 3+3
///   across strands), each yielding one consensus template (R1 + R2) → 4 records.
/// - Codec (`write_codec_input`): 1 molecule (4 overlapping pairs sharing one
///   `MI`, depth 8, split 4+4), one consensus fragment.
fn expected_records(cmd: ConsensusCmd) -> Vec<ExpectedRecord> {
    let rec = |name, seq: &'static str, qual, flags, total_depth, strand_depths| ExpectedRecord {
        name,
        seq,
        qual,
        flags,
        total_depth,
        strand_depths,
    };
    match cmd {
        ConsensusCmd::Simplex => vec![
            rec(":1", "ACGTACGTAC", 45, F_UNMAPPED, 5, None),
            rec(":2", "TTTTAAAAGG", 45, F_UNMAPPED, 4, None),
            rec(":3", "CCCCGGGGTT", 45, F_UNMAPPED, 3, None),
        ],
        ConsensusCmd::Duplex => vec![
            rec(":1", "ACGTACGT", 90, F_DUPLEX_R1, 8, Some((4, 4))),
            rec(":1", "ACGTACGT", 90, F_DUPLEX_R2, 8, Some((4, 4))),
            rec(":2", "TTGGCCAA", 90, F_DUPLEX_R1, 6, Some((3, 3))),
            rec(":2", "TTGGCCAA", 90, F_DUPLEX_R2, 6, Some((3, 3))),
        ],
        ConsensusCmd::Codec => vec![rec(":UMI001", "ACGTACGT", 90, F_UNMAPPED, 8, Some((4, 4)))],
    }
}

/// The *complete* fixture-determined `--stats` TSV row set (`key` → exact
/// `value` string) for each command — every row the tool emits, not a subset.
/// Every raw read is used and none rejected, so all counts are fully determined
/// by the fixture size and every rejection reason is `0`. [`assert_stats_rows`]
/// asserts this map is exactly the emitted key set (a missing OR an unexpected
/// row fails), so a regression that adds, drops, or mis-values any row — a new
/// rejection reason, a renamed key — is caught.
fn expected_stats(cmd: ConsensusCmd) -> &'static [(&'static str, &'static str)] {
    match cmd {
        ConsensusCmd::Simplex => &[
            ("raw_reads_considered", "12"),
            ("raw_reads_rejected", "0"),
            ("raw_reads_used", "12"),
            ("frac_raw_reads_used", "1.000000"),
            ("raw_reads_rejected_for_insufficient_support", "0"),
            ("raw_reads_rejected_for_minority_alignment", "0"),
            ("raw_reads_rejected_for_orphan_consensus", "0"),
            ("raw_reads_rejected_for_zero_bases_post_trimming", "0"),
            ("consensus_reads_emitted", "3"),
        ],
        ConsensusCmd::Duplex => &[
            ("raw_reads_considered", "28"),
            ("raw_reads_rejected", "0"),
            ("raw_reads_used", "28"),
            ("frac_raw_reads_used", "1.000000"),
            ("raw_reads_rejected_for_insufficient_support", "0"),
            ("raw_reads_rejected_for_minority_alignment", "0"),
            ("raw_reads_rejected_for_orphan_consensus", "0"),
            ("raw_reads_rejected_for_zero_bases_post_trimming", "0"),
            ("raw_reads_rejected_for_non_paired_reads", "0"),
            ("raw_reads_rejected_for_single_strand_only", "0"),
            ("raw_reads_rejected_for_potential_umi_collision", "0"),
            ("consensus_reads_emitted", "4"),
        ],
        ConsensusCmd::Codec => &[
            ("raw_reads_considered", "8"),
            ("raw_reads_rejected", "0"),
            ("raw_reads_used", "8"),
            ("frac_raw_reads_used", "1.000000"),
            ("raw_reads_rejected_for_insufficient_support", "0"),
            ("raw_reads_rejected_for_minority_alignment", "0"),
            ("raw_reads_rejected_for_orphan_consensus", "0"),
            ("raw_reads_rejected_for_zero_bases_post_trimming", "0"),
            ("raw_reads_rejected_for_non_paired_reads", "0"),
            ("raw_reads_rejected_for_single_strand_only", "0"),
            ("raw_reads_rejected_for_potential_umi_collision", "0"),
            ("raw_reads_rejected_for_r1_r2_overlap_too_short", "0"),
            ("raw_reads_rejected_for_indel_error_between_strands", "0"),
            ("raw_reads_rejected_for_high_duplex_disagreement", "0"),
            ("raw_reads_rejected_for_clip_overlap_failed", "0"),
            ("raw_reads_rejected_for_not_primary_fr_pair", "0"),
            ("consensus_reads_emitted", "1"),
            ("consensus_reads_rejected_hdd", "0"),
            ("consensus_bases_emitted", "8"),
            ("consensus_duplex_bases_emitted", "8"),
            ("duplex_disagreement_base_count", "0"),
            ("duplex_disagreement_rate", "0.000000"),
        ],
    }
}

/// Reads a data-field integer tag (any SAM integer width) off a `RecordBuf` as
/// `i64`. Panics if the tag is absent or not an integer.
fn int_tag_of(record: &noodles::sam::alignment::RecordBuf, tag: SamTag) -> i64 {
    use noodles::sam::alignment::record_buf::data::field::Value;
    match record.data().get(&tag.to_noodles_tag()) {
        Some(Value::Int8(n)) => i64::from(*n),
        Some(Value::UInt8(n)) => i64::from(*n),
        Some(Value::Int16(n)) => i64::from(*n),
        Some(Value::UInt16(n)) => i64::from(*n),
        Some(Value::Int32(n)) => i64::from(*n),
        Some(Value::UInt32(n)) => i64::from(*n),
        other => panic!("tag {tag:?} must be an integer, got {other:?}"),
    }
}

/// Reads a data-field `Float` tag off a `RecordBuf`. Panics if absent or not a float.
fn float_tag_of(record: &noodles::sam::alignment::RecordBuf, tag: SamTag) -> f32 {
    use noodles::sam::alignment::record_buf::data::field::Value;
    match record.data().get(&tag.to_noodles_tag()) {
        Some(Value::Float(f)) => *f,
        other => panic!("tag {tag:?} must be a float, got {other:?}"),
    }
}

/// Reads a data-field `Int16` array tag off a `RecordBuf`. Panics if absent or
/// not an `Int16` array (the per-base depth/error arrays fgumi emits).
fn int16_array_tag_of(record: &noodles::sam::alignment::RecordBuf, tag: SamTag) -> Vec<i16> {
    use noodles::sam::alignment::record_buf::data::field::Value;
    use noodles::sam::alignment::record_buf::data::field::value::Array;
    match record.data().get(&tag.to_noodles_tag()) {
        Some(Value::Array(Array::Int16(values))) => values.clone(),
        other => panic!("tag {tag:?} must be an Int16 array, got {other:?}"),
    }
}

/// Always-available oracle used when no baseline binary is set: since
/// `FGUMI_BASELINE_BIN` is unset in CI, this fallback is the ONLY oracle CI
/// ever runs for this test. It asserts real content, not just "non-empty" or
/// "SEQ": each record's exact read name, SAM FLAG bits, the unmapped-consensus
/// alignment invariant (no position, MAPQ 0, empty CIGAR, no mate), `SEQ`,
/// per-base qualities, command-specific consensus depth/error tags (scalar
/// `cD`/`cE` plus the per-base `cd` or `ad`/`bd` arrays), and the complete
/// data-tag key set; with `--stats`, the complete fixture-determined stats map
/// (a missing OR unexpected row fails); and with `--rejects`, that the rejects
/// fan-out BAM
/// is header-only (zero records) — these fixtures are all-agreeing and above
/// `--min-reads 2`, so every raw read is consumed and nothing is rejected, so a
/// regression that spuriously fanned reads into rejects would fail here. (Reject
/// *record* parity on reject-producing input is pinned by the per-command suites,
/// e.g. `test_duplex_chain_rejects_parity`.) All expectations are pinned in
/// [`expected_records`] / [`expected_stats`].
fn assert_self_consistent(
    cmd: ConsensusCmd,
    output: &Path,
    stats: Option<&Path>,
    rejects: Option<&Path>,
) {
    let (_, out_records) = read_bam_output(output);

    let expected = expected_records(cmd);

    assert_eq!(
        out_records.len(),
        expected.len(),
        "{} chain output must contain exactly {} consensus record(s) for this fixture; got {} \
         record(s)",
        cmd.name(),
        expected.len(),
        out_records.len()
    );

    for (i, (record, exp)) in out_records.iter().zip(&expected).enumerate() {
        assert_consensus_record(cmd, &format!("{} consensus record {i}", cmd.name()), record, exp);
    }

    if let Some(path) = stats {
        assert_stats_rows(cmd, path);
    }

    if let Some(path) = rejects {
        let (_, rejects_records) = read_bam_output(path);
        assert!(
            rejects_records.is_empty(),
            "{} rejects fan-out BAM must be header-only for an all-agreeing, above-min-reads \
             fixture (nothing is rejected); got {} record(s)",
            cmd.name(),
            rejects_records.len()
        );
    }
}

/// Asserts one emitted consensus record against its fixture-determined
/// expectation across every behaviorally relevant field: read name, exact SAM
/// FLAG bits, the unmapped-consensus alignment invariant (no position, MAPQ 0,
/// empty CIGAR, no mate reference/position), `SEQ`, uniform per-base qualities,
/// the command-specific consensus depth/error tags (scalar `cD`/`cE` plus the
/// per-base `cd` for simplex or `aD`/`bD`/`ad`/`bd` for duplex/codec), and the
/// complete set of data-tag keys (so an added or dropped tag also fails). This
/// keeps the CI-only fallback oracle from passing a regression in a field the
/// value-level checks above do not name.
fn assert_consensus_record(
    cmd: ConsensusCmd,
    label: &str,
    record: &noodles::sam::alignment::RecordBuf,
    exp: &ExpectedRecord,
) {
    // Read name — deterministic `<prefix>:<counter>` (the test prefix is empty).
    let name =
        record.name().map(|n| String::from_utf8_lossy(n.as_ref()).into_owned()).unwrap_or_default();
    assert_eq!(name, exp.name, "{label} must carry the fixture-determined read name");

    // FLAG — exact bits. Consensus reads are unmapped; duplex additionally pins
    // the paired FIRST/LAST_SEGMENT roles. A subset check on individual fields
    // would miss a regression that flipped pairing or the unmapped state.
    assert_eq!(
        u16::from(record.flags()),
        exp.flags,
        "{label} must carry FLAG bits {:#06x}; got {:#06x}",
        exp.flags,
        u16::from(record.flags())
    );

    // Alignment invariant — a consensus read is unmapped (constant across every
    // record this test emits).
    assert_unmapped_alignment_invariant(label, record);

    // SEQ — the family/molecule's uniform sequence (every input read agrees).
    let seq = String::from_utf8_lossy(record.sequence().as_ref()).into_owned();
    assert_eq!(
        seq, exp.seq,
        "{label} must reproduce the fixture's known agreeing sequence (every input read agrees \
         at every base, so the consensus is fully determined)"
    );

    // QUAL — one capped Phred value repeated across the whole read, since the
    // inputs agree at every base. A regression that dropped or rescaled the
    // consensus qualities (which SEQ-only checking cannot see) fails here.
    let quals: Vec<u8> = record.quality_scores().as_ref().to_vec();
    assert_eq!(
        quals,
        vec![exp.qual; exp.seq.len()],
        "{label} must carry uniform Phred {} across all {} bases",
        exp.qual,
        exp.seq.len()
    );

    // Consensus depth/error tags — command-specific and fixture-determined.
    assert_consensus_depth_tags(label, record, exp);

    // Complete tag-key set — the value-level checks above assert a hand-picked
    // subset; this pins the *whole* set of data tags the record carries, so a
    // regression that adds or drops any tag (even one no assertion names) fails.
    let got_keys = record_tag_keys(record);
    let want_keys = {
        let mut keys: Vec<String> =
            expected_tag_keys(cmd).iter().map(|s| (*s).to_owned()).collect();
        keys.sort();
        keys
    };
    assert_eq!(
        got_keys, want_keys,
        "{label} must carry exactly the fixture-determined data tags {want_keys:?}; got {got_keys:?}"
    );
}

/// Asserts the unmapped-consensus alignment invariant: no reference position,
/// MAPQ 0, empty CIGAR, and no mate reference/position. Constant for every
/// consensus record this test emits.
fn assert_unmapped_alignment_invariant(label: &str, record: &noodles::sam::alignment::RecordBuf) {
    assert!(
        record.alignment_start().is_none(),
        "{label} (unmapped consensus) must have no alignment position; got {:?}",
        record.alignment_start()
    );
    assert_eq!(
        record.mapping_quality().map_or(0, |m| m.get()),
        0,
        "{label} (unmapped consensus) must have MAPQ 0"
    );
    assert!(
        record.cigar().as_ref().is_empty(),
        "{label} (unmapped consensus) must have an empty CIGAR; got {:?}",
        record.cigar()
    );
    assert!(
        record.mate_reference_sequence_id().is_none(),
        "{label} (unmapped consensus) must have no mate reference"
    );
    assert!(
        record.mate_alignment_start().is_none(),
        "{label} (unmapped consensus) must have no mate alignment position"
    );
}

/// Asserts the fixture-determined consensus depth/error tags: the scalar total
/// depth `cD` and zero error rate `cE`, plus either the single-strand per-base
/// depth array `cd` (simplex) or the per-strand `aD`/`bD` scalars and `ad`/`bd`
/// per-base arrays (duplex/codec).
fn assert_consensus_depth_tags(
    label: &str,
    record: &noodles::sam::alignment::RecordBuf,
    exp: &ExpectedRecord,
) {
    // `cD` is the total supporting-read depth; `cE` is the per-read error rate,
    // which is exactly 0.0 because every input read agrees.
    assert_eq!(
        int_tag_of(record, SamTag::CD),
        exp.total_depth,
        "{label} must report total consensus depth (cD) == {}",
        exp.total_depth
    );
    let cerr = float_tag_of(record, SamTag::CE);
    assert!(
        cerr.abs() < f32::EPSILON,
        "{label} must report zero consensus error rate (cE) for an all-agreeing family; got {cerr}"
    );

    let base_array =
        |depth: i64| vec![i16::try_from(depth).expect("depth fits i16"); exp.seq.len()];
    match exp.strand_depths {
        // Simplex: a single-strand caller carries the per-base depth array `cd`.
        None => assert_eq!(
            int16_array_tag_of(record, SamTag::CD_BASES),
            base_array(exp.total_depth),
            "{label} per-base depth array (cd) must be [{}; {}]",
            exp.total_depth,
            exp.seq.len()
        ),
        // Duplex/codec: the two single-strand consensuses carry `aD`/`bD` scalar
        // depths and `ad`/`bd` per-base arrays.
        Some((a, b)) => {
            assert_eq!(
                int_tag_of(record, SamTag::AD),
                a,
                "{label} must report A-strand depth (aD) == {a}"
            );
            assert_eq!(
                int_tag_of(record, SamTag::BD),
                b,
                "{label} must report B-strand depth (bD) == {b}"
            );
            assert_eq!(
                int16_array_tag_of(record, SamTag::AD_BASES),
                base_array(a),
                "{label} A-strand per-base depth array (ad) must be [{a}; {}]",
                exp.seq.len()
            );
            assert_eq!(
                int16_array_tag_of(record, SamTag::BD_BASES),
                base_array(b),
                "{label} B-strand per-base depth array (bd) must be [{b}; {}]",
                exp.seq.len()
            );
        }
    }
}

/// The sorted data-tag keys (two-character SAM tags) present on `record`.
fn record_tag_keys(record: &noodles::sam::alignment::RecordBuf) -> Vec<String> {
    let mut keys: Vec<String> = record
        .data()
        .iter()
        .map(|(tag, _)| String::from_utf8_lossy(&<[u8; 2]>::from(tag)).into_owned())
        .collect();
    keys.sort();
    keys
}

/// The complete set of data-tag keys a consensus record carries for each
/// command's fixture, confirmed by reading the chain output back. Simplex
/// records add the raw-UMI `RX` and the single-strand `cd`/`ce` per-base arrays;
/// duplex and codec instead carry the matched per-strand `a*`/`b*` tags. Used to
/// assert the *whole* tag set, not just the values the record-level checks name.
fn expected_tag_keys(cmd: ConsensusCmd) -> &'static [&'static str] {
    match cmd {
        ConsensusCmd::Simplex => &["RG", "MI", "RX", "cD", "cM", "cE", "cd", "ce"],
        // Duplex and codec emit the same 19-tag set (the codec fixture yields a
        // fragment consensus, but with the identical per-strand tag surface).
        ConsensusCmd::Duplex | ConsensusCmd::Codec => &[
            "RG", "MI", "cD", "cM", "cE", "aD", "aM", "aE", "bD", "bM", "bE", "ac", "bc", "ad",
            "bd", "ae", "be", "aq", "bq",
        ],
    }
}

/// Asserts the `--stats` TSV matches the *complete* fixture-determined map in
/// [`expected_stats`]: every expected row is present with its exact value AND no
/// unexpected row appears. Comparing the whole key set — not a subset — means a
/// regression that adds a row (e.g. a new rejection reason), drops one, renames
/// a key, or mis-values any of them fails the oracle rather than slipping past a
/// hand-picked check.
fn assert_stats_rows(cmd: ConsensusCmd, path: &Path) {
    let tsv = std::fs::read_to_string(path).expect("read stats");
    let rows: std::collections::BTreeMap<&str, &str> = tsv
        .lines()
        .skip(1) // header: `key\tvalue\tdescription`
        .filter_map(|line| {
            let mut cols = line.split('\t');
            Some((cols.next()?, cols.next()?))
        })
        .collect();
    let want: std::collections::BTreeMap<&str, &str> =
        expected_stats(cmd).iter().copied().collect();
    assert_eq!(
        rows,
        want,
        "{} --stats TSV must be exactly the fixture-determined row set (a missing OR unexpected \
         row fails); got:\n{tsv}",
        cmd.name()
    );
}
