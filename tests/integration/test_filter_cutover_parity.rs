//! Parity gate for the `filter` command's single-threaded-path retirement (C4):
//! `Filter::execute` no longer has a serial in-process unified-pipeline loop
//! reached when `--threads` is absent — it *always* routes through the
//! declarative chain builder. This test proves that cutover lost nothing
//! user-observable.
//!
//! Two independent things are checked here:
//!
//! 1. **The cutover actually happened** (`filter_no_threads_routes_through_chain`).
//!    A no-`--threads` run now emits the chain-only `"Using pipeline with N
//!    threads"` banner that `ChainBuilder::add_filter` logs and the retired
//!    unified-pipeline tail never did. This is the genuine RED/GREEN
//!    discriminator: before the removal a no-`--threads` run took the serial
//!    path and printed no such line; after it, the chain does.
//!
//! 2. **Output parity with the pre-removal serial path**
//!    (`cutover_matches_baseline`). The current build's `filter` output — kept
//!    records and (when written) the rejects BAM (both byte-identical modulo the
//!    `@PG` line) and the `--stats` TSV — must match the frozen owned-serial
//!    baseline binary. The baseline path comes from `FGUMI_BASELINE_BIN`; when it
//!    is unset (or names a missing file) the case degrades to a self-consistency
//!    oracle (min-reads filtering, base masking, and NM regeneration asserted
//!    directly) rather than skipping — the exact fallback discipline of
//!    `test_sort_cutover_parity.rs`. Cases cover `--min-reads` filtering (records
//!    kept vs rejected), `--ref`-based NM/UQ/MD regeneration, base masking, and
//!    both `filter-by-template` modes, plus `--rejects` and `--stats`.
//!
//! **Not a RED/GREEN gate for the parity half.** Because the removed serial loop
//! and the chain were already output-equivalent (see the in-tree
//! worker-count-independence tests in `test_filter_command.rs`), the baseline
//! byte-parity check passes on both sides of the change — like the sort cutover
//! it guards equivalence, it does not observe a regression the cutover
//! introduces. The chain-banner check in (1) is the part that flips RED→GREEN
//! across the removal.

use std::ffi::OsStr;
use std::path::Path;
use std::process::Command;

use fgumi_raw_bam::{RawRecord, SamBuilder, flags};
use rstest::rstest;
use tempfile::TempDir;

use crate::helpers::bam_generator::{create_minimal_header, create_test_reference, write_bam};
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
fn baseline_bin() -> Option<std::path::PathBuf> {
    let path = std::path::PathBuf::from(std::env::var_os("FGUMI_BASELINE_BIN")?);
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
/// binary-layout regression by normalizing it away. filter regenerates NM/UQ/MD
/// and masks bases in place, so its output record bytes must match the pre-removal
/// path exactly, which is what this comparison pins.
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

/// One mapped consensus read: name, `flags`, position, per-consensus depth (`cD`)
/// and error (`cE`), and per-base depth (`cd`)/error (`ce`) arrays. The sequence
/// is 8 bases; a mismatch to the reference or a masked base changes the
/// regenerated `NM`, and a low per-base depth drives base masking.
fn consensus_read(
    name: &str,
    flags: u16,
    pos: i32,
    depth: i32,
    seq: &[u8],
    per_base_depth: &[u16],
) -> RawRecord {
    let mut b = SamBuilder::new();
    b.read_name(name.as_bytes())
        .flags(flags)
        .ref_id(0)
        .pos(pos)
        .mapq(60)
        .cigar_ops(&[u32::try_from(seq.len()).expect("seq len fits u32") << 4]) // <len>M
        .sequence(seq)
        .qualities(&vec![35u8; seq.len()]);
    b.add_int_tag(SamTag::CD, depth).add_float_tag(SamTag::CE, 0.0_f32);
    b.add_array_u16(SamTag::CD_BASES, per_base_depth)
        .add_array_u16(SamTag::CE_BASES, &vec![0u16; seq.len()]);
    b.build()
}

/// The parity corpus, at reference-aligned positions (`create_test_reference`
/// writes `ACGTACGT` repeated, so an 8-base window at any multiple of 8 is
/// `ACGTACGT`). Exercises every behavior the brief names in one input.
///
/// Four single-end reads (each its own template, so they behave identically in
/// both `--filter-by-template` modes):
/// - `pass`: depth 10, exact match, no masking — kept, `NM` 0.
/// - `low_depth`: depth 1 (< `--min-reads 3`) — rejected by the read-level filter.
/// - `masked`: depth 10 but one per-base depth 1 (< 3) — that base masked to `N`,
///   read still kept (1 `N` of 8 = 0.125 < `--max-no-call-fraction 0.2`), `NM`
///   regenerated over the masked base.
/// - `mismatch`: depth 10, first base `T` vs reference `A` — kept, `NM` 1 after
///   regeneration.
///
/// Plus one paired template (`pair` R1 + R2, adjacent as query-grouped input
/// requires) whose R1 passes (depth 10) but R2 fails (depth 1). This is the
/// record that makes `--filter-by-template` genuinely different from single-read
/// mode: template mode drops the *whole* template (both mates rejected, R1 too,
/// under fgbio's "all primaries must pass" rule), while single-read mode keeps R1
/// and drops only R2. Without it the `template` and `single_read` cases would
/// exercise identical behavior.
fn build_parity_records() -> Vec<RawRecord> {
    let paired_r1 = flags::PAIRED | flags::FIRST_SEGMENT;
    let paired_r2 = flags::PAIRED | flags::LAST_SEGMENT;
    vec![
        consensus_read("pass", 0, 0, 10, b"ACGTACGT", &[10; 8]),
        consensus_read("low_depth", 0, 8, 1, b"ACGTACGT", &[1; 8]),
        consensus_read("masked", 0, 16, 10, b"ACGTACGT", &[10, 10, 10, 1, 10, 10, 10, 10]),
        consensus_read("mismatch", 0, 24, 10, b"TCGTACGT", &[10; 8]),
        consensus_read("pair", paired_r1, 32, 10, b"ACGTACGT", &[10; 8]),
        consensus_read("pair", paired_r2, 40, 1, b"ACGTACGT", &[1; 8]),
    ]
}

/// Writes the parity corpus to `dir/in.bam` (query-grouped header, which filter
/// requires) and returns its path.
fn write_parity_input(dir: &Path) -> std::path::PathBuf {
    let input = dir.join("in.bam");
    write_bam(&input, &create_minimal_header("chr1", 10_000), &build_parity_records());
    input
}

/// Runs `<bin> filter -i <input> -o <output> --ref <ref> --min-reads 3
/// --max-no-call-fraction 0.2 --compression-level 1 --filter-by-template <bool>
/// [--rejects <r>] [--stats <s>]` (no `--threads`, so the current build takes the
/// post-cutover chain path and the baseline takes its serial path) with
/// `RUST_LOG=info`, and returns the process output for stderr assertions.
fn run_filter(
    bin: &Path,
    input: &Path,
    output: &Path,
    ref_path: &Path,
    filter_by_template: bool,
    rejects: Option<&Path>,
    stats: Option<&Path>,
) -> std::process::Output {
    let mut cmd = Command::new(bin);
    cmd.env("RUST_LOG", "info").args([
        OsStr::new("filter"),
        OsStr::new("-i"),
        input.as_os_str(),
        OsStr::new("-o"),
        output.as_os_str(),
        OsStr::new("--ref"),
        ref_path.as_os_str(),
        OsStr::new("--min-reads"),
        OsStr::new("3"),
        OsStr::new("--max-no-call-fraction"),
        OsStr::new("0.2"),
        OsStr::new("--compression-level"),
        OsStr::new("1"),
        OsStr::new("--filter-by-template"),
        OsStr::new(if filter_by_template { "true" } else { "false" }),
    ]);
    if let Some(r) = rejects {
        cmd.args([OsStr::new("--rejects"), r.as_os_str()]);
    }
    if let Some(s) = stats {
        cmd.args([OsStr::new("--stats"), s.as_os_str()]);
    }
    cmd.output().unwrap_or_else(|e| panic!("failed to spawn `{}` filter: {e}", bin.display()))
}

/// A no-`--threads` run now routes through the declarative chain, which logs the
/// `"Using pipeline with N threads"` banner from `ChainBuilder::add_filter`. The
/// retired unified-pipeline tail logged no such line, so this is the RED
/// (pre-removal) → GREEN (post-removal) discriminator for the cutover.
#[test]
fn filter_no_threads_routes_through_chain() {
    let dir = TempDir::new().expect("temp dir");
    let input = write_parity_input(dir.path());
    let ref_path = create_test_reference(dir.path());
    let output = dir.path().join("out.bam");

    let current_bin = Path::new(env!("CARGO_BIN_EXE_fgumi"));
    let out = run_filter(current_bin, &input, &output, &ref_path, true, None, None);
    let stderr = String::from_utf8_lossy(&out.stderr);
    assert!(out.status.success(), "no-threads filter must succeed; stderr:\n{stderr}");
    assert!(
        stderr.contains("Using pipeline with 1 threads"),
        "a no-`--threads` filter must route through the chain (which logs the pipeline banner); \
         the serial path is retired. stderr:\n{stderr}"
    );
    assert!(
        stderr.contains("Starting Filter"),
        "the chain must still emit the `Starting Filter` banner; stderr:\n{stderr}"
    );
}

/// Output parity of the post-cutover chain against the pre-removal serial
/// baseline binary, across `filter-by-template` modes and the `--rejects` /
/// `--stats` outputs — plus the always-available self-consistency oracle when no
/// baseline is set.
///
/// `#[case]` args: a label, `filter_by_template`, `with_rejects`, `with_stats`.
#[rstest]
#[case::template(true, false, false)]
#[case::single_read(false, false, false)]
#[case::template_rejects(true, true, false)]
#[case::template_stats(true, false, true)]
#[case::single_read_rejects(false, true, false)]
#[case::single_read_stats(false, false, true)]
fn cutover_matches_baseline(
    #[case] filter_by_template: bool,
    #[case] with_rejects: bool,
    #[case] with_stats: bool,
) {
    let dir = TempDir::new().expect("temp dir");
    let input = write_parity_input(dir.path());
    let ref_path = create_test_reference(dir.path());

    let current_bin = Path::new(env!("CARGO_BIN_EXE_fgumi"));
    let current_out = dir.path().join("current.bam");
    let current_rejects = dir.path().join("current.rejects.bam");
    let current_stats = dir.path().join("current.stats.txt");
    let cur_rej = with_rejects.then(|| current_rejects.clone());
    let cur_stats = with_stats.then(|| current_stats.clone());
    let current = run_filter(
        current_bin,
        &input,
        &current_out,
        &ref_path,
        filter_by_template,
        cur_rej.as_deref(),
        cur_stats.as_deref(),
    );
    let current_stderr = String::from_utf8_lossy(&current.stderr);
    assert!(current.status.success(), "current filter must succeed; stderr:\n{current_stderr}");

    if let Some(baseline) = baseline_bin() {
        let baseline_out = dir.path().join("baseline.bam");
        let baseline_rejects = dir.path().join("baseline.rejects.bam");
        let baseline_stats = dir.path().join("baseline.stats.txt");
        let base_rej = with_rejects.then(|| baseline_rejects.clone());
        let base_stats = with_stats.then(|| baseline_stats.clone());
        let base = run_filter(
            &baseline,
            &input,
            &baseline_out,
            &ref_path,
            filter_by_template,
            base_rej.as_deref(),
            base_stats.as_deref(),
        );
        assert!(
            base.status.success(),
            "baseline filter failed; stderr:\n{}",
            String::from_utf8_lossy(&base.stderr)
        );

        assert_eq!(
            decompressed_records_without_pg(&current_out),
            decompressed_records_without_pg(&baseline_out),
            "chain filter output diverges from the pre-removal serial baseline binary ({}) \
             after stripping @PG — a real cutover parity bug, not something to relax",
            baseline.display(),
        );
        if with_rejects {
            assert_eq!(
                decompressed_records_without_pg(&current_rejects),
                decompressed_records_without_pg(&baseline_rejects),
                "chain --rejects output diverges from the serial baseline binary"
            );
        }
        if with_stats {
            assert_eq!(
                std::fs::read_to_string(&current_stats).expect("current stats"),
                std::fs::read_to_string(&baseline_stats).expect("baseline stats"),
                "chain --stats TSV diverges from the serial baseline binary"
            );
        }
    } else {
        eprintln!(
            "SKIP baseline half of cutover_matches_baseline[template={filter_by_template}, \
             rejects={with_rejects}, stats={with_stats}]: FGUMI_BASELINE_BIN is unset or does \
             not name an existing file — running self-consistency oracle instead"
        );
        assert_self_consistent(
            filter_by_template,
            &current_out,
            cur_rej.as_deref(),
            cur_stats.as_deref(),
        );
    }
}

/// Always-available oracle used when no baseline binary is set: the chain's
/// output must show exactly the `--min-reads` filtering (including the
/// template-drop that distinguishes the two modes), base masking, and NM
/// regeneration the parity corpus prescribes.
///
/// Kept/rejected sets differ by mode because of the `pair` template (R1 passes,
/// R2 fails at depth 1):
/// - template mode: keeps `pass`, `masked`, `mismatch` (3); rejects `low_depth`
///   and both `pair` mates (the whole template drops even though R1 passed).
/// - single-read mode: additionally keeps `pair` R1 (4 kept); rejects `low_depth`
///   and `pair` R2 (2).
fn assert_self_consistent(
    filter_by_template: bool,
    output: &Path,
    rejects: Option<&Path>,
    stats: Option<&Path>,
) {
    let (total, passed, failed): (usize, usize, usize) =
        if filter_by_template { (6, 3, 3) } else { (6, 4, 2) };

    let (_, kept) = read_bam_output(output);
    assert_eq!(
        kept.len(),
        passed,
        "filter_by_template={filter_by_template}: unexpected kept count"
    );

    let names: Vec<String> = kept.iter().map(record_name).collect();
    assert!(names.contains(&"pass".to_string()), "the passing read must be kept");
    assert!(names.contains(&"masked".to_string()), "the base-masked read must be kept");
    assert!(names.contains(&"mismatch".to_string()), "the mismatch read must be kept");
    assert!(
        !names.contains(&"low_depth".to_string()),
        "the low-depth read must be filtered out by --min-reads"
    );
    // The template-drop discriminator: `pair` R1 survives only in single-read mode.
    let pair_kept = names.iter().filter(|n| *n == "pair").count();
    if filter_by_template {
        assert_eq!(pair_kept, 0, "template mode must drop the whole `pair` template (R1 too)");
    } else {
        assert_eq!(pair_kept, 1, "single-read mode must keep `pair` R1 and drop only R2");
    }

    // Base masking: the `masked` read carries exactly one N, at the base whose
    // per-base depth (index 3) fell below --min-reads — every other base intact.
    let masked = kept.iter().find(|r| record_name(r) == "masked").expect("masked read present");
    let seq: Vec<u8> = masked.sequence().as_ref().to_vec();
    let n_positions: Vec<usize> =
        seq.iter().enumerate().filter(|&(_, &b)| b == b'N').map(|(i, _)| i).collect();
    assert_eq!(n_positions, vec![3], "exactly base index 3 must be masked to N, nothing else");

    // --ref NM/UQ/MD regeneration: exact-match read is NM 0; the single-base
    // mismatch read is NM 1; masking base 3 of the `masked` read (an exact match
    // pre-masking) turns that base into an N/ref difference, so its NM regenerates
    // to 1 as well.
    let pass = kept.iter().find(|r| record_name(r) == "pass").expect("pass read present");
    assert_eq!(nm_tag(pass), Some(0), "exact-match read NM must be regenerated to 0");
    let mismatch =
        kept.iter().find(|r| record_name(r) == "mismatch").expect("mismatch read present");
    assert_eq!(nm_tag(mismatch), Some(1), "mismatch read NM must be regenerated to 1");
    assert_eq!(nm_tag(masked), Some(1), "masked read NM must regenerate over the masked base");

    if let Some(rejects_path) = rejects {
        let (_, rejected) = read_bam_output(rejects_path);
        assert_eq!(
            rejected.len(),
            failed,
            "filter_by_template={filter_by_template}: unexpected rejected count"
        );
        assert!(
            rejected.iter().any(|r| record_name(r) == "low_depth"),
            "low_depth must appear in the rejects output"
        );
    }

    if let Some(stats_path) = stats {
        let tsv = std::fs::read_to_string(stats_path).expect("read stats");
        assert!(tsv.contains(&format!("total_reads\t{total}")), "stats total_reads; got:\n{tsv}");
        assert!(
            tsv.contains(&format!("passed_reads\t{passed}")),
            "stats passed_reads; got:\n{tsv}"
        );
        assert!(
            tsv.contains(&format!("failed_reads\t{failed}")),
            "stats failed_reads; got:\n{tsv}"
        );
    }
}

/// The read name of a parsed `RecordBuf`.
fn record_name(record: &noodles::sam::alignment::RecordBuf) -> String {
    record.name().map(|n| String::from_utf8_lossy(n.as_ref()).into_owned()).unwrap_or_default()
}

/// Reads the integer `NM` aux tag off a parsed `RecordBuf`, or `None` when absent
/// or not an integer.
fn nm_tag(record: &noodles::sam::alignment::RecordBuf) -> Option<i64> {
    use noodles::sam::alignment::record::data::field::Tag;
    let nm: SamTag = "NM".parse().expect("valid NM tag");
    record.data().get(&Tag::from(nm))?.as_int()
}
