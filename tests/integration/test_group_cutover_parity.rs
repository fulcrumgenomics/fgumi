//! Cutover parity gate for the `group` command's legacy-path retirement (C4):
//! `GroupReadsByUmi::execute` no longer has a no-`--threads` hand-rolled
//! streaming engine (`execute_single_threaded`) — it *always* routes through
//! the declarative chain builder (`execute_chain`). This test proves that
//! cutover lost nothing user-observable in the UMI-family-assignment logic.
//!
//! Two independent things are checked here:
//!
//! 1. **The cutover actually happened** (`group_no_threads_routes_through_chain`).
//!    A no-`--threads` run now emits the chain-only "Pipeline `group` ran in"
//!    line from `pipeline::chains::finalize`, and the per-step `Pipeline stats`
//!    table — both come from the pipeline-core runtime, which the retired
//!    single-threaded engine never used. This is the RED (pre-removal) / GREEN
//!    (post-removal) discriminator for the cutover.
//!
//! 2. **Output parity with the pre-removal legacy path**
//!    (`cutover_matches_baseline_*`). The current build's `group` output — BAM
//!    records (byte-identical, modulo the `@PG` line) and `--metrics`/
//!    `--family-size-histogram`/`--grouping-metrics` — must match the frozen
//!    legacy-path baseline binary run with no `--threads` (its hand-rolled
//!    `execute_single_threaded` engine). The baseline path comes from
//!    `FGUMI_BASELINE_BIN`; when unset (or names a missing file) the case
//!    degrades to a self-consistency oracle rather than skipping — the
//!    fallback discipline of `test_sort_cutover_parity.rs`.
//!
//! Coverage spans the shapes the campaign plan calls out as historically
//! divergence-prone between the legacy engine and the chain: the assigner
//! strategies `group` supports (identity / edit / adjacency / paired, the last
//! exercising the `/A`+`/B` duplex-style molecule-id split), paired-end vs
//! single-end (fragment) templates in the same input, the
//! `--index-threshold` indexed-clustering path, and secondary/supplementary
//! reads (which `group`, unlike `dedup`, discards entirely — see #903).

use std::path::Path;
use std::process::Command;

use rstest::rstest;
use tempfile::TempDir;

use fgumi_lib::sam::SamTag;
use fgumi_raw_bam::{RawRecord, SamBuilder, flags};

use crate::helpers::bam_generator::{create_minimal_header, write_bam};
use crate::helpers::read_bam_output;

/// Resolves the saved pre-removal legacy-path baseline binary to compare
/// against.
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
/// The fixtures are built by piping through `fgumi sort` (see
/// `create_sorted_bam`), which stamps its own `@PG` line, so the group input
/// already carries one before `group` adds a second (`PP`-chained) one. The
/// `sort` line is byte-identical between the current build's run and the
/// baseline's — both invoke the *same* `fgumi sort` binary on the same input —
/// so it would compare equal either way; it is the `group`/baseline-binary
/// `@PG` line that differs (`VN`, the git-describe version, and `CL`, which
/// names argv[0] — a different binary path for each side — and the per-run
/// output path). Stripping every `@PG` line rather than only the last one is
/// simply the simpler correct implementation: dropping an already-identical
/// line changes nothing.
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
/// binary-layout regression (e.g. MI-tag injection) by normalizing it away.
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

/// Writes `records` to an unsorted BAM alongside `path`, then runs `fgumi
/// sort --order template-coordinate` to produce the final input at `path`.
///
/// `group` requires *true* template-coordinate order for
/// `RecordPositionGrouper`'s streaming position-group boundaries to line up
/// correctly; writing records in an order that merely satisfies the header's
/// `SO`/`GO`/`SS` tags is not sufficient. Routing through a real `fgumi sort`
/// is the same approach `test_dedup_command.rs::create_sorted_bam_with_header`
/// uses for exactly this reason.
fn create_sorted_bam(path: &Path, header: &noodles::sam::Header, records: &[RawRecord]) {
    let unsorted = path.with_extension("unsorted.bam");
    write_bam(&unsorted, header, records);

    let status = Command::new(env!("CARGO_BIN_EXE_fgumi"))
        .args([
            "sort",
            "-i",
            unsorted.to_str().unwrap(),
            "-o",
            path.to_str().unwrap(),
            "--order",
            "template-coordinate",
            "--key-types",
            "mi",
        ])
        .status()
        .unwrap_or_else(|e| panic!("failed to spawn `fgumi sort` for test input: {e}"));
    assert!(status.success(), "failed to template-coordinate sort test input {}", path.display());
}

//////////////////////////////////////////////////////////////////////////////
// Fixture builders
//////////////////////////////////////////////////////////////////////////////

/// Builds one paired-end template: R1 at `r1_pos` and R2 at `r2_pos` (1-based
/// reference positions), with the given `umi` and strand orientation. Carries
/// `MC` (required by `RecordPositionGrouper` for paired reads).
fn build_pair(
    name: &str,
    umi: &str,
    r1_pos: i32,
    r2_pos: i32,
    r1_reverse: bool,
    base_qual: u8,
) -> (RawRecord, RawRecord) {
    let seq = b"ACGTACGT";
    let cigar_op = (u32::try_from(seq.len()).expect("fits u32")) << 4;
    let r2_reverse = !r1_reverse;
    let tlen = (r2_pos - r1_pos) + i32::try_from(seq.len()).expect("fits i32");

    let mut b1 = SamBuilder::new();
    b1.read_name(name.as_bytes())
        .sequence(seq)
        .qualities(&vec![base_qual; seq.len()])
        .flags(
            flags::PAIRED
                | flags::FIRST_SEGMENT
                | if r1_reverse { flags::REVERSE } else { 0 }
                | if r2_reverse { flags::MATE_REVERSE } else { 0 },
        )
        .ref_id(0)
        .pos(r1_pos - 1)
        .mapq(60)
        .cigar_ops(&[cigar_op])
        .mate_ref_id(0)
        .mate_pos(r2_pos - 1)
        .template_length(tlen)
        .add_string_tag(SamTag::MC, b"8M")
        .add_string_tag(SamTag::RX, umi.as_bytes());

    let mut b2 = SamBuilder::new();
    b2.read_name(name.as_bytes())
        .sequence(seq)
        .qualities(&vec![base_qual; seq.len()])
        .flags(
            flags::PAIRED
                | flags::LAST_SEGMENT
                | if r2_reverse { flags::REVERSE } else { 0 }
                | if r1_reverse { flags::MATE_REVERSE } else { 0 },
        )
        .ref_id(0)
        .pos(r2_pos - 1)
        .mapq(60)
        .cigar_ops(&[cigar_op])
        .mate_ref_id(0)
        .mate_pos(r1_pos - 1)
        .template_length(-tlen)
        .add_string_tag(SamTag::MC, b"8M")
        .add_string_tag(SamTag::RX, umi.as_bytes());

    (b1.build(), b2.build())
}

/// Builds `count` paired templates, all at the same coordinates and UMI.
fn build_family(
    base_name: &str,
    umi: &str,
    count: usize,
    r1_pos: i32,
    r2_pos: i32,
) -> Vec<RawRecord> {
    let mut records = Vec::with_capacity(count * 2);
    for i in 0..count {
        let (r1, r2) = build_pair(&format!("{base_name}_{i}"), umi, r1_pos, r2_pos, false, 30);
        records.push(r1);
        records.push(r2);
    }
    records
}

/// A single-end (unpaired, mapped) read at `pos` with `umi`.
fn build_single_end(name: &str, umi: &str, pos: i32) -> RawRecord {
    let seq = b"ACGTACGT";
    let cigar_op = (u32::try_from(seq.len()).expect("fits u32")) << 4;
    let mut b = SamBuilder::new();
    b.read_name(name.as_bytes())
        .sequence(seq)
        .qualities(&[30; 8])
        .flags(0)
        .ref_id(0)
        .pos(pos - 1)
        .mapq(60)
        .cigar_ops(&[cigar_op])
        .add_string_tag(SamTag::RX, umi.as_bytes());
    b.build()
}

/// A secondary or supplementary read at the same coordinates as a primary
/// pair. `group` (unlike `dedup`) discards secondary/supplementary reads
/// entirely (#903), so unlike `dedup`'s cutover-parity fixture this needs no
/// `tc` tag — it never reaches template building.
fn build_non_primary_at(
    name: &str,
    umi: &str,
    pos: i32,
    mate_pos: i32,
    supplementary: bool,
) -> (RawRecord, RawRecord) {
    let seq = b"ACGTACGT";
    let cigar_op = (u32::try_from(seq.len()).expect("fits u32")) << 4;
    let non_primary_flag = if supplementary { flags::SUPPLEMENTARY } else { flags::SECONDARY };
    let mut b1 = SamBuilder::new();
    b1.read_name(name.as_bytes())
        .sequence(seq)
        .qualities(&[30; 8])
        .flags(flags::PAIRED | flags::FIRST_SEGMENT | flags::MATE_REVERSE | non_primary_flag)
        .ref_id(0)
        .pos(pos - 1)
        .mapq(60)
        .cigar_ops(&[cigar_op])
        .mate_ref_id(0)
        .mate_pos(mate_pos - 1)
        .add_string_tag(SamTag::RX, umi.as_bytes());
    let mut b2 = SamBuilder::new();
    b2.read_name(name.as_bytes())
        .sequence(seq)
        .qualities(&[30; 8])
        .flags(flags::PAIRED | flags::LAST_SEGMENT | flags::REVERSE | non_primary_flag)
        .ref_id(0)
        .pos(mate_pos - 1)
        .mapq(60)
        .cigar_ops(&[cigar_op])
        .mate_ref_id(0)
        .mate_pos(pos - 1)
        .add_string_tag(SamTag::RX, umi.as_bytes());
    (b1.build(), b2.build())
}

/// The main fixture: three independent position groups exercising, in one
/// input, everything a single strategy run needs to prove parity on:
///
/// - `(100, 300)`: a 3-read exact-UMI family (`AAAAAAAA`) plus a 2-read
///   1-mismatch family (`AAAAAAAT`) — under `--strategy identity` these are two
///   separate molecules; under `edit`/`adjacency` they cluster into one. Also
///   carries one secondary and one supplementary pair, which must be dropped
///   entirely from the output regardless of strategy.
/// - `(500, 700)`: a plain 2-read family (`CCCCCCCC`), an ordinary paired
///   control group.
/// - single-end reads at position `900`: 3 unpaired reads sharing UMI
///   `GGGGGGGG` — the fragment shape.
fn build_family_fixture(dir: &Path) -> std::path::PathBuf {
    let input = dir.join("in.bam");
    let header = create_minimal_header("chr1", 10000);

    let mut records = Vec::new();
    records.extend(build_family("exact", "AAAAAAAA", 3, 100, 300));
    records.extend(build_family("near", "AAAAAAAT", 2, 100, 300));
    let (s1, s2) = build_non_primary_at("exact_0_sec", "AAAAAAAA", 100, 300, false);
    records.push(s1);
    records.push(s2);
    let (u1, u2) = build_non_primary_at("exact_0_sup", "AAAAAAAA", 100, 300, true);
    records.push(u1);
    records.push(u2);

    records.extend(build_family("ctrl", "CCCCCCCC", 2, 500, 700));

    for i in 0..3 {
        records.push(build_single_end(&format!("frag_{i}"), "GGGGGGGG", 900));
    }

    create_sorted_bam(&input, &header, &records);
    input
}

/// A duplex-shaped fixture for `--strategy paired`: two templates at the same
/// coordinates carrying the paired-UMI format (`AAAA-CCCC`) in opposite
/// physical orientations (FR vs RF) — the shape that makes the `PairedUmiAssigner`
/// assign the same molecule two distinct `/A`/`/B` suffixes.
fn build_paired_strategy_fixture(dir: &Path) -> std::path::PathBuf {
    let input = dir.join("in.bam");
    let header = create_minimal_header("chr1", 10000);

    let mut records = Vec::new();
    let (r1, r2) = build_pair("duplex_fwd", "AAAA-CCCC", 100, 300, false, 30);
    records.push(r1);
    records.push(r2);
    let (r1, r2) = build_pair("duplex_rev", "AAAA-CCCC", 300, 100, true, 30);
    records.push(r1);
    records.push(r2);

    create_sorted_bam(&input, &header, &records);
    input
}

/// A single large position group of `n` templates with pairwise-distinct UMIs,
/// large enough to cross the default `--index-threshold 100` (adjacency/edit
/// index at `--edits 1`) so both the current build and the baseline exercise
/// the indexed-clustering path.
fn build_index_threshold_fixture(dir: &Path, n: usize) -> std::path::PathBuf {
    let input = dir.join("in.bam");
    let header = create_minimal_header("chr1", 10000);

    let bases = [b'A', b'C', b'G', b'T'];
    let mut records = Vec::with_capacity(n * 2);
    for i in 0..n {
        // Single-parity-check code: the first 7 symbols are the base-4 digits
        // of `i` (4^7 = 16384, far more than `n`), the 8th is their digit-sum
        // mod 4. Flipping exactly one base always changes that sum mod 4, so
        // any two *distinct* codes from this family differ in at least 2
        // positions — never edit distance 1 of each other. Plain sequential
        // base-4 digits (no parity digit) do NOT have this property:
        // consecutive indices differ only in their lowest digit, chaining
        // almost every UMI to a neighbor at edit distance 1 and letting
        // adjacency's directed capture collapse the whole set into one
        // molecule — not the "120 independent molecules" this fixture needs
        // to isolate the indexed-clustering path from UMI-clustering
        // behavior.
        let mut umi = [b'A'; 8];
        let mut v = i;
        let mut parity = 0usize;
        for slot in &mut umi[..7] {
            let digit = v % 4;
            v /= 4;
            parity = (parity + digit) % 4;
            *slot = bases[digit];
        }
        umi[7] = bases[parity];
        let umi = std::str::from_utf8(&umi).expect("ACGT is valid UTF-8");
        let (r1, r2) = build_pair(&format!("idx_{i}"), umi, 100, 300, false, 30);
        records.push(r1);
        records.push(r2);
    }

    create_sorted_bam(&input, &header, &records);
    input
}

//////////////////////////////////////////////////////////////////////////////
// Runner
//////////////////////////////////////////////////////////////////////////////

/// Runs `<bin> group` (no `--threads`, so the current build takes the
/// post-cutover chain path and the baseline takes its legacy engine) with
/// `RUST_LOG=info`, and returns the process output for stderr/exit assertions.
#[allow(clippy::too_many_arguments)]
fn run_group(
    bin: &Path,
    input: &Path,
    output: &Path,
    strategy: &str,
    metrics_prefix: Option<&Path>,
    index_threshold: Option<&str>,
) -> std::process::Output {
    let mut cmd = Command::new(bin);
    cmd.env("RUST_LOG", "info").args([
        "group",
        "-i",
        input.to_str().unwrap(),
        "-o",
        output.to_str().unwrap(),
        "--strategy",
        strategy,
        "--compression-level",
        "1",
    ]);
    if let Some(m) = metrics_prefix {
        cmd.args(["-M", m.to_str().unwrap()]);
    }
    if let Some(t) = index_threshold {
        cmd.args(["--index-threshold", t]);
    }
    cmd.output().unwrap_or_else(|e| panic!("failed to spawn `{}` group: {e}", bin.display()))
}

/// The three files a `--metrics PREFIX` run writes (see
/// `commands::group::with_extension`).
fn metrics_prefix_paths(
    prefix: &Path,
) -> (std::path::PathBuf, std::path::PathBuf, std::path::PathBuf) {
    let ext = |suffix: &str| {
        let mut s = prefix.as_os_str().to_owned();
        s.push(".");
        s.push(suffix);
        std::path::PathBuf::from(s)
    };
    (ext("family_sizes.txt"), ext("grouping_metrics.txt"), ext("position_group_sizes.txt"))
}

//////////////////////////////////////////////////////////////////////////////
// Tests
//////////////////////////////////////////////////////////////////////////////

/// A no-`--threads` run now routes through the declarative chain, which logs
/// "Pipeline `group` ran in" (from `pipeline::chains::finalize`) and the
/// per-step `Pipeline stats` table (from the pipeline-core runtime). The
/// retired single-threaded engine never printed either — a hand-rolled
/// streaming loop with no step-based runtime — so this is the RED
/// (pre-removal) -> GREEN (post-removal) discriminator for the cutover.
#[test]
fn group_no_threads_routes_through_chain() {
    let dir = TempDir::new().expect("temp dir");
    let input = build_family_fixture(dir.path());
    let output = dir.path().join("out.bam");

    let current_bin = Path::new(env!("CARGO_BIN_EXE_fgumi"));
    let out = run_group(current_bin, &input, &output, "identity", None, None);
    let stderr = String::from_utf8_lossy(&out.stderr);
    assert!(out.status.success(), "no-threads group must succeed; stderr:\n{stderr}");
    assert!(
        stderr.contains("Pipeline `group` ran in"),
        "a no-`--threads` group must route through the chain (which logs the pipeline-core \
         finalize line); the single-threaded engine is retired. stderr:\n{stderr}"
    );
    assert!(
        stderr.contains("Starting group"),
        "the chain must still emit the `Starting group` banner; stderr:\n{stderr}"
    );
}

/// Output parity of the post-cutover chain against the pre-removal legacy-path
/// baseline binary (run with no `--threads`, i.e. its single-threaded engine),
/// across the assigner strategies `group` supports — plus the always-available
/// self-consistency oracle when no baseline is set.
///
/// `identity` keeps the exact/near-mismatch families in [`build_family_fixture`]
/// as two distinct molecules; `edit` and `adjacency` cluster them into one. All
/// three exercise the secondary/supplementary-discard and single-end (fragment)
/// position groups identically, since strategy only changes UMI clustering.
#[rstest]
#[case::identity("identity")]
#[case::edit("edit")]
#[case::adjacency("adjacency")]
fn cutover_matches_baseline_by_strategy(#[case] strategy: &str) {
    let dir = TempDir::new().expect("temp dir");
    let input = build_family_fixture(dir.path());

    let current_bin = Path::new(env!("CARGO_BIN_EXE_fgumi"));
    let current_out = dir.path().join("current.bam");
    let current_prefix = dir.path().join("current");
    let current =
        run_group(current_bin, &input, &current_out, strategy, Some(&current_prefix), None);
    let current_stderr = String::from_utf8_lossy(&current.stderr);
    assert!(current.status.success(), "current group must succeed; stderr:\n{current_stderr}");

    if let Some(baseline) = baseline_bin() {
        let baseline_out = dir.path().join("baseline.bam");
        let baseline_prefix = dir.path().join("baseline");
        let base =
            run_group(&baseline, &input, &baseline_out, strategy, Some(&baseline_prefix), None);
        assert!(
            base.status.success(),
            "baseline group failed; stderr:\n{}",
            String::from_utf8_lossy(&base.stderr)
        );

        assert_eq!(
            decompressed_records_without_pg(&current_out),
            decompressed_records_without_pg(&baseline_out),
            "chain group output (--strategy {strategy}) diverges from the pre-removal legacy \
             baseline binary ({}) after stripping @PG — a real cutover parity bug",
            baseline.display(),
        );

        let (cur_fam, cur_grp, cur_pos) = metrics_prefix_paths(&current_prefix);
        let (base_fam, base_grp, base_pos) = metrics_prefix_paths(&baseline_prefix);
        for (label, cur, base) in [
            ("family_sizes", cur_fam, base_fam),
            ("grouping_metrics", cur_grp, base_grp),
            ("position_group_sizes", cur_pos, base_pos),
        ] {
            assert_eq!(
                std::fs::read_to_string(&cur).unwrap_or_else(|e| panic!("current {label}: {e}")),
                std::fs::read_to_string(&base).unwrap_or_else(|e| panic!("baseline {label}: {e}")),
                "chain --metrics {label} TSV (--strategy {strategy}) diverges from the legacy \
                 baseline"
            );
        }
    } else {
        eprintln!(
            "SKIP baseline half of cutover_matches_baseline_by_strategy[{strategy}]: \
             FGUMI_BASELINE_BIN is unset or does not name an existing file — running \
             self-consistency oracle instead"
        );
        assert_self_consistent_family_fixture(&current_out, strategy);
    }
}

/// Output parity for `--strategy paired` on the duplex-shaped (`/A`/`/B`)
/// fixture — the known-divergence shape where a paired UMI read from both
/// physical strands must resolve to one canonical molecule split into two
/// duplex halves.
#[test]
fn cutover_matches_baseline_paired_duplex() {
    let dir = TempDir::new().expect("temp dir");
    let input = build_paired_strategy_fixture(dir.path());

    let current_bin = Path::new(env!("CARGO_BIN_EXE_fgumi"));
    let current_out = dir.path().join("current.bam");
    let current = run_group(current_bin, &input, &current_out, "paired", None, None);
    let current_stderr = String::from_utf8_lossy(&current.stderr);
    assert!(
        current.status.success(),
        "current group (paired) must succeed; stderr:\n{current_stderr}"
    );

    if let Some(baseline) = baseline_bin() {
        let baseline_out = dir.path().join("baseline.bam");
        let base = run_group(&baseline, &input, &baseline_out, "paired", None, None);
        assert!(
            base.status.success(),
            "baseline group (paired) failed; stderr:\n{}",
            String::from_utf8_lossy(&base.stderr)
        );
        assert_eq!(
            decompressed_records_without_pg(&current_out),
            decompressed_records_without_pg(&baseline_out),
            "chain group output (--strategy paired, duplex A/B shape) diverges from the \
             pre-removal legacy baseline binary ({}) after stripping @PG",
            baseline.display(),
        );
    } else {
        eprintln!(
            "SKIP baseline half of cutover_matches_baseline_paired_duplex: FGUMI_BASELINE_BIN \
             is unset or does not name an existing file — running self-consistency oracle instead"
        );
        let (_, records) = read_bam_output(&current_out);
        assert_eq!(records.len(), 4, "both duplex-strand templates (2 reads each) must survive");
        let mi_tags = mi_strings(&records);
        assert_eq!(mi_tags.len(), 4, "every record must carry an MI tag");
        let distinct: std::collections::HashSet<&String> = mi_tags.iter().collect();
        assert_eq!(
            distinct.len(),
            2,
            "the two physical-strand templates must resolve to exactly 2 distinct molecule ids \
             (the /A and /B duplex halves), got: {mi_tags:?}"
        );
        assert!(
            mi_tags.iter().any(|m| m.ends_with("/A")) && mi_tags.iter().any(|m| m.ends_with("/B")),
            "the paired assigner must emit both a /A and a /B suffixed molecule id, got: {mi_tags:?}"
        );
    }
}

/// Output parity across the `--index-threshold`-triggered indexed-clustering
/// path: a single position group of 120 pairwise-distinct UMIs under
/// `--strategy adjacency`, which crosses the default `--index-threshold 100`
/// and forces both binaries onto the N-gram/BK-tree index rather than the
/// brute-force pairwise comparison.
#[test]
fn cutover_matches_baseline_index_threshold() {
    let dir = TempDir::new().expect("temp dir");
    let input = build_index_threshold_fixture(dir.path(), 120);

    let current_bin = Path::new(env!("CARGO_BIN_EXE_fgumi"));
    let current_out = dir.path().join("current.bam");
    let current = run_group(current_bin, &input, &current_out, "adjacency", None, Some("always"));
    let current_stderr = String::from_utf8_lossy(&current.stderr);
    assert!(
        current.status.success(),
        "current group (index-threshold always) must succeed; stderr:\n{current_stderr}"
    );

    if let Some(baseline) = baseline_bin() {
        let baseline_out = dir.path().join("baseline.bam");
        let base = run_group(&baseline, &input, &baseline_out, "adjacency", None, Some("always"));
        assert!(
            base.status.success(),
            "baseline group (index-threshold always) failed; stderr:\n{}",
            String::from_utf8_lossy(&base.stderr)
        );
        assert_eq!(
            decompressed_records_without_pg(&current_out),
            decompressed_records_without_pg(&baseline_out),
            "chain group output (120-UMI indexed adjacency) diverges from the pre-removal \
             legacy baseline binary ({}) after stripping @PG",
            baseline.display(),
        );
    } else {
        eprintln!(
            "SKIP baseline half of cutover_matches_baseline_index_threshold: FGUMI_BASELINE_BIN \
             is unset or does not name an existing file — running self-consistency oracle instead"
        );
        let (_, records) = read_bam_output(&current_out);
        assert_eq!(records.len(), 240, "all 120 templates (2 reads each) must survive");
        let mi_tags = mi_strings(&records);
        let distinct: std::collections::HashSet<&String> = mi_tags.iter().collect();
        assert_eq!(
            distinct.len(),
            120,
            "120 pairwise-distinct UMIs must resolve to 120 distinct molecule ids"
        );
    }
}

//////////////////////////////////////////////////////////////////////////////
// Self-consistency oracle
//////////////////////////////////////////////////////////////////////////////

/// Extracts every record's `MI:Z` tag value (skips records without one).
fn mi_strings(records: &[noodles::sam::alignment::RecordBuf]) -> Vec<String> {
    use noodles::sam::alignment::record::data::field::Tag;
    use noodles::sam::alignment::record_buf::data::field::Value;

    let mi_tag = Tag::from(SamTag::MI);
    records
        .iter()
        .filter_map(|r| match r.data().get(&mi_tag) {
            Some(Value::String(s)) => Some(s.to_string()),
            _ => None,
        })
        .collect()
}

fn mi_of(record: &noodles::sam::alignment::RecordBuf) -> Option<String> {
    use noodles::sam::alignment::record::data::field::Tag;
    use noodles::sam::alignment::record_buf::data::field::Value;

    match record.data().get(&Tag::from(SamTag::MI)) {
        Some(Value::String(s)) => Some(s.to_string()),
        _ => None,
    }
}

/// Always-available oracle used when no baseline binary is set: proves the
/// chain's group output on [`build_family_fixture`] is a real transform, not a
/// passthrough — every kept record carries an `MI` tag, the exact/near-UMI
/// families come out as the strategy-appropriate number of distinct molecules,
/// secondary/supplementary reads are entirely absent from the output, and
/// single-end fragments form their own molecule.
///
/// A no-op run (output MI-family structure identical across an input that
/// should cluster, or secondary/supplementary reads leaking through) fails
/// this — see the per-assertion messages.
fn assert_self_consistent_family_fixture(output: &Path, strategy: &str) {
    let (_, records) = read_bam_output(output);
    // exact(3 pairs=6) + near(2 pairs=4) + ctrl(2 pairs=4) + frag(3) = 17; the
    // 2 secondary + 2 supplementary reads are discarded entirely (#903).
    assert_eq!(
        records.len(),
        17,
        "secondary/supplementary reads must be discarded and every primary kept"
    );
    assert!(
        records.iter().all(|r| !r.flags().is_secondary() && !r.flags().is_supplementary()),
        "no output record may be flagged secondary or supplementary"
    );

    let mi_tags = mi_strings(&records);
    assert_eq!(mi_tags.len(), 17, "every kept record must carry an MI tag");

    // The exact(3)+near(2) position group: under identity, 2 distinct
    // molecules; under edit/adjacency, they cluster into 1.
    let exact_near_mi: std::collections::HashSet<String> = records
        .iter()
        .filter(|r| {
            r.name().is_some_and(|n| {
                let n = String::from_utf8_lossy(n.as_ref());
                n.starts_with("exact_") || n.starts_with("near_")
            })
        })
        .filter_map(mi_of)
        .collect();
    let expected_molecules = if strategy == "identity" { 2 } else { 1 };
    assert_eq!(
        exact_near_mi.len(),
        expected_molecules,
        "--strategy {strategy}: exact(AAAAAAAA)+near(AAAAAAAT) families must resolve to \
         {expected_molecules} distinct molecule id(s), got {exact_near_mi:?}"
    );

    // Control family (2 templates at 500/700): exactly 1 molecule.
    let ctrl_mi: std::collections::HashSet<String> = records
        .iter()
        .filter(|r| {
            r.name().is_some_and(|n| String::from_utf8_lossy(n.as_ref()).starts_with("ctrl_"))
        })
        .filter_map(mi_of)
        .collect();
    assert_eq!(ctrl_mi.len(), 1, "the 2-template control family must resolve to 1 molecule");

    // Single-end fragments (3 unpaired reads, same UMI, same position): one
    // molecule, and none carry the PAIRED flag.
    let frag_records: Vec<_> = records
        .iter()
        .filter(|r| {
            r.name().is_some_and(|n| String::from_utf8_lossy(n.as_ref()).starts_with("frag_"))
        })
        .collect();
    assert_eq!(frag_records.len(), 3, "all 3 single-end fragments must survive");
    assert!(
        frag_records.iter().all(|r| !r.flags().is_segmented()),
        "single-end fragments must not carry the PAIRED flag"
    );
    let frag_mi: std::collections::HashSet<String> =
        frag_records.iter().filter_map(|r| mi_of(r)).collect();
    assert_eq!(frag_mi.len(), 1, "the 3 identical single-end fragments must resolve to 1 molecule");
}
