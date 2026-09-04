//! Parity gate for the `clip` command's single-threaded-path retirement (C4):
//! `Clip::execute` no longer has a serial in-process loop reached when `--threads`
//! is absent — it *always* routes through the declarative chain builder.
//!
//! Two independent things are checked here:
//!
//! 1. **The cutover actually happened** (`clip_no_threads_routes_through_chain`).
//!    A no-`--threads` run now emits the chain-only `"Using pipeline with N
//!    threads"` banner that `ChainBuilder::add_clip` logs and the retired serial
//!    tail never did. This is the genuine RED/GREEN discriminator: before the
//!    removal a no-`--threads` run took the serial path and printed no such line;
//!    after it, the chain does.
//!
//! 2. **Output equivalence** (`cutover_matches_baseline`), checked one of two ways:
//!    - **CI default (no `FGUMI_BASELINE_BIN`): a self-consistency oracle.** The
//!      chain's output is asserted directly — every record survives, the run
//!      *introduces* clipping over the (deliberately pre-clipped) input, and the
//!      `--metrics` counts are non-zero for the shapes exercised. This is what runs
//!      in ordinary CI, so it must fail on a silently no-op'd clip on its own.
//!    - **When `FGUMI_BASELINE_BIN` is set: an independent byte-parity check.** The
//!      current build's `clip` output — records (byte-identical, modulo the `@PG`
//!      line) and the `--metrics` TSV — must match the frozen pre-cutover baseline
//!      binary run WITHOUT `--threads` (its serial engine). This is the stronger
//!      check, but it runs only when the (host-specific, uncommitted) baseline
//!      binary is provided; it never silently skips (the self-consistency oracle
//!      always runs in its place) — the fallback discipline of
//!      `test_sort_cutover_parity`.
//!
//! So this test proves the cutover lost nothing user-observable only to the strength
//! of whichever half ran: the always-on self-consistency oracle in default CI, and
//! the full byte-parity guarantee when a baseline binary is supplied.
//!
//! The corpus deliberately covers clip's historical known-divergence shapes in one
//! query-grouped stream: a lone fragment, a plain overlapping FR pair, a pair with
//! pre-existing SOFT clipping on R1, a pair with pre-existing HARD clipping on R2,
//! and a template carrying a secondary + a supplementary alignment (passed through
//! but with supplementary mate-info repaired, and clipping upgraded template-wide
//! under `--upgrade-clipping`).
//!
//! **Not a RED/GREEN gate for the equivalence half.** Because the removed serial
//! loop and the chain were already output-equivalent, both checks in (2) pass on
//! both sides of the change — they guard equivalence, not a regression the cutover
//! introduces. The chain-banner check in (1) is the part that flips RED→GREEN.

use std::ffi::OsStr;
use std::path::Path;
use std::process::Command;

use fgumi_raw_bam::{RawRecord, SamBuilder, flags};
use rstest::rstest;
use tempfile::TempDir;

use crate::helpers::bam_generator::{create_minimal_header, create_test_reference, write_bam};
use crate::helpers::read_bam_output;

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
/// Both the current build and the baseline binary stamp a single `@PG` line whose
/// `VN` (git-describe version) and `CL` (command line, naming argv[0] — a
/// different binary path for each side — and the per-run output path) necessarily
/// differ between two independent invocations. Stripping the whole line is correct
/// here because neither side emits more than one `@PG` for these inputs (the test
/// BAMs carry no pre-existing `@PG`).
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
/// binary-layout regression by normalizing it away. clip rewrites SEQ/CIGAR and
/// the NM/UQ/MD and mate tags, so that masking risk is exactly what must be
/// avoided.
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

/// Number of records in one copy of the varied corpus block (see
/// [`build_varied_records`]): a fragment (1) + three pairs (2 each) + a
/// primary-pair-plus-secondary-plus-supplementary template (4) = 11.
const RECORDS_PER_BLOCK: usize = 11;

/// How many times the varied block is repeated. Byte parity is a pure
/// record-for-record comparison, so spanning many `GroupBam` batches is not
/// required here (worker-count independence is covered by
/// `test_clip_command`); a handful of copies keeps the input small and fast.
const BLOCKS: usize = 3;

/// Appends one copy of the varied corpus block to `records`, tagged with `i` so
/// every template QNAME is unique. All reads are 8 bases (`ACGTACGT`) so the
/// reference tag regeneration is well-defined; positions are chosen so the FR
/// pairs overlap (`--clip-overlapping-reads` does real work).
///
/// Shapes covered, in order (kept query-grouped: a template's reads are adjacent):
/// 1. a lone fragment (unpaired primary) — the fragment branch;
/// 2. a plain overlapping FR pair, no pre-existing clipping — overlap + fixed;
/// 3. a pair with pre-existing SOFT clipping on R1 (`2S6M`) — the upgrade path;
/// 4. a pair with pre-existing HARD clipping on R2 (`6M2H`);
/// 5. a primary FR pair plus a supplementary R1 and a secondary R2 —
///    secondary/supplementary passthrough and supplementary mate-info repair.
#[allow(clippy::too_many_lines)] // a straight-line, per-shape record builder reads best in one place
fn push_varied_block(records: &mut Vec<RawRecord>, i: usize) {
    // 1. Lone fragment (unpaired primary).
    {
        let mut b = SamBuilder::new();
        b.read_name(format!("frag_{i:05}").as_bytes())
            .sequence(b"ACGTACGT")
            .qualities(&[30; 8])
            .flags(0)
            .ref_id(0)
            .pos(499)
            .mapq(60)
            .cigar_ops(&[8 << 4]); // 8M
        records.push(b.build());
    }

    // 2. Plain overlapping FR pair, no pre-existing clipping.
    {
        let mut r1 = SamBuilder::new();
        r1.read_name(format!("pair1_{i:05}").as_bytes())
            .sequence(b"ACGTACGT")
            .qualities(&[30; 8])
            .flags(flags::PAIRED | flags::FIRST_SEGMENT | flags::MATE_REVERSE)
            .ref_id(0)
            .pos(99)
            .mapq(60)
            .cigar_ops(&[8 << 4]) // 8M
            .mate_ref_id(0)
            .mate_pos(103)
            .template_length(12);
        records.push(r1.build());

        let mut r2 = SamBuilder::new();
        r2.read_name(format!("pair1_{i:05}").as_bytes())
            .sequence(b"ACGTACGT")
            .qualities(&[30; 8])
            .flags(flags::PAIRED | flags::LAST_SEGMENT | flags::REVERSE)
            .ref_id(0)
            .pos(103)
            .mapq(60)
            .cigar_ops(&[8 << 4]) // 8M
            .mate_ref_id(0)
            .mate_pos(99)
            .template_length(-12);
        records.push(r2.build());
    }

    // 3. Pair with pre-existing SOFT clipping on R1 (2S6M).
    {
        let mut r1 = SamBuilder::new();
        r1.read_name(format!("pair2_{i:05}").as_bytes())
            .sequence(b"ACGTACGT")
            .qualities(&[30; 8])
            .flags(flags::PAIRED | flags::FIRST_SEGMENT | flags::MATE_REVERSE)
            .ref_id(0)
            .pos(199)
            .mapq(60)
            .cigar_ops(&[(2 << 4) | 4, 6 << 4]) // 2S6M
            .mate_ref_id(0)
            .mate_pos(203)
            .template_length(12);
        records.push(r1.build());

        let mut r2 = SamBuilder::new();
        r2.read_name(format!("pair2_{i:05}").as_bytes())
            .sequence(b"ACGTACGT")
            .qualities(&[30; 8])
            .flags(flags::PAIRED | flags::LAST_SEGMENT | flags::REVERSE)
            .ref_id(0)
            .pos(203)
            .mapq(60)
            .cigar_ops(&[8 << 4]) // 8M
            .mate_ref_id(0)
            .mate_pos(199)
            .template_length(-12);
        records.push(r2.build());
    }

    // 4. Pair with pre-existing HARD clipping on R2 (6M2H). Hard-clipped bases are
    // absent from SEQ/QUAL, so R2's sequence is only 6 bases long.
    {
        let mut r1 = SamBuilder::new();
        r1.read_name(format!("pair3_{i:05}").as_bytes())
            .sequence(b"ACGTACGT")
            .qualities(&[30; 8])
            .flags(flags::PAIRED | flags::FIRST_SEGMENT | flags::MATE_REVERSE)
            .ref_id(0)
            .pos(299)
            .mapq(60)
            .cigar_ops(&[8 << 4]) // 8M
            .mate_ref_id(0)
            .mate_pos(303)
            .template_length(10);
        records.push(r1.build());

        let mut r2 = SamBuilder::new();
        r2.read_name(format!("pair3_{i:05}").as_bytes())
            .sequence(b"ACGTAC")
            .qualities(&[30; 6])
            .flags(flags::PAIRED | flags::LAST_SEGMENT | flags::REVERSE)
            .ref_id(0)
            .pos(303)
            .mapq(60)
            .cigar_ops(&[6 << 4, (2 << 4) | 5]) // 6M2H
            .mate_ref_id(0)
            .mate_pos(299)
            .template_length(-10);
        records.push(r2.build());
    }

    // 5. A primary FR pair plus a supplementary R1 and a secondary R2. clip passes
    // secondary/supplementary reads through unclipped but repairs supplementary
    // mate info; with --upgrade-clipping it also upgrades clipping template-wide,
    // including on the supplementary.
    {
        let mut r1 = SamBuilder::new();
        r1.read_name(format!("supp_{i:05}").as_bytes())
            .sequence(b"ACGTACGT")
            .qualities(&[30; 8])
            .flags(flags::PAIRED | flags::FIRST_SEGMENT | flags::MATE_REVERSE)
            .ref_id(0)
            .pos(599)
            .mapq(60)
            .cigar_ops(&[8 << 4]) // 8M
            .mate_ref_id(0)
            .mate_pos(603)
            .template_length(12);
        records.push(r1.build());

        let mut r2 = SamBuilder::new();
        r2.read_name(format!("supp_{i:05}").as_bytes())
            .sequence(b"ACGTACGT")
            .qualities(&[30; 8])
            .flags(flags::PAIRED | flags::LAST_SEGMENT | flags::REVERSE)
            .ref_id(0)
            .pos(603)
            .mapq(60)
            .cigar_ops(&[8 << 4]) // 8M
            .mate_ref_id(0)
            .mate_pos(599)
            .template_length(-12);
        records.push(r2.build());

        // Supplementary R1 (a split alignment): carries a pre-existing soft clip so
        // --upgrade-clipping has something to upgrade on a non-primary read.
        let mut supp_r1 = SamBuilder::new();
        supp_r1
            .read_name(format!("supp_{i:05}").as_bytes())
            .sequence(b"ACGTACGT")
            .qualities(&[30; 8])
            .flags(
                flags::PAIRED | flags::FIRST_SEGMENT | flags::SUPPLEMENTARY | flags::MATE_REVERSE,
            )
            .ref_id(0)
            .pos(2999)
            .mapq(60)
            .cigar_ops(&[6 << 4, (2 << 4) | 4]) // 6M2S (soft clips consume 2 query bases)
            .mate_ref_id(0)
            .mate_pos(603);
        records.push(supp_r1.build());

        // Secondary R2.
        let mut sec_r2 = SamBuilder::new();
        sec_r2
            .read_name(format!("supp_{i:05}").as_bytes())
            .sequence(b"ACGTACGT")
            .qualities(&[30; 8])
            .flags(flags::PAIRED | flags::LAST_SEGMENT | flags::SECONDARY | flags::REVERSE)
            .ref_id(0)
            .pos(3999)
            .mapq(60)
            .cigar_ops(&[8 << 4]); // 8M
        records.push(sec_r2.build());
    }
}

/// Builds the varied query-grouped corpus and returns the records.
fn build_varied_records() -> Vec<RawRecord> {
    let mut records = Vec::with_capacity(RECORDS_PER_BLOCK * BLOCKS);
    for i in 0..BLOCKS {
        push_varied_block(&mut records, i);
    }
    assert_eq!(
        records.len(),
        RECORDS_PER_BLOCK * BLOCKS,
        "corpus size must match RECORDS_PER_BLOCK"
    );
    records
}

/// Writes the varied corpus to `dir/in.bam` (query-grouped header) and returns
/// its path.
fn write_varied_input(dir: &Path) -> std::path::PathBuf {
    let input = dir.join("in.bam");
    write_bam(&input, &create_minimal_header("chr1", 10_000), &build_varied_records());
    input
}

/// A representative clip option set that engages every clip shape in the corpus:
/// overlap clipping on the FR pairs, fixed 5'/3' clipping, and the template-wide
/// `--upgrade-clipping` pre-pass (upgrading the pre-existing soft clips to hard,
/// including on the supplementary alignment). `--clip-bases-past-mate` exercises
/// mate-extension clipping.
const CLIP_OPS: &[&str] = &[
    "--clip-overlapping-reads",
    "--clip-bases-past-mate",
    "--upgrade-clipping",
    "--read-one-five-prime",
    "1",
    "--read-two-three-prime",
    "1",
];

/// Runs `<bin> clip -i <input> -o <output> --reference <ref> [-m <metrics>]
/// <ops...>` (no `--threads`, so the current build takes the post-cutover chain
/// path and the baseline takes its serial path) with `RUST_LOG=info`, and returns
/// the process output for stderr assertions.
fn run_clip(
    bin: &Path,
    input: &Path,
    output: &Path,
    reference: &Path,
    metrics: Option<&Path>,
    ops: &[&str],
) -> std::process::Output {
    let mut cmd = Command::new(bin);
    cmd.env("RUST_LOG", "info").args([
        OsStr::new("clip"),
        OsStr::new("-i"),
        input.as_os_str(),
        OsStr::new("-o"),
        output.as_os_str(),
        OsStr::new("--reference"),
        reference.as_os_str(),
    ]);
    if let Some(m) = metrics {
        cmd.args([OsStr::new("-m"), m.as_os_str()]);
    }
    cmd.args(ops.iter().map(OsStr::new));
    cmd.output().unwrap_or_else(|e| panic!("failed to spawn `{}` clip: {e}", bin.display()))
}

/// A no-`--threads` run now routes through the declarative chain, which logs the
/// `"Using pipeline with N threads"` banner from `ChainBuilder::add_clip`. The
/// retired serial tail logged no such line, so this is the RED (pre-removal) →
/// GREEN (post-removal) discriminator for the cutover.
#[test]
fn clip_no_threads_routes_through_chain() {
    let dir = TempDir::new().expect("temp dir");
    let input = write_varied_input(dir.path());
    let reference = create_test_reference(dir.path());
    let output = dir.path().join("out.bam");

    let current_bin = Path::new(env!("CARGO_BIN_EXE_fgumi"));
    let out = run_clip(current_bin, &input, &output, &reference, None, CLIP_OPS);
    let stderr = String::from_utf8_lossy(&out.stderr);
    assert!(out.status.success(), "no-threads clip must succeed; stderr:\n{stderr}");
    assert!(
        stderr.contains("Using pipeline with 1 threads"),
        "a no-`--threads` clip must route through the chain (which logs the pipeline banner); \
         the serial path is retired. stderr:\n{stderr}"
    );
    // Assert a banner line unique to clip's `add_clip` banner block — `"Clip"`
    // alone also matches "Clipping reads" (the timer) and "Clip overlapping reads",
    // so it would not discriminate the actual banner.
    assert!(
        stderr.contains("Clipping mode:"),
        "the chain must still emit the clip banner (the `Clipping mode:` line); stderr:\n{stderr}"
    );
}

/// Output parity of the post-cutover chain against the pre-removal serial baseline
/// binary, plus the always-available self-consistency oracle when no baseline is
/// set.
///
/// `#[case]` args: a label and whether `--metrics` is written. All cases run the
/// full [`CLIP_OPS`] set over the varied corpus, so every known-divergence shape
/// (fragment, pre-clipped pairs, secondary/supplementary, overlapping FR pairs) is
/// exercised in each.
#[rstest]
#[case::no_metrics(false)]
#[case::with_metrics(true)]
fn cutover_matches_baseline(#[case] with_metrics: bool) {
    let dir = TempDir::new().expect("temp dir");
    let input = write_varied_input(dir.path());
    let reference = create_test_reference(dir.path());

    let current_bin = Path::new(env!("CARGO_BIN_EXE_fgumi"));
    let current_out = dir.path().join("current.bam");
    let current_tsv = dir.path().join("current.tsv");
    let current_metrics = with_metrics.then(|| current_tsv.clone());
    let current = run_clip(
        current_bin,
        &input,
        &current_out,
        &reference,
        current_metrics.as_deref(),
        CLIP_OPS,
    );
    let current_stderr = String::from_utf8_lossy(&current.stderr);
    assert!(current.status.success(), "current clip must succeed; stderr:\n{current_stderr}");

    if let Some(baseline) = baseline_bin() {
        let baseline_out = dir.path().join("baseline.bam");
        let baseline_tsv = dir.path().join("baseline.tsv");
        let baseline_metrics = with_metrics.then(|| baseline_tsv.clone());
        let base = run_clip(
            &baseline,
            &input,
            &baseline_out,
            &reference,
            baseline_metrics.as_deref(),
            CLIP_OPS,
        );
        assert!(
            base.status.success(),
            "baseline clip failed; stderr:\n{}",
            String::from_utf8_lossy(&base.stderr)
        );

        assert_eq!(
            decompressed_records_without_pg(&current_out),
            decompressed_records_without_pg(&baseline_out),
            "chain clip output diverges from the pre-removal serial baseline binary ({}) \
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
            "SKIP baseline half of cutover_matches_baseline[with_metrics={with_metrics}]: \
             FGUMI_BASELINE_BIN is unset or does not name an existing file — running \
             self-consistency oracle instead"
        );
        assert_self_consistent(&current_out, current_metrics.as_deref());
    }
}

/// Sum of soft- + hard-clipped bases in a raw input record's CIGAR.
fn clipped_bases_raw(record: &RawRecord) -> usize {
    record
        .cigar_ops_typed()
        .filter(|op| {
            matches!(
                op.kind(),
                fgumi_raw_bam::CigarKind::SoftClip | fgumi_raw_bam::CigarKind::HardClip
            )
        })
        .map(|op| op.len() as usize)
        .sum()
}

/// Sum of soft- + hard-clipped bases in a parsed output record's CIGAR.
fn clipped_bases_record_buf(record: &noodles::sam::alignment::RecordBuf) -> usize {
    use noodles::sam::alignment::record::cigar::op::Kind;
    record
        .cigar()
        .as_ref()
        .iter()
        .filter(|op| matches!(op.kind(), Kind::SoftClip | Kind::HardClip))
        .map(|op| op.len())
        .sum()
}

/// Always-available oracle used when no baseline binary is set. It must fail on a
/// silently no-op'd / passthrough clip on its own (it is the CI-default check).
///
/// The corpus is *deliberately pre-clipped* (soft/hard pre-clips on some pairs), so
/// "an output CIGAR carries a clip op" would pass even on a passthrough — that is
/// exactly the trap this avoids. Instead it proves THIS run introduced clipping:
/// every record survives (clip never adds/drops records), the total clipped-base
/// footprint of the output strictly exceeds the input's (`CLIP_OPS` adds overlap +
/// fixed 5'/3' clipping), and a record that arrived unclipped (pair1's R1) leaves
/// clipped. With `--metrics` it further asserts the recorded counts are non-zero for
/// the shapes exercised (overlap, fixed 5'/3', and the pre-existing `prior`) — not
/// merely that a Fragment row exists.
fn assert_self_consistent(output: &Path, metrics: Option<&Path>) {
    use fgumi_lib::metrics::{ClippingMetrics, ReadType};

    let input_records = build_varied_records();
    let (_, out_records) = read_bam_output(output);
    assert_eq!(
        out_records.len(),
        input_records.len(),
        "clip must keep every record (it never adds or drops records)"
    );

    // clip preserves input order, so records line up by index. A real run must
    // strictly increase the total clipped-base footprint; a passthrough leaves it
    // equal and fails here. This is the check that makes the oracle non-vacuous.
    let input_clipped: usize = input_records.iter().map(clipped_bases_raw).sum();
    let output_clipped: usize = out_records.iter().map(clipped_bases_record_buf).sum();
    assert!(
        output_clipped > input_clipped,
        "clip must introduce new clipping: output clipped bases ({output_clipped}) must exceed \
         input ({input_clipped}); a no-op / passthrough clip would leave them equal"
    );

    // Pinpoint the overlapping pair: pair1's R1 (index 1) arrives as a plain 8M with
    // no clipping, so after fixed 5' + overlap clipping it must carry clip ops. This
    // catches a regression that clipped some records but silently skipped the FR pair.
    assert_eq!(clipped_bases_raw(&input_records[1]), 0, "precondition: pair1 R1 starts unclipped");
    assert!(
        clipped_bases_record_buf(&out_records[1]) > 0,
        "pair1 R1 must be clipped by --read-one-five-prime + --clip-overlapping-reads"
    );

    if let Some(path) = metrics {
        let rows = fgumi_lib::metrics::read_metrics::<_, ClippingMetrics>(path, "clipping")
            .expect("parse clipping metrics");
        let by_type = |t: ReadType| {
            rows.iter()
                .find(|m| m.read_type == t)
                .unwrap_or_else(|| panic!("metrics TSV must carry a {t:?} row"))
        };
        assert_eq!(by_type(ReadType::Fragment).reads, BLOCKS, "one fragment read per corpus block");
        // A real run records non-zero clipping for the shapes CLIP_OPS exercises,
        // not an all-zero table.
        assert!(
            by_type(ReadType::Pair).bases_clipped_overlapping > 0,
            "overlap clipping must be recorded for the FR pairs (rolled up into Pair)"
        );
        assert!(
            by_type(ReadType::ReadOne).bases_clipped_five_prime > 0
                || by_type(ReadType::ReadTwo).bases_clipped_five_prime > 0,
            "--read-one-five-prime must record 5' clipping"
        );
        assert!(
            by_type(ReadType::ReadOne).bases_clipped_three_prime > 0
                || by_type(ReadType::ReadTwo).bases_clipped_three_prime > 0,
            "--read-two-three-prime must record 3' clipping"
        );
        assert!(
            by_type(ReadType::ReadOne).bases_clipped_pre > 0
                || by_type(ReadType::ReadTwo).bases_clipped_pre > 0,
            "pre-existing soft/hard pre-clips must be recorded as `prior` (bases_clipped_pre)"
        );
    }
}
