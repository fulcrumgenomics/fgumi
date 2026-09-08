//! Cutover parity gate for the `dedup` command's legacy-path retirement (C4):
//! `MarkDuplicates::execute` no longer has a no-`--threads` hand-rolled
//! `unified_pipeline` engine — it *always* routes through the declarative chain
//! builder (`execute_chain`). This test proves that cutover lost nothing
//! user-observable in the duplicate-marking / UMI-family-assignment logic.
//!
//! Two independent things are checked here:
//!
//! 1. **The cutover actually happened** (`dedup_no_threads_routes_through_chain`).
//!    A no-`--threads` run now emits the chain-only "Pipeline `dedup` ran in"
//!    line from `pipeline::chains::finalize`, and the per-step `Pipeline stats`
//!    table — both come from the pipeline-core runtime, which the retired
//!    `unified_pipeline` engine never used. This is the RED (pre-removal) /
//!    GREEN (post-removal) discriminator for the cutover.
//!
//! 2. **Output parity with the pre-removal legacy path**
//!    (`cutover_matches_baseline_*`). The current build's `dedup` output — BAM
//!    records (byte-identical, modulo the `@PG` line) and `--metrics` /
//!    `--family-size-histogram` — must match the frozen legacy-path baseline
//!    binary run with no `--threads` (its hand-rolled `unified_pipeline`
//!    engine). The baseline path comes from `FGUMI_BASELINE_BIN`; when unset
//!    (or names a missing file) the case degrades to a self-consistency oracle
//!    rather than skipping — the fallback discipline of
//!    `test_sort_cutover_parity.rs`.
//!
//! Coverage spans the shapes the campaign plan calls out as historically
//! divergence-prone between the legacy engine and the chain: multiple UMI
//! strategies (identity / adjacency / paired, the last exercising the `/A`+`/B`
//! duplex-style molecule-id split), the `--index-threshold` indexed-clustering
//! path (a position group large enough to trigger indexing under the default
//! threshold), single-end (unpaired) templates, and secondary/supplementary
//! reads carrying the `tc` tag.

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
/// `create_sorted_bam`), which stamps its own `@PG` line, so the dedup input
/// already carries one before `dedup` adds a second (`PP`-chained) one. The
/// `sort` line is byte-identical between the current build's run and the
/// baseline's — both invoke the *same* `fgumi sort` binary on the same input —
/// so it would compare equal either way; it is the `dedup`/baseline-binary
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
/// `dedup` requires *true* template-coordinate order, which correctly
/// interleaves `tc`-keyed secondary/supplementary reads at their primary's
/// coordinate regardless of where they were written in the input stream.
/// Writing records in an order that merely satisfies the header's `SO`/`GO`/
/// `SS` tags is not sufficient — an out-of-place secondary/supplementary
/// silently lands in the wrong (or no) position group and is dropped by the
/// `NoPrimaryReads` filter with no error. Routing through a real `fgumi sort`
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

/// Builds `count` paired templates, all at the same coordinates and UMI —
/// a PCR-duplicate family that a duplicate-marking run must collapse to one
/// representative, with the rest flagged `0x400`.
fn build_duplicate_family(
    base_name: &str,
    umi: &str,
    count: usize,
    r1_pos: i32,
    r2_pos: i32,
) -> Vec<RawRecord> {
    let mut records = Vec::with_capacity(count * 2);
    for i in 0..count {
        // Distinct base qualities make the highest-scoring representative
        // unambiguous (no tie whose resolution depends on stream order).
        let qual = 20 + u8::try_from(i).expect("small count") * 5;
        let (r1, r2) = build_pair(&format!("{base_name}_{i}"), umi, r1_pos, r2_pos, false, qual);
        records.push(r1);
        records.push(r2);
    }
    records
}

/// A single-end (unpaired, mapped) read at `pos` with `umi`.
fn build_single_end(name: &str, umi: &str, pos: i32, base_qual: u8) -> RawRecord {
    let seq = b"ACGTACGT";
    let cigar_op = (u32::try_from(seq.len()).expect("fits u32")) << 4;
    let mut b = SamBuilder::new();
    b.read_name(name.as_bytes())
        .sequence(seq)
        .qualities(&vec![base_qual; seq.len()])
        .flags(0)
        .ref_id(0)
        .pos(pos - 1)
        .mapq(60)
        .cigar_ops(&[cigar_op])
        .add_string_tag(SamTag::RX, umi.as_bytes());
    b.build()
}

/// The `tc` aux array `fgumi zipper` stamps onto a secondary/supplementary
/// read: `[tid1, pos1, neg1, tid2, pos2, neg2]`, the *primary* pair's own
/// unclipped-5' reference coordinates and strand. Computed from the actual
/// built R1/R2 records (rather than hand-derived arithmetic) so it can never
/// drift from what `build_pair` actually produced — e.g. a reverse-strand
/// mate's unclipped 5' position is its rightmost aligned base, not `pos`.
fn primary_tc_array(r1: &RawRecord, r2: &RawRecord) -> [i32; 6] {
    use fgumi_raw_bam::RawRecordView;

    let strand_of = |r: &RawRecord| -> i32 {
        i32::from((RawRecordView::new(r.as_ref()).flags() & flags::REVERSE) != 0)
    };
    [
        fgumi_raw_bam::ref_id(r1.as_ref()),
        fgumi_raw_bam::unclipped_5prime_from_raw_bam(r1.as_ref()),
        strand_of(r1),
        fgumi_raw_bam::ref_id(r2.as_ref()),
        fgumi_raw_bam::unclipped_5prime_from_raw_bam(r2.as_ref()),
        strand_of(r2),
    ]
}

/// A secondary or supplementary read carrying the `tc` tag pointing at
/// `primary`'s template coordinate, so it travels in that pair's position
/// group regardless of its own placement — mirroring what `fgumi zipper`
/// stamps in production.
fn build_non_primary(
    name: &str,
    umi: &str,
    primary: (&RawRecord, &RawRecord),
    supplementary: bool,
) -> RawRecord {
    let seq = b"ACGTACGT";
    let cigar_op = (u32::try_from(seq.len()).expect("fits u32")) << 4;
    let non_primary_flag = if supplementary { flags::SUPPLEMENTARY } else { flags::SECONDARY };
    let tc = primary_tc_array(primary.0, primary.1);
    let mut b = SamBuilder::new();
    b.read_name(name.as_bytes())
        .sequence(seq)
        .qualities(&[30; 8])
        .flags(flags::PAIRED | flags::FIRST_SEGMENT | non_primary_flag)
        // Placed far from the primary; grouping must key on `tc`, not this.
        .ref_id(0)
        .pos(tc[1] + 2000 - 1)
        .mapq(60)
        .cigar_ops(&[cigar_op])
        .mate_ref_id(0)
        .mate_pos(tc[4] - 1)
        .add_string_tag(SamTag::RX, umi.as_bytes())
        .add_string_tag(SamTag::MC, b"8M")
        .add_array_i32(SamTag::TC, &tc);
    b.build()
}

/// The main fixture: four independent position groups exercising, in one
/// input, everything a single strategy run needs to prove parity on:
///
/// - `(100, 300)`: a 3-read exact-UMI family (`AAAAAAAA`) plus a 2-read
///   1-mismatch family (`AAAAAAAT`) — under `--strategy identity` these are two
///   separate molecules/duplicate families; under `--strategy adjacency` they
///   cluster into one.
/// - `(500, 700)`: a plain 2-read duplicate family (`CCCCCCCC`), an ordinary
///   control group.
/// - single-end reads at position `900`: 3 unpaired reads sharing UMI
///   `GGGGGGGG` — the single-end/fragment shape.
/// - `(1100, 1300)`: a single template (`TTTTTTTT`) carrying a `tc`-keyed
///   secondary AND a `tc`-keyed supplementary read. Kept as a lone template
///   (not a duplicate family) on its own exclusive position group, because
///   `RecordPositionGrouper` coalesces a tc-keyed non-primary read into its
///   template by QNAME match against the immediately *preceding* accumulated
///   record, not by resolved-key equality alone (see
///   `RecordPositionGrouper::with_secondary_supplementary`'s doc comment) —
///   with two same-coordinate templates, which one a following non-primary
///   read lands next to depends on the sort's coordinate tie-break, not
///   insertion order, so a lone template is what makes the outcome
///   unambiguous. Covers the `secondary_reads`/`supplementary_reads`/
///   `missing_tc_tag` accounting.
fn build_family_fixture(dir: &Path) -> std::path::PathBuf {
    let input = dir.join("in.bam");
    let header = create_minimal_header("chr1", 10000);

    let mut records = Vec::new();
    records.extend(build_duplicate_family("exact", "AAAAAAAA", 3, 100, 300));
    records.extend(build_duplicate_family("near", "AAAAAAAT", 2, 100, 300));
    records.extend(build_duplicate_family("ctrl", "CCCCCCCC", 2, 500, 700));
    for i in 0..3 {
        records.push(build_single_end(&format!("frag_{i}"), "GGGGGGGG", 900, 30));
    }

    // A dedicated 2-template family at its own exclusive coordinate, one
    // template carrying a tc-keyed secondary and the other a tc-keyed
    // supplementary. `RecordPositionGrouper` coalesces a tc-keyed
    // secondary/supplementary into its template by QNAME match against the
    // immediately preceding accumulated record ("coalescing is driven by...
    // stream adjacency", not by the resolved key alone) — so each goes right
    // after its own primary pair, and this family gets its own position
    // group so no *other* template's records can land in between.
    // A single template (not a duplicate family — with two same-coordinate
    // templates, which one a same-position tc-keyed non-primary read lands
    // next to depends on the coordinate-tie-break, not insertion order) so
    // both non-primary reads unambiguously coalesce right after it.
    let (sec_r1, sec_r2) = build_pair("sec_0", "TTTTTTTT", 1100, 1300, false, 30);
    records.push(sec_r1.clone());
    records.push(sec_r2.clone());
    records.push(build_non_primary("sec_0", "TTTTTTTT", (&sec_r1, &sec_r2), false));
    records.push(build_non_primary("sec_0", "TTTTTTTT", (&sec_r1, &sec_r2), true));

    create_sorted_bam(&input, &header, &records);
    input
}

/// A duplex-shaped fixture for `--strategy paired`: two templates at the same
/// coordinates carrying the paired-UMI format (`AAAA-CCCC`) in opposite
/// physical orientations (FR vs RF) — the shape that makes the `PairedUmiAssigner`
/// assign the same molecule two distinct `/A`/`/B` suffixes, exactly like a
/// duplex-sequenced molecule read from both strands.
fn build_paired_strategy_fixture(dir: &Path) -> std::path::PathBuf {
    let input = dir.join("in.bam");
    let header = create_minimal_header("chr1", 10000);

    let mut records = Vec::new();
    // FR: R1 (unclipped 5') at 100 < R2 at 300 -> is_r1_earlier = true.
    let (r1, r2) = build_pair("duplex_fwd", "AAAA-CCCC", 100, 300, false, 30);
    records.push(r1);
    records.push(r2);
    // RF: R1 at 300 > R2 at 100 -> is_r1_earlier = false. Same physical
    // molecule read from the opposite strand.
    let (r1, r2) = build_pair("duplex_rev", "AAAA-CCCC", 300, 100, true, 30);
    records.push(r1);
    records.push(r2);

    create_sorted_bam(&input, &header, &records);
    input
}

/// A single large position group of `n` templates with pairwise-distinct UMIs,
/// large enough to cross the default `--index-threshold 100` (adjacency/edit
/// index at `--edits 1`) so both the current build and the baseline exercise
/// the indexed-clustering path rather than the brute-force one.
fn build_index_threshold_fixture(dir: &Path, n: usize) -> std::path::PathBuf {
    let input = dir.join("in.bam");
    let header = create_minimal_header("chr1", 10000);

    let bases = [b'A', b'C', b'G', b'T'];
    let mut records = Vec::with_capacity(n * 2);
    for i in 0..n {
        // Deterministic 8-mers built as a single-parity-check code: the first
        // 7 symbols are the base-4 digits of `i` (4^7 = 16384, far more than
        // `n`), and the 8th is their digit-sum mod 4. Flipping exactly one
        // base always changes that sum mod 4, so any two *distinct* codes
        // from this family differ in at least 2 positions — never edit
        // distance 1 of each other. Plain sequential base-4 digits (no parity
        // digit) do NOT have this property: consecutive indices differ only
        // in their lowest digit, chaining almost every UMI to a neighbor at
        // edit distance 1 and letting adjacency's directed capture collapse
        // the whole set into one molecule — not the "120 independent
        // molecules" this fixture needs to isolate the indexed-clustering
        // path from UMI-clustering behavior.
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

/// Runs `<bin> dedup` (no `--threads`, so the current build takes the
/// post-cutover chain path and the baseline takes its legacy engine) with
/// `RUST_LOG=info`, and returns the process output for stderr/exit assertions.
#[allow(clippy::too_many_arguments)]
fn run_dedup(
    bin: &Path,
    input: &Path,
    output: &Path,
    strategy: &str,
    metrics: Option<&Path>,
    family_size_histogram: Option<&Path>,
    index_threshold: Option<&str>,
) -> std::process::Output {
    let mut cmd = Command::new(bin);
    cmd.env("RUST_LOG", "info").args([
        "dedup",
        "-i",
        input.to_str().unwrap(),
        "-o",
        output.to_str().unwrap(),
        "--strategy",
        strategy,
        "--compression-level",
        "1",
    ]);
    if let Some(m) = metrics {
        cmd.args(["-m", m.to_str().unwrap()]);
    }
    if let Some(h) = family_size_histogram {
        cmd.args(["-H", h.to_str().unwrap()]);
    }
    if let Some(t) = index_threshold {
        cmd.args(["--index-threshold", t]);
    }
    cmd.output().unwrap_or_else(|e| panic!("failed to spawn `{}` dedup: {e}", bin.display()))
}

//////////////////////////////////////////////////////////////////////////////
// Tests
//////////////////////////////////////////////////////////////////////////////

/// A no-`--threads` run now routes through the declarative chain, which logs
/// "Pipeline `dedup` ran in" (from `pipeline::chains::finalize`) and the
/// per-step `Pipeline stats` table (from the pipeline-core runtime). The
/// retired `unified_pipeline` engine never printed either — a hand-rolled loop
/// with no step-based runtime — so this is the RED (pre-removal) -> GREEN
/// (post-removal) discriminator for the cutover.
#[test]
fn dedup_no_threads_routes_through_chain() {
    let dir = TempDir::new().expect("temp dir");
    let input = build_family_fixture(dir.path());
    let output = dir.path().join("out.bam");

    let current_bin = Path::new(env!("CARGO_BIN_EXE_fgumi"));
    let out = run_dedup(current_bin, &input, &output, "identity", None, None, None);
    let stderr = String::from_utf8_lossy(&out.stderr);
    assert!(out.status.success(), "no-threads dedup must succeed; stderr:\n{stderr}");
    assert!(
        stderr.contains("Pipeline `dedup` ran in"),
        "a no-`--threads` dedup must route through the chain (which logs the pipeline-core \
         finalize line); the unified_pipeline engine is retired. stderr:\n{stderr}"
    );
    assert!(
        stderr.contains("Starting dedup"),
        "the chain must still emit the `Starting dedup` banner; stderr:\n{stderr}"
    );
}

/// Output parity of the post-cutover chain against the pre-removal legacy-path
/// baseline binary (run with no `--threads`, i.e. its `unified_pipeline`
/// engine), across UMI strategies — plus the always-available self-consistency
/// oracle when no baseline is set.
///
/// `identity` keeps the exact/near-mismatch families in [`build_family_fixture`]
/// as two distinct molecules; `adjacency` clusters them into one. Both exercise
/// the secondary/supplementary (`tc`-keyed) and single-end position groups
/// identically, since strategy only changes UMI clustering, not template
/// shape handling.
#[rstest]
#[case::identity("identity")]
#[case::adjacency("adjacency")]
fn cutover_matches_baseline_by_strategy(#[case] strategy: &str) {
    let dir = TempDir::new().expect("temp dir");
    let input = build_family_fixture(dir.path());

    let current_bin = Path::new(env!("CARGO_BIN_EXE_fgumi"));
    let current_out = dir.path().join("current.bam");
    let current_metrics = dir.path().join("current.metrics.txt");
    let current_hist = dir.path().join("current.hist.txt");
    let current = run_dedup(
        current_bin,
        &input,
        &current_out,
        strategy,
        Some(&current_metrics),
        Some(&current_hist),
        None,
    );
    let current_stderr = String::from_utf8_lossy(&current.stderr);
    assert!(current.status.success(), "current dedup must succeed; stderr:\n{current_stderr}");

    if let Some(baseline) = baseline_bin() {
        let baseline_out = dir.path().join("baseline.bam");
        let baseline_metrics = dir.path().join("baseline.metrics.txt");
        let baseline_hist = dir.path().join("baseline.hist.txt");
        let base = run_dedup(
            &baseline,
            &input,
            &baseline_out,
            strategy,
            Some(&baseline_metrics),
            Some(&baseline_hist),
            None,
        );
        assert!(
            base.status.success(),
            "baseline dedup failed; stderr:\n{}",
            String::from_utf8_lossy(&base.stderr)
        );

        assert_eq!(
            decompressed_records_without_pg(&current_out),
            decompressed_records_without_pg(&baseline_out),
            "chain dedup output (--strategy {strategy}) diverges from the pre-removal legacy \
             baseline binary ({}) after stripping @PG — a real cutover parity bug",
            baseline.display(),
        );
        assert_eq!(
            std::fs::read_to_string(&current_metrics).expect("current metrics"),
            std::fs::read_to_string(&baseline_metrics).expect("baseline metrics"),
            "chain --metrics TSV (--strategy {strategy}) diverges from the legacy baseline"
        );
        assert_eq!(
            std::fs::read_to_string(&current_hist).expect("current histogram"),
            std::fs::read_to_string(&baseline_hist).expect("baseline histogram"),
            "chain --family-size-histogram TSV (--strategy {strategy}) diverges from the legacy \
             baseline"
        );
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
    let current = run_dedup(current_bin, &input, &current_out, "paired", None, None, None);
    let current_stderr = String::from_utf8_lossy(&current.stderr);
    assert!(
        current.status.success(),
        "current dedup (paired) must succeed; stderr:\n{current_stderr}"
    );

    if let Some(baseline) = baseline_bin() {
        let baseline_out = dir.path().join("baseline.bam");
        let base = run_dedup(&baseline, &input, &baseline_out, "paired", None, None, None);
        assert!(
            base.status.success(),
            "baseline dedup (paired) failed; stderr:\n{}",
            String::from_utf8_lossy(&base.stderr)
        );
        assert_eq!(
            decompressed_records_without_pg(&current_out),
            decompressed_records_without_pg(&baseline_out),
            "chain dedup output (--strategy paired, duplex A/B shape) diverges from the \
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
        // Two distinct molecule ids (one per duplex strand: `/A` and `/B`), not
        // one (which would mean the strand split collapsed) and not four
        // (which would mean nothing paired up at all).
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
    // `--index-threshold always` is a CLI-level assertion, not a test-side
    // one: `validate_index_threshold` rejects it up front if the resolved
    // strategy/edits combination can never index (see `--index-threshold`'s
    // help), so a successful run here already proves indexing was reachable
    // for `adjacency` at the default `--edits 1`. Neither this test nor the
    // fallback self-consistency oracle below can observe, from the output
    // alone, whether the indexed or brute-force code path actually ran for a
    // given group — both must produce byte-identical clustering by
    // construction, so there is no output-visible difference to assert on.
    let current =
        run_dedup(current_bin, &input, &current_out, "adjacency", None, None, Some("always"));
    let current_stderr = String::from_utf8_lossy(&current.stderr);
    assert!(
        current.status.success(),
        "current dedup (index-threshold always) must succeed; stderr:\n{current_stderr}"
    );

    if let Some(baseline) = baseline_bin() {
        let baseline_out = dir.path().join("baseline.bam");
        let base =
            run_dedup(&baseline, &input, &baseline_out, "adjacency", None, None, Some("always"));
        assert!(
            base.status.success(),
            "baseline dedup (index-threshold always) failed; stderr:\n{}",
            String::from_utf8_lossy(&base.stderr)
        );
        assert_eq!(
            decompressed_records_without_pg(&current_out),
            decompressed_records_without_pg(&baseline_out),
            "chain dedup output (120-UMI indexed adjacency) diverges from the pre-removal \
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
        // Every UMI is pairwise-distinct and none is within edit distance 1 of
        // another under this construction (differing base-4 digit patterns), so
        // every template must remain its own molecule — a passthrough that
        // dropped all clustering (or a bug that over-merged everything into one
        // molecule) both fail this.
        assert_eq!(
            distinct.len(),
            120,
            "120 pairwise-distinct UMIs must resolve to 120 distinct molecule ids"
        );
        // Non-vacuous duplicate-marking check: with one template per molecule,
        // no read should be flagged a duplicate.
        let dup_count = records.iter().filter(|r| r.flags().is_duplicate()).count();
        assert_eq!(dup_count, 0, "singleton families must not be marked duplicate");
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

/// A single record's `MI:Z` tag value, if present.
fn mi_of(record: &noodles::sam::alignment::RecordBuf) -> Option<String> {
    use noodles::sam::alignment::record::data::field::Tag;
    use noodles::sam::alignment::record_buf::data::field::Value;

    match record.data().get(&Tag::from(SamTag::MI)) {
        Some(Value::String(s)) => Some(s.to_string()),
        _ => None,
    }
}

/// Every record whose read name starts with `prefix`.
fn records_named<'a>(
    records: &'a [noodles::sam::alignment::RecordBuf],
    prefix: &str,
) -> Vec<&'a noodles::sam::alignment::RecordBuf> {
    records
        .iter()
        .filter(|r| {
            r.name().is_some_and(|n| String::from_utf8_lossy(n.as_ref()).starts_with(prefix))
        })
        .collect()
}

/// All records matching read name `name` exactly, asserting every one of them
/// agrees on the `0x400` duplicate flag (a template's primaries — and, for
/// `sec_0`, its tc-keyed non-primary records — must always propagate together;
/// a mismatch here is itself a bug independent of which value is expected).
/// Returns that shared flag.
fn is_duplicate_by_name(records: &[noodles::sam::alignment::RecordBuf], name: &str) -> bool {
    let matching: Vec<_> = records
        .iter()
        .filter(|r| r.name().is_some_and(|n| AsRef::<[u8]>::as_ref(n) == name.as_bytes()))
        .collect();
    assert!(!matching.is_empty(), "no output record named {name:?}");
    let flags: std::collections::HashSet<bool> =
        matching.iter().map(|r| r.flags().is_duplicate()).collect();
    assert_eq!(
        flags.len(),
        1,
        "every record named {name:?} must agree on the duplicate flag, got {matching:?}"
    );
    flags.into_iter().next().expect("checked non-empty above")
}

/// Always-available oracle used when no baseline binary is set: proves the
/// chain's dedup output on [`build_family_fixture`] is a real transform, not a
/// passthrough — every kept record carries an `MI` tag, the exact/near-UMI
/// families come out as the strategy-appropriate number of distinct molecules,
/// EVERY template's duplicate flag matches its known representative-selection
/// outcome (pinned per family, per strategy, as exact literals — not just "some
/// duplicate exists somewhere"), the tc-keyed secondary/supplementary reads
/// survive and carry their primary's molecule id, and single-end fragments
/// form their own (unpaired) molecule with the correct representative kept.
///
/// A chain that mismarks duplicates within a family (wrong representative,
/// a family that should merge but doesn't, or vice versa) fails one of the
/// per-template assertions below, not just a weakened aggregate count.
fn assert_self_consistent_family_fixture(output: &Path, strategy: &str) {
    let (_, records) = read_bam_output(output);
    // exact(3 pairs=6) + near(2 pairs=4) + ctrl(2 pairs=4) + frag(3) +
    // sec_0(1 pair=2 primaries + 1 secondary + 1 supplementary=4) = 21.
    assert_eq!(records.len(), 21, "no record may be dropped or added by dedup (mark-only mode)");
    assert_eq!(mi_strings(&records).len(), 21, "every kept record must carry an MI tag");

    // The exact(3)+near(2) position group: under identity, 2 distinct
    // molecules (AAAAAAAA and AAAAAAAT never cluster); under adjacency, they
    // cluster into 1. Distinguishing this is the whole point of running both
    // strategies over the same fixture — a strategy that silently ignored the
    // UMI would collapse to a different count than expected here.
    let mut exact_near = records_named(&records, "exact_");
    exact_near.extend(records_named(&records, "near_"));
    let exact_near_mi: std::collections::HashSet<String> =
        exact_near.iter().filter_map(|r| mi_of(r)).collect();
    let expected_molecules = if strategy == "adjacency" { 1 } else { 2 };
    assert_eq!(
        exact_near_mi.len(),
        expected_molecules,
        "--strategy {strategy}: exact(AAAAAAAA)+near(AAAAAAAT) families must resolve to \
         {expected_molecules} distinct molecule id(s), got {exact_near_mi:?}"
    );

    // Per-template duplicate flags, pinned exactly from `build_duplicate_family`'s
    // known base qualities (20/25/30 for indices 0/1/2 — see
    // `score_template`/`PICARD_MIN_BASE_QUALITY`, which every base here clears):
    // the highest-quality template in a molecule is kept, every other template
    // in that SAME molecule is marked duplicate. `exact_2` (quality 30) is the
    // overall maximum across both `exact` and `near`, so it is the kept
    // representative under EITHER strategy. `near_1` (quality 25) is only its
    // own family's representative when `near` is NOT merged into `exact`
    // (identity) — once adjacency merges the two families, `near_1` loses to
    // `exact_2` and becomes a duplicate too.
    let near_1_kept_alone = strategy != "adjacency";
    for (name, expect_duplicate) in [
        ("exact_0", true),
        ("exact_1", true),
        ("exact_2", false), // highest quality (30) in the merged or unmerged group
        ("near_0", true),
        ("near_1", !near_1_kept_alone),
    ] {
        assert_eq!(
            is_duplicate_by_name(&records, name),
            expect_duplicate,
            "--strategy {strategy}: {name} duplicate flag must be {expect_duplicate}"
        );
    }

    // The dedicated `sec_0` template (its own exclusive position group): both
    // the secondary and the supplementary read must survive, carry the same
    // MI as the primary they were tc-keyed into, and — a singleton family, so
    // never a duplicate — leave the `0x400` flag unset on every one of them
    // (checked via `is_duplicate_by_name`, which also asserts internal
    // agreement across the primary pair + both non-primary records).
    let sec_family = records_named(&records, "sec_0");
    assert_eq!(sec_family.len(), 4, "1 primary pair + 1 secondary + 1 supplementary");
    let sec_mi: std::collections::HashSet<String> =
        sec_family.iter().filter_map(|r| mi_of(r)).collect();
    assert_eq!(sec_mi.len(), 1, "the lone sec_0 template is one molecule");
    let non_primary: Vec<_> = sec_family
        .iter()
        .filter(|r| r.flags().is_secondary() || r.flags().is_supplementary())
        .collect();
    assert_eq!(non_primary.len(), 2, "both the secondary and supplementary read must survive");
    for r in &non_primary {
        assert!(
            mi_of(r).is_some_and(|mi| sec_mi.contains(&mi)),
            "secondary/supplementary reads must carry the tc-keyed template's MI tag"
        );
    }
    assert!(
        !is_duplicate_by_name(&records, "sec_0"),
        "a singleton template (including its tc-keyed secondary/supplementary) \
         must not be marked duplicate"
    );

    // Control family (2 duplicate templates at 500/700, qualities 20/25):
    // exactly 1 molecule, `ctrl_1` (the higher quality) kept, `ctrl_0` marked.
    let ctrl = records_named(&records, "ctrl_");
    let ctrl_mi: std::collections::HashSet<String> = ctrl.iter().filter_map(|r| mi_of(r)).collect();
    assert_eq!(ctrl_mi.len(), 1, "the 2-template control family must resolve to 1 molecule");
    assert!(is_duplicate_by_name(&records, "ctrl_0"), "ctrl_0 (quality 20) must be duplicate");
    assert!(
        !is_duplicate_by_name(&records, "ctrl_1"),
        "ctrl_1 (quality 25, the family's highest) must be kept"
    );

    // Single-end fragments (3 unpaired reads, same UMI, same position, all
    // quality 30 — a full 3-way tie). `mark_duplicates_in_family`'s size-3 fast
    // path resolves a tie by keeping the FIRST record and marking the rest, so
    // `frag_0` is kept and `frag_1`/`frag_2` are duplicates, deterministically
    // (not merely "2 of the 3 are marked, order unspecified").
    let frag = records_named(&records, "frag_");
    assert_eq!(frag.len(), 3, "all 3 single-end fragments must survive");
    assert!(frag.iter().all(|r| !r.flags().is_segmented()), "fragments must not carry PAIRED");
    let frag_mi: std::collections::HashSet<String> = frag.iter().filter_map(|r| mi_of(r)).collect();
    assert_eq!(frag_mi.len(), 1, "the 3 identical single-end fragments must resolve to 1 molecule");
    assert!(!is_duplicate_by_name(&records, "frag_0"), "frag_0 must be the kept tie-break winner");
    assert!(is_duplicate_by_name(&records, "frag_1"), "frag_1 must be marked duplicate");
    assert!(is_duplicate_by_name(&records, "frag_2"), "frag_2 must be marked duplicate");
}
