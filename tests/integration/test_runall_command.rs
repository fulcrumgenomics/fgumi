//! End-to-end CLI parity + error-path tests for `fgumi runall`.
//!
//! `runall` fuses a contiguous slice of the pipeline into one in-memory chain
//! and is spawned here as the real, built `fgumi` binary (never in-process),
//! so these tests exercise the actual CLI surface a user invokes: argument
//! parsing, per-stage `--<stage>::<flag>` wiring, and the fused execution
//! path, end to end.
//!
//! This is a representative subset of the parity coverage the design doc
//! calls for, not a full port of the (now-obsolete) `feat-runall` branch's
//! 4030-line `test_runall_parity.rs` — see
//! `docs/superpowers/specs/2026-09-04-runall-pr-b-command-design.md` for the
//! full contract. It covers:
//!
//! * **Class A** self-pairs (`runall --start-from X --stop-after X` vs
//!   standalone `fgumi X`): `sort`, `group`, and (consensus-gated) `simplex`.
//! * **Class B** compositions (`runall` vs a staged standalone chain piped
//!   through intermediate temp BAMs): `sort→group`, and (consensus-gated)
//!   `sort→simplex`, `group→simplex`, `group→simplex→filter`.
//! * The extract-fusion path (aligner-gated) and the now-dropped
//!   `extract→zipper` "incompatible" guard (spec §6 rule 3).
//! * Determinism, stdin-once, `--help`, and the CLI-level error-path
//!   messages the design doc pins verbatim.
//!
//! Every parity test asserts BOTH record-stream equivalence AND header
//! equivalence ignoring `@PG` (a fused runall chain writes one `@PG`; the
//! staged standalone chain writes one per command) — a header-only
//! regression (dropped `@SQ`/`@RG`, wrong `@HD` `SO`/`GO`/`SS`) would be
//! invisible to a records-only comparison.

use std::ffi::OsStr;
use std::io::Write as _;
use std::path::{Path, PathBuf};
use std::process::{Command, Output, Stdio};

use tempfile::TempDir;

use crate::helpers::bam_generator::{
    create_minimal_header, create_test_reference, create_umi_family_at_pos, write_bam,
};
use crate::helpers::cutover::decompressed_records_without_pg;
use crate::helpers::fastq::write_bgzf_fastq_with_corrupt_last_crc;
use crate::helpers::read_bam_output;
use crate::helpers::{aligner_binary, build_aligner_index, write_gzip_fastq};

// ─────────────────────────── process + assertion helpers ───────────────────────────

/// Spawns the real, built `fgumi` binary with `args` and captures its output.
///
/// Every test in this file drives `runall` (and the standalone commands it is
/// compared against) as a subprocess — never in-process — so these tests pin
/// the actual CLI surface (argument parsing included), not just the library
/// API underneath it.
fn fgumi<I, S>(args: I) -> Output
where
    I: IntoIterator<Item = S>,
    S: AsRef<OsStr>,
{
    Command::new(env!("CARGO_BIN_EXE_fgumi")).args(args).output().expect("failed to spawn fgumi")
}

/// Runs `fgumi` with `args`, panicking with `context` and the captured
/// stderr if it exits non-zero.
fn run_ok<I, S>(args: I, context: &str) -> Output
where
    I: IntoIterator<Item = S>,
    S: AsRef<OsStr>,
{
    let output = fgumi(args);
    assert!(
        output.status.success(),
        "{context} failed (status {:?}): stderr={}",
        output.status,
        String::from_utf8_lossy(&output.stderr)
    );
    output
}

/// Runs `fgumi` with `args` and asserts it exits non-zero with `needle`
/// somewhere in stderr.
fn assert_rejected_with<I, S>(args: I, needle: &str, context: &str)
where
    I: IntoIterator<Item = S>,
    S: AsRef<OsStr>,
{
    let output = fgumi(args);
    let stderr = String::from_utf8_lossy(&output.stderr);
    assert!(!output.status.success(), "{context}: expected failure but the command succeeded");
    assert!(
        stderr.contains(needle),
        "{context}: expected stderr to contain {needle:?}, got:\n{stderr}"
    );
}

fn p(path: &Path) -> &str {
    path.to_str().expect("test path is valid UTF-8")
}

/// Record-stream equivalence: both BAMs are non-empty and their decoded
/// record streams are identical.
///
/// Standalone-command `@PG` provenance is irrelevant to records, so no
/// normalization is needed here (unlike the header comparison below).
fn assert_bams_record_equivalent_nonempty(a: &Path, b: &Path) {
    let (_, records_a) = read_bam_output(a);
    let (_, records_b) = read_bam_output(b);
    assert!(!records_a.is_empty(), "{} produced no records", a.display());
    assert_eq!(
        records_a,
        records_b,
        "record streams differ between {} and {}",
        a.display(),
        b.display()
    );
}

/// Header equivalence ignoring `@PG` provenance (and `@HD` `VN` / `@CO`,
/// which this deliberately does not inspect): compares `@HD` `SO`/`GO`/`SS`,
/// the `@SQ` dictionary, and `@RG` records.
///
/// A fused `runall` chain writes a single `@PG` record while the equivalent
/// staged standalone-command chain writes one per command (with different
/// `ID`/`PN`/`CL` values), so a full-header `==` would fail even when the
/// stage output is otherwise identical — hence the dedicated comparison
/// here rather than reusing `crate::helpers::read_bam_output`'s header
/// (which only normalizes `@PG` `CL`, not the `@PG` record count/identity).
fn assert_bam_headers_equivalent_ignoring_pg(a: &Path, b: &Path) {
    use noodles::sam::header::record::value::map::header::tag::{
        GROUP_ORDER, SORT_ORDER, SUBSORT_ORDER,
    };

    fn read_header(path: &Path) -> noodles::sam::Header {
        let mut reader = noodles::bam::io::Reader::new(std::io::BufReader::new(
            std::fs::File::open(path).expect("open BAM for header comparison"),
        ));
        reader.read_header().expect("read BAM header")
    }

    type SortFields = (Option<Vec<u8>>, Option<Vec<u8>>, Option<Vec<u8>>);

    fn sort_fields(header: &noodles::sam::Header) -> SortFields {
        let Some(hd) = header.header() else { return (None, None, None) };
        let other = hd.other_fields();
        (
            other.get(&SORT_ORDER).map(|v| v.to_vec()),
            other.get(&GROUP_ORDER).map(|v| v.to_vec()),
            other.get(&SUBSORT_ORDER).map(|v| v.to_vec()),
        )
    }

    let header_a = read_header(a);
    let header_b = read_header(b);

    assert_eq!(
        sort_fields(&header_a),
        sort_fields(&header_b),
        "@HD SO/GO/SS differ between {} and {}",
        a.display(),
        b.display()
    );
    assert_eq!(
        header_a.reference_sequences(),
        header_b.reference_sequences(),
        "@SQ differs between {} and {}",
        a.display(),
        b.display()
    );
    assert_eq!(
        header_a.read_groups(),
        header_b.read_groups(),
        "@RG differs between {} and {}",
        a.display(),
        b.display()
    );
}

// ─────────────────────────────────── fixtures ───────────────────────────────────

/// A small mapped, UMI-tagged BAM with three families at distinct positions,
/// deliberately written out of template-coordinate order so `sort` has real
/// work to do.
fn unsorted_bam(dir: &Path) -> PathBuf {
    let header = create_minimal_header("chr1", 10_000);
    let mut records = Vec::new();
    records.extend(create_umi_family_at_pos("ACGT", 3, "fam_a", "ACGTACGTACGT", 30, 800));
    records.extend(create_umi_family_at_pos("TGCA", 3, "fam_b", "TGCATGCATGCA", 30, 200));
    records.extend(create_umi_family_at_pos("GGCC", 3, "fam_c", "GGCCGGCCGGCC", 30, 500));
    let path = dir.join("unsorted.bam");
    write_bam(&path, &header, &records);
    path
}

/// [`unsorted_bam`] piped through a real `fgumi sort` to seed tests that
/// need an already-sorted (template-coordinate) input.
fn sorted_bam(dir: &Path) -> PathBuf {
    let unsorted = unsorted_bam(dir);
    let sorted = dir.join("sorted.bam");
    run_ok(["sort", "-i", p(&unsorted), "-o", p(&sorted)], "fixture setup: fgumi sort");
    sorted
}

/// [`sorted_bam`] piped through a real `fgumi group` (`--edits 0`, given
/// `strategy`) to seed tests that need an already-grouped (MI-tagged) input.
fn grouped_bam(dir: &Path, strategy: &str, tag: &str) -> PathBuf {
    let sorted = sorted_bam(dir);
    let grouped = dir.join(format!("grouped_{tag}.bam"));
    run_ok(
        ["group", "-i", p(&sorted), "-o", p(&grouped), "--strategy", strategy, "--edits", "0"],
        "fixture setup: fgumi group",
    );
    grouped
}

/// Deterministic pseudo-random ACGT sequence (xorshift64), long enough that a
/// short substring is very unlikely to recur elsewhere in it.
///
/// `create_test_reference`'s `ACGTACGT`-repeat sequence is periodic (period
/// 8), so any substring drawn from it realigns ambiguously everywhere and
/// gets MAPQ 0 — which `group`'s MAPQ filter (default `--min-map-q 1`) then
/// rejects. A test that needs its aligned reads to actually survive `group`
/// (not just `align`) needs a non-repetitive reference instead.
fn pseudo_random_sequence(n: usize, seed: u64) -> String {
    let bases = [b'A', b'C', b'G', b'T'];
    let mut state = seed | 1;
    (0..n)
        .map(|_| {
            state ^= state << 13;
            state ^= state >> 7;
            state ^= state << 17;
            bases[usize::try_from(state % 4).expect("state % 4 fits usize")] as char
        })
        .collect()
}

/// Writes a FASTA + `.fai` + `.dict` for a `len`-bp single-contig (`chr1`)
/// reference built from [`pseudo_random_sequence`], returning the FASTA path
/// and the sequence itself (so callers can slice out a substring to use as a
/// read that maps uniquely, unlike a substring of `create_test_reference`'s
/// periodic sequence).
fn write_unique_reference(dir: &Path, len: usize) -> (PathBuf, String) {
    let sequence = pseudo_random_sequence(len, 0x5eed_1234_5678_9abc);
    let ref_path = dir.join("unique_ref.fa");
    let mut fasta = std::fs::File::create(&ref_path).expect("create unique reference fasta");
    writeln!(fasta, ">chr1").expect("write fasta header");
    writeln!(fasta, "{sequence}").expect("write fasta sequence");
    drop(fasta);

    std::fs::write(dir.join("unique_ref.fa.fai"), format!("chr1\t{len}\t6\t{len}\t{}\n", len + 1))
        .expect("write fai");

    let mut dict = std::fs::File::create(dir.join("unique_ref.dict")).expect("create dict");
    writeln!(dict, "@HD\tVN:1.6\tSO:unsorted").expect("write dict @HD");
    writeln!(dict, "@SQ\tSN:chr1\tLN:{len}").expect("write dict @SQ");

    (ref_path, sequence)
}

/// Writes a paired gzip FASTQ pair encoding `molecules.len()` duplex
/// molecules, each with BOTH strands present, over a `write_unique_reference`
/// `sequence`.
///
/// For molecule `(p, r1_umi, r2_umi)`: strand A is a plain FR pair — R1 is
/// the forward-strand `TEMPLATE_LEN`-bp anchor at `p` (paired UMI
/// `r1_umi-r2_umi`), R2 is the anchor `SPAN` bp downstream, written as its
/// reverse complement (so it aligns REVERSE there, standard paired-end FASTQ
/// convention: the raw FASTQ read is the reverse complement of whatever the
/// BAM's reference-forward-oriented `SEQ` field ends up holding). Strand B is
/// the SAME molecule sequenced from the other end: its R1 is the reverse
/// complement of the downstream anchor (so it aligns REVERSE there, at the
/// position strand A's R2 occupied) and its R2 is the upstream anchor as-is
/// (so it aligns FORWARD at strand A's R1 position), with the UMI order
/// swapped (`r2_umi-r1_umi`). `--group::strategy paired` canonicalizes the
/// swapped UMI order and matching swapped orientation into one molecule with
/// `/A` and `/B` reads, matching the convention
/// `bam_generator::create_duplex_grouped_family`'s `strand_reads` helper uses
/// for pre-tagged fixtures — this builds the FASTQ-level equivalent from
/// scratch so a real `extract` + aligner run derives the same shape.
fn write_duplex_umi_fastq(
    dir: &Path,
    sequence: &str,
    molecules: &[(usize, &str, &str)],
) -> (PathBuf, PathBuf) {
    const TEMPLATE_LEN: usize = 36;
    const SPAN: usize = 150;
    const DEPTH: usize = 2;

    let r1 = dir.join("duplex_r1.fq.gz");
    let r2 = dir.join("duplex_r2.fq.gz");
    let qual = "I".repeat(4 + TEMPLATE_LEN);
    let mut r1_records: Vec<(String, String)> = Vec::new();
    let mut r2_records: Vec<(String, String)> = Vec::new();

    for (mi, &(p, r1_umi, r2_umi)) in molecules.iter().enumerate() {
        let anchor_near = &sequence[p..p + TEMPLATE_LEN];
        let anchor_far = &sequence[p + SPAN..p + SPAN + TEMPLATE_LEN];
        let anchor_far_rc = fgumi_dna::reverse_complement_str(anchor_far);

        for i in 0..DEPTH {
            // Strand A: R1 forward at the near anchor, R2 at the far anchor.
            r1_records.push((format!("mol{mi}_A_{i}"), format!("{r1_umi}{anchor_near}")));
            r2_records.push((format!("mol{mi}_A_{i}"), format!("{r2_umi}{anchor_far_rc}")));

            // Strand B: the same molecule read from the other end, UMI order
            // swapped, positions/orientations swapped.
            r1_records.push((format!("mol{mi}_B_{i}"), format!("{r2_umi}{anchor_far_rc}")));
            r2_records.push((format!("mol{mi}_B_{i}"), format!("{r1_umi}{anchor_near}")));
        }
    }

    let r1_slices: Vec<(&str, &str, &str)> =
        r1_records.iter().map(|(n, s)| (n.as_str(), s.as_str(), qual.as_str())).collect();
    let r2_slices: Vec<(&str, &str, &str)> =
        r2_records.iter().map(|(n, s)| (n.as_str(), s.as_str(), qual.as_str())).collect();
    write_gzip_fastq(&r1, &r1_slices);
    write_gzip_fastq(&r2, &r2_slices);
    (r1, r2)
}

/// Converts a paired-end unmapped BAM into an interleaved FASTQ byte stream
/// (R1, R2, R1, R2, ...), suitable for `bwa mem -p <ref> -`.
///
/// This is the staged oracle's hand-rolled stand-in for what `runall`'s fused
/// `Stage::Align` does internally (stream `BamTemplateBatch` records to the
/// aligner subprocess as FASTQ) — there is no standalone `fgumi align`
/// command, so a staged comparison has to reconstruct this step itself
/// rather than delegate to one.
fn bam_to_interleaved_fastq(bam_path: &Path) -> Vec<u8> {
    let mut reader = noodles::bam::io::Reader::new(std::io::BufReader::new(
        std::fs::File::open(bam_path).expect("open BAM for FASTQ conversion"),
    ));
    let header = reader.read_header().expect("read BAM header");
    let mut out = Vec::new();
    for result in reader.record_bufs(&header) {
        let record = result.expect("read BAM record");
        let name = record.name().map(|n| n.to_vec()).unwrap_or_default();
        let seq = record.sequence().as_ref().to_vec();
        let quals: Vec<u8> = record.quality_scores().as_ref().iter().map(|q| q + 33).collect();
        out.extend_from_slice(b"@");
        out.extend_from_slice(&name);
        out.push(b'\n');
        out.extend_from_slice(&seq);
        out.extend_from_slice(b"\n+\n");
        out.extend_from_slice(&quals);
        out.push(b'\n');
    }
    out
}

// ══════════════════════════ Class A: single-stage self-pairs ══════════════════════════

#[test]
fn sort_self_pair_matches_standalone_sort() {
    let tmp = TempDir::new().unwrap();
    let fixture = unsorted_bam(tmp.path());
    let runall_out = tmp.path().join("runall.bam");
    let staged_out = tmp.path().join("staged.bam");

    run_ok(
        [
            "runall",
            "--start-from",
            "sort",
            "--stop-after",
            "sort",
            "-i",
            p(&fixture),
            "-o",
            p(&runall_out),
        ],
        "runall sort->sort",
    );
    run_ok(["sort", "-i", p(&fixture), "-o", p(&staged_out)], "standalone sort");

    assert_bams_record_equivalent_nonempty(&runall_out, &staged_out);
    assert_bam_headers_equivalent_ignoring_pg(&runall_out, &staged_out);
}

#[test]
fn group_self_pair_matches_standalone_group() {
    let tmp = TempDir::new().unwrap();
    let fixture = sorted_bam(tmp.path());
    let runall_out = tmp.path().join("runall.bam");
    let staged_out = tmp.path().join("staged.bam");

    run_ok(
        [
            "runall",
            "--start-from",
            "group",
            "--stop-after",
            "group",
            "-i",
            p(&fixture),
            "-o",
            p(&runall_out),
            "--group::strategy",
            "identity",
            "--group::edits",
            "0",
        ],
        "runall group->group",
    );
    run_ok(
        [
            "group",
            "-i",
            p(&fixture),
            "-o",
            p(&staged_out),
            "--strategy",
            "identity",
            "--edits",
            "0",
        ],
        "standalone group",
    );

    assert_bams_record_equivalent_nonempty(&runall_out, &staged_out);
    assert_bam_headers_equivalent_ignoring_pg(&runall_out, &staged_out);
}

#[cfg(feature = "consensus")]
#[test]
fn simplex_self_pair_matches_standalone_simplex() {
    let tmp = TempDir::new().unwrap();
    let fixture = grouped_bam(tmp.path(), "identity", "simplex_self");
    let runall_out = tmp.path().join("runall.bam");
    let staged_out = tmp.path().join("staged.bam");

    run_ok(
        [
            "runall",
            "--start-from",
            "consensus",
            "--stop-after",
            "consensus",
            "--consensus",
            "simplex",
            "-i",
            p(&fixture),
            "-o",
            p(&runall_out),
            "--simplex::min-reads",
            "1",
        ],
        "runall consensus(simplex)->consensus(simplex)",
    );
    run_ok(
        ["simplex", "-i", p(&fixture), "-o", p(&staged_out), "--min-reads", "1"],
        "standalone simplex",
    );

    assert_bams_record_equivalent_nonempty(&runall_out, &staged_out);
    assert_bam_headers_equivalent_ignoring_pg(&runall_out, &staged_out);
}

/// Spec §7: `--methylation-mode` (+ `--ref`) must actually reach the consensus
/// stage's `#[arg(skip)]` `methylation_mode`/`reference` slots (wired via
/// `resolve_methylation_mode`/`methylation_reference` in
/// `build_stage_options_bag`), not just be accepted and silently dropped.
/// Compares the fused `consensus(simplex)` self-pair against the standalone
/// `fgumi simplex --methylation-mode em-seq --ref ...` oracle — record parity
/// between the two proves the flags were threaded through identically,
/// alongside `simplex_self_pair_matches_standalone_simplex` above proving the
/// non-methylation self-pair.
#[cfg(feature = "consensus")]
#[test]
fn simplex_self_pair_with_methylation_mode_matches_standalone() {
    let tmp = TempDir::new().unwrap();
    let fixture = grouped_bam(tmp.path(), "identity", "simplex_methylation");
    let reference = create_test_reference(tmp.path());
    let runall_out = tmp.path().join("runall.bam");
    let staged_out = tmp.path().join("staged.bam");

    run_ok(
        [
            "runall",
            "--start-from",
            "consensus",
            "--stop-after",
            "consensus",
            "--consensus",
            "simplex",
            "-i",
            p(&fixture),
            "-o",
            p(&runall_out),
            "--simplex::min-reads",
            "1",
            "--methylation-mode",
            "em-seq",
            "--ref",
            p(&reference),
        ],
        "runall consensus(simplex)+methylation-mode",
    );
    run_ok(
        [
            "simplex",
            "-i",
            p(&fixture),
            "-o",
            p(&staged_out),
            "--min-reads",
            "1",
            "--methylation-mode",
            "em-seq",
            "--ref",
            p(&reference),
        ],
        "standalone simplex+methylation-mode",
    );

    assert_bams_record_equivalent_nonempty(&runall_out, &staged_out);
    assert_bam_headers_equivalent_ignoring_pg(&runall_out, &staged_out);
}

// ══════════════════════════ Class B: multi-stage compositions ══════════════════════════

#[test]
fn sort_to_group_matches_staged_chain() {
    let tmp = TempDir::new().unwrap();
    let fixture = unsorted_bam(tmp.path());
    let runall_out = tmp.path().join("runall.bam");
    let staged_sorted = tmp.path().join("staged_sorted.bam");
    let staged_out = tmp.path().join("staged.bam");

    run_ok(
        [
            "runall",
            "--start-from",
            "sort",
            "--stop-after",
            "group",
            "-i",
            p(&fixture),
            "-o",
            p(&runall_out),
            "--group::strategy",
            "identity",
            "--group::edits",
            "0",
        ],
        "runall sort->group",
    );

    run_ok(["sort", "-i", p(&fixture), "-o", p(&staged_sorted)], "staged sort");
    run_ok(
        [
            "group",
            "-i",
            p(&staged_sorted),
            "-o",
            p(&staged_out),
            "--strategy",
            "identity",
            "--edits",
            "0",
        ],
        "staged group",
    );

    assert_bams_record_equivalent_nonempty(&runall_out, &staged_out);
    assert_bam_headers_equivalent_ignoring_pg(&runall_out, &staged_out);
}

#[cfg(feature = "consensus")]
#[test]
fn sort_to_simplex_matches_staged_chain() {
    let tmp = TempDir::new().unwrap();
    let fixture = unsorted_bam(tmp.path());
    let runall_out = tmp.path().join("runall.bam");
    let staged_sorted = tmp.path().join("staged_sorted.bam");
    let staged_grouped = tmp.path().join("staged_grouped.bam");
    let staged_out = tmp.path().join("staged.bam");

    run_ok(
        [
            "runall",
            "--start-from",
            "sort",
            "--stop-after",
            "consensus",
            "--consensus",
            "simplex",
            "-i",
            p(&fixture),
            "-o",
            p(&runall_out),
            "--group::strategy",
            "identity",
            "--group::edits",
            "0",
            "--simplex::min-reads",
            "1",
        ],
        "runall sort->simplex",
    );

    run_ok(["sort", "-i", p(&fixture), "-o", p(&staged_sorted)], "staged sort");
    run_ok(
        [
            "group",
            "-i",
            p(&staged_sorted),
            "-o",
            p(&staged_grouped),
            "--strategy",
            "identity",
            "--edits",
            "0",
        ],
        "staged group",
    );
    run_ok(
        ["simplex", "-i", p(&staged_grouped), "-o", p(&staged_out), "--min-reads", "1"],
        "staged simplex",
    );

    assert_bams_record_equivalent_nonempty(&runall_out, &staged_out);
    assert_bam_headers_equivalent_ignoring_pg(&runall_out, &staged_out);
}

// `codec` requires FR-overlapping paired-end fragments (each molecule's R1/R2
// overlap at the same position) — the single-end UMI-family fixtures used
// elsewhere in this file don't satisfy that, so `--consensus simplex` covers
// the "group -> consensus (one mode)" and "group -> consensus -> filter (one
// mode)" compositions below via those fixtures instead.
// `group_to_codec_matches_staged_chain` below gives codec its own dedicated
// FR-overlapping fixture and record/header parity test; codec's numeric-bounds
// CLI wiring is additionally exercised by `rejects_codec_with_methylation_mode`.

/// One CODEC-shaped read pair: R1 forward, R2 reverse, fully overlapping at
/// the same position (mirrors real CODEC sequencing, where R1/R2 read
/// opposite strands of the same short fragment), sharing an `RX` UMI so the
/// `Group` stage — not a pre-set `MI` — assigns the molecule id. Mirrors
/// `test_runall_chain_transitions.rs`'s `create_codec_umi_pair` (this file's
/// own copy: the two integration-test binaries share fixtures only through
/// the `helpers` module, and this one is specific to the CLI-parity shape
/// here).
fn create_codec_umi_pair(
    name: &str,
    seq: &[u8],
    qual: &[u8],
    ref_start: i32,
    umi: &str,
) -> (fgumi_raw_bam::RawRecord, fgumi_raw_bam::RawRecord) {
    use fgumi_lib::sam::SamTag;
    use fgumi_raw_bam::{SamBuilder, flags};

    let len = seq.len();
    let cigar_op = u32::try_from(len).expect("len fits u32") << 4;
    let template_length = i32::try_from(len).expect("len fits i32");
    let mate_cigar = format!("{len}M");

    let mut b1 = SamBuilder::new();
    b1.read_name(name.as_bytes())
        .sequence(seq)
        .qualities(qual)
        .cigar_ops(&[cigar_op])
        .flags(flags::PAIRED | flags::FIRST_SEGMENT | flags::MATE_REVERSE)
        .ref_id(0)
        .pos(ref_start)
        .mapq(60)
        .mate_ref_id(0)
        .mate_pos(ref_start)
        .template_length(template_length)
        .add_string_tag(SamTag::RX, umi.as_bytes())
        .add_string_tag(SamTag::MC, mate_cigar.as_bytes());

    let mut b2 = SamBuilder::new();
    b2.read_name(name.as_bytes())
        .sequence(seq)
        .qualities(qual)
        .cigar_ops(&[cigar_op])
        .flags(flags::PAIRED | flags::LAST_SEGMENT | flags::REVERSE)
        .ref_id(0)
        .pos(ref_start)
        .mapq(60)
        .mate_ref_id(0)
        .mate_pos(ref_start)
        .template_length(-template_length)
        .add_string_tag(SamTag::RX, umi.as_bytes())
        .add_string_tag(SamTag::MC, mate_cigar.as_bytes());

    (b1.build(), b2.build())
}

/// `Group→Codec` record/header parity: two overlapping CODEC-shaped read
/// pairs sharing one `RX` UMI, grouped by identity, consensus-called by
/// `runall --consensus codec` and compared against the staged standalone
/// `fgumi group | fgumi codec` chain. codec (like simplex) requires a
/// non-`Paired` group strategy (`validate_strategy_for_mode`), hence
/// `--strategy identity` here rather than `paired`.
#[cfg(feature = "consensus")]
#[test]
fn group_to_codec_matches_staged_chain() {
    let tmp = TempDir::new().unwrap();
    let header = create_minimal_header("chr1", 10_000);
    let mut records = Vec::new();
    for i in 0..2 {
        let (r1, r2) =
            create_codec_umi_pair(&format!("pair{i}"), b"ACGTACGTAC", &[30; 10], 500, "ACGT");
        records.push(r1);
        records.push(r2);
    }
    let input = tmp.path().join("codec_input.bam");
    write_bam(&input, &header, &records);

    let runall_out = tmp.path().join("runall.bam");
    let staged_grouped = tmp.path().join("staged_grouped.bam");
    let staged_out = tmp.path().join("staged.bam");

    run_ok(
        [
            "runall",
            "--start-from",
            "group",
            "--stop-after",
            "consensus",
            "--consensus",
            "codec",
            "-i",
            p(&input),
            "-o",
            p(&runall_out),
            "--group::strategy",
            "identity",
            "--group::edits",
            "0",
        ],
        "runall group->codec",
    );

    run_ok(
        [
            "group",
            "-i",
            p(&input),
            "-o",
            p(&staged_grouped),
            "--strategy",
            "identity",
            "--edits",
            "0",
        ],
        "staged group",
    );
    run_ok(["codec", "-i", p(&staged_grouped), "-o", p(&staged_out)], "staged codec");

    assert_bams_record_equivalent_nonempty(&runall_out, &staged_out);
    assert_bam_headers_equivalent_ignoring_pg(&runall_out, &staged_out);
}

#[cfg(feature = "consensus")]
#[test]
fn group_to_simplex_matches_staged_chain() {
    let tmp = TempDir::new().unwrap();
    let fixture = sorted_bam(tmp.path());
    let runall_out = tmp.path().join("runall.bam");
    let staged_grouped = tmp.path().join("staged_grouped.bam");
    let staged_out = tmp.path().join("staged.bam");

    run_ok(
        [
            "runall",
            "--start-from",
            "group",
            "--stop-after",
            "consensus",
            "--consensus",
            "simplex",
            "-i",
            p(&fixture),
            "-o",
            p(&runall_out),
            "--group::strategy",
            "identity",
            "--group::edits",
            "0",
            "--simplex::min-reads",
            "1",
        ],
        "runall group->simplex",
    );

    run_ok(
        [
            "group",
            "-i",
            p(&fixture),
            "-o",
            p(&staged_grouped),
            "--strategy",
            "identity",
            "--edits",
            "0",
        ],
        "staged group",
    );
    run_ok(
        ["simplex", "-i", p(&staged_grouped), "-o", p(&staged_out), "--min-reads", "1"],
        "staged simplex",
    );

    assert_bams_record_equivalent_nonempty(&runall_out, &staged_out);
    assert_bam_headers_equivalent_ignoring_pg(&runall_out, &staged_out);
}

#[cfg(feature = "consensus")]
#[test]
fn group_to_simplex_to_filter_matches_staged_chain() {
    let tmp = TempDir::new().unwrap();
    let fixture = sorted_bam(tmp.path());
    let runall_out = tmp.path().join("runall.bam");
    let staged_grouped = tmp.path().join("staged_grouped.bam");
    let staged_simplex = tmp.path().join("staged_simplex.bam");
    let staged_out = tmp.path().join("staged.bam");

    run_ok(
        [
            "runall",
            "--start-from",
            "group",
            "--stop-after",
            "filter",
            "--consensus",
            "simplex",
            "-i",
            p(&fixture),
            "-o",
            p(&runall_out),
            "--group::strategy",
            "identity",
            "--group::edits",
            "0",
            "--simplex::min-reads",
            "1",
            "--filter::min-reads",
            "1",
        ],
        "runall group->simplex->filter",
    );

    run_ok(
        [
            "group",
            "-i",
            p(&fixture),
            "-o",
            p(&staged_grouped),
            "--strategy",
            "identity",
            "--edits",
            "0",
        ],
        "staged group",
    );
    run_ok(
        ["simplex", "-i", p(&staged_grouped), "-o", p(&staged_simplex), "--min-reads", "1"],
        "staged simplex",
    );
    run_ok(
        ["filter", "-i", p(&staged_simplex), "-o", p(&staged_out), "--min-reads", "1"],
        "staged filter",
    );

    assert_bams_record_equivalent_nonempty(&runall_out, &staged_out);
    assert_bam_headers_equivalent_ignoring_pg(&runall_out, &staged_out);
}

/// Audit A1: the fused filter stage must thread the top-level `--methylation-mode`
/// / `--ref` into filter's methylation-aware options — `--filter::min-conversion-fraction`
/// requires BOTH to be set. Before the wiring fix, the fused stage saw
/// `MethylationMode::Disabled` and rejected the legitimately-set flag ("requires
/// --methylation-mode to be set"), so the methylation filters were unreachable
/// through runall. Compares the fused consensus(simplex)->filter chain against the
/// staged equivalent where standalone simplex and filter each receive
/// `--methylation-mode em-seq --ref`; `--min-conversion-fraction 0.0` exercises the
/// threaded flags (validation requires them) while keeping the output non-empty.
/// The consensus is unmapped, which the methylation filters cannot evaluate, so both
/// chains skip them on every record (covered by
/// `consensus_to_filter_with_methylation_skips_unmapped`).
/// Reuses the same `grouped_bam` + `create_test_reference` fixture as
/// `simplex_self_pair_with_methylation_mode_matches_standalone`, which is known to
/// produce methylation-tagged consensus records.
#[cfg(feature = "consensus")]
#[test]
fn consensus_to_filter_with_methylation_matches_staged_chain() {
    let tmp = TempDir::new().unwrap();
    let fixture = grouped_bam(tmp.path(), "identity", "filter_methylation");
    let reference = create_test_reference(tmp.path());
    let runall_out = tmp.path().join("runall.bam");
    let staged_simplex = tmp.path().join("staged_simplex.bam");
    let staged_out = tmp.path().join("staged.bam");

    run_ok(
        [
            "runall",
            "--start-from",
            "consensus",
            "--stop-after",
            "filter",
            "--consensus",
            "simplex",
            "-i",
            p(&fixture),
            "-o",
            p(&runall_out),
            "--simplex::min-reads",
            "1",
            "--filter::min-reads",
            "1",
            "--filter::min-conversion-fraction",
            "0.0",
            "--methylation-mode",
            "em-seq",
            "--ref",
            p(&reference),
        ],
        "runall consensus(simplex)->filter+methylation",
    );

    run_ok(
        [
            "simplex",
            "-i",
            p(&fixture),
            "-o",
            p(&staged_simplex),
            "--min-reads",
            "1",
            "--methylation-mode",
            "em-seq",
            "--ref",
            p(&reference),
        ],
        "staged simplex+methylation",
    );
    run_ok(
        [
            "filter",
            "-i",
            p(&staged_simplex),
            "-o",
            p(&staged_out),
            "--min-reads",
            "1",
            "--min-conversion-fraction",
            "0.0",
            "--methylation-mode",
            "em-seq",
            "--ref",
            p(&reference),
        ],
        "staged filter+methylation",
    );

    assert_bams_record_equivalent_nonempty(&runall_out, &staged_out);
    assert_bam_headers_equivalent_ignoring_pg(&runall_out, &staged_out);
}

/// The methylation filters need aligned records to find informative positions. A fused
/// consensus->filter chain's consensus is unmapped, so every record is left unfiltered by
/// them, and the run must say so with the count at the end rather than fail or stay silent.
#[cfg(feature = "consensus")]
#[test]
fn consensus_to_filter_with_methylation_skips_unmapped() {
    let tmp = TempDir::new().unwrap();
    let fixture = grouped_bam(tmp.path(), "identity", "filter_methylation_unmapped");
    let reference = create_test_reference(tmp.path());
    let out = tmp.path().join("runall.bam");

    let output = run_ok(
        [
            "runall",
            "--start-from",
            "consensus",
            "--stop-after",
            "filter",
            "--consensus",
            "simplex",
            "-i",
            p(&fixture),
            "-o",
            p(&out),
            "--simplex::min-reads",
            "1",
            "--filter::min-reads",
            "1",
            "--filter::min-conversion-fraction",
            "0.0",
            "--methylation-mode",
            "em-seq",
            "--ref",
            p(&reference),
        ],
        "runall consensus(simplex)->filter+methylation on unmapped consensus",
    );
    let stderr = String::from_utf8_lossy(&output.stderr);
    let (_, records) = read_bam_output(&out);
    let expected =
        format!("none of the {} records were checked by the methylation filters", records.len());
    assert!(!records.is_empty(), "the unmapped consensus must be kept");
    assert!(stderr.contains(&expected), "expected {expected:?} in stderr, got:\n{stderr}");
}

// ══════════════════════════ Extract→Correct (no aligner) ══════════════════════════

/// A small paired gzip FASTQ pair (`r1.fq.gz`, `r2.fq.gz`) with a 4 bp UMI on
/// R1 only (read structures `4M+T` / `+T`) — 2 UMI families x 3 read pairs
/// each. Both R1 UMIs (`ACGT`, `TGCA`) are exact entries in the correct
/// step's own whitelist, so every extracted record is an exact match and
/// `correct` keeps it (no UMI rejects), which is what makes the plain
/// (non-`--rejects`) parity case below meaningful. Returns `(r1_path,
/// r2_path)`.
fn write_extract_correct_fastqs(dir: &Path) -> (PathBuf, PathBuf) {
    let r1 = dir.join("ec_r1.fq.gz");
    let r2 = dir.join("ec_r2.fq.gz");
    let families =
        [("ACGT", "ACGTACGTACGT", "GGTTAACCGGTT"), ("TGCA", "TGCATGCATGCA", "CCAATTGGCCAA")];
    let r1_qual = "I".repeat(4 + 12);
    let r2_qual = "I".repeat(12);
    let mut r1_records: Vec<(String, String)> = Vec::new();
    let mut r2_records: Vec<(String, String)> = Vec::new();
    for (fi, (umi, r1_tmpl, r2_tmpl)) in families.iter().enumerate() {
        for i in 0..3 {
            let name = format!("fam{fi}_{i}");
            r1_records.push((name.clone(), format!("{umi}{r1_tmpl}")));
            r2_records.push((name, (*r2_tmpl).to_string()));
        }
    }
    let r1_slices: Vec<(&str, &str, &str)> =
        r1_records.iter().map(|(n, s)| (n.as_str(), s.as_str(), r1_qual.as_str())).collect();
    let r2_slices: Vec<(&str, &str, &str)> =
        r2_records.iter().map(|(n, s)| (n.as_str(), s.as_str(), r2_qual.as_str())).collect();
    write_gzip_fastq(&r1, &r1_slices);
    write_gzip_fastq(&r2, &r2_slices);
    (r1, r2)
}

/// `runall --start-from extract --stop-after correct` vs the staged standalone
/// `fgumi extract | fgumi correct` chain — no aligner involved, so this closes
/// the "extract→correct builder change only verified by an aligner-gated
/// test" gap (the extract→correct chain-builder wiring itself is exercised
/// in-process by `test_runall_chain_transitions.rs`'s
/// `extract_to_correct_chain_builds_and_runs`; this is the CLI-level parity
/// counterpart).
#[test]
fn extract_to_correct_matches_staged_chain() {
    let tmp = TempDir::new().unwrap();
    let (r1, r2) = write_extract_correct_fastqs(tmp.path());
    let runall_out = tmp.path().join("runall.bam");
    let staged_extracted = tmp.path().join("staged_extracted.bam");
    let staged_out = tmp.path().join("staged.bam");

    run_ok(
        [
            "runall",
            "--start-from",
            "extract",
            "--stop-after",
            "correct",
            "--extract::inputs",
            p(&r1),
            p(&r2),
            "--extract::read-structures",
            "4M+T",
            "+T",
            "--extract::sample",
            "s1",
            "--extract::library",
            "lib1",
            "--correct::umis",
            "ACGT",
            "--correct::umis",
            "TGCA",
            "--correct::min-distance",
            "1",
            "-o",
            p(&runall_out),
        ],
        "runall extract->correct",
    );

    run_ok(
        [
            "extract",
            "--inputs",
            p(&r1),
            p(&r2),
            "--read-structures",
            "4M+T",
            "+T",
            "--sample",
            "s1",
            "--library",
            "lib1",
            "-o",
            p(&staged_extracted),
        ],
        "staged extract",
    );
    run_ok(
        [
            "correct",
            "-i",
            p(&staged_extracted),
            "-o",
            p(&staged_out),
            "--umis",
            "ACGT",
            "--umis",
            "TGCA",
            "--min-distance",
            "1",
        ],
        "staged correct",
    );

    assert_bams_record_equivalent_nonempty(&runall_out, &staged_out);
    assert_bam_headers_equivalent_ignoring_pg(&runall_out, &staged_out);
}

/// Same extract→correct self-pair as above, but with a top-level `--rejects`
/// — exercising the 2-output rejects branch off the `add_correct`
/// `BamTemplateBatch` tail (`build_stage_options_bag`'s self-pair rule: honor
/// top-level `--rejects`, falling back to `--correct::rejects`), previously
/// untested. Every UMI in this fixture is an exact whitelist match (see
/// [`write_extract_correct_fastqs`]), so no UMI is actually rejected — the
/// rejects BAM is expected to be header-only (still non-empty as bytes: a
/// BGZF header + EOF block) rather than record-bearing. What this test
/// actually locks down is that (a) the run succeeds with `--rejects` wired
/// through a self-pair correct stage fed by extract, (b) the rejects file is
/// created, and (c) the kept output is unaffected by rejects tracking being
/// enabled, by comparing it against the same staged oracle used above (run
/// with its own `--rejects`).
#[test]
fn extract_to_correct_with_rejects_matches_staged_chain() {
    let tmp = TempDir::new().unwrap();
    let (r1, r2) = write_extract_correct_fastqs(tmp.path());
    let runall_out = tmp.path().join("runall.bam");
    let runall_rejects = tmp.path().join("runall_rejects.bam");
    let staged_extracted = tmp.path().join("staged_extracted.bam");
    let staged_out = tmp.path().join("staged.bam");
    let staged_rejects = tmp.path().join("staged_rejects.bam");

    run_ok(
        [
            "runall",
            "--start-from",
            "extract",
            "--stop-after",
            "correct",
            "--extract::inputs",
            p(&r1),
            p(&r2),
            "--extract::read-structures",
            "4M+T",
            "+T",
            "--extract::sample",
            "s1",
            "--extract::library",
            "lib1",
            "--correct::umis",
            "ACGT",
            "--correct::umis",
            "TGCA",
            "--correct::min-distance",
            "1",
            "--rejects",
            p(&runall_rejects),
            "-o",
            p(&runall_out),
        ],
        "runall extract->correct --rejects",
    );
    assert!(
        std::fs::metadata(&runall_rejects).is_ok(),
        "runall's --rejects file was not created: {}",
        runall_rejects.display()
    );
    let (_, runall_rejects_records) = read_bam_output(&runall_rejects);
    assert!(
        runall_rejects_records.is_empty(),
        "every UMI in this fixture is an exact whitelist match, so runall's rejects BAM must be \
         header-only, got {} record(s)",
        runall_rejects_records.len()
    );

    run_ok(
        [
            "extract",
            "--inputs",
            p(&r1),
            p(&r2),
            "--read-structures",
            "4M+T",
            "+T",
            "--sample",
            "s1",
            "--library",
            "lib1",
            "-o",
            p(&staged_extracted),
        ],
        "staged extract",
    );
    run_ok(
        [
            "correct",
            "-i",
            p(&staged_extracted),
            "-o",
            p(&staged_out),
            "--umis",
            "ACGT",
            "--umis",
            "TGCA",
            "--min-distance",
            "1",
            "--rejects",
            p(&staged_rejects),
        ],
        "staged correct --rejects",
    );
    assert!(
        std::fs::metadata(&staged_rejects).is_ok(),
        "staged --rejects file was not created: {}",
        staged_rejects.display()
    );
    let (_, staged_rejects_records) = read_bam_output(&staged_rejects);
    assert!(
        staged_rejects_records.is_empty(),
        "every UMI in this fixture is an exact whitelist match, so the staged rejects BAM must be \
         header-only, got {} record(s)",
        staged_rejects_records.len()
    );

    assert_bams_record_equivalent_nonempty(&runall_out, &staged_out);
    assert_bam_headers_equivalent_ignoring_pg(&runall_out, &staged_out);
}

/// Writes a `fgumi simulate aligner` replay BAM standing in for the aligner on
/// an `extract→correct→align` chain over `(r1, r2)`: the standalone
/// `extract | correct` kept output, which is exactly the template stream the
/// fused chain feeds its aligner, in the same order. Its records stay unmapped
/// (valid aligner output), and the header gains the reference's `@SQ` so
/// `AlignAndMerge`'s dict check passes. This keeps the fused-chain tests
/// hermetic: no real aligner, so they run in CI. Returns `(kept, replay)`,
/// `kept` being that `extract | correct` output.
#[cfg(feature = "simulate")]
fn write_correct_replay_bam(
    dir: &Path,
    r1: &Path,
    r2: &Path,
    reference_len: usize,
) -> (PathBuf, PathBuf) {
    use noodles::sam::alignment::io::Write as _;
    use noodles::sam::header::record::value::{Map, map::ReferenceSequence};

    let extracted = dir.join("replay_extracted.bam");
    let kept = dir.join("replay_kept.bam");
    run_ok(
        [
            "extract",
            "--inputs",
            p(r1),
            p(r2),
            "--read-structures",
            "4M+T",
            "4M+T",
            "--sample",
            "s1",
            "--library",
            "lib1",
            "-o",
            p(&extracted),
        ],
        "replay fixture: extract",
    );
    run_ok(
        [
            "correct",
            "-i",
            p(&extracted),
            "-o",
            p(&kept),
            "--umis",
            "AAAA",
            "--umis",
            "CCCC",
            "--min-distance",
            "1",
        ],
        "replay fixture: correct",
    );

    let (mut header, records) = read_bam_output(&kept);
    header.reference_sequences_mut().insert(
        bstr::BString::from("chr1"),
        Map::<ReferenceSequence>::new(std::num::NonZeroUsize::new(reference_len).unwrap()),
    );
    let replay = dir.join("replay.bam");
    let mut writer = noodles::bam::io::Writer::new(std::fs::File::create(&replay).unwrap());
    writer.write_header(&header).unwrap();
    for record in &records {
        writer.write_alignment_record(&header, record).unwrap();
    }
    writer.try_finish().unwrap();
    (kept, replay)
}

/// The shared fixture for the fused `correct` rejects tests: a 4000 bp
/// reference plus a duplex FASTQ pair whose second molecule carries UMIs
/// (`GGGG`/`TTTT`) four mismatches from every whitelist entry (`AAAA`, `CCCC`),
/// so its reads are genuinely rejected. Returns `(reference, r1, r2, aligner
/// command)`, the last replaying [`write_correct_replay_bam`].
#[cfg(feature = "simulate")]
fn correct_rejects_fixture(dir: &Path) -> (PathBuf, PathBuf, PathBuf, String) {
    const REFERENCE_LEN: usize = 4000;
    let (reference, sequence) = write_unique_reference(dir, REFERENCE_LEN);
    let molecules = [(500usize, "AAAA", "CCCC"), (2000usize, "GGGG", "TTTT")];
    let (r1, r2) = write_duplex_umi_fastq(dir, &sequence, &molecules);
    let (_, replay) = write_correct_replay_bam(dir, &r1, &r2, REFERENCE_LEN);
    let aligner_cmd = format!(
        "{} simulate aligner --replay-bam {} {{ref}}",
        env!("CARGO_BIN_EXE_fgumi"),
        replay.display()
    );
    (reference, r1, r2, aligner_cmd)
}

/// `runall --start-from extract` args through `--stop-after <stop_after>` over
/// the [`correct_rejects_fixture`] inputs, with the whitelist wired so the
/// chain includes `correct`. Callers append rejects flags and `-o`.
#[cfg(feature = "simulate")]
fn correct_chain_args<'a>(
    stop_after: &'a str,
    fixture: &'a (PathBuf, PathBuf, PathBuf, String),
) -> Vec<&'a str> {
    let (reference, r1, r2, aligner_cmd) = fixture;
    let mut args = vec![
        "runall",
        "--start-from",
        "extract",
        "--stop-after",
        stop_after,
        "--extract::inputs",
        p(r1),
        p(r2),
        "--extract::read-structures",
        "4M+T",
        "4M+T",
        "--extract::sample",
        "s1",
        "--extract::library",
        "lib1",
        "--correct::umis",
        "AAAA",
        "--correct::umis",
        "CCCC",
        "--correct::min-distance",
        "1",
    ];
    if stop_after != "correct" {
        args.extend(["--ref", p(reference), "--aligner::command", aligner_cmd]);
    }
    args
}

/// Names of every record in `path`, in file order.
#[cfg(feature = "simulate")]
fn record_names(path: &Path) -> Vec<String> {
    let (_, records) = read_bam_output(path);
    records.iter().map(|r| r.name().map(ToString::to_string).unwrap_or_default()).collect()
}

/// `--correct::rejects` on a chain that runs past `correct` (here
/// `extract→zipper`, so correct feeds the fused Align stage) must write the
/// same UMI rejects as the `extract→correct` self-pair. It used to be dropped
/// silently: the run logged the rejected count, exited 0, and wrote no file.
/// Hermetic via the [`write_correct_replay_bam`] replay aligner.
#[cfg(feature = "simulate")]
#[test]
fn extract_to_zipper_writes_correct_rejects_like_self_pair() {
    let tmp = TempDir::new().unwrap();
    let fixture = correct_rejects_fixture(tmp.path());

    let run = |stop_after: &str, out: &Path, rejects: &Path| {
        let mut args = correct_chain_args(stop_after, &fixture);
        args.extend(["--correct::rejects", p(rejects), "-o", p(out)]);
        run_ok(args, &format!("runall extract->{stop_after} --correct::rejects"));
    };

    let self_pair_rejects = tmp.path().join("self_pair_rejects.bam");
    run("correct", &tmp.path().join("self_pair.bam"), &self_pair_rejects);
    let chained_rejects = tmp.path().join("chained_rejects.bam");
    let chained_out = tmp.path().join("chained.bam");
    run("zipper", &chained_out, &chained_rejects);

    assert!(
        chained_rejects.exists(),
        "--correct::rejects was dropped on the extract->zipper chain: {}",
        chained_rejects.display()
    );
    // Molecule 1 (2 strands x 2 read-pairs x 2 records) is exactly what is rejected.
    let rejected = record_names(&self_pair_rejects);
    assert_eq!(rejected.len(), 8, "expected molecule 1's 8 records to be rejected: {rejected:?}");
    assert!(
        rejected.iter().all(|n| n.starts_with("mol1_")),
        "only molecule 1 may be rejected: {rejected:?}"
    );
    assert_eq!(
        decompressed_records_without_pg(&chained_rejects),
        decompressed_records_without_pg(&self_pair_rejects),
        "chained correct's rejects must match the self-pair's byte for byte (@PG aside)"
    );
    // Kept and rejects must partition the input: molecule 0's templates all
    // reach the chained output and none of molecule 1's leak into it.
    let kept = record_names(&chained_out);
    assert!(
        kept.iter().all(|n| n.starts_with("mol0_")),
        "rejected molecule 1 leaked into the chained output: {kept:?}"
    );
    assert_eq!(kept.len(), 8, "every molecule-0 record must reach the chained output: {kept:?}");
}

/// Writes a *mapped* replay BAM for an `extract→correct→align` chain over
/// `(r1, r2)`, plus the standalone `extract | correct` kept BAM it pairs with.
///
/// Every template becomes R1 forward at `100 + 40·t`, R2 reverse at `2000 +
/// 40·t`, and a supplementary copy of R1 at `3000`, each with SEQ copied from
/// `sequence` so the replay is valid aligner output. This gives every
/// `--zipper::*` merge rule something to change: tags on a reverse-strand read
/// for reverse/revcomp and a supplementary read for `tc`. Returns
/// `(kept, replay)`.
#[cfg(feature = "simulate")]
fn write_mapped_replay_bam(dir: &Path, r1: &Path, r2: &Path, sequence: &str) -> (PathBuf, PathBuf) {
    use noodles::core::Position;
    use noodles::sam::alignment::io::Write as _;
    use noodles::sam::alignment::record::Flags;
    use noodles::sam::alignment::record::cigar::op::{Kind, Op};
    use noodles::sam::alignment::record_buf::{Cigar, Sequence};
    use noodles::sam::header::record::value::{Map, map::ReferenceSequence};

    // Reuse the unmapped-replay fixture's kept BAM: it is exactly the stream
    // the fused chain hands its aligner.
    let (kept, _) = write_correct_replay_bam(dir, r1, r2, sequence.len());

    let (mut header, records) = read_bam_output(&kept);
    header.reference_sequences_mut().insert(
        bstr::BString::from("chr1"),
        Map::<ReferenceSequence>::new(std::num::NonZeroUsize::new(sequence.len()).unwrap()),
    );

    let place =
        |record: &noodles::sam::alignment::RecordBuf, pos: usize, mate_pos: usize, flags: Flags| {
            let mut out = record.clone();
            let len = record.sequence().len();
            let bases = sequence.as_bytes()[pos - 1..pos - 1 + len].to_vec();
            *out.flags_mut() = flags;
            *out.reference_sequence_id_mut() = Some(0);
            *out.alignment_start_mut() = Position::new(pos);
            *out.cigar_mut() = Cigar::from(vec![Op::new(Kind::Match, len)]);
            *out.mapping_quality_mut() = noodles::sam::alignment::record::MappingQuality::new(60);
            *out.sequence_mut() = Sequence::from(bases);
            *out.mate_reference_sequence_id_mut() = Some(0);
            *out.mate_alignment_start_mut() = Position::new(mate_pos);
            out
        };

    let replay = dir.join("mapped_replay.bam");
    let mut writer = noodles::bam::io::Writer::new(std::fs::File::create(&replay).unwrap());
    writer.write_header(&header).unwrap();
    for (t, pair) in records.chunks(2).enumerate() {
        let [first, second] = pair else { panic!("kept BAM must hold R1/R2 pairs") };
        let r1_pos = 100 + 40 * t;
        let r2_pos = 2000 + 40 * t;
        let paired = Flags::SEGMENTED | Flags::PROPERLY_SEGMENTED;
        let r1_flags = paired | Flags::FIRST_SEGMENT | Flags::MATE_REVERSE_COMPLEMENTED;
        let r2_flags = paired | Flags::LAST_SEGMENT | Flags::REVERSE_COMPLEMENTED;
        let out_r1 = place(first, r1_pos, r2_pos, r1_flags);
        let out_r2 = place(second, r2_pos, r1_pos, r2_flags);
        let out_supp = place(first, 3000, r2_pos, r1_flags | Flags::SUPPLEMENTARY);
        for record in [&out_r1, &out_r2, &out_supp] {
            writer.write_alignment_record(&header, record).unwrap();
        }
    }
    writer.try_finish().unwrap();
    (kept, replay)
}

/// The fixture for the fused-merge `--zipper::*` tests: a 4000 bp reference, a
/// two-molecule duplex FASTQ pair, and the [`write_mapped_replay_bam`] kept and
/// replay BAMs over it.
#[cfg(feature = "simulate")]
struct MappedReplay {
    reference: PathBuf,
    kept: PathBuf,
    replay: PathBuf,
    /// `(reference, r1, r2, aligner command)` for [`correct_chain_args`].
    chain: (PathBuf, PathBuf, PathBuf, String),
}

#[cfg(feature = "simulate")]
impl MappedReplay {
    fn new(dir: &Path) -> Self {
        const REFERENCE_LEN: usize = 4000;
        let (reference, sequence) = write_unique_reference(dir, REFERENCE_LEN);
        let (r1, r2) = write_duplex_umi_fastq(
            dir,
            &sequence,
            &[(500, "AAAA", "CCCC"), (2000, "GGGG", "TTTT")],
        );
        let (kept, replay) = write_mapped_replay_bam(dir, &r1, &r2, &sequence);
        let aligner_cmd = format!(
            "{} simulate aligner --replay-bam {} {{ref}}",
            env!("CARGO_BIN_EXE_fgumi"),
            replay.display()
        );
        Self { chain: (reference.clone(), r1, r2, aligner_cmd), reference, kept, replay }
    }

    /// Runs the fused chain ending at zipper, starting from `extract` (through
    /// correct) or, when `start_from_align`, from the kept unmapped BAM, with
    /// `flags` given their `--zipper::` prefix. Returns stderr.
    fn run_chain(&self, start_from_align: bool, out: &Path, flags: &[&str]) -> String {
        let output = run_ok(
            self.chain_args(start_from_align, out, flags),
            &format!("runall ->zipper (align start: {start_from_align}) {flags:?}"),
        );
        String::from_utf8_lossy(&output.stderr).into_owned()
    }

    /// The `runall` arguments for [`Self::run_chain`], with `flags` given in their
    /// standalone-zipper spelling and prefixed with `--zipper::` here.
    fn chain_args(&self, start_from_align: bool, out: &Path, flags: &[&str]) -> Vec<String> {
        let prefixed: Vec<String> = flags
            .iter()
            .map(|a| {
                a.strip_prefix("--").map_or_else(|| (*a).to_string(), |f| format!("--zipper::{f}"))
            })
            .collect();
        let (reference, _, _, aligner_cmd) = &self.chain;
        let mut args = if start_from_align {
            vec![
                "runall",
                "--start-from",
                "align",
                "--stop-after",
                "zipper",
                "-i",
                p(&self.kept),
                "--ref",
                p(reference),
                "--aligner::command",
                aligner_cmd,
            ]
        } else {
            correct_chain_args("zipper", &self.chain)
        };
        args.extend(prefixed.iter().map(String::as_str));
        args.extend(["-o", p(out)]);
        args.into_iter().map(str::to_string).collect()
    }

    /// Standalone `fgumi zipper` over the same replay (aligner output) and kept
    /// (unmapped) BAMs — the oracle for the fused merge.
    fn run_standalone(&self, out: &Path, flags: &[&str]) {
        let mut args = vec![
            "zipper",
            "-i",
            p(&self.replay),
            "-u",
            p(&self.kept),
            "-r",
            p(&self.reference),
            "-o",
            p(out),
        ];
        args.extend_from_slice(flags);
        run_ok(args, &format!("standalone zipper {flags:?}"));
    }
}

/// Every `--zipper::*` merge rule must apply on a chain that aligns — from
/// `extract` (correct feeding the fused Align stage) and from `align` (an
/// unmapped BAM) — exactly as standalone `fgumi zipper` applies it to the same
/// aligner output. The fused merge used to be built with empty tag rules,
/// `tc` tags always on, so each flag was
/// parsed and dropped. A no-flag baseline must already match standalone
/// zipper, so a mismatch is the flag, not the fixture; and each flag must
/// change the output, so no case can pass vacuously.
#[cfg(feature = "simulate")]
#[rstest::rstest]
#[case::tags_to_remove(&["--tags-to-remove", "RX"])]
#[case::tags_to_reverse(&["--tags-to-reverse", "RX"])]
#[case::tags_to_revcomp(&["--tags-to-revcomp", "RX"])]
#[case::skip_tc_tags(&["--skip-tc-tags"])]
fn fused_chain_honors_zipper_merge_flags(
    #[case] flags: &[&str],
    #[values(false, true)] start_from_align: bool,
) {
    let tmp = TempDir::new().unwrap();
    let fixture = MappedReplay::new(tmp.path());

    let chain_default = tmp.path().join("chain_default.bam");
    fixture.run_chain(start_from_align, &chain_default, &[]);
    let chain_flagged = tmp.path().join("chain_flagged.bam");
    fixture.run_chain(start_from_align, &chain_flagged, flags);
    let oracle_default = tmp.path().join("oracle_default.bam");
    fixture.run_standalone(&oracle_default, &[]);
    let oracle = tmp.path().join("oracle.bam");
    fixture.run_standalone(&oracle, flags);

    assert_bams_record_equivalent_nonempty(&chain_default, &oracle_default);
    assert_bams_record_equivalent_nonempty(&chain_flagged, &oracle);
    let (_, default_records) = read_bam_output(&chain_default);
    let (_, flagged_records) = read_bam_output(&chain_flagged);
    assert_ne!(
        default_records, flagged_records,
        "{flags:?} must change the fused chain's output on this fixture"
    );
}

/// `--zipper::exclude-missing-reads` has nothing to act on in a fused align
/// stage (a template the aligner drops fails the run). The run must say so —
/// and only when the flag is set — and its output must be record-identical to
/// the same chain without it.
#[cfg(feature = "simulate")]
#[test]
fn fused_chain_warns_exclude_missing_reads_is_inert() {
    const WARNING: &str = "--zipper::exclude-missing-reads has no effect on a fused align stage";
    let tmp = TempDir::new().unwrap();
    let fixture = MappedReplay::new(tmp.path());

    let without = tmp.path().join("without.bam");
    let stderr = fixture.run_chain(false, &without, &[]);
    assert!(!stderr.contains(WARNING), "the warning must not fire without the flag:\n{stderr}");

    let with = tmp.path().join("with.bam");
    let stderr = fixture.run_chain(false, &with, &["--exclude-missing-reads"]);
    assert!(stderr.contains(WARNING), "expected the inert-flag warning, got:\n{stderr}");
    assert_bams_record_equivalent_nonempty(&with, &without);
}

/// `--zipper::restore-unconverted-bases` (v0.7.0 and earlier) was removed: a chain that
/// aligns, from `extract` or from `align`, must reject it with a pointer to the methylation
/// guide rather than run.
#[cfg(feature = "simulate")]
#[rstest::rstest]
fn restore_unconverted_bases_is_rejected(#[values(false, true)] start_from_align: bool) {
    let tmp = TempDir::new().unwrap();
    let fixture = MappedReplay::new(tmp.path());
    let out = tmp.path().join("out.bam");
    assert_rejected_with(
        fixture.chain_args(start_from_align, &out, &["--restore-unconverted-bases"]),
        "--restore-unconverted-bases was removed",
        &format!("runall ->zipper (align start: {start_from_align})"),
    );
}

/// A chained correct that reaches no consensus stage leaves the top-level
/// `--rejects` unconsumed; the run must say so and point at the per-stage flag
/// that captures UMI rejects.
#[cfg(feature = "simulate")]
#[test]
fn chained_correct_top_level_rejects_warns_wired_nowhere() {
    let tmp = TempDir::new().unwrap();
    let fixture = correct_rejects_fixture(tmp.path());
    let rejects = tmp.path().join("rejects.bam");
    let out = tmp.path().join("out.bam");
    let mut args = correct_chain_args("zipper", &fixture);
    args.extend(["--rejects", p(&rejects), "-o", p(&out)]);
    let output = run_ok(args, "runall extract->zipper --rejects");
    let stderr = String::from_utf8_lossy(&output.stderr);
    assert!(
        stderr.contains(
            "--rejects is wired nowhere on this runall chain (it is consumed only by a correct \
             self-pair or a consensus stage). Use --correct::rejects to capture UMI rejects."
        ),
        "expected the wired-nowhere warning with the --correct::rejects hint, got:\n{stderr}"
    );
    assert!(
        !stderr.contains("collects only the consensus stage's rejects"),
        "no consensus stage runs, so the consensus-rejects warning must not fire:\n{stderr}"
    );
    assert!(!rejects.exists(), "nothing consumes --rejects here, so no file may be written");
}

/// On a correct self-pair given both flags, the top-level `--rejects` wins; a
/// different `--correct::rejects` is dead and must be reported, not dropped.
#[cfg(feature = "simulate")]
#[test]
fn self_pair_both_rejects_flags_warns_and_writes_top_level_only() {
    let tmp = TempDir::new().unwrap();
    let fixture = correct_rejects_fixture(tmp.path());
    let top = tmp.path().join("top_rejects.bam");
    let per_stage = tmp.path().join("per_stage_rejects.bam");
    let out = tmp.path().join("out.bam");
    let mut args = correct_chain_args("correct", &fixture);
    args.extend(["--rejects", p(&top), "--correct::rejects", p(&per_stage), "-o", p(&out)]);
    let output = run_ok(args, "runall extract->correct with both rejects flags");
    let stderr = String::from_utf8_lossy(&output.stderr);
    assert!(
        stderr.contains("is ignored: on a correct self-pair the top-level --rejects"),
        "expected the ignored --correct::rejects warning, got:\n{stderr}"
    );
    assert!(!per_stage.exists(), "the ignored --correct::rejects must not be written");
    assert!(
        !stderr.contains("--rejects is wired nowhere"),
        "a self-pair consumes --rejects, so it is not dead:\n{stderr}"
    );
    let rejected = record_names(&top);
    assert_eq!(rejected.len(), 8, "the top-level --rejects must collect molecule 1: {rejected:?}");
    assert!(
        rejected.iter().all(|n| n.starts_with("mol1_")),
        "only molecule 1 may be rejected: {rejected:?}"
    );
    let kept = record_names(&out);
    assert!(
        kept.len() == 8 && kept.iter().all(|n| n.starts_with("mol0_")),
        "the kept output must be exactly molecule 0: {kept:?}"
    );
}

/// Both rejects flags naming one file — even spelled differently, before the
/// file exists — are not a conflict on a self-pair, so there is nothing to
/// warn about.
#[cfg(feature = "simulate")]
#[rstest::rstest]
#[case::same_spelling(false)]
#[case::dot_slash_spelling(true)]
fn self_pair_same_rejects_file_on_both_flags_does_not_warn(#[case] dot_slash: bool) {
    let tmp = TempDir::new().unwrap();
    let fixture = correct_rejects_fixture(tmp.path());
    let rejects = tmp.path().join("rejects.bam");
    let spelled =
        if dot_slash { tmp.path().join(".").join("rejects.bam") } else { rejects.clone() };
    let out = tmp.path().join("out.bam");
    let mut args = correct_chain_args("correct", &fixture);
    args.extend(["--rejects", p(&rejects), "--correct::rejects", p(&spelled), "-o", p(&out)]);
    let output = run_ok(args, "runall extract->correct with one rejects file on both flags");
    let stderr = String::from_utf8_lossy(&output.stderr);
    assert!(!stderr.contains("is ignored"), "one file on both flags is not ignored:\n{stderr}");
    assert_eq!(record_names(&rejects).len(), 8, "molecule 1 must be rejected");
}

/// `--start-from correct` on a chain past correct (from an extracted BAM)
/// writes the same UMI rejects as the correct self-pair on that BAM.
#[cfg(feature = "simulate")]
#[test]
fn correct_start_to_zipper_writes_correct_rejects_like_self_pair() {
    let tmp = TempDir::new().unwrap();
    let fixture = correct_rejects_fixture(tmp.path());
    // `write_correct_replay_bam` leaves the standalone extract output here.
    let extracted = tmp.path().join("replay_extracted.bam");
    let (reference, _, _, aligner_cmd) = &fixture;
    let run = |stop_after: &str, name: &str| -> PathBuf {
        let rejects = tmp.path().join(format!("{name}_rejects.bam"));
        let out = tmp.path().join(format!("{name}.bam"));
        let mut args = vec![
            "runall",
            "--start-from",
            "correct",
            "--stop-after",
            stop_after,
            "-i",
            p(&extracted),
            "--correct::umis",
            "AAAA",
            "--correct::umis",
            "CCCC",
            "--correct::min-distance",
            "1",
            "--correct::rejects",
            p(&rejects),
            "-o",
            p(&out),
        ];
        if stop_after != "correct" {
            args.extend(["--ref", p(reference), "--aligner::command", aligner_cmd]);
        }
        run_ok(args, &format!("runall correct->{stop_after}"));
        rejects
    };
    let self_pair = run("correct", "self_pair");
    let chained = run("zipper", "chained");
    let rejected = record_names(&self_pair);
    assert!(
        rejected.len() == 8 && rejected.iter().all(|n| n.starts_with("mol1_")),
        "the self-pair must reject exactly molecule 1: {rejected:?}"
    );
    assert_eq!(
        decompressed_records_without_pg(&chained),
        decompressed_records_without_pg(&self_pair),
        "chained correct's rejects must match the self-pair's"
    );
}

/// A consensus-reaching chained correct with only the top-level `--rejects`
/// must warn that UMI rejects go uncaptured; adding `--correct::rejects`
/// captures them, so the warning must go away. The warning is emitted before
/// the chain runs, so a bogus aligner command is enough.
#[cfg(feature = "consensus")]
#[test]
fn chained_consensus_rejects_warning_is_gated_on_correct_rejects() {
    let tmp = TempDir::new().unwrap();
    let input = unsorted_bam(tmp.path());
    let (reference, _) = write_unique_reference(tmp.path(), 2000);
    let consensus_rejects = tmp.path().join("consensus_rejects.bam");
    let umi_rejects = tmp.path().join("umi_rejects.bam");
    let out = tmp.path().join("out.bam");
    let base = [
        "runall",
        "--start-from",
        "correct",
        "--stop-after",
        "consensus",
        "--consensus",
        "simplex",
        "-i",
        p(&input),
        "-o",
        p(&out),
        "--correct::umis",
        "ACGT",
        "--correct::min-distance",
        "1",
        "--rejects",
        p(&consensus_rejects),
        "--ref",
        p(&reference),
        "--aligner::command",
        "not-a-real-aligner mem {ref} /dev/stdin",
        "--group::strategy",
        "adjacency",
        "--simplex::min-reads",
        "1",
    ];
    let warning = "collects only the consensus stage's rejects, not correct's UMI rejects";

    let without = fgumi(base);
    let stderr = String::from_utf8_lossy(&without.stderr);
    assert!(stderr.contains(warning), "expected the uncaptured-UMI-rejects warning:\n{stderr}");

    let with = fgumi(base.iter().copied().chain(["--correct::rejects", p(&umi_rejects)]));
    let stderr = String::from_utf8_lossy(&with.stderr);
    assert!(
        stderr.contains("not-a-real-aligner"),
        "the run must reach the aligner, past every warning, for the absence to mean anything:\n{stderr}"
    );
    assert!(!stderr.contains(warning), "--correct::rejects is set, so no warning:\n{stderr}");
}

/// On a chained run that reaches consensus, `--correct::rejects` and the
/// top-level `--rejects` are two different writers (correct's UMI rejects vs
/// the consensus stage's rejects), so one path given to both must be refused,
/// naming both flags. Needs no real aligner: the collision guard runs before
/// the chain is built, and command mode never resolves the aligner binary.
#[cfg(feature = "consensus")]
#[test]
fn chained_correct_rejects_collides_with_consensus_rejects() {
    let tmp = TempDir::new().unwrap();
    let input = unsorted_bam(tmp.path());
    let (reference, _) = write_unique_reference(tmp.path(), 2000);
    let rejects = tmp.path().join("rejects.bam");
    assert_rejected_with(
        [
            "runall",
            "--start-from",
            "correct",
            "--stop-after",
            "consensus",
            "--consensus",
            "simplex",
            "-i",
            p(&input),
            "-o",
            p(&tmp.path().join("out.bam")),
            "--correct::umis",
            "ACGT",
            "--correct::min-distance",
            "1",
            "--correct::rejects",
            p(&rejects),
            "--rejects",
            p(&rejects),
            "--ref",
            p(&reference),
            "--aligner::command",
            "not-a-real-aligner mem {ref} /dev/stdin",
            "--group::strategy",
            "adjacency",
            "--simplex::min-reads",
            "1",
        ],
        "--correct::rejects and --rejects both write to",
        "chained --correct::rejects colliding with consensus --rejects",
    );
}

/// A rejects path that aliases `--input` would truncate the BAM being read;
/// runall must refuse it up front and leave the input intact.
#[test]
fn correct_rejects_aliasing_input_is_refused() {
    let tmp = TempDir::new().unwrap();
    let input = unsorted_bam(tmp.path());
    let before = std::fs::read(&input).unwrap();
    assert_rejected_with(
        [
            "runall",
            "--start-from",
            "correct",
            "--stop-after",
            "correct",
            "-i",
            p(&input),
            "-o",
            p(&tmp.path().join("out.bam")),
            "--correct::umis",
            "ACGT",
            "--correct::min-distance",
            "1",
            "--correct::rejects",
            p(&input),
        ],
        &format!("--correct::rejects '{}' is the same file as --input", input.display()),
        "--correct::rejects aliasing --input",
    );
    assert_eq!(std::fs::read(&input).unwrap(), before, "the input BAM must be left untouched");
}

/// Runs `args`, expecting a refusal whose stderr names `needle`, and checks
/// that `victim` — the file the refused write would have truncated — is intact.
fn assert_alias_refused(args: &[&str], needle: &str, victim: &Path) {
    let before = std::fs::read(victim).unwrap();
    assert_rejected_with(args, needle, needle);
    assert_eq!(
        std::fs::read(victim).unwrap(),
        before,
        "{} must be left untouched",
        victim.display()
    );
}

/// Every file runall reads is guarded against a write target that aliases it,
/// not only `--input`: an extract FASTQ, the zipper `--unmapped` BAM, and a
/// `-o` spelled as the input.
#[test]
fn write_targets_aliasing_other_inputs_are_refused() {
    let tmp = TempDir::new().unwrap();
    let input = unsorted_bam(tmp.path());
    assert_alias_refused(
        &[
            "runall",
            "--start-from",
            "sort",
            "--stop-after",
            "sort",
            "-i",
            p(&input),
            "-o",
            p(&input),
        ],
        &format!("--output '{}' is the same file as --input", input.display()),
        &input,
    );

    let (r1, _) = write_extract_correct_fastqs(tmp.path());
    assert_alias_refused(
        &[
            "runall",
            "--start-from",
            "extract",
            "--stop-after",
            "extract",
            "--extract::inputs",
            p(&r1),
            "--extract::read-structures",
            "4M+T",
            "--extract::sample",
            "s1",
            "--extract::library",
            "lib1",
            "-o",
            p(&r1),
        ],
        &format!("--output '{}' is the same file as --extract::inputs", r1.display()),
        &r1,
    );

    let (reference, _) = write_unique_reference(tmp.path(), 2000);
    let unmapped = tmp.path().join("unmapped.bam");
    std::fs::copy(&input, &unmapped).unwrap();
    assert_alias_refused(
        &[
            "runall",
            "--start-from",
            "zipper",
            "--stop-after",
            "zipper",
            "-i",
            p(&input),
            "--unmapped",
            p(&unmapped),
            "--ref",
            p(&reference),
            "-o",
            p(&unmapped),
        ],
        &format!("--output '{}' is the same file as --unmapped", unmapped.display()),
        &unmapped,
    );
    // The mapped side of a zipper start is `--input`.
    assert_alias_refused(
        &[
            "runall",
            "--start-from",
            "zipper",
            "--stop-after",
            "zipper",
            "-i",
            p(&input),
            "--unmapped",
            p(&unmapped),
            "--ref",
            p(&reference),
            "-o",
            p(&input),
        ],
        &format!("--output '{}' is the same file as --input", input.display()),
        &input,
    );
    // The filter stage's own reference.
    assert_alias_refused(
        &[
            "runall",
            "--start-from",
            "filter",
            "--stop-after",
            "filter",
            "-i",
            p(&input),
            "--filter::min-reads",
            "1",
            "--filter::ref",
            p(&reference),
            "-o",
            p(&reference),
        ],
        &format!("--output '{}' is the same file as --filter::ref", reference.display()),
        &reference,
    );
}

/// A consensus metrics interval list is a read-only input like any other.
#[cfg(feature = "consensus")]
#[test]
fn write_target_aliasing_consensus_intervals_is_refused() {
    let tmp = TempDir::new().unwrap();
    let input = grouped_bam(tmp.path(), "adjacency", "iv");
    let intervals = tmp.path().join("targets.bed");
    std::fs::write(&intervals, "chr1\t0\t1000\n").unwrap();
    let metrics = tmp.path().join("m");
    assert_alias_refused(
        &[
            "runall",
            "--start-from",
            "consensus",
            "--stop-after",
            "consensus",
            "--consensus",
            "simplex",
            "-i",
            p(&input),
            "--simplex::min-reads",
            "1",
            "--simplex::metrics",
            p(&metrics),
            "--simplex::intervals",
            p(&intervals),
            "-o",
            p(&intervals),
        ],
        &format!("--output '{}' is the same file as --simplex::intervals", intervals.display()),
        &intervals,
    );
}

/// `-i -` with an unrelated, already-existing output is not an alias: the
/// stdin identity check must not refuse a legitimate run.
#[cfg(unix)]
#[test]
fn stdin_input_with_an_unrelated_existing_output_is_accepted() {
    let tmp = TempDir::new().unwrap();
    let input = unsorted_bam(tmp.path());
    let out = tmp.path().join("existing_out.bam");
    std::fs::write(&out, b"stale").unwrap();
    let output = Command::new(env!("CARGO_BIN_EXE_fgumi"))
        .args(["runall", "--start-from", "sort", "--stop-after", "sort", "-i", "-", "-o", p(&out)])
        .stdin(std::fs::File::open(&input).unwrap())
        .output()
        .unwrap();
    assert!(
        output.status.success(),
        "an unrelated existing output must be accepted:\n{}",
        String::from_utf8_lossy(&output.stderr)
    );
    let (_, records) = read_bam_output(&out);
    let (_, input_records) = read_bam_output(&input);
    assert_eq!(
        records.len(),
        input_records.len(),
        "the sorted output must hold every input record"
    );
}

/// The dead top-level `--rejects` hint names `--filter::rejects` only when
/// filter runs and that flag is unset.
#[rstest::rstest]
#[case::filter_chain_hints(&["--start-from", "filter", "--stop-after", "filter"], false, true)]
#[case::filter_rejects_already_set(&["--start-from", "filter", "--stop-after", "filter"], true, false)]
#[case::no_filter_no_hint(&["--start-from", "sort", "--stop-after", "sort"], false, false)]
fn dead_rejects_hint_names_filter_rejects_only_when_useful(
    #[case] chain: &[&str],
    #[case] filter_rejects: bool,
    #[case] expect_hint: bool,
) {
    let tmp = TempDir::new().unwrap();
    let input = unsorted_bam(tmp.path());
    let rejects = tmp.path().join("rejects.bam");
    let filter_rejects_path = tmp.path().join("filter_rejects.bam");
    let out = tmp.path().join("out.bam");
    let mut args = vec!["runall"];
    args.extend_from_slice(chain);
    args.extend(["-i", p(&input), "-o", p(&out), "--rejects", p(&rejects)]);
    if chain.contains(&"filter") {
        args.extend(["--filter::min-reads", "1"]);
    }
    if filter_rejects {
        args.extend(["--filter::rejects", p(&filter_rejects_path)]);
    }
    let stderr = String::from_utf8_lossy(&fgumi(&args).stderr).into_owned();
    assert!(stderr.contains("--rejects is wired nowhere on this runall chain"), "{stderr}");
    assert_eq!(
        stderr.contains("Use --filter::rejects to capture filter rejects."),
        expect_hint,
        "hint presence must match the chain:\n{stderr}"
    );
}

/// The per-stage read-only files are guarded too: a rejects path aliasing the
/// UMI whitelist, and an output aliasing the `--ref` FASTA the aligner reads.
#[cfg(feature = "simulate")]
#[test]
fn write_targets_aliasing_umi_files_or_ref_are_refused() {
    let tmp = TempDir::new().unwrap();
    let fixture = correct_rejects_fixture(tmp.path());
    let whitelist = tmp.path().join("umis.txt");
    std::fs::write(&whitelist, "AAAA\nCCCC\n").unwrap();
    let out = tmp.path().join("out.bam");

    let mut args = correct_chain_args("correct", &fixture);
    args.extend(["--correct::umi-files", p(&whitelist), "--correct::rejects", p(&whitelist)]);
    args.extend(["-o", p(&out)]);
    assert_alias_refused(
        &args,
        &format!(
            "--correct::rejects '{}' is the same file as --correct::umi-files",
            whitelist.display()
        ),
        &whitelist,
    );

    let reference = fixture.0.clone();
    let mut args = correct_chain_args("zipper", &fixture);
    args.extend(["-o", p(&reference)]);
    assert_alias_refused(
        &args,
        &format!("--output '{}' is the same file as --ref", reference.display()),
        &reference,
    );
}

/// `-i -` reading a redirected file is guarded against writing that same file:
/// the guard compares against whatever stdin actually is.
#[cfg(unix)]
#[test]
fn write_target_aliasing_redirected_stdin_is_refused() {
    let tmp = TempDir::new().unwrap();
    let input = unsorted_bam(tmp.path());
    let before = std::fs::read(&input).unwrap();
    let output = Command::new(env!("CARGO_BIN_EXE_fgumi"))
        .args([
            "runall",
            "--start-from",
            "sort",
            "--stop-after",
            "sort",
            "-i",
            "-",
            "-o",
            p(&input),
        ])
        .stdin(std::fs::File::open(&input).unwrap())
        .output()
        .unwrap();
    let stderr = String::from_utf8_lossy(&output.stderr);
    assert!(!output.status.success(), "-o aliasing the redirected stdin must be refused");
    assert!(
        stderr.contains(&format!("--output '{}' is the same file as --input '-'", input.display())),
        "expected the aliasing error, got:\n{stderr}"
    );
    assert_eq!(std::fs::read(&input).unwrap(), before, "the input BAM must be left untouched");
}

/// The dead `--stats` hint names `--filter::stats` only when filter runs and
/// that flag is unset.
#[rstest::rstest]
#[case::filter_chain_hints(&["--start-from", "filter", "--stop-after", "filter"], false, true)]
#[case::filter_stats_already_set(&["--start-from", "filter", "--stop-after", "filter"], true, false)]
#[case::no_filter_no_hint(&["--start-from", "sort", "--stop-after", "sort"], false, false)]
fn dead_stats_hint_names_filter_stats_only_when_useful(
    #[case] chain: &[&str],
    #[case] filter_stats: bool,
    #[case] expect_hint: bool,
) {
    let tmp = TempDir::new().unwrap();
    let input = unsorted_bam(tmp.path());
    let stats = tmp.path().join("stats.txt");
    let filter_stats_path = tmp.path().join("filter_stats.txt");
    let out = tmp.path().join("out.bam");
    let mut args = vec!["runall"];
    args.extend_from_slice(chain);
    args.extend(["-i", p(&input), "-o", p(&out), "--stats", p(&stats)]);
    if chain.contains(&"filter") {
        args.extend(["--filter::min-reads", "1"]);
    }
    if filter_stats {
        args.extend(["--filter::stats", p(&filter_stats_path)]);
    }
    let stderr = String::from_utf8_lossy(&fgumi(&args).stderr).into_owned();
    assert!(stderr.contains("--stats is consumed only by the consensus stage"), "{stderr}");
    assert_eq!(
        stderr.contains("Use --filter::stats to capture filter statistics."),
        expect_hint,
        "hint presence must match the chain:\n{stderr}"
    );
}

/// `--correct::rejects` on a chain with no correct stage is dead; the run
/// still succeeds, but must say so rather than silently write nothing.
#[test]
fn correct_rejects_without_correct_stage_warns() {
    let tmp = TempDir::new().unwrap();
    let input = unsorted_bam(tmp.path());
    let rejects = tmp.path().join("rejects.bam");
    let output = run_ok(
        [
            "runall",
            "--start-from",
            "sort",
            "--stop-after",
            "sort",
            "-i",
            p(&input),
            "-o",
            p(&tmp.path().join("out.bam")),
            "--correct::rejects",
            p(&rejects),
        ],
        "runall sort with a dead --correct::rejects",
    );
    let stderr = String::from_utf8_lossy(&output.stderr);
    assert!(
        stderr.contains("--correct::rejects is wired nowhere"),
        "expected the dead-flag warning, got:\n{stderr}"
    );
    assert!(!rejects.exists(), "no correct stage ran, so no rejects file may be written");
}

/// `runall --start-from extract` must honor `--extract::no-check-crc` on an
/// all-BGZF FASTQ file input, as standalone `fgumi extract --no-check-crc`
/// does. The chain's `verify_crc` was hardcoded `true` for a FASTQ source, so
/// the BGZF split decoder aborted with a CRC mismatch even while the run
/// logged `CRC verify: off`. The default (file ⇒ verify) must still reject.
#[test]
fn extract_honors_no_check_crc_on_bgzf_fastq() {
    // Quality-encoding detection samples through a 1 MiB buffer, so the
    // corrupted last block must sit well past the first MiB of decompressed
    // FASTQ for only the split decoder to reach it. Size the fixture by bytes
    // (4x that buffer) with long reads, keeping the record count, and so the
    // debug-build runtime, small.
    const MIN_FASTQ_BYTES: usize = 4 * 1024 * 1024;
    let tmp = TempDir::new().unwrap();
    let fastq = tmp.path().join("reads.fq.gz");
    let bases = "ACGT".repeat(38);
    let quals = "I".repeat(bases.len());
    let mut records = Vec::new();
    let mut num_records = 0;
    while records.len() < MIN_FASTQ_BYTES {
        writeln!(records, "@q{num_records}\n{bases}\n+\n{quals}").unwrap();
        num_records += 1;
    }
    write_bgzf_fastq_with_corrupt_last_crc(&fastq, &records);

    let runall = |out: &Path, extra: &[&str]| {
        let mut args = vec![
            "runall",
            "--start-from",
            "extract",
            "--stop-after",
            "extract",
            "--extract::inputs",
            p(&fastq),
            "--extract::read-structures",
            "+T",
            "--extract::sample",
            "s1",
            "--extract::library",
            "lib1",
            "-o",
            p(out),
        ];
        args.extend_from_slice(extra);
        fgumi(args)
    };

    let skipped = tmp.path().join("skipped.bam");
    let output = runall(&skipped, &["--extract::no-check-crc"]);
    assert!(
        output.status.success(),
        "--extract::no-check-crc must accept a corrupted BGZF CRC32: {}",
        String::from_utf8_lossy(&output.stderr)
    );
    let standalone = tmp.path().join("standalone.bam");
    run_ok(
        [
            "extract",
            "--inputs",
            p(&fastq),
            "--read-structures",
            "+T",
            "--sample",
            "s1",
            "--library",
            "lib1",
            "--no-check-crc",
            "-o",
            p(&standalone),
        ],
        "standalone extract --no-check-crc",
    );
    let (_, records) = read_bam_output(&skipped);
    assert_eq!(records.len(), num_records, "every record must be extracted");
    // The standalone oracle runs the same decoder, so also pin identity and
    // order against the fixture itself — including the corrupted last block.
    for (i, record) in records.iter().enumerate() {
        let name = record.name().map(ToString::to_string).unwrap_or_default();
        assert_eq!(name, format!("q{i}"), "record {i} name/order");
    }
    assert_bams_record_equivalent_nonempty(&skipped, &standalone);

    let verified = tmp.path().join("verified.bam");
    let output = runall(&verified, &[]);
    let stderr = String::from_utf8_lossy(&output.stderr).to_lowercase();
    assert!(!output.status.success(), "the default must reject a corrupted BGZF CRC32");
    // Match the mismatch error itself: a bare "crc" would also match the
    // run's own `CRC verify: on` log line and pass on any unrelated failure.
    assert!(
        stderr.contains("crc32 mismatch") || stderr.contains("checksum mismatch"),
        "the default must fail on the CRC mismatch, not some other error: {stderr}"
    );
    // The corruption sits past quality-encoding detection's window, so only the
    // BGZF split decoder can reach it — pin that it is the step enforcing the
    // default, not some earlier reader.
    assert!(
        stderr.contains("step \"fastqdecompress\" failed"),
        "the BGZF split decoder must be what rejects the corruption: {stderr}"
    );
}

/// `runall --start-from extract` over one FASTQ at `fastq` with `extra` flags;
/// returns the output path and the run's stderr.
fn run_extract_self_pair(
    dir: &Path,
    fastq: &Path,
    name: &str,
    extra: &[&str],
) -> (PathBuf, String) {
    let out = dir.join(format!("{name}.bam"));
    let mut args = vec![
        "runall",
        "--start-from",
        "extract",
        "--stop-after",
        "extract",
        "--extract::inputs",
        p(fastq),
        "--extract::read-structures",
        "+T",
        "--extract::sample",
        "s1",
        "--extract::library",
        "lib1",
        "-o",
        p(&out),
    ];
    args.extend_from_slice(extra);
    let output = run_ok(args, &format!("runall extract ({name})"));
    (out, String::from_utf8_lossy(&output.stderr).into_owned())
}

/// Both FASTQ decode fronts must honor both async-reader flags. The BGZF split
/// read only the top-level `--async-reader` and the fused readers only
/// `--extract::async-reader`, so each flag was a silent no-op on the other
/// path. Pinned through each path's own prefetch log line, plus record
/// identity with the flag off.
#[rstest::rstest]
#[case::split_per_stage_flag(true, "--extract::async-reader", "enabled on the BGZF split")]
#[case::split_top_level_flag(true, "--async-reader", "enabled on the BGZF split")]
#[case::fused_per_stage_flag(false, "--extract::async-reader", "enabled: spawning")]
#[case::fused_top_level_flag(false, "--async-reader", "enabled: spawning")]
fn extract_async_reader_flags_reach_both_fastq_paths(
    #[case] bgzf: bool,
    #[case] flag: &str,
    #[case] log_needle: &str,
) {
    let tmp = TempDir::new().unwrap();
    let fastq = tmp.path().join("reads.fq.gz");
    let records: Vec<(String, String, String)> = (0..200)
        .map(|i| (format!("q{i}"), "ACGTACGTAC".to_string(), "IIIIIIIIII".to_string()))
        .collect();
    if bgzf {
        let mut text = Vec::new();
        for (n, s, q) in &records {
            writeln!(text, "@{n}\n{s}\n+\n{q}").unwrap();
        }
        crate::helpers::fastq::write_bgzf_fastq(&fastq, &text);
    } else {
        let slices: Vec<(&str, &str, &str)> =
            records.iter().map(|(n, s, q)| (n.as_str(), s.as_str(), q.as_str())).collect();
        write_gzip_fastq(&fastq, &slices);
    }

    let (plain, plain_log) = run_extract_self_pair(tmp.path(), &fastq, "plain", &[]);
    assert!(
        !plain_log.contains("async FASTQ reader enabled"),
        "no flag, no prefetch:\n{plain_log}"
    );
    let (prefetched, log) = run_extract_self_pair(tmp.path(), &fastq, "prefetched", &[flag]);
    // Exactly one prefetch per input (one input here): `contains` alone would
    // pass if some other reader spawned the thread, or if it spawned twice.
    assert_eq!(
        log.matches(&format!("async FASTQ reader {log_needle}")).count(),
        1,
        "{flag} must reach the {} path, once per input:\n{log}",
        if bgzf { "BGZF split" } else { "fused reader" }
    );
    if bgzf {
        // On the split, the encoding-detection readers are dropped unread, so
        // they must not spawn a prefetch thread of their own.
        assert!(
            !log.contains("async FASTQ reader enabled: spawning"),
            "the split must be the only prefetch on BGZF input:\n{log}"
        );
    }
    assert_bams_record_equivalent_nonempty(&prefetched, &plain);
}

/// `runall --async-reader` now reaches a FASTQ read from stdin too (through the
/// stdin prefetch wrap), and the records match the same file read directly.
#[test]
fn extract_async_reader_prefetches_stdin_fastq() {
    let tmp = TempDir::new().unwrap();
    let fastq = tmp.path().join("reads.fq.gz");
    let records: Vec<(String, String, String)> = (0..200)
        .map(|i| (format!("q{i}"), "ACGTACGTAC".to_string(), "IIIIIIIIII".to_string()))
        .collect();
    let slices: Vec<(&str, &str, &str)> =
        records.iter().map(|(n, s, q)| (n.as_str(), s.as_str(), q.as_str())).collect();
    write_gzip_fastq(&fastq, &slices);
    let (from_file, _) = run_extract_self_pair(tmp.path(), &fastq, "from_file", &[]);

    let from_stdin = tmp.path().join("from_stdin.bam");
    let output = Command::new(env!("CARGO_BIN_EXE_fgumi"))
        .args([
            "runall",
            "--start-from",
            "extract",
            "--stop-after",
            "extract",
            "--extract::inputs",
            "-",
            "--extract::read-structures",
            "+T",
            "--extract::sample",
            "s1",
            "--extract::library",
            "lib1",
            "--async-reader",
            "-o",
            p(&from_stdin),
        ])
        .stdin(std::fs::File::open(&fastq).unwrap())
        .output()
        .unwrap();
    let log = String::from_utf8_lossy(&output.stderr);
    assert!(output.status.success(), "runall from stdin failed:\n{log}");
    assert!(
        log.contains("async FASTQ reader enabled: spawning fgumi-prefetch thread for stdin"),
        "--async-reader must reach the stdin FASTQ reader:\n{log}"
    );
    assert_eq!(
        log.matches("async FASTQ reader enabled").count(),
        1,
        "stdin must be the only prefetch:\n{log}"
    );
    assert_bams_record_equivalent_nonempty(&from_stdin, &from_file);
}

/// A BAM-source start without `-i` must still be told `--input` is missing,
/// not be misrouted into FASTQ-source (extract) option handling.
#[test]
fn bam_start_without_input_reports_missing_input() {
    let tmp = TempDir::new().unwrap();
    assert_rejected_with(
        [
            "runall",
            "--start-from",
            "sort",
            "--stop-after",
            "sort",
            "-o",
            p(&tmp.path().join("out.bam")),
        ],
        "--input is required with --start-from sort",
        "runall --start-from sort without -i",
    );
}

// ══════════════════════════ Extract→Extract (interleaved, no aligner) ══════════════════════════

/// Validates A2 (`--extract::interleaved`'s interleaved-vs-count check):
/// `runall --start-from extract --stop-after extract --extract::interleaved`
/// vs standalone `fgumi extract --interleaved`. Before A2 this would have
/// failed at the interleaved-vs-count check; it must pass now.
#[test]
fn extract_interleaved_matches_standalone_extract() {
    let tmp = TempDir::new().unwrap();
    let interleaved = tmp.path().join("interleaved.fq.gz");
    // 3 read pairs, interleaved R1,R2,R1,R2,R1,R2: R1 = 4bp UMI + 8bp
    // template (read structure `4M+T`), R2 = 8bp template only (`+T`).
    let pairs = [
        ("read0", "ACGTAAAACCCC", "GGGGTTTT"),
        ("read1", "ACGTAAAACCCC", "GGGGTTTT"),
        ("read2", "ACGTAAAACCCC", "GGGGTTTT"),
    ];
    let r1_qual = "I".repeat(12);
    let r2_qual = "I".repeat(8);
    let mut records: Vec<(&str, &str, &str)> = Vec::new();
    for (name, r1_seq, r2_seq) in &pairs {
        records.push((name, r1_seq, r1_qual.as_str()));
        records.push((name, r2_seq, r2_qual.as_str()));
    }
    write_gzip_fastq(&interleaved, &records);

    let runall_out = tmp.path().join("runall.bam");
    let staged_out = tmp.path().join("staged.bam");

    run_ok(
        [
            "runall",
            "--start-from",
            "extract",
            "--stop-after",
            "extract",
            "--extract::interleaved",
            "--extract::inputs",
            p(&interleaved),
            "--extract::read-structures",
            "4M+T",
            "+T",
            "--extract::sample",
            "s1",
            "--extract::library",
            "lib1",
            "-o",
            p(&runall_out),
        ],
        "runall extract(interleaved)->extract",
    );
    run_ok(
        [
            "extract",
            "--interleaved",
            "--inputs",
            p(&interleaved),
            "--read-structures",
            "4M+T",
            "+T",
            "--sample",
            "s1",
            "--library",
            "lib1",
            "-o",
            p(&staged_out),
        ],
        "standalone extract --interleaved",
    );

    assert_bams_record_equivalent_nonempty(&runall_out, &staged_out);
    assert_bam_headers_equivalent_ignoring_pg(&runall_out, &staged_out);
}

// ══════════════════════════ Extract fusion + the dropped guard ══════════════════════════

/// `extract→zipper` used to be rejected outright ("incompatible: zipper
/// requires `PairedBams`, but extract produces `Fastqs`"). Spec §6 dropped that
/// rule (extract now feeds Align, which is included whenever the chain
/// reaches past `Correct`, so `extract→zipper` builds `[Extract, Align]` and
/// is a normal align-bearing chain). This test pins the guard's absence
/// regardless of whether a real aligner is on `PATH`: with one, the chain
/// runs to completion; without one, it still gets past argument validation
/// and fails later for an unrelated reason (the aligner subprocess not being
/// spawnable) — the assertion that matters either way is that "incompatible"
/// never appears.
#[test]
fn extract_to_zipper_guard_is_dropped() {
    let tmp = TempDir::new().unwrap();
    let r1 = tmp.path().join("r1.fq.gz");
    write_gzip_fastq(&r1, &[("read0", "ACGTACGTACGTACGT", "IIIIIIIIIIIIIIII")]);
    let out = tmp.path().join("out.bam");
    let reference = create_test_reference(tmp.path());

    let (aligner_cmd, have_real_aligner) = if let Some(bin) = aligner_binary() {
        build_aligner_index(&reference, bin);
        (format!("{bin} mem {{ref}} /dev/stdin"), true)
    } else {
        eprintln!(
            "no bwa/bwa-mem3 on PATH: only checking the dropped guard's absence, not running \
             extract->zipper end to end"
        );
        ("definitely-not-a-real-fgumi-test-aligner mem {ref} /dev/stdin".to_string(), false)
    };

    let output = fgumi([
        "runall",
        "--start-from",
        "extract",
        "--stop-after",
        "zipper",
        "--extract::inputs",
        p(&r1),
        "--extract::read-structures",
        "4M+T",
        "--extract::sample",
        "s1",
        "--extract::library",
        "lib1",
        "--ref",
        p(&reference),
        "--aligner::command",
        &aligner_cmd,
        "-o",
        p(&out),
    ]);
    let stderr = String::from_utf8_lossy(&output.stderr);
    assert!(
        !stderr.contains("incompatible"),
        "extract->zipper must not be rejected by the removed incompatibility guard: {stderr}"
    );

    if have_real_aligner {
        assert!(
            output.status.success(),
            "extract->zipper should run to completion when a real aligner is available: {stderr}"
        );
        assert!(
            std::fs::metadata(&out).map(|m| m.len() > 0).unwrap_or(false),
            "expected a non-empty output BAM"
        );
    } else {
        // With the bogus aligner, extract->zipper must still reach the Align
        // stage and fail *there* -- running the non-existent aligner command --
        // rather than pass on an unrelated validation error. Pin that failure
        // boundary: a non-success status plus the bogus aligner's name surfaced
        // from the retained subprocess stderr. Without this, the only assertion
        // on the no-aligner path is the absence of "incompatible", which a
        // spurious earlier error would also satisfy.
        assert!(
            !output.status.success(),
            "extract->zipper with a bogus aligner must fail at the Align stage, not succeed: {stderr}"
        );
        assert!(
            stderr.contains("definitely-not-a-real-fgumi-test-aligner"),
            "expected the bogus aligner's name in the retained stderr, proving the Align stage ran \
             it (not a spurious earlier failure): {stderr}"
        );
    }
}

/// `extract→group` (no UMI): FASTQ in, through the fused `Align` stage (real
/// aligner subprocess + zipper-merge), `Sort`, and `Group` (`--no-umi`), all
/// in one `runall` invocation. Gated on a real aligner being on `PATH`.
///
/// There is no standalone `align` command to build an external staged
/// oracle from (that fusion — aligner subprocess + zipper-merge — is exactly
/// what `Stage::Align` exists to elide), so unlike the Class A/B parity
/// tests above this checks end-to-end success and header/record sanity
/// rather than byte-for-byte parity against a staged chain.
#[test]
fn extract_to_group_no_umi_runs_end_to_end() {
    let Some(bin) = aligner_binary() else {
        eprintln!("no bwa/bwa-mem3 on PATH: skipping extract->group (no-UMI) fusion test");
        return;
    };
    let tmp = TempDir::new().unwrap();
    let (reference, sequence) = write_unique_reference(tmp.path(), 2000);
    build_aligner_index(&reference, bin);
    let r1 = tmp.path().join("r1.fq.gz");
    // A 40bp substring of a 2000bp non-repetitive reference maps uniquely
    // (high MAPQ), unlike a substring of `create_test_reference`'s periodic
    // sequence — see `pseudo_random_sequence`'s doc comment. All three reads
    // share the same 40bp window, so they land in one position group.
    let read = &sequence[500..540];
    let qual = "I".repeat(read.len());
    write_gzip_fastq(
        &r1,
        &[("read0", read, &qual), ("read1", read, &qual), ("read2", read, &qual)],
    );
    let out = tmp.path().join("out.bam");
    let aligner_cmd = format!("{bin} mem {{ref}} /dev/stdin");

    run_ok(
        [
            "runall",
            "--start-from",
            "extract",
            "--stop-after",
            "group",
            "--extract::inputs",
            p(&r1),
            "--extract::read-structures",
            "+T",
            "--extract::sample",
            "s1",
            "--extract::library",
            "lib1",
            "--ref",
            p(&reference),
            "--aligner::command",
            &aligner_cmd,
            "--group::strategy",
            "identity",
            "--group::no-umi",
            "true",
            "-o",
            p(&out),
        ],
        "runall extract(no-umi)->group fused pipeline",
    );

    let (header, records) = read_bam_output(&out);
    // The fixture is fully determined: three single-end reads share one 40bp
    // window, so all three survive the fused chain and land in one position
    // group. Assert the exact count and the MI tag Group is responsible for
    // adding, not just non-emptiness.
    assert_eq!(
        records.len(),
        3,
        "all three input reads should survive the fused extract->group chain"
    );
    assert!(header.header().is_some(), "expected an @HD line in the fused output header");
    let mi_tag = fgumi_lib::sam::SamTag::MI.to_noodles_tag();
    for record in &records {
        assert!(
            record.data().get(&mi_tag).is_some(),
            "group must assign an MI tag to every emitted record"
        );
    }
}

/// `extract→duplex` (with `--correct::umis`, splicing `Correct` between
/// `Extract` and `Align`): the other extract-fusion path the design doc
/// calls for, alongside `extract→group` above. Unlike `extract→group`, this
/// needs a genuinely duplex-shaped input — both strands of each molecule
/// present — for `--group::strategy paired` to assign `/A`/`/B` reads and
/// for `duplex` to have real two-strand input to consensus-call, which
/// [`write_duplex_umi_fastq`] builds. Aligner-gated like `extract→group`.
///
/// Unlike the other Class B/extract-fusion tests, the staged oracle here is
/// built from EVERY individual stage as its own `fgumi` subprocess —
/// `extract | correct | <real aligner> | zipper | sort | group | duplex` —
/// because the brief specifically asks for that fully-staged comparison.
/// There is no standalone `fgumi align` command (that fusion is exactly what
/// `Stage::Align` exists to elide), so the "align" step here is a raw
/// `bwa`/`bwa-mem3` subprocess reading [`bam_to_interleaved_fastq`]'s output,
/// piped straight into `fgumi zipper`.
#[cfg(feature = "consensus")]
#[test]
fn runall_extract_to_duplex_fusion_matches_staged() {
    let Some(bin) = aligner_binary() else {
        eprintln!("no bwa/bwa-mem3 on PATH: skipping extract->duplex fusion test");
        return;
    };
    let tmp = TempDir::new().unwrap();
    let (reference, sequence) = write_unique_reference(tmp.path(), 4000);
    build_aligner_index(&reference, bin);

    // Two duplex molecules (both strands present for each), 2 read-pairs per
    // strand, at well-separated offsets in the 4000bp reference so neither
    // molecule's anchors overlap the other's.
    let molecules = [(500usize, "AAAA", "CCCC"), (2000usize, "GGGG", "TTTT")];
    let (r1, r2) = write_duplex_umi_fastq(tmp.path(), &sequence, &molecules);
    // `-p`: the fused Align stage interleaves paired-end reads into one
    // stdin stream, so the aligner command must tell bwa to treat that
    // stream as paired (interleaved) rather than single-end — unlike the
    // single-end `extract_to_group_no_umi_runs_end_to_end` /
    // `extract_to_zipper_guard_is_dropped` tests, which never hit this
    // because they only ever feed one FASTQ file.
    let aligner_cmd = format!("{bin} mem -p {{ref}} /dev/stdin");

    // ── Fused: runall --start-from extract --stop-after consensus --consensus duplex. ──
    let runall_out = tmp.path().join("runall.bam");
    run_ok(
        [
            "runall",
            "--start-from",
            "extract",
            "--stop-after",
            "consensus",
            "--consensus",
            "duplex",
            "--extract::inputs",
            p(&r1),
            p(&r2),
            "--extract::read-structures",
            "4M+T",
            "4M+T",
            "--extract::sample",
            "s1",
            "--extract::library",
            "lib1",
            "--correct::umis",
            "AAAA",
            "--correct::umis",
            "CCCC",
            "--correct::umis",
            "GGGG",
            "--correct::umis",
            "TTTT",
            "--correct::min-distance",
            "1",
            "--ref",
            p(&reference),
            "--aligner::command",
            &aligner_cmd,
            "--group::strategy",
            "paired",
            "--group::edits",
            "0",
            "-o",
            p(&runall_out),
        ],
        "runall extract->duplex fusion",
    );

    let staged_out = run_staged_extract_to_duplex(tmp.path(), &r1, &r2, &reference, bin);

    assert_bams_record_equivalent_nonempty(&runall_out, &staged_out);
    assert_bam_headers_equivalent_ignoring_pg(&runall_out, &staged_out);
}

/// The fully-staged oracle for [`runall_extract_to_duplex_fusion_matches_staged`]:
/// `extract | correct | <real aligner> | zipper | sort | group | duplex`, each
/// its own `fgumi` subprocess (or, for the aligner step, a raw `bwa`/`bwa-mem3`
/// subprocess — there is no standalone `fgumi align` command). Returns the
/// final duplex-consensus BAM's path.
fn run_staged_extract_to_duplex(
    dir: &Path,
    r1: &Path,
    r2: &Path,
    reference: &Path,
    aligner_bin: &str,
) -> PathBuf {
    let extracted = dir.join("staged_extracted.bam");
    run_ok(
        [
            "extract",
            "--inputs",
            p(r1),
            p(r2),
            "--read-structures",
            "4M+T",
            "4M+T",
            "--sample",
            "s1",
            "--library",
            "lib1",
            "-o",
            p(&extracted),
        ],
        "staged extract",
    );

    let corrected = dir.join("staged_corrected.bam");
    run_ok(
        [
            "correct",
            "-i",
            p(&extracted),
            "-o",
            p(&corrected),
            "--umis",
            "AAAA",
            "--umis",
            "CCCC",
            "--umis",
            "GGGG",
            "--umis",
            "TTTT",
            "--min-distance",
            "1",
        ],
        "staged correct",
    );

    let interleaved = bam_to_interleaved_fastq(&corrected);
    let mut aligner = Command::new(aligner_bin)
        .args(["mem", "-p", p(reference), "-"])
        .stdin(Stdio::piped())
        .stdout(Stdio::piped())
        .stderr(Stdio::null())
        .spawn()
        .expect("failed to spawn staged aligner");
    aligner
        .stdin
        .take()
        .expect("aligner stdin was piped")
        .write_all(&interleaved)
        .expect("write interleaved FASTQ to aligner stdin");
    let aligner_output = aligner.wait_with_output().expect("wait on staged aligner");
    assert!(aligner_output.status.success(), "staged aligner run failed");
    let mapped_sam = dir.join("staged_mapped.sam");
    std::fs::write(&mapped_sam, &aligner_output.stdout).expect("write staged mapped SAM");

    let zipped = dir.join("staged_zipped.bam");
    run_ok(
        [
            "zipper",
            "-i",
            p(&mapped_sam),
            "-u",
            p(&corrected),
            "--reference",
            p(reference),
            "-o",
            p(&zipped),
        ],
        "staged zipper",
    );

    let staged_sorted = dir.join("staged_sorted.bam");
    run_ok(["sort", "-i", p(&zipped), "-o", p(&staged_sorted)], "staged sort");

    let staged_grouped = dir.join("staged_grouped.bam");
    run_ok(
        [
            "group",
            "-i",
            p(&staged_sorted),
            "-o",
            p(&staged_grouped),
            "--strategy",
            "paired",
            "--edits",
            "0",
        ],
        "staged group",
    );

    let staged_out = dir.join("staged.bam");
    run_ok(["duplex", "-i", p(&staged_grouped), "-o", p(&staged_out)], "staged duplex");
    staged_out
}

// ══════════════════════════════════ determinism ══════════════════════════════════

#[cfg(feature = "consensus")]
#[test]
fn runall_record_stream_is_deterministic() {
    let tmp = TempDir::new().unwrap();
    let fixture = sorted_bam(tmp.path());
    let out1 = tmp.path().join("run1.bam");
    let out2 = tmp.path().join("run2.bam");

    for out in [&out1, &out2] {
        run_ok(
            [
                "runall",
                "--start-from",
                "group",
                "--stop-after",
                "consensus",
                "--consensus",
                "simplex",
                "-i",
                p(&fixture),
                "-o",
                p(out),
                "--group::strategy",
                "identity",
                "--group::edits",
                "0",
                "--simplex::min-reads",
                "1",
            ],
            "runall determinism run",
        );
    }

    let (_, records1) = read_bam_output(&out1);
    let (_, records2) = read_bam_output(&out2);
    assert!(!records1.is_empty(), "determinism fixture produced no records");
    assert_eq!(
        records1, records2,
        "runall's record stream is not deterministic across repeated runs"
    );
}

// ══════════════════════════════════ stdin-once ══════════════════════════════════

/// `--input -` must be consumed exactly once: a double-read (or a read that
/// starts partway through) would corrupt the BGZF stream and either fail
/// outright or silently drop/duplicate records, which would show up here as
/// a record-stream mismatch against the equivalent file-input run.
#[test]
fn runall_reads_bam_from_stdin_once() {
    let tmp = TempDir::new().unwrap();
    let fixture = unsorted_bam(tmp.path());
    let file_out = tmp.path().join("from_file.bam");
    let stdin_out = tmp.path().join("from_stdin.bam");

    run_ok(
        [
            "runall",
            "--start-from",
            "sort",
            "--stop-after",
            "sort",
            "-i",
            p(&fixture),
            "-o",
            p(&file_out),
        ],
        "runall via file input",
    );

    let mut cat = Command::new("cat")
        .arg(&fixture)
        .stdout(Stdio::piped())
        .spawn()
        .expect("failed to spawn cat");
    let cat_stdout = cat.stdout.take().expect("cat stdout was piped");
    let status = Command::new(env!("CARGO_BIN_EXE_fgumi"))
        .args([
            "runall",
            "--start-from",
            "sort",
            "--stop-after",
            "sort",
            "-i",
            "-",
            "-o",
            p(&stdin_out),
        ])
        .stdin(Stdio::from(cat_stdout))
        .status()
        .expect("failed to spawn runall with piped stdin");
    let _ = cat.wait();
    assert!(status.success(), "runall via stdin failed");

    assert_bams_record_equivalent_nonempty(&file_out, &stdin_out);
}

// ══════════════════════════════════ --help smoke ══════════════════════════════════

#[test]
fn help_lists_key_flags() {
    let output = fgumi(["runall", "--help"]);
    assert!(output.status.success(), "`fgumi runall --help` should exit 0");
    let stdout = String::from_utf8_lossy(&output.stdout);
    for needle in ["--start-from", "--stop-after", "--sort::max-memory", "--group::strategy"] {
        assert!(stdout.contains(needle), "`runall --help` output is missing {needle:?}:\n{stdout}");
    }
}

// ══════════════════════════════════ error paths ══════════════════════════════════

#[test]
fn rejects_backwards_group_to_sort() {
    assert_rejected_with(
        [
            "runall",
            "--start-from",
            "group",
            "--stop-after",
            "sort",
            "-i",
            "nonexistent-input.bam",
            "-o",
            "out.bam",
        ],
        "comes after --stop-after sort",
        "backwards group->sort",
    );
}

/// `--stop-after align` is rejected at the clap parse layer (it is not a
/// [`StopAfter`] value), before any runtime validation runs.
#[test]
fn rejects_stop_after_align_at_clap_level() {
    let output = fgumi([
        "runall",
        "--start-from",
        "align",
        "--stop-after",
        "align",
        "-i",
        "x",
        "-o",
        "out.bam",
    ]);
    assert!(!output.status.success(), "clap should reject --stop-after align");
    let stderr = String::from_utf8_lossy(&output.stderr);
    assert!(
        stderr.contains("invalid value 'align'"),
        "expected clap's invalid-value error for --stop-after align, got:\n{stderr}"
    );
    // Pin the rejected argument: `--start-from align` is itself valid (it maps to
    // `RunAllStage::AlignAndMerge`, clap name "align"), so the invalid-value error
    // must name `--stop-after` -- otherwise the test could pass on a rejection of
    // the wrong argument.
    assert!(
        stderr.contains("--stop-after"),
        "the invalid-value error must name --stop-after (not --start-from, which accepts 'align'), \
         got:\n{stderr}"
    );
}

#[test]
fn rejects_zipper_without_unmapped() {
    let tmp = TempDir::new().unwrap();
    let input = tmp.path().join("mapped.bam");
    write_bam(&input, &create_minimal_header("chr1", 1000), &[]);

    assert_rejected_with(
        [
            "runall",
            "--start-from",
            "zipper",
            "--stop-after",
            "zipper",
            "-i",
            p(&input),
            "-o",
            "out.bam",
        ],
        "--unmapped is required with --start-from zipper",
        "zipper without --unmapped",
    );
}

#[test]
fn rejects_extract_to_correct_without_umi_source() {
    assert_rejected_with(
        ["runall", "--start-from", "extract", "--stop-after", "correct", "-o", "out.bam"],
        "--correct::umi-files or --correct::umis",
        "extract->correct without a UMI source",
    );
}

#[cfg(feature = "consensus")]
#[test]
fn rejects_codec_with_methylation_mode() {
    let tmp = TempDir::new().unwrap();
    let input = tmp.path().join("in.bam");
    write_bam(&input, &create_minimal_header("chr1", 1000), &[]);

    assert_rejected_with(
        [
            "runall",
            "--start-from",
            "group",
            "--stop-after",
            "consensus",
            "--consensus",
            "codec",
            "--methylation-mode",
            "em-seq",
            "-i",
            p(&input),
            "-o",
            "out.bam",
            "--group::strategy",
            "identity",
        ],
        "not supported with codec consensus",
        "codec + --methylation-mode",
    );
}

/// Spec §7: `--consensus <simplex|duplex|codec>` is required whenever the
/// derived chain reaches the consensus stage. `derive_stages_for` returns this
/// error before the input BAM is ever opened (`nonexistent-input.bam` is never
/// read), mirroring `rejects_duplex_without_paired_strategy` below.
#[cfg(feature = "consensus")]
#[test]
fn rejects_missing_consensus_mode_when_chain_reaches_consensus() {
    assert_rejected_with(
        [
            "runall",
            "--start-from",
            "group",
            "--stop-after",
            "consensus",
            "-i",
            "nonexistent-input.bam",
            "-o",
            "out.bam",
            "--group::strategy",
            "adjacency",
        ],
        "--consensus <simplex|duplex|codec> is required",
        "group->consensus without --consensus",
    );
}

#[cfg(feature = "consensus")]
#[test]
fn rejects_duplex_without_paired_strategy() {
    assert_rejected_with(
        [
            "runall",
            "--start-from",
            "group",
            "--stop-after",
            "consensus",
            "--consensus",
            "duplex",
            "-i",
            "nonexistent-input.bam",
            "-o",
            "out.bam",
            "--group::strategy",
            "identity",
            "--duplex::min-reads",
            "1",
        ],
        "--consensus duplex requires --strategy paired",
        "duplex without --strategy paired",
    );
}

#[test]
fn rejects_ref_without_methylation_mode_on_non_align_chain() {
    assert_rejected_with(
        [
            "runall",
            "--start-from",
            "sort",
            "--stop-after",
            "sort",
            "-i",
            "nonexistent-input.bam",
            "-o",
            "out.bam",
            "--ref",
            "/nonexistent/ref.fa",
        ],
        "--ref requires --methylation-mode to be set",
        "--ref without --methylation-mode on a non-align chain",
    );
}

/// Bonus coverage beyond the required error-path list: `correct→sort`
/// always chains through the align stage (a corrected unmapped BAM has no
/// query-coordinate BAM to sort directly), so it hits the align-stage
/// `--ref` guard rather than the generic "`--ref` requires
/// `--methylation-mode`" one above.
#[test]
fn rejects_correct_to_sort_without_ref() {
    assert_rejected_with(
        [
            "runall",
            "--start-from",
            "correct",
            "--stop-after",
            "sort",
            "-i",
            "nonexistent-input.bam",
            "-o",
            "out.bam",
            "--correct::umis",
            "ACGT",
            "--correct::min-distance",
            "1",
        ],
        "a runall chain that includes align requires --ref",
        "correct->sort without --ref",
    );
}

/// Pins the deliberate, documented absence of `--simplex::allow-unmapped`
/// (spec §12.2 option (b) — see the design doc): `runall` does not expose
/// it, so clap rejects it as an unknown argument rather than silently
/// accepting and ignoring it.
#[cfg(feature = "consensus")]
#[test]
fn simplex_allow_unmapped_is_not_exposed() {
    let output = fgumi([
        "runall",
        "--start-from",
        "group",
        "--stop-after",
        "consensus",
        "--consensus",
        "simplex",
        "--simplex::allow-unmapped",
        "-i",
        "nonexistent-input.bam",
        "-o",
        "out.bam",
        "--group::strategy",
        "identity",
        "--simplex::min-reads",
        "1",
    ]);
    assert!(!output.status.success(), "--simplex::allow-unmapped should be rejected by clap");
    let stderr = String::from_utf8_lossy(&output.stderr);
    assert!(
        stderr.contains("unexpected argument"),
        "expected clap to reject --simplex::allow-unmapped as unknown, got:\n{stderr}"
    );
}
