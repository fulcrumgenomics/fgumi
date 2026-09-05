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
/// `resolve_methylation_mode`/`consensus_reference` in
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
