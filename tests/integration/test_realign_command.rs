//! Integration tests for `fgumi realign`.
//!
//! Most tests need no real aligner: an `awk` command stands in for it, reporting
//! every read unmapped (the `rw` input check, the per-end counts, the intermediate
//! filter), placing a pair at fixed positions (per-base tag orientation, the duplex
//! chain), or listing an in-window `XA` alternative (a swap through the whole
//! pipeline). These run in every CI job.
//!
//! `runall_matches_standalone_chain` runs the same input through
//! `runall --realign` and through that standalone chain and compares the
//! records; `runall_from_extract_realigns_with_the_aligned_header` runs
//! `runall --start-from extract --realign`, where group reads its contig names
//! from the align stage's header. Both are gated like the standalone e2e.
//!
//! `standalone_chain_aligns_consensus` runs the standalone chain
//! `simulate mapped-reads → group --realign-window → simplex → realign → sort`
//! with a real aligner. Like the `runall` align-stage tests, it runs end to end
//! only when `bwa-mem3` is on `PATH` and skips (passes) otherwise; with
//! `FGUMI_BWA_MEM3_REQUIRE_TOOLS` set (as CI's e2e-parity job does) a missing
//! `bwa-mem3` fails instead. The
//! in-process case additionally needs the `aligner-bwa-mem3` feature. Its input
//! and the parity test's include mate-unmapped molecules (carried through group
//! and consensus), and it removes `rw` from some consensus R2s, which realign
//! must align and count as without a window.

use std::path::Path;
use std::process::Command;

use fgumi_lib::sam::SamTag;
use fgumi_raw_bam::{SamBuilder, flags};
use tempfile::TempDir;

use crate::helpers::aligner::{bwa_mem3_or_skip, bwa_mem3_test_ref};
use crate::helpers::bam_generator::{create_minimal_header, create_test_reference, write_bam};
use crate::helpers::cli::utf8_path as p;

fn fgumi() -> Command {
    Command::new(env!("CARGO_BIN_EXE_fgumi"))
}

/// A stand-in aligner for the window-check tests: reads the interleaved FASTQ
/// realign feeds it and reports every read unmapped (flag 77 for R1, 141 for
/// R2), so the run completes without a real aligner or index. Its `@SQ` is that
/// of [`create_test_reference`].
pub(crate) const UNMAPPED_ALIGNER: &str = "awk 'BEGIN { OFS = \"\\t\"; print \"@HD\", \"VN:1.6\"; \
    print \"@SQ\", \"SN:chr1\", \"LN:10000\" } \
    NR % 4 == 1 { name = substr($1, 2) } NR % 4 == 2 { seq = $0 } \
    NR % 4 == 0 { n++; print name, (n % 2 ? 77 : 141), \"*\", 0, 0, \"*\", \"*\", 0, 0, seq, $0 }' \
    # {ref}";

/// An unmapped paired consensus end (with the per-read consensus tags `cD`/`cE`
/// filter requires), with `rw` if given.
fn consensus_read(name: &str, seg: u16, rw: Option<&str>) -> fgumi_raw_bam::RawRecord {
    let mut b = SamBuilder::new();
    b.read_name(name.as_bytes())
        .flags(flags::PAIRED | flags::UNMAPPED | flags::MATE_UNMAPPED | seg)
        .sequence(b"ACGTACGTAC")
        .qualities(&[30; 10]);
    if let Some(rw) = rw {
        b.add_string_tag(SamTag::RW, rw.as_bytes());
    }
    b.add_string_tag(SamTag::MI, name.as_bytes())
        .add_int_tag(SamTag::CD, 5)
        .add_float_tag(SamTag::CE, 0.0);
    b.build()
}

/// `n` consensus pairs `cons:<i>`; `rw(i)` gives each pair's `(R1, R2)` windows.
fn write_consensus(
    dir: &Path,
    n: usize,
    rw: impl Fn(usize) -> (Option<&'static str>, Option<&'static str>),
) -> std::path::PathBuf {
    let header = create_minimal_header("chr1", 10_000);
    let records: Vec<_> = (0..n)
        .flat_map(|i| {
            let name = format!("cons:{i}");
            let (rw1, rw2) = rw(i);
            [
                consensus_read(&name, flags::FIRST_SEGMENT, rw1),
                consensus_read(&name, flags::LAST_SEGMENT, rw2),
            ]
        })
        .collect();
    let input = dir.join("consensus.bam");
    write_bam(&input, &header, &records);
    input
}

/// Runs `fgumi realign` on `input` with [`UNMAPPED_ALIGNER`], returning success
/// and stderr.
fn realign_with_unmapped_aligner(dir: &Path, input: &Path) -> (bool, String) {
    let reference = create_test_reference(dir);
    let output = dir.join("out.bam");
    let result = fgumi()
        .args(["realign", "-i", p(input), "-o", p(&output), "-r", p(&reference)])
        .args(["--aligner::command", UNMAPPED_ALIGNER])
        .env("RUST_LOG", "info")
        .output()
        .expect("spawn fgumi realign");
    (result.status.success(), String::from_utf8_lossy(&result.stderr).into_owned())
}

/// An input that was never window-grouped (no `rw` on any read) fails, naming
/// the first read and how to produce windows: a short input once all of it has
/// been read, a long one once its first 1000 templates have no `rw`. Either way
/// it fails before any read reaches the aligner, so the output holds no
/// records.
#[rstest::rstest]
#[case::short_input(4, "none of the first 4 ")]
#[case::past_the_probe(1500, "none of the first 1000 ")]
fn realign_input_without_rw_fails(#[case] templates: usize, #[case] checked: &str) {
    let tmp = TempDir::new().expect("tempdir");
    let input = write_consensus(tmp.path(), templates, |_| (None, None));
    let (ok, stderr) = realign_with_unmapped_aligner(tmp.path(), &input);
    assert!(!ok, "realign must fail on an input without rw:\n{stderr}");
    assert!(stderr.contains("'cons:0'"), "the error must name the first read:\n{stderr}");
    assert!(stderr.contains(checked), "the error must say how many reads it checked:\n{stderr}");
    assert!(stderr.contains("rw tag"), "the error must name the missing tag:\n{stderr}");
    assert!(stderr.contains("--realign-window"), "the error must say how to fix it:\n{stderr}");
    assert_eq!(output_records(&tmp.path().join("out.bam")), 0, "no read may be written");
}

/// The records in the BAM at `path`; 0 when nothing (or no complete header)
/// was written, as when a run fails before its output header.
fn output_records(path: &Path) -> usize {
    let Ok(file) = std::fs::File::open(path) else {
        return 0;
    };
    let mut reader = noodles::bam::io::Reader::new(std::io::BufReader::new(file));
    if reader.read_header().is_err() {
        return 0;
    }
    reader.records().collect::<std::io::Result<Vec<_>>>().expect("written records decode").len()
}

/// A window-grouped input whose consensus ends lack `rw` where the raw reads
/// never mapped (here: every R2, and one whole pair) realigns; those ends are
/// counted as without a window.
#[test]
fn realign_ends_without_rw_are_counted() {
    let tmp = TempDir::new().expect("tempdir");
    let input = write_consensus(tmp.path(), 3, |i| {
        if i == 1 { (None, None) } else { (Some("chr1:1-100:+"), None) }
    });
    let (ok, stderr) = realign_with_unmapped_aligner(tmp.path(), &input);
    assert!(ok, "realign must accept ends without rw:\n{stderr}");
    assert_eq!(realign_count(&stderr, "Consensus reads realigned:"), 6);
    assert_eq!(realign_count(&stderr, "Unmapped:"), 2, "the windowed R1s");
    assert_eq!(realign_count(&stderr, "Without a window:"), 4);
    // Every end is written, with its consensus tags and without rw.
    let want: Vec<_> = (0..3)
        .flat_map(|i| {
            [(format!("cons:{i}"), false, false, true), (format!("cons:{i}"), true, false, true)]
        })
        .collect();
    let mut want = want;
    want.sort();
    assert_eq!(primary_ends(&tmp.path().join("out.bam")), want);
}

/// A stand-in aligner that places R1 outside its window (chr1:5001, forward)
/// with an in-window `XA` alternative (chr1:101), and R2 inside its window
/// (chr1:121, reverse), so realign keeps R2 and swaps R1. It echoes SEQ/QUAL,
/// so the reads must be their own reverse complement.
const SWAPPING_ALIGNER: &str = "awk 'BEGIN { OFS = \"\\t\"; print \"@HD\", \"VN:1.6\"; \
    print \"@SQ\", \"SN:chr1\", \"LN:10000\" } \
    NR % 4 == 1 { name = substr($1, 2) } NR % 4 == 2 { seq = $0 } \
    NR % 4 == 0 { n++; if (n % 2) print name, 97, \"chr1\", 5001, 60, \"16M\", \"=\", 121, 0, seq, $0, \"NM:i:0\", \"AS:i:16\", \"XA:Z:chr1,+101,16M,0;\"; \
    else print name, 145, \"chr1\", 121, 60, \"16M\", \"=\", 5001, 0, seq, $0, \"NM:i:0\", \"AS:i:16\" }' \
    # {ref}";

/// One output end of [`realign_swaps_to_an_in_window_xa_hit`]: (first segment,
/// position, MAPQ, proper pair, mate position, `XA`, `XS`).
type SwapOutcome = (bool, usize, u8, bool, usize, Option<String>, Option<i64>);

/// Kept and Swapped through the real pipeline with no real aligner: R1, placed
/// outside its window with an in-window `XA` hit, is promoted to it (MAPQ 0, the
/// old placement moved into `XA`, `XS` its old score), R2 stays, and the mate
/// fields agree. bwa had not flagged the pair proper, and a swap never sets the
/// flag (the windows are no insert-size check), so it stays improper.
#[test]
fn realign_swaps_to_an_in_window_xa_hit() {
    let tmp = TempDir::new().expect("tempdir");
    let dir = tmp.path();
    let input = write_swap_input(dir);
    let reference = create_test_reference(dir);
    let output = dir.join("out.bam");
    let result = fgumi()
        .args(["realign", "-i", p(&input), "-o", p(&output), "-r", p(&reference)])
        .args(["--aligner::command", SWAPPING_ALIGNER])
        .env("RUST_LOG", "info")
        .output()
        .expect("spawn fgumi realign");
    let stderr = String::from_utf8_lossy(&result.stderr);
    assert!(result.status.success(), "{stderr}");
    assert!(stderr.contains("Primary already in window: 1"), "{stderr}");
    assert!(stderr.contains("Primary swapped to an in-window XA hit: 1"), "{stderr}");
    assert_swap_outcome(&output);
}

/// How [`realign_reference_modes`] gives realign its reference.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
enum ReferenceMode {
    /// The default, with a `.fai` beside the FASTA: read on demand.
    Indexed,
    /// The default without a `.fai`: falls back to holding it in memory.
    NoFai,
    /// `--in-memory-reference`.
    InMemory,
    /// A BGZF FASTA with the `.fai` (and `.gzi`) `samtools faidx` writes for it;
    /// `true` with `--in-memory-reference`.
    Bgzf(bool),
    /// A plain gzip FASTA (no index can describe it).
    Gzip,
}

/// Compresses the FASTA at `reference` (BGZF, or plain gzip) into
/// `ref.fa.gz` beside it, with its `.dict`. A BGZF one is indexed with
/// `samtools faidx` when `samtools` is on `PATH` (writing `.fai` and `.gzi`);
/// otherwise it gets the plain FASTA's `.fai`, which holds the same
/// uncompressed offsets.
fn compressed_reference(reference: &Path, bgzf: bool) -> std::path::PathBuf {
    use std::io::Write;
    let dir = reference.parent().expect("reference dir");
    let text = std::fs::read(reference).expect("read reference");
    let path = dir.join("ref.fa.gz");
    let file = std::fs::File::create(&path).expect("create");
    if bgzf {
        let mut writer = noodles::bgzf::io::Writer::new(file);
        writer.write_all(&text).expect("write");
        writer.finish().expect("finish");
        let indexed = which::which("samtools").is_ok()
            && Command::new("samtools")
                .args(["faidx", p(&path)])
                .status()
                .is_ok_and(|status| status.success());
        if indexed {
            assert!(dir.join("ref.fa.gz.gzi").exists(), "samtools faidx writes a .gzi");
        } else {
            std::fs::copy(dir.join("ref.fa.fai"), dir.join("ref.fa.gz.fai")).expect("fai");
        }
    } else {
        let mut gz = flate2::write::GzEncoder::new(file, flate2::Compression::fast());
        gz.write_all(&text).expect("write");
        gz.finish().expect("finish");
    }
    std::fs::copy(dir.join("ref.dict"), dir.join("ref.fa.dict")).expect("dict");
    path
}

/// The reference is read on demand by default, falls back to memory with a
/// warning naming why and the memory cost (a compressed FASTA is decompressed,
/// never read at a `.fai`'s uncompressed offsets), and is held in memory with
/// `--in-memory-reference`. Every mode promotes the same end and writes the
/// same records (the promotion regenerates `MD` from the reference).
#[rstest::rstest]
#[case::default_indexed(ReferenceMode::Indexed, "on demand through its .fai index")]
#[case::fallback_without_fai(ReferenceMode::NoFai, "no .fai index found at")]
#[case::opt_out(ReferenceMode::InMemory, "in-memory reference requested")]
#[case::bgzf_with_fai(ReferenceMode::Bgzf(false), "is compressed (gzip or BGZF)")]
#[case::bgzf_in_memory(ReferenceMode::Bgzf(true), "in-memory reference requested")]
#[case::gzip(ReferenceMode::Gzip, "is compressed (gzip or BGZF)")]
fn realign_reference_modes(#[case] mode: ReferenceMode, #[case] logged: &str) {
    let tmp = TempDir::new().expect("tempdir");
    let dir = tmp.path();
    let input = write_swap_input(dir);
    let mut reference = create_test_reference(dir);
    let fai = dir.join("ref.fa.fai");
    let mut extra: Vec<&str> = vec![];
    if mode == ReferenceMode::NoFai {
        std::fs::remove_file(&fai).expect("remove .fai");
    }
    if let ReferenceMode::InMemory | ReferenceMode::Bgzf(true) = mode {
        extra.extend(["--in-memory-reference", "true"]);
    }
    if let ReferenceMode::Bgzf(_) | ReferenceMode::Gzip = mode {
        reference = compressed_reference(&reference, mode != ReferenceMode::Gzip);
    }
    let output = dir.join("out.bam");
    let result = fgumi()
        .args(["realign", "-i", p(&input), "-o", p(&output), "-r", p(&reference)])
        .args(["--aligner::command", SWAPPING_ALIGNER])
        .args(&extra)
        .env("RUST_LOG", "info")
        .output()
        .expect("spawn fgumi realign");
    let stderr = String::from_utf8_lossy(&result.stderr);
    assert!(result.status.success(), "{stderr}");
    assert!(stderr.contains(logged), "{stderr}");
    let warnings = stderr.lines().filter(|l| l.contains("WARN") && l.contains("Realign:")).count();
    match mode {
        ReferenceMode::NoFai => {
            assert_eq!(warnings, 1, "one fallback warning:\n{stderr}");
            assert!(stderr.contains(&format!("{}.fai", reference.display())), "{stderr}");
            assert!(stderr.contains("loaded into memory instead"), "{stderr}");
            assert!(stderr.contains("index it with `samtools faidx`"), "{stderr}");
            assert!(!stderr.contains("decompress it"), "nothing to decompress:\n{stderr}");
            assert!(!stderr.contains("on demand through"), "{stderr}");
        }
        ReferenceMode::Bgzf(false) | ReferenceMode::Gzip => {
            assert_eq!(warnings, 1, "one fallback warning:\n{stderr}");
            assert!(stderr.contains("decompressed and loaded into memory instead"), "{stderr}");
            assert!(stderr.contains("decompress it and index it"), "{stderr}");
            assert!(!stderr.contains("on demand through"), "{stderr}");
        }
        ReferenceMode::Indexed | ReferenceMode::InMemory | ReferenceMode::Bgzf(true) => {
            assert_eq!(warnings, 0, "no fallback warning:\n{stderr}");
        }
    }
    assert!(stderr.contains("Primary swapped to an in-window XA hit: 1"), "{stderr}");
    assert_swap_outcome(&output);
}

/// The single consensus pair of the swap tests: R1 and R2 with the windows
/// [`SWAPPING_ALIGNER`] places them outside and inside.
fn write_swap_input(dir: &Path) -> std::path::PathBuf {
    let records: Vec<_> =
        [(flags::FIRST_SEGMENT, MAPPED_PAIR_RW.0), (flags::LAST_SEGMENT, MAPPED_PAIR_RW.1)]
            .into_iter()
            .map(|(seg, rw)| {
                let mut b = SamBuilder::new();
                b.read_name(b"cons:0")
                    .flags(flags::PAIRED | flags::UNMAPPED | flags::MATE_UNMAPPED | seg)
                    .sequence(b"AAAACCCCGGGGTTTT")
                    .qualities(&[30; 16]);
                b.add_string_tag(SamTag::RW, rw.as_bytes())
                    .add_string_tag(SamTag::MI, b"0")
                    .add_int_tag(SamTag::CD, 5)
                    .add_float_tag(SamTag::CE, 0.0);
                b.build()
            })
            .collect();
    let input = dir.join("consensus.bam");
    write_bam(&input, &create_minimal_header("chr1", 10_000), &records);
    input
}

/// R1 swapped into its window and R2 kept, with the mate fields agreeing, and
/// the promoted R1's `MD` regenerated from the reference.
fn assert_swap_outcome(output: &Path) {
    let mut reader = noodles::bam::io::Reader::new(std::io::BufReader::new(
        std::fs::File::open(output).expect("output"),
    ));
    let header = reader.read_header().expect("header");
    let xa = SamTag::XA.to_noodles_tag();
    let xs = SamTag::XS.to_noodles_tag();
    let mut got: Vec<SwapOutcome> = reader
        .record_bufs(&header)
        .map(|r| {
            let r = r.expect("record");
            let string = |tag| match r.data().get(&tag) {
                Some(noodles::sam::alignment::record_buf::data::field::Value::String(v)) => {
                    Some(v.to_string())
                }
                _ => None,
            };
            (
                r.flags().is_first_segment(),
                r.alignment_start().expect("mapped").get(),
                r.mapping_quality().map_or(255, u8::from),
                r.flags().is_properly_segmented(),
                r.mate_alignment_start().expect("mate mapped").get(),
                string(xa),
                r.data()
                    .get(&xs)
                    .and_then(noodles::sam::alignment::record_buf::data::field::Value::as_int),
            )
        })
        .collect();
    got.sort_by_key(|g| !g.0);
    assert_eq!(
        got,
        vec![
            (true, 101, 0, false, 121, Some("chr1,+5001,16M,0;".to_string()), Some(16)),
            (false, 121, 60, false, 101, None, None),
        ]
    );
    // The promoted R1 (`AAAACCCCGGGGTTTT` at chr1:101, where the reference
    // repeats `ACGT`) gets its MD from the reference.
    let mut reader = noodles::bam::io::Reader::new(std::io::BufReader::new(
        std::fs::File::open(output).expect("output"),
    ));
    let header = reader.read_header().expect("header");
    let md: Vec<String> = reader
        .record_bufs(&header)
        .map(|r| r.expect("record"))
        .filter(|r| r.flags().is_first_segment())
        .filter_map(|r| match r.data().get(&SamTag::MD.to_noodles_tag()) {
            Some(noodles::sam::alignment::record_buf::data::field::Value::String(v)) => {
                Some(v.to_string())
            }
            _ => None,
        })
        .collect();
    assert_eq!(md, ["1C0G0T0A1G0T0A0C1T0A0C0G1"]);
}

/// A stand-in aligner for [`duplex_chain_realigns_into_shared_windows`]: R1 forward
/// at chr1:101, R2 reverse at chr1:201, both `16M`, echoing SEQ/QUAL.
const DUPLEX_PAIR_ALIGNER: &str = "awk 'BEGIN { OFS = \"\\t\"; print \"@HD\", \"VN:1.6\"; \
    print \"@SQ\", \"SN:chr1\", \"LN:10000\" } \
    NR % 4 == 1 { name = substr($1, 2) } NR % 4 == 2 { seq = $0 } \
    NR % 4 == 0 { n++; if (n % 2) print name, 99, \"chr1\", 101, 60, \"16M\", \"=\", 201, 116, seq, $0; \
    else print name, 147, \"chr1\", 201, 60, \"16M\", \"=\", 101, -116, seq, $0 }' \
    # {ref}";

/// One raw read pair of a duplex molecule: the A strand is R1 forward at 101 and
/// R2 reverse at 201; the B strand the reverse (R1 reverse at 201, R2 forward at
/// 101), with the swapped paired UMI. The bases are their own reverse complement.
fn duplex_raw_pair(name: &str, b_strand: bool) -> [fgumi_raw_bam::RawRecord; 2] {
    let (r1_pos, r2_pos) = if b_strand { (201, 101) } else { (101, 201) };
    let umi: &[u8] = if b_strand { b"CCCC-AAAA" } else { b"AAAA-CCCC" };
    let read = |seg: u16, pos: i32, mate_pos: i32| {
        let reverse = pos == 201;
        let mut b = SamBuilder::new();
        b.read_name(name.as_bytes())
            .flags(
                flags::PAIRED
                    | flags::PROPER_PAIR
                    | seg
                    | if reverse { flags::REVERSE } else { flags::MATE_REVERSE },
            )
            .ref_id(0)
            .pos(pos - 1)
            .mapq(60)
            .cigar_ops(&[16 << 4])
            .mate_ref_id(0)
            .mate_pos(mate_pos - 1)
            .template_length(if reverse { -116 } else { 116 })
            .sequence(b"AAAACCCCGGGGTTTT")
            .qualities(&[30; 16]);
        b.add_string_tag(SamTag::RX, umi).add_string_tag(SamTag::MC, b"16M");
        b.build()
    };
    [read(flags::FIRST_SEGMENT, r1_pos, r2_pos), read(flags::LAST_SEGMENT, r2_pos, r1_pos)]
}

/// Duplex end to end with no real aligner: `group --strategy paired
/// --realign-window` gives the A-R1/B-R2 reads one window and the A-R2/B-R1
/// reads another, `duplex` carries each onto its consensus end, and realign
/// matches each end's placement against its window by strand (consensus R1
/// forward at 101, R2 reverse at 201): both are kept in their windows.
#[test]
fn duplex_chain_realigns_into_shared_windows() {
    let tmp = TempDir::new().expect("tempdir");
    let dir = tmp.path();
    let records: Vec<_> = (0..3)
        .flat_map(|i| duplex_raw_pair(&format!("a{i}"), false))
        .chain((0..3).flat_map(|i| duplex_raw_pair(&format!("b{i}"), true)))
        .collect();
    let raw = dir.join("raw.bam");
    write_bam(&raw, &create_minimal_header("chr1", 10_000), &records);
    let reference = create_test_reference(dir);
    let (sorted, grouped, duplex, realigned) = (
        dir.join("sorted.bam"),
        dir.join("grouped.bam"),
        dir.join("duplex.bam"),
        dir.join("out.bam"),
    );
    let run = |args: &[&str]| {
        let result = fgumi().args(args).env("RUST_LOG", "info").output().expect("spawn fgumi");
        let stderr = String::from_utf8_lossy(&result.stderr).into_owned();
        assert!(result.status.success(), "fgumi {args:?} failed:\n{stderr}");
        stderr
    };
    run(&["sort", "-i", p(&raw), "-o", p(&sorted), "--order", "template-coordinate"]);
    run(&[
        "group",
        "-i",
        p(&sorted),
        "-o",
        p(&grouped),
        "--strategy",
        "paired",
        "--realign-window",
    ]);
    run(&["duplex", "-i", p(&grouped), "-o", p(&duplex), "--min-reads", "1"]);
    let stderr = run(&[
        "realign",
        "-i",
        p(&duplex),
        "-o",
        p(&realigned),
        "-r",
        p(&reference),
        "--aligner::command",
        DUPLEX_PAIR_ALIGNER,
    ]);
    assert_eq!(realign_count(&stderr, "Consensus reads realigned:"), 2);
    assert_eq!(realign_count(&stderr, "Primary already in window:"), 2);
    assert_eq!(realign_count(&stderr, "Without a window:"), 0);

    // Each consensus end carries the window of its strand's reads: the A-R1/B-R2
    // reads (forward at 101) for R1, the A-R2/B-R1 reads (reverse at 201) for R2.
    let mut reader = noodles::bam::io::Reader::new(std::io::BufReader::new(
        std::fs::File::open(&duplex).expect("duplex output"),
    ));
    reader.read_header().expect("header");
    let mut windows: Vec<(bool, Option<String>)> = reader
        .records()
        .map(|r| {
            let r = r.expect("record");
            let rw = match r.data().get(&SamTag::RW.to_noodles_tag()) {
                Some(Ok(noodles::sam::alignment::record::data::field::Value::String(v))) => {
                    Some(v.to_string())
                }
                _ => None,
            };
            (r.flags().is_first_segment(), rw)
        })
        .collect();
    windows.sort();
    assert_eq!(
        windows,
        vec![
            (false, Some("chr1:201-216:-".to_string())),
            (true, Some("chr1:101-116:+".to_string())),
        ]
    );

    // Both ends are written in their windows, keeping their consensus tags.
    let mut reader = noodles::bam::io::Reader::new(std::io::BufReader::new(
        std::fs::File::open(&realigned).expect("realign output"),
    ));
    let header = reader.read_header().expect("header");
    let mut placed: Vec<(bool, usize, bool)> = reader
        .record_bufs(&header)
        .map(|r| {
            let r = r.expect("record");
            (
                r.flags().is_first_segment(),
                r.alignment_start().expect("mapped").get(),
                r.flags().is_reverse_complemented(),
            )
        })
        .collect();
    placed.sort_unstable();
    assert_eq!(placed, vec![(false, 201, true), (true, 101, false)]);
    let ends = primary_ends(&realigned);
    assert!(ends.iter().all(|(_, _, rw, tags)| !rw && *tags), "{ends:?}");
}

/// A stand-in aligner for CODEC fragments: places each read forward at
/// chr1:5001, outside its window, with an `XA` alternative forward at chr1:117,
/// both spanning the whole read and echoing SEQ/QUAL.
const CODEC_FRAGMENT_ALIGNER: &str = "awk 'BEGIN { OFS = \"\\t\"; print \"@HD\", \"VN:1.6\"; \
    print \"@SQ\", \"SN:chr1\", \"LN:10000\" } \
    NR % 4 == 1 { name = substr($1, 2) } NR % 4 == 2 { seq = $0 } \
    NR % 4 == 0 { n = length(seq); print name, 0, \"chr1\", 5001, 60, n \"M\", \"*\", 0, 0, seq, $0, \
    \"NM:i:0\", \"AS:i:\" n, \"XA:Z:chr1,+117,\" n \"M,0;\" }' \
    # {ref}";

/// One raw CODEC read pair: R1 forward at chr1:101 and R2 reverse at chr1:105,
/// both `16M`, overlapping on 105-116 and sharing the `RX` UMI.
fn codec_raw_pair(name: &str) -> [fgumi_raw_bam::RawRecord; 2] {
    const MOLECULE: &[u8] = b"ACGTACGTACGTACGTACGT";
    let read = |seg: u16, pos: i32, mate_pos: i32, reverse: bool| {
        let offset = usize::try_from(pos - 101).expect("offset");
        let mut b = SamBuilder::new();
        b.read_name(name.as_bytes())
            .flags(
                flags::PAIRED
                    | flags::PROPER_PAIR
                    | seg
                    | if reverse { flags::REVERSE } else { flags::MATE_REVERSE },
            )
            .ref_id(0)
            .pos(pos - 1)
            .mapq(60)
            .cigar_ops(&[16 << 4])
            .mate_ref_id(0)
            .mate_pos(mate_pos - 1)
            .template_length(if reverse { -20 } else { 20 })
            .sequence(&MOLECULE[offset..offset + 16])
            .qualities(&[30; 16]);
        b.add_string_tag(SamTag::RX, b"ACGT").add_string_tag(SamTag::MC, b"16M");
        b.build()
    };
    [read(flags::FIRST_SEGMENT, 101, 105, false), read(flags::LAST_SEGMENT, 105, 101, true)]
}

/// CODEC end to end with no real aligner: `group --realign-window` gives the R1s
/// the window chr1:101-116:+ and the R2s chr1:105-120:-; `codec` unions them onto
/// its single fragment in R1's orientation (chr1:101-120:+), and realign, given
/// a placement outside it, promotes the `XA` hit at chr1:117, which overlaps only
/// the R2 part of the union and only on R1's strand.
#[test]
fn codec_chain_realigns_into_the_window_union() {
    let tmp = TempDir::new().expect("tempdir");
    let dir = tmp.path();
    let records: Vec<_> = (0..3).flat_map(|i| codec_raw_pair(&format!("c{i}"))).collect();
    let raw = dir.join("raw.bam");
    write_bam(&raw, &create_minimal_header("chr1", 10_000), &records);
    let reference = create_test_reference(dir);
    let (sorted, grouped, codec, realigned) = (
        dir.join("sorted.bam"),
        dir.join("grouped.bam"),
        dir.join("codec.bam"),
        dir.join("out.bam"),
    );
    let run = |args: &[&str]| {
        let result = fgumi().args(args).env("RUST_LOG", "info").output().expect("spawn fgumi");
        let stderr = String::from_utf8_lossy(&result.stderr).into_owned();
        assert!(result.status.success(), "fgumi {args:?} failed:\n{stderr}");
        stderr
    };
    run(&["sort", "-i", p(&raw), "-o", p(&sorted), "--order", "template-coordinate"]);
    run(&[
        "group",
        "-i",
        p(&sorted),
        "-o",
        p(&grouped),
        "--strategy",
        "identity",
        "--edits",
        "0",
        "--realign-window",
    ]);
    run(&[
        "codec",
        "-i",
        p(&grouped),
        "-o",
        p(&codec),
        "--min-reads",
        "1",
        "--min-duplex-length",
        "1",
    ]);
    let stderr = run(&[
        "realign",
        "-i",
        p(&codec),
        "-o",
        p(&realigned),
        "-r",
        p(&reference),
        "--aligner::command",
        CODEC_FRAGMENT_ALIGNER,
    ]);
    assert_eq!(realign_count(&stderr, "Consensus reads realigned:"), 1);
    assert_eq!(realign_count(&stderr, "Primary swapped to an in-window XA hit:"), 1);
    assert_eq!(realign_count(&stderr, "Without a window:"), 0);

    // The one CODEC fragment carries the union of the R1 and R2 windows, on R1's
    // strand.
    let mut reader = noodles::bam::io::Reader::new(std::io::BufReader::new(
        std::fs::File::open(&codec).expect("codec output"),
    ));
    reader.read_header().expect("header");
    let fragments: Vec<(bool, Option<String>)> = reader
        .records()
        .map(|r| {
            let r = r.expect("record");
            let rw = match r.data().get(&SamTag::RW.to_noodles_tag()) {
                Some(Ok(noodles::sam::alignment::record::data::field::Value::String(v))) => {
                    Some(v.to_string())
                }
                _ => None,
            };
            (r.flags().is_segmented(), rw)
        })
        .collect();
    assert_eq!(fragments, vec![(false, Some("chr1:101-120:+".to_string()))]);

    // The fragment is promoted to the in-window hit: forward at 117, MAPQ 0,
    // still a fragment, keeping its consensus tags and losing rw.
    let mut reader = noodles::bam::io::Reader::new(std::io::BufReader::new(
        std::fs::File::open(&realigned).expect("realign output"),
    ));
    let header = reader.read_header().expect("header");
    let placed: Vec<(bool, usize, bool, u8)> = reader
        .record_bufs(&header)
        .map(|r| {
            let r = r.expect("record");
            (
                r.flags().is_segmented(),
                r.alignment_start().expect("mapped").get(),
                r.flags().is_reverse_complemented(),
                r.mapping_quality().map_or(255, |q| q.get()),
            )
        })
        .collect();
    assert_eq!(placed, vec![(false, 117, false, 0)]);
    let ends = primary_ends(&realigned);
    assert!(ends.iter().all(|(_, _, rw, tags)| !rw && *tags), "{ends:?}");
}

/// The awk program of [`UNMAPPED_ALIGNER`], for a stand-in aligner script.
const UNMAPPED_AWK: &str = r#"BEGIN { OFS = "\t"; print "@HD", "VN:1.6"; print "@SQ", "SN:chr1", "LN:10000" } NR % 4 == 1 { name = substr($1, 2) } NR % 4 == 2 { seq = $0 } NR % 4 == 0 { n++; print name, (n % 2 ? 77 : 141), "*", 0, 0, "*", "*", 0, 0, seq, $0 }"#;

/// Writes a stand-in aligner script that records its arguments to the file
/// named by `FGUMI_TEST_ARGV`, then behaves like [`UNMAPPED_ALIGNER`]; returns an
/// `--aligner::command` running it with a `mem` word (`sh <script> mem {ref}`).
fn argv_recording_aligner(dir: &Path) -> String {
    let script = dir.join("aligner.sh");
    std::fs::write(
        &script,
        format!("printf '%s\\n' \"$*\" > \"$FGUMI_TEST_ARGV\"\nexec awk '{UNMAPPED_AWK}'\n"),
    )
    .expect("write aligner script");
    format!("sh {} mem {{ref}}", p(&script))
}

/// bwa's `-h` reaches the realign stage's aligner command only when the user
/// sets it, for `fgumi realign --xa-max-hits` and `runall
/// --realign::xa-max-hits`; unset, the command is run as written.
#[rstest::rstest]
#[case::realign_unset(false, &[], "mem {ref}")]
#[case::realign_set(false, &["--xa-max-hits", "7,9"], "mem -h 7,9 {ref}")]
#[case::runall_unset(true, &[], "mem {ref}")]
#[case::runall_set(true, &["--realign::xa-max-hits", "7,9"], "mem -h 7,9 {ref}")]
fn xa_max_hits_reaches_the_aligner(
    #[case] runall: bool,
    #[case] extra: &[&str],
    #[case] expected_argv: &str,
) {
    let tmp = TempDir::new().expect("tempdir");
    let dir = tmp.path();
    let input = write_consensus(dir, 3, |_| (Some("chr1:1-100:+"), Some("chr1:1-100:-")));
    let reference = create_test_reference(dir);
    let output = dir.join("out.bam");
    let argv = dir.join("argv.txt");
    let aligner = argv_recording_aligner(dir);
    let mut cmd = fgumi();
    if runall {
        cmd.args(["runall", "--start-from", "filter", "--stop-after", "filter", "--realign"])
            .args(["--filter::min-reads", "1", "--ref", p(&reference)]);
    } else {
        cmd.args(["realign", "-r", p(&reference)]);
    }
    let result = cmd
        .args(["-i", p(&input), "-o", p(&output), "--aligner::command", &aligner])
        .args(extra)
        .env("FGUMI_TEST_ARGV", &argv)
        .output()
        .expect("spawn fgumi");
    assert!(result.status.success(), "{}", String::from_utf8_lossy(&result.stderr));
    let recorded = std::fs::read_to_string(&argv).expect("the aligner ran");
    assert_eq!(recorded.trim_end(), expected_argv.replace("{ref}", p(&reference)));
}

/// A `--xa-max-hits` (or runall's `--realign::xa-max-hits`) left without a value
/// before an option that takes one is reported against the flag, not as clap's
/// "unexpected argument" for the option's value; a negative value given as its
/// own word still reaches the non-negative-integer check.
#[rstest::rstest]
#[case::realign_missing(false, "--xa-max-hits", &["--threads", "4"], "a value is required for '--xa-max-hits")]
#[case::realign_negative(false, "--xa-max-hits", &["-1,5"], "'-1' is not a non-negative integer")]
#[case::runall_missing(
    true,
    "--realign::xa-max-hits",
    &["--threads", "4"],
    "a value is required for '--realign::xa-max-hits"
)]
#[case::runall_negative(
    true,
    "--realign::xa-max-hits",
    &["-1,5"],
    "'-1' is not a non-negative integer"
)]
fn xa_max_hits_bad_value_names_the_flag(
    #[case] runall: bool,
    #[case] flag: &str,
    #[case] rest: &[&str],
    #[case] message: &str,
) {
    let mut cmd = fgumi();
    if runall {
        cmd.args(["runall", "--start-from", "filter", "--stop-after", "filter", "--realign"]);
    } else {
        cmd.args(["realign", "-r", "ref.fa"]);
    }
    let result =
        cmd.args(["-i", "in.bam", "-o", "out.bam", flag]).args(rest).output().expect("spawn");
    let stderr = String::from_utf8_lossy(&result.stderr);
    assert_eq!(result.status.code(), Some(2), "{stderr}");
    assert!(stderr.contains(flag), "{stderr}");
    assert!(stderr.contains(message), "{stderr}");
    assert!(!stderr.contains("unexpected argument"), "{stderr}");
}

/// A stand-in aligner that maps every pair FR on `chr1`: R1 forward at 101 and
/// R2 reverse at 121, both `16M` (as [`MAPPED_PAIR_RW`]'s windows expect). It
/// echoes SEQ/QUAL as given, so the reads must be their own reverse complement
/// with uniform qualities.
const MAPPED_PAIR_ALIGNER: &str = "awk 'BEGIN { OFS = \"\\t\"; print \"@HD\", \"VN:1.6\"; \
    print \"@SQ\", \"SN:chr1\", \"LN:10000\" } \
    NR % 4 == 1 { name = substr($1, 2) } NR % 4 == 2 { seq = $0 } \
    NR % 4 == 0 { n++; if (n % 2) print name, 99, \"chr1\", 101, 60, \"16M\", \"=\", 121, 36, seq, $0; \
    else print name, 147, \"chr1\", 121, 60, \"16M\", \"=\", 101, -36, seq, $0 }' \
    # {ref}";

/// The windows [`MAPPED_PAIR_ALIGNER`]'s placements fall in.
const MAPPED_PAIR_RW: (&str, &str) = ("chr1:101-116:+", "chr1:121-136:-");

/// realign merges the aligner's records onto the consensus reads as `zipper
/// --tags-to-reverse Consensus --tags-to-revcomp Consensus` would: a
/// reverse-strand record's per-base consensus arrays are reversed and its
/// per-base bases reverse-complemented; a forward-strand record's are kept.
#[test]
fn realign_reorients_per_base_tags_on_reverse_strand() {
    use noodles::sam::alignment::record_buf::data::field::Value;
    use noodles::sam::alignment::record_buf::data::field::value::Array;

    let tmp = TempDir::new().expect("tempdir");
    let dir = tmp.path();
    // Its own reverse complement, so the stand-in aligner can echo it.
    let seq = b"AAAACCCCGGGGTTTT";
    let depths: Vec<i16> = (0..16).collect();
    let bases = "AAAAAAAACCCCCCCC";
    let records: Vec<_> =
        [(flags::FIRST_SEGMENT, MAPPED_PAIR_RW.0), (flags::LAST_SEGMENT, MAPPED_PAIR_RW.1)]
            .into_iter()
            .map(|(seg, rw)| {
                let mut b = SamBuilder::new();
                b.read_name(b"cons:0")
                    .flags(flags::PAIRED | flags::UNMAPPED | flags::MATE_UNMAPPED | seg)
                    .sequence(seq)
                    .qualities(&[30; 16]);
                b.add_string_tag(SamTag::RW, rw.as_bytes())
                    .add_string_tag(SamTag::MI, b"0")
                    .add_int_tag(SamTag::CD, 5)
                    .add_float_tag(SamTag::CE, 0.0)
                    .add_array_i16(SamTag::CD_BASES, &depths)
                    .add_string_tag(SamTag::AC, bases.as_bytes());
                b.build()
            })
            .collect();
    let input = dir.join("consensus.bam");
    write_bam(&input, &create_minimal_header("chr1", 10_000), &records);
    let reference = create_test_reference(dir);
    let output = dir.join("out.bam");
    let result = fgumi()
        .args(["realign", "-i", p(&input), "-o", p(&output), "-r", p(&reference)])
        .args(["--aligner::command", MAPPED_PAIR_ALIGNER])
        .env("RUST_LOG", "info")
        .output()
        .expect("spawn fgumi realign");
    let stderr = String::from_utf8_lossy(&result.stderr);
    assert!(result.status.success(), "{stderr}");
    assert!(stderr.contains("Primary already in window: 2"), "{stderr}");

    let mut reader = noodles::bam::io::Reader::new(std::io::BufReader::new(
        std::fs::File::open(&output).expect("output"),
    ));
    let header = reader.read_header().expect("header");
    let cd = SamTag::CD_BASES.to_noodles_tag();
    let ac = SamTag::AC.to_noodles_tag();
    let mut seen = Vec::new();
    for record in reader.record_bufs(&header) {
        let record = record.expect("record");
        let reverse = record.flags().is_reverse_complemented();
        let Some(Value::Array(Array::Int16(got_depths))) = record.data().get(&cd) else {
            panic!("cd missing");
        };
        let Some(Value::String(got_bases)) = record.data().get(&ac) else {
            panic!("ac missing");
        };
        let (want_depths, want_bases): (Vec<i16>, &str) = if reverse {
            ((0..16).rev().collect(), "GGGGGGGGTTTTTTTT")
        } else {
            (depths.clone(), bases)
        };
        assert_eq!(got_depths, &want_depths, "cd (reverse = {reverse})");
        assert_eq!(got_bases.to_string(), want_bases, "ac (reverse = {reverse})");
        seen.push(reverse);
    }
    seen.sort_unstable();
    assert_eq!(seen, [false, true], "one forward and one reverse record");
}

/// The name and bases (upper-cased, so soft-masked bases read as any other)
/// of the first contig of the FASTA at `path`.
fn first_contig(path: &Path) -> (String, Vec<u8>) {
    let mut reader =
        noodles::fasta::io::reader::Builder.build_from_path(path).expect("open reference");
    let record = reader.records().next().expect("a contig").expect("read the first contig");
    let name = String::from_utf8_lossy(record.name()).into_owned();
    let bases = record.sequence().as_ref().iter().map(u8::to_ascii_uppercase).collect();
    (name, bases)
}

/// A chimeric consensus read through the real aligner: R1 is 45 bases from one
/// locus of the test reference followed by 105 from another 10 kb away, so bwa
/// makes the 105-base piece primary and the 45-base piece supplementary. R1's
/// window holds only the 45-base piece and R2 maps beside it, as when the raw
/// reads put the molecule there. With no `XA` hit in the window, realign swaps
/// the in-window piece to primary: the primary is that piece with the full read
/// (hard clips turned into soft clips) and its own MAPQ, the old primary becomes
/// a supplementary clipped as bwa clipped the piece, `SA` names each from the
/// other, and the consensus tags stay on both, the per-base ones full-length.
#[rstest::rstest]
#[case::subprocess("bwa-mem3")]
#[cfg_attr(feature = "aligner-bwa-mem3", case::inproc("bwa-mem3-inproc"))]
fn realign_swaps_an_in_window_supplementary_piece_to_primary(#[case] preset: &str) {
    use noodles::sam::alignment::record_buf::data::field::Value;
    use noodles::sam::alignment::record_buf::data::field::value::Array;

    let test = format!("realign_swaps_an_in_window_supplementary_piece_to_primary ({preset})");
    if !bwa_mem3_or_skip(&test) {
        return;
    }
    let Some(reference) = bwa_mem3_test_ref(&test) else {
        return;
    };
    let (contig, bases) = first_contig(&reference);
    assert!(bases.len() >= 15_000, "the test reference is too short");
    let (a, b) = (2_000usize, 12_000usize);
    // R1 in read orientation: the in-window piece, then the out-of-window one.
    let r1: Vec<u8> = [&bases[a..a + 45], &bases[b..b + 105]].concat();
    // R2: reverse strand, 150 bases starting 200 past R1's in-window piece.
    let r2 = fgumi_dna::dna::reverse_complement(&bases[a + 200..a + 350]);
    let quals: Vec<u8> = (0..150u8).map(|i| 20 + i % 20).collect();
    let depths: Vec<i16> = (0..150).collect();
    let rw1 = format!("{contig}:{}-{}:+", a + 1, a + 45);
    let rw2 = format!("{contig}:{}-{}:-", a + 201, a + 350);
    let records: Vec<_> = [(flags::FIRST_SEGMENT, &r1, &rw1), (flags::LAST_SEGMENT, &r2, &rw2)]
        .into_iter()
        .map(|(seg, seq, rw)| {
            let mut builder = SamBuilder::new();
            builder
                .read_name(b"chimera:1")
                .flags(flags::PAIRED | flags::UNMAPPED | flags::MATE_UNMAPPED | seg)
                .sequence(seq)
                .qualities(&quals);
            builder
                .add_string_tag(SamTag::RW, rw.as_bytes())
                .add_string_tag(SamTag::MI, b"1")
                .add_string_tag(SamTag::RX, b"ACGT-TGCA")
                .add_int_tag(SamTag::CD, 9)
                .add_float_tag(SamTag::CE, 0.002)
                .add_array_i16(SamTag::CD_BASES, &depths);
            builder.build()
        })
        .collect();
    let tmp = TempDir::new().expect("tempdir");
    let input = tmp.path().join("consensus.bam");
    write_bam(&input, &create_minimal_header(&contig, bases.len()), &records);
    let output = tmp.path().join("realigned.bam");
    let stderr = crate::helpers::cli::run_args_and_capture_logs(
        &["realign", "-i", p(&input), "-o", p(&output), "-r", p(&reference)]
            .into_iter()
            .chain(["--aligner::preset", preset, "--threads", "2"])
            .map(std::ffi::OsStr::new)
            .collect::<Vec<_>>(),
    );
    assert_eq!(
        realign_count(&stderr, "Primary swapped to an in-window supplementary piece:"),
        1,
        "{stderr}"
    );
    assert_eq!(realign_count(&stderr, "Primary already in window:"), 1, "{stderr}");

    let (_, records) = crate::helpers::read_bam_output(&output);
    let r1_records: Vec<_> = records.iter().filter(|r| r.flags().is_first_segment()).collect();
    assert_eq!(r1_records.len(), 2, "R1's primary and its one piece");
    let primary = r1_records.iter().find(|r| !r.flags().is_supplementary()).expect("primary");
    let piece = r1_records.iter().find(|r| r.flags().is_supplementary()).expect("piece");

    let start = |r: &noodles::sam::alignment::RecordBuf| r.alignment_start().expect("pos").get();
    let cigar = |r: &noodles::sam::alignment::RecordBuf| {
        use noodles::sam::alignment::record::cigar::op::Kind;
        use std::fmt::Write as _;
        let mut cigar = String::new();
        for op in r.cigar().as_ref() {
            let c = match op.kind() {
                Kind::Match => 'M',
                Kind::SoftClip => 'S',
                Kind::HardClip => 'H',
                Kind::Insertion => 'I',
                Kind::Deletion => 'D',
                other => panic!("unexpected op {other:?}"),
            };
            write!(cigar, "{}{c}", op.len()).expect("write");
        }
        cigar
    };
    let string_tag = |r: &noodles::sam::alignment::RecordBuf, tag: SamTag| match r
        .data()
        .get(&tag.to_noodles_tag())
    {
        Some(Value::String(s)) => Some(s.to_string()),
        _ => None,
    };

    // The promoted piece: in the window, forward, the whole read soft-clipped.
    assert_eq!(start(primary), a + 1, "promoted piece position");
    assert!(!primary.flags().is_reverse_complemented());
    assert_eq!(cigar(primary), "45M105S");
    assert!(primary.mapping_quality().is_some_and(|q| q.get() > 0), "the piece's own MAPQ");
    assert_eq!(primary.sequence().as_ref(), r1.as_slice(), "full SEQ");
    assert_eq!(primary.quality_scores().as_ref(), quals.as_slice(), "full QUAL");
    assert!(string_tag(primary, SamTag::XA).is_none(), "no XA on the promoted piece");

    // The demoted primary: out of the window, clipped as bwa clipped the piece.
    assert_eq!(start(piece), b + 1, "demoted primary position");
    let hard = cigar(piece) == "45H105M";
    assert!(hard || cigar(piece) == "45S105M", "demoted CIGAR {}", cigar(piece));
    let want_seq = if hard { &r1[45..] } else { &r1[..] };
    assert_eq!(piece.sequence().as_ref(), want_seq, "demoted SEQ");

    // SA names each from the other, hard clips written as soft clips.
    let mapq =
        |r: &noodles::sam::alignment::RecordBuf| r.mapping_quality().map_or(255, |q| q.get());
    let nm =
        |r: &noodles::sam::alignment::RecordBuf| match r.data().get(&SamTag::NM.to_noodles_tag()) {
            Some(v) => v.as_int().expect("NM"),
            None => panic!("NM missing"),
        };
    assert_eq!(
        string_tag(primary, SamTag::SA),
        Some(format!("{contig},{},+,45S105M,{},{};", b + 1, mapq(piece), nm(piece)))
    );
    assert_eq!(
        string_tag(piece, SamTag::SA),
        Some(format!("{contig},{},+,45M105S,{},{};", a + 1, mapq(primary), nm(primary)))
    );

    // The consensus tags are on both records, the per-base ones full-length
    // (both records are forward, so in input order).
    for record in [primary, piece] {
        assert_eq!(string_tag(record, SamTag::MI).as_deref(), Some("1"), "MI");
        assert_eq!(string_tag(record, SamTag::RX).as_deref(), Some("ACGT-TGCA"), "RX");
        assert!(string_tag(record, SamTag::RW).is_none(), "rw stripped");
        let Some(Value::Array(Array::Int16(cd))) =
            record.data().get(&SamTag::CD_BASES.to_noodles_tag())
        else {
            panic!("cd missing");
        };
        assert_eq!(cd, &depths, "cd");
        assert!(record.data().get(&SamTag::CD.to_noodles_tag()).is_some(), "cD");
        assert!(record.data().get(&SamTag::CE.to_noodles_tag()).is_some(), "cE");
    }

    // R2 kept, and the mates point at each other's primaries.
    let r2_primary = records
        .iter()
        .find(|r| r.flags().is_last_segment() && !r.flags().is_supplementary())
        .expect("R2");
    assert_eq!(start(r2_primary), a + 201, "R2 kept");
    assert_eq!(r2_primary.mate_alignment_start().map(|p| p.get()), Some(a + 1), "R2's mate");
    for record in [primary, piece] {
        assert_eq!(record.mate_alignment_start().map(|p| p.get()), Some(a + 201), "R1's mate");
    }
}

/// The number following `label` on the first logged line that has it, e.g.
/// `Unmapped: 3`.
fn logged(stderr: &str, label: &str) -> u64 {
    stderr
        .lines()
        .find_map(|line| line.split_once(label).map(|(_, n)| n.trim()))
        .unwrap_or_else(|| panic!("no `{label}` line in:\n{stderr}"))
        .parse()
        .expect("count")
}

/// The count realign's summary logs after `label`, read from the `=== Realign
/// Summary ===` block on, so a line with the same label logged by another stage
/// is never read. `Without a window:` is logged only when non-zero, so its
/// absence reads as 0.
fn realign_count(stderr: &str, label: &str) -> u64 {
    let at = stderr
        .find("=== Realign Summary ===")
        .unwrap_or_else(|| panic!("no realign summary in:\n{stderr}"));
    let summary = &stderr[at..];
    if label == "Without a window:" && !summary.contains(label) {
        return 0;
    }
    logged(summary, label)
}

/// `(name, first segment, carries rw, carries cD and cE)` of each primary
/// record in the BAM at `path`.
fn primary_ends(path: &Path) -> Vec<(String, bool, bool, bool)> {
    let mut reader = noodles::bam::io::Reader::new(std::io::BufReader::new(
        std::fs::File::open(path).expect("open output"),
    ));
    reader.read_header().expect("header");
    let has = |record: &noodles::bam::Record, tag: SamTag| {
        record.data().get(&tag.to_noodles_tag()).is_some()
    };
    let mut ends: Vec<_> = reader
        .records()
        .map(|r| r.expect("record"))
        .filter(|r| !r.flags().is_secondary() && !r.flags().is_supplementary())
        .map(|r| {
            (
                r.name().map(ToString::to_string).unwrap_or_default(),
                r.flags().is_first_segment(),
                has(&r, SamTag::RW),
                has(&r, SamTag::CD) && has(&r, SamTag::CE),
            )
        })
        .collect();
    ends.sort();
    ends
}

/// The simulated standalone chain (needs `fgumi simulate`).
#[cfg(feature = "simulate")]
mod chain {
    use std::collections::HashMap;
    use std::fmt::Write as _;
    use std::io::Write;
    use std::path::{Path, PathBuf};

    use fgumi_lib::sam::SamTag;
    use noodles::sam::alignment::record::Cigar as _;
    use noodles::sam::alignment::record::data::field::Value;
    use rstest::rstest;
    use tempfile::TempDir;

    use super::{bwa_mem3_or_skip, fgumi, p, realign_count};
    use crate::helpers::aligner::build_aligner_index;
    use crate::helpers::cli::run_args_and_capture_logs;

    /// Deterministic pseudo-random ACGT sequence (xorshift64).
    fn random_sequence(n: usize, seed: u64) -> String {
        let mut state = seed | 1;
        (0..n)
            .map(|_| {
                state ^= state << 13;
                state ^= state >> 7;
                state ^= state << 17;
                ['A', 'C', 'G', 'T'][usize::try_from(state % 4).expect("fits")]
            })
            .collect()
    }

    /// Writes a two-contig reference (`.fa`, `.fai`, `.dict`). `chr2` repeats a
    /// 2 kb segment of `chr1`, so reads from that segment multi-map and bwa may
    /// place a consensus read on the copy the raw reads did not come from.
    fn write_reference(dir: &Path) -> PathBuf {
        let chr1 = random_sequence(20_000, 0x5eed_0001);
        let chr2 = format!("{}{}", &chr1[5_000..7_000], random_sequence(8_000, 0x5eed_0002));
        let fasta = dir.join("ref.fa");
        let mut fa = std::fs::File::create(&fasta).expect("create fasta");
        let mut fai = String::new();
        let mut dict = String::from("@HD\tVN:1.6\tSO:unsorted\n");
        let mut offset = 0usize;
        for (name, seq) in [("chr1", &chr1), ("chr2", &chr2)] {
            let header = format!(">{name}\n");
            writeln!(fa, "{header}{seq}").expect("write fasta");
            offset += header.len();
            let len = seq.len();
            writeln!(fai, "{name}\t{len}\t{offset}\t{len}\t{}", len + 1).expect("fai");
            writeln!(dict, "@SQ\tSN:{name}\tLN:{len}").expect("dict");
            offset += len + 1;
        }
        std::fs::write(dir.join("ref.fa.fai"), fai).expect("write fai");
        std::fs::write(dir.join("ref.dict"), dict).expect("write dict");
        fasta
    }

    fn run_ok(args: &[&str]) {
        let out = fgumi().args(args).output().expect("spawn fgumi");
        assert!(
            out.status.success(),
            "fgumi {} failed:\n{}",
            args.join(" "),
            String::from_utf8_lossy(&out.stderr)
        );
    }

    /// Every this-many-th simulated molecule has its R2 unmapped (see
    /// [`simulate_with_unmapped_mates`]).
    const UNMAPPED_MATE_EVERY: usize = 25;

    /// `simulate mapped-reads` (200 molecules, 100 bp, `seed`), then every read
    /// pair of each [`UNMAPPED_MATE_EVERY`]-th molecule gets its R2 unmapped
    /// (placed at its mate, as aligners do) and the result is re-sorted
    /// template-coordinate. `group --realign-window` stamps only the mapped
    /// R1s; simplex then calls no R2 consensus from the unmapped R2s and drops
    /// those molecules as orphans, so the chain must carry mate-unmapped input
    /// through group, consensus and realign without tripping the `rw` check.
    fn simulate_with_unmapped_mates(dir: &Path, reference: &Path, seed: u64) -> PathBuf {
        use noodles::sam::alignment::io::Write as _;
        use noodles::sam::alignment::record::Flags;
        use noodles::sam::alignment::record_buf::Cigar;

        let simulated = dir.join("simulated.bam");
        run_ok(&[
            "simulate",
            "mapped-reads",
            "-r",
            p(reference),
            "-o",
            p(&simulated),
            "--truth",
            p(&dir.join("truth.tsv")),
            "--num-molecules",
            "200",
            "--read-length",
            "100",
            "--seed",
            &seed.to_string(),
        ]);
        let edited = dir.join("edited.bam");
        let mut reader =
            noodles::bam::io::reader::Builder.build_from_path(&simulated).expect("open");
        let header = reader.read_header().expect("header");
        let mut writer =
            noodles::bam::io::writer::Builder.build_from_path(&edited).expect("create");
        writer.write_header(&header).expect("write header");
        let mut unmapped = 0usize;
        for result in reader.record_bufs(&header) {
            let mut record = result.expect("record");
            let name = record.name().map(ToString::to_string).unwrap_or_default();
            let molecule: usize = name
                .strip_prefix("mol")
                .and_then(|n| n.split('_').next())
                .and_then(|n| n.parse().ok())
                .unwrap_or_else(|| panic!("unexpected simulated read name {name}"));
            if molecule.is_multiple_of(UNMAPPED_MATE_EVERY) {
                let mut f = record.flags();
                f.remove(Flags::PROPERLY_SEGMENTED);
                if f.is_last_segment() {
                    f.insert(Flags::UNMAPPED);
                    *record.reference_sequence_id_mut() = record.mate_reference_sequence_id();
                    *record.alignment_start_mut() = record.mate_alignment_start();
                    *record.mapping_quality_mut() = None;
                    *record.cigar_mut() = Cigar::default();
                    unmapped += 1;
                } else {
                    // R1's mate is now unmapped; R2 keeps `MC`/`MQ`, which still
                    // describe R1 (template-coordinate keys read its `MC`).
                    f.insert(Flags::MATE_UNMAPPED);
                    *record.mate_reference_sequence_id_mut() = record.reference_sequence_id();
                    *record.mate_alignment_start_mut() = record.alignment_start();
                    record.data_mut().remove(&SamTag::MC.to_noodles_tag());
                    record.data_mut().remove(&SamTag::MQ.to_noodles_tag());
                }
                *record.flags_mut() = f;
                *record.template_length_mut() = 0;
            }
            writer.write_alignment_record(&header, &record).expect("write record");
        }
        writer.try_finish().expect("finish");
        drop(writer);
        assert!(unmapped > 0, "no simulated molecule had its R2 unmapped");
        let mapped = dir.join("mapped.bam");
        run_ok(&["sort", "-i", p(&edited), "-o", p(&mapped), "--order", "template-coordinate"]);
        mapped
    }

    /// Copies the consensus BAM `src` to `dst`, removing `rw` from the R2 of
    /// every consensus pair whose `MI` is a multiple of `every`: stands in for
    /// the consensus ends that carry no window in a window-grouped input (a
    /// duplex strand or CODEC end without mapped raw reads, `--allow-unmapped`
    /// families). Returns how many ends lost `rw`.
    fn strip_some_r2_windows(src: &Path, dst: &Path, every: usize) -> u64 {
        use noodles::sam::alignment::io::Write as _;

        let mut reader = noodles::bam::io::reader::Builder.build_from_path(src).expect("open");
        let header = reader.read_header().expect("header");
        let mut writer = noodles::bam::io::writer::Builder.build_from_path(dst).expect("create");
        writer.write_header(&header).expect("write header");
        let mut stripped = 0;
        for result in reader.record_bufs(&header) {
            let mut record = result.expect("record");
            let mi: usize = match record.data().get(&SamTag::MI.to_noodles_tag()) {
                Some(noodles::sam::alignment::record_buf::data::field::Value::String(v)) => {
                    v.to_string().parse().expect("numeric MI")
                }
                other => panic!("consensus read without a string MI: {other:?}"),
            };
            if record.flags().is_last_segment() && mi.is_multiple_of(every) {
                record.data_mut().remove(&SamTag::RW.to_noodles_tag()).expect("consensus rw");
                stripped += 1;
            }
            writer.write_alignment_record(&header, &record).expect("write record");
        }
        writer.try_finish().expect("finish");
        stripped
    }

    /// Segment key of a record: the read name and which end it is.
    fn segment_key(record: &noodles::bam::Record) -> (String, u8) {
        let f = record.flags();
        let seg = if f.is_last_segment() { 2 } else { u8::from(f.is_first_segment()) };
        (record.name().map(std::string::ToString::to_string).unwrap_or_default(), seg)
    }

    fn string_tag(record: &noodles::bam::Record, tag: SamTag) -> Option<String> {
        match record.data().get(&tag.to_noodles_tag()) {
            Some(Ok(Value::String(v))) => Some(v.to_string()),
            _ => None,
        }
    }

    /// A parsed `rw` value: `contig:start-end:strand` (1-based, inclusive).
    struct Window {
        contig: String,
        start: usize,
        end: usize,
        reverse: bool,
    }

    fn parse_window(rw: &str) -> Window {
        let (rest, strand) = rw.rsplit_once(':').expect("strand");
        let (contig, range) = rest.rsplit_once(':').expect("range");
        let (start, end) = range.split_once('-').expect("start-end");
        Window {
            contig: contig.to_string(),
            start: start.parse().expect("start"),
            end: end.parse().expect("end"),
            reverse: strand == "-",
        }
    }

    /// Upper bounds on the two proper-pair guard counts, from the output.
    /// `ends` maps each mapped consensus primary with a window, keyed by its
    /// template and segment (1 or 2), to whether it is flagged proper and
    /// whether it lies in its window. A guarded end keeps bwa's placement, flag
    /// included, outside its window, so it is one of these ends: kept "with an
    /// in-window mate" only if its mate is in its window, and "whose mate cannot
    /// move" only if it is not. Returns the two counts in that order.
    fn guard_bounds<K: Eq + std::hash::Hash + Clone>(
        ends: &HashMap<(K, u8), (bool, bool)>,
    ) -> (u64, u64) {
        let (mut mate_in, mut mate_out) = (0, 0);
        for ((template, seg), &(proper, inside)) in ends {
            if !proper || inside {
                continue;
            }
            if ends.get(&(template.clone(), 3 - seg)).is_some_and(|&(_, mate_inside)| mate_inside) {
                mate_in += 1;
            } else {
                mate_out += 1;
            }
        }
        (mate_in, mate_out)
    }

    #[derive(Debug, Clone, Copy)]
    enum Backend {
        Subprocess,
        #[cfg_attr(not(feature = "aligner-bwa-mem3"), allow(dead_code))]
        InProcess,
    }

    /// Copies the consensus BAM `src` to `dst`, setting `cD` to 0 on R1 of every
    /// pair whose `MI` is 1 mod 6 and on R2 of every pair whose `MI` is 4 mod 6, so
    /// a per-record `filter --min-reads 1` keeps the other end alone. The lone
    /// mates are spaced out by whole pairs, so the read count between two of
    /// them is odd (a lone mate next to another would restore even parity).
    fn zero_depth_on_one_end(src: &Path, dst: &Path) {
        use noodles::sam::alignment::io::Write as _;
        use noodles::sam::alignment::record_buf::data::field::Value as BufValue;

        let mut reader = noodles::bam::io::reader::Builder.build_from_path(src).expect("open");
        let header = reader.read_header().expect("header");
        let mut writer = noodles::bam::io::writer::Builder.build_from_path(dst).expect("create");
        writer.write_header(&header).expect("write header");
        for result in reader.record_bufs(&header) {
            let mut record = result.expect("record");
            let mi: usize = match record.data().get(&SamTag::MI.to_noodles_tag()) {
                Some(BufValue::String(v)) => v.to_string().parse().expect("numeric MI"),
                other => panic!("consensus read without a string MI: {other:?}"),
            };
            let f = record.flags();
            if (mi % 6 == 1 && f.is_first_segment()) || (mi % 6 == 4 && f.is_last_segment()) {
                record.data_mut().insert(SamTag::CD.to_noodles_tag(), BufValue::Int32(0));
            }
            writer.write_alignment_record(&header, &record).expect("write record");
        }
        writer.try_finish().expect("finish");
    }

    /// A record's tags as `(tag, value)` strings.
    type Tags = Vec<(String, String)>;

    /// A record's `(reference id, 1-based position)`, or its mate's.
    type Placement = (Option<usize>, Option<usize>);

    /// The tags of `record` as `(tag, value)` strings, without `rw`.
    fn tag_values(record: &noodles::bam::Record) -> Tags {
        let mut tags: Vec<(String, String)> = record
            .data()
            .iter()
            .map(|field| {
                let (tag, value) = field.expect("tag");
                (String::from_utf8_lossy(tag.as_ref()).into_owned(), format!("{value:?}"))
            })
            .filter(|(tag, _)| tag != "rw")
            .collect();
        tags.sort();
        tags
    }

    /// A consensus pair whose other end a per-record filter removed reaches
    /// realign as a lone mate (`PAIRED` with its segment bit, mate absent). It
    /// keeps every consensus tag and its R1/R2 identity through realign, uses its
    /// window, and comes out as a pair end whose mate is unmapped: `RNEXT`/`PNEXT`
    /// at its own placement, `TLEN` 0, no `MC`/`MQ`, never proper. Lone R1s and
    /// lone R2s are both covered.
    ///
    /// With a small `--aligner::chunk-size` (15 reads, not a multiple of the
    /// 10-read period of the lone mates) the input crosses many of bwa's `-K`
    /// cuts, which fall only after an even number of reads: the lone mates must
    /// not shift a cut between the two ends of a real pair, so every real pair
    /// still comes out paired, and no read name has more than one primary per
    /// end.
    #[rstest]
    #[case::subprocess(Backend::Subprocess)]
    #[cfg_attr(feature = "aligner-bwa-mem3", case::inproc(Backend::InProcess))]
    fn realign_keeps_a_lone_mates_tags_and_segment(
        #[case] backend: Backend,
        #[values(None, Some("1500"))] chunk_size: Option<&str>,
    ) {
        if !bwa_mem3_or_skip(&format!(
            "realign_keeps_a_lone_mates_tags_and_segment ({backend:?}, {chunk_size:?})"
        )) {
            return;
        }
        let tmp = TempDir::new().expect("tempdir");
        let dir = tmp.path();
        let reference = write_reference(dir);
        build_aligner_index(&reference, "bwa-mem3");
        let mapped = simulate_with_unmapped_mates(dir, &reference, 11);
        let grouped = dir.join("grouped.bam");
        let simplex = dir.join("simplex.bam");
        let lowered = dir.join("lowered.bam");
        let filtered = dir.join("filtered.bam");
        let realigned = dir.join("realigned.bam");
        run_ok(&[
            "group",
            "-i",
            p(&mapped),
            "-o",
            p(&grouped),
            "--strategy",
            "adjacency",
            "--realign-window",
        ]);
        run_ok(&["simplex", "-i", p(&grouped), "-o", p(&simplex), "--min-reads", "1"]);
        zero_depth_on_one_end(&simplex, &lowered);
        run_ok(&[
            "filter",
            "-i",
            p(&lowered),
            "-o",
            p(&filtered),
            "-r",
            p(&reference),
            "--min-reads",
            "1",
            "--filter-by-template",
            "false",
        ]);

        // The lone ends of the filtered input, by name.
        let mut primaries: HashMap<String, Vec<(u8, Tags)>> = HashMap::new();
        let mut reader =
            noodles::bam::io::reader::Builder.build_from_path(&filtered).expect("open");
        reader.read_header().expect("header");
        for record in reader.records() {
            let record = record.expect("record");
            assert!(string_tag(&record, SamTag::RW).is_some(), "every filtered end has rw");
            let (name, seg) = segment_key(&record);
            primaries.entry(name).or_default().push((seg, tag_values(&record)));
        }
        let lone: HashMap<String, (u8, Tags)> = primaries
            .into_iter()
            .filter(|(_, ends)| ends.len() == 1)
            .map(|(name, mut ends)| (name, ends.pop().expect("one end")))
            .collect();
        for seg in [1, 2] {
            assert!(lone.values().any(|(s, _)| *s == seg), "no lone R{seg} in the input");
        }

        let mut args = vec![
            "realign",
            "-i",
            p(&filtered),
            "-o",
            p(&realigned),
            "-r",
            p(&reference),
            "--aligner::preset",
            match backend {
                Backend::Subprocess => "bwa-mem3",
                Backend::InProcess => "bwa-mem3-inproc",
            },
        ];
        if let Some(chunk_size) = chunk_size {
            args.extend(["--aligner::chunk-size", chunk_size]);
        }
        let stderr =
            run_args_and_capture_logs(&args.iter().map(std::ffi::OsStr::new).collect::<Vec<_>>());
        assert!(!stderr.contains("Without a window"), "a lone mate's rw must be used:\n{stderr}");

        // Every read name's primaries: one per end, and a real pair's two ends
        // name each other as mates.
        let mut ends: HashMap<String, Vec<(u8, bool, Placement, Placement)>> = HashMap::new();
        let mut reader =
            noodles::bam::io::reader::Builder.build_from_path(&realigned).expect("open");
        reader.read_header().expect("header");
        for record in reader.records() {
            let record = record.expect("record");
            let f = record.flags();
            if f.is_secondary() || f.is_supplementary() {
                continue;
            }
            let (name, seg) = segment_key(&record);
            let own = (
                record.reference_sequence_id().transpose().expect("tid"),
                record.alignment_start().transpose().expect("pos").map(usize::from),
            );
            let mate = (
                record.mate_reference_sequence_id().transpose().expect("mate tid"),
                record.mate_alignment_start().transpose().expect("mate pos").map(usize::from),
            );
            ends.entry(name).or_default().push((seg, f.is_segmented(), own, mate));
        }
        for (name, primaries) in &ends {
            let mut segs: Vec<u8> = primaries.iter().map(|e| e.0).collect();
            segs.sort_unstable();
            if lone.contains_key(name) {
                assert_eq!(segs.len(), 1, "{name}: one primary for a lone mate, got {primaries:?}");
                continue;
            }
            assert_eq!(segs, [1, 2], "{name}: a real pair must stay paired, got {primaries:?}");
            assert!(primaries.iter().all(|e| e.1), "{name}: lost PAIRED: {primaries:?}");
            let (a, b) = (&primaries[0], &primaries[1]);
            if a.2.1.is_some() && b.2.1.is_some() {
                assert_eq!((a.3, b.3), (b.2, a.2), "{name}: mates must name each other");
            }
        }

        let mut reader =
            noodles::bam::io::reader::Builder.build_from_path(&realigned).expect("open");
        reader.read_header().expect("header");
        let mut checked = 0;
        for record in reader.records() {
            let record = record.expect("record");
            let (name, seg) = segment_key(&record);
            let Some((want_seg, want_tags)) = lone.get(&name) else {
                continue;
            };
            let f = record.flags();
            assert!(f.is_segmented(), "{name}: lost PAIRED ({f:?})");
            assert_eq!(seg, *want_seg, "{name}: lost its R1/R2 identity ({f:?})");
            assert!(f.is_mate_unmapped(), "{name}: its absent mate must be unmapped ({f:?})");
            assert!(!f.is_properly_segmented(), "{name}: a lone mate is never proper");
            assert!(!f.is_mate_reverse_complemented(), "{name}: ({f:?})");
            assert_eq!(record.template_length(), 0, "{name}: TLEN");
            assert!(string_tag(&record, SamTag::MC).is_none(), "{name}: MC");
            assert!(record.data().get(&SamTag::MQ.to_noodles_tag()).is_none(), "{name}: MQ");
            assert!(string_tag(&record, SamTag::RW).is_none(), "{name}: rw must be stripped");
            if f.is_secondary() || f.is_supplementary() {
                continue;
            }
            let own = (
                record.reference_sequence_id().transpose().expect("tid"),
                record.alignment_start().transpose().expect("pos"),
            );
            let mate = (
                record.mate_reference_sequence_id().transpose().expect("mate tid"),
                record.mate_alignment_start().transpose().expect("mate pos"),
            );
            if f.is_unmapped() {
                assert_eq!(mate, (None, None), "{name}: an unplaced lone mate");
            } else {
                assert_eq!(mate, own, "{name}: RNEXT/PNEXT must be the read's own placement");
            }
            // Every consensus tag survives (per-base tags may be re-oriented on the
            // reverse strand, so only their presence is compared for those).
            let got = tag_values(&record);
            for (tag, value) in want_tags {
                let found = got.iter().find(|(t, _)| t == tag);
                assert!(found.is_some(), "{name}: lost {tag}");
                if ["MI", "RX", "RG", "cD", "cM", "cE"].contains(&tag.as_str()) {
                    assert_eq!(&found.expect("found").1, value, "{name}: {tag} changed");
                }
            }
            checked += 1;
        }
        assert_eq!(checked, lone.len(), "one primary per lone mate");
    }

    /// `simulate → group --realign-window → simplex → realign → sort`: every
    /// consensus primary is mapped inside its molecule's window or counted as
    /// outside/unmapped, keeps its consensus tags, and carries no `rw`.
    #[rstest]
    #[case::subprocess(Backend::Subprocess)]
    #[cfg_attr(feature = "aligner-bwa-mem3", case::inproc(Backend::InProcess))]
    fn standalone_chain_aligns_consensus(#[case] backend: Backend) {
        if !bwa_mem3_or_skip(&format!("standalone_chain_aligns_consensus ({backend:?})")) {
            return;
        }
        let tmp = TempDir::new().expect("tempdir");
        let dir = tmp.path();
        let reference = write_reference(dir);
        build_aligner_index(&reference, "bwa-mem3");

        let mapped = simulate_with_unmapped_mates(dir, &reference, 7);
        let grouped = dir.join("grouped.bam");
        let simplex = dir.join("simplex.bam");
        let consensus = dir.join("consensus.bam");
        let realigned = dir.join("realigned.bam");
        let sorted = dir.join("sorted.bam");
        run_ok(&[
            "group",
            "-i",
            p(&mapped),
            "-o",
            p(&grouped),
            "--strategy",
            "adjacency",
            "--edits",
            "1",
            "--realign-window",
        ]);
        run_ok(&["simplex", "-i", p(&grouped), "-o", p(&simplex), "--min-reads", "1"]);
        let stripped = strip_some_r2_windows(&simplex, &consensus, 10);
        assert!(stripped > 0, "no consensus R2 lost its window");

        let preset = match backend {
            Backend::Subprocess => "bwa-mem3",
            Backend::InProcess => "bwa-mem3-inproc",
        };
        let stderr = run_args_and_capture_logs(
            &[
                "realign",
                "-i",
                p(&consensus),
                "-o",
                p(&realigned),
                "-r",
                p(&reference),
                "--aligner::preset",
                preset,
                "--threads",
                "4",
            ]
            .map(std::ffi::OsStr::new),
        );
        run_ok(&["sort", "-i", p(&realigned), "-o", p(&sorted), "--order", "coordinate"]);

        // Each consensus end's window, from the unmapped consensus BAM.
        // `None`: a consensus end whose raw reads never mapped (no window).
        let mut windows: HashMap<(String, u8), Option<String>> = HashMap::new();
        let mut reader =
            noodles::bam::io::reader::Builder.build_from_path(&consensus).expect("open");
        reader.read_header().expect("header");
        for record in reader.records() {
            let record = record.expect("record");
            windows.insert(segment_key(&record), string_tag(&record, SamTag::RW));
        }
        assert!(!windows.is_empty(), "simplex produced no consensus reads");
        let windowless = windows.values().filter(|w| w.is_none()).count() as u64;
        assert_eq!(windowless, stripped, "only the stripped ends lack rw");

        let mut reader = noodles::bam::io::reader::Builder.build_from_path(&sorted).expect("open");
        let header = reader.read_header().expect("header");
        let (mut primaries, mut in_window, mut outside, mut unmapped) = (0u64, 0u64, 0u64, 0u64);
        let mut no_window = 0u64;
        let mut placed: HashMap<(String, u8), (bool, bool)> = HashMap::new();
        for record in reader.records() {
            let record = record.expect("record");
            assert!(string_tag(&record, SamTag::RW).is_none(), "rw must be stripped");
            let f = record.flags();
            if f.is_secondary() || f.is_supplementary() {
                continue;
            }
            primaries += 1;
            let key = segment_key(&record);
            let window = windows.get(&key).unwrap_or_else(|| panic!("{key:?}"));
            assert!(
                record.data().get(&SamTag::CD.to_noodles_tag()).is_some(),
                "{key:?} lost its cD tag"
            );
            assert!(
                record.data().get(&SamTag::CE.to_noodles_tag()).is_some(),
                "{key:?} lost its cE tag"
            );
            let Some(window) = window.as_deref().map(parse_window) else {
                no_window += 1;
                continue;
            };
            if f.is_unmapped() {
                unmapped += 1;
                continue;
            }
            let id = record.reference_sequence_id().expect("tid").expect("tid");
            let (contig, _) = header.reference_sequences().get_index(id).expect("contig");
            let start = record.alignment_start().expect("pos").expect("pos").get();
            let span = record.cigar().alignment_span().expect("cigar");
            let end = start + span.max(1) - 1;
            let inside = contig.as_slice() == window.contig.as_bytes()
                && f.is_reverse_complemented() == window.reverse
                && start <= window.end
                && end >= window.start;
            placed.insert(key, (f.is_properly_segmented(), inside));
            if inside { in_window += 1 } else { outside += 1 }
        }
        assert_eq!(primaries, windows.len() as u64, "one primary per consensus read");
        assert_eq!(realign_count(&stderr, "Consensus reads realigned:"), primaries);
        assert_eq!(realign_count(&stderr, "Unmapped:"), unmapped);
        let guarded_mate_in = realign_count(
            &stderr,
            "Primary kept outside window (proper pair with an in-window mate):",
        );
        let guarded_mate_out = realign_count(
            &stderr,
            "Primary kept outside window (proper pair whose mate cannot move into its window):",
        );
        let (max_mate_in, max_mate_out) = guard_bounds(&placed);
        assert!(guarded_mate_in <= max_mate_in, "{guarded_mate_in} > {max_mate_in}");
        assert!(guarded_mate_out <= max_mate_out, "{guarded_mate_out} > {max_mate_out}");
        let kept_proper_pair = guarded_mate_in + guarded_mate_out;
        assert_eq!(
            realign_count(
                &stderr,
                "Primary outside window (no in-window XA hit or supplementary piece):"
            ) + kept_proper_pair,
            outside,
            "out-of-window primaries are outside or kept by the proper-pair guard"
        );
        assert_eq!(realign_count(&stderr, "Without a window:"), no_window);
        assert_eq!(no_window, windowless, "every consensus end without rw is counted");
        let kept = realign_count(&stderr, "Primary already in window:");
        let swapped = realign_count(&stderr, "Primary swapped to an in-window XA hit:");
        let swapped_supplementary =
            realign_count(&stderr, "Primary swapped to an in-window supplementary piece:");
        assert_eq!(
            kept + swapped + swapped_supplementary,
            in_window,
            "in-window primaries are kept or swapped"
        );
        eprintln!(
            "{backend:?}: {primaries} primaries, {kept} kept, {swapped} swapped, {outside} \
             outside, {unmapped} unmapped, {no_window} without a window"
        );
        assert!(in_window > 0, "no consensus read landed in its window");
        // chr2 repeats a chr1 segment, so some reads bwa places on the copy are
        // swapped back into their window.
        assert!(swapped > 0, "the duplicated segment forced no swap");
    }

    /// Writes into `dir` a copy of the single-line reference at `reference`:
    /// the same contigs, wrapped at 60 bases, with every other 50-base stretch
    /// lowercase (soft-masked), and its own `.fai`, `.dict` and bwa-mem3 index.
    fn soft_masked_wrapped_copy(reference: &Path, dir: &Path) -> PathBuf {
        std::fs::create_dir_all(dir).expect("mkdir");
        let mut reader =
            noodles::fasta::io::reader::Builder.build_from_path(reference).expect("open fasta");
        let mut text = String::new();
        for record in reader.records() {
            let record = record.expect("fasta record");
            let name = std::str::from_utf8(record.name()).expect("utf-8 name");
            writeln!(text, ">{name}").expect("write");
            let masked: Vec<u8> = record
                .sequence()
                .as_ref()
                .iter()
                .enumerate()
                .map(|(i, b)| if (i / 50) % 2 == 1 { b.to_ascii_lowercase() } else { *b })
                .collect();
            for line in masked.chunks(60) {
                writeln!(text, "{}", std::str::from_utf8(line).expect("ascii")).expect("write");
            }
        }
        let fasta = dir.join("ref.fa");
        std::fs::write(&fasta, text).expect("write fasta");
        let index = noodles::fasta::fs::index(&fasta).expect("index fasta");
        noodles::fasta::fai::fs::write(dir.join("ref.fa.fai"), &index).expect("write fai");
        std::fs::copy(reference.with_extension("dict"), dir.join("ref.dict")).expect("copy dict");
        build_aligner_index(&fasta, "bwa-mem3");
        fasta
    }

    /// The inputs of the reference-mode parity tests.
    struct ReferenceModeInputs {
        /// The simulated molecules after `group --realign-window`.
        grouped: PathBuf,
        /// The simulated molecules after `group --realign-window` and `simplex
        /// --min-reads 1`.
        consensus: PathBuf,
        /// A wrapped, soft-masked copy of the reference (see
        /// [`soft_masked_wrapped_copy`]), so the indexed reader must skip line
        /// ends and keep case exactly as the in-memory one does.
        reference: PathBuf,
    }

    fn reference_mode_inputs(dir: &Path) -> ReferenceModeInputs {
        let reference = write_reference(dir);
        build_aligner_index(&reference, "bwa-mem3");
        let mapped = simulate_with_unmapped_mates(dir, &reference, 7);
        let grouped = dir.join("grouped.bam");
        let consensus = dir.join("consensus.bam");
        run_ok(&[
            "group",
            "-i",
            p(&mapped),
            "-o",
            p(&grouped),
            "--strategy",
            "adjacency",
            "--edits",
            "1",
            "--realign-window",
        ]);
        run_ok(&["simplex", "-i", p(&grouped), "-o", p(&consensus), "--min-reads", "1"]);
        let reference = soft_masked_wrapped_copy(&reference, &dir.join("masked"));
        ReferenceModeInputs { grouped, consensus, reference }
    }

    /// By default realign reads the reference regions it needs through the
    /// FASTA's `.fai` index, and writes the same records as with
    /// `--in-memory-reference`. The input forces promotions (the duplicated
    /// segment), so the reference is read.
    #[rstest]
    #[case::subprocess(Backend::Subprocess)]
    #[cfg_attr(feature = "aligner-bwa-mem3", case::inproc(Backend::InProcess))]
    fn indexed_reference_matches_in_memory(#[case] backend: Backend) {
        if !bwa_mem3_or_skip(&format!("indexed_reference_matches_in_memory ({backend:?})")) {
            return;
        }
        let tmp = TempDir::new().expect("tempdir");
        let dir = tmp.path();
        let ReferenceModeInputs { consensus, reference, .. } = reference_mode_inputs(dir);
        let preset = match backend {
            Backend::Subprocess => "bwa-mem3",
            Backend::InProcess => "bwa-mem3-inproc",
        };
        let run = |in_memory: &str| {
            let out = dir.join(format!("realigned.in-memory-{in_memory}.bam"));
            let stderr = run_args_and_capture_logs(
                &[
                    "realign",
                    "-i",
                    p(&consensus),
                    "-o",
                    p(&out),
                    "-r",
                    p(&reference),
                    "--aligner::preset",
                    preset,
                    "--threads",
                    "4",
                    "--in-memory-reference",
                    in_memory,
                ]
                .map(std::ffi::OsStr::new),
            );
            (crate::helpers::read_bam_output(&out).1, stderr)
        };
        let (expected, _) = run("true");
        let (actual, stderr) = run("false");
        assert!(stderr.contains("on demand through its .fai index"), "{stderr}");
        assert!(
            realign_count(&stderr, "Primary swapped to an in-window XA hit:") > 0,
            "no end was promoted, so the reference was never read"
        );
        assert!(!expected.is_empty(), "realign produced no records");
        assert_eq!(actual.len(), expected.len(), "record counts differ");
        for (i, (a, e)) in actual.iter().zip(&expected).enumerate() {
            assert_eq!(a, e, "record {i} differs between the indexed and in-memory reference");
        }
    }

    fn preset(backend: Backend) -> &'static str {
        match backend {
            Backend::Subprocess => "bwa-mem3",
            Backend::InProcess => "bwa-mem3-inproc",
        }
    }

    /// The decoded records of a BAM (no header, so `@PG` is not compared).
    fn records(path: &Path) -> Vec<noodles::sam::alignment::RecordBuf> {
        crate::helpers::read_bam_output(path).1
    }

    /// The same simulated grouped-input through `runall --start-from group
    /// --stop-after consensus --realign` and through the standalone chain
    /// `group --realign-window → simplex → realign → sort --order coordinate`
    /// gives identical records (`@PG` aside: one runall `@PG` vs one per
    /// command).
    #[rstest]
    #[case::subprocess(Backend::Subprocess)]
    #[cfg_attr(feature = "aligner-bwa-mem3", case::inproc(Backend::InProcess))]
    fn runall_matches_standalone_chain(#[case] backend: Backend) {
        if !bwa_mem3_or_skip(&format!("runall_matches_standalone_chain ({backend:?})")) {
            return;
        }
        let tmp = TempDir::new().expect("tempdir");
        let dir = tmp.path();
        let reference = write_reference(dir);
        build_aligner_index(&reference, "bwa-mem3");
        let mapped = simulate_with_unmapped_mates(dir, &reference, 7);

        // Standalone chain.
        let grouped = dir.join("grouped.bam");
        let consensus = dir.join("consensus.bam");
        let realigned = dir.join("realigned.bam");
        let standalone = dir.join("standalone.bam");
        run_ok(&[
            "group",
            "-i",
            p(&mapped),
            "-o",
            p(&grouped),
            "--strategy",
            "adjacency",
            "--edits",
            "1",
            "--realign-window",
        ]);
        run_ok(&["simplex", "-i", p(&grouped), "-o", p(&consensus), "--min-reads", "1"]);
        run_ok(&[
            "realign",
            "-i",
            p(&consensus),
            "-o",
            p(&realigned),
            "-r",
            p(&reference),
            "--aligner::preset",
            preset(backend),
            "--threads",
            "4",
        ]);
        run_ok(&["sort", "-i", p(&realigned), "-o", p(&standalone), "--order", "coordinate"]);

        // runall.
        let fused = dir.join("runall.bam");
        run_ok(&[
            "runall",
            "--start-from",
            "group",
            "--stop-after",
            "consensus",
            "--consensus",
            "simplex",
            "-i",
            p(&mapped),
            "-o",
            p(&fused),
            "--group::strategy",
            "adjacency",
            "--group::edits",
            "1",
            "--simplex::min-reads",
            "1",
            "--ref",
            p(&reference),
            "--aligner::preset",
            preset(backend),
            "--realign",
            "--threads",
            "4",
        ]);

        let expected = records(&standalone);
        assert!(!expected.is_empty(), "the standalone chain produced no records");
        let actual = records(&fused);
        assert_eq!(actual.len(), expected.len(), "record counts differ");
        for (i, (a, e)) in actual.iter().zip(&expected).enumerate() {
            assert_eq!(
                a, e,
                "record {i} differs between runall --realign and the standalone chain"
            );
        }
        let (header, _) = crate::helpers::read_bam_output(&fused);
        let sort_order = header
            .header()
            .and_then(|hd| {
                hd.other_fields()
                    .get(&noodles::sam::header::record::value::map::header::tag::SORT_ORDER)
            })
            .map(ToString::to_string);
        assert_eq!(sort_order.as_deref(), Some("coordinate"));
    }

    /// `runall --start-from consensus --realign`, from the grouped reads, writes
    /// the same records by default, when the realign stage reads the reference
    /// through its `.fai` index (no other stage of this chain holds it), as with
    /// `--realign::in-memory-reference`.
    #[rstest]
    #[case::subprocess(Backend::Subprocess)]
    #[cfg_attr(feature = "aligner-bwa-mem3", case::inproc(Backend::InProcess))]
    fn runall_indexed_reference_matches_in_memory(#[case] backend: Backend) {
        if !bwa_mem3_or_skip(&format!("runall_indexed_reference_matches_in_memory ({backend:?})")) {
            return;
        }
        let tmp = TempDir::new().expect("tempdir");
        let dir = tmp.path();
        let ReferenceModeInputs { grouped, reference, .. } = reference_mode_inputs(dir);
        let run = |in_memory: &str| {
            let out = dir.join(format!("runall.in-memory-{in_memory}.bam"));
            let stderr = run_args_and_capture_logs(
                &[
                    "runall",
                    "--start-from",
                    "consensus",
                    "--stop-after",
                    "consensus",
                    "--consensus",
                    "simplex",
                    "--simplex::min-reads",
                    "1",
                    "-i",
                    p(&grouped),
                    "-o",
                    p(&out),
                    "--ref",
                    p(&reference),
                    "--aligner::preset",
                    preset(backend),
                    "--realign",
                    "--realign::in-memory-reference",
                    in_memory,
                    "--threads",
                    "4",
                ]
                .map(std::ffi::OsStr::new),
            );
            (records(&out), stderr)
        };
        let (expected, _) = run("true");
        let (actual, stderr) = run("false");
        assert!(stderr.contains("on demand through its .fai index"), "{stderr}");
        assert!(
            realign_count(&stderr, "Primary swapped to an in-window XA hit:") > 0,
            "no end was promoted, so the reference was never read"
        );
        assert!(!expected.is_empty(), "runall --realign produced no records");
        assert_eq!(actual.len(), expected.len(), "record counts differ");
        for (i, (a, e)) in actual.iter().zip(&expected).enumerate() {
            assert_eq!(a, e, "record {i} differs between the indexed and in-memory reference");
        }
    }

    /// `runall --start-from extract … --realign`: group runs on the align
    /// stage's output, so its `rw` contig names must come from the aligned
    /// header. Checked twice: with `--stop-after group
    /// --group::realign-window` every raw primary's window names its own
    /// contig, and with `--realign` every consensus primary bwa placed in its
    /// window is counted as kept or swapped by realign (a wrong contig name
    /// would put no window on any contig). The in-process case also runs the
    /// raw-read align and the realign against one shared index.
    #[rstest]
    #[case::subprocess(Backend::Subprocess)]
    #[cfg_attr(feature = "aligner-bwa-mem3", case::inproc(Backend::InProcess))]
    fn runall_from_extract_realigns_with_the_aligned_header(#[case] backend: Backend) {
        if !bwa_mem3_or_skip(&format!(
            "runall_from_extract_realigns_with_the_aligned_header ({backend:?})"
        )) {
            return;
        }
        let tmp = TempDir::new().expect("tempdir");
        let dir = tmp.path();
        let reference = write_reference(dir);
        build_aligner_index(&reference, "bwa-mem3");
        let (r1, r2) = (dir.join("r1.fq.gz"), dir.join("r2.fq.gz"));
        run_ok(&[
            "simulate",
            "fastq-reads",
            "-1",
            p(&r1),
            "-2",
            p(&r2),
            "--truth",
            p(&dir.join("truth.tsv")),
            "-r",
            p(&reference),
            "--num-molecules",
            "200",
            "--read-length",
            "100",
            "--seed",
            "11",
        ]);
        let common = |stop: &str, out: &Path| -> Vec<String> {
            [
                "runall",
                "--start-from",
                "extract",
                "--stop-after",
                stop,
                "--extract::inputs",
                p(&r1),
                p(&r2),
                "--extract::read-structures",
                "8M+T",
                "+T",
                "--extract::sample",
                "s1",
                "--extract::library",
                "lib1",
                "--ref",
                p(&reference),
                "--aligner::preset",
                preset(backend),
                "--group::strategy",
                "adjacency",
                "--group::edits",
                "1",
                "--threads",
                "4",
                "-o",
                p(out),
            ]
            .map(String::from)
            .to_vec()
        };

        // 1. Group's windows name the contig each raw primary is on.
        let grouped = dir.join("grouped.bam");
        let mut args = common("group", &grouped);
        args.push("--group::realign-window".to_string());
        run_ok(&args.iter().map(String::as_str).collect::<Vec<_>>());
        let mut reader = noodles::bam::io::reader::Builder.build_from_path(&grouped).expect("open");
        let header = reader.read_header().expect("header");
        assert!(!header.reference_sequences().is_empty(), "grouped output has no @SQ");
        let mut windows: HashMap<(String, bool), String> = HashMap::new();
        let mut stamped = 0usize;
        for record in reader.records() {
            let record = record.expect("record");
            let f = record.flags();
            if f.is_unmapped() || f.is_secondary() || f.is_supplementary() {
                continue;
            }
            let rw = string_tag(&record, SamTag::RW).expect("raw primaries carry rw");
            let id = record.reference_sequence_id().expect("tid").expect("tid");
            let (contig, _) = header.reference_sequences().get_index(id).expect("contig");
            assert_eq!(
                parse_window(&rw).contig.as_bytes(),
                contig.as_slice(),
                "rw {rw} does not name the record's contig"
            );
            let mi = string_tag(&record, SamTag::MI).expect("MI");
            windows.entry((mi, f.is_first_segment())).or_insert(rw);
            stamped += 1;
        }
        assert!(stamped > 0, "no raw primary carried rw");

        // 2. The full run re-aligns the consensus reads against those windows.
        let realigned = dir.join("realigned.bam");
        let mut args = common("consensus", &realigned);
        args.extend(
            ["--consensus", "simplex", "--simplex::min-reads", "1", "--realign"].map(String::from),
        );
        let stderr =
            run_args_and_capture_logs(&args.iter().map(std::ffi::OsStr::new).collect::<Vec<_>>());
        if matches!(backend, Backend::InProcess) {
            assert_eq!(
                stderr.matches("in-process bwa-mem3: loaded index").count(),
                1,
                "the align and realign stages must load the index once:\n{stderr}"
            );
            assert!(stderr.contains("reusing the index"), "{stderr}");
        }
        let mut reader =
            noodles::bam::io::reader::Builder.build_from_path(&realigned).expect("open");
        let header = reader.read_header().expect("header");
        let (mut primaries, mut in_window, mut outside, mut unmapped) = (0u64, 0u64, 0u64, 0u64);
        let mut last_pos: Option<(usize, usize)> = None;
        let mut placed: HashMap<(String, u8), (bool, bool)> = HashMap::new();
        for record in reader.records() {
            let record = record.expect("record");
            assert!(string_tag(&record, SamTag::RW).is_none(), "rw must be stripped");
            let f = record.flags();
            if let (Some(Ok(id)), Some(Ok(pos))) =
                (record.reference_sequence_id(), record.alignment_start())
            {
                let key = (id, pos.get());
                assert!(last_pos.is_none_or(|last| last <= key), "output is not coordinate sorted");
                last_pos = Some(key);
            }
            if f.is_secondary() || f.is_supplementary() {
                continue;
            }
            primaries += 1;
            if f.is_unmapped() {
                unmapped += 1;
                continue;
            }
            let mi = string_tag(&record, SamTag::MI).expect("consensus MI");
            let Some(rw) = windows.get(&(mi.clone(), f.is_first_segment())) else {
                panic!("no raw window for consensus MI {mi}");
            };
            let window = parse_window(rw);
            let id = record.reference_sequence_id().expect("tid").expect("tid");
            let (contig, _) = header.reference_sequences().get_index(id).expect("contig");
            let start = record.alignment_start().expect("pos").expect("pos").get();
            let span = record.cigar().alignment_span().expect("cigar");
            let end = start + span.max(1) - 1;
            let inside = contig.as_slice() == window.contig.as_bytes()
                && f.is_reverse_complemented() == window.reverse
                && start <= window.end
                && end >= window.start;
            let seg = if f.is_first_segment() { 1 } else { 2 };
            placed.insert((mi, seg), (f.is_properly_segmented(), inside));
            if inside { in_window += 1 } else { outside += 1 }
        }
        assert!(primaries > 0, "no consensus reads");
        assert_eq!(realign_count(&stderr, "Consensus reads realigned:"), primaries);
        assert_eq!(realign_count(&stderr, "Unmapped:"), unmapped);
        let kept = realign_count(&stderr, "Primary already in window:");
        let swapped = realign_count(&stderr, "Primary swapped to an in-window XA hit:");
        let swapped_supplementary =
            realign_count(&stderr, "Primary swapped to an in-window supplementary piece:");
        assert_eq!(
            kept + swapped + swapped_supplementary,
            in_window,
            "in-window primaries are kept or swapped"
        );
        let guarded_mate_in = realign_count(
            &stderr,
            "Primary kept outside window (proper pair with an in-window mate):",
        );
        let guarded_mate_out = realign_count(
            &stderr,
            "Primary kept outside window (proper pair whose mate cannot move into its window):",
        );
        let (max_mate_in, max_mate_out) = guard_bounds(&placed);
        assert!(guarded_mate_in <= max_mate_in, "{guarded_mate_in} > {max_mate_in}");
        assert!(guarded_mate_out <= max_mate_out, "{guarded_mate_out} > {max_mate_out}");
        let kept_proper_pair = guarded_mate_in + guarded_mate_out;
        assert_eq!(
            realign_count(
                &stderr,
                "Primary outside window (no in-window XA hit or supplementary piece):"
            ) + kept_proper_pair,
            outside,
            "out-of-window primaries are outside or kept by the proper-pair guard"
        );
        assert!(
            in_window * 10 >= (primaries - unmapped) * 9,
            "{in_window} of {} mapped consensus primaries in their window",
            primaries - unmapped
        );
    }

    /// [`super::UNMAPPED_ALIGNER`] with the `@SQ` lines of [`write_reference`]'s
    /// contigs: reports every read unmapped, so the realign stage after an
    /// intermediate filter passes the filtered records through without a real
    /// aligner.
    const UNMAPPED_ALIGNER_TWO_CONTIGS: &str = "awk 'BEGIN { OFS = \"\\t\"; \
        print \"@HD\", \"VN:1.6\"; print \"@SQ\", \"SN:chr1\", \"LN:20000\"; \
        print \"@SQ\", \"SN:chr2\", \"LN:10000\" } \
        NR % 4 == 1 { name = substr($1, 2) } NR % 4 == 2 { seq = $0 } \
        NR % 4 == 0 { n++; print name, (n % 2 ? 77 : 141), \"*\", 0, 0, \"*\", \"*\", 0, 0, seq, $0 }' \
        # {ref}";

    /// How the intermediate filter is reached: the run starts at filter (its
    /// input is a consensus BAM), or it calls consensus and continues to filter.
    #[derive(Debug, Clone, Copy)]
    enum FilterForm {
        StartFromFilter,
        StopAfterFilter,
    }

    /// The filter options both sides of a filter comparison use (without the
    /// `--filter::` prefix): enough `--min-reads` to reject some molecules.
    const FILTER_ARGS: [&str; 2] = ["min-reads", "3"];

    /// The window-grouped simplex consensus of `mapped` (`group --realign-window
    /// → simplex --min-reads 1`), written to `dir/consensus.bam`.
    fn window_consensus(dir: &Path, mapped: &Path) -> PathBuf {
        let grouped = dir.join("grouped.bam");
        let consensus = dir.join("consensus.bam");
        run_ok(&[
            "group",
            "-i",
            p(mapped),
            "-o",
            p(&grouped),
            "--strategy",
            "adjacency",
            "--edits",
            "1",
            "--realign-window",
        ]);
        run_ok(&["simplex", "-i", p(&grouped), "-o", p(&consensus), "--min-reads", "1"]);
        consensus
    }

    /// `runall` arguments that run `form` up to and including filter (kept reads
    /// to `out`, rejects to `rejects`, stats to `stats`); with `realign`,
    /// `--realign` and its aligner (`aligner_args`) are added.
    #[allow(clippy::too_many_arguments)]
    fn runall_filter_args(
        form: FilterForm,
        mapped: &Path,
        consensus: &Path,
        out: &Path,
        rejects: &Path,
        stats: &Path,
        by_template: bool,
        realign: Option<&[&str]>,
    ) -> Vec<String> {
        let mut args: Vec<String> = match form {
            FilterForm::StartFromFilter => {
                vec![
                    "runall",
                    "--start-from",
                    "filter",
                    "--stop-after",
                    "filter",
                    "-i",
                    p(consensus),
                ]
            }
            FilterForm::StopAfterFilter => vec![
                "runall",
                "--start-from",
                "group",
                "--stop-after",
                "filter",
                "--consensus",
                "simplex",
                "-i",
                p(mapped),
                "--group::strategy",
                "adjacency",
                "--group::edits",
                "1",
                "--group::realign-window",
                "--simplex::min-reads",
                "1",
            ],
        }
        .into_iter()
        .map(String::from)
        .collect();
        args.extend(
            ["-o", p(out), "--filter::rejects", p(rejects), "--filter::stats", p(stats)]
                .map(String::from),
        );
        args.push(format!("--filter::{}", FILTER_ARGS[0]));
        args.push(FILTER_ARGS[1].to_string());
        args.push("--filter::filter-by-template".to_string());
        args.push(by_template.to_string());
        if let Some(aligner) = realign {
            args.push("--realign".to_string());
            args.extend(aligner.iter().map(|a| (*a).to_string()));
        }
        args
    }

    /// A record as filter leaves it: name, flags, bases, qualities and tags,
    /// minus `rw` (realign strips it). The alignment fields are left out: the
    /// realign stage's aligner sets them.
    type FilteredRecord = (String, u16, Vec<u8>, Vec<u8>, Vec<(String, String)>);

    fn filtered_record(record: &noodles::sam::alignment::RecordBuf) -> FilteredRecord {
        let rw = SamTag::RW.to_noodles_tag();
        let mut tags: Vec<(String, String)> = record
            .data()
            .iter()
            .filter(|(tag, _)| *tag != rw)
            .map(|(tag, value)| {
                (String::from_utf8_lossy(tag.as_ref()).into_owned(), format!("{value:?}"))
            })
            .collect();
        tags.sort();
        (
            record.name().map(ToString::to_string).unwrap_or_default(),
            record.flags().bits(),
            record.sequence().as_ref().to_vec(),
            record.quality_scores().as_ref().to_vec(),
            tags,
        )
    }

    /// The records of `path` as [`filtered_record`]s, sorted (the realigned
    /// output is coordinate sorted, filter's is in input order).
    fn filtered_records(path: &Path) -> Vec<FilteredRecord> {
        let mut records: Vec<_> = records(path).iter().map(filtered_record).collect();
        records.sort();
        records
    }

    /// The intermediate filter keeps exactly the records the terminal filter
    /// keeps, writes the same rejects and the same stats, for both ways of
    /// reaching it and both filter modes. The intermediate filter's output is
    /// observed through the realign stage after it, with an aligner that
    /// reports every read unmapped: realign then hands each filtered record
    /// through (it strips `rw` and takes alignment fields from the aligner),
    /// so the kept records' names, flags, bases, qualities (masking) and tags
    /// are the intermediate filter's. Rejects and stats are written by the
    /// filter stage itself and are compared exactly. The terminal side is
    /// `fgumi filter` on the same consensus input (`--start-from filter`) or
    /// `runall --stop-after filter` without `--realign` (`--stop-after filter`,
    /// with `--group::realign-window` so the consensus reads match). Needs no
    /// real aligner.
    #[rstest]
    fn intermediate_filter_matches_terminal_filter(
        #[values(FilterForm::StartFromFilter, FilterForm::StopAfterFilter)] form: FilterForm,
        #[values(true, false)] by_template: bool,
    ) {
        let tmp = TempDir::new().expect("tempdir");
        let dir = tmp.path();
        let reference = write_reference(dir);
        let mapped = simulate_with_unmapped_mates(dir, &reference, 7);
        let consensus = window_consensus(dir, &mapped);

        let (kept, rejects, stats) =
            (dir.join("terminal.bam"), dir.join("terminal.rej.bam"), dir.join("terminal.stats"));
        match form {
            FilterForm::StartFromFilter => run_ok(&[
                "filter",
                "-i",
                p(&consensus),
                "-o",
                p(&kept),
                "--rejects",
                p(&rejects),
                "--stats",
                p(&stats),
                &format!("--{}", FILTER_ARGS[0]),
                FILTER_ARGS[1],
                "--filter-by-template",
                &by_template.to_string(),
            ]),
            FilterForm::StopAfterFilter => {
                let args = runall_filter_args(
                    form,
                    &mapped,
                    &consensus,
                    &kept,
                    &rejects,
                    &stats,
                    by_template,
                    None,
                );
                run_ok(&args.iter().map(String::as_str).collect::<Vec<_>>());
            }
        }

        let (out, out_rejects, out_stats) =
            (dir.join("realign.bam"), dir.join("realign.rej.bam"), dir.join("realign.stats"));
        let aligner = ["--ref", p(&reference), "--aligner::command", UNMAPPED_ALIGNER_TWO_CONTIGS];
        let args = runall_filter_args(
            form,
            &mapped,
            &consensus,
            &out,
            &out_rejects,
            &out_stats,
            by_template,
            Some(&aligner),
        );
        run_ok(&args.iter().map(String::as_str).collect::<Vec<_>>());

        // The realign stage's trailing sort sorts by coordinate.
        let (header, _) = crate::helpers::read_bam_output(&out);
        let so = header
            .header()
            .and_then(|hd| {
                hd.other_fields()
                    .get(&noodles::sam::header::record::value::map::header::tag::SORT_ORDER)
            })
            .map(ToString::to_string);
        assert_eq!(so.as_deref(), Some("coordinate"), "the --realign output is coordinate sorted");

        let expected = filtered_records(&kept);
        assert!(!expected.is_empty(), "the terminal filter kept nothing");
        let expected_rejects = records(&rejects);
        assert!(!expected_rejects.is_empty(), "the terminal filter rejected nothing");
        assert_eq!(filtered_records(&out), expected, "kept records differ");
        assert_eq!(records(&out_rejects), expected_rejects, "rejects differ");
        assert_eq!(
            std::fs::read_to_string(&out_stats).expect("stats"),
            std::fs::read_to_string(&stats).expect("stats"),
            "filter stats differ"
        );
    }

    /// `runall … filter --realign` runs end to end with a real aligner and
    /// matches the standalone chain `group --realign-window → simplex → filter
    /// → realign → sort --order coordinate` (records identical, `@PG` aside),
    /// for both ways of reaching filter: `--stop-after filter` from the raw
    /// reads, and `--start-from filter` on the standalone chain's consensus.
    /// The output is coordinate sorted, carries no `rw`, every consensus
    /// primary is counted by realign, and the filter's rejects match the
    /// standalone filter's.
    #[rstest]
    #[case::stop_after_filter_subprocess(FilterForm::StopAfterFilter, Backend::Subprocess)]
    #[case::start_from_filter_subprocess(FilterForm::StartFromFilter, Backend::Subprocess)]
    #[cfg_attr(
        feature = "aligner-bwa-mem3",
        case::stop_after_filter_inproc(FilterForm::StopAfterFilter, Backend::InProcess)
    )]
    #[cfg_attr(
        feature = "aligner-bwa-mem3",
        case::start_from_filter_inproc(FilterForm::StartFromFilter, Backend::InProcess)
    )]
    fn runall_filter_realign_matches_standalone_chain(
        #[case] form: FilterForm,
        #[case] backend: Backend,
    ) {
        if !bwa_mem3_or_skip(&format!(
            "runall_filter_realign_matches_standalone_chain ({form:?}, {backend:?})"
        )) {
            return;
        }
        let tmp = TempDir::new().expect("tempdir");
        let dir = tmp.path();
        let reference = write_reference(dir);
        build_aligner_index(&reference, "bwa-mem3");
        let mapped = simulate_with_unmapped_mates(dir, &reference, 7);

        // Standalone chain.
        let consensus = window_consensus(dir, &mapped);
        let (filtered, rejects, stats) =
            (dir.join("filtered.bam"), dir.join("filter.rej.bam"), dir.join("filter.stats"));
        run_ok(&[
            "filter",
            "-i",
            p(&consensus),
            "-o",
            p(&filtered),
            "--rejects",
            p(&rejects),
            "--stats",
            p(&stats),
            &format!("--{}", FILTER_ARGS[0]),
            FILTER_ARGS[1],
        ]);
        let realigned = dir.join("realigned.bam");
        let standalone = dir.join("standalone.bam");
        run_ok(&[
            "realign",
            "-i",
            p(&filtered),
            "-o",
            p(&realigned),
            "-r",
            p(&reference),
            "--aligner::preset",
            preset(backend),
            "--threads",
            "4",
        ]);
        run_ok(&["sort", "-i", p(&realigned), "-o", p(&standalone), "--order", "coordinate"]);

        // runall.
        let (fused, fused_rejects, fused_stats) =
            (dir.join("runall.bam"), dir.join("runall.rej.bam"), dir.join("runall.stats"));
        let aligner =
            ["--ref", p(&reference), "--aligner::preset", preset(backend), "--threads", "4"];
        let args = runall_filter_args(
            form,
            &mapped,
            &consensus,
            &fused,
            &fused_rejects,
            &fused_stats,
            true,
            Some(&aligner),
        );
        let stderr =
            run_args_and_capture_logs(&args.iter().map(std::ffi::OsStr::new).collect::<Vec<_>>());

        let expected = records(&standalone);
        assert!(!expected.is_empty(), "the standalone chain produced no records");
        let actual = records(&fused);
        assert_eq!(actual.len(), expected.len(), "record counts differ");
        for (i, (a, e)) in actual.iter().zip(&expected).enumerate() {
            assert_eq!(a, e, "record {i} differs between runall and the standalone chain");
        }
        let expected_rejects = records(&rejects);
        assert!(!expected_rejects.is_empty(), "the filter rejected nothing");
        assert_eq!(records(&fused_rejects), expected_rejects, "rejects differ");
        assert_eq!(
            std::fs::read_to_string(&fused_stats).expect("stats"),
            std::fs::read_to_string(&stats).expect("stats"),
            "filter stats differ"
        );

        let (header, _) = crate::helpers::read_bam_output(&fused);
        let sort_order = header
            .header()
            .and_then(|hd| {
                hd.other_fields()
                    .get(&noodles::sam::header::record::value::map::header::tag::SORT_ORDER)
            })
            .map(ToString::to_string);
        assert_eq!(sort_order.as_deref(), Some("coordinate"));

        let mut reader = noodles::bam::io::reader::Builder.build_from_path(&fused).expect("open");
        reader.read_header().expect("header");
        let mut primaries = 0u64;
        let mut last_pos: Option<(usize, usize)> = None;
        for record in reader.records() {
            let record = record.expect("record");
            assert!(string_tag(&record, SamTag::RW).is_none(), "rw must be stripped");
            if let (Some(Ok(id)), Some(Ok(pos))) =
                (record.reference_sequence_id(), record.alignment_start())
            {
                let key = (id, pos.get());
                assert!(last_pos.is_none_or(|last| last <= key), "output is not coordinate sorted");
                last_pos = Some(key);
            }
            let f = record.flags();
            if !(f.is_secondary() || f.is_supplementary()) {
                primaries += 1;
            }
        }
        assert_eq!(realign_count(&stderr, "Consensus reads realigned:"), primaries);
        let counted = [
            "Primary already in window:",
            "Primary swapped to an in-window XA hit:",
            "Primary swapped to an in-window supplementary piece:",
            "Primary kept outside window (proper pair with an in-window mate):",
            "Primary kept outside window (proper pair whose mate cannot move into its window):",
            "Primary outside window (no in-window XA hit or supplementary piece):",
            "Unmapped:",
            "Without a window:",
        ]
        .into_iter()
        .map(|label| realign_count(&stderr, label))
        .sum::<u64>();
        assert_eq!(counted, primaries, "every consensus primary is in its window or counted");
    }

    /// The fused `group → simplex` chain carries each segment's raw-read window
    /// onto the consensus reads exactly as the standalone `group
    /// --realign-window → simplex` does: `runall --realign` reads `rw` off the
    /// consensus reads it calls itself. Needs no aligner.
    #[test]
    fn fused_group_to_consensus_carries_rw() {
        let tmp = TempDir::new().expect("tempdir");
        let dir = tmp.path();
        let reference = write_reference(dir);
        let mapped = dir.join("mapped.bam");
        run_ok(&[
            "simulate",
            "mapped-reads",
            "-r",
            p(&reference),
            "-o",
            p(&mapped),
            "--truth",
            p(&dir.join("truth.tsv")),
            "--num-molecules",
            "100",
            "--read-length",
            "100",
            "--seed",
            "3",
        ]);
        let grouped = dir.join("grouped.bam");
        let standalone = dir.join("standalone.bam");
        run_ok(&[
            "group",
            "-i",
            p(&mapped),
            "-o",
            p(&grouped),
            "--strategy",
            "adjacency",
            "--realign-window",
        ]);
        run_ok(&["simplex", "-i", p(&grouped), "-o", p(&standalone), "--min-reads", "1"]);
        let fused = dir.join("fused.bam");
        run_ok(&[
            "runall",
            "--start-from",
            "group",
            "--stop-after",
            "consensus",
            "--consensus",
            "simplex",
            "-i",
            p(&mapped),
            "-o",
            p(&fused),
            "--group::strategy",
            "adjacency",
            "--group::realign-window",
            "--simplex::min-reads",
            "1",
        ]);
        let expected = records(&standalone);
        assert!(!expected.is_empty(), "simplex produced no consensus reads");
        let rw = SamTag::RW.to_noodles_tag();
        assert!(expected.iter().all(|r| r.data().get(&rw).is_some()), "standalone lacks rw");
        assert_eq!(records(&fused), expected);
    }
}

/// `runall --realign` sorts the re-aligned reads by coordinate whatever
/// `--sort::order` says. It warns that the order is ignored only when the user
/// passed one, never on a default run.
#[rstest::rstest]
#[case::default_order(&[], false)]
#[case::explicit_order(&["--sort::order", "queryname"], true)]
fn runall_realign_warns_only_about_an_explicit_sort_order(
    #[case] extra: &[&str],
    #[case] warns: bool,
) {
    let tmp = TempDir::new().expect("tempdir");
    let dir = tmp.path();
    let input = write_consensus(dir, 2, |_| (Some("chr1:1-100:+"), Some("chr1:1-100:-")));
    let reference = create_test_reference(dir);
    let output = dir.join("out.bam");
    let result = fgumi()
        .args(["runall", "--start-from", "filter", "--stop-after", "filter", "--realign"])
        .args(["--filter::min-reads", "1", "--ref", p(&reference)])
        .args(["-i", p(&input), "-o", p(&output), "--aligner::command", UNMAPPED_ALIGNER])
        .args(extra)
        .env("RUST_LOG", "info")
        .output()
        .expect("spawn fgumi");
    let stderr = String::from_utf8_lossy(&result.stderr);
    assert!(result.status.success(), "{stderr}");
    assert_eq!(realign_count(&stderr, "Consensus reads realigned:"), 4, "the chain ran");
    assert_eq!(stderr.contains("--sort::order"), warns, "{stderr}");
}

/// In `runall --consensus simplex --methylation-mode em-seq --realign`, the
/// methylation-aware consensus caller loads `--ref` into memory, and realign
/// shares that copy rather than reading the (indexed) FASTA on demand or
/// loading a second copy.
#[test]
fn runall_realign_shares_the_methylation_callers_reference() {
    let tmp = TempDir::new().expect("tempdir");
    let dir = tmp.path();
    let read = |name: &str, seg: u16, pos: i32, mate_pos: i32, reverse: bool| {
        let mut b = SamBuilder::new();
        b.read_name(name.as_bytes())
            .flags(
                flags::PAIRED
                    | flags::PROPER_PAIR
                    | seg
                    | if reverse { flags::REVERSE } else { flags::MATE_REVERSE },
            )
            .ref_id(0)
            .pos(pos - 1)
            .mapq(60)
            .cigar_ops(&[16 << 4])
            .mate_ref_id(0)
            .mate_pos(mate_pos - 1)
            .template_length(if reverse { -116 } else { 116 })
            .sequence(b"ACGTACGTACGTACGT")
            .qualities(&[40; 16]);
        b.add_string_tag(SamTag::RX, b"AAAA").add_string_tag(SamTag::MC, b"16M");
        b.build()
    };
    let records: Vec<_> = (0..3)
        .flat_map(|i| {
            let name = format!("a{i}");
            [
                read(&name, flags::FIRST_SEGMENT, 101, 201, false),
                read(&name, flags::LAST_SEGMENT, 201, 101, true),
            ]
        })
        .collect();
    let raw = dir.join("raw.bam");
    write_bam(&raw, &create_minimal_header("chr1", 10_000), &records);
    let reference = create_test_reference(dir);
    let (sorted, output) = (dir.join("sorted.bam"), dir.join("out.bam"));
    let sort = fgumi()
        .args(["sort", "-i", p(&raw), "-o", p(&sorted), "--order", "template-coordinate"])
        .output()
        .expect("spawn fgumi sort");
    assert!(sort.status.success(), "{}", String::from_utf8_lossy(&sort.stderr));
    let result = fgumi()
        .args(["runall", "--start-from", "group", "--stop-after", "consensus", "--realign"])
        .args(["--consensus", "simplex", "--methylation-mode", "em-seq", "--ref", p(&reference)])
        .args(["--group::strategy", "adjacency", "--simplex::min-reads", "1"])
        .args(["-i", p(&sorted), "-o", p(&output), "--aligner::command", UNMAPPED_ALIGNER])
        .env("RUST_LOG", "info")
        .output()
        .expect("spawn fgumi");
    let stderr = String::from_utf8_lossy(&result.stderr);
    assert!(result.status.success(), "{stderr}");
    assert_eq!(realign_count(&stderr, "Consensus reads realigned:"), 2, "the chain ran");
    assert!(stderr.contains("Realign: an earlier stage holds"), "{stderr}");
    assert!(stderr.contains("realign shares it"), "{stderr}");
    assert!(!stderr.contains("on demand through"), "{stderr}");
    assert!(!stderr.contains("in-memory reference requested"), "{stderr}");
}
