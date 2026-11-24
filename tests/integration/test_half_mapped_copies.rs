//! Half-mapped pairs that are copies of a fully mapped pair: `sort` places them beside it, `group`
//! gives them its molecule, `dedup` marks them, and consensus uses only their mapped read.
//!
//! Every input is written as SAM and sorted with `fgumi sort --order template-coordinate`.

use std::collections::{BTreeMap, BTreeSet};
use std::ffi::OsStr;
use std::fmt::Write as _;
use std::fs;
use std::path::{Path, PathBuf};

use fgumi_lib::sam::SamTag;
use noodles::sam::alignment::RecordBuf;
use noodles::sam::alignment::record_buf::data::field::Value;
use rstest::rstest;
use tempfile::TempDir;

use crate::helpers::assertions::assert_bam_sorted;
use crate::helpers::bam_generator::write_bam;
use crate::helpers::cli::{run_and_capture_logs, run_args_and_capture_logs, utf8_path as p};
use crate::helpers::cutover::decompressed_records_without_pg;
use crate::helpers::read_bam_output;

const HEADER: &str = "@HD\tVN:1.6\tSO:unsorted\n\
                      @SQ\tSN:chr1\tLN:1000000\n\
                      @SQ\tSN:chr2\tLN:1000000\n\
                      @SQ\tSN:chr3\tLN:1000000\n\
                      @RG\tID:A\tSM:s\tLB:l1\n\
                      @RG\tID:B\tSM:s\tLB:l2\n";

const PAIRED: &[&str] = &["--strategy", "paired"];
const IDENTITY: &[&str] = &["--strategy", "identity"];
const ADJACENCY: &[&str] = &["--strategy", "adjacency"];
const EDIT: &[&str] = &["--strategy", "edit"];
const PARALLEL_PAIRED: &[&str] = &["--strategy", "paired", "--parallel-group-min-templates", "0"];
const PARALLEL_ADJACENCY: &[&str] =
    &["--strategy", "adjacency", "--parallel-group-min-templates", "0"];
const NO_UMI: &[&str] = &["--strategy", "identity", "--no-umi"];

/// Where a read maps: contig, 1-based position, and whether it is on the reverse strand. Reads
/// are 10 bases long, so a reverse read's unclipped 5' position is 9 bases to the right.
type Site = (&'static str, u32, bool);

const LOWER: Site = ("chr1", 501, false);
const UPPER: Site = ("chr1", 701, true);

#[derive(Clone)]
struct Pair {
    name: String,
    r1: Option<Site>,
    r2: Option<Site>,
    umi: String,
    cell: Option<&'static str>,
    read_group: &'static str,
    mapq: u8,
    quality: char,
    extra_flags: u16,
    mate_cigar: bool,
}

impl Pair {
    fn new(name: &str, r1: Option<Site>, r2: Option<Site>, umi: &str) -> Self {
        Self {
            name: name.to_string(),
            r1,
            r2,
            umi: umi.to_string(),
            cell: None,
            read_group: "A",
            mapq: 60,
            quality: 'I',
            extra_flags: 0,
            mate_cigar: true,
        }
    }

    /// A fully mapped pair with R1 at `LOWER` and R2 at `UPPER`.
    fn full(umi: &str) -> Self {
        Self::new("full", Some(LOWER), Some(UPPER), umi)
    }

    /// A fully mapped pair with R1 at `UPPER` and R2 at `LOWER`.
    fn flipped(umi: &str) -> Self {
        Self::new("full", Some(UPPER), Some(LOWER), umi)
    }

    /// A pair with only R1 mapped.
    fn r1_only(site: Site, umi: &str) -> Self {
        Self::new("half", Some(site), None, umi)
    }

    /// A pair with only R2 mapped.
    fn r2_only(site: Site, umi: &str) -> Self {
        Self::new("half", None, Some(site), umi)
    }

    fn named(mut self, name: &str) -> Self {
        self.name = name.to_string();
        self
    }

    fn cell(mut self, cell: &'static str) -> Self {
        self.cell = Some(cell);
        self
    }

    fn read_group(mut self, read_group: &'static str) -> Self {
        self.read_group = read_group;
        self
    }

    fn mapq(mut self, mapq: u8) -> Self {
        self.mapq = mapq;
        self
    }

    fn quality(mut self, quality: char) -> Self {
        self.quality = quality;
        self
    }

    fn extra_flags(mut self, flags: u16) -> Self {
        self.extra_flags = flags;
        self
    }

    fn without_mate_cigar(mut self) -> Self {
        self.mate_cigar = false;
        self
    }

    fn write_sam(&self, sam: &mut String) {
        for (own, mate, segment) in [(self.r1, self.r2, 0x40), (self.r2, self.r1, 0x80)] {
            let mut flags = 0x1 | segment | self.extra_flags;
            if own.is_none() {
                flags |= 0x4;
            }
            if mate.is_none() {
                flags |= 0x8;
            }
            if own.is_some_and(|(_, _, reverse)| reverse) {
                flags |= 0x10;
            }
            if mate.is_some_and(|(_, _, reverse)| reverse) {
                flags |= 0x20;
            }
            let place = own.or(mate).map_or(("*", 0), |(contig, pos, _)| (contig, pos));
            let mate_place = mate.or(own).map_or(("*", 0), |(contig, pos, _)| (contig, pos));
            let (mapq, cigar) = if own.is_some() { (self.mapq, "10M") } else { (0, "*") };
            write!(
                sam,
                "{}\t{flags}\t{}\t{}\t{mapq}\t{cigar}\t{}\t{}\t0\tACGTACGTAC\t{}\tRG:Z:{}\tRX:Z:{}",
                self.name,
                place.0,
                place.1,
                mate_place.0,
                mate_place.1,
                self.quality.to_string().repeat(10),
                self.read_group,
                self.umi,
            )
            .unwrap();
            if let Some(cell) = self.cell {
                write!(sam, "\tCB:Z:{cell}").unwrap();
            }
            if mate.is_some() {
                if self.mate_cigar {
                    sam.push_str("\tMC:Z:10M");
                }
                write!(sam, "\tMQ:i:{}", self.mapq).unwrap();
            }
            sam.push('\n');
        }
    }
}

fn run(args: &[&str]) -> String {
    run_args_and_capture_logs(&args.iter().map(OsStr::new).collect::<Vec<_>>())
}

fn sort(input: &Path, output: &Path, args: &[&str]) -> String {
    let mut all = vec!["sort", "-i", p(input), "-o", p(output), "--order", "template-coordinate"];
    all.extend(args);
    run(&all)
}

fn sort_sam(dir: &TempDir, sam: &str) -> PathBuf {
    let (unsorted, sorted) = (dir.path().join("in.sam"), dir.path().join("in.bam"));
    fs::write(&unsorted, sam).unwrap();
    sort(&unsorted, &sorted, &[]);
    sorted
}

fn sam_of(pairs: &[Pair]) -> String {
    let mut sam = HEADER.to_string();
    for pair in pairs {
        pair.write_sam(&mut sam);
    }
    sam
}

fn run_command(command: &str, input: &Path, output: &Path, args: &[&str]) {
    run_and_capture_logs(command, input, output, &[args, &["--verify"]].concat());
}

fn run_dedup(input: &Path, output: &Path, args: &[&str]) {
    run_command("dedup", input, output, args);
}

fn dedup(sam: &str, args: &[&str]) -> Vec<RecordBuf> {
    output_of("dedup", sam, args)
}

fn group(sam: &str, args: &[&str]) -> Vec<RecordBuf> {
    output_of("group", sam, args)
}

/// Sorts `sam`, runs `command` on it with one and four threads, requires the same bytes from
/// both, and returns the records written.
fn output_of(command: &str, sam: &str, args: &[&str]) -> Vec<RecordBuf> {
    let dir = TempDir::new().unwrap();
    let input = sort_sam(&dir, sam);
    let outputs = ["1", "4"].map(|threads| {
        let output = dir.path().join(format!("out.{threads}.bam"));
        run_command(command, &input, &output, &[args, &["--threads", threads]].concat());
        output
    });
    assert_eq!(
        decompressed_records_without_pg(&outputs[0]),
        decompressed_records_without_pg(&outputs[1]),
        "the output must not depend on --threads"
    );
    assert_bam_sorted(&outputs[0], "template-coordinate", Some("mi"));
    read_bam_output(&outputs[0]).1
}

fn name_of(record: &RecordBuf) -> String {
    record.name().expect("named record").to_string()
}

fn mi_of(record: &RecordBuf) -> Option<String> {
    match record.data().get(&SamTag::MI.to_noodles_tag())? {
        Value::String(mi) => Some(mi.to_string()),
        other => panic!("MI is not a string: {other:?}"),
    }
}

/// The names of the templates marked as duplicates, all records of which must be marked.
fn marked(records: &[RecordBuf]) -> BTreeSet<String> {
    let names = |duplicate: bool| -> BTreeSet<String> {
        records.iter().filter(|r| r.flags().is_duplicate() == duplicate).map(name_of).collect()
    };
    assert!(names(true).is_disjoint(&names(false)), "a template is marked in part");
    names(true)
}

/// Whether the template called `half` is marked, which must be the only one marked.
fn is_marked(pairs: &[Pair], args: &[&str]) -> bool {
    let records = dedup(&sam_of(pairs), args);
    assert!(records.iter().all(|r| mi_of(r).is_some()));
    let marked = marked(&records);
    assert!(marked.iter().all(|name| name == "half"), "{marked:?}");
    !marked.is_empty()
}

#[rstest]
fn copies_of_either_end_are_marked(
    #[values(PAIRED, IDENTITY, EDIT, ADJACENCY, NO_UMI)] args: &[&str],
) {
    let full = Pair::full("AAAA-CCCC");
    assert!(is_marked(&[full.clone(), Pair::r1_only(LOWER, "AAAA-GTGT")], args), "lower end");
    assert!(is_marked(&[full, Pair::r2_only(UPPER, "GTGT-CCCC")], args), "upper end");
}

/// Five full pairs and three half-mapped templates at each end, whose failed reads have random
/// UMIs, leave one template unmarked.
#[rstest]
fn a_deep_family_keeps_one_full_pair(
    #[values(PAIRED, IDENTITY, EDIT, ADJACENCY, NO_UMI)] args: &[&str],
) {
    let mut pairs: Vec<Pair> =
        (0..5).map(|i| Pair::full("ACGTACGTA-TTGCAAGCT").named(&format!("full{i}"))).collect();
    for (i, umi) in ["TACGGACCA", "CTGTTTGTT", "AGATCTCGC"].into_iter().enumerate() {
        let lower = Pair::r1_only(LOWER, &format!("ACGTACGTA-{umi}"));
        let upper = Pair::r2_only(UPPER, &format!("{umi}-TTGCAAGCT"));
        pairs.extend([lower.named(&format!("lower{i}")), upper.named(&format!("upper{i}"))]);
    }
    let records = dedup(&sam_of(&pairs), args);
    assert_eq!(records.len(), 22);
    let unmarked: BTreeSet<String> =
        records.iter().filter(|r| !r.flags().is_duplicate()).map(name_of).collect();
    assert_eq!(unmarked.len(), 1, "{unmarked:?}");
    assert!(unmarked.iter().all(|name| name.starts_with("full")));
}

/// A pair is never marked because of a half-mapped template, whatever their scores.
#[rstest]
fn a_pair_is_kept_over_a_better_half_mapped_template(
    #[values(PAIRED, IDENTITY, NO_UMI)] args: &[&str],
) {
    let pairs = [Pair::full("AAAA-CCCC").quality('#'), Pair::r1_only(LOWER, "AAAA-GTGT")];
    assert!(is_marked(&pairs, args));
}

#[rstest]
#[case::other_umi(Pair::full("AAAA-CCCC"), Pair::r1_only(LOWER, "GGGG-CCCC"))]
#[case::umi_of_the_other_end(Pair::full("AAAA-CCCC"), Pair::r1_only(LOWER, "CCCC-AAAA"))]
#[case::no_umi_on_the_mapped_end(Pair::full("-CCCC"), Pair::r1_only(LOWER, "-GTGT"))]
#[case::other_position(Pair::full("AAAA-CCCC"), Pair::r1_only(("chr1", 502, false), "AAAA-GTGT"))]
#[case::other_strand(Pair::full("AAAA-CCCC"), Pair::r1_only(("chr1", 492, true), "AAAA-GTGT"))]
#[case::other_contig(Pair::full("AAAA-CCCC"), Pair::r1_only(("chr2", 501, false), "AAAA-GTGT"))]
#[case::other_library(
    Pair::full("AAAA-CCCC"),
    Pair::r1_only(LOWER, "AAAA-GTGT").read_group("B")
)]
#[case::other_cell(
    Pair::full("AAAA-CCCC").cell("ACGT"),
    Pair::r1_only(LOWER, "AAAA-GTGT").cell("TTTT")
)]
#[case::no_cell(Pair::full("AAAA-CCCC").cell("ACGT"), Pair::r1_only(LOWER, "AAAA-GTGT"))]
fn a_half_mapped_template_with_a_different_end_is_not_marked(
    #[case] full: Pair,
    #[case] half: Pair,
    #[values(PAIRED, IDENTITY, EDIT, ADJACENCY)] args: &[&str],
) {
    assert!(!is_marked(&[full, half], args));
}

#[rstest]
#[case::library(Pair::read_group, "A", "B")]
#[case::cell(Pair::cell, "AAAA", "ACGT")]
fn copies_are_found_within_a_library_and_cell(
    #[case] set: fn(Pair, &'static str) -> Pair,
    #[case] other: &'static str,
    #[case] own: &'static str,
    #[values(PAIRED, NO_UMI)] args: &[&str],
) {
    let pairs = [
        set(Pair::full("AAAA-CCCC").named("elsewhere"), other),
        set(Pair::full("AAAA-CCCC"), own),
        set(Pair::r1_only(LOWER, "AAAA-GTGT"), own),
    ];
    assert!(is_marked(&pairs, args));
}

#[rstest]
#[case::low_mapping_quality(Pair::full("AAAA-CCCC").mapq(0))]
#[case::qc_fail(Pair::full("AAAA-CCCC").extra_flags(0x200))]
#[case::n_in_umi(Pair::full("AAAA-CCNC"))]
fn a_filtered_pair_marks_nothing(#[case] filtered: Pair) {
    let args = [PAIRED, &["--min-map-q", "1"]].concat();
    assert!(!is_marked(&[filtered, Pair::r1_only(LOWER, "AAAA-GTGT")], &args));
}

#[rstest]
#[case::r1_end(Pair::r1_only(LOWER, "AAAA-GTGT"))]
#[case::r2_end_on_a_later_contig(Pair::r2_only(("chr3", 801, true), "GTGT-CCCC"))]
fn both_ends_of_an_inter_contig_pair_are_known(#[case] half: Pair) {
    let full = Pair::new("full", Some(LOWER), Some(("chr3", 801, true)), "AAAA-CCCC");
    let (lower, upper) = (Some(("chr2", 501, false)), Some(("chr2", 701, true)));
    let between = Pair::new("between", lower, upper, "AAAA-CCCC");
    assert!(is_marked(&[full, between, half], PAIRED));
}

/// A half-mapped R2 is where the pair has its R1, so they are from opposite strands. With UMIs,
/// strands are deduplicated separately whatever the strategy. Without, like Picard, they are not.
#[rstest]
#[case::paired(PAIRED, false)]
#[case::identity(IDENTITY, false)]
#[case::edit(EDIT, false)]
#[case::adjacency(ADJACENCY, false)]
#[case::no_umi(NO_UMI, true)]
fn opposite_strands_are_separate_only_with_umis(
    #[case] args: &[&str],
    #[case] strands_together: bool,
) {
    let other_strand = [Pair::full("AAAA-CCCC"), Pair::r2_only(LOWER, "GTGT-AAAA")];
    assert_eq!(is_marked(&other_strand, args), strands_together);

    let same_strand = [Pair::flipped("CCCC-AAAA"), Pair::r2_only(LOWER, "GTGT-AAAA")];
    assert!(is_marked(&same_strand, args));

    let both_mapped = [Pair::full("AAAA-CCCC"), Pair::flipped("CCCC-AAAA").named("half")];
    assert_eq!(marked(&dedup(&sam_of(&both_mapped), args)).len(), usize::from(strands_together));
}

#[rstest]
#[case::exact("AAAACC", 0, true)]
#[case::one_mismatch_not_allowed("AAAACT", 0, false)]
#[case::one_mismatch_allowed("AAAACT", 1, true)]
#[case::two_mismatches_not_allowed("AAAATT", 1, false)]
#[case::two_mismatches_allowed("AAAATT", 2, true)]
#[case::longer("AAAACCG", 1, false)]
#[case::lower_case("aaaacc", 0, true)]
fn umis_are_compared_within_the_edits(
    #[case] umi: &str,
    #[case] edits: u32,
    #[case] expected: bool,
    #[values("paired", "edit", "adjacency")] strategy: &str,
) {
    let edits = edits.to_string();
    let half = Pair::r1_only(LOWER, &format!("{umi}-GTGTGT"));
    let full = Pair::full("AAAACC-GGGGGG");
    assert_eq!(is_marked(&[full, half], &["--strategy", strategy, "--edits", &edits]), expected);
}

#[rstest]
#[case::same("AAAACCCC", true)]
#[case::lower_case("aaaacccc", true)]
#[case::one_mismatch("AAAACCCT", false)]
#[case::shorter("AAAACCC", false)]
#[case::long("ACGTTACGTTACGTTACGTT", true)]
fn a_single_umi_is_compared_whole(#[case] umi: &str, #[case] expected: bool) {
    let full = if umi.len() > 8 { umi } else { "AAAACCCC" };
    assert_eq!(is_marked(&[Pair::full(full), Pair::r1_only(LOWER, umi)], IDENTITY), expected);
}

#[test]
fn long_umis_are_compared_base_by_base() {
    let pairs = [Pair::full("ACGTTACGTTACGTTACGTT"), Pair::r1_only(LOWER, "ACGTTACGTTACGTTACGTA")];
    assert!(!is_marked(&pairs, IDENTITY));
    assert!(is_marked(&pairs, ADJACENCY));
}

/// With `--min-umi-length`, single UMIs are compared over the length of the shorter one.
#[rstest]
#[case::shorter("ACGTAC", true)]
#[case::longer("ACGTACGTAC", true)]
#[case::shorter_and_different("ACGTCC", false)]
fn single_umis_are_truncated_under_min_umi_length(
    #[case] umi: &str,
    #[case] expected: bool,
    #[values("identity", "adjacency")] strategy: &str,
) {
    let args = ["--strategy", strategy, "--edits", "0", "--min-umi-length", "6"];
    assert_eq!(is_marked(&[Pair::full("ACGTACGT"), Pair::r1_only(LOWER, umi)], &args), expected);
}

/// Matching by way of another half-mapped template would join UMIs two edits apart.
#[test]
fn matches_are_not_chained() {
    let args = ["--strategy", "edit", "--edits", "1"];
    let pairs = [
        Pair::full("AAAAAA-CCCCCC"),
        Pair::r1_only(LOWER, "AAAAAT-GTGTGT").named("near").quality('5'),
        Pair::r1_only(LOWER, "AAAATT-GTGTGT").named("far"),
    ];
    assert_eq!(marked(&dedup(&sam_of(&pairs), &args)), BTreeSet::from(["near".to_string()]));
}

/// Families of half-mapped templates that match no pair keep their best template.
#[rstest]
fn half_mapped_templates_without_a_pair_are_deduplicated_with_each_other(
    #[values(PAIRED, IDENTITY, EDIT, ADJACENCY, NO_UMI)] args: &[&str],
    #[values(true, false)] elsewhere: bool,
) {
    let mut pairs = vec![
        Pair::r1_only(LOWER, "GGGG-TTTT").named("worse").quality('5'),
        Pair::r1_only(LOWER, "GGGG-TTTT").named("better"),
    ];
    if elsewhere {
        let (r1, r2) = (Some(("chr1", 502, false)), Some(UPPER));
        pairs.push(Pair::new("full", r1, r2, "GGGG-TTTT"));
    }
    let records = dedup(&sam_of(&pairs), args);
    assert_eq!(marked(&records), BTreeSet::from(["worse".to_string()]));
    let mis: BTreeSet<_> = records.iter().filter(|r| name_of(r) != "full").map(mi_of).collect();
    assert_eq!(mis.len(), 1, "one family");
}

/// A family that is a copy of a pair is marked whole.
#[test]
fn a_family_of_copies_is_marked_whole() {
    let pairs = [
        Pair::full("AAAA-CCCC"),
        Pair::r1_only(LOWER, "AAAA-GTGT").named("worse").quality('5'),
        Pair::r1_only(LOWER, "AAAA-GTGT").named("better"),
    ];
    let both = BTreeSet::from(["better", "worse"].map(String::from));
    assert_eq!(marked(&dedup(&sam_of(&pairs), IDENTITY)), both);
}

#[test]
fn unpaired_and_unmapped_templates_are_not_marked() {
    let mut sam =
        sam_of(&[Pair::full("AAAA-CCCC"), Pair::new("unmapped", None, None, "AAAA-CCCC")]);
    sam.push_str(
        "single\t0\tchr1\t501\t60\t10M\t*\t0\t0\tACGTACGTAC\tIIIIIIIIII\tRG:Z:A\tRX:Z:AAAA-CCCC\n",
    );
    let records = dedup(&sam, &["--strategy", "identity", "--include-unmapped"]);
    assert_eq!(records.len(), 5);
    assert!(marked(&records).is_empty());
}

#[test]
fn a_pair_without_mc_is_still_an_error() {
    let dir = TempDir::new().unwrap();
    let no_mc = Pair::full("AAAA-CCCC").named("half").without_mate_cigar();
    let input = sort_sam(&dir, &sam_of(&[Pair::full("AAAA-CCCC"), no_mc]));
    let output = dir.path().join("out.bam");
    let result = std::process::Command::new(env!("CARGO_BIN_EXE_fgumi"))
        .args(["dedup", "-i", p(&input), "-o", p(&output), "--strategy", "identity"])
        .output()
        .unwrap();
    assert!(!result.status.success());
    assert!(String::from_utf8_lossy(&result.stderr).contains("requires MC tags"));
}

/// A mapped read whose unmapped mate is absent from the input is still half of a pair.
#[test]
fn a_half_mapped_template_without_its_unmapped_record_is_marked() {
    let mut sam = sam_of(&[Pair::full("AAAA-CCCC")]);
    sam.push_str(
        "half\t73\tchr1\t501\t60\t10M\t=\t501\t0\tACGTACGTAC\tIIIIIIIIII\tRG:Z:A\tRX:Z:AAAA-GTGT\n",
    );
    assert_eq!(marked(&dedup(&sam, PAIRED)).len(), 1);
}

/// The end is the unclipped 5' position, so a soft-clipped copy still matches.
#[test]
fn clipped_ends_match_by_unclipped_position() {
    let mut sam = sam_of(&[Pair::full("AAAA-CCCC")]);
    let tail = "ACGTACGTAC\tIIIIIIIIII\tRG:Z:A\tRX:Z:AAAA-GTGT";
    writeln!(sam, "half\t73\tchr1\t504\t60\t3S7M\t=\t504\t0\t{tail}").unwrap();
    writeln!(sam, "half\t133\tchr1\t504\t0\t*\t=\t504\t0\t{tail}\tMC:Z:3S7M").unwrap();
    assert_eq!(marked(&dedup(&sam, PAIRED)).len(), 1);
}

/// Secondary and supplementary records of a copy are marked and removed with it.
#[test]
fn all_records_of_a_copy_go_with_it() {
    let mut sam = sam_of(&[Pair::full("AAAA-CCCC")]);
    let tail = "ACGTACGTAC\tIIIIIIIIII\tRG:Z:A\tRX:Z:AAAA-GTGT";
    let tc = "tc:B:i,0,501,0,2147483647,2147483647,0";
    writeln!(sam, "half\t73\tchr1\t501\t60\t10M\t=\t501\t0\t{tail}").unwrap();
    writeln!(sam, "half\t133\tchr1\t501\t0\t*\t=\t501\t0\t{tail}\tMC:Z:10M\tMQ:i:60").unwrap();
    writeln!(sam, "half\t2121\tchr2\t9001\t60\t10M\tchr1\t501\t0\t{tail}\t{tc}").unwrap();
    writeln!(sam, "half\t329\tchr3\t9001\t0\t10M\tchr1\t501\t0\t{tail}\t{tc}").unwrap();

    let records = dedup(&sam, PAIRED);
    assert_eq!(records.len(), 6);
    assert_eq!(marked(&records), BTreeSet::from(["half".to_string()]));

    let removed = dedup(&sam, &["--strategy", "paired", "--remove-duplicates"]);
    assert_eq!(removed.len(), 2);
    assert!(removed.iter().all(|r| name_of(r) == "full"));
}

/// Pairs and half-mapped templates at many positions, so that they fall in different batches.
#[test]
fn copies_in_other_batches_are_marked_at_any_thread_count() {
    let mut pairs = Vec::new();
    for i in 0..600u32 {
        let contig = ["chr1", "chr2", "chr3"][(i % 3) as usize];
        let lower = Some((contig, 1_000 + 37 * i, false));
        let upper = Some((contig, 1_300 + 37 * i, true));
        pairs.push(Pair::new(&format!("full{i}"), lower, upper, "ACGTAC-TGCATG"));
        pairs.push(match i % 3 {
            0 => Pair::new(&format!("copy{i}"), lower, None, "ACGTAC-GTGTGT"),
            1 => Pair::new(&format!("copy{i}"), None, upper, "GTGTGT-TGCATG"),
            _ => Pair::new(&format!("other{i}"), lower, None, "GTGTGT-TGCATG"),
        });
    }
    let dir = TempDir::new().unwrap();
    let input = sort_sam(&dir, &sam_of(&pairs));
    let outputs = ["1", "2", "8"].map(|threads| {
        let output = dir.path().join(format!("{threads}.bam"));
        run_dedup(&input, &output, &["--strategy", "paired", "--threads", threads]);
        output
    });
    let bytes = outputs.each_ref().map(|output| decompressed_records_without_pg(output));
    assert_eq!(bytes[0], bytes[1]);
    assert_eq!(bytes[0], bytes[2]);

    let marked = marked(&read_bam_output(&outputs[2]).1);
    assert_eq!(marked.len(), 400);
    assert!(marked.iter().all(|name| name.starts_with("copy")));
}

/// Reads the aggregate metrics row as column name to value.
fn metrics(path: &Path) -> BTreeMap<String, String> {
    let text = fs::read_to_string(path).unwrap();
    let header: Vec<&str> = text.lines().next().unwrap().split('\t').collect();
    let row = text.lines().find(|line| line.contains("All Reads")).expect("aggregate row");
    header.iter().zip(row.split('\t')).map(|(k, v)| ((*k).to_string(), v.to_string())).collect()
}

#[rstest]
fn copies_are_counted_as_duplicate_orphans(
    #[values(false, true)] remove: bool,
    #[values("1", "4")] threads: &str,
) {
    let dir = TempDir::new().unwrap();
    let pairs = [
        Pair::full("AAAA-CCCC"),
        Pair::r1_only(LOWER, "AAAA-GTGT").named("copy"),
        Pair::r1_only(LOWER, "GGGG-GTGT").named("other"),
    ];
    let input = sort_sam(&dir, &sam_of(&pairs));
    let file = |name: &str| dir.path().join(name);
    let (output, counts, families) = (file("out.bam"), file("m.txt"), file("f.txt"));
    let mut args = vec!["--strategy", "identity", "--threads", threads];
    args.extend(["--metrics", p(&counts), "--family-size-histogram", p(&families)]);
    if remove {
        args.push("--remove-duplicates");
    }
    run_dedup(&input, &output, &args);

    let names: BTreeSet<String> = read_bam_output(&output).1.iter().map(name_of).collect();
    assert_eq!(names.contains("copy"), !remove);
    assert!(names.contains("full") && names.contains("other"));

    let counts = metrics(&counts);
    for (column, expected) in [
        ("filtered_templates", "0"),
        ("total_templates", "3"),
        ("unique_templates", "2"),
        ("duplicate_templates", "1"),
        ("total_reads", "6"),
        ("unique_reads", "4"),
        ("duplicate_reads", "2"),
        ("mapped_pairs", "1"),
        ("duplicate_pairs", "0"),
        ("mapped_orphans", "2"),
        ("duplicate_orphans", "1"),
        ("unmapped_pairs", "2"),
    ] {
        assert_eq!(counts[column], expected, "{column}");
    }
    let families = fs::read_to_string(&families).unwrap();
    assert!(families.lines().nth(1).unwrap().starts_with("1\t1\t"), "{families}");
    assert!(families.lines().nth(2).unwrap().starts_with("2\t1\t"), "{families}");
}

/// The molecule identifier of the template called `name`, which all its records must have.
fn mi_of_template(records: &[RecordBuf], name: &str) -> String {
    let mis: BTreeSet<_> =
        records.iter().filter(|r| name_of(r) == name).map(|r| mi_of(r).expect("an MI")).collect();
    assert_eq!(mis.len(), 1, "{name}: {mis:?}");
    mis.into_iter().next().unwrap()
}

/// The templates of each molecule, both strands together, and whether its records are contiguous.
fn molecules(records: &[RecordBuf]) -> BTreeSet<BTreeSet<String>> {
    let base = |r: &RecordBuf| mi_of(r).expect("an MI").split('/').next().unwrap().to_string();
    let mut molecules: BTreeMap<String, BTreeSet<String>> = BTreeMap::new();
    let mut last = None;
    for record in records {
        let base = base(record);
        if last.as_ref() != Some(&base) {
            assert!(!molecules.contains_key(&base), "molecule {base} is not contiguous");
        }
        molecules.entry(base.clone()).or_default().insert(name_of(record));
        last = Some(base);
    }
    molecules.into_values().collect()
}

fn names(names: &[&str]) -> BTreeSet<String> {
    names.iter().map(ToString::to_string).collect()
}

fn deep_family() -> Vec<Pair> {
    let mut pairs: Vec<Pair> =
        (0..5).map(|i| Pair::full("ACGTACGTA-TTGCAAGCT").named(&format!("full{i}"))).collect();
    for (i, umi) in ["TACGGACCA", "CTGTTTGTT", "AGATCTCGC"].into_iter().enumerate() {
        let lower = Pair::r1_only(LOWER, &format!("ACGTACGTA-{umi}"));
        let upper = Pair::r2_only(UPPER, &format!("{umi}-TTGCAAGCT"));
        pairs.extend([lower.named(&format!("lower{i}")), upper.named(&format!("upper{i}"))]);
    }
    pairs
}

#[rstest]
fn group_gives_copies_their_molecule(
    #[values(PAIRED, IDENTITY, EDIT, ADJACENCY, NO_UMI)] args: &[&str],
) {
    let full = Pair::full("AAAA-CCCC");
    for half in [Pair::r1_only(LOWER, "AAAA-GTGT"), Pair::r2_only(UPPER, "GTGT-CCCC")] {
        let records = group(&sam_of(&[full.clone(), half]), args);
        assert_eq!(records.len(), 4, "both records of the copy are kept");
        assert_eq!(mi_of_template(&records, "half"), mi_of_template(&records, "full"));
    }
    assert_eq!(molecules(&group(&sam_of(&deep_family()), args)).len(), 1);
}

#[rstest]
#[case::other_umi(Pair::r1_only(LOWER, "GGGG-CCCC"))]
#[case::other_position(Pair::r1_only(("chr1", 502, false), "AAAA-GTGT"))]
#[case::other_strand(Pair::r1_only(("chr1", 492, true), "AAAA-GTGT"))]
#[case::other_contig(Pair::r1_only(("chr2", 501, false), "AAAA-GTGT"))]
#[case::other_library(Pair::r1_only(LOWER, "AAAA-GTGT").read_group("B"))]
fn group_keeps_other_half_mapped_templates_apart(
    #[case] half: Pair,
    #[values(PAIRED, IDENTITY, EDIT, ADJACENCY)] args: &[&str],
) {
    let records = group(&sam_of(&[Pair::full("AAAA-CCCC"), half]), args);
    assert_eq!(molecules(&records), BTreeSet::from([names(&["full"]), names(&["half"])]));
}

#[test]
fn group_keeps_cells_apart() {
    let full = Pair::full("AAAA-CCCC");
    let half = Pair::r1_only(LOWER, "AAAA-GTGT").cell("TTTT");
    let pairs = [full.clone().cell("ACGT"), half.clone()];
    assert_eq!(molecules(&group(&sam_of(&pairs), PAIRED)).len(), 2);
    assert_eq!(molecules(&group(&sam_of(&[full.cell("TTTT"), half]), PAIRED)).len(), 1);
}

/// Where two pairs with different second ends match, `sort` cannot choose, `group` keeps the
/// template apart and `dedup` marks it.
#[rstest]
fn an_ambiguous_end_stays_on_its_own(#[values(PAIRED, IDENTITY)] args: &[&str]) {
    let longer = Pair::new("longer", Some(LOWER), Some(("chr1", 801, true)), "AAAA-CCCC");
    let sam = sam_of(&[Pair::full("AAAA-CCCC"), longer, Pair::r1_only(LOWER, "AAAA-GTGT")]);
    let records = group(&sam, args);
    let alone = BTreeSet::from([names(&["full"]), names(&["longer"]), names(&["half"])]);
    assert_eq!(molecules(&records), alone);
    assert!(records.iter().all(|r| r.data().get(&SamTag::TC.to_noodles_tag()).is_none()));
    assert_eq!(marked(&dedup(&sam, args)), names(&["half"]));
}

/// As above, but the template is from the other strand, which only `--no-umi` ignores.
#[rstest]
#[case::paired(PAIRED, false)]
#[case::adjacency(ADJACENCY, false)]
#[case::no_umi(NO_UMI, true)]
fn an_ambiguous_end_of_the_other_strand(#[case] args: &[&str], #[case] expected: bool) {
    let longer = Pair::new("longer", Some(LOWER), Some(("chr1", 801, true)), "AAAA-CCCC");
    let pairs = [Pair::full("AAAA-CCCC"), longer, Pair::r2_only(LOWER, "GTGT-AAAA")];
    assert_eq!(is_marked(&pairs, args), expected);
}

/// Two molecules of one position match equally, so `group` joins the template to neither.
#[rstest]
fn two_molecules_are_not_bridged(
    #[values(PAIRED, ADJACENCY, PARALLEL_PAIRED, PARALLEL_ADJACENCY)] args: &[&str],
) {
    let pairs = [
        Pair::full("AAAA-CCCC"),
        Pair::full("GGTT-CCCC").named("other"),
        Pair::r2_only(UPPER, "ACAC-CCCC"),
    ];
    let sam = sam_of(&pairs);
    let alone = BTreeSet::from([names(&["full"]), names(&["other"]), names(&["half"])]);
    assert_eq!(molecules(&group(&sam, args)), alone);
    assert_eq!(marked(&dedup(&sam, &args[..2])), names(&["half"]));
}

/// Copies neither seed molecules nor count towards the assignment of fully mapped pairs.
#[rstest]
fn copies_do_not_change_the_other_molecules(
    #[values("paired", "edit", "adjacency", "identity")] strategy: &str,
) {
    let mut pairs = Vec::new();
    for (i, (umi, count)) in
        [("AAAAAA", 3), ("AAAAAT", 2), ("AAAATT", 1), ("AATTTT", 2)].into_iter().enumerate()
    {
        let of = |j: i32| Pair::full(&format!("{umi}-CCCCCC")).named(&format!("full{i}.{j}"));
        pairs.extend((0..count).map(of));
    }
    let args = ["--strategy", strategy, "--edits", "1"];
    let without = molecules(&group(&sam_of(&pairs), &args));
    for (i, umi) in ["AAAATT"; 9].into_iter().enumerate() {
        pairs.push(Pair::r1_only(LOWER, &format!("{umi}-GTGTGT")).named(&format!("half{i}")));
    }
    let with: BTreeSet<BTreeSet<String>> = molecules(&group(&sam_of(&pairs), &args))
        .into_iter()
        .map(|names| names.into_iter().filter(|name| name.starts_with("full")).collect())
        .filter(|names: &BTreeSet<String>| !names.is_empty())
        .collect();
    assert_eq!(with, without);
}

/// A filtered pair is in no molecule, so the template beside it is grouped on its own.
#[test]
fn a_copy_of_a_filtered_pair_is_kept() {
    let sam = sam_of(&[Pair::full("AAAA-CCCC").mapq(0), Pair::r1_only(LOWER, "AAAA-GTGT")]);
    assert_eq!(molecules(&group(&sam, PAIRED)), BTreeSet::from([names(&["half"])]));
}

/// Families with no fully mapped pair, as at a breakpoint, are grouped as before.
#[rstest]
fn half_mapped_families_are_grouped_alone(
    #[values(PAIRED, IDENTITY, EDIT, ADJACENCY, NO_UMI)] args: &[&str],
) {
    let pairs: Vec<Pair> = (0..3)
        .map(|i| Pair::r1_only(LOWER, "AAAA-CCCC").named(&format!("half{i}")))
        .chain([Pair::new("full", Some(("chr1", 502, false)), Some(UPPER), "AAAA-CCCC")])
        .collect();
    let records = group(&sam_of(&pairs), args);
    let halves = names(&["half0", "half1", "half2"]);
    assert_eq!(molecules(&records), BTreeSet::from([halves, names(&["full"])]));
}

/// Pairs that share a start but not an end are told apart by UMI, and their copies do not
/// interleave with one another's groups.
#[rstest]
fn copies_do_not_split_the_groups_around_them(#[values(PAIRED, ADJACENCY)] args: &[&str]) {
    let mut pairs = Vec::new();
    let mut expected = BTreeSet::new();
    for (end, umi) in ["ACGTAC", "CGTACA", "GTACAC", "TACACG", "AAGGTT"].into_iter().enumerate() {
        let upper = ("chr1", 701 + u32::try_from(end).unwrap(), true);
        let mut molecule = BTreeSet::new();
        for i in 0..4 {
            let (full, copy) = (format!("full{end}.{i}"), format!("copy{end}.{i}"));
            let umi = |other: &str| format!("{umi}-{other}");
            pairs.push(Pair::new(&full, Some(LOWER), Some(upper), &umi("CCCCCC")));
            pairs.push(Pair::r1_only(LOWER, &umi("GTGTGT")).named(&copy));
            molecule.extend([full, copy]);
        }
        expected.insert(molecule);
    }
    assert_eq!(molecules(&group(&sam_of(&pairs), args)), expected);
}

/// A pair in another library or cell does not make an end ambiguous.
#[rstest]
#[case::library(Pair::read_group, "A", "B")]
#[case::cell(Pair::cell, "AAAA", "ACGT")]
fn copies_are_placed_within_a_library_and_cell(
    #[case] set: fn(Pair, &'static str) -> Pair,
    #[case] other: &'static str,
    #[case] own: &'static str,
) {
    let elsewhere = Pair::new("elsewhere", Some(LOWER), Some(("chr1", 801, true)), "AAAA-CCCC");
    let pairs = [
        set(elsewhere, other),
        set(Pair::full("AAAA-CCCC"), own),
        set(Pair::r1_only(LOWER, "AAAA-GTGT"), own),
    ];
    assert_eq!(
        molecules(&group(&sam_of(&pairs), PAIRED)),
        BTreeSet::from([names(&["elsewhere"]), names(&["full", "half"])])
    );
}

/// The strand suffix follows the copy's own orientation, so it can be the other strand's only read.
#[test]
fn an_opposite_strand_copy_is_the_other_strand_of_the_molecule() {
    for (full, half) in [
        (Pair::full("AAAA-CCCC"), Pair::r2_only(LOWER, "GTGT-AAAA")),
        (Pair::full("AAAA-CCCC"), Pair::r1_only(UPPER, "CCCC-GTGT")),
        (Pair::flipped("CCCC-AAAA"), Pair::r1_only(LOWER, "AAAA-GTGT")),
        (Pair::full("CCCC-AAAA"), Pair::r2_only(LOWER, "GTGT-CCCC")),
    ] {
        let alone = mi_of_template(&group(&sam_of(std::slice::from_ref(&full)), PAIRED), "full");
        let both = group(&sam_of(&[full.clone(), half.clone()]), PAIRED);
        assert_eq!(mi_of_template(&both, "full"), alone);
        let other = if alone.ends_with("/A") { "0/B" } else { "0/A" };
        assert_eq!(mi_of_template(&both, "half"), other);

        // What the copy would have got if its other read had mapped.
        let flipped = if full.r1 == Some(LOWER) {
            Pair::flipped(&full.umi.split('-').rev().collect::<Vec<_>>().join("-"))
        } else {
            Pair::full(&full.umi.split('-').rev().collect::<Vec<_>>().join("-"))
        };
        let mapped = group(&sam_of(&[full, flipped.named("half")]), PAIRED);
        assert_eq!(mi_of_template(&mapped, "half"), other);
    }
    for args in [IDENTITY, EDIT, ADJACENCY, NO_UMI] {
        let pairs = [Pair::full("AAAA-CCCC"), Pair::r2_only(LOWER, "GTGT-AAAA")];
        assert_eq!(molecules(&group(&sam_of(&pairs), args)).len(), 2, "{args:?}");
    }
}

fn int_tag(record: &RecordBuf, tag: SamTag) -> Option<i64> {
    record.data().get(&tag.to_noodles_tag()).and_then(Value::as_int)
}

/// Consensus depth counts the mapped reads of copies and nothing else, and the consensus UMI is
/// that of the fully mapped pairs.
#[cfg(feature = "consensus")]
#[rstest]
#[case::deep(deep_family(), [8, 8], "ACGTACGTA-TTGCAAGCT")]
#[case::mini(vec![Pair::full("AAAA-CCCC"), Pair::r1_only(LOWER, "AAAA-GTGT")], [2, 1], "AAAA-CCCC")]
fn consensus_uses_the_mapped_reads_of_copies(
    #[case] pairs: Vec<Pair>,
    #[case] depths: [i64; 2],
    #[case] umi: &str,
    #[values("true", "false")] allow_unmapped: &str,
) {
    let dir = TempDir::new().unwrap();
    let file = |name: &str| dir.path().join(name);
    let (grouped, called) = (file("grouped.bam"), file("called.bam"));
    run_command("group", &sort_sam(&dir, &sam_of(&pairs)), &grouped, IDENTITY_FOR_SIMPLEX);
    let unmapped = format!("--allow-unmapped={allow_unmapped}");
    run(&["simplex", "-i", p(&grouped), "-o", p(&called), "--min-reads", "1", &unmapped]);
    let records = read_bam_output(&called).1;
    assert_eq!(records.len(), 2, "one consensus pair");
    for (record, depth) in records.iter().zip(depths) {
        assert_eq!(int_tag(record, SamTag::CD), Some(depth));
        match record.data().get(&SamTag::RX.to_noodles_tag()) {
            Some(Value::String(rx)) => assert_eq!(rx.to_string(), umi),
            other => panic!("no UMI: {other:?}"),
        }
    }
}

const IDENTITY_FOR_SIMPLEX: &[&str] = &["--strategy", "identity"];

/// Half-mapped templates of the strand that has no fully mapped pair give their end duplex
/// support, and their unmapped reads give the other end nothing.
#[cfg(feature = "consensus")]
#[rstest]
#[case::none(&[], "1,1,0", Some([[2, 0], [2, 0]]))]
#[case::lower_end(&[Pair::r2_only(LOWER, "GTGT-AAAA")], "1,1,0", Some([[2, 1], [2, 0]]))]
#[case::upper_end(&[Pair::r1_only(UPPER, "CCCC-GTGT")], "1,1,0", Some([[2, 0], [2, 1]]))]
#[case::both_ends(
    &[Pair::r2_only(LOWER, "GTGT-AAAA"), Pair::r1_only(UPPER, "CCCC-GTGT").named("half2")],
    "1,1,0",
    Some([[2, 1], [2, 1]])
)]
#[case::none_where_two_strands_are_needed(&[], "2,1,1", None)]
#[case::one_end_where_two_strands_are_needed(
    &[Pair::r2_only(LOWER, "GTGT-AAAA")],
    "2,1,1",
    None
)]
#[case::both_ends_where_two_strands_are_needed(
    &[Pair::r2_only(LOWER, "GTGT-AAAA"), Pair::r1_only(UPPER, "CCCC-GTGT").named("half2")],
    "2,1,1",
    Some([[2, 1], [2, 1]])
)]
fn opposite_strand_copies_complete_a_duplex(
    #[case] copies: &[Pair],
    #[case] min_reads: &str,
    #[case] expected: Option<[[i64; 2]; 2]>,
    #[values("true", "false")] allow_unmapped: &str,
) {
    let dir = TempDir::new().unwrap();
    let file = |name: &str| dir.path().join(name);
    let (grouped, called) = (file("grouped.bam"), file("called.bam"));
    let mut pairs = vec![Pair::full("AAAA-CCCC"), Pair::full("AAAA-CCCC").named("full1")];
    pairs.extend_from_slice(copies);
    run_command("group", &sort_sam(&dir, &sam_of(&pairs)), &grouped, PAIRED);
    let unmapped = format!("--allow-unmapped={allow_unmapped}");
    run(&["duplex", "-i", p(&grouped), "-o", p(&called), "--min-reads", min_reads, &unmapped]);
    let records = read_bam_output(&called).1;
    let depths: Vec<[i64; 2]> = records
        .iter()
        .map(|r| [SamTag::AD, SamTag::BD].map(|tag| int_tag(r, tag).unwrap_or(0)))
        .collect();
    assert_eq!(depths, expected.map_or(Vec::new(), |ends| ends.to_vec()));
}

/// Pairs at many positions with copies, ambiguous copies and lone half-mapped templates.
fn many_templates() -> Vec<Pair> {
    let mut pairs = Vec::new();
    for i in 0..900u32 {
        let contig = ["chr1", "chr2", "chr3"][(i % 3) as usize];
        let lower = Some((contig, 1_000 + 37 * (i / 2), false));
        let upper = Some((contig, 1_300 + 37 * i, true));
        let umi = ["ACGTAC-TGCATG", "CATCAT-GGATCC"][(i / 5 % 2) as usize];
        let (a, b) = umi.split_once('-').unwrap();
        for j in 0..=(i % 2) {
            pairs.push(Pair::new(&format!("full{i}.{j}"), lower, upper, umi));
        }
        match i % 4 {
            0 => pairs.push(Pair::new(&format!("half{i}"), lower, None, &format!("{a}-GTGTGT"))),
            1 => pairs.push(Pair::new(&format!("half{i}"), None, upper, &format!("GTGTGT-{b}"))),
            2 => pairs.push(Pair::new(&format!("half{i}"), lower, None, "GTGTGT-GTGTGT")),
            _ => {}
        }
    }
    pairs
}

/// The same bytes from SAM and BAM input, with and without spilling, at any thread count, and from
/// sorting the output again.
#[test]
fn sort_output_does_not_depend_on_how_it_is_run() {
    let dir = TempDir::new().unwrap();
    let file = |name: &str| dir.path().join(name);
    let (sam, bam) = (file("in.sam"), file("in.bam"));
    fs::write(&sam, sam_of(&many_templates())).unwrap();
    run(&["sort", "-i", p(&sam), "-o", p(&bam), "--order", "queryname"]);

    let sorted = file("sorted.bam");
    let log = sort(&bam, &sorted, &[]);
    assert!(log.contains("half-mapped templates beside a fully mapped pair"), "{log}");
    assert_bam_sorted(&sorted, "template-coordinate", None);
    let expected = decompressed_records_without_pg(&sorted);
    let output = file("out.bam");
    for input in [&sam, &bam, &sorted] {
        for threads in ["1", "4"] {
            for memory in ["1M", "768M"] {
                let args = ["--threads", threads, "--max-memory", memory];
                let log = sort(input, &output, &args);
                let records = |path: &Path| {
                    let mut names: Vec<_> = read_bam_output(path).1.iter().map(name_of).collect();
                    names.sort_unstable();
                    names
                };
                assert_eq!(records(&output), records(&sorted), "{input:?} {args:?}");
                let same = decompressed_records_without_pg(&output);
                if input == &sam {
                    assert_eq!(read_bam_output(&output).1, read_bam_output(&sorted).1);
                } else {
                    assert!(same == expected, "{input:?} {args:?}\n{log}");
                }
            }
        }
    }
}

/// Input without half-mapped templates is not touched.
#[test]
fn sort_adds_no_tag_without_half_mapped_templates() {
    let dir = TempDir::new().unwrap();
    let sorted = sort_sam(
        &dir,
        &sam_of(&[Pair::full("AAAA-CCCC"), Pair::flipped("CCCC-AAAA").named("other")]),
    );
    let records = read_bam_output(&sorted).1;
    assert!(records.iter().all(|r| r.data().get(&SamTag::TC.to_noodles_tag()).is_none()));
}

/// Sorting again without the pair puts its copy back where it was before.
#[test]
fn sorting_again_undoes_a_placement_that_no_longer_holds() {
    let dir = TempDir::new().unwrap();
    let half = Pair::r1_only(LOWER, "AAAA-GTGT");
    let sorted = sort_sam(&dir, &sam_of(&[Pair::full("AAAA-CCCC"), half.clone()]));
    let (header, records) = read_bam_output(&sorted);
    let placed: Vec<_> = records.into_iter().filter(|r| name_of(r) == "half").collect();
    assert!(placed.iter().all(|r| r.data().get(&SamTag::TC.to_noodles_tag()).is_some()));

    let (stale, again) = (dir.path().join("stale.bam"), dir.path().join("again.bam"));
    let raw: Vec<_> = placed
        .iter()
        .map(|r| fgumi_raw_bam::encode_record_buf_to_raw(r, &header).unwrap())
        .collect();
    write_bam(&stale, &header, &raw);
    sort(&stale, &again, &[]);
    let alone = sort_sam(&dir, &sam_of(&[half]));
    assert_eq!(read_bam_output(&again).1, read_bam_output(&alone).1);
}

/// A placement's secondary and supplementary records sort, and are marked, with it.
#[test]
fn all_records_of_a_copy_are_placed() {
    let mut sam = sam_of(&[
        Pair::full("AAAA-CCCC"),
        Pair::new("later", Some(("chr1", 601, false)), Some(UPPER), "AAAA-CCCC"),
    ]);
    let tail = "ACGTACGTAC\tIIIIIIIIII\tRG:Z:A\tRX:Z:AAAA-GTGT";
    let tc = "tc:B:i,0,501,0,2147483647,2147483647,0";
    writeln!(sam, "half\t73\tchr1\t501\t60\t10M\t=\t501\t0\t{tail}").unwrap();
    writeln!(sam, "half\t133\tchr1\t501\t0\t*\t=\t501\t0\t{tail}\tMC:Z:10M\tMQ:i:60").unwrap();
    writeln!(sam, "half\t2121\tchr2\t9001\t60\t10M\tchr1\t501\t0\t{tail}\t{tc}").unwrap();
    let dir = TempDir::new().unwrap();
    let records = read_bam_output(&sort_sam(&dir, &sam)).1;
    let order: Vec<String> = records.iter().map(name_of).collect();
    let halves: Vec<usize> = (0..order.len()).filter(|&i| order[i] == "half").collect();
    assert_eq!(halves.len(), 3);
    assert_eq!(halves[2] - halves[0], 2, "{order:?}");
    assert!(order[..5].contains(&"full".to_string()) && order[5] == "later", "{order:?}");
}

/// A secondary record whose flags do not say that its template is half-mapped goes with it too.
#[test]
fn records_are_held_by_their_tc_tag_too() {
    let mut sam = sam_of(&[
        Pair::full("AAAA-CCCC"),
        Pair::new("later", Some(("chr1", 601, false)), Some(UPPER), "AAAA-CCCC"),
    ]);
    let tail = "ACGTACGTAC\tIIIIIIIIII\tRG:Z:A\tRX:Z:AAAA-GTGT";
    let tc = "tc:B:i,0,501,0,2147483647,2147483647,0";
    writeln!(sam, "half\t73\tchr1\t501\t60\t10M\t=\t501\t0\t{tail}").unwrap();
    writeln!(sam, "half\t133\tchr1\t501\t0\t*\t=\t501\t0\t{tail}\tMC:Z:10M\tMQ:i:60").unwrap();
    writeln!(sam, "half\t321\tchr3\t9001\t0\t10M\tchr1\t501\t0\t{tail}\t{tc}").unwrap();
    let records = dedup(&sam, PAIRED);
    assert_eq!(records.len(), 7);
    assert_eq!(marked(&records), names(&["half"]));
}

/// Merging sorted files keeps copies beside their pairs.
#[test]
fn merged_files_stay_in_order() {
    let dir = TempDir::new().unwrap();
    let file = |name: &str| dir.path().join(name);
    let (odd, even): (Vec<Pair>, Vec<Pair>) = many_templates()
        .into_iter()
        .partition(|pair| pair.name.trim_start_matches(char::is_alphabetic).len() % 2 == 0);
    for (name, pairs) in [("odd", odd), ("even", even)] {
        let sam = file(&format!("{name}.sam"));
        fs::write(&sam, sam_of(&pairs)).unwrap();
        sort(&sam, &file(&format!("{name}.bam")), &[]);
    }
    let merged = file("merged.bam");
    run(&[
        "merge",
        "--order",
        "template-coordinate",
        "-o",
        p(&merged),
        p(&file("odd.bam")),
        p(&file("even.bam")),
    ]);
    assert_bam_sorted(&merged, "template-coordinate", None);
}

/// Input in template-coordinate order from another tool has no placements: `group` keeps the
/// template apart, and `dedup` still marks it.
#[test]
fn input_sorted_by_another_tool() {
    let dir = TempDir::new().unwrap();
    let header = HEADER.replace("SO:unsorted", "SO:unsorted\tGO:query\tSS:template-coordinate");
    let body = sam_of(&[Pair::full("AAAA-CCCC"), Pair::r1_only(LOWER, "AAAA-GTGT")]);
    let input = dir.path().join("other.sam");
    fs::write(&input, body.replace(HEADER, &header)).unwrap();
    let (grouped, deduped) = (dir.path().join("group.bam"), dir.path().join("dedup.bam"));
    run_command("group", &input, &grouped, PAIRED);
    run_command("dedup", &input, &deduped, PAIRED);
    assert_eq!(
        molecules(&read_bam_output(&grouped).1),
        BTreeSet::from([names(&["full"]), names(&["half"])])
    );
    assert_eq!(marked(&read_bam_output(&deduped).1), names(&["half"]));
}

/// A fused `runall` gives the same records as the separate commands.
#[cfg(feature = "consensus")]
#[rstest]
fn runall_matches_separate_commands(#[values("1", "4")] threads: &str) {
    let dir = TempDir::new().unwrap();
    let file = |name: &str| dir.path().join(name);
    let unsorted = file("in.sam");
    fs::write(&unsorted, sam_of(&many_templates())).unwrap();
    let (sorted, grouped, called) = (file("sorted.bam"), file("grouped.bam"), file("called.bam"));
    sort(&unsorted, &sorted, &[]);
    run_command("group", &sorted, &grouped, &["--strategy", "adjacency"]);
    run(&["simplex", "-i", p(&grouped), "-o", p(&called), "--min-reads", "1"]);

    let fused = |stop: &[&str], output: &Path| {
        let mut args = vec!["runall", "--start-from", "sort", "-i", p(&unsorted), "-o", p(output)];
        args.extend(["--threads", threads, "--group::strategy", "adjacency"]);
        args.extend(stop);
        run(&args);
        read_bam_output(output).1
    };
    assert_eq!(
        fused(&["--stop-after", "group"], &file("fused.group.bam")),
        read_bam_output(&grouped).1
    );
    let simplex =
        ["--stop-after", "consensus", "--consensus", "simplex", "--simplex::min-reads", "1"];
    assert_eq!(fused(&simplex, &file("fused.bam")), read_bam_output(&called).1);
}
