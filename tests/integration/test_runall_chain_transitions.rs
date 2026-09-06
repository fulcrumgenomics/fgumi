//! In-process `build_for`+run tests for each multi-stage transition `runall`
//! depends on. These validate the chain plumbing (type-erased hand-offs)
//! BEFORE the CLI is wired, per spec §9.
use std::path::{Path, PathBuf};

use fgumi_lib::aligner::AlignerOptions;
use fgumi_lib::assigner::Strategy;
#[cfg(feature = "consensus")]
use fgumi_lib::commands::codec::CodecOptions;
use fgumi_lib::commands::common::{
    CompressionOptions, QueueMemoryOptions, SchedulerOptions, ThreadingOptions,
};
use fgumi_lib::commands::correct::CorrectOptions;
#[cfg(feature = "consensus")]
use fgumi_lib::commands::duplex::DuplexOptions;
use fgumi_lib::commands::extract::ExtractRunallOptions;
use fgumi_lib::commands::filter::FilterOptions;
use fgumi_lib::commands::group::GroupOptions;
#[cfg(feature = "consensus")]
use fgumi_lib::commands::simplex::SimplexOptions;
use fgumi_lib::commands::sort::SortOptions;
use fgumi_lib::commands::zipper::ZipperOptions;
use fgumi_lib::pipeline::chains::options_bag::AlignOptions;
use fgumi_lib::pipeline::chains::validate::{
    validate_cross_stage_constraints, validate_stage_opts_present, validate_stage_progression,
};
use fgumi_lib::pipeline::chains::{
    ChainSpec, SinkSpec, SourceSpec, Stage, StageOptionsBag, build_for,
};
use fgumi_lib::sam::SamTag;
use fgumi_raw_bam::{RawRecord, SamBuilder, flags};
use tempfile::TempDir;

#[cfg(feature = "consensus")]
use crate::helpers::bam_generator::create_grouped_family;
use crate::helpers::bam_generator::{
    create_minimal_header, create_test_reference, create_umi_family, create_umi_family_at_pos,
    write_bam,
};
use crate::helpers::read_bam_output;
use crate::helpers::{aligner_binary, build_aligner_index, write_gzip_fastq};

/// Assembles a [`ChainSpec`] with the boilerplate fields every test below
/// shares (threading/compression/scheduler/queue-memory defaults, a plain
/// sequential reader, no CRC verification, a nominal `@PG` command line) —
/// mirrors the tail of `extract_to_correct_chain_builds_and_runs`'s literal,
/// factored out so each transition test only spells out what varies:
/// `stages`, `source`, `sink`, and `stage_opts`.
fn base_chain_spec(
    stages: Vec<Stage>,
    source: SourceSpec,
    sink: SinkSpec,
    stage_opts: StageOptionsBag,
) -> ChainSpec {
    ChainSpec {
        stages,
        source,
        sink,
        stage_opts,
        threading: ThreadingOptions { threads: None },
        compression: CompressionOptions::default(),
        scheduler: SchedulerOptions::default(),
        queue_memory: QueueMemoryOptions::default(),
        async_reader: false,
        read_streams: fgumi_bam_io::ReadStreams::Fixed(1),
        verify_crc: false,
        command_line: "test".into(),
    }
}

/// Build a small paired gzip FASTQ (`r1.fq.gz`, `r2.fq.gz`) with a 4 bp UMI on
/// R1 only (read structures `4M+T` / `+T`) — 2 families x 3 read pairs each.
/// The R1 UMI for each family is one of the correct-step's own whitelist
/// entries (`ACGT`, `TGCA`), so every extracted record is an exact match and
/// the correct step keeps it (a mismatched/ambiguous UMI would be dropped,
/// which would make a non-empty-output assertion meaningless). Returns
/// `(r1_path, r2_path)`.
fn write_tiny_umi_fastqs(dir: &Path) -> (PathBuf, PathBuf) {
    let r1 = dir.join("r1.fq.gz");
    let r2 = dir.join("r2.fq.gz");
    // (r1 umi (4bp, matches the correct-step whitelist), r1 template (12bp),
    // r2 template (12bp)).
    let families =
        [("ACGT", "ACGTACGTACGT", "GGTTAACCGGTT"), ("TGCA", "TGCATGCATGCA", "CCAATTGGCCAA")];
    let r1_qual = "I".repeat(4 + 12); // 4 (UMI) + 12 (template)
    let r2_qual = "I".repeat(12);
    let mut r1_owned: Vec<(String, String)> = Vec::new();
    let mut r2_owned: Vec<(String, String)> = Vec::new();
    for (fi, (umi, r1_tmpl, r2_tmpl)) in families.iter().enumerate() {
        for i in 0..3 {
            let name = format!("fam{fi}_{i}");
            r1_owned.push((name.clone(), format!("{umi}{r1_tmpl}"))); // 4 + 12 = 16
            r2_owned.push((name, (*r2_tmpl).to_string())); // 12
        }
    }
    let r1_slices: Vec<(&str, &str, &str)> =
        r1_owned.iter().map(|(n, s)| (n.as_str(), s.as_str(), r1_qual.as_str())).collect();
    let r2_slices: Vec<(&str, &str, &str)> =
        r2_owned.iter().map(|(n, s)| (n.as_str(), s.as_str(), r2_qual.as_str())).collect();
    write_gzip_fastq(&r1, &r1_slices);
    write_gzip_fastq(&r2, &r2_slices);
    (r1, r2)
}

/// Extract→Correct is the one `runall` transition that is unwired at the
/// chain-builder level: `add_correct` used to unconditionally prepend a
/// `GroupByQueryname` step (input `DecodedRecordBatch`) onto whatever tail
/// preceded it. Extract's output tail is `BamTemplateBatch`, so feeding an
/// extract-fed chain into correct was a type-erased runtime panic, not a
/// build-time error. This builds `[Stage::Extract, Stage::Correct]` directly
/// (bypassing the not-yet-existing `RunAll` command) and asserts it runs to
/// completion and produces a non-empty corrected unmapped BAM.
#[test]
fn extract_to_correct_chain_builds_and_runs() {
    let tmp = TempDir::new().unwrap();
    let (r1, r2) = write_tiny_umi_fastqs(tmp.path());
    let out = tmp.path().join("corrected.bam");

    let extract = ExtractRunallOptions {
        inputs: vec![r1.clone(), r2.clone()],
        read_structures: vec!["4M+T".parse().unwrap(), "+T".parse().unwrap()],
        sample: "s1".into(),
        library: "lib1".into(),
        ..ExtractRunallOptions::default()
    };
    let correct = CorrectOptions {
        umis: vec!["ACGT".into(), "TGCA".into()], // 4 bp whitelist, not 2 bp
        min_distance_diff: 1,
        ..CorrectOptions::default()
    };

    let bag = StageOptionsBag {
        extract: Some(extract.to_extract_options()),
        correct: Some(correct),
        ..Default::default()
    };

    let spec = ChainSpec {
        stages: vec![Stage::Extract, Stage::Correct],
        source: SourceSpec::fastqs(
            vec![r1, r2],
            vec!["4M+T".parse().unwrap(), "+T".parse().unwrap()],
        )
        .unwrap(),
        sink: SinkSpec::Bam(out.clone()),
        stage_opts: bag,
        threading: ThreadingOptions { threads: None },
        compression: CompressionOptions::default(),
        scheduler: SchedulerOptions::default(),
        queue_memory: QueueMemoryOptions::default(),
        async_reader: false,
        read_streams: fgumi_bam_io::ReadStreams::Fixed(1),
        verify_crc: false,
        command_line: "test".into(),
    };

    build_for(spec).expect("build").run().expect("run without panic");
    let (_, records) = read_bam_output(&out);
    assert!(!records.is_empty(), "corrected BAM should have records");
}

// ─────────────────────────── Sort→Group, Group→consensus, Consensus→Filter ───────────────────────────

/// `Sort→Group`: an unsorted, mapped, `RX`-tagged BAM (two UMI families at two
/// distinct positions, so sort has real work to do) sorted to
/// template-coordinate order and grouped by UMI adjacency.
#[test]
fn sort_to_group_chain_builds_and_runs() {
    let tmp = TempDir::new().unwrap();
    let header = create_minimal_header("chr1", 10000);
    let mut records = Vec::new();
    records.extend(create_umi_family_at_pos("ACGT", 3, "famA", "ACGTACGTACGT", 30, 500));
    records.extend(create_umi_family_at_pos("TGCA", 3, "famB", "ACGTACGTACGT", 30, 100));
    let input = tmp.path().join("unsorted.bam");
    write_bam(&input, &header, &records);

    let out = tmp.path().join("grouped.bam");
    let sort = SortOptions::default(); // order defaults to TemplateCoordinate
    let mut group = GroupOptions { strategy: Strategy::Adjacency, ..GroupOptions::default() };
    let (s, e) = group.resolve_strategy_and_edits();
    group.effective_strategy = s;
    group.effective_edits = e;

    let bag = StageOptionsBag { sort: Some(sort), group: Some(group), ..Default::default() };
    let spec = base_chain_spec(
        vec![Stage::Sort, Stage::Group],
        SourceSpec::Bam(input),
        SinkSpec::Bam(out.clone()),
        bag,
    );
    build_for(spec).expect("build").run().expect("run");
    let (_, records) = read_bam_output(&out);
    assert!(!records.is_empty(), "grouped BAM should have records");
}

/// `Group→Simplex`: a mapped, `RX`-tagged single-UMI family, grouped by
/// identity (exact-match UMI) directly into a simplex consensus call.
#[cfg(feature = "consensus")]
#[test]
fn group_to_simplex_chain_builds_and_runs() {
    let tmp = TempDir::new().unwrap();
    let header = create_minimal_header("chr1", 10000);
    let records = create_umi_family("ACGT", 3, "fam", "ACGTACGTACGT", 30);
    let input = tmp.path().join("input.bam");
    write_bam(&input, &header, &records);

    let out = tmp.path().join("simplex.bam");
    let mut group = GroupOptions { strategy: Strategy::Identity, ..GroupOptions::default() };
    let (s, e) = group.resolve_strategy_and_edits();
    group.effective_strategy = s;
    group.effective_edits = e;
    let simplex = SimplexOptions { min_reads: 1, ..SimplexOptions::default() };

    let bag = StageOptionsBag { group: Some(group), simplex: Some(simplex), ..Default::default() };
    let spec = base_chain_spec(
        vec![Stage::Group, Stage::Simplex],
        SourceSpec::Bam(input),
        SinkSpec::Bam(out.clone()),
        bag,
    );
    build_for(spec).expect("build").run().expect("run");
    let (_, records) = read_bam_output(&out);
    assert!(!records.is_empty(), "simplex output should have records");
}

/// One paired-UMI (`AAA-TTT`/`TTT-AAA`) read pair for one duplex strand,
/// oriented FR (top strand, `is_b_strand = false`) or RF (bottom strand,
/// `is_b_strand = true`) — the shape `Strategy::Paired` grouping expects so it
/// can assign both strands of one molecule the same base id with `/A`/`/B`
/// suffixes. Mirrors `helpers::bam_generator::create_duplex_grouped_family`'s
/// shape but adds `MC` tags: `RecordPositionGrouper` (the position-based UMI
/// grouper `Group` wires for paired-end input) requires them.
fn duplex_strand_pair(
    name: &str,
    umi: &str,
    sequence: &str,
    quality: u8,
    ref_start: i32,
    is_b_strand: bool,
) -> (RawRecord, RawRecord) {
    let seq = sequence.as_bytes();
    let read_len = seq.len();
    let cigar_op = u32::try_from(read_len).expect("read_len fits u32") << 4;
    let mate_cigar = format!("{read_len}M");
    let span = i32::try_from(read_len + 100).expect("template span fits i32");

    // B strand is the same molecule read the other way round: R1 reverse at
    // the far position, R2 forward at the near one.
    let (r1_start, r2_start, r1_rev, r2_rev) = if is_b_strand {
        (ref_start + 100, ref_start, true, false)
    } else {
        (ref_start, ref_start + 100, false, true)
    };
    let tlen = if is_b_strand { -span } else { span };

    let build = |first: bool| {
        let (pos, mate_pos, rev, mate_rev) = if first {
            (r1_start, r2_start, r1_rev, r2_rev)
        } else {
            (r2_start, r1_start, r2_rev, r1_rev)
        };
        let segment = if first { flags::FIRST_SEGMENT } else { flags::LAST_SEGMENT };
        let mut b = SamBuilder::new();
        b.read_name(name.as_bytes())
            .sequence(seq)
            .qualities(&vec![quality; read_len])
            .flags(
                flags::PAIRED
                    | segment
                    | if rev { flags::REVERSE } else { 0 }
                    | if mate_rev { flags::MATE_REVERSE } else { 0 },
            )
            .ref_id(0)
            .pos(pos - 1)
            .mapq(60)
            .cigar_ops(&[cigar_op])
            .mate_ref_id(0)
            .mate_pos(mate_pos - 1)
            .template_length(if first { tlen } else { -tlen })
            .add_string_tag(SamTag::RX, umi.as_bytes())
            .add_string_tag(SamTag::MC, mate_cigar.as_bytes());
        b.build()
    };

    (build(true), build(false))
}

/// `depth` read pairs on each of the two duplex strands of one molecule.
fn create_duplex_umi_family(
    depth: usize,
    sequence: &str,
    quality: u8,
    ref_start: i32,
) -> Vec<RawRecord> {
    let mut records = Vec::with_capacity(depth * 4);
    for (umi, is_b_strand) in [("AAA-TTT", false), ("TTT-AAA", true)] {
        for i in 0..depth {
            let (r1, r2) = duplex_strand_pair(
                &format!("{umi}_{i}"),
                umi,
                sequence,
                quality,
                ref_start,
                is_b_strand,
            );
            records.push(r1);
            records.push(r2);
        }
    }
    records
}

/// `Group→Duplex`: a paired-UMI (`AAA-TTT`/`TTT-AAA`) two-strand molecule,
/// grouped with `Strategy::Paired` (the strategy duplex consensus requires —
/// see `validate_cross_stage_constraints` Rule 2) directly into a duplex call.
#[cfg(feature = "consensus")]
#[test]
fn group_to_duplex_chain_builds_and_runs() {
    let tmp = TempDir::new().unwrap();
    let header = create_minimal_header("chr1", 10000);
    let records = create_duplex_umi_family(2, "ACGTACGTACGT", 30, 500);
    let input = tmp.path().join("input.bam");
    write_bam(&input, &header, &records);

    let out = tmp.path().join("duplex.bam");
    let mut group = GroupOptions { strategy: Strategy::Paired, ..GroupOptions::default() };
    let (s, e) = group.resolve_strategy_and_edits();
    group.effective_strategy = s;
    group.effective_edits = e;
    let duplex = DuplexOptions::default(); // min_reads = [1]

    let bag = StageOptionsBag { group: Some(group), duplex: Some(duplex), ..Default::default() };
    let spec = base_chain_spec(
        vec![Stage::Group, Stage::Duplex],
        SourceSpec::Bam(input),
        SinkSpec::Bam(out.clone()),
        bag,
    );
    build_for(spec).expect("build").run().expect("run");
    let (_, records) = read_bam_output(&out);
    assert!(!records.is_empty(), "duplex output should have records");
}

/// One CODEC-shaped read pair for `Group→Codec`: R1 forward, R2 reverse,
/// fully overlapping at the same position (mirrors real CODEC sequencing,
/// where R1/R2 read opposite strands of the same short fragment), sharing an
/// `RX` UMI so the upstream `Group` stage — not a pre-set `MI` — assigns the
/// molecule id.
fn create_codec_umi_pair(
    name: &str,
    seq: &[u8],
    qual: &[u8],
    ref_start: i32,
    umi: &str,
) -> (RawRecord, RawRecord) {
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

/// `Group→Codec`: two overlapping CODEC-shaped read pairs sharing one `RX`
/// UMI, grouped by identity into a single molecule and CODEC-consensus-called.
#[cfg(feature = "consensus")]
#[test]
fn group_to_codec_chain_builds_and_runs() {
    let tmp = TempDir::new().unwrap();
    let header = create_minimal_header("chr1", 10000);
    let mut records = Vec::new();
    for i in 0..2 {
        let (r1, r2) =
            create_codec_umi_pair(&format!("pair{i}"), b"ACGTACGTAC", &[30; 10], 500, "ACGT");
        records.push(r1);
        records.push(r2);
    }
    let input = tmp.path().join("input.bam");
    write_bam(&input, &header, &records);

    let out = tmp.path().join("codec.bam");
    let mut group = GroupOptions::default(); // Strategy::Identity, edits 0
    let (s, e) = group.resolve_strategy_and_edits();
    group.effective_strategy = s;
    group.effective_edits = e;
    let codec = CodecOptions::default(); // min_reads = 1, min_duplex_length = 1

    let bag = StageOptionsBag { group: Some(group), codec: Some(codec), ..Default::default() };
    let spec = base_chain_spec(
        vec![Stage::Group, Stage::Codec],
        SourceSpec::Bam(input),
        SinkSpec::Bam(out.clone()),
        bag,
    );
    build_for(spec).expect("build").run().expect("run");
    let (_, records) = read_bam_output(&out);
    assert!(!records.is_empty(), "codec output should have records");
}

/// `Consensus→Filter`: a simplex consensus call immediately followed by
/// `Filter` in the same chain — the `runall consensus → filter` shape that
/// `validate_stage_progression` special-cases (a consensus stage may be
/// followed by exactly one trailing `Filter`).
#[cfg(feature = "consensus")]
#[test]
fn simplex_to_filter_chain_builds_and_runs() {
    let tmp = TempDir::new().unwrap();
    let header = create_minimal_header("chr1", 10000);
    let records = create_grouped_family("ACGT", "1", 3, "fam", "ACGTACGTACGT", 30);
    let input = tmp.path().join("grouped.bam");
    write_bam(&input, &header, &records);

    let out = tmp.path().join("filtered.bam");
    let simplex = SimplexOptions { min_reads: 1, ..SimplexOptions::default() };
    let filter = FilterOptions { min_reads: vec![1], ..FilterOptions::default() };

    let bag =
        StageOptionsBag { simplex: Some(simplex), filter: Some(filter), ..Default::default() };
    let spec = base_chain_spec(
        vec![Stage::Simplex, Stage::Filter],
        SourceSpec::Bam(input),
        SinkSpec::Bam(out.clone()),
        bag,
    );
    build_for(spec).expect("build").run().expect("run");
    let (_, records) = read_bam_output(&out);
    assert!(!records.is_empty(), "filtered output should have records");
}

// ─────────────────────────────────── Zipper→Sort ───────────────────────────────────

/// Matched (unmapped, mapped, reference) triple for the zipper stage: an
/// unmapped BAM with `RX`/`QX` tags and a mapped BAM with the same read
/// names aligned to `chr1`, both in the same (queryname-matching) order.
/// Mirrors `feat-runall`'s `test_runall_parity.rs::zipper_fixture`.
struct ZipperFixture {
    unmapped: PathBuf,
    mapped: PathBuf,
    reference: PathBuf,
}

fn zipper_fixture(dir: &Path) -> ZipperFixture {
    let unmapped_path = dir.join("zipper_unmapped.bam");
    let mapped_path = dir.join("zipper_mapped.bam");
    let reference = create_test_reference(dir);

    let names = ["read1", "read2", "read3", "read4"];

    let unmapped_records: Vec<RawRecord> = names
        .iter()
        .enumerate()
        .map(|(i, name)| {
            let mut b = SamBuilder::new();
            b.read_name(name.as_bytes())
                .sequence(b"ACGTACGT")
                .qualities(&[30; 8])
                .flags(flags::UNMAPPED)
                .add_string_tag(SamTag::RX, format!("AACC{i}").as_bytes())
                .add_string_tag(SamTag::QX, b"IIIIIIII");
            b.build()
        })
        .collect();
    write_bam(&unmapped_path, &noodles::sam::Header::default(), &unmapped_records);

    let mapped_header = create_minimal_header("chr1", 10000);
    let mapped_records: Vec<RawRecord> = names
        .iter()
        .enumerate()
        .map(|(i, name)| {
            let mut b = SamBuilder::new();
            b.read_name(name.as_bytes())
                .ref_id(0)
                .pos(i32::try_from(100 + 50 * i).unwrap_or(0))
                .mapq(60)
                .flags(0)
                .cigar_ops(&[8u32 << 4])
                .sequence(b"ACGTACGT")
                .qualities(&[30; 8]);
            b.build()
        })
        .collect();
    write_bam(&mapped_path, &mapped_header, &mapped_records);

    ZipperFixture { mapped: mapped_path, unmapped: unmapped_path, reference }
}

/// `Zipper→Sort`: zipper merges the matched mapped+unmapped BAMs into
/// `BamTemplateBatch`es and sort (template-coordinate, the default) consumes
/// them directly — no aligner subprocess involved, so this runs live
/// unconditionally.
#[test]
fn zipper_to_sort_chain_builds_and_runs() {
    let tmp = TempDir::new().unwrap();
    let inputs = zipper_fixture(tmp.path());

    let out = tmp.path().join("sorted.bam");
    let zipper = ZipperOptions::default();
    let sort = SortOptions::default(); // order defaults to TemplateCoordinate

    let bag = StageOptionsBag { zipper: Some(zipper), sort: Some(sort), ..Default::default() };
    let spec = base_chain_spec(
        vec![Stage::Zipper, Stage::Sort],
        SourceSpec::PairedBams {
            unmapped: inputs.unmapped,
            mapped: inputs.mapped,
            reference: inputs.reference,
        },
        SinkSpec::Bam(out.clone()),
        bag,
    );
    build_for(spec).expect("build").run().expect("run");
    let (_, records) = read_bam_output(&out);
    assert!(!records.is_empty(), "sorted zipper output should have records");
}

// ─────────────────────────── Align-stage transitions ───────────────────────────
//
// `Stage::Align` spawns a real aligner subprocess. These tests are gated by
// `which::which()` (mirroring `feat-runall`'s `run_aam_parity_test` guard):
// when a `bwa-mem3`/`bwa` binary is on `PATH`, the chain is built AND run
// end-to-end against a tiny real reference; otherwise the assembled
// `ChainSpec` is validated instead (`validate_stage_progression` /
// `validate_stage_opts_present` / `validate_cross_stage_constraints`), so the
// transition's wiring is still exercised at the spec level everywhere.

/// Writes a single-end unmapped BAM: `n` reads, each carrying an `RX` UMI tag
/// and a 12bp template sequence (`ACGTACGTACGT`) that is an exact substring of
/// `create_test_reference`'s repeating `ACGTACGT` `chr1`, so a real aligner
/// maps them (even at low MAPQ, given the repeat). Used as the Align-stage
/// "unmapped BAM" source/tail below.
fn write_unmapped_umi_bam(path: &Path, umi: &str, n: usize) {
    let header = noodles::sam::Header::default();
    let records: Vec<RawRecord> = (0..n)
        .map(|i| {
            let mut b = SamBuilder::new();
            b.read_name(format!("read{i}").as_bytes())
                .sequence(b"ACGTACGTACGT")
                .qualities(&[30; 12])
                .flags(flags::UNMAPPED)
                .add_string_tag(SamTag::RX, umi.as_bytes());
            b.build()
        })
        .collect();
    write_bam(path, &header, &records);
}

/// Asserts the three chain-spec validators all accept `spec` — the
/// no-aligner-available fallback for the `Stage::Align` transition tests.
fn assert_chain_spec_validates(spec: &ChainSpec) {
    validate_stage_progression(spec).expect("stage progression should validate");
    validate_stage_opts_present(spec).expect("stage options should validate");
    validate_cross_stage_constraints(spec).expect("cross-stage constraints should validate");
}

/// `Extract→Align`: a single-end FASTQ (`4M+T`) extracted to an unmapped BAM
/// and fed directly into Align (`add_align` skips `GroupByQueryname` when the
/// upstream tail is already `BamTemplateBatch`).
#[test]
fn extract_to_align_chain_builds_and_runs_or_validates() {
    let tmp = TempDir::new().unwrap();
    let r1 = tmp.path().join("r1.fq.gz");
    // 4bp UMI ("ACGT") + 12bp template exactly matching the test reference.
    write_gzip_fastq(&r1, &[("read0", "ACGTACGTACGTACGT", "IIIIIIIIIIIIIIII")]);

    let extract = ExtractRunallOptions {
        inputs: vec![r1.clone()],
        read_structures: vec!["4M+T".parse().unwrap()],
        sample: "s1".into(),
        library: "lib1".into(),
        ..ExtractRunallOptions::default()
    };
    let source =
        SourceSpec::fastqs(vec![r1], vec!["4M+T".parse().unwrap()]).expect("valid fastq source");
    let out = tmp.path().join("aligned.bam");

    if let Some(bin) = aligner_binary() {
        let reference = create_test_reference(tmp.path());
        build_aligner_index(&reference, bin);
        let align = AlignOptions {
            aligner: AlignerOptions {
                command: Some(format!("{bin} mem {{ref}} /dev/stdin")),
                ..AlignerOptions::default()
            },
            reference,
            aligner_bin: None,
        };
        let bag = StageOptionsBag {
            extract: Some(extract.to_extract_options()),
            aligner: Some(align),
            ..Default::default()
        };
        let spec = base_chain_spec(
            vec![Stage::Extract, Stage::Align],
            source,
            SinkSpec::Bam(out.clone()),
            bag,
        );
        build_for(spec).expect("build").run().expect("run");
        let (_, records) = read_bam_output(&out);
        assert!(!records.is_empty(), "aligned BAM should have records");
    } else {
        eprintln!(
            "no bwa/bwa-mem3 on PATH: validating the Extract→Align spec instead of running it"
        );
        let align = AlignOptions {
            aligner: AlignerOptions {
                command: Some("bwa mem {ref} /dev/stdin".into()),
                ..AlignerOptions::default()
            },
            reference: tmp.path().join("ref.fa"), // need not exist for spec validation
            aligner_bin: None,
        };
        let bag = StageOptionsBag {
            extract: Some(extract.to_extract_options()),
            aligner: Some(align),
            ..Default::default()
        };
        let spec =
            base_chain_spec(vec![Stage::Extract, Stage::Align], source, SinkSpec::Bam(out), bag);
        assert_chain_spec_validates(&spec);
    }
}

/// `Correct→Align`: an unmapped BAM with an exact-whitelist-match `RX` tag,
/// corrected, then fed directly into Align (same `BamTemplateBatch` hand-off
/// as `Extract→Align`).
#[test]
fn correct_to_align_chain_builds_and_runs_or_validates() {
    let tmp = TempDir::new().unwrap();
    let unmapped = tmp.path().join("unmapped.bam");
    write_unmapped_umi_bam(&unmapped, "ACGT", 3);

    let correct = CorrectOptions {
        umis: vec!["ACGT".into(), "TGCA".into()], // 4bp whitelist, not 2bp
        min_distance_diff: 1,
        ..CorrectOptions::default()
    };
    let out = tmp.path().join("aligned.bam");

    if let Some(bin) = aligner_binary() {
        let reference = create_test_reference(tmp.path());
        build_aligner_index(&reference, bin);
        let align = AlignOptions {
            aligner: AlignerOptions {
                command: Some(format!("{bin} mem {{ref}} /dev/stdin")),
                ..AlignerOptions::default()
            },
            reference,
            aligner_bin: None,
        };
        let bag =
            StageOptionsBag { correct: Some(correct), aligner: Some(align), ..Default::default() };
        let spec = base_chain_spec(
            vec![Stage::Correct, Stage::Align],
            SourceSpec::Bam(unmapped),
            SinkSpec::Bam(out.clone()),
            bag,
        );
        build_for(spec).expect("build").run().expect("run");
        let (_, records) = read_bam_output(&out);
        assert!(!records.is_empty(), "aligned BAM should have records");
    } else {
        eprintln!(
            "no bwa/bwa-mem3 on PATH: validating the Correct→Align spec instead of running it"
        );
        let align = AlignOptions {
            aligner: AlignerOptions {
                command: Some("bwa mem {ref} /dev/stdin".into()),
                ..AlignerOptions::default()
            },
            reference: tmp.path().join("ref.fa"),
            aligner_bin: None,
        };
        let bag =
            StageOptionsBag { correct: Some(correct), aligner: Some(align), ..Default::default() };
        let spec = base_chain_spec(
            vec![Stage::Correct, Stage::Align],
            SourceSpec::Bam(unmapped),
            SinkSpec::Bam(out),
            bag,
        );
        assert_chain_spec_validates(&spec);
    }
}

/// `Align→Sort`: Align is the FIRST stage (reads straight from
/// `SourceSpec::Bam`, prepending `GroupByQueryname` itself since the upstream
/// tail is not yet `BamTemplateBatch`), and Sort consumes Align's
/// `BamTemplateBatch` output directly.
#[test]
fn align_to_sort_chain_builds_and_runs_or_validates() {
    let tmp = TempDir::new().unwrap();
    let unmapped = tmp.path().join("unmapped.bam");
    write_unmapped_umi_bam(&unmapped, "ACGT", 3);

    let sort = SortOptions::default(); // order defaults to TemplateCoordinate
    let out = tmp.path().join("sorted.bam");

    if let Some(bin) = aligner_binary() {
        let reference = create_test_reference(tmp.path());
        build_aligner_index(&reference, bin);
        let align = AlignOptions {
            aligner: AlignerOptions {
                command: Some(format!("{bin} mem {{ref}} /dev/stdin")),
                ..AlignerOptions::default()
            },
            reference,
            aligner_bin: None,
        };
        let bag = StageOptionsBag { aligner: Some(align), sort: Some(sort), ..Default::default() };
        let spec = base_chain_spec(
            vec![Stage::Align, Stage::Sort],
            SourceSpec::Bam(unmapped),
            SinkSpec::Bam(out.clone()),
            bag,
        );
        build_for(spec).expect("build").run().expect("run");
        let (_, records) = read_bam_output(&out);
        assert!(!records.is_empty(), "sorted, aligned output should have records");
    } else {
        eprintln!("no bwa/bwa-mem3 on PATH: validating the Align→Sort spec instead of running it");
        let align = AlignOptions {
            aligner: AlignerOptions {
                command: Some("bwa mem {ref} /dev/stdin".into()),
                ..AlignerOptions::default()
            },
            reference: tmp.path().join("ref.fa"),
            aligner_bin: None,
        };
        let bag = StageOptionsBag { aligner: Some(align), sort: Some(sort), ..Default::default() };
        let spec = base_chain_spec(
            vec![Stage::Align, Stage::Sort],
            SourceSpec::Bam(unmapped),
            SinkSpec::Bam(out),
            bag,
        );
        assert_chain_spec_validates(&spec);
    }
}
