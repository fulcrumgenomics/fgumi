//! Fixture + BAM-comparison helpers shared by the bwa-mem3 align-stage parity
//! suites (`align_inproc_parity.rs`, `align_subprocess_split_parity.rs`).
//!
//! Included by both as `mod align_common;` (a directory module, so cargo does not
//! build it as a test target of its own). Each suite uses a subset of the
//! helpers, hence the module-wide `dead_code` allow.
//!
//! The fixture is built to exercise real alignment behaviour rather than only
//! exact matches: a seeded-random reference with a duplicated segment (so reads
//! there multi-map, get MAPQ 0 and an `XA` tag), and reads drawn in nine kinds —
//! exact, substitutions, a small insertion, a small deletion, unmappable random
//! sequence, chimeric (two distant loci, so bwa emits a supplementary),
//! discordant pairs, reads inside the repeat, and zero-length reads (an empty
//! single, or a pair with one or both mates empty). Every unmapped input record
//! carries an `RX` tag and some are QC-fail, so tag and flag transfer through the
//! zipper merge is part of what the suites compare.

#![allow(dead_code)]

use std::collections::{BTreeMap, BTreeSet};
use std::io::Write as _;
use std::path::{Path, PathBuf};

use fgumi_raw_bam::flags::{
    FIRST_SEGMENT, LAST_SEGMENT, PAIRED, QC_FAIL, REVERSE, SECONDARY, SUPPLEMENTARY, UNMAPPED,
};
use fgumi_raw_bam::{RawRecord, RawRecordView, RawTagsView, SamBuilder, SamTag};
use rand::rngs::StdRng;
use rand::{RngExt, SeedableRng};

#[path = "../integration/helpers/aligner.rs"]
mod aligner;

// ---------------------------------------------------------------------------
// Environment gating
// ---------------------------------------------------------------------------

/// Whether missing tools must fail the test rather than skip it. The CI
/// e2e-parity job sets `FGUMI_BWA_MEM3_REQUIRE_TOOLS` so a broken environment is
/// loud instead of a silent no-op; everywhere else the suites skip-as-pass.
pub fn tools_required() -> bool {
    std::env::var_os("FGUMI_BWA_MEM3_REQUIRE_TOOLS").is_some()
}

/// The reference `bwa-mem3` binary from `BWA_MEM3_BIN`, with no `PATH` fallback.
pub fn bwa_mem3_bin_from_env() -> Option<String> {
    std::env::var_os("BWA_MEM3_BIN").map(|bin| bin.to_string_lossy().into_owned())
}

// ---------------------------------------------------------------------------
// Fixture reference
// ---------------------------------------------------------------------------

/// Length of every simulated read before any indel (2×150, the production shape).
pub const READ_LEN: usize = 150;
/// Length of the fixture reference contig.
pub const REF_LEN: usize = 24_000;
/// Start of the segment that is duplicated at [`REPEAT_DST`].
const REPEAT_SRC: usize = 3_000;
/// Length of the duplicated segment: longer than the largest simulated insert,
/// so a whole pair can be drawn from inside it.
const REPEAT_LEN: usize = 1_200;
/// Where the copy of `REPEAT_SRC..REPEAT_SRC + REPEAT_LEN` is placed.
const REPEAT_DST: usize = 15_000;
/// Offset from a chimeric read's first locus to its second, far enough apart
/// that bwa reports the second part as a separate (supplementary) alignment.
const CHIMERA_OFFSET: usize = 9_000;
/// Offset from a discordant pair's R1 locus to its R2 locus (far beyond any
/// plausible insert size).
const DISCORDANT_OFFSET: usize = 6_000;

/// The fixture reference sequence: seeded-random bases with
/// `REPEAT_SRC..REPEAT_SRC + REPEAT_LEN` copied to `REPEAT_DST`.
pub fn fixture_reference_seq() -> Vec<u8> {
    let mut rng = StdRng::seed_from_u64(0x2545_F491_4F6C_DD1D);
    let mut seq: Vec<u8> = (0..REF_LEN).map(|_| random_base(&mut rng)).collect();
    seq.copy_within(REPEAT_SRC..REPEAT_SRC + REPEAT_LEN, REPEAT_DST);
    seq
}

/// A written, indexed fixture reference.
pub struct FixtureRef {
    /// The FASTA path (also the bwa-mem3 index prefix and `runall --ref`).
    pub fasta: PathBuf,
    /// The reference bases, for read simulation.
    pub seq: Vec<u8>,
}

/// Write the fixture reference (`chr1` on one line) plus its `.fai` and `.dict`
/// into `dir`, and build its bwa-mem3 index with `bin`. The index must come from
/// a CLI at the vendored commit so the in-process backend can load it.
pub fn write_fixture_reference(dir: &Path, bin: &str) -> FixtureRef {
    let seq = fixture_reference_seq();
    let fasta = dir.join("ref.fa");
    let len = seq.len();

    let mut file = std::fs::File::create(&fasta).expect("create reference FASTA");
    writeln!(file, ">chr1").expect("write FASTA name");
    file.write_all(&seq).expect("write FASTA seq");
    writeln!(file).expect("write FASTA newline");
    file.flush().expect("flush FASTA");

    std::fs::write(dir.join("ref.fa.fai"), format!("chr1\t{len}\t6\t{len}\t{}\n", len + 1))
        .expect("write .fai");
    std::fs::write(
        dir.join("ref.dict"),
        format!("@HD\tVN:1.6\tSO:unsorted\n@SQ\tSN:chr1\tLN:{len}\n"),
    )
    .expect("write .dict");

    aligner::build_aligner_index(&fasta, bin);
    FixtureRef { fasta, seq }
}

// ---------------------------------------------------------------------------
// Read simulation
// ---------------------------------------------------------------------------

/// The layout of the simulated unmapped input.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum Input {
    /// All paired (FR): R1 forward, R2 reverse-complemented.
    Paired,
    /// All single-end, alternating forward / reverse-complement.
    Single,
    /// Alternating PE and SE templates (PE at even indices), so the `-K` cut can
    /// land between a pair's two reads (see [`mid_pair_chunk`]).
    Mixed,
}

/// What a simulated read (or pair) is designed to exercise in the aligner.
#[derive(Debug, Clone, Copy)]
enum Kind {
    Exact,
    Substitutions,
    Insertion,
    Deletion,
    Unmapped,
    Chimeric,
    Discordant,
    Repeat,
    /// A zero-length SEQ: an empty single, or a pair with one or both mates
    /// empty. bwa-mem3 emits an empty read as unmapped (placed at its mate when
    /// the mate maps).
    Empty,
}

/// Every [`Kind`], cycled through per template shape.
const KINDS: [Kind; 9] = [
    Kind::Exact,
    Kind::Substitutions,
    Kind::Insertion,
    Kind::Deletion,
    Kind::Unmapped,
    Kind::Chimeric,
    Kind::Discordant,
    Kind::Repeat,
    Kind::Empty,
];

/// One simulated template: a name, one (SE) or two (PE) read bodies, and the
/// `RX` / QC-fail state its unmapped records carry.
pub struct Template {
    pub name: String,
    pub reads: Vec<Vec<u8>>,
    pub paired: bool,
    pub rx: String,
    pub qc_fail: bool,
}

fn random_base(rng: &mut StdRng) -> u8 {
    b"ACGT"[rng.random_range(0..4)]
}

fn random_seq(rng: &mut StdRng, len: usize) -> Vec<u8> {
    (0..len).map(|_| random_base(rng)).collect()
}

/// Replace `n` random positions of `read` with a different base.
fn substitute(rng: &mut StdRng, read: &mut [u8], n: usize) {
    for _ in 0..n {
        let at = rng.random_range(0..read.len());
        let old = read[at];
        read[at] = loop {
            let b = random_base(rng);
            if b != old {
                break b;
            }
        };
    }
}

/// A `READ_LEN` read at `start`, drawn with the given kind's mutation (for kinds
/// that mutate a single read in place). `Unmapped`, `Chimeric`, `Discordant` and
/// `Repeat` choose their loci before calling this, so they use it unmutated.
fn draw_read(rng: &mut StdRng, ref_seq: &[u8], start: usize, kind: Kind) -> Vec<u8> {
    match kind {
        Kind::Substitutions => {
            let mut read = ref_seq[start..start + READ_LEN].to_vec();
            substitute(rng, &mut read, 3);
            read
        }
        Kind::Insertion => {
            let mut read = ref_seq[start..start + READ_LEN].to_vec();
            let at = rng.random_range(60..90);
            let ins = random_seq(rng, 2);
            read.splice(at..at, ins);
            read
        }
        Kind::Deletion => {
            // Draw 3 extra bases and delete 3 mid-read, so the read stays READ_LEN.
            let mut read = ref_seq[start..start + READ_LEN + 3].to_vec();
            let at = rng.random_range(60..90);
            read.drain(at..at + 3);
            read
        }
        _ => ref_seq[start..start + READ_LEN].to_vec(),
    }
}

/// Simulate `n` templates of the requested [`Input`] shape from `ref_seq`,
/// cycling through every [`Kind`] for each template shape. Deterministic: the
/// same input BAM feeds every leg of a comparison.
pub fn simulate_templates(ref_seq: &[u8], input: Input, n: usize) -> Vec<Template> {
    assert!(ref_seq.len() >= REF_LEN, "fixture reference is shorter than REF_LEN");
    let mut rng = StdRng::seed_from_u64(0x9E37_79B9_7F4A_7C15);
    let (mut pairs_seen, mut singles_seen) = (0usize, 0usize);
    let mut templates = Vec::with_capacity(n);
    for i in 0..n {
        let paired = match input {
            Input::Paired => true,
            Input::Single => false,
            Input::Mixed => i.is_multiple_of(2),
        };
        let kind = if paired {
            KINDS[pairs_seen % KINDS.len()]
        } else {
            KINDS[singles_seen % KINDS.len()]
        };
        let reads = if paired {
            pairs_seen += 1;
            simulate_pair(&mut rng, ref_seq, kind)
        } else {
            singles_seen += 1;
            let mut read = simulate_single(&mut rng, ref_seq, kind);
            if singles_seen.is_multiple_of(2) {
                read = fgumi_dna::reverse_complement(&read);
            }
            vec![read]
        };
        let rx = String::from_utf8(random_seq(&mut rng, 8)).expect("ACGT is UTF-8");
        let name = if paired { format!("pe{i}") } else { format!("se{i}") };
        templates.push(Template { name, reads, paired, rx, qc_fail: i % 7 == 3 });
    }
    templates
}

fn simulate_pair(rng: &mut StdRng, ref_seq: &[u8], kind: Kind) -> Vec<Vec<u8>> {
    let insert = READ_LEN + 50 + rng.random_range(0..150);
    match kind {
        Kind::Unmapped => vec![random_seq(rng, READ_LEN), random_seq(rng, READ_LEN)],
        Kind::Empty => {
            let start = rng.random_range(0..REF_LEN - insert);
            let r1 = ref_seq[start..start + READ_LEN].to_vec();
            let r2 =
                fgumi_dna::reverse_complement(&ref_seq[start + insert - READ_LEN..start + insert]);
            match rng.random_range(0..3u8) {
                0 => vec![Vec::new(), r2],
                1 => vec![r1, Vec::new()],
                _ => vec![Vec::new(), Vec::new()],
            }
        }
        Kind::Repeat => {
            let start = REPEAT_SRC + rng.random_range(0..=REPEAT_LEN - insert);
            let r1 = ref_seq[start..start + READ_LEN].to_vec();
            let r2 =
                fgumi_dna::reverse_complement(&ref_seq[start + insert - READ_LEN..start + insert]);
            vec![r1, r2]
        }
        Kind::Chimeric => {
            let a = rng.random_range(0..REF_LEN - CHIMERA_OFFSET - READ_LEN);
            let b = a + CHIMERA_OFFSET;
            let mut r1 = ref_seq[a..a + 80].to_vec();
            r1.extend_from_slice(&ref_seq[b..b + 70]);
            let r2 = fgumi_dna::reverse_complement(&ref_seq[a + insert - READ_LEN..a + insert]);
            vec![r1, r2]
        }
        Kind::Discordant => {
            let a = rng.random_range(0..REF_LEN - DISCORDANT_OFFSET - READ_LEN);
            let b = a + DISCORDANT_OFFSET;
            let r1 = ref_seq[a..a + READ_LEN].to_vec();
            let r2 = fgumi_dna::reverse_complement(&ref_seq[b..b + READ_LEN]);
            vec![r1, r2]
        }
        _ => {
            // +3 leaves room for the Deletion kind's extra bases.
            let start = rng.random_range(0..REF_LEN - insert - 3);
            let r1 = draw_read(rng, ref_seq, start, kind);
            let r2_start = start + insert - READ_LEN;
            let r2 = fgumi_dna::reverse_complement(&draw_read(rng, ref_seq, r2_start, kind));
            vec![r1, r2]
        }
    }
}

fn simulate_single(rng: &mut StdRng, ref_seq: &[u8], kind: Kind) -> Vec<u8> {
    match kind {
        Kind::Unmapped => random_seq(rng, READ_LEN),
        Kind::Empty => Vec::new(),
        Kind::Repeat => {
            let start = REPEAT_SRC + rng.random_range(0..=REPEAT_LEN - READ_LEN);
            ref_seq[start..start + READ_LEN].to_vec()
        }
        Kind::Chimeric => {
            let a = rng.random_range(0..REF_LEN - CHIMERA_OFFSET - READ_LEN);
            let mut read = ref_seq[a..a + 80].to_vec();
            read.extend_from_slice(&ref_seq[a + CHIMERA_OFFSET..a + CHIMERA_OFFSET + 70]);
            read
        }
        // A single read has no mate to be discordant with; draw it exact.
        _ => {
            let start = rng.random_range(0..REF_LEN - READ_LEN - 3);
            draw_read(rng, ref_seq, start, kind)
        }
    }
}

/// Write `templates` to a queryname-grouped unmapped BAM (mates adjacent), each
/// record carrying its template's `RX` and, for QC-fail templates, the QC-fail
/// flag — the shape `runall --start-from align` consumes.
pub fn write_unmapped_bam(path: &Path, templates: &[Template]) {
    use noodles::sam::alignment::io::Write as _;

    let header = noodles::sam::Header::default();
    let mut records: Vec<RawRecord> = Vec::new();
    for template in templates {
        for (segment_idx, seq) in template.reads.iter().enumerate() {
            let mut flags = if template.paired {
                let segment = if segment_idx == 0 { FIRST_SEGMENT } else { LAST_SEGMENT };
                UNMAPPED | PAIRED | segment
            } else {
                UNMAPPED
            };
            if template.qc_fail {
                flags |= QC_FAIL;
            }
            let quals = vec![30u8; seq.len()];
            let mut b = SamBuilder::new();
            b.read_name(template.name.as_bytes())
                .flags(flags)
                .sequence(seq)
                .qualities(&quals)
                .add_string_tag(SamTag::RX, template.rx.as_bytes());
            records.push(b.build());
        }
    }

    let mut writer =
        noodles::bam::io::Writer::new(std::fs::File::create(path).expect("create unmapped BAM"));
    writer.write_header(&header).expect("write unmapped header");
    for record in &records {
        let buf = fgumi_raw_bam::raw_record_to_record_buf(record, &header)
            .expect("raw_record_to_record_buf");
        writer.write_alignment_record(&header, &buf).expect("write unmapped record");
    }
    writer.try_finish().expect("finish unmapped BAM");
}

/// A forced mid-pair `-K` cut: the `--aligner::chunk-size` that produces it and
/// the paired template whose two reads it separates.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct MidPairCut {
    pub chunk_size: u64,
    pub template: String,
}

/// A `--aligner::chunk-size` that makes bwa's `-K` cut land between a pair's two
/// reads: the cumulative bases through the first R1 that sits at an even running
/// read count. bwa cuts a chunk at the first read where the bases reach `-K`
/// **and** the read count is even (`bseq_read`), so with this value the first
/// chunk ends on that R1 and its R2 opens the next chunk — bwa-mem3 then aligns
/// the two reads as unpaired singles. Bases only grow along the input, so that R1
/// is the first read to reach the cut and the split pair is exactly
/// [`MidPairCut::template`]. `None` when no such R1 exists (paired-only input
/// keeps every pair inside one chunk).
pub fn mid_pair_chunk(templates: &[Template]) -> Option<MidPairCut> {
    let (mut n_reads, mut bases) = (0u64, 0u64);
    for template in templates {
        for (segment_idx, read) in template.reads.iter().enumerate() {
            n_reads += 1;
            bases += u64::try_from(read.len()).expect("read length fits u64");
            if template.paired && segment_idx == 0 && n_reads.is_multiple_of(2) {
                return Some(MidPairCut { chunk_size: bases, template: template.name.clone() });
            }
        }
    }
    None
}

// ---------------------------------------------------------------------------
// BAM decoding
// ---------------------------------------------------------------------------

/// A decoded BAM: header text, the raw reference list, and each record's body
/// (the bytes after `block_size`).
pub struct Bam {
    pub header_text: String,
    pub refs: Vec<u8>,
    pub records: Vec<Vec<u8>>,
}

fn le_i32(bytes: &[u8], at: usize) -> usize {
    let v = i32::from_le_bytes(bytes[at..at + 4].try_into().expect("4 bytes"));
    usize::try_from(v).expect("BAM length field is non-negative")
}

/// Decompress and split a BAM into header text, reference list and records.
/// BGZF block boundaries differ between writers even for identical content, so
/// the decompressed stream — not the file bytes — is what is compared.
pub fn read_bam(path: &Path) -> Bam {
    let file = std::fs::File::open(path).unwrap_or_else(|e| panic!("open {}: {e}", path.display()));
    let mut reader = noodles::bgzf::io::Reader::new(std::io::BufReader::new(file));
    let mut raw = Vec::new();
    std::io::Read::read_to_end(&mut reader, &mut raw)
        .unwrap_or_else(|e| panic!("decompress {}: {e}", path.display()));
    assert!(raw.len() >= 8 && &raw[0..4] == b"BAM\x01", "{} is not a BAM", path.display());

    let text_end = 8 + le_i32(&raw, 4);
    let header_text = String::from_utf8_lossy(&raw[8..text_end]).into_owned();
    let mut at = text_end;
    let n_ref = le_i32(&raw, at);
    at += 4;
    for _ in 0..n_ref {
        at += 4 + le_i32(&raw, at); // l_name + name
        at += 4; // l_ref
    }
    let refs = raw[text_end..at].to_vec();

    let mut records = Vec::new();
    while at + 4 <= raw.len() {
        let body = at + 4;
        let end = body + le_i32(&raw, at);
        records.push(raw[body..end].to_vec());
        at = end;
    }
    Bam { header_text, refs, records }
}

/// A record body's FLAG.
pub fn flags(body: &[u8]) -> u16 {
    u16::from_le_bytes(body[14..16].try_into().expect("2 bytes"))
}

/// A record body's read name (without the NUL).
pub fn name(body: &[u8]) -> &[u8] {
    let l_read_name = usize::from(body[8]);
    &body[32..32 + l_read_name - 1]
}

/// Whether a record is a primary alignment (neither secondary nor supplementary).
pub fn is_primary(body: &[u8]) -> bool {
    flags(body) & (SECONDARY | SUPPLEMENTARY) == 0
}

/// The value bytes of a string tag on a record, if present.
pub fn string_tag(body: &[u8], tag: SamTag) -> Option<&[u8]> {
    let tag: [u8; 2] = tag.into();
    RawTagsView::new(fgumi_raw_bam::aux_data_slice(body))
        .iter()
        .find(|entry| entry.tag == tag)
        .map(|entry| entry.value_bytes.strip_suffix(b"\0").unwrap_or(entry.value_bytes))
}

// ---------------------------------------------------------------------------
// Normalized comparison
// ---------------------------------------------------------------------------

/// Whether a normalization sorts each record's optional tags.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum TagOrder {
    /// Keep the emitted order (runs of one backend must agree exactly).
    Keep,
    /// Compare tags as a set (the two backends emit tags in different orders).
    Sort,
}

/// One aux tag: identifier, BAM type byte (so an integer's width is compared),
/// and raw value bytes.
pub type AuxEntry = ([u8; 2], u8, Vec<u8>);

/// A record split into its fixed/core bytes (everything before the aux data)
/// and its aux entries.
#[derive(Debug, PartialEq, Eq)]
pub struct NormRecord {
    pub core: Vec<u8>,
    pub aux: Vec<AuxEntry>,
}

/// A BAM normalized for comparison; see [`normalize`].
pub struct NormBam {
    pub header: String,
    pub refs: Vec<u8>,
    pub records: Vec<NormRecord>,
}

/// Normalize `bam` for comparison. Exactly these differences are removed, as
/// they carry no alignment meaning:
///
/// 1. the `CL:` field of every `@PG` header line (the command line differs
///    between presets and records a per-run output path);
/// 2. a git dev suffix on `@PG` `VN:` (see [`without_git_suffix`]): a CLI built
///    from a git checkout reports e.g. `0.12.0-ea58288`, the vendored build
///    `0.12.0`; the release version is still compared;
/// 3. with [`TagOrder::Sort`], the order of each record's optional tags.
///
/// Everything else — the rest of the header, the reference list, every core
/// field and every tag's identifier, type (including integer width) and value —
/// is compared byte for byte.
pub fn normalize(bam: &Bam, order: TagOrder) -> NormBam {
    let header = bam
        .header_text
        .lines()
        .map(|line| {
            if !line.starts_with("@PG") {
                return line.to_owned();
            }
            line.split('\t')
                .filter(|f| !f.starts_with("CL:"))
                // `without_git_suffix` returns a prefix of `vn`, so keep that
                // much of the field after the `VN:` tag.
                .map(|f| {
                    f.strip_prefix("VN:")
                        .map_or(f, |vn| &f[.."VN:".len() + without_git_suffix(vn).len()])
                })
                .collect::<Vec<_>>()
                .join("\t")
        })
        .collect::<Vec<_>>()
        .join("\n");
    let records = bam
        .records
        .iter()
        .map(|body| {
            let offset = fgumi_raw_bam::aux_data_offset_from_record(body)
                .filter(|&o| o <= body.len())
                .expect("record long enough to hold its aux offset");
            let mut aux: Vec<AuxEntry> = RawTagsView::new(&body[offset..])
                .iter()
                .map(|e| (e.tag, e.type_byte, e.value_bytes.to_vec()))
                .collect();
            if order == TagOrder::Sort {
                aux.sort();
            }
            NormRecord { core: body[..offset].to_vec(), aux }
        })
        .collect();
    NormBam { header, refs: bam.refs.clone(), records }
}

/// `version` without the git dev suffix bwa-mem3's `scripts/version.sh` appends
/// off the release tag: a trailing `-dirty`, then a trailing `-<7..=40 hex>`
/// commit.
pub fn without_git_suffix(version: &str) -> &str {
    let version = version.strip_suffix("-dirty").unwrap_or(version);
    match version.rsplit_once('-') {
        Some((base, sha))
            if (7..=40).contains(&sha.len()) && sha.bytes().all(|b| b.is_ascii_hexdigit()) =>
        {
            base
        }
        _ => version,
    }
}

/// Assert two normalized BAMs are identical, naming the first differing record
/// (by index and read name) rather than dumping both files.
pub fn assert_same_bam(expected: &NormBam, actual: &NormBam, label: &str) {
    assert_eq!(expected.header, actual.header, "{label}: headers differ");
    assert_eq!(expected.refs, actual.refs, "{label}: reference lists differ");
    let first_diff = expected.records.iter().zip(&actual.records).position(|(e, a)| e != a);
    if let Some(i) = first_diff {
        let (e, a) = (&expected.records[i], &actual.records[i]);
        panic!(
            "{label}: record {i} differs ('{}' vs '{}'):\n  expected {e:?}\n  actual   {a:?}",
            String::from_utf8_lossy(name(&e.core)),
            String::from_utf8_lossy(name(&a.core)),
        );
    }
    assert_eq!(
        expected.records.len(),
        actual.records.len(),
        "{label}: record counts differ (all shared records match)"
    );
}

// ---------------------------------------------------------------------------
// Fixture sanity + independent expectations
// ---------------------------------------------------------------------------

/// Counts of the alignment shapes the fixture is designed to produce.
#[derive(Debug, Default)]
pub struct Coverage {
    pub mapped: usize,
    pub unmapped: usize,
    pub supplementary: usize,
    pub xa_tagged: usize,
}

/// Count the alignment shapes in `bam`.
pub fn coverage(bam: &Bam) -> Coverage {
    let xa = SamTag::new(b'X', b'A');
    let mut c = Coverage::default();
    for body in &bam.records {
        let f = flags(body);
        if f & UNMAPPED == 0 {
            c.mapped += 1;
        } else {
            c.unmapped += 1;
        }
        if f & SUPPLEMENTARY != 0 {
            c.supplementary += 1;
        }
        if string_tag(body, xa).is_some() {
            c.xa_tagged += 1;
        }
    }
    c
}

/// Assert the output holds every alignment shape the fixture is built to
/// produce — mapped, unmapped, supplementary and `XA`-tagged records — so the
/// comparison can never silently degrade into comparing trivial output.
pub fn assert_fixture_coverage(bam: &Bam, label: &str) {
    let c = coverage(bam);
    assert!(
        c.mapped > 0 && c.unmapped > 0 && c.supplementary > 0 && c.xa_tagged > 0,
        "{label}: fixture no longer exercises every alignment shape: {c:?}"
    );
}

/// Assert the independent contract of the zipper merge against the input: every
/// template yields exactly one primary record per input read, and every primary
/// record carries its template's `RX` and, iff the template was QC-fail, the
/// QC-fail flag. This holds for both halves of a mid-pair split too.
pub fn assert_tags_transferred(bam: &Bam, templates: &[Template], label: &str) {
    let by_name: BTreeMap<&[u8], &Template> =
        templates.iter().map(|t| (t.name.as_bytes(), t)).collect();
    let mut primaries: BTreeMap<&[u8], usize> = BTreeMap::new();
    for body in bam.records.iter().filter(|b| is_primary(b)) {
        let read_name = name(body);
        let shown = String::from_utf8_lossy(read_name);
        let template = by_name.get(read_name).unwrap_or_else(|| {
            panic!("{label}: output record '{shown}' matches no input template")
        });
        assert_eq!(
            string_tag(body, SamTag::RX),
            Some(template.rx.as_bytes()),
            "{label}: primary record of '{shown}' lost (or changed) its RX tag"
        );
        assert_eq!(
            flags(body) & QC_FAIL != 0,
            template.qc_fail,
            "{label}: primary record of '{shown}' has the wrong QC-fail flag"
        );
        *primaries.entry(read_name).or_default() += 1;
    }
    for template in templates {
        assert_eq!(
            primaries.get(template.name.as_bytes()).copied().unwrap_or(0),
            template.reads.len(),
            "{label}: template '{}' should yield one primary record per input read",
            template.name
        );
    }
}

/// Assert every input read survives the merge as exactly one primary record of
/// its own template: the multiset of each template's primary `SEQ`s (reverse
/// complemented back for reverse-strand records) equals its input reads. A merge
/// that paired a mapped record with the wrong template, duplicated one read of a
/// pair, or dropped the other fails here even when the tags still match.
pub fn assert_reads_preserved(bam: &Bam, templates: &[Template], label: &str) {
    let mut seqs: BTreeMap<&[u8], Vec<Vec<u8>>> = BTreeMap::new();
    for body in bam.records.iter().filter(|b| is_primary(b)) {
        let mut seq = RawRecordView::new(body).sequence_vec();
        if flags(body) & REVERSE != 0 {
            seq = fgumi_dna::reverse_complement(&seq);
        }
        seqs.entry(name(body)).or_default().push(seq);
    }
    for template in templates {
        let mut actual = seqs.remove(template.name.as_bytes()).unwrap_or_default();
        actual.sort();
        let mut expected = template.reads.clone();
        expected.sort();
        assert!(
            actual == expected,
            "{label}: template '{}' primary records do not carry its input reads",
            template.name
        );
    }
}

/// Names of paired templates that bwa aligned as a mid-pair split: two or more
/// primary records under one name, none flagged paired.
pub fn split_names(bam: &Bam) -> BTreeSet<Vec<u8>> {
    let mut by_name: BTreeMap<&[u8], (usize, bool)> = BTreeMap::new();
    for body in bam.records.iter().filter(|b| is_primary(b)) {
        let entry = by_name.entry(name(body)).or_default();
        entry.0 += 1;
        entry.1 |= flags(body) & PAIRED != 0;
    }
    by_name
        .into_iter()
        .filter(|(_, (n, any_paired))| *n >= 2 && !any_paired)
        .map(|(n, _)| n.to_vec())
        .collect()
}
