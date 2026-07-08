//! Sort correctness tests for all sort orders and spill configurations.
//!
//! Validates that `fgumi sort` produces correctly ordered output across the
//! full order × spill matrix:
//!
//! - **Orders**: coordinate, queryname-lexicographic, queryname-natural,
//!   template-coordinate.
//! - **Spill regimes**: in-memory (large `-m`), single-spill, many-spill (small
//!   `-m`, 4 threads), and the `--temp-compression 0` uncompressed-bgzf spill
//!   path.
//!
//! ## Independent order oracle (S9b-002)
//!
//! Sort order is NOT proven solely by `fgumi sort --verify`: that path re-derives
//! the key with the SAME extractor functions and comparison operators the sort
//! path uses, so a wrong-but-consistent comparator would sort records into a
//! bogus order that its own verify then blesses. `--verify` is kept here only as
//! a write-corruption / post-write guard (truncation, reordering by the writer,
//! dropped records), demoted from the sole order proof.
//!
//! The actual order assertion is INDEPENDENT of the sort crate:
//! - coordinate and queryname-lexicographic orders are checked against an
//!   in-test expected order computed directly from the fixture records (the test
//!   constructs the records, so the correct `(tid, pos)` / lexicographic name
//!   order is computable without any sort-crate code);
//! - queryname-natural and template-coordinate orders — whose exact rules
//!   (`strnum_cmp` semantics; the packed template cell key) are error-prone to
//!   restate — are cross-checked against `samtools sort -n` /
//!   `samtools sort --template-coordinate`. Those samtools tests RUN in CI
//!   (samtools is installed in the workflow) and are gated only by a runtime
//!   `samtools_available()` check so local dev without samtools skips them
//!   gracefully with a message.
//!
//! ## Fixtures (S9b-001)
//!
//! Input BAMs are built in-process with `SamBuilder` + the noodles BAM writer
//! (no samtools, no committed data files), so the whole matrix runs by default
//! in CI — the spill/streaming-sort path is the branch's headline behaviour and
//! must execute on every PR.

use std::ffi::OsString;
use std::path::Path;
use std::process::Command;

use clap::Parser;
use fgumi_lib::commands::sort::Sort;
use fgumi_lib::sam::SamTag;
use fgumi_raw_bam::{RawRecord, SamBuilder, flags};
use noodles::bam;
use noodles::sam::Header;
use noodles::sam::alignment::io::Write as AlignmentWrite;
use noodles::sam::alignment::record_buf::RecordBuf;
use rstest::rstest;
use tempfile::TempDir;

use crate::helpers::bam_generator::to_record_buf;

// ---------------------------------------------------------------------------
// Fixture construction (in-process, no samtools)
// ---------------------------------------------------------------------------

/// Build a multi-chromosome unsorted SAM header with one read group.
fn unsorted_header() -> Header {
    use bstr::BString;
    use noodles::sam::header::record::value::map::Header as HeaderRecord;
    use noodles::sam::header::record::value::map::header::tag::Tag as HeaderTag;
    use noodles::sam::header::record::value::map::{Map, ReferenceSequence};
    use std::num::NonZeroUsize;

    let HeaderTag::Other(so_tag) = HeaderTag::from(*b"SO") else { unreachable!() };
    let header_map = Map::<HeaderRecord>::builder()
        .insert(so_tag, "unsorted")
        .build()
        .expect("valid header map");

    let mut builder = Header::builder().set_header(header_map);
    for chrom in ["chr1", "chr2", "chr3"] {
        builder = builder.add_reference_sequence(
            BString::from(chrom),
            Map::<ReferenceSequence>::new(NonZeroUsize::new(100_000).expect("non-zero ref len")),
        );
    }
    builder.build()
}

/// Generate a deliberately unsorted set of paired-end records spread across
/// three chromosomes with non-zero-padded read names.
///
/// The records are crafted so the four sort orders are all genuinely distinct:
/// - names are `read_0`, `read_1`, ... `read_N` (no zero padding) so natural and
///   lexicographic orders diverge (`read_2` < `read_10` naturally, but
///   `read_10` < `read_2` lexically);
/// - reads round-robin across chr2/chr1/chr3 and positions count *backwards*
///   within each chromosome, so input order matches no sort order;
/// - each template is assigned to one of four cells (`CB`) so the
///   template-coordinate cell-aware key path is exercised.
///
/// Returns the raw records in (unsorted) emission order.
fn unsorted_records(num_templates: usize) -> Vec<RawRecord> {
    let seq = b"ACGTACGTAC";
    let qual = vec![40u8; seq.len()];
    let cigar = u32::try_from(seq.len()).expect("seq len fits u32") << 4; // 10M (op 0)
    // Mate-cigar (MC) text — every read's mate is also a 10M alignment.
    // `samtools sort --template-coordinate` requires MC on paired-end inputs.
    let mc = format!("{}M", seq.len());

    let mut records = Vec::with_capacity(num_templates * 2 + 20);

    for i in 0..num_templates {
        // ref_id chosen so input order is NOT coordinate order.
        let ref_id: i32 = match i % 3 {
            0 => 1, // chr2
            1 => 0, // chr1
            _ => 2, // chr3
        };
        // Positions decrease as i grows within each residue class -> unsorted.
        let pos: i32 = 50_000 - i32::try_from((i % 1000) * 50).expect("pos offset fits i32");
        let mate_pos = pos + 200;
        let name = format!("read_{i}");
        let cell = format!("cell{}", i % 4);
        let umi = format!("umi_{i}");
        let template_len = 210_i32;

        // R1 (0x63 = PAIRED | PROPER_PAIR | MATE_REVERSE | FIRST_SEGMENT).
        let mut r1 = SamBuilder::new();
        r1.read_name(name.as_bytes())
            .ref_id(ref_id)
            .pos(pos)
            .mapq(60)
            .flags(flags::PAIRED | flags::PROPER_PAIR | flags::MATE_REVERSE | flags::FIRST_SEGMENT)
            .mate_ref_id(ref_id)
            .mate_pos(mate_pos)
            .template_length(template_len)
            .cigar_ops(&[cigar])
            .sequence(seq)
            .qualities(&qual)
            .add_string_tag(SamTag::RG, b"rg1")
            .add_string_tag(SamTag::MI, umi.as_bytes())
            .add_string_tag(SamTag::CB, cell.as_bytes())
            .add_string_tag(SamTag::MC, mc.as_bytes());
        records.push(r1.build());

        // R2 (PAIRED | PROPER_PAIR | REVERSE | LAST_SEGMENT).
        let mut r2 = SamBuilder::new();
        r2.read_name(name.as_bytes())
            .ref_id(ref_id)
            .pos(mate_pos)
            .mapq(60)
            .flags(flags::PAIRED | flags::PROPER_PAIR | flags::REVERSE | flags::LAST_SEGMENT)
            .mate_ref_id(ref_id)
            .mate_pos(pos)
            .template_length(-template_len)
            .cigar_ops(&[cigar])
            .sequence(seq)
            .qualities(&qual)
            .add_string_tag(SamTag::RG, b"rg1")
            .add_string_tag(SamTag::MI, umi.as_bytes())
            .add_string_tag(SamTag::CB, cell.as_bytes())
            .add_string_tag(SamTag::MC, mc.as_bytes());
        records.push(r2.build());
    }

    // A handful of fully-unmapped pairs (tid = -1) to exercise the
    // "no reference sorts last" branch of the coordinate / template key.
    for i in 0..10 {
        let name = format!("unmapped_{i}");
        for &seg in &[flags::FIRST_SEGMENT, flags::LAST_SEGMENT] {
            let mut b = SamBuilder::new();
            b.read_name(name.as_bytes())
                .flags(flags::PAIRED | flags::UNMAPPED | flags::MATE_UNMAPPED | seg)
                .sequence(seq)
                .qualities(&qual)
                .add_string_tag(SamTag::RG, b"rg1");
            records.push(b.build());
        }
    }

    records
}

/// Write `records` under `header` to a BAM at `path` using the noodles writer.
fn write_bam(path: &Path, header: &Header, records: &[RawRecord]) {
    let mut writer =
        bam::io::Writer::new(std::fs::File::create(path).expect("create input BAM file"));
    writer.write_header(header).expect("write header");
    for record in records {
        writer.write_alignment_record(header, &to_record_buf(record)).expect("write record");
    }
    writer.try_finish().expect("finish BAM");
}

/// Build an unsorted input BAM with `num_templates` paired templates and return
/// the temp dir (kept alive) plus the input path.
fn build_unsorted_fixture(num_templates: usize) -> (TempDir, std::path::PathBuf) {
    let dir = TempDir::new().expect("tempdir");
    let input = dir.path().join("unsorted.bam");
    let header = unsorted_header();
    let records = unsorted_records(num_templates);
    write_bam(&input, &header, &records);
    (dir, input)
}

// ---------------------------------------------------------------------------
// fgumi sort invocation
// ---------------------------------------------------------------------------

fn fgumi_sort_in_process(args: &[OsString]) -> anyhow::Result<()> {
    let cmd = Sort::try_parse_from(args.iter().cloned()).expect("failed to parse sort args");
    fgumi_lib::commands::sort::execute_sort_command(&cmd, "fgumi sort")
}

/// Spill regime applied to a sort run. Each variant maps to a concrete memory
/// limit + thread count + temp-compression configuration that forces the
/// targeted number of spills.
#[derive(Clone, Copy, Debug)]
enum Spill {
    /// Large memory budget — everything sorts in memory, no spill files.
    InMemory,
    /// Small memory budget — forces a small number of spill files.
    Single,
    /// Very small memory budget + 4 threads — forces many spill files.
    Many,
    /// Small memory budget with `--temp-compression 0 --temp-codec bgzf` —
    /// exercises the uncompressed-bgzf spill path under spill pressure.
    UncompressedBgzf,
}

impl Spill {
    fn max_memory(self) -> &'static str {
        match self {
            Spill::InMemory => "256M",
            Spill::Single | Spill::UncompressedBgzf => "200K",
            Spill::Many => "64K",
        }
    }

    fn threads(self) -> usize {
        match self {
            Spill::Many => 4,
            _ => 2,
        }
    }

    fn extra_args(self) -> Vec<&'static str> {
        match self {
            Spill::UncompressedBgzf => {
                vec!["--temp-compression", "0", "--temp-codec", "bgzf"]
            }
            _ => vec!["--temp-compression", "1"],
        }
    }
}

/// Sort `input` into `output` with `order` under the given spill regime.
fn run_sort(input: &Path, output: &Path, order: &str, spill: Spill) {
    let mut args: Vec<OsString> = vec![
        "sort".into(),
        "-i".into(),
        input.as_os_str().to_owned(),
        "-o".into(),
        output.as_os_str().to_owned(),
        "--order".into(),
        order.into(),
        "--threads".into(),
        spill.threads().to_string().into(),
        "-m".into(),
        spill.max_memory().into(),
    ];
    for a in spill.extra_args() {
        args.push(a.into());
    }
    fgumi_sort_in_process(&args)
        .unwrap_or_else(|e| panic!("fgumi sort failed (order={order}, spill={spill:?}): {e:#}"));
}

/// Run `fgumi sort --verify` as a write-corruption guard (NOT the order oracle).
fn verify_sorted(bam_path: &Path, order: &str) -> bool {
    let args: Vec<OsString> = vec![
        "sort".into(),
        "--verify".into(),
        "-i".into(),
        bam_path.as_os_str().to_owned(),
        "--order".into(),
        order.into(),
    ];
    fgumi_sort_in_process(&args).is_ok()
}

// ---------------------------------------------------------------------------
// Reading output back (noodles, no samtools)
// ---------------------------------------------------------------------------

/// Read every record of `path` as a `RecordBuf`, in file (sorted) order.
fn read_records(path: &Path) -> Vec<RecordBuf> {
    let mut reader =
        bam::io::Reader::new(std::fs::File::open(path).expect("open output BAM for reading"));
    let header = reader.read_header().expect("read header");
    reader
        .record_bufs(&header)
        .collect::<std::io::Result<Vec<_>>>()
        .expect("read records from output BAM")
}

/// A small, order-relevant projection of a record: `(tid, pos, name, flags)`.
/// `tid == -1` for records with no reference (BAM convention: sort last).
#[derive(Clone, PartialEq, Eq, PartialOrd, Ord, Debug)]
struct RecordKey {
    tid: i32,
    pos: i32,
    name: Vec<u8>,
    flags: u16,
}

fn record_key(r: &RecordBuf) -> RecordKey {
    let tid = r.reference_sequence_id().map_or(-1, |id| i32::try_from(id).expect("tid fits i32"));
    let pos = r
        .alignment_start()
        .map_or(-1, |p| i32::try_from(usize::from(p)).expect("pos fits i32") - 1);
    let name: Vec<u8> = r.name().map(|n| AsRef::<[u8]>::as_ref(n).to_vec()).unwrap_or_default();
    let flags = u16::from(r.flags());
    RecordKey { tid, pos, name, flags }
}

/// Assert exact record-count preservation and no content loss: the multiset of
/// records (by `RecordBuf` equality) must be identical between input and output.
/// Sorting reorders records but must not add, drop, or mutate any.
fn assert_records_preserved(input: &Path, output: &Path) {
    let mut in_recs = read_records(input);
    let mut out_recs = read_records(output);
    assert_eq!(
        in_recs.len(),
        out_recs.len(),
        "record count changed: input {} != output {}",
        in_recs.len(),
        out_recs.len()
    );
    // Sort both by a stable canonical key so the comparison is order-independent;
    // record content must round-trip the sort + (possibly uncompressed-bgzf)
    // write path byte-for-byte.
    let canon = |r: &RecordBuf| record_key(r);
    in_recs.sort_by_key(&canon);
    out_recs.sort_by_key(&canon);
    assert_eq!(in_recs, out_recs, "record content changed across sort (loss/dup/mutation)");
}

// ---------------------------------------------------------------------------
// Independent in-test order oracles (coordinate, queryname-lex)
// ---------------------------------------------------------------------------

/// Coordinate order oracle: `tid` ascending with no-reference (`tid == -1`) last,
/// then `pos` ascending. This is computed independently of the sort crate from
/// the records actually present in the output.
fn assert_coordinate_ordered(records: &[RecordBuf]) {
    // fgumi's coordinate sort key is `(tid << 34) | ((pos+1) << 1) | reverse`
    // (see `RawCoordinateKey`), so the REVERSE strand flag (0x10) is a real,
    // reliably-ordered tertiary component of the key — forward-strand records
    // precede reverse-strand at the same (tid, pos). Including it here catches a
    // comparator/merge regression that reshuffles equal-coordinate records by
    // strand, which a `(tid, pos)`-only projection would miss. The remaining
    // within-(tid, pos, strand) order is the stable input-order tie-break (no
    // name tie-break, matching samtools); it is not asserted here because the
    // parallel/spill merge does not independently guarantee global input order
    // for fully-equal keys, and record integrity is proven separately by
    // `assert_records_preserved`.
    const REVERSE_FLAG: u16 = 0x10;
    let strand = |k: &RecordKey| (k.flags & REVERSE_FLAG) != 0;
    let keys: Vec<RecordKey> = records.iter().map(record_key).collect();
    let mut expected = keys.clone();
    expected.sort_by(|a, b| {
        let a_noref = a.tid < 0;
        let b_noref = b.tid < 0;
        a_noref
            .cmp(&b_noref) // false (has ref) sorts before true (no ref)
            .then(a.tid.cmp(&b.tid))
            .then(a.pos.cmp(&b.pos))
            .then(strand(a).cmp(&strand(b))) // forward (false) before reverse (true)
    });
    let proj = |k: &RecordKey| (k.tid < 0, k.tid, k.pos, strand(k));
    let got: Vec<_> = keys.iter().map(proj).collect();
    let want: Vec<_> = expected.iter().map(proj).collect();
    assert_eq!(
        got, want,
        "output is not coordinate-ordered by (tid, pos, strand) (independent oracle)"
    );
}

/// Queryname-lexicographic order oracle. The lex sort key is the read name
/// compared as raw bytes; the sort is documented stable, so within an
/// equal-name run the records keep their input order — and the fixture emits
/// every template's R1 (`FIRST_SEGMENT`) before its R2 (`LAST_SEGMENT`). The order
/// is therefore FULLY determined, so this oracle uses an order-preserving
/// per-record projection `(name, is_last_segment)` rather than collapsing
/// duplicates by name: a name-only `<=` would pass even if R1/R2 were swapped
/// within a name, but the full projection catches per-record reordering within
/// an equal-name run. Computed independently of the sort crate.
fn assert_queryname_lex_ordered(records: &[RecordBuf]) {
    let proj = |r: &RecordBuf| -> (Vec<u8>, bool) {
        let name = record_key(r).name;
        // FIRST_SEGMENT (0x40) precedes LAST_SEGMENT (0x80) under the stable
        // sort; `is_last_segment` (false < true) encodes that secondary order.
        let is_last = r.flags().is_last_segment();
        (name, is_last)
    };
    let keys: Vec<(Vec<u8>, bool)> = records.iter().map(proj).collect();
    for pair in keys.windows(2) {
        assert!(
            pair[0] <= pair[1],
            "output is not queryname-lexicographic ordered (independent oracle): \
             ({:?}, last={}) > ({:?}, last={})",
            String::from_utf8_lossy(&pair[0].0),
            pair[0].1,
            String::from_utf8_lossy(&pair[1].0),
            pair[1].1,
        );
    }
}

// ---------------------------------------------------------------------------
// samtools cross-check oracle (queryname-natural, template-coordinate)
// ---------------------------------------------------------------------------

fn samtools_available() -> bool {
    which::which("samtools").is_ok()
}

/// Sort `input` with samtools into `output` using `samtools_args` (e.g.
/// `["-n"]` for queryname-natural, `["--template-coordinate"]` for
/// template-coordinate). Panics on failure.
fn samtools_sort(input: &Path, output: &Path, samtools_args: &[&str]) {
    let mut cmd = Command::new("samtools");
    cmd.arg("sort");
    cmd.args(samtools_args);
    cmd.arg("-o");
    cmd.arg(output);
    cmd.arg(input);
    let status = cmd.status().expect("run samtools sort");
    assert!(status.success(), "samtools sort failed");
}

/// Cross-check fgumi's output read-NAME order against samtools'. Used for
/// queryname-natural, where the read name IS the sort key, so the exact name
/// sequence is the contract being proven independently.
fn assert_name_order_matches(fgumi_out: &Path, samtools_out: &Path) {
    let fgumi_names: Vec<Vec<u8>> =
        read_records(fgumi_out).iter().map(|r| record_key(r).name).collect();
    let samtools_names: Vec<Vec<u8>> =
        read_records(samtools_out).iter().map(|r| record_key(r).name).collect();
    assert_eq!(fgumi_names.len(), samtools_names.len(), "fgumi / samtools record count differs");
    assert_eq!(
        fgumi_names, samtools_names,
        "fgumi read-name order disagrees with samtools (independent natural oracle)"
    );
}

/// Cross-check fgumi's template-coordinate POSITION-key sequence against
/// samtools'. The template-coordinate primary key is the template's
/// `(tid, pos)` coordinate; the per-record tie-breaker differs by design (fgumi
/// folds the `CB` cell tag + a name hash into its key, which samtools'
/// `--template-coordinate` does not model), so the read-NAME sequence is a
/// legitimately divergent tie order. The contract that IS shared — and the one
/// proven independently here — is that the ordered sequence of `(tid, pos)` keys
/// matches: both tools place every template at the same point in coordinate
/// order. Comparing the position-key stream is tie-agnostic and so isolates the
/// ordering definition from the unspecified tie order.
fn assert_template_position_order_matches(fgumi_out: &Path, samtools_out: &Path) {
    let proj = |path: &Path| -> Vec<(i32, i32)> {
        read_records(path)
            .iter()
            .map(|r| {
                let k = record_key(r);
                (k.tid, k.pos)
            })
            .collect()
    };
    let fgumi_keys = proj(fgumi_out);
    let samtools_keys = proj(samtools_out);
    assert_eq!(fgumi_keys.len(), samtools_keys.len(), "fgumi / samtools record count differs");
    assert_eq!(
        fgumi_keys, samtools_keys,
        "fgumi template-coordinate (tid,pos) order disagrees with samtools \
         (independent template oracle)"
    );
}

// ---------------------------------------------------------------------------
// Coordinate sort: full spill matrix, independent in-test oracle
// ---------------------------------------------------------------------------

#[rstest]
#[case::in_memory(Spill::InMemory)]
#[case::single_spill(Spill::Single)]
#[case::many_spill(Spill::Many)]
#[case::uncompressed_bgzf_spill(Spill::UncompressedBgzf)]
fn coordinate_sort_matrix(#[case] spill: Spill) {
    let (dir, input) = build_unsorted_fixture(1500);
    let output = dir.path().join("sorted.bam");
    run_sort(&input, &output, "coordinate", spill);

    // Write-corruption guard (NOT the order proof).
    assert!(verify_sorted(&output, "coordinate"), "coordinate --verify guard failed ({spill:?})");

    // Independent order proof + exact count / content preservation.
    let out_recs = read_records(&output);
    assert_coordinate_ordered(&out_recs);
    assert_records_preserved(&input, &output);
}

// ---------------------------------------------------------------------------
// Queryname-lexicographic sort: full spill matrix, independent in-test oracle
// ---------------------------------------------------------------------------

#[rstest]
#[case::in_memory(Spill::InMemory)]
#[case::single_spill(Spill::Single)]
#[case::many_spill(Spill::Many)]
#[case::uncompressed_bgzf_spill(Spill::UncompressedBgzf)]
fn queryname_lex_sort_matrix(#[case] spill: Spill) {
    let (dir, input) = build_unsorted_fixture(1500);
    let output = dir.path().join("sorted.bam");
    run_sort(&input, &output, "queryname", spill);

    assert!(verify_sorted(&output, "queryname"), "queryname-lex --verify guard failed ({spill:?})");

    let out_recs = read_records(&output);
    assert_queryname_lex_ordered(&out_recs);
    assert_records_preserved(&input, &output);
}

// ---------------------------------------------------------------------------
// Queryname-natural sort: full spill matrix, samtools cross-check oracle
// ---------------------------------------------------------------------------

#[rstest]
#[case::in_memory(Spill::InMemory)]
#[case::single_spill(Spill::Single)]
#[case::many_spill(Spill::Many)]
#[case::uncompressed_bgzf_spill(Spill::UncompressedBgzf)]
fn queryname_natural_sort_matrix(#[case] spill: Spill) {
    let (dir, input) = build_unsorted_fixture(1500);
    let output = dir.path().join("sorted.bam");
    run_sort(&input, &output, "queryname::natural", spill);

    // Corruption guard always runs.
    assert!(
        verify_sorted(&output, "queryname::natural"),
        "queryname-natural --verify guard failed ({spill:?})"
    );
    assert_records_preserved(&input, &output);

    // Independent order proof: cross-check the name sequence against samtools.
    if !samtools_available() {
        eprintln!("skipping queryname-natural samtools oracle: samtools not on PATH");
        return;
    }
    let samtools_out = dir.path().join("samtools_natural.bam");
    samtools_sort(&input, &samtools_out, &["-n"]);
    assert_name_order_matches(&output, &samtools_out);
}

// ---------------------------------------------------------------------------
// Template-coordinate sort: full spill matrix, samtools cross-check oracle
// ---------------------------------------------------------------------------

#[rstest]
#[case::in_memory(Spill::InMemory)]
#[case::single_spill(Spill::Single)]
#[case::many_spill(Spill::Many)]
#[case::uncompressed_bgzf_spill(Spill::UncompressedBgzf)]
fn template_coordinate_sort_matrix(#[case] spill: Spill) {
    let (dir, input) = build_unsorted_fixture(1500);
    let output = dir.path().join("sorted.bam");
    run_sort(&input, &output, "template-coordinate", spill);

    assert!(
        verify_sorted(&output, "template-coordinate"),
        "template-coordinate --verify guard failed ({spill:?})"
    );
    assert_records_preserved(&input, &output);

    // Run-to-run determinism: fgumi's template-coordinate tie-break folds the
    // CB cell tag + a read-name hash into its key and legitimately DIVERGES from
    // samtools' per-record tie order, so we do NOT cross-check the full record
    // stream against samtools (that would assert a false equivalence — see
    // `assert_template_position_order_matches`). But fgumi's OWN tie-break is
    // deterministic, so sorting the same input twice must yield the identical
    // full ordered record sequence. This pins the complete record order
    // (including ties) internally, which the position-key samtools cross-check
    // below cannot.
    let output2 = dir.path().join("sorted2.bam");
    run_sort(&input, &output2, "template-coordinate", spill);
    let recs1 = read_records(&output);
    let recs2 = read_records(&output2);
    assert_eq!(
        recs1, recs2,
        "template-coordinate sort is not run-to-run deterministic ({spill:?}): \
         the full ordered record sequence (including tie order) differs between two runs"
    );

    if !samtools_available() {
        eprintln!("skipping template-coordinate samtools oracle: samtools not on PATH");
        return;
    }
    let samtools_out = dir.path().join("samtools_template.bam");
    samtools_sort(&input, &samtools_out, &["--template-coordinate"]);
    // Cross-check only the (tid, pos) position-key stream against samtools — the
    // per-record tie order is intentionally divergent (documented above).
    assert_template_position_order_matches(&output, &samtools_out);
}

// ---------------------------------------------------------------------------
// Template-coordinate: mate-side / strand key coverage (independent oracle)
// ---------------------------------------------------------------------------
//
// The shared `unsorted_records` fixture gives every mapped template the SAME
// orientation (R1 forward, R2 reverse) and a FIXED mate offset (+200), so its
// template-coordinate key collapses to ordinary `(tid, pos)` coordinate order:
// the mate-side (`tid2`/`pos2`) and per-end strand (`neg1`/`neg2`) lanes of the
// samtools-compatible key never independently vary, and a bug confined to them
// would still pass `template_coordinate_sort_matrix`. This fixture instead pins
// several templates to the SAME lower-end `(tid1, pos1)` while varying the
// upper-end chromosome, position, and both ends' strands, so the secondary key
// lanes are the sole discriminator of the correct order.

/// Reference span of every read in the mate-geometry fixture (all `10M`).
const MATE_GEOM_SPAN: i32 = 10;

/// A template-coordinate-orderable projection of an output record, mirroring the
/// canonical key normalization in `fgumi_sort`'s `extract_template_key_inline`:
/// the lower end (by `(tid, unclipped-5')`, ties broken so the reverse end is
/// "upper") becomes `(tid1, pos1, neg1)` and the mate becomes `(tid2, pos2,
/// neg2)`. The returned tuple orders fields exactly as the key compares them —
/// `tid1, tid2, pos1, pos2`, the strands, then the final `is_upper` lane — with
/// each strand encoded as `u8::from(!neg)` so the reverse strand (`neg = true`)
/// sorts first, matching samtools.
///
/// The trailing `is_upper` byte (`0` = lower-of-pair, `1` = upper-of-pair) pins
/// the documented within-template tie-break: both records of a template share
/// the leading six lanes, so without it a deterministic regression that flips
/// the per-record order inside a template would still satisfy the oracle. The
/// production key's penultimate lane is the read-name *hash* (not the literal
/// name, and not lexicographic), which is non-reproducible here; it only decides
/// order between templates that tie on every position/strand lane, and this
/// fixture gives every template a unique position key, so that lane never
/// arbitrates and is intentionally not modeled. Derived purely from record
/// fields (own + mate, plus the `MC` span), independent of the sort crate.
fn template_coord_proj(r: &RecordBuf) -> (i32, i32, i32, i32, u8, u8, u8) {
    let zero_based = |p: noodles::core::Position| i32::try_from(usize::from(p)).expect("pos") - 1;
    let tid = r.reference_sequence_id().map_or(-1, |id| i32::try_from(id).expect("tid fits i32"));
    let pos = zero_based(r.alignment_start().expect("mapped record has alignment start"));
    let is_rev = r.flags().is_reverse_complemented();
    let mate_tid =
        r.mate_reference_sequence_id().map_or(-1, |id| i32::try_from(id).expect("mate tid"));
    let mate_pos = zero_based(r.mate_alignment_start().expect("paired record has mate start"));
    let mate_rev = r.flags().is_mate_reverse_complemented();

    // Unclipped 5' position: alignment start for forward reads, alignment end
    // (start + span - 1) for reverse reads. Every fixture read is a clip-free
    // `10M` alignment, so the span is constant.
    let five_prime = |p: i32, rev: bool| if rev { p + MATE_GEOM_SPAN - 1 } else { p };
    let this5 = five_prime(pos, is_rev);
    let mate5 = five_prime(mate_pos, mate_rev);

    // Canonical normalization (extract_template_key_inline): this read is the
    // "upper" end when its coordinate exceeds the mate's, or ties and it is the
    // reverse strand.
    let is_upper =
        (tid, this5) > (mate_tid, mate5) || ((tid, this5) == (mate_tid, mate5) && is_rev);
    let (tid1, pos1, neg1, tid2, pos2, neg2) = if is_upper {
        (mate_tid, mate5, mate_rev, tid, this5, is_rev)
    } else {
        (tid, this5, is_rev, mate_tid, mate5, mate_rev)
    };
    (tid1, tid2, pos1, pos2, u8::from(!neg1), u8::from(!neg2), u8::from(is_upper))
}

/// Independent template-coordinate order oracle: assert the output record stream
/// is non-decreasing under the full template key projection — the mate-side and
/// strand lanes AND the trailing within-template `is_upper` tie-break. Computed
/// from the records actually present in the output, independently of the sort
/// crate.
fn assert_template_coordinate_ordered(records: &[RecordBuf]) {
    let projs: Vec<_> = records.iter().map(template_coord_proj).collect();
    for pair in projs.windows(2) {
        assert!(
            pair[0] <= pair[1],
            "output is not template-coordinate ordered by the full \
             (tid1, tid2, pos1, pos2, neg1, neg2, is_upper) key (independent oracle): \
             {:?} > {:?}",
            pair[0],
            pair[1],
        );
    }
}

/// Build a paired-end template `name` from its two ends `a` and `b`, each given
/// as `(tid, pos0, reverse)` (0-based `pos`), optionally attaching a cell-barcode
/// (`CB`) tag. Both reads are clip-free `10M` alignments carrying the mate's
/// coordinate/strand and an `MC` tag, so the sort can compute the mate's
/// unclipped 5' position from each record alone.
fn push_template_with_cb(
    records: &mut Vec<RawRecord>,
    name: &str,
    cb: Option<&[u8]>,
    a: (i32, i32, bool),
    b: (i32, i32, bool),
) {
    let seq = b"ACGTACGTAC"; // 10 bases == 10M
    let qual = vec![40u8; seq.len()];
    let cigar = u32::try_from(seq.len()).expect("len fits u32") << 4; // 10M

    let build = |seg_flag: u16, this: (i32, i32, bool), mate: (i32, i32, bool)| {
        let (tid, pos, rev) = this;
        let (mate_tid, mate_pos, mate_rev) = mate;
        let mut f = flags::PAIRED | flags::PROPER_PAIR | seg_flag;
        if rev {
            f |= flags::REVERSE;
        }
        if mate_rev {
            f |= flags::MATE_REVERSE;
        }
        let mut sb = SamBuilder::new();
        sb.read_name(name.as_bytes())
            .ref_id(tid)
            .pos(pos)
            .mapq(60)
            .flags(f)
            .mate_ref_id(mate_tid)
            .mate_pos(mate_pos)
            .template_length(0)
            .cigar_ops(&[cigar])
            .sequence(seq)
            .qualities(&qual)
            .add_string_tag(SamTag::RG, b"rg1")
            .add_string_tag(SamTag::MC, b"10M");
        if let Some(cb) = cb {
            sb.add_string_tag(SamTag::CB, cb);
        }
        sb.build()
    };

    records.push(build(flags::FIRST_SEGMENT, a, b));
    records.push(build(flags::LAST_SEGMENT, b, a));
}

/// Build a paired-end template `name` (no `CB` tag) from its two ends.
fn push_mate_geom_template(
    records: &mut Vec<RawRecord>,
    name: &str,
    a: (i32, i32, bool),
    b: (i32, i32, bool),
) {
    push_template_with_cb(records, name, None, a, b);
}

/// Templates that share a lower-end `(tid1, pos1)` but diverge in the mate-side
/// and strand key lanes, emitted out of template-coordinate order. The correct
/// order is `t_low_rev, t_pos2_lo, t_up_rev, t_pos2_hi, t_tie, t_tid2, t_tid1`;
/// the records are written scrambled so a real sort has work to do and a
/// primary-`(tid1, pos1)`-only key would leave the shared-lower-end group
/// misordered.
fn mate_geometry_records() -> Vec<RawRecord> {
    let mut records = Vec::new();
    // Scrambled emission order. Each end is (tid, pos0, reverse); all five
    // chr1-lower templates share lower end (tid1=0, pos1=1000) except where a
    // reverse lower end is used to vary neg1 at the same 5' coordinate.
    push_mate_geom_template(&mut records, "t_pos2_hi", (0, 1000, false), (0, 3000, false));
    push_mate_geom_template(&mut records, "t_tid2", (0, 1000, false), (1, 500, false));
    // endA forward@5000, endB reverse@4991 → both 5'=5000; reverse end is upper.
    push_mate_geom_template(&mut records, "t_tie", (0, 5000, false), (0, 4991, true));
    push_mate_geom_template(&mut records, "t_up_rev", (0, 1000, false), (0, 2000, true));
    push_mate_geom_template(&mut records, "t_tid1", (1, 100, false), (1, 200, false));
    push_mate_geom_template(&mut records, "t_pos2_lo", (0, 1000, false), (0, 2000, false));
    // reverse lower end at pos 991 → 5'=1000, same pos1 as the forward-lower
    // templates but neg1 differs, so it must sort before t_pos2_lo.
    push_mate_geom_template(&mut records, "t_low_rev", (0, 991, true), (0, 2000, false));
    records
}

// Runs in-memory only: this fixture is a handful of templates crafted to make
// the mate-side/strand key lanes the sole order discriminator, so it would never
// fill the spill thresholds anyway. The comparator is identical on the in-memory
// and spill/merge paths, and the spill/merge path itself is exercised across the
// `64K`/`200K` budgets by the 1500-template `*_sort_matrix` tests above — so a
// spill `#[case]` matrix here would only re-run the same in-memory sorter and
// falsely advertise merge-path coverage it does not provide.
#[test]
fn template_coordinate_mate_geometry_oracle() {
    let dir = TempDir::new().expect("tempdir");
    let input = dir.path().join("unsorted.bam");
    let records = mate_geometry_records();
    write_bam(&input, &unsorted_header(), &records);

    let output = dir.path().join("sorted.bam");
    run_sort(&input, &output, "template-coordinate", Spill::InMemory);

    assert!(
        verify_sorted(&output, "template-coordinate"),
        "template-coordinate --verify guard failed"
    );
    assert_records_preserved(&input, &output);

    // The real contract: the output respects the FULL template key, including
    // the mate-side and strand lanes that the shared-fixture matrix cannot
    // exercise. With several templates pinned to the same lower-end (tid1, pos1),
    // a comparator that dropped pos2/neg1/neg2 would leave that group in input
    // (scrambled) order and break this assertion.
    let out = read_records(&output);
    assert_template_coordinate_ordered(&out);

    // Sanity: the fixture genuinely stresses the secondary lanes — more than one
    // template shares the lower-end (tid1, pos1), so ordering them requires the
    // mate-side/strand key (not just the primary coordinate).
    let mut primary_keys: Vec<(i32, i32)> =
        out.iter().map(template_coord_proj).map(|p| (p.0, p.2)).collect();
    primary_keys.sort_unstable();
    primary_keys.dedup();
    let shared_lower_end = out.len() / 2 > primary_keys.len();
    assert!(
        shared_lower_end,
        "fixture no longer shares a lower-end (tid1, pos1) across templates — \
         the mate-side/strand lanes are no longer exercised"
    );
}

// ---------------------------------------------------------------------------
// Template-coordinate: final CB / read-name-hash tie-break lanes
// ---------------------------------------------------------------------------
//
// The mate-geometry oracle and the samtools `(tid,pos)` cross-check both stop at
// the position/strand lanes, and every template in those fixtures has a unique
// position key — so the production key's trailing CB lane (`cb_hash`) and
// read-name-hash lane never arbitrate, and a regression there (dropped lane,
// reversed comparison, wrong key position) would slip through. This test builds
// template *pairs* that collide on `(tid1, tid2, pos1, pos2, neg1, neg2)` and
// differ in EXACTLY ONE trailing lane, then pins fgumi's order against a FROZEN
// expected order (below). The expectation is a precomputed constant, NOT
// recomputed from `fgumi_sort::cb_hasher` / `LibraryLookup::hash_name` at
// runtime: deriving it from the same hashers the sorter uses would let a bad
// seed or hash-wiring regression move both the output and the expectation
// together and pass silently. Freezing makes such a regression flip the actual
// order against a fixed expectation and fail. The frozen values were captured
// from the fixed-seed hashers for these exact inputs; a deliberate seed change
// is itself a sort-order change and is expected to require updating them.
//
// In-memory only (like `template_coordinate_mate_geometry_oracle`): four
// templates never reach the spill thresholds, the comparator is identical on the
// merge path, and spill/merge is covered by the `*_sort_matrix` tests.
#[test]
fn template_coordinate_cb_and_name_hash_tiebreak() {
    let dir = TempDir::new().expect("tempdir");
    let input = dir.path().join("unsorted.bam");

    // Two colliding pairs at distinct positions so each pair's assertion is
    // independent of the other:
    //  - CB pair  @ chr1:1000(F)/2000(F): identical but for the CB tag
    //    (`cellA`/`cellB`) → the `cb_hash` lane is the sole discriminator.
    //  - name pair @ chr1:5000(F)/6000(F): identical, NO CB (so `cb_hash` is 0
    //    for both), differing only in read name → the name-hash lane decides.
    let mut records = Vec::new();
    push_template_with_cb(&mut records, "cb_a", Some(b"cellA"), (0, 1000, false), (0, 2000, false));
    push_template_with_cb(&mut records, "cb_b", Some(b"cellB"), (0, 1000, false), (0, 2000, false));
    push_template_with_cb(&mut records, "name_aaa", None, (0, 5000, false), (0, 6000, false));
    push_template_with_cb(&mut records, "name_zzz", None, (0, 5000, false), (0, 6000, false));
    write_bam(&input, &unsorted_header(), &records);

    let output = dir.path().join("sorted.bam");
    run_sort(&input, &output, "template-coordinate", Spill::InMemory);
    assert!(verify_sorted(&output, "template-coordinate"), "verify guard failed");
    assert_records_preserved(&input, &output);

    let out = read_records(&output);
    // The ordered `(name, is_first_segment)` slice of just the records whose
    // name is in `names`, in output order. Each colliding template emits two
    // records sharing a name: the FIRST_SEGMENT end is the lower 5' mate
    // (`is_upper = 0`), the LAST_SEGMENT end the upper (`is_upper = 1`). Because
    // the `cb_hash` / `name_hash` lane is MORE significant than the trailing
    // `is_upper` lane, BOTH records of the earlier-hashing template must precede
    // BOTH of the later one — never interleaved. Asserting only which template's
    // first record appears earlier (a bare `first_index` check) would miss an
    // interleaving that still violates the frozen lane order, so pin the full
    // four-record slice.
    let pair_slice = |names: &[&[u8]]| -> Vec<(Vec<u8>, bool)> {
        out.iter()
            .filter(|r| {
                r.name().is_some_and(|n| {
                    let nb = AsRef::<[u8]>::as_ref(n);
                    names.contains(&nb)
                })
            })
            .map(|r| {
                let name = AsRef::<[u8]>::as_ref(&r.name().expect("named record")).to_vec();
                (name, r.flags().is_first_segment())
            })
            .collect()
    };

    // FROZEN expected order (see the header comment for why it is not recomputed
    // from the production hashers). For the fixed-seed `cb_hasher`, `cellA`'s
    // hash sorts before `cellB`'s, so the CB pair comes out cb_a then cb_b. For
    // `LibraryLookup`'s fixed-seed name hasher, `name_zzz`'s hash sorts before
    // `name_aaa`'s, so the name pair comes out name_zzz then name_aaa. Within
    // each template the lower-end (FIRST_SEGMENT) record precedes the upper-end
    // (LAST_SEGMENT) one.
    assert_eq!(
        pair_slice(&[b"cb_a", b"cb_b"]),
        vec![
            (b"cb_a".to_vec(), true),
            (b"cb_a".to_vec(), false),
            (b"cb_b".to_vec(), true),
            (b"cb_b".to_vec(), false),
        ],
        "CB tie-break regressed: expected both records of the cellA template (cb_a) to \
         sort before both of the cellB template (cb_b), lower end first, under the frozen \
         cb_hash order",
    );
    assert_eq!(
        pair_slice(&[b"name_zzz", b"name_aaa"]),
        vec![
            (b"name_zzz".to_vec(), true),
            (b"name_zzz".to_vec(), false),
            (b"name_aaa".to_vec(), true),
            (b"name_aaa".to_vec(), false),
        ],
        "read-name tie-break regressed: expected both records of name_zzz to sort before \
         both of name_aaa, lower end first, under the frozen name_hash order",
    );

    // The tie-break order must also be run-to-run deterministic (the fixed seeds
    // guarantee it; a regression to a per-run seed would break this).
    let output2 = dir.path().join("sorted2.bam");
    run_sort(&input, &output2, "template-coordinate", Spill::InMemory);
    assert_eq!(
        out,
        read_records(&output2),
        "template-coordinate CB/name tie-break is not run-to-run deterministic",
    );
}

// ---------------------------------------------------------------------------
// Thread-count consistency (coordinate, independent oracle)
// ---------------------------------------------------------------------------

/// All thread counts must produce the *same* coordinate-ordered output that
/// preserves every record. The same input is sorted at 1, 2, and 4 threads;
/// each run is self-checked (sortedness + record preservation), and then the
/// order-determining `(tid, pos)` key stream is asserted byte-identical across
/// all thread counts. Self-checking each `#[case]` independently would let a
/// thread-count-dependent ordering regression slip through, so the cross-run
/// comparison is the real contract here. A small `-m` forces spills so the
/// external-merge path is exercised at each thread count.
#[test]
fn coordinate_sort_consistent_across_threads() {
    let (dir, input) = build_unsorted_fixture(1000);

    // The order-determining projection: tid ascending (no-reference `tid == -1`
    // last), then pos. Names within an equal `(tid, pos)` are an unspecified tie
    // (see `assert_coordinate_ordered`), so they are excluded from the
    // cross-thread identity check — only the documented key stream is compared.
    let position_stream = |output: &Path| -> Vec<(bool, i32, i32)> {
        read_records(output).iter().map(record_key).map(|k| (k.tid < 0, k.tid, k.pos)).collect()
    };

    let mut baseline_stream: Option<Vec<(bool, i32, i32)>> = None;
    for threads in [1usize, 2, 4] {
        let output = dir.path().join(format!("sorted_t{threads}.bam"));
        let args: Vec<OsString> = vec![
            "sort".into(),
            "-i".into(),
            input.as_os_str().to_owned(),
            "-o".into(),
            output.as_os_str().to_owned(),
            "--order".into(),
            "coordinate".into(),
            "--threads".into(),
            threads.to_string().into(),
            "-m".into(),
            "200K".into(),
            "--temp-compression".into(),
            "1".into(),
        ];
        fgumi_sort_in_process(&args)
            .unwrap_or_else(|e| panic!("fgumi sort failed (threads={threads}): {e:#}"));

        assert!(
            verify_sorted(&output, "coordinate"),
            "coordinate guard failed (threads={threads})"
        );
        assert_coordinate_ordered(&read_records(&output));
        assert_records_preserved(&input, &output);

        // Cross-thread identity: every thread count must emit the *same*
        // coordinate key stream, not merely a self-consistent one.
        let stream = position_stream(&output);
        match &baseline_stream {
            None => baseline_stream = Some(stream),
            Some(baseline) => assert_eq!(
                &stream, baseline,
                "coordinate key stream differs at threads={threads} vs the 1-thread baseline",
            ),
        }
    }
}
