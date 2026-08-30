//! Smoke test for the sort-command cutover (Task 3): `Sort::execute` runs its
//! coordinate sort through the declarative chain builder (`build_for(spec)?.run()`)
//! instead of the owned `RawExternalSorter` engine, and still produces
//! correctly sorted output.
//!
//! **Not a genuine RED/GREEN gate.** Pre-cutover, `execute_sort`'s owned-engine
//! path already sorts correctly, so this test passes before Task 3's change
//! too — there is no way to observe a real regression here that a broken
//! cutover wouldn't also need to break record count or ordering outright. It
//! exists as a basic end-to-end guard through the cutover, kept green on both
//! sides. The real parity gate — the chain's output is byte-identical to the
//! owned engine's, across sort orders and spill configurations — is Task 4's
//! job, not this test's.

use std::ffi::OsStr;
use std::fs;
use std::path::Path;
use std::process::Command;

use clap::Parser;
use fgumi_lib::commands::command::Command as FgumiCommand;
use fgumi_lib::commands::sort::Sort;
use fgumi_lib::sam::SamTag;
use fgumi_raw_bam::{RawRecord, SamBuilder, flags};
use rstest::rstest;
use tempfile::TempDir;

use crate::helpers::bam_generator::{create_minimal_header, write_bam};
use crate::helpers::read_bam_output;

/// `n` mapped, single-end records on one reference, placed at *descending*
/// positions so the input is deliberately unsorted — the sort stage has real
/// work to do before this test's assertions can distinguish a working
/// cutover from a no-op one.
fn unsorted_records(n: usize) -> Vec<RawRecord> {
    (0..n)
        .map(|i| {
            let pos = i32::try_from((n - i) * 100).expect("pos fits i32");
            let mut b = SamBuilder::new();
            b.read_name(format!("read{i}").as_bytes())
                .ref_id(0)
                .pos(pos)
                .mapq(60)
                .flags(0)
                .cigar_ops(&[4u32 << 4]) // 4M
                .sequence(b"ACGT")
                .qualities(&[30u8; 4]);
            b.build()
        })
        .collect()
}

#[test]
fn sort_command_produces_coordinate_sorted_output_via_chain() {
    let dir = TempDir::new().expect("create temp dir");
    let input_bam = dir.path().join("in.bam");
    let output_bam = dir.path().join("out.bam");

    let header = create_minimal_header("chr1", 1_000_000);
    let records = unsorted_records(500);
    let input_count = records.len();
    write_bam(&input_bam, &header, &records);

    // Pass `&OsStr` directly (not `&str`) so `try_parse_from` doesn't UTF-8
    // round-trip the temp-dir paths, matching the convention in
    // `test_sort_write_index.rs`.
    let cmd = Sort::try_parse_from([
        OsStr::new("sort"),
        OsStr::new("-i"),
        input_bam.as_os_str(),
        OsStr::new("-o"),
        output_bam.as_os_str(),
        OsStr::new("--order"),
        OsStr::new("coordinate"),
    ])
    .expect("failed to parse sort args");

    cmd.execute("fgumi sort").expect("sort command should succeed via the chain");

    let (_, out_records) = read_bam_output(&output_bam);
    assert_eq!(out_records.len(), input_count, "output record count must match input record count");

    let keys: Vec<(usize, usize)> = out_records
        .iter()
        .map(|r| {
            (
                r.reference_sequence_id().unwrap_or(usize::MAX),
                r.alignment_start().map_or(0, usize::from),
            )
        })
        .collect();
    assert!(
        keys.windows(2).all(|w| w[0] <= w[1]),
        "output is not coordinate-sorted (non-decreasing (ref_id, pos) expected): {keys:?}"
    );
}

// ============================================================================
// Task 4: four-order parity vs. the saved baseline binary + samtools
//
// This is the genuine correctness gate for the cutover: the chain builder's
// sort output must be byte-identical to the owned-engine baseline binary
// (modulo the `@PG` header line, which necessarily differs -- different
// `VN` from git-describe, different `CL` naming a different argv[0] and
// output path) across all four sort orders, and order-valid with the same
// record multiset as samtools where samtools has an equivalent order.
//
// A failure here is a real cutover parity bug, not a test to relax.
// ============================================================================

/// Resolves the saved owned-engine baseline binary to compare against.
///
/// Precedence: `FGUMI_BASELINE_BIN` env var, then the known checked-in-CI-host
/// path, then `None` (meaning: skip the baseline half of the gate). A `None`
/// return is only ever a "no oracle available" signal, not a passing result --
/// callers must still require the samtools half (where one exists) or skip the
/// whole case with an explicit `eprintln!`.
fn baseline_bin() -> Option<std::path::PathBuf> {
    if let Ok(path) = std::env::var("FGUMI_BASELINE_BIN") {
        let path = std::path::PathBuf::from(path);
        if path.is_file() {
            return Some(path);
        }
        eprintln!(
            "FGUMI_BASELINE_BIN={} does not name an existing file; baseline oracle unavailable",
            path.display()
        );
        return None;
    }

    let fallback =
        std::path::PathBuf::from("/Users/nhomer/work/git/fgumi/baselines/fgumi-owned-07b4a51b");
    fallback.is_file().then_some(fallback)
}

/// Whether `samtools` is on `PATH` and runnable.
fn samtools_available() -> bool {
    Command::new("samtools").arg("--version").output().map(|o| o.status.success()).unwrap_or(false)
}

/// Runs `<bin> sort -i <input> -o <output> --order <order_flag>` and asserts success.
fn run_fgumi_sort(bin: &Path, input: &Path, output: &Path, order_flag: &str) {
    let status = Command::new(bin)
        .args([
            OsStr::new("sort"),
            OsStr::new("-i"),
            input.as_os_str(),
            OsStr::new("-o"),
            output.as_os_str(),
            OsStr::new("--order"),
            OsStr::new(order_flag),
        ])
        .status()
        .unwrap_or_else(|e| panic!("failed to spawn `{}`: {e}", bin.display()));
    assert!(status.success(), "`{}` sort --order {order_flag} failed", bin.display());
}

/// Runs `<bin> sort --verify -i <path> --order <order_flag>` and returns whether
/// it exits 0 (the file is correctly ordered).
fn fgumi_verify_sorted(bin: &Path, path: &Path, order_flag: &str) -> bool {
    Command::new(bin)
        .args([
            OsStr::new("sort"),
            OsStr::new("--verify"),
            OsStr::new("-i"),
            path.as_os_str(),
            OsStr::new("--order"),
            OsStr::new(order_flag),
        ])
        .status()
        .map(|s| s.success())
        .unwrap_or(false)
}

/// Runs `samtools <args> -o <output> <input>` and asserts success.
fn run_samtools_sort(args: &[&str], input: &Path, output: &Path) {
    let status = Command::new("samtools")
        .args(args)
        .arg("-o")
        .arg(output)
        .arg(input)
        .status()
        .expect("failed to spawn samtools");
    assert!(status.success(), "`samtools {args:?}` failed for {}", input.display());
}

/// Collects a sorted multiset of per-record identities (`QNAME`, `FLAG`,
/// `RNAME`, `POS`) via `samtools view`.
///
/// Mirrors `record_identities` in `test_sort_correctness.rs`: sorting reorders
/// records and fgumi/samtools break ties differently, so identity has to be
/// compared order-independently rather than position-by-position.
fn record_identities(bam_path: &Path) -> Vec<String> {
    let output = Command::new("samtools")
        .args(["view", bam_path.to_str().expect("path is valid UTF-8")])
        .output()
        .expect("samtools view");
    assert!(output.status.success(), "samtools view failed for {}", bam_path.display());
    let mut identities: Vec<String> = String::from_utf8_lossy(&output.stdout)
        .lines()
        .map(|line| {
            let mut fields = line.split('\t');
            let qname = fields.next().unwrap_or_default();
            let flag = fields.next().unwrap_or_default();
            let rname = fields.next().unwrap_or_default();
            let pos = fields.next().unwrap_or_default();
            format!("{qname}\t{flag}\t{rname}\t{pos}")
        })
        .collect();
    identities.sort_unstable();
    identities
}

/// Removes every `@PG` line from a SAM header text blob.
///
/// Both the cutover build and the baseline binary stamp a single `@PG` line
/// whose `VN` (git-describe version) and `CL` (command line, naming argv[0] --
/// a different binary path for each side -- and the per-run output path)
/// necessarily differ between two independent invocations. Stripping the whole
/// line (not just normalizing those two sub-fields) is simplest and correct
/// here because neither side emits more than one `@PG` record for this input
/// (the test BAMs carry no pre-existing `@PG`).
fn strip_pg_lines(text: &str) -> String {
    let mut lines: Vec<&str> = text.split('\n').collect();
    let had_trailing_newline = lines.last() == Some(&"");
    if had_trailing_newline {
        lines.pop();
    }
    lines.retain(|line| !line.starts_with("@PG"));
    let mut out = lines.join("\n");
    if !out.is_empty() || had_trailing_newline {
        out.push('\n');
    }
    out
}

/// Reads a BAM file's raw BGZF stream, decompresses it, strips the `@PG`
/// line(s) from the embedded SAM header text, and returns the resulting bytes
/// (new header + unmodified `n_ref`/reference-list/record bytes).
///
/// This compares the two writers' *decompressed* output rather than
/// `RecordBuf`-parsed records or raw file bytes: BGZF block boundaries and
/// compression levels differ between the two writers even when the logical
/// content is identical, so raw file bytes never match; parsing into
/// `RecordBuf` and re-encoding risks masking a real tag-order or
/// binary-layout regression by normalizing it away. Operating directly on the
/// decompressed BAM binary format keeps everything except the `@PG` line an
/// exact, uninterpreted byte comparison.
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

/// Builds a diverse SAM header for the four-order parity gate: three
/// reference sequences (so coordinate sort has real (ref, pos) work to do)
/// and two read groups on two libraries (so template-coordinate's
/// library-lookup lane is exercised, not just its fallback).
fn diverse_header() -> noodles::sam::Header {
    use bstr::BString;
    use noodles::sam::header::record::value::map::header::tag::Tag as HeaderTag;
    use noodles::sam::header::record::value::map::read_group::tag as rg_tag;
    use noodles::sam::header::record::value::map::{
        Header as HeaderRecord, Map, ReadGroup, ReferenceSequence,
    };
    use std::num::NonZeroUsize;

    let mut header_builder = Map::<HeaderRecord>::builder();
    for &(tag_bytes, value) in
        &[(*b"SO", "unsorted"), (*b"GO", "query"), (*b"SS", "unsorted:template-coordinate")]
    {
        let HeaderTag::Other(tag) = HeaderTag::from(tag_bytes) else { unreachable!() };
        header_builder = header_builder.insert(tag, value);
    }
    let header_map = header_builder.build().expect("valid header map");

    let mut builder = noodles::sam::Header::builder().set_header(header_map);
    for ref_name in ["chr1", "chr2", "chr3"] {
        builder = builder.add_reference_sequence(
            BString::from(ref_name),
            Map::<ReferenceSequence>::new(NonZeroUsize::new(500_000).expect("non-zero")),
        );
    }
    for (id, library) in [("RG1", "libA"), ("RG2", "libB")] {
        let rg = Map::<ReadGroup>::builder()
            .insert(rg_tag::LIBRARY, String::from(library))
            .build()
            .expect("read group builds");
        builder = builder.add_read_group(BString::from(id), rg);
    }
    builder.build()
}

/// Builds `num_pairs` unsorted, paired-end mapped read pairs plus a handful of
/// unmapped pairs, spread across three references, two read groups/libraries,
/// four cell barcodes, and non-zero-padded numeric read names.
///
/// The non-zero-padded names (`read2`, `read10`, ...) are what makes
/// lexicographic and natural queryname order diverge (`"read10" < "read2"`
/// lexicographically, but not numerically) -- without that, the
/// queryname/queryname-natural cases could pass by coincidence even if the
/// comparators were swapped. Positions descend as the pair index increases
/// (mirroring `test_sort_correctness.rs`'s generator) so the input is
/// deliberately unsorted for every order under test, not just coordinate.
fn diverse_unsorted_records(num_pairs: usize) -> Vec<RawRecord> {
    let seq = b"ACGTACGTAC";
    let read_len = u32::try_from(seq.len()).expect("seq len fits u32");
    let cigar_op = read_len << 4; // 10M
    let qual = vec![30u8; seq.len()];

    let mut records = Vec::with_capacity(num_pairs * 2 + 16);
    for i in 0..num_pairs {
        let ref_id = i32::try_from(i % 3).expect("fits i32");
        // Positions cycle downward within each reference, so within-ref order
        // is scrambled rather than merely offset by a constant.
        let pos = i32::try_from(400_000 - (i % 4000) * 100).expect("pos fits i32");
        let mate_pos = pos + 300;
        let tlen = 300 + i32::try_from(seq.len()).expect("fits i32");
        let name = format!("read{i}");
        let rg_id: &[u8] = if i % 2 == 0 { b"RG1" } else { b"RG2" };
        let cell = format!("cell{}", i % 4);
        let mi = format!("umi_{}", i % 37);
        // Both mates share the same match-only CIGAR ("10M"), so each side's
        // `MC` (mate CIGAR) tag is the same literal. `samtools sort
        // --template-coordinate` refuses input without `MC` ("Please run
        // samtools fixmate on file first"); fgumi's template-coordinate sort
        // does not need it, but carrying it keeps one input usable by both
        // oracles instead of needing a `samtools fixmate`-derived copy just
        // for the samtools half.
        let mate_cigar: &[u8] = b"10M";

        let mut r1 = SamBuilder::new();
        r1.read_name(name.as_bytes())
            .ref_id(ref_id)
            .pos(pos)
            .mapq(60)
            .flags(flags::PAIRED | flags::FIRST_SEGMENT | flags::MATE_REVERSE)
            .cigar_ops(&[cigar_op])
            .sequence(seq)
            .qualities(&qual)
            .mate_ref_id(ref_id)
            .mate_pos(mate_pos)
            .template_length(tlen)
            .add_string_tag(SamTag::RG, rg_id)
            .add_string_tag(SamTag::MI, mi.as_bytes())
            .add_string_tag(SamTag::CB, cell.as_bytes())
            .add_string_tag(SamTag::MC, mate_cigar);
        records.push(r1.build());

        let mut r2 = SamBuilder::new();
        r2.read_name(name.as_bytes())
            .ref_id(ref_id)
            .pos(mate_pos)
            .mapq(60)
            .flags(flags::PAIRED | flags::LAST_SEGMENT | flags::REVERSE)
            .cigar_ops(&[cigar_op])
            .sequence(seq)
            .qualities(&qual)
            .mate_ref_id(ref_id)
            .mate_pos(pos)
            .template_length(-tlen)
            .add_string_tag(SamTag::RG, rg_id)
            .add_string_tag(SamTag::MI, mi.as_bytes())
            .add_string_tag(SamTag::CB, cell.as_bytes())
            .add_string_tag(SamTag::MC, mate_cigar);
        records.push(r2.build());
    }

    // A handful of fully unmapped pairs, to exercise the "no coordinate" tail
    // that coordinate and template-coordinate order both have to place.
    for i in 0..8 {
        let name = format!("unmapped{i}");
        for segment_flag in [flags::FIRST_SEGMENT, flags::LAST_SEGMENT] {
            let mut b = SamBuilder::new();
            b.read_name(name.as_bytes())
                .flags(flags::PAIRED | flags::UNMAPPED | flags::MATE_UNMAPPED | segment_flag)
                .sequence(seq)
                .qualities(&qual)
                .add_string_tag(SamTag::RG, b"RG1");
            records.push(b.build());
        }
    }

    records
}

/// Maps an `#[case]` order label to the `--order` flag value `fgumi sort`
/// accepts. The case labels use hyphens (`queryname-natural`) to read cleanly
/// as rstest case names; the CLI's own spelling uses `::` sub-sort syntax.
fn order_flag(order_label: &str) -> &'static str {
    match order_label {
        "coordinate" => "coordinate",
        "queryname" => "queryname",
        "queryname-natural" => "queryname::natural",
        "template-coordinate" => "template-coordinate",
        other => panic!("unhandled order label {other:?} in test case table"),
    }
}

/// Four-order parity gate: the cutover build's sort output must be
/// byte-identical (modulo `@PG`) to the saved owned-engine baseline binary,
/// and -- where samtools has an equivalent sort mode -- order-valid with the
/// same record multiset as samtools' own output.
///
/// **Which oracles run where:** on a host with samtools 1.23.1 on `PATH` but
/// no `FGUMI_BASELINE_BIN` (e.g. plain CI), the coordinate, queryname-natural,
/// and template-coordinate cases still run their samtools cross-check;
/// queryname-lex has no samtools equivalent (samtools' `-n` is natural order,
/// not lexicographic), so it is baseline-only and is *skipped* there. Setting
/// `FGUMI_BASELINE_BIN` (when it names an existing file) adds the baseline half
/// to every case, queryname-lex included.
#[rstest]
#[case::coordinate("coordinate", Some(&["sort"] as &[&str]))]
#[case::queryname_lex("queryname", None)] // no samtools equivalent: samtools -n is natural order
#[case::queryname_natural("queryname-natural", Some(&["sort", "-n"] as &[&str]))]
#[case::template_coordinate(
    "template-coordinate",
    Some(&["sort", "--template-coordinate"] as &[&str])
)]
fn cutover_matches_baseline_and_samtools(
    #[case] order: &str,
    #[case] samtools_args: Option<&[&str]>,
) {
    let flag = order_flag(order);

    let baseline = baseline_bin();
    let run_samtools = samtools_args.is_some() && samtools_available();

    if baseline.is_none() && !run_samtools {
        eprintln!(
            "SKIP cutover_matches_baseline_and_samtools[{order}]: neither FGUMI_BASELINE_BIN \
             (unset or naming a missing file) nor a usable samtools cross-check is available -- \
             no oracle to compare against"
        );
        return;
    }

    let dir = TempDir::new().expect("create temp dir");
    let input_bam = dir.path().join("in.bam");
    write_bam(&input_bam, &diverse_header(), &diverse_unsorted_records(300));

    let current_bin = Path::new(env!("CARGO_BIN_EXE_fgumi"));
    let current_out = dir.path().join("current.bam");
    run_fgumi_sort(current_bin, &input_bam, &current_out, flag);

    match baseline {
        Some(baseline_bin_path) => {
            let baseline_out = dir.path().join("baseline.bam");
            run_fgumi_sort(&baseline_bin_path, &input_bam, &baseline_out, flag);

            let current_bytes = decompressed_records_without_pg(&current_out);
            let baseline_bytes = decompressed_records_without_pg(&baseline_out);
            assert_eq!(
                current_bytes,
                baseline_bytes,
                "cutover output for order={order} diverges from the owned-engine baseline \
                 binary ({}) after stripping the @PG header line -- this is a real cutover \
                 parity bug, not something to relax the assertion for",
                baseline_bin_path.display(),
            );
        }
        None => {
            eprintln!(
                "SKIP baseline half of cutover_matches_baseline_and_samtools[{order}]: \
                 FGUMI_BASELINE_BIN is unset or does not name an existing file"
            );
        }
    }

    match samtools_args {
        Some(args) if run_samtools => {
            let samtools_out = dir.path().join("samtools.bam");
            run_samtools_sort(args, &input_bam, &samtools_out);

            assert!(
                fgumi_verify_sorted(current_bin, &current_out, flag),
                "cutover output for order={order} is not order-valid per `fgumi sort --verify`"
            );

            let current_identities = record_identities(&current_out);
            let samtools_identities = record_identities(&samtools_out);
            assert_eq!(
                current_identities, samtools_identities,
                "cutover output for order={order} is not the same multiset of records \
                 (QNAME/FLAG/RNAME/POS) as samtools' output"
            );
        }
        Some(_) => {
            eprintln!(
                "SKIP samtools half of cutover_matches_baseline_and_samtools[{order}]: \
                 samtools not found on PATH"
            );
        }
        None => {
            eprintln!(
                "no samtools equivalent for order={order} (queryname-lex is fgumi-specific, \
                 samtools' -n is natural order) -- baseline-only oracle"
            );
        }
    }
}

// ============================================================================
// Task 5: --write-index, edge cases, spill identity, env fallback, and guards
//
// Command-level coverage of the pieces the cutover changed or left inert:
// the coordinate-only `--write-index` guard and its inline BAI content, the
// empty/single/already-sorted edges, spill-vs-in-memory byte identity,
// `FGUMI_TMP_DIRS` reaching the chain, and the `--threads 0` /
// accept-but-inert-flag guards. A genuine failure here is a real cutover bug,
// not a test to relax -- see the module-level correctness discipline note.
// ============================================================================

/// `n` mapped, single-end records on one reference, placed at *ascending*
/// positions -- already in coordinate order. Used for the "already sorted"
/// edge case; at `n <= 1` it also serves as the empty/single-record input,
/// where position order is moot.
fn sorted_records(n: usize) -> Vec<RawRecord> {
    (0..n)
        .map(|i| {
            let pos = i32::try_from((i + 1) * 100).expect("pos fits i32");
            let mut b = SamBuilder::new();
            b.read_name(format!("read{i}").as_bytes())
                .ref_id(0)
                .pos(pos)
                .mapq(60)
                .flags(0)
                .cigar_ops(&[4u32 << 4]) // 4M
                .sequence(b"ACGT")
                .qualities(&[30u8; 4]);
            b.build()
        })
        .collect()
}

/// `samtools idxstats <bam>` output: per-reference mapped/unmapped counts
/// read off the `.bai` sidecar currently sitting next to `bam`.
fn idxstats(bam: &Path) -> String {
    let out = Command::new("samtools")
        .args(["idxstats", bam.to_str().expect("path is valid UTF-8")])
        .output()
        .expect("failed to spawn samtools idxstats");
    assert!(
        out.status.success(),
        "samtools idxstats failed for {}: {}",
        bam.display(),
        String::from_utf8_lossy(&out.stderr)
    );
    String::from_utf8_lossy(&out.stdout).into_owned()
}

/// `samtools view <bam> <region>` output: the records a region query returns
/// via the `.bai` sidecar currently sitting next to `bam`.
fn region_records(bam: &Path, region: &str) -> String {
    let out = Command::new("samtools")
        .args(["view", bam.to_str().expect("path is valid UTF-8"), region])
        .output()
        .expect("failed to spawn samtools view");
    assert!(
        out.status.success(),
        "samtools view {region} failed for {}: {}",
        bam.display(),
        String::from_utf8_lossy(&out.stderr)
    );
    String::from_utf8_lossy(&out.stdout).into_owned()
}

/// `fgumi sort --write-index` (coordinate) writes the output BAM plus an
/// inline `.bai` sidecar built by the arena indexer (PR A0). That `.bai` must
/// be equivalent to one `samtools index` builds for the identical output
/// bytes: identical `idxstats`, identical region-query results.
///
/// The inline and samtools-built indexes are compared by staging each, in
/// turn, at the BAM's default `.bai` sidecar path (only one can occupy that
/// path at a time), rather than by asking samtools to look at a differently
/// named index file, which it has no flag for.
#[test]
fn cutover_write_index_matches_samtools_idxstats() {
    if !samtools_available() {
        eprintln!("SKIP cutover_write_index_matches_samtools_idxstats: samtools not on PATH");
        return;
    }

    let dir = TempDir::new().expect("create temp dir");
    let input_bam = dir.path().join("in.bam");
    write_bam(&input_bam, &diverse_header(), &diverse_unsorted_records(1_000));

    let output_bam = dir.path().join("out.bam");
    let bai_sidecar = dir.path().join("out.bam.bai");

    let current_bin = Path::new(env!("CARGO_BIN_EXE_fgumi"));
    let status = Command::new(current_bin)
        .args([
            OsStr::new("sort"),
            OsStr::new("-i"),
            input_bam.as_os_str(),
            OsStr::new("-o"),
            output_bam.as_os_str(),
            OsStr::new("--order"),
            OsStr::new("coordinate"),
            OsStr::new("--write-index"),
            OsStr::new("-m"),
            OsStr::new("64K"),
        ])
        .status()
        .expect("failed to spawn fgumi sort --write-index");
    assert!(status.success(), "fgumi sort --write-index failed");
    assert!(output_bam.exists(), "sorted output BAM was not created");
    assert!(bai_sidecar.exists(), "inline .bai was not created next to the output BAM");

    // Move the inline index aside, then let samtools build its own index over
    // the identical output bytes.
    let fgumi_bai = dir.path().join("fgumi.bai");
    fs::rename(&bai_sidecar, &fgumi_bai).expect("move inline .bai aside");

    let index_status = Command::new("samtools")
        .args(["index", output_bam.to_str().expect("path is valid UTF-8")])
        .status()
        .expect("failed to spawn samtools index");
    assert!(index_status.success(), "samtools index failed");
    assert!(bai_sidecar.exists(), "samtools did not (re)create the .bai sidecar");
    let samtools_bai = dir.path().join("samtools.bai");
    fs::rename(&bai_sidecar, &samtools_bai).expect("move samtools .bai aside");

    let regions = ["chr1", "chr2", "chr3", "chr1:300000-350000", "chr2:350000-400500"];

    fs::copy(&fgumi_bai, &bai_sidecar).expect("stage fgumi .bai");
    let fgumi_idxstats = idxstats(&output_bam);
    let fgumi_regions: Vec<String> =
        regions.iter().map(|r| region_records(&output_bam, r)).collect();
    fs::remove_file(&bai_sidecar).expect("remove staged fgumi .bai");

    fs::copy(&samtools_bai, &bai_sidecar).expect("stage samtools .bai");
    let samtools_idxstats = idxstats(&output_bam);
    let samtools_regions: Vec<String> =
        regions.iter().map(|r| region_records(&output_bam, r)).collect();
    fs::remove_file(&bai_sidecar).expect("remove staged samtools .bai");

    assert_eq!(
        fgumi_idxstats, samtools_idxstats,
        "idxstats via fgumi's inline .bai must match samtools' own index over the same bytes"
    );
    for (region, (fgumi_out, samtools_out)) in
        regions.iter().zip(fgumi_regions.iter().zip(samtools_regions.iter()))
    {
        assert_eq!(
            fgumi_out, samtools_out,
            "region {region}: records retrieved via fgumi's inline .bai differ from samtools' \
             own index over the identical output bytes"
        );
    }
}

/// `--write-index` is coordinate-only; `--order queryname --write-index` must
/// be rejected up front (before any output is opened) and leave no partial
/// output or `.bai` sidecar behind.
#[test]
fn cutover_write_index_rejects_non_coordinate() {
    let dir = TempDir::new().expect("create temp dir");
    let input_bam = dir.path().join("in.bam");
    write_bam(&input_bam, &create_minimal_header("chr1", 1_000_000), &unsorted_records(10));
    let output_bam = dir.path().join("out.bam");

    let cmd = Sort::try_parse_from([
        OsStr::new("sort"),
        OsStr::new("-i"),
        input_bam.as_os_str(),
        OsStr::new("-o"),
        output_bam.as_os_str(),
        OsStr::new("--order"),
        OsStr::new("queryname"),
        OsStr::new("--write-index"),
    ])
    .expect("failed to parse sort args");

    let result = cmd.execute("fgumi sort");
    assert!(result.is_err(), "--order queryname --write-index must be rejected");
    let err_msg = format!("{:#}", result.unwrap_err());
    assert!(
        err_msg.contains("--write-index is only valid for coordinate sort"),
        "unexpected error message for --order queryname --write-index: {err_msg}"
    );

    assert!(!output_bam.exists(), "a rejected --write-index run must leave no partial output BAM");
    let bai_sidecar = dir.path().join("out.bam.bai");
    assert!(!bai_sidecar.exists(), "a rejected --write-index run must leave no partial .bai");
}

/// Empty, single-record, and already-coordinate-sorted inputs must all
/// produce valid, correctly-counted output through the cutover, and --
/// wherever the owned-engine baseline binary is available -- byte-identical
/// output (modulo `@PG`) to it.
#[rstest]
#[case::empty(0)]
#[case::single(1)]
#[case::already_sorted(1000)]
fn cutover_edge_cases_match_baseline(#[case] n: usize) {
    let dir = TempDir::new().expect("create temp dir");
    let input_bam = dir.path().join("in.bam");
    write_bam(&input_bam, &create_minimal_header("chr1", 1_000_000), &sorted_records(n));

    let current_bin = Path::new(env!("CARGO_BIN_EXE_fgumi"));
    let current_out = dir.path().join("current.bam");
    run_fgumi_sort(current_bin, &input_bam, &current_out, "coordinate");

    let (_, out_records) = read_bam_output(&current_out);
    assert_eq!(out_records.len(), n, "output record count must match input record count (n={n})");
    // The baseline half below is skipped whenever FGUMI_BASELINE_BIN is unset, so
    // assert order here too: a count-only assertion cannot see a reordered output.
    assert!(
        fgumi_verify_sorted(current_bin, &current_out, "coordinate"),
        "cutover output for n={n} is not order-valid per `fgumi sort --verify`"
    );

    match baseline_bin() {
        Some(baseline_bin_path) => {
            let baseline_out = dir.path().join("baseline.bam");
            run_fgumi_sort(&baseline_bin_path, &input_bam, &baseline_out, "coordinate");

            assert_eq!(
                decompressed_records_without_pg(&current_out),
                decompressed_records_without_pg(&baseline_out),
                "cutover output for n={n} diverges from the owned-engine baseline binary after \
                 stripping the @PG header line -- this is a real cutover parity bug, not \
                 something to relax the assertion for"
            );
        }
        None => {
            eprintln!(
                "SKIP baseline half of cutover_edge_cases_match_baseline[n={n}]: \
                 FGUMI_BASELINE_BIN is unset or does not name an existing file"
            );
        }
    }
}

/// The same unsorted input, sorted once under a tiny `--max-memory` budget
/// (forcing spill into a `--tmp-dir` tempfile directory) and once under the
/// (large) default budget (fitting entirely in memory), must produce
/// byte-identical output -- spilling must never change what the sort emits,
/// only how it gets there.
///
/// The spill branch runs as a subprocess (mirroring
/// `cutover_honors_fgumi_tmp_dirs_env`) so the test can assert on
/// `SortSummaryFinalizeHook`'s `"Spill runs: {n}"` stderr line -- only emitted
/// when `runs_written > 0` -- rather than merely assuming a 64K budget over
/// 5,000 records spills. Without that check, a future change to the
/// memory-budget math that silently stopped the spill would leave this test
/// green while comparing two in-memory sorts, for the wrong reason.
#[test]
fn cutover_spill_and_nonspill_are_identical() {
    let dir = TempDir::new().expect("create temp dir");
    let input_bam = dir.path().join("in.bam");
    write_bam(&input_bam, &create_minimal_header("chr1", 1_000_000), &unsorted_records(5_000));

    let spill_tmp_dir = TempDir::new().expect("create spill tmp dir");
    let spill_out = dir.path().join("spill.bam");
    let current_bin = Path::new(env!("CARGO_BIN_EXE_fgumi"));
    let spill_output = Command::new(current_bin)
        .env("RUST_LOG", "info")
        .args([
            OsStr::new("sort"),
            OsStr::new("-i"),
            input_bam.as_os_str(),
            OsStr::new("-o"),
            spill_out.as_os_str(),
            OsStr::new("--order"),
            OsStr::new("coordinate"),
            OsStr::new("--max-memory"),
            OsStr::new("64K"),
            OsStr::new("--tmp-dir"),
            spill_tmp_dir.path().as_os_str(),
        ])
        .output()
        .expect("failed to spawn fgumi sort (spill)");
    let spill_stderr = String::from_utf8_lossy(&spill_output.stderr);
    assert!(spill_output.status.success(), "spill sort should succeed:\n{spill_stderr}");
    assert!(
        spill_stderr.contains("Spill runs:"),
        "expected at least one spill run under a 64K budget over 5,000 records -- otherwise \
         this test cannot prove the spill branch actually spilled, only that it resolved a \
         tmp-dir nothing wrote to; stderr:\n{spill_stderr}"
    );

    // Default --max-memory (768M/thread) trivially fits 5,000 tiny records in
    // memory with zero spills. Run this arm as a subprocess too (mirroring the
    // spill arm) so it can assert `"Spill runs:"` is *absent* -- otherwise a
    // future memory-budget change that made the default path spill would leave
    // this test comparing two spill sorts while still named "non-spill".
    let nonspill_out = dir.path().join("nonspill.bam");
    let nonspill_output = Command::new(current_bin)
        .env("RUST_LOG", "info")
        .args([
            OsStr::new("sort"),
            OsStr::new("-i"),
            input_bam.as_os_str(),
            OsStr::new("-o"),
            nonspill_out.as_os_str(),
            OsStr::new("--order"),
            OsStr::new("coordinate"),
        ])
        .output()
        .expect("failed to spawn fgumi sort (non-spill)");
    let nonspill_stderr = String::from_utf8_lossy(&nonspill_output.stderr);
    assert!(nonspill_output.status.success(), "non-spill sort should succeed:\n{nonspill_stderr}");
    assert!(
        !nonspill_stderr.contains("Spill runs:"),
        "the default budget must keep this sort in memory, otherwise this test compares two \
         spill sorts; stderr:\n{nonspill_stderr}"
    );

    // Assert order against an independent oracle before the byte-identity check:
    // comparing two fgumi runs alone cannot catch a spill bug that misorders both
    // the same way.
    assert!(
        fgumi_verify_sorted(current_bin, &spill_out, "coordinate"),
        "spill output is not order-valid per `fgumi sort --verify`"
    );
    assert_eq!(
        decompressed_records_without_pg(&spill_out),
        decompressed_records_without_pg(&nonspill_out),
        "spill and non-spill sorts of the same input must produce byte-identical output"
    );
}

/// `FGUMI_TMP_DIRS` must reach `SortOptions.tmp_dirs` via `execute_sort` when
/// no `--tmp-dir` flag is given, and the sort must actually spill into it.
///
/// This has to run as a subprocess: `std::env::set_var` is `unsafe` under
/// edition 2024 and this crate forbids unsafe code, so the only way to set an
/// environment variable for a run without touching the current process's
/// environment is to hand it to a child's `Command::env`.
///
/// Evidence gathered from the child's stderr (default `info`-level logging,
/// pinned explicitly via `RUST_LOG` in case the ambient environment overrides
/// it):
/// - `execute_sort`'s own `"Temp directories: {joined}"` banner names exactly
///   the resolved `tmp_dirs`, so its presence with our env tempdir's path
///   proves the env var reached `SortOptions.tmp_dirs`.
/// - `SortSummaryFinalizeHook`'s `"Spill runs: {n}"` line is only emitted
///   when `runs_written > 0`, so its presence proves the tiny `--max-memory`
///   budget actually forced the sort to spill (into that same resolved
///   directory) rather than merely resolving a path nothing used.
#[test]
fn cutover_honors_fgumi_tmp_dirs_env() {
    let dir = TempDir::new().expect("create temp dir");
    let input_bam = dir.path().join("in.bam");
    let records = unsorted_records(8_000);
    let input_count = records.len();
    write_bam(&input_bam, &create_minimal_header("chr1", 1_000_000), &records);
    let output_bam = dir.path().join("out.bam");

    let env_tmp_dir = TempDir::new().expect("create env tmp dir");

    let output = Command::new(env!("CARGO_BIN_EXE_fgumi"))
        .env("RUST_LOG", "info")
        .env("FGUMI_TMP_DIRS", env_tmp_dir.path())
        .args([
            OsStr::new("sort"),
            OsStr::new("-i"),
            input_bam.as_os_str(),
            OsStr::new("-o"),
            output_bam.as_os_str(),
            OsStr::new("--order"),
            OsStr::new("coordinate"),
            OsStr::new("--max-memory"),
            OsStr::new("64K"),
        ])
        // --tmp-dir deliberately left unset, so tmp_dirs can only have come
        // from the FGUMI_TMP_DIRS env var above.
        .output()
        .expect("failed to spawn fgumi sort subprocess");

    let stderr = String::from_utf8_lossy(&output.stderr);
    assert!(output.status.success(), "fgumi sort (FGUMI_TMP_DIRS) failed:\n{stderr}");

    let expected_banner = format!("Temp directories: {}", env_tmp_dir.path().display());
    assert!(
        stderr.contains(&expected_banner),
        "expected the resolved-tmp-dirs banner to name the FGUMI_TMP_DIRS path \
         ({expected_banner:?}); stderr:\n{stderr}"
    );
    assert!(
        stderr.contains("Spill runs:"),
        "expected at least one spill run under a 64K budget over {input_count} records -- \
         otherwise this test cannot prove the env-resolved directory was actually used, only \
         that it was resolved; stderr:\n{stderr}"
    );

    let (_, out_records) = read_bam_output(&output_bam);
    assert_eq!(out_records.len(), input_count, "output record count must match input");
    // A record count cannot see a reordering; assert order against the independent
    // `--verify` oracle so a spill that misordered records fails here.
    assert!(
        fgumi_verify_sorted(Path::new(env!("CARGO_BIN_EXE_fgumi")), &output_bam, "coordinate"),
        "FGUMI_TMP_DIRS spill output is not order-valid per `fgumi sort --verify`"
    );
}

/// `fgumi sort --threads 0` must be rejected cleanly, before any output is
/// created.
///
/// In practice this bails out of `resolve_memory_budget` inside
/// `execute_sort`'s pre-chain setup (`self.memory_budget_threads()` returns 0
/// verbatim for `--threads 0`, and `resolve_memory_budget` rejects a
/// zero-thread budget with this exact message) -- the chain is never built
/// for this case, so `ChainBuilder::new`'s own identical `num_threads == 0`
/// guard is not what fires here. Both call sites bail with the same message
/// ("--threads must be at least 1"), so this test's assertion holds
/// regardless of which one is reached; it also cross-checks against the
/// owned-engine baseline binary where available, since the baseline's
/// `execute_sort` hits the identical `resolve_memory_budget` guard (the
/// chain-cutover in Task 3 did not touch this code path).
#[test]
fn cutover_threads_zero_is_rejected() {
    let dir = TempDir::new().expect("create temp dir");
    let input_bam = dir.path().join("in.bam");
    write_bam(&input_bam, &create_minimal_header("chr1", 1_000_000), &unsorted_records(10));
    let output_bam = dir.path().join("out.bam");

    let cmd = Sort::try_parse_from([
        OsStr::new("sort"),
        OsStr::new("-i"),
        input_bam.as_os_str(),
        OsStr::new("-o"),
        output_bam.as_os_str(),
        OsStr::new("--order"),
        OsStr::new("coordinate"),
        OsStr::new("--threads"),
        OsStr::new("0"),
    ])
    .expect("failed to parse sort args");

    let result = cmd.execute("fgumi sort");
    assert!(result.is_err(), "--threads 0 must be rejected");
    let err_msg = format!("{:#}", result.unwrap_err());
    assert!(
        err_msg.contains("--threads must be at least 1"),
        "unexpected error message for --threads 0: {err_msg}"
    );
    assert!(!output_bam.exists(), "a rejected --threads 0 run must leave no partial output");

    match baseline_bin() {
        Some(baseline_bin_path) => {
            let baseline_output = Command::new(&baseline_bin_path)
                .args([
                    OsStr::new("sort"),
                    OsStr::new("-i"),
                    input_bam.as_os_str(),
                    OsStr::new("-o"),
                    output_bam.as_os_str(),
                    OsStr::new("--order"),
                    OsStr::new("coordinate"),
                    OsStr::new("--threads"),
                    OsStr::new("0"),
                ])
                .output()
                .expect("failed to spawn baseline binary");
            assert!(
                !baseline_output.status.success(),
                "baseline binary's --threads 0 unexpectedly succeeded"
            );
            let baseline_err = String::from_utf8_lossy(&baseline_output.stderr);
            assert!(
                baseline_err.contains("--threads must be at least 1"),
                "baseline binary's --threads 0 error message differs from the cutover's \
                 (\"--threads must be at least 1\"): {baseline_err}"
            );
        }
        None => {
            eprintln!(
                "SKIP baseline comparison in cutover_threads_zero_is_rejected: \
                 FGUMI_BASELINE_BIN is unset or does not name an existing file"
            );
        }
    }
}

/// `--sort-stats` and `--read-streams` are inert on the chain: neither field
/// even exists on the chain-facing `SortOptions` struct (`to_sort_options`
/// drops both), so nothing downstream reads them. Restoring
/// `--read-streams`'s real behavior is tracked as R7b. Passing them must be
/// accept-but-inert -- not a hard error -- and the sort must still produce
/// correct, coordinate-sorted output.
#[test]
fn cutover_inert_flags_do_not_error() {
    let dir = TempDir::new().expect("create temp dir");
    let input_bam = dir.path().join("in.bam");
    let records = unsorted_records(200);
    let input_count = records.len();
    write_bam(&input_bam, &create_minimal_header("chr1", 1_000_000), &records);
    let output_bam = dir.path().join("out.bam");

    // Run as a subprocess (with `RUST_LOG=info` pinned so the ambient environment
    // cannot filter the notices out) so the test can assert on stderr: the ignored
    // flags are `warn`-level and each must announce that it is inert, otherwise a
    // silently-dropped flag would leave the user thinking it took effect.
    let current_bin = Path::new(env!("CARGO_BIN_EXE_fgumi"));
    let output = Command::new(current_bin)
        .env("RUST_LOG", "info")
        .args([
            OsStr::new("sort"),
            OsStr::new("-i"),
            input_bam.as_os_str(),
            OsStr::new("-o"),
            output_bam.as_os_str(),
            OsStr::new("--order"),
            OsStr::new("coordinate"),
            OsStr::new("--sort-stats"),
            OsStr::new("--read-streams"),
            OsStr::new("4"),
        ])
        .output()
        .expect("failed to spawn fgumi sort (inert flags)");
    let stderr = String::from_utf8_lossy(&output.stderr);
    assert!(
        output.status.success(),
        "--sort-stats/--read-streams are accept-but-inert on the chain and must not error:\n{stderr}"
    );
    assert!(
        stderr.contains("--read-streams=4 is ignored by the sort chain"),
        "the ignored `--read-streams` notice must be visible on stderr; stderr:\n{stderr}"
    );
    assert!(
        stderr.contains("--sort-stats is ignored by the sort chain"),
        "the ignored `--sort-stats` notice must be visible on stderr; stderr:\n{stderr}"
    );

    let (_, out_records) = read_bam_output(&output_bam);
    assert_eq!(out_records.len(), input_count, "output record count must match input");
    let keys: Vec<(usize, usize)> = out_records
        .iter()
        .map(|r| {
            (
                r.reference_sequence_id().unwrap_or(usize::MAX),
                r.alignment_start().map_or(0, usize::from),
            )
        })
        .collect();
    assert!(
        keys.windows(2).all(|w| w[0] <= w[1]),
        "output is not coordinate-sorted with --sort-stats/--read-streams set: {keys:?}"
    );
}

/// The owned engine always verified BGZF CRC32, including on stdin input
/// (`decompress_block` has no CRC opt-out at all). Standalone `fgumi sort`'s
/// input decode must reject a corrupted block on stdin too, not silently sort
/// past it.
///
/// Only the block's footer CRC32 is corrupted (one bit flipped, mirroring
/// `decompress_opts_skips_crc_on_stored_block` in `fgumi-bgzf`): the
/// compressed payload is left untouched and decodes normally to the same
/// bytes, so this corruption is detectable *only* by comparing the
/// decompressed output's CRC32 against the (now-wrong) footer value. A
/// structural decode failure could never catch it, so a pass here proves a
/// live CRC check ran -- not just that some unrelated error-detection caught
/// the file.
///
/// 20,000 records (rather than a handful) is deliberate: it forces the
/// uncompressed SAM header + record stream well past a single 64 KiB BGZF
/// block, so the *last* real block is a pure record (body) block, never the
/// one `fgumi_bam_io::read_header_and_replay` decompresses in full while
/// parsing the header via noodles. Corrupting the last block therefore
/// exercises stdin sort's own record-decode path (the arena ingest
/// `InflateToArena` uses for the sole-`[Stage::Sort]` chain), not just the
/// separate header-parse tee.
///
/// **Not a RED/GREEN gate for the `verify_crc: true` change in
/// `execute_sort`.** This test passes identically with or without that flag:
/// standalone sort's `[Stage::Sort]`-only chain always takes the arena-ingest
/// path (`InflateToArena` -> `fgumi_bgzf::decompress_into_slice`), which has
/// no CRC opt-out and checks unconditionally regardless of `ChainSpec.verify_crc`
/// -- the flag is only read by `build_bam_decode_preamble`'s `BgzfDecompress`,
/// a path standalone sort never takes. It exists as a regression guard on its
/// own terms (stdin corruption must be rejected end-to-end, through whichever
/// layer catches it), and as a tripwire: if a future chain-topology change
/// ever routes standalone sort through `BgzfDecompress` instead, making
/// `verify_crc` load-bearing, a stale `effective_check_crc()`-style stdin skip
/// would fail this test rather than silently reintroducing the regression.
#[test]
fn cutover_stdin_input_detects_corrupt_crc() {
    let dir = TempDir::new().expect("create temp dir");
    let input_bam = dir.path().join("in.bam");
    write_bam(&input_bam, &create_minimal_header("chr1", 2_100_000), &unsorted_records(20_000));

    let raw = fs::read(&input_bam).expect("read seed BAM");
    let mut reader = std::io::Cursor::new(raw.as_slice());
    let mut blocks = fgumi_bgzf::read_raw_blocks(&mut reader, 4_096).expect("parse BGZF blocks");
    let target_idx = blocks
        .iter()
        .enumerate()
        .filter(|(_, b)| !b.is_eof() && !b.is_empty())
        .map(|(i, _)| i)
        .next_back()
        .expect("expected at least one real (non-EOF) BGZF block to corrupt");
    assert!(
        target_idx > 0,
        "expected 20,000 records to span more than one real BGZF block (got only 1) -- \
         corrupting it would land in the header-parse block instead of a pure body block"
    );
    let target = &mut blocks[target_idx];
    let crc_off = target.data.len() - fgumi_bgzf::BGZF_FOOTER_SIZE;
    target.data[crc_off] ^= 0x01;

    let mut corrupted = Vec::with_capacity(raw.len());
    for block in &blocks {
        corrupted.extend_from_slice(&block.data);
    }
    // `read_raw_blocks` silently drops BGZF EOF-marker blocks it reads (see
    // its own doc comment), so re-append the standard marker for a
    // well-formed, correctly-terminated stream.
    corrupted.extend_from_slice(&fgumi_bgzf::BGZF_EOF);
    let corrupted_bam = dir.path().join("corrupted.bam");
    fs::write(&corrupted_bam, &corrupted).expect("write corrupted BAM");

    let output_bam = dir.path().join("out.bam");
    let current_bin = Path::new(env!("CARGO_BIN_EXE_fgumi"));
    let output = Command::new(current_bin)
        .args([
            OsStr::new("sort"),
            OsStr::new("-i"),
            OsStr::new("-"),
            OsStr::new("-o"),
            output_bam.as_os_str(),
            OsStr::new("--order"),
            OsStr::new("coordinate"),
        ])
        .stdin(std::process::Stdio::from(
            fs::File::open(&corrupted_bam).expect("open corrupted BAM to pipe"),
        ))
        .output()
        .expect("failed to spawn fgumi sort on corrupted stdin input");

    assert!(
        !output.status.success(),
        "fgumi sort must reject a corrupted-CRC BGZF block on stdin, not silently sort past it"
    );
    // Wording varies by which layer catches it (fgumi-bgzf's own "CRC32
    // mismatch" vs. noodles' "checksum mismatch" in the header-parse tee), so
    // accept either rather than pinning one literal message.
    let stderr = String::from_utf8_lossy(&output.stderr).to_lowercase();
    assert!(
        stderr.contains("crc") || stderr.contains("checksum"),
        "expected the failure to name CRC/checksum verification as the cause; stderr:\n{stderr}"
    );
    // Unlike the upfront `--write-index`/`--threads 0` guards, this failure
    // surfaces mid-stream (well after the output file was opened for
    // writing), so a truncated partial output file is expected here -- not
    // asserted against.
}
