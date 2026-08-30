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
/// `FGUMI_BASELINE_BIN` (or having the default baseline path present) adds the
/// baseline half to every case, queryname-lex included.
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
             (nor the default baseline path) nor a usable samtools cross-check is available -- \
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
                 FGUMI_BASELINE_BIN not set and the default baseline path was not found"
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
