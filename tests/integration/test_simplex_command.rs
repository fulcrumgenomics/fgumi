//! End-to-end CLI tests for the simplex command.
//!
//! These tests invoke `Simplex::execute()` in-process and validate:
//! 1. Basic simplex consensus calling from grouped reads
//! 2. Statistics output
//! 3. Rejected reads output

use bstr::BString;
use clap::Parser;
use fgumi_dna::reverse_complement;
use fgumi_lib::commands::command::Command;
use fgumi_lib::commands::simplex::Simplex;
use fgumi_lib::sam::SamTag;
use noodles::bam;
use noodles::sam::Header;
use noodles::sam::alignment::io::Write as AlignmentWrite;
use noodles::sam::alignment::record_buf::RecordBuf;
use noodles::sam::header::record::value::Map;
use noodles::sam::header::record::value::map::ReadGroup;
use noodles::sam::header::record::value::map::read_group::tag as rg_tag;
use std::fs;
use std::num::NonZeroUsize;
use std::path::Path;
use tempfile::TempDir;

use crate::helpers::bam_generator::{create_minimal_header, create_umi_family, to_record_buf};

/// Write grouped BAM file (reads grouped by MI tag).
fn create_grouped_bam(path: &Path, families: Vec<(&str, Vec<fgumi_raw_bam::RawRecord>)>) {
    let header = create_minimal_header("chr1", 10000);
    let mut writer =
        bam::io::Writer::new(fs::File::create(path).expect("Failed to create BAM file"));
    writer.write_header(&header).expect("Failed to write header");

    for (mi, records) in families {
        for raw in &records {
            // Convert to RecordBuf, add MI tag, write
            use noodles::sam::alignment::record::data::field::Tag;
            use noodles::sam::alignment::record_buf::data::field::Value;
            let mut record = to_record_buf(raw);
            let mi_tag = Tag::from(fgumi_lib::sam::SamTag::MI);
            record.data_mut().insert(mi_tag, Value::from(mi));
            writer.write_alignment_record(&header, &record).expect("Failed to write record");
        }
    }
    writer.try_finish().expect("Failed to finish BAM");
}

/// Test basic simplex consensus calling.
#[test]
fn test_simplex_command_basic_consensus() {
    let temp_dir = TempDir::new().unwrap();
    let input_bam = temp_dir.path().join("input.bam");
    let output_bam = temp_dir.path().join("output.bam");

    // Create two families: 5 reads each
    let family1 = create_umi_family("ACGT", 5, "fam1", "ACGTACGT", 30);
    let family2 = create_umi_family("TGCA", 5, "fam2", "TTTTAAAA", 30);
    create_grouped_bam(&input_bam, vec![("1", family1), ("2", family2)]);

    let cmd = Simplex::try_parse_from([
        "simplex",
        "--input",
        input_bam.to_str().unwrap(),
        "--output",
        output_bam.to_str().unwrap(),
        "--min-reads",
        "2",
        "--compression-level",
        "1",
    ])
    .expect("failed to parse simplex args");
    cmd.execute("fgumi simplex").expect("Simplex command failed");
    assert!(output_bam.exists(), "Output BAM not created");

    // Read output and verify consensus reads were produced
    let mut reader = bam::io::Reader::new(fs::File::open(&output_bam).unwrap());
    let _header = reader.read_header().unwrap();
    let mut count = 0;
    for result in reader.records() {
        let record = result.expect("Failed to read record");
        // Verify consensus tags exist
        let cd_tag = SamTag::CD.to_noodles_tag();
        assert!(record.data().get(&cd_tag).is_some(), "Consensus should have cD tag");
        count += 1;
    }
    assert!(count > 0, "Should have produced consensus reads");
}

/// Consensus-unify regression: `simplex` with `--threads` unset must produce
/// output identical to `--threads 1`. Both now route through the chain (unset
/// resolves to a one-worker chain), so this pins that the fold of the former
/// single-threaded path onto the chain is output-preserving.
#[test]
fn simplex_no_threads_matches_threads_1() {
    let temp_dir = TempDir::new().unwrap();
    let input_bam = temp_dir.path().join("input.bam");
    let family1 = create_umi_family("ACGT", 5, "fam1", "ACGTACGT", 30);
    let family2 = create_umi_family("TGCA", 5, "fam2", "TTTTAAAA", 30);
    create_grouped_bam(&input_bam, vec![("1", family1), ("2", family2)]);

    let run = |out: &Path, stats: &Path, threads: Option<&str>| {
        let mut args: Vec<String> = vec![
            "simplex".into(),
            "--input".into(),
            input_bam.to_str().unwrap().into(),
            "--output".into(),
            out.to_str().unwrap().into(),
            "--stats".into(),
            stats.to_str().unwrap().into(),
            "--min-reads".into(),
            "2".into(),
            "--compression-level".into(),
            "1".into(),
        ];
        if let Some(t) = threads {
            args.push("--threads".into());
            args.push(t.into());
        }
        Simplex::try_parse_from(args).expect("parse").execute("fgumi simplex").expect("simplex");
    };
    let out_none = temp_dir.path().join("none.bam");
    let out_t1 = temp_dir.path().join("t1.bam");
    let stats_none = temp_dir.path().join("none.stats.txt");
    let stats_t1 = temp_dir.path().join("t1.stats.txt");
    run(&out_none, &stats_none, None);
    run(&out_t1, &stats_t1, Some("1"));

    let read = |p: &Path| {
        let mut r = bam::io::Reader::new(fs::File::open(p).unwrap());
        let h = r.read_header().unwrap();
        let recs: Vec<RecordBuf> = r.record_bufs(&h).map(|x| x.unwrap()).collect();
        (h, recs)
    };
    let (h_none, recs_none) = read(&out_none);
    let (h_t1, recs_t1) = read(&out_t1);
    assert_eq!(
        recs_none, recs_t1,
        "simplex --threads unset must match --threads 1 (consensus records)"
    );
    assert_eq!(h_none, h_t1, "simplex --threads unset must match --threads 1 (output header)");
    assert_eq!(
        fs::read_to_string(&stats_none).unwrap(),
        fs::read_to_string(&stats_t1).unwrap(),
        "simplex --threads unset must match --threads 1 (--stats file)"
    );

    // --- Independent oracle (NOT a snapshot) -------------------------------
    // The parity asserts above only prove the two --threads paths AGREE; a
    // shared bug in the folded chain (wrong consensus math / wrong record
    // count) would corrupt BOTH identically and still pass. Derive the
    // expected output from the input and pin it against `recs_none` (the
    // --threads-unset path that PR 549 folds onto the chain).
    //
    // Input: two clean, error-free single-end UMI families (5 identical reads
    // each), `--min-reads 2`. The consensus of N identical error-free reads is
    // the input template, and each family yields exactly one single-end
    // consensus read → 2 consensus records, one per family. Templates are
    // `ACGTACGT` (fam1) and `TTTTAAAA` (fam2). Consensus reads are emitted
    // unmapped in read (forward) orientation, so no record is
    // reverse-complemented; each stored SEQ equals its template forward. We
    // still revcomp defensively for any reverse-strand record so the oracle
    // stays correct if orientation ever changes.
    let expected_consensus_reads: usize = 2;
    assert_eq!(
        recs_none.len(),
        expected_consensus_reads,
        "two clean families → two single-end consensus reads",
    );
    // Forward-oriented SEQ of each consensus record (revcomp iff reverse-strand).
    let mut forward_seqs: Vec<Vec<u8>> = recs_none
        .iter()
        .map(|rec| {
            let seq = rec.sequence().as_ref().to_vec();
            if rec.flags().is_reverse_complemented() { reverse_complement(&seq) } else { seq }
        })
        .collect();
    forward_seqs.sort();
    // Non-vacuous: a single-base change to either expected template (e.g.
    // `ACGTACGA`) fails this equality.
    assert_eq!(
        forward_seqs,
        vec![b"ACGTACGT".to_vec(), b"TTTTAAAA".to_vec()],
        "clean families → consensus SEQ equals the input templates (forward-oriented)",
    );
    // Second independent oracle: the --stats `consensus_reads_emitted` value
    // (vertical `key<TAB>value<TAB>description` layout) must equal the count.
    let emitted = parse_stats_kv_usize(&stats_none, "consensus_reads_emitted");
    assert_eq!(
        emitted, expected_consensus_reads,
        "--stats consensus_reads_emitted must equal the emitted consensus count",
    );
}

/// Parse a `usize` value from a vertical `key<TAB>value<TAB>description` stats
/// file (the fgbio-style consensus stats layout), by key.
fn parse_stats_kv_usize(path: &Path, key: &str) -> usize {
    let text = fs::read_to_string(path).expect("read stats file");
    let prefix = format!("{key}\t");
    text.lines()
        .find(|line| line.starts_with(&prefix))
        .and_then(|line| line.split('\t').nth(1))
        .unwrap_or_else(|| panic!("stats file missing key {key}"))
        .parse()
        .expect("stats value parses as usize")
}

/// Test simplex command with statistics output.
#[test]
fn test_simplex_command_with_stats() {
    let temp_dir = TempDir::new().unwrap();
    let input_bam = temp_dir.path().join("input.bam");
    let output_bam = temp_dir.path().join("output.bam");
    let stats_path = temp_dir.path().join("stats.txt");

    let family = create_umi_family("ACGT", 5, "fam1", "ACGTACGT", 30);
    create_grouped_bam(&input_bam, vec![("1", family)]);

    let cmd = Simplex::try_parse_from([
        "simplex",
        "--input",
        input_bam.to_str().unwrap(),
        "--output",
        output_bam.to_str().unwrap(),
        "--min-reads",
        "1",
        "--stats",
        stats_path.to_str().unwrap(),
        "--compression-level",
        "1",
    ])
    .expect("failed to parse simplex args");
    cmd.execute("fgumi simplex").expect("Simplex command with stats failed");
    assert!(stats_path.exists(), "Stats file not created");
}

/// Test simplex command with rejects output.
#[test]
fn test_simplex_command_with_rejects() {
    let temp_dir = TempDir::new().unwrap();
    let input_bam = temp_dir.path().join("input.bam");
    let output_bam = temp_dir.path().join("output.bam");
    let rejects_bam = temp_dir.path().join("rejects.bam");

    // Create one family with 5 reads (passes min-reads=2) and one with 1 read (fails)
    let family1 = create_umi_family("ACGT", 5, "fam1", "ACGTACGT", 30);
    let family2 = create_umi_family("TGCA", 1, "fam2", "TTTTAAAA", 30);
    create_grouped_bam(&input_bam, vec![("1", family1), ("2", family2)]);

    let cmd = Simplex::try_parse_from([
        "simplex",
        "--input",
        input_bam.to_str().unwrap(),
        "--output",
        output_bam.to_str().unwrap(),
        "--min-reads",
        "2",
        "--rejects",
        rejects_bam.to_str().unwrap(),
        "--compression-level",
        "1",
    ])
    .expect("failed to parse simplex args");
    cmd.execute("fgumi simplex").expect("Simplex command with rejects failed");
    assert!(rejects_bam.exists(), "Rejects BAM not created");
}

/// Regression test for X5-001 at the simplex `--rejects` fan-out site: a
/// fully-clean batch must not wedge the rejects branch. The rejects fan-out's
/// branch B uses a `ByItemOrdinal` reorder stage, so every batch serial must
/// produce a branch-B item; a batch with zero rejects previously pushed nothing
/// (`Process2Output::only_a`), leaving a permanent serial gap that stalls
/// `try_pop_in_order` forever. The fix emits an empty (zero-byte) rejects block
/// for clean batches so serials stay dense.
///
/// Construction: the consensus step batches its input by a byte budget
/// (`per_step_byte_limit`, 4 MiB), so the clean prefix must exceed that budget
/// to guarantee the leading emitted batch(es) carry no rejects. Many clean
/// families — each its own MI group (depth 2, passing `--min-reads 2`) —
/// together clear 4 MiB (4000 families × 2 reads × 600 bp), so the batcher
/// flushes all-clean leading batches before a singleton family that
/// `--min-reads 2` rejects in a later batch. A *single* giant family would be
/// one MI group the batcher co-locates with the reject, so no reject-free batch
/// would form — see the construction comment in the body.
#[test]
fn test_simplex_command_rejects_all_clean_batch_does_not_wedge() {
    let temp_dir = TempDir::new().unwrap();
    let input_bam = temp_dir.path().join("input.bam");
    let output_bam = temp_dir.path().join("output.bam");
    let rejects_bam = temp_dir.path().join("rejects.bam");

    let long_seq: String = "ACGT".repeat(150); // 600 bp
    // MANY clean families — each its own MI group, depth 2 (passes --min-reads
    // 2). A single giant family is ONE MI group that the byte-budget batcher
    // co-locates with the reject, so no reject-free batch ever forms (the
    // all-clean-batch X5-001 trigger). Many families let the batcher flush
    // all-clean leading batches before the reject's batch. 4000 families × 2
    // reads × 600 bp clears the 4 MiB `per_step_byte_limit` several times over.
    let family_count = 4000;
    let mi_tags: Vec<String> = (0..family_count).map(|i| format!("mi{i:04}")).collect();
    let mut families: Vec<(&str, Vec<fgumi_raw_bam::RawRecord>)> =
        Vec::with_capacity(family_count + 1);
    for (i, mi) in mi_tags.iter().enumerate() {
        families.push((
            mi.as_str(),
            create_umi_family("ACGTACGT", 2, &format!("clean{i}"), &long_seq, 30),
        ));
    }
    // Singleton family rejected by --min-reads 2, written last so it lands in
    // the final batch — after at least one all-clean leading batch.
    families.push(("zzz_reject", create_umi_family("TTTTTTTT", 1, "rej", &long_seq, 30)));
    create_grouped_bam(&input_bam, families);

    let cmd = Simplex::try_parse_from([
        "simplex",
        "--input",
        input_bam.to_str().unwrap(),
        "--output",
        output_bam.to_str().unwrap(),
        "--rejects",
        rejects_bam.to_str().unwrap(),
        "--min-reads",
        "2",
        "--threads",
        "4",
        // Short deadlock window: a regression wedges and is failed in ~18s
        // (3s × the 6× fatal multiplier) instead of the ~60s default; the fixed
        // path finishes in well under a second.
        "--deadlock-timeout",
        "3",
        "--compression-level",
        "1",
    ])
    .expect("failed to parse simplex args");
    // Pre-fix this wedges on the rejects sink; post-fix it returns.
    cmd.execute("fgumi simplex")
        .expect("simplex --rejects must complete even with an all-clean leading batch");

    assert!(rejects_bam.exists(), "Rejects BAM not created");
    crate::helpers::assertions::assert_has_bgzf_eof(&rejects_bam);

    // The singleton "rej_0" read must reach --rejects byte-identical to its
    // input record (the reject path passes source reads through unchanged).
    // Compare full `RecordBuf`s against the record as written to the input BAM
    // — that captures every field/tag (incl. the MI added at grouping time),
    // so a mutated sequence/quality/tag on the passed-through reject is caught,
    // not just a wrong QNAME.
    let expected: Vec<RecordBuf> = crate::helpers::assertions::read_record_bufs(&input_bam)
        .into_iter()
        .filter(|r| r.name().is_some_and(|n| n == "rej_0"))
        .collect();
    assert_eq!(expected.len(), 1, "fixture sanity: exactly one reject input record");

    let observed = crate::helpers::assertions::read_record_bufs(&rejects_bam);
    assert_eq!(
        observed, expected,
        "the singleton reject must reach --rejects byte-identical to its input record"
    );
}

/// Creates a header with multiple read groups, each with SM, LB, PL tags.
fn create_header_with_read_groups(ref_name: &str, ref_len: usize) -> Header {
    use noodles::sam::header::record::value::map::Header as HeaderRecord;
    use noodles::sam::header::record::value::map::Map as HeaderRecordMap;
    use noodles::sam::header::record::value::map::ReferenceSequence;
    use noodles::sam::header::record::value::map::header::tag::Tag as HeaderTag;

    let HeaderTag::Other(sort_order_tag) = HeaderTag::from([b'S', b'O']) else { unreachable!() };
    let HeaderTag::Other(group_order_tag) = HeaderTag::from([b'G', b'O']) else { unreachable!() };
    let HeaderTag::Other(sub_sort_tag) = HeaderTag::from([b'S', b'S']) else { unreachable!() };

    let header_map = HeaderRecordMap::<HeaderRecord>::builder()
        .insert(sort_order_tag, "unsorted")
        .insert(group_order_tag, "query")
        .insert(sub_sort_tag, "template-coordinate")
        .build()
        .expect("valid header map");

    let reference_sequence = Map::<ReferenceSequence>::new(
        NonZeroUsize::new(ref_len).expect("reference length must be non-zero"),
    );

    let rg1 = Map::<ReadGroup>::builder()
        .insert(rg_tag::SAMPLE, String::from("SampleA"))
        .insert(rg_tag::LIBRARY, String::from("LibA"))
        .insert(rg_tag::PLATFORM, String::from("illumina"))
        .build()
        .expect("valid RG1");

    let rg2 = Map::<ReadGroup>::builder()
        .insert(rg_tag::SAMPLE, String::from("SampleA"))
        .insert(rg_tag::LIBRARY, String::from("LibB"))
        .insert(rg_tag::PLATFORM, String::from("ILLUMINA"))
        .build()
        .expect("valid RG2");

    Header::builder()
        .set_header(header_map)
        .add_reference_sequence(BString::from(ref_name), reference_sequence)
        .add_read_group(BString::from("RG1"), rg1)
        .add_read_group(BString::from("RG2"), rg2)
        .build()
}

/// Write grouped BAM file with a custom header that has multiple read groups.
fn create_grouped_bam_with_header(
    path: &Path,
    header: &Header,
    families: Vec<(&str, Vec<fgumi_raw_bam::RawRecord>)>,
) {
    let mut writer =
        bam::io::Writer::new(fs::File::create(path).expect("Failed to create BAM file"));
    writer.write_header(header).expect("Failed to write header");

    for (mi, records) in families {
        for raw in &records {
            use noodles::sam::alignment::record::data::field::Tag;
            use noodles::sam::alignment::record_buf::data::field::Value;
            let mut record = to_record_buf(raw);
            let mi_tag = Tag::from(fgumi_lib::sam::SamTag::MI);
            record.data_mut().insert(mi_tag, Value::from(mi));
            writer.write_alignment_record(header, &record).expect("Failed to write record");
        }
    }
    writer.finish(header).expect("Failed to finish BAM");
}

/// Test that simplex output header collapses read group attributes from input.
#[test]
fn test_simplex_command_collapses_read_group_attributes() {
    let temp_dir = TempDir::new().unwrap();
    let input_bam = temp_dir.path().join("input.bam");
    let output_bam = temp_dir.path().join("output.bam");

    // Create input BAM with two read groups having different LB but same SM
    let header = create_header_with_read_groups("chr1", 10000);
    let family = create_umi_family("ACGT", 5, "fam1", "ACGTACGT", 30);
    create_grouped_bam_with_header(&input_bam, &header, vec![("1", family)]);

    let cmd = Simplex::try_parse_from([
        "simplex",
        "--input",
        input_bam.to_str().unwrap(),
        "--output",
        output_bam.to_str().unwrap(),
        "--min-reads",
        "2",
        "--compression-level",
        "1",
    ])
    .expect("failed to parse simplex args");
    cmd.execute("fgumi simplex").expect("Simplex command failed");

    // Read the output header and verify collapsed read group attributes
    let mut reader = bam::io::Reader::new(fs::File::open(&output_bam).unwrap());
    let output_header = reader.read_header().unwrap();
    let read_groups = output_header.read_groups();

    assert_eq!(read_groups.len(), 1, "Should have exactly one output read group");

    let rg = read_groups.get(&BString::from("A")).expect("Output RG 'A' not found");

    // SM: both input RGs have "SampleA" -> deduplicated to single value
    assert_eq!(
        rg.other_fields().get(&rg_tag::SAMPLE).map(std::string::ToString::to_string),
        Some("SampleA".to_string()),
        "SM tag should be 'SampleA'",
    );

    // LB: "LibA" and "LibB" -> comma-joined
    assert_eq!(
        rg.other_fields().get(&rg_tag::LIBRARY).map(std::string::ToString::to_string),
        Some("LibA,LibB".to_string()),
        "LB tag should be comma-joined distinct values",
    );

    // PL: "illumina" and "ILLUMINA" -> uppercased and deduplicated
    assert_eq!(
        rg.other_fields().get(&rg_tag::PLATFORM).map(std::string::ToString::to_string),
        Some("ILLUMINA".to_string()),
        "PL tag should be uppercased and deduplicated",
    );
}

/// Verifies that the rejects BAM advertises the **input** header (RGs included),
/// while the primary output BAM advertises the **consensus** header (RGs collapsed).
///
/// This is the end-to-end check for the `secondary_output_header` parameter on
/// [`run_bam_pipeline_from_reader_with_secondary`][rbprws]. A regression where a
/// caller accidentally passes the primary `output_header` (consensus header) for
/// the secondary would collapse rejects' RGs to "A" — this test catches that.
///
/// [rbprws]: fgumi_lib::unified_pipeline::run_bam_pipeline_from_reader_with_secondary
#[test]
fn test_simplex_command_rejects_inherits_input_read_groups() {
    let temp_dir = TempDir::new().unwrap();
    let input_bam = temp_dir.path().join("input.bam");
    let output_bam = temp_dir.path().join("output.bam");
    let rejects_bam = temp_dir.path().join("rejects.bam");

    // Build an input with two distinct read groups (RG1, RG2) that the
    // consensus header would collapse to a single RG "A".
    let header = create_header_with_read_groups("chr1", 10000);
    // One family that passes consensus (kept), one singleton that doesn't (rejected).
    let kept = create_umi_family("ACGT", 5, "kept", "ACGTACGT", 30);
    let rejected = create_umi_family("TGCA", 1, "reject_singleton", "TTTTAAAA", 30);
    create_grouped_bam_with_header(&input_bam, &header, vec![("1", kept), ("2", rejected)]);

    let cmd = Simplex::try_parse_from([
        "simplex",
        "--input",
        input_bam.to_str().unwrap(),
        "--output",
        output_bam.to_str().unwrap(),
        "--rejects",
        rejects_bam.to_str().unwrap(),
        "--min-reads",
        "2",
        "--threads",
        "2",
        "--compression-level",
        "1",
    ])
    .expect("failed to parse simplex args");
    cmd.execute("fgumi simplex").expect("Simplex command with rejects failed");

    // Primary output BAM: consensus header collapses RG1+RG2 -> single "A".
    let primary_header =
        bam::io::Reader::new(fs::File::open(&output_bam).unwrap()).read_header().unwrap();
    let primary_rgs: Vec<BString> = primary_header.read_groups().keys().cloned().collect();
    assert_eq!(
        primary_rgs,
        vec![BString::from("A")],
        "primary output should have the collapsed consensus RG 'A'",
    );

    // Rejects BAM: input header preserved verbatim — both RG1 and RG2 present.
    let mut reader = bam::io::Reader::new(fs::File::open(&rejects_bam).unwrap());
    let rejects_header = reader.read_header().unwrap();
    let mut rejects_rgs: Vec<BString> = rejects_header.read_groups().keys().cloned().collect();
    rejects_rgs.sort();
    assert_eq!(
        rejects_rgs,
        vec![BString::from("RG1"), BString::from("RG2")],
        "rejects BAM should inherit the input header's read groups verbatim",
    );

    // Sanity-check that the rejection path actually fired: the singleton
    // input record (TGCA family, depth=1) must land in the rejects BAM.
    // Assert the read *name* (not just the count) so a regression that
    // rejects the wrong family (e.g. the kept "ACGT" family instead of the
    // singleton "TGCA" family) is caught — a count-only check would still
    // pass on a swap.
    let reject_names: Vec<String> = reader
        .records()
        .map(|result| {
            result
                .expect("Failed to read rejects record")
                .name()
                .expect("reject record missing read name")
                .to_string()
        })
        .collect();
    assert_eq!(
        reject_names.len(),
        1,
        "rejects BAM should contain the singleton record dropped by --min-reads",
    );
    assert!(
        reject_names[0].starts_with("reject_singleton"),
        "expected the singleton family in rejects, got {reject_names:?}",
    );
}
