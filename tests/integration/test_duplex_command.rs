//! End-to-end CLI tests for the duplex command.
//!
//! These tests invoke `Duplex::execute()` in-process and validate:
//! 1. Basic duplex consensus calling from paired-UMI grouped reads
//! 2. Statistics output
//! 3. Rejected reads output

use clap::Parser;
use fgumi_lib::commands::command::Command;
use fgumi_lib::commands::duplex::Duplex;
use fgumi_lib::sam::SamTag;
use fgumi_raw_bam::{RawRecord, RawRecordView, SamBuilder, flags};
use noodles::bam;
use noodles::sam::alignment::io::Write as AlignmentWrite;
use rstest::rstest;
use std::collections::HashMap;
use std::fs;
use std::path::{Path, PathBuf};
use tempfile::TempDir;

use crate::helpers::bam_generator::{create_minimal_header, to_record_buf};

/// Create a paired-end read pair for duplex consensus testing.
///
/// Duplex consensus requires reads grouped by MI tag with /A and /B strand suffixes.
/// For /A strand: R1 forward at pos, R2 reverse at pos+100 (FR orientation)
/// For /B strand: R1 reverse at pos+100, R2 forward at pos (RF orientation)
fn create_duplex_read_pair(
    name: &str,
    mi_tag: &str,
    sequence: &str,
    quality: u8,
    ref_start: i32,
    is_b_strand: bool,
) -> (RawRecord, RawRecord) {
    create_duplex_read_pair_with_sequences(
        name,
        mi_tag,
        sequence,
        sequence,
        quality,
        ref_start,
        is_b_strand,
    )
}

/// As [`create_duplex_read_pair`], but with distinct (equal-length) R1 and R2 sequences so
/// tests can tell which input read a consensus read was built from.
fn create_duplex_read_pair_with_sequences(
    name: &str,
    mi_tag: &str,
    r1_sequence: &str,
    r2_sequence: &str,
    quality: u8,
    ref_start: i32,
    is_b_strand: bool,
) -> (RawRecord, RawRecord) {
    let cigar_ops = [u32::try_from(r1_sequence.len()).expect("read_len fits u32") << 4];
    create_duplex_read_pair_with_cigar(
        name,
        mi_tag,
        r1_sequence,
        r2_sequence,
        quality,
        ref_start,
        is_b_strand,
        &cigar_ops,
    )
}

/// As [`create_duplex_read_pair_with_sequences`], but with a caller-supplied CIGAR on both
/// reads. The other constructors derive an all-`M` CIGAR from the sequence length, so they
/// cannot express the differing alignment patterns the minority-alignment filter keys on.
#[expect(clippy::too_many_arguments, reason = "test fixture builder mirrors the BAM fields")]
fn create_duplex_read_pair_with_cigar(
    name: &str,
    mi_tag: &str,
    r1_sequence: &str,
    r2_sequence: &str,
    quality: u8,
    ref_start: i32,
    is_b_strand: bool,
    cigar_ops: &[u32],
) -> (RawRecord, RawRecord) {
    assert_eq!(r1_sequence.len(), r2_sequence.len(), "R1 and R2 must be the same length");
    let seq = r1_sequence.as_bytes();
    let read_len = seq.len();

    let (r1_start, r2_start, r1_rev, r2_rev) = if is_b_strand {
        // B strand: RF orientation — R1 reverse at far position, R2 forward at near
        (ref_start + 100, ref_start, true, false)
    } else {
        // A strand: FR orientation — R1 forward at near position, R2 reverse at far
        (ref_start, ref_start + 100, false, true)
    };

    #[expect(
        clippy::cast_possible_truncation,
        clippy::cast_possible_wrap,
        reason = "test data with known small values"
    )]
    let tlen: i32 = if is_b_strand { -((read_len + 100) as i32) } else { (read_len + 100) as i32 };

    let r1_flags = flags::PAIRED
        | flags::FIRST_SEGMENT
        | if r1_rev { flags::REVERSE } else { 0 }
        | if r2_rev { flags::MATE_REVERSE } else { 0 };

    let r2_flags = flags::PAIRED
        | flags::LAST_SEGMENT
        | if r2_rev { flags::REVERSE } else { 0 }
        | if r1_rev { flags::MATE_REVERSE } else { 0 };

    let r1 = {
        let mut b = SamBuilder::new();
        b.read_name(name.as_bytes())
            .sequence(seq)
            .qualities(&vec![quality; read_len])
            .flags(r1_flags)
            .ref_id(0)
            .pos(r1_start - 1)
            .mapq(60)
            .cigar_ops(cigar_ops)
            .mate_ref_id(0)
            .mate_pos(r2_start - 1)
            .template_length(tlen)
            .add_string_tag(SamTag::MI, mi_tag.as_bytes());
        b.build()
    };

    let r2 = {
        let mut b = SamBuilder::new();
        b.read_name(name.as_bytes())
            .sequence(r2_sequence.as_bytes())
            .qualities(&vec![quality; read_len])
            .flags(r2_flags)
            .ref_id(0)
            .pos(r2_start - 1)
            .mapq(60)
            .cigar_ops(cigar_ops)
            .mate_ref_id(0)
            .mate_pos(r1_start - 1)
            .template_length(-tlen)
            .add_string_tag(SamTag::MI, mi_tag.as_bytes());
        b.build()
    };

    (r1, r2)
}

/// Create a BAM with duplex-grouped reads (MI tags with /A and /B strand suffixes).
pub(crate) fn create_duplex_bam(path: &Path, molecules: Vec<Vec<(RawRecord, RawRecord)>>) {
    create_duplex_bam_with_header(path, &create_minimal_header("chr1", 10000), molecules);
}

/// Create a BAM with duplex-grouped reads using a caller-supplied header, so tests can
/// pin behavior that depends on header content (e.g. `@RG` collapsing).
fn create_duplex_bam_with_header(
    path: &Path,
    header: &noodles::sam::Header,
    molecules: Vec<Vec<(RawRecord, RawRecord)>>,
) {
    let mut writer =
        bam::io::Writer::new(fs::File::create(path).expect("Failed to create BAM file"));
    writer.write_header(header).expect("Failed to write header");

    for pairs in molecules {
        for (r1, r2) in pairs {
            writer.write_alignment_record(header, &to_record_buf(&r1)).expect("Failed to write R1");
            writer.write_alignment_record(header, &to_record_buf(&r2)).expect("Failed to write R2");
        }
    }
    writer.try_finish().expect("Failed to finish BAM");
}

/// Create a single duplex molecule with `depth` read pairs on each strand.
pub(crate) fn create_duplex_molecule(
    mi_id: &str,
    sequence: &str,
    quality: u8,
    ref_start: i32,
    depth: usize,
) -> Vec<(RawRecord, RawRecord)> {
    let mut molecule = Vec::new();
    for i in 0..depth {
        let (r1, r2) = create_duplex_read_pair(
            &format!("ab_{i}"),
            &format!("{mi_id}/A"),
            sequence,
            quality,
            ref_start,
            false,
        );
        molecule.push((r1, r2));
    }
    for i in 0..depth {
        let (r1, r2) = create_duplex_read_pair(
            &format!("ba_{i}"),
            &format!("{mi_id}/B"),
            sequence,
            quality,
            ref_start,
            true,
        );
        molecule.push((r1, r2));
    }
    molecule
}

/// Test basic duplex consensus calling.
#[test]
fn test_duplex_command_basic() {
    let temp_dir = TempDir::new().unwrap();
    let input_bam = temp_dir.path().join("input.bam");
    let output_bam = temp_dir.path().join("output.bam");

    let molecule = create_duplex_molecule("1", "ACGTACGT", 30, 100, 3);
    create_duplex_bam(&input_bam, vec![molecule]);

    let cmd = Duplex::try_parse_from([
        "duplex",
        "--input",
        input_bam.to_str().unwrap(),
        "--output",
        output_bam.to_str().unwrap(),
        "--min-reads",
        "1",
        "--compression-level",
        "1",
    ])
    .expect("failed to parse duplex args");
    cmd.execute("fgumi duplex").expect("Duplex command failed");
    assert!(output_bam.exists(), "Output BAM not created");

    // Verify consensus reads were produced
    let mut reader = bam::io::Reader::new(fs::File::open(&output_bam).unwrap());
    let _header = reader.read_header().unwrap();
    let mut count = 0;
    for result in reader.records() {
        let record = result.expect("Failed to read record");
        let cd_tag = SamTag::CD.to_noodles_tag();
        assert!(record.data().get(&cd_tag).is_some(), "Duplex consensus should have cD tag");
        count += 1;
    }
    assert!(count > 0, "Should have produced duplex consensus reads");
}

/// Test duplex command with rejects output.
///
/// Runs the multi-threaded pipeline with `--min-reads 2` and a singleton /A
/// template that fails the single-strand `min-reads` check. Verifies that:
///
/// 1. The rejects BAM is created and its `@HD` sort fields (`SO`/`GO`/`SS`)
///    match the input BAM, because rejects flow through the unified
///    pipeline's secondary output in batch/input order.
/// 2. The singleton template's raw records are actually streamed to the
///    rejects BAM rather than being silently dropped.
#[test]
fn test_duplex_command_with_rejects() {
    let temp_dir = TempDir::new().unwrap();
    let input_bam = temp_dir.path().join("input.bam");
    let output_bam = temp_dir.path().join("output.bam");
    let rejects_bam = temp_dir.path().join("rejects.bam");

    // Kept molecule: 3 /A pairs + 3 /B pairs
    let kept = create_duplex_molecule("1", "ACGTACGT", 30, 100, 3);
    // Singleton /A pair (fails the /A strand's min-reads=2 check)
    let singleton = vec![create_duplex_read_pair("solo", "2/A", "ACGTACGT", 30, 500, false)];
    create_duplex_bam(&input_bam, vec![kept, singleton]);

    let cmd = Duplex::try_parse_from([
        "duplex",
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
    .expect("failed to parse duplex args");
    cmd.execute("fgumi duplex").expect("Duplex command with rejects failed");
    assert!(rejects_bam.exists(), "Rejects BAM not created");

    crate::helpers::assertions::assert_rejects_header_matches_input(&rejects_bam, &input_bam);

    // The singleton /A template (R1 + R2) should be streamed to the rejects BAM.
    let mut reader = bam::io::Reader::new(fs::File::open(&rejects_bam).unwrap());
    let _header = reader.read_header().unwrap();
    let mut reject_count = 0;
    for result in reader.records() {
        result.expect("Failed to read reject record");
        reject_count += 1;
    }
    assert_eq!(
        reject_count, 2,
        "Singleton paired-end /A template should be streamed to rejects BAM (R1 + R2)"
    );
}

/// Regression test: each rejected input record must appear in the rejects BAM
/// exactly once.
///
/// The duplex caller accumulates rejects from two sources: records dropped at
/// the single-strand (ss) layer, and records dropped at the duplex layer. When
/// the duplex layer rejects a group it returns every raw input record — which
/// already includes anything the ss layer would have rejected. Appending both
/// sources naively would emit the overlapping records twice.
///
/// The molecules below all reject at the duplex layer, so this pins the
/// deduplication of the *whole-group* source across several independent
/// rejection paths: the rejects BAM must contain each `(read_name, flags)`
/// tuple at most once, with the total record count matching the expected
/// input-record count. The single-strand source, and its reconciliation with
/// `--stats`, is covered by
/// `test_duplex_rejects_bam_reconciles_with_stats`.
#[test]
fn test_duplex_command_rejects_contain_no_duplicates() {
    let temp_dir = TempDir::new().unwrap();
    let input_bam = temp_dir.path().join("input.bam");
    let output_bam = temp_dir.path().join("output.bam");
    let rejects_bam = temp_dir.path().join("rejects.bam");

    // Kept molecule: produces a consensus, contributes zero rejects.
    let kept = create_duplex_molecule("1", "ACGTACGT", 30, 100, 3);
    // Four singleton molecules (2 records each) — two on each strand — that
    // reject via the insufficient-reads path. Using molecule IDs in a single
    // `singletons` vec avoids the similar_names clippy lint.
    let singletons: Vec<Vec<(RawRecord, RawRecord)>> = vec![
        vec![create_duplex_read_pair("soloA1", "2/A", "ACGTACGT", 30, 500, false)],
        vec![create_duplex_read_pair("soloA2", "3/A", "ACGTACGT", 30, 800, false)],
        vec![create_duplex_read_pair("soloB1", "4/B", "ACGTACGT", 30, 1100, true)],
        vec![create_duplex_read_pair("soloB2", "5/B", "ACGTACGT", 30, 1400, true)],
    ];

    // Expected rejects = 4 singletons × 2 records = 8. The kept molecule
    // contributes 0 rejects.
    let expected_rejects = u32::try_from(singletons.len() * 2).expect("rejects count fits in u32");

    let mut molecules = vec![kept];
    molecules.extend(singletons);
    create_duplex_bam(&input_bam, molecules);

    let cmd = Duplex::try_parse_from([
        "duplex",
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
    .expect("failed to parse duplex args");
    cmd.execute("fgumi duplex").expect("Duplex command with rejects failed");
    assert!(rejects_bam.exists(), "Rejects BAM not created");

    crate::helpers::assertions::assert_has_bgzf_eof(&rejects_bam);
    crate::helpers::assertions::assert_rejects_header_matches_input(&rejects_bam, &input_bam);

    let mut reader = bam::io::Reader::new(fs::File::open(&rejects_bam).unwrap());
    let _header = reader.read_header().unwrap();
    let mut seen: std::collections::HashMap<(String, u16), u32> = std::collections::HashMap::new();
    let mut total = 0u32;
    for result in reader.records() {
        let record = result.expect("Failed to read reject record");
        let name = record.name().expect("reject record missing read name").to_string();
        let flags = u16::from(record.flags());
        *seen.entry((name, flags)).or_insert(0) += 1;
        total += 1;
    }

    assert_eq!(
        total, expected_rejects,
        "rejects BAM record count should match the expected rejected-input count",
    );
    for ((name, flags), count) in &seen {
        assert_eq!(
            *count, 1,
            "record ({name}, flags={flags}) appears {count} times in the rejects BAM",
        );
    }
}

/// #757: the rejects BAM must reconcile with `--stats`.
///
/// The molecule here emits a duplex consensus *and* drops reads: `minority` carries a
/// minority indel pattern on both ends, so the single-strand layer drops it from both
/// combined alignment groups (AB-R1 + BA-R2 and AB-R2 + BA-R1) while the three majority
/// templates per strand still consense. Those drops used to be recorded on the composed
/// single-strand caller, whose statistics the duplex caller never merged and whose rejected
/// records it never received, so they reached neither output — `raw_reads_rejected`
/// under-reported and the rejects BAM was empty for a molecule that plainly rejected reads.
///
/// Asserted in both the single-threaded and pipeline paths, since they drain the caller's
/// rejects and statistics through different code.
#[rstest]
#[case::single_threaded(None)]
#[case::threaded(Some("2"))]
fn test_duplex_rejects_bam_reconciles_with_stats(#[case] threads: Option<&str>) {
    let temp_dir = TempDir::new().unwrap();
    let input_bam = temp_dir.path().join("input.bam");
    let output_bam = temp_dir.path().join("output.bam");
    let rejects_bam = temp_dir.path().join("rejects.bam");
    let stats_path = temp_dir.path().join("stats.txt");

    // 3M1D5M against the majority 8M: same query length, different alignment pattern.
    let gapped = [3u32 << 4, (1u32 << 4) | 2, 5u32 << 4];
    let mut molecule = create_duplex_molecule("1", "ACGTACGT", 30, 100, 3);
    molecule.push(create_duplex_read_pair_with_cigar(
        "minority", "1/A", "ACGTACGT", "ACGTACGT", 30, 100, false, &gapped,
    ));
    create_duplex_bam(&input_bam, vec![molecule]);

    let mut args = vec![
        "duplex",
        "--input",
        input_bam.to_str().unwrap(),
        "--output",
        output_bam.to_str().unwrap(),
        "--rejects",
        rejects_bam.to_str().unwrap(),
        "--stats",
        stats_path.to_str().unwrap(),
        "--min-reads",
        "1",
        "--compression-level",
        "1",
    ];
    if let Some(threads) = threads {
        args.extend_from_slice(&["--threads", threads]);
    }
    Duplex::try_parse_from(args)
        .expect("failed to parse duplex args")
        .execute("fgumi duplex")
        .expect("Duplex command failed");

    let mut output_reader = bam::io::Reader::new(fs::File::open(&output_bam).unwrap());
    output_reader.read_header().expect("read output header");
    assert_eq!(
        output_reader.records().count(),
        2,
        "the majority-alignment reads must still emit a duplex R1/R2 pair"
    );

    // Key the input records by (name, flags) so each reject can be compared against the
    // exact bytes it came from: a rejects BAM that dropped a tag or rewrote a CIGAR would
    // still satisfy a name-and-count assertion.
    let expected_rejects: HashMap<(Vec<u8>, u16), Vec<u8>> = read_raw_records(&input_bam)
        .into_iter()
        .filter(|record| fgumi_raw_bam::read_name(record) == b"minority")
        .map(|record| {
            let view = RawRecordView::new(&record);
            ((view.read_name().to_vec(), view.flags()), record)
        })
        .collect();
    assert_eq!(expected_rejects.len(), 2, "the minority template contributes two records");

    let mut seen: Vec<(Vec<u8>, u16)> = Vec::new();
    for record in read_raw_records(&rejects_bam) {
        let view = RawRecordView::new(&record);
        let key = (view.read_name().to_vec(), view.flags());
        let expected = expected_rejects.get(&key).unwrap_or_else(|| {
            panic!("unexpected reject record {}", String::from_utf8_lossy(&key.0))
        });
        assert_eq!(
            &record, expected,
            "each reject must be byte-for-byte identical to its input record — every field \
             and tag preserved on the --rejects path"
        );
        assert!(!seen.contains(&key), "a reject record was written more than once");
        seen.push(key);
    }

    let stats = fs::read_to_string(&stats_path).expect("read stats");
    let stat = |key: &str| -> usize {
        stats
            .lines()
            .find_map(|line| line.strip_prefix(&format!("{key}\t")))
            .and_then(|rest| rest.split('\t').next())
            .unwrap_or_else(|| panic!("stats must contain a {key} row"))
            .parse()
            .unwrap_or_else(|_| panic!("{key} must be an integer"))
    };
    assert_eq!(
        seen.len(),
        expected_rejects.len(),
        "the rejects BAM must hold both ends of the minority template"
    );
    assert_eq!(
        stat("raw_reads_rejected_for_minority_alignment"),
        2,
        "both ends of the minority template must be counted as minority-alignment"
    );
    assert_eq!(
        seen.len(),
        stat("raw_reads_rejected"),
        "the rejects BAM record count must equal raw_reads_rejected"
    );
}

/// #792: reads that trim to zero length must reach both the rejects BAM and `--stats`.
///
/// The duplex path built its source reads with a bare `filter_map`, silently discarding every
/// read `create_source_read` rejected (a read that quality-trims to zero) — so it reached neither
/// output. fgbio rejects such a read as `ZeroPostAfterTrimming` on the duplex caller's own writer
/// and counter (`UmiConsensusCaller.toSourceRead`), one layer, both outputs; the vanilla path
/// already matches that, and the duplex path now must too.
///
/// The `trimmed` template's two reads are entirely below `--min-input-base-quality`, so with
/// `--trim` they trim to zero length while the three majority templates per strand still consense.
/// Asserted single-threaded and threaded, since they drain rejects and statistics through
/// different code.
#[rstest]
#[case::single_threaded(None)]
#[case::threaded(Some("2"))]
fn test_duplex_zero_length_rejects_reconcile_with_stats(#[case] threads: Option<&str>) {
    let temp_dir = TempDir::new().unwrap();
    let input_bam = temp_dir.path().join("input.bam");
    let output_bam = temp_dir.path().join("output.bam");
    let rejects_bam = temp_dir.path().join("rejects.bam");
    let stats_path = temp_dir.path().join("stats.txt");

    let mut molecule = create_duplex_molecule("1", "ACGTACGT", 30, 100, 3);
    // One extra /A template whose bases are all below --min-input-base-quality (10), so with
    // --trim it trims to zero length and create_source_read rejects it.
    molecule.push(create_duplex_read_pair_with_cigar(
        "trimmed",
        "1/A",
        "ACGTACGT",
        "ACGTACGT",
        2,
        100,
        false,
        &[8u32 << 4],
    ));
    create_duplex_bam(&input_bam, vec![molecule]);

    let mut args = vec![
        "duplex",
        "--input",
        input_bam.to_str().unwrap(),
        "--output",
        output_bam.to_str().unwrap(),
        "--rejects",
        rejects_bam.to_str().unwrap(),
        "--stats",
        stats_path.to_str().unwrap(),
        "--min-reads",
        "1",
        "--trim",
        "--compression-level",
        "1",
    ];
    if let Some(threads) = threads {
        args.extend_from_slice(&["--threads", threads]);
    }
    Duplex::try_parse_from(args)
        .expect("failed to parse duplex args")
        .execute("fgumi duplex")
        .expect("Duplex command failed");

    let mut output_reader = bam::io::Reader::new(fs::File::open(&output_bam).unwrap());
    output_reader.read_header().expect("read output header");
    assert_eq!(
        output_reader.records().count(),
        2,
        "the majority-alignment reads must still emit a duplex R1/R2 pair"
    );

    // Key input records by (name, flags) so each reject is compared byte-for-byte against the
    // exact bytes it came from.
    let expected_rejects: HashMap<(Vec<u8>, u16), Vec<u8>> = read_raw_records(&input_bam)
        .into_iter()
        .filter(|record| fgumi_raw_bam::read_name(record) == b"trimmed")
        .map(|record| {
            let view = RawRecordView::new(&record);
            ((view.read_name().to_vec(), view.flags()), record)
        })
        .collect();
    assert_eq!(expected_rejects.len(), 2, "the trimmed template contributes two records");

    let mut seen: Vec<(Vec<u8>, u16)> = Vec::new();
    for record in read_raw_records(&rejects_bam) {
        let view = RawRecordView::new(&record);
        let key = (view.read_name().to_vec(), view.flags());
        let expected = expected_rejects.get(&key).unwrap_or_else(|| {
            panic!("unexpected reject record {}", String::from_utf8_lossy(&key.0))
        });
        assert_eq!(
            &record, expected,
            "each reject must be byte-for-byte identical to its input record"
        );
        assert!(!seen.contains(&key), "a reject record was written more than once");
        seen.push(key);
    }

    let stats = fs::read_to_string(&stats_path).expect("read stats");
    let stat = |key: &str| -> usize {
        stats
            .lines()
            .find_map(|line| line.strip_prefix(&format!("{key}\t")))
            .and_then(|rest| rest.split('\t').next())
            .unwrap_or_else(|| panic!("stats must contain a {key} row"))
            .parse()
            .unwrap_or_else(|_| panic!("{key} must be an integer"))
    };
    assert_eq!(
        seen.len(),
        expected_rejects.len(),
        "the rejects BAM must hold both ends of the trimmed template"
    );
    assert_eq!(
        stat("raw_reads_rejected_for_zero_bases_post_trimming"),
        2,
        "both ends of the trimmed template must be counted as zero-bases-post-trimming"
    );
    assert_eq!(
        seen.len(),
        stat("raw_reads_rejected"),
        "the rejects BAM record count must equal raw_reads_rejected"
    );
}

/// #792 (follow-up): zero-length rejects must reach the rejects BAM in input order.
///
/// A read that trims to zero length is dropped by `create_source_read` while its X/Y source
/// vectors are built, before alignment filtering. For a `/B` template that means R1 lands in the
/// Y partition and R2 in the X partition, so collecting the dropped raws X-partition-first would
/// write R2 ahead of R1 — the reverse of input order, which wrote R1 then R2. The zero-length
/// path must merge both partitions back into the group's canonical order, exactly as the
/// alignment-filter reject path does.
///
/// Mirrors `test_duplex_rejects_preserve_input_order_for_b_strand_minority`, but the minority
/// template is rejected for trimming to zero length rather than for a divergent alignment, so it
/// exercises the separate zero-length collection path. Asserted single-threaded and threaded,
/// since they drain rejects through different code.
#[rstest]
#[case::single_threaded(None)]
#[case::threaded(Some("2"))]
fn test_duplex_zero_length_rejects_preserve_input_order_for_b_strand(
    #[case] threads: Option<&str>,
) {
    let temp_dir = TempDir::new().unwrap();
    let input_bam = temp_dir.path().join("input.bam");
    let output_bam = temp_dir.path().join("output.bam");
    let rejects_bam = temp_dir.path().join("rejects.bam");

    let mut molecule = create_duplex_molecule("1", "ACGTACGT", 30, 100, 3);
    // Minority template on the *B* strand whose bases are all below --min-input-base-quality (10):
    // with --trim both ends trim to zero length and `create_source_read` rejects them. Its R1
    // becomes BA-R1 (lands in Y) and its R2 becomes BA-R2 (lands in X), so the two ends are
    // dropped by separate partitions — exactly the split that can invert their order.
    molecule.push(create_duplex_read_pair_with_cigar(
        "trimmed",
        "1/B",
        "ACGTACGT",
        "ACGTACGT",
        2,
        100,
        true,
        &[8u32 << 4],
    ));
    create_duplex_bam(&input_bam, vec![molecule]);

    let mut args = vec![
        "duplex",
        "--input",
        input_bam.to_str().unwrap(),
        "--output",
        output_bam.to_str().unwrap(),
        "--rejects",
        rejects_bam.to_str().unwrap(),
        "--min-reads",
        "1",
        "--trim",
        "--compression-level",
        "1",
    ];
    if let Some(threads) = threads {
        args.extend_from_slice(&["--threads", threads]);
    }
    Duplex::try_parse_from(args)
        .expect("failed to parse duplex args")
        .execute("fgumi duplex")
        .expect("Duplex command failed");

    // The trimmed template's two input records, in the exact order the input wrote them (R1 then
    // R2). The majority templates all consense, so they contribute no rejects and the rejects BAM
    // must reproduce this sequence and nothing else.
    let expected_sequence: Vec<Vec<u8>> = read_raw_records(&input_bam)
        .into_iter()
        .filter(|record| fgumi_raw_bam::read_name(record) == b"trimmed")
        .collect();
    assert_eq!(expected_sequence.len(), 2, "the trimmed template contributes two records");

    let reject_sequence = read_raw_records(&rejects_bam);
    assert_eq!(
        reject_sequence, expected_sequence,
        "the /B trimmed template's zero-length rejects must be serialized in input order (R1 then \
         R2), not X-before-Y order (R2 then R1)"
    );
}

/// #757 (follow-up): single-strand rejects must reach the rejects BAM in input order.
///
/// The two combined alignment groups split a template's ends across strands: X = AB-R1 +
/// BA-R2 and Y = AB-R2 + BA-R1. For a `/B` template that means X carries its R2 while Y
/// carries its R1, so recording X's rejects before Y's would write R2 ahead of R1 — the
/// reverse of input order, which wrote R1 then R2. The reconciliation must merge both
/// strands' rejects back into the group's canonical order before appending them.
///
/// Mirrors `test_duplex_rejects_bam_reconciles_with_stats`, but puts the minority template on
/// the `/B` strand and asserts the serialized reject *sequence* equals the input sequence,
/// which a name-and-count check (as the reconcile test uses) cannot catch.
#[rstest]
#[case::single_threaded(None)]
#[case::threaded(Some("2"))]
fn test_duplex_rejects_preserve_input_order_for_b_strand_minority(#[case] threads: Option<&str>) {
    let temp_dir = TempDir::new().unwrap();
    let input_bam = temp_dir.path().join("input.bam");
    let output_bam = temp_dir.path().join("output.bam");
    let rejects_bam = temp_dir.path().join("rejects.bam");

    // 3M1D5M against the majority 8M: same query length, minority alignment pattern.
    let gapped = [3u32 << 4, (1u32 << 4) | 2, 5u32 << 4];
    let mut molecule = create_duplex_molecule("1", "ACGTACGT", 30, 100, 3);
    // Minority template on the *B* strand: its R1 becomes BA-R1 (lands in Y) and its R2
    // becomes BA-R2 (lands in X), so the two ends are dropped by separate
    // `filter_by_alignment` calls — exactly the split that can invert their order.
    molecule.push(create_duplex_read_pair_with_cigar(
        "minority", "1/B", "ACGTACGT", "ACGTACGT", 30, 100, true, &gapped,
    ));
    create_duplex_bam(&input_bam, vec![molecule]);

    let mut args = vec![
        "duplex",
        "--input",
        input_bam.to_str().unwrap(),
        "--output",
        output_bam.to_str().unwrap(),
        "--rejects",
        rejects_bam.to_str().unwrap(),
        "--min-reads",
        "1",
        "--compression-level",
        "1",
    ];
    if let Some(threads) = threads {
        args.extend_from_slice(&["--threads", threads]);
    }
    Duplex::try_parse_from(args)
        .expect("failed to parse duplex args")
        .execute("fgumi duplex")
        .expect("Duplex command failed");

    // The minority template's two input records, in the exact order the input wrote them
    // (R1 then R2). The majority templates all consense, so they contribute no rejects and
    // the rejects BAM must reproduce this sequence and nothing else.
    let expected_sequence: Vec<Vec<u8>> = read_raw_records(&input_bam)
        .into_iter()
        .filter(|record| fgumi_raw_bam::read_name(record) == b"minority")
        .collect();
    assert_eq!(expected_sequence.len(), 2, "the minority template contributes two records");

    let reject_sequence = read_raw_records(&rejects_bam);
    assert_eq!(
        reject_sequence, expected_sequence,
        "the /B minority template's rejects must be serialized in input order (R1 then R2), \
         not X-before-Y order (R2 then R1)"
    );
}

/// Reads every record of a BAM back as raw bytes, so tests can compare whole serialized
/// records rather than only the fields noodles decodes ergonomically.
fn read_raw_records(path: &Path) -> Vec<Vec<u8>> {
    let (mut reader, _header) =
        fgumi_bam_io::create_raw_bam_reader(path, 1).expect("open BAM for raw reading");
    let mut records = Vec::new();
    let mut record = RawRecord::new();
    while reader.read_record(&mut record).expect("read raw record") != 0 {
        records.push(record.as_ref().to_vec());
    }
    records
}

/// Test duplex command with statistics output.
#[test]
fn test_duplex_command_with_stats() {
    let temp_dir = TempDir::new().unwrap();
    let input_bam = temp_dir.path().join("input.bam");
    let output_bam = temp_dir.path().join("output.bam");
    let stats_path = temp_dir.path().join("stats.txt");

    let molecule = create_duplex_molecule("1", "ACGTACGT", 30, 100, 3);
    create_duplex_bam(&input_bam, vec![molecule]);

    let cmd = Duplex::try_parse_from([
        "duplex",
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
    .expect("failed to parse duplex args");
    cmd.execute("fgumi duplex").expect("Duplex command with stats failed");
    assert!(stats_path.exists(), "Stats file not created");

    // Single-threaded duplex must emit the same fgbio seeded key-value format as the
    // multi-threaded path (thread-independent metrics): the always-seeded `usedByDuplex`
    // rejection rows must be present even at zero. The wide (per-field) format the
    // single-threaded path previously wrote would not contain these KV keys.
    let stats = std::fs::read_to_string(&stats_path).expect("read stats");
    for key in [
        "raw_reads_rejected_for_non_paired_reads",
        "raw_reads_rejected_for_single_strand_only",
        "raw_reads_rejected_for_potential_umi_collision",
    ] {
        assert!(
            stats.contains(key),
            "single-threaded duplex stats must emit the seeded KV row {key}:\n{stats}"
        );
    }
}

/// Build a molecule observed on a single strand only: `depth` read pairs all carrying
/// the same `/A` or `/B` MI suffix.
fn create_single_strand_molecule(
    mi_id: &str,
    sequence: &str,
    quality: u8,
    ref_start: i32,
    depth: usize,
    is_b_strand: bool,
) -> Vec<(RawRecord, RawRecord)> {
    let strand = if is_b_strand { 'B' } else { 'A' };
    (0..depth)
        .map(|i| {
            create_duplex_read_pair(
                &format!("{strand}_{i}"),
                &format!("{mi_id}/{strand}"),
                sequence,
                quality,
                ref_start,
                is_b_strand,
            )
        })
        .collect()
}

/// A third `--min-reads` value of 0 is fgbio's single-strand mode: molecules seen on only
/// one strand still yield a consensus. With a non-zero third value they are rejected.
///
/// The first two `--min-reads` values still apply in single-strand mode, and they count
/// TEMPLATES, not read records: `has_minimum_number_of_reads` counts only the R1 of each
/// pair, matching fgbio's `x.count(r => r.paired && r.firstOfPair)`. The threshold cases
/// below straddle the molecule's template count to pin that — see `MOLECULE_DEPTH`.
#[rstest]
fn test_duplex_single_strand_mode_emits_one_strand_molecules(
    #[values(None, Some("2"))] threads: Option<&str>,
    #[values(true, false)] is_b_strand: bool,
) {
    /// Read pairs on the molecule's single strand: 3 templates, i.e. 6 read records.
    const MOLECULE_DEPTH: usize = 3;

    let temp_dir = TempDir::new().unwrap();
    let input_bam = temp_dir.path().join("input.bam");

    let molecule =
        create_single_strand_molecule("1", "ACGTACGT", 30, 100, MOLECULE_DEPTH, is_b_strand);
    create_duplex_bam(&input_bam, vec![molecule]);

    // Identity of one consensus record: `(is_r1, molecule id from the read name, bases)`. A
    // bare record count cannot tell a correct R1/R2 pair from two empty, duplicated, or
    // wrong-molecule records, so the assertions below work off this instead.
    let run = |min_reads: &str, output_bam: &Path| -> Vec<(bool, String, String)> {
        let mut args = vec![
            "duplex",
            "--input",
            input_bam.to_str().unwrap(),
            "--output",
            output_bam.to_str().unwrap(),
            "--min-reads",
            min_reads,
            "--compression-level",
            "1",
        ];
        if let Some(threads) = threads {
            args.extend_from_slice(&["--threads", threads]);
        }
        Duplex::try_parse_from(args)
            .expect("failed to parse duplex args")
            .execute("fgumi duplex")
            .expect("Duplex command failed");
        let mut reader = bam::io::Reader::new(fs::File::open(output_bam).unwrap());
        let _header = reader.read_header().unwrap();
        reader
            .records()
            .map(|result| {
                let record = result.expect("failed to read consensus record");
                let name = String::from_utf8(
                    record.name().expect("consensus record must have a read name").to_vec(),
                )
                .expect("consensus read name must be UTF-8");
                // Consensus reads are named `<read-group-prefix>:<base MI>`, so the molecule
                // identity survives in the name even though the `/A`-`/B` suffix does not.
                let (_prefix, molecule) = name.rsplit_once(':').unwrap_or_else(|| {
                    panic!("consensus read name '{name}' must be '<prefix>:<MI>'")
                });
                let flags = record.flags();
                assert!(flags.is_segmented(), "consensus record '{name}' must be paired");
                assert!(
                    flags.is_first_segment() != flags.is_last_segment(),
                    "consensus record '{name}' must be exactly one of R1/R2"
                );
                let bases: String = record.sequence().iter().map(char::from).collect();
                assert_eq!(
                    record.quality_scores().as_ref().len(),
                    bases.len(),
                    "consensus record '{name}' must carry one quality per base"
                );
                (flags.is_first_segment(), molecule.to_string(), bases)
            })
            .collect()
    };

    // The molecule's reads are `sequence` long, and a single-strand consensus over identical
    // Q30 reads reproduces those bases (reverse-complemented on whichever end the strand puts
    // in reverse orientation), so neither end may come back empty or N-masked.
    let sequence = "ACGTACGT";
    let expected_bases = [sequence.to_string(), reverse_complement(sequence)];

    // The whole contract for an emitted pair: two records, exactly one R1, both carrying the
    // input molecule's id and reproducing its bases. Every accepting case below asserts this,
    // so none of them can pass on a bare record count.
    let assert_is_the_molecules_pair = |reads: &[(bool, String, String)], context: &str| {
        assert_eq!(reads.len(), 2, "{context}: expected an R1/R2 consensus pair, got {reads:?}");
        assert_eq!(
            reads.iter().filter(|(is_r1, _, _)| *is_r1).count(),
            1,
            "{context}: the pair must be exactly one consensus R1 and one consensus R2, got \
             {reads:?}"
        );
        for (_, molecule, bases) in reads {
            assert_eq!(
                molecule, "1",
                "{context}: consensus must carry the input molecule's id, got {molecule}"
            );
            assert!(
                expected_bases.contains(bases),
                "{context}: consensus bases {bases:?} must reproduce the molecule's reads (one \
                 of {expected_bases:?})"
            );
        }
    };

    assert_is_the_molecules_pair(
        &run("1,1,0", &temp_dir.path().join("single-strand.bam")),
        "a one-strand molecule with a third --min-reads value of 0",
    );

    assert!(
        run("1", &temp_dir.path().join("both-strands.bam")).is_empty(),
        "a one-strand molecule must still be rejected when both strands are required"
    );

    // Template-count boundary. The molecule is MOLECULE_DEPTH templates but twice that many
    // read records, so a threshold in between separates the two readings: at exactly
    // MOLECULE_DEPTH the molecule must still be emitted, and one above it must be rejected.
    // An implementation counting records instead of templates would emit in both cases.
    let at_depth = format!("{MOLECULE_DEPTH},{MOLECULE_DEPTH},0");
    assert_is_the_molecules_pair(
        &run(&at_depth, &temp_dir.path().join("at-template-count.bam")),
        &format!("--min-reads {at_depth} equals the molecule's {MOLECULE_DEPTH} templates"),
    );
    let above_depth = format!("{},{},0", MOLECULE_DEPTH + 1, MOLECULE_DEPTH + 1);
    assert!(
        run(&above_depth, &temp_dir.path().join("above-template-count.bam")).is_empty(),
        "--min-reads {above_depth} exceeds the molecule's {MOLECULE_DEPTH} templates, so no \
         consensus may be emitted -- the {} read records must not count toward the threshold",
        MOLECULE_DEPTH * 2
    );
}

/// An out-of-order `--min-reads` is an argument error, not a mid-run failure.
///
/// Both execution paths must reject it before touching the output. That is only
/// automatic single-threaded, where the caller is constructed up front; `--threads`
/// mode builds its caller inside the pipeline's per-batch closure, i.e. after the
/// output BAM has been created, so without an up-front check the same bad value
/// would surface as a wrapped worker I/O error over a partial file.
#[rstest]
fn test_duplex_rejects_out_of_order_min_reads_before_writing_output(
    #[values(None, Some("2"))] threads: Option<&str>,
    #[values("1,2,3", "3,1,2")] min_reads: &str,
) {
    let temp_dir = TempDir::new().unwrap();
    let input_bam = temp_dir.path().join("input.bam");
    let output_bam = temp_dir.path().join("output.bam");

    create_duplex_bam(
        &input_bam,
        vec![create_single_strand_molecule("1", "ACGTACGT", 30, 100, 3, false)],
    );

    let mut args = vec![
        "duplex",
        "--input",
        input_bam.to_str().unwrap(),
        "--output",
        output_bam.to_str().unwrap(),
        "--min-reads",
        min_reads,
        "--compression-level",
        "1",
    ];
    if let Some(threads) = threads {
        args.extend_from_slice(&["--threads", threads]);
    }

    let error = Duplex::try_parse_from(args)
        .expect("failed to parse duplex args")
        .execute("fgumi duplex")
        .expect_err(&format!("--min-reads {min_reads} must be rejected"));
    let message = format!("{error:#}");
    assert!(
        message.contains("min-reads values must be specified high to low"),
        "--min-reads {min_reads} must fail as an argument error, got: {message}"
    );
    assert!(
        !output_bam.exists(),
        "--min-reads {min_reads} must be rejected before the output BAM is created, but \
         {} was left behind",
        output_bam.display()
    );
}

/// Reverse complement, for deriving the source-molecule orientation of a reverse-strand read.
fn reverse_complement(sequence: &str) -> String {
    sequence
        .chars()
        .rev()
        .map(|base| match base {
            'A' => 'T',
            'C' => 'G',
            'G' => 'C',
            'T' => 'A',
            other => other,
        })
        .collect()
}

/// In single-strand mode a consensus read must be built from the input read of the same
/// end, for a `/B`-only molecule as much as for a `/A`-only one: consensus R1 from the
/// molecule's R1s and consensus R2 from its R2s, each in source-molecule orientation.
///
/// This matches fgbio, which treats whichever strand group is present as the "AB" group.
#[rstest]
fn test_duplex_single_strand_consensus_reads_follow_input_read_ends(
    #[values(true, false)] is_b_strand: bool,
    #[values(None, Some("2"))] threads: Option<&str>,
) {
    let temp_dir = TempDir::new().unwrap();
    let input_bam = temp_dir.path().join("input.bam");
    let output_bam = temp_dir.path().join("output.bam");

    // Deliberately not reverse complements of each other, so that "built from R1" and
    // "built from R2" cannot produce the same bases.
    let (r1_seq, r2_seq) = ("AAAACCCC", "ACACACAC");
    let strand = if is_b_strand { 'B' } else { 'A' };
    let molecule = (0..3)
        .map(|i| {
            create_duplex_read_pair_with_sequences(
                &format!("{strand}_{i}"),
                &format!("1/{strand}"),
                r1_seq,
                r2_seq,
                30,
                100,
                is_b_strand,
            )
        })
        .collect();
    create_duplex_bam(&input_bam, vec![molecule]);

    let mut args = vec![
        "duplex",
        "--input",
        input_bam.to_str().unwrap(),
        "--output",
        output_bam.to_str().unwrap(),
        "--min-reads",
        "1,1,0",
        "--compression-level",
        "1",
    ];
    if let Some(threads) = threads {
        args.extend_from_slice(&["--threads", threads]);
    }
    Duplex::try_parse_from(args)
        .expect("failed to parse duplex args")
        .execute("fgumi duplex")
        .expect("Duplex command failed");

    // `/A` reads are R1-forward/R2-reverse and `/B` reads the reverse, so the source-molecule
    // orientation of each end flips with the strand.
    let (expected_r1, expected_r2) = if is_b_strand {
        (reverse_complement(r1_seq), r2_seq.to_string())
    } else {
        (r1_seq.to_string(), reverse_complement(r2_seq))
    };

    let mut reader = bam::io::Reader::new(fs::File::open(&output_bam).unwrap());
    let _header = reader.read_header().unwrap();
    let mut seen_r1 = None;
    let mut seen_r2 = None;
    for result in reader.records() {
        let record = result.expect("Failed to read record");
        let bases: String = record.sequence().iter().map(char::from).collect();
        if record.flags().is_first_segment() {
            seen_r1 = Some(bases);
        } else {
            seen_r2 = Some(bases);
        }
    }

    assert_eq!(
        seen_r1.as_deref(),
        Some(expected_r1.as_str()),
        "consensus R1 must be built from the molecule's R1 reads"
    );
    assert_eq!(
        seen_r2.as_deref(),
        Some(expected_r2.as_str()),
        "consensus R2 must be built from the molecule's R2 reads"
    );
}

/// Build a `/A` read pair whose mates align over exactly the same reference bases, so the
/// overlapping-bases consensus has something to correct.
fn create_fully_overlapping_pair(
    name: &str,
    mi_tag: &str,
    r1_sequence: &str,
    r2_sequence: &str,
    quality: u8,
    ref_start: i32,
) -> (RawRecord, RawRecord) {
    assert_eq!(r1_sequence.len(), r2_sequence.len(), "R1 and R2 must be the same length");
    let read_len = r1_sequence.len();
    let cigar_op = u32::try_from(read_len).expect("read_len fits u32") << 4;
    let build = |sequence: &str, read_flags: u16| {
        let mut b = SamBuilder::new();
        b.read_name(name.as_bytes())
            .sequence(sequence.as_bytes())
            .qualities(&vec![quality; read_len])
            .flags(read_flags)
            .ref_id(0)
            .pos(ref_start - 1)
            .mapq(60)
            .cigar_ops(&[cigar_op])
            .mate_ref_id(0)
            .mate_pos(ref_start - 1)
            .template_length(0)
            .add_string_tag(SamTag::MI, mi_tag.as_bytes());
        b.build()
    };
    (
        build(r1_sequence, flags::PAIRED | flags::FIRST_SEGMENT | flags::MATE_REVERSE),
        build(r2_sequence, flags::PAIRED | flags::LAST_SEGMENT | flags::REVERSE),
    )
}

/// R1's stored bases in the overlap fixture below.
const OVERLAP_R1_BASES: &str = "ACGTACGT";
/// R2's stored bases in the overlap fixture below. BAM `SEQ` is stored in reference
/// orientation, so these ARE R2's reference bases: they disagree with `OVERLAP_R1_BASES`
/// at the LAST reference base (offset 7), which is the base the overlapping consensus
/// masks. R2 is `REVERSE`, so its sequencing orientation is the reverse complement
/// `TCGTACGT` — what the uncorrected consensus emits for R2.
const OVERLAP_R2_BASES: &str = "ACGTACGA";

/// Write the single-strand overlapping-mates fixture: one `/A` pair whose mates align over
/// exactly the same reference bases at equal quality but disagree at one of them.
///
/// One BAM builder serves both the assertion test and the `#[ignore]`d regen test below,
/// so the bytes the offline fgbio run consumes are the bytes the command under test reads.
fn write_overlapping_single_strand_fixture(path: &Path) {
    let molecule = vec![create_fully_overlapping_pair(
        "ab_0",
        "1/A",
        OVERLAP_R1_BASES,
        OVERLAP_R2_BASES,
        30,
        100,
    )];
    create_duplex_bam(path, vec![molecule]);
}

/// Re-emit the single-strand overlapping-mates fixture BAM for a fresh fgbio capture.
///
/// Mirrors `fgumi_consensus`'s `regen_write` helper — writes to `FGUMI_ORACLE_BAM_OUT` when
/// set, otherwise to a logged temp path — which that crate keeps `pub(crate)`, so it cannot
/// be shared with this integration test.
#[test]
#[ignore = "regen: writes the fixture BAM (set FGUMI_ORACLE_BAM_OUT), then run fgbio on it"]
fn regen_duplex_single_strand_overlap_fixture() {
    let path = std::env::var_os("FGUMI_ORACLE_BAM_OUT").map_or_else(
        || {
            let (_file, path) = tempfile::NamedTempFile::new()
                .expect("temp fixture BAM")
                .keep()
                .expect("persist temp fixture BAM");
            path
        },
        PathBuf::from,
    );
    write_overlapping_single_strand_fixture(&path);
    eprintln!("wrote fixture BAM to {}", path.display());
}

/// The overlapping-bases consensus is applied to single-strand molecules too.
///
/// It is skipped for groups without both strands on the grounds that no duplex consensus
/// can be built from them — which stops being true once the third `--min-reads` value is 0.
///
/// Both the corrected and uncorrected outputs are pinned to a real fgbio run on this exact
/// fixture rather than to a "contains an N" predicate, so a masking rule that differs from
/// fgbio in WHICH base it masks, or in what it leaves behind, is caught. Oracle provenance
/// — fgbio `4.1.0`, on the BAM `regen_duplex_single_strand_overlap_fixture` emits:
///
/// ```text
/// fgbio CallDuplexConsensusReads -i <fixture>.bam -o out.bam -M 1 1 0 \
///   --consensus-call-overlapping-bases {true,false}
/// samtools view out.bam   # SEQ of the R1 (0x40) and R2 (0x80) records
/// ```
#[rstest]
fn test_duplex_single_strand_mode_applies_overlapping_consensus(
    #[values(None, Some("2"))] threads: Option<&str>,
) {
    let temp_dir = TempDir::new().unwrap();
    let input_bam = temp_dir.path().join("input.bam");

    write_overlapping_single_strand_fixture(&input_bam);

    // Each consensus record as `(is_r1, bases, qualities)`. Qualities are part of the pinned
    // identity: the overlapping consensus combines the mates' evidence where they AGREE,
    // which raises the consensus quality, so bases alone would not catch a correction that
    // masked the right base but mis-weighted the rest.
    let consensus_reads = |overlapping: &str, output_bam: &Path| -> Vec<(bool, String, String)> {
        let mut args = vec![
            "duplex",
            "--input",
            input_bam.to_str().unwrap(),
            "--output",
            output_bam.to_str().unwrap(),
            "--min-reads",
            "1,1,0",
            "--consensus-call-overlapping-bases",
            overlapping,
            "--compression-level",
            "1",
        ];
        if let Some(threads) = threads {
            args.extend_from_slice(&["--threads", threads]);
        }
        Duplex::try_parse_from(args)
            .expect("failed to parse duplex args")
            .execute("fgumi duplex")
            .expect("Duplex command failed");

        let mut reader = bam::io::Reader::new(fs::File::open(output_bam).unwrap());
        let _header = reader.read_header().unwrap();
        let mut reads: Vec<(bool, String, String)> = reader
            .records()
            .map(|result| {
                let record = result.expect("Failed to read record");
                let bases: String = record.sequence().iter().map(char::from).collect();
                // Render qualities as the printable SAM characters the captured fgbio
                // `samtools view` output shows, so expected and actual read alike.
                let quals: String =
                    record.quality_scores().as_ref().iter().map(|&q| char::from(q + 33)).collect();
                (record.flags().is_first_segment(), bases, quals)
            })
            .collect();
        reads.sort_by_key(|(is_r1, _, _)| !*is_r1);
        reads
    };

    // fgbio's captured output. The mates disagree at the LAST reference base, which the
    // overlapping consensus masks to N at qual 2 in both mates. In sequencing orientation
    // that lands at the end of R1 (forward) and at the start of R2 (reverse), and fgbio
    // treats the two ends differently: the leading no-call survives as an `N` with depth 0,
    // while the trailing one is trimmed, leaving R1 one base shorter than R2. Confirmed by
    // running the same fgbio command on a fixture that disagrees at the FIRST reference
    // base instead, which swaps the two records' shapes exactly.
    let corrected = consensus_reads("true", &temp_dir.path().join("corrected.bam"));
    assert_eq!(
        corrected,
        vec![
            (true, "ACGTACG".to_string(), "HHHHHHH".to_string()),
            (false, "NCGTACGT".to_string(), "#HHHHHHH".to_string()),
        ],
        "corrected single-strand output must match the captured fgbio 4.1.0 run"
    );

    // Without the correction each mate is consensus-called on its own, so both keep all
    // eight bases -- R2 reverse-complemented into sequencing orientation -- and the quality
    // stays at what one Q30 read supports rather than the combined value above.
    let uncorrected = consensus_reads("false", &temp_dir.path().join("uncorrected.bam"));
    assert_eq!(
        uncorrected,
        vec![
            (true, "ACGTACGT".to_string(), ">>>>>>>>".to_string()),
            (false, "TCGTACGT".to_string(), ">>>>>>>>".to_string()),
        ],
        "uncorrected single-strand output must match the captured fgbio 4.1.0 run"
    );
}

/// The reported `consensus_reads_emitted` must equal the number of consensus records
/// actually written, in both the single-threaded and pipeline (`--threads N`) paths.
///
/// The duplex caller emits R1 and R2 per molecule, so a metric that counts templates
/// rather than reads reports half the records in the output BAM.
#[rstest]
#[case::single_threaded(None)]
#[case::threaded(Some("2"))]
fn test_duplex_stats_consensus_reads_matches_records_written(#[case] threads: Option<&str>) {
    let temp_dir = TempDir::new().unwrap();
    let input_bam = temp_dir.path().join("input.bam");
    let output_bam = temp_dir.path().join("output.bam");
    let stats_path = temp_dir.path().join("stats.txt");

    let molecule_ids: Vec<String> = (1..=3).map(|i: i32| i.to_string()).collect();
    let molecules = (1..=3)
        .map(|i| create_duplex_molecule(&i.to_string(), "ACGTACGT", 30, 100 + i * 1000, 3))
        .collect();
    create_duplex_bam(&input_bam, molecules);

    let mut args = vec![
        "duplex",
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
    ];
    if let Some(threads) = threads {
        args.extend_from_slice(&["--threads", threads]);
    }
    Duplex::try_parse_from(args)
        .expect("failed to parse duplex args")
        .execute("fgumi duplex")
        .expect("Duplex command failed");

    // Independent oracle: every input molecule must yield exactly one consensus R1 and one
    // consensus R2, so the output holds 2 records per molecule and each molecule's identity
    // survives. A bare `records_written > 0` would accept a dropped mate, a dropped molecule,
    // or a duplicated record whenever `consensus_reads_emitted` is wrong in the same way.
    let mut reader = bam::io::Reader::new(fs::File::open(&output_bam).unwrap());
    let _header = reader.read_header().unwrap();
    let mut records_written = 0_usize;
    // Per molecule: (R1 count, R2 count), keyed by the molecule ID trailing the consensus
    // read name (`<read-group-prefix>:<base MI>`).
    let mut ends_by_molecule: HashMap<String, (usize, usize)> = HashMap::new();
    for result in reader.records() {
        let record = result.expect("failed to read consensus record");
        records_written += 1;
        let name = String::from_utf8(
            record.name().expect("consensus record must have a read name").to_vec(),
        )
        .expect("consensus read name must be UTF-8");
        let (_prefix, molecule) = name
            .rsplit_once(':')
            .unwrap_or_else(|| panic!("consensus read name '{name}' must be '<prefix>:<MI>'"));
        let flags = record.flags();
        assert!(flags.is_segmented(), "consensus record '{name}' must be paired");
        let ends = ends_by_molecule.entry(molecule.to_string()).or_default();
        if flags.is_first_segment() {
            ends.0 += 1;
        } else if flags.is_last_segment() {
            ends.1 += 1;
        } else {
            panic!("consensus record '{name}' is neither R1 nor R2");
        }
    }

    assert_eq!(
        records_written,
        2 * molecule_ids.len(),
        "each of the {} input molecules must emit exactly one consensus R1 and one R2",
        molecule_ids.len()
    );
    for mi in &molecule_ids {
        assert_eq!(
            ends_by_molecule.get(mi.as_str()).copied(),
            Some((1, 1)),
            "molecule {mi} must emit exactly one consensus R1 and one consensus R2, \
             got {:?}",
            ends_by_molecule.get(mi.as_str())
        );
    }

    let stats = fs::read_to_string(&stats_path).expect("read stats");
    let emitted: usize = stats
        .lines()
        .find_map(|line| line.strip_prefix("consensus_reads_emitted\t"))
        .and_then(|rest| rest.split('\t').next())
        .expect("stats must contain a consensus_reads_emitted row")
        .parse()
        .expect("consensus_reads_emitted must be an integer");

    assert_eq!(
        emitted, records_written,
        "consensus_reads_emitted must count consensus reads written to the output BAM"
    );
}

/// Verifies that the duplex consensus output BAM advertises an *unmapped consensus*
/// header, identically in the single-threaded (`--threads` unset) and pipeline
/// (`--threads N`) paths.
///
/// fgbio's `UmiConsensusCaller.outputHeader` builds a brand-new `SAMFileHeader`:
/// no sequence dictionary, the single collapsed `@RG`, `SO:unsorted` + `GO:query`,
/// and `SamOrder.Unsorted.applyTo` explicitly *clears* `SS`. Consensus reads are
/// unmapped, so carrying the input's `@SQ` dictionary and its
/// `SS:unsorted:template-coordinate` sub-sort forward is wrong twice over: the
/// records reference no contig, and the header falsely advertises a
/// template-coordinate order that downstream sort-order preconditions
/// (`check_consensus_sort_order`) accept without warning.
///
/// `simplex` and `codec` build the consensus header in both threading paths; this
/// test pins `duplex` to the same contract — @SQ, @HD, @RG, and the exact fgbio
/// provenance @CO — and checks the emitted records honour it (unmapped, carrying the
/// declared @RG).
#[rstest]
#[case::single_threaded(None)]
#[case::pipeline(Some("2"))]
fn test_duplex_command_emits_unmapped_consensus_header(#[case] threads: Option<&str>) {
    use bstr::BString;
    use noodles::sam::alignment::record::data::field::Tag;
    use noodles::sam::alignment::record_buf::data::field::Value;
    use noodles::sam::header::record::value::Map;
    use noodles::sam::header::record::value::map::ReadGroup;
    use noodles::sam::header::record::value::map::header::tag as header_tag;
    use noodles::sam::header::record::value::map::read_group::tag as rg_tag;

    // Two input read groups with distinct SM/LB/PL/PU/CN/DS values. fgbio collapses each
    // tag's distinct values into the single output @RG (comma-joined in read-group order,
    // with PL upper-cased); the expectations below are fixed and independent of the
    // implementation so a collapse that drops, reorders, or mangles a tag fails here.
    let input_read_groups: &[(&str, &[(_, &str)])] = &[
        (
            "rgA",
            &[
                (rg_tag::SAMPLE, "sampleA"),
                (rg_tag::LIBRARY, "libA"),
                (rg_tag::PLATFORM, "illumina"),
                (rg_tag::PLATFORM_UNIT, "unitA"),
                (rg_tag::SEQUENCING_CENTER, "centerA"),
                (rg_tag::DESCRIPTION, "descA"),
            ],
        ),
        (
            "rgB",
            &[
                (rg_tag::SAMPLE, "sampleB"),
                (rg_tag::LIBRARY, "libB"),
                (rg_tag::PLATFORM, "iontorrent"),
                (rg_tag::PLATFORM_UNIT, "unitB"),
                (rg_tag::SEQUENCING_CENTER, "centerB"),
                (rg_tag::DESCRIPTION, "descB"),
            ],
        ),
    ];
    let expected_collapsed_rg: &[(_, &str)] = &[
        (rg_tag::SAMPLE, "sampleA,sampleB"),
        (rg_tag::LIBRARY, "libA,libB"),
        (rg_tag::PLATFORM, "ILLUMINA,IONTORRENT"),
        (rg_tag::PLATFORM_UNIT, "unitA,unitB"),
        (rg_tag::SEQUENCING_CENTER, "centerA,centerB"),
        (rg_tag::DESCRIPTION, "descA,descB"),
    ];

    let temp_dir = TempDir::new().unwrap();
    let input_bam = temp_dir.path().join("input.bam");
    let output_bam = temp_dir.path().join("output.bam");

    let mut input_header = create_minimal_header("chr1", 10000);
    for (id, tags) in input_read_groups {
        let mut rg = Map::<ReadGroup>::builder();
        for (tag, value) in *tags {
            rg = rg.insert(*tag, (*value).to_string());
        }
        input_header
            .read_groups_mut()
            .insert(BString::from(*id), rg.build().expect("valid read group"));
    }

    let molecule = create_duplex_molecule("1", "ACGTACGT", 30, 100, 3);
    create_duplex_bam_with_header(&input_bam, &input_header, vec![molecule]);

    let mut args = vec![
        "duplex",
        "--input",
        input_bam.to_str().unwrap(),
        "--output",
        output_bam.to_str().unwrap(),
        "--min-reads",
        "1",
        "--compression-level",
        "1",
    ];
    if let Some(threads) = threads {
        args.extend_from_slice(&["--threads", threads]);
    }
    let cmd = Duplex::try_parse_from(args).expect("failed to parse duplex args");
    let read_group_id = cmd.read_group.read_group_id.clone();
    cmd.execute("fgumi duplex").expect("Duplex command failed");

    let mut input_reader = bam::io::Reader::new(fs::File::open(&input_bam).unwrap());
    let input_read_group_count = input_reader.read_header().unwrap().read_groups().len();

    let mut reader = bam::io::Reader::new(fs::File::open(&output_bam).unwrap());
    let output_header = reader.read_header().unwrap();

    // Consensus reads are unmapped: the input's sequence dictionary must be dropped.
    assert!(
        output_header.reference_sequences().is_empty(),
        "consensus output header must not carry @SQ lines, found {:?}",
        output_header.reference_sequences().keys().collect::<Vec<_>>(),
    );

    let hd = output_header.header().expect("output header must have an @HD line");
    assert_eq!(
        hd.other_fields().get(&header_tag::SORT_ORDER).map(std::string::ToString::to_string),
        Some("unsorted".to_string()),
        "consensus output must be SO:unsorted",
    );
    assert_eq!(
        hd.other_fields().get(&header_tag::GROUP_ORDER).map(std::string::ToString::to_string),
        Some("query".to_string()),
        "consensus output must be GO:query",
    );
    assert_eq!(
        hd.other_fields().get(&header_tag::SUBSORT_ORDER).map(std::string::ToString::to_string),
        None,
        "consensus output must not advertise an SS sub-sort (fgbio clears it)",
    );

    // The records carry RG:Z:<read_group_id>, so the header must declare exactly that
    // read group and no other, with the collapsed SM/LB/PL/PU/CN/DS attributes fgbio
    // derives from the input read groups.
    let read_groups = output_header.read_groups();
    assert_eq!(read_groups.len(), 1, "consensus output must declare exactly one @RG");
    let output_rg = read_groups
        .get(&BString::from(read_group_id.as_str()))
        .expect("consensus output must declare the @RG the consensus records reference");
    for (tag, expected) in expected_collapsed_rg {
        assert_eq!(
            output_rg.other_fields().get(tag).map(std::string::ToString::to_string).as_deref(),
            Some(*expected),
            "collapsed @RG must carry the expected {tag} value",
        );
    }
    // No stray attributes beyond the collapsed set fgbio produces.
    assert_eq!(
        output_rg.other_fields().len(),
        expected_collapsed_rg.len(),
        "collapsed @RG must carry exactly the collapsed SM/LB/PL/PU/CN/DS attributes, found {:?}",
        output_rg.other_fields().keys().collect::<Vec<_>>(),
    );

    // fgbio records provenance as a single @CO naming the collapsed read group and the
    // number of input read groups it was built from. Assert the whole comment set against
    // an expectation derived from the input header, so a malformed comment or a wrong
    // input-@RG count fails here rather than passing a substring check.
    let expected_comment = format!(
        "Read group {read_group_id} contains consensus reads generated from \
         {input_read_group_count} input read groups."
    );
    assert_eq!(
        output_header.comments().iter().map(std::string::ToString::to_string).collect::<Vec<_>>(),
        vec![expected_comment],
        "consensus output must carry exactly the fgbio-style @CO provenance comment",
    );

    // The header contract only holds if the records honour it: an @SQ-free header cannot
    // describe a mapped record, and the lone @RG must be the one the records reference.
    let records: Vec<_> = reader
        .record_bufs(&output_header)
        .map(|result| result.expect("failed to read consensus record"))
        .collect();
    assert!(!records.is_empty(), "duplex must emit at least one consensus record");
    for record in &records {
        assert!(
            record.flags().is_unmapped(),
            "consensus records must be unmapped, found flags {:?}",
            record.flags(),
        );
        let rg = record
            .data()
            .get(&Tag::from(SamTag::RG))
            .expect("consensus record must carry an RG tag");
        match rg {
            Value::String(id) => assert_eq!(
                String::from_utf8_lossy(id),
                read_group_id,
                "consensus records must reference the @RG the header declares",
            ),
            other => panic!("RG tag must be a string, found {other:?}"),
        }
    }
}

/// Builds an A-strand pair whose R2 is unmapped, matching the geometry of
/// [`create_duplex_read_pair_with_sequences`] so it lands in the same molecule.
///
/// The R2 keeps its REVERSE flag so strand grouping is unaffected; only the mapping is removed.
fn create_half_mapped_a_strand_pair(
    name: &str,
    mi_tag: &str,
    r1_sequence: &str,
    r2_sequence: &str,
    quality: u8,
    ref_start: i32,
) -> (RawRecord, RawRecord) {
    let read_len = r1_sequence.len();
    let cigar_op = u32::try_from(read_len).expect("read_len fits u32") << 4;
    #[expect(
        clippy::cast_possible_truncation,
        clippy::cast_possible_wrap,
        reason = "test data with known small values"
    )]
    let tlen: i32 = (read_len + 100) as i32;

    let r1 = {
        let mut b = SamBuilder::new();
        b.read_name(name.as_bytes())
            .sequence(r1_sequence.as_bytes())
            .qualities(&vec![quality; read_len])
            .flags(flags::PAIRED | flags::FIRST_SEGMENT | flags::MATE_REVERSE)
            .ref_id(0)
            .pos(ref_start - 1)
            .mapq(60)
            .cigar_ops(&[cigar_op])
            .mate_ref_id(0)
            .mate_pos(ref_start + 100 - 1)
            .template_length(tlen)
            .add_string_tag(SamTag::MI, mi_tag.as_bytes());
        b.build()
    };

    let r2 = {
        let mut b = SamBuilder::new();
        b.read_name(name.as_bytes())
            .sequence(r2_sequence.as_bytes())
            .qualities(&vec![quality; read_len])
            .flags(flags::PAIRED | flags::LAST_SEGMENT | flags::REVERSE | flags::UNMAPPED)
            .ref_id(0)
            .pos(ref_start + 100 - 1)
            .mapq(0)
            .mate_ref_id(0)
            .mate_pos(ref_start - 1)
            .template_length(-tlen)
            .add_string_tag(SamTag::MI, mi_tag.as_bytes());
        b.build()
    };

    (r1, r2)
}

/// The unmapped end of a half-mapped pair must not contribute to a single-strand consensus when
/// mapped reads are present on that strand. fgumi always prefers mapped reads.
#[test]
fn test_duplex_ignores_unmapped_end_when_mapped_reads_present() {
    let temp_dir = TempDir::new().unwrap();
    let input_bam = temp_dir.path().join("input.bam");
    let output_bam = temp_dir.path().join("output.bam");

    let molecule = vec![
        create_duplex_read_pair_with_sequences(
            "ab1",
            "1/A",
            "ACGTACGTAC",
            "GGGGGGGGGG",
            30,
            100,
            false,
        ),
        create_duplex_read_pair_with_sequences(
            "ab2",
            "1/A",
            "ACGTACGTAC",
            "GGGGGGGGGG",
            30,
            100,
            false,
        ),
        create_half_mapped_a_strand_pair("ab3", "1/A", "ACGTACGTAC", "TTTTTTTTTT", 30, 100),
        create_duplex_read_pair_with_sequences(
            "ba1",
            "1/B",
            "ACGTACGTAC",
            "GGGGGGGGGG",
            30,
            100,
            true,
        ),
        create_duplex_read_pair_with_sequences(
            "ba2",
            "1/B",
            "ACGTACGTAC",
            "GGGGGGGGGG",
            30,
            100,
            true,
        ),
    ];
    create_duplex_bam(&input_bam, vec![molecule]);

    let cmd = Duplex::try_parse_from([
        "duplex",
        "--input",
        input_bam.to_str().unwrap(),
        "--output",
        output_bam.to_str().unwrap(),
        "--min-reads",
        "1",
        "--compression-level",
        "1",
    ])
    .expect("failed to parse duplex args");
    cmd.execute("fgumi duplex").expect("Duplex command failed");

    let mut reader = bam::io::Reader::new(fs::File::open(&output_bam).unwrap());
    let _header = reader.read_header().unwrap();
    let records: Vec<bam::Record> =
        reader.records().map(|r| r.expect("failed to read record")).collect();
    assert_eq!(records.len(), 2, "expected one duplex consensus read pair");

    let depth = |record: &bam::Record, tag: SamTag| -> i64 {
        match record.data().get(&tag.to_noodles_tag()) {
            Some(Ok(value)) => value.as_int().expect("depth tag must be an integer"),
            other => panic!("record must carry an integer {tag:?} tag, got {other:?}"),
        }
    };

    let r1 = records
        .iter()
        .find(|r| r.flags().bits() & flags::FIRST_SEGMENT != 0)
        .expect("missing R1 consensus");
    let r2 = records
        .iter()
        .find(|r| r.flags().bits() & flags::LAST_SEGMENT != 0)
        .expect("missing R2 consensus");

    // All three AB-R1s are mapped, so all three contribute to R1's AB single-strand consensus.
    assert_eq!(depth(r1, SamTag::AD), 3, "R1 should use all three mapped AB reads");

    // Only the two mapped AB-R2s may contribute to R2's AB single-strand consensus. AD is the
    // oracle here rather than the emitted bases: the record's SEQ is the *duplex* consensus of
    // the AB and BA strands, so it does not expose the AB single-strand bases this test is about.
    // (The simplex counterpart can assert bases because it has only one strand.)
    assert_eq!(depth(r2, SamTag::AD), 2, "R2 must use only the two mapped AB reads");
}
