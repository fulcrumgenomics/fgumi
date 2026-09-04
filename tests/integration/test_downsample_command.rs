//! Integration tests for the downsample command.
//!
//! These tests invoke the downsample `Command::execute()` in-process rather than spawning
//! the actual `fgumi downsample` binary.

use clap::Parser;
use fgumi_lib::commands::command::Command;
use fgumi_lib::commands::downsample::Downsample;
use fgumi_lib::sam::SamTag;
use fgumi_raw_bam::SamBuilder;
use noodles::bam;
use noodles::sam::alignment::io::Write as AlignmentWrite;
use noodles::sam::alignment::record::data::field::Tag;
use noodles::sam::alignment::record_buf::data::field::Value;
use std::collections::{HashMap, HashSet};
use std::fs;
use std::path::PathBuf;
use tempfile::TempDir;

use crate::helpers::bam_generator::create_minimal_header;

/// MI tag constant
const MI_TAG: Tag = SamTag::MI.to_noodles_tag();

/// Create a grouped BAM file with MI tags (simulating output from group).
fn create_grouped_bam(path: &PathBuf, families: Vec<(&str, usize)>) {
    let header = create_minimal_header("chr1", 10000);

    let mut writer =
        bam::io::Writer::new(fs::File::create(path).expect("Failed to create BAM file"));

    writer.write_header(&header).expect("Failed to write header");

    // Create records grouped by MI tag
    let mut read_idx = 0;
    for (mi, count) in families {
        for _ in 0..count {
            let raw = {
                let mut b = SamBuilder::new();
                b.read_name(format!("read_{read_idx}").as_bytes())
                    .sequence(b"ACGT")
                    .qualities(&[30; 4])
                    .ref_id(0)
                    .pos(99) // 0-based (alignment_start 100)
                    .mapq(60)
                    .cigar_ops(&[4u32 << 4]) // 4M
                    .add_string_tag(SamTag::MI, mi.as_bytes());
                b.build()
            };
            let record =
                fgumi_raw_bam::raw_record_to_record_buf(&raw, &noodles::sam::Header::default())
                    .expect("raw_record_to_record_buf failed in test");

            writer.write_alignment_record(&header, &record).expect("Failed to write record");
            read_idx += 1;
        }
    }

    writer.try_finish().expect("Failed to finish BAM");
}

/// Flip a byte in the last BGZF block's CRC32 footer, so decoding that block
/// fails only when CRC verification is on. Adapted from PR2's dedup integration
/// helper. Requires the file to span at least two blocks so the corrupted block
/// comes *after* reader construction succeeds: the single-threaded raw reader
/// (`FgumiBgzfReader`) applies `verify_crc` uniformly to every block, so
/// corrupting block 0 (the header's block) with verification on would fail while
/// building the reader, before the intended read/count assertion can run.
fn corrupt_last_block_crc(path: &std::path::Path) {
    let mut bytes = fs::read(path).expect("read bam for corruption");
    let mut cursor: &[u8] = &bytes;
    let blocks = fgumi_lib::bgzf_reader::read_raw_blocks(&mut cursor, 100_000)
        .expect("read bgzf blocks from test bam");
    assert!(
        blocks.len() >= 2,
        "test input must span >= 2 BGZF blocks so the corrupted block isn't the header's; \
         got {} -- generate more records",
        blocks.len()
    );
    let offset: usize =
        blocks[..blocks.len() - 1].iter().map(fgumi_lib::bgzf_reader::RawBgzfBlock::len).sum();
    let last = blocks.last().expect("checked len >= 2 above");
    // `read_raw_blocks` drops every BGZF EOF marker, so summing the returned
    // (real) block lengths yields the last block's on-disk offset only when no
    // marker sits *between* real blocks. Guard that: everything past the last
    // framed block must be whole trailing EOF markers (a writer may emit more
    // than one). An intermediate marker would leave real data here instead and
    // shift `crc_off` onto an unrelated byte, which the `>= 2` guard cannot
    // detect.
    let eof = &fgumi_lib::bgzf_reader::BGZF_EOF;
    let tail = &bytes[offset + last.len()..];
    assert!(
        tail.len().is_multiple_of(eof.len())
            && tail.chunks_exact(eof.len()).all(|chunk| chunk == &eof[..]),
        "bytes after the last framed block must be only trailing BGZF EOF markers; \
         an intermediate marker would invalidate the CRC offset"
    );
    let crc_off = offset + last.len() - fgumi_lib::bgzf_reader::BGZF_FOOTER_SIZE;
    bytes[crc_off] ^= 0x01;
    fs::write(path, bytes).expect("write corrupted bam");
}

/// `--no-check-crc` must let downsample's single-threaded raw reader accept a
/// corrupted BGZF CRC32 (it decodes through fgumi-bgzf, honoring the flag, #800).
#[test]
fn test_downsample_no_check_crc_accepts_corrupted_crc() {
    let temp_dir = TempDir::new().expect("Failed to create temp dir");
    let input_bam = temp_dir.path().join("input.bam");
    let output_bam = temp_dir.path().join("output.bam");

    // Enough records to span multiple BGZF blocks (see corrupt_last_block_crc).
    create_grouped_bam(&input_bam, vec![("MI1", 3000)]);
    corrupt_last_block_crc(&input_bam);

    let cmd = Downsample::try_parse_from([
        "downsample",
        "-i",
        input_bam.to_str().unwrap(),
        "-o",
        output_bam.to_str().unwrap(),
        "-f",
        "1.0",
        "--seed",
        "42",
        "--no-check-crc",
    ])
    .expect("failed to parse downsample args");
    cmd.execute("fgumi downsample")
        .expect("--no-check-crc must accept a corrupted BGZF CRC32 and complete");

    // Assert the corrupted block was decoded, not silently dropped: the input is
    // one MI family of 3000 records at -f 1.0 (every record kept), so all 3000
    // must survive the --no-check-crc decode. Pin the complete read-name set,
    // not just the count — a bare `records().count()` can tally error items and
    // passes even if records were dropped and replaced or duplicated.
    let mut reader =
        bam::io::reader::Builder.build_from_path(&output_bam).expect("open output BAM");
    reader.read_header().expect("read output BAM header");
    let records: Vec<_> = reader
        .records()
        .collect::<std::io::Result<Vec<_>>>()
        .expect("every output record must decode under --no-check-crc");
    let mut names: Vec<String> = records
        .iter()
        .map(|r| {
            let name = r.name().expect("record has a name");
            String::from_utf8_lossy(AsRef::<[u8]>::as_ref(&name)).into_owned()
        })
        .collect();
    names.sort();
    let mut expected: Vec<String> = (0..3000).map(|i| format!("read_{i}")).collect();
    expected.sort();
    assert_eq!(
        names, expected,
        "read_0..read_2999 must all survive the --no-check-crc decode of the corrupted block, \
         by identity"
    );
}

/// Default (verify-on for file input) must reject the same corrupted CRC32.
#[test]
fn test_downsample_rejects_corrupted_crc_by_default() {
    let temp_dir = TempDir::new().expect("Failed to create temp dir");
    let input_bam = temp_dir.path().join("input.bam");
    let output_bam = temp_dir.path().join("output.bam");

    create_grouped_bam(&input_bam, vec![("MI1", 3000)]);
    corrupt_last_block_crc(&input_bam);

    let cmd = Downsample::try_parse_from([
        "downsample",
        "-i",
        input_bam.to_str().unwrap(),
        "-o",
        output_bam.to_str().unwrap(),
        "-f",
        "1.0",
        "--seed",
        "42",
    ])
    .expect("failed to parse downsample args");
    let err = cmd
        .execute("fgumi downsample")
        .expect_err("default (verify-on for file input) must reject a corrupted BGZF CRC32");
    let message = format!("{err:#}");
    assert!(message.to_uppercase().contains("CRC32"), "error should mention CRC32: {message}");
}

/// `--check-crc` must also reject the corrupted CRC32 (forces verification on).
#[test]
fn test_downsample_check_crc_rejects_corrupted_crc() {
    let temp_dir = TempDir::new().expect("Failed to create temp dir");
    let input_bam = temp_dir.path().join("input.bam");
    let output_bam = temp_dir.path().join("output.bam");

    create_grouped_bam(&input_bam, vec![("MI1", 3000)]);
    corrupt_last_block_crc(&input_bam);

    let cmd = Downsample::try_parse_from([
        "downsample",
        "-i",
        input_bam.to_str().unwrap(),
        "-o",
        output_bam.to_str().unwrap(),
        "-f",
        "1.0",
        "--seed",
        "42",
        "--check-crc",
    ])
    .expect("failed to parse downsample args");
    let err = cmd
        .execute("fgumi downsample")
        .expect_err("--check-crc must reject a corrupted BGZF CRC32");
    let message = format!("{err:#}");
    assert!(message.to_uppercase().contains("CRC32"), "error should mention CRC32: {message}");
}

/// Read records from a BAM file.
fn read_bam_records(path: &PathBuf) -> Vec<noodles::sam::alignment::RecordBuf> {
    let mut reader = bam::io::reader::Builder.build_from_path(path).expect("Failed to open BAM");
    let header = reader.read_header().expect("Failed to read header");

    reader.record_bufs(&header).map(|r| r.expect("Failed to read record")).collect()
}

/// Count the number of unique MI values in a BAM file.
fn count_unique_mis(path: &PathBuf) -> usize {
    let records = read_bam_records(path);
    let mis: std::collections::HashSet<String> = records
        .iter()
        .filter_map(|r| {
            r.data().get(&MI_TAG).map(|v| {
                if let Value::String(s) = v {
                    s.to_string()
                } else {
                    panic!("MI tag is not a string")
                }
            })
        })
        .collect();
    mis.len()
}

/// Test basic downsampling functionality.
#[test]
fn test_downsample_basic() {
    let temp_dir = TempDir::new().expect("Failed to create temp dir");
    let input_bam = temp_dir.path().join("input.bam");
    let output_bam = temp_dir.path().join("output.bam");

    // Create input BAM with 10 families of 5 reads each
    let families: Vec<(&str, usize)> =
        (0..10).map(|i| (Box::leak(format!("{i}").into_boxed_str()) as &str, 5)).collect();
    create_grouped_bam(&input_bam, families);

    // Run downsample with fraction=0.5 and seed for reproducibility
    let cmd = Downsample::try_parse_from([
        "downsample",
        "-i",
        input_bam.to_str().unwrap(),
        "-o",
        output_bam.to_str().unwrap(),
        "-f",
        "0.5",
        "--seed",
        "42",
        "--compression-level",
        "1",
    ])
    .expect("failed to parse downsample args");
    cmd.execute("fgumi downsample").expect("Downsample command failed");
    assert!(output_bam.exists(), "Output BAM not created");

    // Verify that some families were kept (with seed 42 and fraction 0.5, should be around 5)
    let output_records = read_bam_records(&output_bam);
    assert!(!output_records.is_empty(), "Output should have some records");
    assert!(output_records.len() < 50, "Output should have fewer records than input (50)");

    let output_families = count_unique_mis(&output_bam);
    assert!(output_families > 0, "Should have kept some families");
    assert!(output_families < 10, "Should have fewer than 10 families");
}

/// Test downsampling with rejects output.
#[test]
fn test_downsample_with_rejects() {
    let temp_dir = TempDir::new().expect("Failed to create temp dir");
    let input_bam = temp_dir.path().join("input.bam");
    let output_bam = temp_dir.path().join("output.bam");
    let rejects_bam = temp_dir.path().join("rejects.bam");

    // Create input BAM with 10 families
    let families: Vec<(&str, usize)> =
        (0..10).map(|i| (Box::leak(format!("{i}").into_boxed_str()) as &str, 5)).collect();
    create_grouped_bam(&input_bam, families);

    let cmd = Downsample::try_parse_from([
        "downsample",
        "-i",
        input_bam.to_str().unwrap(),
        "-o",
        output_bam.to_str().unwrap(),
        "-f",
        "0.5",
        "--seed",
        "42",
        "--rejects",
        rejects_bam.to_str().unwrap(),
        "--compression-level",
        "1",
    ])
    .expect("failed to parse downsample args");
    cmd.execute("fgumi downsample").expect("Downsample command failed");
    assert!(output_bam.exists(), "Output BAM not created");
    assert!(rejects_bam.exists(), "Rejects BAM not created");

    // Kept + rejected should equal input
    let output_count = read_bam_records(&output_bam).len();
    let rejects_count = read_bam_records(&rejects_bam).len();

    assert_eq!(output_count + rejects_count, 50, "Total records should be preserved");
}

/// Test determinism with same seed.
#[test]
fn test_downsample_deterministic() {
    let temp_dir = TempDir::new().expect("Failed to create temp dir");
    let input_bam = temp_dir.path().join("input.bam");
    let output1_bam = temp_dir.path().join("output1.bam");
    let output2_bam = temp_dir.path().join("output2.bam");

    // Create input BAM
    let families: Vec<(&str, usize)> =
        (0..20).map(|i| (Box::leak(format!("{i}").into_boxed_str()) as &str, 3)).collect();
    create_grouped_bam(&input_bam, families);

    // Run twice with same seed
    for output in [&output1_bam, &output2_bam] {
        let cmd = Downsample::try_parse_from([
            "downsample",
            "-i",
            input_bam.to_str().unwrap(),
            "-o",
            output.to_str().unwrap(),
            "-f",
            "0.3",
            "--seed",
            "12345",
            "--compression-level",
            "1",
        ])
        .expect("failed to parse downsample args");
        cmd.execute("fgumi downsample").expect("Downsample command failed");
    }

    // Both outputs should be identical
    let records1 = read_bam_records(&output1_bam);
    let records2 = read_bam_records(&output2_bam);

    assert_eq!(records1.len(), records2.len(), "Same seed should produce same count");

    // Check that the same families were selected
    let mi_set1: std::collections::HashSet<_> = records1
        .iter()
        .filter_map(|r| {
            r.data().get(&MI_TAG).map(|v| {
                if let Value::String(s) = v {
                    s.to_string()
                } else {
                    panic!("MI tag is not a string")
                }
            })
        })
        .collect();

    let mi_set2: std::collections::HashSet<_> = records2
        .iter()
        .filter_map(|r| {
            r.data().get(&MI_TAG).map(|v| {
                if let Value::String(s) = v {
                    s.to_string()
                } else {
                    panic!("MI tag is not a string")
                }
            })
        })
        .collect();

    assert_eq!(mi_set1, mi_set2, "Same seed should select same families");
}

/// Test histogram output.
#[test]
fn test_downsample_histograms() {
    let temp_dir = TempDir::new().expect("Failed to create temp dir");
    let input_bam = temp_dir.path().join("input.bam");
    let output_bam = temp_dir.path().join("output.bam");
    let hist_kept = temp_dir.path().join("hist_kept.txt");
    let hist_rejected = temp_dir.path().join("hist_rejected.txt");

    // Create input BAM with varied family sizes
    let families = vec![("0", 1), ("1", 2), ("2", 3), ("3", 4), ("4", 5)];
    create_grouped_bam(&input_bam, families);

    let cmd = Downsample::try_parse_from([
        "downsample",
        "-i",
        input_bam.to_str().unwrap(),
        "-o",
        output_bam.to_str().unwrap(),
        "-f",
        "0.5",
        "--seed",
        "42",
        "--histogram-kept",
        hist_kept.to_str().unwrap(),
        "--histogram-rejected",
        hist_rejected.to_str().unwrap(),
        "--compression-level",
        "1",
    ])
    .expect("failed to parse downsample args");
    cmd.execute("fgumi downsample").expect("Downsample command failed");
    assert!(hist_kept.exists(), "Kept histogram not created");
    assert!(hist_rejected.exists(), "Rejected histogram not created");

    // Check histogram format
    let kept_contents = fs::read_to_string(&hist_kept).unwrap();
    assert!(kept_contents.contains("family_size\tcount"), "Histogram should have header");

    let rejected_contents = fs::read_to_string(&hist_rejected).unwrap();
    assert!(rejected_contents.contains("family_size\tcount"), "Histogram should have header");
}

/// Test with fraction=1.0 (keep all).
#[test]
fn test_downsample_keep_all() {
    let temp_dir = TempDir::new().expect("Failed to create temp dir");
    let input_bam = temp_dir.path().join("input.bam");
    let output_bam = temp_dir.path().join("output.bam");

    let families = vec![("0", 5), ("1", 3), ("2", 7)];
    create_grouped_bam(&input_bam, families);

    let cmd = Downsample::try_parse_from([
        "downsample",
        "-i",
        input_bam.to_str().unwrap(),
        "-o",
        output_bam.to_str().unwrap(),
        "-f",
        "1.0",
        "--compression-level",
        "1",
    ])
    .expect("failed to parse downsample args");
    cmd.execute("fgumi downsample").expect("Downsample command failed");

    // All records should be kept
    let output_records = read_bam_records(&output_bam);
    assert_eq!(output_records.len(), 15, "All 15 records should be kept with fraction=1.0");
}

/// Collect `(read_name, MI)` for every record in a BAM, in file order. Used to
/// assert source-record IDENTITY (which named reads survived), not just counts.
fn read_bam_name_mi_pairs(path: &PathBuf) -> Vec<(String, String)> {
    read_bam_records(path)
        .iter()
        .map(|r| {
            let name = r.name().map(ToString::to_string).unwrap_or_default();
            let mi = match r.data().get(&MI_TAG) {
                Some(Value::String(s)) => s.to_string(),
                Some(_) => panic!("MI tag is not a string"),
                None => panic!("record is missing its MI tag"),
            };
            (name, mi)
        })
        .collect()
}

/// The molecule base of an MI: everything before the last `/` (mirrors
/// `extract_mi_base`), so `100/A` and `100/B` share the base `100`.
fn mi_base(mi: &str) -> &str {
    match mi.rfind('/') {
        Some(i) if i > 0 => &mi[..i],
        _ => mi,
    }
}

/// End-to-end: by default (molecule grouping), the two duplex strands of a
/// molecule are sampled as ONE family — kept or dropped together — so no
/// surviving molecule base ever appears with only one of its strands. This holds
/// for ANY seed, which is what makes the assertion robust: it is the atomicity
/// guarantee, not a specific RNG outcome. The complementary `--per-strand`
/// behavior (strands sampled independently) is covered by
/// `test_downsample_per_strand_splits_duplex_strands` and the `FamilyIterator`
/// unit tests.
#[test]
fn test_downsample_default_keeps_duplex_strands_together() {
    let temp_dir = TempDir::new().expect("Failed to create temp dir");
    let input_bam = temp_dir.path().join("input.bam");
    let output_bam = temp_dir.path().join("output.bam");

    // 12 duplex molecules, each with an /A strand (3 reads) and a /B strand
    // (2 reads) laid out consecutively as group emits them, plus one simplex
    // family with no strand suffix. Enough molecules that a fraction of 0.5
    // keeps some and drops some, exercising the all-or-nothing path on both.
    // create_grouped_bam names reads `read_{idx}` sequentially in family order,
    // so molecule i owns read indices [i*5, i*5+5) — the exact set we assert
    // survives as a unit below.
    let mut families: Vec<(String, usize)> = Vec::new();
    for i in 0..12 {
        families.push((format!("{i}/A"), 3));
        families.push((format!("{i}/B"), 2));
    }
    families.push(("999".to_string(), 4)); // simplex family, no suffix
    let families_ref: Vec<(&str, usize)> =
        families.iter().map(|(mi, n)| (mi.as_str(), *n)).collect();
    create_grouped_bam(&input_bam, families_ref);

    // No flag: molecule grouping is the default.
    let cmd = Downsample::try_parse_from([
        "downsample",
        "-i",
        input_bam.to_str().unwrap(),
        "-o",
        output_bam.to_str().unwrap(),
        "-f",
        "0.5",
        "--seed",
        "7",
        "--compression-level",
        "1",
    ])
    .expect("failed to parse downsample args");
    cmd.execute("fgumi downsample").expect("Downsample command failed");

    // Group the surviving reads by molecule base: which named reads and which
    // strand suffixes survived under each base.
    let mut names_by_base: HashMap<String, HashSet<String>> = HashMap::new();
    let mut strands_by_base: HashMap<String, HashSet<String>> = HashMap::new();
    for (name, mi) in read_bam_name_mi_pairs(&output_bam) {
        let base = mi_base(&mi).to_string();
        let suffix = mi.strip_prefix(base.as_str()).unwrap_or("").to_string();
        names_by_base.entry(base.clone()).or_default().insert(name);
        strands_by_base.entry(base).or_default().insert(suffix);
    }

    // At least one duplex molecule survived and at least one was dropped, so the
    // test actually exercises both branches (not a vacuous pass on an empty or
    // full output).
    let duplex_survivors = (0..12).filter(|i| names_by_base.contains_key(&i.to_string())).count();
    assert!(duplex_survivors > 0, "expected some duplex molecules to survive at f=0.5");
    assert!(duplex_survivors < 12, "expected some duplex molecules to be dropped at f=0.5");

    // The core invariant, asserted by source-record IDENTITY (not just counts):
    // every surviving duplex molecule kept BOTH strands and EXACTLY its own five
    // reads (`read_{i*5}`..`read_{i*5+4}`: 3 from /A + 2 from /B). A molecule with
    // only one strand — or with a read that belongs to a different molecule —
    // would mean the strands were sampled independently or mis-grouped.
    for i in 0..12 {
        let base = i.to_string();
        if let Some(names) = names_by_base.get(&base) {
            let expected_names: HashSet<String> =
                (i * 5..i * 5 + 5).map(|idx| format!("read_{idx}")).collect();
            assert_eq!(
                names, &expected_names,
                "molecule {base} survived with reads {names:?}, expected exactly {expected_names:?}"
            );
            assert_eq!(
                strands_by_base.get(&base),
                Some(&HashSet::from(["/A".to_string(), "/B".to_string()])),
                "molecule {base} survived without both /A and /B strands present"
            );
        }
    }
}

/// By default at `f=1.0`, a simplex family (an MI with no `/A`/`/B` strand
/// suffix) is preserved in full: exactly its own records survive, each still
/// carrying its original MI. Deterministic complement to
/// `keeps_duplex_strands_together`, which samples a mixed input at `f=0.5` where
/// any single family's survival is not guaranteed. As per CLAUDE.md path
/// instructions, which require coverage of simplex-record preservation by
/// identity, not just its shape.
#[test]
fn test_downsample_default_preserves_simplex_family() {
    let temp_dir = TempDir::new().expect("Failed to create temp dir");
    let input_bam = temp_dir.path().join("input.bam");
    let output_bam = temp_dir.path().join("output.bam");

    // A single simplex family: MI `999`, four reads (read_0..read_3), no suffix.
    create_grouped_bam(&input_bam, vec![("999", 4)]);

    // No flag: molecule grouping is the default (a no-op for a simplex MI,
    // which has no `/` to strip).
    let cmd = Downsample::try_parse_from([
        "downsample",
        "-i",
        input_bam.to_str().unwrap(),
        "-o",
        output_bam.to_str().unwrap(),
        "-f",
        "1.0",
        "--compression-level",
        "1",
    ])
    .expect("failed to parse downsample args");
    cmd.execute("fgumi downsample").expect("Downsample command failed");

    // f=1.0 keeps everything: exactly the four input reads survive, each with MI
    // `999`. Assert identity (which named reads), not just the count/MI shape.
    let pairs = read_bam_name_mi_pairs(&output_bam);
    let got: HashSet<(String, String)> = pairs.iter().cloned().collect();
    let expected: HashSet<(String, String)> =
        (0..4).map(|idx| (format!("read_{idx}"), "999".to_string())).collect();
    assert_eq!(
        got, expected,
        "all four simplex reads (read_0..read_3) must be kept with MI 999 at f=1.0, got {pairs:?}"
    );
}

/// With `--per-strand`, the legacy behavior is restored: a molecule's two strands
/// are sampled INDEPENDENTLY, so a surviving molecule base may have only one of
/// its strands. This is the opt-out inverse of the default atomicity invariant.
/// Seed-pinned so the input is large enough that some molecule splits under
/// independent draws (the whole reason `--per-strand` collapses duplex families).
#[test]
fn test_downsample_per_strand_splits_duplex_strands() {
    let temp_dir = TempDir::new().expect("Failed to create temp dir");
    let input_bam = temp_dir.path().join("input.bam");
    let output_bam = temp_dir.path().join("output.bam");

    // 12 duplex molecules laid out as group emits them (/A then /B per molecule).
    let mut families: Vec<(String, usize)> = Vec::new();
    for i in 0..12 {
        families.push((format!("{i}/A"), 3));
        families.push((format!("{i}/B"), 2));
    }
    let families_ref: Vec<(&str, usize)> =
        families.iter().map(|(mi, n)| (mi.as_str(), *n)).collect();
    create_grouped_bam(&input_bam, families_ref);

    let cmd = Downsample::try_parse_from([
        "downsample",
        "-i",
        input_bam.to_str().unwrap(),
        "-o",
        output_bam.to_str().unwrap(),
        "-f",
        "0.5",
        "--seed",
        "7",
        "--per-strand",
        "--compression-level",
        "1",
    ])
    .expect("failed to parse downsample args");
    cmd.execute("fgumi downsample").expect("Downsample command failed");

    // Each raw MI (`i/A`, `i/B`) was its own family with its own draw. Collect the
    // surviving read NAMES per raw MI, and the surviving strand suffixes per base.
    // create_grouped_bam names reads `read_{idx}` sequentially in family order, so
    // molecule i occupies [i*5, i*5+5): `i/A` = read_{i*5..i*5+3}, `i/B` = the next 2.
    let mut names_by_raw_mi: HashMap<String, HashSet<String>> = HashMap::new();
    let mut strands_by_base: HashMap<String, HashSet<String>> = HashMap::new();
    for (name, mi) in read_bam_name_mi_pairs(&output_bam) {
        let base = mi_base(&mi).to_string();
        let suffix = mi.strip_prefix(base.as_str()).unwrap_or("").to_string();
        names_by_raw_mi.entry(mi).or_default().insert(name);
        strands_by_base.entry(base).or_default().insert(suffix);
    }

    // Per-strand keeps or drops each raw MI family ATOMICALLY (it is still a
    // family-atomic sampler — only the family key differs). So every RETAINED raw
    // MI must contain exactly its own source reads — all three of an `/A`, both of
    // a `/B` — never a partial family. This is the raw-MI identity contract, the
    // per-strand analogue of the by-molecule identity check.
    for (mi, names) in &names_by_raw_mi {
        let i: usize = mi_base(mi).parse().expect("base MI is an integer in this fixture");
        let expected: HashSet<String> = if mi.ends_with("/A") {
            (i * 5..i * 5 + 3).map(|idx| format!("read_{idx}")).collect()
        } else {
            (i * 5 + 3..i * 5 + 5).map(|idx| format!("read_{idx}")).collect()
        };
        assert_eq!(
            names, &expected,
            "raw MI {mi} survived with reads {names:?}, expected exactly its own family {expected:?}"
        );
    }

    // The defining property of per-strand sampling: at least one molecule survived
    // with only ONE of its two strands. (Under the default molecule grouping this
    // can never happen — see keeps_duplex_strands_together.)
    let split = strands_by_base.values().any(|strands| strands.len() == 1);
    assert!(
        split,
        "expected at least one molecule to survive with a single strand under --per-strand, \
         got {strands_by_base:?}"
    );
}

/// Test error handling for invalid fractions. Each case writes to a
/// distinct output path so an artifact left by one case can never make a
/// later assertion pass by accident.
#[rstest::rstest]
#[case::zero("0.0", "output_zero.bam")]
#[case::greater_than_one("1.5", "output_gt1.bam")]
fn test_downsample_invalid_fraction(#[case] fraction: &str, #[case] output_name: &str) {
    let temp_dir = TempDir::new().expect("Failed to create temp dir");
    let input_bam = temp_dir.path().join("input.bam");
    let output_bam = temp_dir.path().join(output_name);

    create_grouped_bam(&input_bam, vec![("0", 5)]);

    let cmd = Downsample::try_parse_from([
        "downsample",
        "-i",
        input_bam.to_str().unwrap(),
        "-o",
        output_bam.to_str().unwrap(),
        "-f",
        fraction,
        "--compression-level",
        "1",
    ])
    .expect("failed to parse downsample args");
    assert!(cmd.execute("fgumi downsample").is_err(), "Fraction={fraction} should fail");
}
