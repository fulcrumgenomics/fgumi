//! Custom assertion helpers for integration tests.
//!
//! These helpers provide reusable assertions for verifying BAM record properties
//! in integration tests.

#![allow(dead_code)]
#![allow(clippy::cast_precision_loss)]

use fgumi_lib::sam::SamTag;
use noodles::sam::alignment::record::data::field::Tag;
use noodles::sam::alignment::record_buf::RecordBuf;
use noodles::sam::alignment::record_buf::data::field::Value;

/// The standard 28-byte BGZF EOF marker block.
const BGZF_EOF: [u8; 28] = [
    0x1f, 0x8b, 0x08, 0x04, 0x00, 0x00, 0x00, 0x00, 0x00, 0xff, 0x06, 0x00, 0x42, 0x43, 0x02, 0x00,
    0x1b, 0x00, 0x03, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00,
];

/// Reads a required string tag off an emitted BAM record.
///
/// The `bam::Record` counterpart to [`assert_mi_tag`] and friends, for tests that
/// compare a whole emitted record against an expected one rather than checking a
/// single field. Panics rather than returning an `Option`: a missing tag is a
/// defect in the command under test, not a case the caller should handle.
///
/// # Panics
///
/// Panics if the tag is absent, unreadable, or not a string.
pub fn string_tag(record: &noodles::bam::Record, tag: SamTag) -> String {
    match record.data().get(&tag.to_noodles_tag()) {
        Some(Ok(noodles::sam::alignment::record::data::field::Value::String(value))) => {
            value.to_string()
        }
        other => panic!("record must carry a string {tag:?} tag, got {other:?}"),
    }
}

/// Reads a required integer tag off an emitted BAM record.
///
/// Widens every SAM integer width to `i64`, so a caller's expectation does not
/// have to predict whether a small count was encoded as `Int8` or `Int32`.
///
/// # Panics
///
/// Panics if the tag is absent, unreadable, or not an integer.
pub fn int_tag(record: &noodles::bam::Record, tag: SamTag) -> i64 {
    match record.data().get(&tag.to_noodles_tag()) {
        Some(Ok(value)) => {
            value.as_int().unwrap_or_else(|| panic!("{tag:?} must be an integer, got {value:?}"))
        }
        other => panic!("record must carry an integer {tag:?} tag, got {other:?}"),
    }
}

/// Asserts that an in-memory BGZF stream ends with the standard 28-byte EOF block.
///
/// The byte-slice counterpart to [`assert_has_bgzf_eof`], for output captured
/// from a pipe rather than written to a file. `label` names the run in the
/// failure, since callers check several spellings of the same output.
///
/// # Panics
///
/// Panics if `bytes` is too small or is missing the BGZF EOF block.
pub fn assert_bytes_have_bgzf_eof(bytes: &[u8], label: &str) {
    let eof_len = BGZF_EOF.len();
    assert!(
        bytes.len() >= eof_len,
        "{label}: too small to contain BGZF EOF ({} bytes)",
        bytes.len()
    );
    assert_eq!(
        &bytes[bytes.len() - eof_len..],
        &BGZF_EOF,
        "{label}: missing the BGZF EOF block — the stream reads as truncated \
         (`samtools quickcheck` fails) even though its records parse"
    );
}

/// Asserts that a file ends with the standard 28-byte BGZF EOF marker block.
///
/// # Panics
///
/// Panics if the file cannot be read, is too small, or is missing the BGZF EOF block.
pub fn assert_has_bgzf_eof(path: &std::path::Path) {
    let data = std::fs::read(path).expect("Failed to read file for EOF check");
    assert_bytes_have_bgzf_eof(&data, &path.display().to_string());
}

/// Asserts that a rejects BAM's `@HD` sort-order fields (`SO`, `GO`, `SS`)
/// match those of the input BAM.
///
/// Pipeline commands route rejects through the unified pipeline's first-class
/// secondary output, which emits rejects in batch-input order (a subset of an
/// SO-X stream is still SO-X). The rejects BAM therefore inherits the input
/// header verbatim — including whatever sort-order claim the input made (or
/// the absence of one).
///
/// # Panics
///
/// Panics if either file cannot be opened, the header cannot be read, or the
/// sort-order fields differ.
pub fn assert_rejects_header_matches_input(
    rejects_path: &std::path::Path,
    input_path: &std::path::Path,
) {
    #[derive(Debug, Eq, PartialEq)]
    struct SortFields {
        so: Option<Vec<u8>>,
        go: Option<Vec<u8>>,
        ss: Option<Vec<u8>>,
    }

    fn read_sort_fields(path: &std::path::Path) -> SortFields {
        use noodles::sam::header::record::value::map::header::tag::{
            GROUP_ORDER, SORT_ORDER, SUBSORT_ORDER,
        };

        let mut reader = noodles::bam::io::Reader::new(
            std::fs::File::open(path).expect("Failed to open BAM for header check"),
        );
        let header = reader.read_header().expect("Failed to read BAM header");
        // `@HD` is optional in noodles' Header — treat its absence as "no fields".
        let Some(hdr_map) = header.header() else {
            return SortFields { so: None, go: None, ss: None };
        };
        // Use typed `Other` tag constants rather than raw byte slices so this
        // keeps working if noodles ever promotes these to first-class fields
        // on `Map<Header>` (currently 0.82 keeps them in `other_fields`).
        let other = hdr_map.other_fields();
        SortFields {
            so: other.get(&SORT_ORDER).map(|v| v.to_vec()),
            go: other.get(&GROUP_ORDER).map(|v| v.to_vec()),
            ss: other.get(&SUBSORT_ORDER).map(|v| v.to_vec()),
        }
    }

    let actual = read_sort_fields(rejects_path);
    let expected = read_sort_fields(input_path);
    assert_eq!(
        actual,
        expected,
        "rejects @HD sort fields should match input ({} vs {})",
        rejects_path.display(),
        input_path.display()
    );
}

/// Asserts that a record has a specific MI (molecule ID) tag value.
///
/// # Panics
///
/// Panics if the MI tag is missing or has unexpected value.
pub fn assert_mi_tag(record: &RecordBuf, expected: &str) {
    let mi_tag = Tag::from(SamTag::MI);
    let mi_value = record.data().get(&mi_tag).expect("Record should have MI tag");

    match mi_value {
        Value::String(s) => {
            let s_bytes: &[u8] = s.as_ref();
            assert_eq!(
                s_bytes,
                expected.as_bytes(),
                "MI tag mismatch for record {:?}",
                record.name()
            );
        }
        _ => panic!("MI tag should be a string"),
    }
}

/// Asserts that a record has a specific RX (UMI) tag value.
///
/// # Panics
///
/// Panics if the RX tag is missing or has unexpected value.
pub fn assert_rx_tag(record: &RecordBuf, expected: &str) {
    let rx_tag = Tag::from(SamTag::RX);
    let rx_value = record.data().get(&rx_tag).expect("Record should have RX tag");

    match rx_value {
        Value::String(s) => {
            let s_bytes: &[u8] = s.as_ref();
            assert_eq!(
                s_bytes,
                expected.as_bytes(),
                "RX tag mismatch for record {:?}",
                record.name()
            );
        }
        _ => panic!("RX tag should be a string"),
    }
}

/// Asserts that consensus quality improved compared to input reads.
///
/// Checks that the consensus read has:
/// - Higher minimum quality score than any input read
/// - Mean quality >= mean of input reads
///
/// # Panics
///
/// Panics if consensus quality is not improved.
pub fn assert_consensus_quality_improved(input_reads: &[RecordBuf], consensus: &RecordBuf) {
    // Calculate mean quality of input reads
    let input_mean_quality: f64 = input_reads
        .iter()
        .flat_map(|r| r.quality_scores().as_ref().iter())
        .map(|&q| f64::from(q))
        .sum::<f64>()
        / input_reads.iter().map(|r| r.quality_scores().len()).sum::<usize>() as f64;

    // Calculate mean quality of consensus
    let consensus_mean_quality: f64 =
        consensus.quality_scores().as_ref().iter().map(|&q| f64::from(q)).sum::<f64>()
            / consensus.quality_scores().len() as f64;

    assert!(
        consensus_mean_quality >= input_mean_quality,
        "Consensus mean quality ({consensus_mean_quality}) should be >= input mean quality ({input_mean_quality})"
    );
}

/// Asserts that all records in a family have the same MI tag.
///
/// # Panics
///
/// Panics if records have different MI tags or any record is missing MI tag.
pub fn assert_same_molecule_id(records: &[RecordBuf]) {
    if records.is_empty() {
        return;
    }

    let mi_tag = Tag::from(SamTag::MI);
    let first_mi = records[0].data().get(&mi_tag).expect("First record should have MI tag");

    for record in records.iter().skip(1) {
        let mi = record.data().get(&mi_tag).expect("Record should have MI tag");
        assert_eq!(mi, first_mi, "All records in family should have same MI tag");
    }
}

/// Asserts that reads are properly paired (R1 and R2 flags set correctly).
///
/// # Panics
///
/// Panics if pairing is incorrect.
pub fn assert_proper_pairing(r1: &RecordBuf, r2: &RecordBuf) {
    assert!(r1.flags().is_segmented(), "R1 should have SEGMENTED flag");
    assert!(r2.flags().is_segmented(), "R2 should have SEGMENTED flag");

    assert!(r1.flags().is_first_segment(), "R1 should have FIRST_SEGMENT flag");
    assert!(!r2.flags().is_first_segment(), "R2 should not have FIRST_SEGMENT flag");

    assert!(r2.flags().is_last_segment(), "R2 should have LAST_SEGMENT flag");
    assert!(!r1.flags().is_last_segment(), "R1 should not have LAST_SEGMENT flag");
}

/// Asserts that two UMI families have different molecule IDs.
///
/// # Panics
///
/// Panics if the families have the same molecule ID.
pub fn assert_different_molecule_ids(family1: &[RecordBuf], family2: &[RecordBuf]) {
    if family1.is_empty() || family2.is_empty() {
        return;
    }

    let mi_tag = Tag::from(SamTag::MI);

    let mi1 = family1[0].data().get(&mi_tag).expect("Family 1 should have MI tag");
    let mi2 = family2[0].data().get(&mi_tag).expect("Family 2 should have MI tag");

    assert_ne!(mi1, mi2, "Different UMI families should have different molecule IDs");
}

/// Asserts that `bam` is actually sorted in `order`, using fgumi's own canonical
/// sort-order verifier (`fgumi sort --verify`).
///
/// A command's `@HD` `SO`/`GO`/`SS` tags are only a *promise* that the output is in
/// that order; this checks the promise against the bytes the command wrote. Pass
/// `key_types = Some("mi")` for template-coordinate BAMs that carry MI tags (e.g.
/// `group`/`dedup` output), mirroring `fgumi sort --key-types`.
///
/// # Panics
///
/// Panics if the `fgumi` binary cannot be spawned or reports any sort-order violation.
pub fn assert_bam_sorted(bam: &std::path::Path, order: &str, key_types: Option<&str>) {
    let bam_path = bam.to_str().expect("BAM path is valid UTF-8");
    let mut cmd = std::process::Command::new(env!("CARGO_BIN_EXE_fgumi"));
    cmd.args(["sort", "-i", bam_path, "--verify", "--order", order]);
    if let Some(kt) = key_types {
        cmd.args(["--key-types", kt]);
    }
    let status = cmd.status().expect("failed to spawn `fgumi sort --verify`");
    assert!(
        status.success(),
        "{bam_path} is not sorted by {order}: `fgumi sort --verify` reported violations"
    );
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::helpers::bam_generator::to_record_buf;
    use fgumi_raw_bam::{SamBuilder, flags};

    /// Locks down the noodles behavior that `assert_rejects_header_matches_input`
    /// relies on: SO/GO/SS are accessible via the typed `Other` tag constants on
    /// `Map<Header>::other_fields()`. If noodles ever moves these to first-class
    /// fields on `Map<Header>`, this test will fail and the helper above must be
    /// updated to use the new accessors instead.
    #[test]
    fn sort_fields_are_readable_from_other_fields() {
        use noodles::sam::header::record::value::map::header::tag::{
            GROUP_ORDER, SORT_ORDER, SUBSORT_ORDER,
        };

        let header_text =
            "@HD\tVN:1.6\tSO:queryname\tGO:none\tSS:queryname:natural\n@PG\tID:test\tPN:test\n";
        let header: noodles::sam::Header = header_text.parse().unwrap();
        let hdr_map = header.header().expect("@HD should be present");
        let other = hdr_map.other_fields();
        assert_eq!(other.get(&SORT_ORDER).map(|v| v.to_vec()), Some(b"queryname".to_vec()),);
        assert_eq!(other.get(&GROUP_ORDER).map(|v| v.to_vec()), Some(b"none".to_vec()));
        assert_eq!(
            other.get(&SUBSORT_ORDER).map(|v| v.to_vec()),
            Some(b"queryname:natural".to_vec()),
        );
    }

    #[test]
    fn test_assert_mi_tag() {
        let raw = {
            let mut b = SamBuilder::new();
            b.read_name(b"test")
                .sequence(b"ACGT")
                .qualities(&[30; 4])
                .add_string_tag(SamTag::MI, b"molecule_123");
            b.build()
        };
        let record = to_record_buf(&raw);
        assert_mi_tag(&record, "molecule_123");
    }

    #[test]
    #[should_panic(expected = "MI tag mismatch")]
    fn test_assert_mi_tag_mismatch() {
        let raw = {
            let mut b = SamBuilder::new();
            b.read_name(b"test")
                .sequence(b"ACGT")
                .qualities(&[30; 4])
                .add_string_tag(SamTag::MI, b"molecule_123");
            b.build()
        };
        let record = to_record_buf(&raw);
        assert_mi_tag(&record, "wrong_id");
    }

    #[test]
    fn test_assert_proper_pairing() {
        let raw_r1 = {
            let mut b = SamBuilder::new();
            b.read_name(b"pair1")
                .sequence(b"ACGT")
                .qualities(&[30; 4])
                .flags(flags::PAIRED | flags::FIRST_SEGMENT);
            b.build()
        };
        let raw_r2 = {
            let mut b = SamBuilder::new();
            b.read_name(b"pair1")
                .sequence(b"TGCA")
                .qualities(&[30; 4])
                .flags(flags::PAIRED | flags::LAST_SEGMENT);
            b.build()
        };
        let r1 = to_record_buf(&raw_r1);
        let r2 = to_record_buf(&raw_r2);
        assert_proper_pairing(&r1, &r2);
    }
}
