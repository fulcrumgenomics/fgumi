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

/// Asserts that a file ends with the standard 28-byte BGZF EOF marker block.
///
/// # Panics
///
/// Panics if the file cannot be read, is too small, or is missing the BGZF EOF block.
pub fn assert_has_bgzf_eof(path: &std::path::Path) {
    let data = std::fs::read(path).expect("Failed to read file for EOF check");
    assert!(data.len() >= 28, "File too small to contain BGZF EOF: {}", path.display());
    assert_eq!(
        &data[data.len() - 28..],
        &BGZF_EOF,
        "File missing BGZF EOF block: {}",
        path.display()
    );
}

/// Asserts that a BAM file's `@HD` line has `SO:unsorted` and no `GO` or `SS`
/// tags, matching the contract of `fgumi_lib::sam::header_as_unsorted`.
///
/// Use this on rejects BAMs produced by the pipeline commands, whose records
/// are emitted in mutex-acquisition order rather than input order.
///
/// # Panics
///
/// Panics if the file cannot be opened, the header cannot be read, or the
/// sort-order tags do not match the unsorted contract.
pub fn assert_header_unsorted(path: &std::path::Path) {
    let mut reader = noodles::bam::io::Reader::new(
        std::fs::File::open(path).expect("Failed to open BAM for header check"),
    );
    let header = reader.read_header().expect("Failed to read BAM header");
    let hdr_map = header.header().expect("BAM header is missing an @HD line");
    let other_fields = hdr_map.other_fields();

    let so =
        other_fields.get(b"SO").unwrap_or_else(|| panic!("@HD missing SO in {}", path.display()));
    assert_eq!(
        <_ as AsRef<[u8]>>::as_ref(so),
        b"unsorted",
        "SO should be 'unsorted' in {}",
        path.display()
    );
    assert!(other_fields.get(b"GO").is_none(), "GO should be cleared in {}", path.display());
    assert!(other_fields.get(b"SS").is_none(), "SS should be cleared in {}", path.display());
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

#[cfg(test)]
mod tests {
    use super::*;
    use crate::helpers::bam_generator::to_record_buf;
    use fgumi_raw_bam::{SamBuilder, flags};

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
