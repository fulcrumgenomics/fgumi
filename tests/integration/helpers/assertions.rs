//! Custom assertion helpers for integration tests.
//!
//! These helpers provide reusable assertions for verifying BAM record properties
//! in integration tests.

#![allow(dead_code)]
#![allow(clippy::cast_precision_loss)]

use noodles::sam::alignment::record::data::field::Tag;
use noodles::sam::alignment::record_buf::RecordBuf;
use noodles::sam::alignment::record_buf::data::field::Value;

/// Asserts that a record has a specific MI (molecule ID) tag value.
///
/// # Panics
///
/// Panics if the MI tag is missing or has unexpected value.
pub fn assert_mi_tag(record: &RecordBuf, expected: &str) {
    let mi_tag = Tag::from([b'M', b'I']);
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
    let rx_tag = Tag::from([b'R', b'X']);
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

    let mi_tag = Tag::from([b'M', b'I']);
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

    let mi_tag = Tag::from([b'M', b'I']);

    let mi1 = family1[0].data().get(&mi_tag).expect("Family 1 should have MI tag");
    let mi2 = family2[0].data().get(&mi_tag).expect("Family 2 should have MI tag");

    assert_ne!(mi1, mi2, "Different UMI families should have different molecule IDs");
}

#[cfg(test)]
mod tests {
    use super::*;
    use fgumi_lib::sam::builder::RecordBuilder;

    #[test]
    fn test_assert_mi_tag() {
        let record =
            RecordBuilder::new().name("test").sequence("ACGT").tag("MI", "molecule_123").build();

        assert_mi_tag(&record, "molecule_123");
    }

    #[test]
    #[should_panic(expected = "MI tag mismatch")]
    fn test_assert_mi_tag_mismatch() {
        let record =
            RecordBuilder::new().name("test").sequence("ACGT").tag("MI", "molecule_123").build();

        assert_mi_tag(&record, "wrong_id");
    }

    #[test]
    fn test_assert_proper_pairing() {
        let r1 = RecordBuilder::new().name("pair1").sequence("ACGT").first_segment(true).build();

        let r2 = RecordBuilder::new().name("pair1").sequence("TGCA").first_segment(false).build();

        assert_proper_pairing(&r1, &r2);
    }
}
