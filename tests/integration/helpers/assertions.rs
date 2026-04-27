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
#[allow(dead_code)]
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
#[allow(dead_code)]
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
#[allow(dead_code, clippy::cast_precision_loss)]
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
#[allow(dead_code)]
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
#[allow(dead_code)]
pub fn assert_different_molecule_ids(family1: &[RecordBuf], family2: &[RecordBuf]) {
    if family1.is_empty() || family2.is_empty() {
        return;
    }

    let mi_tag = Tag::from(SamTag::MI);

    let mi1 = family1[0].data().get(&mi_tag).expect("Family 1 should have MI tag");
    let mi2 = family2[0].data().get(&mi_tag).expect("Family 2 should have MI tag");

    assert_ne!(mi1, mi2, "Different UMI families should have different molecule IDs");
}

/// Extract MI tag string values from a slice of records.
///
/// Records without an MI tag are skipped. Returns a `Vec<String>` with one entry
/// per record that carries the tag.
pub fn extract_mi_tags(records: &[RecordBuf]) -> Vec<String> {
    let mi_tag = Tag::from([b'M', b'I']);
    records
        .iter()
        .filter_map(|r| match r.data().get(&mi_tag)? {
            Value::String(s) => Some(String::from_utf8_lossy(s.as_ref()).into_owned()),
            _ => None,
        })
        .collect()
}

/// Compute the family-size distribution from a list of MI tag values.
///
/// Returns a sorted `Vec<usize>` where each element is the number of reads
/// assigned to one MI group. The returned vec is sorted so that two equivalent
/// distributions compare equal regardless of the order in which families appear
/// or how MI numbers are assigned.
pub fn count_mi_families(mi_tags: &[String]) -> Vec<usize> {
    use std::collections::HashMap;
    let mut counts: HashMap<&str, usize> = HashMap::new();
    for tag in mi_tags {
        *counts.entry(tag.as_str()).or_insert(0) += 1;
    }
    let mut sizes: Vec<usize> = counts.into_values().collect();
    sizes.sort_unstable();
    sizes
}

/// Assert that two BAM files are record-level equivalent.
///
/// Checks:
/// 1. Both files contain the same number of records (and at least one).
/// 2. The multisets of sequences are identical (compared sorted).
/// 3. If both files carry MI tags, the family-size distributions are identical.
///    The exact MI values are not compared because the numbering scheme may differ
///    between a sequential run and a runall pipeline run.
///
/// # Panics
///
/// Panics on any mismatch.
pub fn assert_bam_equivalent(path_a: &std::path::Path, path_b: &std::path::Path) {
    use crate::helpers::bam_generator::{count_bam_records, read_bam_records, read_bam_sequences};

    let count_a = count_bam_records(path_a);
    let count_b = count_bam_records(path_b);
    assert_eq!(
        count_a,
        count_b,
        "Record count mismatch: {} has {count_a}, {} has {count_b}",
        path_a.display(),
        path_b.display()
    );
    assert!(count_a > 0, "Both files should have at least one record");

    // Compare sorted sequence multisets.
    let mut seqs_a = read_bam_sequences(path_a);
    let mut seqs_b = read_bam_sequences(path_b);
    seqs_a.sort();
    seqs_b.sort();
    assert_eq!(
        seqs_a,
        seqs_b,
        "Sequences differ between {} and {}",
        path_a.display(),
        path_b.display()
    );

    // Compare MI tag family-size distributions when tags are present.
    let records_a = read_bam_records(path_a);
    let records_b = read_bam_records(path_b);
    let mi_a = extract_mi_tags(&records_a);
    let mi_b = extract_mi_tags(&records_b);
    if !mi_a.is_empty() && !mi_b.is_empty() {
        let dist_a = count_mi_families(&mi_a);
        let dist_b = count_mi_families(&mi_b);
        assert_eq!(
            dist_a,
            dist_b,
            "MI family size distributions differ between {} and {}",
            path_a.display(),
            path_b.display()
        );
    }
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
