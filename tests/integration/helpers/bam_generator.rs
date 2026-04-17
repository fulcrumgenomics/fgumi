//! Utilities for generating test BAM data programmatically.

use fgumi_raw_bam::{RawRecord, SamBuilder, flags};
use noodles::sam::Header;
use noodles::sam::alignment::record_buf::RecordBuf;

/// Convert a `RawRecord` to a noodles `RecordBuf` using the default (empty) header.
///
/// Used to bridge the raw-byte builder API with the noodles BAM writer for test
/// file creation.
pub fn to_record_buf(raw: &RawRecord) -> RecordBuf {
    fgumi_raw_bam::raw_record_to_record_buf(raw, &noodles::sam::Header::default())
        .expect("raw_record_to_record_buf should succeed in test")
}

/// Creates a UMI family with specified parameters.
///
/// All reads in the family are mapped to reference sequence 0 at position 100
/// with a simple match CIGAR. This ensures they pass the group command's
/// unmapped filter.
///
/// # Arguments
///
/// * `umi` - The UMI sequence to assign
/// * `depth` - Number of reads in the family
/// * `base_name` - Base name for reads (will be suffixed with index)
/// * `sequence` - The read sequence (all reads will have this sequence)
/// * `quality` - Base quality score for all bases
///
/// # Returns
///
/// Vector of `RawRecord` representing the UMI family
pub fn create_umi_family(
    umi: &str,
    depth: usize,
    base_name: &str,
    sequence: &str,
    quality: u8,
) -> Vec<RawRecord> {
    let seq = sequence.as_bytes();
    let cigar_op = u32::try_from(seq.len()).expect("seq.len() fits u32") << 4; // NM
    (0..depth)
        .map(|i| {
            let name = format!("{base_name}_{i}");
            let mut b = SamBuilder::new();
            b.read_name(name.as_bytes())
                .ref_id(0)
                .pos(99)
                .mapq(60)
                .flags(0)
                .cigar_ops(&[cigar_op])
                .sequence(seq)
                .qualities(&vec![quality; seq.len()])
                .add_string_tag(b"RX", umi.as_bytes());
            b.build()
        })
        .collect()
}

/// Creates a paired-end UMI family.
///
/// All reads are mapped to reference sequence 0. R1 is mapped at position 100,
/// R2 at position 200 (to simulate insert size).
///
/// # Arguments
///
/// * `umi` - The UMI sequence (for paired UMIs, use "AAAA-CCCC" format)
/// * `depth` - Number of read pairs in the family
/// * `base_name` - Base name for reads
/// * `r1_sequence` - R1 sequence
/// * `r2_sequence` - R2 sequence
/// * `quality` - Base quality score
///
/// # Returns
///
/// Vector of `RawRecord` with R1 and R2 reads properly flagged
pub fn create_paired_umi_family(
    umi: &str,
    depth: usize,
    base_name: &str,
    r1_sequence: &str,
    r2_sequence: &str,
    quality: u8,
) -> Vec<RawRecord> {
    let r1_seq = r1_sequence.as_bytes();
    let r2_seq = r2_sequence.as_bytes();
    let r1_cigar = u32::try_from(r1_seq.len()).expect("r1_seq.len() fits u32") << 4;
    let r2_cigar = u32::try_from(r2_seq.len()).expect("r2_seq.len() fits u32") << 4;

    // R1 at pos 99 (0-based); R2 at pos 199. Template spans from R1 start through end
    // of R2 — 100bp gap + R2 length — and R2 carries the negated length.
    let template_len = i32::try_from(100 + r2_seq.len()).expect("template length fits i32");

    let mut records = Vec::new();

    for i in 0..depth {
        let read_name = format!("{base_name}_{i}");

        // R1: paired + first segment
        let mut b1 = SamBuilder::new();
        b1.read_name(read_name.as_bytes())
            .ref_id(0)
            .pos(99)
            .mapq(60)
            .flags(flags::PAIRED | flags::FIRST_SEGMENT)
            .mate_ref_id(0)
            .mate_pos(199)
            .template_length(template_len)
            .cigar_ops(&[r1_cigar])
            .sequence(r1_seq)
            .qualities(&vec![quality; r1_seq.len()])
            .add_string_tag(b"RX", umi.as_bytes());
        records.push(b1.build());

        // R2: paired + last segment
        let mut b2 = SamBuilder::new();
        b2.read_name(read_name.as_bytes())
            .ref_id(0)
            .pos(199)
            .mapq(60)
            .flags(flags::PAIRED | flags::LAST_SEGMENT)
            .mate_ref_id(0)
            .mate_pos(99)
            .template_length(-template_len)
            .cigar_ops(&[r2_cigar])
            .sequence(r2_seq)
            .qualities(&vec![quality; r2_seq.len()])
            .add_string_tag(b"RX", umi.as_bytes());
        records.push(b2.build());
    }

    records
}

/// Creates a consensus read with specified metrics.
///
/// # Arguments
///
/// * `name` - Read name
/// * `sequence` - Consensus sequence
/// * `depth_max` - Maximum depth (cD tag)
/// * `depth_min` - Minimum depth (cM tag)
/// * `error_rate` - Error rate (cE tag)
/// * `mean_quality` - Mean base quality score
///
/// # Returns
///
/// `RawRecord` representing a consensus read
pub fn create_consensus_read(
    name: &str,
    sequence: &str,
    depth_max: i32,
    depth_min: i32,
    error_rate: f32,
    mean_quality: u8,
) -> RawRecord {
    let seq = sequence.as_bytes();
    let mut b = SamBuilder::new();
    b.read_name(name.as_bytes())
        .sequence(seq)
        .qualities(&vec![mean_quality; seq.len()])
        .add_int_tag(b"cD", depth_max)
        .add_int_tag(b"cM", depth_min)
        .add_float_tag(b"cE", error_rate);
    b.build()
}

/// Builds a SAM header with the given header tags and one reference sequence.
///
/// This is a shared helper used by the public header constructors to avoid
/// duplicating the noodles header building boilerplate.
fn build_header_with_tags(ref_name: &str, ref_len: usize, tags: &[([u8; 2], &str)]) -> Header {
    use bstr::BString;
    use noodles::sam::header::record::value::map::Map as HeaderRecordMap;
    use noodles::sam::header::record::value::map::header::tag::Tag as HeaderTag;
    use noodles::sam::header::record::value::{
        Map, map::Header as HeaderRecord, map::ReferenceSequence,
    };
    use std::num::NonZeroUsize;

    let mut builder = HeaderRecordMap::<HeaderRecord>::builder();
    for &(tag_bytes, value) in tags {
        let HeaderTag::Other(tag) = HeaderTag::from(tag_bytes) else { unreachable!() };
        builder = builder.insert(tag, value);
    }
    let header_map = builder.build().expect("valid header map");

    let reference_sequence = Map::<ReferenceSequence>::new(
        NonZeroUsize::new(ref_len).expect("reference length must be non-zero"),
    );

    Header::builder()
        .set_header(header_map)
        .add_reference_sequence(BString::from(ref_name), reference_sequence)
        .build()
}

/// Creates a minimal SAM header with one reference sequence.
///
/// The header is configured with template-coordinate sort order tags
/// (SO:unsorted, GO:query, SS:template-coordinate) to be compatible with
/// the group command.
///
/// # Arguments
///
/// * `ref_name` - Reference sequence name (e.g., "chr1")
/// * `ref_len` - Reference sequence length
///
/// # Returns
///
/// Configured SAM Header with template-coordinate sort order
pub fn create_minimal_header(ref_name: &str, ref_len: usize) -> Header {
    build_header_with_tags(
        ref_name,
        ref_len,
        &[(*b"SO", "unsorted"), (*b"GO", "query"), (*b"SS", "template-coordinate")],
    )
}

/// Creates a coordinate-sorted SAM header with one reference sequence.
///
/// # Arguments
///
/// * `ref_name` - Reference sequence name (e.g., "chr1")
/// * `ref_len` - Reference sequence length
///
/// # Returns
///
/// Configured SAM Header with SO:coordinate sort order
pub fn create_coordinate_sorted_header(ref_name: &str, ref_len: usize) -> Header {
    build_header_with_tags(ref_name, ref_len, &[(*b"SO", "coordinate")])
}

/// Creates a UMI family with intentional sequencing errors.
///
/// # Arguments
///
/// * `base_umi` - The "true" UMI sequence
/// * `error_umi` - The error variant of the UMI
/// * `base_depth` - Number of reads with `base_umi`
/// * `error_depth` - Number of reads with `error_umi`
/// * `sequence` - Read sequence
///
/// # Returns
///
/// Combined vector of reads with both UMI variants
pub fn create_umi_family_with_errors(
    base_umi: &str,
    error_umi: &str,
    base_depth: usize,
    error_depth: usize,
    sequence: &str,
) -> Vec<RawRecord> {
    let mut records = Vec::new();

    // Add base UMI reads
    records.extend(create_umi_family(
        base_umi,
        base_depth,
        &format!("base_{base_umi}"),
        sequence,
        30,
    ));

    // Add error UMI reads
    records.extend(create_umi_family(
        error_umi,
        error_depth,
        &format!("error_{error_umi}"),
        sequence,
        30,
    ));

    records
}

/// Create a reference FASTA + FAI + sequence dictionary that matches the test header
/// (chr1, 10000bp).
///
/// Writes `ref.fa`, `ref.fa.fai`, and `ref.dict` into the given directory.
/// Returns the path to `ref.fa`.
pub fn create_test_reference(dir: &std::path::Path) -> std::path::PathBuf {
    use std::io::Write;

    let ref_path = dir.join("ref.fa");
    let fai_path = dir.join("ref.fa.fai");
    let dict_path = dir.join("ref.dict");

    let ref_seq = "ACGTACGT".repeat(1250); // 10000 bases
    let mut fasta = std::fs::File::create(&ref_path).unwrap();
    writeln!(fasta, ">chr1").unwrap();
    writeln!(fasta, "{ref_seq}").unwrap();
    fasta.flush().unwrap();

    let fai_content = "chr1\t10000\t6\t10000\t10001\n";
    std::fs::write(&fai_path, fai_content).unwrap();

    let mut dict = std::fs::File::create(&dict_path).unwrap();
    writeln!(dict, "@HD\tVN:1.6\tSO:unsorted").unwrap();
    writeln!(dict, "@SQ\tSN:chr1\tLN:10000").unwrap();
    dict.flush().unwrap();

    ref_path
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_create_umi_family() {
        let family = create_umi_family("ACGTACGT", 5, "test", "AAAA", 30);
        assert_eq!(family.len(), 5);

        for (i, record) in family.iter().enumerate() {
            let buf = to_record_buf(record);
            assert_eq!(
                buf.name().map(std::convert::AsRef::as_ref),
                Some(format!("test_{i}").as_bytes())
            );
            assert_eq!(buf.sequence().as_ref(), b"AAAA");
        }
    }

    #[test]
    fn test_create_paired_umi_family() {
        let family = create_paired_umi_family("ACGT-TGCA", 3, "pair", "AAAA", "TTTT", 30);

        // Should have 6 records (3 pairs)
        assert_eq!(family.len(), 6);

        // Check R1/R2 flags
        let r1_flags = family[0].flags();
        let r2_flags = family[1].flags();
        assert_ne!(r1_flags & flags::FIRST_SEGMENT, 0, "R1 should have FIRST_SEGMENT flag");
        assert_eq!(r2_flags & flags::FIRST_SEGMENT, 0, "R2 should not have FIRST_SEGMENT flag");
        assert_ne!(r2_flags & flags::LAST_SEGMENT, 0, "R2 should have LAST_SEGMENT flag");
    }

    #[test]
    fn test_create_consensus_read() {
        let consensus = create_consensus_read("cons1", "ACGT", 10, 5, 0.01, 35);
        let buf = to_record_buf(&consensus);

        assert_eq!(buf.name().map(std::convert::AsRef::as_ref), Some(b"cons1".as_ref()));
        assert_eq!(buf.sequence().as_ref(), b"ACGT");
        assert_eq!(buf.quality_scores().as_ref(), &[35, 35, 35, 35]);
    }

    #[test]
    fn test_create_minimal_header() {
        let header = create_minimal_header("chr1", 1000);
        assert_eq!(header.reference_sequences().len(), 1);
    }
}
