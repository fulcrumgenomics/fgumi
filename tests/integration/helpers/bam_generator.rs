//! Utilities for generating test BAM data programmatically.

use fgumi_lib::sam::builder::{ConsensusTagsBuilder, RecordBuilder};
use noodles::sam::Header;
use noodles::sam::alignment::record_buf::RecordBuf;

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
/// Vector of `RecordBuf` representing the UMI family
pub fn create_umi_family(
    umi: &str,
    depth: usize,
    base_name: &str,
    sequence: &str,
    quality: u8,
) -> Vec<RecordBuf> {
    (0..depth)
        .map(|i| {
            RecordBuilder::new()
                .name(&format!("{base_name}_{i}"))
                .sequence(sequence)
                .qualities(&vec![quality; sequence.len()])
                .reference_sequence_id(0)
                .alignment_start(100)
                .mapping_quality(60)
                .tag("RX", umi)
                .build()
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
/// Vector of `RecordBuf` with R1 and R2 reads properly flagged
pub fn create_paired_umi_family(
    umi: &str,
    depth: usize,
    base_name: &str,
    r1_sequence: &str,
    r2_sequence: &str,
    quality: u8,
) -> Vec<RecordBuf> {
    let mut records = Vec::new();

    for i in 0..depth {
        let read_name = format!("{base_name}_{i}");

        // R1
        let r1 = RecordBuilder::new()
            .name(&read_name)
            .sequence(r1_sequence)
            .qualities(&vec![quality; r1_sequence.len()])
            .paired(true)
            .first_segment(true)
            .reference_sequence_id(0)
            .alignment_start(100)
            .mapping_quality(60)
            .tag("RX", umi)
            .build();
        records.push(r1);

        // R2
        let r2 = RecordBuilder::new()
            .name(&read_name)
            .sequence(r2_sequence)
            .qualities(&vec![quality; r2_sequence.len()])
            .paired(true)
            .first_segment(false)
            .reference_sequence_id(0)
            .alignment_start(200)
            .mapping_quality(60)
            .tag("RX", umi)
            .build();
        records.push(r2);
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
/// `RecordBuf` representing a consensus read
pub fn create_consensus_read(
    name: &str,
    sequence: &str,
    depth_max: i32,
    depth_min: i32,
    error_rate: f32,
    mean_quality: u8,
) -> RecordBuf {
    RecordBuilder::new()
        .name(name)
        .sequence(sequence)
        .qualities(&vec![mean_quality; sequence.len()])
        .consensus_tags(
            ConsensusTagsBuilder::new()
                .depth_max(depth_max)
                .depth_min(depth_min)
                .error_rate(error_rate),
        )
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
    use bstr::BString;
    use noodles::sam::header::record::value::map::Map as HeaderRecordMap;
    use noodles::sam::header::record::value::map::header::tag::Tag as HeaderTag;
    use noodles::sam::header::record::value::{
        Map, map::Header as HeaderRecord, map::ReferenceSequence,
    };
    use std::num::NonZeroUsize;

    // Create header with template-coordinate sort order tags
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

    Header::builder()
        .set_header(header_map)
        .add_reference_sequence(BString::from(ref_name), reference_sequence)
        .build()
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
) -> Vec<RecordBuf> {
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

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_create_umi_family() {
        let family = create_umi_family("ACGTACGT", 5, "test", "AAAA", 30);
        assert_eq!(family.len(), 5);

        for (i, record) in family.iter().enumerate() {
            assert_eq!(
                record.name().map(std::convert::AsRef::as_ref),
                Some(format!("test_{i}").as_bytes())
            );
            assert_eq!(record.sequence().as_ref(), b"AAAA");
        }
    }

    #[test]
    fn test_create_paired_umi_family() {
        let family = create_paired_umi_family("ACGT-TGCA", 3, "pair", "AAAA", "TTTT", 30);

        // Should have 6 records (3 pairs)
        assert_eq!(family.len(), 6);

        // Check R1/R2 flags
        assert!(family[0].flags().is_first_segment());
        assert!(!family[1].flags().is_first_segment());
        assert!(family[1].flags().is_last_segment());
    }

    #[test]
    fn test_create_consensus_read() {
        let consensus = create_consensus_read("cons1", "ACGT", 10, 5, 0.01, 35);

        assert_eq!(consensus.name().map(std::convert::AsRef::as_ref), Some(b"cons1".as_ref()));
        assert_eq!(consensus.sequence().as_ref(), b"ACGT");
        assert_eq!(consensus.quality_scores().as_ref(), &[35, 35, 35, 35]);
    }

    #[test]
    fn test_create_minimal_header() {
        let header = create_minimal_header("chr1", 1000);
        assert_eq!(header.reference_sequences().len(), 1);
    }
}
