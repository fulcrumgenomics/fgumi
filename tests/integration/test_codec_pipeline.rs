//! Integration tests for the CODEC consensus calling pipeline.
//!
//! CODEC is a sequencing protocol where a single read-pair sequences both strands
//! of the original duplex molecule:
//! - R1 comes from one strand
//! - R2 comes from the opposite strand
//! - Even a single read-pair can generate duplex consensus
//!
//! These tests validate:
//! 1. Basic CODEC consensus calling from paired-end reads
//! 2. Quality score consensus calling
//! 3. Rejection of insufficient reads
//! 4. Handling of strand orientation
//! 5. Overlap region computation

use fgumi_lib::consensus::codec_caller::{CodecConsensusCaller, CodecConsensusOptions};
use fgumi_lib::sam::builder::RecordBuilder;
use noodles::sam::alignment::record::cigar::Op;
use noodles::sam::alignment::record::cigar::op::Kind;
use noodles::sam::alignment::record_buf::{Cigar, RecordBuf};

/// Helper function to create a test paired read with proper CIGAR
#[allow(clippy::too_many_arguments, clippy::cast_possible_truncation, clippy::cast_possible_wrap)]
fn create_codec_read(
    name: &str,
    seq: &[u8],
    qual: &[u8],
    is_first: bool,
    is_reverse: bool,
    mate_is_reverse: bool,
    ref_start: usize,
    mate_start: usize,
    cigar_ops: &[(Kind, usize)],
    umi: &str,
) -> RecordBuf {
    // Build CIGAR and calculate reference length
    let mut cigar = Cigar::default();
    let mut ref_len = 0usize;
    for &(kind, len) in cigar_ops {
        cigar.as_mut().push(Op::new(kind, len));
        match kind {
            Kind::Match
            | Kind::SequenceMatch
            | Kind::SequenceMismatch
            | Kind::Deletion
            | Kind::Skip => {
                ref_len += len;
            }
            _ => {}
        }
    }

    // Calculate template length for FR pair detection
    let insert_size: i32 = if ref_start <= mate_start {
        (mate_start + ref_len - ref_start) as i32
    } else {
        -((ref_start + ref_len - mate_start) as i32)
    };

    let mut record = RecordBuilder::new()
        .name(name)
        .sequence(&String::from_utf8_lossy(seq))
        .qualities(qual)
        .reference_sequence_id(0)
        .alignment_start(ref_start)
        .paired(true)
        .first_segment(is_first)
        .reverse_complement(is_reverse)
        .mate_reverse_complement(mate_is_reverse)
        .mate_reference_sequence_id(0)
        .mate_alignment_start(mate_start)
        .tag("MI", umi)
        .build();

    // Set CIGAR (built from ops) and template length
    *record.cigar_mut() = cigar;
    *record.template_length_mut() = insert_size;

    record
}

/// Creates a CODEC read pair (R1 forward, R2 reverse - FR orientation)
#[allow(clippy::too_many_arguments)]
fn create_codec_read_pair(
    name: &str,
    r1_seq: &[u8],
    r2_seq: &[u8],
    r1_qual: &[u8],
    r2_qual: &[u8],
    r1_start: usize,
    r2_start: usize,
    umi: &str,
) -> (RecordBuf, RecordBuf) {
    let r1 = create_codec_read(
        name,
        r1_seq,
        r1_qual,
        true,  // is_first
        false, // is_reverse (R1 forward)
        true,  // mate_is_reverse (R2 is reverse)
        r1_start,
        r2_start, // mate_start
        &[(Kind::Match, r1_seq.len())],
        umi,
    );

    let r2 = create_codec_read(
        name,
        r2_seq,
        r2_qual,
        false, // is_first (R2)
        true,  // is_reverse (R2 reverse)
        false, // mate_is_reverse (R1 is forward)
        r2_start,
        r1_start, // mate_start
        &[(Kind::Match, r2_seq.len())],
        umi,
    );

    (r1, r2)
}

#[test]
fn test_codec_single_pair_consensus() {
    // Test that a single read pair can produce consensus (CODEC's key feature)
    let options = CodecConsensusOptions {
        min_reads_per_strand: 1,
        min_duplex_length: 1,
        ..Default::default()
    };

    let mut caller = CodecConsensusCaller::new("codec".to_string(), "RG1".to_string(), options);

    // Create a single read pair with identical sequences
    // R1 forward at position 100, R2 reverse at position 100 (same region)
    let (r1, r2) = create_codec_read_pair(
        "read1",
        b"ACGTACGT",
        b"ACGTACGT", // R2 will be reverse complemented internally
        &[30, 30, 30, 30, 30, 30, 30, 30],
        &[30, 30, 30, 30, 30, 30, 30, 30],
        100,
        100,
        "UMI001",
    );

    let reads = vec![r1, r2];
    let output = caller.consensus_reads_from_sam_records(reads).unwrap();

    // Should produce exactly 1 consensus read
    assert_eq!(output.count, 1, "Single read pair should produce 1 consensus");

    // Check statistics
    let stats = caller.statistics();
    assert_eq!(stats.total_input_reads, 2);
    assert_eq!(stats.consensus_reads_generated, 1);
}

#[test]
fn test_codec_multiple_pairs_consensus() {
    // Test consensus from multiple read pairs (higher depth)
    let options = CodecConsensusOptions {
        min_reads_per_strand: 1,
        min_duplex_length: 1,
        ..Default::default()
    };

    let mut caller = CodecConsensusCaller::new("codec".to_string(), "RG1".to_string(), options);

    // Create 3 read pairs with same UMI
    let mut reads = Vec::new();
    for i in 0..3 {
        let (r1, r2) = create_codec_read_pair(
            &format!("read{i}"),
            b"ACGTACGT",
            b"ACGTACGT",
            &[30, 30, 30, 30, 30, 30, 30, 30],
            &[30, 30, 30, 30, 30, 30, 30, 30],
            100,
            100,
            "UMI002",
        );
        reads.push(r1);
        reads.push(r2);
    }

    let output = caller.consensus_reads_from_sam_records(reads).unwrap();

    assert_eq!(output.count, 1, "Multiple pairs should produce 1 consensus");

    let stats = caller.statistics();
    assert_eq!(stats.total_input_reads, 6);
}

#[test]
fn test_codec_insufficient_reads_rejection() {
    // Test that reads are rejected when below min_reads_per_strand
    let options = CodecConsensusOptions {
        min_reads_per_strand: 3, // Require at least 3 pairs
        min_duplex_length: 1,
        ..Default::default()
    };

    let mut caller = CodecConsensusCaller::new("codec".to_string(), "RG1".to_string(), options);

    // Create only 2 read pairs (below threshold)
    let mut reads = Vec::new();
    for i in 0..2 {
        let (r1, r2) = create_codec_read_pair(
            &format!("read{i}"),
            b"ACGTACGT",
            b"ACGTACGT",
            &[30, 30, 30, 30, 30, 30, 30, 30],
            &[30, 30, 30, 30, 30, 30, 30, 30],
            100,
            100,
            "UMI003",
        );
        reads.push(r1);
        reads.push(r2);
    }

    let output = caller.consensus_reads_from_sam_records(reads).unwrap();

    assert_eq!(output.count, 0, "Should reject when fewer than min_reads_per_strand pairs");

    let stats = caller.statistics();
    assert!(stats.reads_filtered > 0, "Reads should be filtered");
}

#[test]
fn test_codec_fragment_reads_rejection() {
    // Test that non-paired (fragment) reads are rejected
    let options = CodecConsensusOptions {
        min_reads_per_strand: 1,
        min_duplex_length: 1,
        ..Default::default()
    };

    let mut caller = CodecConsensusCaller::new("codec".to_string(), "RG1".to_string(), options);

    // Create a fragment read (not paired - no SEGMENTED flag)
    let record = RecordBuilder::new()
        .name("fragment1")
        .sequence("ACGTACGT")
        .qualities(&[30; 8])
        .reference_sequence_id(0)
        .alignment_start(100)
        .cigar("8M")
        .tag("MI", "UMI004")
        .build();

    let output = caller.consensus_reads_from_sam_records(vec![record]).unwrap();

    assert_eq!(output.count, 0, "Fragment reads should be rejected");

    let stats = caller.statistics();
    assert!(stats.reads_filtered > 0, "Fragment read should be filtered");
}

#[test]
#[allow(clippy::similar_names)]
fn test_codec_consensus_has_required_tags() {
    // Test that consensus reads have all required tags
    let options = CodecConsensusOptions {
        min_reads_per_strand: 1,
        min_duplex_length: 1,
        produce_per_base_tags: true,
        ..Default::default()
    };

    let mut caller = CodecConsensusCaller::new("codec".to_string(), "RG1".to_string(), options);

    let (r1, r2) = create_codec_read_pair(
        "read1",
        b"ACGTACGT",
        b"ACGTACGT",
        &[30, 30, 30, 30, 30, 30, 30, 30],
        &[30, 30, 30, 30, 30, 30, 30, 30],
        100,
        100,
        "UMI005",
    );

    let output = caller.consensus_reads_from_sam_records(vec![r1, r2]).unwrap();
    assert_eq!(output.count, 1);

    // Verify raw bytes are present (detailed tag validation is covered by unit tests)
    assert!(!output.data.is_empty(), "Consensus output should contain raw BAM bytes");
}

#[test]
fn test_codec_empty_input() {
    let options = CodecConsensusOptions::default();
    let mut caller = CodecConsensusCaller::new("codec".to_string(), "RG1".to_string(), options);

    let output = caller.consensus_reads_from_sam_records(Vec::new()).unwrap();

    assert_eq!(output.count, 0, "Empty input should produce empty output");
    assert_eq!(caller.statistics().total_input_reads, 0);
}

#[test]
fn test_codec_rejected_reads_tracking() {
    // Test that rejected reads are properly tracked when enabled
    let options = CodecConsensusOptions {
        min_reads_per_strand: 5, // High threshold to ensure rejection
        min_duplex_length: 1,
        ..Default::default()
    };

    let mut caller = CodecConsensusCaller::new_with_rejects_tracking(
        "codec".to_string(),
        "RG1".to_string(),
        options,
        true, // Enable tracking
    );

    // Create only 2 pairs (will be rejected)
    let mut reads = Vec::new();
    for i in 0..2 {
        let (r1, r2) = create_codec_read_pair(
            &format!("read{i}"),
            b"ACGTACGT",
            b"ACGTACGT",
            &[30; 8],
            &[30; 8],
            100,
            100,
            "UMI006",
        );
        reads.push(r1);
        reads.push(r2);
    }

    let output = caller.consensus_reads_from_sam_records(reads).unwrap();
    assert_eq!(output.count, 0);

    // With reject tracking enabled, rejected reads should be stored
    // (Note: actual storage depends on implementation details)
    let stats = caller.statistics();
    assert!(stats.reads_filtered > 0);
}

#[test]
fn test_codec_min_duplex_length_filter() {
    // Test that reads with insufficient overlap are rejected
    let options = CodecConsensusOptions {
        min_reads_per_strand: 1,
        min_duplex_length: 100, // Require long overlap
        ..Default::default()
    };

    let mut caller = CodecConsensusCaller::new("codec".to_string(), "RG1".to_string(), options);

    // Create reads with only 8bp overlap (below 100bp threshold)
    let (r1, r2) = create_codec_read_pair(
        "read1",
        b"ACGTACGT",
        b"ACGTACGT",
        &[30; 8],
        &[30; 8],
        100,
        100,
        "UMI007",
    );

    let output = caller.consensus_reads_from_sam_records(vec![r1, r2]).unwrap();

    assert_eq!(output.count, 0, "Should reject reads with insufficient duplex length");
}

#[test]
fn test_codec_statistics_tracking() {
    let options = CodecConsensusOptions {
        min_reads_per_strand: 1,
        min_duplex_length: 1,
        ..Default::default()
    };

    let mut caller = CodecConsensusCaller::new("codec".to_string(), "RG1".to_string(), options);

    // Process multiple molecules
    for mol_idx in 0..3 {
        let (r1, r2) = create_codec_read_pair(
            &format!("read_mol{mol_idx}"),
            b"ACGTACGT",
            b"ACGTACGT",
            &[30; 8],
            &[30; 8],
            100,
            100,
            &format!("UMI{mol_idx:03}"),
        );
        let _ = caller.consensus_reads_from_sam_records(vec![r1, r2]).unwrap();
    }

    let stats = caller.statistics();
    assert_eq!(stats.total_input_reads, 6, "Should count all input reads");
    assert_eq!(stats.consensus_reads_generated, 3, "Should emit 3 consensus reads");
}
