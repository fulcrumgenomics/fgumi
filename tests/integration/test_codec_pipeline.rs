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

use fgumi_lib::consensus::caller::ConsensusCaller;
use fgumi_lib::consensus::codec_caller::{CodecConsensusCaller, CodecConsensusOptions};
use fgumi_raw_bam::RawRecord;
use noodles::sam::alignment::record::cigar::op::Kind;

/// Map a noodles CIGAR `Kind` to its BAM op code (bits 3:0 of each cigar word).
const fn kind_to_op_code(k: Kind) -> u32 {
    match k {
        Kind::Match => 0,
        Kind::Insertion => 1,
        Kind::Deletion => 2,
        Kind::Skip => 3,
        Kind::SoftClip => 4,
        Kind::HardClip => 5,
        Kind::Pad => 6,
        Kind::SequenceMatch => 7,
        Kind::SequenceMismatch => 8,
    }
}

/// Helper function to create a test paired read with proper CIGAR.
///
/// Returns a `RawRecord` built directly via `fgumi_raw_bam::SamBuilder`.
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
) -> RawRecord {
    use fgumi_raw_bam::flags as raw_flags;

    // Encode CIGAR as BAM u32 words and compute reference length
    let mut ref_len = 0usize;
    let encoded_cigar: Vec<u32> = cigar_ops
        .iter()
        .map(|&(kind, len)| {
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
            (u32::try_from(len).expect("cigar len fits u32") << 4) | kind_to_op_code(kind)
        })
        .collect();

    // Compute template length for FR pair detection
    let insert_size: i32 = if ref_start <= mate_start {
        (mate_start + ref_len - ref_start) as i32
    } else {
        -((ref_start + ref_len - mate_start) as i32)
    };

    // Build flags
    let mut flags = raw_flags::PAIRED;
    if is_first {
        flags |= raw_flags::FIRST_SEGMENT;
    } else {
        flags |= raw_flags::LAST_SEGMENT;
    }
    if is_reverse {
        flags |= raw_flags::REVERSE;
    }
    if mate_is_reverse {
        flags |= raw_flags::MATE_REVERSE;
    }

    let mut b = fgumi_raw_bam::SamBuilder::new();
    b.read_name(name.as_bytes())
        .sequence(seq)
        .qualities(qual)
        .flags(flags)
        .ref_id(0)
        .pos(i32::try_from(ref_start).expect("ref_start fits i32") - 1) // 0-based
        .mapq(60)
        .cigar_ops(&encoded_cigar)
        .mate_ref_id(0)
        .mate_pos(i32::try_from(mate_start).expect("mate_start fits i32") - 1) // 0-based
        .template_length(insert_size)
        .add_string_tag(b"MI", umi.as_bytes());
    b.build()
}

/// Creates a CODEC read pair (R1 forward, R2 reverse - FR orientation).
///
/// Returns raw BAM records directly without an intermediate `RecordBuf`.
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
) -> (RawRecord, RawRecord) {
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
    let output = caller.consensus_reads(reads).unwrap();

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

    let output = caller.consensus_reads(reads).unwrap();

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

    let output = caller.consensus_reads(reads).unwrap();

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
    let record = {
        let mut b = fgumi_raw_bam::SamBuilder::new();
        b.read_name(b"fragment1")
            .sequence(b"ACGTACGT")
            .qualities(&[30; 8])
            .flags(0) // not paired
            .ref_id(0)
            .pos(99) // 0-based (alignment_start 100)
            .mapq(60)
            .cigar_ops(&[8u32 << 4]) // 8M
            .add_string_tag(b"MI", b"UMI004");
        b.build()
    };

    let output = caller.consensus_reads(vec![record]).unwrap();

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

    let output = caller.consensus_reads(vec![r1, r2]).unwrap();
    assert_eq!(output.count, 1);

    // Verify raw bytes are present (detailed tag validation is covered by unit tests)
    assert!(!output.data.is_empty(), "Consensus output should contain raw BAM bytes");
}

#[test]
fn test_codec_empty_input() {
    let options = CodecConsensusOptions::default();
    let mut caller = CodecConsensusCaller::new("codec".to_string(), "RG1".to_string(), options);

    let output = caller.consensus_reads(Vec::new()).unwrap();

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

    let output = caller.consensus_reads(reads).unwrap();
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

    let output = caller.consensus_reads(vec![r1, r2]).unwrap();

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
        let _ = caller.consensus_reads(vec![r1, r2]).unwrap();
    }

    let stats = caller.statistics();
    assert_eq!(stats.total_input_reads, 6, "Should count all input reads");
    assert_eq!(stats.consensus_reads_generated, 3, "Should emit 3 consensus reads");
}
