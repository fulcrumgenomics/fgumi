//! Error path integration tests.
//!
//! These tests verify that error conditions are handled correctly,
//! including validation failures, missing files, and invalid inputs.

use fgumi_lib::consensus::{ConsensusSequence, FilterConfig, FilterThresholds};

// ==================== ConsensusSequence Error Paths ====================

#[test]
fn test_consensus_sequence_empty_operations() {
    let seq = ConsensusSequence::new();

    // Operations on empty sequence should be safe
    assert_eq!(seq.max_depth(), 0);
    assert_eq!(seq.min_depth(), 0);
    assert!((seq.error_rate() - 0.0).abs() < f32::EPSILON);
    assert!(seq.bases().is_empty());
    assert!(seq.quals().is_empty());
}

#[test]
#[should_panic(expected = "new_length")]
fn test_consensus_sequence_padded_shorter_panics() {
    let seq = ConsensusSequence::from_vecs(
        vec![b'A', b'C', b'G', b'T'],
        vec![30, 30, 30, 30],
        vec![10, 10, 10, 10],
        vec![0, 0, 0, 0],
    );

    // Trying to pad to shorter length should panic
    let _ = seq.padded(2, false, b'N', 2);
}

// ==================== FilterConfig Error Paths ====================

#[test]
#[should_panic(expected = "min-reads values must be specified high to low")]
fn test_filter_config_for_duplex_invalid_strand_greater_than_duplex() {
    let duplex =
        FilterThresholds { min_reads: 5, max_read_error_rate: 0.1, max_base_error_rate: 0.2 };
    let strand =
        FilterThresholds { min_reads: 10, max_read_error_rate: 0.1, max_base_error_rate: 0.2 };

    // strand.min_reads (10) > duplex.min_reads (5) should panic
    let _ = FilterConfig::for_duplex(duplex, strand, Some(13), None, 0.1);
}

#[test]
#[should_panic(expected = "max-read-error-rate for AB")]
fn test_filter_config_for_duplex_asymmetric_invalid_error_rate() {
    let duplex =
        FilterThresholds { min_reads: 10, max_read_error_rate: 0.05, max_base_error_rate: 0.1 };
    let ab = FilterThresholds { min_reads: 5, max_read_error_rate: 0.2, max_base_error_rate: 0.1 };
    let ba = FilterThresholds { min_reads: 3, max_read_error_rate: 0.1, max_base_error_rate: 0.2 };

    // ab.max_read_error_rate (0.2) > ba.max_read_error_rate (0.1) should panic
    let _ = FilterConfig::for_duplex_asymmetric(duplex, ab, ba, Some(13), None, 0.1);
}

// Note: ConsensusCallingOptions validation tests are in the binary crate's
// unit tests (src/commands/common.rs) since they require access to the
// command module which is not exposed from the library.
