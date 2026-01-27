//! Integration tests for fgumi.
//!
//! Run with: `cargo test --test integration_tests`
//!
//! These tests validate end-to-end workflows spanning multiple modules.

#![allow(clippy::cast_precision_loss)]

use fgumi_lib::logging::{format_duration, format_percent, format_rate};
use fgumi_lib::metrics::{ConsensusMetrics, FamilySizeMetrics, UmiGroupingMetrics};
use fgumi_lib::umi::assigner::{AdjacencyUmiAssigner, IdentityUmiAssigner, UmiAssigner};
use std::collections::HashSet;
use std::time::Duration;

/// Helper to create a simple test UMI vector
fn create_test_umis_simple() -> Vec<String> {
    vec![
        "ACGTACGT".to_string(),
        "ACGTACGT".to_string(),
        "ACGTACGT".to_string(),
        "TGCATGCA".to_string(),
        "TGCATGCA".to_string(),
        "GGGGGGGG".to_string(),
    ]
}

#[test]
fn test_identity_assigner_basic_workflow() {
    // Create test data: 3 unique UMIs
    let umis = create_test_umis_simple();

    // Apply identity UMI assigner
    let assigner = IdentityUmiAssigner::default();
    let assignments = assigner.assign(&umis);

    // Helper to find the MoleculeId for a given UMI string
    let get_assignment = |umi: &str| -> Option<&fgumi_lib::template::MoleculeId> {
        umis.iter().position(|u| u == umi).map(|idx| &assignments[idx])
    };

    // Verify correct number of molecule IDs
    let unique_ids: HashSet<_> = assignments.iter().collect();
    assert_eq!(unique_ids.len(), 3, "Should have 3 unique molecule IDs for 3 different UMIs");

    // Verify all reads with same UMI got same molecule ID
    assert_eq!(
        get_assignment("ACGTACGT"),
        get_assignment("ACGTACGT"),
        "Same UMI should get same molecule ID"
    );

    // Verify different UMIs got different molecule IDs
    assert_ne!(
        get_assignment("ACGTACGT"),
        get_assignment("TGCATGCA"),
        "Different UMIs should get different molecule IDs"
    );
}

#[test]
fn test_adjacency_assigner_with_error_correction() {
    // Create test data with intentional sequencing errors
    let mut umis = Vec::new();

    // Base UMI: AAAAAAAA (50 reads)
    for _ in 0..50 {
        umis.push("AAAAAAAA".to_string());
    }

    // Error UMI: AAAAAAAC (5 reads, 1bp error)
    for _ in 0..5 {
        umis.push("AAAAAAAC".to_string());
    }

    // Distinct UMI: TTTTTTTT (30 reads)
    for _ in 0..30 {
        umis.push("TTTTTTTT".to_string());
    }

    // Apply adjacency assigner with 1 mismatch allowed
    let assigner = AdjacencyUmiAssigner::new(1, 1, 100);
    let assignments = assigner.assign(&umis);

    // Helper to find the MoleculeId for a given UMI string
    let get_assignment = |umi: &str| -> Option<&fgumi_lib::template::MoleculeId> {
        umis.iter().position(|u| u == umi).map(|idx| &assignments[idx])
    };

    // Verify error correction: AAAAAAAA and AAAAAAAC should be merged
    // because 5 < 50/2 + 1 = 26 (count gradient rule)
    assert_eq!(
        get_assignment("AAAAAAAA"),
        get_assignment("AAAAAAAC"),
        "Error UMI should be captured by high-abundance base UMI"
    );

    // Verify TTTTTTTT gets separate molecule ID
    assert_ne!(
        get_assignment("AAAAAAAA"),
        get_assignment("TTTTTTTT"),
        "Distinct UMI family should not be merged"
    );

    // Should have 2 unique molecule IDs total
    let unique_ids: HashSet<_> = assignments.iter().collect();
    assert_eq!(
        unique_ids.len(),
        2,
        "Should have 2 molecule IDs: merged AAAA* family and TTTT* family"
    );
}

#[test]
fn test_adjacency_respects_count_gradient() {
    // Test that adjacency method does NOT merge UMIs when count gradient prevents it
    let mut umis = Vec::new();

    // Base UMI: GGGGGGGG (30 reads)
    for _ in 0..30 {
        umis.push("GGGGGGGG".to_string());
    }

    // Similar UMI: GGGGGGGC (20 reads, 1bp error)
    // These should NOT be merged because 20 >= 30/2 + 1 = 16
    for _ in 0..20 {
        umis.push("GGGGGGGC".to_string());
    }

    // Apply adjacency assigner
    let assigner = AdjacencyUmiAssigner::new(1, 1, 100);
    let assignments = assigner.assign(&umis);

    // Helper to find the MoleculeId for a given UMI string
    let get_assignment = |umi: &str| -> Option<&fgumi_lib::template::MoleculeId> {
        umis.iter().position(|u| u == umi).map(|idx| &assignments[idx])
    };

    // Verify they are NOT merged (count gradient prevents it)
    assert_ne!(
        get_assignment("GGGGGGGG"),
        get_assignment("GGGGGGGC"),
        "Similar-abundance UMIs should not be merged due to count gradient"
    );

    // Should have 2 unique molecule IDs
    let unique_ids: HashSet<_> = assignments.iter().collect();
    assert_eq!(unique_ids.len(), 2, "Should have 2 separate molecule IDs due to count gradient");
}

// Metrics Integration Tests

#[test]
fn test_umi_grouping_metrics_workflow() {
    let mut metrics = UmiGroupingMetrics::new();

    // Simulate UMI grouping workflow
    metrics.total_records = 10_000;
    metrics.accepted_records = 9_500;
    metrics.discarded_non_pf = 200;
    metrics.discarded_ns_in_umi = 150;
    metrics.discarded_umi_too_short = 150;
    metrics.unique_molecule_ids = 1_000;
    metrics.total_families = 1_000;
    metrics.avg_reads_per_molecule = 9.5;

    assert_eq!(metrics.total_records, 10_000);
    assert_eq!(metrics.unique_molecule_ids, 1_000);
    assert!((metrics.avg_reads_per_molecule - 9.5).abs() < f64::EPSILON);
}

#[test]
fn test_family_size_metrics_distribution() {
    let mut size_1 = FamilySizeMetrics::new(1);
    let mut size_2 = FamilySizeMetrics::new(2);
    let mut size_5 = FamilySizeMetrics::new(5);

    size_1.count = 100;
    size_1.fraction = 0.5;
    size_1.fraction_gt_or_eq_family_size = 1.0;

    size_2.count = 50;
    size_2.fraction = 0.25;
    size_2.fraction_gt_or_eq_family_size = 0.5;

    size_5.count = 50;
    size_5.fraction = 0.25;
    size_5.fraction_gt_or_eq_family_size = 0.25;

    assert_eq!(size_1.family_size, 1);
    assert_eq!(size_2.family_size, 2);
    assert_eq!(size_5.family_size, 5);
}

#[test]
fn test_metrics_basic_fields() {
    // Test that metrics fields can be set and read
    let mut metrics = ConsensusMetrics::new();
    metrics.total_input_reads = 1000;
    metrics.consensus_reads = 800;
    metrics.filtered_reads = 200;

    assert_eq!(metrics.total_input_reads, 1000);
    assert_eq!(metrics.consensus_reads, 800);
    assert_eq!(metrics.filtered_reads, 200);
}

// Logging Integration Tests

#[test]
fn test_format_percent_integration() {
    // Test percentage formatting in realistic scenarios
    let pass_rate = 0.9543;
    assert_eq!(format_percent(pass_rate, 2), "95.43%");

    let low_rate = 0.0123;
    assert_eq!(format_percent(low_rate, 2), "1.23%");

    let perfect_rate = 1.0;
    assert_eq!(format_percent(perfect_rate, 1), "100.0%");
}

#[test]
fn test_format_duration_realistic() {
    // Test duration formatting for typical processing times
    let short_job = Duration::from_secs(45);
    assert_eq!(format_duration(short_job), "45s");

    let medium_job = Duration::from_secs(125);
    assert_eq!(format_duration(medium_job), "2m 5s");

    let long_job = Duration::from_secs(7200); // Exactly 2 hours
    assert_eq!(format_duration(long_job), "2h");
}

#[test]
fn test_format_rate_with_realistic_data() {
    // Simulate processing 100k reads in 10 seconds
    let count = 100_000;
    let duration = Duration::from_secs(10);
    let rate = format_rate(count, duration);
    assert!(rate.contains("10,000 items/s"));

    // Simulate slow processing
    let slow_count = 50;
    let slow_duration = Duration::from_secs(60);
    let slow_rate = format_rate(slow_count, slow_duration);
    assert!(slow_rate.contains("items/min"));
}
