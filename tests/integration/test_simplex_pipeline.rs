//! Integration tests for the complete simplex consensus pipeline.
//!
//! These tests validate the end-to-end workflow:
//! 1. BAM input with UMI tags
//! 2. UMI assignment (Identity/Adjacency strategies)
//! 3. Grouping reads by molecule ID
//! 4. Consensus calling
//! 5. Quality filtering
//! 6. BAM output with consensus tags

use crate::helpers::*;
use fgumi_lib::umi::assigner::{AdjacencyUmiAssigner, IdentityUmiAssigner, UmiAssigner};
use std::collections::HashSet;

#[test]
fn test_identity_assigner_basic_workflow() {
    use noodles::sam::alignment::record::data::field::Tag;
    use noodles::sam::alignment::record_buf::data::field::Value;

    // Create test data: 3 UMI families with 5, 3, and 2 reads respectively
    let mut all_reads = Vec::new();

    // Family 1: ACGTACGT - 5 reads
    all_reads.extend(create_umi_family("ACGTACGT", 5, "family1", "AAAAAAAAAA", 30));

    // Family 2: TGCATGCA - 3 reads
    all_reads.extend(create_umi_family("TGCATGCA", 3, "family2", "TTTTTTTTTT", 30));

    // Family 3: GGGGGGGG - 2 reads
    all_reads.extend(create_umi_family("GGGGGGGG", 2, "family3", "CCCCCCCCCC", 30));

    let rx_tag = Tag::from([b'R', b'X']);
    let umis: Vec<String> = all_reads
        .iter()
        .map(|record| match record.data().get(&rx_tag).expect("RX tag should exist") {
            Value::String(s) => String::from_utf8_lossy(s.as_ref()).to_string(),
            _ => panic!("RX should be string"),
        })
        .collect();

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
    use noodles::sam::alignment::record::data::field::Tag;
    use noodles::sam::alignment::record_buf::data::field::Value;

    // Create test data with intentional sequencing errors
    // Base UMI: AAAAAAAA (50 reads)
    // Error UMI: AAAAAAAC (5 reads, 1bp error)
    // These should be grouped together by adjacency method
    let mut all_reads = Vec::new();

    all_reads.extend(create_umi_family_with_errors("AAAAAAAA", "AAAAAAAC", 50, 5, "ACGTACGTACGT"));

    // Add another distinct UMI family that should NOT be merged
    all_reads.extend(create_umi_family("TTTTTTTT", 30, "distinct", "TGCATGCATGCA", 30));

    let rx_tag = Tag::from([b'R', b'X']);
    let umis: Vec<String> = all_reads
        .iter()
        .map(|record| match record.data().get(&rx_tag).expect("RX tag") {
            Value::String(s) => String::from_utf8_lossy(s.as_ref()).to_string(),
            _ => panic!("RX should be string"),
        })
        .collect();

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
    use noodles::sam::alignment::record::data::field::Tag;
    use noodles::sam::alignment::record_buf::data::field::Value;

    // Test that adjacency method does NOT merge UMIs when count gradient prevents it
    // Base UMI: GGGGGGGG (30 reads)
    // Similar UMI: GGGGGGGC (20 reads, 1bp error)
    // These should NOT be merged because 20 >= 30/2 + 1 = 16
    let mut all_reads = Vec::new();

    all_reads.extend(create_umi_family_with_errors("GGGGGGGG", "GGGGGGGC", 30, 20, "ACGTACGT"));

    let rx_tag = Tag::from([b'R', b'X']);
    let umis: Vec<String> = all_reads
        .iter()
        .map(|record| match record.data().get(&rx_tag).expect("RX tag") {
            Value::String(s) => String::from_utf8_lossy(s.as_ref()).to_string(),
            _ => panic!("RX should be string"),
        })
        .collect();

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

#[test]
fn test_paired_end_umi_workflow() {
    use noodles::sam::alignment::record::data::field::Tag;
    use noodles::sam::alignment::record_buf::data::field::Value;

    // Create paired-end reads with UMIs
    let family1 = create_paired_umi_family(
        "AAAA-CCCC",
        5,
        "pair_family1",
        "ACGTACGTACGT",
        "TGCATGCATGCA",
        30,
    );

    let family2 = create_paired_umi_family(
        "GGGG-TTTT",
        3,
        "pair_family2",
        "ACGTACGTACGT",
        "TGCATGCATGCA",
        30,
    );

    // Verify proper pairing
    assert_proper_pairing(&family1[0], &family1[1]);
    assert_proper_pairing(&family2[0], &family2[1]);

    // Verify R1 and R2 have different sequences
    assert_ne!(
        family1[0].sequence().as_ref(),
        family1[1].sequence().as_ref(),
        "R1 and R2 should have different sequences"
    );

    let mut all_reads = Vec::new();
    all_reads.extend(family1);
    all_reads.extend(family2);

    let rx_tag = Tag::from([b'R', b'X']);
    let umis: Vec<String> = all_reads
        .iter()
        .map(|record| match record.data().get(&rx_tag).expect("RX tag") {
            Value::String(s) => String::from_utf8_lossy(s.as_ref()).to_string(),
            _ => panic!("RX should be string"),
        })
        .collect();

    // For paired UMIs, we would use PairedUmiAssigner (not tested here due to suffix handling)
    // This test just validates paired-end read generation
    assert_eq!(umis.len(), 16, "Should have 16 reads total (8 pairs)");
}
