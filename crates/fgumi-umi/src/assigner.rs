// Copyright (c) 2025, Fulcrum Genomics LLC
// All rights reserved.

//! # UMI Assignment Strategies
//!
//! This module provides algorithms for grouping sequencing reads by their Unique Molecular
//! Identifiers (UMIs). UMIs are short random sequences added during library preparation to
//! tag individual DNA/RNA molecules, enabling accurate removal of PCR duplicates and
//! quantification of original molecules.
//!
//! ## Overview
//!
//! UMI assignment is the process of determining which reads originate from the same source
//! molecule by comparing their UMI sequences. This is complicated by sequencing errors in
//! the UMIs themselves, which can cause reads from the same molecule to appear to have
//! different UMIs. The strategies in this module balance error correction against
//! over-correction (incorrectly merging distinct molecules).
//!
//! ## Algorithms
//!
//! Four complementary strategies are provided, each with different error tolerance and
//! computational characteristics:
//!
//! ### Identity Strategy
//!
//! Groups only reads with identical UMI sequences (case-insensitive). This is the fastest
//! strategy but provides no error correction. Each unique UMI sequence receives a distinct
//! molecule ID.
//!
//! **Algorithm**: Direct hash-based grouping after uppercasing.
//!
//! **Complexity**: O(n) time, O(u) space where n = read count, u = unique UMI count.
//!
//! **When to use**: High-quality UMIs (Q30+), maximum specificity required, or as a
//! baseline for comparison.
//!
//! ### Edit Distance Strategy
//!
//! Groups UMIs that are within a specified Hamming distance (mismatches). Uses iterative
//! clustering where UMIs matching any member of a group are added to that group. This can
//! lead to transitive clustering: if UMI A matches B and B matches C, all three are grouped
//! even if A and C differ by more than the threshold.
//!
//! **Algorithm**:
//! 1. For each UMI, find all existing groups containing at least one UMI within threshold
//! 2. If no matches, create a new group
//! 3. If one match, add UMI to that group
//! 4. If multiple matches, merge all matching groups
//!
//! **Complexity**: O(n * u * m) worst case, where m = UMI length.
//!
//! **When to use**: Moderate error rates, smaller datasets, or when transitive clustering
//! is acceptable.
//!
//! ### Adjacency Strategy
//!
//! Implements the directed adjacency graph method from Smith et al. (2017). This method
//! prevents over-correction by using a count gradient: a UMI can only be "captured" by
//! another UMI if (a) they are within edit distance and (b) the child UMI has count
//! strictly less than `parent_count / 2 + 1`. This ensures that two similarly abundant
//! UMIs are never merged.
//!
//! **Algorithm**:
//! 1. Count UMIs and sort by decreasing abundance
//! 2. Build directed graph where edges represent "captures":
//!    - Process UMIs in order of abundance (most abundant first)
//!    - Each unassigned UMI becomes a root node
//!    - Search remaining unassigned UMIs for potential children
//!    - Add directed edge if distance <= threshold AND `child_count` < `parent_count/2` + 1
//!    - Recursively process children as potential parents
//! 3. Assign same molecule ID to each root and all its descendants
//!
//! **Mathematical basis**: The count gradient condition approximates a 2:1 evidence ratio.
//! If a true UMI appears 100 times and an error appears 45 times, they will merge (45 < 51).
//! But if both appear 50 times, they stay separate, preventing false merges of distinct
//! molecules that happen to have similar UMIs.
//!
//! **Complexity**: O(u^2 * m) worst case, O(u log u + u * k * m) typical case where
//! k << u is the average number of candidates within count threshold. Multi-threaded
//! implementation available for better performance.
//!
//! **When to use**: High-depth sequencing (>10x per molecule), moderate to high error
//! rates, or when conservative grouping is critical (e.g., variant calling).
//!
//! ### Paired Strategy
//!
//! Extends the adjacency method for paired UMIs in duplex sequencing. Paired UMIs have
//! the form "AAAA-CCCC" where the two parts come from opposite ends of the original
//! molecule. Critically, reads from opposite strands have reversed UMI pairs: "AAAA-CCCC"
//! and "CCCC-AAAA" represent the same molecule.
//!
//! **Algorithm**:
//! 1. Canonicalize UMIs (use lexically lower of A-B and B-A) for counting
//! 2. Build adjacency graph on canonical forms (as in Adjacency strategy)
//! 3. For each root and its descendants:
//!    - Compare each UMI to root to determine strand assignment
//!    - If closer to root, assign /A suffix; if closer to reverse, assign /B suffix
//! 4. Map both orientations to appropriate strand-specific molecule IDs
//!
//! **Matching**: When checking if UMIs are adjacent, considers both A-B and B-A forms
//! (up to 2x edits allowed across the pair).
//!
//! **Complexity**: Similar to Adjacency, with 2x cost for reverse comparisons.
//!
//! **When to use**: Duplex sequencing libraries where both strands are sequenced
//! independently, enabling strand-aware consensus calling and error correction.
//!
//! ## Algorithm References
//!
//! The Adjacency and Paired strategies implement the directed adjacency graph method
//! described in:
//!
//! - Smith, T., Heger, A., & Sudbery, I. (2017). "UMI-tools: modeling sequencing errors
//!   in Unique Molecular Identifiers to improve quantification accuracy."
//!   *Genome Research*, 27(3), 491-499.
//!   DOI: [10.1101/gr.209601.116](https://doi.org/10.1101/gr.209601.116)
//!
//! The count gradient approach (`child_count < parent_count / 2 + 1`) is key to
//! preventing over-correction while still allowing error correction.
//!
//! ## Usage Examples
//!
//! ### Basic Usage - Identity Strategy
//!
//! ```
//! use fgumi_lib::umi::assigner::{Strategy, UmiAssigner};
//!
//! let strategy = Strategy::Identity;
//! let assigner = strategy.new_assigner(0); // No mismatches allowed
//!
//! let umis = vec![
//!     "ACGT".to_string(),
//!     "ACGT".to_string(),  // Same as first
//!     "ACTT".to_string(),  // Different - will not group with first
//! ];
//!
//! let assignments = assigner.assign(&umis);
//! // Returns a HashMap mapping each input UMI to its molecule ID
//! // Two ACGT UMIs map to one ID, one ACTT UMI maps to a different ID
//! use std::collections::HashSet;
//! let unique_ids: HashSet<_> = assignments.iter().collect();
//! assert_eq!(unique_ids.len(), 2); // Two unique molecule IDs
//! ```
//!
//! ### Error-Tolerant Assignment - Adjacency Strategy
//!
//! ```
//! use fgumi_lib::umi::assigner::{Strategy, UmiAssigner};
//!
//! let strategy = Strategy::Adjacency;
//! let assigner = strategy.new_assigner(1); // Allow 1 mismatch
//!
//! let umis = vec![
//!     "ACGT".to_string(),
//!     "ACGT".to_string(),
//!     "ACGT".to_string(),
//!     "ACTT".to_string(),  // 1 mismatch from ACGT, low abundance
//! ];
//!
//! let assignments = assigner.assign(&umis);
//! // ACGT appears 3 times, ACTT appears 1 time
//! // Since 1 < 3/2 + 1 = 2.5, ACTT will be captured by ACGT
//! // All reads get the same molecule ID
//! ```
//!
//! ### Paired UMI Assignment - Duplex Sequencing
//!
//! ```
//! use fgumi_lib::umi::assigner::{Strategy, UmiAssigner, TOP_STRAND_DUPLEX, BOTTOM_STRAND_DUPLEX};
//!
//! let strategy = Strategy::Paired;
//! let assigner = strategy.new_assigner(1); // Allow 1 mismatch per UMI half
//!
//! let umis = vec![
//!     "AAAA-CCCC".to_string(),  // Top strand read
//!     "CCCC-AAAA".to_string(),  // Bottom strand read (same molecule!)
//! ];
//!
//! let assignments = assigner.assign(&umis);
//! // Both UMIs get same base molecule ID but different suffixes:
//! // "AAAA-CCCC" -> "0/A"
//! // "CCCC-AAAA" -> "0/B"
//! // This enables separate consensus calling per strand before duplex consensus
//! ```
//!
//! ## Performance Characteristics
//!
//! | Strategy   | Time Complexity      | Space | Best For                    |
//! |------------|---------------------|-------|------------------------------|
//! | Identity   | O(n)                | O(u)  | Q30+ UMIs, maximum speed     |
//! | Edit       | O(n * u * m)        | O(u)  | Small datasets, simple logic |
//! | Adjacency  | O(u^2 * m) worst    | O(u)  | High depth, conservative     |
//! |            | O(u log u) typical  |       | grouping                     |
//! | Paired     | O(u^2 * m) worst    | O(u)  | Duplex sequencing           |
//!
//! Where: n = total reads, u = unique UMIs, m = UMI length
//!
//! **Multi-threading**: Adjacency and Paired strategies support multi-threading via
//! `new_assigner_with_threads()`, which parallelizes the candidate matching step.
//!
//! ## Choosing a Strategy
//!
//! - **Use Identity** if: UMI quality is very high (Q35+), you want maximum speed, or
//!   you need to exactly match reference implementations.
//!
//! - **Use Edit** if: You have moderate error rates, smaller datasets (<100K UMIs per
//!   position), and transitive clustering is acceptable for your application.
//!
//! - **Use Adjacency** if: You have high sequencing depth (10-100x+ per molecule),
//!   need conservative error correction, or are doing variant calling where
//!   over-correction could introduce false negatives.
//!
//! - **Use Paired** if: You are processing duplex sequencing data where both strands
//!   are tagged with dual UMIs and you need strand-aware consensus calling.

use ahash::AHashMap;
use anyhow::Result;
use rayon::ThreadPoolBuilder;
use rayon::prelude::*;
use std::any::Any;
use std::collections::{BTreeMap, BTreeSet, HashMap, HashSet, VecDeque};
use std::sync::atomic::{AtomicU64, Ordering};

use fgumi_dna::BitEnc;
use crate::MoleculeId;

/// Default minimum number of unique UMIs to use N-gram/BK-tree index for candidate search.
/// Below this threshold, linear scan is used. Based on benchmarks: index provides ~10%
/// speedup even at 100 UMIs/position with negligible construction overhead.
pub const DEFAULT_INDEX_THRESHOLD: usize = 100;

/// A node in the BK-tree.
///
/// Each node stores a UMI's `BitEnc` encoding and its index in the original list.
/// Children are organized by their Hamming distance to this node's UMI.
struct BkNode {
    /// `BitEnc` encoding for fast Hamming distance
    encoded: BitEnc,
    /// Index in the original UMI list
    index: usize,
    /// Children indexed by Hamming distance
    children: BTreeMap<u32, Box<BkNode>>,
}

/// BK-tree for efficient Hamming distance queries on UMIs.
///
/// A BK-tree (Burkhard-Keller tree) is a metric tree that exploits the triangle
/// inequality to efficiently find all strings within a given edit distance.
/// For Hamming distance queries with threshold k, it achieves O(k * log n)
/// average case complexity instead of O(n) for linear scan.
///
/// This implementation uses `BitEnc` for O(1) Hamming distance computation.
///
/// # Algorithm
///
/// The tree is built by inserting UMIs one by one:
/// 1. The first UMI becomes the root
/// 2. For each subsequent UMI, traverse from root following edges labeled
///    with the Hamming distance to each node's UMI
/// 3. Insert as a new child when no edge exists for that distance
///
/// To query for all UMIs within distance k of a query:
/// 1. Compute distance d from query to current node
/// 2. If d <= k, include this node in results
/// 3. Recursively search children with edge labels in range [d-k, d+k]
///    (triangle inequality guarantees other children can't contain matches)
pub struct BkTree {
    root: Option<Box<BkNode>>,
    size: usize,
}

impl BkTree {
    /// Create a new empty BK-tree.
    #[must_use]
    pub fn new() -> Self {
        Self { root: None, size: 0 }
    }

    /// Build a BK-tree from a list of UMIs with their indices.
    ///
    /// Returns None if any UMI contains non-ACGT bases (cannot be encoded).
    #[must_use]
    pub fn from_umis(umis: &[(usize, &str)]) -> Option<Self> {
        let mut tree = Self::new();
        for &(idx, umi) in umis {
            let encoded = BitEnc::from_bytes(umi.as_bytes())?;
            tree.insert(encoded, idx);
        }
        Some(tree)
    }

    /// Insert a UMI with its index into the tree.
    fn insert(&mut self, encoded: BitEnc, index: usize) {
        self.size += 1;
        match self.root {
            Some(ref mut root) => Self::insert_into(root, encoded, index),
            None => {
                self.root =
                    Some(Box::new(BkNode { encoded, index, children: BTreeMap::new() }));
            }
        }
    }

    /// Recursively insert a UMI into the subtree rooted at `node`.
    fn insert_into(node: &mut BkNode, encoded: BitEnc, index: usize) {
        let dist = node.encoded.hamming_distance(&encoded);
        match node.children.get_mut(&dist) {
            Some(child) => Self::insert_into(child, encoded, index),
            None => {
                node.children
                    .insert(dist, Box::new(BkNode { encoded, index, children: BTreeMap::new() }));
            }
        }
    }

    /// Find all UMIs within the given Hamming distance of the query.
    ///
    /// Returns a vector of (index, distance) pairs for all matching UMIs.
    /// Returns empty if query cannot be encoded (contains non-ACGT bases).
    #[must_use]
    pub fn find_within(&self, query: &str, max_distance: u32) -> Vec<(usize, u32)> {
        let Some(query_enc) = BitEnc::from_bytes(query.as_bytes()) else {
            return Vec::new();
        };

        let mut results = Vec::new();
        if let Some(ref root) = self.root {
            Self::find_within_recursive(root, &query_enc, max_distance, &mut results);
        }
        results
    }

    fn find_within_recursive(
        node: &BkNode,
        query: &BitEnc,
        max_distance: u32,
        results: &mut Vec<(usize, u32)>,
    ) {
        let dist = node.encoded.hamming_distance(query);

        if dist <= max_distance {
            results.push((node.index, dist));
        }

        // Triangle inequality: only search children with distances in [dist-k, dist+k]
        let min_child_dist = dist.saturating_sub(max_distance);
        let max_child_dist = dist + max_distance;

        for (&child_dist, child) in &node.children {
            if child_dist >= min_child_dist && child_dist <= max_child_dist {
                Self::find_within_recursive(child, query, max_distance, results);
            }
        }
    }

    /// Returns the number of UMIs in the tree.
    #[must_use]
    pub fn len(&self) -> usize {
        self.size
    }

    /// Returns true if the tree is empty.
    #[must_use]
    pub fn is_empty(&self) -> bool {
        self.size == 0
    }
}

impl Default for BkTree {
    fn default() -> Self {
        Self::new()
    }
}

/// N-gram index for fast candidate filtering using the pigeonhole principle.
///
/// For UMIs of length M and max mismatches k, if two UMIs differ by at most k positions,
/// at least one of their k+1 equal-length partitions must match exactly. This allows
/// using hash table lookups instead of tree traversal for candidate generation.
///
/// # Algorithm
///
/// 1. Split each UMI into k+1 non-overlapping partitions of length ⌊M/(k+1)⌋
/// 2. Build k+1 hash tables mapping partition values to UMI indices
/// 3. To query: look up query's partitions in each hash table, union candidates
/// 4. Verify each candidate with actual Hamming distance
///
/// # Complexity
///
/// - Build: O(n) where n = number of UMIs
/// - Query: O(k + c) where k = max mismatches, c = candidates to verify
///
/// For typical parameters (12bp UMIs, k=1), this gives O(1) expected query time
/// due to hash table lookups, compared to O(log n) for BK-tree.
pub struct NgramIndex {
    /// Number of partitions (k+1)
    num_partitions: usize,
    /// Length of each partition in bases
    partition_len: usize,
    /// For each partition position: `partition_bits` → `Vec<(index, BitEnc)>`
    /// Uses `AHashMap` for faster lookups than std `HashMap`
    partitions: Vec<AHashMap<u32, Vec<(usize, BitEnc)>>>,
}

impl NgramIndex {
    /// Build an N-gram index from UMIs with their indices.
    ///
    /// Returns None if:
    /// - The UMI list is empty
    /// - Any UMI contains non-ACGT bases
    /// - UMIs have inconsistent lengths
    ///
    /// # Arguments
    ///
    /// * `umis` - Slice of (index, UMI string) pairs
    /// * `max_mismatches` - Maximum allowed mismatches (determines partition count)
    #[must_use]
    pub fn new(umis: &[(usize, &str)], max_mismatches: u32) -> Option<Self> {
        let first_umi = umis.first()?.1;
        let umi_len = first_umi.len();
        let num_partitions = (max_mismatches + 1) as usize;
        let partition_len = umi_len / num_partitions;

        // Need at least 1 base per partition for meaningful filtering
        if partition_len == 0 {
            return None;
        }

        let mut partitions: Vec<AHashMap<u32, Vec<(usize, BitEnc)>>> =
            (0..num_partitions).map(|_| AHashMap::new()).collect();

        for &(idx, umi) in umis {
            // Verify consistent length
            if umi.len() != umi_len {
                return None;
            }

            let full_enc = BitEnc::from_bytes(umi.as_bytes())?;

            for (p, partition_map) in partitions.iter_mut().enumerate() {
                let start_base = p * partition_len;
                let partition_bits = full_enc.extract_bits(start_base, partition_len);

                partition_map.entry(partition_bits).or_default().push((idx, full_enc));
            }
        }

        Some(Self { num_partitions, partition_len, partitions })
    }

    /// Find all UMIs within the given Hamming distance of the query.
    ///
    /// Returns a vector of (index, distance) pairs for all matching UMIs.
    /// Returns empty if query cannot be encoded or has wrong length.
    #[must_use]
    pub fn find_within(&self, query: &str, max_distance: u32) -> Vec<(usize, u32)> {
        let Some(query_enc) = BitEnc::from_bytes(query.as_bytes()) else {
            return Vec::new();
        };

        // Verify query length matches index
        let expected_len = self.num_partitions * self.partition_len;
        if query.len() < expected_len {
            return Vec::new();
        }

        let mut seen = HashSet::new();
        let mut results = Vec::new();

        // Check each partition - pigeonhole guarantees at least one will match
        for (p, partition_map) in self.partitions.iter().enumerate() {
            let start_base = p * self.partition_len;
            let query_partition = query_enc.extract_bits(start_base, self.partition_len);

            if let Some(candidates) = partition_map.get(&query_partition) {
                for &(idx, ref full_enc) in candidates {
                    // Deduplicate candidates found via multiple partitions
                    if seen.insert(idx) {
                        let dist = query_enc.hamming_distance(full_enc);
                        if dist <= max_distance {
                            results.push((idx, dist));
                        }
                    }
                }
            }
        }

        results
    }
}

/// Type representing a UMI sequence
pub type Umi = String;

// Note: MoleculeId enum is imported from crate::template

/// Suffix for top strand duplex molecules
///
/// Used by the `Paired` strategy to distinguish reads originating from the top strand
/// of duplex molecules. Top strand reads have their molecule ID suffixed with "/A".
pub const TOP_STRAND_DUPLEX: &str = "/A";

/// Suffix for bottom strand duplex molecules  
///
/// Used by the `Paired` strategy to distinguish reads originating from the bottom strand
/// of duplex molecules. Bottom strand reads have their molecule ID suffixed with "/B".
pub const BOTTOM_STRAND_DUPLEX: &str = "/B";

/// The UMI assignment strategy
///
/// Determines how reads are grouped based on their UMI sequences. Each strategy makes
/// different tradeoffs between speed, error tolerance, and grouping behavior.
#[derive(Debug, Clone, Copy)]
#[cfg_attr(feature = "cli", derive(clap::ValueEnum))]
pub enum Strategy {
    /// Only reads with identical UMI sequences are grouped together
    ///
    /// The fastest strategy but provides no error correction. Two UMIs that differ by even
    /// a single base will be assigned to different molecules. Best used when UMI quality
    /// is very high or when maximum specificity is required.
    Identity,

    /// Cluster reads based on mismatches between UMIs
    ///
    /// Groups UMIs that are within a specified edit distance. Uses a simple clustering
    /// algorithm that can merge groups transitively. For example, if A matches B and B
    /// matches C, all three will be grouped together even if A and C differ by more than
    /// the threshold. More permissive than adjacency method.
    Edit,

    /// Directed adjacency method from `umi_tools` that allows errors with count gradient
    ///
    /// Implements the directed adjacency graph method where high-abundance UMIs can "capture"
    /// lower-abundance UMIs that are within edit distance. A UMI can only be assigned to
    /// another UMI if its count is less than half of the parent's count. This prevents
    /// over-correction and is more conservative than the simple edit distance method.
    Adjacency,

    /// Similar to adjacency but for paired UMIs (A-B related to B-A)
    ///
    /// Extends the adjacency method for paired UMIs (e.g., "ACGT-TGCA"). Treats A-B and
    /// B-A as related molecules from opposite strands of duplex sequencing. Assigns
    /// strand-specific molecule IDs with /A and /B suffixes to enable duplex consensus
    /// calling.
    Paired,
}

impl Strategy {
    /// Create a new UMI assigner for this strategy
    ///
    /// # Arguments
    ///
    /// * `edits` - Maximum number of mismatches allowed between UMIs. Must be 0 for Identity strategy.
    ///
    /// # Returns
    ///
    /// A boxed trait object implementing the `UmiAssigner` trait.
    ///
    /// # Panics
    ///
    /// Panics if `edits` is non-zero when using the Identity strategy.
    #[must_use]
    pub fn new_assigner(&self, edits: u32) -> Box<dyn UmiAssigner> {
        self.new_assigner_full(edits, 1, DEFAULT_INDEX_THRESHOLD)
    }

    /// Create a new UMI assigner with thread count specified
    ///
    /// # Arguments
    ///
    /// * `edits` - Maximum number of mismatches allowed between UMIs
    /// * `threads` - Number of threads to use (only applies to Adjacency and Paired strategies)
    ///
    /// # Returns
    ///
    /// A boxed trait object implementing the `UmiAssigner` trait
    #[must_use]
    pub fn new_assigner_with_threads(&self, edits: u32, threads: usize) -> Box<dyn UmiAssigner> {
        self.new_assigner_full(edits, threads, DEFAULT_INDEX_THRESHOLD)
    }

    /// Create a new UMI assigner with all parameters specified
    ///
    /// # Arguments
    ///
    /// * `edits` - Maximum number of mismatches allowed between UMIs
    /// * `threads` - Number of threads to use (only applies to Adjacency and Paired strategies)
    /// * `index_threshold` - Minimum UMIs per position to use N-gram/BK-tree index
    ///
    /// # Returns
    ///
    /// A boxed trait object implementing the `UmiAssigner` trait
    ///
    /// # Panics
    ///
    /// Panics if `edits` is non-zero with the `Identity` strategy.
    #[must_use]
    pub fn new_assigner_full(
        &self,
        edits: u32,
        threads: usize,
        index_threshold: usize,
    ) -> Box<dyn UmiAssigner> {
        match self {
            Strategy::Identity => {
                assert_eq!(edits, 0, "Edits should be zero when using the identity UMI assigner");
                Box::new(IdentityUmiAssigner::new())
            }
            Strategy::Edit => Box::new(SimpleErrorUmiAssigner::new(edits)),
            Strategy::Adjacency => {
                Box::new(AdjacencyUmiAssigner::new(edits, threads, index_threshold))
            }
            Strategy::Paired => {
                Box::new(PairedUmiAssigner::new_with_threads(edits, threads, index_threshold))
            }
        }
    }
}

/// Count mismatches between two sequences
///
/// Compares two sequences position by position and counts the number of positions where
/// the characters differ. If the sequences have different lengths, returns `usize::MAX`.
///
/// # Arguments
///
/// * `a` - First sequence
/// * `b` - Second sequence
///
/// # Returns
///
/// Number of mismatched positions, or `usize::MAX` if lengths differ.
///
/// # Examples
///
/// ```
/// use fgumi_lib::umi::assigner::count_mismatches;
///
/// assert_eq!(count_mismatches("ACGT", "ACGT"), 0);
/// assert_eq!(count_mismatches("ACGT", "ACTT"), 1);
/// assert_eq!(count_mismatches("ACGT", "TGCA"), 4);
/// assert_eq!(count_mismatches("ACG", "ACGT"), usize::MAX);
/// ```
#[must_use]
pub fn count_mismatches(a: &str, b: &str) -> usize {
    if a.len() != b.len() {
        return usize::MAX;
    }
    a.chars().zip(b.chars()).filter(|(x, y)| x != y).count()
}

/// Check if two UMI sequences match within a maximum mismatch threshold.
///
/// This is an optimized version that exits early when the mismatch count exceeds
/// the threshold, avoiding unnecessary comparisons. Use this instead of
/// `count_mismatches(a, b) <= max` when you only need a boolean result.
///
/// # Arguments
///
/// * `a` - First UMI sequence
/// * `b` - Second UMI sequence
/// * `max_mismatches` - Maximum number of allowed mismatches
///
/// # Returns
///
/// `true` if the sequences are the same length and differ by at most `max_mismatches`
/// positions, `false` otherwise.
///
/// # Examples
///
/// ```
/// use fgumi_lib::umi::assigner::matches_within_threshold;
///
/// assert!(matches_within_threshold("ACGT", "ACGT", 1));
/// assert!(matches_within_threshold("ACGT", "ACTT", 1));
/// assert!(!matches_within_threshold("ACGT", "AATT", 1));
/// assert!(!matches_within_threshold("ACG", "ACGT", 1)); // different lengths
/// ```
///
/// # Panics
///
/// Does not panic in practice; internal slice indexing is guaranteed by length checks.
#[must_use]
pub fn matches_within_threshold(a: &str, b: &str, max_mismatches: usize) -> bool {
    if a.len() != b.len() {
        return false;
    }

    let a_bytes = a.as_bytes();
    let b_bytes = b.as_bytes();
    let len = a_bytes.len();
    let too_many = max_mismatches + 1;
    let mut mismatches = 0;

    // Process 8 bytes at a time using u64 XOR (finds differing bytes quickly)
    let chunks = len / 8;
    for i in 0..chunks {
        let offset = i * 8;
        // Safe: we know offset + 8 <= len because i < chunks = len / 8
        let a_chunk = u64::from_ne_bytes(a_bytes[offset..offset + 8].try_into().unwrap());
        let b_chunk = u64::from_ne_bytes(b_bytes[offset..offset + 8].try_into().unwrap());
        let diff = a_chunk ^ b_chunk;
        if diff != 0 {
            // Count differing bytes in this chunk
            // Each differing byte contributes at least 1 to popcount of byte-masked diff
            for j in 0..8 {
                if (diff >> (j * 8)) & 0xFF != 0 {
                    mismatches += 1;
                    if mismatches >= too_many {
                        return false;
                    }
                }
            }
        }
    }

    // Handle remaining bytes
    for i in (chunks * 8)..len {
        if a_bytes[i] != b_bytes[i] {
            mismatches += 1;
            if mismatches >= too_many {
                return false;
            }
        }
    }

    true
}

/// Trait for UMI assignment strategies
///
/// Defines the interface for assigning molecular identifiers to reads based on their UMI
/// sequences. Implementations provide different strategies for grouping UMIs that may
/// differ due to sequencing errors.
///
/// All implementations must be thread-safe (`Send + Sync`) to enable parallel processing.
pub trait UmiAssigner: Send + Sync {
    /// Assign molecule IDs to UMIs
    ///
    /// Takes a slice of raw UMI sequences and returns a Vec where `result[i]` is the
    /// `MoleculeId` assigned to `raw_umis[i]`. UMIs that should be grouped together will
    /// receive the same molecule ID.
    ///
    /// # Arguments
    ///
    /// * `raw_umis` - Slice of UMI sequences to assign. May contain duplicates.
    ///
    /// # Returns
    ///
    /// `Vec<MoleculeId>` where `result[i]` is the assignment for `raw_umis[i]`.
    fn assign(&self, raw_umis: &[Umi]) -> Vec<MoleculeId>;

    /// Check if two UMIs are the same
    ///
    /// Default implementation uses simple string equality. Override for custom
    /// matching logic (e.g., treating A-B equal to B-A in paired UMIs).
    ///
    /// # Arguments
    ///
    /// * `a` - First UMI
    /// * `b` - Second UMI
    ///
    /// # Returns
    ///
    /// `true` if the UMIs should be considered the same, `false` otherwise.
    fn is_same_umi(&self, a: &str, b: &str) -> bool {
        a == b
    }

    /// Get canonical form of UMI
    ///
    /// Returns a canonical representation of the UMI for comparison and deduplication.
    /// Default implementation returns the UMI unchanged. Override to implement custom
    /// canonicalization (e.g., lexical ordering for paired UMIs).
    ///
    /// # Arguments
    ///
    /// * `umi` - UMI sequence to canonicalize
    ///
    /// # Returns
    ///
    /// Canonical form of the UMI
    fn canonicalize(&self, umi: &str) -> String {
        umi.to_string()
    }

    /// Whether to split templates by pair orientation
    ///
    /// Returns `true` if templates should be split based on whether read pairs are in
    /// F1R2 or F2R1 orientation. Returns `false` for paired UMI strategies where both
    /// orientations are intentionally combined.
    ///
    /// # Returns
    ///
    /// `true` if templates should be split by orientation, `false` otherwise.
    fn split_templates_by_pair_orientation(&self) -> bool {
        true
    }

    /// Downcast to concrete type (for pattern matching)
    ///
    /// Returns a reference to the trait object as `Any`, enabling runtime type checking
    /// and downcasting to concrete types.
    ///
    /// # Returns
    ///
    /// Reference to self as `&dyn Any`
    fn as_any(&self) -> &dyn Any;
}

/// Identity UMI assigner - only exact matches
///
/// Groups reads only when their UMI sequences are identical (after converting to uppercase).
/// This is the fastest strategy but provides no error correction. Each unique UMI sequence
/// is assigned a separate molecule ID.
///
/// # Thread Safety
///
/// Uses atomic counter for ID generation, making it safe to use across threads.
pub struct IdentityUmiAssigner {
    counter: AtomicU64,
}

impl IdentityUmiAssigner {
    /// Create a new identity UMI assigner
    ///
    /// # Returns
    ///
    /// A new `IdentityUmiAssigner` instance with counter initialized to 0.
    #[must_use]
    pub fn new() -> Self {
        Self { counter: AtomicU64::new(0) }
    }

    /// Generate the next unique molecule ID
    ///
    /// Atomically increments the counter and returns the previous value.
    ///
    /// # Returns
    ///
    /// A unique molecule ID as u64
    fn next_id(&self) -> u64 {
        self.counter.fetch_add(1, Ordering::SeqCst)
    }
}

impl Default for IdentityUmiAssigner {
    fn default() -> Self {
        Self::new()
    }
}

impl UmiAssigner for IdentityUmiAssigner {
    fn assign(&self, raw_umis: &[Umi]) -> Vec<MoleculeId> {
        if raw_umis.is_empty() {
            return Vec::new();
        }

        // Map each UMI to its canonical (uppercase) form and track unique canonicals
        let mut canonical_to_id: HashMap<String, MoleculeId> = HashMap::new();
        let mut unique_canonicals: Vec<String> = Vec::new();

        for umi in raw_umis {
            let canonical = umi.to_uppercase();
            if !canonical_to_id.contains_key(&canonical) {
                unique_canonicals.push(canonical);
            }
        }

        // Sort for deterministic ID assignment
        unique_canonicals.sort();

        // Assign IDs in sorted order
        for canonical in unique_canonicals {
            canonical_to_id.insert(canonical, MoleculeId::Single(self.next_id()));
        }

        // Build result Vec indexed by input position
        raw_umis.iter().map(|umi| canonical_to_id[&umi.to_uppercase()]).collect()
    }

    fn as_any(&self) -> &dyn Any {
        self
    }
}

/// Simple error-tolerant UMI assigner
///
/// Groups UMIs that are within a specified edit distance using a simple clustering algorithm.
/// When a new UMI is encountered, it is compared against existing groups. If it matches any
/// UMI in a group, it is added to that group. If it matches multiple groups, those groups
/// are merged. This can lead to transitive clustering where UMIs are grouped together even
/// if they differ by more than the threshold.
///
/// # Algorithm
///
/// 1. For each UMI, find all existing groups that contain at least one UMI within threshold
/// 2. If no matches, create a new group
/// 3. If one match, add to that group
/// 4. If multiple matches, merge all matching groups
///
/// # Thread Safety
///
/// Uses atomic counter for ID generation, making it safe to use across threads.
pub struct SimpleErrorUmiAssigner {
    max_mismatches: u32,
    counter: AtomicU64,
}

impl SimpleErrorUmiAssigner {
    /// Create a new simple error-tolerant UMI assigner
    ///
    /// # Arguments
    ///
    /// * `max_mismatches` - Maximum number of mismatches allowed between UMIs in the same group
    ///
    /// # Returns
    ///
    /// A new `SimpleErrorUmiAssigner` instance
    #[must_use]
    pub fn new(max_mismatches: u32) -> Self {
        Self { max_mismatches, counter: AtomicU64::new(0) }
    }

    /// Generate the next unique molecule ID
    ///
    /// Atomically increments the counter and returns the previous value.
    ///
    /// # Returns
    ///
    /// A unique molecule ID as u64
    fn next_id(&self) -> u64 {
        self.counter.fetch_add(1, Ordering::SeqCst)
    }
}

impl UmiAssigner for SimpleErrorUmiAssigner {
    fn assign(&self, raw_umis: &[Umi]) -> Vec<MoleculeId> {
        if raw_umis.is_empty() {
            return Vec::new();
        }

        let mut umi_sets: Vec<BTreeSet<Umi>> = Vec::new();
        // Track seen UMIs to skip duplicate processing
        let mut seen: std::collections::HashSet<&str> =
            std::collections::HashSet::with_capacity(raw_umis.len());

        for umi in raw_umis {
            // Skip duplicate UMIs - they're already in a set
            if !seen.insert(umi.as_str()) {
                continue;
            }

            let mut matched_sets = Vec::new();

            // Find all sets that match this UMI
            for (idx, set) in umi_sets.iter().enumerate() {
                if set
                    .iter()
                    .any(|other| matches_within_threshold(other, umi, self.max_mismatches as usize))
                {
                    matched_sets.push(idx);
                }
            }

            match matched_sets.len() {
                0 => {
                    // No matches - create new set
                    let mut new_set = BTreeSet::new();
                    new_set.insert(umi.clone());
                    umi_sets.push(new_set);
                }
                1 => {
                    // One match - add to that set
                    umi_sets[matched_sets[0]].insert(umi.clone());
                }
                _ => {
                    // Collect UMIs from all matched sets
                    let mut merged = BTreeSet::new();
                    merged.insert(umi.clone());

                    for &idx in matched_sets.iter().rev() {
                        let set = umi_sets.remove(idx);
                        merged.extend(set);
                    }

                    umi_sets.push(merged);
                }
            }
        }

        // Sort UMI sets by their smallest UMI for deterministic ID assignment
        umi_sets.sort_by(|a, b| {
            let a_min = a.iter().next();
            let b_min = b.iter().next();
            a_min.cmp(&b_min)
        });

        // Build map from UMI string to MoleculeId
        let mut umi_to_id: AHashMap<String, MoleculeId> = AHashMap::new();
        for set in umi_sets {
            let id = MoleculeId::Single(self.next_id());
            for umi in set {
                umi_to_id.insert(umi, id);
            }
        }

        // Build result Vec indexed by input position
        raw_umis.iter().map(|umi| umi_to_id[umi]).collect()
    }

    fn as_any(&self) -> &dyn Any {
        self
    }
}

/// Node in the adjacency graph
///
/// Represents a UMI in the directed adjacency graph. Each node tracks its
/// count (number of times observed), child nodes, and assignment status.
/// The UMI string is not stored; use the node index to look up the UMI
/// from the original `umi_counts` slice.
#[derive(Debug)]
struct Node {
    /// Number of times this UMI was observed
    count: usize,
    /// Indices of child nodes in the graph
    children: Vec<usize>,
    /// Whether this node has been assigned to a molecule
    assigned: bool,
}

/// Adjacency UMI assigner using directed graph method
///
/// Implements the directed adjacency graph method from UMI-tools. This method builds a
/// directed graph where high-abundance UMIs can "capture" lower-abundance UMIs that are
/// within edit distance. A key constraint is that a UMI can only be a child of another
/// UMI if its count is less than half of the parent's count (`child_count < parent_count / 2 + 1`).
///
/// This prevents over-correction by ensuring that two similarly abundant UMIs are not
/// incorrectly merged, even if they are within edit distance. The method is more
/// conservative than simple edit distance clustering.
///
/// # Algorithm
///
/// 1. Count UMIs and sort by decreasing abundance
/// 2. Process UMIs in order of abundance
/// 3. For each unassigned UMI (becomes a root):
///    - Search for unassigned UMIs with `count < root_count / 2 + 1`
///    - If within edit distance, add as child and continue recursively
/// 4. Assign same molecule ID to root and all its descendants
///
/// # Multi-threading
///
/// When threads > 1, the matching step is parallelized using rayon for better performance
/// with high UMI counts at the same genomic position.
///
/// # Thread Safety
///
/// Uses atomic counter for ID generation, making it safe to use across threads.
pub struct AdjacencyUmiAssigner {
    max_mismatches: u32,
    threads: usize,
    index_threshold: usize,
    counter: AtomicU64,
}

impl AdjacencyUmiAssigner {
    /// Create a new adjacency UMI assigner
    ///
    /// # Arguments
    ///
    /// * `max_mismatches` - Maximum number of mismatches allowed for UMIs to be adjacent
    /// * `threads` - Number of threads to use for matching (default: 1)
    /// * `index_threshold` - Minimum UMIs per position to use N-gram/BK-tree index (default: 1000)
    ///
    /// # Returns
    ///
    /// A new `AdjacencyUmiAssigner` instance
    #[must_use]
    pub fn new(max_mismatches: u32, threads: usize, index_threshold: usize) -> Self {
        Self { max_mismatches, threads, index_threshold, counter: AtomicU64::new(0) }
    }

    /// Generate the next unique molecule ID
    ///
    /// Atomically increments the counter and returns the previous value.
    ///
    /// # Returns
    ///
    /// A unique molecule ID as u64
    pub(crate) fn next_id(&self) -> u64 {
        self.counter.fetch_add(1, Ordering::SeqCst)
    }

    /// Assign IDs to nodes using `BitEnc`-based `umi_counts`.
    ///
    /// Returns a `Vec<MoleculeId>` where `result[i]` is the `MoleculeId` for `umi_counts[i]`.
    fn assign_ids_to_nodes_bitenc(
        &self,
        nodes: &[Node],
        roots: &[usize],
        umi_counts: &[(BitEnc, usize, usize)], // (encoded, count, first_raw_index)
    ) -> Vec<MoleculeId> {
        // Result Vec indexed by unique UMI index (same order as umi_counts)
        let mut result = vec![MoleculeId::None; umi_counts.len()];

        for &root_idx in roots {
            let id = MoleculeId::Single(self.next_id());
            let root = &nodes[root_idx];

            // Assign ID to root
            result[root_idx] = id;

            // Assign ID to all descendants
            let mut stack = root.children.clone();
            while let Some(child_idx) = stack.pop() {
                let child = &nodes[child_idx];
                result[child_idx] = id;
                stack.extend(&child.children);
            }
        }

        result
    }

    /// Build adjacency graph using `BitEnc` for efficient matching.
    fn build_adjacency_graph_bitenc(
        &self,
        umi_counts: &[(BitEnc, usize, usize)], // (encoded, count, first_raw_index)
        raw_umis: &[Umi],
    ) -> (Vec<Node>, Vec<usize>) {
        // Build nodes
        let mut nodes: Vec<Node> = umi_counts
            .iter()
            .map(|(_, count, _)| Node { count: *count, children: Vec::new(), assigned: false })
            .collect();

        // Build count index lookup for optimization
        let mut count_index: Vec<(usize, usize)> = Vec::new();
        let mut last_count = nodes[0].count;
        count_index.push((last_count, 0));

        for (idx, node) in nodes.iter().enumerate().skip(1) {
            if node.count != last_count {
                count_index.push((node.count, idx));
                last_count = node.count;
            }
        }

        let mut roots = Vec::new();
        let mut queue = VecDeque::new();

        // Create thread pool once for parallel matching (if threads > 1)
        let pool = if self.threads > 1 {
            Some(
                ThreadPoolBuilder::new()
                    .num_threads(self.threads)
                    .build()
                    .expect("Failed to create thread pool"),
            )
        } else {
            None
        };

        // Try to build index for fast candidate search (if enough UMIs)
        let ngram_index = if nodes.len() >= self.index_threshold && self.max_mismatches == 1 {
            // For BitEnc-based matching, we need to use the original strings for NgramIndex
            // since it internally converts to BitEnc
            let umis: Vec<(usize, &str)> = umi_counts
                .iter()
                .enumerate()
                .map(|(i, (_, _, raw_idx))| (i, raw_umis[*raw_idx].as_str()))
                .collect();
            NgramIndex::new(&umis, self.max_mismatches)
        } else {
            None
        };

        // For k>1, we could use BK-tree, but for now use linear scan with BitEnc
        // (BK-tree uses strings internally, so we'd need to refactor it too)

        // Build adjacency graph
        for root_idx in 0..nodes.len() {
            if nodes[root_idx].assigned {
                continue;
            }

            roots.push(root_idx);
            queue.push_back(root_idx);

            while let Some(idx) = queue.pop_front() {
                nodes[idx].assigned = true;

                let max_child_count = nodes[idx].count / 2 + 1;

                // Find where to start searching (for count constraint)
                let search_from = count_index
                    .iter()
                    .find(|(count, _)| *count <= max_child_count)
                    .map_or(nodes.len(), |(_, idx)| *idx);

                if search_from < nodes.len() {
                    let parent_enc = &umi_counts[idx].0;

                    // Find matching children using index or linear scan
                    let matching_indices: Vec<usize> = if let Some(ref index) = ngram_index {
                        // N-gram accelerated search
                        let parent_str = &raw_umis[umi_counts[idx].2];
                        let candidates = index.find_within(parent_str, self.max_mismatches);
                        candidates
                            .into_iter()
                            .filter_map(|(child_idx, _dist)| {
                                if !nodes[child_idx].assigned
                                    && nodes[child_idx].count <= max_child_count
                                {
                                    Some(child_idx)
                                } else {
                                    None
                                }
                            })
                            .collect()
                    } else if let Some(ref pool) = pool {
                        // Parallel linear search using BitEnc hamming distance
                        let start_idx = search_from;
                        let max_mm = self.max_mismatches;
                        pool.install(|| {
                            (start_idx..nodes.len())
                                .into_par_iter()
                                .filter(|&child_idx| {
                                    !nodes[child_idx].assigned
                                        && nodes[child_idx].count <= max_child_count
                                        && parent_enc.hamming_distance(&umi_counts[child_idx].0)
                                            <= max_mm
                                })
                                .collect()
                        })
                    } else {
                        // Sequential linear search using BitEnc hamming distance
                        let start_idx = search_from;
                        (start_idx..nodes.len())
                            .filter(|&child_idx| {
                                !nodes[child_idx].assigned
                                    && nodes[child_idx].count <= max_child_count
                                    && parent_enc.hamming_distance(&umi_counts[child_idx].0)
                                        <= self.max_mismatches
                            })
                            .collect()
                    };

                    // Add matching children
                    for child_idx in matching_indices {
                        nodes[idx].children.push(child_idx);
                        queue.push_back(child_idx);
                        nodes[child_idx].assigned = true;
                    }
                }
            }
        }

        (nodes, roots)
    }

    /// Build adjacency graph using BFS with a custom matcher function
    ///
    /// This is the core graph-building logic shared by both `AdjacencyUmiAssigner`
    /// and `PairedUmiAssigner`. The matcher function determines how UMIs are compared
    /// for adjacency.
    ///
    /// # Arguments
    ///
    /// * `umi_counts` - Pre-sorted vector of (UMI, count) tuples (sorted by descending count)
    /// * `matcher` - Function that returns true if two UMIs should be considered adjacent
    ///
    /// # Returns
    ///
    /// Tuple of (nodes, roots) where nodes is the built graph and roots are the root indices
    #[expect(clippy::too_many_lines, reason = "complex graph-building algorithm with profiling")]
    fn build_adjacency_graph<F>(
        &self,
        umi_counts: &[(Umi, usize)],
        matcher: F,
    ) -> (Vec<Node>, Vec<usize>)
    where
        F: Fn(&str, &str) -> bool + Sync,
    {
        // Build nodes (UMI is not stored; use node index to look up from umi_counts)
        let mut nodes: Vec<Node> = umi_counts
            .iter()
            .map(|(_, count)| Node { count: *count, children: Vec::new(), assigned: false })
            .collect();

        // Build count index lookup for optimization
        let mut count_index: Vec<(usize, usize)> = Vec::new();
        let mut last_count = nodes[0].count;
        count_index.push((last_count, 0));

        for (idx, node) in nodes.iter().enumerate().skip(1) {
            if node.count != last_count {
                count_index.push((node.count, idx));
                last_count = node.count;
            }
        }

        let mut roots = Vec::new();
        let mut queue = VecDeque::new();

        // Create thread pool once for parallel matching (if threads > 1)
        let pool = if self.threads > 1 {
            Some(
                ThreadPoolBuilder::new()
                    .num_threads(self.threads)
                    .build()
                    .expect("Failed to create thread pool"),
            )
        } else {
            None
        };

        // Try to build index for fast candidate search (if enough UMIs)
        // Use N-gram index for k=1 (most efficient), BK-tree for k>1
        let (ngram_index, bktree) = if nodes.len() >= self.index_threshold {
            let umis: Vec<(usize, &str)> =
                umi_counts.iter().enumerate().map(|(i, (umi, _))| (i, umi.as_str())).collect();

            if self.max_mismatches == 1 {
                // N-gram is most efficient for k=1
                (NgramIndex::new(&umis, self.max_mismatches), None)
            } else {
                // BK-tree for k>1
                (None, BkTree::from_umis(&umis))
            }
        } else {
            (None, None)
        };

        // Build adjacency graph
        for root_idx in 0..nodes.len() {
            if nodes[root_idx].assigned {
                continue;
            }

            roots.push(root_idx);
            queue.push_back(root_idx);

            while let Some(idx) = queue.pop_front() {
                nodes[idx].assigned = true;

                let max_child_count = nodes[idx].count / 2 + 1;

                // Find where to start searching (for count constraint)
                let search_from = count_index
                    .iter()
                    .find(|(count, _)| *count <= max_child_count)
                    .map_or(nodes.len(), |(_, idx)| *idx);

                if search_from < nodes.len() {
                    let parent_umi = &umi_counts[idx].0;

                    // Find matching children using index (N-gram or BK-tree) or linear scan
                    let matching_indices: Vec<usize> = if let Some(ref index) = ngram_index {
                        // N-gram accelerated search: candidates already verified by distance
                        // No need to call matcher() again - find_within() already checked Hamming distance
                        let candidates = index.find_within(parent_umi, self.max_mismatches);
                        candidates
                            .into_iter()
                            .filter_map(|(child_idx, _dist)| {
                                if !nodes[child_idx].assigned
                                    && nodes[child_idx].count <= max_child_count
                                {
                                    Some(child_idx)
                                } else {
                                    None
                                }
                            })
                            .collect()
                    } else if let Some(ref tree) = bktree {
                        // BK-tree accelerated search: candidates already verified by distance
                        // No need to call matcher() again - find_within() already checked Hamming distance
                        let candidates = tree.find_within(parent_umi, self.max_mismatches);
                        candidates
                            .into_iter()
                            .filter_map(|(child_idx, _dist)| {
                                if !nodes[child_idx].assigned
                                    && nodes[child_idx].count <= max_child_count
                                {
                                    Some(child_idx)
                                } else {
                                    None
                                }
                            })
                            .collect()
                    } else if let Some(ref pool) = pool {
                        // Parallel linear search using the custom pool
                        let start_idx = search_from;
                        pool.install(|| {
                            (start_idx..nodes.len())
                                .into_par_iter()
                                .filter(|&child_idx| {
                                    // Check cheap conditions first (assigned, count) before expensive matcher
                                    !nodes[child_idx].assigned
                                        && nodes[child_idx].count <= max_child_count
                                        && matcher(parent_umi, &umi_counts[child_idx].0)
                                })
                                .collect()
                        })
                    } else {
                        // Sequential linear search
                        let start_idx = search_from;
                        (start_idx..nodes.len())
                            .filter(|&child_idx| {
                                // Check cheap conditions first (assigned, count) before expensive matcher
                                !nodes[child_idx].assigned
                                    && nodes[child_idx].count <= max_child_count
                                    && matcher(parent_umi, &umi_counts[child_idx].0)
                            })
                            .collect()
                    };

                    // Add matching children
                    for child_idx in matching_indices {
                        nodes[idx].children.push(child_idx);
                        queue.push_back(child_idx);
                        nodes[child_idx].assigned = true;
                    }
                }
            }
        }

        (nodes, roots)
    }
}

impl UmiAssigner for AdjacencyUmiAssigner {
    fn assign(&self, raw_umis: &[Umi]) -> Vec<MoleculeId> {
        if raw_umis.is_empty() {
            return Vec::new();
        }

        // Profiling: track time spent in each phase
        #[cfg(feature = "profile-adjacency")]
        let t_start = std::time::Instant::now();

        // Count UMIs using compact BitEnc representation (saves ~60-70% memory)
        // We store index into raw_umis to avoid cloning strings
        // Also build a map from BitEnc -> unique_index for result lookup
        let mut umi_counts: Vec<(BitEnc, usize, usize)> = {
            // Map from BitEnc -> (count, first_raw_index)
            let mut counts: AHashMap<BitEnc, (usize, usize)> =
                AHashMap::with_capacity(raw_umis.len() / 4); // Estimate ~4 reads per UMI

            for (i, umi) in raw_umis.iter().enumerate() {
                // Convert to compact representation (skips dashes in paired UMIs)
                let Some(encoded) = BitEnc::from_umi_str(umi) else {
                    // UMI contains invalid characters - skip it
                    continue;
                };

                if let Some((count, _idx)) = counts.get_mut(&encoded) {
                    *count += 1;
                } else {
                    // Store index into raw_umis instead of cloning
                    counts.insert(encoded, (1, i));
                }
            }
            // Convert to Vec<(BitEnc, count, first_raw_index)>
            counts.into_iter().map(|(enc, (count, idx))| (enc, count, idx)).collect()
        };

        // Handle case where all UMIs had invalid characters
        if umi_counts.is_empty() {
            return vec![MoleculeId::None; raw_umis.len()];
        }

        #[cfg(feature = "profile-adjacency")]
        let t_after_count = std::time::Instant::now();

        // Sort by descending count, then by original string for determinism
        umi_counts.sort_by(|a, b| b.1.cmp(&a.1).then_with(|| raw_umis[a.2].cmp(&raw_umis[b.2])));

        // Build lookup from BitEnc to sorted index (after sorting)
        let bitenc_to_unique_idx: AHashMap<BitEnc, usize> =
            umi_counts.iter().enumerate().map(|(idx, (enc, _, _))| (*enc, idx)).collect();

        #[cfg(feature = "profile-adjacency")]
        let t_after_sort = std::time::Instant::now();

        if umi_counts.len() == 1 {
            let id = MoleculeId::Single(self.next_id());
            return vec![id; raw_umis.len()];
        }

        // Check all UMIs have same length (in bases, not string length)
        let umi_base_length = umi_counts[0].0.len();
        for (enc, _, _) in &umi_counts {
            assert!(
                enc.len() == umi_base_length,
                "Multiple UMI lengths: {} vs {}",
                enc.len(),
                umi_base_length
            );
        }

        // Build adjacency graph using BitEnc-based matching
        let (nodes, roots) = self.build_adjacency_graph_bitenc(&umi_counts, raw_umis);

        #[cfg(feature = "profile-adjacency")]
        let t_after_graph = std::time::Instant::now();

        // Assign IDs to nodes and build result Vec
        let unique_assignments = self.assign_ids_to_nodes_bitenc(&nodes, &roots, &umi_counts);

        // Map each input UMI to its MoleculeId using BitEnc lookup
        let result: Vec<MoleculeId> = raw_umis
            .iter()
            .map(|umi| {
                if let Some(enc) = BitEnc::from_umi_str(umi) {
                    if let Some(&unique_idx) = bitenc_to_unique_idx.get(&enc) {
                        return unique_assignments[unique_idx];
                    }
                }
                MoleculeId::None // Invalid UMI
            })
            .collect();

        #[cfg(feature = "profile-adjacency")]
        {
            let t_end = std::time::Instant::now();
            let num_unique = umi_counts.len();
            let num_reads = raw_umis.len();

            // Bucket by unique UMI count
            let bucket = match num_unique {
                0..=10 => "1-10",
                11..=50 => "11-50",
                51..=100 => "51-100",
                101..=500 => "101-500",
                501..=1000 => "501-1000",
                _ => "1000+",
            };

            eprintln!(
                "PROFILE bucket={} unique={} reads={} L={} count={:?} sort={:?} graph={:?} assign={:?} total={:?}",
                bucket,
                num_unique,
                num_reads,
                umi_base_length,
                t_after_count.duration_since(t_start),
                t_after_sort.duration_since(t_after_count),
                t_after_graph.duration_since(t_after_sort),
                t_end.duration_since(t_after_graph),
                t_end.duration_since(t_start),
            );
        }

        result
    }

    fn as_any(&self) -> &dyn Any {
        self
    }
}

/// Paired UMI assigner for UMIs of form A-B
///
/// Extends the adjacency method for paired UMIs used in duplex sequencing. Paired UMIs have
/// the form "ACGT-TGCA" where the two parts come from opposite ends of the original molecule.
/// Importantly, A-B and B-A represent reads from opposite strands of the same molecule.
///
/// This assigner:
/// - Treats A-B and B-A as related (counts them together)
/// - Builds adjacency graph on canonical (lexically lower) forms
/// - Assigns strand-specific molecule IDs: `base_id/A` for top strand, `base_id/B` for bottom strand
/// - Determines strand assignment based on similarity to root UMI
///
/// # Algorithm
///
/// 1. Canonicalize paired UMIs (use lexically lower of A-B and B-A)
/// 2. Count canonical UMIs and build adjacency graph (similar to `AdjacencyUmiAssigner`)
/// 3. For each root and its descendants:
///    - Compare to root to determine strand (closer to root or its reverse)
///    - Assign /A suffix if closer to root, /B if closer to reverse
/// 4. Map both A-B and B-A forms to appropriate strand IDs
///
/// # Thread Safety
///
/// Uses atomic counter for ID generation, making it safe to use across threads.
pub struct PairedUmiAssigner {
    adjacency: AdjacencyUmiAssigner,
    lower_prefix: String,
    higher_prefix: String,
}

impl PairedUmiAssigner {
    /// Create a new paired UMI assigner
    ///
    /// # Arguments
    ///
    /// * `max_mismatches` - Maximum number of mismatches allowed for UMIs to be adjacent
    ///
    /// # Returns
    ///
    /// A new `PairedUmiAssigner` instance
    #[must_use]
    pub fn new(max_mismatches: u32) -> Self {
        Self::new_with_threads(max_mismatches, 1, DEFAULT_INDEX_THRESHOLD)
    }

    /// Create a new paired UMI assigner with specified thread count
    ///
    /// # Arguments
    ///
    /// * `max_mismatches` - Maximum number of mismatches allowed for UMIs to be adjacent
    /// * `threads` - Number of threads to use for matching
    /// * `index_threshold` - Minimum UMIs per position to use N-gram/BK-tree index
    ///
    /// # Returns
    ///
    /// A new `PairedUmiAssigner` instance
    #[must_use]
    pub fn new_with_threads(max_mismatches: u32, threads: usize, index_threshold: usize) -> Self {
        let prefix_len = max_mismatches as usize + 1;
        Self {
            adjacency: AdjacencyUmiAssigner::new(max_mismatches, threads, index_threshold),
            lower_prefix: "a".repeat(prefix_len),
            higher_prefix: "b".repeat(prefix_len),
        }
    }

    /// Get the prefix for lower-ordered reads in a pair
    ///
    /// # Returns
    ///
    /// String slice containing the lower prefix
    pub fn lower_read_umi_prefix(&self) -> &str {
        &self.lower_prefix
    }

    /// Get the prefix for higher-ordered reads in a pair
    ///
    /// # Returns
    ///
    /// String slice containing the higher prefix
    pub fn higher_read_umi_prefix(&self) -> &str {
        &self.higher_prefix
    }

    /// Split paired UMI into two parts
    ///
    /// Splits a paired UMI like "ACGT-TGCA" into its two components.
    ///
    /// # Arguments
    ///
    /// * `umi` - Paired UMI string (must contain exactly one '-')
    ///
    /// # Returns
    ///
    /// Result containing tuple of (first part, second part) as string slices
    ///
    /// # Errors
    ///
    /// Returns error if UMI doesn't contain exactly one '-'
    fn split(umi: &str) -> Result<(&str, &str)> {
        // Use byte operations for faster lookup since UMIs are ASCII
        let bytes = umi.as_bytes();
        if let Some(pos) = bytes.iter().position(|&b| b == b'-') {
            // Check for additional dashes after first
            if bytes[pos + 1..].contains(&b'-') {
                return Err(anyhow::anyhow!(
                    "UMI '{umi}' is not a valid paired UMI (expected format: 'A-B', found multiple '-')"
                ));
            }
            // SAFETY: '-' is a single byte ASCII, so splitting at its position is valid UTF-8
            Ok((&umi[..pos], &umi[pos + 1..]))
        } else {
            Err(anyhow::anyhow!("UMI '{umi}' is not a valid paired UMI (expected format: 'A-B')"))
        }
    }

    /// Reverse a paired UMI (A-B -> B-A)
    ///
    /// Swaps the two parts of a paired UMI, representing the read from the opposite strand.
    ///
    /// # Arguments
    ///
    /// * `umi` - Paired UMI to reverse
    ///
    /// # Returns
    ///
    /// Result containing the reversed UMI
    ///
    /// # Errors
    ///
    /// Returns error if UMI is not a valid paired UMI
    fn reverse(umi: &str) -> Result<String> {
        let (first, second) = Self::split(umi)?;
        // Pre-allocate exact capacity: second.len() + 1 (for '-') + first.len()
        let mut result = String::with_capacity(umi.len());
        result.push_str(second);
        result.push('-');
        result.push_str(first);
        Ok(result)
    }

    /// Get canonical form (lexically lower of A-B and B-A)
    ///
    /// Returns the lexically lower of the two orientations. This ensures that A-B and B-A
    /// map to the same canonical form for counting and grouping.
    ///
    /// # Arguments
    ///
    /// * `umi` - Paired UMI to canonicalize
    ///
    /// # Returns
    ///
    /// Result containing the canonical (lexically lower) form
    ///
    /// # Errors
    ///
    /// Returns error if UMI is not a valid paired UMI
    fn canonicalize_paired(umi: &str) -> Result<String> {
        let (first, second) = Self::split(umi)?;
        // Compare parts directly to avoid allocating reversed string when not needed
        // If first <= second, original "first-second" is canonical
        // If first > second, reversed "second-first" is canonical
        if first <= second {
            Ok(umi.to_string())
        } else {
            // Only allocate when we need the reversed form
            let mut result = String::with_capacity(umi.len());
            result.push_str(second);
            result.push('-');
            result.push_str(first);
            Ok(result)
        }
    }

    /// Count UMIs, treating A-B and B-A as the same
    ///
    /// Counts canonical forms of paired UMIs, so A-B and B-A contribute to the same count.
    ///
    /// # Arguments
    ///
    /// * `umis` - Slice of paired UMI sequences
    ///
    /// # Returns
    ///
    /// Vector of (canonical UMI, count) tuples
    fn count_paired(umis: &[Umi]) -> Vec<(Umi, usize)> {
        let mut counts: AHashMap<Umi, usize> = AHashMap::new();
        for umi in umis {
            let canonical = Self::canonicalize_paired(umi)
                .expect("UMI should be valid paired format (validated in assign())");
            *counts.entry(canonical).or_insert(0) += 1;
        }
        counts.into_iter().collect()
    }

    /// Check if UMIs match (including reversed)
    ///
    /// Returns true if the two UMIs are within edit distance, considering both the
    /// forward and reverse orientations. For example, "ACGT-TGCA" matches "TGCA-ACGT"
    /// with 0 edits even though they are different strings.
    ///
    /// # Arguments
    ///
    /// * `lhs` - First UMI
    /// * `rhs` - Second UMI
    ///
    /// # Returns
    ///
    /// `true` if UMIs match (forward or reversed) within threshold
    fn matches_paired(&self, lhs: &str, rhs: &str) -> bool {
        let max_mismatches = self.adjacency.max_mismatches as usize;
        if matches_within_threshold(lhs, rhs, max_mismatches) {
            return true;
        }
        if let Ok(lhs_rev) = Self::reverse(lhs) {
            matches_within_threshold(&lhs_rev, rhs, max_mismatches)
        } else {
            false
        }
    }
}

impl UmiAssigner for PairedUmiAssigner {
    fn assign(&self, raw_umis: &[Umi]) -> Vec<MoleculeId> {
        if raw_umis.is_empty() {
            return Vec::new();
        }

        // Validate all are paired UMIs
        for umi in raw_umis {
            assert!((umi.split('-').count() == 2), "UMI {umi} is not a paired UMI");
        }

        // Count with A-B and B-A combined
        let mut umi_counts = Self::count_paired(raw_umis);
        umi_counts.sort_by(|a, b| b.1.cmp(&a.1).then_with(|| a.0.cmp(&b.0)));

        if umi_counts.len() == 1 {
            let id = self.adjacency.next_id();
            let ab = MoleculeId::PairedA(id);
            let ba = MoleculeId::PairedB(id);

            return raw_umis
                .iter()
                .map(|umi| {
                    let reversed = Self::reverse(umi)
                        .expect("UMI should be valid paired format (validated above)");
                    if *umi < reversed { ab } else { ba }
                })
                .collect();
        }

        // Build adjacency graph using the shared helper with our paired matcher
        let (nodes, roots) = self
            .adjacency
            .build_adjacency_graph(&umi_counts, |lhs, rhs| self.matches_paired(lhs, rhs));

        // Assign IDs with A/B variants - build a map from UMI string to MoleculeId
        // Use owned strings since we need to store reversed UMIs
        let mut umi_to_id: AHashMap<String, MoleculeId> = AHashMap::new();

        for &root_idx in &roots {
            let id = self.adjacency.next_id();
            let ab = MoleculeId::PairedA(id);
            let ba = MoleculeId::PairedB(id);

            let root = &nodes[root_idx];
            let root_umi = &umi_counts[root_idx].0;
            let root_rev = Self::reverse(root_umi)
                .expect("UMI should be valid paired format (validated above)");

            umi_to_id.insert(root_umi.clone(), ab);
            umi_to_id.insert(root_rev, ba);

            // Process descendants
            let mut stack = root.children.clone();
            while let Some(child_idx) = stack.pop() {
                let child = &nodes[child_idx];
                let child_umi = &umi_counts[child_idx].0;
                let child_rev = Self::reverse(child_umi)
                    .expect("UMI should be valid paired format (validated above)");

                // Determine which strand based on similarity to root
                if count_mismatches(root_umi, child_umi) < count_mismatches(root_umi, &child_rev) {
                    umi_to_id.insert(child_umi.clone(), ab);
                    umi_to_id.insert(child_rev, ba);
                } else {
                    umi_to_id.insert(child_umi.clone(), ba);
                    umi_to_id.insert(child_rev, ab);
                }

                stack.extend(&child.children);
            }
        }

        // Build result Vec indexed by input position
        raw_umis.iter().map(|umi| umi_to_id.get(umi).copied().unwrap_or(MoleculeId::None)).collect()
    }

    fn is_same_umi(&self, a: &str, b: &str) -> bool {
        if a == b {
            return true;
        }
        // Check if A-B equals B-A
        if let (Ok((a1, a2)), Ok((b1, b2))) = (Self::split(a), Self::split(b)) {
            a1 == b2 && a2 == b1
        } else {
            false
        }
    }

    fn canonicalize(&self, umi: &str) -> String {
        Self::canonicalize_paired(umi).unwrap_or_else(|_| umi.to_string())
    }

    fn split_templates_by_pair_orientation(&self) -> bool {
        false
    }

    fn as_any(&self) -> &dyn Any {
        self
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use rstest::rstest;
    use std::collections::{HashMap, HashSet};

    // ========================================================================
    // Helper Functions
    // ========================================================================

    /// Helper function to repeat an element n times
    fn repeat_n<T: Clone>(item: T, n: usize) -> Vec<T> {
        vec![item; n]
    }

    /// Helper to convert assignments to groups for easier testing
    /// Returns a vector of sets of UMIs that are assigned the same ID
    fn group_assignments(
        umis: &[Umi],
        assignments: &[MoleculeId],
        strip_suffix: bool,
    ) -> Vec<HashSet<String>> {
        let mut groups: HashMap<String, HashSet<String>> = HashMap::new();

        for (umi, id) in umis.iter().zip(assignments.iter()) {
            let group_key = if strip_suffix {
                // Strip /A or /B suffix
                id.base_id_string()
            } else {
                id.to_string()
            };
            groups.entry(group_key).or_default().insert(umi.clone());
        }

        groups.into_values().collect()
    }

    /// Helper to build a lookup map from UMI string to `MoleculeId` for test assertions
    fn build_assignment_map(
        umis: &[Umi],
        assignments: &[MoleculeId],
    ) -> HashMap<String, MoleculeId> {
        umis.iter().zip(assignments.iter()).map(|(umi, id)| (umi.clone(), *id)).collect()
    }

    /// Helper to check if two group collections are equivalent
    fn assert_groups_eq(groups: &[HashSet<String>], expected: &[HashSet<String>]) {
        assert_eq!(groups.len(), expected.len(), "Different number of groups");
        for exp_group in expected {
            assert!(
                groups.iter().any(|g| g == exp_group),
                "Expected group {exp_group:?} not found in actual groups"
            );
        }
    }

    // ========================================================================
    // Basic Function Tests
    // ========================================================================

    #[rstest]
    #[case("AAAA", "AAAA", 0, "identical sequences")]
    #[case("AAAA", "AAAC", 1, "one mismatch at end")]
    #[case("AAAA", "CCCC", 4, "all positions differ")]
    #[case("AAAA", "AACC", 2, "two mismatches")]
    #[case("AAA", "AAAA", usize::MAX, "different lengths")]
    #[case("", "", 0, "empty sequences")]
    #[case("A", "A", 0, "single character match")]
    #[case("A", "T", 1, "single character mismatch")]
    fn test_count_mismatches(
        #[case] a: &str,
        #[case] b: &str,
        #[case] expected: usize,
        #[case] description: &str,
    ) {
        assert_eq!(count_mismatches(a, b), expected, "Failed for: {description}");
    }

    // ========================================================================
    // IdentityUmiAssigner Tests
    // ========================================================================

    #[test]
    fn test_identity_assigner_basic() {
        let assigner = IdentityUmiAssigner::new();
        let umis = vec![
            "AAAAAA".to_string(),
            "AAAAAT".to_string(),
            "ACGTAC".to_string(),
            "ACGTAC".to_string(),
            "AAAAAA".to_string(),
        ];

        let assignments = assigner.assign(&umis);

        // Count unique IDs
        let unique_ids: HashSet<_> = assignments.iter().collect();
        assert_eq!(unique_ids.len(), 3);
    }

    #[test]
    fn test_identity_assigner_groups_exact_matches() {
        let assigner = IdentityUmiAssigner::new();
        let umis = vec![
            "AAAAAA".to_string(),
            "AAAAAT".to_string(),
            "ACGTAC".to_string(),
            "ACGTAC".to_string(),
            "AAAAAA".to_string(),
        ];

        let assignments = assigner.assign(&umis);
        let groups = group_assignments(&umis, &assignments, false);

        // Should create 3 groups: {AAAAAA}, {AAAAAT}, {ACGTAC}
        assert_eq!(groups.len(), 3);

        let expected_groups: Vec<HashSet<String>> = vec![
            vec!["AAAAAA".to_string()].into_iter().collect(),
            vec!["AAAAAT".to_string()].into_iter().collect(),
            vec!["ACGTAC".to_string()].into_iter().collect(),
        ];

        assert_groups_eq(&groups, &expected_groups);
    }

    // ========================================================================
    // SimpleErrorUmiAssigner Tests
    // ========================================================================

    #[test]
    fn test_simple_error_assigner_basic() {
        let assigner = SimpleErrorUmiAssigner::new(1);
        let umis = vec![
            "AAAAAA".to_string(),
            "AAAATT".to_string(),
            "AAAATA".to_string(),
            "TGCACC".to_string(),
            "TGCACG".to_string(),
            "GGCGGC".to_string(),
        ];

        let assignments = assigner.assign(&umis);

        // Should create 3 groups
        let unique_ids: HashSet<_> = assignments.iter().collect();
        assert_eq!(unique_ids.len(), 3);
    }

    #[test]
    #[expect(clippy::similar_names, reason = "umi variable names intentionally similar")]
    fn test_simple_error_assigner_groups_by_mismatches() {
        let assigner = SimpleErrorUmiAssigner::new(1);
        let umis = vec![
            "AAAAAA".to_string(),
            "AAAATT".to_string(),
            "AAAATA".to_string(),
            "TGCACC".to_string(),
            "TGCACG".to_string(),
            "GGCGGC".to_string(),
            "GGCGGC".to_string(),
            "GGCGGC".to_string(),
        ];

        let assignments = assigner.assign(&umis);
        let groups = group_assignments(&umis, &assignments, false);

        // Should create 3 groups
        assert_eq!(groups.len(), 3);

        // Verify the specific groupings
        let group1: HashSet<String> =
            vec!["AAAAAA".to_string(), "AAAATT".to_string(), "AAAATA".to_string()]
                .into_iter()
                .collect();
        let group2: HashSet<String> = vec!["GGCGGC".to_string()].into_iter().collect();
        let group3: HashSet<String> =
            vec!["TGCACC".to_string(), "TGCACG".to_string()].into_iter().collect();

        assert!(groups.iter().any(|g| g == &group1));
        assert!(groups.iter().any(|g| g == &group2));
        assert!(groups.iter().any(|g| g == &group3));
    }

    #[test]
    fn test_simple_error_assigner_groups_everything_with_high_edits() {
        let assigner = SimpleErrorUmiAssigner::new(6);
        let umis = vec![
            "AAAAAA".to_string(),
            "AAAATT".to_string(),
            "AAAATA".to_string(),
            "TGCACC".to_string(),
            "TGCACG".to_string(),
            "GGCGGC".to_string(),
            "GGCGGC".to_string(),
            "GGCGGC".to_string(),
        ];

        let assignments = assigner.assign(&umis);
        let groups = group_assignments(&umis, &assignments, false);

        // With edits=6, everything should be in one group
        assert_eq!(groups.len(), 1);
    }

    #[test]
    fn test_simple_error_with_zero_edits() {
        let assigner = SimpleErrorUmiAssigner::new(0);
        let umis = vec!["AAAAAA".to_string(), "AAAAAA".to_string(), "AAAAAT".to_string()];
        let assignments = assigner.assign(&umis);
        let groups = group_assignments(&umis, &assignments, false);

        // With 0 edits, should behave like identity
        assert_eq!(groups.len(), 2);
    }

    #[test]
    fn test_adjacency_with_zero_edits() {
        let assigner = AdjacencyUmiAssigner::new(0, 1, DEFAULT_INDEX_THRESHOLD);
        let umis = vec!["AAAAAA".to_string(), "AAAAAA".to_string(), "AAAAAT".to_string()];
        let assignments = assigner.assign(&umis);
        let groups = group_assignments(&umis, &assignments, false);

        // With 0 edits, only exact matches group together
        assert_eq!(groups.len(), 2);
    }

    /// Test for the backward capture bug fix.
    ///
    /// This tests a scenario where during BFS traversal, a child node should be able
    /// to capture an unassigned node at an earlier index. The bug was that
    /// `start_idx = search_from.max(idx + 1)` prevented searching earlier indices.
    ///
    /// Scenario with edits=1:
    /// - Node 0: count=100, UMI="AAAA"
    /// - Node 1: count=2,   UMI="CCAA" (2 edits from AAAA, 1 edit from ACAA)
    /// - Node 2: count=2,   UMI="GGGG" (4 edits from everything)
    /// - Node 3: count=2,   UMI="ACAA" (1 edit from AAAA, 1 edit from CCAA)
    ///
    /// Expected behavior:
    /// - AAAA captures ACAA (1 edit)
    /// - ACAA captures CCAA (1 edit) via backward search
    /// - Result: AAAA, ACAA, CCAA in same group; GGGG separate
    #[test]
    fn test_adjacency_backward_capture() {
        let assigner = AdjacencyUmiAssigner::new(1, 1, DEFAULT_INDEX_THRESHOLD);

        // Construct input with counts: AAAA(100), CCAA(1), GGGG(2), ACAA(2)
        // Note: CCAA count=1 so it satisfies count constraint: 1 < ACAA_count/2+1 = 2
        let mut umis = Vec::new();
        umis.extend(std::iter::repeat_n("AAAA".to_string(), 100));
        umis.extend(std::iter::repeat_n("CCAA".to_string(), 1));
        umis.extend(std::iter::repeat_n("GGGG".to_string(), 2));
        umis.extend(std::iter::repeat_n("ACAA".to_string(), 2));

        let assignments = assigner.assign(&umis);
        let map = build_assignment_map(&umis, &assignments);

        // AAAA, CCAA, ACAA should all have same ID
        // AAAA captures ACAA (1 edit), ACAA captures CCAA (1 edit via backward search)
        assert_eq!(
            map.get("AAAA"),
            map.get("CCAA"),
            "CCAA should be captured via ACAA intermediate"
        );
        assert_eq!(map.get("AAAA"), map.get("ACAA"), "ACAA should be captured by AAAA");

        // GGGG should be separate (4 edits from everything)
        assert_ne!(map.get("AAAA"), map.get("GGGG"), "GGGG should be separate");
    }

    // ========================================================================
    // PairedUmiAssigner Tests
    // ========================================================================

    #[test]
    fn test_paired_assigner_basic() {
        let assigner = PairedUmiAssigner::new(1);
        let umis = vec!["AAAA-CCCC".to_string(), "CCCC-AAAA".to_string()];

        let assignments = assigner.assign(&umis);

        // Should have different suffixes but same base ID
        let ids: Vec<_> = assignments.iter().collect();
        assert_eq!(ids.len(), 2);

        // Extract base IDs (before the slash)
        let base_ids: HashSet<_> = ids.iter().map(|id| id.base_id_string()).collect();
        assert_eq!(base_ids.len(), 1);
    }

    #[test]
    fn test_paired_assigns_ab_and_ba_with_different_suffix() {
        let umis = vec!["AAAA-CCCC".to_string(), "CCCC-AAAA".to_string()];
        let assigner = PairedUmiAssigner::new(1);
        let assignments = assigner.assign(&umis);

        // With strip_suffix=false, should be different groups
        let groups_no_strip = group_assignments(&umis, &assignments, false);
        assert_eq!(groups_no_strip.len(), 2);

        // With strip_suffix=true, should be same base group
        let groups_strip = group_assignments(&umis, &assignments, true);
        assert_eq!(groups_strip.len(), 1);
    }

    #[test]
    fn test_paired_assigns_ab_and_ba_with_errors() {
        // Original UMIs have count=3, error variants have count=1
        // Count constraint: child < parent/2+1 means 1 < 3/2+1=2 ✓
        let mut umis = Vec::new();
        umis.extend(std::iter::repeat_n("AAAA-CCCC".to_string(), 3));
        umis.extend(std::iter::repeat_n("CCCC-CAAA".to_string(), 1));
        umis.extend(std::iter::repeat_n("AAAA-GGGG".to_string(), 3));
        umis.extend(std::iter::repeat_n("GGGG-AAGA".to_string(), 1));

        let assigner = PairedUmiAssigner::new(1);
        let assignments = assigner.assign(&umis);

        // With strip_suffix=false, should be 4 groups (each UMI string is distinct)
        let groups_no_strip = group_assignments(&umis, &assignments, false);
        assert_eq!(groups_no_strip.len(), 4);

        // With strip_suffix=true, should be 2 base groups
        // Group 1: AAAA-CCCC and CCCC-CAAA (canonically CAAA-CCCC, 1 edit from AAAA-CCCC)
        // Group 2: AAAA-GGGG and GGGG-AAGA (canonically AAGA-GGGG, 1 edit from AAAA-GGGG)
        let groups_strip = group_assignments(&umis, &assignments, true);
        assert_eq!(groups_strip.len(), 2);
    }

    #[test]
    fn test_paired_handles_errors_in_first_base_changing_lexical_ordering() {
        let mut umis = Vec::new();
        umis.extend(repeat_n("GTGT-ACAC".to_string(), 500));
        umis.extend(repeat_n("ACAC-GTGT".to_string(), 460));
        umis.extend(repeat_n("GTGT-TCAC".to_string(), 6));
        umis.extend(repeat_n("TCAC-GTGT".to_string(), 6));
        umis.extend(repeat_n("GTGT-TGAC".to_string(), 1));

        let assigner = PairedUmiAssigner::new(1);
        let assignments = assigner.assign(&umis);

        // With strip_suffix=false
        let groups_no_strip = group_assignments(&umis, &assignments, false);
        assert_eq!(groups_no_strip.len(), 2);

        // Verify specific groupings
        let group1: HashSet<String> =
            vec!["GTGT-ACAC".to_string(), "GTGT-TCAC".to_string(), "GTGT-TGAC".to_string()]
                .into_iter()
                .collect();
        let group2: HashSet<String> =
            vec!["TCAC-GTGT".to_string(), "ACAC-GTGT".to_string()].into_iter().collect();

        assert!(groups_no_strip.contains(&group1));
        assert!(groups_no_strip.contains(&group2));

        // With strip_suffix=true, should all be in one base group
        let groups_strip = group_assignments(&umis, &assignments, true);
        assert_eq!(groups_strip.len(), 1);
    }

    #[test]
    fn test_paired_counts_ab_and_ba_together_when_constructing_graph() {
        // Since the graph only creates nodes where count(child) <= count(parent) / 2 + 1,
        // it should group everything together in the first set, but not in the second set.
        let mut umis1 = Vec::new();
        umis1.extend(repeat_n("AAAA-CCCC".to_string(), 256));
        umis1.extend(repeat_n("AAAA-CCCG".to_string(), 64));
        umis1.extend(repeat_n("CCCG-AAAA".to_string(), 64));

        let mut umis2 = Vec::new();
        umis2.extend(repeat_n("AAAA-CCCC".to_string(), 256));
        umis2.extend(repeat_n("AAAA-CCCG".to_string(), 128));
        umis2.extend(repeat_n("CCCG-AAAA".to_string(), 128));

        let assigner = PairedUmiAssigner::new(1);
        let assignments1 = assigner.assign(&umis1);
        let assignments2 = assigner.assign(&umis2);

        let groups1_strip = group_assignments(&umis1, &assignments1, true);
        let groups2_strip = group_assignments(&umis2, &assignments2, true);

        // umis1: combined count of AAAA-CCCG/CCCG-AAAA is 128, which is <= 256/2 + 1, so should merge
        assert_eq!(groups1_strip.len(), 1);

        // umis2: combined count of AAAA-CCCG/CCCG-AAAA is 256, which is > 256/2 + 1, so should NOT merge
        assert_eq!(groups2_strip.len(), 2);
    }

    #[test]
    #[should_panic(expected = "is not a paired UMI")]
    fn test_paired_fails_if_supplied_non_paired_umis() {
        let umis = vec!["AAAAAAAA".to_string(), "GGGGGGGG".to_string()];
        let assigner = PairedUmiAssigner::new(1);
        let _ = assigner.assign(&umis);
    }

    // ========================================================================
    // PairedUmiAssigner Helper Function Tests
    // ========================================================================
    #[rstest]
    #[case("ACGT-TGCA", true, Some(("ACGT", "TGCA")), "normal paired UMI")]
    #[case("ACT-", true, Some(("ACT", "")), "empty second part")]
    #[case("-ACT", true, Some(("", "ACT")), "empty first part")]
    #[case("ACTACT", false, None, "no delimiter")]
    #[case("A-B-C", false, None, "multiple delimiters")]
    fn test_paired_split(
        #[case] umi: &str,
        #[case] should_succeed: bool,
        #[case] expected: Option<(&str, &str)>,
        #[case] description: &str,
    ) {
        let result = PairedUmiAssigner::split(umi);
        if should_succeed {
            assert!(result.is_ok(), "Failed for: {description}");
            assert_eq!(result.unwrap(), expected.unwrap(), "Failed for: {description}");
        } else {
            assert!(result.is_err(), "Should have failed for: {description}");
        }
    }

    #[rstest]
    #[case("ACGT-TGCA", "TGCA-ACGT", "normal paired UMI")]
    #[case("ACT-", "-ACT", "empty second part")]
    #[case("-ACT", "ACT-", "empty first part")]
    #[case("AAAA-TTTT", "TTTT-AAAA", "symmetric UMI")]
    fn test_paired_reverse(#[case] input: &str, #[case] expected: &str, #[case] description: &str) {
        let result = PairedUmiAssigner::reverse(input);
        assert!(result.is_ok(), "Failed for: {description}");
        assert_eq!(result.unwrap(), expected, "Failed for: {description}");
    }

    #[rstest]
    #[case("ACGT-TGCA", "ACGT-TGCA", "already canonical")]
    #[case("TGCA-ACGT", "ACGT-TGCA", "needs reversal")]
    #[case("ACT-", "-ACT", "empty second canonicalizes to empty first")]
    #[case("-ACT", "-ACT", "empty first is canonical")]
    #[case("AAAA-TTTT", "AAAA-TTTT", "symmetric but first is lexically lower")]
    fn test_paired_canonicalize(
        #[case] input: &str,
        #[case] expected: &str,
        #[case] description: &str,
    ) {
        let result = PairedUmiAssigner::canonicalize_paired(input);
        assert!(result.is_ok(), "Failed for: {description}");
        assert_eq!(result.unwrap(), expected, "Failed for: {description}");
    }

    #[test]
    fn test_paired_count() {
        let umis = vec![
            "ACGT-TGCA".to_string(),
            "TGCA-ACGT".to_string(),
            "AAAA-TTTT".to_string(),
            "ACGT-TGCA".to_string(),
        ];

        let counts = PairedUmiAssigner::count_paired(&umis);

        // ACGT-TGCA and TGCA-ACGT should be counted together (canonicalized to ACGT-TGCA)
        // So we should have ACGT-TGCA: 3, AAAA-TTTT: 1
        assert_eq!(counts.len(), 2);

        let mut count_map: HashMap<String, usize> = HashMap::new();
        for (umi, count) in counts {
            count_map.insert(umi, count);
        }

        assert_eq!(count_map.get("ACGT-TGCA"), Some(&3));
        assert_eq!(count_map.get("AAAA-TTTT"), Some(&1));
    }

    #[test]
    fn test_paired_is_same_umi() {
        let assigner = PairedUmiAssigner::new(1);

        // Same UMI
        assert!(assigner.is_same_umi("ACGT-TGCA", "ACGT-TGCA"));

        // Reversed UMI
        assert!(assigner.is_same_umi("ACGT-TGCA", "TGCA-ACGT"));

        // Different UMI
        assert!(!assigner.is_same_umi("ACGT-TGCA", "AAAA-TTTT"));
    }

    #[test]
    fn test_paired_matches() {
        let assigner = PairedUmiAssigner::new(1);

        // Exact match
        assert!(assigner.matches_paired("ACGT-TGCA", "ACGT-TGCA"));

        // Reversed match
        assert!(assigner.matches_paired("ACGT-TGCA", "TGCA-ACGT"));

        // Within 1 mismatch
        assert!(assigner.matches_paired("ACGT-TGCA", "ACTT-TGCA"));

        // Within 1 mismatch (reversed)
        assert!(assigner.matches_paired("ACGT-TGCA", "TGCA-ACTT"));

        // Too many mismatches
        assert!(!assigner.matches_paired("ACGT-TGCA", "AAAA-TTTT"));
    }

    /// Test for the backward capture bug fix in `PairedUmiAssigner`.
    ///
    /// This tests the same scenario as `test_adjacency_backward_capture` but with
    /// paired UMIs. During BFS traversal, a child node should be able to capture
    /// an unassigned node at an earlier index.
    ///
    /// Scenario with edits=1:
    /// - Root: "AAAA-XXXX" (count=100)
    /// - Target: "CCAA-XXXX" (count=1) - 2 edits from root, won't be captured by root
    /// - Intermediate: "ACAA-XXXX" (count=2) - 1 edit from root, 1 edit from target
    ///
    /// Note: Target count=1 satisfies count constraint: `1 < intermediate_count/2+1 = 2`
    /// Expected: Root captures intermediate, intermediate captures target via backward search
    #[test]
    fn test_paired_backward_capture() {
        let assigner = PairedUmiAssigner::new(1);

        // Construct paired UMIs with the triangle inequality scenario
        let mut umis = Vec::new();

        // High count root: "AAAA-XXXX"
        umis.extend(std::iter::repeat_n("AAAA-XXXX".to_string(), 100));

        // Low count target: "CCAA-XXXX" (count=1) - 2 edits from root, won't be captured by root
        umis.extend(std::iter::repeat_n("CCAA-XXXX".to_string(), 1));

        // Low count intermediate: "ACAA-XXXX" - 1 edit from root AND 1 edit from target
        umis.extend(std::iter::repeat_n("ACAA-XXXX".to_string(), 2));

        let assignments = assigner.assign(&umis);
        let map = build_assignment_map(&umis, &assignments);

        // Helper to extract base ID (before /A or /B suffix)
        let get_base_id = |umi: &str| -> String { map.get(umi).unwrap().base_id_string() };

        // All three should have the same base ID
        // Root captures intermediate (1 edit), intermediate captures target (1 edit via backward search)
        assert_eq!(
            get_base_id("AAAA-XXXX"),
            get_base_id("CCAA-XXXX"),
            "CCAA-XXXX should be captured via ACAA-XXXX intermediate"
        );
        assert_eq!(
            get_base_id("AAAA-XXXX"),
            get_base_id("ACAA-XXXX"),
            "ACAA-XXXX should be captured by AAAA-XXXX"
        );
    }

    // ========================================================================
    // Edge Case Tests
    // ========================================================================

    #[test]
    fn test_empty_umi_list() {
        let assigner = IdentityUmiAssigner::new();
        let umis: Vec<String> = vec![];
        let assignments = assigner.assign(&umis);
        assert_eq!(assignments.len(), 0);
    }

    #[test]
    fn test_single_umi() {
        let assigner = IdentityUmiAssigner::new();
        let umis = vec!["ACGTACGT".to_string()];
        let assignments = assigner.assign(&umis);
        assert_eq!(assignments.len(), 1);
        // Every UMI gets an assignment (by index), verify it's not None
        assert_ne!(assignments[0], MoleculeId::None);
    }

    #[test]
    fn test_all_identical_umis() {
        let assigner = IdentityUmiAssigner::new();
        let umis = repeat_n("ACGTACGT".to_string(), 100);
        let assignments = assigner.assign(&umis);

        let unique_ids: HashSet<_> = assignments.iter().collect();
        assert_eq!(unique_ids.len(), 1);
    }

    #[test]
    fn test_umis_with_empty_parts_in_paired() {
        let umis = vec!["ACT-".to_string(), "-ACT".to_string(), "ACT-".to_string()];
        let assigner = PairedUmiAssigner::new(1);
        let assignments = assigner.assign(&umis);

        // ACT- and -ACT are reverses of each other
        let groups_strip = group_assignments(&umis, &assignments, true);
        assert_eq!(groups_strip.len(), 1);
    }

    // ========================================================================
    // Property-Based Tests
    // ========================================================================

    use proptest::prelude::*;
    use proptest::strategy::Strategy as PropStrategy;

    /// Generate a random UMI sequence of length 4-12
    fn umi_strategy() -> impl PropStrategy<Value = String> {
        "[ACGT]{4,12}".prop_map(|s| s.clone())
    }

    /// Generate a paired UMI (A-B format)
    fn paired_umi_strategy() -> impl PropStrategy<Value = String> {
        ("[ACGT]{4,8}", "[ACGT]{4,8}").prop_map(|(a, b)| format!("{a}-{b}"))
    }

    proptest! {
        /// Property: count_mismatches is symmetric
        #[test]
        fn prop_count_mismatches_symmetric(a in umi_strategy(), b in umi_strategy()) {
            assert_eq!(count_mismatches(&a, &b), count_mismatches(&b, &a));
        }

        /// Property: count_mismatches with self is always 0
        #[test]
        fn prop_count_mismatches_self_is_zero(umi in umi_strategy()) {
            assert_eq!(count_mismatches(&umi, &umi), 0);
        }

        /// Property: Identity assigner assigns same ID to duplicate UMIs
        #[test]
        fn prop_identity_assigns_duplicates_to_same_id(umi in umi_strategy()) {
            let assigner = IdentityUmiAssigner::new();
            let umis = vec![umi.clone(), umi.clone(), umi.clone()];
            let assignments = assigner.assign(&umis);

            // All three should have the same ID
            let ids: HashSet<_> = assignments.iter().collect();
            assert_eq!(ids.len(), 1, "Duplicate UMIs should have same ID");
        }

        /// Property: Every input UMI gets an assignment
        #[test]
        fn prop_every_umi_gets_assignment(umis in prop::collection::vec(umi_strategy(), 1..20)) {
            let assigner = IdentityUmiAssigner::new();
            let assignments = assigner.assign(&umis);

            // Every UMI should be in the assignments (assignments Vec has same length as input)
            assert_eq!(assignments.len(), umis.len(), "Assignments Vec should match input length");
            // Every assignment should be valid (not None)
            for (i, assignment) in assignments.iter().enumerate() {
                assert_ne!(*assignment, MoleculeId::None, "UMI at index {i} should have valid assignment");
            }
        }

        /// Property: Paired UMI reverse is reversible
        #[test]
        fn prop_paired_reverse_is_reversible(umi in paired_umi_strategy()) {
            let reversed = PairedUmiAssigner::reverse(&umi).unwrap();
            let double_reversed = PairedUmiAssigner::reverse(&reversed).unwrap();
            assert_eq!(umi, double_reversed, "Reversing twice should return original");
        }

        /// Property: Paired UMI canonicalize is idempotent
        #[test]
        fn prop_paired_canonicalize_is_idempotent(umi in paired_umi_strategy()) {
            let canonical1 = PairedUmiAssigner::canonicalize_paired(&umi).unwrap();
            let canonical2 = PairedUmiAssigner::canonicalize_paired(&canonical1).unwrap();
            assert_eq!(canonical1, canonical2, "Canonicalizing twice should give same result");
        }

        /// Property: Paired UMI and its reverse have same canonical form
        #[test]
        fn prop_paired_canonicalize_matches_reverse(umi in paired_umi_strategy()) {
            let reversed = PairedUmiAssigner::reverse(&umi).unwrap();
            let canonical_forward = PairedUmiAssigner::canonicalize_paired(&umi).unwrap();
            let canonical_reverse = PairedUmiAssigner::canonicalize_paired(&reversed).unwrap();
            assert_eq!(
                canonical_forward,
                canonical_reverse,
                "UMI and its reverse should have same canonical form"
            );
        }

        /// Property: Adjacency assigner with 0 edits behaves like identity (fixed-length UMIs)
        #[test]
        fn prop_adjacency_zero_edits_like_identity(umis in prop::collection::vec("[ACGT]{8}", 1..10)) {
            let identity_assigner = IdentityUmiAssigner::new();
            let adjacency_assigner = AdjacencyUmiAssigner::new(0, 1, DEFAULT_INDEX_THRESHOLD);

            let identity_assignments = identity_assigner.assign(&umis);
            let adjacency_assignments = adjacency_assigner.assign(&umis);

            // Should have same number of unique molecule IDs
            let identity_ids: HashSet<_> = identity_assignments.iter().collect();
            let adjacency_ids: HashSet<_> = adjacency_assignments.iter().collect();
            assert_eq!(
                identity_ids.len(),
                adjacency_ids.len(),
                "With 0 edits, adjacency should create same number of groups as identity"
            );
        }

        /// Property: matches_within_threshold is consistent with count_mismatches
        #[test]
        fn prop_matches_within_threshold_consistent(a in "[ACGT]{8}", b in "[ACGT]{8}", max in 0usize..5) {
            let mismatches = count_mismatches(&a, &b);
            let matches = matches_within_threshold(&a, &b, max);
            assert_eq!(matches, mismatches <= max, "matches_within_threshold should match count_mismatches <= threshold");
        }

        /// Property: matches_within_threshold is symmetric
        #[test]
        fn prop_matches_within_threshold_symmetric(a in "[ACGT]{8}", b in "[ACGT]{8}", max in 0usize..5) {
            assert_eq!(
                matches_within_threshold(&a, &b, max),
                matches_within_threshold(&b, &a, max),
                "matches_within_threshold should be symmetric"
            );
        }
    }

    // ========================================================================
    // matches_within_threshold Unit Tests
    // ========================================================================

    #[test]
    fn test_matches_within_threshold_exact_match() {
        assert!(matches_within_threshold("ACGT", "ACGT", 0));
        assert!(matches_within_threshold("ACGT", "ACGT", 1));
        assert!(matches_within_threshold("ACGT", "ACGT", 5));
    }

    #[test]
    fn test_matches_within_threshold_one_mismatch() {
        // One mismatch at position 3
        assert!(!matches_within_threshold("ACGT", "ACGC", 0));
        assert!(matches_within_threshold("ACGT", "ACGC", 1));
        assert!(matches_within_threshold("ACGT", "ACGC", 2));
    }

    #[test]
    fn test_matches_within_threshold_multiple_mismatches() {
        // Two mismatches
        assert!(!matches_within_threshold("ACGT", "AATT", 0));
        assert!(!matches_within_threshold("ACGT", "AATT", 1));
        assert!(matches_within_threshold("ACGT", "AATT", 2));
        assert!(matches_within_threshold("ACGT", "AATT", 3));
    }

    #[test]
    fn test_matches_within_threshold_different_lengths() {
        assert!(!matches_within_threshold("ACGT", "ACG", 5));
        assert!(!matches_within_threshold("ACG", "ACGT", 5));
        assert!(!matches_within_threshold("", "ACGT", 10));
    }

    #[test]
    fn test_matches_within_threshold_empty_strings() {
        assert!(matches_within_threshold("", "", 0));
        assert!(matches_within_threshold("", "", 5));
    }

    #[test]
    fn test_matches_within_threshold_all_mismatches() {
        // All four positions are mismatches
        assert!(!matches_within_threshold("AAAA", "TTTT", 3));
        assert!(matches_within_threshold("AAAA", "TTTT", 4));
        assert!(matches_within_threshold("AAAA", "TTTT", 5));
    }

    #[test]
    fn test_matches_within_threshold_early_exit_behavior() {
        // Test that early exit works correctly by using long strings
        // With max_mismatches=1, should return false quickly when finding 2nd mismatch
        let a = "AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA";
        let b = "TTAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA"; // 2 mismatches at start

        // Should return false without checking the whole string
        assert!(!matches_within_threshold(a, b, 1));

        // Should return true with threshold of 2
        assert!(matches_within_threshold(a, b, 2));
    }

    // ========================================================================
    // Determinism Tests
    // ========================================================================

    /// Helper to verify determinism by checking that the grouping pattern is consistent.
    /// Two assignments are structurally equal if the same UMIs get the same MI within each run.
    fn assignments_structurally_equal(
        umis_a: &[Umi],
        a: &[MoleculeId],
        umis_b: &[Umi],
        b: &[MoleculeId],
    ) -> bool {
        if a.len() != b.len() {
            return false;
        }
        let map_a = build_assignment_map(umis_a, a);
        let map_b = build_assignment_map(umis_b, b);
        if map_a.len() != map_b.len() {
            return false;
        }
        // For each pair of UMIs, check if they have the same MI in both assignments
        let keys: Vec<_> = map_a.keys().collect();
        for i in 0..keys.len() {
            for j in i + 1..keys.len() {
                let same_in_a = map_a.get(keys[i]) == map_a.get(keys[j]);
                let same_in_b = map_b.get(keys[i]) == map_b.get(keys[j]);
                if same_in_a != same_in_b {
                    return false;
                }
            }
        }
        true
    }

    #[test]
    fn test_identity_assigner_deterministic() {
        // Verify that identity assigner produces consistent groupings on repeated runs
        // (the exact MI values may differ due to internal counter, but the structure should match)
        let assigner = IdentityUmiAssigner::new();
        let umis = vec![
            "GGGGGG".to_string(),
            "AAAAAA".to_string(),
            "TTTTTT".to_string(),
            "CCCCCC".to_string(),
            "AAAAAA".to_string(),
            "GGGGGG".to_string(),
        ];

        // Run multiple times and ensure grouping patterns are identical
        let assignments1 = assigner.assign(&umis);
        let assignments2 = assigner.assign(&umis);
        let assignments3 = assigner.assign(&umis);

        assert!(
            assignments_structurally_equal(&umis, &assignments1, &umis, &assignments2),
            "Assignment grouping should be deterministic"
        );
        assert!(
            assignments_structurally_equal(&umis, &assignments2, &umis, &assignments3),
            "Assignment grouping should be deterministic"
        );

        // Also verify that duplicate UMIs get the same MI within each run
        let map1 = build_assignment_map(&umis, &assignments1);
        assert_eq!(map1.get("AAAAAA"), map1.get("AAAAAA"), "Same UMI should get same MI");
    }

    #[test]
    fn test_simple_error_assigner_deterministic() {
        // Verify that edit distance assigner produces consistent groupings on repeated runs
        let assigner = SimpleErrorUmiAssigner::new(1);
        let umis = vec![
            "GGGGGG".to_string(),
            "GGGGGC".to_string(), // 1 edit from GGGGGG
            "AAAAAA".to_string(),
            "AAAAAC".to_string(), // 1 edit from AAAAAA
            "TTTTTT".to_string(),
        ];

        let assignments1 = assigner.assign(&umis);
        let assignments2 = assigner.assign(&umis);
        let assignments3 = assigner.assign(&umis);

        assert!(
            assignments_structurally_equal(&umis, &assignments1, &umis, &assignments2),
            "Edit assigner grouping should be deterministic"
        );
        assert!(
            assignments_structurally_equal(&umis, &assignments2, &umis, &assignments3),
            "Edit assigner grouping should be deterministic"
        );

        // Verify UMIs within 1 edit get same MI
        let map1 = build_assignment_map(&umis, &assignments1);
        assert_eq!(
            map1.get("GGGGGG"),
            map1.get("GGGGGC"),
            "UMIs within edit distance should get same MI"
        );
        assert_eq!(
            map1.get("AAAAAA"),
            map1.get("AAAAAC"),
            "UMIs within edit distance should get same MI"
        );
    }

    #[test]
    fn test_adjacency_assigner_deterministic() {
        // Verify that adjacency assigner produces consistent groupings on repeated runs
        let assigner = AdjacencyUmiAssigner::new(1, 1, DEFAULT_INDEX_THRESHOLD);
        let umis = vec![
            "GGGGGG".to_string(),
            "GGGGGC".to_string(),
            "GGGGGG".to_string(),
            "AAAAAA".to_string(),
            "AAAAAC".to_string(),
            "TTTTTT".to_string(),
            "TTTTTT".to_string(),
            "TTTTTT".to_string(),
        ];

        let assignments1 = assigner.assign(&umis);
        let assignments2 = assigner.assign(&umis);
        let assignments3 = assigner.assign(&umis);

        assert!(
            assignments_structurally_equal(&umis, &assignments1, &umis, &assignments2),
            "Adjacency assigner grouping should be deterministic"
        );
        assert!(
            assignments_structurally_equal(&umis, &assignments2, &umis, &assignments3),
            "Adjacency assigner grouping should be deterministic"
        );
    }

    #[test]
    fn test_adjacency_equal_count_deterministic() {
        // Regression test for deterministic tie-breaking when UMIs have equal counts.
        // Without lexicographic tie-breaking in the sort, the order of equal-count UMIs
        // depends on HashMap iteration order, causing non-deterministic groupings.
        //
        // Setup: AAAAAA and AAAGAC both have count=2, AAAAAC has count=1
        // Both AAAAAA and AAAGAC are within 1 edit of AAAAAC (can capture it)
        // AAAGAC is 2 edits from AAAAAA (positions 3 and 5), so they won't be grouped together
        // With deterministic sorting, AAAAAA (lexicographically first) should capture AAAAAC
        let assigner = AdjacencyUmiAssigner::new(1, 1, DEFAULT_INDEX_THRESHOLD);
        let umis = vec![
            "AAAAAA".to_string(),
            "AAAAAA".to_string(),
            "AAAGAC".to_string(),
            "AAAGAC".to_string(),
            "AAAAAC".to_string(), // count=1, within 1 edit of both AAAAAA and AAAGAC
        ];

        let assignments = assigner.assign(&umis);
        let map = build_assignment_map(&umis, &assignments);

        // AAAAAA should capture AAAAAC (both get same MI) because AAAAAA < AAAGAC lexicographically
        assert_eq!(
            map.get("AAAAAA"),
            map.get("AAAAAC"),
            "AAAAAA should capture AAAAAC (lexicographically first with equal count)"
        );

        // AAAGAC should be in its own group
        assert_ne!(map.get("AAAAAA"), map.get("AAAGAC"), "AAAGAC should be in a separate group");

        // Verify determinism across multiple runs
        let umis2 = vec![
            "AAAAAA".to_string(),
            "AAAAAA".to_string(),
            "AAAGAC".to_string(),
            "AAAGAC".to_string(),
            "AAAAAC".to_string(),
        ];
        let assignments2 = AdjacencyUmiAssigner::new(1, 1, DEFAULT_INDEX_THRESHOLD).assign(&umis2);
        assert!(
            assignments_structurally_equal(&umis, &assignments, &umis2, &assignments2),
            "Equal-count adjacency grouping should be deterministic across runs"
        );
    }

    // ==================== Index Threshold Tests ====================

    #[test]
    fn test_default_index_threshold_value() {
        // Verify the default threshold constant has the expected value
        assert_eq!(DEFAULT_INDEX_THRESHOLD, 100);
    }

    #[test]
    fn test_adjacency_index_vs_linear_produces_same_results() {
        // Generate enough UMIs to trigger index usage
        let mut umis = Vec::new();

        // Add 200 unique UMIs to exceed the default threshold
        for i in 0..200 {
            let umi = format!("{:0>8}", format!("{i:b}").replace('0', "A").replace('1', "C"));
            umis.push(umi.chars().take(8).collect::<String>());
        }

        // With index (threshold = 100, so 200 UMIs will trigger indexing)
        let indexed_assignments = AdjacencyUmiAssigner::new(1, 1, 100).assign(&umis);

        // Without index (threshold = 1000, so 200 UMIs will use linear scan)
        let linear_assignments = AdjacencyUmiAssigner::new(1, 1, 1000).assign(&umis);

        // Both should produce structurally equivalent results
        assert!(
            assignments_structurally_equal(&umis, &indexed_assignments, &umis, &linear_assignments),
            "Index-based and linear scan should produce equivalent groupings"
        );
    }

    #[test]
    fn test_adjacency_threshold_zero_always_uses_linear() {
        // threshold=0 should never use the index (always linear scan)
        // Need enough duplicates for count gradient: child < parent/2 + 1
        let mut umis = vec!["AAAAAAAA".to_string(); 10]; // Parent with 10 reads
        umis.push("AAAAAAAC".to_string()); // Child with 1 read (1 < 10/2 + 1 = 6)
        umis.extend(vec!["TTTTTTTT".to_string(); 5]); // Distinct UMI

        let assignments = AdjacencyUmiAssigner::new(1, 1, 0).assign(&umis);
        let map = build_assignment_map(&umis, &assignments);

        // Should still produce correct results
        assert_eq!(
            map.get("AAAAAAAA"),
            map.get("AAAAAAAC"),
            "1-mismatch UMIs should be grouped together due to count gradient"
        );
        assert_ne!(
            map.get("AAAAAAAA"),
            map.get("TTTTTTTT"),
            "Distinct UMIs should be in separate groups"
        );
    }

    #[test]
    fn test_strategy_new_assigner_full_passes_threshold() {
        // Verify that Strategy::new_assigner_full creates assigners with the correct threshold
        use super::Strategy as UmiStrategy;

        // Need count gradient for adjacency merging: child < parent/2 + 1
        let mut umis = vec!["AAAAAA".to_string(); 10]; // Parent with 10 reads
        umis.push("AAAAAC".to_string()); // Child with 1 read
        umis.extend(vec!["TTTTTT".to_string(); 5]); // Distinct UMI

        // Using Adjacency strategy with custom threshold
        let assigner = UmiStrategy::Adjacency.new_assigner_full(1, 1, 100);
        let assignments = assigner.assign(&umis);
        let map = build_assignment_map(&umis, &assignments);

        // Verify correct grouping behavior
        assert_eq!(
            map.get("AAAAAA"),
            map.get("AAAAAC"),
            "1-mismatch UMIs should be grouped with Adjacency strategy"
        );
        assert_ne!(map.get("AAAAAA"), map.get("TTTTTT"), "Distinct UMIs should be separate");
    }

    #[test]
    fn test_strategy_new_assigner_uses_default_threshold() {
        // Verify that new_assigner and new_assigner_threaded use DEFAULT_INDEX_THRESHOLD
        use super::Strategy as UmiStrategy;

        // Use simple UMIs to verify structural equality between assigners
        let mut umis = vec!["AAAAAA".to_string(); 10];
        umis.push("AAAAAC".to_string());
        umis.extend(vec!["TTTTTT".to_string(); 5]);

        // new_assigner should use DEFAULT_INDEX_THRESHOLD
        let assigner1 = UmiStrategy::Adjacency.new_assigner(1);
        let assignments1 = assigner1.assign(&umis);

        // new_assigner_full with explicit DEFAULT_INDEX_THRESHOLD should produce same result
        let assigner2 = UmiStrategy::Adjacency.new_assigner_full(1, 1, DEFAULT_INDEX_THRESHOLD);
        let assignments2 = assigner2.assign(&umis);

        assert!(
            assignments_structurally_equal(&umis, &assignments1, &umis, &assignments2),
            "new_assigner should use DEFAULT_INDEX_THRESHOLD"
        );
    }

    #[test]
    fn test_paired_assigner_respects_index_threshold() {
        // Verify PairedUmiAssigner passes threshold to underlying AdjacencyUmiAssigner
        // Need count gradient for adjacency merging
        let mut umis = vec!["AAAA-CCCC".to_string(); 10]; // Parent with 10 reads
        umis.push("AAAC-CCCC".to_string()); // Child with 1 read (1 mismatch in first half)
        umis.extend(vec!["TTTT-GGGG".to_string(); 5]); // Distinct UMI

        let assigner = PairedUmiAssigner::new_with_threads(1, 1, 100);
        let assignments = assigner.assign(&umis);
        let assignment_map = build_assignment_map(&umis, &assignments);

        // Should group the similar paired UMIs
        assert_eq!(
            assignment_map.get("AAAA-CCCC"),
            assignment_map.get("AAAC-CCCC"),
            "1-mismatch paired UMIs should be grouped due to count gradient"
        );
        assert_ne!(
            assignment_map.get("AAAA-CCCC"),
            assignment_map.get("TTTT-GGGG"),
            "Distinct paired UMIs should be separate"
        );
    }

    #[test]
    fn test_adjacency_large_dataset_with_indexing() {
        // Test with a dataset large enough to benefit from indexing
        let mut umis = Vec::new();

        // Create 150 unique UMIs - enough to trigger indexing at threshold=100
        for i in 0..150 {
            // Generate unique 8-character UMIs using base conversion
            let bases = ["A", "C", "G", "T"];
            let mut umi = String::new();
            let mut n = i;
            for _ in 0..8 {
                umi.push_str(bases[n % 4]);
                n /= 4;
            }
            umis.push(umi);
        }

        // Run with indexing (threshold=100, we have 150 UMIs)
        let assignments_indexed = AdjacencyUmiAssigner::new(1, 1, 100).assign(&umis);

        // Run without indexing (threshold > 150)
        let assignments_linear = AdjacencyUmiAssigner::new(1, 1, 200).assign(&umis);

        // Verify indexing produces same result as linear scan
        // (the actual number of groups depends on UMI edit distances)
        let groups_indexed: std::collections::HashSet<_> = assignments_indexed.iter().collect();
        let groups_linear: std::collections::HashSet<_> = assignments_linear.iter().collect();
        assert_eq!(
            groups_indexed.len(),
            groups_linear.len(),
            "Indexed and linear should produce same number of groups"
        );

        // Verify all UMIs got assigned
        assert_eq!(assignments_indexed.len(), 150, "All UMIs should be assigned");
    }

    // ==================== BkTree Unit Tests ====================

    #[test]
    fn test_bktree_new_empty() {
        let tree = BkTree::new();
        assert!(tree.is_empty());
        assert_eq!(tree.len(), 0);
    }

    #[test]
    fn test_bktree_from_umis_single() {
        let umis = vec![(0, "AAAAAAAA")];
        let tree = BkTree::from_umis(&umis).expect("Should build tree");
        assert_eq!(tree.len(), 1);
        assert!(!tree.is_empty());
    }

    #[test]
    fn test_bktree_from_umis_multiple() {
        let umis = vec![(0, "AAAAAAAA"), (1, "CCCCCCCC"), (2, "GGGGGGGG"), (3, "TTTTTTTT")];
        let tree = BkTree::from_umis(&umis).expect("Should build tree");
        assert_eq!(tree.len(), 4);
    }

    #[test]
    fn test_bktree_from_umis_invalid_base() {
        let umis = vec![(0, "AAAANAA")]; // N is not valid
        assert!(BkTree::from_umis(&umis).is_none());
    }

    #[test]
    fn test_bktree_find_exact_match() {
        let umis = vec![(0, "AAAAAAAA"), (1, "CCCCCCCC"), (2, "GGGGGGGG")];
        let tree = BkTree::from_umis(&umis).expect("Should build tree");

        let results = tree.find_within("AAAAAAAA", 0);
        assert_eq!(results.len(), 1);
        assert_eq!(results[0], (0, 0)); // index 0, distance 0
    }

    #[test]
    fn test_bktree_find_within_distance_1() {
        let umis = vec![
            (0, "AAAAAAAA"),
            (1, "AAAAAAAC"), // 1 mismatch from AAAAAAAA
            (2, "CCCCCCCC"), // 8 mismatches from AAAAAAAA
        ];
        let tree = BkTree::from_umis(&umis).expect("Should build tree");

        let results = tree.find_within("AAAAAAAA", 1);
        assert_eq!(results.len(), 2);

        let indices: HashSet<_> = results.iter().map(|(i, _)| *i).collect();
        assert!(indices.contains(&0));
        assert!(indices.contains(&1));
        assert!(!indices.contains(&2));
    }

    #[test]
    fn test_bktree_find_no_match() {
        let umis = vec![(0, "AAAAAAAA"), (1, "CCCCCCCC")];
        let tree = BkTree::from_umis(&umis).expect("Should build tree");

        // Query with very different UMI, distance 0
        let results = tree.find_within("TTTTTTTT", 0);
        assert!(results.is_empty());
    }

    #[test]
    fn test_bktree_find_invalid_query() {
        let umis = vec![(0, "AAAAAAAA")];
        let tree = BkTree::from_umis(&umis).expect("Should build tree");

        // Query with invalid base returns empty
        let results = tree.find_within("AAANAAAA", 1);
        assert!(results.is_empty());
    }

    #[test]
    fn test_bktree_triangle_inequality() {
        // Test that BK-tree correctly uses triangle inequality for pruning
        let umis = vec![
            (0, "AAAAAAAA"), // root
            (1, "AAAAAAAC"), // dist 1 from root
            (2, "AAAAAACC"), // dist 2 from root
            (3, "AAAAAACG"), // dist 2 from root
            (4, "CCCCCCCC"), // dist 8 from root
        ];
        let tree = BkTree::from_umis(&umis).expect("Should build tree");

        // Find all within distance 1 of AAAAAAAC
        let results = tree.find_within("AAAAAAAC", 1);

        // Should find: AAAAAAAC (dist 0), AAAAAAAA (dist 1), AAAAAACC (dist 1)
        assert_eq!(results.len(), 3);

        let result_map: std::collections::HashMap<_, _> = results.into_iter().collect();
        assert_eq!(result_map.get(&0), Some(&1)); // AAAAAAAA, dist 1
        assert_eq!(result_map.get(&1), Some(&0)); // AAAAAAAC, dist 0
        assert_eq!(result_map.get(&2), Some(&1)); // AAAAAACC, dist 1
    }

    // ==================== NgramIndex Unit Tests ====================

    #[test]
    fn test_ngram_index_new_single() {
        let umis = vec![(0, "AAAAAAAA")];
        let index = NgramIndex::new(&umis, 1);
        assert!(index.is_some());
    }

    #[test]
    fn test_ngram_index_new_empty() {
        let umis: Vec<(usize, &str)> = vec![];
        let index = NgramIndex::new(&umis, 1);
        assert!(index.is_none());
    }

    #[test]
    fn test_ngram_index_new_invalid_base() {
        let umis = vec![(0, "AAAANAA")];
        let index = NgramIndex::new(&umis, 1);
        assert!(index.is_none());
    }

    #[test]
    fn test_ngram_index_inconsistent_lengths() {
        let umis = vec![
            (0, "AAAAAAAA"),
            (1, "CCCCCC"), // Different length
        ];
        let index = NgramIndex::new(&umis, 1);
        assert!(index.is_none());
    }

    #[test]
    fn test_ngram_index_partition_too_short() {
        // 3bp UMI with 3 mismatches = 4 partitions = 0bp each
        // partition_len = 3 / 4 = 0, which returns None
        let umis = vec![(0, "AAA")];
        let index = NgramIndex::new(&umis, 3);
        assert!(index.is_none()); // partition_len is 0
    }

    #[test]
    fn test_ngram_index_find_exact_match() {
        let umis = vec![(0, "AAAAAAAA"), (1, "CCCCCCCC"), (2, "GGGGGGGG")];
        let index = NgramIndex::new(&umis, 1).expect("Should build index");

        let results = index.find_within("AAAAAAAA", 0);
        assert_eq!(results.len(), 1);
        assert_eq!(results[0], (0, 0));
    }

    #[test]
    fn test_ngram_index_find_within_distance_1() {
        let umis = vec![
            (0, "AAAAAAAA"),
            (1, "AAAAAAAC"), // 1 mismatch
            (2, "CCCCCCCC"), // 8 mismatches
        ];
        let index = NgramIndex::new(&umis, 1).expect("Should build index");

        let results = index.find_within("AAAAAAAA", 1);
        assert_eq!(results.len(), 2);

        let indices: HashSet<_> = results.iter().map(|(i, _)| *i).collect();
        assert!(indices.contains(&0));
        assert!(indices.contains(&1));
    }

    #[test]
    fn test_ngram_index_pigeonhole_principle() {
        // For k=1 mismatches, we have 2 partitions
        // If UMIs differ by 1, at least one partition must match
        let umis = vec![
            (0, "AAAACCCC"),
            (1, "AAAACCCG"), // 1 mismatch in second half
            (2, "GAAACCCC"), // 1 mismatch in first half
        ];
        let index = NgramIndex::new(&umis, 1).expect("Should build index");

        // Query AAAACCCC should find both 1-mismatch UMIs
        let results = index.find_within("AAAACCCC", 1);
        assert_eq!(results.len(), 3);

        let indices: HashSet<_> = results.iter().map(|(i, _)| *i).collect();
        assert!(indices.contains(&0)); // exact match
        assert!(indices.contains(&1)); // first partition matches (AAAA)
        assert!(indices.contains(&2)); // second partition matches (CCCC)
    }

    #[test]
    fn test_ngram_index_find_invalid_query() {
        let umis = vec![(0, "AAAAAAAA")];
        let index = NgramIndex::new(&umis, 1).expect("Should build index");

        // Query with invalid base returns empty
        let results = index.find_within("AAANAAAA", 1);
        assert!(results.is_empty());
    }

    #[test]
    fn test_ngram_index_find_wrong_length() {
        let umis = vec![(0, "AAAAAAAA")];
        let index = NgramIndex::new(&umis, 1).expect("Should build index");

        // Query with wrong length returns empty
        let results = index.find_within("AAAA", 1);
        assert!(results.is_empty());
    }

    #[test]
    fn test_bktree_vs_ngram_same_results() {
        // Both data structures should return the same results
        let umis = vec![
            (0, "AAAAAAAA"),
            (1, "AAAAAAAC"),
            (2, "AAAAAACT"),
            (3, "CCCCCCCC"),
            (4, "CCCCCCCT"),
            (5, "GGGGGGGG"),
        ];

        let tree = BkTree::from_umis(&umis).expect("Should build tree");
        let index = NgramIndex::new(&umis, 1).expect("Should build index");

        for (_, query) in &umis {
            let tree_results = tree.find_within(query, 1);
            let index_results = index.find_within(query, 1);

            let tree_set: HashSet<_> = tree_results.into_iter().collect();
            let index_set: HashSet<_> = index_results.into_iter().collect();

            assert_eq!(
                tree_set, index_set,
                "BkTree and NgramIndex should return same results for query {query}"
            );
        }
    }
}
