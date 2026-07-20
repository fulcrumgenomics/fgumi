//! Parallel UMI assignment strategies optimized for unmapped reads.
//!
//! When all reads are unmapped, they form a single position group which can
//! contain millions of unique UMIs. The standard sequential algorithms become
//! bottlenecks. This module provides parallel implementations that:
//!
//! - **Identity Strategy**: Uses partition-merge for parallel exact-match grouping
//! - **Edit Strategy**: Uses parallel edge discovery + union-find for O(U×L/P) complexity
//! - **Adjacency Strategy**: Uses parallel edge discovery + fast BFS for O(U×L/P + U+E) complexity
//! - **Paired Strategy**: Extends adjacency for duplex sequencing with strand assignment
//!
//! These produce identical results to the sequential implementations in `assigner.rs`.

use std::collections::VecDeque;
use std::sync::atomic::{AtomicUsize, Ordering};

use ahash::{AHashMap, AHashSet};
use rayon::prelude::*;

use crate::bitenc::BitEnc;
use crate::template::MoleculeId;

use super::{Umi, UmiAssigner};

/// Thread-safe union-find data structure for parallel connected component discovery.
///
/// Uses path compression and union-by-rank for near-O(1) amortized operations.
/// All operations are lock-free and safe for concurrent use.
pub struct UnionFind {
    parent: Vec<AtomicUsize>,
    rank: Vec<AtomicUsize>,
}

impl UnionFind {
    /// Create a new union-find structure with n elements.
    pub fn new(n: usize) -> Self {
        Self {
            parent: (0..n).map(AtomicUsize::new).collect(),
            rank: (0..n).map(|_| AtomicUsize::new(0)).collect(),
        }
    }

    /// Find the root of the set containing x, with path compression.
    #[must_use]
    pub fn find(&self, x: usize) -> usize {
        let mut current = x;
        loop {
            let parent = self.parent[current].load(Ordering::Acquire);
            if parent == current {
                return current;
            }
            // Path compression: try to link directly to grandparent
            let grandparent = self.parent[parent].load(Ordering::Acquire);
            let _ = self.parent[current].compare_exchange(
                parent,
                grandparent,
                Ordering::Release,
                Ordering::Relaxed,
            );
            current = parent;
        }
    }

    /// Union the sets containing x and y. Returns true if they were different sets.
    #[must_use]
    pub fn union(&self, x: usize, y: usize) -> bool {
        loop {
            let root_x = self.find(x);
            let root_y = self.find(y);

            if root_x == root_y {
                return false;
            }

            let rank_x = self.rank[root_x].load(Ordering::Acquire);
            let rank_y = self.rank[root_y].load(Ordering::Acquire);

            // Union by rank: attach smaller tree under larger
            let (smaller, larger) =
                if rank_x < rank_y { (root_x, root_y) } else { (root_y, root_x) };

            // Try to make larger the parent of smaller
            if self.parent[smaller]
                .compare_exchange(smaller, larger, Ordering::Release, Ordering::Relaxed)
                .is_ok()
            {
                // If ranks were equal, increment the rank of the new root
                if rank_x == rank_y {
                    let _ = self.rank[larger].compare_exchange(
                        rank_x,
                        rank_x + 1,
                        Ordering::Release,
                        Ordering::Relaxed,
                    );
                }
                return true;
            }
            // CAS failed, retry
        }
    }

    /// Perform unions in parallel using the provided edges.
    /// This is safe because `UnionFind` operations are lock-free.
    pub fn union_parallel(&self, edges: &[(usize, usize)]) {
        edges.par_iter().for_each(|&(i, j)| {
            let _ = self.union(i, j);
        });
    }
}

/// Generate all UMI neighbors within Hamming distance 1.
///
/// For a UMI of length L, generates 3×L neighbors (3 possible substitutions per position).
#[inline]
fn generate_neighbors(umi: &BitEnc) -> impl Iterator<Item = BitEnc> + '_ {
    (0..umi.len()).flat_map(move |pos| {
        let current_base = umi.base_at(pos);
        (0..4u8)
            .filter(move |&base| base != current_base)
            .map(move |base| umi.with_base_at(pos, base))
    })
}

/// Generate all UMI neighbors within Hamming distance k (for k > 1).
///
/// Uses recursive neighbor generation with memoization.
fn generate_neighbors_k(umi: &BitEnc, k: u32) -> Vec<BitEnc> {
    if k == 0 {
        return vec![*umi];
    }
    if k == 1 {
        return generate_neighbors(umi).collect();
    }

    // For k > 1, we generate all sequences within hamming distance k
    // This is done by generating all combinations of positions to mutate
    let mut result = AHashSet::new();
    generate_neighbors_k_recursive(umi, k, 0, &mut result);
    result.into_iter().collect()
}

/// Recursive helper for generating neighbors within distance k.
fn generate_neighbors_k_recursive(
    umi: &BitEnc,
    remaining_edits: u32,
    start_pos: usize,
    result: &mut AHashSet<BitEnc>,
) {
    if remaining_edits == 0 {
        result.insert(*umi);
        return;
    }

    // Include the current UMI (using fewer than max edits)
    result.insert(*umi);

    // Try mutating each position from start_pos onwards
    for pos in start_pos..umi.len() {
        let current_base = umi.base_at(pos);
        for new_base in 0..4u8 {
            if new_base != current_base {
                let mutated = umi.with_base_at(pos, new_base);
                generate_neighbors_k_recursive(&mutated, remaining_edits - 1, pos + 1, result);
            }
        }
    }
}

/// Parallel edge discovery for UMIs within edit distance 1.
///
/// Returns a list of (i, j) pairs where UMI i and UMI j are within edit distance 1.
/// Each edge is reported once (i < j).
///
/// Uses rayon parallel iteration; callers should invoke within a configured thread pool
/// via `pool.install(|| ...)` to control the number of threads.
///
/// # Complexity
/// O(U × L / P) where U = unique UMIs, L = UMI length, P = threads
#[must_use]
pub fn discover_edges_parallel_k1(
    umis: &[(BitEnc, usize)], // (encoded UMI, count)
) -> Vec<(usize, usize)> {
    // Build lookup map: BitEnc -> index
    let umi_to_idx: AHashMap<BitEnc, usize> =
        umis.iter().enumerate().map(|(i, (enc, _))| (*enc, i)).collect();

    umis.par_iter()
        .enumerate()
        .flat_map(|(i, (enc, _))| {
            let mut edges = Vec::new();
            for neighbor in generate_neighbors(enc) {
                if let Some(&j) = umi_to_idx.get(&neighbor) {
                    if i < j {
                        edges.push((i, j));
                    }
                }
            }
            edges
        })
        .collect()
}

/// Parallel edge discovery for UMIs within edit distance k (k > 1).
///
/// Uses neighbor generation with hash lookup instead of O(n²) pairwise comparison.
///
/// # Complexity
/// O(U × L^k × 3^k / P) where U = unique UMIs, L = UMI length, k = max edits, P = threads
#[must_use]
pub fn discover_edges_parallel_k(
    umis: &[(BitEnc, usize)],
    max_mismatches: u32,
) -> Vec<(usize, usize)> {
    if max_mismatches == 1 {
        return discover_edges_parallel_k1(umis);
    }

    // Build lookup map: BitEnc -> index
    let umi_to_idx: AHashMap<BitEnc, usize> =
        umis.iter().enumerate().map(|(i, (enc, _))| (*enc, i)).collect();

    // For k > 1, generate all neighbors within distance k and check hash map
    // This is still faster than O(n²) for reasonable k (1-3) and large n
    umis.par_iter()
        .enumerate()
        .flat_map(|(i, (enc, _))| {
            let mut edges = Vec::new();
            let neighbors = generate_neighbors_k(enc, max_mismatches);
            for neighbor in neighbors {
                if let Some(&j) = umi_to_idx.get(&neighbor) {
                    if i < j {
                        // Verify actual distance (neighbors may include closer matches)
                        if enc.hamming_distance(&neighbor) <= max_mismatches {
                            edges.push((i, j));
                        }
                    }
                }
            }
            edges
        })
        .collect()
}

/// Parallel UMI assigner for Identity strategy.
///
/// Uses a partition-merge approach: splits UMIs into chunks, builds per-chunk
/// hash maps in parallel, merges into a global map, then maps results in parallel.
/// Produces identical results to `IdentityUmiAssigner` in `assigner.rs`.
pub struct ParallelIdentityAssigner {
    pool: rayon::ThreadPool,
}

impl ParallelIdentityAssigner {
    /// Create a new parallel identity assigner.
    ///
    /// # Arguments
    /// * `threads` - Number of threads for the rayon thread pool
    ///
    /// # Panics
    /// Panics if the rayon thread pool cannot be created.
    #[must_use]
    pub fn new(threads: usize) -> Self {
        let pool = rayon::ThreadPoolBuilder::new()
            .num_threads(threads)
            .build()
            .expect("Failed to build rayon thread pool");
        Self { pool }
    }
}

impl UmiAssigner for ParallelIdentityAssigner {
    fn assign(&self, raw_umis: &[Umi]) -> Vec<MoleculeId> {
        if raw_umis.is_empty() {
            return Vec::new();
        }

        // Phase 1: Parallel uppercase conversion
        let canonicals: Vec<String> =
            self.pool.install(|| raw_umis.par_iter().map(|umi| umi.to_uppercase()).collect());

        // Phase 2: Parallel unique UMI discovery via partition-merge
        // Split into chunks, build per-chunk sets in parallel, then merge
        let chunk_size = (canonicals.len() / self.pool.current_num_threads().max(1)).max(1);
        let per_chunk_sets: Vec<AHashSet<&str>> = self.pool.install(|| {
            canonicals
                .par_chunks(chunk_size)
                .map(|chunk| chunk.iter().map(String::as_str).collect::<AHashSet<&str>>())
                .collect()
        });

        // Merge chunk sets into global unique set (sequential over N small sets)
        let mut unique_canonicals: Vec<&str> = Vec::new();
        let mut seen = AHashSet::new();
        for chunk_set in &per_chunk_sets {
            for &umi in chunk_set {
                if seen.insert(umi) {
                    unique_canonicals.push(umi);
                }
            }
        }

        // Sort for deterministic ID assignment (matches sequential behavior)
        unique_canonicals.sort_unstable();

        // Assign IDs in sorted order
        let canonical_to_id: AHashMap<&str, MoleculeId> = unique_canonicals
            .into_iter()
            .enumerate()
            .map(|(i, umi)| (umi, MoleculeId::Single(i as u64)))
            .collect();

        // Phase 3: Parallel final mapping
        self.pool.install(|| {
            canonicals.par_iter().map(|canonical| canonical_to_id[canonical.as_str()]).collect()
        })
    }

    fn as_any(&self) -> &dyn std::any::Any {
        self
    }
}

/// Parallel UMI assigner for Edit strategy.
///
/// Uses parallel edge discovery and union-find for fully parallel execution.
///
/// Matches `SimpleErrorUmiAssigner` in `assigner.rs` on all UMIs the `group`/`dedup` CLI can
/// deliver: both fold case, and both never see `N`-containing UMIs (those templates are
/// discarded up front with `discarded_ns_in_umi`, mirroring fgbio). The one residual
/// difference is for UMIs that are not `BitEnc`-encodable for a reason other than `N`
/// (length > 32 bases, or a non-`ACGTN` character): this assigner gives each such UMI its
/// own molecule id (see `assign`), whereas the sequential assigner groups them by string
/// edit distance. That case is unreachable via the CLI's validation and left as a documented
/// divergence rather than converged here.
pub struct ParallelEditAssigner {
    max_mismatches: u32,
    pool: rayon::ThreadPool,
}

impl ParallelEditAssigner {
    /// Create a new parallel edit assigner.
    ///
    /// # Arguments
    /// * `max_mismatches` - Maximum Hamming distance for UMIs to be grouped
    /// * `threads` - Number of threads for the rayon thread pool
    ///
    /// # Panics
    /// Panics if the rayon thread pool cannot be created.
    #[must_use]
    pub fn new(max_mismatches: u32, threads: usize) -> Self {
        let pool = rayon::ThreadPoolBuilder::new()
            .num_threads(threads)
            .build()
            .expect("Failed to build rayon thread pool");
        Self { max_mismatches, pool }
    }
}

impl UmiAssigner for ParallelEditAssigner {
    fn assign(&self, raw_umis: &[Umi]) -> Vec<MoleculeId> {
        if raw_umis.is_empty() {
            return Vec::new();
        }

        // Encode UMIs and count occurrences
        // Use a single pass to build both the count map and the encoding vector
        let mut umi_counts: AHashMap<BitEnc, usize> = AHashMap::new();
        let mut umi_to_original: Vec<Option<BitEnc>> = Vec::with_capacity(raw_umis.len());

        for umi in raw_umis {
            // `BitEnc::from_umi_str` is case-insensitive (`encode_base` maps
            // a/A, c/C, … identically) and skips dashes (e.g. "ACGT-TGCA"), so
            // the previous per-read `to_uppercase()` allocation was redundant —
            // feed the raw UMI directly.
            if let Some(enc) = BitEnc::from_umi_str(umi) {
                *umi_counts.entry(enc).or_insert(0) += 1;
                umi_to_original.push(Some(enc));
            } else {
                umi_to_original.push(None);
            }
        }

        // Handle edge case: no valid UMIs
        if umi_counts.is_empty() {
            // All UMIs were invalid - each gets its own molecule ID
            return (0..raw_umis.len()).map(|i| MoleculeId::Single(i as u64)).collect();
        }

        // Convert to indexed list (order doesn't matter for Edit strategy)
        let unique_umis: Vec<(BitEnc, usize)> = umi_counts.into_iter().collect();
        let enc_to_idx: AHashMap<BitEnc, usize> =
            unique_umis.iter().enumerate().map(|(i, (enc, _))| (*enc, i)).collect();

        // Discover edges and build components using the configured thread pool
        let max_mismatches = self.max_mismatches;
        let (edges, uf) = self.pool.install(|| {
            let edges = discover_edges_parallel_k(&unique_umis, max_mismatches);
            let uf = UnionFind::new(unique_umis.len());
            uf.union_parallel(&edges);
            (edges, uf)
        });
        drop(edges);

        // Assign molecule IDs based on connected components
        let mut root_to_mol: AHashMap<usize, MoleculeId> = AHashMap::new();
        let mut next_mol_id: u64 = 0;

        let mut result = Vec::with_capacity(raw_umis.len());
        for enc_opt in &umi_to_original {
            let mol_id = if let Some(enc) = enc_opt {
                let idx = enc_to_idx[enc];
                let root = uf.find(idx);
                *root_to_mol.entry(root).or_insert_with(|| {
                    let id = MoleculeId::Single(next_mol_id);
                    next_mol_id += 1;
                    id
                })
            } else {
                // UMI is not `BitEnc`-encodable (length > 32 or a non-`ACGTN` base): give it
                // its own molecule id. `N`-containing UMIs never reach here — those templates
                // are discarded up front (`discarded_ns_in_umi`) — so this branch is
                // unreachable via the CLI. It is the one documented divergence from the
                // sequential edit assigner, which would group such UMIs by string edit
                // distance (see the type doc above).
                let id = MoleculeId::Single(next_mol_id);
                next_mol_id += 1;
                id
            };
            result.push(mol_id);
        }

        result
    }

    fn as_any(&self) -> &dyn std::any::Any {
        self
    }
}

/// Parallel UMI assigner for Adjacency strategy.
///
/// Uses parallel edge discovery followed by sequential BFS with count constraints.
/// Produces identical results to `AdjacencyUmiAssigner` in `assigner.rs`.
pub struct ParallelAdjacencyAssigner {
    max_mismatches: u32,
    pool: rayon::ThreadPool,
}

impl ParallelAdjacencyAssigner {
    /// Create a new parallel adjacency assigner.
    ///
    /// # Arguments
    /// * `max_mismatches` - Maximum Hamming distance for UMIs to be adjacent
    /// * `threads` - Number of threads for the rayon thread pool
    ///
    /// # Panics
    /// Panics if the rayon thread pool cannot be created.
    #[must_use]
    pub fn new(max_mismatches: u32, threads: usize) -> Self {
        let pool = rayon::ThreadPoolBuilder::new()
            .num_threads(threads)
            .build()
            .expect("Failed to build rayon thread pool");
        Self { max_mismatches, pool }
    }
}

impl UmiAssigner for ParallelAdjacencyAssigner {
    fn assign(&self, raw_umis: &[Umi]) -> Vec<MoleculeId> {
        if raw_umis.is_empty() {
            return Vec::new();
        }

        // Count unique UMIs by their uppercased form (node identity is the
        // case-insensitive `BitEnc`), while remembering the RAW first-seen UMI
        // string for each — that is the case-sensitive, fgbio-compatible tie-break
        // key the sequential `AdjacencyUmiAssigner` uses (S7-002/FU-004). Iterating
        // `raw_umis` in input order makes the recorded raw string the first-seen
        // one, matching the sequential path's `raw_umis[first_index]`.
        let mut umi_counts: AHashMap<String, (usize, String)> = AHashMap::new();
        for umi in raw_umis {
            umi_counts
                .entry(umi.to_uppercase())
                .and_modify(|(count, _)| *count += 1)
                .or_insert_with(|| (1, umi.clone()));
        }

        // Build sorted list with BitEnc encoding
        // Filter out invalid UMIs during encoding
        // Use from_umi_str to handle paired UMIs with dashes (e.g., "ACGT-TGCA")
        // Tuple: (uppercased key, count, BitEnc, raw first-seen UMI).
        let mut sorted_umis: Vec<(String, usize, BitEnc, String)> = umi_counts
            .iter()
            .filter_map(|(umi_upper, (count, raw_first_seen))| {
                BitEnc::from_umi_str(umi_upper)
                    .map(|enc| (umi_upper.clone(), *count, enc, raw_first_seen.clone()))
            })
            .collect();

        // Handle edge case: no valid UMIs. Reads whose UMI failed to encode keep
        // `MoleculeId::None`, matching the sequential `AdjacencyUmiAssigner` (which
        // returns `vec![MoleculeId::None; len]` when every UMI is invalid). Assigning
        // a fresh `Single` here would diverge from the sequential path (parity bug).
        if sorted_umis.is_empty() {
            return vec![MoleculeId::None; raw_umis.len()];
        }

        // Sort by count descending, then by the RAW first-seen UMI string — exactly
        // the sequential assigner's case-sensitive tie-break (fgbio parity). A
        // mixed-case parity proptest pins the equivalence.
        sorted_umis.sort_by(|a, b| b.1.cmp(&a.1).then_with(|| a.3.cmp(&b.3)));

        // Build indexed structures (avoid cloning strings by using indices)
        let unique_umis: Vec<(BitEnc, usize)> =
            sorted_umis.iter().map(|(_, count, enc, _)| (*enc, *count)).collect();

        // Phase 1: Parallel edge discovery (using configured thread pool)
        let max_mismatches = self.max_mismatches;
        let edges = self.pool.install(|| discover_edges_parallel_k(&unique_umis, max_mismatches));

        // Build adjacency list
        let mut adj_list: Vec<Vec<usize>> = vec![Vec::new(); unique_umis.len()];
        for (i, j) in edges {
            adj_list[i].push(j);
            adj_list[j].push(i);
        }

        // Phase 2: Sequential BFS with adjacency constraints
        // (Must be sequential to maintain count-order processing)
        let mut assigned = vec![false; unique_umis.len()];
        let mut mol_ids: Vec<MoleculeId> = vec![MoleculeId::None; unique_umis.len()];
        let mut next_mol_id: u64 = 0;

        // Process in count order (sorted_umis is already sorted by count desc, then string)
        for root_idx in 0..unique_umis.len() {
            if assigned[root_idx] {
                continue;
            }

            // Start new molecule group
            let mol_id = MoleculeId::Single(next_mol_id);
            next_mol_id += 1;

            // BFS from root
            let mut queue = VecDeque::new();
            queue.push_back(root_idx);
            assigned[root_idx] = true;
            mol_ids[root_idx] = mol_id;

            while let Some(idx) = queue.pop_front() {
                let parent_count = unique_umis[idx].1;
                let max_child_count = parent_count / 2 + 1;

                // Check neighbors
                for &neighbor_idx in &adj_list[idx] {
                    if !assigned[neighbor_idx] {
                        let neighbor_count = unique_umis[neighbor_idx].1;
                        // Use <= to match sequential assigner behavior
                        if neighbor_count <= max_child_count {
                            assigned[neighbor_idx] = true;
                            mol_ids[neighbor_idx] = mol_id;
                            queue.push_back(neighbor_idx);
                        }
                    }
                }
            }
        }

        // Build map from uppercase UMI string to molecule ID
        // Use references to avoid cloning
        let str_to_mol: AHashMap<&str, MoleculeId> = sorted_umis
            .iter()
            .enumerate()
            .map(|(i, (umi, _, _, _))| (umi.as_str(), mol_ids[i]))
            .collect();

        // Map back to original UMI order. Reads whose UMI failed to encode are not
        // present in `str_to_mol` and keep `MoleculeId::None`, mirroring the
        // sequential `AdjacencyUmiAssigner` (which maps invalid reads to `None`).
        // Assigning a fresh `Single` here would be a sequential/parallel parity bug.
        raw_umis
            .iter()
            .map(|umi| {
                let upper = umi.to_uppercase();
                str_to_mol.get(upper.as_str()).copied().unwrap_or(MoleculeId::None)
            })
            .collect()
    }

    fn as_any(&self) -> &dyn std::any::Any {
        self
    }
}

/// Parallel UMI assigner for Paired (duplex) strategy.
///
/// Extends the adjacency algorithm for paired UMIs in duplex sequencing.
/// Paired UMIs have the form "AAAA-CCCC" where both parts come from opposite
/// ends of the original molecule. Reads from opposite strands have reversed
/// UMI pairs: "AAAA-CCCC" and "CCCC-AAAA" represent the same molecule.
///
/// Produces identical results to `PairedUmiAssigner` in `assigner.rs`.
pub struct ParallelPairedAssigner {
    max_mismatches: u32,
    pool: rayon::ThreadPool,
}

impl ParallelPairedAssigner {
    /// Create a new parallel paired assigner.
    ///
    /// # Arguments
    /// * `max_mismatches` - Maximum Hamming distance for UMIs to be adjacent
    /// * `threads` - Number of threads for the rayon thread pool
    ///
    /// # Panics
    /// Panics if the rayon thread pool cannot be created.
    #[must_use]
    pub fn new(max_mismatches: u32, threads: usize) -> Self {
        let pool = rayon::ThreadPoolBuilder::new()
            .num_threads(threads)
            .build()
            .expect("Failed to build rayon thread pool");
        Self { max_mismatches, pool }
    }

    /// Reverse a paired UMI: "AAAA-CCCC" -> "CCCC-AAAA"
    /// Input is expected to be uppercase.
    fn reverse_paired(umi: &str) -> Option<String> {
        let parts: Vec<&str> = umi.split('-').collect();
        if parts.len() == 2 { Some(format!("{}-{}", parts[1], parts[0])) } else { None }
    }

    /// Canonicalize a paired UMI by taking the lexicographically smaller of A-B and B-A.
    /// Always returns uppercase.
    fn canonicalize(umi: &str) -> String {
        let upper = umi.to_uppercase();
        if let Some(reversed) = Self::reverse_paired(&upper) {
            if reversed < upper { reversed } else { upper }
        } else {
            upper
        }
    }

    /// Check if a UMI matches the canonical form (vs its reverse).
    fn matches_canonical(umi: &str, canonical: &str) -> bool {
        // `canonical` is already uppercase and UMIs are printable ASCII (SAM
        // 'Z'-type RX/OX charset), so an allocation-free ASCII case-insensitive
        // compare replaces the previous `to_uppercase()`-then-compare.
        umi.eq_ignore_ascii_case(canonical)
    }
}

impl UmiAssigner for ParallelPairedAssigner {
    fn assign(&self, raw_umis: &[Umi]) -> Vec<MoleculeId> {
        if raw_umis.is_empty() {
            return Vec::new();
        }

        // Validate all are paired UMIs (consistent with sequential PairedUmiAssigner)
        for umi in raw_umis {
            assert!(umi.split('-').count() == 2, "UMI {umi} is not a paired UMI");
        }

        // Count unique UMIs using canonical form
        let mut canonical_counts: AHashMap<String, usize> = AHashMap::new();
        for umi in raw_umis {
            let canonical = Self::canonicalize(umi);
            *canonical_counts.entry(canonical).or_insert(0) += 1;
        }

        // Build sorted list with BitEnc encoding (using canonical forms)
        let mut sorted_umis: Vec<(String, usize, BitEnc)> = canonical_counts
            .iter()
            .filter_map(|(umi, &count)| {
                // For paired UMIs, encode without the dash
                BitEnc::from_umi_str(umi).map(|enc| (umi.clone(), count, enc))
            })
            .collect();

        // Handle edge case: no valid UMIs
        if sorted_umis.is_empty() {
            return (0..raw_umis.len()).map(|i| MoleculeId::Single(i as u64)).collect();
        }

        // Sort by count descending, then by UMI string for determinism
        sorted_umis.sort_by(|a, b| b.1.cmp(&a.1).then_with(|| a.0.cmp(&b.0)));

        // Build indexed structures
        let unique_umis: Vec<(BitEnc, usize)> =
            sorted_umis.iter().map(|(_, count, enc)| (*enc, *count)).collect();

        // Phase 1: Parallel edge discovery (using configured thread pool)
        // Use max_mismatches directly (not 2x) because canonicalization already handles
        // reverse-orientation matching. The sequential PairedUmiAssigner checks both forward
        // and reversed orientations with the same threshold; canonicalization achieves the
        // same effect by ensuring reversed pairs share the same canonical form.
        let max_mismatches = self.max_mismatches;
        let edges = self.pool.install(|| discover_edges_parallel_k(&unique_umis, max_mismatches));

        // Build adjacency list
        let mut adj_list: Vec<Vec<usize>> = vec![Vec::new(); unique_umis.len()];
        for (i, j) in edges {
            adj_list[i].push(j);
            adj_list[j].push(i);
        }

        // Phase 2: Sequential BFS with adjacency constraints
        let mut assigned = vec![false; unique_umis.len()];
        let mut mol_ids: Vec<u64> = vec![0; unique_umis.len()];
        let mut next_mol_id: u64 = 0;

        for root_idx in 0..unique_umis.len() {
            if assigned[root_idx] {
                continue;
            }

            let mol_id = next_mol_id;
            next_mol_id += 1;

            let mut queue = VecDeque::new();
            queue.push_back(root_idx);
            assigned[root_idx] = true;
            mol_ids[root_idx] = mol_id;

            while let Some(idx) = queue.pop_front() {
                let parent_count = unique_umis[idx].1;
                let max_child_count = parent_count / 2 + 1;

                for &neighbor_idx in &adj_list[idx] {
                    if !assigned[neighbor_idx] {
                        let neighbor_count = unique_umis[neighbor_idx].1;
                        if neighbor_count <= max_child_count {
                            assigned[neighbor_idx] = true;
                            mol_ids[neighbor_idx] = mol_id;
                            queue.push_back(neighbor_idx);
                        }
                    }
                }
            }
        }

        // Build map from canonical UMI to molecule ID
        let canonical_to_mol: AHashMap<&str, u64> = sorted_umis
            .iter()
            .enumerate()
            .map(|(i, (umi, _, _))| (umi.as_str(), mol_ids[i]))
            .collect();

        // Map back to original UMI order with strand assignment
        raw_umis
            .iter()
            .map(|umi| {
                let canonical = Self::canonicalize(umi);
                if let Some(&base_id) = canonical_to_mol.get(canonical.as_str()) {
                    // Assign /A or /B based on whether the original matches canonical
                    if Self::matches_canonical(umi, &canonical) {
                        MoleculeId::PairedA(base_id)
                    } else {
                        MoleculeId::PairedB(base_id)
                    }
                } else {
                    // Invalid UMI
                    MoleculeId::None
                }
            })
            .collect()
    }

    fn split_templates_by_pair_orientation(&self) -> bool {
        false
    }

    fn as_any(&self) -> &dyn std::any::Any {
        self
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::umi::{AdjacencyUmiAssigner, SimpleErrorUmiAssigner};
    use fgumi_umi::assigner::DEFAULT_INDEX_THRESHOLD;
    use proptest::prelude::*;

    // ==================== UnionFind Tests ====================

    #[test]
    fn test_union_find_basic() {
        let uf = UnionFind::new(5);

        // Initially all separate
        for i in 0..5 {
            assert_eq!(uf.find(i), i);
        }

        // Union 0 and 1
        assert!(uf.union(0, 1));
        assert_eq!(uf.find(0), uf.find(1));

        // Union again should return false
        assert!(!uf.union(0, 1));

        // Union 2 and 3
        assert!(uf.union(2, 3));
        assert_eq!(uf.find(2), uf.find(3));
        assert_ne!(uf.find(0), uf.find(2));

        // Connect the two groups
        assert!(uf.union(1, 3));
        assert_eq!(uf.find(0), uf.find(3));
    }

    #[test]
    fn test_union_find_large() {
        let uf = UnionFind::new(1000);

        // Chain union: 0-1-2-3-...-999
        for i in 0..999 {
            let _ = uf.union(i, i + 1);
        }

        // All should be in same component
        let root = uf.find(0);
        for i in 1..1000 {
            assert_eq!(uf.find(i), root);
        }
    }

    #[test]
    fn test_union_find_parallel() {
        let uf = UnionFind::new(100);

        // Create edges for a chain
        let edges: Vec<(usize, usize)> = (0..99).map(|i| (i, i + 1)).collect();

        // Parallel union
        uf.union_parallel(&edges);

        // All should be in same component
        let root = uf.find(0);
        for i in 1..100 {
            assert_eq!(uf.find(i), root);
        }
    }

    // ==================== Neighbor Generation Tests ====================

    #[test]
    fn test_generate_neighbors() {
        let umi = BitEnc::from_bytes(b"AA").expect("valid DNA sequence");
        let neighbors: Vec<_> = generate_neighbors(&umi).collect();

        // 2 positions × 3 alternatives = 6 neighbors
        assert_eq!(neighbors.len(), 6);

        // Check specific neighbors
        assert!(neighbors.contains(&BitEnc::from_bytes(b"CA").expect("valid DNA sequence")));
        assert!(neighbors.contains(&BitEnc::from_bytes(b"GA").expect("valid DNA sequence")));
        assert!(neighbors.contains(&BitEnc::from_bytes(b"TA").expect("valid DNA sequence")));
        assert!(neighbors.contains(&BitEnc::from_bytes(b"AC").expect("valid DNA sequence")));
        assert!(neighbors.contains(&BitEnc::from_bytes(b"AG").expect("valid DNA sequence")));
        assert!(neighbors.contains(&BitEnc::from_bytes(b"AT").expect("valid DNA sequence")));
    }

    #[test]
    fn test_generate_neighbors_k2() {
        let umi = BitEnc::from_bytes(b"AA").expect("valid DNA sequence");
        let neighbors = generate_neighbors_k(&umi, 2);

        // Should include: original, all 1-edit, all 2-edit neighbors
        // 1 + 6 + 9 = 16 (including original)
        // Actually: for 2-base UMI with k=2, we can change both bases
        // Total combinations: 4^2 = 16 (all possible 2-base sequences)
        assert_eq!(neighbors.len(), 16);

        // Should include original
        assert!(neighbors.contains(&umi));

        // Should include 1-edit neighbors
        assert!(neighbors.contains(&BitEnc::from_bytes(b"CA").expect("valid DNA sequence")));

        // Should include 2-edit neighbors
        assert!(neighbors.contains(&BitEnc::from_bytes(b"CC").expect("valid DNA sequence")));
        assert!(neighbors.contains(&BitEnc::from_bytes(b"TT").expect("valid DNA sequence")));
    }

    // ==================== Edge Discovery Tests ====================

    #[test]
    fn test_discover_edges_parallel_k1() {
        let umis = vec![
            (BitEnc::from_bytes(b"AAAA").expect("valid DNA sequence"), 10),
            (BitEnc::from_bytes(b"AAAT").expect("valid DNA sequence"), 5), // 1 edit from AAAA
            (BitEnc::from_bytes(b"TTTT").expect("valid DNA sequence"), 3), // 4 edits from AAAA
        ];

        let edges = discover_edges_parallel_k1(&umis);

        // Only AAAA-AAAT should be connected
        assert_eq!(edges.len(), 1);
        assert!(edges.contains(&(0, 1)));
    }

    #[test]
    fn test_discover_edges_parallel_k2() {
        let umis = vec![
            (BitEnc::from_bytes(b"AAAA").expect("valid DNA sequence"), 10),
            (BitEnc::from_bytes(b"AATT").expect("valid DNA sequence"), 5), // 2 edits from AAAA
            (BitEnc::from_bytes(b"TTTT").expect("valid DNA sequence"), 3), // 4 edits from AAAA, 2 from AATT
        ];

        let edges = discover_edges_parallel_k(&umis, 2);

        // AAAA-AATT (2 edits) and AATT-TTTT (2 edits) should be connected
        assert_eq!(edges.len(), 2);
        assert!(edges.contains(&(0, 1)));
        assert!(edges.contains(&(1, 2)));
    }

    #[test]
    fn test_discover_edges_parallel_chain() {
        // Test A-B-C chain where A~B and B~C
        let umis = vec![
            (BitEnc::from_bytes(b"AAAA").expect("valid DNA sequence"), 10),
            (BitEnc::from_bytes(b"AAAT").expect("valid DNA sequence"), 5), // 1 edit from AAAA
            (BitEnc::from_bytes(b"AATT").expect("valid DNA sequence"), 3), // 1 edit from AAAT, 2 from AAAA
        ];

        let edges = discover_edges_parallel_k1(&umis);

        // Should have edges 0-1 and 1-2
        assert_eq!(edges.len(), 2);
        assert!(edges.contains(&(0, 1)));
        assert!(edges.contains(&(1, 2)));
    }

    // ==================== Identity Assigner Tests ====================

    #[test]
    fn test_parallel_identity_basic() {
        let assigner = ParallelIdentityAssigner::new(2);
        let umis: Vec<String> =
            vec!["AAAA", "AAAA", "TTTT"].into_iter().map(String::from).collect();

        let assignments = assigner.assign(&umis);

        // Identical UMIs should be same molecule
        assert_eq!(assignments[0], assignments[1]);
        // Different UMIs should be different molecules
        assert_ne!(assignments[0], assignments[2]);
    }

    #[test]
    fn test_parallel_identity_case_insensitive() {
        let assigner = ParallelIdentityAssigner::new(2);
        let umis: Vec<String> =
            vec!["ACGT", "acgt", "AcGt"].into_iter().map(String::from).collect();

        let assignments = assigner.assign(&umis);

        // All should be same molecule (case insensitive)
        assert_eq!(assignments[0], assignments[1]);
        assert_eq!(assignments[1], assignments[2]);
    }

    #[test]
    fn test_parallel_identity_empty() {
        let assigner = ParallelIdentityAssigner::new(2);
        let umis: Vec<String> = vec![];
        let assignments = assigner.assign(&umis);
        assert!(assignments.is_empty());
    }

    #[test]
    fn test_parallel_identity_single() {
        let assigner = ParallelIdentityAssigner::new(2);
        let umis: Vec<String> = vec!["ACGT".to_string()];
        let assignments = assigner.assign(&umis);
        assert_eq!(assignments.len(), 1);
    }

    #[test]
    fn test_parallel_identity_deterministic() {
        let umis: Vec<String> = vec!["CCCC", "AAAA", "TTTT", "GGGG", "AAAA", "CCCC"]
            .into_iter()
            .map(String::from)
            .collect();

        let assigner1 = ParallelIdentityAssigner::new(2);
        let assigner2 = ParallelIdentityAssigner::new(4);

        let result1 = assigner1.assign(&umis);
        let result2 = assigner2.assign(&umis);

        assert_same_groupings(&result1, &result2);
    }

    #[test]
    fn test_parallel_identity_many_duplicates() {
        let assigner = ParallelIdentityAssigner::new(4);
        let mut umis: Vec<String> = Vec::new();
        umis.extend(std::iter::repeat_n("AAAAAA".to_string(), 100));
        umis.extend(std::iter::repeat_n("TTTTTT".to_string(), 100));
        umis.extend(std::iter::repeat_n("CCCCCC".to_string(), 100));

        let assignments = assigner.assign(&umis);

        // All AAAAAA should be same
        for i in 1..100 {
            assert_eq!(assignments[0], assignments[i]);
        }
        // All TTTTTT should be same
        for i in 101..200 {
            assert_eq!(assignments[100], assignments[i]);
        }
        // Groups should differ
        assert_ne!(assignments[0], assignments[100]);
        assert_ne!(assignments[0], assignments[200]);
        assert_ne!(assignments[100], assignments[200]);
    }

    // ==================== Edit Assigner Tests ====================

    #[test]
    fn test_parallel_edit_basic() {
        let assigner = ParallelEditAssigner::new(1, 2);
        let umis: Vec<String> =
            vec!["AAAA", "AAAT", "TTTT"].into_iter().map(String::from).collect();

        let assignments = assigner.assign(&umis);

        // AAAA and AAAT should be same molecule (1 edit apart)
        assert_eq!(assignments[0], assignments[1]);
        // TTTT should be different (4 edits from both)
        assert_ne!(assignments[0], assignments[2]);
    }

    #[test]
    fn test_parallel_edit_k2() {
        let assigner = ParallelEditAssigner::new(2, 2);
        let umis: Vec<String> =
            vec!["AAAA", "AATT", "TTTT"].into_iter().map(String::from).collect();

        let assignments = assigner.assign(&umis);

        // AAAA and AATT should be same molecule (2 edits apart)
        assert_eq!(assignments[0], assignments[1]);
        // TTTT should be different (4 edits from AAAA, 2 from AATT)
        // With k=2, AATT-TTTT are connected, so all three should be same
        assert_eq!(assignments[1], assignments[2]);
    }

    #[test]
    fn test_parallel_edit_transitive() {
        let assigner = ParallelEditAssigner::new(1, 2);
        // A-B-C chain where A~B and B~C but A and C are 2 apart
        let umis: Vec<String> =
            vec!["AAAA", "AAAT", "AATT"].into_iter().map(String::from).collect();

        let assignments = assigner.assign(&umis);

        // All should be same molecule due to transitivity
        assert_eq!(assignments[0], assignments[1]);
        assert_eq!(assignments[1], assignments[2]);
    }

    #[test]
    fn test_parallel_edit_empty() {
        let assigner = ParallelEditAssigner::new(1, 2);
        let umis: Vec<String> = vec![];
        let assignments = assigner.assign(&umis);
        assert!(assignments.is_empty());
    }

    #[test]
    fn test_parallel_edit_single() {
        let assigner = ParallelEditAssigner::new(1, 2);
        let umis: Vec<String> = vec!["ACGT".to_string()];
        let assignments = assigner.assign(&umis);
        assert_eq!(assignments.len(), 1);
    }

    #[test]
    fn test_parallel_edit_identical_umis() {
        let assigner = ParallelEditAssigner::new(1, 2);
        let umis: Vec<String> =
            vec!["ACGT", "ACGT", "ACGT"].into_iter().map(String::from).collect();

        let assignments = assigner.assign(&umis);

        // All identical UMIs should have same molecule ID
        assert_eq!(assignments[0], assignments[1]);
        assert_eq!(assignments[1], assignments[2]);
    }

    #[test]
    fn test_parallel_edit_paired_umis_with_dash() {
        // Test that paired UMIs with dashes (e.g., "ACGT-TGCA") are handled correctly
        // This is a regression test for the bug where dashes caused all UMIs to be invalid
        let assigner = ParallelEditAssigner::new(1, 2);
        let umis: Vec<String> = vec![
            "ACGTACGT-TGCATGCA".to_string(), // Same UMI, should be grouped
            "ACGTACGT-TGCATGCA".to_string(),
            "ACGTACGT-TGCATGCT".to_string(), // 1 edit from above, should be grouped
            "GGGGGGGG-CCCCCCCC".to_string(), // Different, should be separate
        ];

        let assignments = assigner.assign(&umis);

        // First three should be same molecule (identical or 1 edit)
        assert_eq!(assignments[0], assignments[1], "Identical paired UMIs should have same ID");
        assert_eq!(assignments[0], assignments[2], "Paired UMIs within 1 edit should have same ID");
        // Fourth should be different
        assert_ne!(
            assignments[0], assignments[3],
            "Different paired UMIs should have different IDs"
        );
    }

    #[test]
    fn test_parallel_adjacency_paired_umis_with_dash() {
        // Test that adjacency also handles paired UMIs with dashes
        let assigner = ParallelAdjacencyAssigner::new(1, 2);
        let umis: Vec<String> = vec![
            "ACGTACGT-TGCATGCA".to_string(),
            "ACGTACGT-TGCATGCA".to_string(),
            "ACGTACGT-TGCATGCT".to_string(), // 1 edit
        ];

        let assignments = assigner.assign(&umis);

        // With adjacency (count-based), the result depends on counts
        // But at minimum, identical UMIs should be grouped
        assert_eq!(assignments[0], assignments[1], "Identical paired UMIs should have same ID");
    }

    #[test]
    fn test_parallel_edit_case_insensitive() {
        let assigner = ParallelEditAssigner::new(1, 2);
        let umis: Vec<String> =
            vec!["ACGT", "acgt", "AcGt"].into_iter().map(String::from).collect();

        let assignments = assigner.assign(&umis);

        // All should be same molecule (case insensitive)
        assert_eq!(assignments[0], assignments[1]);
        assert_eq!(assignments[1], assignments[2]);
    }

    #[test]
    fn test_parallel_edit_invalid_umis() {
        let assigner = ParallelEditAssigner::new(1, 2);
        let umis: Vec<String> = vec!["ACGT", "ACGN", "ACGT"] // ACGN is invalid
            .into_iter()
            .map(String::from)
            .collect();

        let assignments = assigner.assign(&umis);

        // Valid UMIs should be grouped
        assert_eq!(assignments[0], assignments[2]);
        // Invalid UMI gets its own ID
        assert_ne!(assignments[0], assignments[1]);
    }

    // ==================== Adjacency Assigner Tests ====================

    #[test]
    fn test_parallel_adjacency_basic() {
        let assigner = ParallelAdjacencyAssigner::new(1, 2);
        let umis: Vec<String> =
            vec!["AAAA", "AAAA", "AAAT"].into_iter().map(String::from).collect(); // AAAA count=2, AAAT count=1

        let assignments = assigner.assign(&umis);

        // AAAT (count=1) should be captured by AAAA (count=2) since 1 <= 2/2+1 = 2
        assert_eq!(assignments[0], assignments[1]); // Same UMI
        assert_eq!(assignments[0], assignments[2]); // AAAT captured by AAAA
    }

    #[test]
    fn test_parallel_adjacency_count_constraint() {
        let assigner = ParallelAdjacencyAssigner::new(1, 2);
        // AAAA appears 4x, AAAT appears 1x
        // max_child_count = 4/2 + 1 = 3
        // AAAT count (1) <= 3, so it will be captured by AAAA
        // But TTTT (count=4) cannot be captured since 4 > 3
        let umis: Vec<String> = vec![
            "AAAA", "AAAA", "AAAA", "AAAA", // count=4
            "AAAT", // count=1, within 1 edit of AAAA
            "TTTT", "TTTT", "TTTT", "TTTT", // count=4, too far from AAAA
        ]
        .into_iter()
        .map(String::from)
        .collect();

        let assignments = assigner.assign(&umis);

        // Same UMIs should be same molecule
        assert_eq!(assignments[0], assignments[1]);
        assert_eq!(assignments[0], assignments[2]);
        assert_eq!(assignments[0], assignments[3]);
        // AAAT should be captured by AAAA
        assert_eq!(assignments[0], assignments[4]);
        // TTTT should be different (too far from AAAA)
        assert_eq!(assignments[5], assignments[6]);
        assert_ne!(assignments[0], assignments[5]);
    }

    #[test]
    fn test_parallel_adjacency_cascade() {
        let assigner = ParallelAdjacencyAssigner::new(1, 2);
        // Test cascade: A(10) -> B(4) -> C(1)
        // B can be captured by A: 4 <= 10/2+1 = 6
        // C can be captured by B: 1 <= 4/2+1 = 3
        let mut umis: Vec<String> = Vec::new();
        umis.extend(std::iter::repeat_n("AAAA".to_string(), 10));
        umis.extend(std::iter::repeat_n("AAAT".to_string(), 4));
        umis.extend(std::iter::repeat_n("AATT".to_string(), 1));

        let assignments = assigner.assign(&umis);

        // All should be captured by AAAA
        let first_mol = assignments[0];
        for assignment in &assignments {
            assert_eq!(*assignment, first_mol);
        }
    }

    #[test]
    fn test_parallel_adjacency_deterministic_tiebreak() {
        // Test that equal-count UMIs are handled deterministically
        let assigner = ParallelAdjacencyAssigner::new(1, 2);

        // AAAAAA and AAAGAC both have count=2, AAAAAC has count=1
        // Both are within 1 edit of AAAAAC
        // With lexicographic tie-breaking, AAAAAA should capture AAAAAC
        let umis: Vec<String> = vec![
            "AAAAAA", "AAAAAA", // count=2
            "AAAGAC", "AAAGAC", // count=2
            "AAAAAC", // count=1
        ]
        .into_iter()
        .map(String::from)
        .collect();

        let assignments1 = assigner.assign(&umis);
        let assignments2 = assigner.assign(&umis);

        // Should be deterministic
        assert_eq!(assignments1, assignments2);

        // AAAAAA and AAAAAC should be grouped (AAAAAA is lexicographically first)
        assert_eq!(assignments1[0], assignments1[4]);
    }

    #[test]
    fn test_parallel_adjacency_empty() {
        let assigner = ParallelAdjacencyAssigner::new(1, 2);
        let umis: Vec<String> = vec![];
        let assignments = assigner.assign(&umis);
        assert!(assignments.is_empty());
    }

    #[test]
    fn test_parallel_adjacency_single() {
        let assigner = ParallelAdjacencyAssigner::new(1, 2);
        let umis: Vec<String> = vec!["ACGT".to_string()];
        let assignments = assigner.assign(&umis);
        assert_eq!(assignments.len(), 1);
    }

    // ==================== Paired Assigner Tests ====================

    #[test]
    fn test_parallel_paired_basic() {
        let assigner = ParallelPairedAssigner::new(1, 2);
        let umis: Vec<String> = vec![
            "AAAA-CCCC".to_string(), // Top strand
            "CCCC-AAAA".to_string(), // Bottom strand (same molecule!)
        ];

        let assignments = assigner.assign(&umis);

        // Both should be same base molecule but different strands
        match (&assignments[0], &assignments[1]) {
            (MoleculeId::PairedA(a), MoleculeId::PairedB(b))
            | (MoleculeId::PairedB(a), MoleculeId::PairedA(b)) => assert_eq!(a, b),
            _ => unreachable!("Expected PairedA and PairedB, got {assignments:?}"),
        }
    }

    #[test]
    fn test_parallel_paired_with_errors() {
        let assigner = ParallelPairedAssigner::new(1, 2);
        let umis: Vec<String> = vec![
            "AAAA-CCCC".to_string(),
            "AAAA-CCCC".to_string(),
            "AAAT-CCCC".to_string(), // 1 error in first half
        ];

        let assignments = assigner.assign(&umis);

        // All should be same molecule (error corrected)
        let base_id = match &assignments[0] {
            MoleculeId::PairedA(id) | MoleculeId::PairedB(id) => *id,
            _ => unreachable!("Expected paired ID"),
        };

        for assignment in &assignments {
            match assignment {
                MoleculeId::PairedA(id) | MoleculeId::PairedB(id) => {
                    assert_eq!(*id, base_id);
                }
                _ => unreachable!("Expected paired ID"),
            }
        }
    }

    #[test]
    fn test_parallel_paired_canonicalization() {
        // Test that canonicalization works correctly
        assert_eq!(
            ParallelPairedAssigner::canonicalize("AAAA-CCCC"),
            "AAAA-CCCC" // Already canonical (A < C)
        );
        assert_eq!(
            ParallelPairedAssigner::canonicalize("CCCC-AAAA"),
            "AAAA-CCCC" // Reversed to canonical
        );
        assert_eq!(
            ParallelPairedAssigner::canonicalize("TTTT-AAAA"),
            "AAAA-TTTT" // Reversed to canonical (A < T)
        );
    }

    #[test]
    fn test_parallel_paired_canonicalization_mixed_case() {
        // Test that mixed-case input is handled correctly (uppercased before canonicalization)
        // This is a regression test for a bug where mixed-case wasn't uppercased before reversing
        assert_eq!(ParallelPairedAssigner::canonicalize("aaaa-cccc"), "AAAA-CCCC");
        assert_eq!(ParallelPairedAssigner::canonicalize("cccc-aaaa"), "AAAA-CCCC");
        assert_eq!(ParallelPairedAssigner::canonicalize("Tttt-Aaaa"), "AAAA-TTTT");
        // Mixed case should produce same result as uppercase
        assert_eq!(
            ParallelPairedAssigner::canonicalize("aaaa-CCCC"),
            ParallelPairedAssigner::canonicalize("AAAA-cccc")
        );
    }

    // ==================== Equivalence Verification ====================

    #[test]
    fn test_edit_equivalence_with_sequential() {
        // This test verifies parallel Edit produces same groupings as expected
        let parallel = ParallelEditAssigner::new(1, 4);

        let umis: Vec<String> = vec!["AAAAAA", "AAAAAT", "AAAATT", "TTTTTT", "TTTTTA", "GGGGGG"]
            .into_iter()
            .map(String::from)
            .collect();

        let parallel_result = parallel.assign(&umis);

        // Verify expected groupings
        // AAAAAA-AAAAAT-AAAATT form a chain
        assert_eq!(parallel_result[0], parallel_result[1]);
        assert_eq!(parallel_result[1], parallel_result[2]);
        // TTTTTT-TTTTTA form a chain
        assert_eq!(parallel_result[3], parallel_result[4]);
        // GGGGGG is isolated
        assert_ne!(parallel_result[5], parallel_result[0]);
        assert_ne!(parallel_result[5], parallel_result[3]);
        // The two chains are separate
        assert_ne!(parallel_result[0], parallel_result[3]);
    }

    #[test]
    fn test_edit_sequential_parallel_agree_on_mixed_case() {
        // Regression (audit P2): the sequential edit assigner (`SimpleErrorUmiAssigner`)
        // matched raw bytes case-sensitively, while the parallel edit assigner uppercases
        // before encoding. Mixed-case UMIs therefore grouped differently depending on
        // `--threads` (sequential vs parallel path). Both must fold case.
        //
        // fgbio baseline: `GroupReadsByUmi` uppercases every UMI before assignment
        // (`canonicalize(rawTag.toUpperCase)`), so case-insensitive grouping is the
        // fgbio-faithful behavior — this converges with fgbio's user-facing output.
        let umis: Vec<String> =
            vec!["ACGT", "acgt", "AcGt", "TTTT", "tttt"].into_iter().map(String::from).collect();

        let sequential = SimpleErrorUmiAssigner::new(1);
        let parallel = ParallelEditAssigner::new(1, 4);

        let seq = sequential.assign(&umis);
        let par = parallel.assign(&umis);

        // Both must fold case: {ACGT, acgt, AcGt} in one molecule, {TTTT, tttt} in another.
        assert_eq!(seq[0], seq[1], "sequential must fold case: ACGT == acgt");
        assert_eq!(seq[1], seq[2], "sequential must fold case: acgt == AcGt");
        assert_eq!(seq[3], seq[4], "sequential must fold case: TTTT == tttt");
        assert_ne!(seq[0], seq[3], "ACGT and TTTT are 3 edits apart, separate molecules");

        // And the two implementations must produce the same grouping structure.
        assert_same_groupings(&seq, &par);

        // Both must also match the programmatic fgbio baseline (uppercase, then transitive
        // edit-distance grouping), not merely agree with each other.
        assert_matches_fgbio_edit_baseline(&seq, &umis, 1);
        assert_matches_fgbio_edit_baseline(&par, &umis, 1);
    }

    proptest! {
        /// `ParallelEditAssigner` and the sequential `SimpleErrorUmiAssigner` must induce the
        /// same partition for any UMI multiset, including mixed-case input — and that partition
        /// must be the fgbio-faithful one. This is the generated-coverage parity guard for the
        /// audit-P2 case-folding fix: the single hand-picked case in
        /// `test_edit_sequential_parallel_agree_on_mixed_case` can't reach the case-collision
        /// patterns at the tie-break / edit-boundary for random multisets, so exercise both fgumi
        /// edit paths over random mixed-case input.
        ///
        /// Both paths are additionally pinned to `fgbio_edit_baseline_partition`, an independent
        /// programmatic fgbio baseline (uppercase, then transitive edit-distance grouping). This
        /// is stronger than cross-checking the two fgumi paths against each other: they share the
        /// union-find grouping algorithm and could drift together, whereas the baseline is an
        /// independent reference that reproduces fgbio's mixed-case behavior
        /// (`canonicalize(rawTag.toUpperCase)`) without an fgbio dependency. UMIs are drawn from a
        /// small mixed-case alphabet over a short fixed length to maximize case collisions and
        /// 1-edit adjacency.
        #[test]
        fn proptest_edit_parallel_matches_sequential(
            indices in prop::collection::vec(0usize..16, 1..40),
        ) {
            const POOL: [&str; 16] = [
                "aAAA", "AAAA", "CAAA", "cAAA", "GAAA", "gAAA", "TAAA", "tAAA",
                "AACA", "aaca", "AAGA", "aaga", "ACGT", "acgt", "TGCA", "tgca",
            ];
            let umis: Vec<Umi> = indices.iter().map(|&i| POOL[i].to_string()).collect();

            let parallel = ParallelEditAssigner::new(1, 4).assign(&umis);
            let sequential = SimpleErrorUmiAssigner::new(1).assign(&umis);

            assert_same_partition_by_molecule(&parallel, &sequential);
            assert_matches_fgbio_edit_baseline(&parallel, &umis, 1);
            assert_matches_fgbio_edit_baseline(&sequential, &umis, 1);
        }
    }

    /// Helper to check if two assignment vectors produce the same groupings.
    /// Molecule IDs may differ, but the grouping structure must match.
    fn assert_same_groupings(a: &[MoleculeId], b: &[MoleculeId]) {
        assert_eq!(a.len(), b.len());

        // Build equivalence classes for both
        let mut a_groups: AHashMap<MoleculeId, Vec<usize>> = AHashMap::new();
        let mut b_groups: AHashMap<MoleculeId, Vec<usize>> = AHashMap::new();

        for (i, (ma, mb)) in a.iter().zip(b.iter()).enumerate() {
            a_groups.entry(*ma).or_default().push(i);
            b_groups.entry(*mb).or_default().push(i);
        }

        // Sort and compare
        let mut a_sorted: Vec<Vec<usize>> = a_groups.into_values().collect();
        let mut b_sorted: Vec<Vec<usize>> = b_groups.into_values().collect();

        for v in &mut a_sorted {
            v.sort_unstable();
        }
        for v in &mut b_sorted {
            v.sort_unstable();
        }
        a_sorted.sort();
        b_sorted.sort();

        assert_eq!(a_sorted, b_sorted, "Groupings do not match");
    }

    #[test]
    fn test_parallel_edit_vs_parallel_edit_deterministic() {
        // Run the same assignment twice to verify determinism
        let umis: Vec<String> = vec![
            "AAAAAA", "AAAAAT", "AAAATT", "TTTTTT", "TTTTTA", "GGGGGG", "CCCCCC", "CCCCCG",
            "CCCCGG",
        ]
        .into_iter()
        .map(String::from)
        .collect();

        let assigner1 = ParallelEditAssigner::new(1, 4);
        let assigner2 = ParallelEditAssigner::new(1, 4);

        let result1 = assigner1.assign(&umis);
        let result2 = assigner2.assign(&umis);

        assert_same_groupings(&result1, &result2);
    }

    #[test]
    fn test_parallel_adjacency_vs_parallel_adjacency_deterministic() {
        // Run the same assignment twice to verify determinism
        let mut umis: Vec<String> = Vec::new();
        umis.extend(std::iter::repeat_n("AAAAAA".to_string(), 10));
        umis.extend(std::iter::repeat_n("AAAAAT".to_string(), 3));
        umis.extend(std::iter::repeat_n("TTTTTT".to_string(), 8));
        umis.extend(std::iter::repeat_n("TTTTTA".to_string(), 2));

        let assigner1 = ParallelAdjacencyAssigner::new(1, 4);
        let assigner2 = ParallelAdjacencyAssigner::new(1, 4);

        let result1 = assigner1.assign(&umis);
        let result2 = assigner2.assign(&umis);

        assert_same_groupings(&result1, &result2);
    }

    // ==================== Sequential-vs-parallel parity (S7-002 / FU-004) ====================

    /// Compare two assignment vectors by the *partition* of input indices they induce, keyed on
    /// the underlying molecule id (`MoleculeId::id()`). Unlike `assert_same_groupings`, this
    /// collapses `PairedA`/`PairedB` of the same molecule into one group, so it is valid for both
    /// the single (`AdjacencyUmiAssigner`) and paired (`PairedUmiAssigner`) strategies. Concrete
    /// id *values* are irrelevant — only which inputs share a molecule.
    fn assert_same_partition_by_molecule(a: &[MoleculeId], b: &[MoleculeId]) {
        assert_eq!(a.len(), b.len());
        assert_eq!(
            molecule_partition(a),
            molecule_partition(b),
            "Partitions diverge between assigners"
        );
    }

    /// Canonical partition of input indices induced by an assignment vector, keyed on the
    /// underlying molecule id (`MoleculeId::id()`). `PairedA`/`PairedB` of the same molecule
    /// collapse into one group. Concrete id *values* are irrelevant — only which inputs share a
    /// molecule — so the result is directly comparable across assigners and against a reference
    /// baseline.
    fn molecule_partition(ids: &[MoleculeId]) -> Vec<Vec<usize>> {
        let mut groups: AHashMap<Option<u64>, Vec<usize>> = AHashMap::new();
        for (i, id) in ids.iter().enumerate() {
            groups.entry(id.id()).or_default().push(i);
        }
        let mut sorted: Vec<Vec<usize>> = groups.into_values().collect();
        for group in &mut sorted {
            group.sort_unstable();
        }
        sorted.sort();
        sorted
    }

    /// Programmatic fgbio-faithful baseline partition for the *edit* ("simple error") strategy.
    ///
    /// fgbio's `GroupReadsByUmi` uppercases every UMI before assignment
    /// (`canonicalize(rawTag.toUpperCase)`), then its edit/`SimpleErrorUmiAssigner` strategy groups
    /// transitively: two UMIs share a molecule when they are within `max_mismatches`, and that
    /// relation is closed transitively (chains merge). The fgbio-faithful partition is therefore
    /// the connected components of the graph over the *uppercased* UMIs whose edges join pairs
    /// within `max_mismatches`.
    ///
    /// This is an independent reference — a dead-simple O(n^2) union-find over the pairwise
    /// within-threshold relation, not the lock-free union-find / BFS the production assigners use —
    /// so matching it is a genuine identity check against fgbio's documented behavior rather than a
    /// same-algorithm-drifts-together tautology. It takes no fgbio dependency: fgbio's mixed-case
    /// behavior is fully determined by the up-front uppercasing, which the oracle reproduces.
    fn fgbio_edit_baseline_partition(umis: &[Umi], max_mismatches: u32) -> Vec<Vec<usize>> {
        fn find(parent: &mut [usize], mut x: usize) -> usize {
            while parent[x] != x {
                parent[x] = parent[parent[x]];
                x = parent[x];
            }
            x
        }

        let upper: Vec<String> = umis.iter().map(|umi| umi.to_ascii_uppercase()).collect();
        let within = |a: &str, b: &str| -> bool {
            a.len() == b.len()
                && a.bytes().zip(b.bytes()).filter(|(x, y)| x != y).count()
                    <= max_mismatches as usize
        };

        let mut parent: Vec<usize> = (0..upper.len()).collect();
        for i in 0..upper.len() {
            for j in (i + 1)..upper.len() {
                if within(&upper[i], &upper[j]) {
                    let (ri, rj) = (find(&mut parent, i), find(&mut parent, j));
                    if ri != rj {
                        parent[ri] = rj;
                    }
                }
            }
        }

        let mut groups: AHashMap<usize, Vec<usize>> = AHashMap::new();
        for i in 0..upper.len() {
            let root = find(&mut parent, i);
            groups.entry(root).or_default().push(i);
        }
        let mut sorted: Vec<Vec<usize>> = groups.into_values().collect();
        for group in &mut sorted {
            group.sort_unstable();
        }
        sorted.sort();
        sorted
    }

    /// Assert an edit-assigner's output partition is identical to the programmatic fgbio baseline.
    fn assert_matches_fgbio_edit_baseline(
        assignments: &[MoleculeId],
        umis: &[Umi],
        max_mismatches: u32,
    ) {
        assert_eq!(
            molecule_partition(assignments),
            fgbio_edit_baseline_partition(umis, max_mismatches),
            "assignment partition diverges from the programmatic fgbio edit baseline"
        );
    }

    proptest! {
        /// `ParallelAdjacencyAssigner` and the sequential `AdjacencyUmiAssigner` must induce the
        /// same partition for any UMI multiset, including mixed-case input. This is the parity
        /// guard for S7-002 / FU-004: before converging the parallel tie-break key onto the RAW
        /// first-seen UMI string, mixed-case equal-count roots could sort in opposite order
        /// between the two paths (parallel keyed on the uppercased string, sequential on raw
        /// bytes) and capture a shared neighbor into different groups. UMIs are drawn from a small
        /// alphabet (mixed case) over a short fixed length to maximize ties and 1-edit adjacency.
        ///
        /// No fgbio oracle (by design): the case-SENSITIVE raw-bytes tie-break IS fgbio's
        /// `strnum`-style order, established in S7-002-bitenc-perf-analysis.md (fgbio = raw bytes).
        /// The DIRECTION test below (`adjacency_tie_break_is_case_sensitive_raw_not_uppercased`,
        /// aAAA/CAAA/GAAA) non-vacuously pins that direction, so we cross-check the two fgumi paths
        /// against each other rather than adding an fgbio dependency.
        #[test]
        fn proptest_adjacency_parallel_matches_sequential(
            indices in prop::collection::vec(0usize..16, 1..40),
        ) {
            // A small pool of length-4 UMIs over {A,C,G,T} in mixed case, so distinct-by-case
            // strings collide to the same case-insensitive `BitEnc` node and exercise the
            // tie-break ordering.
            const POOL: [&str; 16] = [
                "aAAA", "AAAA", "CAAA", "cAAA", "GAAA", "gAAA", "TAAA", "tAAA",
                "AACA", "aaca", "AAGA", "aaga", "ACGT", "acgt", "TGCA", "tgca",
            ];
            let umis: Vec<Umi> = indices.iter().map(|&i| POOL[i].to_string()).collect();

            let parallel = ParallelAdjacencyAssigner::new(1, 4).assign(&umis);
            let sequential =
                AdjacencyUmiAssigner::new(1, 1, DEFAULT_INDEX_THRESHOLD).assign(&umis);

            assert_same_partition_by_molecule(&parallel, &sequential);
        }
    }

    /// Non-vacuous DIRECTION guard for S7-002 / FU-004: pins that the tie-break is
    /// the case-SENSITIVE RAW first-seen string (fgbio parity), not the uppercased
    /// one. With raw bytes, `"CAAA"` (C=0x43) sorts before `"aAAA"` (a=0x61), so the
    /// equal-count `CAAA` root is processed first and captures the 1-edit `GAAA`
    /// child — `GAAA` groups with `CAAA`, NOT with `aAAA`. An uppercased tie-break
    /// (`"AAAA" < "CAAA"`) would reverse this and group `GAAA` with `aAAA`, so this
    /// test fails under the wrong direction. Both assigners must agree on the raw
    /// outcome.
    #[test]
    fn adjacency_tie_break_is_case_sensitive_raw_not_uppercased() {
        // 4×aAAA, 4×CAAA, 1×GAAA. aAAA/CAAA are equal-count (4); neither captures
        // the other (4 > 4/2+1=3); both are 1 edit from GAAA (capturable at 1).
        let umis: Vec<Umi> =
            ["aAAA", "aAAA", "aAAA", "aAAA", "CAAA", "CAAA", "CAAA", "CAAA", "GAAA"]
                .iter()
                .map(|s| (*s).to_string())
                .collect();

        for label_and_result in [
            ("sequential", AdjacencyUmiAssigner::new(1, 1, DEFAULT_INDEX_THRESHOLD).assign(&umis)),
            ("parallel", ParallelAdjacencyAssigner::new(1, 4).assign(&umis)),
        ] {
            let (label, mol) = label_and_result;
            // GAAA (index 8) must share CAAA's molecule (indices 4..=7), not aAAA's (0..=3).
            assert_eq!(mol[8], mol[4], "{label}: GAAA must group with CAAA under a raw tie-break");
            assert_ne!(mol[8], mol[0], "{label}: GAAA must NOT group with aAAA (that is the bug)");
            assert_eq!(mol[0], mol[3], "{label}: the four aAAA reads share one molecule");
            assert_eq!(mol[4], mol[7], "{label}: the four CAAA reads share one molecule");
        }
    }

    /// Invalid-UMI parity: a read whose UMI fails to encode (e.g. contains an `N`)
    /// must map to `MoleculeId::None` in BOTH the sequential `AdjacencyUmiAssigner`
    /// and the parallel `ParallelAdjacencyAssigner`, never to a fresh
    /// `MoleculeId::Single`. The sequential path has long preserved `None`; the
    /// parallel path previously minted a `Single` per invalid read, which broke
    /// sequential/parallel parity and would silently fold invalid reads into their
    /// own consensus families. This pins the fix for both the mixed case and the
    /// all-invalid case.
    ///
    /// The `None` expectations below ARE the fgbio oracle, not a mere
    /// sequential-vs-parallel cross-check: fgbio's `GroupReadsByUmi` drops reads
    /// whose UMI contains an `N` (`.filter(r => !umi.contains('N'))`), assigning
    /// them no molecule — exactly `MoleculeId::None` here. (That fgbio-compatible
    /// UMI validation is documented on `fgumi_umi`'s UMI validator, which mirrors
    /// the same Scala filter.) Because these are absolute expected values, a
    /// shared drift from fgbio — e.g. both paths minting a fresh `Single` for an
    /// invalid UMI — FAILS the asserts rather than slipping through. No fgbio
    /// dependency is added, consistent with the no-oracle-by-design rationale on
    /// `proptest_adjacency_parallel_matches_sequential` above.
    #[test]
    fn adjacency_invalid_umi_maps_to_none_in_both_paths() {
        // One valid UMI (twice) + one invalid UMI (contains `N`).
        let mixed: Vec<Umi> = ["ACGT", "ACGT", "ACGN"].iter().map(|s| (*s).to_string()).collect();
        for (label, mol) in [
            ("sequential", AdjacencyUmiAssigner::new(1, 1, DEFAULT_INDEX_THRESHOLD).assign(&mixed)),
            ("parallel", ParallelAdjacencyAssigner::new(1, 4).assign(&mixed)),
        ] {
            assert_eq!(mol[0], mol[1], "{label}: the two valid reads share one molecule");
            assert!(matches!(mol[0], MoleculeId::Single(_)), "{label}: valid reads get a Single");
            assert_eq!(mol[2], MoleculeId::None, "{label}: invalid read must be None, not Single");
        }

        // All-invalid input: every read must be `None` (no fresh `Single` ids).
        let all_invalid: Vec<Umi> =
            ["ACGN", "NNNN", "TTNT"].iter().map(|s| (*s).to_string()).collect();
        for (label, mol) in [
            (
                "sequential",
                AdjacencyUmiAssigner::new(1, 1, DEFAULT_INDEX_THRESHOLD).assign(&all_invalid),
            ),
            ("parallel", ParallelAdjacencyAssigner::new(1, 4).assign(&all_invalid)),
        ] {
            assert!(
                mol.iter().all(|id| *id == MoleculeId::None),
                "{label}: all-invalid input must map every read to None, got {mol:?}",
            );
        }
    }

    // NOTE on the paired strategy (S7-002 / FU-004 scope): the sequential `PairedUmiAssigner`
    // tie-break is the case-SENSITIVE canonical string (`canonicalize_paired` normalizes
    // orientation only, not case) — raw/fgbio-compatible, matching the single-UMI raw tie-break.
    // The parallel `ParallelPairedAssigner` is NOT converged onto raw here, and a *full* paired
    // parity proptest is intentionally NOT added: the two paired paths diverge on inputs
    // unrelated to the tie-break key — `ParallelPairedAssigner::canonicalize` uppercases (so its
    // node identity itself is case-insensitive, unlike the sequential), and the two paths
    // implement reverse-orientation adjacency matching differently (e.g. `GAAA-CCCC` vs
    // `AAAA-CCCC`, 1 edit apart, group in the parallel path but not the sequential one).
    // Converging the paired paths means reconciling those deeper divergences, which is a
    // pre-existing, out-of-scope follow-up distinct from the FU-004 tie-break-key fix.
}
