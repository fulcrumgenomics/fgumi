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

use fgumi_umi::assigner::assert_uniform_umi_length;

use super::{Umi, UmiAssigner};

/// Assign each distinct invalid (non-encodable) UMI its own molecule, with identical uppercased
/// strings sharing -- mirroring fgbio's `Map[Umi, MoleculeId]` (every UMI it sees is assigned a
/// molecule) and the sequential assigners. Used for the "all UMIs invalid" edge case shared by
/// every parallel assigner.
fn all_invalid_molecule_ids(raw_umis: &[Umi]) -> Vec<MoleculeId> {
    let mut invalid_to_id: AHashMap<String, MoleculeId> = AHashMap::new();
    let mut next_mol_id: u64 = 0;
    raw_umis
        .iter()
        .map(|umi| {
            *invalid_to_id.entry(umi.to_uppercase()).or_insert_with(|| {
                let id = MoleculeId::Single(next_mol_id);
                next_mol_id += 1;
                id
            })
        })
        .collect()
}

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

/// Parallel reverse-orientation edge discovery for the paired assigner.
///
/// Returns undirected edges `(i, j)` (with `i < j`) where the REVERSE of canonical
/// form `i` is within `max_mismatches` of canonical form `j`. This is the half of
/// fgbio's `matches = within(l, r) OR within(reverse(l), r)` that canonicalization
/// alone misses: a single mismatch can flip which half is lexically smaller, so two
/// reads of one molecule canonicalize to forms that are far apart forward yet within
/// threshold when one is reversed (GRP3-01). `reverse_encs[i]` must be the encoding
/// of the halves-swapped `i`th canonical form.
#[must_use]
fn discover_paired_reverse_edges(
    forward: &[(BitEnc, usize)],
    reverse_encs: &[BitEnc],
    max_mismatches: u32,
) -> Vec<(usize, usize)> {
    let umi_to_idx: AHashMap<BitEnc, usize> =
        forward.iter().enumerate().map(|(i, (enc, _))| (*enc, i)).collect();

    reverse_encs
        .par_iter()
        .enumerate()
        .flat_map(|(i, rev_enc)| {
            let mut edges = Vec::new();
            let neighbors: Vec<BitEnc> = if max_mismatches == 1 {
                generate_neighbors(rev_enc).collect()
            } else {
                generate_neighbors_k(rev_enc, max_mismatches)
            };
            for neighbor in neighbors {
                if let Some(&j) = umi_to_idx.get(&neighbor) {
                    // Distinct canonical forms are never exact reverses of each other,
                    // so `i == j` can only arise for a palindromic form; skip it.
                    if i != j && rev_enc.hamming_distance(&neighbor) <= max_mismatches {
                        edges.push((i.min(j), i.max(j)));
                    }
                }
            }
            edges
        })
        .collect()
}

/// Strand (`true` = A) for reads spelled exactly as each canonical form, relative to
/// that molecule's root, mirroring the sequential `PairedUmiAssigner`.
///
/// The root's canonical is strand A; another member's canonical is strand A iff it is
/// strictly closer (fewer mismatches) to the root's canonical than its reverse is.
/// Reads in the reverse spelling take the opposite strand. Without this, reads reached
/// via a reverse-orientation edge (which represent the OPPOSITE strand) would all be
/// mislabeled strand A.
///
/// Palindromic canonical forms (identical halves) are their own reverse, so their
/// forward and reverse spellings are the same string; the sequential assigner inserts
/// the reverse spelling LAST for them. For a palindromic ROOT that reverse spelling is
/// `B` (overwriting the `A` it would otherwise get); for a palindromic non-root child
/// it is `A`.
fn paired_canonical_strands(
    unique_umis: &[(BitEnc, usize)],
    reverse_encs: &[BitEnc],
    mol_ids: &[u64],
    cluster_root: &[usize],
) -> Vec<bool> {
    (0..unique_umis.len())
        .map(|i| {
            let mol_id = usize::try_from(mol_ids[i]).expect("molecule id fits in usize");
            let root_i = cluster_root[mol_id];
            let is_palindrome = unique_umis[i].0 == reverse_encs[i];
            if i == root_i {
                return !is_palindrome; // root: A, unless palindrome -> B
            }
            if is_palindrome {
                return true; // non-root palindrome child -> A
            }
            let root_enc = &unique_umis[root_i].0;
            let forward = root_enc.hamming_distance(&unique_umis[i].0);
            let reversed = root_enc.hamming_distance(&reverse_encs[i]);
            forward < reversed
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
/// Produces identical results to `SimpleErrorUmiAssigner` in `assigner.rs`.
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
            let umi_upper = umi.to_uppercase();
            // Use from_umi_str to handle paired UMIs with dashes (e.g., "ACGT-TGCA")
            if let Some(enc) = BitEnc::from_umi_str(&umi_upper) {
                *umi_counts.entry(enc).or_insert(0) += 1;
                umi_to_original.push(Some(enc));
            } else {
                umi_to_original.push(None);
            }
        }

        // Handle edge case: no valid UMIs. Give each distinct invalid string its own molecule
        // (identical strings share), mirroring fgbio's per-string assignment and the sequential
        // edit assigner.
        if umi_counts.is_empty() {
            return all_invalid_molecule_ids(raw_umis);
        }

        // Convert to indexed list (order doesn't matter for Edit strategy)
        let unique_umis: Vec<(BitEnc, usize)> = umi_counts.into_iter().collect();

        // Reject differing-length UMIs the way fgbio (and the sequential assigner) does.
        assert_uniform_umi_length(unique_umis.iter().map(|(enc, _)| enc.len()));

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
        let mut invalid_to_id: AHashMap<String, MoleculeId> = AHashMap::new();
        let mut next_mol_id: u64 = 0;

        let mut result = Vec::with_capacity(raw_umis.len());
        for (raw, enc_opt) in raw_umis.iter().zip(&umi_to_original) {
            let mol_id = if let Some(enc) = enc_opt {
                let idx = enc_to_idx[enc];
                let root = uf.find(idx);
                *root_to_mol.entry(root).or_insert_with(|| {
                    let id = MoleculeId::Single(next_mol_id);
                    next_mol_id += 1;
                    id
                })
            } else {
                // Each distinct invalid (non-encodable) UMI gets its own molecule (identical
                // strings share), mirroring fgbio's per-string assignment and the sequential
                // edit assigner. Invalid UMIs never join a valid molecule. See the
                // cross-assigner parity note in the tests module.
                *invalid_to_id.entry(raw.to_uppercase()).or_insert_with(|| {
                    let id = MoleculeId::Single(next_mol_id);
                    next_mol_id += 1;
                    id
                })
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

        // Count unique UMIs
        let mut umi_counts: AHashMap<String, usize> = AHashMap::new();
        for umi in raw_umis {
            *umi_counts.entry(umi.to_uppercase()).or_insert(0) += 1;
        }

        // Build sorted list with BitEnc encoding
        // Filter out invalid UMIs during encoding
        // Use from_umi_str to handle paired UMIs with dashes (e.g., "ACGT-TGCA")
        let mut sorted_umis: Vec<(String, usize, BitEnc)> = umi_counts
            .iter()
            .filter_map(|(umi, &count)| {
                BitEnc::from_umi_str(umi).map(|enc| (umi.clone(), count, enc))
            })
            .collect();

        // Handle edge case: no valid UMIs. Give each distinct invalid string its own molecule
        // (identical strings share), mirroring fgbio's per-string assignment and the sequential
        // adjacency assigner.
        if sorted_umis.is_empty() {
            return all_invalid_molecule_ids(raw_umis);
        }

        // Sort by count descending, then by the case-folded UMI string for determinism.
        // `a.0` is the uppercased UMI (the count-map key), so this tie-break is
        // case-insensitive and matches the sequential assigner, which folds case in its
        // own tie-break. Keeping both case-folded is what makes the two implementations
        // group mixed-case inputs identically.
        sorted_umis.sort_by(|a, b| b.1.cmp(&a.1).then_with(|| a.0.cmp(&b.0)));

        // Build indexed structures (avoid cloning strings by using indices)
        let unique_umis: Vec<(BitEnc, usize)> =
            sorted_umis.iter().map(|(_, count, enc)| (*enc, *count)).collect();

        // Reject differing-length UMIs the way fgbio (and the sequential assigner) does.
        assert_uniform_umi_length(unique_umis.iter().map(|(enc, _)| enc.len()));

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
            .map(|(i, (umi, _, _))| (umi.as_str(), mol_ids[i]))
            .collect();

        // Map back to original UMI order. Each distinct invalid (non-encodable) UMI gets its own
        // molecule (identical strings share, continuing the `next_mol_id` sequence), mirroring
        // fgbio's per-string assignment and the sequential adjacency assigner. Invalid UMIs never
        // join a valid molecule. See the cross-assigner parity note in the tests module.
        let mut invalid_to_id: AHashMap<String, MoleculeId> = AHashMap::new();
        raw_umis
            .iter()
            .map(|umi| {
                let upper = umi.to_uppercase();
                str_to_mol.get(upper.as_str()).copied().unwrap_or_else(|| {
                    *invalid_to_id.entry(upper).or_insert_with(|| {
                        let id = MoleculeId::Single(next_mol_id);
                        next_mol_id += 1;
                        id
                    })
                })
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
        umi.to_uppercase() == canonical
    }
}

impl UmiAssigner for ParallelPairedAssigner {
    // One cohesive algorithm (encode → build reverse-aware adjacency → BFS clusters →
    // canonical strands → final per-UMI molecule/strand map); splitting it would hurt
    // readability more than the length costs. Matches other assigner/pipeline sites.
    #[allow(clippy::too_many_lines)]
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

        // Handle edge case: no valid UMIs. Give each distinct invalid string its own molecule
        // (identical strings share), mirroring fgbio's per-string assignment, the main path
        // below, and the sequential paired assigner.
        if sorted_umis.is_empty() {
            return all_invalid_molecule_ids(raw_umis);
        }

        // Sort by count descending, then by UMI string for determinism
        sorted_umis.sort_by(|a, b| b.1.cmp(&a.1).then_with(|| a.0.cmp(&b.0)));

        // Build indexed structures
        let unique_umis: Vec<(BitEnc, usize)> =
            sorted_umis.iter().map(|(_, count, enc)| (*enc, *count)).collect();

        // Reject differing-length paired UMIs the way fgbio (and the sequential assigner)
        // does. `BitEnc::len` excludes the dash, so this compares base length.
        assert_uniform_umi_length(unique_umis.iter().map(|(enc, _)| enc.len()));

        // Guard the pre-existing `BitEnc`-drops-dash limitation. `BitEnc::from_umi_str`
        // discards the `-`, so the parallel path's forward AND reverse edit distances are
        // dash-blind. fgbio (`GroupReadsByUmi.PairedUmiAssigner`) and the sequential
        // `PairedUmiAssigner` instead compare the dash-delimited string byte-for-byte
        // (dash-sensitive, both orientations). The two agree exactly when every UMI has
        // SYMMETRIC halves (both sides the same length): the dash then sits at a fixed
        // position, so dropping it loses no information and swapping halves is a clean
        // reverse. With ASYMMETRIC halves the dash position varies, and the dash-blind path
        // diverges in more than one way -- distinct forms with the same concatenated bases
        // but different splits collide to one encoding (`AC-GTA` / `ACG-TA`), and the
        // halves-swapped reverse encoding can forge a within-threshold reverse edge that the
        // dash-sensitive reference never draws (`A-AC` reverse `ACA` is Hamming-1 from
        // `A-CT`'s `ACT`). A collision-only check misses the latter.
        //
        // A full split-aware rework of the hot path is the complete fix (tracked in #586);
        // until then, delegate any pool that contains an asymmetric-halves UMI to the
        // sequential assigner, which is byte-for-byte faithful to fgbio. Symmetric pools
        // (the common case) take the fast parallel path unchanged. Checked over `sorted_umis`
        // -- the encodable canonical forms that actually feed edge discovery.
        {
            let has_asymmetric_halves = sorted_umis.iter().any(|(umi, _, _)| {
                let (left, right) = umi.split_once('-').expect("validated paired UMI");
                left.len() != right.len()
            });
            if has_asymmetric_halves {
                return crate::umi::PairedUmiAssigner::new(self.max_mismatches).assign(raw_umis);
            }
        }

        // Encode the REVERSE of each canonical form (halves swapped) alongside it, so
        // the reverse-orientation edge pass and the strand assignment can both use it.
        let reverse_encs: Vec<BitEnc> = sorted_umis
            .iter()
            .map(|(umi, _, _)| {
                let reversed = Self::reverse_paired(umi).unwrap_or_else(|| umi.clone());
                BitEnc::from_umi_str(&reversed)
                    .expect("canonical paired UMI reverses to a valid UMI")
            })
            .collect();

        // Phase 1: Parallel edge discovery (using configured thread pool).
        //
        // Two reads map to the same molecule when they are within `max_mismatches`
        // in EITHER orientation, exactly as the sequential `PairedUmiAssigner`'s
        // `matches_paired` tests `within(l, r) OR within(reverse(l), r)`.
        // Canonicalization alone does NOT capture the reverse case: a single mismatch
        // can flip which half is lexically smaller, so two reads of one molecule
        // canonicalize to forms that are far apart forward yet within threshold when
        // one is reversed (GRP3-01). So we union the forward edges with a
        // reverse-orientation edge pass.
        let max_mismatches = self.max_mismatches;
        let edges = self.pool.install(|| {
            let forward = discover_edges_parallel_k(&unique_umis, max_mismatches);
            let reverse =
                discover_paired_reverse_edges(&unique_umis, &reverse_encs, max_mismatches);
            let mut set: AHashSet<(usize, usize)> = forward.into_iter().collect();
            set.extend(reverse);
            set
        });

        // Build adjacency list. `edges` is an `AHashSet` whose iteration order is seeded
        // per process at runtime (unlike the unpaired path, which feeds a deterministically
        // ordered rayon `Vec` straight in), so the order neighbors land in each bucket is
        // nondeterministic. The BFS below visits neighbors in bucket order and count-gates
        // them (`neighbor_count <= max_child_count`), so a threshold-boundary neighbor
        // reachable from more than one node can be claimed differently depending on that
        // order -- a run-to-run nondeterministic partition. Sort each bucket after
        // construction to pin a deterministic neighbor order (ascending index == count
        // descending, then UMI string, matching the sorted order the BFS processes roots in).
        let mut adj_list: Vec<Vec<usize>> = vec![Vec::new(); unique_umis.len()];
        for (i, j) in edges {
            adj_list[i].push(j);
            adj_list[j].push(i);
        }
        for neighbors in &mut adj_list {
            neighbors.sort_unstable();
        }

        // Phase 2: Sequential BFS with adjacency constraints. Track each molecule's
        // root canonical (the highest-count member, first in sorted order) so strand
        // can be assigned relative to it, matching the sequential assigner.
        let mut assigned = vec![false; unique_umis.len()];
        let mut mol_ids: Vec<u64> = vec![0; unique_umis.len()];
        let mut cluster_root: Vec<usize> = Vec::new();
        let mut next_mol_id: u64 = 0;

        for root_idx in 0..unique_umis.len() {
            if assigned[root_idx] {
                continue;
            }

            let mol_id = next_mol_id;
            next_mol_id += 1;
            cluster_root.push(root_idx);

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

        let canonical_strand =
            paired_canonical_strands(&unique_umis, &reverse_encs, &mol_ids, &cluster_root);

        // Final pass: map each raw UMI back to its molecule + strand via `canonical_to_idx`;
        // each distinct invalid (non-encodable) UMI instead gets its own `Single` molecule
        // keyed by its raw uppercase string (no valid strand so never `PairedA`/`PairedB`; see
        // the tests-module parity note). Keying by the raw uppercase string -- not the canonical
        // form -- matches the sequential paired assigner (`assign_with_invalid_fallback`) and the
        // `all_invalid_molecule_ids` branch, which both key by `to_uppercase()`. Two invalid UMIs
        // that are reverses of each other (`ACGN-TTTT` / `TTTT-ACGN`) canonicalize to the same
        // form; keying by canonical would merge them into one molecule here while the sequential
        // and all-invalid paths keep them distinct -- a `--threads`-dependent divergence. fgumi's
        // chosen model is one molecule per distinct raw UMI string (an intentional fgbio
        // divergence, since `BitEnc` cannot represent the invalid base).
        let canonical_to_idx: AHashMap<&str, usize> =
            sorted_umis.iter().enumerate().map(|(i, (umi, _, _))| (umi.as_str(), i)).collect();
        let mut invalid_to_id: AHashMap<String, MoleculeId> = AHashMap::new();
        raw_umis
            .iter()
            .map(|umi| {
                let canonical = Self::canonicalize(umi);
                if let Some(&i) = canonical_to_idx.get(canonical.as_str()) {
                    let base_id = mol_ids[i];
                    // Reads spelled as the canonical form take `canonical_strand[i]`;
                    // reads in the reverse spelling take the opposite strand.
                    let is_canonical_spelling = Self::matches_canonical(umi, &canonical);
                    if canonical_strand[i] == is_canonical_spelling {
                        MoleculeId::PairedA(base_id)
                    } else {
                        MoleculeId::PairedB(base_id)
                    }
                } else {
                    *invalid_to_id.entry(umi.to_uppercase()).or_insert_with(|| {
                        let id = MoleculeId::Single(next_mol_id);
                        next_mol_id += 1;
                        id
                    })
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
    use crate::umi::{AdjacencyUmiAssigner, PairedUmiAssigner, SimpleErrorUmiAssigner};
    use rstest::rstest;

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

    // fgbio's GroupReadsByUmi throws on UMIs of differing length; every parallel
    // mismatch-based assigner must reject them the same way rather than silently
    // under-group (tracker GRP-01), matching the sequential assigners.
    #[rstest]
    #[case::edit(|| Box::new(ParallelEditAssigner::new(1, 2)) as Box<dyn UmiAssigner>, &["AAAA", "AAA"])]
    #[case::adjacency(
        || Box::new(ParallelAdjacencyAssigner::new(1, 2)) as Box<dyn UmiAssigner>,
        &["AAAA", "AAA"]
    )]
    #[case::paired(
        || Box::new(ParallelPairedAssigner::new(1, 2)) as Box<dyn UmiAssigner>,
        &["ACT-ACT", "ACT-AC"]
    )]
    #[should_panic(expected = "Multiple UMI lengths")]
    fn test_parallel_assigner_rejects_differing_umi_lengths(
        #[case] make: fn() -> Box<dyn UmiAssigner>,
        #[case] umis: &[&str],
    ) {
        let assigner = make();
        let umis: Vec<Umi> = umis.iter().map(|s| (*s).to_string()).collect();
        let _ = assigner.assign(&umis);
    }

    // Complement to `test_parallel_assigner_rejects_differing_umi_lengths`: a non-encodable
    // (invalid-base) UMI whose length differs from the valid UMIs must NOT trip the
    // differing-length guard. The parallel assigners never feed non-encodable UMIs into the
    // length check (they are dropped during graph construction), so the mixed valid/invalid
    // GRP-01 parity scenario groups the valid UMIs together and gives the invalid one its own
    // molecule instead of panicking. This mirrors the sequential
    // `test_sequential_assigner_excludes_invalid_umis_from_length_guard` and guards against a
    // future refactor moving the guard ahead of the encoding filter and silently reintroducing
    // the parity break. Cross-assigner grouping parity for invalid UMIs is enforced separately
    // by `test_sequential_and_parallel_assigners_induce_same_partition`.
    //
    // fgbio divergence here is *intentional*, so these assertions encode fgumi's chosen behavior
    // rather than fgbio parity: fgbio never drops non-encodable UMIs, so on this exact mixed input
    // its assigners reject the whole family instead of grouping the valid subset. The Edit strategy
    // (SimpleErrorUmiAssigner) throws via Sequences.countMismatches's require(s1.length == s2.length)
    // (fgbio Sequences.scala:99), and Paired (PairedUmiAssigner extends AdjacencyUmiAssigner) throws
    // via the umiLength require (fgbio GroupReadsByUmi.scala:283) -- both see the length-3 NNN /
    // length-5 NN-NN beside the length-4 AAAA / length-7 ACT-ACT. fgumi deliberately tolerates
    // invalid-base UMIs of a differing length by excluding them from the guard (GRP-01) and grouping
    // only the encodable UMIs. Differing-length *valid* UMIs still match fgbio by rejecting (see
    // test_parallel_assigner_rejects_differing_umi_lengths).
    #[rstest]
    #[case::edit(
        || Box::new(ParallelEditAssigner::new(1, 2)) as Box<dyn UmiAssigner>,
        &["AAAA", "AAAA", "NNN"]
    )]
    #[case::paired(
        || Box::new(ParallelPairedAssigner::new(1, 2)) as Box<dyn UmiAssigner>,
        &["ACT-ACT", "ACT-ACT", "NN-NN"]
    )]
    fn test_parallel_assigner_excludes_invalid_umis_from_length_guard(
        #[case] make: fn() -> Box<dyn UmiAssigner>,
        #[case] umis: &[&str],
    ) {
        let assigner = make();
        let umis: Vec<Umi> = umis.iter().map(|s| (*s).to_string()).collect();
        // Must not panic despite the invalid UMI having a different length.
        let ids = assigner.assign(&umis);
        assert_eq!(ids.len(), umis.len());
        assert_eq!(ids[0], ids[1], "identical valid UMIs must share a molecule");
        assert_ne!(ids[0], ids[2], "invalid-base UMI must not join the valid molecule");
        assert_ne!(
            ids[2],
            MoleculeId::None,
            "invalid-base UMI gets its own molecule (fgbio assigns every UMI a molecule)"
        );
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

    // Run on BOTH the sequential and parallel edit assigners: an invalid (non-encodable) UMI must
    // get its own molecule -- distinct from any valid molecule -- by each, while the surrounding
    // valid UMIs still group. This mirrors fgbio's per-string molecule assignment (every UMI it
    // sees gets a molecule). `ACGN` is a same-length 1-mismatch neighbour of `ACGT`; the
    // sequential assigner used to string-match and *group* it (fgbio does group N-neighbours,
    // but the parallel BitEnc path structurally cannot), so this pins the corrected,
    // implementation-independent behavior: `ACGN` is isolated, not merged into `ACGT`.
    #[rstest]
    #[case::sequential(|| Box::new(SimpleErrorUmiAssigner::new(1)) as Box<dyn UmiAssigner>)]
    #[case::parallel(|| Box::new(ParallelEditAssigner::new(1, 2)) as Box<dyn UmiAssigner>)]
    fn test_edit_invalid_umis_get_own_molecule(#[case] make: fn() -> Box<dyn UmiAssigner>) {
        let assigner = make();
        let umis: Vec<String> = vec!["ACGT", "ACGN", "ACGT"] // ACGN is invalid (non-encodable)
            .into_iter()
            .map(String::from)
            .collect();

        let assignments = assigner.assign(&umis);

        // Valid UMIs group together.
        assert_eq!(assignments[0], assignments[2], "identical valid UMIs must share a molecule");
        // The invalid UMI gets its own molecule, distinct from the valid one (never merged into
        // it), and it is a real assignment rather than left unassigned.
        assert_ne!(assignments[1], assignments[0], "invalid UMI must not join the valid molecule");
        assert_ne!(
            assignments[1],
            MoleculeId::None,
            "invalid UMI gets its own molecule (fgbio assigns every UMI a molecule)"
        );
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

    /// Returns true iff two assignment vectors induce the same partition of input
    /// indices (i.e. identical grouping structure), ignoring the concrete
    /// `MoleculeId` labels. Strand-sensitive: reads on opposite strands of the same
    /// molecule are treated as different groups (their `MoleculeId`s differ).
    fn same_partition(a: &[MoleculeId], b: &[MoleculeId]) -> bool {
        a.len() == b.len()
            && (0..a.len()).all(|i| (0..a.len()).all(|j| (a[i] == a[j]) == (b[i] == b[j])))
    }

    /// Returns true iff two assignment vectors induce the same *base-molecule*
    /// partition (ignoring the `/A` `/B` strand suffix). This is what detects
    /// molecule over-splitting: two reads of one molecule on opposite strands share
    /// a base id but not a full `MoleculeId`, so [`same_partition`] alone misses a
    /// split that lands them in different base molecules (GRP3-01).
    fn same_base_partition(a: &[MoleculeId], b: &[MoleculeId]) -> bool {
        let base_same = |xs: &[MoleculeId], i: usize, j: usize| {
            xs[i].base_id_string() == xs[j].base_id_string()
        };
        a.len() == b.len()
            && (0..a.len()).all(|i| (0..a.len()).all(|j| base_same(a, i, j) == base_same(b, i, j)))
    }

    /// Full parallel-vs-sequential equivalence: identical base-molecule partition,
    /// identical strand partition, AND identical ABSOLUTE strand orientation
    /// (`PairedA` vs `PairedB` vs `Single`) at every read position.
    ///
    /// The base-molecule and strand partitions are compared up to molecule-id
    /// relabeling, but the strand is pinned absolutely. `same_partition` alone still
    /// passes under a per-molecule A<->B swap, so an inversion in the root-relative
    /// strand logic (palindrome roots/children especially) would slip through every
    /// caller. Requiring the [`MoleculeId`] enum discriminants to match
    /// position-by-position closes that gap: any caller asserting
    /// `assignments_equivalent(sequential, parallel)` is now asserting the parallel
    /// assigner reproduces the sequential assigner's exact `/A` `/B` labeling, not merely
    /// a partition that agrees up to strand relabeling.
    fn assignments_equivalent(a: &[MoleculeId], b: &[MoleculeId]) -> bool {
        same_base_partition(a, b)
            && same_partition(a, b)
            && a.len() == b.len()
            && a.iter()
                .zip(b)
                .all(|(left, right)| std::mem::discriminant(left) == std::mem::discriminant(right))
    }

    // ==================== Mixed-case parity (sequential vs parallel) ====================
    //
    // fgbio baseline for the three parity tests below. fgbio's `GroupReadsByUmi` CLI
    // uppercases every UMI before assignment (`canonicalize(rawTag.toUpperCase)`), so at the
    // user-facing level fgbio groups mixed-case input case-insensitively. fgumi folds case
    // inside the assigner instead, producing the SAME grouping as fgbio's CLI — these tests
    // pin the concrete grouping that fgbio's CLI emits, so they are a real fgbio-derived
    // expected output, not just a sequential-vs-parallel cross-check.
    //
    // The expected groupings below were captured by driving fgbio's actual assigner classes
    // (fgbio @ src/main/scala/com/fulcrumgenomics/umi/GroupReadsByUmi.scala) on these exact
    // inputs, uppercased as the CLI does:
    //   * edit  [ACGT, acgt]                        -> 1 molecule {ACGT, acgt}
    //   * adjacency [AAAaAA x2, AAAGAC x2, AAAAAC]  -> {AAAAAA(<-AAAaAA), AAAAAC}, {AAAGAC}
    //   * paired [ACGT-TGCA, acgt-tgca, TGCA-ACGT]  -> one molecule, fwd spellings on strand
    //                                                  A, reverse spelling on strand B
    // On UPPERCASE input case folding is a no-op, so fgumi == fgbio there trivially too. The
    // only place fgumi differs from fgbio is fgbio's raw assigner *classes* (case-sensitive
    // counting/tie-break), which fgbio never reaches with mixed case because its CLI
    // pre-uppercases; fgumi moving the fold into the assigner is the safer, equivalent choice.

    #[test]
    fn test_sequential_and_parallel_adjacency_agree_on_mixed_case() {
        // The sequential and parallel adjacency assigners are documented to produce
        // identical groupings. Counting is case-insensitive in both (BitEnc folds case),
        // so their equal-count tie-break must also be case-insensitive to agree.
        //
        // This is the same topology as `test_parallel_adjacency_deterministic_tiebreak`
        // (AAAAAA & AAAGAC both count=2, AAAAAC count=1 within 1 edit of both), but the
        // first node's reads are spelled "AAAaAA". Raw-string ordering sorts "AAAaAA"
        // AFTER "AAAGAC" ('a' > 'G'), while case-folded ordering sorts "AAAAAA" BEFORE
        // "AAAGAC". A raw tie-break would make AAAGAC the root (capturing AAAAAC); a
        // case-folded tie-break makes AAAAAA the root. The two assigners must agree.
        let umis: Vec<String> = vec!["AAAaAA", "AAAaAA", "AAAGAC", "AAAGAC", "AAAAAC"]
            .into_iter()
            .map(String::from)
            .collect();

        let sequential = crate::umi::AdjacencyUmiAssigner::new(1, 1, 100).assign(&umis);
        let parallel = ParallelAdjacencyAssigner::new(1, 2).assign(&umis);

        assert!(
            same_partition(&sequential, &parallel),
            "sequential and parallel adjacency must group mixed-case UMIs identically;\n  \
             sequential: {sequential:?}\n  parallel:   {parallel:?}"
        );

        // Pin the full fgbio-CLI grouping: {AAAAAA(<-AAAaAA), AAAAAC} and {AAAGAC} separate.
        // The case-folded root (AAAAAA, from "AAAaAA") captures AAAAAC, and AAAGAC stays in
        // its own molecule, in both implementations.
        assert_eq!(sequential[0], sequential[4], "sequential: AAAaAA should capture AAAAAC");
        assert_eq!(parallel[0], parallel[4], "parallel: AAAaAA should capture AAAAAC");
        assert_ne!(sequential[0], sequential[2], "sequential: AAAGAC is a separate molecule");
        assert_ne!(parallel[0], parallel[2], "parallel: AAAGAC is a separate molecule");
    }

    #[test]
    fn test_sequential_and_parallel_edit_agree_on_mixed_case() {
        // The parallel edit assigner uppercases before matching, so "ACGT" and "acgt"
        // are the same molecule. The sequential edit assigner must agree: matching is
        // case-insensitive. Without case folding, "ACGT" and "acgt" differ at all four
        // positions (case-sensitive), so the sequential assigner would split them.
        let umis: Vec<String> = vec!["ACGT", "acgt"].into_iter().map(String::from).collect();

        let sequential = crate::umi::SimpleErrorUmiAssigner::new(1).assign(&umis);
        let parallel = ParallelEditAssigner::new(1, 2).assign(&umis);

        assert!(
            same_partition(&sequential, &parallel),
            "sequential and parallel edit must group mixed-case UMIs identically;\n  \
             sequential: {sequential:?}\n  parallel:   {parallel:?}"
        );
        assert_eq!(sequential[0], sequential[1], "ACGT and acgt are the same molecule");
    }

    #[test]
    fn test_sequential_and_parallel_paired_agree_on_mixed_case() {
        // The parallel paired assigner canonicalizes to uppercase, so forward,
        // reversed, and lowercase spellings of the same molecule group together with
        // consistent strand (A/B) assignment. The sequential paired assigner must agree.
        // Without case folding it counts/matches on the raw string, splitting "acgt-tgca"
        // from "ACGT-TGCA".
        //
        // All three spell the same molecule: "ACGT-TGCA" (fwd), "acgt-tgca" (fwd, lower),
        // "TGCA-ACGT" (reverse strand). Expected: indices 0,1 share a strand, index 2 the
        // opposite strand.
        let umis: Vec<String> =
            vec!["ACGT-TGCA", "acgt-tgca", "TGCA-ACGT"].into_iter().map(String::from).collect();

        let sequential = crate::umi::PairedUmiAssigner::new(1).assign(&umis);
        let parallel = ParallelPairedAssigner::new(1, 2).assign(&umis);

        assert!(
            same_partition(&sequential, &parallel),
            "sequential and parallel paired must group mixed-case UMIs identically \
             (including strand);\n  sequential: {sequential:?}\n  parallel:   {parallel:?}"
        );
        // 0 and 1 (both forward) share a strand; 2 (reverse) is the opposite strand of the
        // SAME molecule. Assert the full expected grouping (not just sequential/parallel
        // parity): all three indices must share one base molecule ID, with index 2 on the
        // opposite A/B strand. `same_partition` alone would still pass if index 2 were
        // assigned to an unrelated molecule, so pin the base-molecule and strand
        // relationships explicitly for both implementations.
        //
        // Pin the ABSOLUTE A/B orientation, not just "same"/"opposite": the canonical form of
        // this molecule is "ACGT-TGCA" (the lexicographically smaller of ACGT-TGCA / TGCA-ACGT,
        // uppercased), so the forward spellings (0, 1) match the canonical and land on strand A
        // while the reverse spelling (2) lands on strand B. A global A/B inversion would still
        // satisfy the same/opposite checks below but would diverge from this fgbio-derived
        // expectation, so assert the concrete PairedA/PairedB variants.
        assert!(
            matches!(sequential[0], MoleculeId::PairedA(_))
                && matches!(sequential[1], MoleculeId::PairedA(_))
                && matches!(sequential[2], MoleculeId::PairedB(_)),
            "sequential: expected fwd,fwd,rev => A,A,B; got {sequential:?}"
        );
        assert!(
            matches!(parallel[0], MoleculeId::PairedA(_))
                && matches!(parallel[1], MoleculeId::PairedA(_))
                && matches!(parallel[2], MoleculeId::PairedB(_)),
            "parallel: expected fwd,fwd,rev => A,A,B; got {parallel:?}"
        );
        assert_eq!(sequential[0], sequential[1], "forward spellings share a strand");
        assert_eq!(parallel[0], parallel[1], "parallel: forward spellings share a strand");
        assert_eq!(
            sequential[0].base_id_string(),
            sequential[2].base_id_string(),
            "sequential: reverse spelling should share the base molecule"
        );
        assert_eq!(
            parallel[0].base_id_string(),
            parallel[2].base_id_string(),
            "parallel: reverse spelling should share the base molecule"
        );
        assert_ne!(sequential[0], sequential[2], "reverse spelling is the opposite strand");
        assert_ne!(parallel[0], parallel[2], "parallel: reverse spelling is the opposite strand");
    }

    // ==================== Cross-assigner parity (sequential vs parallel) ====================
    //
    // Production selects between the sequential and parallel assigners purely on the worker
    // thread count (`create_umi_assigner` in `commands/group.rs`: `use_parallel && threads > 1`
    // picks the parallel struct, otherwise the sequential one). The two MUST induce the same
    // partition of the input UMIs for every strategy, or the same reads group differently
    // depending only on `--threads` -- a silent, user-visible correctness bug.
    //
    // Invalid (non-ACGT-encodable, e.g. `N`-containing) UMIs are the sharp edge. fgbio's
    // assigner returns a `Map[Umi, MoleculeId]` keyed by the UMI string, so *every* UMI it sees
    // gets a molecule: identical strings always share, distinct strings never merge. fgumi
    // mirrors that model: an invalid UMI forms its own molecule keyed by its uppercased string --
    // distinct from every valid molecule and from every other distinct invalid string, with
    // identical invalid strings sharing. The one thing fgumi cannot mirror is fgbio grouping a
    // near (within-threshold) `N`-neighbour -- fgbio treats `N` as a mismatch character, but
    // `BitEnc` cannot represent `N`, so fgumi isolates such a UMI instead. Both fgumi assigners
    // agree on that isolation, which is what this harness pins.
    //
    // These inputs are rare but reachable: fgumi (like fgbio, GroupReadsByUmi.scala:673) filters
    // only *uppercase* `N` UMIs upstream (`validate_umi` -> `discarded_ns_in_umi`), so a UMI with
    // a lowercase `n` or a non-ACGT IUPAC base passes the filter (as it does in fgbio) and reaches
    // `assign()` as non-encodable. The two implementations must therefore agree on it.
    //
    // `same_partition` compares grouping structure only, ignoring concrete `MoleculeId`
    // labels, so it is robust to the two assigners numbering molecules in different orders.
    #[rstest]
    // --- Edit strategy ---
    #[case::edit_valid_grouping(
        || Box::new(SimpleErrorUmiAssigner::new(1)) as Box<dyn UmiAssigner>,
        || Box::new(ParallelEditAssigner::new(1, 2)) as Box<dyn UmiAssigner>,
        &["ACGT", "ACGT", "ACGA", "TTTT"]
    )]
    #[case::edit_duplicate_invalid(
        || Box::new(SimpleErrorUmiAssigner::new(1)) as Box<dyn UmiAssigner>,
        || Box::new(ParallelEditAssigner::new(1, 2)) as Box<dyn UmiAssigner>,
        &["AAAA", "NNN", "NNN"]
    )]
    #[case::edit_distinct_invalid(
        || Box::new(SimpleErrorUmiAssigner::new(1)) as Box<dyn UmiAssigner>,
        || Box::new(ParallelEditAssigner::new(1, 2)) as Box<dyn UmiAssigner>,
        &["AAAA", "NNNN", "XXXX"]
    )]
    #[case::edit_invalid_within_threshold(
        || Box::new(SimpleErrorUmiAssigner::new(1)) as Box<dyn UmiAssigner>,
        || Box::new(ParallelEditAssigner::new(1, 2)) as Box<dyn UmiAssigner>,
        &["ACGT", "ACGN", "ACGT"]
    )]
    // --- Adjacency strategy ---
    #[case::adjacency_valid_grouping(
        || Box::new(AdjacencyUmiAssigner::new(1, 1, 100)) as Box<dyn UmiAssigner>,
        || Box::new(ParallelAdjacencyAssigner::new(1, 2)) as Box<dyn UmiAssigner>,
        &["AAAA", "AAAA", "AAAT"]
    )]
    #[case::adjacency_single_invalid(
        || Box::new(AdjacencyUmiAssigner::new(1, 1, 100)) as Box<dyn UmiAssigner>,
        || Box::new(ParallelAdjacencyAssigner::new(1, 2)) as Box<dyn UmiAssigner>,
        &["AAAA", "NNN"]
    )]
    #[case::adjacency_duplicate_invalid_fast_path(
        || Box::new(AdjacencyUmiAssigner::new(1, 1, 100)) as Box<dyn UmiAssigner>,
        || Box::new(ParallelAdjacencyAssigner::new(1, 2)) as Box<dyn UmiAssigner>,
        &["AAAA", "AAAA", "NNN", "NNN"]
    )]
    #[case::adjacency_duplicate_invalid_normal_path(
        || Box::new(AdjacencyUmiAssigner::new(1, 1, 100)) as Box<dyn UmiAssigner>,
        || Box::new(ParallelAdjacencyAssigner::new(1, 2)) as Box<dyn UmiAssigner>,
        &["AAAA", "GGGG", "NNN", "NNN"]
    )]
    // --- Paired strategy ---
    #[case::paired_valid_grouping(
        || Box::new(PairedUmiAssigner::new(1)) as Box<dyn UmiAssigner>,
        || Box::new(ParallelPairedAssigner::new(1, 2)) as Box<dyn UmiAssigner>,
        &["ACT-ACT", "ACT-ACT", "ACT-TGA"]
    )]
    #[case::paired_duplicate_invalid(
        || Box::new(PairedUmiAssigner::new(1)) as Box<dyn UmiAssigner>,
        || Box::new(ParallelPairedAssigner::new(1, 2)) as Box<dyn UmiAssigner>,
        &["ACT-ACT", "NN-NN", "NN-NN"]
    )]
    // Same-length invalid neighbour: `ACGN-TTTT` is a 1-mismatch neighbour of the valid
    // `ACGT-TTTT` (and, unlike `NN-NN` vs `ACT-ACT` above, the SAME length, so it can string-
    // match). The parallel paired assigner drops it before grouping (not BitEnc-encodable) and
    // gives it its own `Single`; the sequential paired assigner must agree. Before the paired-path
    // GRP-01 fix the sequential assigner counted invalid UMIs into `umi_counts`, so this invalid
    // neighbour string-matched into the valid molecule (index 1 joined index 0) -- diverging from
    // the parallel path depending only on `--threads`. With one valid UMI it exercises the
    // `umi_counts.len() == 1` fast path.
    #[case::paired_invalid_neighbor_fast_path(
        || Box::new(PairedUmiAssigner::new(1)) as Box<dyn UmiAssigner>,
        || Box::new(ParallelPairedAssigner::new(1, 2)) as Box<dyn UmiAssigner>,
        &["ACGT-TTTT", "ACGN-TTTT"]
    )]
    // As above but with two distinct valid molecules, exercising the adjacency-graph path
    // (`umi_counts.len() > 1`) rather than the fast path. `ACGN-TTTT` is a 1-mismatch neighbour
    // of `ACGT-TTTT`; it must still be isolated as its own `Single` in both assigners while the
    // two valid UMIs (>1 mismatch apart) stay in separate molecules.
    #[case::paired_invalid_neighbor_graph_path(
        || Box::new(PairedUmiAssigner::new(1)) as Box<dyn UmiAssigner>,
        || Box::new(ParallelPairedAssigner::new(1, 2)) as Box<dyn UmiAssigner>,
        &["ACGT-TTTT", "GGGG-TTTT", "ACGN-TTTT"]
    )]
    // Reversed invalid UMIs beside a valid one: `ACGN-TTTT` and `TTTT-ACGN` are reverses of each
    // other and both non-encodable. fgumi keys each distinct invalid UMI by its raw uppercase
    // string (not its canonical form), so the two reverses land in DISTINCT molecules -- matching
    // the sequential assigner and the all-invalid branch. Before the fix the parallel main-path
    // fallback keyed by canonical form, merging the two reverses into one molecule whenever a
    // valid UMI was also present (so the all-invalid fast path was not taken) -- a
    // `--threads`-dependent divergence. A valid `ACGT-ACGT` keeps that main path live.
    #[case::paired_reversed_invalid(
        || Box::new(PairedUmiAssigner::new(1)) as Box<dyn UmiAssigner>,
        || Box::new(ParallelPairedAssigner::new(1, 2)) as Box<dyn UmiAssigner>,
        &["ACGT-ACGT", "ACGN-TTTT", "TTTT-ACGN"]
    )]
    fn test_sequential_and_parallel_assigners_induce_same_partition(
        #[case] make_sequential: fn() -> Box<dyn UmiAssigner>,
        #[case] make_parallel: fn() -> Box<dyn UmiAssigner>,
        #[case] umis: &[&str],
    ) {
        let umis: Vec<Umi> = umis.iter().map(|s| (*s).to_string()).collect();
        let sequential = make_sequential().assign(&umis);
        let parallel = make_parallel().assign(&umis);

        assert_eq!(sequential.len(), umis.len());
        assert_eq!(parallel.len(), umis.len());
        assert!(
            same_partition(&sequential, &parallel),
            "sequential and parallel assigners must induce the same partition;\n  \
             umis:       {umis:?}\n  sequential: {sequential:?}\n  parallel:   {parallel:?}"
        );
    }

    /// GRP3-01: canonicalization is NOT sufficient to capture reverse-orientation
    /// adjacency. A single mismatch can flip which half is lexically smaller, so two
    /// reads of the same molecule land in canonical forms that are far apart FORWARD
    /// but 1 apart REVERSED. "CAAA-GTTT" and "ATTT-CAAA" are the same molecule (the
    /// reverse of "CAAA-GTTT" is "GTTT-CAAA", which is Hamming-1 from "ATTT-CAAA"),
    /// but their canonical forms "CAAA-GTTT" and "ATTT-CAAA" differ at every base
    /// forward. The sequential `PairedUmiAssigner` groups them via its reverse check;
    /// the parallel assigner must too, or it over-splits the molecule.
    #[test]
    fn test_parallel_paired_groups_reverse_orientation_edge() {
        let umis: Vec<String> =
            vec!["CAAA-GTTT", "ATTT-CAAA"].into_iter().map(String::from).collect();

        let sequential = crate::umi::PairedUmiAssigner::new(1).assign(&umis);
        let parallel = ParallelPairedAssigner::new(1, 2).assign(&umis);

        // Sanity: the sequential assigner groups them into one molecule.
        assert_eq!(
            sequential[0].base_id_string(),
            sequential[1].base_id_string(),
            "sequential should group the reverse-orientation pair; got {sequential:?}"
        );
        // The parallel assigner must match the sequential base-molecule partition
        // (they are one molecule) AND the strand partition (opposite strands).
        assert!(
            assignments_equivalent(&sequential, &parallel),
            "parallel paired must group the reverse-orientation pair like sequential;\n  \
             sequential: {sequential:?}\n  parallel:   {parallel:?}"
        );
    }

    /// A paired UMI's `-` split point is dropped by `BitEnc::from_umi_str`, so two
    /// DISTINCT canonical forms with the same concatenated bases but different split
    /// points encode to the same `BitEnc`. This is only possible with ASYMMETRIC halves
    /// (read-1 UMI length != read-2 UMI length); with symmetric halves the dash sits at a
    /// fixed position, so distinct forms can never share concatenated bases.
    ///
    /// fgbio (`GroupReadsByUmi.PairedUmiAssigner`) and fgumi's sequential
    /// `PairedUmiAssigner` compare the dash-delimited string byte-for-byte (dash-sensitive,
    /// both orientations), so they keep such forms in separate molecules. The parallel
    /// path's `BitEnc`-keyed edge discovery cannot see the dash: an error-neighbour of the
    /// shared encoding bridges the two forms in `BitEnc` space, over-merging (or, via the
    /// `umi_to_idx` collision losing an index, mis-partitioning) versus the sequential
    /// reference -- a `--threads`-dependent divergence.
    ///
    /// Here `AC-GTT` (x4) is a dash-sensitive 1-neighbour of `AC-GTA` but NOT of `ACG-TA`
    /// (their dashes are at different positions -> >=2 mismatches). All three encode with
    /// `AC-GTA` and `ACG-TA` colliding on `ACGTA`. The sequential/fgbio partition is
    /// `{AC-GTT, AC-GTA}` + `{ACG-TA}`; the un-guarded parallel path instead isolates
    /// `AC-GTA` and groups `ACG-TA` with `AC-GTT`. The collision guard detects the shared
    /// encoding and delegates the pool to the sequential assigner, restoring parity.
    #[test]
    fn test_parallel_paired_mixed_half_length_collision_matches_sequential() {
        let umis: Vec<String> = vec!["AC-GTT", "AC-GTT", "AC-GTT", "AC-GTT", "AC-GTA", "ACG-TA"]
            .into_iter()
            .map(String::from)
            .collect();

        let sequential = crate::umi::PairedUmiAssigner::new(1).assign(&umis);
        let parallel = ParallelPairedAssigner::new(1, 2).assign(&umis);

        // Sanity: the sequential (fgbio-faithful) partition keeps `ACG-TA` (index 5)
        // separate from the `AC-GTT`/`AC-GTA` molecule (indices 0-4).
        assert_ne!(
            sequential[4].base_id_string(),
            sequential[5].base_id_string(),
            "sequential should keep the different-split-point form separate; got {sequential:?}"
        );
        assert_eq!(
            sequential[0].base_id_string(),
            sequential[4].base_id_string(),
            "sequential should group AC-GTT with its 1-neighbour AC-GTA; got {sequential:?}"
        );

        // The parallel assigner must induce the same base-molecule partition despite the
        // `BitEnc` collision between `AC-GTA` and `ACG-TA`.
        assert!(
            same_partition(&sequential, &parallel),
            "parallel paired must match the sequential partition on mixed-half-length \
             collisions;\n  sequential: {sequential:?}\n  parallel:   {parallel:?}"
        );
    }

    /// Asymmetric halves can diverge WITHOUT a forward-encoding collision. Dropping the
    /// `-` moves the split, so the halves-swapped reverse encoding used by
    /// [`discover_paired_reverse_edges`] no longer corresponds to a dash-sensitive reverse:
    /// it can spuriously fall within threshold of another form's forward encoding, forging a
    /// reverse edge that fgbio and the sequential `PairedUmiAssigner` (both dash-sensitive in
    /// each orientation) never draw.
    ///
    /// At `max_mismatches = 1`, `A-AC` (x2) and `A-CT` are distinct molecules sequentially:
    /// `matches_paired` compares `A-AC`/`A-CT` (2 mismatches) and `AC-A`/`A-CT` (2
    /// mismatches), both above threshold. But their forward encodings `AAC` and `ACT` are
    /// distinct (no collision to catch), while the reverse encoding of `A-AC` is `ACA`, which
    /// is Hamming-1 from `ACT` -> the parallel reverse pass forges an `A-AC`--`A-CT` edge and
    /// over-merges. A collision-only guard misses this; delegating every asymmetric-halves
    /// pool to the sequential assigner covers it.
    #[test]
    fn test_parallel_paired_asymmetric_noncollision_matches_sequential() {
        let umis: Vec<String> =
            vec!["A-AC", "A-AC", "A-CT"].into_iter().map(String::from).collect();

        let sequential = crate::umi::PairedUmiAssigner::new(1).assign(&umis);
        let parallel = ParallelPairedAssigner::new(1, 2).assign(&umis);

        // Sanity: sequentially `A-CT` (index 2) is its own molecule, distinct from the
        // `A-AC` pair (indices 0-1) -- no forward-encoding collision is involved.
        assert_ne!(
            sequential[0].base_id_string(),
            sequential[2].base_id_string(),
            "sequential should keep A-AC and A-CT separate; got {sequential:?}"
        );

        assert!(
            assignments_equivalent(&sequential, &parallel),
            "parallel paired must match sequential on asymmetric-halves pools even without \
             a forward-encoding collision;\n  sequential: {sequential:?}\n  parallel:   {parallel:?}"
        );
    }

    /// Assert `actual` matches a hand-derived fgbio oracle.
    ///
    /// `expected[i] = (group, strand)`, where `group` is an arbitrary label shared by
    /// every read fgbio places in one *base* molecule and `strand` is fgbio's absolute
    /// duplex suffix (`'A'` for `/A` / [`MoleculeId::PairedA`], `'B'` for `/B` /
    /// [`MoleculeId::PairedB`]). This checks the base-molecule partition and each read's
    /// absolute strand, but not the concrete numeric molecule id (a relabeling).
    ///
    /// The oracle values are derived by hand from fgbio's `PairedUmiAssigner`
    /// (`GroupReadsByUmi.assignIdsToNodes`): the canonical (lexically-smaller-half-first)
    /// spelling of a cluster's root maps to `/A` and its reverse to `/B`; a descendant
    /// takes `/A` for the orientation closer to the root and `/B` for the other; and
    /// fgbio's `Map` last-write-wins gives the palindrome quirks (a palindrome root
    /// canonical spelling resolves to `/B`, a palindrome descendant to `/A`). Because the
    /// oracle is independent of BOTH fgumi assigners, a shared sequential/parallel defect
    /// cannot satisfy it.
    fn assert_matches_fgbio_oracle(actual: &[MoleculeId], expected: &[(u32, char)]) {
        assert_eq!(
            actual.len(),
            expected.len(),
            "length mismatch: actual={actual:?} expected={expected:?}"
        );
        // Base-molecule partition: reads i and j share a base id IFF the oracle groups match.
        for i in 0..actual.len() {
            for j in 0..actual.len() {
                let actual_same = actual[i].base_id_string() == actual[j].base_id_string();
                let expected_same = expected[i].0 == expected[j].0;
                assert_eq!(
                    actual_same, expected_same,
                    "base-partition mismatch at ({i},{j}); actual={actual:?} expected={expected:?}"
                );
            }
        }
        // Absolute strand per read.
        for (i, (id, &(_, strand))) in actual.iter().zip(expected).enumerate() {
            let actual_strand = match id {
                MoleculeId::PairedA(_) => 'A',
                MoleculeId::PairedB(_) => 'B',
                other => panic!("read {i}: expected a paired strand, got {other:?}"),
            };
            assert_eq!(
                actual_strand, strand,
                "strand mismatch at read {i}; actual={actual:?} expected={expected:?}"
            );
        }
    }

    /// GRP3-01 fgbio oracle: pin the fgbio-derived base-molecule partition AND absolute
    /// `/A` `/B` strand for the reverse-aware paired cases, and assert BOTH the sequential
    /// `PairedUmiAssigner` and the parallel `ParallelPairedAssigner` reproduce it.
    ///
    /// The other reverse-orientation suites compare parallel against the sequential
    /// assigner only, so a defect shared by both would pass. This table's expectations
    /// are hand-derived from fgbio's `PairedUmiAssigner` (see
    /// [`assert_matches_fgbio_oracle`]), giving an oracle independent of either fgumi
    /// implementation. Each case carries its own `max_mismatches` threshold. Covers reverse
    /// edges, palindrome roots, palindrome children, equal-count tie-breaking, the child-count
    /// threshold boundary (inclusion and exclusion) at `max_mismatches = 1`, and a reverse-only
    /// edge at `max_mismatches = 2`, per the grouping path's fgbio-parity requirement.
    #[rstest]
    // Reverse-orientation edge that canonicalization alone misses; equal counts, so the
    // root is the lexically-smaller canonical "ATTT-CAAA" (=> /A) and "CAAA-GTTT" the
    // reversed descendant (=> /B). The tie-break choosing the root is what fixes strand
    // here: rooting on the other member would invert both labels.
    #[case::reverse_edge_equal_count_tiebreak(
        &["CAAA-GTTT", "ATTT-CAAA"],
        &[(0, 'B'), (0, 'A')],
        1,
    )]
    // Palindrome root (halves equal => umi == reverse(umi)): fgbio maps the root canonical
    // spelling then its (identical) reverse, so last-write-wins lands the palindrome root
    // on /B. Its non-palindrome descendant "ACGT-ACGA", reached in the reverse spelling,
    // takes /A.
    #[case::palindrome_root(
        &["ACGT-ACGT", "ACGT-ACGT", "ACGT-ACGA"],
        &[(0, 'B'), (0, 'B'), (0, 'A')],
        1,
    )]
    // Palindrome descendant "AAAA-AAAA" under a non-palindrome root "AAAA-AAAC": fgbio's
    // else-branch maps the palindrome child to /B then its identical reverse to /A, so
    // last-write-wins lands the palindrome child on /A. The reverse spelling of the root
    // ("AAAC-AAAA", index 2) is the strand-/B counterpart, so not every read is /A.
    #[case::palindrome_child(
        &["AAAA-AAAC", "AAAA-AAAC", "AAAC-AAAA", "AAAA-AAAA"],
        &[(0, 'A'), (0, 'A'), (0, 'B'), (0, 'A')],
        1,
    )]
    // Threshold boundary, INCLUDED: descendant count 3 == root_count(5)/2 + 1 == 3, so the
    // "AAAA-CCCG" node joins the "AAAA-CCCC" molecule (one group).
    #[case::threshold_boundary_included(
        &[
            "AAAA-CCCC", "AAAA-CCCC", "AAAA-CCCC", "AAAA-CCCC", "AAAA-CCCC",
            "AAAA-CCCG", "AAAA-CCCG", "AAAA-CCCG",
        ],
        &[(0, 'A'), (0, 'A'), (0, 'A'), (0, 'A'), (0, 'A'), (0, 'A'), (0, 'A'), (0, 'A')],
        1,
    )]
    // Threshold boundary, EXCLUDED: descendant count 4 > root_count(5)/2 + 1 == 3, so the
    // "AAAA-CCCG" node splits off into its own molecule (group 1), each a /A root.
    #[case::threshold_boundary_excluded(
        &[
            "AAAA-CCCC", "AAAA-CCCC", "AAAA-CCCC", "AAAA-CCCC", "AAAA-CCCC",
            "AAAA-CCCG", "AAAA-CCCG", "AAAA-CCCG", "AAAA-CCCG",
        ],
        &[(0, 'A'), (0, 'A'), (0, 'A'), (0, 'A'), (0, 'A'), (1, 'A'), (1, 'A'), (1, 'A'), (1, 'A')],
        1,
    )]
    // Reverse-only edge that needs BOTH orientation-awareness AND a distance-2 threshold: the
    // two canonical forms "ATTT-CAAA" and "CAAT-GTTT" are 7 mismatches apart forward, so they
    // never draw a forward edge, but reverse("ATTT-CAAA") = "CAAA-ATTT" is exactly 2 mismatches
    // from "CAAT-GTTT" (positions 4 and 6). At max_mismatches = 1 they are separate molecules;
    // at max_mismatches = 2 the reverse edge merges them. Equal counts, so the root is the
    // lexically-smaller canonical "ATTT-CAAA" (=> /A); its reverse-orientation descendant
    // "CAAT-GTTT" is closer to the root's reverse than to the root, so its canonical spelling
    // takes /B. This exercises the k=2 fgbio oracle path, which the other reverse cases (all
    // k=1) leave to sequential-vs-parallel parity only.
    #[case::reverse_only_distance_two(
        &["ATTT-CAAA", "CAAT-GTTT"],
        &[(0, 'A'), (0, 'B')],
        2,
    )]
    fn test_paired_matches_fgbio_oracle(
        #[case] umis: &[&str],
        #[case] expected: &[(u32, char)],
        #[case] max_mismatches: u32,
    ) {
        let umis: Vec<Umi> = umis.iter().map(|s| (*s).to_string()).collect();

        // Both the sequential assigner AND the parallel assigner (at multiple thread
        // counts) must reproduce the independent fgbio oracle at the case's threshold.
        let sequential = crate::umi::PairedUmiAssigner::new(max_mismatches).assign(&umis);
        assert_matches_fgbio_oracle(&sequential, expected);
        for threads in [1usize, 4, 16] {
            let parallel = ParallelPairedAssigner::new(max_mismatches, threads).assign(&umis);
            assert_matches_fgbio_oracle(&parallel, expected);
        }
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

    /// GRP3-T2: parallel-vs-sequential PAIRED parity property test.
    ///
    /// Builds a pool of paired UMIs from a handful of random molecules, emitting each
    /// molecule in both orientations (`L-R` and `R-L`) and occasionally mutating a
    /// single base — the exact shape that produces reverse-orientation adjacency
    /// edges (GRP3-01). The parallel `ParallelPairedAssigner` must produce a grouping
    /// equivalent to the sequential `PairedUmiAssigner` (identical base-molecule AND
    /// strand partition) at every thread count, which also pins thread-determinism
    /// (threads 1 ≡ 4 ≡ 16).
    mod paired_parity_proptest {
        use super::*;
        use proptest::prelude::*;

        const BASES: [char; 4] = ['A', 'C', 'G', 'T'];

        /// Build the UMI pool from random molecules and per-read (molecule, orientation,
        /// optional single-base mutation) descriptors.
        fn build_pool(
            molecules: &[(Vec<u8>, Vec<u8>)],
            reads: &[(usize, bool, Option<usize>)],
        ) -> Vec<String> {
            let half = |v: &[u8]| -> String { v.iter().map(|&b| BASES[b as usize]).collect() };
            let mols: Vec<(String, String)> =
                molecules.iter().map(|(l, r)| (half(l), half(r))).collect();
            reads
                .iter()
                .map(|&(mi, reversed, mutate)| {
                    let (l, r) = &mols[mi % mols.len()];
                    let mut chars: Vec<char> = if reversed {
                        format!("{r}-{l}").chars().collect()
                    } else {
                        format!("{l}-{r}").chars().collect()
                    };
                    if let Some(pos) = mutate {
                        // Map a base index 0..8 to a string index, skipping the dash at 4.
                        let idx = if pos < 4 { pos } else { pos + 1 };
                        if let Some(cur) = chars.get(idx).copied() {
                            if cur != '-' {
                                let ci = BASES.iter().position(|&c| c == cur).unwrap_or(0);
                                chars[idx] = BASES[(ci + 1) % 4];
                            }
                        }
                    }
                    chars.into_iter().collect()
                })
                .collect()
        }

        proptest! {
            #![proptest_config(ProptestConfig::with_cases(300))]
            #[test]
            fn prop_parallel_paired_matches_sequential(
                molecules in prop::collection::vec(
                    (prop::collection::vec(0u8..4, 4), prop::collection::vec(0u8..4, 4)),
                    1..=5,
                ),
                reads in prop::collection::vec(
                    (0usize..5, any::<bool>(), prop::option::of(0usize..8)),
                    1..=24,
                ),
                // Exercise the k=2 reverse-edge path (`generate_neighbors_k`, k>1) too,
                // not just the k=1 fast path. These pools stay well under the sequential
                // assigner's index threshold, so the sequential also uses its matcher
                // (reverse-aware) path and remains the correct oracle at k=2.
                max_mismatches in 1u32..=2,
            ) {
                let umis = build_pool(&molecules, &reads);
                let sequential = crate::umi::PairedUmiAssigner::new(max_mismatches).assign(&umis);
                for threads in [1usize, 4, 16] {
                    let parallel = ParallelPairedAssigner::new(max_mismatches, threads).assign(&umis);
                    prop_assert!(
                        assignments_equivalent(&sequential, &parallel),
                        "max_mismatches={max_mismatches} threads={threads}\n  umis={umis:?}\n  \
                         seq={sequential:?}\n  par={parallel:?}"
                    );
                    // Same thread count, run twice: the two runs build independent
                    // edge-set `AHashMap`s/`AHashSet`s (seeded per process at runtime), so
                    // this is the direct probe for the adjacency-order nondeterminism the
                    // adjacency-bucket sort in `ParallelPairedAssigner::assign` guards
                    // against: without that sort, a threshold-boundary assignment could
                    // differ between these two runs.
                    let parallel_again =
                        ParallelPairedAssigner::new(max_mismatches, threads).assign(&umis);
                    prop_assert!(
                        assignments_equivalent(&parallel, &parallel_again),
                        "same-thread nondeterminism: max_mismatches={max_mismatches} \
                         threads={threads}\n  umis={umis:?}\n  run1={parallel:?}\n  \
                         run2={parallel_again:?}"
                    );
                }
            }
        }
    }
}
