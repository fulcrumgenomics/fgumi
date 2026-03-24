//! Loser tree (tournament tree) for efficient k-way merge.
//!
//! A loser tree supports k-way merge with `log2(k)` comparisons per element,
//! compared to `2·log2(k)` for a binary heap. Each internal node stores the
//! *loser* of the comparison at that level, and updating after replacing the
//! winner only requires walking leaf-to-root with one comparison per level.
//!
//! # Algorithm
//!
//! Uses the incremental insertion approach (no power-of-2 padding required):
//! - Array of `k` entries: `losers[0]` = overall winner, `losers[1..k]` = internal losers
//! - Leaf mapping: source `i` maps to internal node `(k + i) / 2`
//! - Parent navigation: `parent(node) = node / 2`
//! - Build: insert sources `0..k` one at a time, each sifting up from leaf to root
//!
//! Based on the approach used in Apache `DataFusion` (`arrow-rs`).
//!
//! # References
//!
//! - Knuth, *The Art of Computer Programming*, Vol. 3, Section 5.4.1
//! - Apache `DataFusion` `SortPreservingMergeStream` (PR #4301)
//! - Grafana blog: "The Loser Tree Data Structure"

use std::cmp::Ordering;

/// Sentinel value indicating an unoccupied internal node during build.
const EMPTY: usize = usize::MAX;

/// A loser tree for k-way merge of sorted sequences.
///
/// The tree has `k` leaves (one per source). Internal nodes store the source
/// index of the loser at each tournament level. `losers[0]` holds the overall
/// winner (the source with the current minimum key).
pub struct LoserTree<K> {
    /// `losers[0]` = winner, `losers[1..k]` = internal loser nodes.
    losers: Vec<usize>,
    /// Current key for each source.
    keys: Vec<K>,
    /// Number of sources.
    k: usize,
    /// Per-source active flag. Inactive sources always lose.
    active: Vec<bool>,
    /// Number of sources still active.
    num_active: usize,
}

impl<K: Ord> LoserTree<K> {
    /// Create a new loser tree from initial keys (one per source).
    ///
    /// # Panics
    ///
    /// Panics if `initial_keys` is empty.
    #[must_use]
    pub fn new(initial_keys: Vec<K>) -> Self {
        let k = initial_keys.len();
        assert!(k > 0, "loser tree requires at least one source");

        let mut tree = Self {
            losers: vec![EMPTY; k],
            keys: initial_keys,
            k,
            active: vec![true; k],
            num_active: k,
        };
        tree.build();
        tree
    }

    /// Build the tournament by inserting each source from `0` to `k-1`.
    ///
    /// Each source sifts up from its leaf node toward the root. At occupied
    /// internal nodes, the loser stays and the winner continues up. At empty
    /// nodes (value = `EMPTY`), the winner is written and sifting stops.
    fn build(&mut self) {
        for i in 0..self.k {
            let mut winner = i;
            let mut node = self.leaf_to_node(i);

            while node > 0 && self.losers[node] != EMPTY {
                let stored = self.losers[node];
                if self.is_greater(winner, stored) {
                    self.losers[node] = winner;
                    winner = stored;
                }
                node >>= 1;
            }

            self.losers[node] = winner;
        }
    }

    /// Sift source `source` up from its leaf, replaying the tournament.
    /// Used after replacing or removing the winner. Cost: `log2(k)` comparisons.
    fn replay(&mut self, source: usize) {
        let mut winner = source;
        let mut node = self.leaf_to_node(source);

        while node > 0 {
            let challenger = self.losers[node];
            if self.is_greater(winner, challenger) {
                self.losers[node] = winner;
                winner = challenger;
            }
            node >>= 1;
        }

        self.losers[0] = winner;
    }

    /// Map source index to its internal node index.
    #[inline]
    fn leaf_to_node(&self, source: usize) -> usize {
        (self.k + source) >> 1
    }

    /// Returns true if source `a` is greater than source `b` (a loses to b).
    /// Inactive sources are always greater. Ties broken by source index.
    #[inline]
    fn is_greater(&self, a: usize, b: usize) -> bool {
        match (self.active[a], self.active[b]) {
            (true, true) => {
                matches!(self.keys[a].cmp(&self.keys[b]).then_with(|| a.cmp(&b)), Ordering::Greater)
            }
            (true, false) => false,
            (false, true) => true,
            (false, false) => a > b,
        }
    }

    /// Get the index of the current winner (source with minimum key).
    #[inline]
    #[must_use]
    pub fn winner(&self) -> usize {
        self.losers[0]
    }

    /// Check if the winner is still an active source.
    #[inline]
    #[must_use]
    pub fn winner_is_active(&self) -> bool {
        self.num_active > 0 && self.active[self.losers[0]]
    }

    /// Replace the winner's key with a new key and replay the tournament.
    /// Costs `log2(k)` comparisons.
    #[inline]
    pub fn replace_winner(&mut self, new_key: K) {
        let w = self.losers[0];
        self.keys[w] = new_key;
        self.replay(w);
    }

    /// Mark the winner's source as exhausted and replay the tournament.
    /// Costs `log2(k)` comparisons.
    pub fn remove_winner(&mut self) {
        let w = self.losers[0];
        self.active[w] = false;
        self.num_active -= 1;
        self.replay(w);
    }

    /// Get a reference to the winner's current key.
    #[inline]
    #[must_use]
    pub fn winner_key(&self) -> &K {
        &self.keys[self.losers[0]]
    }

    /// Number of sources in the tree (including exhausted ones).
    #[inline]
    #[must_use]
    pub fn len(&self) -> usize {
        self.k
    }

    /// Whether the tree has no sources.
    ///
    /// Always returns `false` since the constructor requires at least one source.
    /// Present for clippy's `len_without_is_empty` lint.
    #[inline]
    #[must_use]
    pub fn is_empty(&self) -> bool {
        self.k == 0
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_single_source() {
        let mut tree = LoserTree::new(vec![5]);
        assert_eq!(tree.winner(), 0);
        assert_eq!(*tree.winner_key(), 5);
        assert!(tree.winner_is_active());

        tree.replace_winner(10);
        assert_eq!(tree.winner(), 0);
        assert_eq!(*tree.winner_key(), 10);

        tree.remove_winner();
        assert!(!tree.winner_is_active());
    }

    #[test]
    fn test_two_sources() {
        let mut tree = LoserTree::new(vec![3, 1]);
        assert_eq!(tree.winner(), 1);
        assert_eq!(*tree.winner_key(), 1);

        tree.replace_winner(5);
        assert_eq!(tree.winner(), 0);
        assert_eq!(*tree.winner_key(), 3);
    }

    #[test]
    fn test_three_sources() {
        let tree = LoserTree::new(vec![30, 10, 20]);
        assert_eq!(tree.winner(), 1);
    }

    #[test]
    fn test_merge_three_sorted_sequences() {
        let sources = [vec![1, 4, 7], vec![2, 5, 8], vec![3, 6, 9]];
        let mut indices = [0usize; 3];
        let initial_keys: Vec<i32> = sources.iter().map(|s| s[0]).collect();
        let mut tree = LoserTree::new(initial_keys);

        let mut result = Vec::new();
        while tree.winner_is_active() {
            let w = tree.winner();
            result.push(*tree.winner_key());
            indices[w] += 1;
            if indices[w] < sources[w].len() {
                tree.replace_winner(sources[w][indices[w]]);
            } else {
                tree.remove_winner();
            }
        }

        assert_eq!(result, vec![1, 2, 3, 4, 5, 6, 7, 8, 9]);
    }

    #[test]
    fn test_merge_with_duplicates_stable() {
        let sources = [vec![1, 3, 3], vec![1, 2, 3]];
        let mut indices = [0usize; 2];
        let initial_keys: Vec<i32> = sources.iter().map(|s| s[0]).collect();
        let mut tree = LoserTree::new(initial_keys);

        let mut result = Vec::new();
        let mut winner_sources = Vec::new();
        while tree.winner_is_active() {
            let w = tree.winner();
            result.push(*tree.winner_key());
            winner_sources.push(w);
            indices[w] += 1;
            if indices[w] < sources[w].len() {
                tree.replace_winner(sources[w][indices[w]]);
            } else {
                tree.remove_winner();
            }
        }

        assert_eq!(result, vec![1, 1, 2, 3, 3, 3]);
        assert_eq!(winner_sources[0], 0);
        assert_eq!(winner_sources[1], 1);
    }

    #[test]
    fn test_many_sources_descending_keys() {
        let keys: Vec<i32> = (0..16).rev().collect();
        let mut tree = LoserTree::new(keys);

        let mut result = Vec::new();
        while tree.winner_is_active() {
            result.push(*tree.winner_key());
            tree.remove_winner();
        }

        assert_eq!(result, (0..16).collect::<Vec<i32>>());
    }

    #[test]
    fn test_non_power_of_two() {
        let keys = vec![10, 30, 20, 50, 40];
        let mut tree = LoserTree::new(keys);

        let mut result = Vec::new();
        while tree.winner_is_active() {
            result.push(*tree.winner_key());
            tree.remove_winner();
        }

        assert_eq!(result, vec![10, 20, 30, 40, 50]);
    }

    #[test]
    fn test_merge_longer_sequences() {
        let sources = [vec![1, 5, 9, 13], vec![2, 6, 10], vec![3, 7, 11, 14, 15], vec![4, 8, 12]];
        let mut indices = [0usize; 4];
        let initial_keys: Vec<i32> = sources.iter().map(|s| s[0]).collect();
        let mut tree = LoserTree::new(initial_keys);

        let mut result = Vec::new();
        while tree.winner_is_active() {
            let w = tree.winner();
            result.push(*tree.winner_key());
            indices[w] += 1;
            if indices[w] < sources[w].len() {
                tree.replace_winner(sources[w][indices[w]]);
            } else {
                tree.remove_winner();
            }
        }

        assert_eq!(result, (1..=15).collect::<Vec<i32>>());
    }

    #[test]
    fn test_large_fan_in() {
        let n: i32 = 64;
        let sources: Vec<Vec<i32>> = (0..n).map(|i| vec![i]).collect();
        let initial_keys: Vec<i32> = sources.iter().map(|s| s[0]).collect();
        let mut tree = LoserTree::new(initial_keys);

        let mut result = Vec::new();
        while tree.winner_is_active() {
            result.push(*tree.winner_key());
            tree.remove_winner();
        }

        assert_eq!(result, (0..n).collect::<Vec<i32>>());
    }

    #[test]
    fn test_all_same_keys_stable() {
        let keys = vec![42, 42, 42, 42];
        let mut tree = LoserTree::new(keys);

        let mut winners = Vec::new();
        while tree.winner_is_active() {
            winners.push(tree.winner());
            tree.remove_winner();
        }

        assert_eq!(winners, vec![0, 1, 2, 3]);
    }

    #[test]
    fn test_k_equals_7() {
        let keys = vec![70, 10, 50, 30, 60, 20, 40];
        let mut tree = LoserTree::new(keys);

        let mut result = Vec::new();
        while tree.winner_is_active() {
            result.push(*tree.winner_key());
            tree.remove_winner();
        }

        assert_eq!(result, vec![10, 20, 30, 40, 50, 60, 70]);
    }

    #[test]
    fn test_interleaved_merge() {
        let sources = [vec![1, 3, 5, 7, 9], vec![2, 4, 6, 8, 10]];
        let mut indices = [0usize; 2];
        let initial_keys: Vec<i32> = sources.iter().map(|s| s[0]).collect();
        let mut tree = LoserTree::new(initial_keys);

        let mut result = Vec::new();
        while tree.winner_is_active() {
            let w = tree.winner();
            result.push(*tree.winner_key());
            indices[w] += 1;
            if indices[w] < sources[w].len() {
                tree.replace_winner(sources[w][indices[w]]);
            } else {
                tree.remove_winner();
            }
        }

        assert_eq!(result, (1..=10).collect::<Vec<i32>>());
    }
}
