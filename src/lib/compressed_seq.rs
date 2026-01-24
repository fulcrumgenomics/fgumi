//! Compressed DNA sequence storage for reference genomes.
//!
//! This module provides memory-efficient storage of DNA reference sequences
//! using 2-bit encoding with support for N bases and case preservation.
//!
//! Memory usage comparison for a 3GB human genome:
//! - Uncompressed (1 byte/base): ~3 GB
//! - Compressed: ~0.375 bytes/base = ~1.1 GB
//!   - 2-bit sequence: 0.25 bytes/base
//!   - Lowercase bitmap: 0.125 bytes/base
//!   - N positions: negligible (sparse `HashSet`)

use std::collections::HashSet;

/// A memory-efficient representation of a DNA sequence.
///
/// Uses 2-bit encoding (A=0, C=1, G=2, T=3) stored in Vec<u64>,
/// with separate tracking for N bases and lowercase positions.
#[derive(Clone, Debug)]
pub struct CompressedSequence {
    /// 2-bit encoded bases, 32 bases per u64.
    /// N bases are encoded as A (0) but tracked separately.
    bits: Vec<u64>,
    /// Bitmap for lowercase positions (including lowercase 'n'), 64 positions per u64.
    lowercase: Vec<u64>,
    /// Set of positions (0-indexed) containing N bases (both 'N' and 'n').
    n_positions: HashSet<usize>,
    /// Total number of bases in the sequence.
    len: usize,
}

impl CompressedSequence {
    /// Create a compressed sequence from a byte slice.
    ///
    /// Supports A, C, G, T, N (and lowercase variants).
    #[must_use]
    pub fn from_bytes(seq: &[u8]) -> Self {
        let len = seq.len();

        // Calculate required capacity
        let bits_len = len.div_ceil(32); // 32 bases per u64
        let lowercase_len = len.div_ceil(64); // 64 positions per u64

        let mut bits = vec![0u64; bits_len];
        let mut lowercase = vec![0u64; lowercase_len];
        let mut n_positions = HashSet::new();

        for (i, &base) in seq.iter().enumerate() {
            // Track lowercase
            if base.is_ascii_lowercase() {
                let word_idx = i / 64;
                let bit_idx = i % 64;
                lowercase[word_idx] |= 1u64 << bit_idx;
            }

            // Encode base (uppercase for encoding)
            let encoded = match base.to_ascii_uppercase() {
                b'A' => 0u64,
                b'C' => 1u64,
                b'G' => 2u64,
                b'T' => 3u64,
                b'N' => {
                    n_positions.insert(i);
                    0u64 // Store as A, track in n_positions
                }
                _ => 0u64, // Unknown bases treated as A
            };

            let word_idx = i / 32;
            let bit_idx = (i % 32) * 2;
            bits[word_idx] |= encoded << bit_idx;
        }

        Self { bits, lowercase, n_positions, len }
    }

    /// Get the number of bases in the sequence.
    #[inline]
    #[must_use]
    pub const fn len(&self) -> usize {
        self.len
    }

    /// Check if the sequence is empty.
    #[inline]
    #[must_use]
    pub const fn is_empty(&self) -> bool {
        self.len == 0
    }

    /// Get the base at a 0-indexed position.
    ///
    /// Returns the base as a byte (ASCII), preserving case.
    /// Returns None if the position is out of bounds.
    #[inline]
    #[must_use]
    pub fn base_at(&self, pos: usize) -> Option<u8> {
        if pos >= self.len {
            return None;
        }

        // Check lowercase
        let lc_word_idx = pos / 64;
        let lc_bit_idx = pos % 64;
        let is_lowercase = (self.lowercase[lc_word_idx] >> lc_bit_idx) & 1 == 1;

        // Check if it's an N
        if self.n_positions.contains(&pos) {
            return Some(if is_lowercase { b'n' } else { b'N' });
        }

        // Extract the 2-bit encoded value
        let word_idx = pos / 32;
        let bit_idx = (pos % 32) * 2;
        let encoded = (self.bits[word_idx] >> bit_idx) & 0b11;

        // Decode to base
        let base = match encoded {
            0 => b'A',
            1 => b'C',
            2 => b'G',
            3 => b'T',
            _ => unreachable!(),
        };

        if is_lowercase { Some(base.to_ascii_lowercase()) } else { Some(base) }
    }

    /// Fetch a subsequence as a byte vector.
    ///
    /// Parameters are 0-indexed, inclusive start, exclusive end.
    /// Returns None if the range is invalid.
    #[must_use]
    pub fn fetch(&self, start: usize, end: usize) -> Option<Vec<u8>> {
        if start > end || end > self.len {
            return None;
        }

        let result: Vec<u8> = (start..end).filter_map(|i| self.base_at(i)).collect();

        if result.len() == end - start { Some(result) } else { None }
    }

    /// Calculate the memory usage in bytes.
    #[must_use]
    pub fn memory_usage(&self) -> usize {
        let bits_bytes = self.bits.capacity() * std::mem::size_of::<u64>();
        let lowercase_bytes = self.lowercase.capacity() * std::mem::size_of::<u64>();
        // HashSet overhead is approximately 24 bytes per entry plus base size
        let n_positions_bytes = self.n_positions.len() * 32;
        let struct_overhead = std::mem::size_of::<Self>();

        bits_bytes + lowercase_bytes + n_positions_bytes + struct_overhead
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_basic_encoding() {
        let seq = CompressedSequence::from_bytes(b"ACGT");
        assert_eq!(seq.len(), 4);
        assert_eq!(seq.base_at(0), Some(b'A'));
        assert_eq!(seq.base_at(1), Some(b'C'));
        assert_eq!(seq.base_at(2), Some(b'G'));
        assert_eq!(seq.base_at(3), Some(b'T'));
        assert_eq!(seq.base_at(4), None);
    }

    #[test]
    fn test_lowercase_preserved() {
        let seq = CompressedSequence::from_bytes(b"AcGt");
        assert_eq!(seq.base_at(0), Some(b'A'));
        assert_eq!(seq.base_at(1), Some(b'c'));
        assert_eq!(seq.base_at(2), Some(b'G'));
        assert_eq!(seq.base_at(3), Some(b't'));
    }

    #[test]
    fn test_n_bases() {
        let seq = CompressedSequence::from_bytes(b"ACNGT");
        assert_eq!(seq.base_at(0), Some(b'A'));
        assert_eq!(seq.base_at(1), Some(b'C'));
        assert_eq!(seq.base_at(2), Some(b'N'));
        assert_eq!(seq.base_at(3), Some(b'G'));
        assert_eq!(seq.base_at(4), Some(b'T'));
    }

    #[test]
    fn test_fetch_subsequence() {
        let seq = CompressedSequence::from_bytes(b"ACGTNNACGT");
        assert_eq!(seq.fetch(0, 4), Some(b"ACGT".to_vec()));
        assert_eq!(seq.fetch(4, 6), Some(b"NN".to_vec()));
        assert_eq!(seq.fetch(6, 10), Some(b"ACGT".to_vec()));
        assert_eq!(seq.fetch(2, 8), Some(b"GTNNAC".to_vec()));
    }

    #[test]
    fn test_fetch_boundary_cases() {
        let seq = CompressedSequence::from_bytes(b"ACGT");
        assert_eq!(seq.fetch(0, 0), Some(vec![])); // Empty range
        assert_eq!(seq.fetch(0, 4), Some(b"ACGT".to_vec())); // Full sequence
        assert_eq!(seq.fetch(4, 4), Some(vec![])); // Empty at end
        assert_eq!(seq.fetch(0, 5), None); // Beyond end
        assert_eq!(seq.fetch(5, 5), None); // Start beyond end
    }

    #[test]
    fn test_long_sequence() {
        // Test crossing u64 boundaries (32 bases per word for bits, 64 for lowercase)
        let seq_str = "ACGT".repeat(25); // 100 bases
        let seq = CompressedSequence::from_bytes(seq_str.as_bytes());

        assert_eq!(seq.len(), 100);

        // Check bases at u64 boundaries
        for i in 0..100 {
            let expected = match i % 4 {
                0 => b'A',
                1 => b'C',
                2 => b'G',
                3 => b'T',
                _ => unreachable!(),
            };
            assert_eq!(seq.base_at(i), Some(expected), "Mismatch at position {i}");
        }
    }

    #[test]
    fn test_long_sequence_with_lowercase() {
        // Test lowercase tracking across u64 boundaries
        let mut seq_bytes = vec![b'A'; 100];
        // Set positions 63 and 64 to lowercase (crosses u64 boundary for lowercase bitmap)
        seq_bytes[63] = b'c';
        seq_bytes[64] = b'g';

        let seq = CompressedSequence::from_bytes(&seq_bytes);

        assert_eq!(seq.base_at(62), Some(b'A'));
        assert_eq!(seq.base_at(63), Some(b'c'));
        assert_eq!(seq.base_at(64), Some(b'g'));
        assert_eq!(seq.base_at(65), Some(b'A'));
    }

    #[test]
    fn test_all_n_sequence() {
        let seq = CompressedSequence::from_bytes(b"NNNNNN");
        assert_eq!(seq.len(), 6);
        for i in 0..6 {
            assert_eq!(seq.base_at(i), Some(b'N'));
        }
    }

    #[test]
    fn test_empty_sequence() {
        let seq = CompressedSequence::from_bytes(b"");
        assert_eq!(seq.len(), 0);
        assert!(seq.is_empty());
        assert_eq!(seq.base_at(0), None);
        assert_eq!(seq.fetch(0, 0), Some(vec![]));
    }

    #[test]
    fn test_memory_savings() {
        // For a 1MB sequence, verify memory is less than uncompressed
        let size = 1_000_000;
        let seq_bytes: Vec<u8> = (0..size)
            .map(|i| match i % 4 {
                0 => b'A',
                1 => b'C',
                2 => b'G',
                _ => b'T',
            })
            .collect();

        let seq = CompressedSequence::from_bytes(&seq_bytes);

        let memory = seq.memory_usage();
        let uncompressed = size;

        // Should use less than 50% of uncompressed size
        assert!(
            memory < uncompressed / 2,
            "Memory usage {memory} should be less than {uncompressed}/2"
        );
    }

    #[test]
    fn test_mixed_case_all_bases() {
        let seq = CompressedSequence::from_bytes(b"AaCcGgTtNn");
        assert_eq!(seq.base_at(0), Some(b'A'));
        assert_eq!(seq.base_at(1), Some(b'a'));
        assert_eq!(seq.base_at(2), Some(b'C'));
        assert_eq!(seq.base_at(3), Some(b'c'));
        assert_eq!(seq.base_at(4), Some(b'G'));
        assert_eq!(seq.base_at(5), Some(b'g'));
        assert_eq!(seq.base_at(6), Some(b'T'));
        assert_eq!(seq.base_at(7), Some(b't'));
        assert_eq!(seq.base_at(8), Some(b'N'));
        // Lowercase 'n' is preserved
        assert_eq!(seq.base_at(9), Some(b'n'));
    }
}
