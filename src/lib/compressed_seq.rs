//! Compressed DNA sequence storage for reference genomes.
//!
//! This module provides memory-efficient storage of DNA reference sequences
//! using 2-bit encoding with support for N bases and case preservation.
//!
//! Memory usage comparison for a 3GB human genome:
//! - Uncompressed (1 byte/base): ~3 GB
//! - Compressed: ~0.5 bytes/base = ~1.5 GB
//!   - 2-bit sequence: 0.25 bytes/base
//!   - Lowercase bitmap: 0.125 bytes/base
//!   - N positions bitmap: 0.125 bytes/base

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
    /// Bitmap for N positions (both 'N' and 'n'), 64 positions per u64.
    n_mask: Vec<u64>,
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
        let bitmap_len = len.div_ceil(64); // 64 positions per u64

        let mut bits = vec![0u64; bits_len];
        let mut lowercase = vec![0u64; bitmap_len];
        let mut n_mask = vec![0u64; bitmap_len];

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
                    // Track N in bitmap (same pattern as lowercase)
                    let word_idx = i / 64;
                    let bit_idx = i % 64;
                    n_mask[word_idx] |= 1u64 << bit_idx;
                    0u64 // Store as A, track in n_mask
                }
                _ => 0u64, // Unknown bases treated as A
            };

            let word_idx = i / 32;
            let bit_idx = (i % 32) * 2;
            bits[word_idx] |= encoded << bit_idx;
        }

        Self { bits, lowercase, n_mask, len }
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

        // Check lowercase and N using bitmaps (same pattern, same word index)
        let bitmap_word_idx = pos / 64;
        let bitmap_bit_idx = pos % 64;
        let is_lowercase = (self.lowercase[bitmap_word_idx] >> bitmap_bit_idx) & 1 == 1;
        let is_n = (self.n_mask[bitmap_word_idx] >> bitmap_bit_idx) & 1 == 1;

        // Check if it's an N
        if is_n {
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
    ///
    /// This method is optimized for batch decoding - it processes multiple bases
    /// per iteration by working word-by-word through the underlying bit vectors,
    /// rather than calling `base_at()` for each position individually.
    #[must_use]
    pub fn fetch(&self, start: usize, end: usize) -> Option<Vec<u8>> {
        if start > end || end > self.len {
            return None;
        }

        let len = end - start;
        if len == 0 {
            return Some(Vec::new());
        }

        let mut result = Vec::with_capacity(len);

        // Lookup table for 2-bit to base conversion
        const DECODE_TABLE: [u8; 4] = [b'A', b'C', b'G', b'T'];

        // Process each position, but batch the word lookups
        let mut pos = start;
        while pos < end {
            // Calculate word indices for current position
            let bits_word_idx = pos / 32;
            let bitmap_word_idx = pos / 64;

            // Load words once for this chunk
            let bits_word = self.bits[bits_word_idx];
            let n_word = self.n_mask[bitmap_word_idx];
            let lc_word = self.lowercase[bitmap_word_idx];

            // Calculate how many bases we can process from current words
            // For bits: up to end of current 32-base word
            // For bitmaps: up to end of current 64-position word
            let bits_word_end = ((bits_word_idx + 1) * 32).min(end);
            let bitmap_word_end = ((bitmap_word_idx + 1) * 64).min(end);
            let chunk_end = bits_word_end.min(bitmap_word_end);

            // Process all bases in this chunk using the cached words
            while pos < chunk_end {
                let bitmap_bit_idx = pos % 64;
                let is_n = (n_word >> bitmap_bit_idx) & 1 == 1;
                let is_lowercase = (lc_word >> bitmap_bit_idx) & 1 == 1;

                let base = if is_n {
                    if is_lowercase { b'n' } else { b'N' }
                } else {
                    let bits_bit_idx = (pos % 32) * 2;
                    let encoded = ((bits_word >> bits_bit_idx) & 0b11) as usize;
                    let b = DECODE_TABLE[encoded];
                    if is_lowercase { b.to_ascii_lowercase() } else { b }
                };

                result.push(base);
                pos += 1;
            }
        }

        Some(result)
    }

    /// Calculate the memory usage in bytes.
    #[must_use]
    pub fn memory_usage(&self) -> usize {
        let bits_bytes = self.bits.capacity() * std::mem::size_of::<u64>();
        let lowercase_bytes = self.lowercase.capacity() * std::mem::size_of::<u64>();
        let n_mask_bytes = self.n_mask.capacity() * std::mem::size_of::<u64>();
        let struct_overhead = std::mem::size_of::<Self>();

        bits_bytes + lowercase_bytes + n_mask_bytes + struct_overhead
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
        // With n_mask bitmap: 2-bit seq (0.25) + lowercase (0.125) + n_mask (0.125) = 0.5 bytes/base
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

        // Should use less than 60% of uncompressed size (actual is ~50%)
        assert!(
            memory < (uncompressed * 6 / 10),
            "Memory usage {memory} should be less than {}/10 * 6",
            uncompressed
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

    // === Word boundary tests (critical for batch decode optimization) ===

    #[test]
    fn test_fetch_within_single_word() {
        // 32 bases per u64 word for bits
        let seq = CompressedSequence::from_bytes(b"ACGTACGTACGTACGTACGTACGTACGTACGT"); // 32 bases
        assert_eq!(seq.fetch(0, 32), Some(b"ACGTACGTACGTACGTACGTACGTACGTACGT".to_vec()));
        assert_eq!(seq.fetch(4, 8), Some(b"ACGT".to_vec()));
        assert_eq!(seq.fetch(28, 32), Some(b"ACGT".to_vec()));
    }

    #[test]
    fn test_fetch_crossing_bits_word_boundary() {
        // Fetch spanning positions 30-34 crosses the 32-base bits word boundary
        let seq = CompressedSequence::from_bytes(&b"ACGT".repeat(20)); // 80 bases
        assert_eq!(seq.fetch(30, 34), Some(b"GTAC".to_vec()));
    }

    #[test]
    fn test_fetch_crossing_bitmap_word_boundary() {
        // 64 positions per u64 word for bitmaps (lowercase and n_mask)
        let seq = CompressedSequence::from_bytes(&b"ACGT".repeat(20)); // 80 bases
        // Fetch spanning positions 62-66 crosses the 64-position bitmap boundary
        assert_eq!(seq.fetch(62, 66), Some(b"GTAC".to_vec()));
    }

    #[test]
    fn test_fetch_spanning_multiple_words() {
        let seq = CompressedSequence::from_bytes(&b"ACGT".repeat(30)); // 120 bases = ~4 words
        assert_eq!(seq.fetch(0, 120), Some(b"ACGT".repeat(30)));
        assert_eq!(seq.fetch(10, 100), Some(b"ACGT".repeat(30)[10..100].to_vec()));
    }

    // === N base tests (bitmap) ===

    #[test]
    fn test_fetch_with_n_bases() {
        let seq = CompressedSequence::from_bytes(b"ACNGT");
        assert_eq!(seq.fetch(0, 5), Some(b"ACNGT".to_vec()));
        assert_eq!(seq.fetch(2, 3), Some(b"N".to_vec()));
    }

    #[test]
    fn test_fetch_all_n_sequence() {
        let seq = CompressedSequence::from_bytes(b"NNNNNNNN");
        assert_eq!(seq.fetch(0, 8), Some(b"NNNNNNNN".to_vec()));
    }

    #[test]
    fn test_fetch_n_at_word_boundary() {
        // N at position 31 and 32 (crossing bits word boundary)
        let mut bytes = b"ACGT".repeat(8); // 32 bases: ends with ...G(30)T(31)
        bytes[31] = b'N'; // position 31 becomes N
        bytes.extend_from_slice(b"NACGT"); // positions 32-36: N(32)A(33)C(34)G(35)T(36)
        let seq = CompressedSequence::from_bytes(&bytes);
        // Positions 30-34: G(30), N(31), N(32), A(33), C(34)
        assert_eq!(seq.fetch(30, 35), Some(b"GNNAC".to_vec()));
    }

    #[test]
    fn test_fetch_n_at_bitmap_boundary() {
        // N at positions 63 and 64 (crossing 64-bit bitmap boundary)
        let mut bytes = vec![b'A'; 70];
        bytes[63] = b'N';
        bytes[64] = b'N';
        let seq = CompressedSequence::from_bytes(&bytes);
        assert_eq!(seq.fetch(62, 66), Some(b"ANNA".to_vec()));
    }

    // === Lowercase tests (bitmap) ===

    #[test]
    fn test_fetch_with_lowercase() {
        let seq = CompressedSequence::from_bytes(b"AcGt");
        assert_eq!(seq.fetch(0, 4), Some(b"AcGt".to_vec()));
    }

    #[test]
    fn test_fetch_lowercase_n() {
        let seq = CompressedSequence::from_bytes(b"ACnGT");
        assert_eq!(seq.fetch(0, 5), Some(b"ACnGT".to_vec()));
        assert_eq!(seq.fetch(2, 3), Some(b"n".to_vec()));
    }

    #[test]
    fn test_fetch_lowercase_at_bitmap_boundary() {
        // lowercase at position 63 and 64 (crossing 64-bit lowercase bitmap boundary)
        let mut bytes = vec![b'A'; 70];
        bytes[63] = b'c';
        bytes[64] = b'g';
        let seq = CompressedSequence::from_bytes(&bytes);
        assert_eq!(seq.fetch(62, 66), Some(b"AcgA".to_vec()));
    }

    // === Large sequence tests ===

    #[test]
    fn test_fetch_large_sequence() {
        // 10,000 bases - exercises multiple words
        let bases: Vec<u8> = (0..10_000)
            .map(|i| match i % 4 {
                0 => b'A',
                1 => b'C',
                2 => b'G',
                _ => b'T',
            })
            .collect();
        let seq = CompressedSequence::from_bytes(&bases);

        // Fetch various ranges
        assert_eq!(seq.fetch(0, 100).unwrap().len(), 100);
        assert_eq!(seq.fetch(5000, 5100).unwrap().len(), 100);
        assert_eq!(seq.fetch(9900, 10000).unwrap().len(), 100);

        // Verify content
        assert_eq!(seq.fetch(0, 4), Some(b"ACGT".to_vec()));
        assert_eq!(seq.fetch(9996, 10000), Some(b"ACGT".to_vec()));
    }

    #[test]
    fn test_fetch_with_scattered_n_bases() {
        // N bases at various positions including word boundaries
        let mut bytes: Vec<u8> = (0..1000)
            .map(|i| match i % 4 {
                0 => b'A',
                1 => b'C',
                2 => b'G',
                _ => b'T',
            })
            .collect();
        // Add N's at word boundaries and random positions
        for pos in [0, 31, 32, 63, 64, 100, 500, 999] {
            bytes[pos] = b'N';
        }
        let seq = CompressedSequence::from_bytes(&bytes);

        // Verify N positions are correct
        let fetched = seq.fetch(0, 1000).unwrap();
        for pos in [0, 31, 32, 63, 64, 100, 500, 999] {
            assert_eq!(fetched[pos], b'N', "N not found at position {pos}");
        }
    }

    // === Regression tests (ensure fetch matches base_at) ===

    #[test]
    fn test_fetch_matches_base_at() {
        // Ensure batch fetch produces same results as individual base_at calls
        let seq = CompressedSequence::from_bytes(b"AcNgTACGTacgtNNNN");
        let fetched = seq.fetch(0, 17).unwrap();
        for (i, &expected) in fetched.iter().enumerate() {
            assert_eq!(seq.base_at(i), Some(expected), "Mismatch at position {i}");
        }
    }

    #[test]
    fn test_fetch_consistency_random_ranges() {
        let seq = CompressedSequence::from_bytes(&b"ACGTNacgtn".repeat(100));
        // Test various range combinations
        for start in (0..900).step_by(37) {
            for len in [1, 10, 33, 64, 100] {
                let end = (start + len).min(1000);
                let fetched = seq.fetch(start, end).unwrap();
                for (i, &expected) in fetched.iter().enumerate() {
                    assert_eq!(
                        seq.base_at(start + i),
                        Some(expected),
                        "Mismatch at position {} (range {start}..{end})",
                        start + i
                    );
                }
            }
        }
    }

    #[test]
    fn test_fetch_single_base_all_positions() {
        // Test fetching single bases at each position matches base_at
        let seq = CompressedSequence::from_bytes(&b"ACGTNacgtn".repeat(10));
        for i in 0..100 {
            let fetched = seq.fetch(i, i + 1).unwrap();
            assert_eq!(fetched.len(), 1, "Single base fetch should return 1 element");
            assert_eq!(
                Some(fetched[0]),
                seq.base_at(i),
                "Single base fetch mismatch at position {i}"
            );
        }
    }

    #[test]
    fn test_fetch_invalid_range() {
        let seq = CompressedSequence::from_bytes(b"ACGT");
        assert_eq!(seq.fetch(5, 6), None); // start beyond end
        assert_eq!(seq.fetch(0, 5), None); // end beyond length
        assert_eq!(seq.fetch(3, 2), None); // start > end
    }
}
