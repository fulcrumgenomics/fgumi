//! A 2-bit DNA encoding for fast UMI comparison.
//!
//! This module provides efficient storage and comparison of DNA sequences
//! using 2 bits per base. This is particularly useful for UMI correction
//! where each observed UMI must be compared against thousands of expected UMIs.
//!
//! # Example
//!
//! ```
//! use fgumi_dna::bitenc::BitEnc;
//!
//! let umi1 = BitEnc::from_bytes(b"ACGT").unwrap();
//! let umi2 = BitEnc::from_bytes(b"ACTT").unwrap();
//! assert_eq!(umi1.hamming_distance(&umi2), 1);
//! ```

/// A 2-bit encoded DNA sequence stored in a u64.
///
/// Supports sequences up to 32 bases (64 bits).
/// Each base is encoded as: A=0, C=1, G=2, T=3.
#[derive(Clone, Copy, Debug, PartialEq, Eq, Hash)]
pub struct BitEnc {
    /// The encoded sequence, with bases packed from LSB.
    bits: u64,
    /// Number of bases in the sequence.
    len: u8,
}

impl BitEnc {
    /// Encode a single DNA base to 2 bits.
    #[inline]
    const fn encode_base(base: u8) -> Option<u64> {
        match base {
            b'A' | b'a' => Some(0),
            b'C' | b'c' => Some(1),
            b'G' | b'g' => Some(2),
            b'T' | b't' => Some(3),
            _ => None,
        }
    }

    /// Create a `BitEnc` from a byte slice.
    ///
    /// Returns None if the sequence contains non-ACGT bases or exceeds 32 bases.
    #[inline]
    #[must_use]
    pub fn from_bytes(seq: &[u8]) -> Option<Self> {
        if seq.len() > 32 {
            return None;
        }

        let mut bits: u64 = 0;
        for (i, &base) in seq.iter().enumerate() {
            let encoded = Self::encode_base(base)?;
            bits |= encoded << (i * 2);
        }

        let len = u8::try_from(seq.len()).ok()?;
        Some(Self { bits, len })
    }

    /// Create a `BitEnc` from a UMI string, skipping non-ACGT characters (e.g., dashes in paired UMIs).
    ///
    /// This is useful for paired UMIs like "ACGT-TGCA" where the dash should be ignored.
    /// Returns None if the sequence contains invalid bases (not ACGT or dash) or exceeds 32 bases.
    #[inline]
    #[must_use]
    pub fn from_umi_str(umi: &str) -> Option<Self> {
        let mut bits: u64 = 0;
        let mut base_count: usize = 0;

        for &byte in umi.as_bytes() {
            if let Some(encoded) = Self::encode_base(byte) {
                if base_count >= 32 {
                    return None;
                }
                bits |= encoded << (base_count * 2);
                base_count += 1;
            } else if byte != b'-' {
                // Invalid character (not ACGT and not dash)
                return None;
            }
            // Dash is silently skipped
        }

        let len = u8::try_from(base_count).ok()?;
        Some(Self { bits, len })
    }

    /// Get the number of bases in this sequence.
    #[inline]
    #[must_use]
    pub const fn len(&self) -> usize {
        self.len as usize
    }

    /// Check if the sequence is empty.
    #[inline]
    #[must_use]
    pub const fn is_empty(&self) -> bool {
        self.len == 0
    }

    /// Compute the Hamming distance between two encoded sequences.
    ///
    /// Both sequences must have the same length (debug assertion).
    #[inline]
    #[must_use]
    pub fn hamming_distance(&self, other: &Self) -> u32 {
        debug_assert_eq!(self.len, other.len, "Sequences must have equal length");

        // XOR to find differing bits
        let diff = self.bits ^ other.bits;

        // For 2-bit encoding, a position differs if either of its 2 bits differ
        let odd_bits = diff & 0xAAAA_AAAA_AAAA_AAAA;
        let even_bits = diff & 0x5555_5555_5555_5555;

        // Combine: a position differs if odd OR even bit is set
        let differs = (odd_bits >> 1) | even_bits;

        differs.count_ones()
    }

    /// Extract bits for bases `[start_base, start_base + len)` as a u32.
    ///
    /// Each base is 2 bits, so this can extract up to 16 bases into a u32.
    /// Used for N-gram partitioning in similarity search.
    #[inline]
    #[must_use]
    pub fn extract_bits(&self, start_base: usize, num_bases: usize) -> u32 {
        debug_assert!(num_bases <= 16, "Can only extract up to 16 bases into u32");
        debug_assert!(
            start_base + num_bases <= self.len as usize,
            "Extraction range exceeds sequence length"
        );
        let start_bit = start_base * 2;
        let mask = (1u64 << (num_bases * 2)) - 1;
        #[expect(clippy::cast_possible_truncation, reason = "masked to at most 32 bits")]
        let result = ((self.bits >> start_bit) & mask) as u32;
        result
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_from_bytes() {
        // Valid sequences
        let enc = BitEnc::from_bytes(b"ACGT").unwrap();
        assert_eq!(enc.len(), 4);
        assert_eq!(enc.bits, 0b11_10_01_00); // T=3, G=2, C=1, A=0

        // Lowercase works
        assert!(BitEnc::from_bytes(b"acgt").is_some());

        // Invalid: non-ACGT bases
        assert!(BitEnc::from_bytes(b"ACGN").is_none());

        // Invalid: too long (>32 bases)
        assert!(BitEnc::from_bytes(&[b'A'; 33]).is_none());
        assert!(BitEnc::from_bytes(&[b'A'; 32]).is_some());
    }

    #[test]
    fn test_hamming_distance() {
        // Identical sequences
        let seq1 = BitEnc::from_bytes(b"ACGTACGT").unwrap();
        let seq2 = BitEnc::from_bytes(b"ACGTACGT").unwrap();
        assert_eq!(seq1.hamming_distance(&seq2), 0);

        // One difference
        let seq3 = BitEnc::from_bytes(b"ACGTACTT").unwrap();
        assert_eq!(seq1.hamming_distance(&seq3), 1);

        // All different
        let all_a = BitEnc::from_bytes(b"AAAA").unwrap();
        let all_t = BitEnc::from_bytes(b"TTTT").unwrap();
        assert_eq!(all_a.hamming_distance(&all_t), 4);

        // Typical 18bp UMI
        let umi1 = BitEnc::from_bytes(b"AACAACACATCTACCTTC").unwrap();
        let umi2 = BitEnc::from_bytes(b"AACAACACATCTACCTTA").unwrap();
        assert_eq!(umi1.hamming_distance(&umi2), 1);
    }

    #[test]
    fn test_extract_bits() {
        let enc = BitEnc::from_bytes(b"ACGTACGT").unwrap();

        // First and last 4 bases are identical (ACGT)
        assert_eq!(enc.extract_bits(0, 4), enc.extract_bits(4, 4));

        // Single base extraction
        assert_eq!(enc.extract_bits(0, 1), 0); // A=0
        assert_eq!(enc.extract_bits(1, 1), 1); // C=1
        assert_eq!(enc.extract_bits(2, 1), 2); // G=2
        assert_eq!(enc.extract_bits(3, 1), 3); // T=3

        // Middle extraction (GTAC)
        assert_eq!(enc.extract_bits(2, 4), 0b01_00_11_10);
    }

    #[test]
    fn test_from_umi_str() {
        // Simple UMI without dash
        let enc = BitEnc::from_umi_str("ACGT").unwrap();
        assert_eq!(enc.len(), 4);
        assert_eq!(enc, BitEnc::from_bytes(b"ACGT").unwrap());

        // Paired UMI with dash - dash should be skipped
        let paired = BitEnc::from_umi_str("ACGT-TGCA").unwrap();
        assert_eq!(paired.len(), 8);
        assert_eq!(paired, BitEnc::from_bytes(b"ACGTTGCA").unwrap());

        // Real paired UMI from test data
        let real = BitEnc::from_umi_str("GTCTGAGATC-AATCTTTAAT").unwrap();
        assert_eq!(real.len(), 20);

        // Lowercase works
        let lower = BitEnc::from_umi_str("acgt-tgca").unwrap();
        assert_eq!(lower, paired);

        // Invalid character (not ACGT, not dash)
        assert!(BitEnc::from_umi_str("ACGT-NGCA").is_none());

        // Multiple dashes work
        let multi = BitEnc::from_umi_str("AC-GT-TG").unwrap();
        assert_eq!(multi.len(), 6);
        assert_eq!(multi, BitEnc::from_bytes(b"ACGTTG").unwrap());
    }

    #[test]
    fn test_paired_umi_hamming() {
        // Two paired UMIs with 1 mismatch
        let umi1 = BitEnc::from_umi_str("GTCTGAGATC-AATCTTTAAT").unwrap();
        let umi2 = BitEnc::from_umi_str("GTCTGAGATC-AATCTTTAAC").unwrap(); // T->C at end
        assert_eq!(umi1.hamming_distance(&umi2), 1);

        // Two identical paired UMIs
        let umi3 = BitEnc::from_umi_str("AAAGCGATGC-CCAGTTAACC").unwrap();
        let umi4 = BitEnc::from_umi_str("AAAGCGATGC-CCAGTTAACC").unwrap();
        assert_eq!(umi3.hamming_distance(&umi4), 0);
    }
}
