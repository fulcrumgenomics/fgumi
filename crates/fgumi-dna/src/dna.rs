//! DNA sequence utilities.
//!
//! This module provides common DNA sequence operations like reverse complement.

use crate::{NO_CALL_BASE, NO_CALL_BASE_LOWER};

/// Complements a single DNA base, normalizing to uppercase except for N/n.
///
/// Returns the Watson-Crick complement: A<->T, C<->G.
/// Both uppercase and lowercase A/T/C/G input are normalized to uppercase output.
///
/// # Case preservation for N/n
///
/// Uppercase 'N' and lowercase 'n' are returned unchanged (case preserved).
/// This is required for fgbio compatibility in CODEC consensus calling:
///
/// - **CODEC** uses lowercase 'n' for padding when R1/R2 have different lengths
///   (see `codec_caller::pad_consensus`). These padded bases are then reverse
///   complemented for ac/bc tags (see `codec_caller::reverse_complement_ss`).
///   fgbio preserves lowercase 'n' through this operation.
///
/// - **Simplex/Duplex** use only uppercase 'N' for no-call bases in production
///   (see `vanilla_caller` lines 640, 1024, 1078 and `duplex_caller` `NO_CALL` const).
///   The `padded_default` function with lowercase 'n' exists but is only used in tests.
///
/// - **`Zipper/tag_reversal`** do not create 'n' bases, only reverse complement
///   existing tag values which come from upstream processing.
#[inline]
#[must_use]
pub const fn complement_base(base: u8) -> u8 {
    match base {
        b'A' | b'a' => b'T',
        b'T' | b't' => b'A',
        b'C' | b'c' => b'G',
        b'G' | b'g' => b'C',
        NO_CALL_BASE => NO_CALL_BASE,             // 'N' stays 'N'
        NO_CALL_BASE_LOWER => NO_CALL_BASE_LOWER, // 'n' stays 'n' (fgbio compat)
        _ => base,
    }
}

/// Reverse complements a DNA sequence.
///
/// Returns the reverse complement of the input sequence, normalizing to uppercase.
/// A<->T, C<->G, N and other bases are unchanged.
///
/// # Examples
///
/// ```
/// use fgumi_dna::dna::reverse_complement;
///
/// assert_eq!(reverse_complement(b"ACGT"), b"ACGT".to_vec());
/// assert_eq!(reverse_complement(b"AAAA"), b"TTTT".to_vec());
/// assert_eq!(reverse_complement(b"ACGTN"), b"NACGT".to_vec());
/// assert_eq!(reverse_complement(b"acgt"), b"ACGT".to_vec());  // normalized to uppercase
/// ```
#[must_use]
pub fn reverse_complement(seq: &[u8]) -> Vec<u8> {
    seq.iter().rev().map(|&base| complement_base(base)).collect()
}

/// Complements a single DNA base, preserving case.
///
/// Returns the Watson-Crick complement: A<->T, C<->G.
/// Unlike [`complement_base`], this preserves the case of the input:
/// uppercase input produces uppercase output, lowercase produces lowercase.
///
/// # Examples
///
/// ```
/// use fgumi_dna::dna::complement_base_preserve_case;
///
/// assert_eq!(complement_base_preserve_case(b'A'), b'T');
/// assert_eq!(complement_base_preserve_case(b'a'), b't');
/// ```
#[inline]
#[must_use]
pub const fn complement_base_preserve_case(base: u8) -> u8 {
    match base {
        b'A' => b'T',
        b'T' => b'A',
        b'C' => b'G',
        b'G' => b'C',
        b'a' => b't',
        b't' => b'a',
        b'c' => b'g',
        b'g' => b'c',
        other => other,
    }
}

/// Reverse complements a DNA string, preserving case.
///
/// Returns the reverse complement of the input string.
/// Unlike [`reverse_complement`], this preserves the case of each base.
/// Non-ASCII characters are passed through unchanged.
///
/// # Examples
///
/// ```
/// use fgumi_dna::dna::reverse_complement_str;
///
/// assert_eq!(reverse_complement_str("ACGT"), "ACGT");
/// assert_eq!(reverse_complement_str("AAAA"), "TTTT");
/// assert_eq!(reverse_complement_str("acgt"), "acgt");
/// ```
#[must_use]
pub fn reverse_complement_str(seq: &str) -> String {
    seq.chars()
        .rev()
        .map(|c| {
            if c.is_ascii() {
                complement_base_preserve_case(c as u8) as char
            } else {
                c
            }
        })
        .collect()
}

#[cfg(test)]
mod tests {
    use rstest::rstest;

    use super::*;

    #[test]
    fn test_complement_base() {
        // Uppercase bases
        assert_eq!(complement_base(b'A'), b'T');
        assert_eq!(complement_base(b'T'), b'A');
        assert_eq!(complement_base(b'C'), b'G');
        assert_eq!(complement_base(b'G'), b'C');

        // Lowercase normalized to uppercase
        assert_eq!(complement_base(b'a'), b'T');
        assert_eq!(complement_base(b't'), b'A');
        assert_eq!(complement_base(b'c'), b'G');
        assert_eq!(complement_base(b'g'), b'C');

        // N bases preserve case (fgbio uses lowercase 'n' for padding)
        assert_eq!(complement_base(b'N'), b'N');
        assert_eq!(complement_base(b'n'), b'n');

        // IUPAC ambiguity codes unchanged
        for code in [b'R', b'Y', b'S', b'W', b'K', b'M', b'B', b'D', b'H', b'V'] {
            assert_eq!(complement_base(code), code);
        }

        // Special characters unchanged
        for c in [b'.', b'-', b'*', b'0', b'X'] {
            assert_eq!(complement_base(c), c);
        }
    }

    #[test]
    fn test_reverse_complement() {
        // Basic cases
        assert_eq!(reverse_complement(b""), b"".to_vec());
        assert_eq!(reverse_complement(b"A"), b"T".to_vec());
        assert_eq!(reverse_complement(b"AAAA"), b"TTTT".to_vec());
        assert_eq!(reverse_complement(b"TTTT"), b"AAAA".to_vec());
        assert_eq!(reverse_complement(b"CCCC"), b"GGGG".to_vec());
        assert_eq!(reverse_complement(b"GGGG"), b"CCCC".to_vec());
        assert_eq!(reverse_complement(b"ACGTN"), b"NACGT".to_vec());

        // Lowercase normalized to uppercase
        assert_eq!(reverse_complement(b"acgt"), b"ACGT".to_vec());
        assert_eq!(reverse_complement(b"AcGt"), b"ACGT".to_vec());

        // Palindromic sequences
        assert_eq!(reverse_complement(b"ACGT"), b"ACGT".to_vec());
        assert_eq!(reverse_complement(b"GAATTC"), b"GAATTC".to_vec());

        // Double operation returns original (uppercase)
        let seq = b"ACGTACGT";
        assert_eq!(reverse_complement(&reverse_complement(seq)), seq.to_vec());
    }

    #[rstest]
    #[case(b'A', b'T')]
    #[case(b'T', b'A')]
    #[case(b'C', b'G')]
    #[case(b'G', b'C')]
    #[case(b'a', b't')]
    #[case(b't', b'a')]
    #[case(b'c', b'g')]
    #[case(b'g', b'c')]
    #[case(b'N', b'N')]
    #[case(b'n', b'n')]
    #[case(b'.', b'.')]
    #[case(b'-', b'-')]
    #[case(b'*', b'*')]
    #[case(b'0', b'0')]
    #[case(b'X', b'X')]
    fn test_complement_base_preserve_case(#[case] input: u8, #[case] expected: u8) {
        assert_eq!(complement_base_preserve_case(input), expected);
    }

    #[test]
    fn test_reverse_complement_str() {
        // Basic cases
        assert_eq!(reverse_complement_str(""), "");
        assert_eq!(reverse_complement_str("ACGT"), "ACGT");
        assert_eq!(reverse_complement_str("AAAA"), "TTTT");

        // Case preserved
        assert_eq!(reverse_complement_str("acgt"), "acgt");
        assert_eq!(reverse_complement_str("AcGt"), "aCgT");

        // Special characters unchanged
        assert_eq!(reverse_complement_str("A-C-G-T"), "A-C-G-T");

        // Double operation preserves original
        for seq in ["ACGT", "acgt", "AcGt"] {
            assert_eq!(reverse_complement_str(&reverse_complement_str(seq)), seq);
        }
    }

    #[test]
    fn test_reverse_complement_str_non_ascii() {
        // Non-ASCII characters (e.g. multi-byte UTF-8) should pass through unchanged
        let input = "A\u{00e9}T"; // 'é' is a 2-byte UTF-8 character
        let result = reverse_complement_str(input);
        assert_eq!(result, "A\u{00e9}T");
    }
}
