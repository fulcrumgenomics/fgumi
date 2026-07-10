//! DNA sequence utilities.
//!
//! This module provides common DNA sequence operations like reverse complement.

/// Lookup table mapping each byte to its Watson–Crick / IUPAC complement.
///
/// Single source of truth for base complementation, shared by [`complement_base`]
/// (and thus every consensus reverse-complement path) and by `fgumi-raw-bam`'s
/// per-base tag reverse-complement. It mirrors fgbio's `Sequences.complement` /
/// htsjdk `SequenceUtil`:
///
/// - **Full IUPAC**: `A<->T`, `C<->G`, `U->A`, `R<->Y`, `K<->M`, `B<->V`,
///   `D<->H`, and the self-complements `S`, `W`, `N`.
/// - **Case-preserving**: lowercase input yields lowercase output (`a->t`),
///   matching fgbio — e.g. CODEC's lowercase `n` padding survives round trips.
/// - **Unknown bytes pass through unchanged** (identity default). fgbio *throws*
///   on invalid bytes; fgumi keeps the non-panicking passthrough because such
///   bytes are unreachable on real ACGTN(+IUPAC) sequence/tag data and the
///   `simulate` command relies on passthrough (`fastq_reads.rs`). Strict
///   validation, if ever needed, belongs at a dedicated entry point, not here.
///
/// Indexed by a `u8` into a 256-entry table, so the bounds check is elided:
/// complementation is a single branch-free load.
pub const COMPLEMENT: [u8; 256] = {
    let mut table = [0u8; 256];
    // Identity default so unknown bytes pass through unchanged. Counter is a `u8`
    // (widening `as usize` index, no truncating cast) with a break at 255 to
    // cover the full 0..=255 range without overflowing the increment.
    let mut b: u8 = 0;
    loop {
        table[b as usize] = b;
        if b == u8::MAX {
            break;
        }
        b += 1;
    }
    // Uppercase Watson–Crick + IUPAC.
    table[b'A' as usize] = b'T';
    table[b'T' as usize] = b'A';
    table[b'C' as usize] = b'G';
    table[b'G' as usize] = b'C';
    table[b'U' as usize] = b'A';
    table[b'R' as usize] = b'Y';
    table[b'Y' as usize] = b'R';
    table[b'S' as usize] = b'S';
    table[b'W' as usize] = b'W';
    table[b'K' as usize] = b'M';
    table[b'M' as usize] = b'K';
    table[b'B' as usize] = b'V';
    table[b'V' as usize] = b'B';
    table[b'D' as usize] = b'H';
    table[b'H' as usize] = b'D';
    table[b'N' as usize] = b'N';
    // Lowercase (case preserved).
    table[b'a' as usize] = b't';
    table[b't' as usize] = b'a';
    table[b'c' as usize] = b'g';
    table[b'g' as usize] = b'c';
    table[b'u' as usize] = b'a';
    table[b'r' as usize] = b'y';
    table[b'y' as usize] = b'r';
    table[b's' as usize] = b's';
    table[b'w' as usize] = b'w';
    table[b'k' as usize] = b'm';
    table[b'm' as usize] = b'k';
    table[b'b' as usize] = b'v';
    table[b'v' as usize] = b'b';
    table[b'd' as usize] = b'h';
    table[b'h' as usize] = b'd';
    table[b'n' as usize] = b'n';
    table
};

/// Complements a single DNA base (IUPAC-aware, case-preserving).
///
/// Watson–Crick + IUPAC complement via [`COMPLEMENT`]: `A<->T`, `C<->G`,
/// `R<->Y`, etc. Case is preserved (`a->t`) and unknown bytes pass through
/// unchanged. See [`COMPLEMENT`] for the exact contract and rationale.
#[inline]
#[must_use]
pub const fn complement_base(base: u8) -> u8 {
    COMPLEMENT[base as usize]
}

/// Reverse complements a DNA sequence (IUPAC-aware, case-preserving).
///
/// Reverses the sequence and complements each base via [`complement_base`]:
/// full IUPAC (`A<->T`, `C<->G`, `R<->Y`, …), case preserved, unknown bytes
/// unchanged.
///
/// # Examples
///
/// ```
/// use fgumi_dna::dna::reverse_complement;
///
/// assert_eq!(reverse_complement(b"ACGT"), b"ACGT".to_vec());
/// assert_eq!(reverse_complement(b"AAAA"), b"TTTT".to_vec());
/// assert_eq!(reverse_complement(b"ACGTN"), b"NACGT".to_vec());
/// assert_eq!(reverse_complement(b"NRG"), b"CYN".to_vec()); // IUPAC-aware
/// assert_eq!(reverse_complement(b"acgt"), b"acgt".to_vec()); // case preserved
/// ```
#[must_use]
pub fn reverse_complement(seq: &[u8]) -> Vec<u8> {
    seq.iter().rev().map(|&base| complement_base(base)).collect()
}

/// Reverse complements a DNA string (IUPAC-aware, case-preserving).
///
/// Reverses the string and complements each ASCII byte via [`complement_base`]
/// (case preserved, IUPAC-aware, unknown bytes unchanged). Non-ASCII characters
/// pass through unchanged.
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
        .map(|c| if c.is_ascii() { complement_base(c as u8) as char } else { c })
        .collect()
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_complement_base_iupac_and_case() {
        // fgbio Sequences.complement: full IUPAC table, case-preserving (R2-SEQ-01/02).
        assert_eq!(complement_base(b'R'), b'Y');
        assert_eq!(complement_base(b'Y'), b'R');
        assert_eq!(complement_base(b'S'), b'S');
        assert_eq!(complement_base(b'W'), b'W');
        assert_eq!(complement_base(b'K'), b'M');
        assert_eq!(complement_base(b'M'), b'K');
        assert_eq!(complement_base(b'B'), b'V');
        assert_eq!(complement_base(b'V'), b'B');
        assert_eq!(complement_base(b'D'), b'H');
        assert_eq!(complement_base(b'H'), b'D');
        assert_eq!(complement_base(b'U'), b'A');
        // Lowercase IUPAC/ACGT preserve case.
        assert_eq!(complement_base(b'a'), b't');
        assert_eq!(complement_base(b'u'), b'a'); // RNA uracil, case preserved
        assert_eq!(complement_base(b'r'), b'y');
        assert_eq!(complement_base(b'k'), b'm');
        // Unknown bytes pass through unchanged (R2-SEQ-03: documented passthrough).
        assert_eq!(complement_base(b'.'), b'.');
        assert_eq!(complement_base(b'X'), b'X');
    }

    #[test]
    fn test_reverse_complement_fgbio_parity() {
        // fgbio SequencesTest.scala:160,161 — the headline parity cases (R2-SEQ-01/02).
        assert_eq!(reverse_complement(b"NRG"), b"CYN".to_vec());
        assert_eq!(reverse_complement(b"acacNNNN"), b"NNNNgtgt".to_vec());
        // RNA uracil complements to A (matches fgbio revcomp of 'r'-prefixed UMIs).
        assert_eq!(reverse_complement(b"AAU"), b"ATT".to_vec());
        assert_eq!(reverse_complement(b"U"), b"A".to_vec());
    }

    #[test]
    fn test_complement_base() {
        // Uppercase bases
        assert_eq!(complement_base(b'A'), b'T');
        assert_eq!(complement_base(b'T'), b'A');
        assert_eq!(complement_base(b'C'), b'G');
        assert_eq!(complement_base(b'G'), b'C');

        // Lowercase preserves case (matches fgbio)
        assert_eq!(complement_base(b'a'), b't');
        assert_eq!(complement_base(b't'), b'a');
        assert_eq!(complement_base(b'c'), b'g');
        assert_eq!(complement_base(b'g'), b'c');

        // N bases preserve case (fgbio uses lowercase 'n' for padding)
        assert_eq!(complement_base(b'N'), b'N');
        assert_eq!(complement_base(b'n'), b'n');

        // Non-base bytes pass through unchanged
        for c in [b'.', b'-', b'*', b'0'] {
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

        // Case preserved (matches fgbio)
        assert_eq!(reverse_complement(b"acgt"), b"acgt".to_vec());
        assert_eq!(reverse_complement(b"AcGt"), b"aCgT".to_vec());

        // Palindromic sequences
        assert_eq!(reverse_complement(b"ACGT"), b"ACGT".to_vec());
        assert_eq!(reverse_complement(b"GAATTC"), b"GAATTC".to_vec());

        // Double operation returns original (uppercase)
        let seq = b"ACGTACGT";
        assert_eq!(reverse_complement(&reverse_complement(seq)), seq.to_vec());
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
