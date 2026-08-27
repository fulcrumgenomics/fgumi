//! Read-name UMI normalization shared across commands.
//!
//! Both `extract` (FASTQ read names) and `copy-umi` (BAM read names) take the raw
//! UMI text from a read-name field and normalize it into the canonical `RX`-tag
//! form. This module holds that normalization so the two commands stay in lock-step
//! and byte-for-byte faithful to fgbio's `Umis.extractUmisFromReadName`
//! (`Umis.scala:85-126`). Locating the raw UMI within the read name (which field,
//! how strict about field count) is the caller's job and differs between commands;
//! only the shared normalization lives here.

use anyhow::{Result, bail};
use fgumi_dna::dna::reverse_complement;

/// Prefix marking a read-name UMI segment that must be reverse-complemented,
/// matching fgbio `Umis.RevcompPrefix` (`Umis.scala:33`).
pub(crate) const REVCOMP_PREFIX: u8 = b'r';

/// Delimiter separating dual UMIs in a read name, translated to `-` on extraction,
/// matching fgbio's `umiDelimiter` default (`Umis.scala:88`).
pub(crate) const UMI_DELIMITER: u8 = b'+';

/// Returns true if `b` is a valid SAM UMI character (`ACGTN-`), matching fgbio's
/// `Umis.isValidUmiCharacter` (`Umis.scala:128-130`).
///
/// This is intentionally distinct from [`fgumi_umi::validate_umi`], which
/// implements the more permissive `GroupReadsByUmi` counting rule (only rejects
/// upper-case `N`); the read-name extraction path uses fgbio's strict alphabet.
fn is_valid_umi_char(b: u8) -> bool {
    matches!(b, b'A' | b'C' | b'G' | b'T' | b'N' | b'-')
}

/// Normalizes a UMI extracted from the last field of a read name, matching
/// fgbio's `Umis.extractUmisFromReadName` (strict mode, as called by
/// `FastqToBam`; `Umis.scala:85-126`).
///
/// With `reverseComplementPrefixedUmis` defaulting to true, fgbio:
/// 1. reverse-complements any `r`-prefixed segment (EXT-04); when the UMI is
///    `+`-delimited, only the `r`-prefixed segments are reverse-complemented,
/// 2. translates the `+` dual-UMI delimiter to `-`,
/// 3. upper-cases the result (EXT-03), and
/// 4. throws if any character is outside `ACGTN-` (EXT-01, strict mode).
///
/// Reverse-complementing reuses [`fgumi_dna::dna::reverse_complement`], which
/// complements `ACGT`, maps RNA `U`→`A`, and preserves `N` — matching fgbio's
/// `Sequences.complement` over the bases a read-name UMI may contain. Any other
/// character is passed through and then rejected by the `ACGTN-` validation
/// below, mirroring fgbio's strict-mode rejection.
///
/// # Errors
/// Returns an error if the normalized UMI contains a character outside `ACGTN-`,
/// mirroring fgbio's strict-mode `IllegalArgumentException`.
pub(crate) fn normalize_read_name_umi(raw: &[u8]) -> Result<Vec<u8>> {
    // fgbio keys the transform on whether an 'r' appears anywhere and whether a
    // '+' delimiter appears at an index > 0 (a leading '+' is not a delimiter).
    let has_revcomp_prefix = raw.contains(&REVCOMP_PREFIX);
    let has_delimiter = raw.iter().position(|&b| b == UMI_DELIMITER).is_some_and(|idx| idx > 0);

    let transformed: Vec<u8> = match (has_revcomp_prefix, has_delimiter) {
        // Reverse-complement each `r`-prefixed segment, join with '-'.
        (true, true) => {
            let mut out = Vec::with_capacity(raw.len());
            for (i, segment) in raw.split(|&b| b == UMI_DELIMITER).enumerate() {
                if i > 0 {
                    out.push(b'-');
                }
                match segment.split_first() {
                    Some((&REVCOMP_PREFIX, rest)) => out.extend(reverse_complement(rest)),
                    // `stripPrefix('r')` is a no-op when the segment is not
                    // prefixed with 'r'; the segment is copied verbatim.
                    _ => out.extend_from_slice(segment),
                }
            }
            out
        }
        // Single UMI: reverse-complement after stripping a leading 'r'.
        (true, false) => {
            let stripped = raw.strip_prefix(&[REVCOMP_PREFIX][..]).unwrap_or(raw);
            reverse_complement(stripped)
        }
        // No revcomp: just translate the '+' delimiter to '-'.
        (false, true) => raw.iter().map(|&b| if b == UMI_DELIMITER { b'-' } else { b }).collect(),
        (false, false) => raw.to_vec(),
    };

    // fgbio upper-cases the UMI after any reverse-complementing.
    let umi: Vec<u8> = transformed.iter().map(u8::to_ascii_uppercase).collect();

    // Strict mode: reject any character outside the SAM UMI alphabet.
    if let Some(&bad) = umi.iter().find(|&&b| !is_valid_umi_char(b)) {
        bail!(
            "Invalid UMI '{}' extracted from read name (illegal character '{}')",
            String::from_utf8_lossy(&umi),
            bad as char,
        );
    }

    Ok(umi)
}

#[cfg(test)]
mod tests {
    use super::*;
    use rstest::rstest;

    #[rstest]
    #[case::single(b"ACGT", b"ACGT")]
    #[case::lowercase_upcased(b"acgt", b"ACGT")]
    #[case::plus_to_hyphen(b"ACGT+CAGA", b"ACGT-CAGA")]
    #[case::already_hyphen(b"ACGT-CAGA", b"ACGT-CAGA")]
    #[case::r_single_revcomped(b"rAAAA", b"TTTT")]
    #[case::r_lowercase_revcomped(b"racgt", b"ACGT")]
    #[case::r_only_prefixed_segment(b"rAAAA+CCCC", b"TTTT-CCCC")]
    #[case::r_both_segments(b"rAAAA+rCCCC", b"TTTT-GGGG")]
    #[case::uracil_maps_to_a(b"rUACG", b"CGTA")]
    #[case::trailing_plus(b"AAAA+", b"AAAA-")]
    fn normalizes(#[case] raw: &[u8], #[case] expected: &[u8]) {
        assert_eq!(normalize_read_name_umi(raw).unwrap(), expected);
    }

    #[rstest]
    #[case::illegal_k(b"CCKC")]
    #[case::illegal_x(b"XYZ")]
    #[case::illegal_in_dual(b"ACGT+CCKC")]
    fn rejects_illegal_characters(#[case] raw: &[u8]) {
        assert!(normalize_read_name_umi(raw).is_err());
    }

    #[rstest]
    #[case::a_is_valid(b'A', true)]
    #[case::n_is_valid(b'N', true)]
    #[case::hyphen_is_valid(b'-', true)]
    #[case::lowercase_invalid(b'a', false)]
    #[case::plus_invalid(b'+', false)]
    #[case::k_invalid(b'K', false)]
    fn validates_umi_chars(#[case] b: u8, #[case] expected: bool) {
        assert_eq!(is_valid_umi_char(b), expected);
    }
}
