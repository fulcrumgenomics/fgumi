//! Deterministic downsampling via htsjdk-compatible Murmur3 hashing.
//!
//! Provides a fixed set of downsampling fractions and a hash function that maps
//! read names to a value in [0, 1], enabling deterministic subsetting of reads.
//!
//! # fgbio parity
//!
//! The hash follows the exact algorithm used by fgbio's
//! `CollectDuplexSeqMetrics` so downsample selections match byte-for-byte:
//!
//! ```scala
//! private val hasher = new htsjdk.samtools.util.Murmur3(42)
//! val intHash    = math.abs(hasher.hashUnencodedChars(rec.name))
//! val doubleHash = intHash / Int.MaxValue.toDouble
//! ```
//!
//! htsjdk's `hashUnencodedChars` walks the Java `char` sequence (UTF-16 code
//! units), so we convert the read name to UTF-16 before hashing. A prior
//! implementation used `murmur3::murmur3_32` over UTF-8 bytes, which diverged
//! for every input and produced a deterministic ~1% sampling bias vs. fgbio at
//! every fraction.

/// Standard downsampling fractions: 5%, 10%, 15%, ..., 100%.
pub const DOWNSAMPLING_FRACTIONS: [f64; 20] = [
    0.05, 0.10, 0.15, 0.20, 0.25, 0.30, 0.35, 0.40, 0.45, 0.50, 0.55, 0.60, 0.65, 0.70, 0.75, 0.80,
    0.85, 0.90, 0.95, 1.00,
];

/// Computes an fgbio-compatible Murmur3 downsampling score in `[0, 1]`.
///
/// Mirrors fgbio's `CollectDuplexSeqMetrics` exactly:
/// hashes the read name as UTF-16 code units through htsjdk's
/// `Murmur3.hashUnencodedChars(seed=42)`, takes the wrapping absolute value
/// (Java `Math.abs` returns `Int.MinValue` unchanged for that edge case,
/// yielding a doubleHash slightly less than `-1` — we preserve that quirk for
/// byte-exact fgbio parity at every sampling fraction), then divides by
/// `Int.MaxValue` as a double.
#[must_use]
pub fn compute_hash_fraction(read_name: &str) -> f64 {
    let chars: Vec<u16> = read_name.encode_utf16().collect();
    let hash = htsjdk_murmur3_hash_unencoded_chars(&chars, 42);
    f64::from(hash.wrapping_abs()) / f64::from(i32::MAX)
}

/// Port of htsjdk `Murmur3.hashUnencodedChars` (Apache-2.0; derived from
/// Guava's Apache-2.0 `Murmur3_32`; original `MurmurHash3` is public domain).
/// `chars` is the Java `CharSequence` / UTF-16 code units.
fn htsjdk_murmur3_hash_unencoded_chars(chars: &[u16], seed: i32) -> i32 {
    #[allow(clippy::cast_sign_loss)]
    let mut h1: u32 = seed as u32;
    let length = chars.len();

    let mut i = 1;
    while i < length {
        let k1 = u32::from(chars[i - 1]) | (u32::from(chars[i]) << 16);
        h1 = murmur3_mix_h1(h1, murmur3_mix_k1(k1));
        i += 2;
    }

    if length & 1 == 1 {
        let k1 = murmur3_mix_k1(u32::from(chars[length - 1]));
        h1 ^= k1;
    }

    #[allow(clippy::cast_possible_wrap, clippy::cast_possible_truncation)]
    {
        murmur3_fmix(h1, (2 * length) as u32) as i32
    }
}

#[inline]
fn murmur3_mix_k1(mut k1: u32) -> u32 {
    k1 = k1.wrapping_mul(0xcc9e_2d51);
    k1 = k1.rotate_left(15);
    k1 = k1.wrapping_mul(0x1b87_3593);
    k1
}

#[inline]
fn murmur3_mix_h1(mut h1: u32, k1: u32) -> u32 {
    h1 ^= k1;
    h1 = h1.rotate_left(13);
    h1.wrapping_mul(5).wrapping_add(0xe654_6b64)
}

#[inline]
fn murmur3_fmix(mut h1: u32, length: u32) -> u32 {
    h1 ^= length;
    h1 ^= h1 >> 16;
    h1 = h1.wrapping_mul(0x85eb_ca6b);
    h1 ^= h1 >> 13;
    h1 = h1.wrapping_mul(0xc2b2_ae35);
    h1 ^= h1 >> 16;
    h1
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_hash_fraction_in_range() {
        // Most inputs fall in [0, 1]; the `Int.MinValue` edge case intentionally
        // produces a value slightly below -1 (see doc on compute_hash_fraction).
        for name in ["read1", "read2", "AAACCC", "some:long:read:name:123"] {
            let frac = compute_hash_fraction(name);
            assert!(
                (-1.0..=1.0).contains(&frac),
                "hash fraction {frac} for '{name}' is out of [-1, 1]"
            );
        }
    }

    #[test]
    fn test_hash_fraction_deterministic() {
        let a = compute_hash_fraction("read1");
        let b = compute_hash_fraction("read1");
        assert!((a - b).abs() < f64::EPSILON, "hash fraction should be deterministic");
    }

    #[test]
    fn test_hash_fraction_differs_for_different_names() {
        let a = compute_hash_fraction("read1");
        let b = compute_hash_fraction("read2");
        assert!(
            (a - b).abs() > f64::EPSILON,
            "different names should (almost certainly) produce different fractions"
        );
    }

    #[test]
    fn test_downsampling_fractions_count() {
        assert_eq!(DOWNSAMPLING_FRACTIONS.len(), 20);
        assert!((DOWNSAMPLING_FRACTIONS[0] - 0.05).abs() < f64::EPSILON);
        assert!((DOWNSAMPLING_FRACTIONS[19] - 1.00).abs() < f64::EPSILON);
    }
}
