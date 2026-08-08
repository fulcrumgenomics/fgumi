//! fgbio-parity hashing of BAM read names.

/// The Murmur3-32 body htsjdk applies to a Java `CharSequence`, over `len` UTF-16
/// code units read through `at`. Indexed rather than iterator-based because the
/// algorithm consumes code units pairwise.
#[expect(
    clippy::cast_sign_loss,
    clippy::cast_possible_wrap,
    clippy::cast_possible_truncation,
    reason = "lossless bit reinterpretations required for Java Murmur3 parity: the seed and \
              final hash are Java `int`s reinterpreted as `u32` for the mixing arithmetic, and \
              `len` is a read name / char count that never approaches u32::MAX"
)]
fn murmur3_hash_indexed(len: usize, at: impl Fn(usize) -> u16, seed: i32) -> i32 {
    let mut h1: u32 = seed as u32;

    let mut i = 1;
    while i < len {
        let k1 = u32::from(at(i - 1)) | (u32::from(at(i)) << 16);
        h1 = murmur3_mix_h1(h1, murmur3_mix_k1(k1));
        i += 2;
    }

    if len & 1 == 1 {
        let k1 = murmur3_mix_k1(u32::from(at(len - 1)));
        h1 ^= k1;
    }

    murmur3_fmix(h1, (2 * len) as u32) as i32
}

/// Port of htsjdk `Murmur3.hashUnencodedChars` (Apache-2.0; derived from
/// Guava's Apache-2.0 `Murmur3_32`; original `MurmurHash3` is public domain).
/// `chars` is the Java `CharSequence` / UTF-16 code units.
///
/// The returned hash is **signed**, and must be compared as such: htsjdk hashes
/// by Java `int`, so about half of all inputs hash negative. Comparing as `u32`
/// silently reorders every one of them.
#[must_use]
pub fn htsjdk_murmur3_hash_unencoded_chars(chars: &[u16], seed: i32) -> i32 {
    murmur3_hash_indexed(chars.len(), |i| chars[i], seed)
}

/// Ranks a read name for fgbio-compatible downsampling.
///
/// Mirrors fgbio's `VanillaUmiConsensusCaller.downsampleRank`: the rank is
/// `new Murmur3(42).hashUnencodedChars(readName)`, and downsampling keeps the
/// lowest-ranking reads. Because both ends of a template share a read name they
/// share a rank, so a cap retains or discards a template as a unit.
///
/// The returned rank is **signed**, and must be compared as such: fgbio ranks by
/// Scala `Int`, so about half of all read names rank negative. Comparing as
/// `u32` silently reorders every one of them.
///
/// htsjdk hashes the UTF-16 code units of a Java `String`. The SAM spec restricts
/// read names to printable ASCII (`[!-?A-~]{1,254}`), for which widening each
/// byte to a `u16` is exactly that code-unit sequence — so this avoids a UTF-8
/// decode on a per-read hot path. A non-ASCII byte would be a spec violation;
/// it still hashes deterministically, just not identically to htsjdk on the same
/// bytes.
#[must_use]
pub fn fgbio_read_name_rank(name: &[u8]) -> i32 {
    murmur3_hash_indexed(name.len(), |i| u16::from(name[i]), 42)
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
    use rstest::rstest;

    /// Reference values captured from htsjdk's `Murmur3(42).hashUnencodedChars(s)`
    /// by running htsjdk itself, not fgumi's port. These are the independent
    /// oracle: if the port drifts, these fail.
    #[rstest]
    #[case::q0("q0", -593_808_727)]
    #[case::q2("q2", -974_105_965)]
    #[case::q3("q3", -1_135_430_185)]
    #[case::q9("q9", 98_042_550)]
    #[case::read0("read0", 916_970_908)]
    #[case::read5("read5", -1_573_193_749)]
    fn rank_matches_htsjdk_reference_vectors(#[case] name: &str, #[case] expected: i32) {
        assert_eq!(fgbio_read_name_rank(name.as_bytes()), expected, "rank mismatch on {name:?}");
    }

    /// The ASCII fast path must agree with decoding to UTF-16 first, which is
    /// what htsjdk does to a Java `String`.
    #[rstest]
    #[case::short("A")]
    #[case::typical("H0164ALXX140820:2:1101:10003:23260")]
    #[case::odd_length("abc")]
    #[case::even_length("abcd")]
    fn ascii_fast_path_matches_utf16_widening(#[case] name: &str) {
        let widened: Vec<u16> = name.encode_utf16().collect();
        assert_eq!(
            fgbio_read_name_rank(name.as_bytes()),
            htsjdk_murmur3_hash_unencoded_chars(&widened, 42)
        );
    }

    /// Mates share a read name and therefore a rank — this is what makes both
    /// ends of a template survive or die together.
    #[test]
    fn identical_names_rank_identically() {
        assert_eq!(fgbio_read_name_rank(b"frag:1"), fgbio_read_name_rank(b"frag:1"));
    }
}
