use crate::fields::{l_seq, qual_offset, seq_offset};

/// BAM 4-bit base encoding -> ASCII lookup table.
///
/// Index is the 4-bit nibble value: 0=`=`, 1=`A`, 2=`C`, 4=`G`, 8=`T`, 15=`N`.
pub const BAM_BASE_TO_ASCII: [u8; 16] = [
    b'=', b'A', b'C', b'M', b'G', b'R', b'S', b'V', b'T', b'W', b'Y', b'H', b'K', b'D', b'B', b'N',
];

/// Lookup table mapping ASCII base characters to 4-bit BAM codes.
///
/// Duplicated from `vendored/bam_codec/encoder/sequence.rs` to keep
/// `bam_fields` self-contained. The table is `const` so there is zero
/// runtime cost.
const SEQ_CODES: [u8; 256] = build_seq_codes();

const fn build_seq_codes() -> [u8; 256] {
    // SAM spec §4.2.3: =ACMGRSVTWYHKDBN -> 0..15
    // BASES has exactly 16 elements so `i` always fits in u8.
    const BASES: [u8; 16] = *b"=ACMGRSVTWYHKDBN";
    const N: u8 = 0x0F;
    let mut codes = [N; 256];
    let mut i: u8 = 0;
    while (i as usize) < BASES.len() {
        let base = BASES[i as usize];
        let code = i;
        codes[base as usize] = code;
        codes[base.to_ascii_lowercase() as usize] = code;
        i += 1;
    }
    codes
}

/// Decode table: 4-bit BAM base code -> ASCII base.
const BASE_DECODE: [u8; 16] = *b"=ACMGRSVTWYHKDBN";

/// Extract a 4-bit base from packed sequence data.
///
/// BAM packs two bases per byte: high nibble = even index, low nibble = odd index.
/// Returns the 4-bit encoded base (1=A, 2=C, 4=G, 8=T, 15=N).
#[inline]
#[must_use]
pub fn get_base(bam: &[u8], seq_off: usize, position: usize) -> u8 {
    let byte = bam[seq_off + position / 2];
    if position.is_multiple_of(2) { byte >> 4 } else { byte & 0xF }
}

/// Set a single base at a position in packed 4-bit BAM sequence data.
///
/// Converts the given ASCII base to a 4-bit BAM code and writes it in place.
/// Unknown bases are encoded as N (0xF).
#[inline]
pub fn set_base(bam: &mut [u8], seq_off: usize, position: usize, base: u8) {
    let encoded = SEQ_CODES[base as usize];
    let byte_idx = seq_off + position / 2;
    if position.is_multiple_of(2) {
        bam[byte_idx] = (encoded << 4) | (bam[byte_idx] & 0x0F);
    } else {
        bam[byte_idx] = (bam[byte_idx] & 0xF0) | encoded;
    }
}

/// Set the 4-bit base at a position to N (0xF) in packed sequence data.
#[inline]
pub fn mask_base(bam: &mut [u8], seq_off: usize, position: usize) {
    let byte_idx = seq_off + position / 2;
    if position.is_multiple_of(2) {
        bam[byte_idx] = (bam[byte_idx] & 0x0F) | 0xF0; // high nibble
    } else {
        bam[byte_idx] = (bam[byte_idx] & 0xF0) | 0x0F; // low nibble
    }
}

/// Check if a 4-bit encoded base at position is N (0xF).
#[inline]
#[must_use]
pub fn is_base_n(bam: &[u8], seq_off: usize, position: usize) -> bool {
    get_base(bam, seq_off, position) == 0xF
}

/// Bulk-extract the full sequence from a BAM record as ASCII bases.
///
/// Decodes the packed 4-bit sequence data into a `Vec<u8>` of ASCII bases.
#[must_use]
pub fn extract_sequence(bam: &[u8]) -> Vec<u8> {
    let l = l_seq(bam) as usize;
    let off = seq_offset(bam);
    let mut bases = Vec::with_capacity(l);
    for i in 0..l {
        let code = get_base(bam, off, i);
        bases.push(BASE_DECODE[code as usize]);
    }
    bases
}

/// Pack ASCII bases into BAM 4-bit-per-base format, appending to `dst`.
///
/// Uses 16-base chunked processing for cache efficiency, matching the
/// htslib/vendored encoder strategy. When `l_seq` is odd the bottom
/// 4 bits of the last byte are zero-padded per the SAM spec.
const PACK_CHUNK: usize = 16;

#[inline]
pub fn pack_sequence_into(dst: &mut Vec<u8>, bases: &[u8]) {
    if bases.is_empty() {
        return;
    }
    let packed_len = bases.len().div_ceil(2);
    dst.reserve(packed_len);
    let mut chunks = bases.chunks_exact(PACK_CHUNK);
    for chunk in chunks.by_ref() {
        dst.push((SEQ_CODES[chunk[0] as usize] << 4) | SEQ_CODES[chunk[1] as usize]);
        dst.push((SEQ_CODES[chunk[2] as usize] << 4) | SEQ_CODES[chunk[3] as usize]);
        dst.push((SEQ_CODES[chunk[4] as usize] << 4) | SEQ_CODES[chunk[5] as usize]);
        dst.push((SEQ_CODES[chunk[6] as usize] << 4) | SEQ_CODES[chunk[7] as usize]);
        dst.push((SEQ_CODES[chunk[8] as usize] << 4) | SEQ_CODES[chunk[9] as usize]);
        dst.push((SEQ_CODES[chunk[10] as usize] << 4) | SEQ_CODES[chunk[11] as usize]);
        dst.push((SEQ_CODES[chunk[12] as usize] << 4) | SEQ_CODES[chunk[13] as usize]);
        dst.push((SEQ_CODES[chunk[14] as usize] << 4) | SEQ_CODES[chunk[15] as usize]);
    }

    let remainder = chunks.remainder();
    let mut pairs = remainder.chunks_exact(2);
    for pair in pairs.by_ref() {
        dst.push((SEQ_CODES[pair[0] as usize] << 4) | SEQ_CODES[pair[1] as usize]);
    }
    if let Some(&last) = pairs.remainder().first() {
        dst.push(SEQ_CODES[last as usize] << 4);
    }
}

/// Extract a quality score at a given position.
#[inline]
#[must_use]
pub fn get_qual(bam: &[u8], qual_off: usize, position: usize) -> u8 {
    bam[qual_off + position]
}

/// Set the quality score at a given position.
#[inline]
pub fn set_qual(bam: &mut [u8], qual_off: usize, position: usize, value: u8) {
    bam[qual_off + position] = value;
}

/// Zero-copy access to quality scores in a BAM record.
///
/// Returns a slice of the raw Phred quality scores (not Phred+33).
#[inline]
#[must_use]
pub fn quality_scores_slice(bam: &[u8]) -> &[u8] {
    let l = l_seq(bam) as usize;
    let off = qual_offset(bam);
    &bam[off..off + l]
}

/// Mutable zero-copy access to quality scores in a BAM record.
#[inline]
pub fn quality_scores_slice_mut(bam: &mut [u8]) -> &mut [u8] {
    let l = l_seq(bam) as usize;
    let off = qual_offset(bam);
    &mut bam[off..off + l]
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::testutil::*;

    // ========================================================================
    // get_base and mask_base tests
    // ========================================================================

    #[test]
    fn test_get_base_and_mask() {
        // Build a record with seq_len=4, fill seq bytes manually
        let mut rec = make_bam_bytes(0, 0, 0, b"rea", &[], 4, -1, -1, &[]);
        let so = seq_offset(&rec);
        // Packed: A=1, C=2 -> byte = 0x12; G=4, T=8 -> byte = 0x48
        rec[so] = 0x12; // bases 0,1 = A,C
        rec[so + 1] = 0x48; // bases 2,3 = G,T

        assert_eq!(get_base(&rec, so, 0), 1); // A
        assert_eq!(get_base(&rec, so, 1), 2); // C
        assert_eq!(get_base(&rec, so, 2), 4); // G
        assert_eq!(get_base(&rec, so, 3), 8); // T

        // Mask base 1 (C -> N)
        mask_base(&mut rec, so, 1);
        assert_eq!(get_base(&rec, so, 1), 0xF); // N
        assert_eq!(get_base(&rec, so, 0), 1); // A unchanged

        // Mask base 2 (G -> N, high nibble)
        mask_base(&mut rec, so, 2);
        assert_eq!(get_base(&rec, so, 2), 0xF); // N
        assert_eq!(get_base(&rec, so, 3), 8); // T unchanged
    }

    #[test]
    fn test_mask_base_even_position() {
        let mut rec = make_bam_bytes(0, 0, 0, b"rea", &[], 4, -1, -1, &[]);
        let so = seq_offset(&rec);
        rec[so] = 0x12; // A(1),C(2)
        rec[so + 1] = 0x48; // G(4),T(8)

        // Mask position 0 (even, high nibble)
        mask_base(&mut rec, so, 0);
        assert_eq!(get_base(&rec, so, 0), 0xF); // N
        assert_eq!(get_base(&rec, so, 1), 2); // C preserved

        // Mask position 3 (odd, low nibble)
        mask_base(&mut rec, so, 3);
        assert_eq!(get_base(&rec, so, 3), 0xF); // N
        assert_eq!(get_base(&rec, so, 2), 4); // G preserved
    }

    #[test]
    fn test_mask_base_sets_n() {
        let mut rec = make_bam_bytes(0, 0, 0, b"rd", &[], 4, -1, -1, &[]);
        let so = seq_offset(&rec);
        // Set all bases to ACGT
        set_base(&mut rec, so, 0, b'A');
        set_base(&mut rec, so, 1, b'C');
        set_base(&mut rec, so, 2, b'G');
        set_base(&mut rec, so, 3, b'T');
        // Mask positions 1 and 3
        mask_base(&mut rec, so, 1);
        mask_base(&mut rec, so, 3);
        assert_eq!(get_base(&rec, so, 0), 1); // A preserved
        assert_eq!(get_base(&rec, so, 1), 15); // N
        assert_eq!(get_base(&rec, so, 2), 4); // G preserved
        assert_eq!(get_base(&rec, so, 3), 15); // N
    }

    // ========================================================================
    // get_qual and set_qual tests
    // ========================================================================

    #[test]
    fn test_get_qual_and_set_qual() {
        let mut rec = make_bam_bytes(0, 0, 0, b"rea", &[], 4, -1, -1, &[]);
        let qo = qual_offset(&rec);
        rec[qo] = 30;
        rec[qo + 1] = 40;
        assert_eq!(get_qual(&rec, qo, 0), 30);
        assert_eq!(get_qual(&rec, qo, 1), 40);

        set_qual(&mut rec, qo, 0, 10);
        assert_eq!(get_qual(&rec, qo, 0), 10);
        assert_eq!(get_qual(&rec, qo, 1), 40); // unchanged
    }

    #[test]
    fn test_set_qual_boundary_values() {
        let mut rec = make_bam_bytes(0, 0, 0, b"rd", &[], 3, -1, -1, &[]);
        let qo = qual_offset(&rec);
        set_qual(&mut rec, qo, 0, 0); // min
        set_qual(&mut rec, qo, 1, 93); // max typical
        set_qual(&mut rec, qo, 2, 255); // max possible
        assert_eq!(get_qual(&rec, qo, 0), 0);
        assert_eq!(get_qual(&rec, qo, 1), 93);
        assert_eq!(get_qual(&rec, qo, 2), 255);
    }

    // ========================================================================
    // extract_sequence tests
    // ========================================================================

    #[test]
    fn test_extract_sequence() {
        let mut rec = make_bam_bytes(0, 0, 0, b"rea", &[], 4, -1, -1, &[]);
        let so = seq_offset(&rec);
        rec[so] = 0x12; // A, C
        rec[so + 1] = 0x48; // G, T
        let seq = extract_sequence(&rec);
        assert_eq!(seq, b"ACGT");
    }

    #[test]
    fn test_extract_sequence_odd_length() {
        // Odd-length seq: the last nibble in the packed byte is in the high nibble
        let mut rec = make_bam_bytes(0, 0, 0, b"rd", &[], 3, -1, -1, &[]);
        let so = seq_offset(&rec);
        rec[so] = 0x12; // A(1), C(2)
        rec[so + 1] = 0x40; // G(4), padding(0)
        let seq = extract_sequence(&rec);
        assert_eq!(seq, b"ACG");
    }

    #[test]
    fn test_extract_sequence_all_n() {
        let mut rec = make_bam_bytes(0, 0, 0, b"rd", &[], 4, -1, -1, &[]);
        let so = seq_offset(&rec);
        rec[so] = 0xFF; // N(15), N(15)
        rec[so + 1] = 0xFF; // N(15), N(15)
        let seq = extract_sequence(&rec);
        assert_eq!(seq, b"NNNN");
    }

    #[test]
    fn test_extract_sequence_single_base() {
        let mut rec = make_bam_bytes(0, 0, 0, b"rd", &[], 1, -1, -1, &[]);
        let so = seq_offset(&rec);
        rec[so] = 0x80; // T(8)
        let seq = extract_sequence(&rec);
        assert_eq!(seq, b"T");
    }

    // ========================================================================
    // quality_scores_slice tests
    // ========================================================================

    #[test]
    fn test_quality_scores_slice() {
        let mut rec = make_bam_bytes(0, 0, 0, b"rea", &[], 4, -1, -1, &[]);
        let qo = qual_offset(&rec);
        rec[qo] = 30;
        rec[qo + 1] = 25;
        rec[qo + 2] = 20;
        rec[qo + 3] = 15;
        let quals = quality_scores_slice(&rec);
        assert_eq!(quals, &[30, 25, 20, 15]);
    }

    #[test]
    fn test_quality_scores_slice_high_qualities() {
        let mut rec = make_bam_bytes(0, 0, 0, b"rd", &[], 3, -1, -1, &[]);
        let qo = qual_offset(&rec);
        rec[qo] = 40;
        rec[qo + 1] = 93; // max Phred
        rec[qo + 2] = 0; // min Phred
        let quals = quality_scores_slice(&rec);
        assert_eq!(quals, &[40, 93, 0]);
    }

    #[test]
    fn test_quality_scores_slice_mut() {
        let mut rec = make_bam_bytes(0, 0, 0, b"rea", &[], 4, -1, -1, &[]);
        let quals = quality_scores_slice_mut(&mut rec);
        quals[0] = 30;
        quals[1] = 25;
        assert_eq!(quality_scores_slice(&rec), &[30, 25, 0, 0]);
    }

    #[test]
    fn test_quality_scores_slice_mut_modify_all() {
        let mut rec = make_bam_bytes(0, 0, 0, b"rd", &[], 4, -1, -1, &[]);
        let quals = quality_scores_slice_mut(&mut rec);
        for (i, q) in quals.iter_mut().enumerate() {
            *q = u8::try_from(i * 10).unwrap();
        }
        assert_eq!(quality_scores_slice(&rec), &[0, 10, 20, 30]);
    }

    // ========================================================================
    // set_base tests
    // ========================================================================

    #[test]
    fn test_set_base() {
        let mut rec = make_bam_bytes(0, 0, 0, b"rea", &[], 4, -1, -1, &[]);
        let so = seq_offset(&rec);
        // Start with AAAA (all 0x11 pairs)
        rec[so] = 0x11; // A=1, A=1
        rec[so + 1] = 0x11;

        // set_base takes ASCII bases, not codes
        // Set position 0 to T (code 8)
        set_base(&mut rec, so, 0, b'T');
        assert_eq!(get_base(&rec, so, 0), 8); // T
        assert_eq!(get_base(&rec, so, 1), 1); // A unchanged

        // Set position 1 to C (code 2)
        set_base(&mut rec, so, 1, b'C');
        assert_eq!(get_base(&rec, so, 0), 8); // T unchanged
        assert_eq!(get_base(&rec, so, 1), 2); // C

        // Set position 3 to G (code 4) - odd position
        set_base(&mut rec, so, 3, b'G');
        assert_eq!(get_base(&rec, so, 2), 1); // A unchanged
        assert_eq!(get_base(&rec, so, 3), 4); // G
    }

    #[test]
    fn test_set_base_all_bases() {
        // Verify set_base/get_base roundtrip for all standard bases
        let mut rec = make_bam_bytes(0, 0, 0, b"rd", &[], 8, -1, -1, &[]);
        let so = seq_offset(&rec);
        let bases = [b'A', b'C', b'G', b'T', b'N', b'R', b'Y', b'M'];
        let codes: [u8; 8] = [1, 2, 4, 8, 15, 5, 10, 3];
        for (i, &base) in bases.iter().enumerate() {
            set_base(&mut rec, so, i, base);
        }
        for (i, &expected_code) in codes.iter().enumerate() {
            assert_eq!(
                get_base(&rec, so, i),
                expected_code,
                "base {} at pos {i}",
                bases[i] as char
            );
        }
    }

    #[test]
    fn test_set_base_preserves_neighbor() {
        // Ensure setting an even position doesn't corrupt the odd neighbor and vice versa
        let mut rec = make_bam_bytes(0, 0, 0, b"rd", &[], 2, -1, -1, &[]);
        let so = seq_offset(&rec);
        set_base(&mut rec, so, 0, b'A');
        set_base(&mut rec, so, 1, b'T');
        assert_eq!(get_base(&rec, so, 0), 1); // A
        assert_eq!(get_base(&rec, so, 1), 8); // T
        // Now change only position 0
        set_base(&mut rec, so, 0, b'G');
        assert_eq!(get_base(&rec, so, 0), 4); // G
        assert_eq!(get_base(&rec, so, 1), 8); // T still intact
    }

    // ========================================================================
    // pack_sequence_into tests
    // ========================================================================

    #[test]
    fn test_pack_sequence_into_empty() {
        let mut dst = Vec::new();
        pack_sequence_into(&mut dst, b"");
        assert!(dst.is_empty());
    }

    #[test]
    fn test_pack_sequence_into_even() {
        let mut dst = Vec::new();
        pack_sequence_into(&mut dst, b"ACGT");
        // A=1, C=2, G=4, T=8
        assert_eq!(dst, [0x12, 0x48]);
    }

    #[test]
    fn test_pack_sequence_into_odd() {
        let mut dst = Vec::new();
        pack_sequence_into(&mut dst, b"ACG");
        // A=1, C=2, G=4, pad=0
        assert_eq!(dst, [0x12, 0x40]);
    }

    #[test]
    fn test_pack_sequence_into_single_base() {
        let mut dst = Vec::new();
        pack_sequence_into(&mut dst, b"T");
        // T=8, pad=0
        assert_eq!(dst, [0x80]);
    }

    #[test]
    fn test_pack_sequence_into_17_bases() {
        // 17 bases: exercises 16-base chunked path + 1-base remainder
        let mut dst = Vec::new();
        pack_sequence_into(&mut dst, b"ACGTACGTACGTACGTA");
        assert_eq!(dst.len(), 9); // ceil(17/2)
        // First 8 bytes (16 bases chunked)
        assert_eq!(dst[0], 0x12); // AC
        assert_eq!(dst[1], 0x48); // GT
        assert_eq!(dst[2], 0x12); // AC
        assert_eq!(dst[3], 0x48); // GT
        assert_eq!(dst[4], 0x12); // AC
        assert_eq!(dst[5], 0x48); // GT
        assert_eq!(dst[6], 0x12); // AC
        assert_eq!(dst[7], 0x48); // GT
        // Last byte (1 base + padding)
        assert_eq!(dst[8], 0x10); // A=1, pad=0
    }

    #[test]
    fn test_pack_sequence_into_n_bases() {
        let mut dst = Vec::new();
        pack_sequence_into(&mut dst, b"NN");
        assert_eq!(dst, [0xFF]); // N=15, N=15
    }
}
