//! SIMD lexer: classifies 64-byte blocks into newline and ACGT bitmasks
//! with 2-bit DNA encoding.
//!
//! Three implementations are provided:
//! - NEON (aarch64): primary target for Apple Silicon / Graviton
//! - AVX2 (`x86_64`): for Intel/AMD servers
//! - Scalar fallback: for testing and unsupported architectures

use crate::bitmask::FastqBitmask;

// ============================================================================
// Public dispatch
// ============================================================================

/// Lex a 64-byte block, producing only the newline bitmask.
///
/// This is the fast path for record boundary detection. Use [`lex_block_full`]
/// when ACGT classification and 2-bit encoding are also needed.
#[inline]
pub fn lex_block(block: &[u8; 64]) -> u64 {
    #[cfg(target_arch = "aarch64")]
    {
        unsafe { lex_block_newlines_neon(block) }
    }

    #[cfg(target_arch = "x86_64")]
    {
        if is_x86_feature_detected!("avx2") {
            unsafe { lex_block_newlines_avx2(block) }
        } else {
            lex_block_newlines_scalar(block)
        }
    }

    #[cfg(not(any(target_arch = "aarch64", target_arch = "x86_64")))]
    {
        lex_block_newlines_scalar(block)
    }
}

/// Lex a 64-byte block, producing newline bitmask, ACGT bitmask, and 2-bit encoding.
///
/// This is the full classification path for when 2-bit DNA encoding is needed
/// (e.g., for UMI extraction in a fused pipeline).
#[inline]
#[must_use]
pub fn lex_block_full(block: &[u8; 64]) -> FastqBitmask {
    #[cfg(target_arch = "aarch64")]
    {
        unsafe { lex_block_neon(block) }
    }

    #[cfg(target_arch = "x86_64")]
    {
        if is_x86_feature_detected!("avx2") {
            unsafe { lex_block_avx2(block) }
        } else {
            lex_block_scalar(block)
        }
    }

    #[cfg(not(any(target_arch = "aarch64", target_arch = "x86_64")))]
    {
        lex_block_scalar(block)
    }
}

// ============================================================================
// Shared constants
// ============================================================================

/// 2-bit encoding table for fgumi's `BitEnc` mapping: A=0, C=1, G=2, T=3.
///
/// Indexed by `(ascii_byte >> 1) & 0x0F`. Only entries for A/C/G/T (and lowercase)
/// produce valid encodings; all others are garbage (but masked out by `is_acgt`).
///
/// ASCII values:
/// - A=0x41, a=0x61 → (>>1)&0xF = 0x0 → encode as 0
/// - C=0x43, c=0x63 → (>>1)&0xF = 0x1 → encode as 1
/// - G=0x47, g=0x67 → (>>1)&0xF = 0x3 → encode as 2
/// - T=0x54, t=0x74 → (>>1)&0xF = 0xA → encode as 3
const ENCODE_LUT: [u8; 16] = [
    0, // 0x0: A/a
    1, // 0x1: C/c
    0, // 0x2: (unused)
    2, // 0x3: G/g
    0, // 0x4: (unused)
    0, // 0x5: (unused)
    0, // 0x6: (unused)
    0, // 0x7: (unused)
    0, // 0x8: (unused)
    0, // 0x9: (unused)
    3, // 0xA: T/t
    0, // 0xB: (unused)
    0, // 0xC: (unused)
    0, // 0xD: (unused)
    0, // 0xE: (unused)
    0, // 0xF: (unused)
];

// ============================================================================
// Scalar fallback
// ============================================================================

/// Scalar newline-only detection. Used as a fallback on `x86_64` without AVX2
/// and on unsupported architectures, plus as a reference in tests.
#[allow(dead_code)]
#[inline]
fn lex_block_newlines_scalar(block: &[u8; 64]) -> u64 {
    let mut newlines: u64 = 0;
    for (i, &byte) in block.iter().enumerate() {
        if byte == b'\n' {
            newlines |= 1u64 << i;
        }
    }
    newlines
}

/// Scalar full classification. Used as a fallback on `x86_64` without AVX2,
/// on unsupported architectures, and as a reference in tests.
#[allow(dead_code)]
#[inline]
pub fn lex_block_scalar(block: &[u8; 64]) -> FastqBitmask {
    let mut newlines: u64 = 0;
    let mut is_acgt: u64 = 0;
    let mut two_bits: u128 = 0;

    for (i, &byte) in block.iter().enumerate() {
        if byte == b'\n' {
            newlines |= 1u64 << i;
        }
        let upper = byte & 0xDF; // case-insensitive
        if upper == b'A' || upper == b'C' || upper == b'G' || upper == b'T' {
            is_acgt |= 1u64 << i;
            let enc = u128::from(ENCODE_LUT[((byte >> 1) & 0x0F) as usize]);
            two_bits |= enc << (i * 2);
        }
    }
    FastqBitmask { newlines, is_acgt, two_bits }
}

// ============================================================================
// NEON (aarch64)
// ============================================================================

#[cfg(target_arch = "aarch64")]
use std::arch::aarch64::uint8x16_t;

/// Extract a 16-bit mask from a NEON comparison result (each byte is 0x00 or 0xFF).
/// Returns the mask in the low 16 bits of a `u64`.
#[cfg(target_arch = "aarch64")]
#[target_feature(enable = "neon")]
#[inline]
unsafe fn neon_movemask(cmp: uint8x16_t) -> u64 {
    use std::arch::aarch64::{
        vaddv_u8, vcreate_u8, vget_high_u8, vget_low_u8, vmul_u8, vshrq_n_u8,
    };

    // NEON intrinsics are safe to call inside a #[target_feature(enable = "neon")] function
    let shifted = vshrq_n_u8(cmp, 7);
    let weights = vcreate_u8(u64::from_le_bytes([1, 2, 4, 8, 16, 32, 64, 128]));
    let low = vget_low_u8(shifted);
    let high = vget_high_u8(shifted);
    let low_sum = u64::from(vaddv_u8(vmul_u8(low, weights)));
    let high_sum = u64::from(vaddv_u8(vmul_u8(high, weights)));
    low_sum | (high_sum << 8)
}

#[cfg(target_arch = "aarch64")]
/// NEON newline-only detection (fast path).
///
/// # Safety
///
/// Requires NEON support (mandatory on aarch64).
#[target_feature(enable = "neon")]
unsafe fn lex_block_newlines_neon(block: &[u8; 64]) -> u64 {
    use std::arch::aarch64::{vceqq_u8, vdupq_n_u8, vld1q_u8};

    unsafe {
        let newline = vdupq_n_u8(b'\n');
        let mut result: u64 = 0;

        for chunk_idx in 0..4 {
            let v = vld1q_u8(block.as_ptr().add(chunk_idx * 16));
            result |= neon_movemask(vceqq_u8(v, newline)) << (chunk_idx * 16);
        }

        result
    }
}

#[cfg(target_arch = "aarch64")]
/// NEON full classification: newlines, ACGT, and 2-bit encoding.
///
/// # Safety
///
/// Requires NEON support (mandatory on aarch64).
#[target_feature(enable = "neon")]
unsafe fn lex_block_neon(block: &[u8; 64]) -> FastqBitmask {
    use std::arch::aarch64::{
        vandq_u8, vceqq_u8, vdupq_n_u8, vld1q_u8, vorrq_u8, vqtbl1q_u8, vshrq_n_u8,
    };

    /// Extract 2-bit encoding from 16 bytes, returning a u32 (16 bases * 2 bits = 32 bits).
    #[inline]
    #[allow(clippy::cast_possible_truncation)]
    unsafe fn neon_two_bits(v: uint8x16_t, lut: uint8x16_t) -> u32 {
        unsafe {
            let idx = vandq_u8(vshrq_n_u8::<1>(v), vdupq_n_u8(0x0F));
            let encoded = vqtbl1q_u8(lut, idx);

            let bit0_mask = vdupq_n_u8(0x01);
            let bit0 = vandq_u8(encoded, bit0_mask);
            let bit1 = vandq_u8(vshrq_n_u8::<1>(encoded), bit0_mask);

            let bit0_packed = neon_movemask(vceqq_u8(bit0, bit0_mask)) as u32;
            let bit1_packed = neon_movemask(vceqq_u8(bit1, bit0_mask)) as u32;

            interleave_bits(bit0_packed as u16, bit1_packed as u16)
        }
    }

    unsafe {
        let newline_vec = vdupq_n_u8(b'\n');
        let case_mask = vdupq_n_u8(0xDF);
        let a_vec = vdupq_n_u8(b'A');
        let c_vec = vdupq_n_u8(b'C');
        let g_vec = vdupq_n_u8(b'G');
        let t_vec = vdupq_n_u8(b'T');
        let lut = vld1q_u8(ENCODE_LUT.as_ptr());

        let mut newlines: u64 = 0;
        let mut is_acgt: u64 = 0;
        let mut two_bits: u128 = 0;

        for chunk_idx in 0..4 {
            let v = vld1q_u8(block.as_ptr().add(chunk_idx * 16));

            let nl_cmp = vceqq_u8(v, newline_vec);
            newlines |= neon_movemask(nl_cmp) << (chunk_idx * 16);

            let upper = vandq_u8(v, case_mask);
            let acgt_cmp = vorrq_u8(
                vorrq_u8(vceqq_u8(upper, a_vec), vceqq_u8(upper, c_vec)),
                vorrq_u8(vceqq_u8(upper, g_vec), vceqq_u8(upper, t_vec)),
            );
            is_acgt |= neon_movemask(acgt_cmp) << (chunk_idx * 16);

            let chunk_two_bits = u128::from(neon_two_bits(v, lut));
            two_bits |= chunk_two_bits << (chunk_idx * 32);
        }

        FastqBitmask { newlines, is_acgt, two_bits }
    }
}

// ============================================================================
// AVX2 (x86_64)
// ============================================================================

/// Reinterpret the low 32 bits of an `i32` movemask result as `u32`.
/// `_mm256_movemask_epi8` returns `i32` but the value is a 32-bit bitmask.
#[cfg(target_arch = "x86_64")]
#[inline]
const fn movemask_u32(mask: i32) -> u32 {
    mask.cast_unsigned()
}

/// Convert a `u8` constant to `i8` for `_mm256_set1_epi8` without wrapping lint.
#[cfg(target_arch = "x86_64")]
#[inline]
const fn byte_as_i8(b: u8) -> i8 {
    i8::from_ne_bytes([b])
}

#[cfg(target_arch = "x86_64")]
/// AVX2 newline-only detection (fast path).
///
/// # Safety
///
/// Requires AVX2 support.
#[target_feature(enable = "avx2")]
unsafe fn lex_block_newlines_avx2(block: &[u8; 64]) -> u64 {
    use std::arch::x86_64::{
        _mm256_cmpeq_epi8, _mm256_loadu_si256, _mm256_movemask_epi8, _mm256_set1_epi8,
    };

    unsafe {
        let newline = _mm256_set1_epi8(i8::from_ne_bytes([b'\n']));
        let v1 = _mm256_loadu_si256(block.as_ptr().cast());
        let v2 = _mm256_loadu_si256(block.as_ptr().add(32).cast());
        let nl1 = movemask_u32(_mm256_movemask_epi8(_mm256_cmpeq_epi8(v1, newline)));
        let nl2 = movemask_u32(_mm256_movemask_epi8(_mm256_cmpeq_epi8(v2, newline)));
        u64::from(nl1) | (u64::from(nl2) << 32)
    }
}

#[cfg(target_arch = "x86_64")]
/// AVX2 implementation: classifies newlines, ACGT bases, and produces 2-bit encoding.
///
/// # Safety
///
/// Requires AVX2 support (checked by caller via `is_x86_feature_detected!`).
#[target_feature(enable = "avx2")]
unsafe fn lex_block_avx2(block: &[u8; 64]) -> FastqBitmask {
    use std::arch::x86_64::{
        __m256i, _mm256_and_si256, _mm256_cmpeq_epi8, _mm256_loadu_si256, _mm256_movemask_epi8,
        _mm256_or_si256, _mm256_set1_epi8, _mm256_shuffle_epi8, _mm256_srli_epi16,
    };

    /// Extract 2-bit encoding from 32 bytes, returning a u64 (32 bases * 2 bits).
    #[inline]
    unsafe fn avx2_two_bits(v: __m256i, lut: __m256i, one_mask: __m256i) -> u64 {
        unsafe {
            let idx = _mm256_and_si256(_mm256_srli_epi16(v, 1), _mm256_set1_epi8(0x0F));
            let encoded = _mm256_shuffle_epi8(lut, idx);

            let bit0 = _mm256_and_si256(encoded, one_mask);
            let bit1 = _mm256_and_si256(_mm256_srli_epi16(encoded, 1), one_mask);

            let bit0_mask = movemask_u32(_mm256_movemask_epi8(_mm256_cmpeq_epi8(bit0, one_mask)));
            let bit1_mask = movemask_u32(_mm256_movemask_epi8(_mm256_cmpeq_epi8(bit1, one_mask)));

            interleave_bits_32(bit0_mask, bit1_mask)
        }
    }

    unsafe {
        let newline = _mm256_set1_epi8(byte_as_i8(b'\n'));
        let case_mask = _mm256_set1_epi8(byte_as_i8(0xDF));
        let a_vec = _mm256_set1_epi8(byte_as_i8(b'A'));
        let c_vec = _mm256_set1_epi8(byte_as_i8(b'C'));
        let g_vec = _mm256_set1_epi8(byte_as_i8(b'G'));
        let t_vec = _mm256_set1_epi8(byte_as_i8(b'T'));

        // Build LUT for both 128-bit lanes (AVX2 shuffle operates per-lane)
        let lut = _mm256_loadu_si256([ENCODE_LUT, ENCODE_LUT].as_ptr().cast());
        let one_mask = _mm256_set1_epi8(0x01);

        let v1 = _mm256_loadu_si256(block.as_ptr().cast());
        let v2 = _mm256_loadu_si256(block.as_ptr().add(32).cast());

        // Newlines
        let nl1 = movemask_u32(_mm256_movemask_epi8(_mm256_cmpeq_epi8(v1, newline)));
        let nl2 = movemask_u32(_mm256_movemask_epi8(_mm256_cmpeq_epi8(v2, newline)));
        let newlines = u64::from(nl1) | (u64::from(nl2) << 32);

        // ACGT classification
        let upper1 = _mm256_and_si256(v1, case_mask);
        let upper2 = _mm256_and_si256(v2, case_mask);
        let acgt1 = _mm256_or_si256(
            _mm256_or_si256(_mm256_cmpeq_epi8(upper1, a_vec), _mm256_cmpeq_epi8(upper1, c_vec)),
            _mm256_or_si256(_mm256_cmpeq_epi8(upper1, g_vec), _mm256_cmpeq_epi8(upper1, t_vec)),
        );
        let acgt2 = _mm256_or_si256(
            _mm256_or_si256(_mm256_cmpeq_epi8(upper2, a_vec), _mm256_cmpeq_epi8(upper2, c_vec)),
            _mm256_or_si256(_mm256_cmpeq_epi8(upper2, g_vec), _mm256_cmpeq_epi8(upper2, t_vec)),
        );
        let acgt_mask1 = movemask_u32(_mm256_movemask_epi8(acgt1));
        let acgt_mask2 = movemask_u32(_mm256_movemask_epi8(acgt2));
        let is_acgt = u64::from(acgt_mask1) | (u64::from(acgt_mask2) << 32);

        // 2-bit encoding
        let tb1 = avx2_two_bits(v1, lut, one_mask);
        let tb2 = avx2_two_bits(v2, lut, one_mask);
        let two_bits = u128::from(tb1) | (u128::from(tb2) << 64);

        FastqBitmask { newlines, is_acgt, two_bits }
    }
}

// ============================================================================
// Bit interleaving helpers
// ============================================================================

/// Interleave bits from `lo` and `hi`: bit `i` of `lo` goes to position `2*i`,
/// bit `i` of `hi` goes to position `2*i+1`. Processes `N` bits from each input.
#[cfg(any(target_arch = "aarch64", target_arch = "x86_64"))]
#[inline]
const fn interleave_bits_n<const N: usize>(lo: u64, hi: u64) -> u64 {
    let mut result: u64 = 0;
    let mut i = 0;
    while i < N {
        result |= ((lo >> i) & 1) << (i * 2);
        result |= ((hi >> i) & 1) << (i * 2 + 1);
        i += 1;
    }
    result
}

/// Interleave two 16-bit values into a 32-bit result.
#[cfg(target_arch = "aarch64")]
#[expect(clippy::cast_possible_truncation, reason = "16 interleaved bits fit in u32")]
#[inline]
const fn interleave_bits(lo: u16, hi: u16) -> u32 {
    interleave_bits_n::<16>(lo as u64, hi as u64) as u32
}

/// Interleave two 32-bit values into a 64-bit result.
#[cfg(target_arch = "x86_64")]
#[inline]
const fn interleave_bits_32(lo: u32, hi: u32) -> u64 {
    interleave_bits_n::<32>(lo as u64, hi as u64)
}

// ============================================================================
// Tests
// ============================================================================

#[cfg(test)]
mod tests {
    use super::*;

    fn make_block(data: &[u8]) -> [u8; 64] {
        let mut block = [0u8; 64];
        let len = data.len().min(64);
        block[..len].copy_from_slice(&data[..len]);
        block
    }

    // --- Newline tests (unchanged) ---

    #[test]
    fn test_scalar_no_newlines() {
        let block = make_block(b"ACGTACGTACGTACGTACGTACGTACGTACGT");
        let bm = lex_block_scalar(&block);
        assert_eq!(bm.newlines, 0);
    }

    #[test]
    fn test_scalar_all_newlines() {
        let block = [b'\n'; 64];
        let bm = lex_block_scalar(&block);
        assert_eq!(bm.newlines, u64::MAX);
    }

    #[test]
    fn test_scalar_known_positions() {
        let mut block = [b'A'; 64];
        block[0] = b'\n';
        block[3] = b'\n';
        block[63] = b'\n';
        let bm = lex_block_scalar(&block);
        assert_eq!(bm.newlines, (1 << 0) | (1 << 3) | (1 << 63));
    }

    #[test]
    fn test_scalar_fastq_record() {
        let block = make_block(b"@r1\nACGT\n+\nIIII\n");
        let bm = lex_block_scalar(&block);
        assert_eq!(bm.newlines, (1 << 3) | (1 << 8) | (1 << 10) | (1 << 15));
    }

    // --- ACGT classification tests ---

    #[test]
    fn test_scalar_acgt_classification() {
        let block = make_block(b"ACGTacgt");
        let bm = lex_block_scalar(&block);
        // All 8 bases should be classified as ACGT
        assert_eq!(bm.is_acgt & 0xFF, 0xFF);
        // Padding zeros are not ACGT
        assert_eq!(bm.is_acgt >> 8, 0);
    }

    #[test]
    fn test_scalar_non_acgt() {
        let block = make_block(b"@+N\nIIII");
        let bm = lex_block_scalar(&block);
        // None of these are ACGT
        assert_eq!(bm.is_acgt & 0xFF, 0);
    }

    #[test]
    fn test_scalar_mixed_acgt() {
        let block = make_block(b"@ACNGT\n");
        let bm = lex_block_scalar(&block);
        // Positions: @(0)=no, A(1)=yes, C(2)=yes, N(3)=no, G(4)=yes, T(5)=yes, \n(6)=no
        assert_eq!(bm.is_acgt & 0x7F, 0b011_0110);
    }

    // --- 2-bit encoding tests ---

    #[test]
    fn test_scalar_two_bit_encoding() {
        // A=0, C=1, G=2, T=3 (fgumi BitEnc mapping)
        let block = make_block(b"ACGT");
        let bm = lex_block_scalar(&block);

        // Position 0: A → 0b00
        assert_eq!(bm.two_bits & 0x3, 0);
        // Position 1: C → 0b01
        assert_eq!((bm.two_bits >> 2) & 0x3, 1);
        // Position 2: G → 0b10
        assert_eq!((bm.two_bits >> 4) & 0x3, 2);
        // Position 3: T → 0b11
        assert_eq!((bm.two_bits >> 6) & 0x3, 3);
    }

    #[test]
    fn test_scalar_two_bit_case_insensitive() {
        let upper = make_block(b"ACGT");
        let lower = make_block(b"acgt");
        let bm_upper = lex_block_scalar(&upper);
        let bm_lower = lex_block_scalar(&lower);
        // First 8 bits of two_bits should match (4 bases * 2 bits each)
        assert_eq!(bm_upper.two_bits & 0xFF, bm_lower.two_bits & 0xFF);
    }

    // --- SIMD vs scalar consistency (newlines only) ---

    #[test]
    fn test_simd_matches_scalar_newlines() {
        let mut block = [b'A'; 64];
        for (i, byte) in block.iter_mut().enumerate() {
            if i % 7 == 0 {
                *byte = b'\n';
            }
        }
        let simd_result = lex_block(&block);
        let scalar_result = lex_block_scalar(&block).newlines;
        assert_eq!(simd_result, scalar_result, "Newline bitmask mismatch");
    }

    #[test]
    fn test_simd_newlines_every_position() {
        for pos in 0..64 {
            let mut block = [b'X'; 64];
            block[pos] = b'\n';
            let simd_result = lex_block(&block);
            let scalar_result = lex_block_scalar(&block).newlines;
            assert_eq!(simd_result, scalar_result, "Mismatch at position {pos}");
        }
    }

    // --- SIMD vs scalar consistency (full classification) ---

    #[test]
    fn test_simd_full_matches_scalar_all_fields() {
        let block = make_block(b"@read1\nACGTACGTacgtNNNN\n+\nIIIIIIIIIIIIIIII\n");
        let simd_result = lex_block_full(&block);
        let scalar_result = lex_block_scalar(&block);
        assert_eq!(simd_result.newlines, scalar_result.newlines, "newlines mismatch");
        assert_eq!(simd_result.is_acgt, scalar_result.is_acgt, "is_acgt mismatch");
        assert_eq!(simd_result.two_bits, scalar_result.two_bits, "two_bits mismatch");
    }

    #[test]
    fn test_simd_full_matches_scalar_all_acgt() {
        let block = {
            let mut b = [0u8; 64];
            for (i, byte) in b.iter_mut().enumerate() {
                *byte = b"ACGTacgt"[i % 8];
            }
            b
        };
        let simd_result = lex_block_full(&block);
        let scalar_result = lex_block_scalar(&block);
        assert_eq!(simd_result.is_acgt, scalar_result.is_acgt, "is_acgt mismatch");
        assert_eq!(simd_result.two_bits, scalar_result.two_bits, "two_bits mismatch");
    }

    #[test]
    fn test_simd_full_matches_scalar_every_position() {
        for base in [b'A', b'C', b'G', b'T', b'a', b'c', b'g', b't'] {
            for pos in 0..64 {
                let mut block = [b'N'; 64];
                block[pos] = base;
                let simd = lex_block_full(&block);
                let scalar = lex_block_scalar(&block);
                assert_eq!(
                    simd.is_acgt, scalar.is_acgt,
                    "is_acgt mismatch for {base} at pos {pos}"
                );
                let mask = 0x3u128 << (pos * 2);
                assert_eq!(
                    simd.two_bits & mask,
                    scalar.two_bits & mask,
                    "two_bits mismatch for {} at pos {pos}",
                    base as char,
                );
            }
        }
    }
}
