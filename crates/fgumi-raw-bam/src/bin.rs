//! BAM bin field (`reg2bin`) computation.
//!
//! The BAM binary record carries a 16-bit `bin` field (bytes 10–11) that indexes
//! the alignment into the UCSC binning scheme (SAM spec §5.3). The bin is a pure
//! function of the alignment's reference span — its 0-based start (`POS`) and the
//! reference length consumed by its CIGAR — so any operation that moves `POS` or
//! rewrites the CIGAR (clipping, mate relocation, unmapping) invalidates it and
//! must recompute it.
//!
//! The noodles BAM *encoder* recomputes the bin from `POS`+CIGAR whenever a
//! `RecordBuf` is serialized, and htsjdk does the same on write — but fgumi's
//! raw-byte pipelines mutate the encoded record in place and emit the bytes
//! verbatim, so there is no encode step to refresh the field. These helpers are
//! the raw-byte equivalent: call [`set_bin_from_raw_bam`] (or the
//! `RawRecord::recompute_bin` convenience wrapper) after any in-place `POS`/CIGAR
//! change so the emitted record carries a correct bin.

use crate::cigar::reference_length_from_raw_bam;
use crate::fields::pos;

/// The BAM "unmapped bin" — `reg2bin(-1, 0)` per SAM spec §4.2.1. Records with no
/// alignment position (`POS < 0`) use this value.
pub const UNMAPPED_BIN: u16 = 4680;

/// Compute the BAM bin for a reference region using the SAM spec §5.3 `reg2bin`
/// reference implementation.
///
/// `beg` is the 0-based leftmost aligned position; `end` is the 0-based position
/// one past the rightmost aligned base (i.e. a half-open `[beg, end)` interval,
/// which equals `POS + reference_length`). Callers with an unmapped region should
/// use [`UNMAPPED_BIN`] directly rather than passing negative coordinates here.
///
/// The layer constants `((1 << k) - 1) / 7` are written verbatim from the SAM
/// specification's reference C code; clippy's `eq_op` lint flags the smallest
/// layer (`((1 << 3) - 1) / 7`) as a tautology, but the form is kept to stay
/// textually identical to the spec.
#[allow(clippy::eq_op, clippy::cast_possible_truncation, clippy::cast_sign_loss)]
#[must_use]
pub fn reg2bin(beg: i32, end: i32) -> u16 {
    // `end` arrives half-open; the spec reference code works with the 0-based
    // position of the last aligned base, so decrement once up front.
    let end = end - 1;
    let bin = if beg >> 14 == end >> 14 {
        ((1 << 15) - 1) / 7 + (beg >> 14)
    } else if beg >> 17 == end >> 17 {
        ((1 << 12) - 1) / 7 + (beg >> 17)
    } else if beg >> 20 == end >> 20 {
        ((1 << 9) - 1) / 7 + (beg >> 20)
    } else if beg >> 23 == end >> 23 {
        ((1 << 6) - 1) / 7 + (beg >> 23)
    } else if beg >> 26 == end >> 26 {
        ((1 << 3) - 1) / 7 + (beg >> 26)
    } else {
        0
    };
    // Truncate overflowing bin IDs, matching noodles/htsjdk behaviour.
    bin as u16
}

/// Compute the BAM bin for a raw BAM record from its `POS` and CIGAR reference
/// span. Returns [`UNMAPPED_BIN`] when the record has no alignment position
/// (`POS < 0`), matching htsjdk `SAMRecord.computeIndexingBin` (which keys off
/// the alignment start, not the unmapped flag — a flag-unmapped but *placed*
/// read still gets its position's bin).
#[must_use]
pub fn bin_from_raw_bam(bam: &[u8]) -> u16 {
    let start = pos(bam);
    if start < 0 {
        return UNMAPPED_BIN;
    }
    // A placed read with a CIGAR that consumes no reference (e.g. all soft-clip)
    // still occupies a single reference position; clamp to 1 so its bin reflects
    // its start rather than degenerating (matches htsjdk `alignmentEnd = start`).
    let ref_len = reference_length_from_raw_bam(bam).max(1);
    reg2bin(start, start + ref_len)
}

/// Recompute and write the BAM bin field (bytes 10–11, little-endian) of a raw
/// BAM record in place, from its current `POS` and CIGAR. Call after any in-place
/// mutation of `POS` or the CIGAR.
pub fn set_bin_from_raw_bam(bam: &mut [u8]) {
    let value = bin_from_raw_bam(bam);
    bam[10..12].copy_from_slice(&value.to_le_bytes());
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::SamBuilder;
    use crate::testutil::encode_op;
    use rstest::rstest;

    // CIGAR op codes (BAM encoding): M=0, S=4.
    const M: u32 = 0;
    const S: u32 = 4;

    // Authoritative `reg2bin` vectors ported verbatim from htsjdk
    // `GenomicIndexUtilTest.testRegionToBinDataProvider` (regionToBin uses the
    // identical `[beg, end)` 0-based convention). Exercises every binning level:
    // level 5 (16 kb), level 4 (128 kb), level 3 (1 Mb), level 2 (8 Mb),
    // level 1 (64 Mb), and level 0 (whole-window bin 0). Note htsjdk's `1<<26+1`
    // is `1<<27` under Java/Rust shift-vs-`+` precedence.
    #[rstest]
    #[case::empty_at_zero(0, 0, 0)]
    #[case::base_at_pos1(1, 1, 4681)]
    #[case::first_16kb_window(0, 1 << 14, 4681)]
    #[case::spans_into_second_16kb(0, (1 << 14) + 1, 585)]
    #[case::empty_at_16kb_boundary(1 << 14, 1 << 14, 585)]
    #[case::base_past_16kb((1 << 14) + 1, (1 << 14) + 1, 4682)]
    #[case::range_within_first_128kb(1 << 14, 1 << 17, 585)]
    #[case::spans_into_second_128kb(1 << 14, (1 << 17) + 1, 73)]
    #[case::empty_at_128kb_boundary(1 << 17, 1 << 17, 73)]
    #[case::base_past_128kb((1 << 17) + 1, (1 << 17) + 1, 4689)]
    #[case::range_within_first_1mb(1 << 17, 1 << 20, 73)]
    #[case::spans_into_second_1mb(1 << 17, (1 << 20) + 1, 9)]
    #[case::empty_at_1mb_boundary(1 << 20, 1 << 20, 9)]
    #[case::base_past_1mb((1 << 20) + 1, (1 << 20) + 1, 4745)]
    #[case::range_within_first_8mb(1 << 20, 1 << 23, 9)]
    #[case::spans_into_second_8mb(1 << 20, (1 << 23) + 1, 1)]
    #[case::empty_at_8mb_boundary(1 << 23, 1 << 23, 1)]
    #[case::base_past_8mb((1 << 23) + 1, (1 << 23) + 1, 5193)]
    #[case::range_within_first_64mb(1 << 23, 1 << 26, 1)]
    #[case::spans_into_second_64mb(1 << 23, (1 << 26) + 1, 0)]
    #[case::empty_at_64mb_boundary(1 << 26, 1 << 26, 0)]
    #[case::base_past_64mb((1 << 26) + 1, (1 << 26) + 1, 8777)]
    #[case::range_within_first_512mb(1 << 26, 1 << 27, 2)]
    fn reg2bin_matches_htsjdk(#[case] beg: i32, #[case] end: i32, #[case] expected: u16) {
        assert_eq!(reg2bin(beg, end), expected);
    }

    /// Build a minimal mapped record at `pos` with a single-op CIGAR, then read
    /// its stored bin bytes back after recompute.
    fn build_and_recompute(ref_id: i32, pos: i32, op: u32, op_len: usize) -> u16 {
        let mut b = SamBuilder::new();
        let mut rec = b
            .ref_id(ref_id)
            .pos(pos)
            .read_name(b"r")
            .cigar_ops(&[encode_op(op, op_len)])
            .sequence(&vec![b'A'; op_len])
            .qualities(&vec![30u8; op_len])
            .build();
        set_bin_from_raw_bam(rec.as_mut_vec());
        rec.bin()
    }

    #[test]
    fn set_bin_from_raw_bam_writes_position_bin_for_mapped_read() {
        // 84M at pos 16416 -> [16416, 16500) -> level-5 bin 4682.
        assert_eq!(build_and_recompute(0, 16416, M, 84), 4682);
    }

    #[test]
    fn set_bin_from_raw_bam_writes_level4_bin_when_straddling_boundary() {
        // 200M at pos 16300 -> [16300, 16500) straddles the 16 kb boundary -> 585.
        assert_eq!(build_and_recompute(0, 16300, M, 200), 585);
    }

    #[test]
    fn set_bin_from_raw_bam_uses_unmapped_bin_when_pos_negative() {
        // pos = -1 -> UNMAPPED_BIN regardless of CIGAR bytes present.
        assert_eq!(build_and_recompute(-1, -1, S, 100), UNMAPPED_BIN);
    }

    #[test]
    fn set_bin_from_raw_bam_clamps_zero_reference_span_placed_read() {
        // A placed read (pos >= 0) whose CIGAR consumes no reference (all soft-clip)
        // still gets its start position's bin via the `ref_len.max(1)` clamp:
        // 100S at pos 16416 -> reg2bin(16416, 16417) -> 4682, not the unmapped bin.
        assert_eq!(build_and_recompute(0, 16416, S, 100), 4682);
    }

    #[test]
    fn bin_from_raw_bam_ignores_flag_uses_alignment_start() {
        // A placed read (pos >= 0) gets its position bin from bin_from_raw_bam.
        let mut b = SamBuilder::new();
        let rec = b
            .ref_id(0)
            .pos(16416)
            .read_name(b"r")
            .cigar_ops(&[encode_op(M, 84)])
            .sequence(&[b'A'; 84])
            .qualities(&[30u8; 84])
            .build();
        assert_eq!(bin_from_raw_bam(rec.as_ref()), 4682);
    }
}
