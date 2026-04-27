//! Numeric-equivalence regression test for `unclipped_five_prime_position`.
//!
//! Both `fgumi duplex-metrics` and `fgumi simplex-metrics` depend on
//! [`fgumi_lib::commands::shared_metrics::unclipped_five_prime_position`] to
//! compute the reverse-strand `unclippedEnd = alignmentEnd + trailing clips`
//! coordinate used in [`ReadInfoKey`]. During the raw-byte migration the
//! implementation was swapped from a noodles-`RecordBuf`-based algorithm to
//! the zero-allocation `fgumi_raw_bam::unclipped_5prime_from_raw_bam` helper.
//! If the two diverge on any CIGAR shape, `ReadInfoKey` hashes change, which
//! silently shifts the downsampling cohort and corrupts every metrics TSV
//! emitted by either command.
//!
//! Each case here is computed independently via the htsjdk / fgbio reference
//! formula (walking the CIGAR ops directly) and compared against the
//! `shared_metrics::unclipped_five_prime_position` result on a `RawRecord`
//! carrying the same CIGAR. Cases cover both forward and reverse strands with
//! combinations of soft (`S`) and hard (`H`) clips — the surface that matters
//! for reverse-strand metrics correctness.
//!
//! CIGAR op codes in BAM: `0=M 1=I 2=D 3=N 4=S 5=H 6=P 7== 8=X`.

#![cfg(feature = "simplex")]

use fgumi_lib::commands::shared_metrics::unclipped_five_prime_position;
use fgumi_raw_bam::{RawRecord, SamBuilder, flags as raw_flags};
use rstest::rstest;

/// Reference implementation: fgbio / htsjdk's `SAMRecord.getUnclippedStart/End`.
///
/// - Forward: `alignment_start_1based - (leading S) - (leading H)`
/// - Reverse: `alignment_end_1based   + (trailing S) + (trailing H)`
///
/// where `alignment_end_1based = alignment_start_1based + ref_consumed - 1`,
/// and "leading" / "trailing" clip counts include both soft (4) and hard (5)
/// clip ops from the start / end of the CIGAR.
fn reference_unclipped_5prime(
    alignment_start_1based: i32,
    is_reverse: bool,
    cigar_ops: &[(u32, u32)], // (op_code, op_length)
) -> i32 {
    let consumes_ref = |op: u32| matches!(op, 0 | 2 | 3 | 7 | 8);
    let is_clip = |op: u32| matches!(op, 4 | 5);

    if is_reverse {
        let ref_len: i32 = cigar_ops
            .iter()
            .filter(|(op, _)| consumes_ref(*op))
            .map(|(_, len)| (*len).cast_signed())
            .sum();
        let alignment_end = alignment_start_1based + ref_len - 1;
        let trailing_clips: i32 = cigar_ops
            .iter()
            .rev()
            .take_while(|(op, _)| is_clip(*op))
            .map(|(_, len)| (*len).cast_signed())
            .sum();
        alignment_end + trailing_clips
    } else {
        let leading_clips: i32 = cigar_ops
            .iter()
            .take_while(|(op, _)| is_clip(*op))
            .map(|(_, len)| (*len).cast_signed())
            .sum();
        alignment_start_1based - leading_clips
    }
}

/// Build a minimal `RawRecord` with the given strand / position / CIGAR.
fn build_raw(alignment_start_0based: i32, is_reverse: bool, cigar_ops: &[(u32, u32)]) -> RawRecord {
    let encoded: Vec<u32> = cigar_ops.iter().map(|(op, len)| (*len << 4) | (*op & 0xF)).collect();
    let mut flags: u16 = 0;
    if is_reverse {
        flags |= raw_flags::REVERSE;
    }
    // Compute consumed-query length so sequence/qualities length matches the CIGAR.
    let consumes_query = |op: u32| matches!(op, 0 | 1 | 4 | 7 | 8);
    let query_len: usize =
        cigar_ops.iter().filter(|(op, _)| consumes_query(*op)).map(|(_, len)| *len as usize).sum();
    let seq = vec![b'A'; query_len];
    let quals = vec![30u8; query_len];
    let mut b = SamBuilder::new();
    b.read_name(b"read1")
        .flags(flags)
        .ref_id(0)
        .pos(alignment_start_0based)
        .mapq(60)
        .cigar_ops(&encoded)
        .sequence(&seq)
        .qualities(&quals);
    b.build()
}

#[rstest]
// --- Reverse strand (primary regression target for the raw-byte migration) ---
// pos_0based=100 → pos_1based=101; alignment_end = 101 + ref_len - 1;
// unclipped_5' = alignment_end + trailing_soft_plus_hard_clips.
#[case::reverse_all_match(100, true, &[(0, 50)], 150)]
#[case::reverse_trailing_soft(100, true, &[(0, 96), (4, 4)], 200)]
#[case::reverse_leading_and_trailing_soft(100, true, &[(4, 4), (0, 92), (4, 4)], 196)]
#[case::reverse_trailing_hard(100, true, &[(0, 50), (5, 3)], 153)]
#[case::reverse_trailing_soft_then_hard(100, true, &[(0, 50), (4, 4), (5, 2)], 156)]
#[case::reverse_indel_middle(100, true, &[(0, 40), (1, 5), (0, 10), (4, 3)], 153)]
#[case::reverse_deletion_middle(100, true, &[(0, 40), (2, 5), (0, 10), (4, 3)], 158)]
// --- Forward strand (sanity — guards against future breakage). ---
#[case::forward_all_match(100, false, &[(0, 50)], 101)]
#[case::forward_leading_soft(100, false, &[(4, 4), (0, 96)], 97)]
#[case::forward_leading_hard_and_soft(100, false, &[(5, 3), (4, 4), (0, 93)], 94)]
fn unclipped_5prime_matches_fgbio_reference(
    #[case] alignment_start_0based: i32,
    #[case] is_reverse: bool,
    #[case] cigar_ops: &[(u32, u32)],
    #[case] expected: i32,
) {
    let rec = build_raw(alignment_start_0based, is_reverse, cigar_ops);
    let from_helper = unclipped_five_prime_position(&rec).expect("mapped read with a CIGAR");

    // Double-check the expected value actually matches the inline reference
    // formula, so a typo in the #[case] literal can't silently pass.
    let alignment_start_1based = alignment_start_0based + 1;
    let from_reference = reference_unclipped_5prime(alignment_start_1based, is_reverse, cigar_ops);
    assert_eq!(
        from_reference, expected,
        "reference formula disagrees with hand-computed expected — fix the #[case] literal"
    );

    assert_eq!(
        from_helper,
        expected,
        "unclipped_five_prime_position diverges from fgbio reference for \
         CIGAR {:?} on the {} strand",
        cigar_ops,
        if is_reverse { "reverse" } else { "forward" },
    );
}

/// Unmapped reads must return `None` — the downsampling cohort logic skips them.
#[test]
fn unclipped_5prime_returns_none_for_unmapped() {
    let mut b = SamBuilder::new();
    b.read_name(b"read1")
        .flags(raw_flags::UNMAPPED)
        .ref_id(-1)
        .pos(-1)
        .mapq(0)
        .cigar_ops(&[])
        .sequence(b"ACGT")
        .qualities(&[30, 30, 30, 30]);
    let rec = b.build();
    assert_eq!(unclipped_five_prime_position(&rec), None);
}
