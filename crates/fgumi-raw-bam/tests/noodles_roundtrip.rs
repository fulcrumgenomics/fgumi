#![cfg(feature = "noodles")]

use fgumi_raw_bam::testutil::{encode_op, make_bam_bytes};
use fgumi_raw_bam::{RawRecord, raw_records_to_record_bufs};
use noodles::sam::alignment::record::Cigar;

/// Decode a single [`RawRecord`] through noodles and return the parsed record.
///
/// Panics if noodles cannot parse the bytes — that failure is the test signal.
fn decode_with_noodles(rec: &RawRecord) -> noodles::sam::alignment::RecordBuf {
    let buffers = vec![rec.as_ref().to_vec()];
    let mut bufs = raw_records_to_record_bufs(&buffers).expect("noodles decode failed");
    bufs.pop().expect("expected at least one decoded record")
}

/// Build a minimal valid [`RawRecord`] for use as proptest input.
///
/// Constructs a record with the given name, a single Match CIGAR op covering
/// `l_seq` bases, and no aux data.  The sequence bytes are left as all-zeros
/// (the BAM encoding for all-N) which is valid.
fn build_record(name: &[u8], l_seq: usize) -> RawRecord {
    let cigar = if l_seq > 0 { vec![encode_op(0, l_seq)] } else { vec![] };
    let bytes = make_bam_bytes(0, 0, 0, name, &cigar, l_seq, -1, -1, b"");
    RawRecord::from(bytes)
}

proptest::proptest! {
    /// Verify that [`RawRecord::set_read_name`] produces bytes that noodles decodes
    /// to the same name, across many random name strings.
    #[test]
    fn set_read_name_roundtrips(name in "[A-Za-z0-9_]{1,200}") {
        let mut rec = build_record(b"placeholder", 4);
        rec.set_read_name(name.as_bytes());

        // Internal view must reflect the new name.
        proptest::prop_assert_eq!(rec.read_name(), name.as_bytes());

        // Noodles must decode the same name.
        let buf = decode_with_noodles(&rec);
        let decoded = buf.name().map(std::convert::AsRef::as_ref);
        proptest::prop_assert_eq!(decoded, Some(name.as_bytes()));
    }

    /// Verify that [`RawRecord::set_cigar_ops`] produces bytes that noodles decodes
    /// to the same number of CIGAR operations, across many random op sequences.
    #[test]
    fn set_cigar_ops_roundtrips(
        ops in proptest::collection::vec(1u32..=100u32, 1..=10),
    ) {
        // Build all-Match ops so the sequence length equals the sum of op lengths.
        let total: usize = ops.iter().map(|&l| l as usize).sum();
        // CIGAR word encodes length in upper bits, op type (M=0) in lower 4 bits.
        let cigar_ops: Vec<u32> = ops.iter().map(|&len| len << 4).collect();

        let mut rec = build_record(b"r", total);
        rec.set_cigar_ops(&cigar_ops);

        // Internal view must reflect the new CIGAR.
        proptest::prop_assert_eq!(rec.cigar_ops_vec(), cigar_ops.clone());

        // Noodles must decode the same number of ops.
        let buf = decode_with_noodles(&rec);
        let decoded_n_ops = buf.cigar().iter().count();
        proptest::prop_assert_eq!(decoded_n_ops, ops.len());
    }

    /// Verify that [`RawRecord::set_sequence_and_qualities`] produces bytes that
    /// noodles can parse without error, and that the internal view matches across
    /// many random sequence lengths.
    #[test]
    fn set_sequence_and_qualities_roundtrips(
        n in 1usize..=200,
    ) {
        // Build deterministic bases + quals from the length.
        let bases: Vec<u8> = (0..n).map(|i| b"ACGTN"[i % 5]).collect();
        // i % 40 is always 0..=39, so the truncation is safe.
        #[allow(clippy::cast_possible_truncation)]
        let quals: Vec<u8> = (0..n).map(|i| (i % 40) as u8).collect();

        let mut rec = build_record(b"r", 1);
        // Update the CIGAR to match the new sequence length before writing seq+qual.
        // n <= 200 and u32::MAX >> 4 is well above 200, so the cast is safe.
        #[allow(clippy::cast_possible_truncation)]
        rec.set_cigar_ops(&[(n as u32) << 4]);
        rec.set_sequence_and_qualities(&bases, &quals);

        // Internal view must reflect the new sequence and qualities.
        proptest::prop_assert_eq!(rec.sequence_vec(), bases.clone());
        proptest::prop_assert_eq!(rec.quality_scores().to_vec(), quals.clone());

        // Noodles must be able to parse the record without error.
        decode_with_noodles(&rec);
    }
}
