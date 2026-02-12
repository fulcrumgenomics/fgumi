/// Simplify CIGAR operations from raw BAM u32 ops.
///
/// Same logic as `cigar_utils::simplify_cigar` but operates on raw BAM CIGAR
/// u32 words instead of noodles `Cigar`. Converts S, =, X, H operations to M
/// and coalesces adjacent operations of the same type.
///
/// Uses `noodles::sam::alignment::record::cigar::op::Kind` for the output
/// representation to stay compatible with the existing `SimplifiedCigar` type.
#[must_use]
pub fn simplify_cigar_from_raw(
    cigar_ops: &[u32],
) -> Vec<(noodles::sam::alignment::record::cigar::op::Kind, usize)> {
    use noodles::sam::alignment::record::cigar::op::Kind;

    let mut simplified = Vec::new();

    for &raw_op in cigar_ops {
        let op_len = (raw_op >> 4) as usize;
        let op_type = raw_op & 0xF;

        // Map BAM CIGAR op to Kind
        let kind = match op_type {
            0 => Kind::Match,
            1 => Kind::Insertion,
            2 => Kind::Deletion,
            3 => Kind::Skip,
            4 => Kind::SoftClip,
            5 => Kind::HardClip,
            6 => Kind::Pad,
            7 => Kind::SequenceMatch,
            8 => Kind::SequenceMismatch,
            _ => continue,
        };

        // Simplify: convert S, =, X, H to M
        let new_kind = match kind {
            Kind::SoftClip | Kind::SequenceMatch | Kind::SequenceMismatch | Kind::HardClip => {
                Kind::Match
            }
            _ => kind,
        };

        // Coalesce adjacent operations of the same type
        if let Some((last_kind, last_len)) = simplified.last_mut() {
            if *last_kind == new_kind {
                *last_len += op_len;
                continue;
            }
        }

        simplified.push((new_kind, op_len));
    }

    simplified
}

/// Decode raw BAM byte records to noodles `RecordBuf`.
///
/// This is a temporary bridge for consensus callers (duplex, codec) that still
/// use `RecordBuf` internally. It constructs a minimal BAM stream and uses
/// noodles reader to parse records.
///
/// # Errors
///
/// Returns an error if the BAM header cannot be read or if any record fails to
/// parse (e.g., due to malformed binary data).
///
/// # Panics
///
/// Panics if any individual record byte length exceeds `u32::MAX`.
pub fn raw_records_to_record_bufs(
    records: &[Vec<u8>],
) -> anyhow::Result<Vec<noodles::sam::alignment::RecordBuf>> {
    use std::io::Cursor;

    let header = noodles::sam::Header::default();
    let mut bam_data: Vec<u8> = Vec::new();

    // BAM magic
    bam_data.extend_from_slice(b"BAM\x01");
    // Header text length = 0
    bam_data.extend_from_slice(&0u32.to_le_bytes());
    // n_ref = 0
    bam_data.extend_from_slice(&0u32.to_le_bytes());

    for raw in records {
        let block_size = u32::try_from(raw.len()).expect("record too large for BAM block_size");
        bam_data.extend_from_slice(&block_size.to_le_bytes());
        bam_data.extend_from_slice(raw);
    }

    let cursor = Cursor::new(bam_data);
    let mut reader = noodles::bam::io::Reader::from(cursor);
    let _ = reader.read_header()?;

    let mut result = Vec::with_capacity(records.len());
    for record_result in reader.record_bufs(&header) {
        result.push(record_result?);
    }

    Ok(result)
}

#[cfg(test)]
#[allow(clippy::identity_op)]
mod tests {
    use super::*;
    use crate::builder::*;
    use crate::testutil::*;

    // ========================================================================
    // simplify_cigar_from_raw tests
    // ========================================================================

    #[test]
    fn test_simplify_cigar_from_raw_basic() {
        use noodles::sam::alignment::record::cigar::op::Kind;
        // 5S10M3I5M4S => all S/= become M, so: 5M+10M+3I+5M+4M => 15M+3I+9M
        let cigar = &[
            encode_op(4, 5),  // 5S
            encode_op(0, 10), // 10M
            encode_op(1, 3),  // 3I
            encode_op(0, 5),  // 5M
            encode_op(4, 4),  // 4S
        ];
        let result = simplify_cigar_from_raw(cigar);
        assert_eq!(
            result,
            vec![
                (Kind::Match, 15), // 5S+10M coalesced
                (Kind::Insertion, 3),
                (Kind::Match, 9), // 5M+4S coalesced
            ]
        );
    }

    #[test]
    fn test_simplify_cigar_from_raw_eq_and_x() {
        use noodles::sam::alignment::record::cigar::op::Kind;
        // 5=3X2D4= => all =, X become M: 5M+3M+2D+4M => 8M+2D+4M
        let cigar = &[
            encode_op(7, 5), // 5=
            encode_op(8, 3), // 3X
            encode_op(2, 2), // 2D
            encode_op(7, 4), // 4=
        ];
        let result = simplify_cigar_from_raw(cigar);
        assert_eq!(result, vec![(Kind::Match, 8), (Kind::Deletion, 2), (Kind::Match, 4),]);
    }

    #[test]
    fn test_simplify_cigar_from_raw_empty() {
        let result = simplify_cigar_from_raw(&[]);
        assert!(result.is_empty());
    }

    #[test]
    fn test_simplify_cigar_from_raw_hard_clips() {
        use noodles::sam::alignment::record::cigar::op::Kind;
        // 5H10M5H -> H becomes M, coalesced: 20M
        let cigar = &[encode_op(5, 5), encode_op(0, 10), encode_op(5, 5)];
        let result = simplify_cigar_from_raw(cigar);
        assert_eq!(result, vec![(Kind::Match, 20)]);
    }

    #[test]
    fn test_simplify_cigar_from_raw_mixed_operations() {
        use noodles::sam::alignment::record::cigar::op::Kind;
        // 2H3S5M2I3M1D4M4S1H -> H,S,M coalesce to M; I,D stay
        // 2M+3M+5M = 10M, 2I, 3M, 1D, 4M+4M+1M = 9M
        let cigar = &[
            encode_op(5, 2), // 2H
            encode_op(4, 3), // 3S
            encode_op(0, 5), // 5M
            encode_op(1, 2), // 2I
            encode_op(0, 3), // 3M
            encode_op(2, 1), // 1D
            encode_op(0, 4), // 4M
            encode_op(4, 4), // 4S
            encode_op(5, 1), // 1H
        ];
        let result = simplify_cigar_from_raw(cigar);
        assert_eq!(
            result,
            vec![
                (Kind::Match, 10),
                (Kind::Insertion, 2),
                (Kind::Match, 3),
                (Kind::Deletion, 1),
                (Kind::Match, 9),
            ]
        );
    }

    #[test]
    fn test_simplify_cigar_from_raw_skip_and_pad() {
        use noodles::sam::alignment::record::cigar::op::Kind;
        // 5M3N5M -> N (skip) preserved, not converted to M
        let cigar = &[encode_op(0, 5), encode_op(3, 3), encode_op(0, 5)];
        let result = simplify_cigar_from_raw(cigar);
        assert_eq!(result, vec![(Kind::Match, 5), (Kind::Skip, 3), (Kind::Match, 5)]);
    }

    #[test]
    fn test_simplify_cigar_from_raw_all_soft_clips() {
        use noodles::sam::alignment::record::cigar::op::Kind;
        // 10S -> becomes 10M
        let cigar = &[encode_op(4, 10)];
        let result = simplify_cigar_from_raw(cigar);
        assert_eq!(result, vec![(Kind::Match, 10)]);
    }

    // ========================================================================
    // raw_records_to_record_bufs tests
    // ========================================================================

    #[test]
    fn test_raw_records_to_record_bufs_single() {
        let mut builder = UnmappedBamRecordBuilder::new();
        builder.build_record(b"read1", 0, b"ACGT", &[30, 25, 20, 15]);
        builder.append_string_tag(b"MI", b"1");
        let raw = builder.as_bytes().to_vec();

        let bufs = raw_records_to_record_bufs(&[raw]).unwrap();
        assert_eq!(bufs.len(), 1);
        assert_eq!(bufs[0].name().map(std::convert::AsRef::as_ref), Some(b"read1".as_ref()));
    }

    #[test]
    fn test_raw_records_to_record_bufs_multiple() {
        let mut records = Vec::new();
        for i in 0..3 {
            let mut builder = UnmappedBamRecordBuilder::new();
            let name = format!("read{i}");
            builder.build_record(name.as_bytes(), 0, b"AC", &[30, 25]);
            records.push(builder.as_bytes().to_vec());
        }

        let bufs = raw_records_to_record_bufs(&records).unwrap();
        assert_eq!(bufs.len(), 3);
        for (i, buf) in bufs.iter().enumerate() {
            let expected = format!("read{i}");
            assert_eq!(buf.name().map(<_ as AsRef<[u8]>>::as_ref), Some(expected.as_bytes()));
        }
    }

    #[test]
    fn test_raw_records_to_record_bufs_empty() {
        let bufs = raw_records_to_record_bufs(&[]).unwrap();
        assert!(bufs.is_empty());
    }
}
