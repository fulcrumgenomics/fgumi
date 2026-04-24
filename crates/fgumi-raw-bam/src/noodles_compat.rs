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

/// Length of the 4-byte little-endian `block_size` prefix in the BAM framing protocol.
///
/// When `noodles::bam::io::Writer<Vec<u8>>` writes a record it prepends a `u32` `block_size`
/// field before the record payload.  This constant names that prefix so it can be stripped.
const BLOCK_SIZE_PREFIX_LEN: usize = 4;

/// Encode a single `RecordBuf` to raw BAM bytes and wrap in a [`RawRecord`].
///
/// Uses noodles' BAM writer machinery writing to an in-memory buffer to produce
/// the record's binary representation. The 4-byte `block_size` prefix written
/// by the BAM framing protocol is stripped so the returned `RawRecord` contains
/// only the record payload bytes (matching the layout of records produced by
/// [`raw_records_to_record_bufs`] and consumed by [`RawRecord`] field accessors).
///
/// # Errors
///
/// Returns an error if encoding fails (e.g., the record contains an invalid
/// reference-sequence ID or alignment position for the given `header`).
pub fn encode_record_buf_to_raw(
    buf: &noodles::sam::alignment::RecordBuf,
    header: &noodles::sam::Header,
) -> anyhow::Result<crate::RawRecord> {
    use noodles::sam::alignment::io::Write as _;

    // Writer<Vec<u8>> is uncompressed — writes 4-byte block_size + record bytes.
    let mut writer = noodles::bam::io::Writer::from(Vec::new());
    writer
        .write_alignment_record(header, buf)
        .map_err(|e| anyhow::anyhow!("encode_record_buf_to_raw: {e}"))?;
    let raw = writer.into_inner();

    // Skip the 4-byte little-endian block_size prefix to get the record bytes.
    anyhow::ensure!(
        raw.len() >= BLOCK_SIZE_PREFIX_LEN,
        "encode_record_buf_to_raw: encoded output too short ({} bytes)",
        raw.len()
    );
    Ok(crate::RawRecord::from(raw[BLOCK_SIZE_PREFIX_LEN..].to_vec()))
}

/// Decode a single raw BAM record to a noodles `RecordBuf`, using the supplied
/// SAM header for the reference-sequence dictionary.
///
/// Delegates to [`raw_records_to_record_bufs_with_header`] for the single-record
/// case. Passing a non-empty header is required for round-tripping mapped
/// records — without the matching `@SQ` dictionary, reference IDs do not
/// resolve back to names on re-encode.
///
/// # Errors
///
/// Returns an error if the record bytes cannot be decoded (e.g., malformed
/// binary data) or the header cannot be serialized into the decode stream.
pub fn raw_record_to_record_buf(
    rec: &crate::RawRecord,
    header: &noodles::sam::Header,
) -> anyhow::Result<noodles::sam::alignment::RecordBuf> {
    let mut bufs = raw_records_to_record_bufs_with_header(&[rec.as_ref().to_vec()], header)?;
    bufs.pop().ok_or_else(|| anyhow::anyhow!("raw_record_to_record_buf: no record decoded"))
}

/// Amortized `RecordBuf` → `RawRecord` encoder.
///
/// Holds a reference to the SAM header and an internal scratch `Vec<u8>` so
/// that per-call allocations for the encoded output are avoided.  The scratch
/// buffer is reused on every call to [`encode`] and [`encode_into`].
///
/// # Examples
///
/// ```rust,ignore
/// let mut enc = RecordBufEncoder::new(&header);
/// for record in records {
///     let raw = enc.encode(&record)?;
///     // use raw …
/// }
/// ```
pub struct RecordBufEncoder<'h> {
    header: &'h noodles::sam::Header,
    /// Scratch buffer reused across encode calls to amortize allocations.
    scratch: Vec<u8>,
}

impl<'h> RecordBufEncoder<'h> {
    /// Creates a new encoder that encodes against `header`.
    #[must_use]
    pub fn new(header: &'h noodles::sam::Header) -> Self {
        Self { header, scratch: Vec::with_capacity(512) }
    }

    /// Encode `buf` into a freshly-allocated [`RawRecord`].
    ///
    /// The internal scratch buffer is reused for the encoding step, so no
    /// allocation is needed there.  The returned [`RawRecord`] owns its bytes.
    ///
    /// # Errors
    ///
    /// Returns an error if encoding fails.
    pub fn encode(
        &mut self,
        buf: &noodles::sam::alignment::RecordBuf,
    ) -> anyhow::Result<crate::RawRecord> {
        use noodles::sam::alignment::io::Write as _;

        self.scratch.clear();
        let mut writer = noodles::bam::io::Writer::from(&mut self.scratch);
        writer
            .write_alignment_record(self.header, buf)
            .map_err(|e| anyhow::anyhow!("RecordBufEncoder::encode: {e}"))?;

        // Skip the 4-byte block_size prefix.
        anyhow::ensure!(
            self.scratch.len() >= BLOCK_SIZE_PREFIX_LEN,
            "RecordBufEncoder::encode: encoded output too short ({} bytes)",
            self.scratch.len()
        );
        Ok(crate::RawRecord::from(self.scratch[BLOCK_SIZE_PREFIX_LEN..].to_vec()))
    }

    /// Encode `buf` into the caller-supplied `out` [`RawRecord`], replacing its
    /// contents.
    ///
    /// Reuses both the internal scratch buffer and the storage inside `out` to
    /// minimise allocations per call.
    ///
    /// # Errors
    ///
    /// Returns an error if encoding fails.
    pub fn encode_into(
        &mut self,
        buf: &noodles::sam::alignment::RecordBuf,
        out: &mut crate::RawRecord,
    ) -> anyhow::Result<()> {
        use noodles::sam::alignment::io::Write as _;

        self.scratch.clear();
        let mut writer = noodles::bam::io::Writer::from(&mut self.scratch);
        writer
            .write_alignment_record(self.header, buf)
            .map_err(|e| anyhow::anyhow!("RecordBufEncoder::encode_into: {e}"))?;

        // Skip the 4-byte block_size prefix.
        anyhow::ensure!(
            self.scratch.len() >= BLOCK_SIZE_PREFIX_LEN,
            "RecordBufEncoder::encode_into: encoded output too short ({} bytes)",
            self.scratch.len()
        );
        let out_vec = out.as_mut_vec();
        out_vec.clear();
        out_vec.extend_from_slice(&self.scratch[BLOCK_SIZE_PREFIX_LEN..]);
        Ok(())
    }
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
    let header = noodles::sam::Header::default();
    raw_records_to_record_bufs_with_header(records, &header)
}

/// Decode raw BAM byte records to noodles `RecordBuf` using the supplied header.
///
/// Serializes `header` into the synthetic BAM stream so mapped records decode
/// against the correct reference-sequence dictionary. Use this variant whenever
/// records may reference `@SQ` entries — e.g. to round-trip mapped records
/// through noodles. For unmapped-only streams, the header-less
/// [`raw_records_to_record_bufs`] is equivalent.
///
/// # Errors
///
/// Returns an error if the header cannot be serialized or if any record fails
/// to parse.
///
/// # Panics
///
/// Panics if any individual record byte length exceeds `u32::MAX`.
pub fn raw_records_to_record_bufs_with_header(
    records: &[Vec<u8>],
    header: &noodles::sam::Header,
) -> anyhow::Result<Vec<noodles::sam::alignment::RecordBuf>> {
    use anyhow::Context as _;
    use std::io::Cursor;

    // Ask noodles to emit a full BAM prelude (magic + header text + ref dict)
    // for the supplied header. Writing zero records and dropping the writer
    // yields just the prelude bytes.
    let mut prelude: Vec<u8> = Vec::new();
    {
        let mut writer = noodles::bam::io::Writer::from(&mut prelude);
        writer.write_header(header).context("writing synthetic BAM header")?;
    }

    let mut bam_data = prelude;
    for raw in records {
        let block_size = u32::try_from(raw.len()).expect("record too large for BAM block_size");
        bam_data.extend_from_slice(&block_size.to_le_bytes());
        bam_data.extend_from_slice(raw);
    }

    let cursor = Cursor::new(bam_data);
    let mut reader = noodles::bam::io::Reader::from(cursor);
    let _ = reader.read_header().context("reading synthetic BAM header")?;

    let mut result = Vec::with_capacity(records.len());
    for record_result in reader.record_bufs(header) {
        result.push(record_result?);
    }

    Ok(result)
}

#[cfg(test)]
#[allow(clippy::identity_op)]
mod tests {
    use super::*;
    use crate::SamTag;
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
        let mut builder = UnmappedSamBuilder::new();
        builder.build_record(b"read1", 0, b"ACGT", &[30, 25, 20, 15]);
        builder.append_string_tag(SamTag::MI, b"1");
        let raw = builder.as_bytes().to_vec();

        let bufs = raw_records_to_record_bufs(&[raw])
            .expect("converting single raw record to RecordBuf should succeed");
        assert_eq!(bufs.len(), 1);
        assert_eq!(bufs[0].name().map(std::convert::AsRef::as_ref), Some(b"read1".as_ref()));
    }

    #[test]
    fn test_raw_records_to_record_bufs_multiple() {
        let mut records = Vec::new();
        for i in 0..3 {
            let mut builder = UnmappedSamBuilder::new();
            let name = format!("read{i}");
            builder.build_record(name.as_bytes(), 0, b"AC", &[30, 25]);
            records.push(builder.as_bytes().to_vec());
        }

        let bufs = raw_records_to_record_bufs(&records)
            .expect("converting multiple raw records to RecordBufs should succeed");
        assert_eq!(bufs.len(), 3);
        for (i, buf) in bufs.iter().enumerate() {
            let expected = format!("read{i}");
            assert_eq!(buf.name().map(<_ as AsRef<[u8]>>::as_ref), Some(expected.as_bytes()));
        }
    }

    #[test]
    fn test_raw_records_to_record_bufs_empty() {
        let bufs =
            raw_records_to_record_bufs(&[]).expect("converting empty record list should succeed");
        assert!(bufs.is_empty());
    }

    // ========================================================================
    // encode_record_buf_to_raw / raw_record_to_record_buf round-trip tests
    // ========================================================================

    #[test]
    fn test_encode_record_buf_to_raw_single() {
        let mut builder = UnmappedSamBuilder::new();
        builder.build_record(b"read1", 0, b"ACGT", &[30, 25, 20, 15]);
        builder.append_string_tag(SamTag::MI, b"1");
        let original_raw = crate::RawRecord::from(builder.as_bytes().to_vec());

        // Decode to RecordBuf then re-encode and compare field-by-field.
        let header = noodles::sam::Header::default();
        let buf = raw_records_to_record_bufs(&[original_raw.as_ref().to_vec()])
            .expect("decode should succeed")
            .pop()
            .expect("should have one record");

        let reencoded = encode_record_buf_to_raw(&buf, &header)
            .expect("encode_record_buf_to_raw should succeed");
        assert_eq!(reencoded.as_ref(), original_raw.as_ref());
    }

    #[test]
    fn test_raw_record_to_record_buf_single() {
        let mut builder = UnmappedSamBuilder::new();
        builder.build_record(b"myread", 0, b"GCTA", &[40, 30, 20, 10]);
        let raw = crate::RawRecord::from(builder.as_bytes().to_vec());

        let header = noodles::sam::Header::default();
        let buf = raw_record_to_record_buf(&raw, &header)
            .expect("raw_record_to_record_buf should succeed");
        assert_eq!(buf.name().map(std::convert::AsRef::as_ref), Some(b"myread".as_ref()));
    }

    // ========================================================================
    // RecordBufEncoder tests
    // ========================================================================

    #[test]
    fn test_record_buf_encoder_encode() {
        let mut builder = UnmappedSamBuilder::new();
        builder.build_record(b"enc1", 0, b"TTTT", &[20, 20, 20, 20]);
        let original_raw = crate::RawRecord::from(builder.as_bytes().to_vec());

        let header = noodles::sam::Header::default();
        let buf = raw_records_to_record_bufs(&[original_raw.as_ref().to_vec()])
            .expect("decode should succeed")
            .pop()
            .expect("one record");

        let mut enc = RecordBufEncoder::new(&header);
        let result = enc.encode(&buf).expect("encode should succeed");
        assert_eq!(result.as_ref(), original_raw.as_ref());
    }

    #[test]
    fn test_record_buf_encoder_encode_into() {
        let mut builder = UnmappedSamBuilder::new();
        builder.build_record(b"enc2", 0, b"CCCC", &[10, 10, 10, 10]);
        let original_raw = crate::RawRecord::from(builder.as_bytes().to_vec());

        let header = noodles::sam::Header::default();
        let buf = raw_records_to_record_bufs(&[original_raw.as_ref().to_vec()])
            .expect("decode should succeed")
            .pop()
            .expect("one record");

        let mut encoder = RecordBufEncoder::new(&header);
        let mut out = crate::RawRecord::new();
        encoder.encode_into(&buf, &mut out).expect("encode_into should succeed");
        assert_eq!(out.as_ref(), original_raw.as_ref());
    }

    #[test]
    fn test_record_buf_encoder_encode_into_reuses_buffer() {
        // Encodes two records into the same `out` RawRecord to verify reuse.
        let header = noodles::sam::Header::default();
        let mut encoder = RecordBufEncoder::new(&header);
        let mut out = crate::RawRecord::new();

        for name in &[b"aa" as &[u8], b"bb"] {
            let mut builder = UnmappedSamBuilder::new();
            builder.build_record(name, 0, b"AC", &[30, 25]);
            let raw = crate::RawRecord::from(builder.as_bytes().to_vec());
            let buf = raw_records_to_record_bufs(&[raw.as_ref().to_vec()])
                .expect("decode should succeed")
                .pop()
                .expect("one record");
            encoder.encode_into(&buf, &mut out).expect("encode_into should succeed");
            // Verify round-trip: re-decode the encoded bytes and check name.
            let decoded = raw_records_to_record_bufs(&[out.as_ref().to_vec()])
                .expect("re-decode should succeed")
                .pop()
                .expect("one record");
            assert_eq!(decoded.name().map(AsRef::as_ref), Some(*name));
        }
    }
}
