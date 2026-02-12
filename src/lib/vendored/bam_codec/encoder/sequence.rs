use std::io;

use noodles::sam::alignment::record::Sequence;

use super::num::{write_u8, write_u32_le};

pub(super) fn write_length(dst: &mut Vec<u8>, base_count: usize) -> io::Result<()> {
    let n =
        u32::try_from(base_count).map_err(|e| io::Error::new(io::ErrorKind::InvalidInput, e))?;

    write_u32_le(dst, n);

    Ok(())
}

#[allow(clippy::needless_pass_by_value)]
pub(super) fn write_sequence<S>(
    dst: &mut Vec<u8>,
    read_length: usize,
    sequence: S,
) -> io::Result<()>
where
    S: Sequence,
{
    const EQ: u8 = b'=';

    if sequence.is_empty() {
        return Ok(());
    }

    // § 1.4.10 "`SEQ`" (2022-08-22): "If not a '*', the length of the sequence must equal the sum
    // of lengths of `M`/`I`/`S`/`=`/`X` operations in `CIGAR`."
    if read_length > 0 && sequence.len() != read_length {
        return Err(io::Error::new(
            io::ErrorKind::InvalidInput,
            "read length-sequence length mismatch",
        ));
    }

    let mut bases = sequence.iter();

    while let Some(l) = bases.next() {
        // § 4.2.3 "SEQ and QUAL encoding" (2021-06-03): "When `l_seq` is odd the bottom 4 bits of
        // the last byte are undefined, but we recommend writing these as zero."
        let r = bases.next().unwrap_or(EQ);
        let n = (encode_base(l) << 4) | encode_base(r);
        write_u8(dst, n);
    }

    Ok(())
}

/// Fast path: encode sequence from a slice (no dynamic dispatch).
///
/// Uses 16-base chunked processing (matching htslib's NN=16 strategy in sam.c)
/// for better cache locality and reduced loop overhead.
#[inline]
pub(super) fn write_sequence_from_slice(
    dst: &mut Vec<u8>,
    read_length: usize,
    bases: &[u8],
) -> io::Result<()> {
    // Process 16 bases (8 output bytes) at a time for cache efficiency
    // This matches htslib's NN=16 chunking strategy (sam.c:621-636)
    const CHUNK_SIZE: usize = 16;

    if bases.is_empty() {
        return Ok(());
    }

    if read_length > 0 && bases.len() != read_length {
        return Err(io::Error::new(
            io::ErrorKind::InvalidInput,
            "read length-sequence length mismatch",
        ));
    }

    // Reserve space for packed output
    let packed_len = bases.len().div_ceil(2);
    dst.reserve(packed_len);

    let mut chunks = bases.chunks_exact(CHUNK_SIZE);

    for chunk in chunks.by_ref() {
        // Unrolled loop processes 16 bases -> 8 bytes
        // Helps CPU pipelining and may enable auto-vectorization
        dst.push((CODES[chunk[0] as usize] << 4) | CODES[chunk[1] as usize]);
        dst.push((CODES[chunk[2] as usize] << 4) | CODES[chunk[3] as usize]);
        dst.push((CODES[chunk[4] as usize] << 4) | CODES[chunk[5] as usize]);
        dst.push((CODES[chunk[6] as usize] << 4) | CODES[chunk[7] as usize]);
        dst.push((CODES[chunk[8] as usize] << 4) | CODES[chunk[9] as usize]);
        dst.push((CODES[chunk[10] as usize] << 4) | CODES[chunk[11] as usize]);
        dst.push((CODES[chunk[12] as usize] << 4) | CODES[chunk[13] as usize]);
        dst.push((CODES[chunk[14] as usize] << 4) | CODES[chunk[15] as usize]);
    }

    // Handle remainder (< 16 bases) with 2-base pairs
    let remainder = chunks.remainder();
    let mut pairs = remainder.chunks_exact(2);
    for pair in pairs.by_ref() {
        let l = CODES[pair[0] as usize];
        let r = CODES[pair[1] as usize];
        dst.push((l << 4) | r);
    }

    // Handle final odd base if present
    if let Some(&last) = pairs.remainder().first() {
        let l = CODES[last as usize];
        dst.push(l << 4); // Lower 4 bits are zero (padding)
    }

    Ok(())
}

/// Lookup table for base encoding (compile-time generated).
const CODES: [u8; 256] = build_codes();

// § 4.2.3 "SEQ and QUAL encoding" (2023-11-16): "The case-insensitive base codes [...] are mapped
// to [0, 15] respectively with all other characters mapping to 'N' (value 15)".
#[inline]
fn encode_base(n: u8) -> u8 {
    CODES[usize::from(n)]
}

#[allow(clippy::cast_possible_truncation)]
const fn build_codes() -> [u8; 256] {
    // § 4.2.3 "SEQ and QUAL encoding" (2024-11-06)
    const BASES: [u8; 16] = *b"=ACMGRSVTWYHKDBN";
    const N: u8 = 0x0f;

    let mut i = 0;
    let mut codes = [N; 256];

    while i < BASES.len() {
        let base = BASES[i];

        // SAFETY: `i < BASES.len() <= u8::MAX`.
        let code = i as u8;

        // SAFETY: `base <= codes.len() <= u8::MAX`.
        codes[base as usize] = code;
        codes[base.to_ascii_lowercase() as usize] = code;

        i += 1;
    }

    codes
}

#[cfg(test)]
mod tests {
    use noodles::sam::alignment::record_buf::Sequence as SequenceBuf;

    use super::*;

    #[test]
    fn test_write_length() -> io::Result<()> {
        let mut buf = Vec::new();
        write_length(&mut buf, 8)?;
        assert_eq!(buf, [0x08, 0x00, 0x00, 0x00]);
        Ok(())
    }

    #[cfg(not(target_pointer_width = "16"))]
    #[test]
    fn test_write_length_with_out_of_range_base_count() {
        let mut buf = Vec::new();

        assert!(matches!(
            write_length(&mut buf, usize::MAX),
            Err(e) if e.kind() == io::ErrorKind::InvalidInput
        ));
    }

    #[test]
    fn test_write_sequence() -> Result<(), Box<dyn std::error::Error>> {
        fn t(buf: &mut Vec<u8>, sequence: &SequenceBuf, expected: &[u8]) -> io::Result<()> {
            buf.clear();
            write_sequence(buf, sequence.len(), sequence)?;
            assert_eq!(buf, expected);
            Ok(())
        }

        let mut buf = Vec::new();

        t(&mut buf, &SequenceBuf::default(), &[])?;
        t(&mut buf, &SequenceBuf::from(b"ACG"), &[0x12, 0x40])?;
        t(&mut buf, &SequenceBuf::from(b"ACGT"), &[0x12, 0x48])?;

        buf.clear();
        write_sequence(&mut buf, 2, &SequenceBuf::default())?;
        assert!(buf.is_empty());

        buf.clear();
        let sequence = SequenceBuf::from(b"A");
        assert!(matches!(
            write_sequence(&mut buf, 2, &sequence),
            Err(e) if e.kind() == io::ErrorKind::InvalidInput,
        ));

        Ok(())
    }

    #[test]
    fn test_encode_base() {
        const BASES: [u8; 16] = *b"=ACMGRSVTWYHKDBN";

        for (i, b) in (0..).zip(BASES) {
            assert_eq!(encode_base(b), i);
            assert_eq!(encode_base(b.to_ascii_lowercase()), i);
        }

        assert_eq!(encode_base(b'X'), 15);
        assert_eq!(encode_base(b'x'), 15);
    }

    #[test]
    fn test_write_sequence_from_slice() {
        let mut buf = Vec::new();

        // Empty sequence
        write_sequence_from_slice(&mut buf, 0, b"").unwrap();
        assert!(buf.is_empty());

        // Even length (4 bases -> 2 bytes)
        buf.clear();
        write_sequence_from_slice(&mut buf, 4, b"ACGT").unwrap();
        assert_eq!(buf, [0x12, 0x48]); // A=1, C=2, G=4, T=8

        // Odd length (3 bases -> 2 bytes, last nibble padded)
        buf.clear();
        write_sequence_from_slice(&mut buf, 3, b"ACG").unwrap();
        assert_eq!(buf, [0x12, 0x40]); // A=1, C=2, G=4, pad=0

        // Single base
        buf.clear();
        write_sequence_from_slice(&mut buf, 1, b"A").unwrap();
        assert_eq!(buf, [0x10]); // A=1, pad=0

        // Length mismatch error
        buf.clear();
        assert!(write_sequence_from_slice(&mut buf, 2, b"A").is_err());

        // Verify matches trait-based version
        buf.clear();
        write_sequence_from_slice(&mut buf, 4, b"ACGT").unwrap();
        let mut buf2 = Vec::new();
        write_sequence(&mut buf2, 4, &SequenceBuf::from(b"ACGT")).unwrap();
        assert_eq!(buf, buf2);
    }
}
