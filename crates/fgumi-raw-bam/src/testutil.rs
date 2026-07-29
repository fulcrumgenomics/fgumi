/// Parsed BAM record for test assertions.
///
/// Provides convenient access to fields extracted from raw BAM bytes
/// produced by `UnmappedSamBuilder`. Used only in tests.
pub struct ParsedBamRecord {
    pub name: Vec<u8>,
    pub flag: u16,
    pub bases: Vec<u8>,
    pub quals: Vec<u8>,
    pub aux_data: Vec<u8>,
}

impl ParsedBamRecord {
    /// Parse a single record from raw bytes (without `block_size` prefix).
    #[must_use]
    pub fn from_bytes(data: &[u8]) -> Self {
        let l_read_name = data[8] as usize;
        let n_cigar_op = u16::from_le_bytes([data[12], data[13]]) as usize;
        let flag = u16::from_le_bytes([data[14], data[15]]);
        let l_seq = u32::from_le_bytes([data[16], data[17], data[18], data[19]]) as usize;

        let name_start = 32;
        let name_end = name_start + l_read_name - 1; // exclude NUL
        let name = data[name_start..name_end].to_vec();

        let cigar_start = name_start + l_read_name;
        let seq_start = cigar_start + n_cigar_op * 4;
        let packed_seq_len = l_seq.div_ceil(2);
        let qual_start = seq_start + packed_seq_len;
        let aux_start = qual_start + l_seq;

        // Unpack bases from 4-bit encoding
        let packed = &data[seq_start..seq_start + packed_seq_len];
        let bases = unpack_sequence_for_test(packed, l_seq);

        let quals = data[qual_start..qual_start + l_seq].to_vec();
        let aux_data = data[aux_start..].to_vec();

        Self { name, flag, bases, quals, aux_data }
    }

    /// Parse all records from a `ConsensusOutput` (`block_size`-prefixed concatenation).
    #[must_use]
    pub fn parse_all(data: &[u8]) -> Vec<Self> {
        let mut records = Vec::new();
        let mut offset = 0;
        while offset + 4 <= data.len() {
            let block_size = u32::from_le_bytes([
                data[offset],
                data[offset + 1],
                data[offset + 2],
                data[offset + 3],
            ]) as usize;
            offset += 4;
            records.push(Self::from_bytes(&data[offset..offset + block_size]));
            offset += block_size;
        }
        records
    }

    /// Find a Z-type string tag value by tag name.
    #[must_use]
    pub fn get_string_tag(&self, tag: impl fgumi_tag::AsTagBytes) -> Option<Vec<u8>> {
        let tag = tag.as_tag_bytes();
        find_z_tag_in_aux(&self.aux_data, *tag)
    }

    /// Find an integer tag value (c/s/i type) by tag name.
    #[must_use]
    pub fn get_int_tag(&self, tag: impl fgumi_tag::AsTagBytes) -> Option<i32> {
        let tag = tag.as_tag_bytes();
        find_int_tag_in_aux(&self.aux_data, *tag)
    }

    /// Find a float tag value by tag name.
    #[must_use]
    pub fn get_float_tag(&self, tag: impl fgumi_tag::AsTagBytes) -> Option<f32> {
        let tag = tag.as_tag_bytes();
        find_float_tag_in_aux(&self.aux_data, *tag)
    }

    /// Find a B:s (i16 array) tag value by tag name.
    #[must_use]
    pub fn get_i16_array_tag(&self, tag: impl fgumi_tag::AsTagBytes) -> Option<Vec<i16>> {
        let tag = tag.as_tag_bytes();
        find_i16_array_tag_in_aux(&self.aux_data, *tag)
    }

    /// Find a B:C (u8 array) tag value by tag name.
    #[must_use]
    pub fn get_u8_array_tag(&self, tag: impl fgumi_tag::AsTagBytes) -> Option<Vec<u8>> {
        let tag = tag.as_tag_bytes();
        find_u8_array_tag_in_aux(&self.aux_data, *tag)
    }
}

fn unpack_sequence_for_test(packed: &[u8], l_seq: usize) -> Vec<u8> {
    use crate::sequence::BAM_BASE_TO_ASCII;
    let mut bases = Vec::with_capacity(l_seq);
    for i in 0..l_seq {
        let byte = packed[i / 2];
        let code = if i % 2 == 0 { byte >> 4 } else { byte & 0x0F };
        bases.push(BAM_BASE_TO_ASCII[code as usize]);
    }
    bases
}

fn find_z_tag_in_aux(aux: &[u8], tag: [u8; 2]) -> Option<Vec<u8>> {
    crate::tags::find_string_tag(aux, tag).map(<[u8]>::to_vec)
}

#[expect(clippy::cast_possible_truncation, reason = "test values always fit in i32")]
fn find_int_tag_in_aux(aux: &[u8], tag: [u8; 2]) -> Option<i32> {
    crate::tags::find_int_tag(aux, tag).map(|v| v as i32)
}

fn find_float_tag_in_aux(aux: &[u8], tag: [u8; 2]) -> Option<f32> {
    crate::tags::find_float_tag(aux, tag)
}

fn find_u8_array_tag_in_aux(aux: &[u8], tag: [u8; 2]) -> Option<Vec<u8>> {
    let arr = crate::tags::find_array_tag(aux, tag)?;
    (arr.elem_type == b'C').then(|| arr.data.to_vec())
}

/// Advance `i` by `n` bytes, or `None` if that would run past `len`.
///
/// Every payload step in [`find_i16_array_tag_in_aux`] goes through this so a
/// truncated aux block returns `None` instead of panicking on an out-of-range
/// index.
fn advance_within(i: usize, n: usize, len: usize) -> Option<usize> {
    let end = i.checked_add(n)?;
    (end <= len).then_some(end)
}

/// Scan an aux block for a `B:s` (int16 array) tag.
///
/// Fails closed on malformed input: a truncated scalar, an unterminated `Z`/`H`
/// string, a `B` header or payload that runs past the end of `aux`, an unknown
/// `B` subtype (whose element width is unknowable, so the cursor cannot be
/// advanced safely), or an unknown type code all yield `None` rather than
/// indexing out of bounds.
fn find_i16_array_tag_in_aux(aux: &[u8], tag: [u8; 2]) -> Option<Vec<i16>> {
    let mut i = 0;
    while i + 3 <= aux.len() {
        let t = [aux[i], aux[i + 1]];
        let typ = aux[i + 2];
        i += 3;
        match typ {
            b'B' => {
                // Subtype byte + 4-byte little-endian element count.
                let header_end = advance_within(i, 5, aux.len())?;
                let sub = aux[i];
                let count =
                    u32::from_le_bytes([aux[i + 1], aux[i + 2], aux[i + 3], aux[i + 4]]) as usize;
                i = header_end;

                let elem_size = match sub {
                    b'c' | b'C' => 1,
                    b's' | b'S' => 2,
                    b'i' | b'I' | b'f' => 4,
                    // Not a BAM array subtype: the payload width is unknown, so
                    // skipping it would resync the cursor onto garbage.
                    _ => return None,
                };
                let payload_end = advance_within(i, count.checked_mul(elem_size)?, aux.len())?;

                if t == tag && sub == b's' {
                    return Some(
                        aux[i..payload_end]
                            .chunks_exact(2)
                            .map(|c| i16::from_le_bytes([c[0], c[1]]))
                            .collect(),
                    );
                }
                i = payload_end;
            }
            b'A' | b'c' | b'C' => i = advance_within(i, 1, aux.len())?,
            b's' | b'S' => i = advance_within(i, 2, aux.len())?,
            b'i' | b'I' | b'f' => i = advance_within(i, 4, aux.len())?,
            b'Z' | b'H' => {
                // NUL-terminated; a missing terminator is malformed, not an
                // implicit end-of-block.
                let nul = aux.get(i..)?.iter().position(|&b| b == 0)?;
                i = advance_within(i, nul + 1, aux.len())?;
            }
            _ => return None,
        }
    }
    None
}

/// Construct a raw BAM byte array for testing.
///
/// IMPORTANT: `name` length + 1 (for null terminator) should be divisible
/// by 4 to maintain alignment for CIGAR ops.  Use names like b"rea" (3+1=4)
/// or b"readABC" (7+1=8).
///
/// # Panics
///
/// Panics if `name` length exceeds 254 bytes, `cigar_ops` length exceeds `u16::MAX`,
/// or `seq_len` exceeds `u32::MAX`.
#[must_use]
#[allow(clippy::too_many_arguments)]
pub fn make_bam_bytes(
    tid: i32,
    pos: i32,
    flag: u16,
    name: &[u8],
    cigar_ops: &[u32],
    seq_len: usize,
    mate_tid: i32,
    mate_pos: i32,
    aux_data: &[u8],
) -> Vec<u8> {
    let l_read_name = u8::try_from(name.len() + 1).expect("name length + 1 must fit in u8"); // +1 for null terminator
    let n_cigar_op = u16::try_from(cigar_ops.len()).expect("cigar_ops length must fit in u16");
    let seq_bytes = seq_len.div_ceil(2);
    let total =
        32 + l_read_name as usize + cigar_ops.len() * 4 + seq_bytes + seq_len + aux_data.len();
    let mut buf = vec![0u8; total];

    // Fixed header fields
    buf[0..4].copy_from_slice(&tid.to_le_bytes());
    buf[4..8].copy_from_slice(&pos.to_le_bytes());
    buf[8] = l_read_name;
    buf[9] = 0; // mapq
    buf[10..12].copy_from_slice(&0u16.to_le_bytes()); // bin
    buf[12..14].copy_from_slice(&n_cigar_op.to_le_bytes());
    buf[14..16].copy_from_slice(&flag.to_le_bytes());
    buf[16..20]
        .copy_from_slice(&u32::try_from(seq_len).expect("seq_len must fit in u32").to_le_bytes());
    buf[20..24].copy_from_slice(&mate_tid.to_le_bytes());
    buf[24..28].copy_from_slice(&mate_pos.to_le_bytes());
    buf[28..32].copy_from_slice(&0i32.to_le_bytes()); // tlen

    // Read name + null terminator
    let name_start = 32;
    buf[name_start..name_start + name.len()].copy_from_slice(name);
    buf[name_start + name.len()] = 0; // null terminator

    // CIGAR ops
    let cigar_start = name_start + l_read_name as usize;
    for (i, &op) in cigar_ops.iter().enumerate() {
        let offset = cigar_start + i * 4;
        buf[offset..offset + 4].copy_from_slice(&op.to_le_bytes());
    }

    // Sequence bytes (all zeros) and quality bytes (all zeros) are already zero

    // Aux data
    let aux_start = cigar_start + cigar_ops.len() * 4 + seq_bytes + seq_len;
    buf[aux_start..aux_start + aux_data.len()].copy_from_slice(aux_data);

    buf
}

/// Build a `make_bam_bytes` record with a custom template length (tlen).
#[must_use]
#[allow(clippy::too_many_arguments)]
pub fn make_bam_bytes_with_tlen(
    tid: i32,
    pos: i32,
    flag: u16,
    name: &[u8],
    cigar_ops: &[u32],
    seq_len: usize,
    mate_tid: i32,
    mate_pos: i32,
    tlen: i32,
    aux_data: &[u8],
) -> Vec<u8> {
    let mut rec =
        make_bam_bytes(tid, pos, flag, name, cigar_ops, seq_len, mate_tid, mate_pos, aux_data);
    rec[28..32].copy_from_slice(&tlen.to_le_bytes());
    rec
}

/// Encode a single CIGAR op.  `op_type`: M=0, I=1, D=2, N=3, S=4, H=5, P=6, `=7`, X=8.
///
/// # Panics
///
/// Panics if `len` exceeds `u32::MAX`.
#[must_use]
pub fn encode_op(op_type: u32, len: usize) -> u32 {
    (u32::try_from(len).expect("CIGAR op length must fit in u32") << 4) | op_type
}

/// Helper: build raw aux bytes for a B-type array tag.
#[must_use]
pub fn make_b_array_tag(tag: [u8; 2], elem_type: u8, count: u32, elements: &[u8]) -> Vec<u8> {
    let mut aux = vec![tag[0], tag[1], b'B', elem_type];
    aux.extend_from_slice(&count.to_le_bytes());
    aux.extend_from_slice(elements);
    aux
}

/// # Panics
///
/// Panics if `values` length exceeds `u32::MAX`.
#[must_use]
pub fn make_b_int_array_tag(tag: [u8; 2], values: &[i32]) -> Vec<u8> {
    let bytes: Vec<u8> = values.iter().flat_map(|v| v.to_le_bytes()).collect();
    make_b_array_tag(
        tag,
        b'i',
        u32::try_from(values.len()).expect("array length must fit in u32"),
        &bytes,
    )
}

/// # Panics
///
/// Panics if `values` length exceeds `u32::MAX`.
#[must_use]
pub fn make_b_float_array_tag(tag: [u8; 2], values: &[f32]) -> Vec<u8> {
    let bytes: Vec<u8> = values.iter().flat_map(|v| v.to_le_bytes()).collect();
    make_b_array_tag(
        tag,
        b'f',
        u32::try_from(values.len()).expect("array length must fit in u32"),
        &bytes,
    )
}

/// # Panics
///
/// Panics if `values` length exceeds `u32::MAX`.
#[must_use]
pub fn make_b_uint8_array_tag(tag: [u8; 2], values: &[u8]) -> Vec<u8> {
    make_b_array_tag(
        tag,
        b'C',
        u32::try_from(values.len()).expect("array length must fit in u32"),
        values,
    )
}

/// # Panics
///
/// Panics if `values` length exceeds `u32::MAX`.
#[must_use]
pub fn make_b_int8_array_tag(tag: [u8; 2], values: &[i8]) -> Vec<u8> {
    let bytes: Vec<u8> = values.iter().map(|&v| v.cast_unsigned()).collect();
    make_b_array_tag(
        tag,
        b'c',
        u32::try_from(values.len()).expect("array length must fit in u32"),
        &bytes,
    )
}

/// # Panics
///
/// Panics if `values` length exceeds `u32::MAX`.
#[must_use]
pub fn make_b_int16_array_tag(tag: [u8; 2], values: &[i16]) -> Vec<u8> {
    let bytes: Vec<u8> = values.iter().flat_map(|v| v.to_le_bytes()).collect();
    make_b_array_tag(
        tag,
        b's',
        u32::try_from(values.len()).expect("array length must fit in u32"),
        &bytes,
    )
}

/// # Panics
///
/// Panics if `values` length exceeds `u32::MAX`.
#[must_use]
pub fn make_b_uint16_array_tag(tag: [u8; 2], values: &[u16]) -> Vec<u8> {
    let bytes: Vec<u8> = values.iter().flat_map(|v| v.to_le_bytes()).collect();
    make_b_array_tag(
        tag,
        b'S',
        u32::try_from(values.len()).expect("array length must fit in u32"),
        &bytes,
    )
}

/// # Panics
///
/// Panics if `values` length exceeds `u32::MAX`.
#[must_use]
pub fn make_b_uint32_array_tag(tag: [u8; 2], values: &[u32]) -> Vec<u8> {
    let bytes: Vec<u8> = values.iter().flat_map(|v| v.to_le_bytes()).collect();
    make_b_array_tag(
        tag,
        b'I',
        u32::try_from(values.len()).expect("array length must fit in u32"),
        &bytes,
    )
}

#[cfg(test)]
mod tests {
    use super::*;
    use proptest::prelude::*;
    use rstest::rstest;

    /// A well-formed `B:s` array is still found after the bounds hardening.
    #[test]
    fn finds_a_well_formed_i16_array() {
        // MI:i:7 then CD:B:s,[1,2]
        let mut aux = vec![b'M', b'I', b'i'];
        aux.extend_from_slice(&7i32.to_le_bytes());
        aux.extend_from_slice(b"CDBs");
        aux.extend_from_slice(&2u32.to_le_bytes());
        aux.extend_from_slice(&1i16.to_le_bytes());
        aux.extend_from_slice(&2i16.to_le_bytes());

        assert_eq!(find_i16_array_tag_in_aux(&aux, [b'C', b'D']), Some(vec![1, 2]));
        assert_eq!(find_i16_array_tag_in_aux(&aux, [b'X', b'X']), None);
    }

    /// Malformed aux data must return `None`, never panic.
    ///
    /// The scanner walks caller-supplied bytes, so every payload step has to be
    /// bounds-checked; before hardening each of these indexed past the end of the
    /// slice and panicked.
    #[rstest]
    #[case::truncated_b_header(vec![b'C', b'D', b'B'])]
    #[case::truncated_b_count(vec![b'C', b'D', b'B', b's', 0x02, 0x00])]
    #[case::b_payload_shorter_than_count({
        let mut v = vec![b'C', b'D', b'B', b's'];
        v.extend_from_slice(&9u32.to_le_bytes()); // claims 9 elements
        v.extend_from_slice(&1i16.to_le_bytes()); // supplies 1
        v
    })]
    #[case::b_count_overflows_usize({
        let mut v = vec![b'C', b'D', b'B', b'i'];
        v.extend_from_slice(&u32::MAX.to_le_bytes());
        v
    })]
    #[case::unknown_b_subtype({
        let mut v = vec![b'C', b'D', b'B', b'?'];
        v.extend_from_slice(&1u32.to_le_bytes());
        v.push(0);
        v
    })]
    #[case::truncated_scalar(vec![b'N', b'M', b'i', 0x01, 0x02])]
    #[case::unterminated_z_string(vec![b'R', b'X', b'Z', b'A', b'C', b'G'])]
    #[case::unknown_type_code(vec![b'Z', b'Z', b'?', 0x00])]
    fn malformed_aux_returns_none_instead_of_panicking(#[case] aux: Vec<u8>) {
        assert_eq!(find_i16_array_tag_in_aux(&aux, [b'C', b'D']), None);
    }

    proptest! {
        /// The "fail closed, never read past the buffer" contract, over arbitrary bytes.
        ///
        /// The case table above pins the truncation shapes worth naming; this covers the
        /// ones nobody thought of. The scanner walks caller-supplied bytes with a cursor
        /// it advances by type-dependent widths, which is exactly the shape where a
        /// hand-picked corpus misses a case — so the property is simply that no input
        /// panics, and that a returned array is self-consistent.
        #[test]
        fn arbitrary_aux_bytes_never_panic(aux in proptest::collection::vec(any::<u8>(), 0..256)) {
            // Must not panic for any tag, present or absent.
            let found = find_i16_array_tag_in_aux(&aux, [b'C', b'D']);
            let _ = find_i16_array_tag_in_aux(&aux, [b'M', b'I']);

            // Any array we do return must fit inside the block it was parsed from:
            // two bytes per element cannot exceed the input length.
            if let Some(vals) = found {
                prop_assert!(
                    vals.len() * 2 <= aux.len(),
                    "returned {} i16 values from a {}-byte aux block",
                    vals.len(),
                    aux.len(),
                );
            }
        }

        /// A well-formed `B:s` array is always recovered exactly, whatever the values.
        ///
        /// Guards the hardening against the opposite failure: bounds checks that are too
        /// strict and start rejecting valid input.
        #[test]
        fn well_formed_i16_arrays_round_trip(vals in proptest::collection::vec(any::<i16>(), 0..64)) {
            let mut aux: Vec<u8> = vec![b'C', b'D', b'B', b's'];
            aux.extend_from_slice(&u32::try_from(vals.len()).unwrap().to_le_bytes());
            for v in &vals {
                aux.extend_from_slice(&v.to_le_bytes());
            }
            prop_assert_eq!(find_i16_array_tag_in_aux(&aux, [b'C', b'D']), Some(vals));
        }
    }
}
