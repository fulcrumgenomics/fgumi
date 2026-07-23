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

fn find_i16_array_tag_in_aux(aux: &[u8], tag: [u8; 2]) -> Option<Vec<i16>> {
    let mut i = 0;
    while i + 3 <= aux.len() {
        let t = [aux[i], aux[i + 1]];
        let typ = aux[i + 2];
        i += 3;
        match typ {
            b'B' => {
                // Fail closed on a truncated or unrecognized array rather than
                // indexing past the end, mirroring `tags::parse_array_tag_at`.
                let header = aux.get(i..i + 5)?;
                let sub = header[0];
                let count =
                    u32::from_le_bytes([header[1], header[2], header[3], header[4]]) as usize;
                i += 5;
                let elem_size = match sub {
                    b'c' | b'C' => 1,
                    b's' | b'S' => 2,
                    b'i' | b'I' | b'f' => 4,
                    // An unknown subtype has no known width, so the walker cannot
                    // resynchronize past it.
                    _ => return None,
                };
                let payload = aux.get(i..i.checked_add(count.checked_mul(elem_size)?)?)?;
                i += payload.len();
                if t == tag && sub == b's' {
                    return Some(
                        payload
                            .chunks_exact(2)
                            .map(|elem| i16::from_le_bytes([elem[0], elem[1]]))
                            .collect(),
                    );
                }
            }
            b'A' | b'c' | b'C' => {
                i += 1;
            }
            b's' | b'S' => {
                i += 2;
            }
            b'i' | b'I' | b'f' => {
                i += 4;
            }
            b'Z' | b'H' => {
                while i < aux.len() && aux[i] != 0 {
                    i += 1;
                }
                i += 1;
            }
            _ => break,
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

    /// Append a tag of the given type code to an aux block.
    fn push_tag(aux: &mut Vec<u8>, tag: [u8; 2], typ: u8, payload: &[u8]) {
        aux.extend_from_slice(&tag);
        aux.push(typ);
        aux.extend_from_slice(payload);
    }

    /// Append a `B:s` (int16 array) tag.
    fn push_i16_array(aux: &mut Vec<u8>, tag: [u8; 2], values: &[i16]) {
        aux.extend_from_slice(&tag);
        aux.push(b'B');
        aux.push(b's');
        aux.extend_from_slice(&u32::try_from(values.len()).expect("count fits u32").to_le_bytes());
        for v in values {
            aux.extend_from_slice(&v.to_le_bytes());
        }
    }

    /// The aux walker has to skip a tag of every BAM type code to reach a later
    /// tag. If it advances by the wrong width for any of them it desynchronizes
    /// and reads garbage as tag identifiers -- which, in a test oracle, means
    /// real tag assertions can silently pass against the wrong bytes. Each case
    /// places one preceding tag type in front of the target array.
    #[rstest::rstest]
    #[case::char_a(b'A', b"x")]
    #[case::int8_c(b'c', b"\xFF")]
    #[case::uint8_upper_c(b'C', b"\xC8")]
    #[case::int16_s(b's', &[0x34, 0x12])]
    #[case::uint16_upper_s(b'S', &[0xFF, 0xFF])]
    #[case::int32_i(b'i', &[1, 0, 0, 0])]
    #[case::uint32_upper_i(b'I', &[0xFF, 0xFF, 0xFF, 0xFF])]
    #[case::float_f(b'f', &[0, 0, 0x80, 0x3F])]
    #[case::string_z(b'Z', b"hello\0")]
    #[case::hex_h(b'H', b"1AE301\0")]
    fn test_aux_walker_skips_every_type_code(#[case] typ: u8, #[case] payload: &[u8]) {
        let expected = vec![10i16, -20, 30];

        let mut aux = Vec::new();
        push_tag(&mut aux, *fgumi_tag::SamTag::MI, typ, payload);
        push_i16_array(&mut aux, *fgumi_tag::SamTag::AD_BASES, &expected);

        assert_eq!(
            find_i16_array_tag_in_aux(&aux, *fgumi_tag::SamTag::AD_BASES),
            Some(expected),
            "walker desynchronized while skipping a '{}' tag",
            typ as char,
        );
    }

    /// Arrays of every element type must also be skipped at the right width.
    #[rstest::rstest]
    #[case::int8_array(b'c', 1)]
    #[case::uint8_array(b'C', 1)]
    #[case::int16_array(b's', 2)]
    #[case::uint16_array(b'S', 2)]
    #[case::int32_array(b'i', 4)]
    #[case::uint32_array(b'I', 4)]
    #[case::float_array(b'f', 4)]
    fn test_aux_walker_skips_arrays_of_every_element_type(
        #[case] subtype: u8,
        #[case] elem_size: usize,
    ) {
        let expected = vec![7i16, 8];
        let count = 3usize;

        let mut aux = Vec::new();
        // A non-matching array tag of the given subtype, then the target.
        aux.extend_from_slice(&*fgumi_tag::SamTag::MI);
        aux.push(b'B');
        aux.push(subtype);
        aux.extend_from_slice(&u32::try_from(count).expect("fits").to_le_bytes());
        aux.extend(std::iter::repeat_n(0u8, count * elem_size));
        push_i16_array(&mut aux, *fgumi_tag::SamTag::AD_BASES, &expected);

        assert_eq!(
            find_i16_array_tag_in_aux(&aux, *fgumi_tag::SamTag::AD_BASES),
            Some(expected),
            "walker desynchronized while skipping a 'B:{}' array",
            subtype as char,
        );
    }

    /// A same-subtype array whose tag does not match must be skipped rather than
    /// returned -- the `B:s` branch returns early only on a tag match.
    #[test]
    fn test_aux_walker_skips_non_matching_i16_array() {
        let expected = vec![1i16, 2, 3];
        let mut aux = Vec::new();
        push_i16_array(&mut aux, *fgumi_tag::SamTag::BD_BASES, &[99, 98]);
        push_i16_array(&mut aux, *fgumi_tag::SamTag::AD_BASES, &expected);

        assert_eq!(find_i16_array_tag_in_aux(&aux, *fgumi_tag::SamTag::AD_BASES), Some(expected));
    }

    #[test]
    fn test_aux_walker_returns_none_for_absent_tag() {
        let mut aux = Vec::new();
        push_i16_array(&mut aux, *fgumi_tag::SamTag::BD_BASES, &[1, 2]);
        assert_eq!(find_i16_array_tag_in_aux(&aux, *fgumi_tag::SamTag::AD_BASES), None);
    }

    /// Malformed `B` encodings must fail closed rather than index past the end
    /// of the aux block. Each case truncates the target array after a different
    /// part of its encoding.
    #[rstest::rstest]
    #[case::truncated_before_subtype(0)]
    #[case::truncated_mid_count(3)]
    #[case::truncated_before_payload(5)]
    #[case::truncated_mid_payload(8)]
    fn test_aux_walker_fails_closed_on_truncated_array(#[case] kept_after_type_byte: usize) {
        let mut aux = Vec::new();
        push_i16_array(&mut aux, *fgumi_tag::SamTag::AD_BASES, &[1i16, 2, 3]);
        // Keep the tag name and `B` type byte, then cut the encoding short.
        aux.truncate(3 + kept_after_type_byte);

        assert_eq!(find_i16_array_tag_in_aux(&aux, *fgumi_tag::SamTag::AD_BASES), None);
    }

    #[test]
    fn test_aux_walker_fails_closed_on_unknown_array_subtype() {
        let mut aux = Vec::new();
        aux.extend_from_slice(&*fgumi_tag::SamTag::MI);
        aux.push(b'B');
        aux.push(b'?');
        aux.extend_from_slice(&2u32.to_le_bytes());
        aux.extend_from_slice(&[0, 0]);
        // A well-formed target after the bad array is unreachable: the walker
        // cannot know how wide the unknown subtype is, so it must stop.
        push_i16_array(&mut aux, *fgumi_tag::SamTag::AD_BASES, &[1i16, 2]);

        assert_eq!(find_i16_array_tag_in_aux(&aux, *fgumi_tag::SamTag::AD_BASES), None);
    }

    /// An unterminated `Z`/`H` payload has no NUL to stop the scan on, so the
    /// walk must end at the end of the block rather than running off it.
    #[rstest::rstest]
    #[case::unterminated_string(b'Z', b"nonul")]
    #[case::unterminated_hex(b'H', b"1AE3")]
    fn test_aux_walker_fails_closed_on_malformed_scalar(#[case] typ: u8, #[case] payload: &[u8]) {
        let mut aux = Vec::new();
        push_tag(&mut aux, *fgumi_tag::SamTag::MI, typ, payload);

        assert_eq!(find_i16_array_tag_in_aux(&aux, *fgumi_tag::SamTag::AD_BASES), None);
    }

    /// An unknown scalar type code has no width to advance by, so the walk must
    /// stop there rather than guess one. An aux block that simply ends after the
    /// bad tag would return `None` either way; a well-formed target placed after
    /// it is what proves the walker stopped instead of resynchronizing onto it.
    #[test]
    fn test_aux_walker_stops_at_unknown_scalar_type_code() {
        let mut aux = Vec::new();
        push_tag(&mut aux, *fgumi_tag::SamTag::MI, b'?', &[0, 0, 0, 0]);
        push_i16_array(&mut aux, *fgumi_tag::SamTag::AD_BASES, &[1i16, 2]);

        assert_eq!(find_i16_array_tag_in_aux(&aux, *fgumi_tag::SamTag::AD_BASES), None);
    }

    /// Trailing bytes too short to hold a tag must end the walk without being
    /// read as one -- the target found earlier still comes back. Querying an
    /// absent tag is what actually reaches the remainder: a match returns from
    /// inside the `B` arm before the outer bounds guard is re-tested, so only a
    /// miss carries the walk into the partial header.
    #[test]
    fn test_aux_walker_ignores_trailing_partial_tag() {
        let expected = vec![5i16];
        let mut aux = Vec::new();
        push_i16_array(&mut aux, *fgumi_tag::SamTag::AD_BASES, &expected);
        aux.extend_from_slice(&*fgumi_tag::SamTag::MI);

        assert_eq!(find_i16_array_tag_in_aux(&aux, *fgumi_tag::SamTag::AD_BASES), Some(expected));
        assert_eq!(find_i16_array_tag_in_aux(&aux, *fgumi_tag::SamTag::BD_BASES), None);
    }

    /// A count that overflows the aux block must be rejected on the *matching*
    /// tag too -- that is the branch that reads the payload out.
    #[test]
    fn test_aux_walker_fails_closed_on_overlong_count() {
        let mut aux = Vec::new();
        aux.extend_from_slice(&*fgumi_tag::SamTag::AD_BASES);
        aux.push(b'B');
        aux.push(b's');
        aux.extend_from_slice(&u32::MAX.to_le_bytes());
        aux.extend_from_slice(&[1, 0, 2, 0]);

        assert_eq!(find_i16_array_tag_in_aux(&aux, *fgumi_tag::SamTag::AD_BASES), None);
    }

    // ========================================================================
    // proptest: the aux walker fails closed on arbitrary and truncated bytes
    //
    // The cases above pin hand-picked malformed encodings; random buffers cover
    // the offsets nobody thought to enumerate (partial trailing tags of every
    // length, truncations behind an arbitrary run of preceding tags).
    // ========================================================================

    /// Fixed-width scalar type codes and their payload widths.
    const SCALAR_TYPES: [(u8, usize); 8] =
        [(b'A', 1), (b'c', 1), (b'C', 1), (b's', 2), (b'S', 2), (b'i', 4), (b'I', 4), (b'f', 4)];

    /// `B` array element types and their element widths.
    const ARRAY_ELEMENT_TYPES: [(u8, usize); 7] =
        [(b'c', 1), (b'C', 1), (b's', 2), (b'S', 2), (b'i', 4), (b'I', 4), (b'f', 4)];

    /// Every tag encoding [`encode_tag`] can emit: one per scalar type code, one
    /// per array element type, then `Z` and `H`. Derived from the tables so a
    /// new type code widens the generated space instead of going untested.
    const TAG_KIND_COUNT: usize = SCALAR_TYPES.len() + ARRAY_ELEMENT_TYPES.len() + 2;

    /// Payload bytes supplied to [`encode_tag`], sized for the widest tag it
    /// emits: a two-element `B:i`/`B:f` array.
    const FILLER_LEN: usize = 8;

    /// Encode one well-formed aux tag of the `kind`-th type code, drawing its
    /// payload from `filler`.
    ///
    /// No byte this emits is `b'a'` or `b'd'`, so no window of the encoding can
    /// spell the `ad` tag name. That is what lets the truncation property below
    /// assert the walker returns *the planted array or nothing*: a walk that
    /// desynchronizes into these bytes cannot resynchronize onto a header that
    /// looks like the target.
    fn encode_tag(name: [u8; 2], kind: usize, filler: [u8; FILLER_LEN]) -> Vec<u8> {
        let mut out = name.to_vec();
        if let Some(&(typ, width)) = SCALAR_TYPES.get(kind) {
            out.push(typ);
            out.extend_from_slice(&filler[..width]);
        } else if let Some(&(subtype, width)) = ARRAY_ELEMENT_TYPES.get(kind - SCALAR_TYPES.len()) {
            let count = usize::from(filler[0]) % 3;
            out.push(b'B');
            out.push(subtype);
            out.extend_from_slice(&u32::try_from(count).expect("count fits u32").to_le_bytes());
            out.extend_from_slice(&filler[..count * width]);
        } else {
            let terminated = if kind == TAG_KIND_COUNT - 2 { b'Z' } else { b'H' };
            out.push(terminated);
            out.extend_from_slice(&filler[..usize::from(filler[0]) % 5]);
            out.push(0);
        }
        out
    }

    /// Bytes that make random buffers actually look like aux data: mostly
    /// structural bytes the walker branches on, with arbitrary bytes mixed in.
    fn aux_ish_byte() -> impl proptest::strategy::Strategy<Value = u8> {
        proptest::prop_oneof![
            4 => proptest::strategy::Just(b'B'),
            4 => proptest::strategy::Just(b's'),
            2 => proptest::strategy::Just(b'a'),
            2 => proptest::strategy::Just(b'd'),
            2 => proptest::strategy::Just(b'Z'),
            2 => proptest::strategy::Just(0u8),
            8 => proptest::prelude::any::<u8>(),
        ]
    }

    proptest::proptest! {
        /// Aux bytes reach this walker straight off disk, so arbitrary input has
        /// to fail closed rather than panic. Anything it does hand back must
        /// have come from inside the block: an array of `n` `i16`s cannot fit in
        /// fewer than `2n` bytes, so a longer result means it read past the end.
        #[test]
        fn test_aux_walker_fails_closed_on_arbitrary_bytes(
            aux in proptest::collection::vec(aux_ish_byte(), 0..64),
            query in proptest::array::uniform2(proptest::prelude::any::<u8>()),
        ) {
            for tag in [query, *fgumi_tag::SamTag::AD_BASES] {
                if let Some(found) = find_i16_array_tag_in_aux(&aux, tag) {
                    proptest::prop_assert!(found.len() * 2 <= aux.len());
                }
            }
        }

        /// A planted array must be reachable behind any run of well-formed tags
        /// -- each has to be skipped at exactly the right width -- and every
        /// truncation of that block must yield the planted array or nothing.
        #[test]
        fn test_aux_walker_finds_planted_array_behind_arbitrary_tags(
            prefix in proptest::collection::vec(
                (
                    (b'A'..=b'Z', b'0'..=b'9'),
                    0..TAG_KIND_COUNT,
                    proptest::array::uniform8(0x01u8..=0x39),
                ),
                0..6,
            ),
            // Values are built from two low bytes rather than arbitrary `i16`s
            // so the payload cannot spell an `ad` tag name either.
            value_bytes in proptest::collection::vec((0x01u8..=0x39, 0x01u8..=0x39), 0..8),
        ) {
            let values: Vec<i16> =
                value_bytes.iter().map(|&(lo, hi)| i16::from_le_bytes([lo, hi])).collect();

            let mut aux = Vec::new();
            for &((first, second), kind, filler) in &prefix {
                aux.extend_from_slice(&encode_tag([first, second], kind, filler));
            }
            push_i16_array(&mut aux, *fgumi_tag::SamTag::AD_BASES, &values);

            proptest::prop_assert_eq!(
                find_i16_array_tag_in_aux(&aux, *fgumi_tag::SamTag::AD_BASES),
                Some(values.clone()),
            );

            for cut in 0..aux.len() {
                let found = find_i16_array_tag_in_aux(&aux[..cut], *fgumi_tag::SamTag::AD_BASES);
                proptest::prop_assert!(
                    found.is_none() || found.as_ref() == Some(&values),
                    "truncating to {cut} bytes produced {found:?}, not the planted {values:?}",
                );
            }
        }
    }
}
