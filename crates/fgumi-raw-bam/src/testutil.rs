/// Parsed BAM record for test assertions.
///
/// Provides convenient access to fields extracted from raw BAM bytes
/// produced by `UnmappedBamRecordBuilder`. Used only in tests.
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
    pub fn get_string_tag(&self, tag: &[u8; 2]) -> Option<Vec<u8>> {
        find_z_tag_in_aux(&self.aux_data, *tag)
    }

    /// Find an integer tag value (c/s/i type) by tag name.
    #[must_use]
    pub fn get_int_tag(&self, tag: &[u8; 2]) -> Option<i32> {
        find_int_tag_in_aux(&self.aux_data, *tag)
    }

    /// Find a float tag value by tag name.
    #[must_use]
    pub fn get_float_tag(&self, tag: &[u8; 2]) -> Option<f32> {
        find_float_tag_in_aux(&self.aux_data, *tag)
    }

    /// Find a B:s (i16 array) tag value by tag name.
    #[must_use]
    pub fn get_i16_array_tag(&self, tag: &[u8; 2]) -> Option<Vec<i16>> {
        find_i16_array_tag_in_aux(&self.aux_data, *tag)
    }
}

fn unpack_sequence_for_test(packed: &[u8], l_seq: usize) -> Vec<u8> {
    const DECODE: [u8; 16] = *b"=ACMGRSVTWYHKDBN";
    let mut bases = Vec::with_capacity(l_seq);
    for i in 0..l_seq {
        let byte = packed[i / 2];
        let code = if i % 2 == 0 { byte >> 4 } else { byte & 0x0F };
        bases.push(DECODE[code as usize]);
    }
    bases
}

fn find_z_tag_in_aux(aux: &[u8], tag: [u8; 2]) -> Option<Vec<u8>> {
    let mut i = 0;
    while i + 3 <= aux.len() {
        let t = [aux[i], aux[i + 1]];
        let typ = aux[i + 2];
        i += 3;
        match typ {
            b'Z' => {
                let start = i;
                while i < aux.len() && aux[i] != 0 {
                    i += 1;
                }
                if t == tag {
                    return Some(aux[start..i].to_vec());
                }
                i += 1; // skip NUL
            }
            b'c' => {
                i += 1;
            }
            b's' => {
                i += 2;
            }
            b'i' | b'f' => {
                i += 4;
            }
            b'B' => {
                let sub = aux[i];
                i += 1;
                let count =
                    u32::from_le_bytes([aux[i], aux[i + 1], aux[i + 2], aux[i + 3]]) as usize;
                i += 4;
                let elem_size = match sub {
                    b's' | b'S' => 2,
                    b'i' | b'I' | b'f' => 4,
                    _ => 1,
                };
                i += count * elem_size;
            }
            _ => break,
        }
    }
    None
}

fn find_int_tag_in_aux(aux: &[u8], tag: [u8; 2]) -> Option<i32> {
    let mut i = 0;
    while i + 3 <= aux.len() {
        let t = [aux[i], aux[i + 1]];
        let typ = aux[i + 2];
        i += 3;
        match typ {
            b'c' => {
                let v = i32::from(aux[i].cast_signed());
                if t == tag {
                    return Some(v);
                }
                i += 1;
            }
            b's' => {
                let v = i32::from(i16::from_le_bytes([aux[i], aux[i + 1]]));
                if t == tag {
                    return Some(v);
                }
                i += 2;
            }
            b'i' => {
                let v = i32::from_le_bytes([aux[i], aux[i + 1], aux[i + 2], aux[i + 3]]);
                if t == tag {
                    return Some(v);
                }
                i += 4;
            }
            b'f' => {
                i += 4;
            }
            b'Z' => {
                while i < aux.len() && aux[i] != 0 {
                    i += 1;
                }
                i += 1;
            }
            b'B' => {
                let sub = aux[i];
                i += 1;
                let count =
                    u32::from_le_bytes([aux[i], aux[i + 1], aux[i + 2], aux[i + 3]]) as usize;
                i += 4;
                let elem_size = match sub {
                    b's' | b'S' => 2,
                    b'i' | b'I' | b'f' => 4,
                    _ => 1,
                };
                i += count * elem_size;
            }
            _ => break,
        }
    }
    None
}

fn find_float_tag_in_aux(aux: &[u8], tag: [u8; 2]) -> Option<f32> {
    let mut i = 0;
    while i + 3 <= aux.len() {
        let t = [aux[i], aux[i + 1]];
        let typ = aux[i + 2];
        i += 3;
        match typ {
            b'f' => {
                let v = f32::from_le_bytes([aux[i], aux[i + 1], aux[i + 2], aux[i + 3]]);
                if t == tag {
                    return Some(v);
                }
                i += 4;
            }
            b'c' => {
                i += 1;
            }
            b's' => {
                i += 2;
            }
            b'i' => {
                i += 4;
            }
            b'Z' => {
                while i < aux.len() && aux[i] != 0 {
                    i += 1;
                }
                i += 1;
            }
            b'B' => {
                let sub = aux[i];
                i += 1;
                let count =
                    u32::from_le_bytes([aux[i], aux[i + 1], aux[i + 2], aux[i + 3]]) as usize;
                i += 4;
                let elem_size = match sub {
                    b's' | b'S' => 2,
                    b'i' | b'I' | b'f' => 4,
                    _ => 1,
                };
                i += count * elem_size;
            }
            _ => break,
        }
    }
    None
}

fn find_i16_array_tag_in_aux(aux: &[u8], tag: [u8; 2]) -> Option<Vec<i16>> {
    let mut i = 0;
    while i + 3 <= aux.len() {
        let t = [aux[i], aux[i + 1]];
        let typ = aux[i + 2];
        i += 3;
        match typ {
            b'B' => {
                let sub = aux[i];
                i += 1;
                let count =
                    u32::from_le_bytes([aux[i], aux[i + 1], aux[i + 2], aux[i + 3]]) as usize;
                i += 4;
                if t == tag && sub == b's' {
                    let mut vals = Vec::with_capacity(count);
                    for _ in 0..count {
                        vals.push(i16::from_le_bytes([aux[i], aux[i + 1]]));
                        i += 2;
                    }
                    return Some(vals);
                }
                let elem_size = match sub {
                    b's' | b'S' => 2,
                    b'i' | b'I' | b'f' => 4,
                    _ => 1,
                };
                i += count * elem_size;
            }
            b'c' => {
                i += 1;
            }
            b's' => {
                i += 2;
            }
            b'i' | b'f' => {
                i += 4;
            }
            b'Z' => {
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
    let l_read_name = u8::try_from(name.len() + 1).unwrap(); // +1 for null terminator
    let n_cigar_op = u16::try_from(cigar_ops.len()).unwrap();
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
    buf[16..20].copy_from_slice(&u32::try_from(seq_len).unwrap().to_le_bytes());
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
    (u32::try_from(len).unwrap() << 4) | op_type
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
    make_b_array_tag(tag, b'i', u32::try_from(values.len()).unwrap(), &bytes)
}

/// # Panics
///
/// Panics if `values` length exceeds `u32::MAX`.
#[must_use]
pub fn make_b_float_array_tag(tag: [u8; 2], values: &[f32]) -> Vec<u8> {
    let bytes: Vec<u8> = values.iter().flat_map(|v| v.to_le_bytes()).collect();
    make_b_array_tag(tag, b'f', u32::try_from(values.len()).unwrap(), &bytes)
}

/// # Panics
///
/// Panics if `values` length exceeds `u32::MAX`.
#[must_use]
pub fn make_b_uint8_array_tag(tag: [u8; 2], values: &[u8]) -> Vec<u8> {
    make_b_array_tag(tag, b'C', u32::try_from(values.len()).unwrap(), values)
}

/// # Panics
///
/// Panics if `values` length exceeds `u32::MAX`.
#[must_use]
pub fn make_b_int8_array_tag(tag: [u8; 2], values: &[i8]) -> Vec<u8> {
    let bytes: Vec<u8> = values.iter().map(|&v| v.cast_unsigned()).collect();
    make_b_array_tag(tag, b'c', u32::try_from(values.len()).unwrap(), &bytes)
}

/// # Panics
///
/// Panics if `values` length exceeds `u32::MAX`.
#[must_use]
pub fn make_b_int16_array_tag(tag: [u8; 2], values: &[i16]) -> Vec<u8> {
    let bytes: Vec<u8> = values.iter().flat_map(|v| v.to_le_bytes()).collect();
    make_b_array_tag(tag, b's', u32::try_from(values.len()).unwrap(), &bytes)
}

/// # Panics
///
/// Panics if `values` length exceeds `u32::MAX`.
#[must_use]
pub fn make_b_uint16_array_tag(tag: [u8; 2], values: &[u16]) -> Vec<u8> {
    let bytes: Vec<u8> = values.iter().flat_map(|v| v.to_le_bytes()).collect();
    make_b_array_tag(tag, b'S', u32::try_from(values.len()).unwrap(), &bytes)
}

/// # Panics
///
/// Panics if `values` length exceeds `u32::MAX`.
#[must_use]
pub fn make_b_uint32_array_tag(tag: [u8; 2], values: &[u32]) -> Vec<u8> {
    let bytes: Vec<u8> = values.iter().flat_map(|v| v.to_le_bytes()).collect();
    make_b_array_tag(tag, b'I', u32::try_from(values.len()).unwrap(), &bytes)
}
