use crate::fields::{TAG_FIXED_SIZES, aux_data_offset, aux_data_offset_from_record, aux_data_slice, tag_value_size};

/// Find a string (Z-type) tag in auxiliary data, returning value bytes without null terminator.
#[must_use]
pub fn find_string_tag<'a>(aux_data: &'a [u8], tag: &[u8; 2]) -> Option<&'a [u8]> {
    let mut p = 0;
    while p + 3 <= aux_data.len() {
        let t = &aux_data[p..p + 2];
        let val_type = aux_data[p + 2];

        if t == tag {
            return match val_type {
                b'Z' => {
                    let start = p + 3;
                    let end = aux_data[start..].iter().position(|&b| b == 0)?;
                    Some(&aux_data[start..start + end])
                }
                _ => None,
            };
        }

        if let Some(size) = tag_value_size(val_type, &aux_data[p + 3..]) {
            p += 3 + size;
        } else {
            break;
        }
    }
    None
}

/// Check whether a tag exists in auxiliary data, returning its type byte if found.
///
/// Returns `Some(type_byte)` (e.g. `b'Z'`, `b'C'`, `b'i'`) if the tag is present,
/// `None` if the tag is absent.
#[must_use]
pub fn find_tag_type(aux_data: &[u8], tag: &[u8; 2]) -> Option<u8> {
    let mut p = 0;
    while p + 3 <= aux_data.len() {
        let t = &aux_data[p..p + 2];
        let val_type = aux_data[p + 2];

        if t == tag {
            return Some(val_type);
        }

        if let Some(size) = tag_value_size(val_type, &aux_data[p + 3..]) {
            p += 3 + size;
        } else {
            break;
        }
    }
    None
}

/// Find a string tag in a complete BAM record.
#[must_use]
pub fn find_string_tag_in_record<'a>(bam: &'a [u8], tag: &[u8; 2]) -> Option<&'a [u8]> {
    let aux = aux_data_slice(bam);
    if aux.is_empty() {
        return None;
    }
    find_string_tag(aux, tag)
}

/// Find the byte range `[start, end)` of an entire tag entry (tag+type+value) in aux data.
///
/// Returns offsets relative to the start of `aux_data`.
/// Returns `None` if the tag is not found.
#[must_use]
pub fn find_tag_bounds(aux_data: &[u8], tag: &[u8; 2]) -> Option<(usize, usize)> {
    let mut p = 0;
    while p + 3 <= aux_data.len() {
        let t = &aux_data[p..p + 2];
        let val_type = aux_data[p + 2];

        if let Some(size) = tag_value_size(val_type, &aux_data[p + 3..]) {
            let entry_end = p + 3 + size;
            if t == tag {
                return Some((p, entry_end));
            }
            p = entry_end;
        } else {
            break;
        }
    }
    None
}

/// Find a uint8 (C-type) tag value in auxiliary data.
#[must_use]
pub fn find_uint8_tag(aux_data: &[u8], tag: &[u8; 2]) -> Option<u8> {
    let mut p = 0;
    while p + 3 <= aux_data.len() {
        let t = &aux_data[p..p + 2];
        let val_type = aux_data[p + 2];

        if t == tag && val_type == b'C' && p + 4 <= aux_data.len() {
            return Some(aux_data[p + 3]);
        }

        if let Some(size) = tag_value_size(val_type, &aux_data[p + 3..]) {
            p += 3 + size;
        } else {
            break;
        }
    }
    None
}

/// Find a float (f-type) tag value in auxiliary data.
#[must_use]
pub fn find_float_tag(aux_data: &[u8], tag: &[u8; 2]) -> Option<f32> {
    let mut p = 0;
    while p + 3 <= aux_data.len() {
        let t = &aux_data[p..p + 2];
        let val_type = aux_data[p + 2];

        if t == tag && val_type == b'f' && p + 7 <= aux_data.len() {
            return Some(f32::from_le_bytes([
                aux_data[p + 3],
                aux_data[p + 4],
                aux_data[p + 5],
                aux_data[p + 6],
            ]));
        }

        if let Some(size) = tag_value_size(val_type, &aux_data[p + 3..]) {
            p += 3 + size;
        } else {
            break;
        }
    }
    None
}

/// Find an integer tag value in auxiliary data.
///
/// Supports signed/unsigned byte, short, and int types (c/C/s/S/i/I).
#[must_use]
pub fn find_int_tag(aux_data: &[u8], tag: &[u8; 2]) -> Option<i64> {
    let mut p = 0;
    while p + 3 <= aux_data.len() {
        let t = &aux_data[p..p + 2];
        let val_type = aux_data[p + 2];

        if t == tag {
            return match val_type {
                b'c' if p + 4 <= aux_data.len() => Some(i64::from(aux_data[p + 3].cast_signed())),
                b'C' if p + 4 <= aux_data.len() => Some(i64::from(aux_data[p + 3])),
                b's' if p + 5 <= aux_data.len() => {
                    Some(i64::from(i16::from_le_bytes([aux_data[p + 3], aux_data[p + 4]])))
                }
                b'S' if p + 5 <= aux_data.len() => {
                    Some(i64::from(u16::from_le_bytes([aux_data[p + 3], aux_data[p + 4]])))
                }
                b'i' if p + 7 <= aux_data.len() => Some(i64::from(i32::from_le_bytes([
                    aux_data[p + 3],
                    aux_data[p + 4],
                    aux_data[p + 5],
                    aux_data[p + 6],
                ]))),
                b'I' if p + 7 <= aux_data.len() => Some(i64::from(u32::from_le_bytes([
                    aux_data[p + 3],
                    aux_data[p + 4],
                    aux_data[p + 5],
                    aux_data[p + 6],
                ]))),
                _ => None,
            };
        }

        if let Some(size) = tag_value_size(val_type, &aux_data[p + 3..]) {
            p += 3 + size;
        } else {
            break;
        }
    }
    None
}

/// Find MI (Molecular Identifier) tag in auxiliary data.
///
/// Returns the value as `(integer_value, is_A_suffix)`.
/// - For string values like "12345" or "12345/A" or "12345/B"
/// - For integer values, returns `(value, true)`
/// - Returns `None` if MI tag not found.
#[must_use]
pub fn find_mi_tag(aux_data: &[u8]) -> Option<(u64, bool)> {
    let mut pos = 0;
    while pos + 3 <= aux_data.len() {
        let tag = &aux_data[pos..pos + 2];
        let val_type = aux_data[pos + 2];

        if tag == b"MI" {
            return match val_type {
                // String type - parse "12345" or "12345/A" or "12345/B"
                b'Z' => {
                    let start = pos + 3;
                    let end = aux_data[start..].iter().position(|&b| b == 0)?;
                    let s = &aux_data[start..start + end];
                    parse_mi_bytes(s)
                }
                // Integer types (with per-type bounds checks)
                // Signed types: return None for negative values (MI must be non-negative)
                b'c' if pos + 4 <= aux_data.len() => {
                    let v = i8::from_le_bytes([aux_data[pos + 3]]);
                    if v >= 0 { Some((u64::from(v.cast_unsigned()), true)) } else { None }
                }
                b'C' if pos + 4 <= aux_data.len() => Some((u64::from(aux_data[pos + 3]), true)),
                b's' if pos + 5 <= aux_data.len() => {
                    let v = i16::from_le_bytes([aux_data[pos + 3], aux_data[pos + 4]]);
                    if v >= 0 { Some((u64::from(v.cast_unsigned()), true)) } else { None }
                }
                b'S' if pos + 5 <= aux_data.len() => Some((
                    u64::from(u16::from_le_bytes([aux_data[pos + 3], aux_data[pos + 4]])),
                    true,
                )),
                b'i' if pos + 7 <= aux_data.len() => {
                    let v = i32::from_le_bytes([
                        aux_data[pos + 3],
                        aux_data[pos + 4],
                        aux_data[pos + 5],
                        aux_data[pos + 6],
                    ]);
                    if v >= 0 { Some((u64::from(v.cast_unsigned()), true)) } else { None }
                }
                b'I' if pos + 7 <= aux_data.len() => Some((
                    u64::from(u32::from_le_bytes([
                        aux_data[pos + 3],
                        aux_data[pos + 4],
                        aux_data[pos + 5],
                        aux_data[pos + 6],
                    ])),
                    true,
                )),
                _ => None,
            };
        }

        if let Some(size) = tag_value_size(val_type, &aux_data[pos + 3..]) {
            pos += 3 + size;
        } else {
            break;
        }
    }
    None
}

/// Parse MI tag bytes to `(integer, is_A_suffix)`.
///
/// Format: "12345" or "12345/A" or "12345/B"
/// Returns `None` if the string contains non-digit characters in the numeric portion.
#[inline]
fn parse_mi_bytes(s: &[u8]) -> Option<(u64, bool)> {
    if s.is_empty() {
        return None;
    }

    // Check for /A or /B suffix
    let (num_part, is_a) = if s.len() >= 2 && s[s.len() - 2] == b'/' {
        let suffix = s[s.len() - 1];
        (&s[..s.len() - 2], suffix != b'B')
    } else {
        (s, true)
    };

    if num_part.is_empty() {
        return None;
    }

    // Parse the numeric part
    let mut value = 0u64;
    for &c in num_part {
        if c.is_ascii_digit() {
            value = value.saturating_mul(10).saturating_add(u64::from(c - b'0'));
        } else {
            return None; // Invalid format
        }
    }

    Some((value, is_a))
}

/// Find MI tag in a complete BAM record.
///
/// Returns `(integer_value, is_A_suffix)` or `(0, true)` if not found.
#[must_use]
pub fn find_mi_tag_in_record(bam: &[u8]) -> (u64, bool) {
    let l_read_name = bam[8] as usize;
    let n_cigar_op = u16::from_le_bytes([bam[12], bam[13]]) as usize;
    let l_seq = u32::from_le_bytes([bam[16], bam[17], bam[18], bam[19]]) as usize;

    let aux_start = aux_data_offset(l_read_name, n_cigar_op, l_seq);

    if aux_start >= bam.len() {
        return (0, true);
    }

    find_mi_tag(&bam[aux_start..]).unwrap_or((0, true))
}

/// Find the MC (mate CIGAR) tag in auxiliary data.
///
/// Returns the CIGAR string, or None if not found.
#[must_use]
pub fn find_mc_tag(aux_data: &[u8]) -> Option<&str> {
    let mut pos = 0;
    while pos + 3 <= aux_data.len() {
        let tag = &aux_data[pos..pos + 2];
        let val_type = aux_data[pos + 2];

        if tag == b"MC" {
            return match val_type {
                b'Z' => {
                    let start = pos + 3;
                    let end = aux_data[start..].iter().position(|&b| b == 0)?;
                    std::str::from_utf8(&aux_data[start..start + end]).ok()
                }
                _ => None,
            };
        }

        if let Some(size) = tag_value_size(val_type, &aux_data[pos + 3..]) {
            pos += 3 + size;
        } else {
            break;
        }
    }
    None
}

/// Find MC tag in a complete BAM record.
#[must_use]
pub fn find_mc_tag_in_record(bam: &[u8]) -> Option<&str> {
    let l_read_name = bam[8] as usize;
    let n_cigar_op = u16::from_le_bytes([bam[12], bam[13]]) as usize;
    let l_seq = u32::from_le_bytes([bam[16], bam[17], bam[18], bam[19]]) as usize;

    let aux_start = aux_data_offset(l_read_name, n_cigar_op, l_seq);

    if aux_start >= bam.len() {
        return None;
    }

    find_mc_tag(&bam[aux_start..])
}

/// Result of a single-pass extraction of multiple string tags from aux data.
/// Used by `compute_group_key_from_raw` to avoid scanning aux data 3 times.
pub struct AuxStringTags<'a> {
    pub rg: Option<&'a [u8]>,
    pub cell: Option<&'a [u8]>,
    pub mc: Option<&'a str>,
}

/// Extract RG, cell barcode, and MC tags in a single pass over the aux data.
#[must_use]
pub fn extract_aux_string_tags<'a>(aux_data: &'a [u8], cell_tag: &[u8; 2]) -> AuxStringTags<'a> {
    let mut result = AuxStringTags { rg: None, cell: None, mc: None };
    let mut found = 0u8; // bit 0=RG, bit 1=cell, bit 2=MC
    let mut p = 0;
    while p + 3 <= aux_data.len() {
        let t = [aux_data[p], aux_data[p + 1]];
        let val_type = aux_data[p + 2];

        if val_type == b'Z' {
            let start = p + 3;
            if let Some(end) = aux_data[start..].iter().position(|&b| b == 0) {
                let value = &aux_data[start..start + end];
                if t == *b"RG" {
                    result.rg = Some(value);
                    found |= 1;
                } else if t == *cell_tag {
                    result.cell = Some(value);
                    found |= 2;
                } else if t == *b"MC" {
                    result.mc = std::str::from_utf8(value).ok();
                    found |= 4;
                }
                if found == 7 {
                    return result;
                }
                p = start + end + 1;
            } else {
                break;
            }
            continue;
        }

        if let Some(size) = tag_value_size(val_type, &aux_data[p + 3..]) {
            p += 3 + size;
        } else {
            break;
        }
    }
    result
}

/// Zero-allocation reference to a B-type array tag in aux data.
pub struct ArrayTagRef<'a> {
    /// Element bytes (raw, little-endian).
    pub data: &'a [u8],
    /// Sub-type byte (`b'c'`, `b'C'`, `b's'`, `b'S'`, `b'i'`, `b'I'`, `b'f'`).
    pub elem_type: u8,
    /// Number of elements.
    pub count: usize,
    /// Size of each element in bytes (1, 2, or 4).
    pub elem_size: usize,
}

/// Find a B-type (array) tag in auxiliary data, returning a zero-allocation reference.
///
/// Scans aux linearly (same pattern as `find_string_tag`). Returns `None` if
/// the tag is absent or is not of type `B`.
#[must_use]
pub fn find_array_tag<'a>(aux_data: &'a [u8], tag: &[u8; 2]) -> Option<ArrayTagRef<'a>> {
    let mut p = 0;
    while p + 3 <= aux_data.len() {
        let t = &aux_data[p..p + 2];
        let val_type = aux_data[p + 2];

        if t == tag && val_type == b'B' {
            let data_start = p + 3;
            if data_start + 5 > aux_data.len() {
                return None;
            }
            let elem_type = aux_data[data_start];
            let count = u32::from_le_bytes([
                aux_data[data_start + 1],
                aux_data[data_start + 2],
                aux_data[data_start + 3],
                aux_data[data_start + 4],
            ]) as usize;
            let elem_size = TAG_FIXED_SIZES[elem_type as usize] as usize;
            if elem_size == 0 {
                return None;
            }
            let elements_start = data_start + 5;
            let total_bytes = count.checked_mul(elem_size)?;
            let elements_end = elements_start.checked_add(total_bytes)?;
            if elements_end > aux_data.len() {
                return None;
            }
            return Some(ArrayTagRef {
                data: &aux_data[elements_start..elements_end],
                elem_type,
                count,
                elem_size,
            });
        }

        if let Some(size) = tag_value_size(val_type, &aux_data[p + 3..]) {
            p += 3 + size;
        } else {
            break;
        }
    }
    None
}

/// Read one element from an `ArrayTagRef` as `u16`.
///
/// Handles `C` (u8), `S` (u16), `s` (i16, clamped to 0), `c` (i8, clamped to 0) sub-types.
/// Other sub-types return 0.
#[inline]
#[must_use]
pub fn array_tag_element_u16(tag_ref: &ArrayTagRef, index: usize) -> u16 {
    if index >= tag_ref.count {
        return 0;
    }
    let off = index * tag_ref.elem_size;
    match tag_ref.elem_type {
        b'C' => u16::from(tag_ref.data[off]),
        b'S' => u16::from_le_bytes([tag_ref.data[off], tag_ref.data[off + 1]]),
        b's' => {
            let v = i16::from_le_bytes([tag_ref.data[off], tag_ref.data[off + 1]]);
            v.max(0).cast_unsigned()
        }
        b'c' => {
            let v = tag_ref.data[off].cast_signed();
            u16::from(v.max(0).cast_unsigned())
        }
        _ => 0,
    }
}

/// Read all elements from an `ArrayTagRef` as a `Vec<u16>`.
///
/// This is useful when you need to release the immutable borrow on the record
/// before mutating it (e.g., masking bases while reading per-base depth/error tags).
#[must_use]
pub fn array_tag_to_vec_u16(tag_ref: &ArrayTagRef) -> Vec<u16> {
    (0..tag_ref.count).map(|i| array_tag_element_u16(tag_ref, i)).collect()
}

/// Append a string (Z-type) tag to a BAM record.
///
/// The tag is appended at the end of the record: `[tag_byte_1, tag_byte_2, 'Z', value..., NUL]`.
pub fn append_string_tag(record: &mut Vec<u8>, tag: &[u8; 2], value: &[u8]) {
    record.push(tag[0]);
    record.push(tag[1]);
    record.push(b'Z');
    record.extend_from_slice(value);
    record.push(0); // null terminator
}

/// Append an integer tag using the smallest signed type that fits.
///
/// Encodes as:
/// - `i8` (type `'c'`): if value in `[-128, 127]`
/// - `i16` (type `'s'`): if value in `[-32768, 32767]`
/// - `i32` (type `'i'`): otherwise
///
/// This matches the behavior of [`crate::sam::to_smallest_signed_int`].
pub fn append_int_tag(record: &mut Vec<u8>, tag: &[u8; 2], value: i32) {
    record.push(tag[0]);
    record.push(tag[1]);
    if let Ok(v) = i8::try_from(value) {
        record.push(b'c');
        record.push(v.cast_unsigned());
    } else if let Ok(v) = i16::try_from(value) {
        record.push(b's');
        record.extend_from_slice(&v.to_le_bytes());
    } else {
        record.push(b'i');
        record.extend_from_slice(&value.to_le_bytes());
    }
}

/// Append a float (`f`-type) tag to a BAM record.
pub fn append_float_tag(record: &mut Vec<u8>, tag: &[u8; 2], value: f32) {
    record.push(tag[0]);
    record.push(tag[1]);
    record.push(b'f');
    record.extend_from_slice(&value.to_le_bytes());
}

/// Append an `i16` array (`B:s`-type) tag to a BAM record.
///
/// Format: `[tag0, tag1, 'B', 's', count_u32_le, values_i16_le...]`
///
/// # Panics
///
/// Panics if `values.len()` exceeds `u32::MAX`.
pub fn append_i16_array_tag(record: &mut Vec<u8>, tag: &[u8; 2], values: &[i16]) {
    record.push(tag[0]);
    record.push(tag[1]);
    record.push(b'B');
    record.push(b's');
    record.extend_from_slice(&u32::try_from(values.len()).expect("array length exceeds u32").to_le_bytes());
    for &v in values {
        record.extend_from_slice(&v.to_le_bytes());
    }
}

/// Append a Phred+33 encoded quality string (`Z`-type) tag.
///
/// Converts raw Phred scores (0-93) to ASCII (Phred+33) and writes
/// directly as a null-terminated string tag. Avoids intermediate String allocation.
pub fn append_phred33_string_tag(record: &mut Vec<u8>, tag: &[u8; 2], quals: &[u8]) {
    record.push(tag[0]);
    record.push(tag[1]);
    record.push(b'Z');
    for &q in quals {
        debug_assert!(q <= 93, "Phred score out of range: {q}");
        record.push(q.saturating_add(33));
    }
    record.push(0); // null terminator
}

/// Remove a tag from a BAM record. No-op if the tag is not found.
pub fn remove_tag(record: &mut Vec<u8>, tag: &[u8; 2]) {
    let Some(aux_start) = aux_data_offset_from_record(record) else {
        return;
    };
    if aux_start >= record.len() {
        return;
    }
    if let Some((start, end)) = find_tag_bounds(&record[aux_start..], tag) {
        let abs_start = aux_start + start;
        let abs_end = aux_start + end;
        record.drain(abs_start..abs_end);
    }
}

/// Update the value of a string (Z-type) tag in a BAM record.
///
/// If the tag exists, its value is replaced with `new_value`.
/// If the tag does not exist, it is appended.
///
/// **Note:** `find_tag_bounds` is type-agnostic, so if the tag exists with a
/// non-Z type (e.g., integer), this will replace it with a Z-type string tag.
/// This is intentional for MI/OX tags which should always be strings.
pub fn update_string_tag(record: &mut Vec<u8>, tag: &[u8; 2], new_value: &[u8]) {
    let aux_start = aux_data_offset_from_record(record).unwrap_or(record.len());
    if aux_start < record.len() {
        if let Some((start, end)) = find_tag_bounds(&record[aux_start..], tag) {
            let abs_start = aux_start + start;
            let abs_end = aux_start + end;
            let old_value_len = end - start - 4; // subtract tag(2) + type(1) + NUL(1)
            if old_value_len == new_value.len() {
                // Same length: overwrite value bytes in-place (no memmove)
                let value_start = abs_start + 3; // skip tag(2) + type(1)
                record[value_start..value_start + new_value.len()].copy_from_slice(new_value);
            } else {
                // Different length: splice replacement
                let mut replacement = Vec::with_capacity(3 + new_value.len() + 1);
                replacement.push(tag[0]);
                replacement.push(tag[1]);
                replacement.push(b'Z');
                replacement.extend_from_slice(new_value);
                replacement.push(0);
                record.splice(abs_start..abs_end, replacement);
            }
            return;
        }
    }
    // Tag not found — append
    append_string_tag(record, tag, new_value);
}

/// Update an existing int tag in-place, or append if absent.
///
/// If the existing tag is a 4-byte `i`/`I` type, the value is overwritten in-place.
/// Otherwise the tag is removed and re-appended using `append_int_tag`
/// (smallest signed type that fits).
pub fn update_int_tag(record: &mut Vec<u8>, tag: &[u8; 2], value: i32) {
    let aux_start = aux_data_offset_from_record(record).unwrap_or(record.len());
    if aux_start < record.len() {
        if let Some((start, end)) = find_tag_bounds(&record[aux_start..], tag) {
            let abs_start = aux_start + start;
            let abs_end = aux_start + end;
            let val_type = record[abs_start + 2];
            // If 4-byte integer type, overwrite in-place
            if matches!(val_type, b'i' | b'I') && (abs_end - abs_start) == 7 {
                record[abs_start + 3..abs_start + 7].copy_from_slice(&value.to_le_bytes());
                return;
            }
            // Different size — remove and re-append
            record.drain(abs_start..abs_end);
            append_int_tag(record, tag, value);
            return;
        }
    }
    // Tag not found — append
    append_int_tag(record, tag, value);
}

/// Reverse elements of a B-type array tag in place. No-op if tag not found.
pub fn reverse_array_tag_in_place(record: &mut [u8], aux_offset: usize, tag: &[u8; 2]) {
    if aux_offset >= record.len() {
        return;
    }
    let aux_data = &record[aux_offset..];
    // First find the tag to get its offset and element info
    let mut p = 0;
    while p + 3 <= aux_data.len() {
        let t = &aux_data[p..p + 2];
        let val_type = aux_data[p + 2];

        if t == tag && val_type == b'B' {
            let data_start = p + 3;
            if data_start + 5 > aux_data.len() {
                return;
            }
            let elem_type = aux_data[data_start];
            let count = u32::from_le_bytes([
                aux_data[data_start + 1],
                aux_data[data_start + 2],
                aux_data[data_start + 3],
                aux_data[data_start + 4],
            ]) as usize;
            let elem_size = TAG_FIXED_SIZES[elem_type as usize] as usize;
            if elem_size == 0 || count == 0 {
                return;
            }
            let Some(total_bytes) = count.checked_mul(elem_size) else {
                return;
            };
            let elements_start = aux_offset + data_start + 5;
            let Some(elements_end) = elements_start.checked_add(total_bytes) else {
                return;
            };
            if elements_end > record.len() {
                return;
            }
            // Reverse elements by swapping elem_size-byte chunks
            let mut i = 0;
            let mut j = count - 1;
            while i < j {
                let off_i = elements_start + i * elem_size;
                let off_j = elements_start + j * elem_size;
                for k in 0..elem_size {
                    record.swap(off_i + k, off_j + k);
                }
                i += 1;
                j -= 1;
            }
            return;
        }

        if let Some(size) = tag_value_size(val_type, &aux_data[p + 3..]) {
            p += 3 + size;
        } else {
            break;
        }
    }
}

/// Reverse bytes of a Z-type string tag value in place. No-op if tag not found.
pub fn reverse_string_tag_in_place(record: &mut [u8], aux_offset: usize, tag: &[u8; 2]) {
    if aux_offset >= record.len() {
        return;
    }
    let aux_data = &record[aux_offset..];
    let mut p = 0;
    while p + 3 <= aux_data.len() {
        let t = &aux_data[p..p + 2];
        let val_type = aux_data[p + 2];

        if t == tag && val_type == b'Z' {
            let start = aux_offset + p + 3;
            if let Some(nul_off) = record[start..].iter().position(|&b| b == 0) {
                let end = start + nul_off;
                if end > start {
                    record[start..end].reverse();
                }
            }
            return;
        }

        if let Some(size) = tag_value_size(val_type, &aux_data[p + 3..]) {
            p += 3 + size;
        } else {
            break;
        }
    }
}

/// Reverse-complement a Z-type string tag value in place (A<->T, C<->G).
pub fn reverse_complement_string_tag_in_place(record: &mut [u8], aux_offset: usize, tag: &[u8; 2]) {
    if aux_offset >= record.len() {
        return;
    }
    let aux_data = &record[aux_offset..];
    let mut p = 0;
    while p + 3 <= aux_data.len() {
        let t = &aux_data[p..p + 2];
        let val_type = aux_data[p + 2];

        if t == tag && val_type == b'Z' {
            let start = aux_offset + p + 3;
            if let Some(nul_off) = record[start..].iter().position(|&b| b == 0) {
                let end = start + nul_off;
                if end > start {
                    // Reverse
                    record[start..end].reverse();
                    // Complement each base
                    for b in &mut record[start..end] {
                        *b = match *b {
                            b'A' | b'a' => b'T',
                            b'T' | b't' => b'A',
                            b'C' | b'c' => b'G',
                            b'G' | b'g' => b'C',
                            other => other,
                        };
                    }
                }
            }
            return;
        }

        if let Some(size) = tag_value_size(val_type, &aux_data[p + 3..]) {
            p += 3 + size;
        } else {
            break;
        }
    }
}

#[cfg(test)]
#[allow(clippy::identity_op)]
mod tests {
    use super::*;
    use crate::testutil::*;
    use rstest::rstest;

    // ========================================================================
    // find_mi_tag tests
    // ========================================================================

    #[test]
    fn test_find_mi_tag() {
        // Test numeric MI tag: MI:Z:12345
        let aux_data = b"MIZ12345\x00";
        let result = find_mi_tag(aux_data);
        assert_eq!(result, Some((12345, true)));

        // Test MI tag with /A suffix: MI:Z:12345/A
        let aux_data = b"MIZ12345/A\x00";
        let result = find_mi_tag(aux_data);
        assert_eq!(result, Some((12345, true)));

        // Test MI tag with /B suffix: MI:Z:12345/B
        let aux_data = b"MIZ12345/B\x00";
        let result = find_mi_tag(aux_data);
        assert_eq!(result, Some((12345, false)));

        // Test integer MI tag: MI:i:42
        let mut aux_data = [b'M', b'I', b'i', 0, 0, 0, 0];
        aux_data[3..7].copy_from_slice(&42i32.to_le_bytes());
        let result = find_mi_tag(&aux_data);
        assert_eq!(result, Some((42, true)));
    }

    #[test]
    fn test_find_mi_tag_empty_aux() {
        let aux_data: &[u8] = &[];
        assert_eq!(find_mi_tag(aux_data), None);
    }

    #[rstest]
    #[case::unsigned_byte(b'C', &[200], Some((200, true)))]
    #[case::unsigned_short(b'S', &5000u16.to_le_bytes(), Some((5000, true)))]
    #[case::unsigned_int(b'I', &100_000u32.to_le_bytes(), Some((100_000, true)))]
    #[case::positive_signed_byte(b'c', &[42], Some((42, true)))]
    #[case::positive_signed_short(b's', &100i16.to_le_bytes(), Some((100, true)))]
    #[case::positive_signed_int(b'i', &12345i32.to_le_bytes(), Some((12345, true)))]
    #[case::negative_signed_byte(b'c', &[(-1i8).cast_unsigned()], None)]
    #[case::negative_signed_short(b's', &(-1i16).to_le_bytes(), None)]
    #[case::negative_signed_int(b'i', &(-1i32).to_le_bytes(), None)]
    #[case::float_type(b'f', &1.0f32.to_le_bytes(), None)]
    fn test_find_mi_tag_by_type(
        #[case] type_byte: u8,
        #[case] value_bytes: &[u8],
        #[case] expected: Option<(u64, bool)>,
    ) {
        let mut aux = vec![b'M', b'I', type_byte];
        aux.extend_from_slice(value_bytes);
        assert_eq!(find_mi_tag(&aux), expected);
    }

    #[test]
    fn test_find_mi_tag_after_other_tags() {
        // Put another tag before MI
        // XY:C:42 then MI:Z:99\0
        let mut aux = Vec::new();
        aux.extend_from_slice(b"XYC"); // tag XY, type C
        aux.push(42); // value
        aux.extend_from_slice(b"MIZ99\x00"); // MI:Z:99
        assert_eq!(find_mi_tag(&aux), Some((99, true)));
    }

    #[test]
    fn test_find_mi_tag_invalid_string() {
        // MI:Z:abc\0 -> non-numeric chars return None
        let aux = b"MIZabc\x00";
        assert_eq!(find_mi_tag(aux), None);
    }

    #[test]
    fn test_find_mi_tag_empty_string_value() {
        // MI:Z: with empty string -> MI:Z:\0
        let aux = b"MIZ\x00";
        // parse_mi_bytes gets an empty slice, returns None
        assert_eq!(find_mi_tag(aux), None);
    }

    #[test]
    fn test_find_mi_tag_suffix_other_than_ab() {
        // "12345/C" -> suffix != 'B', so is_a=true
        let aux = b"MIZ12345/C\x00";
        assert_eq!(find_mi_tag(aux), Some((12345, true)));
    }

    #[test]
    fn test_find_mi_tag_large_number() {
        // Large number: "9999999999"
        let aux = b"MIZ9999999999\x00";
        assert_eq!(find_mi_tag(aux), Some((9_999_999_999, true)));
    }

    #[rstest]
    #[case::truncated_signed_byte(b"MIc" as &[u8])]
    #[case::truncated_signed_short(&[b'M', b'I', b's', 42] as &[u8])]
    #[case::truncated_signed_int(&[b'M', b'I', b'i', 42, 0] as &[u8])]
    fn test_find_mi_tag_truncated_returns_none(#[case] aux: &[u8]) {
        // Truncated tag data returns None instead of panicking
        assert_eq!(find_mi_tag(aux), None);
    }

    #[test]
    fn test_find_mi_tag_just_slash_suffix() {
        // MI:Z:/B -> num_part is empty after stripping suffix, should return None
        let aux = b"MIZ/B\x00";
        assert_eq!(find_mi_tag(aux), None);
    }

    // ========================================================================
    // find_mi_tag_in_record tests
    // ========================================================================

    #[test]
    fn test_find_mi_tag_in_record_no_aux() {
        let rec = make_bam_bytes(0, 0, 0, b"rea", &[], 0, -1, -1, &[]);
        assert_eq!(find_mi_tag_in_record(&rec), (0, true));
    }

    // ========================================================================
    // find_string_tag tests
    // ========================================================================

    #[test]
    fn test_find_string_tag_present() {
        let aux = b"RGZsample1\x00";
        assert_eq!(find_string_tag(aux, b"RG"), Some(b"sample1".as_ref()));
    }

    #[test]
    fn test_find_string_tag_absent() {
        let aux = b"RGZsample1\x00";
        assert_eq!(find_string_tag(aux, b"RX"), None);
    }

    #[test]
    fn test_find_string_tag_after_other_tags() {
        // NM:C:5 then RX:Z:ACGT
        let mut aux = Vec::new();
        aux.extend_from_slice(b"NMC");
        aux.push(5);
        aux.extend_from_slice(b"RXZACGT\x00");
        assert_eq!(find_string_tag(&aux, b"RX"), Some(b"ACGT".as_ref()));
    }

    #[test]
    fn test_find_string_tag_in_record() {
        let aux = b"RXZhello\x00";
        let rec = make_bam_bytes(0, 0, 0, b"rea", &[], 0, -1, -1, aux);
        assert_eq!(find_string_tag_in_record(&rec, b"RX"), Some(b"hello".as_ref()));
    }

    #[test]
    fn test_find_string_tag_non_z_type_returns_none() {
        // Tag matches but type is 'C' not 'Z'
        let aux = [b'R', b'X', b'C', 42];
        assert_eq!(find_string_tag(&aux, b"RX"), None);
    }

    #[test]
    fn test_find_string_tag_in_record_no_aux() {
        // Record with no aux data (aux_start >= len)
        let rec = make_bam_bytes(0, 0, 0, b"rea", &[], 0, -1, -1, &[]);
        assert_eq!(find_string_tag_in_record(&rec, b"RX"), None);
    }

    #[test]
    fn test_find_string_tag_truncated_z_value() {
        // RX:Z:hello but no null terminator — should return None
        let aux = b"RXZhello";
        assert_eq!(find_string_tag(aux, b"RX"), None);
    }

    // --- Per-type find_string_tag tests ---

    // --- B:c array (int8 array) ---

    #[test]
    fn test_find_string_tag_cannot_find_b_int8_array() {
        let aux = make_b_int8_array_tag(*b"Xc", &[-1, 0, 1]);
        assert_eq!(find_string_tag(&aux, b"Xc"), None);
    }

    // --- B:s array (int16 array) ---

    #[test]
    fn test_find_string_tag_cannot_find_b_int16_array() {
        let aux = make_b_int16_array_tag(*b"Xs", &[-100, 0, 200]);
        assert_eq!(find_string_tag(&aux, b"Xs"), None);
    }

    // --- B:S array (uint16 array) ---

    #[test]
    fn test_find_string_tag_cannot_find_b_uint16_array() {
        let aux = make_b_uint16_array_tag(*b"XS", &[100, 200, 300]);
        assert_eq!(find_string_tag(&aux, b"XS"), None);
    }

    // --- B:I array (uint32 array) ---

    #[test]
    fn test_find_string_tag_cannot_find_b_uint32_array() {
        let aux = make_b_uint32_array_tag(*b"XI", &[1000, 2000, 3000]);
        assert_eq!(find_string_tag(&aux, b"XI"), None);
    }

    // --- B:i array (the pa tag type) ---

    #[test]
    fn test_find_string_tag_cannot_find_b_int_array() {
        let aux = make_b_int_array_tag(*b"pa", &[0, 27_056_961, 0, 207, 60005, 1]);
        // BUG: find_string_tag returns None for B:i tags
        assert_eq!(find_string_tag(&aux, b"pa"), None);
    }

    // --- B:f array ---

    #[test]
    fn test_find_string_tag_cannot_find_b_float_array() {
        let aux = make_b_float_array_tag(*b"XF", &[1.0, 2.5, 3.0]);
        assert_eq!(find_string_tag(&aux, b"XF"), None);
    }

    // --- B:C array (uint8 array) ---

    #[test]
    fn test_find_string_tag_cannot_find_b_uint8_array() {
        let aux = make_b_uint8_array_tag(*b"XC", &[10, 20, 30]);
        assert_eq!(find_string_tag(&aux, b"XC"), None);
    }

    // --- H (hex string) type ---

    #[test]
    fn test_find_string_tag_cannot_find_h_type() {
        // H-type hex string: tag, 'H', hex bytes, NUL
        let aux: &[u8] = b"XHH1A2B\x00";
        // find_string_tag only handles Z, not H
        assert_eq!(find_string_tag(aux, b"XH"), None);
    }

    // --- Z (NUL-terminated string) type ---

    #[test]
    fn test_find_string_tag_finds_z_type() {
        let aux: &[u8] = b"RXZhello\x00";
        assert_eq!(find_string_tag(aux, b"RX"), Some(b"hello".as_ref()));
    }

    // --- c (int8) type ---

    #[test]
    fn test_find_string_tag_cannot_find_c_type() {
        let aux: &[u8] = &[b'X', b'c', b'c', 0xFE]; // Xc:c:-2
        assert_eq!(find_string_tag(aux, b"Xc"), None);
    }

    // --- C (uint8) type ---

    #[test]
    fn test_find_string_tag_cannot_find_upper_c_type() {
        let aux: &[u8] = &[b'X', b'C', b'C', 200];
        assert_eq!(find_string_tag(aux, b"XC"), None);
    }

    // --- s (int16) type ---

    #[test]
    fn test_find_string_tag_cannot_find_s_type() {
        let val = 300i16.to_le_bytes();
        let aux: &[u8] = &[b'X', b's', b's', val[0], val[1]];
        assert_eq!(find_string_tag(aux, b"Xs"), None);
    }

    // --- S (uint16) type ---

    #[test]
    fn test_find_string_tag_cannot_find_upper_s_type() {
        let val = 50_000_u16.to_le_bytes();
        let aux: &[u8] = &[b'X', b'S', b'S', val[0], val[1]];
        assert_eq!(find_string_tag(aux, b"XS"), None);
    }

    // --- i (int32) type ---

    #[test]
    fn test_find_string_tag_cannot_find_i_type() {
        let val = 100_000_i32.to_le_bytes();
        let aux: &[u8] = &[b'X', b'i', b'i', val[0], val[1], val[2], val[3]];
        assert_eq!(find_string_tag(aux, b"Xi"), None);
    }

    // --- I (uint32) type ---

    #[test]
    fn test_find_string_tag_cannot_find_upper_i_type() {
        let val = 3_000_000_000_u32.to_le_bytes();
        let aux: &[u8] = &[b'X', b'I', b'I', val[0], val[1], val[2], val[3]];
        assert_eq!(find_string_tag(aux, b"XI"), None);
    }

    // --- A (single char) type ---

    #[test]
    fn test_find_string_tag_cannot_find_a_type() {
        let aux = [b'X', b'A', b'A', b'G']; // XA:A:G
        assert_eq!(find_string_tag(&aux, b"XA"), None);
    }

    // --- f (float) type ---

    #[test]
    fn test_find_string_tag_cannot_find_f_type() {
        let mut aux = vec![b'X', b'F', b'f'];
        aux.extend_from_slice(&1.5f32.to_le_bytes());
        assert_eq!(find_string_tag(&aux, b"XF"), None);
    }

    // --- Tag after a B-array tag (verifies traversal works) ---

    #[test]
    fn test_find_string_tag_after_b_array() {
        let mut aux = make_b_int_array_tag(*b"pa", &[1, 2, 3]);
        aux.extend_from_slice(b"RXZhello\x00");
        // Can find the Z-tag after the B-array
        assert_eq!(find_string_tag(&aux, b"RX"), Some(b"hello".as_ref()));
    }

    // ========================================================================
    // find_tag_bounds tests
    // ========================================================================

    #[test]
    fn test_find_tag_bounds() {
        // NM:C:5 (4 bytes) then RX:Z:ACG\0 (7 bytes)
        let mut aux = Vec::new();
        aux.extend_from_slice(b"NMC");
        aux.push(5);
        aux.extend_from_slice(b"RXZACG\x00");
        assert_eq!(find_tag_bounds(&aux, b"NM"), Some((0, 4)));
        assert_eq!(find_tag_bounds(&aux, b"RX"), Some((4, 11)));
        assert_eq!(find_tag_bounds(&aux, b"XX"), None);
    }

    #[test]
    fn test_find_tag_bounds_empty() {
        assert_eq!(find_tag_bounds(&[], b"RX"), None);
    }

    #[test]
    fn test_find_tag_bounds_truncated() {
        // Only 2 bytes — not enough for a tag entry
        assert_eq!(find_tag_bounds(b"RX", b"RX"), None);
    }

    #[test]
    fn test_find_tag_bounds_unknown_type_after_target() {
        // First tag: AA:C:1 (findable), second tag: BB with unknown type
        let mut aux = Vec::new();
        aux.extend_from_slice(&[b'A', b'A', b'C', 1]); // AA:C:1
        aux.extend_from_slice(&[b'B', b'B', b'?', 0]); // unknown type
        assert_eq!(find_tag_bounds(&aux, b"AA"), Some((0, 4)));
        // BB has unknown type '?' → tag_value_size returns None → loop breaks → None
        assert_eq!(find_tag_bounds(&aux, b"BB"), None);
    }

    // ========================================================================
    // find_int_tag tests
    // ========================================================================

    #[rstest]
    #[case::signed_int(b'i', &42i32.to_le_bytes(), Some(42))]
    #[case::signed_byte(b'c', &[(-5i8).cast_unsigned()], Some(-5))]
    #[case::unsigned_byte(b'C', &[10], Some(10))]
    #[case::signed_short(b's', &(-123i16).to_le_bytes(), Some(-123))]
    #[case::unsigned_short(b'S', &500u16.to_le_bytes(), Some(500))]
    #[case::unsigned_int(b'I', &100_000u32.to_le_bytes(), Some(100_000))]
    #[case::signed_int_negative(b'i', &(-99999i32).to_le_bytes(), Some(-99999))]
    #[case::float_type_returns_none(b'f', &1.0f32.to_le_bytes(), None)]
    fn test_find_int_tag_by_type(
        #[case] type_byte: u8,
        #[case] value_bytes: &[u8],
        #[case] expected: Option<i64>,
    ) {
        let mut aux = vec![b'X', b'Y', type_byte];
        aux.extend_from_slice(value_bytes);
        assert_eq!(find_int_tag(&aux, b"XY"), expected);
    }

    #[test]
    fn test_find_int_tag_not_found() {
        let aux = [b'X', b'Y', b'C', 10];
        assert_eq!(find_int_tag(&aux, b"ZZ"), None);
    }

    // --- Per-type find_int_tag tests ---

    #[test]
    fn test_find_int_tag_cannot_find_b_int8_array() {
        let aux = make_b_int8_array_tag(*b"Xc", &[-1, 0, 1]);
        assert_eq!(find_int_tag(&aux, b"Xc"), None);
    }

    #[test]
    fn test_find_int_tag_cannot_find_b_int16_array() {
        let aux = make_b_int16_array_tag(*b"Xs", &[-100, 0, 200]);
        assert_eq!(find_int_tag(&aux, b"Xs"), None);
    }

    #[test]
    fn test_find_int_tag_cannot_find_b_uint16_array() {
        let aux = make_b_uint16_array_tag(*b"XS", &[100, 200, 300]);
        assert_eq!(find_int_tag(&aux, b"XS"), None);
    }

    #[test]
    fn test_find_int_tag_cannot_find_b_uint32_array() {
        let aux = make_b_uint32_array_tag(*b"XI", &[1000, 2000, 3000]);
        assert_eq!(find_int_tag(&aux, b"XI"), None);
    }

    #[test]
    fn test_find_int_tag_cannot_find_b_int_array() {
        let aux = make_b_int_array_tag(*b"pa", &[0, 27_056_961, 0, 207, 60005, 1]);
        // BUG: find_int_tag returns None for B:i tags
        assert_eq!(find_int_tag(&aux, b"pa"), None);
    }

    #[test]
    fn test_find_int_tag_cannot_find_b_float_array() {
        let aux = make_b_float_array_tag(*b"XF", &[1.0, 2.5, 3.0]);
        assert_eq!(find_int_tag(&aux, b"XF"), None);
    }

    #[test]
    fn test_find_int_tag_cannot_find_b_uint8_array() {
        let aux = make_b_uint8_array_tag(*b"XC", &[10, 20, 30]);
        assert_eq!(find_int_tag(&aux, b"XC"), None);
    }

    #[test]
    fn test_find_int_tag_cannot_find_h_type() {
        let aux: &[u8] = b"XHH1A2B\x00";
        assert_eq!(find_int_tag(aux, b"XH"), None);
    }

    #[test]
    fn test_find_int_tag_cannot_find_z_type() {
        let aux: &[u8] = b"RXZhello\x00";
        assert_eq!(find_int_tag(aux, b"RX"), None);
    }

    #[test]
    fn test_find_int_tag_finds_c_type() {
        let aux: &[u8] = &[b'X', b'c', b'c', 5];
        assert_eq!(find_int_tag(aux, b"Xc"), Some(5));
    }

    #[test]
    fn test_find_int_tag_finds_upper_c_type() {
        let aux: &[u8] = &[b'N', b'M', b'C', 42];
        assert_eq!(find_int_tag(aux, b"NM"), Some(42));
    }

    #[test]
    fn test_find_int_tag_finds_s_type() {
        let val = 300i16.to_le_bytes();
        let aux: &[u8] = &[b'X', b's', b's', val[0], val[1]];
        assert_eq!(find_int_tag(aux, b"Xs"), Some(300));
    }

    #[test]
    fn test_find_int_tag_finds_upper_s_type() {
        let val = 50_000_u16.to_le_bytes();
        let aux: &[u8] = &[b'X', b'S', b'S', val[0], val[1]];
        assert_eq!(find_int_tag(aux, b"XS"), Some(50_000));
    }

    #[test]
    fn test_find_int_tag_finds_i_type() {
        let val = 100_000_i32.to_le_bytes();
        let aux: &[u8] = &[b'X', b'i', b'i', val[0], val[1], val[2], val[3]];
        assert_eq!(find_int_tag(aux, b"Xi"), Some(100_000));
    }

    #[test]
    fn test_find_int_tag_finds_upper_i_type() {
        let val = 3_000_000_000_u32.to_le_bytes();
        let aux: &[u8] = &[b'X', b'I', b'I', val[0], val[1], val[2], val[3]];
        // find_int_tag returns i64, so this should work for large uint32 values
        assert_eq!(find_int_tag(aux, b"XI"), Some(3_000_000_000));
    }

    #[test]
    fn test_find_int_tag_cannot_find_a_type() {
        let aux = [b'X', b'A', b'A', b'G'];
        assert_eq!(find_int_tag(&aux, b"XA"), None);
    }

    #[test]
    fn test_find_int_tag_cannot_find_f_type() {
        let mut aux = vec![b'X', b'F', b'f'];
        aux.extend_from_slice(&1.5f32.to_le_bytes());
        assert_eq!(find_int_tag(&aux, b"XF"), None);
    }

    #[test]
    fn test_find_int_tag_after_b_array() {
        let mut aux = make_b_int_array_tag(*b"pa", &[1, 2, 3]);
        aux.extend_from_slice(&[b'N', b'M', b'C', 5]); // NM:C:5
        assert_eq!(find_int_tag(&aux, b"NM"), Some(5));
    }

    // ========================================================================
    // find_uint8_tag tests
    // ========================================================================

    #[test]
    fn test_find_uint8_tag() {
        let aux = [b'M', b'Q', b'C', 30];
        assert_eq!(find_uint8_tag(&aux, b"MQ"), Some(30));
    }

    #[rstest]
    #[case::unsigned_byte(b'C', &[30], Some(30), Some(30))]
    #[case::signed_byte(b'c', &[30], None, Some(30))]
    #[case::signed_short(b's', &(255i16).to_le_bytes(), None, Some(255))]
    fn test_mq_tag_type_dispatch(
        #[case] type_byte: u8,
        #[case] value_bytes: &[u8],
        #[case] expected_uint8: Option<u8>,
        #[case] expected_int: Option<i64>,
    ) {
        let mut aux = vec![b'M', b'Q', type_byte];
        aux.extend_from_slice(value_bytes);
        assert_eq!(find_uint8_tag(&aux, b"MQ"), expected_uint8);
        assert_eq!(find_int_tag(&aux, b"MQ"), expected_int);
    }

    #[test]
    fn test_find_uint8_tag_wrong_type() {
        // Tag exists as type 'i' not 'C' — should not match
        let mut aux = vec![b'M', b'Q', b'i'];
        aux.extend_from_slice(&42i32.to_le_bytes());
        assert_eq!(find_uint8_tag(&aux, b"MQ"), None);
    }

    #[test]
    fn test_find_uint8_tag_not_found() {
        let aux = [b'A', b'B', b'C', 10];
        assert_eq!(find_uint8_tag(&aux, b"XY"), None);
    }

    // ========================================================================
    // find_float_tag tests
    // ========================================================================

    #[test]
    fn test_find_float_tag() {
        let val: f32 = 1.234;
        let mut aux = vec![b'X', b'F', b'f'];
        aux.extend_from_slice(&val.to_le_bytes());
        let result = find_float_tag(&aux, b"XF").unwrap();
        assert!((result - val).abs() < 0.001);
    }

    #[test]
    fn test_find_float_tag_not_found() {
        let aux = [b'A', b'B', b'C', 10];
        assert_eq!(find_float_tag(&aux, b"XY"), None);
    }

    #[test]
    fn test_find_float_tag_wrong_type() {
        // Tag exists as 'C' (uint8), not 'f' — should not match
        let aux = [b'X', b'F', b'C', 42];
        assert_eq!(find_float_tag(&aux, b"XF"), None);
    }

    // ========================================================================
    // find_mc_tag tests
    // ========================================================================

    #[test]
    fn test_find_mc_tag_present() {
        // MC:Z:10M5S\0
        let aux = b"MCZ10M5S\x00";
        assert_eq!(find_mc_tag(aux), Some("10M5S"));
    }

    #[test]
    fn test_find_mc_tag_absent() {
        // Some other tag but no MC
        let aux = b"NMC\x05\x00\x00\x00"; // NM:i:5
        assert_eq!(find_mc_tag(aux), None);
    }

    #[test]
    fn test_find_mc_tag_after_other_tags() {
        // NM:C:5 (1 byte) then MC:Z:15M\0
        // NM tag: b'N', b'M', b'C', 5
        // MC tag: b'M', b'C', b'Z', b'1', b'5', b'M', 0
        let mut aux = Vec::new();
        aux.extend_from_slice(b"NMC"); // tag NM, type C (unsigned byte)
        aux.push(5); // value
        aux.extend_from_slice(b"MCZ15M\x00"); // MC:Z:15M
        assert_eq!(find_mc_tag(&aux), Some("15M"));
    }

    #[test]
    fn test_find_mc_tag_non_z_type() {
        // MC tag with non-Z type → None
        let aux = [b'M', b'C', b'C', 42]; // MC:C:42
        assert_eq!(find_mc_tag(&aux), None);
    }

    #[test]
    fn test_find_mc_tag_empty() {
        assert_eq!(find_mc_tag(&[]), None);
    }

    // ========================================================================
    // find_mc_tag_in_record tests
    // ========================================================================

    #[test]
    fn test_find_mc_tag_in_record_no_aux() {
        let rec = make_bam_bytes(0, 0, 0, b"rea", &[], 0, -1, -1, &[]);
        assert_eq!(find_mc_tag_in_record(&rec), None);
    }

    // ========================================================================
    // extract_aux_string_tags tests
    // ========================================================================

    #[test]
    fn test_extract_aux_string_tags_all_found() {
        // Build aux data with RG, CB (cell barcode), and MC tags
        let mut aux = Vec::new();
        aux.extend_from_slice(b"RGZsample1\x00");
        aux.extend_from_slice(b"CBZcell42\x00");
        aux.extend_from_slice(b"MCZ10M5S\x00");
        let result = extract_aux_string_tags(&aux, b"CB");
        assert_eq!(result.rg, Some(b"sample1".as_ref()));
        assert_eq!(result.cell, Some(b"cell42".as_ref()));
        assert_eq!(result.mc, Some("10M5S"));
    }

    #[test]
    fn test_extract_aux_string_tags_early_exit() {
        // All three tags present, should stop scanning after finding all
        let mut aux = Vec::new();
        aux.extend_from_slice(b"RGZrg1\x00");
        aux.extend_from_slice(b"CBZbc1\x00");
        aux.extend_from_slice(b"MCZ5M\x00");
        // Add extra tags that should never be reached
        aux.extend_from_slice(&[b'X', b'Y', b'C', 99]);
        let result = extract_aux_string_tags(&aux, b"CB");
        assert_eq!(result.rg, Some(b"rg1".as_ref()));
        assert_eq!(result.cell, Some(b"bc1".as_ref()));
        assert_eq!(result.mc, Some("5M"));
    }

    #[test]
    fn test_extract_aux_string_tags_partial() {
        // Only RG present, others missing
        let aux = b"RGZsample\x00";
        let result = extract_aux_string_tags(aux, b"CB");
        assert_eq!(result.rg, Some(b"sample".as_ref()));
        assert!(result.cell.is_none());
        assert!(result.mc.is_none());
    }

    #[test]
    fn test_extract_aux_string_tags_with_non_string_tags() {
        // Mix of string and non-string tags
        let mut aux = Vec::new();
        aux.extend_from_slice(&[b'N', b'M', b'C', 5]); // NM:C:5
        aux.extend_from_slice(b"RGZlib1\x00"); // RG:Z:lib1
        aux.extend_from_slice(&[b'A', b'S', b'C', 30]); // AS:C:30
        aux.extend_from_slice(b"MCZ20M\x00"); // MC:Z:20M
        let result = extract_aux_string_tags(&aux, b"CB");
        assert_eq!(result.rg, Some(b"lib1".as_ref()));
        assert!(result.cell.is_none());
        assert_eq!(result.mc, Some("20M"));
    }

    #[test]
    fn test_extract_aux_string_tags_empty() {
        let result = extract_aux_string_tags(&[], b"CB");
        assert!(result.rg.is_none());
        assert!(result.cell.is_none());
        assert!(result.mc.is_none());
    }

    #[test]
    fn test_extract_aux_string_tags_mc_invalid_utf8() {
        // MC tag with invalid UTF-8 bytes → mc should be None
        let mut aux = Vec::new();
        aux.extend_from_slice(b"MCZ");
        aux.extend_from_slice(&[0xFF, 0xFE, 0xFD]); // invalid UTF-8
        aux.push(0); // null terminator
        let result = extract_aux_string_tags(&aux, b"CB");
        assert!(result.mc.is_none()); // from_utf8 fails
    }

    #[test]
    fn test_extract_aux_string_tags_truncated_z() {
        // RG:Z:lib1 with no null → should break out, no tags found
        let aux = b"RGZlib1";
        let result = extract_aux_string_tags(aux, b"CB");
        assert!(result.rg.is_none());
    }

    // ========================================================================
    // find_tag_type tests (per-type)
    // ========================================================================

    #[test]
    fn test_find_tag_type_finds_b_int8_array() {
        let aux = make_b_int8_array_tag(*b"Xc", &[-1, 0, 1]);
        assert_eq!(find_tag_type(&aux, b"Xc"), Some(b'B'));
    }

    #[test]
    fn test_find_tag_type_finds_b_int16_array() {
        let aux = make_b_int16_array_tag(*b"Xs", &[-100, 0, 200]);
        assert_eq!(find_tag_type(&aux, b"Xs"), Some(b'B'));
    }

    #[test]
    fn test_find_tag_type_finds_b_uint16_array() {
        let aux = make_b_uint16_array_tag(*b"XS", &[100, 200, 300]);
        assert_eq!(find_tag_type(&aux, b"XS"), Some(b'B'));
    }

    #[test]
    fn test_find_tag_type_finds_b_uint32_array() {
        let aux = make_b_uint32_array_tag(*b"XI", &[1000, 2000, 3000]);
        assert_eq!(find_tag_type(&aux, b"XI"), Some(b'B'));
    }

    #[test]
    fn test_find_tag_type_finds_b_int_array() {
        let aux = make_b_int_array_tag(*b"pa", &[0, 27_056_961, 0, 207, 60005, 1]);
        // find_tag_type correctly finds B-array tags
        assert_eq!(find_tag_type(&aux, b"pa"), Some(b'B'));
    }

    #[test]
    fn test_find_tag_type_finds_b_float_array() {
        let aux = make_b_float_array_tag(*b"XF", &[1.0, 2.5, 3.0]);
        assert_eq!(find_tag_type(&aux, b"XF"), Some(b'B'));
    }

    #[test]
    fn test_find_tag_type_finds_b_uint8_array() {
        let aux = make_b_uint8_array_tag(*b"XC", &[10, 20, 30]);
        assert_eq!(find_tag_type(&aux, b"XC"), Some(b'B'));
    }

    #[test]
    fn test_find_tag_type_finds_h_type() {
        let aux: &[u8] = b"XHH1A2B\x00";
        assert_eq!(find_tag_type(aux, b"XH"), Some(b'H'));
    }

    #[test]
    fn test_find_tag_type_finds_z_type() {
        let aux: &[u8] = b"RXZhello\x00";
        assert_eq!(find_tag_type(aux, b"RX"), Some(b'Z'));
    }

    #[test]
    fn test_find_tag_type_finds_c_type() {
        let aux: &[u8] = &[b'X', b'c', b'c', 5];
        assert_eq!(find_tag_type(aux, b"Xc"), Some(b'c'));
    }

    #[test]
    fn test_find_tag_type_finds_upper_c_type() {
        let aux: &[u8] = &[b'N', b'M', b'C', 42];
        assert_eq!(find_tag_type(aux, b"NM"), Some(b'C'));
    }

    #[test]
    fn test_find_tag_type_finds_s_type() {
        let val = 300i16.to_le_bytes();
        let aux: &[u8] = &[b'X', b's', b's', val[0], val[1]];
        assert_eq!(find_tag_type(aux, b"Xs"), Some(b's'));
    }

    #[test]
    fn test_find_tag_type_finds_upper_s_type() {
        let val = 50_000_u16.to_le_bytes();
        let aux: &[u8] = &[b'X', b'S', b'S', val[0], val[1]];
        assert_eq!(find_tag_type(aux, b"XS"), Some(b'S'));
    }

    #[test]
    fn test_find_tag_type_finds_i_type() {
        let val = 100_000_i32.to_le_bytes();
        let aux: &[u8] = &[b'X', b'i', b'i', val[0], val[1], val[2], val[3]];
        assert_eq!(find_tag_type(aux, b"Xi"), Some(b'i'));
    }

    #[test]
    fn test_find_tag_type_finds_upper_i_type() {
        let val = 3_000_000_000_u32.to_le_bytes();
        let aux: &[u8] = &[b'X', b'I', b'I', val[0], val[1], val[2], val[3]];
        assert_eq!(find_tag_type(aux, b"XI"), Some(b'I'));
    }

    #[test]
    fn test_find_tag_type_finds_a_type() {
        let aux = [b'X', b'A', b'A', b'G'];
        assert_eq!(find_tag_type(&aux, b"XA"), Some(b'A'));
    }

    #[test]
    fn test_find_tag_type_finds_f_type() {
        let mut aux = vec![b'X', b'F', b'f'];
        aux.extend_from_slice(&1.5f32.to_le_bytes());
        assert_eq!(find_tag_type(&aux, b"XF"), Some(b'f'));
    }

    #[test]
    fn test_find_tag_type_after_b_array() {
        let mut aux = make_b_int_array_tag(*b"pa", &[1, 2, 3]);
        aux.extend_from_slice(b"RXZhello\x00");
        assert_eq!(find_tag_type(&aux, b"RX"), Some(b'Z'));
    }

    // ========================================================================
    // append_string_tag tests
    // ========================================================================

    #[test]
    fn test_append_string_tag() {
        let mut rec = make_bam_bytes(0, 0, 0, b"rea", &[], 0, -1, -1, &[]);
        let orig_len = rec.len();
        append_string_tag(&mut rec, b"MI", b"12345");
        assert_eq!(rec.len(), orig_len + 2 + 1 + 5 + 1); // tag(2) + type(1) + value(5) + NUL(1)
        // Verify we can find it back
        let aux_start = aux_data_offset_from_record(&rec).unwrap();
        assert_eq!(find_string_tag(&rec[aux_start..], b"MI"), Some(b"12345".as_ref()));
    }

    // ========================================================================
    // remove_tag tests
    // ========================================================================

    #[test]
    fn test_remove_tag_present() {
        let aux = b"MIZ42\x00";
        let mut rec = make_bam_bytes(0, 0, 0, b"rea", &[], 0, -1, -1, aux);
        assert!(find_string_tag_in_record(&rec, b"MI").is_some());
        remove_tag(&mut rec, b"MI");
        assert!(find_string_tag_in_record(&rec, b"MI").is_none());
    }

    #[test]
    fn test_remove_tag_absent() {
        let mut rec = make_bam_bytes(0, 0, 0, b"rea", &[], 0, -1, -1, &[]);
        let orig_len = rec.len();
        remove_tag(&mut rec, b"MI"); // should be no-op
        assert_eq!(rec.len(), orig_len);
    }

    #[test]
    fn test_remove_tag_between_tags() {
        // Three tags: AA:C:1, BB:Z:hi\0, CC:C:2
        let mut aux = Vec::new();
        aux.extend_from_slice(&[b'A', b'A', b'C', 1]); // AA:C:1
        aux.extend_from_slice(b"BBZhi\x00"); // BB:Z:hi
        aux.extend_from_slice(&[b'C', b'C', b'C', 2]); // CC:C:2
        let mut rec = make_bam_bytes(0, 0, 0, b"rea", &[], 0, -1, -1, &aux);
        remove_tag(&mut rec, b"BB");
        // AA and CC should still be findable
        let aux_start = aux_data_offset_from_record(&rec).unwrap();
        assert_eq!(find_uint8_tag(&rec[aux_start..], b"AA"), Some(1));
        assert_eq!(find_uint8_tag(&rec[aux_start..], b"CC"), Some(2));
        assert!(find_string_tag(&rec[aux_start..], b"BB").is_none());
    }

    #[test]
    fn test_remove_tag_no_aux_data() {
        // Record where aux_data_offset >= record length
        let mut rec = make_bam_bytes(0, 0, 0, b"rea", &[], 0, -1, -1, &[]);
        let orig_len = rec.len();
        remove_tag(&mut rec, b"MI");
        assert_eq!(rec.len(), orig_len); // no-op
    }

    // ========================================================================
    // update_string_tag tests
    // ========================================================================

    #[test]
    fn test_update_string_tag_existing() {
        let aux = b"RXZold\x00";
        let mut rec = make_bam_bytes(0, 0, 0, b"rea", &[], 0, -1, -1, aux);
        update_string_tag(&mut rec, b"RX", b"newvalue");
        assert_eq!(find_string_tag_in_record(&rec, b"RX"), Some(b"newvalue".as_ref()));
    }

    #[test]
    fn test_update_string_tag_new() {
        let mut rec = make_bam_bytes(0, 0, 0, b"rea", &[], 0, -1, -1, &[]);
        update_string_tag(&mut rec, b"RX", b"added");
        assert_eq!(find_string_tag_in_record(&rec, b"RX"), Some(b"added".as_ref()));
    }

    #[test]
    fn test_update_string_tag_same_length() {
        // Test same-length fast path (copy_from_slice, no splice)
        let aux = b"RXZold1\x00";
        let mut rec = make_bam_bytes(0, 0, 0, b"rea", &[], 0, -1, -1, aux);
        let orig_len = rec.len();
        update_string_tag(&mut rec, b"RX", b"new2");
        // Record length should not change for same-length update
        assert_eq!(rec.len(), orig_len);
        assert_eq!(find_string_tag_in_record(&rec, b"RX"), Some(b"new2".as_ref()));
    }

    #[test]
    fn test_update_string_tag_different_length_shorter() {
        // Test splice path: old value longer than new value
        let aux = b"RXZlongvalue\x00";
        let mut rec = make_bam_bytes(0, 0, 0, b"rea", &[], 0, -1, -1, aux);
        update_string_tag(&mut rec, b"RX", b"hi");
        assert_eq!(find_string_tag_in_record(&rec, b"RX"), Some(b"hi".as_ref()));
    }

    #[test]
    fn test_update_string_tag_different_length_longer() {
        // Test splice path: old value shorter than new value
        let aux = b"RXZhi\x00";
        let mut rec = make_bam_bytes(0, 0, 0, b"rea", &[], 0, -1, -1, aux);
        update_string_tag(&mut rec, b"RX", b"longvalue");
        assert_eq!(find_string_tag_in_record(&rec, b"RX"), Some(b"longvalue".as_ref()));
    }

    #[test]
    fn test_update_string_tag_preserves_other_tags() {
        // Test that updating one tag doesn't corrupt adjacent tags
        let mut aux = Vec::new();
        aux.extend_from_slice(&[b'A', b'A', b'C', 1]); // AA:C:1
        aux.extend_from_slice(b"RXZold\x00"); // RX:Z:old
        aux.extend_from_slice(&[b'C', b'C', b'C', 2]); // CC:C:2
        let mut rec = make_bam_bytes(0, 0, 0, b"rea", &[], 0, -1, -1, &aux);
        update_string_tag(&mut rec, b"RX", b"newval");
        let aux_start = aux_data_offset_from_record(&rec).unwrap();
        assert_eq!(find_uint8_tag(&rec[aux_start..], b"AA"), Some(1));
        assert_eq!(find_string_tag(&rec[aux_start..], b"RX"), Some(b"newval".as_ref()));
        assert_eq!(find_uint8_tag(&rec[aux_start..], b"CC"), Some(2));
    }

    #[test]
    fn test_update_string_tag_no_aux_appends() {
        let mut rec = make_bam_bytes(0, 0, 0, b"rea", &[], 0, -1, -1, &[]);
        update_string_tag(&mut rec, b"MI", b"42");
        assert_eq!(find_string_tag_in_record(&rec, b"MI"), Some(b"42".as_ref()));
    }

    // ========================================================================
    // append_int_tag tests
    // ========================================================================

    #[test]
    fn test_append_int_tag_i8() {
        let mut rec = Vec::new();
        append_int_tag(&mut rec, b"cD", 42);
        assert_eq!(rec, [b'c', b'D', b'c', 42]);
    }

    #[test]
    fn test_append_int_tag_negative_i8() {
        let mut rec = Vec::new();
        append_int_tag(&mut rec, b"cM", -5);
        assert_eq!(rec, [b'c', b'M', b'c', (-5i8).cast_unsigned()]);
    }

    #[test]
    fn test_append_int_tag_i16() {
        let mut rec = Vec::new();
        append_int_tag(&mut rec, b"cD", 200);
        let v = 200i16.to_le_bytes();
        assert_eq!(rec, [b'c', b'D', b's', v[0], v[1]]);
    }

    #[test]
    fn test_append_int_tag_i32() {
        let mut rec = Vec::new();
        append_int_tag(&mut rec, b"cD", 100_000);
        let v = 100_000i32.to_le_bytes();
        assert_eq!(rec, [b'c', b'D', b'i', v[0], v[1], v[2], v[3]]);
    }

    // ========================================================================
    // append_float_tag tests
    // ========================================================================

    #[test]
    fn test_append_float_tag() {
        let mut rec = Vec::new();
        append_float_tag(&mut rec, b"cE", 0.05);
        let v = 0.05f32.to_le_bytes();
        assert_eq!(rec, [b'c', b'E', b'f', v[0], v[1], v[2], v[3]]);
    }

    // ========================================================================
    // append_i16_array_tag tests
    // ========================================================================

    #[test]
    fn test_append_i16_array_tag_empty() {
        let mut rec = Vec::new();
        append_i16_array_tag(&mut rec, b"cd", &[]);
        assert_eq!(rec, [b'c', b'd', b'B', b's', 0, 0, 0, 0]);
    }

    #[test]
    fn test_append_i16_array_tag_values() {
        let mut rec = Vec::new();
        append_i16_array_tag(&mut rec, b"cd", &[10, 20, 5]);
        let mut expected = vec![b'c', b'd', b'B', b's'];
        expected.extend_from_slice(&3u32.to_le_bytes());
        expected.extend_from_slice(&10i16.to_le_bytes());
        expected.extend_from_slice(&20i16.to_le_bytes());
        expected.extend_from_slice(&5i16.to_le_bytes());
        assert_eq!(rec, expected);
    }

    // ========================================================================
    // append_phred33_string_tag tests
    // ========================================================================

    #[test]
    fn test_append_phred33_string_tag() {
        let mut buf = Vec::new();
        append_phred33_string_tag(&mut buf, b"aq", &[0, 10, 30, 40]);
        assert_eq!(buf[0], b'a');
        assert_eq!(buf[1], b'q');
        assert_eq!(buf[2], b'Z');
        assert_eq!(buf[3], b'!'); // 0 + 33
        assert_eq!(buf[4], b'+'); // 10 + 33
        assert_eq!(buf[5], b'?'); // 30 + 33
        assert_eq!(buf[6], b'I'); // 40 + 33
        assert_eq!(buf[7], 0); // NUL
    }

    // ========================================================================
    // Dedup pa-tag validation bug reproduction
    // ========================================================================

    #[test]
    fn test_dedup_pa_tag_check_fails_on_b_array() {
        // Simulate the exact pa tag as produced by fgumi zipper: pa:B:i,0,27_056_961,0,207,60005,1
        let aux = make_b_int_array_tag(*b"pa", &[0, 27_056_961, 0, 207, 60005, 1]);
        let pa_tag_bytes: [u8; 2] = *b"pa";

        // This is the exact check from dedup.rs:932-934
        let found_by_dedup = find_string_tag(&aux, &pa_tag_bytes).is_some()
            || find_int_tag(&aux, &pa_tag_bytes).is_some();

        // BUG: dedup thinks the pa tag is missing even though it's present
        assert!(!found_by_dedup, "dedup check should fail to find B:i pa tag");

        // But find_tag_type correctly finds it
        assert!(find_tag_type(&aux, &pa_tag_bytes).is_some());
    }
}
