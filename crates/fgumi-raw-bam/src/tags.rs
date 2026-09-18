use fgumi_tag::{AsTagBytes, SamTag};

use crate::fields::{
    RawRecordMut, RawRecordView, TAG_FIXED_SIZES, aux_data_offset_from_record, aux_data_slice,
    nul_offset, tag_value_size,
};

/// Find a tag's position and type byte in auxiliary data.
///
/// Returns `(offset, type_byte)` where `offset` is the position of the tag entry
/// (the first byte of the 2-byte tag identifier). Returns `None` if not found.
#[must_use]
fn find_tag_position(aux_data: &[u8], tag: [u8; 2]) -> Option<(usize, u8)> {
    // Compare as u16 — a single integer compare instead of a slice compare.
    // LE byte order is preserved: both values are built from the same two bytes
    // in the same order, so equality is equivalent to the byte-slice compare.
    let tag_u16 = u16::from_le_bytes(tag);
    let mut p = 0;
    while p + 3 <= aux_data.len() {
        let entry_u16 = u16::from_le_bytes([aux_data[p], aux_data[p + 1]]);
        let val_type = aux_data[p + 2];

        if entry_u16 == tag_u16 {
            return Some((p, val_type));
        }

        if let Some(size) = tag_value_size(val_type, &aux_data[p + 3..]) {
            p += 3 + size;
        } else {
            break;
        }
    }
    None
}

/// Find a string (Z-type) tag in auxiliary data, returning value bytes without null terminator.
#[inline]
#[must_use]
pub fn find_string_tag(aux_data: &[u8], tag: impl AsTagBytes) -> Option<&[u8]> {
    let tag = tag.as_tag_bytes();
    let (p, val_type) = find_tag_position(aux_data, *tag)?;
    if val_type != b'Z' {
        return None;
    }
    let start = p + 3;
    let end = nul_offset(&aux_data[start..])?;
    Some(&aux_data[start..start + end])
}

/// Find a string (Z-type) tag and return its value position as `(offset, len)` —
/// where `offset` is the offset into `aux_data` of the first value byte (after the
/// 2-byte tag identifier and 1-byte type byte) and `len` is the value length in
/// bytes, excluding the NUL terminator.
///
/// Used to cache a Z-tag's position so a later caller can slice the value without
/// re-scanning aux data.
///
/// Returns `None` if the tag is absent, is not Z-type, has no NUL terminator, or
/// would not fit in the returned types (`u32` offset, `u16` length). BAM records
/// are bounded well under 4 GiB and aux tag values well under 64 KiB, so the
/// width-fit checks are defensive only.
#[inline]
#[must_use]
pub fn find_string_tag_position(aux_data: &[u8], tag: impl AsTagBytes) -> Option<(u32, u16)> {
    let tag = tag.as_tag_bytes();
    let (p, val_type) = find_tag_position(aux_data, *tag)?;
    if val_type != b'Z' {
        return None;
    }
    let start = p + 3;
    let len = nul_offset(&aux_data[start..])?;
    Some((u32::try_from(start).ok()?, u16::try_from(len).ok()?))
}

/// Check whether a tag exists in auxiliary data, returning its type byte if found.
///
/// Returns `Some(type_byte)` (e.g. `b'Z'`, `b'C'`, `b'i'`) if the tag is present,
/// `None` if the tag is absent.
#[inline]
#[must_use]
pub fn find_tag_type(aux_data: &[u8], tag: impl AsTagBytes) -> Option<u8> {
    let tag = tag.as_tag_bytes();
    find_tag_position(aux_data, *tag).map(|(_, val_type)| val_type)
}

/// Find a string tag in a complete BAM record.
#[inline]
#[must_use]
pub fn find_string_tag_in_record(bam: &[u8], tag: impl AsTagBytes) -> Option<&[u8]> {
    let tag = tag.as_tag_bytes();
    let aux = aux_data_slice(bam);
    if aux.is_empty() {
        return None;
    }
    find_string_tag(aux, tag)
}

/// Find two string (Z-type) tags in a complete BAM record in a **single** aux-data
/// walk, returning each value's bytes (without the NUL terminator) or `None` if that
/// tag is absent or not Z-type.
///
/// Equivalent to calling [`find_string_tag_in_record`] twice, but it walks the aux
/// block only once — relevant on hot paths that need both an MI tag and a cell-barcode
/// tag per record. The returned slices borrow the record's aux data.
#[inline]
#[must_use]
pub fn find_two_string_tags_in_record(
    bam: &[u8],
    first: impl AsTagBytes,
    second: impl AsTagBytes,
) -> (Option<&[u8]>, Option<&[u8]>) {
    let first = u16::from_le_bytes(*first.as_tag_bytes());
    let second = u16::from_le_bytes(*second.as_tag_bytes());
    let aux = aux_data_slice(bam);

    let mut first_done = false;
    let mut second_done = false;
    let mut first_val: Option<&[u8]> = None;
    let mut second_val: Option<&[u8]> = None;
    let mut p = 0;
    while p + 3 <= aux.len() {
        let entry_u16 = u16::from_le_bytes([aux[p], aux[p + 1]]);
        let val_type = aux[p + 2];

        // Each tag resolves on its FIRST matching entry, exactly as
        // `find_string_tag_in_record` does (it returns the first tag match and
        // yields a value only if that match is a NUL-terminated `Z` entry). A
        // non-`Z` (or unterminated) first match resolves the tag to `None`; we
        // must NOT skip it and pick up a later `Z` duplicate, or the two-tag
        // walk would disagree with two single-tag lookups on malformed aux.
        let matches_first = entry_u16 == first && !first_done;
        let matches_second = entry_u16 == second && !second_done;
        if matches_first || matches_second {
            let value = if val_type == b'Z' {
                let start = p + 3;
                aux[start..].iter().position(|&b| b == 0).map(|len| &aux[start..start + len])
            } else {
                None
            };
            // Independent (not `else`) so a single entry resolves both slots
            // when `first == second`, preserving the "call twice" contract.
            if matches_first {
                first_done = true;
                first_val = value;
            }
            if matches_second {
                second_done = true;
                second_val = value;
            }
        }

        // Stop once both tags have been resolved (each to a value or `None`).
        if first_done && second_done {
            break;
        }

        match tag_value_size(val_type, &aux[p + 3..]) {
            Some(size) => p += 3 + size,
            None => break,
        }
    }
    (first_val, second_val)
}

/// Find the byte range `[start, end)` of an entire tag entry (tag+type+value) in aux data.
///
/// Returns offsets relative to the start of `aux_data`.
/// Returns `None` if the tag is not found.
#[inline]
#[must_use]
pub(crate) fn find_tag_bounds(aux_data: &[u8], tag: impl AsTagBytes) -> Option<(usize, usize)> {
    let tag = tag.as_tag_bytes();
    let (p, val_type) = find_tag_position(aux_data, *tag)?;
    let size = tag_value_size(val_type, &aux_data[p + 3..])?;
    Some((p, p + 3 + size))
}

/// Find a uint8 (C-type) tag value in auxiliary data.
#[inline]
#[must_use]
pub fn find_uint8_tag(aux_data: &[u8], tag: impl AsTagBytes) -> Option<u8> {
    let tag = tag.as_tag_bytes();
    let (p, val_type) = find_tag_position(aux_data, *tag)?;
    if val_type == b'C' && p + 4 <= aux_data.len() { Some(aux_data[p + 3]) } else { None }
}

/// Find a float (f-type) tag value in auxiliary data.
#[inline]
#[must_use]
pub fn find_float_tag(aux_data: &[u8], tag: impl AsTagBytes) -> Option<f32> {
    let tag = tag.as_tag_bytes();
    let (p, val_type) = find_tag_position(aux_data, *tag)?;
    if val_type == b'f' && p + 7 <= aux_data.len() {
        Some(f32::from_le_bytes([
            aux_data[p + 3],
            aux_data[p + 4],
            aux_data[p + 5],
            aux_data[p + 6],
        ]))
    } else {
        None
    }
}

/// Find an integer tag value in auxiliary data.
///
/// Supports signed/unsigned byte, short, and int types (c/C/s/S/i/I).
#[inline]
#[must_use]
pub fn find_int_tag(aux_data: &[u8], tag: impl AsTagBytes) -> Option<i64> {
    let tag = tag.as_tag_bytes();
    let (p, val_type) = find_tag_position(aux_data, *tag)?;
    extract_int_value(aux_data, p, val_type)
}

/// Extract an integer value at position `p` with the given type byte.
///
/// Shared by [`find_int_tag`] and `find_mi_tag`, and by callers that walk the
/// aux block themselves and already hold the tag's position and type byte — they
/// get the decode ladder without paying for a second `find_tag_position` scan.
///
/// Returns `None` for a non-integer type byte, and for any `p` whose value bytes
/// do not lie wholly within `aux_data`.
#[must_use]
pub fn extract_int_value(aux_data: &[u8], p: usize, val_type: u8) -> Option<i64> {
    // The value bytes start after the 2-byte tag id and 1-byte type at `p`.
    // `checked_add` rejects a near-`usize::MAX` `p` (which can reach here from a
    // caller rather than a `find_tag_position` scan) instead of wrapping.
    decode_int_value(val_type, aux_data.get(p.checked_add(3)?..)?)
}

/// Decode an integer aux value from its raw value bytes (the bytes *after* the
/// two-byte tag id and one-byte type), given the BAM type byte.
///
/// This is the decode ladder shared by [`extract_int_value`] (which slices the
/// value out of a full aux block first) and by callers that already hold a tag's
/// value bytes — e.g. from an [`AuxTagsIter`]/[`TagEntry`] walk — so they need
/// not re-scan for the position. Returns `None` for a non-integer type byte or a
/// value slice too short for the type's width.
#[must_use]
pub fn decode_int_value(val_type: u8, value_bytes: &[u8]) -> Option<i64> {
    Some(match val_type {
        b'c' => i64::from(value_bytes.first()?.cast_signed()),
        b'C' => i64::from(*value_bytes.first()?),
        b's' => i64::from(i16::from_le_bytes(value_bytes.get(0..2)?.try_into().ok()?)),
        b'S' => i64::from(u16::from_le_bytes(value_bytes.get(0..2)?.try_into().ok()?)),
        b'i' => i64::from(i32::from_le_bytes(value_bytes.get(0..4)?.try_into().ok()?)),
        b'I' => i64::from(u32::from_le_bytes(value_bytes.get(0..4)?.try_into().ok()?)),
        _ => return None,
    })
}

/// Decode an integer aux value directly from its type byte and value bytes.
///
/// The value-bytes form of [`extract_int_value`], for callers that already hold
/// a tag entry's `(type_byte, value_bytes)` (e.g. from an [`AuxTagsIter`] walk)
/// and want the decoded integer without re-deriving the position. Returns `None`
/// for a non-integer type byte, or if `value_bytes` is too short for the type.
#[must_use]
fn int_from_value_bytes(val_type: u8, value_bytes: &[u8]) -> Option<i64> {
    match val_type {
        b'c' => value_bytes.first().map(|&b| i64::from(b.cast_signed())),
        b'C' => value_bytes.first().map(|&b| i64::from(b)),
        b's' => value_bytes.get(..2).map(|b| i64::from(i16::from_le_bytes([b[0], b[1]]))),
        b'S' => value_bytes.get(..2).map(|b| i64::from(u16::from_le_bytes([b[0], b[1]]))),
        b'i' => {
            value_bytes.get(..4).map(|b| i64::from(i32::from_le_bytes([b[0], b[1], b[2], b[3]])))
        }
        b'I' => {
            value_bytes.get(..4).map(|b| i64::from(u32::from_le_bytes([b[0], b[1], b[2], b[3]])))
        }
        _ => None,
    }
}

/// Find MI (Molecular Identifier) tag in auxiliary data.
///
/// Returns the value as `(integer_value, is_A_suffix)`.
/// - For string values like "12345" or "12345/A" or "12345/B"
/// - For integer values, returns `(value, true)`
/// - Returns `None` if MI tag not found.
#[must_use]
pub(crate) fn find_mi_tag(aux_data: &[u8]) -> Option<(u64, bool)> {
    let (pos, val_type) = find_tag_position(aux_data, SamTag::MI.into())?;
    if val_type == b'Z' {
        // String type - parse "12345" or "12345/A" or "12345/B"
        let start = pos + 3;
        let end = nul_offset(&aux_data[start..])?;
        parse_mi_bytes(&aux_data[start..start + end])
    } else {
        // Integer types: delegate to shared extractor, reject negative values
        let v = extract_int_value(aux_data, pos, val_type)?;
        if v >= 0 {
            #[expect(clippy::cast_sign_loss, reason = "guarded by v >= 0")]
            Some((v as u64, true))
        } else {
            None
        }
    }
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

/// Find the MC (mate CIGAR) tag in auxiliary data.
///
/// Returns the CIGAR bytes, or `None` if not found.
///
/// No UTF-8 validation: `MC` is a CIGAR that every consumer parses as bytes,
/// and this accessor is the byte-level sibling of the batch extractors
/// ([`extract_aux_string_tags`], [`extract_template_aux_tags`]). Those two
/// hand the raw value over, so gating here would make the two paths disagree
/// on a non-UTF-8 value -- and `src/lib/grouper.rs` cross-checks them against
/// each other, with duplicate `MC` entries documented as the *only* case where
/// they may differ.
#[must_use]
pub(crate) fn find_mc_tag(aux_data: &[u8]) -> Option<&[u8]> {
    find_string_tag(aux_data, SamTag::MC)
}

/// Find MC tag in a complete BAM record.
#[must_use]
pub fn find_mc_tag_in_record(bam: &[u8]) -> Option<&[u8]> {
    find_mc_tag(aux_data_slice(bam))
}

/// Result of a single-pass extraction of multiple string tags from aux data.
/// Used by `compute_group_key_from_raw` to avoid scanning aux data 3 times.
pub struct AuxStringTags<'a> {
    pub rg: Option<&'a [u8]>,
    pub cell: Option<&'a [u8]>,
    pub mc: Option<&'a [u8]>,
    /// Aux-relative `(value_offset, value_len)` of the UMI tag's value bytes
    /// (excluding NUL terminator), when the caller supplied a `umi_tag` to look
    /// for. `None` when no `umi_tag` was supplied, or when the tag is absent
    /// or non-Z-typed. Aux-relative — to convert to record-relative for use
    /// with `update_string_tag_at_position` etc., add the record's aux offset.
    pub umi_position: Option<(u32, u16)>,
}

/// Extract RG, cell barcode, MC, and optionally UMI tag position in a single
/// pass over the aux data.
#[inline]
#[must_use]
pub fn extract_aux_string_tags(
    aux_data: &[u8],
    cell_tag: impl AsTagBytes,
    umi_tag: Option<[u8; 2]>,
) -> AuxStringTags<'_> {
    let cell_tag = cell_tag.as_tag_bytes();
    let mut result = AuxStringTags { rg: None, cell: None, mc: None, umi_position: None };
    let mut found = 0u8; // bit 0=RG, bit 1=cell, bit 2=MC, bit 3=UMI
    // Early-exit mask: only require UMI if the caller asked for it.
    let target = 7u8 | if umi_tag.is_some() { 8 } else { 0 };
    let mut p = 0;
    while p + 3 <= aux_data.len() {
        let t = [aux_data[p], aux_data[p + 1]];
        let val_type = aux_data[p + 2];

        if val_type == b'Z' {
            let start = p + 3;
            if let Some(end) = nul_offset(&aux_data[start..]) {
                let value = &aux_data[start..start + end];
                // Independent `if` checks (not a chained `else if`) so a caller
                // that configures overlapping selected tags -- e.g. `cell_tag ==
                // MC`, or `umi_tag` equal to RG/cell/MC -- still captures every
                // selected field from the one matching aux entry. A chained
                // `else if` would populate only the first-matching field. This
                // mirrors the sibling `extract_template_aux_tags`.
                if t == SamTag::RG {
                    result.rg = Some(value);
                    found |= 1;
                }
                if t == *cell_tag {
                    result.cell = Some(value);
                    found |= 2;
                }
                if t == SamTag::MC {
                    result.mc = Some(value);
                    found |= 4;
                }
                if matches!(umi_tag, Some(ut) if t == ut) {
                    // Width checks are defensive only (BAM records < 4 GiB,
                    // aux values < 64 KiB); silently skip if they fail.
                    if let (Ok(off_u32), Ok(len_u16)) = (u32::try_from(start), u16::try_from(end)) {
                        result.umi_position = Some((off_u32, len_u16));
                        found |= 8;
                    }
                }
                if found == target {
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

/// Result of a single-pass extraction of all tags needed for template-coordinate sorting.
/// Extracts MI (integer or string), RG, cell barcode, and MC in one scan of aux data.
pub struct TemplateAuxTags<'a> {
    /// MI tag value: (`molecular_id`, `is_A_suffix`). Defaults to `(0, true)` if not found.
    pub mi: (u64, bool),
    /// RG (read group) tag value.
    pub rg: Option<&'a [u8]>,
    /// Cell barcode tag value.
    pub cell: Option<&'a [u8]>,
    /// MC (mate CIGAR) tag value.
    pub mc: Option<&'a [u8]>,
}

/// Extract MI, RG, cell barcode, and MC tags in a single pass over aux data.
///
/// This replaces 4 separate linear scans with one, reducing the cost of aux tag
/// extraction from O(4n) to O(n) where n is the aux data length.
#[must_use]
pub fn extract_template_aux_tags(bam: &[u8], cell_tag: Option<SamTag>) -> TemplateAuxTags<'_> {
    let aux_data = aux_data_slice(bam);
    let mut result = TemplateAuxTags { mi: (0, true), rg: None, cell: None, mc: None };
    // Bits: 0=MI, 1=RG, 2=cell, 3=MC
    let target_bits: u8 = if cell_tag.is_some() { 0xF } else { 0b1011 };
    let mut found = 0u8;
    let mut p = 0;

    while p + 3 <= aux_data.len() {
        let t = [aux_data[p], aux_data[p + 1]];
        let val_type = aux_data[p + 2];

        if t == SamTag::MI {
            if val_type == b'Z' {
                let start = p + 3;
                if let Some(end) = nul_offset(&aux_data[start..]) {
                    result.mi = parse_mi_bytes(&aux_data[start..start + end]).unwrap_or((0, true));
                    p = start + end + 1;
                } else {
                    break;
                }
            } else if let Some(v) = extract_int_value(aux_data, p, val_type) {
                if v >= 0 {
                    #[expect(clippy::cast_sign_loss, reason = "guarded by v >= 0")]
                    {
                        result.mi = (v as u64, true);
                    }
                }
                if let Some(size) = tag_value_size(val_type, &aux_data[p + 3..]) {
                    p += 3 + size;
                } else {
                    break;
                }
            } else if let Some(size) = tag_value_size(val_type, &aux_data[p + 3..]) {
                p += 3 + size;
            } else {
                break;
            }
            found |= 1;
            if found & target_bits == target_bits {
                return result;
            }
            continue;
        }

        if val_type == b'Z' {
            let start = p + 3;
            if let Some(end) = nul_offset(&aux_data[start..]) {
                let value = &aux_data[start..start + end];
                if t == SamTag::RG {
                    result.rg = Some(value);
                    found |= 2;
                }
                if cell_tag.is_some_and(|ct| t == *ct) {
                    result.cell = Some(value);
                    found |= 4;
                }
                if t == SamTag::MC {
                    result.mc = Some(value);
                    found |= 8;
                }
                if found & target_bits == target_bits {
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
#[derive(Copy, Clone, Debug, PartialEq)]
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
/// Returns `None` if the tag is absent or is not of type `B`.
#[inline]
#[must_use]
pub fn find_array_tag(aux_data: &[u8], tag: impl AsTagBytes) -> Option<ArrayTagRef<'_>> {
    let tag = tag.as_tag_bytes();
    let (p, val_type) = find_tag_position(aux_data, *tag)?;
    if val_type != b'B' {
        return None;
    }
    parse_array_tag_at(aux_data, p + 3)
}

/// Parse B-type array tag data starting at the given offset (after the type byte).
///
/// Shared by [`find_array_tag`] and [`reverse_array_tag_in_place`].
fn parse_array_tag_at(aux_data: &[u8], data_start: usize) -> Option<ArrayTagRef<'_>> {
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
    Some(ArrayTagRef { data: &aux_data[elements_start..elements_end], elem_type, count, elem_size })
}

/// Read the `tc` template-coordinate tag from aux data.
///
/// `fgumi zipper` stamps the exact template coordinate of the primary pair onto
/// secondary/supplementary reads as a 6-element `B:i` array
/// `[tid1, pos1, neg1, tid2, pos2, neg2]` (canonically ordered `tid1 <= tid2`).
/// This lets the template-coordinate sort and dedup place a supplementary read at
/// its primary pair's coordinate — which a supplementary record cannot otherwise
/// reconstruct, since it carries its own and its mate's position but not its own
/// primary's.
///
/// Returns `None` when the tag is absent or is not a 6-element `B:i` array.
#[must_use]
pub fn read_tc_template_coordinate(aux_data: &[u8]) -> Option<[i32; 6]> {
    let arr = find_array_tag(aux_data, [b't', b'c'])?;
    if arr.elem_type != b'i' || arr.count != 6 {
        return None;
    }
    let mut out = [0i32; 6];
    for (i, slot) in out.iter_mut().enumerate() {
        let off = i * 4;
        *slot = i32::from_le_bytes([
            arr.data[off],
            arr.data[off + 1],
            arr.data[off + 2],
            arr.data[off + 3],
        ]);
    }
    Some(out)
}

/// Zero-copy typed BAM aux tag value.
///
/// Borrows directly into the raw record bytes — no allocation, no decode beyond the tag header.
/// Integer variants are always widened to `i64` regardless of the on-disk type byte.
#[derive(Copy, Clone, Debug, PartialEq)]
pub enum TagValue<'a> {
    /// `A` — single printable ASCII character.
    Char(u8),
    /// `c`/`C`/`s`/`S`/`i`/`I` — any integer type, widened to `i64`.
    Int(i64),
    /// `f` — 32-bit IEEE float.
    Float(f32),
    /// `Z` — NUL-terminated string; bytes exclude the NUL.
    String(&'a [u8]),
    /// `H` — NUL-terminated hex string; bytes exclude the NUL.
    Hex(&'a [u8]),
    /// `B` — typed array; zero-allocation reference into the record.
    Array(ArrayTagRef<'a>),
}

impl TagValue<'_> {
    /// Returns the canonical BAM aux type byte for this variant.
    ///
    /// Integer variants always return `b'i'` (the widened representation);
    /// callers that need the exact on-disk width should use [`RawTagsView::iter`] directly.
    /// Array variants always return `b'B'` — the element subtype is available via
    /// [`ArrayTagRef::elem_type`].
    #[must_use]
    pub fn type_byte(&self) -> u8 {
        match self {
            Self::Char(_) => b'A',
            Self::Int(_) => b'i',
            Self::Float(_) => b'f',
            Self::String(_) => b'Z',
            Self::Hex(_) => b'H',
            Self::Array(_) => b'B',
        }
    }
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

/// Append a pre-encoded aux tag entry to a BAM record, verbatim.
///
/// Writes `[tag_byte_1, tag_byte_2, type_byte, value_bytes...]` to the end of the
/// record's aux section. `type_byte` and `value_bytes` are the exact on-disk BAM
/// encoding of the value (e.g. a `Z` value's trailing NUL, or a `B` array's
/// subtype + count + elements), as yielded by [`TagEntry`]. The value is copied
/// with no decode/re-encode, so the destination entry is byte-identical to the
/// source apart from the two-byte tag identifier.
///
/// This is the shared primitive behind copying/renaming a tag: pair it with
/// [`remove_tag`] on the destination first if an existing value must be
/// overwritten. The caller is responsible for supplying a well-formed
/// `(type_byte, value_bytes)`; no validation is performed.
#[inline]
pub fn append_raw_tag(
    record: &mut Vec<u8>,
    tag: impl AsTagBytes,
    type_byte: u8,
    value_bytes: &[u8],
) {
    let tag = tag.as_tag_bytes();
    record.push(tag[0]);
    record.push(tag[1]);
    record.push(type_byte);
    record.extend_from_slice(value_bytes);
}

/// Append a string (Z-type) tag to a BAM record.
///
/// The tag is appended at the end of the record: `[tag_byte_1, tag_byte_2, 'Z', value..., NUL]`.
#[inline]
pub(crate) fn append_string_tag(record: &mut Vec<u8>, tag: impl AsTagBytes, value: &[u8]) {
    let tag = tag.as_tag_bytes();
    record.push(tag[0]);
    record.push(tag[1]);
    record.push(b'Z');
    record.extend_from_slice(value);
    record.push(0); // null terminator
}

/// Append an integer tag using the smallest type that fits the value.
///
/// For non-negative values, unsigned types are preferred over signed types to minimize the number
/// of bytes written, matching the behavior of htsjdk and samtools.
///
/// Encodes as:
/// - `i8`  (type `'c'`): if value in `[-128, 127]`
/// - `u8`  (type `'C'`): if value in `[128, 255]`
/// - `i16` (type `'s'`): if value in `[-32768, -129]`
/// - `u16` (type `'S'`): if value in `[256, 65535]`
/// - `i32` (type `'i'`): otherwise
#[inline]
pub(crate) fn append_int_tag(record: &mut Vec<u8>, tag: impl AsTagBytes, value: i32) {
    let tag = tag.as_tag_bytes();
    record.push(tag[0]);
    record.push(tag[1]);
    if let Ok(v) = i8::try_from(value) {
        record.push(b'c');
        record.push(v.cast_unsigned());
    } else if let Ok(v) = u8::try_from(value) {
        record.push(b'C');
        record.push(v);
    } else if let Ok(v) = u16::try_from(value) {
        record.push(b'S');
        record.extend_from_slice(&v.to_le_bytes());
    } else if let Ok(v) = i16::try_from(value) {
        record.push(b's');
        record.extend_from_slice(&v.to_le_bytes());
    } else {
        record.push(b'i');
        record.extend_from_slice(&value.to_le_bytes());
    }
}

/// Append a float (`f`-type) tag to a BAM record.
#[inline]
pub(crate) fn append_float_tag(record: &mut Vec<u8>, tag: impl AsTagBytes, value: f32) {
    let tag = tag.as_tag_bytes();
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
#[inline]
pub(crate) fn append_i16_array_tag(record: &mut Vec<u8>, tag: impl AsTagBytes, values: &[i16]) {
    let tag = tag.as_tag_bytes();
    record.push(tag[0]);
    record.push(tag[1]);
    record.push(b'B');
    record.push(b's');
    record.extend_from_slice(
        &u32::try_from(values.len()).expect("array length exceeds u32").to_le_bytes(),
    );
    for &v in values {
        record.extend_from_slice(&v.to_le_bytes());
    }
}

/// Append a `u8` array (`B:C`-type) tag to a BAM record.
///
/// Format: `[tag0, tag1, 'B', 'C', count_u32_le, values_u8...]`
///
/// # Panics
///
/// Panics if `values.len()` exceeds `u32::MAX`.
#[inline]
pub(crate) fn append_u8_array_tag(record: &mut Vec<u8>, tag: impl AsTagBytes, values: &[u8]) {
    let tag = tag.as_tag_bytes();
    record.push(tag[0]);
    record.push(tag[1]);
    record.push(b'B');
    record.push(b'C');
    record.extend_from_slice(
        &u32::try_from(values.len()).expect("array length exceeds u32").to_le_bytes(),
    );
    record.extend_from_slice(values);
}

/// Append an `i8` array (`B:c`-type) tag to a BAM record.
///
/// Format: `[tag0, tag1, 'B', 'c', count_u32_le, values_i8...]`
///
/// # Panics
///
/// Panics if `values.len()` exceeds `u32::MAX`.
#[inline]
pub(crate) fn append_i8_array_tag(record: &mut Vec<u8>, tag: impl AsTagBytes, values: &[i8]) {
    let tag = tag.as_tag_bytes();
    record.push(tag[0]);
    record.push(tag[1]);
    record.push(b'B');
    record.push(b'c');
    record.extend_from_slice(
        &u32::try_from(values.len()).expect("array length exceeds u32").to_le_bytes(),
    );
    for &v in values {
        record.push(v.cast_unsigned());
    }
}

/// Append a `u32` array (`B:I`-type) tag to a BAM record.
///
/// Format: `[tag0, tag1, 'B', 'I', count_u32_le, values_u32_le...]`
///
/// # Panics
///
/// Panics if `values.len()` exceeds `u32::MAX`.
#[inline]
pub(crate) fn append_u32_array_tag(record: &mut Vec<u8>, tag: impl AsTagBytes, values: &[u32]) {
    let tag = tag.as_tag_bytes();
    record.push(tag[0]);
    record.push(tag[1]);
    record.push(b'B');
    record.push(b'I');
    record.extend_from_slice(
        &u32::try_from(values.len()).expect("array length exceeds u32").to_le_bytes(),
    );
    for &v in values {
        record.extend_from_slice(&v.to_le_bytes());
    }
}

/// Append a Phred+33 encoded quality string (`Z`-type) tag.
///
/// Converts raw Phred scores (0-93) to ASCII (Phred+33) and writes
/// directly as a null-terminated string tag. Avoids intermediate String allocation.
#[inline]
pub(crate) fn append_phred33_string_tag(record: &mut Vec<u8>, tag: impl AsTagBytes, quals: &[u8]) {
    let tag = tag.as_tag_bytes();
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
#[inline]
pub fn remove_tag(record: &mut Vec<u8>, tag: impl AsTagBytes) {
    let tag = tag.as_tag_bytes();
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
#[inline]
pub fn update_string_tag(record: &mut Vec<u8>, tag: impl AsTagBytes, new_value: &[u8]) {
    let tag = tag.as_tag_bytes();
    let aux_start = aux_data_offset_from_record(record).unwrap_or(record.len());
    if aux_start < record.len()
        && let Some((start, end)) = find_tag_bounds(&record[aux_start..], tag)
    {
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
    // Tag not found — append
    append_string_tag(record, tag, new_value);
}

/// Update an existing int tag in-place, or append if absent.
///
/// If the existing tag is a 4-byte `i`/`I` type, the value is overwritten in-place.
/// Otherwise the tag is removed and re-appended using `append_int_tag`
/// (smallest type that fits).
#[inline]
pub fn update_int_tag(record: &mut Vec<u8>, tag: impl AsTagBytes, value: i32) {
    let tag = tag.as_tag_bytes();
    let aux_start = aux_data_offset_from_record(record).unwrap_or(record.len());
    if aux_start < record.len()
        && let Some((start, end)) = find_tag_bounds(&record[aux_start..], tag)
    {
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
    // Tag not found — append
    append_int_tag(record, tag, value);
}

/// Reverse elements of a B-type array tag in place. No-op if tag not found.
#[inline]
pub fn reverse_array_tag_in_place(record: &mut [u8], aux_offset: usize, tag: impl AsTagBytes) {
    let tag = tag.as_tag_bytes();
    if aux_offset >= record.len() {
        return;
    }
    let Some((p, val_type)) = find_tag_position(&record[aux_offset..], *tag) else {
        return;
    };
    if val_type != b'B' {
        return;
    }
    // Extract array metadata before taking mutable borrow
    let Some(arr) = parse_array_tag_at(&record[aux_offset..], p + 3) else {
        return;
    };
    let count = arr.count;
    let elem_size = arr.elem_size;
    if count == 0 {
        return;
    }
    let elements_start = aux_offset + p + 3 + 5;
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
}

/// Find the byte range of a Z-type string tag value within a mutable record.
///
/// Returns `Some((start, end))` where `start..end` is the value range (excluding NUL).
/// Returns `None` if the tag is not found or is not Z-type.
fn find_string_tag_range(record: &[u8], aux_offset: usize, tag: [u8; 2]) -> Option<(usize, usize)> {
    if aux_offset >= record.len() {
        return None;
    }
    let (p, val_type) = find_tag_position(&record[aux_offset..], tag)?;
    if val_type != b'Z' {
        return None;
    }
    let start = aux_offset + p + 3;
    let nul_off = nul_offset(&record[start..])?;
    let end = start + nul_off;
    if end > start { Some((start, end)) } else { None }
}

/// Reverse bytes of a Z-type string tag value in place. No-op if tag not found.
#[inline]
pub fn reverse_string_tag_in_place(record: &mut [u8], aux_offset: usize, tag: impl AsTagBytes) {
    let tag = tag.as_tag_bytes();
    if let Some((start, end)) = find_string_tag_range(record, aux_offset, *tag) {
        record[start..end].reverse();
    }
}

/// Reverse-complement a Z-type string tag value in place.
///
/// Uses `fgumi_dna::COMPLEMENT`, the shared IUPAC-aware, case-preserving
/// complement table (unknown bytes pass through), so tag reverse-complement
/// matches the consensus base paths and fgbio/htsjdk (`R<->Y`, etc.).
#[inline]
pub fn reverse_complement_string_tag_in_place(
    record: &mut [u8],
    aux_offset: usize,
    tag: impl AsTagBytes,
) {
    let tag = tag.as_tag_bytes();
    if let Some((start, end)) = find_string_tag_range(record, aux_offset, *tag) {
        record[start..end].reverse();
        for b in &mut record[start..end] {
            *b = fgumi_dna::COMPLEMENT[*b as usize];
        }
    }
}

/// Append an `i32` array (`B:i`-type) tag to a BAM record.
///
/// Format: `[tag0, tag1, 'B', 'i', count_u32_le, values_i32_le...]`
///
/// # Panics
///
/// Panics if `values.len()` exceeds `u32::MAX`.
#[inline]
pub fn append_i32_array_tag(record: &mut Vec<u8>, tag: impl AsTagBytes, values: &[i32]) {
    let tag = tag.as_tag_bytes();
    record.push(tag[0]);
    record.push(tag[1]);
    record.push(b'B');
    record.push(b'i');
    record.extend_from_slice(
        &u32::try_from(values.len()).expect("array length exceeds u32").to_le_bytes(),
    );
    for &v in values {
        record.extend_from_slice(&v.to_le_bytes());
    }
}

/// Append a `u16` array (`B:S`-type) tag to a BAM record.
///
/// Format: `[tag0, tag1, 'B', 'S', count_u32_le, values_u16_le...]`
///
/// # Panics
///
/// Panics if `values.len()` exceeds `u32::MAX`.
#[inline]
pub(crate) fn append_u16_array_tag(record: &mut Vec<u8>, tag: impl AsTagBytes, values: &[u16]) {
    let tag = tag.as_tag_bytes();
    record.push(tag[0]);
    record.push(tag[1]);
    record.push(b'B');
    record.push(b'S');
    record.extend_from_slice(
        &u32::try_from(values.len()).expect("array length exceeds u32").to_le_bytes(),
    );
    for &v in values {
        record.extend_from_slice(&v.to_le_bytes());
    }
}

/// Append an `f32` array (`B:f`-type) tag to a BAM record.
///
/// Format: `[tag0, tag1, 'B', 'f', count_u32_le, values_f32_le...]`
///
/// # Panics
///
/// Panics if `values.len()` exceeds `u32::MAX`.
#[inline]
pub(crate) fn append_f32_array_tag(record: &mut Vec<u8>, tag: impl AsTagBytes, values: &[f32]) {
    let tag = tag.as_tag_bytes();
    record.push(tag[0]);
    record.push(tag[1]);
    record.push(b'B');
    record.push(b'f');
    record.extend_from_slice(
        &u32::try_from(values.len()).expect("array length exceeds u32").to_le_bytes(),
    );
    for &v in values {
        record.extend_from_slice(&v.to_le_bytes());
    }
}

/// Normalize an integer tag to the smallest signed type that fits its value.
///
/// If the tag exists and is an integer type, it is re-encoded using
/// [`append_signed_int_tag`] semantics (signed types only: i8 → i16 → i32).
/// No-op if the tag is not found or is not an integer type.
#[inline]
pub fn normalize_int_tag_to_smallest_signed(record: &mut Vec<u8>, tag: impl AsTagBytes) {
    let tag = tag.as_tag_bytes();
    let Some(aux_start) = aux_data_offset_from_record(record) else {
        return;
    };
    if aux_start >= record.len() {
        return;
    }
    let aux_data = &record[aux_start..];
    let Some(value) = find_int_tag(aux_data, tag) else {
        return;
    };
    let Ok(value_i32) = i32::try_from(value) else {
        return;
    };
    // Remove and re-append with smallest *signed* encoding
    // (i8 -> i16 -> i32, matching fgbio's to_smallest_signed_int)
    remove_tag(record, tag);
    append_signed_int_tag(record, tag, value_i32);
}

/// Append an integer tag using the smallest *signed* type that fits.
///
/// Unlike `append_int_tag` (which prefers unsigned types), this uses only
/// signed types: `i8` (type `'c'`) -> `i16` (type `'s'`) -> `i32` (type `'i'`).
/// This matches fgbio's `to_smallest_signed_int` encoding.
#[inline]
pub fn append_signed_int_tag(record: &mut Vec<u8>, tag: impl AsTagBytes, value: i32) {
    let tag = tag.as_tag_bytes();
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

/// Copy auxiliary tags from source aux data to a destination record, skipping specified tags.
///
/// Iterates all tags in `src_aux` and appends each tag entry (tag + type + value bytes)
/// to `dest`, unless the tag's two-byte key is in `skip_tags`.
pub(crate) fn copy_aux_tags(src_aux: &[u8], dest: &mut Vec<u8>, skip_tags: &[SamTag]) {
    let mut offset = 0;
    while offset + 3 <= src_aux.len() {
        let tag_key = [src_aux[offset], src_aux[offset + 1]];
        let val_type = src_aux[offset + 2];
        let value_start = offset + 3;

        let Some(value_size) = tag_value_size(val_type, &src_aux[value_start..]) else {
            break;
        };
        let entry_end = value_start + value_size;
        if entry_end > src_aux.len() {
            break;
        }

        // Copy unless this tag should be skipped
        if !skip_tags.iter().any(|t| *t == tag_key) {
            dest.extend_from_slice(&src_aux[offset..entry_end]);
        }

        offset = entry_end;
    }
}

/// Zero-allocation entry yielded by tag iterators.
#[derive(Copy, Clone, Debug)]
pub struct TagEntry<'a> {
    pub tag: [u8; 2],
    pub type_byte: u8,
    pub value_bytes: &'a [u8],
}

/// Borrowed read-only view over a BAM record's auxiliary tag section.
#[derive(Copy, Clone, Debug)]
pub struct RawTagsView<'a>(&'a [u8]);

impl<'a> RawTagsView<'a> {
    /// Wraps a raw auxiliary-data slice.
    #[inline]
    #[must_use]
    pub const fn new(aux: &'a [u8]) -> Self {
        Self(aux)
    }

    /// Returns the underlying aux-data bytes.
    #[inline]
    #[must_use]
    pub const fn as_bytes(&self) -> &'a [u8] {
        self.0
    }

    /// Returns the number of bytes in the aux section.
    #[inline]
    #[must_use]
    pub const fn len(&self) -> usize {
        self.0.len()
    }

    /// Returns `true` if the aux section is empty (no tags present).
    #[inline]
    #[must_use]
    pub const fn is_empty(&self) -> bool {
        self.0.is_empty()
    }

    /// Returns `true` if the named tag is present in the aux section.
    #[inline]
    #[must_use]
    pub fn contains(&self, tag: impl AsTagBytes) -> bool {
        find_tag_type(self.0, tag).is_some()
    }

    /// Returns the value of a `Z`-type (string) tag, without the null terminator.
    ///
    /// Returns `None` if the tag is absent or is not `Z`-typed.
    #[inline]
    #[must_use]
    pub fn find_string(&self, tag: impl AsTagBytes) -> Option<&'a [u8]> {
        find_string_tag(self.0, tag)
    }

    /// Returns the value of any integer tag (`c/C/s/S/i/I`) widened to `i64`.
    ///
    /// Returns `None` if the tag is absent or is not an integer type.
    #[inline]
    #[must_use]
    pub fn find_int(&self, tag: impl AsTagBytes) -> Option<i64> {
        find_int_tag(self.0, tag)
    }

    /// Returns the value of a `C`-type (unsigned byte) tag.
    ///
    /// Returns `None` if the tag is absent or is not `C`-typed.
    #[inline]
    #[must_use]
    pub fn find_uint8(&self, tag: impl AsTagBytes) -> Option<u8> {
        find_uint8_tag(self.0, tag)
    }

    /// Returns the value of an `f`-type (32-bit float) tag.
    ///
    /// Returns `None` if the tag is absent or is not `f`-typed.
    #[inline]
    #[must_use]
    pub fn find_float(&self, tag: impl AsTagBytes) -> Option<f32> {
        find_float_tag(self.0, tag)
    }

    /// Returns a zero-allocation reference to a `B`-type (array) tag.
    ///
    /// Returns `None` if the tag is absent or is not `B`-typed.
    #[inline]
    #[must_use]
    pub fn find_array(&self, tag: impl AsTagBytes) -> Option<ArrayTagRef<'a>> {
        find_array_tag(self.0, tag)
    }

    /// Returns the MI (Molecular Identifier) tag as `(value, is_A_suffix)`.
    ///
    /// Returns `None` if the tag is absent or cannot be parsed.
    #[inline]
    #[must_use]
    pub fn find_mi(&self) -> Option<(u64, bool)> {
        find_mi_tag(self.0)
    }

    /// Returns the MC (mate CIGAR) tag as raw CIGAR bytes.
    ///
    /// Returns `None` if the tag is absent or is not `Z`-typed.
    #[inline]
    #[must_use]
    pub fn find_mc(&self) -> Option<&'a [u8]> {
        find_mc_tag(self.0)
    }

    /// Returns the typed value of a tag as a [`TagValue`], borrowing into the record bytes.
    ///
    /// Handles all BAM aux types: `A`, `c`/`C`/`s`/`S`/`i`/`I`, `f`, `Z`, `H`, `B`.
    /// Integer types are widened to `i64`. Returns `None` if the tag is absent or malformed.
    #[inline]
    #[must_use]
    pub fn get(&self, tag: impl AsTagBytes) -> Option<TagValue<'a>> {
        let tag = tag.as_tag_bytes();
        let (p, val_type) = find_tag_position(self.0, *tag)?;
        let start = p + 3;
        let aux = self.0;
        match val_type {
            b'A' => aux.get(start).copied().map(TagValue::Char),
            b'c' | b'C' | b's' | b'S' | b'i' | b'I' => {
                extract_int_value(aux, p, val_type).map(TagValue::Int)
            }
            b'f' => {
                if start + 4 > aux.len() {
                    return None;
                }
                Some(TagValue::Float(f32::from_le_bytes([
                    aux[start],
                    aux[start + 1],
                    aux[start + 2],
                    aux[start + 3],
                ])))
            }
            b'Z' => {
                let end = nul_offset(&aux[start..])?;
                Some(TagValue::String(&aux[start..start + end]))
            }
            b'H' => {
                let end = nul_offset(&aux[start..])?;
                Some(TagValue::Hex(&aux[start..start + end]))
            }
            b'B' => parse_array_tag_at(aux, start).map(TagValue::Array),
            _ => None,
        }
    }

    /// Iterate over all tags as `(tag, TagValue)` pairs, borrowing into the record bytes.
    ///
    /// Skips any tag whose value cannot be decoded (malformed entries stop iteration).
    /// Decodes each entry directly from the iterator's yielded bytes, so
    /// iteration is O(n) in the aux size and duplicate tag keys each yield their
    /// own value.
    pub fn iter_typed(&self) -> impl Iterator<Item = ([u8; 2], TagValue<'a>)> + 'a {
        RawTagsView::new(self.0)
            .iter()
            .filter_map(|entry| decode_tag_entry(&entry).map(|v| (entry.tag, v)))
    }
}

/// Decode a [`TagEntry`] into a [`TagValue`], borrowing into the same aux bytes
/// that backed the entry. Returns `None` for malformed payloads. Used by
/// [`RawTagsView::iter_typed`] to avoid rescanning.
fn decode_tag_entry<'a>(entry: &TagEntry<'a>) -> Option<TagValue<'a>> {
    let bytes = entry.value_bytes;
    match entry.type_byte {
        b'A' => bytes.first().copied().map(TagValue::Char),
        b'c' => bytes.first().map(|&b| TagValue::Int(i64::from(b.cast_signed()))),
        b'C' => bytes.first().map(|&b| TagValue::Int(i64::from(b))),
        b's' if bytes.len() >= 2 => {
            Some(TagValue::Int(i64::from(i16::from_le_bytes([bytes[0], bytes[1]]))))
        }
        b'S' if bytes.len() >= 2 => {
            Some(TagValue::Int(i64::from(u16::from_le_bytes([bytes[0], bytes[1]]))))
        }
        b'i' if bytes.len() >= 4 => Some(TagValue::Int(i64::from(i32::from_le_bytes([
            bytes[0], bytes[1], bytes[2], bytes[3],
        ])))),
        b'I' if bytes.len() >= 4 => Some(TagValue::Int(i64::from(u32::from_le_bytes([
            bytes[0], bytes[1], bytes[2], bytes[3],
        ])))),
        b'f' if bytes.len() >= 4 => {
            Some(TagValue::Float(f32::from_le_bytes([bytes[0], bytes[1], bytes[2], bytes[3]])))
        }
        // Z/H payloads include the trailing NUL; strip it for TagValue.
        b'Z' => bytes
            .split_last()
            .and_then(|(&last, rest)| (last == 0).then_some(rest))
            .map(TagValue::String),
        b'H' => bytes
            .split_last()
            .and_then(|(&last, rest)| (last == 0).then_some(rest))
            .map(TagValue::Hex),
        b'B' => parse_array_tag_at(bytes, 0).map(TagValue::Array),
        _ => None,
    }
}

/// Forward-only iterator over aux tag entries.
///
/// Yields `TagEntry` values containing the 2-byte tag, the type byte, and
/// the raw value bytes (whose length is determined by the type per BAM spec).
/// Stops at the first malformed entry.
pub struct AuxTagsIter<'a> {
    aux: &'a [u8],
    pos: usize,
}

impl<'a> Iterator for AuxTagsIter<'a> {
    type Item = TagEntry<'a>;
    fn next(&mut self) -> Option<Self::Item> {
        if self.pos + 3 > self.aux.len() {
            return None;
        }
        let tag = [self.aux[self.pos], self.aux[self.pos + 1]];
        let type_byte = self.aux[self.pos + 2];
        let value_start = self.pos + 3;
        let size = tag_value_size(type_byte, &self.aux[value_start..])?;
        let end = value_start.checked_add(size)?;
        if end > self.aux.len() {
            return None;
        }
        let entry = TagEntry { tag, type_byte, value_bytes: &self.aux[value_start..end] };
        self.pos = end;
        Some(entry)
    }
}

impl<'a> RawTagsView<'a> {
    /// Iterate over aux tag entries without allocating.
    #[inline]
    #[must_use]
    pub fn iter(&self) -> AuxTagsIter<'a> {
        AuxTagsIter { aux: self.0, pos: 0 }
    }
}

impl<'a> IntoIterator for &RawTagsView<'a> {
    type Item = TagEntry<'a>;
    type IntoIter = AuxTagsIter<'a>;
    fn into_iter(self) -> Self::IntoIter {
        self.iter()
    }
}

impl<'a> RawTagsView<'a> {
    /// Single-pass extraction of (RG, cell, MC) string tags from this aux
    /// section, optionally also capturing the UMI tag's position.
    #[inline]
    #[must_use]
    pub fn extract_string_batch(
        &self,
        cell_tag: impl AsTagBytes,
        umi_tag: Option<[u8; 2]>,
    ) -> AuxStringTags<'a> {
        extract_aux_string_tags(self.0, cell_tag, umi_tag)
    }
}

impl<'a> RawRecordView<'a> {
    /// Returns a read-only view over this record's auxiliary tag section.
    #[inline]
    #[must_use]
    pub fn tags(&self) -> RawTagsView<'a> {
        RawTagsView::new(aux_data_slice(self.as_bytes()))
    }

    /// Single-pass extraction of (MI, RG, cell, MC) from this record's aux data.
    #[inline]
    #[must_use]
    pub fn template_aux_tags(&self, cell_tag: Option<SamTag>) -> TemplateAuxTags<'a> {
        extract_template_aux_tags(self.as_bytes(), cell_tag)
    }
}

/// Borrowed fixed-length mutable view over a BAM record's aux section.
///
/// Cannot change the byte count of the aux section. Hosts in-place writers
/// for B-array elements, Z-tag byte flips, and same-type tag overwrites.
pub struct RawTagsMut<'a>(&'a mut [u8]);

impl<'a> RawTagsMut<'a> {
    /// Wraps a raw auxiliary-data slice mutably.
    #[inline]
    #[must_use]
    pub fn new(aux: &'a mut [u8]) -> Self {
        Self(aux)
    }

    /// Returns a read-only view over the aux section.
    #[inline]
    #[must_use]
    pub fn view(&self) -> RawTagsView<'_> {
        RawTagsView::new(self.0)
    }

    /// Overwrite a single u8 element of a B:C array tag in place.
    ///
    /// No-op on tag absence, type mismatch (must be `B` with sub-type `C`), or out-of-range index.
    #[inline]
    pub fn set_array_element_u8(&mut self, tag: impl AsTagBytes, index: usize, value: u8) {
        Self::set_array_element_le(self.0, *tag.as_tag_bytes(), b'C', 1, index, &[value]);
    }

    /// Overwrite a single u16 element of a B:S array tag in place.
    ///
    /// No-op on tag absence, type mismatch (must be `B` with sub-type `S`), or out-of-range index.
    #[inline]
    pub fn set_array_element_u16(&mut self, tag: impl AsTagBytes, index: usize, value: u16) {
        Self::set_array_element_le(
            self.0,
            *tag.as_tag_bytes(),
            b'S',
            2,
            index,
            &value.to_le_bytes(),
        );
    }

    /// Overwrite a single i16 element of a B:s array tag in place.
    ///
    /// No-op on tag absence, type mismatch (must be `B` with sub-type `s`), or out-of-range index.
    #[inline]
    pub fn set_array_element_i16(&mut self, tag: impl AsTagBytes, index: usize, value: i16) {
        Self::set_array_element_le(
            self.0,
            *tag.as_tag_bytes(),
            b's',
            2,
            index,
            &value.to_le_bytes(),
        );
    }

    /// Overwrite a single i32 element of a B:i array tag in place.
    ///
    /// No-op on tag absence, type mismatch (must be `B` with sub-type `i`), or out-of-range index.
    #[inline]
    pub fn set_array_element_i32(&mut self, tag: impl AsTagBytes, index: usize, value: i32) {
        Self::set_array_element_le(
            self.0,
            *tag.as_tag_bytes(),
            b'i',
            4,
            index,
            &value.to_le_bytes(),
        );
    }

    /// Overwrite a single f32 element of a B:f array tag in place.
    ///
    /// No-op on tag absence, type mismatch (must be `B` with sub-type `f`), or out-of-range index.
    #[inline]
    pub fn set_array_element_f32(&mut self, tag: impl AsTagBytes, index: usize, value: f32) {
        Self::set_array_element_le(
            self.0,
            *tag.as_tag_bytes(),
            b'f',
            4,
            index,
            &value.to_le_bytes(),
        );
    }

    /// Internal helper: locate a B-type array of the expected sub-type and overwrite
    /// element `index` with the provided little-endian bytes.
    ///
    /// The layout of a B-type tag entry in the aux section is:
    /// `[tag0, tag1, 'B', sub_type, count(4 bytes LE), elem0..elemN]`
    fn set_array_element_le(
        aux: &mut [u8],
        tag: [u8; 2],
        expected_sub: u8,
        elem_size: usize,
        index: usize,
        le_bytes: &[u8],
    ) {
        let Some((p, val_type)) = find_tag_position(aux, tag) else { return };
        if val_type != b'B' || p + 8 > aux.len() {
            return;
        }
        if aux[p + 3] != expected_sub {
            return;
        }
        let count = u32::from_le_bytes([aux[p + 4], aux[p + 5], aux[p + 6], aux[p + 7]]) as usize;
        if index >= count {
            return;
        }
        let off = p + 8 + index * elem_size;
        if off + elem_size > aux.len() {
            return;
        }
        aux[off..off + elem_size].copy_from_slice(le_bytes);
    }

    /// Reverse element bytes of a B-type array tag in place.
    #[inline]
    pub fn reverse_array(&mut self, tag: impl AsTagBytes) {
        reverse_array_tag_in_place(self.0, 0, tag);
    }

    /// Reverse the bytes of a Z-type string tag value in place.
    #[inline]
    pub fn reverse_string(&mut self, tag: impl AsTagBytes) {
        reverse_string_tag_in_place(self.0, 0, tag);
    }

    /// Reverse-complement a Z-type string tag value (A<->T, C<->G).
    #[inline]
    pub fn reverse_complement_string(&mut self, tag: impl AsTagBytes) {
        reverse_complement_string_tag_in_place(self.0, 0, tag);
    }

    /// Overwrite an existing Z-type tag if the new value matches its current length.
    ///
    /// Returns `false` (no-op) if absent, wrong type, or different length.
    pub fn set_string_in_place(&mut self, tag: impl AsTagBytes, value: &[u8]) -> bool {
        let tag = tag.as_tag_bytes();
        let Some((p, val_type)) = find_tag_position(self.0, *tag) else {
            return false;
        };
        if val_type != b'Z' {
            return false;
        }
        let start = p + 3;
        let Some(nul_off) = nul_offset(&self.0[start..]) else {
            return false;
        };
        if nul_off != value.len() {
            return false;
        }
        self.0[start..start + value.len()].copy_from_slice(value);
        true
    }

    /// Overwrite an existing integer tag if the new value fits its current type byte.
    ///
    /// Returns `false` if absent or value doesn't fit the existing type.
    pub fn set_int_in_place(&mut self, tag: impl AsTagBytes, value: i64) -> bool {
        let tag = tag.as_tag_bytes();
        let Some((p, val_type)) = find_tag_position(self.0, *tag) else {
            return false;
        };
        // `find_tag_position` short-circuits on the first matching key, so the
        // payload hasn't been length-validated. Guard every write: malformed
        // / truncated aux bytes must return false, not panic.
        let (need, data_start) = match val_type {
            b'c' | b'C' => (1, p + 3),
            b's' | b'S' => (2, p + 3),
            b'i' | b'I' => (4, p + 3),
            _ => return false,
        };
        let Some(target) = self.0.get_mut(data_start..data_start + need) else {
            return false;
        };
        match val_type {
            b'c' => {
                if let Ok(v) = i8::try_from(value) {
                    target[0] = v.cast_unsigned();
                    true
                } else {
                    false
                }
            }
            b'C' => {
                if let Ok(v) = u8::try_from(value) {
                    target[0] = v;
                    true
                } else {
                    false
                }
            }
            b's' => {
                if let Ok(v) = i16::try_from(value) {
                    target.copy_from_slice(&v.to_le_bytes());
                    true
                } else {
                    false
                }
            }
            b'S' => {
                if let Ok(v) = u16::try_from(value) {
                    target.copy_from_slice(&v.to_le_bytes());
                    true
                } else {
                    false
                }
            }
            b'i' => {
                if let Ok(v) = i32::try_from(value) {
                    target.copy_from_slice(&v.to_le_bytes());
                    true
                } else {
                    false
                }
            }
            b'I' => {
                if let Ok(v) = u32::try_from(value) {
                    target.copy_from_slice(&v.to_le_bytes());
                    true
                } else {
                    false
                }
            }
            _ => unreachable!("val_type was filtered by the match above"),
        }
    }

    /// Overwrite an existing `f`-type tag in place.
    ///
    /// Returns `false` if absent, wrong type, or the payload is truncated.
    #[inline]
    pub fn set_float_in_place(&mut self, tag: impl AsTagBytes, value: f32) -> bool {
        let tag = tag.as_tag_bytes();
        let Some((p, val_type)) = find_tag_position(self.0, *tag) else {
            return false;
        };
        if val_type != b'f' {
            return false;
        }
        // Guard against truncated payloads (see `set_int_in_place` for context).
        let Some(target) = self.0.get_mut(p + 3..p + 7) else {
            return false;
        };
        target.copy_from_slice(&value.to_le_bytes());
        true
    }
}

impl RawRecordMut<'_> {
    /// Borrow a fixed-length mutable view over this record's aux section.
    #[inline]
    #[must_use]
    pub fn tags_mut(&mut self) -> RawTagsMut<'_> {
        // Compute offset and length via immutable borrow; both immutable borrows end
        // before the mutable borrow on the last line (NLL two-phase borrows).
        let len = self.as_bytes().len();
        let off = aux_data_offset_from_record(self.as_bytes()).unwrap_or(len);
        RawTagsMut::new(&mut self.as_bytes_mut()[off..])
    }
}

/// A set membership test over two-byte SAM tag keys, used by
/// [`RawTagsEditor::rebuild_with`] to decide which existing tags to drop.
///
/// Implemented for a plain `[[u8; 2]]` slice (cheap for the small fixed sets a
/// caller like `regenerate_alignment_tags_raw` removes — `NM`/`UQ`/`MD`) and
/// for [`TagBitset`] (O(1) probes for the large per-run sets a caller like
/// `zipper` builds once and reuses).
pub trait TagKeySet {
    /// Returns `true` if `tag` is in the set.
    fn contains_tag(&self, tag: [u8; 2]) -> bool;
}

impl TagKeySet for [[u8; 2]] {
    #[inline]
    fn contains_tag(&self, tag: [u8; 2]) -> bool {
        self.contains(&tag)
    }
}

impl<const N: usize> TagKeySet for [[u8; 2]; N] {
    #[inline]
    fn contains_tag(&self, tag: [u8; 2]) -> bool {
        self.contains(&tag)
    }
}

/// A membership set over two-byte SAM tag names, backed by a 256×256 bit table
/// (65536 bits in 1024 `u64` words, indexed by `(byte0 << 8) | byte1`).
///
/// A direct bit test on the two raw tag bytes — no hashing, no UTF-8 conversion.
/// Built once per run from the tag names a caller wants to match, then probed
/// once per tag per record on a hot path (e.g. `zipper`'s tag-copy loop, where
/// it measurably beats a `HashSet`). All SAM tags are exactly two bytes, so only
/// two-byte names are representable; a longer or shorter name can never equal a
/// real tag and is dropped at construction.
///
/// Membership is tested against the tag's **raw two bytes** — no UTF-8
/// conversion. This is deliberate and differs from a `HashSet<String>` keyed on
/// `str::from_utf8(bytes).unwrap_or("")`: that path aliased a non-UTF-8 tag pair
/// (reachable only from a malformed BAM) to the empty string, so it could match
/// an empty filter, whereas here the raw bytes are matched directly and a
/// non-two-byte filter is never inserted. Do not reintroduce the UTF-8
/// conversion.
#[derive(Debug, Clone)]
pub struct TagBitset {
    words: Box<[u64; 1024]>,
}

impl TagBitset {
    /// An empty set.
    #[must_use]
    pub fn new() -> Self {
        Self { words: Box::new([0u64; 1024]) }
    }

    #[inline]
    fn bit_index(tag: [u8; 2]) -> usize {
        (usize::from(tag[0]) << 8) | usize::from(tag[1])
    }

    /// Insert a two-byte tag key.
    #[inline]
    pub fn insert(&mut self, tag: [u8; 2]) {
        let i = Self::bit_index(tag);
        self.words[i >> 6] |= 1u64 << (i & 63);
    }

    /// Returns `true` if the two-byte tag key is present.
    #[inline]
    #[must_use]
    pub fn contains(&self, tag: [u8; 2]) -> bool {
        let i = Self::bit_index(tag);
        (self.words[i >> 6] >> (i & 63)) & 1 != 0
    }

    /// Returns `true` if no tag key is present (every bit clear). Cheap enough
    /// to let a caller skip a per-record operation entirely when the set is
    /// empty (e.g. `zipper` skipping its remove pass when no tags are removed).
    #[inline]
    #[must_use]
    pub fn is_empty(&self) -> bool {
        self.words.iter().all(|&w| w == 0)
    }

    /// Build a set from tag names, keeping only the exactly-two-byte names (any
    /// other length can never equal a real SAM tag and is skipped). Accepts any
    /// string-like item (`&str`, `String`, `&String`, …) so callers holding a
    /// `&[String]`, a `&[&str]`, or an iterator of either are not forced to
    /// allocate.
    #[must_use]
    pub fn from_names<S: AsRef<str>>(names: impl IntoIterator<Item = S>) -> Self {
        let mut set = Self::new();
        for name in names {
            if let [b0, b1] = *name.as_ref().as_bytes() {
                set.insert([b0, b1]);
            }
        }
        set
    }
}

impl Default for TagBitset {
    fn default() -> Self {
        Self::new()
    }
}

impl TagKeySet for TagBitset {
    #[inline]
    fn contains_tag(&self, tag: [u8; 2]) -> bool {
        self.contains(tag)
    }
}

/// Length-changing tag editor.
///
/// Borrows a full BAM record byte buffer plus the cached aux offset (so update/
/// append/remove ops don't repeat the header scan). Splices the aux section
/// in place via `Vec::splice` for length-changing operations.
pub struct RawTagsEditor<'a> {
    record: &'a mut Vec<u8>,
    aux_offset: usize,
}

impl<'a> RawTagsEditor<'a> {
    /// Borrow a full BAM record byte buffer.
    ///
    /// Computes and caches `aux_offset` once. If the record is too short for
    /// a valid header, `aux_offset` falls back to `record.len()` so all tag
    /// ops behave as if the aux section is empty.
    #[inline]
    pub fn from_vec(record: &'a mut Vec<u8>) -> Self {
        let aux_offset = aux_data_offset_from_record(record).unwrap_or(record.len());
        Self { record, aux_offset }
    }

    /// Returns the cached aux offset (start of tag section).
    #[inline]
    #[must_use]
    pub fn aux_offset(&self) -> usize {
        self.aux_offset
    }

    /// Borrow the aux section as a read-only view.
    #[inline]
    #[must_use]
    pub fn view(&self) -> RawTagsView<'_> {
        let off = self.aux_offset.min(self.record.len());
        RawTagsView::new(&self.record[off..])
    }

    // -- append (always grow) --

    /// Append a string (Z-type) tag to the record.
    #[inline]
    pub fn append_string(&mut self, tag: impl AsTagBytes, value: &[u8]) {
        append_string_tag(self.record, tag, value);
    }

    /// Append an integer tag, choosing the smallest fitting unsigned/signed type.
    #[inline]
    pub fn append_int(&mut self, tag: impl AsTagBytes, value: i32) {
        append_int_tag(self.record, tag, value);
    }

    /// Append a signed integer tag, choosing the smallest fitting signed type.
    #[inline]
    pub fn append_signed_int(&mut self, tag: impl AsTagBytes, value: i32) {
        append_signed_int_tag(self.record, tag, value);
    }

    /// Append a float (`f`-type) tag to the record.
    #[inline]
    pub fn append_float(&mut self, tag: impl AsTagBytes, value: f32) {
        append_float_tag(self.record, tag, value);
    }

    /// Append a Phred+33-encoded quality string tag.
    ///
    /// Converts raw Phred scores (0–93) to ASCII (Phred+33) and writes as a Z-type tag.
    #[inline]
    pub fn append_phred33_string(&mut self, tag: impl AsTagBytes, quals: &[u8]) {
        append_phred33_string_tag(self.record, tag, quals);
    }

    /// Append a `B:C` (u8 array) tag to the record.
    #[inline]
    pub fn append_array_u8(&mut self, tag: impl AsTagBytes, values: &[u8]) {
        append_u8_array_tag(self.record, tag, values);
    }

    /// Append a `B:s` (i16 array) tag to the record.
    #[inline]
    pub fn append_array_i16(&mut self, tag: impl AsTagBytes, values: &[i16]) {
        append_i16_array_tag(self.record, tag, values);
    }

    /// Append a `B:i` (i32 array) tag to the record.
    #[inline]
    pub fn append_array_i32(&mut self, tag: impl AsTagBytes, values: &[i32]) {
        append_i32_array_tag(self.record, tag, values);
    }

    /// Append a `B:c` (i8 array) tag to the record.
    #[inline]
    pub fn append_array_i8(&mut self, tag: impl AsTagBytes, values: &[i8]) {
        append_i8_array_tag(self.record, tag, values);
    }

    /// Append a `B:I` (u32 array) tag to the record.
    #[inline]
    pub fn append_array_u32(&mut self, tag: impl AsTagBytes, values: &[u32]) {
        append_u32_array_tag(self.record, tag, values);
    }

    // -- update / remove / normalize --

    /// Update an existing integer tag in place, or append it if absent.
    #[inline]
    pub fn update_int(&mut self, tag: impl AsTagBytes, value: i32) {
        update_int_tag(self.record, tag, value);
    }

    /// Update an existing string tag in place, or append it if absent.
    #[inline]
    pub fn update_string(&mut self, tag: impl AsTagBytes, value: &[u8]) {
        update_string_tag(self.record, tag, value);
    }

    /// Update an existing `f`-type tag in place, or remove + append if the existing tag has a
    /// different type or is absent.
    pub fn update_float(&mut self, tag: impl AsTagBytes, value: f32) {
        let tag = tag.as_tag_bytes();
        let off = self.aux_offset.min(self.record.len());
        if let Some((p, val_type)) = find_tag_position(&self.record[off..], *tag) {
            let abs = off + p;
            if val_type == b'f' && abs + 7 <= self.record.len() {
                self.record[abs + 3..abs + 7].copy_from_slice(&value.to_le_bytes());
                return;
            }
            if let Some(size) = tag_value_size(val_type, &self.record[abs + 3..]) {
                self.record.drain(abs..abs + 3 + size);
            }
        }
        append_float_tag(self.record, tag, value);
    }

    /// Update an existing `B:C` (u8 array) tag in place if the length matches, else remove and
    /// append.
    pub fn update_array_u8(&mut self, tag: impl AsTagBytes, values: &[u8]) {
        let tag = tag.as_tag_bytes();
        let off = self.aux_offset.min(self.record.len());
        if let Some((p, val_type)) = find_tag_position(&self.record[off..], *tag) {
            let abs = off + p;
            if val_type == b'B' && abs + 8 <= self.record.len() && self.record[abs + 3] == b'C' {
                let count = u32::from_le_bytes([
                    self.record[abs + 4],
                    self.record[abs + 5],
                    self.record[abs + 6],
                    self.record[abs + 7],
                ]) as usize;
                if count == values.len() {
                    let body = abs + 8;
                    self.record[body..body + values.len()].copy_from_slice(values);
                    return;
                }
            }
            if let Some(size) = tag_value_size(val_type, &self.record[abs + 3..]) {
                self.record.drain(abs..abs + 3 + size);
            }
        }
        append_u8_array_tag(self.record, tag, values);
    }

    /// Update an existing `B:S` (u16 array) tag in place if the length matches, else remove and
    /// append.
    pub fn update_array_u16(&mut self, tag: impl AsTagBytes, values: &[u16]) {
        let tag = tag.as_tag_bytes();
        let off = self.aux_offset.min(self.record.len());
        if let Some((p, val_type)) = find_tag_position(&self.record[off..], *tag) {
            let abs = off + p;
            if val_type == b'B' && abs + 8 <= self.record.len() && self.record[abs + 3] == b'S' {
                let count = u32::from_le_bytes([
                    self.record[abs + 4],
                    self.record[abs + 5],
                    self.record[abs + 6],
                    self.record[abs + 7],
                ]) as usize;
                if count == values.len() {
                    let body = abs + 8;
                    for (i, &v) in values.iter().enumerate() {
                        let dst = body + i * 2;
                        self.record[dst..dst + 2].copy_from_slice(&v.to_le_bytes());
                    }
                    return;
                }
            }
            if let Some(size) = tag_value_size(val_type, &self.record[abs + 3..]) {
                self.record.drain(abs..abs + 3 + size);
            }
        }
        append_u16_array_tag(self.record, tag, values);
    }

    /// Update an existing `B:s` (i16 array) tag in place if the length matches, else remove and
    /// append.
    pub fn update_array_i16(&mut self, tag: impl AsTagBytes, values: &[i16]) {
        let tag = tag.as_tag_bytes();
        let off = self.aux_offset.min(self.record.len());
        if let Some((p, val_type)) = find_tag_position(&self.record[off..], *tag) {
            let abs = off + p;
            if val_type == b'B' && abs + 8 <= self.record.len() && self.record[abs + 3] == b's' {
                let count = u32::from_le_bytes([
                    self.record[abs + 4],
                    self.record[abs + 5],
                    self.record[abs + 6],
                    self.record[abs + 7],
                ]) as usize;
                if count == values.len() {
                    let body = abs + 8;
                    for (i, &v) in values.iter().enumerate() {
                        let dst = body + i * 2;
                        self.record[dst..dst + 2].copy_from_slice(&v.to_le_bytes());
                    }
                    return;
                }
            }
            if let Some(size) = tag_value_size(val_type, &self.record[abs + 3..]) {
                self.record.drain(abs..abs + 3 + size);
            }
        }
        append_i16_array_tag(self.record, tag, values);
    }

    /// Update an existing `B:i` (i32 array) tag in place if the length matches, else remove and
    /// append.
    pub fn update_array_i32(&mut self, tag: impl AsTagBytes, values: &[i32]) {
        let tag = tag.as_tag_bytes();
        let off = self.aux_offset.min(self.record.len());
        if let Some((p, val_type)) = find_tag_position(&self.record[off..], *tag) {
            let abs = off + p;
            if val_type == b'B' && abs + 8 <= self.record.len() && self.record[abs + 3] == b'i' {
                let count = u32::from_le_bytes([
                    self.record[abs + 4],
                    self.record[abs + 5],
                    self.record[abs + 6],
                    self.record[abs + 7],
                ]) as usize;
                if count == values.len() {
                    let body = abs + 8;
                    for (i, &v) in values.iter().enumerate() {
                        let dst = body + i * 4;
                        self.record[dst..dst + 4].copy_from_slice(&v.to_le_bytes());
                    }
                    return;
                }
            }
            if let Some(size) = tag_value_size(val_type, &self.record[abs + 3..]) {
                self.record.drain(abs..abs + 3 + size);
            }
        }
        append_i32_array_tag(self.record, tag, values);
    }

    /// Update an existing `B:c` (i8 array) tag in place if the length matches, else remove and
    /// append.
    pub fn update_array_i8(&mut self, tag: impl AsTagBytes, values: &[i8]) {
        let tag = tag.as_tag_bytes();
        let off = self.aux_offset.min(self.record.len());
        if let Some((p, val_type)) = find_tag_position(&self.record[off..], *tag) {
            let abs = off + p;
            if val_type == b'B' && abs + 8 <= self.record.len() && self.record[abs + 3] == b'c' {
                let count = u32::from_le_bytes([
                    self.record[abs + 4],
                    self.record[abs + 5],
                    self.record[abs + 6],
                    self.record[abs + 7],
                ]) as usize;
                if count == values.len() {
                    let body = abs + 8;
                    for (i, &v) in values.iter().enumerate() {
                        self.record[body + i] = v.cast_unsigned();
                    }
                    return;
                }
            }
            if let Some(size) = tag_value_size(val_type, &self.record[abs + 3..]) {
                self.record.drain(abs..abs + 3 + size);
            }
        }
        append_i8_array_tag(self.record, tag, values);
    }

    /// Update an existing `B:I` (u32 array) tag in place if the length matches, else remove and
    /// append.
    pub fn update_array_u32(&mut self, tag: impl AsTagBytes, values: &[u32]) {
        let tag = tag.as_tag_bytes();
        let off = self.aux_offset.min(self.record.len());
        if let Some((p, val_type)) = find_tag_position(&self.record[off..], *tag) {
            let abs = off + p;
            if val_type == b'B' && abs + 8 <= self.record.len() && self.record[abs + 3] == b'I' {
                let count = u32::from_le_bytes([
                    self.record[abs + 4],
                    self.record[abs + 5],
                    self.record[abs + 6],
                    self.record[abs + 7],
                ]) as usize;
                if count == values.len() {
                    let body = abs + 8;
                    for (i, &v) in values.iter().enumerate() {
                        let dst = body + i * 4;
                        self.record[dst..dst + 4].copy_from_slice(&v.to_le_bytes());
                    }
                    return;
                }
            }
            if let Some(size) = tag_value_size(val_type, &self.record[abs + 3..]) {
                self.record.drain(abs..abs + 3 + size);
            }
        }
        append_u32_array_tag(self.record, tag, values);
    }

    /// Update an existing `B:f` (f32 array) tag in place if the length matches, else remove and
    /// append.
    pub fn update_array_f32(&mut self, tag: impl AsTagBytes, values: &[f32]) {
        let tag = tag.as_tag_bytes();
        let off = self.aux_offset.min(self.record.len());
        if let Some((p, val_type)) = find_tag_position(&self.record[off..], *tag) {
            let abs = off + p;
            if val_type == b'B' && abs + 8 <= self.record.len() && self.record[abs + 3] == b'f' {
                let count = u32::from_le_bytes([
                    self.record[abs + 4],
                    self.record[abs + 5],
                    self.record[abs + 6],
                    self.record[abs + 7],
                ]) as usize;
                if count == values.len() {
                    let body = abs + 8;
                    for (i, &v) in values.iter().enumerate() {
                        let dst = body + i * 4;
                        self.record[dst..dst + 4].copy_from_slice(&v.to_le_bytes());
                    }
                    return;
                }
            }
            if let Some(size) = tag_value_size(val_type, &self.record[abs + 3..]) {
                self.record.drain(abs..abs + 3 + size);
            }
        }
        append_f32_array_tag(self.record, tag, values);
    }

    /// Remove a tag from the record. No-op if the tag is not found.
    #[inline]
    pub fn remove(&mut self, tag: impl AsTagBytes) {
        remove_tag(self.record, tag);
    }

    /// Re-encode an existing integer tag using the smallest fitting signed type.
    ///
    /// No-op if the tag is not found or is not an integer type.
    #[inline]
    pub fn normalize_int_to_smallest_signed(&mut self, tag: impl AsTagBytes) {
        normalize_int_tag_to_smallest_signed(self.record, tag);
    }

    /// Copy entries from a source aux view to this record, optionally skipping tags by key.
    ///
    /// All tags in `src` are appended to the record unless their two-byte key appears in `skip`.
    #[inline]
    pub fn copy_from(&mut self, src: RawTagsView<'_>, skip: &[SamTag]) {
        copy_aux_tags(src.as_bytes(), self.record, skip);
    }

    /// Rebuild the aux block in a single pass: drop every existing tag whose
    /// two-byte key is in `remove`, keep the rest in their original order, then
    /// append `adds` (in slice order).
    ///
    /// Upsert semantics: an added tag whose key also exists among the survivors
    /// is dropped from the survivors first, so the appended value wins with no
    /// duplicate — matching the `remove_tag(tag); append(tag, ..)` idiom this
    /// replaces. If `adds` itself carries the same key more than once, only the
    /// last occurrence is emitted (again matching the idiom, where each
    /// `remove_tag` wipes the prior append), so the output never contains a
    /// duplicate key.
    ///
    /// Cost: a single walk of the aux with one allocation and one splice back,
    /// versus N independent O(aux) `remove_tag` scans (each its own
    /// `drain`/`extend` splice) for N per-tag operations. The per-survivor
    /// `adds`-membership check makes the walk O(aux · |adds|) in the strict
    /// sense, but `adds` is expected to be tiny (a handful of tags), so in
    /// practice it is O(aux); the win over the old idiom is the collapsed
    /// allocations and splices, not the asymptotic key-compare count.
    ///
    /// It does NOT reverse/revcomp: a caller transforming negative-strand tags
    /// runs `reverse_*_in_place` as a separate pass over the rebuilt record.
    ///
    /// A record too short to hold a valid header (aux offset past the end) is
    /// treated as having an empty aux block — `rebuild_with` then appends
    /// `adds` only. A malformed entry mid-aux stops the walk (same tolerance as
    /// [`AuxTagsIter`]); bytes past it are dropped.
    pub fn rebuild_with<M: TagKeySet + ?Sized>(&mut self, remove: &M, adds: &[TagEntry<'_>]) {
        let off = self.aux_offset.min(self.record.len());

        // A key that will be (re)appended from `adds` is dropped from the
        // survivors so the appended value wins with no duplicate. `adds` is
        // tiny (a handful of tags), so a linear membership check per survivor is
        // cheaper than allocating a second set.
        let dropped_by_add = |tag: [u8; 2]| adds.iter().any(|a| a.tag == tag);

        // Size the rebuilt aux from the current aux length plus the adds, so the
        // buffer is allocated once.
        let adds_len: usize = adds.iter().map(|a| 3 + a.value_bytes.len()).sum();
        let mut new_aux = Vec::with_capacity((self.record.len() - off) + adds_len);

        // Single pass over the existing aux: copy survivors verbatim.
        for entry in &RawTagsView::new(&self.record[off..]) {
            if remove.contains_tag(entry.tag) || dropped_by_add(entry.tag) {
                continue;
            }
            new_aux.push(entry.tag[0]);
            new_aux.push(entry.tag[1]);
            new_aux.push(entry.type_byte);
            new_aux.extend_from_slice(entry.value_bytes);
        }

        // Append the new entries in slice order. If a key appears more than once
        // in `adds`, only its last occurrence is emitted, so the output carries
        // no duplicate key — matching the `remove_tag(tag); append(tag, ..)`
        // idiom, where each per-tag `remove_tag` wipes the prior append.
        for (idx, a) in adds.iter().enumerate() {
            if adds[idx + 1..].iter().any(|b| b.tag == a.tag) {
                continue;
            }
            new_aux.push(a.tag[0]);
            new_aux.push(a.tag[1]);
            new_aux.push(a.type_byte);
            new_aux.extend_from_slice(a.value_bytes);
        }

        // One splice back over the old aux region (off..end).
        self.record.splice(off.., new_aux);
    }

    /// Like [`Self::rebuild_with`], but during the same single walk it also normalizes
    /// each integer tag in `normalize` to the smallest *signed* width that fits
    /// (i8 → i16 → i32, matching fgbio's `to_smallest_signed_int`), relocating it
    /// to the tail of the rebuilt aux — reproducing byte-for-byte what
    /// `rebuild_with(remove, adds)` followed by
    /// `normalize_int_tag_to_smallest_signed(t)` for each `t` in `normalize`
    /// produces today, but without the extra `find_tag_position` scans and
    /// `Vec::drain` memmoves that the separate per-tag normalization pass costs.
    pub fn rebuild_with_int_normalized<M: TagKeySet + ?Sized>(
        &mut self,
        remove: &M,
        adds: &[TagEntry<'_>],
        normalize: &[[u8; 2]],
        scratch: &mut Vec<u8>,
    ) {
        /// Upper bound on `normalize.len()` so the captured-values array lives on
        /// the stack (no per-record allocation). zipper passes 2 (AS, XS); the
        /// cap is generous. Overflow degrades gracefully to a heap `Vec`.
        const INLINE_NORMALIZE: usize = 8;

        // Nothing to normalize: identical to `rebuild_with`.
        if normalize.is_empty() {
            self.rebuild_with(remove, adds);
            return;
        }

        // Duplicate keys in `normalize` fall back to the literal sequential
        // oracle. With a repeated key, `normalize_int_tag_to_smallest_signed`
        // peels off and relocates the FIRST remaining occurrence on each call,
        // so N occurrences of a key normalize N distinct aux entries — which the
        // single-pass fast path (one capture slot per distinct key) cannot
        // reproduce. A `normalize` list with repeated keys is degenerate and
        // never happens on the hot path (zipper passes [AS, XS]), so paying the
        // extra scans there keeps the fast path simple and provably correct.
        if normalize.iter().enumerate().any(|(i, tag)| normalize[i + 1..].contains(tag)) {
            self.rebuild_with(remove, adds);
            for &tag in normalize {
                normalize_int_tag_to_smallest_signed(self.record, tag);
            }
            return;
        }

        let off = self.aux_offset.min(self.record.len());
        let dropped_by_add = |tag: [u8; 2]| adds.iter().any(|a| a.tag == tag);

        // For each `normalize` tag, the i32 value to relocate-and-normalize at the
        // tail, or `None` to leave it untouched. Seeded from `adds` (an added tag
        // wins the upsert, so it is the value that gets normalized — matching
        // Step 3 running before Step 5 today); the survivor walk below fills in
        // the value for a tag that `adds` does not carry. A tag whose winning
        // value is not an integer, or does not fit i32, is NOT captured here —
        // it stays verbatim wherever it already sits, exactly as
        // `normalize_int_tag_to_smallest_signed`'s early-return leaves it.
        //
        // Held in a small stack array for the common case (≤ INLINE_NORMALIZE
        // tags), spilling to the heap only if a caller passes more.
        let mut captured_inline = [None::<i32>; INLINE_NORMALIZE];
        let mut captured_spill: Vec<Option<i32>> = if normalize.len() > INLINE_NORMALIZE {
            vec![None; normalize.len()]
        } else {
            Vec::new()
        };
        let captured: &mut [Option<i32>] = if normalize.len() > INLINE_NORMALIZE {
            &mut captured_spill
        } else {
            &mut captured_inline[..normalize.len()]
        };
        for (i, &t) in normalize.iter().enumerate() {
            captured[i] = adds
                .iter()
                .rev()
                .find(|a| a.tag == t)
                .and_then(|a| int_from_value_bytes(a.type_byte, a.value_bytes))
                .and_then(|v| i32::try_from(v).ok());
        }

        // Tracks the FIRST survivor occurrence of each normalize tag, so a
        // duplicate-key aux (out of BAM spec but not impossible) matches the
        // sequential oracle exactly: `find_int_tag`/`remove_tag` act on the first
        // key match only, so only the first occurrence is captured/relocated and
        // any later duplicate is passed through verbatim — including the case
        // where the first occurrence is non-integer (which must then leave a
        // later integer duplicate untouched, not capture it). One flag per tag,
        // sized to `normalize.len()` (inline for the common case, spilling to the
        // heap alongside `captured` when a caller passes more than the cap) so the
        // seen-state scales to any list length rather than silently dropping
        // entries past a fixed bitmask width.
        let mut seen_inline = [false; INLINE_NORMALIZE];
        let mut seen_spill: Vec<bool> = if normalize.len() > INLINE_NORMALIZE {
            vec![false; normalize.len()]
        } else {
            Vec::new()
        };
        let seen: &mut [bool] = if normalize.len() > INLINE_NORMALIZE {
            &mut seen_spill
        } else {
            &mut seen_inline[..normalize.len()]
        };

        // Build the rebuilt aux into the caller's reusable scratch buffer (its
        // allocation is reused across records on the serial merge thread), then
        // splice it back via `drain` so the scratch keeps its capacity.
        scratch.clear();
        let adds_len: usize = adds.iter().map(|a| 3 + a.value_bytes.len()).sum();
        scratch.reserve((self.record.len() - off) + adds_len);

        // Single pass over the existing aux: copy survivors, capturing any
        // normalize tag that is the upsert winner and re-encodable to i32.
        for entry in &RawTagsView::new(&self.record[off..]) {
            if remove.contains_tag(entry.tag) || dropped_by_add(entry.tag) {
                continue;
            }
            // A normalize tag reaching here is the winner (`adds` does not carry
            // it, else it was dropped above). Only the FIRST occurrence is
            // considered (first-key-match semantics of `find_int_tag`): capture
            // and drop it iff it is an integer that fits i32; otherwise, and for
            // every later duplicate occurrence, leave it verbatim in place.
            if let Some(i) = normalize.iter().position(|&n| n == entry.tag)
                && !seen[i]
            {
                seen[i] = true;
                if let Some(v) = int_from_value_bytes(entry.type_byte, entry.value_bytes)
                    .and_then(|v| i32::try_from(v).ok())
                {
                    captured[i] = Some(v);
                    continue;
                }
            }
            scratch.push(entry.tag[0]);
            scratch.push(entry.tag[1]);
            scratch.push(entry.type_byte);
            scratch.extend_from_slice(entry.value_bytes);
        }

        // Append the adds (dedup last-wins, as `rebuild_with`), skipping any add
        // that is being relocated as a normalized tag (captured from `adds`).
        for (idx, a) in adds.iter().enumerate() {
            if adds[idx + 1..].iter().any(|b| b.tag == a.tag) {
                continue;
            }
            let relocated =
                normalize.iter().position(|&n| n == a.tag).is_some_and(|i| captured[i].is_some());
            if relocated {
                continue;
            }
            scratch.push(a.tag[0]);
            scratch.push(a.tag[1]);
            scratch.push(a.type_byte);
            scratch.extend_from_slice(a.value_bytes);
        }

        // Append the normalized tags at the tail, in `normalize` order, using the
        // same smallest-signed encoder the sequential normalize pass used — so
        // the tail bytes and order are byte-identical to today. `normalize` has
        // distinct keys here (a duplicate-key list took the sequential fallback
        // above), so each captured slot maps to a unique tag emitted once.
        for (i, &t) in normalize.iter().enumerate() {
            if let Some(v) = captured[i] {
                append_signed_int_tag(scratch, t, v);
            }
        }

        self.record.splice(off.., scratch.drain(..));
    }
}

#[cfg(test)]
#[allow(clippy::identity_op)]
mod tests {
    use super::*;
    use crate::testutil::*;
    use fgumi_tag::SamTag;
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
    // find_string_tag tests
    // ========================================================================

    #[test]
    fn test_find_string_tag_present() {
        let aux = b"RGZsample1\x00";
        assert_eq!(find_string_tag(aux, SamTag::RG), Some(b"sample1".as_ref()));
    }

    #[test]
    fn test_find_string_tag_absent() {
        let aux = b"RGZsample1\x00";
        assert_eq!(find_string_tag(aux, SamTag::RX), None);
    }

    #[test]
    fn test_find_string_tag_after_other_tags() {
        // NM:C:5 then RX:Z:ACGT
        let mut aux = Vec::new();
        aux.extend_from_slice(b"NMC");
        aux.push(5);
        aux.extend_from_slice(b"RXZACGT\x00");
        assert_eq!(find_string_tag(&aux, SamTag::RX), Some(b"ACGT".as_ref()));
    }

    #[test]
    fn test_find_string_tag_in_record() {
        let aux = b"RXZhello\x00";
        let rec = make_bam_bytes(0, 0, 0, b"rea", &[], 0, -1, -1, aux);
        assert_eq!(find_string_tag_in_record(&rec, SamTag::RX), Some(b"hello".as_ref()));
    }

    #[test]
    fn test_find_string_tag_non_z_type_returns_none() {
        // Tag matches but type is 'C' not 'Z'
        let aux = [b'R', b'X', b'C', 42];
        assert_eq!(find_string_tag(&aux, SamTag::RX), None);
    }

    #[test]
    fn test_find_string_tag_in_record_no_aux() {
        // Record with no aux data (aux_start >= len)
        let rec = make_bam_bytes(0, 0, 0, b"rea", &[], 0, -1, -1, &[]);
        assert_eq!(find_string_tag_in_record(&rec, SamTag::RX), None);
    }

    #[test]
    fn test_find_string_tag_truncated_z_value() {
        // RX:Z:hello but no null terminator — should return None
        let aux = b"RXZhello";
        assert_eq!(find_string_tag(aux, SamTag::RX), None);
    }

    // ========================================================================
    // find_string_tag_position tests
    // ========================================================================

    #[test]
    fn test_find_string_tag_position_first_tag() {
        // "RXZACGT\0" — tag header occupies bytes 0..3; value occupies 3..7; NUL at 7
        let aux = b"RXZACGT\x00";
        let (offset, len) =
            find_string_tag_position(aux, SamTag::RX).expect("UMI tag should be present");
        assert_eq!(offset, 3);
        assert_eq!(len, 4);
        let bytes = &aux[offset as usize..offset as usize + len as usize];
        assert_eq!(bytes, b"ACGT");
        assert_eq!(bytes, find_string_tag(aux, SamTag::RX).unwrap());
    }

    #[test]
    fn test_find_string_tag_position_after_other_tags() {
        // NM:C:5 then RX:Z:ACGT — the UMI value should be located after the NM entry.
        let mut aux = Vec::new();
        aux.extend_from_slice(b"NMC");
        aux.push(5);
        aux.extend_from_slice(b"RXZACGT\x00");
        let (offset, len) = find_string_tag_position(&aux, SamTag::RX).unwrap();
        assert_eq!(&aux[offset as usize..offset as usize + len as usize], b"ACGT");
        assert_eq!(len, 4);
    }

    #[test]
    fn test_find_string_tag_position_empty_value() {
        // RX:Z: (empty UMI) — len should be 0 with a valid offset.
        let aux = b"RXZ\x00";
        let (offset, len) = find_string_tag_position(aux, SamTag::RX).unwrap();
        assert_eq!(offset, 3);
        assert_eq!(len, 0);
    }

    #[test]
    fn test_find_string_tag_position_absent_returns_none() {
        let aux = b"RGZsample1\x00";
        assert_eq!(find_string_tag_position(aux, SamTag::RX), None);
    }

    #[test]
    fn test_find_string_tag_position_non_z_returns_none() {
        // Tag matches but type is 'C' (uint8), not 'Z'.
        let aux = [b'R', b'X', b'C', 42];
        assert_eq!(find_string_tag_position(&aux, SamTag::RX), None);
    }

    #[test]
    fn test_find_string_tag_position_truncated_no_nul() {
        // RX:Z:hello but no NUL terminator — should return None.
        let aux = b"RXZhello";
        assert_eq!(find_string_tag_position(aux, SamTag::RX), None);
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
        let aux = make_b_uint16_array_tag(SamTag::XS.into(), &[100, 200, 300]);
        assert_eq!(find_string_tag(&aux, SamTag::XS), None);
    }

    // --- B:I array (uint32 array) ---

    #[test]
    fn test_find_string_tag_cannot_find_b_uint32_array() {
        let aux = make_b_uint32_array_tag(*b"XI", &[1000, 2000, 3000]);
        assert_eq!(find_string_tag(&aux, b"XI"), None);
    }

    // --- B:i array (the type the `tc` template-coordinate tag uses) ---

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
        assert_eq!(find_string_tag(aux, SamTag::RX), Some(b"hello".as_ref()));
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
        assert_eq!(find_string_tag(aux, SamTag::XS), None);
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
        assert_eq!(find_string_tag(&aux, SamTag::RX), Some(b"hello".as_ref()));
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
        assert_eq!(find_tag_bounds(&aux, SamTag::NM), Some((0, 4)));
        assert_eq!(find_tag_bounds(&aux, SamTag::RX), Some((4, 11)));
        assert_eq!(find_tag_bounds(&aux, b"XX"), None);
    }

    #[test]
    fn test_find_tag_bounds_empty() {
        assert_eq!(find_tag_bounds(&[], SamTag::RX), None);
    }

    #[test]
    fn test_find_tag_bounds_truncated() {
        // Only 2 bytes — not enough for a tag entry
        assert_eq!(find_tag_bounds(b"RX", SamTag::RX), None);
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

    // ========================================================================
    // extract_int_value bounds tests
    // ========================================================================

    /// `extract_int_value` is public and takes `p` from the caller rather than from a
    /// `find_tag_position` scan, so a position near `usize::MAX` reaches the decode
    /// ladder. With unchecked `p + 3` / `p + 4` / `p + 7` the guard wraps to a small
    /// number, passes the length check, and the ladder then decodes bytes from a
    /// wrapped offset — the wrong value, silently, in a release build.
    #[rstest]
    #[case::signed_byte(b'c')]
    #[case::unsigned_byte(b'C')]
    #[case::signed_short(b's')]
    #[case::unsigned_short(b'S')]
    #[case::signed_int(b'i')]
    #[case::unsigned_int(b'I')]
    fn test_extract_int_value_rejects_positions_that_overflow(#[case] type_byte: u8) {
        let aux = [b'X', b'Y', type_byte, 1, 2, 3, 4, 5];
        // Every `p` within a value width of the top of the address space: `p + 3`,
        // `p + 4`, and `p + 7` all wrap for at least one of these.
        for delta in 0..8usize {
            let p = usize::MAX - delta;
            assert_eq!(
                extract_int_value(&aux, p, type_byte),
                None,
                "type {} at p = usize::MAX - {delta} must not decode",
                type_byte as char
            );
        }
    }

    /// The in-bounds ladder is unchanged: a value that ends exactly at the end of the
    /// aux block still decodes, and one byte short of that still does not.
    #[rstest]
    #[case::signed_byte(b'c', 1)]
    #[case::unsigned_byte(b'C', 1)]
    #[case::signed_short(b's', 2)]
    #[case::unsigned_short(b'S', 2)]
    #[case::signed_int(b'i', 4)]
    #[case::unsigned_int(b'I', 4)]
    fn test_extract_int_value_boundary(#[case] type_byte: u8, #[case] width: usize) {
        let mut exact = vec![b'X', b'Y', type_byte];
        exact.extend(std::iter::repeat_n(0u8, width));
        assert_eq!(extract_int_value(&exact, 0, type_byte), Some(0));

        let truncated = &exact[..exact.len() - 1];
        assert_eq!(extract_int_value(truncated, 0, type_byte), None);
    }

    /// `decode_int_value` decodes each integer type from bare value bytes (no
    /// tag/type header), returns `None` for a non-integer type, and `None` when
    /// the slice is one byte short of the type's width.
    #[rstest]
    #[case::signed_byte_neg(b'c', &[(-5i8).cast_unsigned()], Some(-5))]
    #[case::unsigned_byte(b'C', &[200u8], Some(200))]
    #[case::signed_short_neg(b's', &(-200i16).to_le_bytes(), Some(-200))]
    #[case::unsigned_short(b'S', &50_000u16.to_le_bytes(), Some(50_000))]
    #[case::signed_int(b'i', &(-100_000i32).to_le_bytes(), Some(-100_000))]
    #[case::unsigned_int(b'I', &3_000_000_000u32.to_le_bytes(), Some(3_000_000_000))]
    #[case::float_type_is_not_int(b'f', &1.0f32.to_le_bytes(), None)]
    #[case::string_type_is_not_int(b'Z', b"5\x00", None)]
    fn test_decode_int_value(
        #[case] type_byte: u8,
        #[case] value: &[u8],
        #[case] expected: Option<i64>,
    ) {
        assert_eq!(decode_int_value(type_byte, value), expected);
    }

    #[rstest]
    #[case::signed_byte(b'c', 1)]
    #[case::signed_short(b's', 2)]
    #[case::signed_int(b'i', 4)]
    fn test_decode_int_value_short_slice_is_none(#[case] type_byte: u8, #[case] width: usize) {
        let short = vec![0u8; width - 1];
        assert_eq!(decode_int_value(type_byte, &short), None);
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
        let aux = make_b_uint16_array_tag(SamTag::XS.into(), &[100, 200, 300]);
        assert_eq!(find_int_tag(&aux, SamTag::XS), None);
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
        assert_eq!(find_int_tag(aux, SamTag::RX), None);
    }

    #[test]
    fn test_find_int_tag_finds_c_type() {
        let aux: &[u8] = &[b'X', b'c', b'c', 5];
        assert_eq!(find_int_tag(aux, b"Xc"), Some(5));
    }

    #[test]
    fn test_find_int_tag_finds_upper_c_type() {
        let aux: &[u8] = &[b'N', b'M', b'C', 42];
        assert_eq!(find_int_tag(aux, SamTag::NM), Some(42));
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
        assert_eq!(find_int_tag(aux, SamTag::XS), Some(50_000));
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
        assert_eq!(find_int_tag(&aux, SamTag::NM), Some(5));
    }

    // ========================================================================
    // find_uint8_tag tests
    // ========================================================================

    #[test]
    fn test_find_uint8_tag() {
        let aux = [b'M', b'Q', b'C', 30];
        assert_eq!(find_uint8_tag(&aux, SamTag::MQ), Some(30));
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
        assert_eq!(find_uint8_tag(&aux, SamTag::MQ), expected_uint8);
        assert_eq!(find_int_tag(&aux, SamTag::MQ), expected_int);
    }

    #[test]
    fn test_find_uint8_tag_wrong_type() {
        // Tag exists as type 'i' not 'C' — should not match
        let mut aux = vec![b'M', b'Q', b'i'];
        aux.extend_from_slice(&42i32.to_le_bytes());
        assert_eq!(find_uint8_tag(&aux, SamTag::MQ), None);
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
        let result = find_float_tag(&aux, b"XF").expect("XF float tag should be found");
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
        assert_eq!(find_mc_tag(aux), Some(b"10M5S".as_slice()));
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
        assert_eq!(find_mc_tag(&aux), Some(b"15M".as_slice()));
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
        let result = extract_aux_string_tags(&aux, SamTag::CB, None);
        assert_eq!(result.rg, Some(b"sample1".as_ref()));
        assert_eq!(result.cell, Some(b"cell42".as_ref()));
        assert_eq!(result.mc, Some(b"10M5S".as_slice()));
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
        let result = extract_aux_string_tags(&aux, SamTag::CB, None);
        assert_eq!(result.rg, Some(b"rg1".as_ref()));
        assert_eq!(result.cell, Some(b"bc1".as_ref()));
        assert_eq!(result.mc, Some(b"5M".as_slice()));
    }

    #[test]
    fn test_extract_aux_string_tags_partial() {
        // Only RG present, others missing
        let aux = b"RGZsample\x00";
        let result = extract_aux_string_tags(aux, SamTag::CB, None);
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
        let result = extract_aux_string_tags(&aux, SamTag::CB, None);
        assert_eq!(result.rg, Some(b"lib1".as_ref()));
        assert!(result.cell.is_none());
        assert_eq!(result.mc, Some(b"20M".as_slice()));
    }

    #[test]
    fn test_extract_aux_string_tags_empty() {
        let result = extract_aux_string_tags(&[], SamTag::CB, None);
        assert!(result.rg.is_none());
        assert!(result.cell.is_none());
        assert!(result.mc.is_none());
    }

    /// A non-UTF-8 `MC` is handed over as bytes rather than discarded.
    ///
    /// This is a deliberate behavior change. Extraction used to run
    /// `str::from_utf8(..).ok()` on the value, so a tag with any non-UTF-8 byte
    /// became `None` and the mate position was used unadjusted. `MC` is a CIGAR
    /// that every consumer immediately re-borrows as bytes, and validating it
    /// per record cost 2.1% of `fgumi sort`'s serial Phase 1 thread, so the gate
    /// is gone.
    ///
    /// The observable difference is narrow: the clip parsers read a leading run
    /// of ASCII digits and stop at the first byte they cannot interpret, so a
    /// value that is garbage from the start still yields zero clips -- the same
    /// answer discarding it gave. Only a value with a *valid CIGAR prefix*
    /// followed by invalid bytes now changes, and it changes toward parsing the
    /// prefix, which is what samtools does.
    #[test]
    fn test_extract_aux_string_tags_keeps_a_non_utf8_mc_as_bytes() {
        let mut aux = Vec::new();
        aux.extend_from_slice(b"MCZ");
        aux.extend_from_slice(&[0xFF, 0xFE, 0xFD]); // not valid UTF-8
        aux.push(0); // NUL terminator
        let result = extract_aux_string_tags(&aux, SamTag::CB, None);
        assert_eq!(result.mc, Some([0xFF, 0xFE, 0xFD].as_slice()));
        // Garbage from the first byte still parses to no clips, so this value
        // places the mate exactly where discarding the tag used to.
        assert_eq!(crate::cigar::mate_unclipped_5prime(100, false, result.mc.expect("mc")), 100);
    }

    #[test]
    fn test_extract_aux_string_tags_truncated_z() {
        // RG:Z:lib1 with no null → should break out, no tags found
        let aux = b"RGZlib1";
        let result = extract_aux_string_tags(aux, SamTag::CB, None);
        assert!(result.rg.is_none());
    }

    #[test]
    fn test_extract_aux_string_tags_captures_umi_position() {
        // Build aux: NM:i:5  RX:Z:ACGT  RG:Z:lib1
        let mut aux = Vec::new();
        aux.extend_from_slice(b"NMC");
        aux.push(5);
        aux.extend_from_slice(b"RXZACGT\x00");
        aux.extend_from_slice(b"RGZlib1\x00");
        let result = extract_aux_string_tags(&aux, SamTag::CB, Some(*SamTag::RX));
        assert_eq!(result.rg, Some(b"lib1".as_ref()));
        let (off, len) = result.umi_position.expect("UMI position should be captured");
        assert_eq!(&aux[off as usize..off as usize + len as usize], b"ACGT");
    }

    #[test]
    fn test_extract_aux_string_tags_umi_none_when_tag_absent() {
        // Aux has RG but no RX.
        let aux = b"RGZlib1\x00";
        let result = extract_aux_string_tags(aux, SamTag::CB, Some(*SamTag::RX));
        assert!(result.rg.is_some());
        assert!(result.umi_position.is_none());
    }

    #[test]
    fn test_extract_aux_string_tags_umi_skipped_when_caller_disinterested() {
        // RX present but caller passed umi_tag=None — must not capture.
        let aux = b"RXZACGT\x00";
        let result = extract_aux_string_tags(aux, SamTag::CB, None);
        assert!(result.umi_position.is_none());
    }

    /// Overlapping selection: when `cell_tag == MC`, the single `MC` aux entry
    /// must populate BOTH the cell field and the mc field. The chained
    /// `else if` this replaced captured only the first-matching field (cell) and
    /// left `mc` = None, which on the grouping/dedup key path silently dropped
    /// the mate-CIGAR contribution to the group key. Parity with the sibling
    /// `extract_template_aux_tags`, which already uses independent `if` checks.
    #[test]
    fn test_extract_aux_string_tags_captures_overlapping_cell_and_mc() {
        let mut aux = Vec::new();
        aux.extend_from_slice(b"RGZlib1\x00");
        aux.extend_from_slice(b"MCZ5M\x00");
        // cell_tag == MC: the one MC entry is selected as both cell and mc.
        let result = extract_aux_string_tags(&aux, SamTag::MC, None);
        assert_eq!(result.rg, Some(b"lib1".as_ref()));
        assert_eq!(result.cell, Some(b"5M".as_ref()), "cell must capture the MC entry");
        assert_eq!(result.mc, Some(b"5M".as_slice()), "mc must still capture the same entry");
    }

    /// Overlapping selection on the UMI axis: when `umi_tag == RG`, the single
    /// `RG` entry must populate BOTH the rg field and the UMI position.
    #[test]
    fn test_extract_aux_string_tags_captures_overlapping_rg_and_umi() {
        let aux = b"RGZlib1\x00";
        let result = extract_aux_string_tags(aux, SamTag::CB, Some(*SamTag::RG));
        assert_eq!(result.rg, Some(b"lib1".as_ref()), "rg must still capture");
        let (off, len) = result.umi_position.expect("UMI position must also capture the RG entry");
        assert_eq!(&aux[off as usize..off as usize + len as usize], b"lib1");
    }

    /// Non-overlap pin: with disjoint selected tags (the default `cell_tag ==
    /// CB`), each field is sourced from its own entry and nothing cross-fills.
    #[test]
    fn test_extract_aux_string_tags_non_overlapping_tags_stay_separate() {
        let mut aux = Vec::new();
        aux.extend_from_slice(b"RGZlib1\x00");
        aux.extend_from_slice(b"CBZbc1\x00");
        aux.extend_from_slice(b"MCZ5M\x00");
        let result = extract_aux_string_tags(&aux, SamTag::CB, None);
        assert_eq!(result.rg, Some(b"lib1".as_ref()));
        assert_eq!(result.cell, Some(b"bc1".as_ref()));
        assert_eq!(result.mc, Some(b"5M".as_slice()));
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
        let aux = make_b_uint16_array_tag(SamTag::XS.into(), &[100, 200, 300]);
        assert_eq!(find_tag_type(&aux, SamTag::XS), Some(b'B'));
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
        assert_eq!(find_tag_type(aux, SamTag::RX), Some(b'Z'));
    }

    #[test]
    fn test_find_tag_type_finds_c_type() {
        let aux: &[u8] = &[b'X', b'c', b'c', 5];
        assert_eq!(find_tag_type(aux, b"Xc"), Some(b'c'));
    }

    #[test]
    fn test_find_tag_type_finds_upper_c_type() {
        let aux: &[u8] = &[b'N', b'M', b'C', 42];
        assert_eq!(find_tag_type(aux, SamTag::NM), Some(b'C'));
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
        assert_eq!(find_tag_type(aux, SamTag::XS), Some(b'S'));
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
        assert_eq!(find_tag_type(&aux, SamTag::RX), Some(b'Z'));
    }

    // ========================================================================
    // append_string_tag tests
    // ========================================================================

    #[test]
    fn test_append_string_tag() {
        let mut rec = make_bam_bytes(0, 0, 0, b"rea", &[], 0, -1, -1, &[]);
        let orig_len = rec.len();
        append_string_tag(&mut rec, SamTag::MI, b"12345");
        assert_eq!(rec.len(), orig_len + 2 + 1 + 5 + 1); // tag(2) + type(1) + value(5) + NUL(1)
        // Verify we can find it back
        let aux_start = aux_data_offset_from_record(&rec)
            .expect("record should have valid header for aux offset");
        assert_eq!(find_string_tag(&rec[aux_start..], SamTag::MI), Some(b"12345".as_ref()));
    }

    // ========================================================================
    // remove_tag tests
    // ========================================================================

    #[test]
    fn test_remove_tag_present() {
        let aux = b"MIZ42\x00";
        let mut rec = make_bam_bytes(0, 0, 0, b"rea", &[], 0, -1, -1, aux);
        assert!(find_string_tag_in_record(&rec, SamTag::MI).is_some());
        remove_tag(&mut rec, SamTag::MI);
        assert!(find_string_tag_in_record(&rec, SamTag::MI).is_none());
    }

    #[test]
    fn test_remove_tag_absent() {
        let mut rec = make_bam_bytes(0, 0, 0, b"rea", &[], 0, -1, -1, &[]);
        let orig_len = rec.len();
        remove_tag(&mut rec, SamTag::MI); // should be no-op
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
        let aux_start = aux_data_offset_from_record(&rec)
            .expect("record should have valid header for aux offset");
        assert_eq!(find_uint8_tag(&rec[aux_start..], b"AA"), Some(1));
        assert_eq!(find_uint8_tag(&rec[aux_start..], b"CC"), Some(2));
        assert!(find_string_tag(&rec[aux_start..], b"BB").is_none());
    }

    #[test]
    fn test_remove_tag_no_aux_data() {
        // Record where aux_data_offset >= record length
        let mut rec = make_bam_bytes(0, 0, 0, b"rea", &[], 0, -1, -1, &[]);
        let orig_len = rec.len();
        remove_tag(&mut rec, SamTag::MI);
        assert_eq!(rec.len(), orig_len); // no-op
    }

    // ========================================================================
    // update_string_tag tests
    // ========================================================================

    #[test]
    fn test_update_string_tag_existing() {
        let aux = b"RXZold\x00";
        let mut rec = make_bam_bytes(0, 0, 0, b"rea", &[], 0, -1, -1, aux);
        update_string_tag(&mut rec, SamTag::RX, b"newvalue");
        assert_eq!(find_string_tag_in_record(&rec, SamTag::RX), Some(b"newvalue".as_ref()));
    }

    #[test]
    fn test_update_string_tag_new() {
        let mut rec = make_bam_bytes(0, 0, 0, b"rea", &[], 0, -1, -1, &[]);
        update_string_tag(&mut rec, SamTag::RX, b"added");
        assert_eq!(find_string_tag_in_record(&rec, SamTag::RX), Some(b"added".as_ref()));
    }

    #[test]
    fn test_update_string_tag_same_length() {
        // Test same-length fast path (copy_from_slice, no splice)
        let aux = b"RXZold1\x00";
        let mut rec = make_bam_bytes(0, 0, 0, b"rea", &[], 0, -1, -1, aux);
        let orig_len = rec.len();
        update_string_tag(&mut rec, SamTag::RX, b"new2");
        // Record length should not change for same-length update
        assert_eq!(rec.len(), orig_len);
        assert_eq!(find_string_tag_in_record(&rec, SamTag::RX), Some(b"new2".as_ref()));
    }

    #[test]
    fn test_update_string_tag_different_length_shorter() {
        // Test splice path: old value longer than new value
        let aux = b"RXZlongvalue\x00";
        let mut rec = make_bam_bytes(0, 0, 0, b"rea", &[], 0, -1, -1, aux);
        update_string_tag(&mut rec, SamTag::RX, b"hi");
        assert_eq!(find_string_tag_in_record(&rec, SamTag::RX), Some(b"hi".as_ref()));
    }

    #[test]
    fn test_update_string_tag_different_length_longer() {
        // Test splice path: old value shorter than new value
        let aux = b"RXZhi\x00";
        let mut rec = make_bam_bytes(0, 0, 0, b"rea", &[], 0, -1, -1, aux);
        update_string_tag(&mut rec, SamTag::RX, b"longvalue");
        assert_eq!(find_string_tag_in_record(&rec, SamTag::RX), Some(b"longvalue".as_ref()));
    }

    #[test]
    fn test_update_string_tag_preserves_other_tags() {
        // Test that updating one tag doesn't corrupt adjacent tags
        let mut aux = Vec::new();
        aux.extend_from_slice(&[b'A', b'A', b'C', 1]); // AA:C:1
        aux.extend_from_slice(b"RXZold\x00"); // RX:Z:old
        aux.extend_from_slice(&[b'C', b'C', b'C', 2]); // CC:C:2
        let mut rec = make_bam_bytes(0, 0, 0, b"rea", &[], 0, -1, -1, &aux);
        update_string_tag(&mut rec, SamTag::RX, b"newval");
        let aux_start = aux_data_offset_from_record(&rec)
            .expect("record should have valid header for aux offset");
        assert_eq!(find_uint8_tag(&rec[aux_start..], b"AA"), Some(1));
        assert_eq!(find_string_tag(&rec[aux_start..], SamTag::RX), Some(b"newval".as_ref()));
        assert_eq!(find_uint8_tag(&rec[aux_start..], b"CC"), Some(2));
    }

    #[test]
    fn test_update_string_tag_no_aux_appends() {
        let mut rec = make_bam_bytes(0, 0, 0, b"rea", &[], 0, -1, -1, &[]);
        update_string_tag(&mut rec, SamTag::MI, b"42");
        assert_eq!(find_string_tag_in_record(&rec, SamTag::MI), Some(b"42".as_ref()));
    }

    // ========================================================================
    // append_int_tag tests
    // ========================================================================

    #[test]
    fn test_append_int_tag_i8() {
        let mut rec = Vec::new();
        append_int_tag(&mut rec, SamTag::CD, 42);
        assert_eq!(rec, [b'c', b'D', b'c', 42]);
    }

    #[test]
    fn test_append_int_tag_negative_i8() {
        let mut rec = Vec::new();
        append_int_tag(&mut rec, SamTag::CM, -5);
        assert_eq!(rec, [b'c', b'M', b'c', (-5i8).cast_unsigned()]);
    }

    #[test]
    fn test_append_int_tag_u8() {
        let mut rec = Vec::new();
        append_int_tag(&mut rec, SamTag::CD, 200);
        assert_eq!(rec, [b'c', b'D', b'C', 200]);
    }

    #[test]
    fn test_append_int_tag_negative_i16() {
        let mut rec = Vec::new();
        append_int_tag(&mut rec, SamTag::CD, -200);
        let v = (-200i16).to_le_bytes();
        assert_eq!(rec, [b'c', b'D', b's', v[0], v[1]]);
    }

    #[test]
    fn test_append_int_tag_u16() {
        let mut rec = Vec::new();
        append_int_tag(&mut rec, SamTag::CD, 1000);
        let v = 1000u16.to_le_bytes();
        assert_eq!(rec, [b'c', b'D', b'S', v[0], v[1]]);
    }

    #[test]
    fn test_append_int_tag_i32() {
        let mut rec = Vec::new();
        append_int_tag(&mut rec, SamTag::CD, 100_000);
        let v = 100_000i32.to_le_bytes();
        assert_eq!(rec, [b'c', b'D', b'i', v[0], v[1], v[2], v[3]]);
    }

    #[rstest]
    #[case::max_i8(127, b'c')]
    #[case::min_u8(128, b'C')]
    #[case::max_u8(255, b'C')]
    #[case::min_u16(256, b'S')]
    #[case::max_u16(65535, b'S')]
    #[case::min_i32(65536, b'i')]
    #[case::max_neg_i8(-128, b'c')]
    #[case::min_neg_i16(-129, b's')]
    #[case::min_i16(-32768, b's')]
    #[case::neg_i32(-32769, b'i')]
    #[case::i32_min(i32::MIN, b'i')]
    #[case::i32_max(i32::MAX, b'i')]
    fn test_append_int_tag_boundaries(#[case] value: i32, #[case] expected_type: u8) {
        let mut rec = Vec::new();
        append_int_tag(&mut rec, b"XX", value);
        assert_eq!(rec[2], expected_type);
    }

    // ========================================================================
    // append_float_tag tests
    // ========================================================================

    #[test]
    fn test_append_float_tag() {
        let mut rec = Vec::new();
        append_float_tag(&mut rec, SamTag::CE, 0.05);
        let v = 0.05f32.to_le_bytes();
        assert_eq!(rec, [b'c', b'E', b'f', v[0], v[1], v[2], v[3]]);
    }

    // ========================================================================
    // append_i16_array_tag tests
    // ========================================================================

    #[test]
    fn test_append_i16_array_tag_empty() {
        let mut rec = Vec::new();
        append_i16_array_tag(&mut rec, SamTag::CD_BASES, &[]);
        assert_eq!(rec, [b'c', b'd', b'B', b's', 0, 0, 0, 0]);
    }

    #[test]
    fn test_append_i16_array_tag_values() {
        let mut rec = Vec::new();
        append_i16_array_tag(&mut rec, SamTag::CD_BASES, &[10, 20, 5]);
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
        append_phred33_string_tag(&mut buf, SamTag::AQ, &[0, 10, 30, 40]);
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
    // find_array_tag tests
    // ========================================================================

    #[test]
    fn test_find_array_tag_int_array() {
        let aux = make_b_int_array_tag(SamTag::CD_BASES.into(), &[10, 20, 30]);
        let tag_ref = find_array_tag(&aux, SamTag::CD_BASES).expect("cd array tag should be found");
        assert_eq!(tag_ref.elem_type, b'i');
        assert_eq!(tag_ref.count, 3);
        assert_eq!(tag_ref.elem_size, 4);
        assert_eq!(tag_ref.data.len(), 12);
    }

    #[test]
    fn test_find_array_tag_uint8_array() {
        let aux = make_b_uint8_array_tag(*b"XC", &[1, 2, 3, 4]);
        let tag_ref = find_array_tag(&aux, b"XC").expect("XC array tag should be found");
        assert_eq!(tag_ref.elem_type, b'C');
        assert_eq!(tag_ref.count, 4);
        assert_eq!(tag_ref.elem_size, 1);
    }

    #[test]
    fn test_find_array_tag_not_found() {
        let aux = make_b_int_array_tag(SamTag::CD_BASES.into(), &[10]);
        assert!(find_array_tag(&aux, b"ZZ").is_none());
    }

    #[test]
    fn test_find_array_tag_not_b_type() {
        // Tag exists but as Z type, not B
        let aux = b"cdZhello\x00";
        assert!(find_array_tag(aux.as_ref(), SamTag::CD_BASES).is_none());
    }

    #[test]
    fn test_find_array_tag_after_other_tags() {
        let mut aux = Vec::new();
        aux.extend_from_slice(&[b'N', b'M', b'C', 5]); // NM:C:5
        aux.extend_from_slice(&make_b_int16_array_tag(SamTag::CD_BASES.into(), &[10, 20]));
        let tag_ref = find_array_tag(&aux, SamTag::CD_BASES).expect("cd array tag should be found");
        assert_eq!(tag_ref.elem_type, b's');
        assert_eq!(tag_ref.count, 2);
    }

    // ========================================================================
    // array_tag_element_u16 tests
    // ========================================================================

    #[test]
    fn test_array_tag_element_u16_uint8() {
        let aux = make_b_uint8_array_tag(SamTag::CD_BASES.into(), &[10, 200, 0]);
        let tag_ref = find_array_tag(&aux, SamTag::CD_BASES).expect("cd array tag should be found");
        assert_eq!(array_tag_element_u16(&tag_ref, 0), 10);
        assert_eq!(array_tag_element_u16(&tag_ref, 1), 200);
        assert_eq!(array_tag_element_u16(&tag_ref, 2), 0);
    }

    #[test]
    fn test_array_tag_element_u16_uint16() {
        let aux = make_b_uint16_array_tag(SamTag::CD_BASES.into(), &[100, 60000, 0]);
        let tag_ref = find_array_tag(&aux, SamTag::CD_BASES).expect("cd array tag should be found");
        assert_eq!(array_tag_element_u16(&tag_ref, 0), 100);
        assert_eq!(array_tag_element_u16(&tag_ref, 1), 60000);
        assert_eq!(array_tag_element_u16(&tag_ref, 2), 0);
    }

    #[test]
    fn test_array_tag_element_u16_int16_clamped() {
        let aux = make_b_int16_array_tag(SamTag::CD_BASES.into(), &[-5, 100, 0]);
        let tag_ref = find_array_tag(&aux, SamTag::CD_BASES).expect("cd array tag should be found");
        // Negative values clamped to 0
        assert_eq!(array_tag_element_u16(&tag_ref, 0), 0);
        assert_eq!(array_tag_element_u16(&tag_ref, 1), 100);
    }

    #[test]
    fn test_array_tag_element_u16_int8_clamped() {
        let aux = make_b_int8_array_tag(SamTag::CD_BASES.into(), &[-1, 42, 0]);
        let tag_ref = find_array_tag(&aux, SamTag::CD_BASES).expect("cd array tag should be found");
        assert_eq!(array_tag_element_u16(&tag_ref, 0), 0); // clamped
        assert_eq!(array_tag_element_u16(&tag_ref, 1), 42);
    }

    #[test]
    fn test_array_tag_element_u16_out_of_bounds() {
        let aux = make_b_uint8_array_tag(SamTag::CD_BASES.into(), &[10]);
        let tag_ref = find_array_tag(&aux, SamTag::CD_BASES).expect("cd array tag should be found");
        assert_eq!(array_tag_element_u16(&tag_ref, 1), 0); // out of bounds
    }

    #[test]
    fn test_array_tag_element_u16_unsupported_type() {
        // i32 array: not directly supported by array_tag_element_u16
        let aux = make_b_int_array_tag(SamTag::CD_BASES.into(), &[42]);
        let tag_ref = find_array_tag(&aux, SamTag::CD_BASES).expect("cd array tag should be found");
        assert_eq!(array_tag_element_u16(&tag_ref, 0), 0); // unsupported returns 0
    }

    // ========================================================================
    // array_tag_to_vec_u16 tests
    // ========================================================================

    #[test]
    fn test_array_tag_to_vec_u16() {
        let aux = make_b_uint8_array_tag(SamTag::CD_BASES.into(), &[10, 20, 30]);
        let tag_ref = find_array_tag(&aux, SamTag::CD_BASES).expect("cd array tag should be found");
        assert_eq!(array_tag_to_vec_u16(&tag_ref), vec![10u16, 20, 30]);
    }

    // ========================================================================
    // update_int_tag tests
    // ========================================================================

    #[test]
    fn test_update_int_tag_existing_i32_in_place() {
        // Create a record with NM:i:42 (4-byte int, in-place update)
        let mut rec = make_bam_bytes(0, 0, 0, b"rea", &[], 0, -1, -1, &[]);
        // Manually append an i32 tag
        rec.extend_from_slice(b"NMi");
        rec.extend_from_slice(&42i32.to_le_bytes());
        update_int_tag(&mut rec, SamTag::NM, 99);
        let aux = aux_data_slice(&rec);
        assert_eq!(find_int_tag(aux, SamTag::NM), Some(99));
    }

    #[test]
    fn test_update_int_tag_existing_different_size() {
        // Create a record with NM:c:5 (1-byte int), update to a value needing 2 bytes
        let mut rec = make_bam_bytes(0, 0, 0, b"rea", &[], 0, -1, -1, &[]);
        rec.extend_from_slice(&[b'N', b'M', b'c', 5]);
        update_int_tag(&mut rec, SamTag::NM, 300);
        let aux = aux_data_slice(&rec);
        assert_eq!(find_int_tag(aux, SamTag::NM), Some(300));
    }

    #[test]
    fn test_update_int_tag_absent() {
        let mut rec = make_bam_bytes(0, 0, 0, b"rea", &[], 0, -1, -1, &[]);
        update_int_tag(&mut rec, SamTag::NM, 42);
        let aux = aux_data_slice(&rec);
        assert_eq!(find_int_tag(aux, SamTag::NM), Some(42));
    }

    // ========================================================================
    // reverse_array_tag_in_place tests
    // ========================================================================

    #[test]
    fn test_reverse_array_tag_in_place_i16() {
        let mut aux = Vec::new();
        aux.extend_from_slice(&make_b_int16_array_tag(SamTag::CD_BASES.into(), &[10, 20, 30]));
        let mut rec = make_bam_bytes(0, 0, 0, b"rea", &[], 0, -1, -1, &aux);
        let aux_offset = aux_data_offset_from_record(&rec)
            .expect("record should have valid header for aux offset");
        reverse_array_tag_in_place(&mut rec, aux_offset, SamTag::CD_BASES);
        let tag_ref = find_array_tag(&rec[aux_offset..], SamTag::CD_BASES)
            .expect("cd array tag should be found in record");
        let values: Vec<i16> = (0..tag_ref.count)
            .map(|i| {
                let off = i * tag_ref.elem_size;
                i16::from_le_bytes([tag_ref.data[off], tag_ref.data[off + 1]])
            })
            .collect();
        assert_eq!(values, vec![30, 20, 10]);
    }

    #[test]
    fn test_reverse_array_tag_in_place_uint8() {
        let mut aux = Vec::new();
        aux.extend_from_slice(&make_b_uint8_array_tag(SamTag::CD_BASES.into(), &[1, 2, 3, 4]));
        let mut rec = make_bam_bytes(0, 0, 0, b"rea", &[], 0, -1, -1, &aux);
        let aux_offset = aux_data_offset_from_record(&rec)
            .expect("record should have valid header for aux offset");
        reverse_array_tag_in_place(&mut rec, aux_offset, SamTag::CD_BASES);
        let tag_ref = find_array_tag(&rec[aux_offset..], SamTag::CD_BASES)
            .expect("cd array tag should be found in record");
        assert_eq!(tag_ref.data, &[4, 3, 2, 1]);
    }

    #[test]
    fn test_reverse_array_tag_in_place_not_found() {
        let mut rec = make_bam_bytes(0, 0, 0, b"rea", &[], 0, -1, -1, &[]);
        let aux_offset = aux_data_offset_from_record(&rec)
            .expect("record should have valid header for aux offset");
        // Should be a no-op
        reverse_array_tag_in_place(&mut rec, aux_offset, SamTag::CD_BASES);
    }

    #[test]
    fn test_reverse_array_tag_in_place_offset_past_end() {
        let mut rec = make_bam_bytes(0, 0, 0, b"rea", &[], 0, -1, -1, &[]);
        // aux_offset >= record.len() should be a no-op
        let offset = rec.len() + 10;
        reverse_array_tag_in_place(&mut rec, offset, SamTag::CD_BASES);
    }

    // ========================================================================
    // reverse_string_tag_in_place tests
    // ========================================================================

    #[test]
    fn test_reverse_string_tag_in_place() {
        let aux = b"RXZhello\x00";
        let mut rec = make_bam_bytes(0, 0, 0, b"rea", &[], 0, -1, -1, aux);
        let aux_offset = aux_data_offset_from_record(&rec)
            .expect("record should have valid header for aux offset");
        reverse_string_tag_in_place(&mut rec, aux_offset, SamTag::RX);
        assert_eq!(find_string_tag(&rec[aux_offset..], SamTag::RX), Some(b"olleh".as_ref()));
    }

    #[test]
    fn test_reverse_string_tag_in_place_single_char() {
        let aux = b"RXZa\x00";
        let mut rec = make_bam_bytes(0, 0, 0, b"rea", &[], 0, -1, -1, aux);
        let aux_offset = aux_data_offset_from_record(&rec)
            .expect("record should have valid header for aux offset");
        reverse_string_tag_in_place(&mut rec, aux_offset, SamTag::RX);
        assert_eq!(find_string_tag(&rec[aux_offset..], SamTag::RX), Some(b"a".as_ref()));
    }

    #[test]
    fn test_reverse_string_tag_in_place_not_found() {
        let mut rec = make_bam_bytes(0, 0, 0, b"rea", &[], 0, -1, -1, &[]);
        let aux_offset = aux_data_offset_from_record(&rec)
            .expect("record should have valid header for aux offset");
        // No-op
        reverse_string_tag_in_place(&mut rec, aux_offset, SamTag::RX);
    }

    #[test]
    fn test_reverse_string_tag_in_place_offset_past_end() {
        let mut rec = make_bam_bytes(0, 0, 0, b"rea", &[], 0, -1, -1, &[]);
        let offset = rec.len() + 10;
        reverse_string_tag_in_place(&mut rec, offset, SamTag::RX);
    }

    // ========================================================================
    // reverse_complement_string_tag_in_place tests
    // ========================================================================

    #[test]
    fn test_reverse_complement_string_tag_in_place() {
        let aux = b"RXZACGT\x00";
        let mut rec = make_bam_bytes(0, 0, 0, b"rea", &[], 0, -1, -1, aux);
        let aux_offset = aux_data_offset_from_record(&rec)
            .expect("record should have valid header for aux offset");
        reverse_complement_string_tag_in_place(&mut rec, aux_offset, SamTag::RX);
        assert_eq!(find_string_tag(&rec[aux_offset..], SamTag::RX), Some(b"ACGT".as_ref()));
    }

    #[test]
    fn test_reverse_complement_string_tag_in_place_lowercase() {
        let aux = b"RXZacgt\x00";
        let mut rec = make_bam_bytes(0, 0, 0, b"rea", &[], 0, -1, -1, aux);
        let aux_offset = aux_data_offset_from_record(&rec)
            .expect("record should have valid header for aux offset");
        reverse_complement_string_tag_in_place(&mut rec, aux_offset, SamTag::RX);
        // case-preserving: acgt -> reverse tgca -> complement acgt
        assert_eq!(find_string_tag(&rec[aux_offset..], SamTag::RX), Some(b"acgt".as_ref()));
    }

    #[test]
    fn test_reverse_complement_string_tag_in_place_iupac_and_case() {
        // IUPAC-aware + case-preserving, sharing fgumi-dna's COMPLEMENT table (R2-SEQ-04).
        // "NRG" -> reverse "GRN" -> complement "CYN"
        let aux = b"RXZNRG\x00";
        let mut rec = make_bam_bytes(0, 0, 0, b"rea", &[], 0, -1, -1, aux);
        let aux_offset = aux_data_offset_from_record(&rec)
            .expect("record should have valid header for aux offset");
        reverse_complement_string_tag_in_place(&mut rec, aux_offset, SamTag::RX);
        assert_eq!(find_string_tag(&rec[aux_offset..], SamTag::RX), Some(b"CYN".as_ref()));
    }

    #[test]
    fn test_reverse_complement_string_tag_in_place_with_n() {
        let aux = b"RXZANGT\x00";
        let mut rec = make_bam_bytes(0, 0, 0, b"rea", &[], 0, -1, -1, aux);
        let aux_offset = aux_data_offset_from_record(&rec)
            .expect("record should have valid header for aux offset");
        reverse_complement_string_tag_in_place(&mut rec, aux_offset, SamTag::RX);
        // ANGT -> reverse TGNA -> complement ACNT
        assert_eq!(find_string_tag(&rec[aux_offset..], SamTag::RX), Some(b"ACNT".as_ref()));
    }

    #[test]
    fn test_reverse_complement_string_tag_in_place_not_found() {
        let mut rec = make_bam_bytes(0, 0, 0, b"rea", &[], 0, -1, -1, &[]);
        let aux_offset = aux_data_offset_from_record(&rec)
            .expect("record should have valid header for aux offset");
        reverse_complement_string_tag_in_place(&mut rec, aux_offset, SamTag::RX);
    }

    #[test]
    fn test_reverse_complement_string_tag_in_place_offset_past_end() {
        let mut rec = make_bam_bytes(0, 0, 0, b"rea", &[], 0, -1, -1, &[]);
        let offset = rec.len() + 10;
        reverse_complement_string_tag_in_place(&mut rec, offset, SamTag::RX);
    }

    // ========================================================================
    // extract_template_aux_tags tests
    // ========================================================================

    #[test]
    fn test_extract_template_aux_tags_all_found() {
        // Build a BAM record with MI, RG, CB, and MC tags
        let mut aux = Vec::new();
        aux.extend_from_slice(b"MIZ42\x00");
        aux.extend_from_slice(b"RGZsample1\x00");
        aux.extend_from_slice(b"CBZcell99\x00");
        aux.extend_from_slice(b"MCZ10M5S\x00");
        let rec = make_bam_bytes(0, 100, 0, b"read1", &[], 4, -1, -1, &aux);

        let result = extract_template_aux_tags(&rec, Some(SamTag::CB));
        assert_eq!(result.mi, (42, true));
        assert_eq!(result.rg, Some(b"sample1".as_ref()));
        assert_eq!(result.cell, Some(b"cell99".as_ref()));
        assert_eq!(result.mc, Some(b"10M5S".as_slice()));
    }

    #[test]
    fn test_extract_template_aux_tags_no_cell_tag() {
        let mut aux = Vec::new();
        aux.extend_from_slice(b"MIZ7/B\x00");
        aux.extend_from_slice(b"RGZlib1\x00");
        aux.extend_from_slice(b"MCZ20M\x00");
        let rec = make_bam_bytes(0, 0, 0, b"r1", &[], 4, -1, -1, &aux);

        let result = extract_template_aux_tags(&rec, None);
        assert_eq!(result.mi, (7, false));
        assert_eq!(result.rg, Some(b"lib1".as_ref()));
        assert!(result.cell.is_none());
        assert_eq!(result.mc, Some(b"20M".as_slice()));
    }

    #[test]
    fn test_extract_template_aux_tags_mi_integer_type() {
        // MI as integer tag (C type = u8)
        let mut aux = Vec::new();
        aux.extend_from_slice(&[b'M', b'I', b'C', 99]);
        aux.extend_from_slice(b"RGZrg0\x00");
        let rec = make_bam_bytes(0, 0, 0, b"r1", &[], 4, -1, -1, &aux);

        let result = extract_template_aux_tags(&rec, None);
        assert_eq!(result.mi, (99, true));
        assert_eq!(result.rg, Some(b"rg0".as_ref()));
    }

    #[test]
    fn test_extract_template_aux_tags_partial() {
        // Only RG present
        let aux = b"RGZsample\x00";
        let rec = make_bam_bytes(0, 0, 0, b"r1", &[], 4, -1, -1, aux);

        let result = extract_template_aux_tags(&rec, Some(SamTag::CB));
        assert_eq!(result.mi, (0, true)); // default
        assert_eq!(result.rg, Some(b"sample".as_ref()));
        assert!(result.cell.is_none());
        assert!(result.mc.is_none());
    }

    #[test]
    fn test_extract_template_aux_tags_empty_aux() {
        let rec = make_bam_bytes(0, 0, 0, b"r1", &[], 4, -1, -1, &[]);

        let result = extract_template_aux_tags(&rec, Some(SamTag::CB));
        assert_eq!(result.mi, (0, true));
        assert!(result.rg.is_none());
        assert!(result.cell.is_none());
        assert!(result.mc.is_none());
    }

    #[test]
    fn test_extract_template_aux_tags_with_non_string_tags() {
        // Non-string tags interspersed with target tags
        let mut aux = Vec::new();
        aux.extend_from_slice(&[b'N', b'M', b'C', 5]); // NM:C:5
        aux.extend_from_slice(b"MIZ100\x00");
        aux.extend_from_slice(&[b'A', b'S', b'C', 30]); // AS:C:30
        aux.extend_from_slice(b"RGZlib2\x00");
        aux.extend_from_slice(b"MCZ5M\x00");
        let rec = make_bam_bytes(0, 0, 0, b"r1", &[], 4, -1, -1, &aux);

        let result = extract_template_aux_tags(&rec, None);
        assert_eq!(result.mi, (100, true));
        assert_eq!(result.rg, Some(b"lib2".as_ref()));
        assert_eq!(result.mc, Some(b"5M".as_slice()));
    }

    // ========================================================================
    // Dedup pa-tag validation bug reproduction
    // ========================================================================

    #[test]
    fn test_dedup_pa_tag_check_fails_on_b_array() {
        // The legacy `pa:B:i,0,27_056_961,0,207,60005,1` sort key that `fgumi zipper` wrote
        // before 0.2.0, when the tag was renamed to `tc`. The name is incidental to the bug —
        // what defeats the dedup check is the `B:i` array type, which `tc` still uses.
        let aux = make_b_int_array_tag(*b"pa", &[0, 27_056_961, 0, 207, 60005, 1]);
        let pa_tag_bytes: [u8; 2] = *b"pa";

        // This is the exact check from dedup.rs:932-934
        let found_by_dedup = find_string_tag(&aux, pa_tag_bytes).is_some()
            || find_int_tag(&aux, pa_tag_bytes).is_some();

        // BUG: dedup thinks the pa tag is missing even though it's present
        assert!(!found_by_dedup, "dedup check should fail to find B:i pa tag");

        // But find_tag_type correctly finds it
        assert!(find_tag_type(&aux, pa_tag_bytes).is_some());
    }

    // ========================================================================
    // append_i32_array_tag tests
    // ========================================================================

    #[test]
    fn test_append_i32_array_tag() {
        let mut rec = make_bam_bytes(0, 0, 0, b"r1", &[], 4, -1, -1, &[]);
        let values = [1i32, -200, 300_000];
        append_i32_array_tag(&mut rec, b"pa", &values);

        let aux = aux_data_slice(&rec);
        let tag_type = find_tag_type(aux, b"pa");
        assert_eq!(tag_type, Some(b'B'));

        // Verify the array contents
        let arr = find_array_tag(aux, b"pa").unwrap();
        assert_eq!(arr.count, 3);
        assert_eq!(arr.elem_type, b'i');
    }

    #[test]
    fn test_append_i32_array_tag_empty() {
        let mut rec = make_bam_bytes(0, 0, 0, b"r1", &[], 4, -1, -1, &[]);
        append_i32_array_tag(&mut rec, b"pa", &[]);

        let aux = aux_data_slice(&rec);
        let arr = find_array_tag(aux, b"pa").unwrap();
        assert_eq!(arr.count, 0);
    }

    // ========================================================================
    // normalize_int_tag_to_smallest_signed tests
    // ========================================================================

    #[test]
    fn test_normalize_int_tag_i32_to_i8() {
        // AS:i:77 (stored as i32) should normalize to AS:c:77 (i8)
        let mut aux = Vec::new();
        aux.extend_from_slice(b"ASi");
        aux.extend_from_slice(&77i32.to_le_bytes());
        let mut rec = make_bam_bytes(0, 0, 0, b"r1", &[], 4, -1, -1, &aux);

        normalize_int_tag_to_smallest_signed(&mut rec, SamTag::AS);

        let aux_data = aux_data_slice(&rec);
        let (_, val_type) = find_tag_position(aux_data, SamTag::AS.into()).unwrap();
        assert_eq!(val_type, b'c', "77 should fit in i8 (type 'c')");
        assert_eq!(find_int_tag(aux_data, SamTag::AS), Some(77));
    }

    #[test]
    fn test_normalize_int_tag_i32_to_i16() {
        // AS:i:200 (stored as i32) should normalize to AS:s:200 (i16, signed-only encoding)
        let mut aux = Vec::new();
        aux.extend_from_slice(b"ASi");
        aux.extend_from_slice(&200i32.to_le_bytes());
        let mut rec = make_bam_bytes(0, 0, 0, b"r1", &[], 4, -1, -1, &aux);

        normalize_int_tag_to_smallest_signed(&mut rec, SamTag::AS);

        let aux_data = aux_data_slice(&rec);
        let (_, val_type) = find_tag_position(aux_data, SamTag::AS.into()).unwrap();
        assert_eq!(val_type, b's', "200 should be i16 (type 's') with signed-only encoding");
        assert_eq!(find_int_tag(aux_data, SamTag::AS), Some(200));
    }

    #[test]
    fn test_normalize_int_tag_preserves_large_value() {
        // AS:i:100000 should remain i32
        let mut aux = Vec::new();
        aux.extend_from_slice(b"ASi");
        aux.extend_from_slice(&100_000i32.to_le_bytes());
        let mut rec = make_bam_bytes(0, 0, 0, b"r1", &[], 4, -1, -1, &aux);

        normalize_int_tag_to_smallest_signed(&mut rec, SamTag::AS);

        let aux_data = aux_data_slice(&rec);
        let (_, val_type) = find_tag_position(aux_data, SamTag::AS.into()).unwrap();
        assert_eq!(val_type, b'i', "100000 requires i32 (type 'i')");
        assert_eq!(find_int_tag(aux_data, SamTag::AS), Some(100_000));
    }

    #[test]
    fn test_normalize_int_tag_missing_is_noop() {
        let mut rec = make_bam_bytes(0, 0, 0, b"r1", &[], 4, -1, -1, &[]);
        let original = rec.clone();
        normalize_int_tag_to_smallest_signed(&mut rec, SamTag::AS);
        assert_eq!(rec, original);
    }

    // ========================================================================
    // copy_aux_tags tests
    // ========================================================================

    #[test]
    fn test_copy_aux_tags_all() {
        // Source has RX:Z:ACGT and NM:C:5
        let mut src_aux = Vec::new();
        src_aux.extend_from_slice(b"RXZ\x41\x43\x47\x54\x00"); // RX:Z:ACGT
        src_aux.extend_from_slice(&[b'N', b'M', b'C', 5]); // NM:C:5

        let mut dest = Vec::new();
        copy_aux_tags(&src_aux, &mut dest, &[]);

        assert_eq!(dest, src_aux, "All tags should be copied when no skip list");
    }

    #[test]
    fn test_copy_aux_tags_with_skip() {
        // Source has RX:Z:ACGT and NM:C:5
        let mut src_aux = Vec::new();
        src_aux.extend_from_slice(b"RXZ\x41\x43\x47\x54\x00"); // RX:Z:ACGT
        src_aux.extend_from_slice(&[b'N', b'M', b'C', 5]); // NM:C:5

        let mut dest = Vec::new();
        copy_aux_tags(&src_aux, &mut dest, &[SamTag::NM]);

        // Only RX should be copied
        assert_eq!(dest, b"RXZ\x41\x43\x47\x54\x00");
    }

    #[test]
    fn test_copy_aux_tags_skip_all() {
        let mut src_aux = Vec::new();
        src_aux.extend_from_slice(b"RXZ\x41\x43\x47\x54\x00");
        src_aux.extend_from_slice(&[b'N', b'M', b'C', 5]);

        let mut dest = Vec::new();
        copy_aux_tags(&src_aux, &mut dest, &[SamTag::RX, SamTag::NM]);

        assert!(dest.is_empty(), "All tags skipped should produce empty dest");
    }

    #[test]
    fn test_copy_aux_tags_empty_source() {
        let mut dest = Vec::new();
        copy_aux_tags(&[], &mut dest, &[]);
        assert!(dest.is_empty());
    }

    #[rstest]
    #[case::empty(&[])]
    #[case::single(&[42u8])]
    #[case::multiple(&[0u8, 128, 255])]
    fn test_append_u8_array_tag_round_trip(#[case] values: &[u8]) {
        let tag = SamTag::ML;
        let mut record = Vec::new();
        append_u8_array_tag(&mut record, tag, values);

        // Verify wire format: [tag0, tag1, 'B', 'C', count_u32_le, values...]
        assert_eq!(record[0], b'M');
        assert_eq!(record[1], b'L');
        assert_eq!(record[2], b'B');
        assert_eq!(record[3], b'C');
        let count = u32::from_le_bytes([record[4], record[5], record[6], record[7]]) as usize;
        assert_eq!(count, values.len());
        assert_eq!(&record[8..], values);

        // Round-trip through find_array_tag
        let arr = find_array_tag(&record, tag).expect("tag should be found");
        assert_eq!(arr.elem_type, b'C');
        assert_eq!(arr.count, values.len());
        assert_eq!(arr.data, values);
    }

    #[test]
    fn test_raw_tags_view_construction_and_find_string() {
        use crate::fields::RawRecordView;
        let aux = b"RGZmysample\0";
        let rec = make_bam_bytes(0, 0, 0, b"r", &[], 0, -1, -1, aux);
        let v = RawRecordView::new(&rec);
        let tags = v.tags();
        assert_eq!(tags.find_string(SamTag::RG), Some(b"mysample".as_slice()));
        assert!(!tags.is_empty());
        assert!(tags.contains(SamTag::RG));
        assert!(!tags.contains(SamTag::NM));
    }

    #[test]
    fn test_raw_tags_iter_basic() {
        use crate::fields::RawRecordView;
        // RG:Z:rg1\0 NM:i:5 (4 bytes LE)
        let mut aux: Vec<u8> = b"RGZrg1\0".to_vec();
        aux.extend_from_slice(b"NMi");
        aux.extend_from_slice(&5i32.to_le_bytes());
        let rec = make_bam_bytes(0, 0, 0, b"r", &[], 0, -1, -1, &aux);

        let v = RawRecordView::new(&rec);
        let mut entries: Vec<(String, u8, Vec<u8>)> = v
            .tags()
            .iter()
            .map(|e| {
                (String::from_utf8(e.tag.to_vec()).unwrap(), e.type_byte, e.value_bytes.to_vec())
            })
            .collect();
        entries.sort();
        assert_eq!(
            entries,
            vec![
                ("NM".into(), b'i', 5i32.to_le_bytes().to_vec()),
                ("RG".into(), b'Z', b"rg1\0".to_vec()),
            ]
        );
    }

    #[test]
    fn test_raw_tags_extract_string_batch() {
        use crate::fields::RawRecordView;
        let aux = b"RGZmygrp\0BCZACGT\0MCZ50M\0";
        let rec = make_bam_bytes(0, 0, 0, b"r", &[], 0, -1, -1, aux);
        let s = RawRecordView::new(&rec).tags().extract_string_batch(SamTag::BC, None);
        assert_eq!(s.rg, Some(b"mygrp".as_slice()));
        assert_eq!(s.cell, Some(b"ACGT".as_slice()));
        assert_eq!(s.mc, Some(b"50M".as_slice()));
    }

    #[test]
    fn test_mc_is_returned_as_bytes_without_a_utf8_gate() {
        // `MC` is only ever consumed as CIGAR bytes, so extraction hands the raw
        // value over rather than validating it as UTF-8 first. A value with a
        // valid CIGAR prefix and a non-UTF-8 byte after it therefore reaches the
        // parser, which reads the prefix and stops -- where the earlier
        // `str::from_utf8(..).ok()` gate discarded the whole tag and silently
        // left the mate position unadjusted.
        let aux = b"MCZ10S40M\xff\0RGZrg1\0";
        let rec = make_bam_bytes(0, 0, 0, b"r", &[], 0, -1, -1, aux);
        let tags = extract_template_aux_tags(&rec, None);
        assert_eq!(tags.mc, Some(b"10S40M\xff".as_slice()));
        assert_eq!(tags.rg, Some(b"rg1".as_slice()), "the scan continues past MC");
    }

    #[test]
    fn test_mc_bytes_round_trip_through_the_clip_parser() {
        // The point of dropping the gate is that the bytes still parse, so pin
        // the value the sort key actually uses rather than only the field.
        let aux = b"MCZ10S40M\0";
        let rec = make_bam_bytes(0, 0, 0, b"r", &[], 0, -1, -1, aux);
        let mc = extract_template_aux_tags(&rec, None).mc.expect("MC present");
        assert_eq!(crate::cigar::mate_unclipped_5prime(100, false, mc), 90);
    }

    #[test]
    fn test_raw_tags_mut_set_array_element() {
        use crate::fields::{RawRecordMut, RawRecordView};
        // Build a record with B:S array tag "bq" of [10, 20, 30, 40]
        let mut aux = Vec::new();
        aux.extend_from_slice(b"bqBS");
        aux.extend_from_slice(&4u32.to_le_bytes());
        for v in [10u16, 20, 30, 40] {
            aux.extend_from_slice(&v.to_le_bytes());
        }
        let mut rec = make_bam_bytes(0, 0, 0, b"r", &[], 0, -1, -1, &aux);

        {
            let mut m = RawRecordMut::new(&mut rec);
            let mut tm = m.tags_mut();
            tm.set_array_element_u16(SamTag::BQ, 2, 99);
        }
        let v = RawRecordView::new(&rec);
        let arr = v.tags().find_array(SamTag::BQ).expect("array tag");
        assert_eq!(arr.count, 4);
        assert_eq!(array_tag_element_u16(&arr, 2), 99);
    }

    // ========================================================================
    // RawTagsMut reverse + in-place setter tests
    // ========================================================================

    #[test]
    fn test_raw_tags_mut_reverse_string_in_place() {
        use crate::fields::{RawRecordMut, RawRecordView};
        let aux = b"BCZACGT\0";
        let mut rec = make_bam_bytes(0, 0, 0, b"r", &[], 0, -1, -1, aux);
        {
            let mut m = RawRecordMut::new(&mut rec);
            m.tags_mut().reverse_string(SamTag::BC);
        }
        assert_eq!(
            RawRecordView::new(&rec).tags().find_string(SamTag::BC),
            Some(b"TGCA".as_slice())
        );
    }

    #[test]
    fn test_raw_tags_mut_set_string_in_place_same_length() {
        use crate::fields::{RawRecordMut, RawRecordView};
        let aux = b"BCZACGT\0";
        let mut rec = make_bam_bytes(0, 0, 0, b"r", &[], 0, -1, -1, aux);
        let ok = {
            let mut m = RawRecordMut::new(&mut rec);
            m.tags_mut().set_string_in_place(SamTag::BC, b"TTTT")
        };
        assert!(ok);
        assert_eq!(
            RawRecordView::new(&rec).tags().find_string(SamTag::BC),
            Some(b"TTTT".as_slice())
        );
    }

    #[test]
    fn test_raw_tags_mut_set_string_in_place_different_length_returns_false() {
        use crate::fields::RawRecordMut;
        let aux = b"BCZACGT\0";
        let mut rec = make_bam_bytes(0, 0, 0, b"r", &[], 0, -1, -1, aux);
        let ok = {
            let mut m = RawRecordMut::new(&mut rec);
            m.tags_mut().set_string_in_place(SamTag::BC, b"AC")
        };
        assert!(!ok);
    }

    #[test]
    fn test_raw_tags_mut_set_int_in_place_fits() {
        use crate::fields::{RawRecordMut, RawRecordView};
        // NM:i:5 (4-byte signed int)
        let mut aux = Vec::from(b"NMi".as_slice());
        aux.extend_from_slice(&5i32.to_le_bytes());
        let mut rec = make_bam_bytes(0, 0, 0, b"r", &[], 0, -1, -1, &aux);
        let ok = {
            let mut m = RawRecordMut::new(&mut rec);
            m.tags_mut().set_int_in_place(SamTag::NM, 100_000)
        };
        assert!(ok);
        assert_eq!(RawRecordView::new(&rec).tags().find_int(SamTag::NM), Some(100_000));
    }

    #[test]
    fn test_raw_tags_mut_set_int_in_place_doesnt_fit_returns_false() {
        use crate::fields::RawRecordMut;
        // NM:c:5 (1-byte signed)
        let aux = b"NMc\x05";
        let mut rec = make_bam_bytes(0, 0, 0, b"r", &[], 0, -1, -1, aux);
        let ok = {
            let mut m = RawRecordMut::new(&mut rec);
            m.tags_mut().set_int_in_place(SamTag::NM, 100_000)
        };
        assert!(!ok);
    }

    #[test]
    fn test_raw_tags_mut_set_int_in_place_truncated_returns_false_not_panic() {
        use crate::fields::RawRecordMut;
        // Aux declares NM as type 'i' (needs 4 bytes) but only carries 2 bytes
        // of payload. The in-place setter must return false and NOT panic on
        // the out-of-bounds slice index.
        let aux = b"NMi\x01\x02"; // truncated i32
        let mut rec = make_bam_bytes(0, 0, 0, b"r", &[], 0, -1, -1, aux);
        let ok = {
            let mut m = RawRecordMut::new(&mut rec);
            m.tags_mut().set_int_in_place(SamTag::NM, 7)
        };
        assert!(!ok, "truncated int tag must return false");
    }

    #[test]
    fn test_raw_tags_mut_set_float_in_place_truncated_returns_false_not_panic() {
        use crate::fields::RawRecordMut;
        // Aux declares AS as type 'f' but only carries 2 bytes of payload.
        let aux = b"ASf\x01\x02";
        let mut rec = make_bam_bytes(0, 0, 0, b"r", &[], 0, -1, -1, aux);
        let ok = {
            let mut m = RawRecordMut::new(&mut rec);
            m.tags_mut().set_float_in_place(SamTag::AS, 1.0)
        };
        assert!(!ok, "truncated float tag must return false");
    }

    #[test]
    fn test_raw_tags_mut_set_float_in_place() {
        use crate::fields::{RawRecordMut, RawRecordView};
        let mut aux = Vec::from(b"ASf".as_slice());
        aux.extend_from_slice(&12.5f32.to_le_bytes());
        let mut rec = make_bam_bytes(0, 0, 0, b"r", &[], 0, -1, -1, &aux);
        let ok = {
            let mut m = RawRecordMut::new(&mut rec);
            m.tags_mut().set_float_in_place(SamTag::AS, 99.25)
        };
        assert!(ok);
        let got = RawRecordView::new(&rec).tags().find_float(SamTag::AS).unwrap();
        assert!((got - 99.25).abs() < 1e-6);
    }

    // ========================================================================
    // RawTagsEditor tests
    // ========================================================================

    #[test]
    fn test_raw_tags_editor_from_vec_caches_aux_offset() {
        let aux = b"RGZmygrp\0";
        let mut rec = make_bam_bytes(0, 0, 0, b"r", &[], 0, -1, -1, aux);
        let expected_off = aux_data_offset_from_record(&rec).unwrap();
        let editor = RawTagsEditor::from_vec(&mut rec);
        assert_eq!(editor.aux_offset(), expected_off);
        assert_eq!(editor.view().find_string(SamTag::RG), Some(b"mygrp".as_slice()));
    }

    #[test]
    fn test_editor_append_remove_update_int_string_roundtrip() {
        use crate::fields::RawRecordView;
        let mut rec = make_bam_bytes(0, 0, 0, b"r", &[], 0, -1, -1, &[]);
        {
            let mut ed = RawTagsEditor::from_vec(&mut rec);
            ed.append_string(SamTag::RG, b"mygrp");
            ed.append_int(SamTag::NM, 5);
            ed.update_int(SamTag::NM, 7);
            ed.update_string(SamTag::RG, b"newgrp");
        }
        let v = RawRecordView::new(&rec);
        assert_eq!(v.tags().find_string(SamTag::RG), Some(b"newgrp".as_slice()));
        assert_eq!(v.tags().find_int(SamTag::NM), Some(7));

        {
            let mut ed = RawTagsEditor::from_vec(&mut rec);
            ed.remove(SamTag::NM);
        }
        assert_eq!(RawRecordView::new(&rec).tags().find_int(SamTag::NM), None);
    }

    #[test]
    fn test_editor_normalize_int_to_smallest_signed() {
        let mut rec = make_bam_bytes(0, 0, 0, b"r", &[], 0, -1, -1, &[]);
        {
            let mut ed = RawTagsEditor::from_vec(&mut rec);
            ed.append_int(SamTag::NM, 100_000); // Will encode as 'i' (i32) — fits 100,000
            ed.normalize_int_to_smallest_signed(SamTag::NM);
        }
        let aux = aux_data_slice(&rec);
        let (p, t) = find_tag_position(aux, SamTag::NM.into()).unwrap();
        assert_eq!(t, b'i');
        assert_eq!(i32::from_le_bytes([aux[p + 3], aux[p + 4], aux[p + 5], aux[p + 6]]), 100_000);
    }

    #[test]
    fn test_editor_append_phred33_string() {
        use crate::fields::RawRecordView;
        let mut rec = make_bam_bytes(0, 0, 0, b"r", &[], 0, -1, -1, &[]);
        {
            let mut ed = RawTagsEditor::from_vec(&mut rec);
            ed.append_phred33_string(SamTag::OQ, &[30, 31, 32, 33]);
        }
        // Phred+33: 30 -> '?', 31 -> '@', 32 -> 'A', 33 -> 'B'
        assert_eq!(
            RawRecordView::new(&rec).tags().find_string(SamTag::OQ),
            Some(b"?@AB".as_slice())
        );
    }

    #[test]
    fn test_editor_append_array_i16() {
        use crate::fields::RawRecordView;
        let mut rec = make_bam_bytes(0, 0, 0, b"r", &[], 0, -1, -1, &[]);
        {
            let mut ed = RawTagsEditor::from_vec(&mut rec);
            ed.append_array_i16(b"sq", &[-100, 0, 100]);
        }
        let arr = RawRecordView::new(&rec).tags().find_array(b"sq").expect("present");
        assert_eq!(arr.count, 3);
        assert_eq!(arr.elem_type, b's');
    }

    #[test]
    fn test_editor_copy_from_skip() {
        use crate::fields::RawRecordView;
        let src_aux = b"RGZmygrp\0NMc\x05ASc\x0a";
        let src_rec = make_bam_bytes(0, 0, 0, b"r", &[], 0, -1, -1, src_aux);
        let mut dst_rec = make_bam_bytes(0, 0, 0, b"r", &[], 0, -1, -1, &[]);
        {
            let src_view = RawRecordView::new(&src_rec);
            let mut ed = RawTagsEditor::from_vec(&mut dst_rec);
            ed.copy_from(src_view.tags(), &[SamTag::NM]);
        }
        let dst = RawRecordView::new(&dst_rec);
        assert_eq!(dst.tags().find_string(SamTag::RG), Some(b"mygrp".as_slice()));
        assert_eq!(dst.tags().find_int(SamTag::NM), None); // skipped
        assert_eq!(dst.tags().find_int(SamTag::AS), Some(10));
    }

    // ========================================================================
    // update_float + update_array_* tests
    // ========================================================================

    #[test]
    fn test_editor_update_float() {
        use crate::fields::RawRecordView;
        let mut rec = make_bam_bytes(0, 0, 0, b"r", &[], 0, -1, -1, &[]);
        {
            let mut ed = RawTagsEditor::from_vec(&mut rec);
            ed.append_float(SamTag::AS, 12.5);
            ed.update_float(SamTag::AS, 99.25);
        }
        assert!(
            (RawRecordView::new(&rec).tags().find_float(SamTag::AS).unwrap() - 99.25).abs() < 1e-6
        );

        // Updating a non-existent float tag inserts it
        let mut rec2 = make_bam_bytes(0, 0, 0, b"r", &[], 0, -1, -1, &[]);
        {
            let mut ed = RawTagsEditor::from_vec(&mut rec2);
            ed.update_float(b"BB", 1.5);
        }
        assert!((RawRecordView::new(&rec2).tags().find_float(b"BB").unwrap() - 1.5).abs() < 1e-6);
    }

    #[test]
    fn test_editor_update_array_u16_same_length() {
        use crate::fields::RawRecordView;
        let mut rec = make_bam_bytes(0, 0, 0, b"r", &[], 0, -1, -1, &[]);
        {
            let mut ed = RawTagsEditor::from_vec(&mut rec);
            ed.update_array_u16(SamTag::BQ, &[0u16, 1, 2, 3]); // not present -> append
            ed.update_array_u16(SamTag::BQ, &[10, 20, 30, 40]); // same length -> in-place
        }
        let arr = RawRecordView::new(&rec).tags().find_array(SamTag::BQ).expect("present");
        let vals = array_tag_to_vec_u16(&arr);
        assert_eq!(vals, vec![10, 20, 30, 40]);
    }

    #[test]
    fn test_editor_update_array_u16_different_length_splices() {
        use crate::fields::RawRecordView;
        let mut rec = make_bam_bytes(0, 0, 0, b"r", &[], 0, -1, -1, &[]);
        {
            let mut ed = RawTagsEditor::from_vec(&mut rec);
            ed.update_array_u16(SamTag::BQ, &[1u16, 2, 3]); // not present -> append
            ed.update_array_u16(SamTag::BQ, &[7u16, 8, 9, 10, 11]); // grows
        }
        let arr = RawRecordView::new(&rec).tags().find_array(SamTag::BQ).expect("present");
        let vals = array_tag_to_vec_u16(&arr);
        assert_eq!(vals, vec![7, 8, 9, 10, 11]);
    }

    #[test]
    fn test_editor_update_array_i32_grow_then_in_place() {
        use crate::fields::RawRecordView;
        let mut rec = make_bam_bytes(0, 0, 0, b"r", &[], 0, -1, -1, &[]);
        {
            let mut ed = RawTagsEditor::from_vec(&mut rec);
            ed.update_array_i32(b"sc", &[100, 200, 300]);
            // Same-length in-place
            ed.update_array_i32(b"sc", &[-1, -2, -3]);
        }
        let arr = RawRecordView::new(&rec).tags().find_array(b"sc").expect("present");
        assert_eq!(arr.elem_type, b'i');
        assert_eq!(arr.count, 3);
    }

    #[test]
    fn test_editor_update_array_f32() {
        use crate::fields::RawRecordView;
        let mut rec = make_bam_bytes(0, 0, 0, b"r", &[], 0, -1, -1, &[]);
        {
            let mut ed = RawTagsEditor::from_vec(&mut rec);
            ed.update_array_f32(b"fa", &[1.0, 2.0, 3.0]);
        }
        let arr = RawRecordView::new(&rec).tags().find_array(b"fa").expect("present");
        assert_eq!(arr.elem_type, b'f');
        assert_eq!(arr.count, 3);
    }

    #[test]
    fn test_editor_update_array_u8_same_length() {
        use crate::fields::RawRecordView;
        let mut rec = make_bam_bytes(0, 0, 0, b"r", &[], 0, -1, -1, &[]);
        {
            let mut ed = RawTagsEditor::from_vec(&mut rec);
            ed.update_array_u8(SamTag::ML, &[10u8, 20, 30]); // not present -> append
            ed.update_array_u8(SamTag::ML, &[40u8, 50, 60]); // same length -> in-place
        }
        let arr = RawRecordView::new(&rec).tags().find_array(SamTag::ML).expect("present");
        assert_eq!(arr.elem_type, b'C');
        assert_eq!(arr.count, 3);
        assert_eq!(arr.data, &[40u8, 50, 60]);
    }

    #[test]
    fn test_editor_update_array_i16_same_length() {
        use crate::fields::RawRecordView;
        let mut rec = make_bam_bytes(0, 0, 0, b"r", &[], 0, -1, -1, &[]);
        {
            let mut ed = RawTagsEditor::from_vec(&mut rec);
            ed.update_array_i16(b"sq", &[-100i16, 0, 100]); // not present -> append
            ed.update_array_i16(b"sq", &[-200i16, 0, 200]); // same length -> in-place
        }
        let arr = RawRecordView::new(&rec).tags().find_array(b"sq").expect("present");
        assert_eq!(arr.elem_type, b's');
        assert_eq!(arr.count, 3);
    }

    #[test]
    fn test_editor_update_array_i8_append_then_in_place_then_grow() {
        use crate::fields::RawRecordView;
        let mut rec = make_bam_bytes(0, 0, 0, b"r", &[], 0, -1, -1, &[]);
        {
            let mut ed = RawTagsEditor::from_vec(&mut rec);
            ed.update_array_i8(b"sb", &[-1i8, 0, 1]); // append (not present)
            ed.update_array_i8(b"sb", &[-42i8, 7, -7]); // in-place (same length)
            ed.update_array_i8(b"sb", &[10i8, 20, 30, 40, 50]); // grow
        }
        let arr = RawRecordView::new(&rec).tags().find_array(b"sb").expect("present");
        assert_eq!(arr.elem_type, b'c');
        assert_eq!(arr.count, 5);
        assert_eq!(arr.elem_size, 1);
        let decoded: Vec<i8> = arr.data.iter().map(|&b| b.cast_signed()).collect();
        assert_eq!(decoded, vec![10i8, 20, 30, 40, 50]);
    }

    #[test]
    fn test_editor_update_array_u32_append_then_in_place_then_grow() {
        use crate::fields::RawRecordView;
        let mut rec = make_bam_bytes(0, 0, 0, b"r", &[], 0, -1, -1, &[]);
        {
            let mut ed = RawTagsEditor::from_vec(&mut rec);
            ed.update_array_u32(b"uI", &[1u32, 2, 3]); // append (not present)
            ed.update_array_u32(b"uI", &[100u32, 200, 300]); // in-place (same length)
            ed.update_array_u32(b"uI", &[u32::MAX, 0, u32::MAX / 2, 42]); // grow
        }
        let arr = RawRecordView::new(&rec).tags().find_array(b"uI").expect("present");
        assert_eq!(arr.elem_type, b'I');
        assert_eq!(arr.count, 4);
        assert_eq!(arr.elem_size, 4);
        let decoded: Vec<u32> = arr
            .data
            .chunks_exact(4)
            .map(|c| u32::from_le_bytes([c[0], c[1], c[2], c[3]]))
            .collect();
        assert_eq!(decoded, vec![u32::MAX, 0, u32::MAX / 2, 42]);
    }

    // ========================================================================
    // TagValue / RawTagsView::get / iter_typed tests
    // ========================================================================

    #[test]
    fn test_tag_value_get_char_a() {
        // Build aux with XA:A:! (single ASCII char 0x21)
        let aux = [b'X', b'A', b'A', b'!'];
        let view = RawTagsView::new(&aux);
        assert_eq!(view.get(b"XA"), Some(TagValue::Char(b'!')));
    }

    #[rstest]
    #[case::signed_byte(b'c', vec![42u8], 42i64)]
    #[case::unsigned_byte(b'C', vec![200u8], 200i64)]
    #[case::signed_short(b's', (-300i16).to_le_bytes().to_vec(), -300i64)]
    #[case::unsigned_short(b'S', 50_000u16.to_le_bytes().to_vec(), 50_000i64)]
    #[case::signed_int(b'i', (-100_000i32).to_le_bytes().to_vec(), -100_000i64)]
    #[case::unsigned_int(b'I', 3_000_000_000u32.to_le_bytes().to_vec(), 3_000_000_000i64)]
    fn test_tag_value_get_int_variants(
        #[case] type_byte: u8,
        #[case] value_bytes: Vec<u8>,
        #[case] expected: i64,
    ) {
        let mut aux = vec![b'X', b'Y', type_byte];
        aux.extend_from_slice(&value_bytes);
        let view = RawTagsView::new(&aux);
        assert_eq!(view.get(b"XY"), Some(TagValue::Int(expected)));
    }

    #[test]
    fn test_tag_value_get_float() {
        let value: f32 = 1.5;
        let mut aux = vec![b'A', b'S', b'f'];
        aux.extend_from_slice(&value.to_le_bytes());
        let view = RawTagsView::new(&aux);
        match view.get(SamTag::AS) {
            Some(TagValue::Float(got)) => assert_eq!(got.to_bits(), value.to_bits()),
            other => panic!("expected TagValue::Float, got {other:?}"),
        }
    }

    #[test]
    fn test_tag_value_get_string_z() {
        // RX:Z:hello\0 — get should return TagValue::String without NUL
        let aux = b"RXZhello\x00";
        let view = RawTagsView::new(aux.as_ref());
        assert_eq!(view.get(SamTag::RX), Some(TagValue::String(b"hello")));
    }

    #[test]
    fn test_tag_value_get_hex_h() {
        // XH:H:1A2B\0 — H-type is distinct from Z-type
        let aux = b"XHH1A2B\x00";
        let view = RawTagsView::new(aux.as_ref());
        assert_eq!(view.get(b"XH"), Some(TagValue::Hex(b"1A2B")));
        // Confirm it does NOT return TagValue::String for H tags
        assert_ne!(view.get(b"XH"), Some(TagValue::String(b"1A2B")));
    }

    #[test]
    fn test_tag_value_get_array_b() {
        // Build B:i array [10, 20, 30]
        let aux = make_b_int_array_tag(*b"pa", &[10, 20, 30]);
        let view = RawTagsView::new(&aux);
        match view.get(b"pa") {
            Some(TagValue::Array(arr)) => {
                assert_eq!(arr.elem_type, b'i');
                assert_eq!(arr.count, 3);
            }
            other => panic!("expected TagValue::Array, got {other:?}"),
        }
    }

    #[test]
    fn test_tag_value_get_truncated_a_returns_none() {
        // Aux bytes "XA" + 'A' type byte but NO value byte. `get` must return
        // None rather than panic on out-of-bounds indexing.
        let aux = vec![b'X', b'A', b'A'];
        let view = RawTagsView::new(&aux);
        assert_eq!(view.get(b"XA"), None);
    }

    #[test]
    fn test_tag_value_type_byte() {
        // Scalar variants return their canonical BAM aux type byte.
        assert_eq!(TagValue::Char(b'A').type_byte(), b'A');
        assert_eq!(TagValue::Int(42).type_byte(), b'i');
        assert_eq!(TagValue::Float(1.0).type_byte(), b'f');
        assert_eq!(TagValue::String(b"xy").type_byte(), b'Z');
        assert_eq!(TagValue::Hex(b"AF").type_byte(), b'H');

        // Array variants must return b'B' (the BAM aux type byte), not the
        // element subtype — otherwise `TagValue::Array(B:i)` would look like
        // a scalar `i` tag to callers.
        let aux = make_b_int_array_tag(*b"pa", &[1, 2, 3]);
        let view = RawTagsView::new(&aux);
        let Some(TagValue::Array(arr)) = view.get(b"pa") else {
            panic!("expected Array");
        };
        assert_eq!(arr.elem_type, b'i', "element subtype preserved on ArrayTagRef");
        assert_eq!(TagValue::Array(arr).type_byte(), b'B', "aux type byte for array is B");
    }

    #[test]
    fn test_tag_value_iter_typed_yields_duplicate_keys_distinctly() {
        // Two entries with the same tag but different values. The old
        // implementation rescanned via `get()` and would yield the first
        // value twice; the direct-decode path must yield each on its own.
        let mut aux = Vec::new();
        aux.extend_from_slice(&[b'X', b'N', b'C', 7u8]); // XN:C:7
        aux.extend_from_slice(&[b'X', b'N', b'C', 42u8]); // XN:C:42 (duplicate key)

        let view = RawTagsView::new(&aux);
        let pairs: Vec<([u8; 2], TagValue<'_>)> = view.iter_typed().collect();
        assert_eq!(pairs.len(), 2);
        assert_eq!(pairs[0], (*b"XN", TagValue::Int(7)));
        assert_eq!(pairs[1], (*b"XN", TagValue::Int(42)));
    }

    #[test]
    fn test_tag_value_iter_typed_roundtrip() {
        // Record with 3 tags of different types: NM:C:5, RG:Z:lib1, AS:f:1.5
        let mut aux = Vec::new();
        aux.extend_from_slice(&[b'N', b'M', b'C', 5u8]); // NM:C:5
        aux.extend_from_slice(b"RGZlib1\x00"); // RG:Z:lib1
        aux.extend_from_slice(b"ASf"); // AS:f:1.5
        aux.extend_from_slice(&1.5f32.to_le_bytes());

        let view = RawTagsView::new(&aux);
        let pairs: Vec<([u8; 2], TagValue<'_>)> = view.iter_typed().collect();

        assert_eq!(pairs.len(), 3);
        assert_eq!(pairs[0], (<[u8; 2]>::from(SamTag::NM), TagValue::Int(5)));
        assert_eq!(pairs[1], (<[u8; 2]>::from(SamTag::RG), TagValue::String(b"lib1")));
        let (tag2, val2) = pairs[2];
        assert_eq!(tag2, <[u8; 2]>::from(SamTag::AS));
        match val2 {
            TagValue::Float(f) => assert_eq!(f.to_bits(), 1.5f32.to_bits()),
            other => panic!("expected TagValue::Float(1.5), got {other:?}"),
        }
    }

    // ========================================================================
    // AsTagBytes smoke test — SamTag and &[u8; 2] both accepted
    // ========================================================================

    #[test]
    fn test_find_string_tag_accepts_sam_tag() {
        // Build a tiny aux block: RX:Z:ACGT\0
        let mut aux = Vec::new();
        aux.extend_from_slice(b"RXZACGT\0");
        // Accept SamTag.
        let v = find_string_tag(&aux, SamTag::RX);
        assert_eq!(v, Some(b"ACGT".as_ref()));
        // And still accept &[u8; 2] for backward compat.
        let v2 = find_string_tag(&aux, b"RX");
        assert_eq!(v2, Some(b"ACGT".as_ref()));
    }

    #[test]
    fn append_raw_tag_writes_tag_type_and_value_verbatim() {
        // Start from an aux block holding RX:Z:ACGT and copy it verbatim under CB.
        let mut record = b"RXZACGT\0".to_vec();
        let entry = RawTagsView::new(&record).iter().next().expect("one source tag");
        let (type_byte, value_bytes) = (entry.type_byte, entry.value_bytes.to_vec());

        append_raw_tag(&mut record, SamTag::CB, type_byte, &value_bytes);

        // The appended entry is byte-identical to the source apart from the tag id,
        // so the destination value (NUL terminator included) round-trips.
        assert_eq!(find_string_tag(&record, SamTag::CB), Some(b"ACGT".as_ref()));
        assert_eq!(&record[8..], b"CBZACGT\0", "verbatim tag+type+value bytes");
    }

    #[test]
    fn test_find_two_string_tags_both_present_either_order() {
        // MI before CB.
        let aux = b"MIZ7\x00CBZACGT\x00";
        let rec = make_bam_bytes(0, 0, 0, b"rea", &[], 0, -1, -1, aux);
        let (mi, cb) = find_two_string_tags_in_record(&rec, SamTag::MI, SamTag::CB);
        assert_eq!(mi, Some(b"7".as_ref()));
        assert_eq!(cb, Some(b"ACGT".as_ref()));

        // CB before MI: result is independent of aux ordering and of arg order.
        let aux = b"CBZACGT\x00MIZ7\x00";
        let rec = make_bam_bytes(0, 0, 0, b"rea", &[], 0, -1, -1, aux);
        let (mi, cb) = find_two_string_tags_in_record(&rec, SamTag::MI, SamTag::CB);
        assert_eq!(mi, Some(b"7".as_ref()));
        assert_eq!(cb, Some(b"ACGT".as_ref()));
    }

    #[test]
    fn test_find_two_string_tags_agrees_with_single_lookups() {
        // Interleave a non-Z tag to make the walk non-trivial.
        let mut aux = Vec::new();
        aux.extend_from_slice(b"NMC");
        aux.push(3);
        aux.extend_from_slice(b"MIZ42\x00");
        aux.extend_from_slice(b"CBZTGCA\x00");
        let rec = make_bam_bytes(0, 0, 0, b"rea", &[], 0, -1, -1, &aux);
        let (mi, cb) = find_two_string_tags_in_record(&rec, SamTag::MI, SamTag::CB);
        assert_eq!(mi, find_string_tag_in_record(&rec, SamTag::MI));
        assert_eq!(cb, find_string_tag_in_record(&rec, SamTag::CB));
    }

    #[test]
    fn test_find_two_string_tags_first_match_wins_on_malformed_dup() {
        // A non-Z MI entry precedes a Z MI duplicate. `find_string_tag_in_record`
        // resolves the FIRST MI match (non-Z → None) and never reaches the later
        // Z duplicate; the two-tag walk must match that rather than skipping the
        // non-Z entry and returning the later Z value.
        let mut aux = Vec::new();
        aux.extend_from_slice(b"MIC"); // MI, type 'C' (uint8) ...
        aux.push(9); //                ... 1-byte value (non-Z → MI resolves None)
        aux.extend_from_slice(b"MIZ77\x00"); // later Z duplicate — must be ignored
        aux.extend_from_slice(b"CBZACGT\x00");
        let rec = make_bam_bytes(0, 0, 0, b"rea", &[], 0, -1, -1, &aux);

        let (mi, cb) = find_two_string_tags_in_record(&rec, SamTag::MI, SamTag::CB);
        assert_eq!(mi, None, "first (non-Z) MI match resolves the tag to None");
        assert_eq!(cb, Some(b"ACGT".as_ref()));
        // Must agree with two independent single-tag lookups.
        assert_eq!(mi, find_string_tag_in_record(&rec, SamTag::MI));
        assert_eq!(cb, find_string_tag_in_record(&rec, SamTag::CB));
    }

    /// A `Z` entry with no terminating NUL, in the two shapes that behave
    /// differently.
    ///
    /// The doc contract says a first match that is non-`Z` *or unterminated*
    /// resolves the tag to `None`. Only the non-`Z` half was pinned, and the
    /// unterminated half is subtler than it reads: the value scan looks for the
    /// next NUL anywhere in the remaining aux block, not just within the entry.
    /// So "unterminated" resolves to `None` only when there is no NUL left at
    /// all — if a later entry supplies one, the value runs *through* that entry
    /// and swallows it.
    ///
    /// Both shapes are pinned because they fail differently: the first stops the
    /// walk via `tag_value_size` returning `None`, the second consumes the
    /// remaining bytes as one oversized value. Either way the two-tag walk must
    /// agree with two independent single-tag lookups, which is the contract that
    /// makes `find_two_string_tags_in_record` a safe substitution for calling
    /// `find_string_tag_in_record` twice.
    #[rstest]
    // No NUL anywhere after the tag: the value scan finds nothing, so MI is
    // `None`, and the walk stops rather than advancing past a sizeless entry.
    #[case::no_nul_anywhere(b"MIZ77".as_ref(), None, None)]
    // A later entry supplies the NUL, so MI's value runs through `CBZACGT` and
    // consumes it — CB is never seen as an entry of its own.
    #[case::nul_supplied_by_a_later_entry(
        b"MIZ77CBZACGT\x00".as_ref(),
        Some(b"77CBZACGT".as_ref()),
        None
    )]
    fn test_find_two_string_tags_unterminated_z_first_match(
        #[case] aux: &[u8],
        #[case] expected_mi: Option<&[u8]>,
        #[case] expected_cb: Option<&[u8]>,
    ) {
        let rec = make_bam_bytes(0, 0, 0, b"rea", &[], 0, -1, -1, aux);

        let (mi, cb) = find_two_string_tags_in_record(&rec, SamTag::MI, SamTag::CB);
        assert_eq!(mi, expected_mi);
        assert_eq!(cb, expected_cb);
        // The parity contract: identical to two independent single-tag lookups.
        assert_eq!(
            mi,
            find_string_tag_in_record(&rec, SamTag::MI),
            "MI must match a single lookup"
        );
        assert_eq!(
            cb,
            find_string_tag_in_record(&rec, SamTag::CB),
            "CB must match a single lookup"
        );
    }

    #[test]
    fn test_find_two_string_tags_same_tag_fills_both() {
        // When the same tag is requested for both first and second, a single
        // matching aux entry must populate both outputs (not just one).
        let aux = b"MIZ7\x00CBZACGT\x00";
        let rec = make_bam_bytes(0, 0, 0, b"rea", &[], 0, -1, -1, aux);
        let (a, b) = find_two_string_tags_in_record(&rec, SamTag::MI, SamTag::MI);
        assert_eq!(a, Some(b"7".as_ref()));
        assert_eq!(b, Some(b"7".as_ref()));
        // Both outputs must agree with an independent single-tag lookup, proving
        // the same-tag case behaves like two separate `find_string_tag_in_record`
        // calls (the contract the `mi_group.rs` production path relies on).
        let mi = find_string_tag_in_record(&rec, SamTag::MI);
        assert_eq!(a, mi, "first output must match an independent MI lookup");
        assert_eq!(b, mi, "second output must match an independent MI lookup");

        // Same-tag MISS path: when the (identical) requested tag is absent, both
        // outputs must be None — the zero-match branch of the same-tag case,
        // matching two independent single-tag lookups that each miss.
        let rec_without_mi = make_bam_bytes(0, 0, 0, b"rea", &[], 0, -1, -1, b"CBZACGT\x00");
        let (a_missing, b_missing) =
            find_two_string_tags_in_record(&rec_without_mi, SamTag::MI, SamTag::MI);
        assert_eq!(a_missing, None);
        assert_eq!(b_missing, None);
        assert_eq!(a_missing, find_string_tag_in_record(&rec_without_mi, SamTag::MI));
    }

    // ========================================================================
    // RawTagsEditor::rebuild_with — single-pass aux rebuild
    // ========================================================================

    /// Build a record whose aux block is the concatenation of `entries`
    /// (each already in on-disk `[tag0, tag1, type, value...]` form).
    fn rec_with_aux(entries: &[&[u8]]) -> Vec<u8> {
        let aux: Vec<u8> = entries.iter().flat_map(|e| e.iter().copied()).collect();
        make_bam_bytes(0, 0, 0, b"r", &[], 0, -1, -1, &aux)
    }

    /// The naive oracle: `remove_tag` for every key in `remove`, then
    /// `remove_tag(add.tag); append_raw_tag(add)` for every add — the exact
    /// per-tag loop `rebuild_with` replaces. Returns the full record bytes.
    fn oracle_rebuild(base: &[u8], remove: &[[u8; 2]], adds: &[TagEntry<'_>]) -> Vec<u8> {
        let mut rec = base.to_vec();
        for &k in remove {
            remove_tag(&mut rec, k);
        }
        for a in adds {
            remove_tag(&mut rec, a.tag);
            append_raw_tag(&mut rec, a.tag, a.type_byte, a.value_bytes);
        }
        rec
    }

    /// Collect the aux entries (tag, type, value bytes) of a record, in order.
    fn aux_entries(rec: &[u8]) -> Vec<([u8; 2], u8, Vec<u8>)> {
        RawTagsView::new(aux_data_slice(rec))
            .iter()
            .map(|e| (e.tag, e.type_byte, e.value_bytes.to_vec()))
            .collect()
    }

    /// `SamTag` -> raw two-byte key, so tests use the tag constants (not bare
    /// byte literals) for tag identifiers. `BX` has no `SamTag` constant and is
    /// path-scoped-allowlisted in the tag-literal check as an opaque fixture.
    fn key(t: SamTag) -> [u8; 2] {
        t.into()
    }
    const BX: [u8; 2] = *b"BX";

    #[test]
    fn rebuild_with_drops_removed_tags_and_keeps_order() {
        // Aux: RG(Z), NM(C=1B), MD(Z). Remove NM; keep RG, MD in order.
        let mut rec = rec_with_aux(&[b"RGZg\x00", b"NMC\x05", b"MDZ10\x00"]);
        {
            let mut ed = RawTagsEditor::from_vec(&mut rec);
            ed.rebuild_with(&[key(SamTag::NM)], &[]);
        }
        let got: Vec<[u8; 2]> = aux_entries(&rec).into_iter().map(|(t, _, _)| t).collect();
        assert_eq!(got, vec![key(SamTag::RG), key(SamTag::MD)], "NM dropped; RG,MD kept in order");
    }

    #[test]
    fn rebuild_with_appends_adds_verbatim_after_survivors() {
        let mut rec = rec_with_aux(&[b"RGZg\x00"]);
        let adds = [
            TagEntry { tag: BX, type_byte: b'Z', value_bytes: b"ACGT\x00" },
            TagEntry { tag: key(SamTag::NM), type_byte: b'C', value_bytes: &[7] },
        ];
        {
            let mut ed = RawTagsEditor::from_vec(&mut rec);
            ed.rebuild_with(&[] as &[[u8; 2]], &adds);
        }
        let got = aux_entries(&rec);
        assert_eq!(got[0].0, key(SamTag::RG), "survivor first");
        assert_eq!((got[1].0, got[1].1, got[1].2.as_slice()), (BX, b'Z', b"ACGT\x00".as_slice()));
        assert_eq!(
            (got[2].0, got[2].1, got[2].2.as_slice()),
            (key(SamTag::NM), b'C', [7].as_slice())
        );
    }

    #[test]
    fn rebuild_with_upserts_an_added_key_that_already_exists() {
        // NM exists as C=1; adding NM as i=999 must replace it (no duplicate),
        // and the new value must appear appended (after other survivors).
        let mut rec = rec_with_aux(&[b"NMC\x01", b"RGZg\x00"]);
        let nm_i = {
            let mut v = vec![b'N', b'M', b'i'];
            v.extend_from_slice(&999i32.to_le_bytes());
            v
        };
        let adds = [TagEntry { tag: key(SamTag::NM), type_byte: b'i', value_bytes: &nm_i[3..] }];
        {
            let mut ed = RawTagsEditor::from_vec(&mut rec);
            ed.rebuild_with(&[] as &[[u8; 2]], &adds);
        }
        let entries = aux_entries(&rec);
        let nm: Vec<_> = entries.iter().filter(|(t, _, _)| *t == key(SamTag::NM)).collect();
        assert_eq!(nm.len(), 1, "exactly one NM after upsert");
        assert_eq!(nm[0].1, b'i', "new NM type wins");
        let v = RawTagsView::new(aux_data_slice(&rec)).find_int(SamTag::NM);
        assert_eq!(v, Some(999));
        // Survivor RG still present.
        assert!(entries.iter().any(|(t, _, _)| *t == key(SamTag::RG)));
    }

    #[test]
    fn rebuild_with_on_empty_aux_appends_only() {
        let mut rec = make_bam_bytes(0, 0, 0, b"r", &[], 0, -1, -1, &[]);
        let adds = [TagEntry { tag: key(SamTag::RG), type_byte: b'Z', value_bytes: b"g\x00" }];
        {
            let mut ed = RawTagsEditor::from_vec(&mut rec);
            ed.rebuild_with(&[key(SamTag::NM)], &adds); // remove of an absent tag is a no-op
        }
        assert_eq!(
            RawTagsView::new(aux_data_slice(&rec)).find_string(SamTag::RG),
            Some(b"g".as_slice())
        );
    }

    #[test]
    fn rebuild_with_matches_the_naive_remove_append_oracle() {
        // A representative mixed aux block: Z, C, Z, f.
        let mut f_tag = vec![b'A', b'S', b'f'];
        f_tag.extend_from_slice(&3.5f32.to_le_bytes());
        let base = rec_with_aux(&[b"RGZgrp\x00", b"NMC\x03", b"MDZ50\x00", &f_tag]);

        let bx = TagEntry { tag: BX, type_byte: b'Z', value_bytes: b"ACGTAC\x00" };
        // RG is both removed and re-added — upsert should land it once, appended.
        let rg = TagEntry { tag: key(SamTag::RG), type_byte: b'Z', value_bytes: b"newgrp\x00" };
        let adds = [bx, rg];
        let remove = [key(SamTag::NM)];

        let expected = oracle_rebuild(&base, &remove, &adds);

        let mut got = base.clone();
        {
            let mut ed = RawTagsEditor::from_vec(&mut got);
            ed.rebuild_with(&remove, &adds);
        }
        assert_eq!(
            aux_entries(&got),
            aux_entries(&expected),
            "rebuild_with must produce the same aux entries as the naive remove+append loop"
        );
    }

    #[test]
    fn rebuild_with_readds_a_key_that_is_also_in_the_remove_set() {
        // AS is both in `remove` and in `adds` (the normalize idiom: drop the
        // old encoding, append the canonical one). The result must carry exactly
        // one AS, with the added value — same as remove_tag(AS) then append(AS).
        let mut rec = rec_with_aux(&[b"ASC\x05", b"RGZg\x00"]);
        let as_add = TagEntry { tag: key(SamTag::AS), type_byte: b'c', value_bytes: &[5] };
        {
            let mut ed = RawTagsEditor::from_vec(&mut rec);
            ed.rebuild_with(&[key(SamTag::AS)], &[as_add]);
        }
        let entries = aux_entries(&rec);
        let as_count = entries.iter().filter(|(t, _, _)| *t == key(SamTag::AS)).count();
        assert_eq!(as_count, 1, "exactly one AS after remove+readd");
        assert_eq!(
            entries.iter().find(|(t, _, _)| *t == key(SamTag::AS)).unwrap().1,
            b'c',
            "readded type"
        );
        assert!(entries.iter().any(|(t, _, _)| *t == key(SamTag::RG)), "RG survivor kept");
    }

    #[test]
    fn rebuild_with_tagbitset_removes_the_configured_set() {
        let mut set = TagBitset::new();
        set.insert(key(SamTag::NM));
        set.insert(key(SamTag::MD));
        let mut rec = rec_with_aux(&[b"RGZg\x00", b"NMC\x05", b"MDZ10\x00"]);
        {
            let mut ed = RawTagsEditor::from_vec(&mut rec);
            ed.rebuild_with(&set, &[]);
        }
        let got: Vec<[u8; 2]> = aux_entries(&rec).into_iter().map(|(t, _, _)| t).collect();
        assert_eq!(got, vec![key(SamTag::RG)], "TagBitset set of NM and MD dropped");
    }

    #[test]
    fn rebuild_with_dedups_duplicate_keys_within_adds_last_wins() {
        // `adds` carries NM twice (C:1 then i:999). The naive idiom does
        // `remove_tag(NM); append(NM:C:1); remove_tag(NM); append(NM:i:999)`,
        // so only the last NM survives — rebuild_with must match that.
        let mut rec = rec_with_aux(&[b"RGZg\x00"]);
        let nm_i = 999i32.to_le_bytes();
        let adds = [
            TagEntry { tag: key(SamTag::NM), type_byte: b'C', value_bytes: &[1] },
            TagEntry { tag: key(SamTag::NM), type_byte: b'i', value_bytes: &nm_i },
        ];
        {
            let mut ed = RawTagsEditor::from_vec(&mut rec);
            ed.rebuild_with(&[] as &[[u8; 2]], &adds);
        }

        let entries = aux_entries(&rec);
        let nm: Vec<_> = entries.iter().filter(|(t, _, _)| *t == key(SamTag::NM)).collect();
        assert_eq!(nm.len(), 1, "exactly one NM after dedup within adds");
        assert_eq!(nm[0].1, b'i', "the last NM (i:999) wins");
        assert_eq!(RawTagsView::new(aux_data_slice(&rec)).find_int(SamTag::NM), Some(999));
        assert!(entries.iter().any(|(t, _, _)| *t == key(SamTag::RG)), "RG survivor kept");

        // Byte-identical to the naive remove+append oracle for the same inputs.
        let base = rec_with_aux(&[b"RGZg\x00"]);
        let expected = oracle_rebuild(&base, &[], &adds);
        assert_eq!(aux_entries(&rec), aux_entries(&expected));
    }

    #[test]
    fn rebuild_with_clamps_an_aux_offset_past_the_end_of_a_truncated_record() {
        // A 20-byte record whose header declares l_seq = 1000, so the computed
        // aux offset (32 + l_read_name + n_cigar*4 + ceil(l_seq/2) + l_seq) lands
        // far past the 20-byte buffer. `from_vec` caches that out-of-range offset;
        // `rebuild_with` must clamp it to `record.len()` (the empty-aux path) and
        // append `adds` without panicking or reading out of bounds.
        let mut rec = vec![0u8; 20];
        rec[8] = 1; // l_read_name (n_cigar at 12..14 = 0)
        rec[16..20].copy_from_slice(&1000u32.to_le_bytes()); // l_seq
        assert!(
            aux_data_offset_from_record(&rec).unwrap() > rec.len(),
            "test precondition: header must declare an out-of-range aux offset"
        );

        let adds = [TagEntry { tag: key(SamTag::RG), type_byte: b'Z', value_bytes: b"g\x00" }];
        {
            let mut ed = RawTagsEditor::from_vec(&mut rec);
            ed.rebuild_with(&[] as &[[u8; 2]], &adds); // must not panic
        }
        // The appended entry lands verbatim at the clamped offset (== old len).
        assert_eq!(rec.len(), 20 + 5, "RG:Z:g\\0 is 5 bytes, appended after the 20-byte header");
        assert_eq!(&rec[20..], b"RGZg\x00");
    }

    #[test]
    fn rebuild_with_stops_at_a_malformed_aux_entry() {
        // Aux: RG(Z), then a truncated tag ("XN" with type 'i' but no value
        // bytes) — the walk must stop at the malformed entry and drop it, per
        // the documented AuxTagsIter tolerance, leaving only the RG survivor.
        let mut rec = rec_with_aux(&[b"RGZg\x00", b"XNi"]);
        {
            let mut ed = RawTagsEditor::from_vec(&mut rec);
            ed.rebuild_with(&[] as &[[u8; 2]], &[]);
        }
        let got: Vec<[u8; 2]> = aux_entries(&rec).into_iter().map(|(t, _, _)| t).collect();
        assert_eq!(got, vec![key(SamTag::RG)], "walk stops at malformed entry; RG kept");
    }

    // ========================================================================
    // RawTagsEditor::rebuild_with_int_normalized — fold AS/XS int normalization
    // into the single rebuild walk. Its contract is defined by an oracle: it
    // must produce byte-identical output to `rebuild_with(remove, adds)` followed
    // by `normalize_int_tag_to_smallest_signed(t)` for each `t` in `normalize`.
    // ========================================================================

    /// Build raw entry bytes `[tag0, tag1, type, value...]`.
    fn ent(tag: [u8; 2], ty: u8, val: &[u8]) -> Vec<u8> {
        let mut v = vec![tag[0], tag[1], ty];
        v.extend_from_slice(val);
        v
    }

    fn as_xs() -> [[u8; 2]; 2] {
        [key(SamTag::AS), key(SamTag::XS)]
    }

    /// Today's behavior, the oracle: `rebuild_with`, then normalize each tag in order.
    fn oracle_rebuild_then_normalize(
        base: &[u8],
        remove: &[[u8; 2]],
        adds: &[TagEntry<'_>],
        normalize: &[[u8; 2]],
    ) -> Vec<u8> {
        let mut rec = base.to_vec();
        {
            let mut ed = RawTagsEditor::from_vec(&mut rec);
            ed.rebuild_with(remove, adds);
        }
        for &t in normalize {
            normalize_int_tag_to_smallest_signed(&mut rec, t);
        }
        rec
    }

    fn folded_rebuild_normalize(
        base: &[u8],
        remove: &[[u8; 2]],
        adds: &[TagEntry<'_>],
        normalize: &[[u8; 2]],
    ) -> Vec<u8> {
        let mut rec = base.to_vec();
        {
            let mut scratch = Vec::new();
            let mut ed = RawTagsEditor::from_vec(&mut rec);
            ed.rebuild_with_int_normalized(remove, adds, normalize, &mut scratch);
        }
        rec
    }

    fn assert_folded_matches_oracle(
        base: &[u8],
        remove: &[[u8; 2]],
        adds: &[TagEntry<'_>],
        normalize: &[[u8; 2]],
        label: &str,
    ) {
        let expected = oracle_rebuild_then_normalize(base, remove, adds, normalize);
        let got = folded_rebuild_normalize(base, remove, adds, normalize);
        assert_eq!(
            aux_entries(&got),
            aux_entries(&expected),
            "{label}: folded aux entries must equal the rebuild+normalize oracle"
        );
        assert_eq!(
            got, expected,
            "{label}: folded record bytes must be byte-identical to the oracle"
        );
    }

    #[test]
    fn folded_normalize_common_i32_as_xs() {
        // bwa-mem3 shape: AS:i:139 and XS:i:0 (int32), interleaved with other
        // tags, plus a copied-in BX. Both must end normalized (AS->i16, XS->i8)
        // and relocated to the tail, after the BX add.
        let as_i = ent(key(SamTag::AS), b'i', &139i32.to_le_bytes());
        let xs_i = ent(key(SamTag::XS), b'i', &0i32.to_le_bytes());
        let base = rec_with_aux(&[b"RGZg\x00", &as_i, b"NMC\x02", &xs_i]);
        let bx = TagEntry { tag: BX, type_byte: b'Z', value_bytes: b"ACGT\x00" };
        assert_folded_matches_oracle(&base, &[] as &[[u8; 2]], &[bx], &as_xs(), "common i32 AS/XS");
    }

    #[test]
    fn folded_normalize_only_as_present() {
        let as_i = ent(key(SamTag::AS), b'i', &200i32.to_le_bytes());
        let base = rec_with_aux(&[b"RGZg\x00", &as_i]);
        assert_folded_matches_oracle(&base, &[] as &[[u8; 2]], &[], &as_xs(), "only AS present");
    }

    #[test]
    fn folded_normalize_both_absent_is_noop() {
        let base = rec_with_aux(&[b"RGZg\x00", b"NMC\x02"]);
        assert_folded_matches_oracle(&base, &[] as &[[u8; 2]], &[], &as_xs(), "AS/XS absent");
    }

    #[test]
    fn folded_normalize_as_already_smallest_width() {
        // AS already 'c' (fits i8). The oracle still remove+re-appends it (moving
        // it to the tail as 'c'); the folded pass must reproduce that relocation.
        let as_c = ent(key(SamTag::AS), b'c', &[5]);
        let base = rec_with_aux(&[&as_c, b"RGZg\x00"]);
        assert_folded_matches_oracle(&base, &[] as &[[u8; 2]], &[], &as_xs(), "AS already i8");
    }

    #[test]
    fn folded_normalize_non_int_as_left_in_place() {
        // AS is a Z-string, not an integer: normalize is a no-op, so AS must stay
        // verbatim at its original position (NOT relocated).
        let as_z = ent(key(SamTag::AS), b'Z', b"hello\x00");
        let base = rec_with_aux(&[&as_z, b"RGZg\x00"]);
        assert_folded_matches_oracle(
            &base,
            &[] as &[[u8; 2]],
            &[],
            &as_xs(),
            "non-int AS left in place",
        );
    }

    #[test]
    fn folded_normalize_value_exceeds_i32_left_in_place() {
        // AS:I (u32) with a value above i32::MAX: normalize early-returns without
        // touching it, so AS stays verbatim in place.
        let as_big = ent(key(SamTag::AS), b'I', &3_000_000_000u32.to_le_bytes());
        let base = rec_with_aux(&[&as_big, b"RGZg\x00"]);
        assert_folded_matches_oracle(
            &base,
            &[] as &[[u8; 2]],
            &[],
            &as_xs(),
            "AS value exceeds i32",
        );
    }

    #[test]
    fn folded_normalize_as_supplied_via_adds_upsert() {
        // A survivor AS:i:139 is upserted by an adds AS:i:200; normalization must
        // apply to the winning (post-upsert) value.
        let as_i = ent(key(SamTag::AS), b'i', &139i32.to_le_bytes());
        let base = rec_with_aux(&[&as_i, b"RGZg\x00"]);
        let as_add =
            TagEntry { tag: key(SamTag::AS), type_byte: b'i', value_bytes: &200i32.to_le_bytes() };
        assert_folded_matches_oracle(
            &base,
            &[] as &[[u8; 2]],
            &[as_add],
            &as_xs(),
            "AS supplied via adds",
        );
    }

    #[test]
    fn folded_normalize_with_remove_and_copy_set() {
        // Full merge_raw_with shape: a remove-set, a copied BX, AS/XS to normalize.
        let as_i = ent(key(SamTag::AS), b'i', &1000i32.to_le_bytes());
        let xs_i = ent(key(SamTag::XS), b'i', &50i32.to_le_bytes());
        let base = rec_with_aux(&[b"RGZg\x00", b"NMC\x03", &as_i, b"MDZ50\x00", &xs_i]);
        let bx = TagEntry { tag: BX, type_byte: b'Z', value_bytes: b"ACGTAC\x00" };
        assert_folded_matches_oracle(
            &base,
            &[key(SamTag::NM)],
            &[bx],
            &as_xs(),
            "remove + copy + normalize",
        );
    }

    #[test]
    fn folded_normalize_empty_aux_with_adds() {
        let base = make_bam_bytes(0, 0, 0, b"r", &[], 0, -1, -1, &[]);
        let bx = TagEntry { tag: BX, type_byte: b'Z', value_bytes: b"AC\x00" };
        assert_folded_matches_oracle(&base, &[] as &[[u8; 2]], &[bx], &as_xs(), "empty aux + adds");
    }

    /// A random integer aux value in any of the six BAM int encodings — including
    /// `i`/`I` values above `i32::MAX`, which exercise the "left in place" path.
    fn arb_int_bytes() -> impl proptest::strategy::Strategy<Value = (u8, Vec<u8>)> {
        use proptest::prelude::*;
        prop_oneof![
            any::<i8>().prop_map(|v| (b'c', vec![v.cast_unsigned()])),
            any::<u8>().prop_map(|v| (b'C', vec![v])),
            any::<i16>().prop_map(|v| (b's', v.to_le_bytes().to_vec())),
            any::<u16>().prop_map(|v| (b'S', v.to_le_bytes().to_vec())),
            any::<i32>().prop_map(|v| (b'i', v.to_le_bytes().to_vec())),
            any::<u32>().prop_map(|v| (b'I', v.to_le_bytes().to_vec())),
        ]
    }

    proptest::proptest! {
        /// The folded pass must be byte-identical to `rebuild_with` + sequential
        /// `normalize_int_tag_to_smallest_signed(AS)`/`(XS)` for arbitrary AS/XS
        /// encodings, presence, an optional AS upsert, a copied BX, and a removed NM.
        #[test]
        fn folded_normalize_matches_oracle_prop(
            as_present in proptest::prelude::any::<bool>(),
            as_enc in arb_int_bytes(),
            xs_present in proptest::prelude::any::<bool>(),
            xs_enc in arb_int_bytes(),
            add_bx in proptest::prelude::any::<bool>(),
            add_as in proptest::option::of(proptest::prelude::any::<i32>()),
            remove_nm in proptest::prelude::any::<bool>(),
        ) {
            let mut entries: Vec<Vec<u8>> = vec![ent(key(SamTag::RG), b'Z', b"g\x00")];
            if as_present { entries.push(ent(key(SamTag::AS), as_enc.0, &as_enc.1)); }
            entries.push(ent(key(SamTag::NM), b'C', &[3]));
            entries.push(ent(key(SamTag::MD), b'Z', b"50\x00"));
            if xs_present { entries.push(ent(key(SamTag::XS), xs_enc.0, &xs_enc.1)); }
            let refs: Vec<&[u8]> = entries.iter().map(Vec::as_slice).collect();
            let base = rec_with_aux(&refs);

            let as_add_bytes = add_as.map(i32::to_le_bytes);
            let mut adds: Vec<TagEntry<'_>> = Vec::new();
            if add_bx {
                adds.push(TagEntry { tag: BX, type_byte: b'Z', value_bytes: b"AC\x00" });
            }
            if let Some(ref b) = as_add_bytes {
                adds.push(TagEntry { tag: key(SamTag::AS), type_byte: b'i', value_bytes: b });
            }
            let remove: Vec<[u8; 2]> = if remove_nm { vec![key(SamTag::NM)] } else { vec![] };

            let expected = oracle_rebuild_then_normalize(&base, &remove, &adds, &as_xs());
            let got = folded_rebuild_normalize(&base, &remove, &adds, &as_xs());
            proptest::prop_assert_eq!(got, expected);
        }
    }

    #[test]
    fn folded_normalize_empty_normalize_list_equals_rebuild_with() {
        // The empty-`normalize` fast path must be identical to plain `rebuild_with`.
        let base = rec_with_aux(&[b"RGZg\x00", b"NMC\x02"]);
        let bx = TagEntry { tag: BX, type_byte: b'Z', value_bytes: b"AC\x00" };
        let mut got = base.clone();
        {
            let mut sc = Vec::new();
            let mut ed = RawTagsEditor::from_vec(&mut got);
            ed.rebuild_with_int_normalized(&[] as &[[u8; 2]], &[bx], &[], &mut sc);
        }
        let mut expected = base.clone();
        {
            let mut ed = RawTagsEditor::from_vec(&mut expected);
            ed.rebuild_with(&[] as &[[u8; 2]], &[bx]);
        }
        assert_eq!(got, expected, "empty normalize must equal rebuild_with");
    }

    #[test]
    fn folded_normalize_spills_to_heap_beyond_inline_cap() {
        // Exercise the heap-spill branch: pass more normalize tags than the
        // inline stack array holds (> INLINE_NORMALIZE = 8). Tag keys are built
        // programmatically (`Zx`) to avoid literal-tag churn; each is an int tag
        // wide enough that normalization actually shrinks it, so the spill
        // indexing is genuinely exercised, not a no-op.
        let n: u8 = 12; // > INLINE_NORMALIZE
        let keys: Vec<[u8; 2]> = (0..n).map(|i| [b'Z', b'0' + i]).collect();
        let entries: Vec<Vec<u8>> = keys
            .iter()
            .enumerate()
            .map(|(i, k)| ent(*k, b'i', &(300 + i32::try_from(i).unwrap()).to_le_bytes()))
            .collect();
        let refs: Vec<&[u8]> = entries.iter().map(Vec::as_slice).collect();
        let base = rec_with_aux(&refs);
        // Oracle: rebuild_with(no-op) then normalize each key sequentially.
        let expected = oracle_rebuild_then_normalize(&base, &[] as &[[u8; 2]], &[], &keys);
        let got = folded_rebuild_normalize(&base, &[] as &[[u8; 2]], &[], &keys);
        assert_eq!(
            aux_entries(&got),
            aux_entries(&expected),
            "heap-spill normalize path must match the sequential oracle"
        );
        assert_eq!(got, expected, "heap-spill normalize path must be byte-identical");
    }

    #[test]
    fn folded_normalize_scales_past_former_bitmask_width() {
        // Regression: the seen-state that tracks each normalize tag's first
        // survivor occurrence was a u64 bitmask capped at 64 entries — a debug
        // panic and, in release, a silent failure to capture (leaving entries at
        // indices >= 64 unnormalized, diverging from the sequential oracle). Pass
        // more than 64 distinct int tags (each wide enough that normalization
        // shrinks it) so the tags beyond the old cap are genuinely normalized.
        let n: u16 = 70; // > former 64-bit cap and > INLINE_NORMALIZE
        let keys: Vec<[u8; 2]> = (0..n).map(|i| [b'Z', 32 + u8::try_from(i).unwrap()]).collect();
        let entries: Vec<Vec<u8>> = keys
            .iter()
            .enumerate()
            .map(|(i, k)| ent(*k, b'i', &(300 + i32::try_from(i).unwrap()).to_le_bytes()))
            .collect();
        let refs: Vec<&[u8]> = entries.iter().map(Vec::as_slice).collect();
        let base = rec_with_aux(&refs);
        let expected = oracle_rebuild_then_normalize(&base, &[] as &[[u8; 2]], &[], &keys);
        let got = folded_rebuild_normalize(&base, &[] as &[[u8; 2]], &[], &keys);
        assert_eq!(
            aux_entries(&got),
            aux_entries(&expected),
            "normalize list beyond 64 entries must match the sequential oracle"
        );
        assert_eq!(got, expected, "normalize list beyond 64 entries must be byte-identical");
    }

    #[test]
    fn folded_normalize_repeated_normalize_key_matches_sequential_oracle() {
        // A `normalize` list with a repeated key ([AS, XS, AS]) is degenerate, so
        // the folded path falls back to the literal sequential oracle. Assert the
        // fallback reproduces it for both a survivor-sourced value (AS in the aux)
        // and an adds-sourced value (AS supplied via adds).
        let normalize = [key(SamTag::AS), key(SamTag::XS), key(SamTag::AS)];

        // Survivor-sourced: AS/XS live in the aux.
        let as_i = ent(key(SamTag::AS), b'i', &139i32.to_le_bytes());
        let xs_i = ent(key(SamTag::XS), b'i', &0i32.to_le_bytes());
        let base = rec_with_aux(&[b"RGZg\x00", &as_i, &xs_i]);
        assert_folded_matches_oracle(
            &base,
            &[] as &[[u8; 2]],
            &[],
            &normalize,
            "repeated AS in normalize, survivor value",
        );

        // Adds-sourced: AS supplied via adds (upsert winner).
        let base2 = rec_with_aux(&[b"RGZg\x00", &xs_i]);
        let as_add =
            TagEntry { tag: key(SamTag::AS), type_byte: b'i', value_bytes: &139i32.to_le_bytes() };
        assert_folded_matches_oracle(
            &base2,
            &[] as &[[u8; 2]],
            &[as_add],
            &normalize,
            "repeated AS in normalize, adds value",
        );
    }

    #[test]
    fn folded_normalize_repeated_key_with_duplicate_aux_entries_normalizes_both() {
        // Doubly-degenerate: TWO integer AS entries in the aux AND a repeated key
        // in `normalize` ([AS, AS]). The sequential oracle calls
        // normalize_int_tag_to_smallest_signed(AS) twice; each call peels off and
        // relocates the first remaining AS, so BOTH aux entries are normalized and
        // moved to the tail (AS=300 fits i16, AS=50 fits i8). The single-pass fast
        // path captures only the first occurrence per key, so this MUST take the
        // sequential fallback — assert byte parity with the oracle.
        let as1 = ent(key(SamTag::AS), b'i', &300i32.to_le_bytes());
        let as2 = ent(key(SamTag::AS), b'i', &50i32.to_le_bytes());
        let base = rec_with_aux(&[&as1, b"RGZg\x00", &as2]);
        assert_folded_matches_oracle(
            &base,
            &[] as &[[u8; 2]],
            &[],
            &[key(SamTag::AS), key(SamTag::AS)],
            "duplicate AS aux entries with repeated AS normalize key",
        );
    }

    #[test]
    fn folded_normalize_duplicate_key_both_int_keeps_first_match_semantics() {
        // Degenerate (out-of-BAM-spec) aux with AS twice. The oracle
        // (find_int_tag/remove_tag = first key match) normalizes the FIRST AS and
        // relocates it to the tail, leaving the second AS verbatim in place. The
        // folded path must reproduce that exactly.
        let as1 = ent(key(SamTag::AS), b'i', &300i32.to_le_bytes());
        let as2 = ent(key(SamTag::AS), b'i', &50i32.to_le_bytes());
        let base = rec_with_aux(&[&as1, b"RGZg\x00", &as2]);
        assert_folded_matches_oracle(
            &base,
            &[] as &[[u8; 2]],
            &[],
            &as_xs(),
            "duplicate AS, both int",
        );
    }

    #[test]
    fn folded_normalize_duplicate_key_non_int_first_left_verbatim() {
        // Doubly-degenerate: AS appears twice, the first occurrence non-integer.
        // find_int_tag matches the first key (non-int) and returns None, so the
        // oracle normalizes nothing and leaves BOTH entries verbatim — the folded
        // path must NOT reach past the non-int first occurrence to the int second.
        let as_z = ent(key(SamTag::AS), b'Z', b"x\x00");
        let as_i = ent(key(SamTag::AS), b'i', &50i32.to_le_bytes());
        let base = rec_with_aux(&[&as_z, &as_i, b"RGZg\x00"]);
        assert_folded_matches_oracle(
            &base,
            &[] as &[[u8; 2]],
            &[],
            &as_xs(),
            "duplicate AS, non-int first",
        );
    }

    #[test]
    fn folded_normalize_adds_supplies_non_int_and_overflow_tags() {
        // `adds` carries AS as a non-integer and XS as a u32 above i32::MAX: both
        // are non-normalizable, so each must be appended verbatim (not relocated),
        // matching rebuild_with + a no-op normalize.
        let base = rec_with_aux(&[b"RGZg\x00"]);
        let xs_big = 3_000_000_000u32.to_le_bytes();
        let adds = [
            TagEntry { tag: key(SamTag::AS), type_byte: b'Z', value_bytes: b"hi\x00" },
            TagEntry { tag: key(SamTag::XS), type_byte: b'I', value_bytes: &xs_big },
        ];
        assert_folded_matches_oracle(
            &base,
            &[] as &[[u8; 2]],
            &adds,
            &as_xs(),
            "adds non-int/overflow",
        );
    }
}
