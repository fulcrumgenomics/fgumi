//! Direct field extraction from BAM binary records.
//!
//! This module provides functions to extract fields directly from BAM records
//! in their binary format, bypassing noodles' parsing for improved performance.
//!
//! # BAM Record Binary Layout
//!
//! ```text
//! Offset  Size  Field
//! ------  ----  -----
//! 0-3     4     refID (i32) - reference sequence ID
//! 4-7     4     pos (i32) - 0-based leftmost position
//! 8       1     l_read_name (u8) - length of read name + NUL
//! 9       1     mapq (u8) - mapping quality
//! 10-11   2     bin (u16) - BAM bin
//! 12-13   2     n_cigar_op (u16) - number of CIGAR operations
//! 14-15   2     flag (u16) - bitwise flags
//! 16-19   4     l_seq (u32) - length of sequence
//! 20-23   4     next_refID (i32) - mate reference sequence ID
//! 24-27   4     next_pos (i32) - mate 0-based position
//! 28-31   4     tlen (i32) - template length
//! 32+     var   read_name (l_read_name bytes, null-terminated)
//! ```

/// Minimum length of a valid BAM record (the 32-byte fixed header).
/// All read-primitive functions that access fixed-offset fields (`flags`, `mapq`,
/// `ref_id`, `pos`, etc.) require at least this many bytes.
pub const MIN_BAM_HEADER_LEN: usize = 32;

/// BAM flag bits.
pub mod flags {
    /// Read is paired in sequencing.
    pub const PAIRED: u16 = 0x1;
    /// Read is unmapped.
    pub const UNMAPPED: u16 = 0x4;
    /// Mate is unmapped.
    pub const MATE_UNMAPPED: u16 = 0x8;
    /// Read is reverse complemented.
    pub const REVERSE: u16 = 0x10;
    /// Mate is reverse complemented.
    pub const MATE_REVERSE: u16 = 0x20;
    /// First segment in template (R1).
    pub const FIRST_SEGMENT: u16 = 0x40;
    /// Last segment in template (R2).
    pub const LAST_SEGMENT: u16 = 0x80;
    /// Secondary alignment.
    pub const SECONDARY: u16 = 0x100;
    /// Not passing quality controls.
    pub const QC_FAIL: u16 = 0x200;
    /// PCR or optical duplicate.
    pub const DUPLICATE: u16 = 0x400;
    /// Supplementary alignment.
    pub const SUPPLEMENTARY: u16 = 0x800;
}

/// Extract coordinate key fields directly from BAM bytes.
///
/// Returns (tid, pos, reverse, name).
#[inline]
#[must_use]
pub fn extract_coordinate_fields(bam_bytes: &[u8]) -> (i32, i32, bool, &[u8]) {
    // Read fixed fields at known offsets
    let tid = i32::from_le_bytes([bam_bytes[0], bam_bytes[1], bam_bytes[2], bam_bytes[3]]);
    let pos = i32::from_le_bytes([bam_bytes[4], bam_bytes[5], bam_bytes[6], bam_bytes[7]]);
    let l_read_name = bam_bytes[8] as usize;
    let flag = u16::from_le_bytes([bam_bytes[14], bam_bytes[15]]);

    let reverse = (flag & flags::REVERSE) != 0;
    let unmapped = (flag & flags::UNMAPPED) != 0;

    // Read name starts at byte 32, exclude null terminator
    let name = if l_read_name > 1 { &bam_bytes[32..32 + l_read_name - 1] } else { &[] };

    if unmapped { (i32::MAX, i32::MAX, false, name) } else { (tid, pos, reverse, name) }
}

/// Fields extracted for template-coordinate sorting.
pub struct TemplateCoordFields<'a> {
    pub tid: i32,
    pub pos: i32,
    pub mate_tid: i32,
    pub mate_pos: i32,
    pub reverse: bool,
    pub mate_reverse: bool,
    pub unmapped: bool,
    pub mate_unmapped: bool,
    pub paired: bool,
    pub name: &'a [u8],
    pub l_read_name: usize,
    pub n_cigar_op: usize,
    pub l_seq: usize,
}

/// Extract template-coordinate key fields directly from BAM bytes.
#[inline]
#[must_use]
pub fn extract_template_coordinate_fields(bam_bytes: &[u8]) -> TemplateCoordFields<'_> {
    let tid = i32::from_le_bytes([bam_bytes[0], bam_bytes[1], bam_bytes[2], bam_bytes[3]]);
    let pos = i32::from_le_bytes([bam_bytes[4], bam_bytes[5], bam_bytes[6], bam_bytes[7]]);
    let l_read_name = bam_bytes[8] as usize;
    let n_cigar_op = u16::from_le_bytes([bam_bytes[12], bam_bytes[13]]) as usize;
    let flag = u16::from_le_bytes([bam_bytes[14], bam_bytes[15]]);
    let l_seq =
        u32::from_le_bytes([bam_bytes[16], bam_bytes[17], bam_bytes[18], bam_bytes[19]]) as usize;
    let mate_tid = i32::from_le_bytes([bam_bytes[20], bam_bytes[21], bam_bytes[22], bam_bytes[23]]);
    let mate_pos = i32::from_le_bytes([bam_bytes[24], bam_bytes[25], bam_bytes[26], bam_bytes[27]]);

    // Read name starts at byte 32, exclude null terminator
    let name = if l_read_name > 1 { &bam_bytes[32..32 + l_read_name - 1] } else { &[] };

    // Extract flags
    let unmapped = (flag & flags::UNMAPPED) != 0;
    let mate_unmapped = (flag & flags::MATE_UNMAPPED) != 0;
    let reverse = (flag & flags::REVERSE) != 0;
    let mate_reverse = (flag & flags::MATE_REVERSE) != 0;
    let paired = (flag & flags::PAIRED) != 0;

    TemplateCoordFields {
        tid,
        pos,
        mate_tid,
        mate_pos,
        reverse,
        mate_reverse,
        unmapped,
        mate_unmapped,
        paired,
        name,
        l_read_name,
        n_cigar_op,
        l_seq,
    }
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
                    if v >= 0 { Some((v as u64, true)) } else { None }
                }
                b'C' if pos + 4 <= aux_data.len() => Some((u64::from(aux_data[pos + 3]), true)),
                b's' if pos + 5 <= aux_data.len() => {
                    let v = i16::from_le_bytes([aux_data[pos + 3], aux_data[pos + 4]]);
                    if v >= 0 { Some((v as u64, true)) } else { None }
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
                    if v >= 0 { Some((v as u64, true)) } else { None }
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

/// Calculate the offset to auxiliary data in a BAM record.
///
/// `aux_offset = 32 + l_read_name + n_cigar_op*4 + (l_seq+1)/2 + l_seq`
#[inline]
#[must_use]
pub fn aux_data_offset(l_read_name: usize, n_cigar_op: usize, l_seq: usize) -> usize {
    32 + l_read_name + n_cigar_op * 4 + l_seq.div_ceil(2) + l_seq
}

/// Fixed-size tag value sizes indexed by type byte.
/// 0 = variable/unknown (needs special handling).
const TAG_FIXED_SIZES: [u8; 256] = {
    let mut table = [0u8; 256];
    table[b'A' as usize] = 1;
    table[b'c' as usize] = 1;
    table[b'C' as usize] = 1;
    table[b's' as usize] = 2;
    table[b'S' as usize] = 2;
    table[b'i' as usize] = 4;
    table[b'I' as usize] = 4;
    table[b'f' as usize] = 4;
    table
};

/// Calculate the size of a tag value based on its type.
#[inline]
#[must_use]
pub fn tag_value_size(val_type: u8, data: &[u8]) -> Option<usize> {
    let fixed = TAG_FIXED_SIZES[val_type as usize];
    if fixed > 0 {
        return Some(fixed as usize);
    }
    match val_type {
        b'Z' | b'H' => Some(data.iter().position(|&b| b == 0)? + 1),
        b'B' => {
            if data.len() < 5 {
                return None;
            }
            let elem_type = data[0];
            let count = u32::from_le_bytes([data[1], data[2], data[3], data[4]]) as usize;
            let elem_size = TAG_FIXED_SIZES[elem_type as usize] as usize;
            if elem_size == 0 {
                return None;
            }
            Some(5 + count * elem_size)
        }
        _ => None,
    }
}

// ============================================================================
// Read Primitives
// ============================================================================
//
// All read primitives below assume `bam.len() >= 32` (the BAM fixed-header
// size).  Passing a shorter slice will panic on out-of-bounds indexing.
// Callers must validate record length before invoking these functions.

/// Extract flags (u16) from a BAM record.
#[inline]
#[must_use]
pub fn flags(bam: &[u8]) -> u16 {
    u16::from_le_bytes([bam[14], bam[15]])
}

/// Extract mapping quality from a BAM record.
#[inline]
#[must_use]
pub fn mapq(bam: &[u8]) -> u8 {
    bam[9]
}

/// Extract reference sequence ID from a BAM record.
#[inline]
#[must_use]
pub fn ref_id(bam: &[u8]) -> i32 {
    i32::from_le_bytes([bam[0], bam[1], bam[2], bam[3]])
}

/// Extract 0-based leftmost position from a BAM record.
#[inline]
#[must_use]
pub fn pos(bam: &[u8]) -> i32 {
    i32::from_le_bytes([bam[4], bam[5], bam[6], bam[7]])
}

/// Extract `l_read_name` (length of read name + NUL) from a BAM record.
#[inline]
#[must_use]
pub fn l_read_name(bam: &[u8]) -> u8 {
    bam[8]
}

/// Extract number of CIGAR operations from a BAM record.
#[inline]
#[must_use]
pub fn n_cigar_op(bam: &[u8]) -> u16 {
    u16::from_le_bytes([bam[12], bam[13]])
}

/// Extract sequence length from a BAM record.
#[inline]
#[must_use]
pub fn l_seq(bam: &[u8]) -> u32 {
    u32::from_le_bytes([bam[16], bam[17], bam[18], bam[19]])
}

/// Extract mate reference sequence ID from a BAM record.
#[inline]
#[must_use]
pub fn mate_ref_id(bam: &[u8]) -> i32 {
    i32::from_le_bytes([bam[20], bam[21], bam[22], bam[23]])
}

/// Extract mate 0-based position from a BAM record.
#[inline]
#[must_use]
pub fn mate_pos(bam: &[u8]) -> i32 {
    i32::from_le_bytes([bam[24], bam[25], bam[26], bam[27]])
}

/// Extract read name (without null terminator) from a BAM record.
#[inline]
#[must_use]
pub fn read_name(bam: &[u8]) -> &[u8] {
    let l = bam[8] as usize;
    if l > 1 { &bam[32..32 + l - 1] } else { &[] }
}

/// Calculate the offset to auxiliary data for a complete BAM record.
///
/// Returns `None` if the record is too short to read the required header
/// fields (needs at least 20 bytes for `l_read_name`, `n_cigar_op`, `l_seq`).
#[inline]
#[must_use]
pub fn aux_data_offset_from_record(bam: &[u8]) -> Option<usize> {
    if bam.len() < 20 {
        return None;
    }
    let l_rn = bam[8] as usize;
    let n_co = u16::from_le_bytes([bam[12], bam[13]]) as usize;
    let l_s = u32::from_le_bytes([bam[16], bam[17], bam[18], bam[19]]) as usize;
    Some(aux_data_offset(l_rn, n_co, l_s))
}

/// Get auxiliary data as a byte slice from a complete BAM record.
///
/// Returns an empty slice for truncated records or records with no aux data.
#[inline]
#[must_use]
pub fn aux_data_slice(bam: &[u8]) -> &[u8] {
    match aux_data_offset_from_record(bam) {
        Some(offset) if offset <= bam.len() => &bam[offset..],
        _ => &[],
    }
}

/// Calculate the offset to sequence data in a BAM record.
#[inline]
#[must_use]
pub fn seq_offset(bam: &[u8]) -> usize {
    let l_rn = bam[8] as usize;
    let n_co = u16::from_le_bytes([bam[12], bam[13]]) as usize;
    32 + l_rn + n_co * 4
}

/// Calculate the offset to quality data in a BAM record.
#[inline]
#[must_use]
pub fn qual_offset(bam: &[u8]) -> usize {
    let l_rn = bam[8] as usize;
    let n_co = u16::from_le_bytes([bam[12], bam[13]]) as usize;
    let l_s = u32::from_le_bytes([bam[16], bam[17], bam[18], bam[19]]) as usize;
    32 + l_rn + n_co * 4 + l_s.div_ceil(2)
}

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

/// Extract a 4-bit base from packed sequence data.
///
/// BAM packs two bases per byte: high nibble = even index, low nibble = odd index.
/// Returns the 4-bit encoded base (1=A, 2=C, 4=G, 8=T, 15=N).
#[inline]
#[must_use]
pub fn get_base(bam: &[u8], seq_off: usize, position: usize) -> u8 {
    let byte = bam[seq_off + position / 2];
    if position.is_multiple_of(2) { byte >> 4 } else { byte & 0xF }
}

/// Extract a quality score at a given position.
#[inline]
#[must_use]
pub fn get_qual(bam: &[u8], qual_off: usize, position: usize) -> u8 {
    bam[qual_off + position]
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
                b'c' if p + 4 <= aux_data.len() => Some(i64::from(aux_data[p + 3] as i8)),
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

// ============================================================================
// Write Primitives
// ============================================================================

/// Set flags (u16) in a BAM record.
#[inline]
pub fn set_flags(bam: &mut [u8], new_flags: u16) {
    bam[14..16].copy_from_slice(&new_flags.to_le_bytes());
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
        record.push(v as u8);
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
pub fn append_i16_array_tag(record: &mut Vec<u8>, tag: &[u8; 2], values: &[i16]) {
    record.push(tag[0]);
    record.push(tag[1]);
    record.push(b'B');
    record.push(b's');
    record.extend_from_slice(&(values.len() as u32).to_le_bytes());
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

/// Set the 4-bit base at a position to N (0xF) in packed sequence data.
#[inline]
pub fn mask_base(bam: &mut [u8], seq_off: usize, position: usize) {
    let byte_idx = seq_off + position / 2;
    if position.is_multiple_of(2) {
        bam[byte_idx] = (bam[byte_idx] & 0x0F) | 0xF0; // high nibble
    } else {
        bam[byte_idx] = (bam[byte_idx] & 0xF0) | 0x0F; // low nibble
    }
}

/// Set the quality score at a given position.
#[inline]
pub fn set_qual(bam: &mut [u8], qual_off: usize, position: usize, value: u8) {
    bam[qual_off + position] = value;
}

// ============================================================================
// Sequence Packing (ASCII bases -> 4-bit packed BAM format)
// ============================================================================

/// Lookup table mapping ASCII base characters to 4-bit BAM codes.
///
/// Duplicated from `vendored/bam_codec/encoder/sequence.rs` to keep
/// `bam_fields` self-contained. The table is `const` so there is zero
/// runtime cost.
const SEQ_CODES: [u8; 256] = build_seq_codes();

const fn build_seq_codes() -> [u8; 256] {
    // SAM spec §4.2.3: =ACMGRSVTWYHKDBN -> 0..15
    const BASES: [u8; 16] = *b"=ACMGRSVTWYHKDBN";
    const N: u8 = 0x0F;
    let mut codes = [N; 256];
    let mut i = 0;
    while i < BASES.len() {
        let base = BASES[i];
        let code = i as u8;
        codes[base as usize] = code;
        codes[base.to_ascii_lowercase() as usize] = code;
        i += 1;
    }
    codes
}

/// Pack ASCII bases into BAM 4-bit-per-base format, appending to `dst`.
///
/// Uses 16-base chunked processing for cache efficiency, matching the
/// htslib/vendored encoder strategy. When `l_seq` is odd the bottom
/// 4 bits of the last byte are zero-padded per the SAM spec.
#[inline]
pub fn pack_sequence_into(dst: &mut Vec<u8>, bases: &[u8]) {
    if bases.is_empty() {
        return;
    }
    let packed_len = bases.len().div_ceil(2);
    dst.reserve(packed_len);

    const CHUNK: usize = 16;
    let mut chunks = bases.chunks_exact(CHUNK);
    for chunk in chunks.by_ref() {
        dst.push((SEQ_CODES[chunk[0] as usize] << 4) | SEQ_CODES[chunk[1] as usize]);
        dst.push((SEQ_CODES[chunk[2] as usize] << 4) | SEQ_CODES[chunk[3] as usize]);
        dst.push((SEQ_CODES[chunk[4] as usize] << 4) | SEQ_CODES[chunk[5] as usize]);
        dst.push((SEQ_CODES[chunk[6] as usize] << 4) | SEQ_CODES[chunk[7] as usize]);
        dst.push((SEQ_CODES[chunk[8] as usize] << 4) | SEQ_CODES[chunk[9] as usize]);
        dst.push((SEQ_CODES[chunk[10] as usize] << 4) | SEQ_CODES[chunk[11] as usize]);
        dst.push((SEQ_CODES[chunk[12] as usize] << 4) | SEQ_CODES[chunk[13] as usize]);
        dst.push((SEQ_CODES[chunk[14] as usize] << 4) | SEQ_CODES[chunk[15] as usize]);
    }

    let remainder = chunks.remainder();
    let mut pairs = remainder.chunks_exact(2);
    for pair in pairs.by_ref() {
        dst.push((SEQ_CODES[pair[0] as usize] << 4) | SEQ_CODES[pair[1] as usize]);
    }
    if let Some(&last) = pairs.remainder().first() {
        dst.push(SEQ_CODES[last as usize] << 4);
    }
}

// ============================================================================
// Unmapped BAM Record Builder
// ============================================================================

/// Reusable builder for constructing unmapped BAM records as raw bytes.
///
/// Produces BAM record bytes directly, bypassing noodles `RecordBuf`
/// serialization. Designed for the consensus output hot path where
/// records are always unmapped (`ref_id=-1`, `pos=-1`, no CIGAR).
///
/// The internal buffer is reused across calls to [`Self::build_record`] +
/// [`Self::clear`] to avoid repeated allocation.
///
/// # Usage
///
/// ```rust,ignore
/// let mut builder = UnmappedBamRecordBuilder::new();
///
/// builder.build_record(b"cons:1:ACG-TCG", my_flags, &bases, &quals);
/// builder.append_string_tag(b"RG", b"sample1");
/// builder.append_int_tag(b"cD", max_depth);
/// builder.append_float_tag(b"cE", error_rate);
/// builder.append_i16_array_tag(b"cd", &depth_array);
/// builder.write_with_block_size(&mut output);
///
/// builder.clear();
/// // ... build next record ...
/// ```
pub struct UnmappedBamRecordBuilder {
    buf: Vec<u8>,
    sealed: bool,
}

/// Unmapped BAM bin (SAM spec §4.2.1: `reg2bin(-1, 0)` = 4680).
const UNMAPPED_BIN: u16 = 4680;

impl UnmappedBamRecordBuilder {
    /// Create a new builder with default capacity (512 bytes).
    #[must_use]
    pub fn new() -> Self {
        Self::with_capacity(512)
    }

    /// Create a new builder with the given initial capacity.
    #[must_use]
    pub fn with_capacity(capacity: usize) -> Self {
        Self { buf: Vec::with_capacity(capacity), sealed: false }
    }

    /// Clear the builder for reuse without deallocating.
    pub fn clear(&mut self) {
        self.buf.clear();
        self.sealed = false;
    }

    /// Write the 32-byte fixed header, read name, packed sequence, and
    /// quality scores for an unmapped record.
    ///
    /// After calling this, append tags via the `append_*_tag` methods,
    /// then call [`Self::write_with_block_size`] or [`Self::as_bytes`].
    ///
    /// Hardcoded unmapped fields: `ref_id=-1`, `pos=-1`, `mapq=0`,
    /// `bin=4680`, `n_cigar_op=0`, `next_ref_id=-1`, `next_pos=-1`,
    /// `tlen=0`.
    ///
    /// # Arguments
    ///
    /// * `name`  — read name **without** null terminator
    /// * `flag`  — BAM flags (`u16`)
    /// * `bases` — ASCII sequence (e.g. `b"ACGT"`)
    /// * `quals` — raw Phred quality scores (0-93, **not** Phred+33)
    ///
    /// # Panics
    ///
    /// Panics if `quals` is non-empty and `bases.len() != quals.len()`.
    pub fn build_record(&mut self, name: &[u8], flag: u16, bases: &[u8], quals: &[u8]) {
        assert!(
            bases.len() == quals.len() || quals.is_empty(),
            "bases.len() ({}) != quals.len() ({})",
            bases.len(),
            quals.len(),
        );

        self.buf.clear();
        self.sealed = false;

        assert!(name.len() < 255, "read name too long ({} bytes, max 254)", name.len());
        let l_read_name = (name.len() + 1) as u8; // +1 for NUL
        let l_seq = bases.len() as u32;
        let packed_seq_len = bases.len().div_ceil(2);

        // Pre-reserve for header + name + seq + qual + typical tag overhead
        self.buf.reserve(32 + l_read_name as usize + packed_seq_len + bases.len() + 100);

        // === Fixed 32-byte header ===
        self.buf.extend_from_slice(&(-1i32).to_le_bytes()); // ref_id = -1
        self.buf.extend_from_slice(&(-1i32).to_le_bytes()); // pos = -1
        self.buf.push(l_read_name); // l_read_name
        self.buf.push(0); // mapq = 0
        self.buf.extend_from_slice(&UNMAPPED_BIN.to_le_bytes()); // bin = 4680
        self.buf.extend_from_slice(&0u16.to_le_bytes()); // n_cigar_op = 0
        self.buf.extend_from_slice(&flag.to_le_bytes()); // flags
        self.buf.extend_from_slice(&l_seq.to_le_bytes()); // l_seq
        self.buf.extend_from_slice(&(-1i32).to_le_bytes()); // next_ref_id = -1
        self.buf.extend_from_slice(&(-1i32).to_le_bytes()); // next_pos = -1
        self.buf.extend_from_slice(&0i32.to_le_bytes()); // tlen = 0

        // === Read name + NUL ===
        self.buf.extend_from_slice(name);
        self.buf.push(0);

        // === Packed sequence ===
        pack_sequence_into(&mut self.buf, bases);

        // === Quality scores ===
        if quals.is_empty() && !bases.is_empty() {
            // Missing quality: 0xFF per SAM spec
            self.buf.resize(self.buf.len() + bases.len(), 0xFF);
        } else {
            self.buf.extend_from_slice(quals);
        }

        self.sealed = true;
    }

    /// Append a string (`Z`-type) tag.
    #[inline]
    pub fn append_string_tag(&mut self, tag: &[u8; 2], value: &[u8]) {
        debug_assert!(self.sealed, "must call build_record before appending tags");
        append_string_tag(&mut self.buf, tag, value);
    }

    /// Append an integer tag using the smallest signed type that fits.
    #[inline]
    pub fn append_int_tag(&mut self, tag: &[u8; 2], value: i32) {
        debug_assert!(self.sealed, "must call build_record before appending tags");
        append_int_tag(&mut self.buf, tag, value);
    }

    /// Append a float (`f`-type) tag.
    #[inline]
    pub fn append_float_tag(&mut self, tag: &[u8; 2], value: f32) {
        debug_assert!(self.sealed, "must call build_record before appending tags");
        append_float_tag(&mut self.buf, tag, value);
    }

    /// Append an `i16` array (`B:s`-type) tag.
    #[inline]
    pub fn append_i16_array_tag(&mut self, tag: &[u8; 2], values: &[i16]) {
        debug_assert!(self.sealed, "must call build_record before appending tags");
        append_i16_array_tag(&mut self.buf, tag, values);
    }

    /// Append a Phred+33 encoded quality string (`Z`-type) tag.
    #[inline]
    pub fn append_phred33_string_tag(&mut self, tag: &[u8; 2], quals: &[u8]) {
        debug_assert!(self.sealed, "must call build_record before appending tags");
        append_phred33_string_tag(&mut self.buf, tag, quals);
    }

    /// Get the completed record bytes (**without** the 4-byte `block_size` prefix).
    #[inline]
    #[must_use]
    pub fn as_bytes(&self) -> &[u8] {
        debug_assert!(self.sealed, "must call build_record first");
        &self.buf
    }

    /// Write the record with a 4-byte little-endian `block_size` prefix to `output`.
    ///
    /// This matches the output format expected by the pipeline serialize
    /// functions and the BGZF compression step.
    #[inline]
    pub fn write_with_block_size(&self, output: &mut Vec<u8>) {
        debug_assert!(self.sealed, "must call build_record first");
        let block_size = self.buf.len() as u32;
        output.extend_from_slice(&block_size.to_le_bytes());
        output.extend_from_slice(&self.buf);
    }
}

impl Default for UnmappedBamRecordBuilder {
    fn default() -> Self {
        Self::new()
    }
}

/// Parsed BAM record for test assertions.
///
/// Provides convenient access to fields extracted from raw BAM bytes
/// produced by `UnmappedBamRecordBuilder`. Used only in tests.
#[cfg(test)]
pub struct ParsedBamRecord {
    pub name: Vec<u8>,
    pub flag: u16,
    pub bases: Vec<u8>,
    pub quals: Vec<u8>,
    pub aux_data: Vec<u8>,
}

#[cfg(test)]
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

#[cfg(test)]
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

#[cfg(test)]
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
                    b'c' | b'C' => 1,
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

#[cfg(test)]
fn find_int_tag_in_aux(aux: &[u8], tag: [u8; 2]) -> Option<i32> {
    let mut i = 0;
    while i + 3 <= aux.len() {
        let t = [aux[i], aux[i + 1]];
        let typ = aux[i + 2];
        i += 3;
        match typ {
            b'c' => {
                let v = i32::from(aux[i] as i8);
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
                    b'c' | b'C' => 1,
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

#[cfg(test)]
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
                    b'c' | b'C' => 1,
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

#[cfg(test)]
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
                    b'c' | b'C' => 1,
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

// ============================================================================
// Coordinate Helpers (1-based, matching noodles)
// ============================================================================

/// Calculate unclipped 5' coordinate using 1-based positions (matching noodles).
///
/// Adds 1 to the 0-based BAM position before applying clip adjustment,
/// producing values identical to `record_utils::unclipped_five_prime_position()`.
///
/// Returns 0 for unmapped reads (matching `get_unclipped_position_for_groupkey`).
/// Returns `i32::MAX` for mapped reads with no CIGAR operations.
#[inline]
#[must_use]
pub fn unclipped_5prime_1based(
    pos_0based: i32,
    reverse: bool,
    unmapped: bool,
    cigar_ops: &[u32],
) -> i32 {
    if unmapped {
        return 0;
    }
    if cigar_ops.is_empty() {
        return i32::MAX;
    }
    // Convert to 1-based, then apply clip adjustment
    let pos_1based = pos_0based + 1;
    if reverse {
        unclipped_end_from_cigar(pos_1based, cigar_ops)
    } else {
        unclipped_start_from_cigar(pos_1based, cigar_ops)
    }
}

/// Calculate mate's unclipped 5' coordinate using 1-based positions (matching noodles).
///
/// Adds 1 to the 0-based BAM mate position before applying clip adjustment.
#[inline]
#[must_use]
pub fn mate_unclipped_5prime_1based(
    mate_pos_0based: i32,
    mate_reverse: bool,
    mc_cigar: &str,
) -> i32 {
    let mate_pos_1based = mate_pos_0based + 1;
    if mate_reverse {
        unclipped_other_end(mate_pos_1based, mc_cigar)
    } else {
        unclipped_other_start(mate_pos_1based, mc_cigar)
    }
}

// ============================================================================
// Raw Bytes Comparison Functions (Optimization 5)
// ============================================================================

use std::cmp::Ordering;

/// Compare two BAM records for coordinate ordering using raw bytes.
///
/// This is faster than extracting keys into structs because it avoids
/// allocations and directly reads from the BAM binary format.
///
/// For coordinate sorting (following samtools behavior):
/// - Reads with valid tid (>= 0) are sorted by (tid, pos, reverse, name)
/// - Reads with tid = -1 (no reference) sort at the end by name
/// - Unmapped reads with a valid tid use that tid for sorting (typically mate's position)
#[inline]
#[must_use]
pub fn compare_coordinate_raw(a: &[u8], b: &[u8]) -> Ordering {
    // Extract fields at fixed offsets
    let a_tid = i32::from_le_bytes([a[0], a[1], a[2], a[3]]);
    let b_tid = i32::from_le_bytes([b[0], b[1], b[2], b[3]]);

    let a_pos = i32::from_le_bytes([a[4], a[5], a[6], a[7]]);
    let b_pos = i32::from_le_bytes([b[4], b[5], b[6], b[7]]);

    let a_flag = u16::from_le_bytes([a[14], a[15]]);
    let b_flag = u16::from_le_bytes([b[14], b[15]]);

    // Handle reads with no reference (tid = -1) - sort last
    // Unmapped reads with a valid tid (mate's position) sort by that position
    let a_no_ref = a_tid < 0;
    let b_no_ref = b_tid < 0;

    match (a_no_ref, b_no_ref) {
        (true, false) => return Ordering::Greater,
        (false, true) => return Ordering::Less,
        (true, true) => return compare_names_raw(a, b),
        _ => {}
    }

    // Compare tid, pos, reverse, then name
    a_tid
        .cmp(&b_tid)
        .then_with(|| a_pos.cmp(&b_pos))
        .then_with(|| {
            let a_rev = (a_flag & flags::REVERSE) != 0;
            let b_rev = (b_flag & flags::REVERSE) != 0;
            a_rev.cmp(&b_rev)
        })
        .then_with(|| compare_names_raw(a, b))
}

/// Compare read names directly from BAM bytes.
#[inline]
#[must_use]
pub fn compare_names_raw(a: &[u8], b: &[u8]) -> Ordering {
    let a_name_len = a[8] as usize;
    let b_name_len = b[8] as usize;

    // Exclude null terminator
    let a_len = if a_name_len > 0 { a_name_len - 1 } else { 0 };
    let b_len = if b_name_len > 0 { b_name_len - 1 } else { 0 };

    let a_name = &a[32..32 + a_len];
    let b_name = &b[32..32 + b_len];

    a_name.cmp(b_name)
}

/// Zero-allocation unclipped 5' position from raw BAM bytes (0-based, for sorting).
///
/// Unlike [`unclipped_5prime_from_raw_bam`] which adds 1 for 1-based group keys,
/// this uses the raw 0-based BAM position directly — matching what the sort
/// comparator needs.
#[inline]
fn unclipped_5prime_sort(bam: &[u8], pos_0based: i32, reverse: bool) -> i32 {
    let n_cigar_op = u16::from_le_bytes([bam[12], bam[13]]) as usize;
    if n_cigar_op == 0 {
        return pos_0based;
    }
    let l_read_name = bam[8] as usize;
    let cigar_start = 32 + l_read_name;
    if cigar_start + n_cigar_op * 4 > bam.len() {
        return pos_0based;
    }
    if reverse {
        unclipped_end_from_raw_cigar(pos_0based, bam, cigar_start, n_cigar_op)
    } else {
        unclipped_start_from_raw_cigar(pos_0based, bam, cigar_start, n_cigar_op)
    }
}

/// Compare for template-coordinate ordering using raw bytes.
///
/// This matches samtools' template-coordinate sorting which uses unclipped 5' positions.
#[inline]
#[must_use]
pub fn compare_template_coordinate_raw(a: &[u8], b: &[u8]) -> Ordering {
    // Extract all needed fields from both records
    let a_tid = i32::from_le_bytes([a[0], a[1], a[2], a[3]]);
    let a_pos = i32::from_le_bytes([a[4], a[5], a[6], a[7]]);
    let a_flag = u16::from_le_bytes([a[14], a[15]]);
    let a_mate_tid = i32::from_le_bytes([a[20], a[21], a[22], a[23]]);
    let a_mate_pos = i32::from_le_bytes([a[24], a[25], a[26], a[27]]);

    let b_tid = i32::from_le_bytes([b[0], b[1], b[2], b[3]]);
    let b_pos = i32::from_le_bytes([b[4], b[5], b[6], b[7]]);
    let b_flag = u16::from_le_bytes([b[14], b[15]]);
    let b_mate_tid = i32::from_le_bytes([b[20], b[21], b[22], b[23]]);
    let b_mate_pos = i32::from_le_bytes([b[24], b[25], b[26], b[27]]);

    // Extract strand information
    let a_reverse = (a_flag & flags::REVERSE) != 0;
    let b_reverse = (b_flag & flags::REVERSE) != 0;
    let a_mate_reverse = (a_flag & flags::MATE_REVERSE) != 0;
    let b_mate_reverse = (b_flag & flags::MATE_REVERSE) != 0;

    // Get MC tags for mate's unclipped position
    let a_mc = find_mc_tag_in_record(a);
    let b_mc = find_mc_tag_in_record(b);

    // Calculate unclipped 5' positions (zero-allocation: reads CIGAR from raw bytes)
    let a_unclipped_pos = unclipped_5prime_sort(a, a_pos, a_reverse);
    let b_unclipped_pos = unclipped_5prime_sort(b, b_pos, b_reverse);

    let a_mate_unclipped_pos =
        a_mc.map(|mc| mate_unclipped_5prime(a_mate_pos, a_mate_reverse, mc)).unwrap_or(a_mate_pos);
    let b_mate_unclipped_pos =
        b_mc.map(|mc| mate_unclipped_5prime(b_mate_pos, b_mate_reverse, mc)).unwrap_or(b_mate_pos);

    // Compute canonical positions using unclipped coordinates
    let (a_tid1, a_tid2, a_pos1, a_pos2, a_neg1, a_neg2, a_upper) =
        canonical_template_pos_unclipped(
            a_tid,
            a_unclipped_pos,
            a_mate_tid,
            a_mate_unclipped_pos,
            a_flag,
            a_reverse,
            a_mate_reverse,
        );
    let (b_tid1, b_tid2, b_pos1, b_pos2, b_neg1, b_neg2, b_upper) =
        canonical_template_pos_unclipped(
            b_tid,
            b_unclipped_pos,
            b_mate_tid,
            b_mate_unclipped_pos,
            b_flag,
            b_reverse,
            b_mate_reverse,
        );

    // Compare canonical positions (samtools order)
    a_tid1
        .cmp(&b_tid1)
        .then_with(|| a_tid2.cmp(&b_tid2))
        .then_with(|| a_pos1.cmp(&b_pos1))
        .then_with(|| a_pos2.cmp(&b_pos2))
        .then_with(|| {
            // Reverse strand sorts before forward (samtools convention)
            match (a_neg1, b_neg1) {
                (true, false) => Ordering::Less,
                (false, true) => Ordering::Greater,
                _ => Ordering::Equal,
            }
        })
        .then_with(|| match (a_neg2, b_neg2) {
            (true, false) => Ordering::Less,
            (false, true) => Ordering::Greater,
            _ => Ordering::Equal,
        })
        .then_with(|| compare_mi_tags_raw(a, b))
        .then_with(|| compare_names_raw(a, b))
        .then_with(|| a_upper.cmp(&b_upper))
}

/// Compute canonical template position using pre-computed unclipped positions.
///
/// This is used for template-coordinate sorting where positions are already
/// computed as unclipped 5' coordinates.
///
/// Returns `(tid1, tid2, pos1, pos2, neg1, neg2, is_upper)`.
#[inline]
fn canonical_template_pos_unclipped(
    tid: i32,
    unclipped_pos: i32,
    mate_tid: i32,
    mate_unclipped_pos: i32,
    flag: u16,
    reverse: bool,
    mate_reverse: bool,
) -> (i32, i32, i32, i32, bool, bool, bool) {
    let unmapped = (flag & flags::UNMAPPED) != 0;
    let mate_unmapped = (flag & flags::MATE_UNMAPPED) != 0;
    let paired = (flag & flags::PAIRED) != 0;

    if unmapped && (!paired || mate_unmapped) {
        // Completely unmapped (no mapped mate) - sort to end
        (i32::MAX, i32::MAX, i32::MAX, i32::MAX, false, false, false)
    } else if unmapped {
        // Unmapped read with mapped mate - use mate's position as primary key
        // This keeps unmapped reads with their mapped mates (samtools behavior)
        (mate_tid, i32::MAX, mate_unclipped_pos, i32::MAX, mate_reverse, false, true)
    } else if !paired || mate_unmapped {
        // Mapped read with unmapped mate - use MAX for tid2/pos2 (samtools behavior)
        (tid, i32::MAX, unclipped_pos, i32::MAX, reverse, false, false)
    } else {
        // Both reads mapped - canonical ordering
        // Samtools logic: is_upper if tid > mate_tid, or pos > mate_pos, or (pos equal and reverse)
        let is_upper = (tid, unclipped_pos) > (mate_tid, mate_unclipped_pos)
            || ((tid, unclipped_pos) == (mate_tid, mate_unclipped_pos) && reverse);

        if is_upper {
            // Swap: mate's position comes first
            (mate_tid, tid, mate_unclipped_pos, unclipped_pos, mate_reverse, reverse, true)
        } else {
            // No swap: this read's position comes first
            (tid, mate_tid, unclipped_pos, mate_unclipped_pos, reverse, mate_reverse, false)
        }
    }
}

/// Compare MI tags directly from raw BAM bytes.
///
/// Uses the pre-parsed integer values from `find_mi_tag_in_record`.
#[inline]
fn compare_mi_tags_raw(a: &[u8], b: &[u8]) -> Ordering {
    let (a_mi, _) = find_mi_tag_in_record(a);
    let (b_mi, _) = find_mi_tag_in_record(b);

    // Compare as integers (0 means no MI tag found)
    a_mi.cmp(&b_mi)
}

/// Compare for queryname ordering using raw bytes.
#[inline]
#[must_use]
pub fn compare_queryname_raw(a: &[u8], b: &[u8]) -> Ordering {
    compare_names_raw(a, b).then_with(|| {
        let a_flag = u16::from_le_bytes([a[14], a[15]]);
        let b_flag = u16::from_le_bytes([b[14], b[15]]);
        a_flag.cmp(&b_flag)
    })
}

// ============================================================================
// Unclipped Position Calculation (for template-coordinate sorting)
// ============================================================================

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

/// Parse leading clips (S or H) from a CIGAR string.
///
/// Returns the total length of leading soft/hard clips.
#[inline]
fn parse_leading_clips(cigar: &str) -> i32 {
    let mut clipped = 0i32;
    let mut num_start = 0;
    let bytes = cigar.as_bytes();

    for (i, &c) in bytes.iter().enumerate() {
        if c.is_ascii_digit() {
            continue;
        }
        // Parse the number before this operation
        let num: i32 = cigar[num_start..i].parse().unwrap_or(0);

        if c == b'S' || c == b'H' {
            clipped += num;
            num_start = i + 1;
        } else {
            // Non-clip operation, stop processing
            break;
        }
    }

    clipped
}

/// Parse reference length and trailing clips from a CIGAR string.
///
/// Returns `(reference_length, trailing_clips)`.
/// Reference length is the sum of M/D/N/=/X operations.
/// Trailing clips are S/H operations after the last reference-consuming op.
#[inline]
fn parse_ref_len_and_trailing_clips(cigar: &str) -> (i32, i32) {
    let mut ref_len = 0i32;
    let mut trailing_clips = 0i32;
    let mut num_start = 0;
    let mut saw_ref_op = false;
    let bytes = cigar.as_bytes();

    for (i, &c) in bytes.iter().enumerate() {
        if c.is_ascii_digit() {
            continue;
        }

        let num: i32 = cigar[num_start..i].parse().unwrap_or(0);
        num_start = i + 1;

        match c {
            b'M' | b'D' | b'N' | b'=' | b'X' => {
                ref_len += num;
                trailing_clips = 0; // Reset trailing clips
                saw_ref_op = true;
            }
            b'S' | b'H' if saw_ref_op => {
                trailing_clips += num;
            }
            _ => {}
        }
    }

    (ref_len, trailing_clips)
}

/// Calculate unclipped start position from BAM record.
///
/// For forward strand reads, this is the 5' position.
/// Returns: `pos - leading_clips` (0-based like the input position).
#[inline]
#[must_use]
pub fn unclipped_start_from_cigar(pos: i32, cigar_ops: &[u32]) -> i32 {
    let mut clipped = 0i32;

    for &op in cigar_ops {
        let op_len = (op >> 4) as i32;
        let op_type = (op & 0xF) as u8;

        match op_type {
            4 | 5 => clipped += op_len, // S (4) or H (5)
            _ => break,                 // Non-clip operation
        }
    }

    pos - clipped
}

/// Calculate unclipped end position from BAM record.
///
/// For reverse strand reads, this is the 5' position.
/// Returns: `pos + ref_length + trailing_clips - 1` (0-based).
#[inline]
#[must_use]
pub fn unclipped_end_from_cigar(pos: i32, cigar_ops: &[u32]) -> i32 {
    let mut ref_len = 0i32;
    let mut trailing_clips = 0i32;
    let mut saw_ref_op = false;

    for &op in cigar_ops {
        let op_len = (op >> 4) as i32;
        let op_type = (op & 0xF) as u8;

        match op_type {
            0 | 2 | 3 | 7 | 8 => {
                // M (0), D (2), N (3), = (7), X (8)
                ref_len += op_len;
                trailing_clips = 0;
                saw_ref_op = true;
            }
            4 | 5 if saw_ref_op => {
                // S (4) or H (5) after ref-consuming ops
                trailing_clips += op_len;
            }
            _ => {}
        }
    }

    pos + ref_len + trailing_clips - 1
}

/// Calculate mate's unclipped start from MC tag CIGAR string.
///
/// For forward strand mates, this is the 5' position.
#[inline]
#[must_use]
pub fn unclipped_other_start(mate_pos: i32, mc_cigar: &str) -> i32 {
    mate_pos - parse_leading_clips(mc_cigar)
}

/// Calculate mate's unclipped end from MC tag CIGAR string.
///
/// For reverse strand mates, this is the 5' position.
#[inline]
#[must_use]
pub fn unclipped_other_end(mate_pos: i32, mc_cigar: &str) -> i32 {
    let (ref_len, trailing_clips) = parse_ref_len_and_trailing_clips(mc_cigar);
    mate_pos + ref_len + trailing_clips - 1
}

/// Extract CIGAR operations from BAM record.
#[inline]
#[must_use]
pub fn get_cigar_ops(bam: &[u8]) -> Vec<u32> {
    let l_read_name = bam[8] as usize;
    let n_cigar_op = u16::from_le_bytes([bam[12], bam[13]]) as usize;

    if n_cigar_op == 0 {
        return Vec::new();
    }

    let cigar_start = 32 + l_read_name;
    let cigar_end = cigar_start + n_cigar_op * 4;

    if cigar_end > bam.len() {
        return Vec::new();
    }

    // Read CIGAR ops bytewise to avoid alignment issues — the CIGAR data offset
    // (32 + l_read_name) is not guaranteed to be 4-byte aligned.
    let cigar_bytes = &bam[cigar_start..cigar_end];
    cigar_bytes.chunks_exact(4).map(|c| u32::from_le_bytes([c[0], c[1], c[2], c[3]])).collect()
}

/// Read a single CIGAR op directly from raw BAM bytes at the given byte offset.
#[inline]
fn cigar_op_at(bam: &[u8], offset: usize) -> u32 {
    u32::from_le_bytes([bam[offset], bam[offset + 1], bam[offset + 2], bam[offset + 3]])
}

/// Compute unclipped start by iterating raw CIGAR bytes directly (zero allocation).
///
/// Equivalent to `unclipped_start_from_cigar(pos, &get_cigar_ops(bam))` but avoids
/// the `Vec<u32>` allocation.
#[inline]
fn unclipped_start_from_raw_cigar(pos: i32, bam: &[u8], cigar_start: usize, n_ops: usize) -> i32 {
    let mut clipped = 0i32;
    for i in 0..n_ops {
        let op = cigar_op_at(bam, cigar_start + i * 4);
        let op_type = (op & 0xF) as u8;
        match op_type {
            4 | 5 => clipped += (op >> 4) as i32, // S or H
            _ => break,
        }
    }
    pos - clipped
}

/// Compute unclipped end by iterating raw CIGAR bytes directly (zero allocation).
///
/// Equivalent to `unclipped_end_from_cigar(pos, &get_cigar_ops(bam))` but avoids
/// the `Vec<u32>` allocation.
#[inline]
fn unclipped_end_from_raw_cigar(pos: i32, bam: &[u8], cigar_start: usize, n_ops: usize) -> i32 {
    let mut ref_len = 0i32;
    let mut trailing_clips = 0i32;
    let mut saw_ref_op = false;

    for i in 0..n_ops {
        let op = cigar_op_at(bam, cigar_start + i * 4);
        let op_len = (op >> 4) as i32;
        let op_type = (op & 0xF) as u8;

        match op_type {
            0 | 2 | 3 | 7 | 8 => {
                ref_len += op_len;
                trailing_clips = 0;
                saw_ref_op = true;
            }
            4 | 5 if saw_ref_op => {
                trailing_clips += op_len;
            }
            _ => {}
        }
    }

    pos + ref_len + trailing_clips - 1
}

/// Compute unclipped 5' position directly from raw BAM bytes (zero allocation).
///
/// This is the primary entry point for raw-byte callers, replacing the pattern:
/// `let cigar = get_cigar_ops(bam); unclipped_5prime_1based(pos, reverse, unmapped, &cigar)`
#[inline]
#[must_use]
pub fn unclipped_5prime_from_raw_bam(bam: &[u8]) -> i32 {
    let flg = flags(bam);
    let unmapped = (flg & flags::UNMAPPED) != 0;
    if unmapped {
        return 0;
    }

    let n_cigar_op = u16::from_le_bytes([bam[12], bam[13]]) as usize;
    if n_cigar_op == 0 {
        return i32::MAX;
    }

    let l_read_name = bam[8] as usize;
    let cigar_start = 32 + l_read_name;
    let cigar_end = cigar_start + n_cigar_op * 4;
    if cigar_end > bam.len() {
        return i32::MAX;
    }

    let reverse = (flg & flags::REVERSE) != 0;
    let pos_1based = pos(bam) + 1;

    if reverse {
        unclipped_end_from_raw_cigar(pos_1based, bam, cigar_start, n_cigar_op)
    } else {
        unclipped_start_from_raw_cigar(pos_1based, bam, cigar_start, n_cigar_op)
    }
}

/// Compute reference length directly from raw CIGAR bytes (zero allocation).
#[inline]
#[must_use]
pub fn reference_length_from_raw_bam(bam: &[u8]) -> i32 {
    let n_cigar_op = u16::from_le_bytes([bam[12], bam[13]]) as usize;
    if n_cigar_op == 0 {
        return 0;
    }
    let l_read_name = bam[8] as usize;
    let cigar_start = 32 + l_read_name;
    let cigar_end = cigar_start + n_cigar_op * 4;
    if cigar_end > bam.len() {
        return 0;
    }

    let mut ref_len = 0i32;
    for i in 0..n_cigar_op {
        let op = cigar_op_at(bam, cigar_start + i * 4);
        let op_type = op & 0xF;
        if matches!(op_type, 0 | 2 | 3 | 7 | 8) {
            ref_len += (op >> 4) as i32;
        }
    }
    ref_len
}

/// Calculate reference-consuming length from CIGAR operations.
///
/// This is the sum of M/D/N/=/X operations, which represents how many
/// reference bases the alignment spans. Used for BAM index generation.
#[inline]
#[must_use]
pub fn reference_length_from_cigar(cigar_ops: &[u32]) -> i32 {
    let mut ref_len = 0i32;

    for &op in cigar_ops {
        let op_len = (op >> 4) as i32;
        let op_type = op & 0xF;

        // M (0), D (2), N (3), = (7), X (8) consume reference bases
        if matches!(op_type, 0 | 2 | 3 | 7 | 8) {
            ref_len += op_len;
        }
    }

    ref_len
}

/// Calculate unclipped 5' coordinate for this read.
///
/// For forward strand: `unclipped_start` (leftmost including clips)
/// For reverse strand: `unclipped_end` (rightmost including clips)
#[inline]
#[must_use]
pub fn unclipped_5prime(pos: i32, reverse: bool, cigar_ops: &[u32]) -> i32 {
    if reverse {
        unclipped_end_from_cigar(pos, cigar_ops)
    } else {
        unclipped_start_from_cigar(pos, cigar_ops)
    }
}

/// Calculate mate's unclipped 5' coordinate from MC tag.
///
/// For forward strand mate: `unclipped_start`
/// For reverse strand mate: `unclipped_end`
#[inline]
#[must_use]
pub fn mate_unclipped_5prime(mate_pos: i32, mate_reverse: bool, mc_cigar: &str) -> i32 {
    if mate_reverse {
        unclipped_other_end(mate_pos, mc_cigar)
    } else {
        unclipped_other_start(mate_pos, mc_cigar)
    }
}

// ============================================================================
// Consensus-oriented read primitives
// ============================================================================

/// Decode table: 4-bit BAM base code → ASCII base.
const BASE_DECODE: [u8; 16] = *b"=ACMGRSVTWYHKDBN";

/// Set a single base at a position in packed 4-bit BAM sequence data.
///
/// Converts the given ASCII base to a 4-bit BAM code and writes it in place.
/// Unknown bases are encoded as N (0xF).
#[inline]
pub fn set_base(bam: &mut [u8], seq_off: usize, position: usize, base: u8) {
    let encoded = SEQ_CODES[base as usize];
    let byte_idx = seq_off + position / 2;
    if position.is_multiple_of(2) {
        bam[byte_idx] = (encoded << 4) | (bam[byte_idx] & 0x0F);
    } else {
        bam[byte_idx] = (bam[byte_idx] & 0xF0) | encoded;
    }
}

/// Bulk-extract the full sequence from a BAM record as ASCII bases.
///
/// Decodes the packed 4-bit sequence data into a `Vec<u8>` of ASCII bases.
#[must_use]
pub fn extract_sequence(bam: &[u8]) -> Vec<u8> {
    let l = l_seq(bam) as usize;
    let off = seq_offset(bam);
    let mut bases = Vec::with_capacity(l);
    for i in 0..l {
        let code = get_base(bam, off, i);
        bases.push(BASE_DECODE[code as usize]);
    }
    bases
}

/// Zero-copy access to quality scores in a BAM record.
///
/// Returns a slice of the raw Phred quality scores (not Phred+33).
#[inline]
#[must_use]
pub fn quality_scores_slice(bam: &[u8]) -> &[u8] {
    let l = l_seq(bam) as usize;
    let off = qual_offset(bam);
    &bam[off..off + l]
}

/// Mutable zero-copy access to quality scores in a BAM record.
#[inline]
pub fn quality_scores_slice_mut(bam: &mut [u8]) -> &mut [u8] {
    let l = l_seq(bam) as usize;
    let off = qual_offset(bam);
    &mut bam[off..off + l]
}

/// Compute 1-based alignment end position from raw BAM bytes.
///
/// Returns `pos + ref_len` (1-based inclusive end, matching noodles convention).
/// Returns `None` if the record is unmapped or has no CIGAR.
#[inline]
#[must_use]
pub fn alignment_end_from_raw(bam: &[u8]) -> Option<usize> {
    let p = pos(bam);
    if p < 0 {
        return None;
    }
    let ref_len = reference_length_from_raw_bam(bam);
    if ref_len == 0 {
        return None;
    }
    // pos is 0-based, convert to 1-based and add ref_len - 1 for inclusive end
    Some((p + ref_len) as usize)
}

/// Compute 1-based alignment start position from raw BAM bytes.
///
/// Returns `None` if `pos < 0` (unmapped).
#[inline]
#[must_use]
pub fn alignment_start_from_raw(bam: &[u8]) -> Option<usize> {
    let p = pos(bam);
    if p < 0 { None } else { Some((p + 1) as usize) }
}

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

/// Extract template length (tlen) from a BAM record.
#[inline]
#[must_use]
pub fn template_length(bam: &[u8]) -> i32 {
    i32::from_le_bytes([bam[28], bam[29], bam[30], bam[31]])
}

/// Check if a single read is part of an FR (forward-reverse) pair using raw BAM bytes.
///
/// This is the raw-byte equivalent of `record_utils::is_fr_pair_from_tags`.
/// Returns `true` if the read is paired, both read and mate are mapped,
/// on the same reference, and in FR orientation (positive strand 5' < negative strand 5').
#[must_use]
pub fn is_fr_pair_raw(bam: &[u8]) -> bool {
    let flg = flags(bam);

    // Must be paired
    if flg & flags::PAIRED == 0 {
        return false;
    }

    // Both read and mate must be mapped
    if flg & flags::UNMAPPED != 0 || flg & flags::MATE_UNMAPPED != 0 {
        return false;
    }

    // Must be on the same reference
    let this_ref_id = ref_id(bam);
    let m_ref_id = mate_ref_id(bam);
    if this_ref_id != m_ref_id {
        return false;
    }

    // Must be on opposite strands for FR or RF
    let is_reverse = flg & flags::REVERSE != 0;
    let mate_is_reverse = flg & flags::MATE_REVERSE != 0;
    if is_reverse == mate_is_reverse {
        return false;
    }

    // Determine if FR or RF using htsjdk's logic:
    // positiveStrandFivePrimePos = readIsOnReverseStrand ? mateStart : alignmentStart
    // negativeStrandFivePrimePos = readIsOnReverseStrand ? alignmentEnd : alignmentStart + insertSize
    let alignment_start = (pos(bam) + 1) as usize; // 1-based
    let m_start = (mate_pos(bam) + 1) as usize; // 1-based
    let insert_size = template_length(bam);

    let (positive_five_prime, negative_five_prime) = if is_reverse {
        // This read is on reverse strand, mate is on positive strand
        let ref_len = reference_length_from_raw_bam(bam) as usize;
        let end = alignment_start + ref_len.saturating_sub(1);
        (m_start as i32, end as i32)
    } else {
        // This read is on positive strand, mate is on reverse strand
        (alignment_start as i32, alignment_start as i32 + insert_size)
    };

    // FR if positive strand 5' < negative strand 5'
    positive_five_prime < negative_five_prime
}

/// Virtual CIGAR clipping on raw u32 CIGAR ops.
///
/// Computes clipped CIGAR ops and the number of reference bases consumed by clipping,
/// without modifying any record. This is the raw-byte equivalent of
/// `SamRecordClipper::clip_start_of_read` / `clip_end_of_read` with Hard clipping mode.
///
/// Like the Clipper, this first accounts for existing H+S clips at the relevant end.
/// If `clip_amount` <= existing clips, only soft clips are upgraded to hard clips
/// (no alignment change). Otherwise, `clip_amount - existing_clips` bases are clipped
/// from the alignment.
///
/// Returns `(new_cigar_ops, ref_bases_consumed)`.
/// `ref_bases_consumed` is used to adjust `alignment_start` for start-clipping.
#[must_use]
pub fn clip_cigar_ops_raw(
    cigar_ops: &[u32],
    clip_amount: usize,
    from_start: bool,
) -> (Vec<u32>, usize) {
    if clip_amount == 0 || cigar_ops.is_empty() {
        return (cigar_ops.to_vec(), 0);
    }

    // Helper to encode a CIGAR op as raw u32
    let encode_op = |op_type: u32, len: usize| -> u32 { ((len as u32) << 4) | op_type };

    // Count existing H+S clips at the relevant end (matching clip_*_of_read)
    let existing_clip: usize = if from_start {
        cigar_ops
            .iter()
            .take_while(|&&op| matches!(op & 0xF, 4 | 5))
            .map(|&op| (op >> 4) as usize)
            .sum()
    } else {
        cigar_ops
            .iter()
            .rev()
            .take_while(|&&op| matches!(op & 0xF, 4 | 5))
            .map(|&op| (op >> 4) as usize)
            .sum()
    };

    if clip_amount <= existing_clip {
        // Just upgrade soft clips to hard clips (no alignment change)
        upgrade_clipping_raw(cigar_ops, clip_amount, from_start, encode_op)
    } else {
        // Clip into alignment: clip (clip_amount - existing_clip) additional bases
        let alignment_clip = clip_amount - existing_clip;
        if from_start {
            clip_cigar_start_raw(cigar_ops, alignment_clip, encode_op)
        } else {
            clip_cigar_end_raw(cigar_ops, alignment_clip, encode_op)
        }
    }
}

/// Upgrade existing soft clips to hard clips without changing alignment.
///
/// Matches `SamRecordClipper::upgrade_clipping` in Hard mode.
/// Used when `clip_amount` <= existing H+S clips.
fn upgrade_clipping_raw(
    cigar_ops: &[u32],
    clip_amount: usize,
    from_start: bool,
    encode_op: impl Fn(u32, usize) -> u32,
) -> (Vec<u32>, usize) {
    if from_start {
        let mut existing_hard = 0usize;
        let mut existing_soft = 0usize;
        let mut skip_count = 0usize;

        for &op in cigar_ops {
            if op & 0xF == 5 {
                existing_hard += (op >> 4) as usize;
                skip_count += 1;
            } else {
                break;
            }
        }
        for &op in &cigar_ops[skip_count..] {
            if op & 0xF == 4 {
                existing_soft += (op >> 4) as usize;
                skip_count += 1;
            } else {
                break;
            }
        }

        let length_to_upgrade = existing_soft.min(clip_amount.saturating_sub(existing_hard));
        let new_hard = existing_hard + length_to_upgrade;
        let remaining_soft = existing_soft - length_to_upgrade;

        let mut result = Vec::new();
        result.push(encode_op(5, new_hard));
        if remaining_soft > 0 {
            result.push(encode_op(4, remaining_soft));
        }
        result.extend_from_slice(&cigar_ops[skip_count..]);

        (result, 0)
    } else {
        let mut existing_hard = 0usize;
        let mut existing_soft = 0usize;
        let mut skip_count = 0usize;

        for &op in cigar_ops.iter().rev() {
            if op & 0xF == 5 {
                existing_hard += (op >> 4) as usize;
                skip_count += 1;
            } else {
                break;
            }
        }
        let end_idx = cigar_ops.len() - skip_count;
        for &op in cigar_ops[..end_idx].iter().rev() {
            if op & 0xF == 4 {
                existing_soft += (op >> 4) as usize;
                skip_count += 1;
            } else {
                break;
            }
        }

        let length_to_upgrade = existing_soft.min(clip_amount.saturating_sub(existing_hard));
        let new_hard = existing_hard + length_to_upgrade;
        let remaining_soft = existing_soft - length_to_upgrade;

        let end_content = cigar_ops.len() - skip_count;
        let mut result = cigar_ops[..end_content].to_vec();
        if remaining_soft > 0 {
            result.push(encode_op(4, remaining_soft));
        }
        result.push(encode_op(5, new_hard));

        (result, 0)
    }
}

/// Clip from the start of alignment (hard clip mode).
fn clip_cigar_start_raw(
    cigar_ops: &[u32],
    clip_amount: usize,
    encode_op: impl Fn(u32, usize) -> u32,
) -> (Vec<u32>, usize) {
    // Extract existing hard and soft clips from the start
    let mut existing_hard_clip = 0usize;
    let mut existing_soft_clip = 0usize;
    let mut skip_count = 0usize;

    for &op in cigar_ops {
        let op_type = op & 0xF;
        let op_len = (op >> 4) as usize;
        if op_type == 5 {
            // H
            existing_hard_clip += op_len;
            skip_count += 1;
        } else {
            break;
        }
    }
    for &op in &cigar_ops[skip_count..] {
        let op_type = op & 0xF;
        let op_len = (op >> 4) as usize;
        if op_type == 4 {
            // S
            existing_soft_clip += op_len;
            skip_count += 1;
        } else {
            break;
        }
    }

    let post_clip_ops = &cigar_ops[skip_count..];

    let mut read_bases_clipped = 0usize;
    let mut ref_bases_clipped = 0usize;
    let mut new_ops: Vec<u32> = Vec::new();
    let mut idx = 0;

    // Clip operations from the start
    while idx < post_clip_ops.len() {
        let op = post_clip_ops[idx];
        let op_type = op & 0xF;
        let op_len = (op >> 4) as usize;

        // Check if we've clipped enough but need to skip trailing deletions
        if read_bases_clipped == clip_amount && new_ops.is_empty() && op_type == 2 {
            // Deletion: skip it (ref-consuming only)
            ref_bases_clipped += op_len;
            idx += 1;
            continue;
        }

        if read_bases_clipped >= clip_amount {
            break;
        }

        // Note: S (soft clip) is NOT read-consuming here, matching the Clipper.
        // After stripping leading H+S, remaining S ops (e.g. trailing) should not
        // count toward read bases clipped.
        let consumes_read = matches!(op_type, 0 | 1 | 7 | 8); // M, I, =, X
        let consumes_ref = matches!(op_type, 0 | 2 | 3 | 7 | 8); // M, D, N, =, X

        if consumes_read && op_len > (clip_amount - read_bases_clipped) {
            if op_type == 1 {
                // Insertion: consume entire at clip boundary
                read_bases_clipped += op_len;
            } else {
                // Split the operation
                let remaining_clip = clip_amount - read_bases_clipped;
                let remaining_length = op_len - remaining_clip;
                read_bases_clipped += remaining_clip;
                if consumes_ref {
                    ref_bases_clipped += remaining_clip;
                }
                new_ops.push(encode_op(op_type, remaining_length));
            }
        } else {
            if consumes_read {
                read_bases_clipped += op_len;
            }
            if consumes_ref {
                ref_bases_clipped += op_len;
            }
        }

        idx += 1;
    }

    // Add remaining operations
    new_ops.extend_from_slice(&post_clip_ops[idx..]);

    // Hard clip mode: convert all existing soft clips to hard clips
    let added_hard_clip = existing_soft_clip + read_bases_clipped;
    let total_hard_clip = existing_hard_clip + added_hard_clip;
    let mut result = Vec::with_capacity(1 + new_ops.len());
    result.push(encode_op(5, total_hard_clip)); // H
    result.extend(new_ops);

    (result, ref_bases_clipped)
}

/// Clip from the end of alignment (hard clip mode).
fn clip_cigar_end_raw(
    cigar_ops: &[u32],
    clip_amount: usize,
    encode_op: impl Fn(u32, usize) -> u32,
) -> (Vec<u32>, usize) {
    // Extract existing hard and soft clips from the end
    let mut existing_hard_clip = 0usize;
    let mut existing_soft_clip = 0usize;
    let mut skip_count = 0usize;

    for &op in cigar_ops.iter().rev() {
        let op_type = op & 0xF;
        let op_len = (op >> 4) as usize;
        if op_type == 5 {
            // H
            existing_hard_clip += op_len;
            skip_count += 1;
        } else {
            break;
        }
    }
    let end_idx = cigar_ops.len() - skip_count;
    for &op in cigar_ops[..end_idx].iter().rev() {
        let op_type = op & 0xF;
        let op_len = (op >> 4) as usize;
        if op_type == 4 {
            // S
            existing_soft_clip += op_len;
            skip_count += 1;
        } else {
            break;
        }
    }

    let post_clip_end = cigar_ops.len() - skip_count;
    let post_clip_ops = &cigar_ops[..post_clip_end];

    let mut read_bases_clipped = 0usize;
    let mut new_ops: Vec<u32> = Vec::new();
    let mut idx = post_clip_ops.len();

    // Clip operations from the end (working backwards)
    while idx > 0 {
        let op = post_clip_ops[idx - 1];
        let op_type = op & 0xF;
        let op_len = (op >> 4) as usize;

        // Check if we've clipped enough but need to skip adjacent deletions
        if read_bases_clipped == clip_amount && new_ops.is_empty() && op_type == 2 {
            // Deletion: skip it
            idx -= 1;
            continue;
        }

        if read_bases_clipped >= clip_amount {
            break;
        }

        // Note: S (soft clip) is NOT read-consuming here, matching the Clipper.
        let consumes_read = matches!(op_type, 0 | 1 | 7 | 8); // M, I, =, X

        if consumes_read && op_len > (clip_amount - read_bases_clipped) {
            if op_type == 1 {
                // Insertion: consume entire at clip boundary
                read_bases_clipped += op_len;
            } else {
                // Split the operation
                let remaining_clip = clip_amount - read_bases_clipped;
                let remaining_length = op_len - remaining_clip;
                read_bases_clipped += remaining_clip;
                new_ops.push(encode_op(op_type, remaining_length));
            }
        } else if consumes_read {
            read_bases_clipped += op_len;
        }

        idx -= 1;
    }

    // Add remaining operations (new_ops is in reverse, remaining are forward)
    let remaining: Vec<u32> = post_clip_ops[..idx].to_vec();
    let mut result = remaining;
    // new_ops collected in reverse order, need to reverse
    new_ops.reverse();
    result.extend(new_ops);

    // Hard clip mode: convert all existing soft clips to hard clips
    let added_hard_clip = existing_soft_clip + read_bases_clipped;
    let total_hard_clip = existing_hard_clip + added_hard_clip;
    result.push(encode_op(5, total_hard_clip)); // H

    (result, 0) // ref_bases_consumed is 0 for end clipping (no position adjustment)
}

/// Returns the query position (1-based) at a given reference position, from raw CIGAR ops.
///
/// This is the raw-byte equivalent of `CodecConsensusCaller::read_pos_at_ref_pos`.
///
/// # Arguments
/// * `cigar_ops` - Raw u32 CIGAR operations
/// * `alignment_start` - 1-based alignment start position
/// * `ref_pos` - 1-based reference position to query
/// * `return_last_base_if_deleted` - If true, returns last query position before deletion
///
/// Returns `None` if the position falls outside the alignment or in a deletion
/// (when `return_last_base_if_deleted` is false).
#[must_use]
pub fn read_pos_at_ref_pos_raw(
    cigar_ops: &[u32],
    alignment_start: usize,
    ref_pos: usize,
    return_last_base_if_deleted: bool,
) -> Option<usize> {
    if ref_pos < alignment_start {
        return None;
    }

    let mut ref_offset = 0usize;
    let mut query_offset = 0usize;

    for &op in cigar_ops {
        let op_type = op & 0xF;
        let op_len = (op >> 4) as usize;

        let consumes_ref = matches!(op_type, 0 | 2 | 3 | 7 | 8); // M, D, N, =, X
        let consumes_query = matches!(op_type, 0 | 1 | 4 | 7 | 8); // M, I, S, =, X

        let op_ref_start = alignment_start + ref_offset;

        if consumes_ref {
            let op_ref_end = op_ref_start + op_len - 1;

            if ref_pos >= op_ref_start && ref_pos <= op_ref_end {
                if consumes_query {
                    // M, =, X: we have a base at this position
                    let offset_in_op = ref_pos - op_ref_start;
                    return Some(query_offset + offset_in_op + 1); // 1-based
                }
                // D, N: position falls in a deletion
                if return_last_base_if_deleted {
                    return Some(if query_offset > 0 { query_offset } else { 1 });
                }
                return None;
            }
        }

        if consumes_ref {
            ref_offset += op_len;
        }
        if consumes_query {
            query_offset += op_len;
        }
    }

    None
}

/// Compute the query-consuming length of CIGAR operations (the "read length").
///
/// This is the sum of M/I/S/=/X operations.
#[inline]
#[must_use]
pub fn query_length_from_cigar(cigar_ops: &[u32]) -> usize {
    let mut len = 0usize;
    for &op in cigar_ops {
        let op_type = op & 0xF;
        let op_len = (op >> 4) as usize;
        if matches!(op_type, 0 | 1 | 4 | 7 | 8) {
            len += op_len;
        }
    }
    len
}

/// Compute mate overlap clip count directly from raw BAM bytes.
///
/// This is the raw-byte equivalent of `record_utils::num_bases_extending_past_mate`.
/// Uses MC tag, mate position, and flags from raw bytes to determine how many
/// bases extend past the mate's unclipped boundary.
///
/// Returns 0 if not an FR pair or if required information is missing.
#[must_use]
pub fn num_bases_extending_past_mate_raw(bam: &[u8]) -> usize {
    let flg = flags(bam);

    // Only applies to paired, mapped reads with mapped mates
    if flg & flags::PAIRED == 0 || flg & flags::UNMAPPED != 0 || flg & flags::MATE_UNMAPPED != 0 {
        return 0;
    }

    // Check FR pair: read and mate on same reference, opposite strands
    let is_reverse = flg & flags::REVERSE != 0;
    let mate_is_reverse = flg & flags::MATE_REVERSE != 0;

    // FR pair requires opposite strand orientations
    if is_reverse == mate_is_reverse {
        return 0;
    }

    let this_ref_id = ref_id(bam);
    let m_ref_id = mate_ref_id(bam);
    if this_ref_id != m_ref_id {
        return 0;
    }

    // Need MC tag for mate CIGAR information
    let aux = aux_data_slice(bam);
    let Some(mc_bytes) = find_string_tag(aux, b"MC") else {
        return 0;
    };
    let Ok(mc_cigar) = std::str::from_utf8(mc_bytes) else {
        return 0;
    };

    let this_pos_0based = pos(bam);
    let m_pos_0based = mate_pos(bam);
    // Convert to 1-based for coordinate calculations
    let this_pos = this_pos_0based + 1;
    let m_pos = m_pos_0based + 1;

    // Calculate read length from CIGAR (query-consuming ops)
    let cigar_ops = get_cigar_ops(bam);
    let read_length: usize = cigar_ops
        .iter()
        .map(|&op| {
            let op_type = op & 0xF;
            let op_len = (op >> 4) as usize;
            // M(0), I(1), S(4), =(7), X(8) consume query
            if matches!(op_type, 0 | 1 | 4 | 7 | 8) { op_len } else { 0 }
        })
        .sum();

    if is_reverse {
        // Negative strand: check if read extends before mate's unclipped start
        let mate_us = unclipped_other_start(m_pos, mc_cigar) as usize;

        if (this_pos as usize) <= mate_us {
            compute_bases_before_ref_pos(bam, &cigar_ops, this_pos, mate_us)
        } else {
            // Only clip excess soft-clipped bases at the start
            let leading_sc = leading_soft_clip_from_ops(&cigar_ops);
            let gap = this_pos as usize - mate_us;
            leading_sc.saturating_sub(gap)
        }
    } else {
        // Positive strand: check if read extends past mate's unclipped end
        let ref_len = reference_length_from_cigar(&cigar_ops);
        let alignment_end = (this_pos + ref_len - 1) as usize;
        let mate_ue = unclipped_other_end(m_pos, mc_cigar) as usize;

        if alignment_end >= mate_ue {
            // Compute read position at mate's unclipped end
            // Simplified: use reference position mapping
            let bases_past = compute_bases_past_ref_pos(bam, &cigar_ops, this_pos, mate_ue);
            read_length.saturating_sub(bases_past)
        } else {
            // Only clip excess soft-clipped bases
            let trailing_sc = trailing_soft_clip_from_ops(&cigar_ops);
            let gap = mate_ue - alignment_end;
            trailing_sc.saturating_sub(gap)
        }
    }
}

/// Compute number of read bases at or past a reference position (for positive strand).
///
/// Returns the 1-based read position at the given reference position,
/// or 0 if the position falls in a deletion or outside the alignment.
fn compute_bases_past_ref_pos(
    _bam: &[u8],
    cigar_ops: &[u32],
    alignment_start_1based: i32,
    target_ref_pos: usize,
) -> usize {
    let mut ref_pos = alignment_start_1based as usize;
    let mut read_pos: usize = 0;

    for &op in cigar_ops {
        let op_type = op & 0xF;
        let op_len = (op >> 4) as usize;

        match op_type {
            0 | 7 | 8 => {
                // M, =, X: consume both query and reference
                for _ in 0..op_len {
                    read_pos += 1;
                    if ref_pos == target_ref_pos {
                        return read_pos;
                    }
                    ref_pos += 1;
                }
            }
            1 => {
                // I: consume query only
                read_pos += op_len;
            }
            4 => {
                // S: consume query only
                read_pos += op_len;
            }
            2 | 3 => {
                // D, N: consume reference only
                for _ in 0..op_len {
                    if ref_pos == target_ref_pos {
                        return 0; // Position in deletion
                    }
                    ref_pos += 1;
                }
            }
            _ => {}
        }
    }

    0
}

/// Compute number of read bases before a reference position (for negative strand).
fn compute_bases_before_ref_pos(
    _bam: &[u8],
    cigar_ops: &[u32],
    alignment_start_1based: i32,
    target_ref_pos: usize,
) -> usize {
    let mut ref_pos = alignment_start_1based as usize;
    let mut read_pos: usize = 0;

    for &op in cigar_ops {
        let op_type = op & 0xF;
        let op_len = (op >> 4) as usize;

        match op_type {
            0 | 7 | 8 => {
                // M, =, X: consume both query and reference
                for _ in 0..op_len {
                    read_pos += 1;
                    if ref_pos == target_ref_pos {
                        return read_pos.saturating_sub(1);
                    }
                    ref_pos += 1;
                }
            }
            1 => {
                // I: consume query only
                read_pos += op_len;
            }
            4 => {
                // S: consume query only
                read_pos += op_len;
            }
            2 | 3 => {
                // D, N: consume reference only
                for _ in 0..op_len {
                    if ref_pos == target_ref_pos {
                        return 0;
                    }
                    ref_pos += 1;
                }
            }
            _ => {}
        }
    }

    0
}

/// Count trailing soft clips from CIGAR ops.
fn trailing_soft_clip_from_ops(cigar_ops: &[u32]) -> usize {
    let mut trailing = 0usize;
    for &op in cigar_ops.iter().rev() {
        let op_type = op & 0xF;
        let op_len = (op >> 4) as usize;
        match op_type {
            4 => trailing += op_len, // S
            5 => {}                  // H - skip
            _ => break,
        }
    }
    trailing
}

/// Count leading soft clips from CIGAR ops.
fn leading_soft_clip_from_ops(cigar_ops: &[u32]) -> usize {
    let mut leading = 0usize;
    for &op in cigar_ops {
        let op_type = op & 0xF;
        let op_len = (op >> 4) as usize;
        match op_type {
            4 => leading += op_len, // S
            5 => {}                 // H - skip
            _ => break,
        }
    }
    leading
}

/// Decode raw BAM byte records to noodles `RecordBuf`.
///
/// This is a temporary bridge for consensus callers (duplex, codec) that still
/// use `RecordBuf` internally. It constructs a minimal BAM stream and uses
/// noodles reader to parse records.
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
        let block_size = raw.len() as u32;
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
    use rstest::rstest;

    #[test]
    fn test_extract_coordinate_fields_basic() {
        // Create a minimal BAM record
        let mut bam_bytes = vec![0u8; 48];
        // tid = 1
        bam_bytes[0..4].copy_from_slice(&1i32.to_le_bytes());
        // pos = 100
        bam_bytes[4..8].copy_from_slice(&100i32.to_le_bytes());
        // l_read_name = 5 (4 chars + null)
        bam_bytes[8] = 5;
        // flag = 0 (mapped, forward)
        bam_bytes[14..16].copy_from_slice(&0u16.to_le_bytes());
        // read name = "test"
        bam_bytes[32..36].copy_from_slice(b"test");

        let (tid, pos, reverse, name) = extract_coordinate_fields(&bam_bytes);
        assert_eq!(tid, 1);
        assert_eq!(pos, 100);
        assert!(!reverse);
        assert_eq!(name, b"test");
    }

    #[test]
    fn test_extract_coordinate_fields_unmapped() {
        let mut bam_bytes = vec![0u8; 48];
        bam_bytes[14..16].copy_from_slice(&flags::UNMAPPED.to_le_bytes());
        bam_bytes[8] = 1; // empty name

        let (tid, pos, reverse, _name) = extract_coordinate_fields(&bam_bytes);
        assert_eq!(tid, i32::MAX);
        assert_eq!(pos, i32::MAX);
        assert!(!reverse);
    }

    #[test]
    fn test_extract_coordinate_fields_reverse() {
        let mut bam_bytes = vec![0u8; 48];
        bam_bytes[14..16].copy_from_slice(&flags::REVERSE.to_le_bytes());
        bam_bytes[8] = 1;

        let (_tid, _pos, reverse, _name) = extract_coordinate_fields(&bam_bytes);
        assert!(reverse);
    }

    #[test]
    fn test_aux_data_offset() {
        // l_read_name=5, n_cigar_op=2, l_seq=50
        // offset = 32 + 5 + 2*4 + (50+1)/2 + 50 = 32 + 5 + 8 + 25 + 50 = 120
        assert_eq!(aux_data_offset(5, 2, 50), 120);
    }

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

    // ========================================================================
    // Helper: make_bam_bytes
    // ========================================================================

    /// Construct a raw BAM byte array for testing.
    ///
    /// IMPORTANT: `name` length + 1 (for null terminator) should be divisible
    /// by 4 to maintain alignment for CIGAR ops.  Use names like b"rea" (3+1=4)
    /// or b"readABC" (7+1=8).
    #[allow(clippy::too_many_arguments)]
    fn make_bam_bytes(
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
        let l_read_name = (name.len() + 1) as u8; // +1 for null terminator
        let n_cigar_op = cigar_ops.len() as u16;
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
        buf[16..20].copy_from_slice(&(seq_len as u32).to_le_bytes());
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

    // ========================================================================
    // compare_coordinate_raw tests
    // ========================================================================

    #[test]
    fn test_compare_coordinate_raw_same_records() {
        let rec = make_bam_bytes(1, 100, 0, b"rea", &[], 0, -1, -1, &[]);
        assert_eq!(compare_coordinate_raw(&rec, &rec), Ordering::Equal);
    }

    #[test]
    fn test_compare_coordinate_raw_different_tid() {
        let a = make_bam_bytes(0, 100, 0, b"rea", &[], 0, -1, -1, &[]);
        let b = make_bam_bytes(2, 100, 0, b"rea", &[], 0, -1, -1, &[]);
        assert_eq!(compare_coordinate_raw(&a, &b), Ordering::Less);
        assert_eq!(compare_coordinate_raw(&b, &a), Ordering::Greater);
    }

    #[test]
    fn test_compare_coordinate_raw_different_pos() {
        let a = make_bam_bytes(1, 50, 0, b"rea", &[], 0, -1, -1, &[]);
        let b = make_bam_bytes(1, 200, 0, b"rea", &[], 0, -1, -1, &[]);
        assert_eq!(compare_coordinate_raw(&a, &b), Ordering::Less);
        assert_eq!(compare_coordinate_raw(&b, &a), Ordering::Greater);
    }

    #[test]
    fn test_compare_coordinate_raw_reverse_strand() {
        // Forward strand (flag=0) should sort before reverse (flag=REVERSE)
        let fwd = make_bam_bytes(1, 100, 0, b"rea", &[], 0, -1, -1, &[]);
        let rev = make_bam_bytes(1, 100, flags::REVERSE, b"rea", &[], 0, -1, -1, &[]);
        assert_eq!(compare_coordinate_raw(&fwd, &rev), Ordering::Less);
        assert_eq!(compare_coordinate_raw(&rev, &fwd), Ordering::Greater);
    }

    #[test]
    fn test_compare_coordinate_raw_name_tiebreak() {
        let a = make_bam_bytes(1, 100, 0, b"aaa", &[], 0, -1, -1, &[]);
        let b = make_bam_bytes(1, 100, 0, b"zzz", &[], 0, -1, -1, &[]);
        assert_eq!(compare_coordinate_raw(&a, &b), Ordering::Less);
        assert_eq!(compare_coordinate_raw(&b, &a), Ordering::Greater);
    }

    #[test]
    fn test_compare_coordinate_raw_no_ref_sorts_last() {
        // tid=-1 (no reference) should sort after mapped records
        let mapped = make_bam_bytes(1, 100, 0, b"rea", &[], 0, -1, -1, &[]);
        let no_ref = make_bam_bytes(-1, -1, 0, b"rea", &[], 0, -1, -1, &[]);
        assert_eq!(compare_coordinate_raw(&mapped, &no_ref), Ordering::Less);
        assert_eq!(compare_coordinate_raw(&no_ref, &mapped), Ordering::Greater);
    }

    #[test]
    fn test_compare_coordinate_raw_both_no_ref() {
        // Both tid=-1, compare by name
        let a = make_bam_bytes(-1, -1, 0, b"aaa", &[], 0, -1, -1, &[]);
        let b = make_bam_bytes(-1, -1, 0, b"zzz", &[], 0, -1, -1, &[]);
        assert_eq!(compare_coordinate_raw(&a, &b), Ordering::Less);
        assert_eq!(compare_coordinate_raw(&b, &a), Ordering::Greater);
    }

    #[test]
    fn test_compare_coordinate_raw_one_no_ref_each_direction() {
        let mapped = make_bam_bytes(5, 999, 0, b"rea", &[], 0, -1, -1, &[]);
        let no_ref = make_bam_bytes(-1, 0, 0, b"rea", &[], 0, -1, -1, &[]);
        // Mapped < no_ref
        assert_eq!(compare_coordinate_raw(&mapped, &no_ref), Ordering::Less);
        // no_ref > mapped
        assert_eq!(compare_coordinate_raw(&no_ref, &mapped), Ordering::Greater);
    }

    // ========================================================================
    // compare_queryname_raw tests
    // ========================================================================

    #[test]
    fn test_compare_queryname_raw_name_ordering() {
        let a = make_bam_bytes(0, 0, 0, b"abc", &[], 0, -1, -1, &[]);
        let b = make_bam_bytes(0, 0, 0, b"xyz", &[], 0, -1, -1, &[]);
        assert_eq!(compare_queryname_raw(&a, &b), Ordering::Less);
        assert_eq!(compare_queryname_raw(&b, &a), Ordering::Greater);
    }

    #[test]
    fn test_compare_queryname_raw_same_name_flag_tiebreak() {
        // Same name, different flags -> flag decides
        let a = make_bam_bytes(0, 0, 0x40, b"rea", &[], 0, -1, -1, &[]); // first in pair
        let b = make_bam_bytes(0, 0, 0x80, b"rea", &[], 0, -1, -1, &[]); // second in pair
        assert_eq!(compare_queryname_raw(&a, &b), Ordering::Less);
        assert_eq!(compare_queryname_raw(&b, &a), Ordering::Greater);
    }

    #[test]
    fn test_compare_queryname_raw_equal() {
        let rec = make_bam_bytes(1, 100, 0x40, b"rea", &[], 0, -1, -1, &[]);
        assert_eq!(compare_queryname_raw(&rec, &rec), Ordering::Equal);
    }

    // ========================================================================
    // Unclipped position tests
    // ========================================================================

    #[test]
    fn test_unclipped_start_from_cigar_no_clips() {
        // 10M = (10 << 4) | 0
        let cigar = &[(10 << 4) | 0];
        assert_eq!(unclipped_start_from_cigar(100, cigar), 100);
    }

    #[test]
    fn test_unclipped_start_from_cigar_soft_clip() {
        // 5S10M: soft clip (op type 4), then match
        let cigar = &[(5 << 4) | 4, (10 << 4) | 0];
        assert_eq!(unclipped_start_from_cigar(100, cigar), 95);
    }

    #[test]
    fn test_unclipped_start_from_cigar_hard_clip() {
        // 3H10M: hard clip (op type 5), then match
        let cigar = &[(3 << 4) | 5, (10 << 4) | 0];
        assert_eq!(unclipped_start_from_cigar(100, cigar), 97);
    }

    #[test]
    fn test_unclipped_end_from_cigar_no_clips() {
        // 10M: end = pos + 10 - 1 = 109
        let cigar = &[(10 << 4) | 0];
        assert_eq!(unclipped_end_from_cigar(100, cigar), 109);
    }

    #[test]
    fn test_unclipped_end_from_cigar_trailing_clips() {
        // 10M5S3H: end = 100 + 10 + 5 + 3 - 1 = 117
        let cigar = &[(10 << 4) | 0, (5 << 4) | 4, (3 << 4) | 5];
        assert_eq!(unclipped_end_from_cigar(100, cigar), 117);
    }

    #[test]
    fn test_unclipped_5prime_forward_vs_reverse() {
        // 5S10M3S: forward uses start, reverse uses end
        let cigar = &[(5 << 4) | 4, (10 << 4) | 0, (3 << 4) | 4];
        // Forward: unclipped_start = 100 - 5 = 95
        assert_eq!(unclipped_5prime(100, false, cigar), 95);
        // Reverse: unclipped_end = 100 + 10 + 3 - 1 = 112
        assert_eq!(unclipped_5prime(100, true, cigar), 112);
    }

    // ========================================================================
    // CIGAR string parsing tests (parse_leading_clips, parse_ref_len_and_trailing_clips)
    // ========================================================================

    #[test]
    fn test_parse_leading_clips_soft() {
        // "5S10M" -> leading clips = 5
        assert_eq!(parse_leading_clips("5S10M"), 5);
    }

    #[test]
    fn test_parse_leading_clips_hard_and_soft() {
        // "3H5S10M" -> leading clips = 3 + 5 = 8
        assert_eq!(parse_leading_clips("3H5S10M"), 8);
    }

    #[test]
    fn test_parse_leading_clips_no_clips() {
        // "10M" -> leading clips = 0
        assert_eq!(parse_leading_clips("10M"), 0);
    }

    #[test]
    fn test_parse_ref_len_and_trailing_clips_basic() {
        // "10M5S" -> ref_len=10, trailing_clips=5
        assert_eq!(parse_ref_len_and_trailing_clips("10M5S"), (10, 5));
    }

    #[test]
    fn test_parse_ref_len_and_trailing_clips_complex() {
        // "5S10M2I3D5M3S2H"
        // ref consuming: 10M + 3D + 5M = 18
        // trailing clips: 3S + 2H = 5
        assert_eq!(parse_ref_len_and_trailing_clips("5S10M2I3D5M3S2H"), (18, 5));
    }

    // ========================================================================
    // reference_length_from_cigar tests
    // ========================================================================

    #[test]
    fn test_reference_length_simple_match() {
        // 50M
        let cigar = &[(50 << 4) | 0];
        assert_eq!(reference_length_from_cigar(cigar), 50);
    }

    #[test]
    fn test_reference_length_with_deletions() {
        // 10M3D5M2N8M
        // M (0), D (2), N (3) all consume reference
        // ref_len = 10 + 3 + 5 + 2 + 8 = 28
        let cigar = &[
            (10 << 4) | 0, // 10M
            (3 << 4) | 2,  // 3D
            (5 << 4) | 0,  // 5M
            (2 << 4) | 3,  // 2N
            (8 << 4) | 0,  // 8M
        ];
        assert_eq!(reference_length_from_cigar(cigar), 28);
    }

    #[test]
    fn test_reference_length_with_insertions() {
        // 10M5I10M: insertions don't consume reference
        // ref_len = 10 + 10 = 20
        let cigar = &[
            (10 << 4) | 0, // 10M
            (5 << 4) | 1,  // 5I
            (10 << 4) | 0, // 10M
        ];
        assert_eq!(reference_length_from_cigar(cigar), 20);
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

    // ========================================================================
    // tag_value_size tests
    // ========================================================================

    #[test]
    fn test_tag_value_size_fixed() {
        assert_eq!(tag_value_size(b'A', &[0]), Some(1));
        assert_eq!(tag_value_size(b'c', &[0]), Some(1));
        assert_eq!(tag_value_size(b'C', &[0]), Some(1));
        assert_eq!(tag_value_size(b's', &[0, 0]), Some(2));
        assert_eq!(tag_value_size(b'S', &[0, 0]), Some(2));
        assert_eq!(tag_value_size(b'i', &[0, 0, 0, 0]), Some(4));
        assert_eq!(tag_value_size(b'I', &[0, 0, 0, 0]), Some(4));
        assert_eq!(tag_value_size(b'f', &[0, 0, 0, 0]), Some(4));
    }

    #[test]
    fn test_tag_value_size_string() {
        // Z type: null-terminated string "hello\0" -> size = 6
        let data = b"hello\x00";
        assert_eq!(tag_value_size(b'Z', data), Some(6));

        // H type: hex-encoded string "ABCD\0" -> size = 5
        let data = b"ABCD\x00";
        assert_eq!(tag_value_size(b'H', data), Some(5));
    }

    #[test]
    fn test_tag_value_size_array() {
        // B type: array header is elem_type (1 byte) + count (4 bytes) + data
        // B:i:3 -> 3 int32 elements = 5 + 3*4 = 17
        let mut data = vec![b'i']; // elem_type
        data.extend_from_slice(&3u32.to_le_bytes()); // count = 3
        data.extend_from_slice(&[0; 12]); // 3 * 4 bytes of data
        assert_eq!(tag_value_size(b'B', &data), Some(17));
    }

    // ========================================================================
    // find_mi_tag edge cases
    // ========================================================================

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
    #[case::negative_signed_byte(b'c', &[(-1i8) as u8], None)]
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

    // ========================================================================
    // parse_mi_bytes tests (tested through find_mi_tag with Z-type values)
    // ========================================================================

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

    // ========================================================================
    // compare_template_coordinate_raw tests
    // ========================================================================

    #[test]
    fn test_compare_template_coordinate_raw_equal() {
        // Two identical unmapped, unpaired records -> should be Equal
        let rec = make_bam_bytes(
            -1,
            -1,
            flags::UNMAPPED, // unmapped, not paired
            b"rea",
            &[],
            0,
            -1,
            -1,
            &[],
        );
        assert_eq!(compare_template_coordinate_raw(&rec, &rec), Ordering::Equal);
    }

    #[test]
    fn test_compare_template_coordinate_raw_different_tid() {
        // Two mapped, unpaired records with different tids
        // 10M cigar for proper alignment
        let cigar = &[(10 << 4) | 0]; // 10M
        let a = make_bam_bytes(0, 100, 0, b"rea", cigar, 10, -1, -1, &[]);
        let b = make_bam_bytes(2, 100, 0, b"rea", cigar, 10, -1, -1, &[]);
        assert_eq!(compare_template_coordinate_raw(&a, &b), Ordering::Less);
        assert_eq!(compare_template_coordinate_raw(&b, &a), Ordering::Greater);
    }

    // ========================================================================
    // Read primitive tests
    // ========================================================================

    #[test]
    fn test_read_primitives() {
        let rec = make_bam_bytes(
            3,                              // tid
            200,                            // pos
            flags::PAIRED | flags::REVERSE, // flags
            b"rea",                         // name
            &[(10 << 4) | 0],               // 10M cigar
            10,                             // seq_len
            5,                              // mate_tid
            300,                            // mate_pos
            &[],
        );
        // Set mapq
        let mut rec = rec;
        rec[9] = 42;

        assert_eq!(ref_id(&rec), 3);
        assert_eq!(pos(&rec), 200);
        assert_eq!(flags(&rec), flags::PAIRED | flags::REVERSE);
        assert_eq!(mapq(&rec), 42);
        assert_eq!(l_read_name(&rec), 4); // "rea" + NUL
        assert_eq!(n_cigar_op(&rec), 1);
        assert_eq!(l_seq(&rec), 10);
        assert_eq!(mate_ref_id(&rec), 5);
        assert_eq!(mate_pos(&rec), 300);
        assert_eq!(read_name(&rec), b"rea");
    }

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

    #[rstest]
    #[case::signed_int(b'i', &42i32.to_le_bytes(), Some(42))]
    #[case::signed_byte(b'c', &[(-5i8) as u8], Some(-5))]
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
    fn test_find_uint8_tag() {
        let aux = [b'M', b'Q', b'C', 30];
        assert_eq!(find_uint8_tag(&aux, b"MQ"), Some(30));
    }

    #[test]
    fn test_find_float_tag() {
        let val: f32 = 1.234;
        let mut aux = vec![b'X', b'F', b'f'];
        aux.extend_from_slice(&val.to_le_bytes());
        let result = find_float_tag(&aux, b"XF").unwrap();
        assert!((result - val).abs() < 0.001);
    }

    // ========================================================================
    // Write primitive tests
    // ========================================================================

    #[test]
    fn test_set_flags() {
        let mut rec = make_bam_bytes(0, 0, 0, b"rea", &[], 0, -1, -1, &[]);
        assert_eq!(flags(&rec), 0);
        set_flags(&mut rec, flags::DUPLICATE | flags::PAIRED);
        assert_eq!(flags(&rec), flags::DUPLICATE | flags::PAIRED);
    }

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

    // ========================================================================
    // Sequence and quality tests
    // ========================================================================

    #[test]
    fn test_get_base_and_mask() {
        // Build a record with seq_len=4, fill seq bytes manually
        let mut rec = make_bam_bytes(0, 0, 0, b"rea", &[], 4, -1, -1, &[]);
        let so = seq_offset(&rec);
        // Packed: A=1, C=2 -> byte = 0x12; G=4, T=8 -> byte = 0x48
        rec[so] = 0x12; // bases 0,1 = A,C
        rec[so + 1] = 0x48; // bases 2,3 = G,T

        assert_eq!(get_base(&rec, so, 0), 1); // A
        assert_eq!(get_base(&rec, so, 1), 2); // C
        assert_eq!(get_base(&rec, so, 2), 4); // G
        assert_eq!(get_base(&rec, so, 3), 8); // T

        // Mask base 1 (C -> N)
        mask_base(&mut rec, so, 1);
        assert_eq!(get_base(&rec, so, 1), 0xF); // N
        assert_eq!(get_base(&rec, so, 0), 1); // A unchanged

        // Mask base 2 (G -> N, high nibble)
        mask_base(&mut rec, so, 2);
        assert_eq!(get_base(&rec, so, 2), 0xF); // N
        assert_eq!(get_base(&rec, so, 3), 8); // T unchanged
    }

    #[test]
    fn test_get_qual_and_set_qual() {
        let mut rec = make_bam_bytes(0, 0, 0, b"rea", &[], 4, -1, -1, &[]);
        let qo = qual_offset(&rec);
        rec[qo] = 30;
        rec[qo + 1] = 40;
        assert_eq!(get_qual(&rec, qo, 0), 30);
        assert_eq!(get_qual(&rec, qo, 1), 40);

        set_qual(&mut rec, qo, 0, 10);
        assert_eq!(get_qual(&rec, qo, 0), 10);
        assert_eq!(get_qual(&rec, qo, 1), 40); // unchanged
    }

    // ========================================================================
    // Coordinate helper tests (1-based)
    // ========================================================================

    #[test]
    fn test_unclipped_5prime_1based_unmapped() {
        assert_eq!(unclipped_5prime_1based(100, false, true, &[(10 << 4) | 0]), 0);
    }

    #[test]
    fn test_unclipped_5prime_1based_no_cigar() {
        assert_eq!(unclipped_5prime_1based(100, false, false, &[]), i32::MAX);
    }

    #[test]
    fn test_unclipped_5prime_1based_forward() {
        // 5S10M: forward 5' = pos+1 - 5 = 96
        let cigar = &[(5 << 4) | 4, (10 << 4) | 0];
        assert_eq!(unclipped_5prime_1based(100, false, false, cigar), 96);
    }

    #[test]
    fn test_unclipped_5prime_1based_reverse() {
        // 10M5S: reverse 5' = pos+1 + 10 + 5 - 1 = 115
        let cigar = &[(10 << 4) | 0, (5 << 4) | 4];
        assert_eq!(unclipped_5prime_1based(100, true, false, cigar), 115);
    }

    #[test]
    fn test_mate_unclipped_5prime_1based_forward() {
        // MC=5S10M: forward 5' = pos+1 - 5 = 96
        assert_eq!(mate_unclipped_5prime_1based(100, false, "5S10M"), 96);
    }

    #[test]
    fn test_mate_unclipped_5prime_1based_reverse() {
        // MC=10M5S: reverse 5' = pos+1 + 10 + 5 - 1 = 115
        assert_eq!(mate_unclipped_5prime_1based(100, true, "10M5S"), 115);
    }

    // ========================================================================
    // Tests for MQ tag type mismatch bug (raw-byte pipeline issue)
    // find_uint8_tag only matches type 'C'; callers should use find_int_tag for signed types
    // ========================================================================

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

    // ========================================================================
    // Tests for bounds check panics (raw-byte pipeline issue)
    // ========================================================================

    #[rstest]
    #[case::truncated_signed_byte(b"MIc" as &[u8])]
    #[case::truncated_signed_short(&[b'M', b'I', b's', 42] as &[u8])]
    #[case::truncated_signed_int(&[b'M', b'I', b'i', 42, 0] as &[u8])]
    fn test_find_mi_tag_truncated_returns_none(#[case] aux: &[u8]) {
        // Truncated tag data returns None instead of panicking
        assert_eq!(find_mi_tag(aux), None);
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

    // ========================================================================
    // unclipped_5prime_from_raw_bam tests
    // ========================================================================

    #[test]
    fn test_unclipped_5prime_from_raw_bam_unmapped() {
        let rec =
            make_bam_bytes(0, 100, flags::UNMAPPED, b"rea", &[(10 << 4) | 0], 10, -1, -1, &[]);
        assert_eq!(unclipped_5prime_from_raw_bam(&rec), 0);
    }

    #[test]
    fn test_unclipped_5prime_from_raw_bam_no_cigar() {
        let rec = make_bam_bytes(0, 100, 0, b"rea", &[], 0, -1, -1, &[]);
        assert_eq!(unclipped_5prime_from_raw_bam(&rec), i32::MAX);
    }

    #[test]
    fn test_unclipped_5prime_from_raw_bam_forward() {
        // 5S10M: forward, pos_0based=100 → 1-based=101, unclipped = 101 - 5 = 96
        let cigar = &[(5 << 4) | 4, (10 << 4) | 0]; // 5S10M
        let rec = make_bam_bytes(0, 100, 0, b"rea", cigar, 15, -1, -1, &[]);
        assert_eq!(unclipped_5prime_from_raw_bam(&rec), 96);
    }

    #[test]
    fn test_unclipped_5prime_from_raw_bam_reverse() {
        // 10M5S: reverse, pos_0based=100 → 1-based=101, end = 101 + 10 + 5 - 1 = 115
        let cigar = &[(10 << 4) | 0, (5 << 4) | 4]; // 10M5S
        let rec = make_bam_bytes(0, 100, flags::REVERSE, b"rea", cigar, 15, -1, -1, &[]);
        assert_eq!(unclipped_5prime_from_raw_bam(&rec), 115);
    }

    // ========================================================================
    // reference_length_from_raw_bam tests
    // ========================================================================

    #[test]
    fn test_reference_length_from_raw_bam_no_cigar() {
        let rec = make_bam_bytes(0, 100, 0, b"rea", &[], 0, -1, -1, &[]);
        assert_eq!(reference_length_from_raw_bam(&rec), 0);
    }

    #[test]
    fn test_reference_length_from_raw_bam_simple() {
        // 50M
        let cigar = &[(50 << 4) | 0];
        let rec = make_bam_bytes(0, 100, 0, b"rea", cigar, 50, -1, -1, &[]);
        assert_eq!(reference_length_from_raw_bam(&rec), 50);
    }

    #[test]
    fn test_reference_length_from_raw_bam_with_insertions() {
        // 10M5I10M: insertions don't consume reference, ref_len = 20
        let cigar = &[(10 << 4) | 0, (5 << 4) | 1, (10 << 4) | 0];
        let rec = make_bam_bytes(0, 100, 0, b"rea", cigar, 25, -1, -1, &[]);
        assert_eq!(reference_length_from_raw_bam(&rec), 20);
    }

    #[test]
    fn test_find_int_tag_not_found() {
        let aux = [b'X', b'Y', b'C', 10];
        assert_eq!(find_int_tag(&aux, b"ZZ"), None);
    }

    // (find_mi_tag integer type tests consolidated into test_find_mi_tag_by_type above)

    // ========================================================================
    // parse_mi_bytes edge cases
    // ========================================================================

    #[test]
    fn test_find_mi_tag_just_slash_suffix() {
        // MI:Z:/B -> num_part is empty after stripping suffix, should return None
        let aux = b"MIZ/B\x00";
        assert_eq!(find_mi_tag(aux), None);
    }

    // ========================================================================
    // find_string_tag non-Z type
    // ========================================================================

    #[test]
    fn test_find_string_tag_non_z_type_returns_none() {
        // Tag matches but type is 'C' not 'Z'
        let aux = [b'R', b'X', b'C', 42];
        assert_eq!(find_string_tag(&aux, b"RX"), None);
    }

    #[test]
    fn test_aux_data_offset_from_record_truncated_no_panic() {
        // aux_data_offset_from_record returns None for records too short to read header
        assert_eq!(aux_data_offset_from_record(&[0u8; 19]), None);
        assert!(aux_data_offset_from_record(&[0u8; 20]).is_some());
    }

    #[test]
    fn test_aux_data_slice_out_of_bounds_returns_empty() {
        // When aux_data_offset exceeds record length, aux_data_slice returns empty
        let mut bam = make_bam_bytes(0, 0, 0, b"rea", &[], 10, -1, -1, &[]);

        // Modify l_seq to unreasonably large value — offset will exceed record length
        bam[16..20].copy_from_slice(&1000u32.to_le_bytes());

        let offset = aux_data_offset_from_record(&bam).unwrap();
        assert!(offset > bam.len());

        // aux_data_slice returns empty instead of panicking
        assert!(aux_data_slice(&bam).is_empty());
    }

    // ========================================================================
    // tag_value_size edge cases
    // ========================================================================

    #[test]
    fn test_tag_value_size_unknown_type() {
        // Unknown type byte → None
        assert_eq!(tag_value_size(b'?', &[0; 10]), None);
    }

    #[test]
    fn test_tag_value_size_b_array_invalid_elem_type() {
        // B-array with elem_type 'Z' (not a fixed-size type) → None
        let mut data = vec![b'Z']; // invalid elem_type
        data.extend_from_slice(&3u32.to_le_bytes()); // count
        data.extend_from_slice(&[0; 12]); // data
        assert_eq!(tag_value_size(b'B', &data), None);
    }

    #[test]
    fn test_tag_value_size_b_array_truncated_header() {
        // B-array with less than 5 bytes of header → None
        let data = [b'i', 0, 0]; // only 3 bytes, need 5
        assert_eq!(tag_value_size(b'B', &data), None);
    }

    #[test]
    fn test_tag_value_size_string_no_null() {
        // Z-type string without null terminator → None
        let data = b"hello";
        assert_eq!(tag_value_size(b'Z', data), None);
    }

    // ========================================================================
    // find_string_tag_in_record edge cases
    // ========================================================================

    #[test]
    fn test_find_string_tag_in_record_no_aux() {
        // Record with no aux data (aux_start >= len)
        let rec = make_bam_bytes(0, 0, 0, b"rea", &[], 0, -1, -1, &[]);
        assert_eq!(find_string_tag_in_record(&rec, b"RX"), None);
    }

    // ========================================================================
    // get_cigar_ops edge cases
    // ========================================================================

    #[test]
    fn test_get_cigar_ops_no_cigar() {
        let rec = make_bam_bytes(0, 100, 0, b"rea", &[], 0, -1, -1, &[]);
        assert!(get_cigar_ops(&rec).is_empty());
    }

    #[test]
    fn test_get_cigar_ops_truncated_record() {
        // Record claims n_cigar_op=5 but is too short to hold that many
        let mut rec = make_bam_bytes(0, 100, 0, b"rea", &[(10 << 4) | 0], 10, -1, -1, &[]);
        // Override n_cigar_op to a large value
        rec[12..14].copy_from_slice(&100u16.to_le_bytes());
        assert!(get_cigar_ops(&rec).is_empty());
    }

    #[test]
    fn test_get_cigar_ops_multiple_ops() {
        let cigar = &[(5 << 4) | 4, (10 << 4) | 0, (3 << 4) | 4]; // 5S10M3S
        let rec = make_bam_bytes(0, 100, 0, b"rea", cigar, 18, -1, -1, &[]);
        let ops = get_cigar_ops(&rec);
        assert_eq!(ops.len(), 3);
        assert_eq!(ops[0], (5 << 4) | 4);
        assert_eq!(ops[1], (10 << 4) | 0);
        assert_eq!(ops[2], (3 << 4) | 4);
    }

    // ========================================================================
    // unclipped_5prime_from_raw_bam truncated record
    // ========================================================================

    #[test]
    fn test_unclipped_5prime_from_raw_bam_truncated_cigar() {
        // Record with n_cigar_op > 0 but record too short to hold CIGAR data
        let mut rec = make_bam_bytes(0, 100, 0, b"rea", &[(10 << 4) | 0], 10, -1, -1, &[]);
        rec[12..14].copy_from_slice(&100u16.to_le_bytes()); // claim 100 CIGAR ops
        assert_eq!(unclipped_5prime_from_raw_bam(&rec), i32::MAX);
    }

    // ========================================================================
    // reference_length_from_raw_bam truncated record
    // ========================================================================

    #[test]
    fn test_reference_length_from_raw_bam_truncated_cigar() {
        let mut rec = make_bam_bytes(0, 100, 0, b"rea", &[(10 << 4) | 0], 10, -1, -1, &[]);
        rec[12..14].copy_from_slice(&100u16.to_le_bytes());
        assert_eq!(reference_length_from_raw_bam(&rec), 0);
    }

    // (find_int_tag type tests consolidated into test_find_int_tag_by_type above)

    // ========================================================================
    // find_uint8_tag: tag found but wrong type
    // ========================================================================

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
    // find_float_tag: not found, wrong type
    // ========================================================================

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
    // find_tag_bounds edge cases
    // ========================================================================

    #[test]
    fn test_find_tag_bounds_empty() {
        assert_eq!(find_tag_bounds(&[], b"RX"), None);
    }

    #[test]
    fn test_find_tag_bounds_truncated() {
        // Only 2 bytes — not enough for a tag entry
        assert_eq!(find_tag_bounds(b"RX", b"RX"), None);
    }

    // ========================================================================
    // extract_aux_string_tags: MC tag with invalid UTF-8
    // ========================================================================

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

    // (find_mi_tag positive signed & unknown type tests consolidated into test_find_mi_tag_by_type above)

    // ========================================================================
    // find_mi_tag_in_record: no aux data
    // ========================================================================

    #[test]
    fn test_find_mi_tag_in_record_no_aux() {
        let rec = make_bam_bytes(0, 0, 0, b"rea", &[], 0, -1, -1, &[]);
        assert_eq!(find_mi_tag_in_record(&rec), (0, true));
    }

    #[test]
    fn test_find_mc_tag_in_record_no_aux() {
        let rec = make_bam_bytes(0, 0, 0, b"rea", &[], 0, -1, -1, &[]);
        assert_eq!(find_mc_tag_in_record(&rec), None);
    }

    // ========================================================================
    // read_name edge cases
    // ========================================================================

    #[test]
    fn test_read_name_empty() {
        // l_read_name = 1 means just the null terminator → empty name
        let mut rec = vec![0u8; 34];
        rec[8] = 1; // l_read_name = 1
        rec[32] = 0; // null terminator only
        assert_eq!(read_name(&rec), b"");
    }

    #[test]
    fn test_read_name_zero_length() {
        // l_read_name = 0 → empty name (edge case)
        let mut rec = vec![0u8; 33];
        rec[8] = 0;
        assert_eq!(read_name(&rec), b"");
    }

    // ========================================================================
    // mask_base: odd position already tested, verify even position nibble isolation
    // ========================================================================

    #[test]
    fn test_mask_base_even_position() {
        let mut rec = make_bam_bytes(0, 0, 0, b"rea", &[], 4, -1, -1, &[]);
        let so = seq_offset(&rec);
        rec[so] = 0x12; // A(1),C(2)
        rec[so + 1] = 0x48; // G(4),T(8)

        // Mask position 0 (even, high nibble)
        mask_base(&mut rec, so, 0);
        assert_eq!(get_base(&rec, so, 0), 0xF); // N
        assert_eq!(get_base(&rec, so, 1), 2); // C preserved

        // Mask position 3 (odd, low nibble)
        mask_base(&mut rec, so, 3);
        assert_eq!(get_base(&rec, so, 3), 0xF); // N
        assert_eq!(get_base(&rec, so, 2), 4); // G preserved
    }

    // ========================================================================
    // canonical_template_pos_unclipped coverage (via compare_template_coordinate_raw)
    // ========================================================================

    #[test]
    fn test_compare_template_coordinate_raw_unmapped_unpaired() {
        // Unmapped, unpaired → (MAX, MAX, MAX, MAX, false, false, false)
        let a = make_bam_bytes(-1, -1, flags::UNMAPPED, b"aaa", &[], 0, -1, -1, &[]);
        let b = make_bam_bytes(-1, -1, flags::UNMAPPED, b"zzz", &[], 0, -1, -1, &[]);
        // Both fully unmapped, compare by name
        assert_eq!(compare_template_coordinate_raw(&a, &b), Ordering::Less);
    }

    #[test]
    fn test_compare_template_coordinate_raw_unmapped_with_mapped_mate() {
        // Unmapped read with mapped mate: uses mate's position
        let cigar = &[(10 << 4) | 0]; // 10M
        let a = make_bam_bytes(0, -1, flags::UNMAPPED | flags::PAIRED, b"rea", &[], 0, 0, 100, &[]);
        let b = make_bam_bytes(0, 200, flags::PAIRED, b"rea", cigar, 10, 0, 100, &[]);
        // a is unmapped with mate at tid=0,pos=100; b is mapped at tid=0,pos=200
        // a sorts by mate position (100) which is < b's position (200)
        let cmp = compare_template_coordinate_raw(&a, &b);
        assert_ne!(cmp, Ordering::Equal);
    }

    #[test]
    fn test_compare_template_coordinate_raw_mapped_with_unmapped_mate() {
        // Mapped read with unmapped mate: uses MAX for mate position
        let cigar = &[(10 << 4) | 0]; // 10M
        let a = make_bam_bytes(
            0,
            100,
            flags::PAIRED | flags::MATE_UNMAPPED,
            b"rea",
            cigar,
            10,
            -1,
            -1,
            &[],
        );
        let b = make_bam_bytes(
            0,
            200,
            flags::PAIRED | flags::MATE_UNMAPPED,
            b"rea",
            cigar,
            10,
            -1,
            -1,
            &[],
        );
        // Both mapped, mates unmapped. a at pos 100 < b at pos 200
        assert_eq!(compare_template_coordinate_raw(&a, &b), Ordering::Less);
    }

    #[test]
    fn test_compare_template_coordinate_raw_both_mapped_canonical_swap() {
        // Both mapped, paired. Test the canonical swap (is_upper branch)
        let cigar = &[(10 << 4) | 0]; // 10M
        // a: tid=0,pos=200 with mate at tid=0,pos=100 → is_upper=true → swap
        let a = make_bam_bytes(0, 200, flags::PAIRED, b"rea", cigar, 10, 0, 100, &[]);
        // b: tid=0,pos=100 with mate at tid=0,pos=200 → is_upper=false → no swap
        let b = make_bam_bytes(0, 100, flags::PAIRED, b"rea", cigar, 10, 0, 200, &[]);
        // After canonical ordering, both should produce same (tid1=0,tid2=0,pos1=~100,pos2=~200)
        // so they compare equal on positions, differentiated by is_upper
        let cmp = compare_template_coordinate_raw(&a, &b);
        // a is_upper=true, b is_upper=false → a > b (true > false)
        assert_eq!(cmp, Ordering::Greater);
    }

    #[test]
    fn test_compare_template_coordinate_raw_neg_strand_tiebreak() {
        // Test the neg1/neg2 tiebreak in compare_template_coordinate_raw
        let cigar = &[(10 << 4) | 0]; // 10M
        // a: forward strand
        let a = make_bam_bytes(
            0,
            100,
            flags::PAIRED | flags::MATE_UNMAPPED,
            b"rea",
            cigar,
            10,
            -1,
            -1,
            &[],
        );
        // b: reverse strand at same position
        let b = make_bam_bytes(
            0,
            100,
            flags::PAIRED | flags::MATE_UNMAPPED | flags::REVERSE,
            b"rea",
            cigar,
            10,
            -1,
            -1,
            &[],
        );
        // Reverse strand sorts before forward in samtools convention (neg1=true < neg1=false)
        let cmp = compare_template_coordinate_raw(&a, &b);
        assert_ne!(cmp, Ordering::Equal);
    }

    // ========================================================================
    // compare_mi_tags_raw (exercised through compare_template_coordinate_raw)
    // ========================================================================

    #[test]
    fn test_compare_template_coordinate_raw_mi_tag_tiebreak() {
        // Two identical records except for MI tag value
        let cigar = &[(10 << 4) | 0]; // 10M
        let mut aux_a = Vec::new();
        aux_a.extend_from_slice(b"MIZ10\x00");
        let a = make_bam_bytes(
            0,
            100,
            flags::PAIRED | flags::MATE_UNMAPPED,
            b"rea",
            cigar,
            10,
            -1,
            -1,
            &aux_a,
        );
        let mut aux_b = Vec::new();
        aux_b.extend_from_slice(b"MIZ20\x00");
        let b = make_bam_bytes(
            0,
            100,
            flags::PAIRED | flags::MATE_UNMAPPED,
            b"rea",
            cigar,
            10,
            -1,
            -1,
            &aux_b,
        );
        // MI 10 < MI 20
        assert_eq!(compare_template_coordinate_raw(&a, &b), Ordering::Less);
    }

    // ========================================================================
    // find_mc_tag edge cases
    // ========================================================================

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
    // unclipped_start/end from_raw_cigar: leading S+H, trailing mixed
    // ========================================================================

    #[test]
    fn test_unclipped_5prime_from_raw_bam_hard_and_soft_clip() {
        // 3H5S10M4S2H: forward → start = 101 - 3 - 5 = 93
        let cigar = &[
            (3 << 4) | 5,  // 3H
            (5 << 4) | 4,  // 5S
            (10 << 4) | 0, // 10M
            (4 << 4) | 4,  // 4S
            (2 << 4) | 5,  // 2H
        ];
        let rec = make_bam_bytes(0, 100, 0, b"rea", cigar, 15, -1, -1, &[]);
        // Forward: pos_1based=101, leading clips=3+5=8, unclipped_start = 101 - 8 = 93
        assert_eq!(unclipped_5prime_from_raw_bam(&rec), 93);
    }

    #[test]
    fn test_unclipped_5prime_from_raw_bam_reverse_hard_and_soft_clip() {
        // 3H5S10M4S2H: reverse → end = 101 + 10 + 4 + 2 - 1 = 116
        let cigar = &[
            (3 << 4) | 5,  // 3H
            (5 << 4) | 4,  // 5S
            (10 << 4) | 0, // 10M
            (4 << 4) | 4,  // 4S
            (2 << 4) | 5,  // 2H
        ];
        let rec = make_bam_bytes(0, 100, flags::REVERSE, b"rea", cigar, 15, -1, -1, &[]);
        // Reverse: pos_1based=101, ref_len=10, trailing_clips=4+2=6, end = 101 + 10 + 6 - 1 = 116
        assert_eq!(unclipped_5prime_from_raw_bam(&rec), 116);
    }

    // ========================================================================
    // reference_length_from_raw_bam: mixed ops
    // ========================================================================

    #[test]
    fn test_reference_length_from_raw_bam_mixed_ops() {
        // 5S10M2D3N5M3I8X2H: ref consuming = 10M+2D+3N+5M+8X = 28
        let cigar = &[
            (5 << 4) | 4,  // 5S
            (10 << 4) | 0, // 10M
            (2 << 4) | 2,  // 2D
            (3 << 4) | 3,  // 3N
            (5 << 4) | 0,  // 5M
            (3 << 4) | 1,  // 3I
            (8 << 4) | 8,  // 8X (= mismatch, consumes ref)
            (2 << 4) | 5,  // 2H
        ];
        let rec = make_bam_bytes(0, 100, 0, b"rea", cigar, 31, -1, -1, &[]);
        assert_eq!(reference_length_from_raw_bam(&rec), 28);
    }

    // ========================================================================
    // find_string_tag: truncated Z string (no null terminator found)
    // ========================================================================

    #[test]
    fn test_find_string_tag_truncated_z_value() {
        // RX:Z:hello but no null terminator — should return None
        let aux = b"RXZhello";
        assert_eq!(find_string_tag(aux, b"RX"), None);
    }

    // ========================================================================
    // extract_aux_string_tags: truncated Z value
    // ========================================================================

    #[test]
    fn test_extract_aux_string_tags_truncated_z() {
        // RG:Z:lib1 with no null → should break out, no tags found
        let aux = b"RGZlib1";
        let result = extract_aux_string_tags(aux, b"CB");
        assert!(result.rg.is_none());
    }

    // ========================================================================
    // find_tag_bounds: with unknown tag type (breaks out of loop)
    // ========================================================================

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
    // remove_tag: record with no aux data
    // ========================================================================

    #[test]
    fn test_remove_tag_no_aux_data() {
        // Record where aux_data_offset >= record length
        let mut rec = make_bam_bytes(0, 0, 0, b"rea", &[], 0, -1, -1, &[]);
        let orig_len = rec.len();
        remove_tag(&mut rec, b"MI");
        assert_eq!(rec.len(), orig_len); // no-op
    }

    // ========================================================================
    // update_string_tag: record with no aux data (append path)
    // ========================================================================

    #[test]
    fn test_update_string_tag_no_aux_appends() {
        let mut rec = make_bam_bytes(0, 0, 0, b"rea", &[], 0, -1, -1, &[]);
        update_string_tag(&mut rec, b"MI", b"42");
        assert_eq!(find_string_tag_in_record(&rec, b"MI"), Some(b"42".as_ref()));
    }

    // ========================================================================
    // Mate unclipped positions: string CIGAR helpers
    // ========================================================================

    #[test]
    fn test_unclipped_other_start_no_clips() {
        assert_eq!(unclipped_other_start(100, "10M"), 100);
    }

    #[test]
    fn test_unclipped_other_end_no_trailing_clips() {
        // 10M: end = 100 + 10 + 0 - 1 = 109
        assert_eq!(unclipped_other_end(100, "10M"), 109);
    }

    #[test]
    fn test_unclipped_other_start_complex() {
        // 3H5S10M: leading = 3+5=8, start = 100 - 8 = 92
        assert_eq!(unclipped_other_start(100, "3H5S10M"), 92);
    }

    #[test]
    fn test_unclipped_other_end_complex() {
        // 10M5S3H: ref=10, trailing=5+3=8, end = 100 + 10 + 8 - 1 = 117
        assert_eq!(unclipped_other_end(100, "10M5S3H"), 117);
    }

    // ========================================================================
    // parse_leading_clips: only clips (no non-clip ops)
    // ========================================================================

    #[test]
    fn test_parse_leading_clips_all_clips() {
        // "5S3H" — all clips, no ref ops
        assert_eq!(parse_leading_clips("5S3H"), 8);
    }

    // ========================================================================
    // parse_ref_len_and_trailing_clips: no ref ops
    // ========================================================================

    #[test]
    fn test_parse_ref_len_and_trailing_clips_no_ref_ops() {
        // "5S3H" — no ref-consuming ops
        assert_eq!(parse_ref_len_and_trailing_clips("5S3H"), (0, 0));
    }

    #[test]
    fn test_parse_ref_len_and_trailing_clips_eq_and_x() {
        // "10=3X" → ref_len=13, trailing_clips=0
        assert_eq!(parse_ref_len_and_trailing_clips("10=3X"), (13, 0));
    }

    // ========================================================================
    // compare_names_raw: zero-length names
    // ========================================================================

    #[test]
    fn test_compare_names_raw_empty_names() {
        let mut a = vec![0u8; 33];
        a[8] = 0; // l_read_name = 0
        let mut b = vec![0u8; 33];
        b[8] = 0;
        assert_eq!(compare_names_raw(&a, &b), Ordering::Equal);
    }

    // ========================================================================
    // Tag lookup coverage for all BAM aux tag types
    //
    // The pa tag (used in dedup for secondary/supplementary tracking) is a
    // B:i array.  find_string_tag only handles Z-type, find_int_tag only
    // handles scalar int types (c/C/s/S/i/I).  Neither can find B-array
    // tags, so the dedup pa-tag validation falsely reports them as missing.
    //
    // find_tag_type handles all types correctly since it only checks for
    // tag presence regardless of type.
    // ========================================================================

    /// Helper: build raw aux bytes for a B-type array tag.
    fn make_b_array_tag(tag: [u8; 2], elem_type: u8, count: u32, elements: &[u8]) -> Vec<u8> {
        let mut aux = vec![tag[0], tag[1], b'B', elem_type];
        aux.extend_from_slice(&count.to_le_bytes());
        aux.extend_from_slice(elements);
        aux
    }

    fn make_b_int_array_tag(tag: [u8; 2], values: &[i32]) -> Vec<u8> {
        let bytes: Vec<u8> = values.iter().flat_map(|v| v.to_le_bytes()).collect();
        make_b_array_tag(tag, b'i', values.len() as u32, &bytes)
    }

    fn make_b_float_array_tag(tag: [u8; 2], values: &[f32]) -> Vec<u8> {
        let bytes: Vec<u8> = values.iter().flat_map(|v| v.to_le_bytes()).collect();
        make_b_array_tag(tag, b'f', values.len() as u32, &bytes)
    }

    fn make_b_uint8_array_tag(tag: [u8; 2], values: &[u8]) -> Vec<u8> {
        make_b_array_tag(tag, b'C', values.len() as u32, values)
    }

    fn make_b_int8_array_tag(tag: [u8; 2], values: &[i8]) -> Vec<u8> {
        let bytes: Vec<u8> = values.iter().map(|&v| v as u8).collect();
        make_b_array_tag(tag, b'c', values.len() as u32, &bytes)
    }

    fn make_b_int16_array_tag(tag: [u8; 2], values: &[i16]) -> Vec<u8> {
        let bytes: Vec<u8> = values.iter().flat_map(|v| v.to_le_bytes()).collect();
        make_b_array_tag(tag, b's', values.len() as u32, &bytes)
    }

    fn make_b_uint16_array_tag(tag: [u8; 2], values: &[u16]) -> Vec<u8> {
        let bytes: Vec<u8> = values.iter().flat_map(|v| v.to_le_bytes()).collect();
        make_b_array_tag(tag, b'S', values.len() as u32, &bytes)
    }

    fn make_b_uint32_array_tag(tag: [u8; 2], values: &[u32]) -> Vec<u8> {
        let bytes: Vec<u8> = values.iter().flat_map(|v| v.to_le_bytes()).collect();
        make_b_array_tag(tag, b'I', values.len() as u32, &bytes)
    }

    // --- B:c array (int8 array) ---

    #[test]
    fn test_find_string_tag_cannot_find_b_int8_array() {
        let aux = make_b_int8_array_tag(*b"Xc", &[-1, 0, 1]);
        assert_eq!(find_string_tag(&aux, b"Xc"), None);
    }

    #[test]
    fn test_find_int_tag_cannot_find_b_int8_array() {
        let aux = make_b_int8_array_tag(*b"Xc", &[-1, 0, 1]);
        assert_eq!(find_int_tag(&aux, b"Xc"), None);
    }

    #[test]
    fn test_find_tag_type_finds_b_int8_array() {
        let aux = make_b_int8_array_tag(*b"Xc", &[-1, 0, 1]);
        assert_eq!(find_tag_type(&aux, b"Xc"), Some(b'B'));
    }

    // --- B:s array (int16 array) ---

    #[test]
    fn test_find_string_tag_cannot_find_b_int16_array() {
        let aux = make_b_int16_array_tag(*b"Xs", &[-100, 0, 200]);
        assert_eq!(find_string_tag(&aux, b"Xs"), None);
    }

    #[test]
    fn test_find_int_tag_cannot_find_b_int16_array() {
        let aux = make_b_int16_array_tag(*b"Xs", &[-100, 0, 200]);
        assert_eq!(find_int_tag(&aux, b"Xs"), None);
    }

    #[test]
    fn test_find_tag_type_finds_b_int16_array() {
        let aux = make_b_int16_array_tag(*b"Xs", &[-100, 0, 200]);
        assert_eq!(find_tag_type(&aux, b"Xs"), Some(b'B'));
    }

    // --- B:S array (uint16 array) ---

    #[test]
    fn test_find_string_tag_cannot_find_b_uint16_array() {
        let aux = make_b_uint16_array_tag(*b"XS", &[100, 200, 300]);
        assert_eq!(find_string_tag(&aux, b"XS"), None);
    }

    #[test]
    fn test_find_int_tag_cannot_find_b_uint16_array() {
        let aux = make_b_uint16_array_tag(*b"XS", &[100, 200, 300]);
        assert_eq!(find_int_tag(&aux, b"XS"), None);
    }

    #[test]
    fn test_find_tag_type_finds_b_uint16_array() {
        let aux = make_b_uint16_array_tag(*b"XS", &[100, 200, 300]);
        assert_eq!(find_tag_type(&aux, b"XS"), Some(b'B'));
    }

    // --- B:I array (uint32 array) ---

    #[test]
    fn test_find_string_tag_cannot_find_b_uint32_array() {
        let aux = make_b_uint32_array_tag(*b"XI", &[1000, 2000, 3000]);
        assert_eq!(find_string_tag(&aux, b"XI"), None);
    }

    #[test]
    fn test_find_int_tag_cannot_find_b_uint32_array() {
        let aux = make_b_uint32_array_tag(*b"XI", &[1000, 2000, 3000]);
        assert_eq!(find_int_tag(&aux, b"XI"), None);
    }

    #[test]
    fn test_find_tag_type_finds_b_uint32_array() {
        let aux = make_b_uint32_array_tag(*b"XI", &[1000, 2000, 3000]);
        assert_eq!(find_tag_type(&aux, b"XI"), Some(b'B'));
    }

    // --- B:i array (the pa tag type) ---

    #[test]
    fn test_find_string_tag_cannot_find_b_int_array() {
        let aux = make_b_int_array_tag(*b"pa", &[0, 27_056_961, 0, 207, 60005, 1]);
        // BUG: find_string_tag returns None for B:i tags
        assert_eq!(find_string_tag(&aux, b"pa"), None);
    }

    #[test]
    fn test_find_int_tag_cannot_find_b_int_array() {
        let aux = make_b_int_array_tag(*b"pa", &[0, 27_056_961, 0, 207, 60005, 1]);
        // BUG: find_int_tag returns None for B:i tags
        assert_eq!(find_int_tag(&aux, b"pa"), None);
    }

    #[test]
    fn test_find_tag_type_finds_b_int_array() {
        let aux = make_b_int_array_tag(*b"pa", &[0, 27_056_961, 0, 207, 60005, 1]);
        // find_tag_type correctly finds B-array tags
        assert_eq!(find_tag_type(&aux, b"pa"), Some(b'B'));
    }

    // --- B:f array ---

    #[test]
    fn test_find_string_tag_cannot_find_b_float_array() {
        let aux = make_b_float_array_tag(*b"XF", &[1.0, 2.5, 3.0]);
        assert_eq!(find_string_tag(&aux, b"XF"), None);
    }

    #[test]
    fn test_find_int_tag_cannot_find_b_float_array() {
        let aux = make_b_float_array_tag(*b"XF", &[1.0, 2.5, 3.0]);
        assert_eq!(find_int_tag(&aux, b"XF"), None);
    }

    #[test]
    fn test_find_tag_type_finds_b_float_array() {
        let aux = make_b_float_array_tag(*b"XF", &[1.0, 2.5, 3.0]);
        assert_eq!(find_tag_type(&aux, b"XF"), Some(b'B'));
    }

    // --- B:C array (uint8 array) ---

    #[test]
    fn test_find_string_tag_cannot_find_b_uint8_array() {
        let aux = make_b_uint8_array_tag(*b"XC", &[10, 20, 30]);
        assert_eq!(find_string_tag(&aux, b"XC"), None);
    }

    #[test]
    fn test_find_int_tag_cannot_find_b_uint8_array() {
        let aux = make_b_uint8_array_tag(*b"XC", &[10, 20, 30]);
        assert_eq!(find_int_tag(&aux, b"XC"), None);
    }

    #[test]
    fn test_find_tag_type_finds_b_uint8_array() {
        let aux = make_b_uint8_array_tag(*b"XC", &[10, 20, 30]);
        assert_eq!(find_tag_type(&aux, b"XC"), Some(b'B'));
    }

    // --- H (hex string) type ---

    #[test]
    fn test_find_string_tag_cannot_find_h_type() {
        // H-type hex string: tag, 'H', hex bytes, NUL
        let aux: &[u8] = b"XHH1A2B\x00";
        // find_string_tag only handles Z, not H
        assert_eq!(find_string_tag(aux, b"XH"), None);
    }

    #[test]
    fn test_find_int_tag_cannot_find_h_type() {
        let aux: &[u8] = b"XHH1A2B\x00";
        assert_eq!(find_int_tag(aux, b"XH"), None);
    }

    #[test]
    fn test_find_tag_type_finds_h_type() {
        let aux: &[u8] = b"XHH1A2B\x00";
        assert_eq!(find_tag_type(aux, b"XH"), Some(b'H'));
    }

    // --- Z (NUL-terminated string) type ---

    #[test]
    fn test_find_string_tag_finds_z_type() {
        let aux: &[u8] = b"RXZhello\x00";
        assert_eq!(find_string_tag(aux, b"RX"), Some(b"hello".as_ref()));
    }

    #[test]
    fn test_find_int_tag_cannot_find_z_type() {
        let aux: &[u8] = b"RXZhello\x00";
        assert_eq!(find_int_tag(aux, b"RX"), None);
    }

    #[test]
    fn test_find_tag_type_finds_z_type() {
        let aux: &[u8] = b"RXZhello\x00";
        assert_eq!(find_tag_type(aux, b"RX"), Some(b'Z'));
    }

    // --- c (int8) type ---

    #[test]
    fn test_find_string_tag_cannot_find_c_type() {
        let aux: &[u8] = &[b'X', b'c', b'c', 0xFE]; // Xc:c:-2
        assert_eq!(find_string_tag(aux, b"Xc"), None);
    }

    #[test]
    fn test_find_int_tag_finds_c_type() {
        let aux: &[u8] = &[b'X', b'c', b'c', 5];
        assert_eq!(find_int_tag(aux, b"Xc"), Some(5));
    }

    #[test]
    fn test_find_tag_type_finds_c_type() {
        let aux: &[u8] = &[b'X', b'c', b'c', 5];
        assert_eq!(find_tag_type(aux, b"Xc"), Some(b'c'));
    }

    // --- C (uint8) type ---

    #[test]
    fn test_find_string_tag_cannot_find_upper_c_type() {
        let aux: &[u8] = &[b'X', b'C', b'C', 200];
        assert_eq!(find_string_tag(aux, b"XC"), None);
    }

    #[test]
    fn test_find_int_tag_finds_upper_c_type() {
        let aux: &[u8] = &[b'N', b'M', b'C', 42];
        assert_eq!(find_int_tag(aux, b"NM"), Some(42));
    }

    #[test]
    fn test_find_tag_type_finds_upper_c_type() {
        let aux: &[u8] = &[b'N', b'M', b'C', 42];
        assert_eq!(find_tag_type(aux, b"NM"), Some(b'C'));
    }

    // --- s (int16) type ---

    #[test]
    fn test_find_string_tag_cannot_find_s_type() {
        let val = 300i16.to_le_bytes();
        let aux: &[u8] = &[b'X', b's', b's', val[0], val[1]];
        assert_eq!(find_string_tag(aux, b"Xs"), None);
    }

    #[test]
    fn test_find_int_tag_finds_s_type() {
        let val = 300i16.to_le_bytes();
        let aux: &[u8] = &[b'X', b's', b's', val[0], val[1]];
        assert_eq!(find_int_tag(aux, b"Xs"), Some(300));
    }

    #[test]
    fn test_find_tag_type_finds_s_type() {
        let val = 300i16.to_le_bytes();
        let aux: &[u8] = &[b'X', b's', b's', val[0], val[1]];
        assert_eq!(find_tag_type(aux, b"Xs"), Some(b's'));
    }

    // --- S (uint16) type ---

    #[test]
    fn test_find_string_tag_cannot_find_upper_s_type() {
        let val = 50_000_u16.to_le_bytes();
        let aux: &[u8] = &[b'X', b'S', b'S', val[0], val[1]];
        assert_eq!(find_string_tag(aux, b"XS"), None);
    }

    #[test]
    fn test_find_int_tag_finds_upper_s_type() {
        let val = 50_000_u16.to_le_bytes();
        let aux: &[u8] = &[b'X', b'S', b'S', val[0], val[1]];
        assert_eq!(find_int_tag(aux, b"XS"), Some(50_000));
    }

    #[test]
    fn test_find_tag_type_finds_upper_s_type() {
        let val = 50_000_u16.to_le_bytes();
        let aux: &[u8] = &[b'X', b'S', b'S', val[0], val[1]];
        assert_eq!(find_tag_type(aux, b"XS"), Some(b'S'));
    }

    // --- i (int32) type ---

    #[test]
    fn test_find_string_tag_cannot_find_i_type() {
        let val = 100_000_i32.to_le_bytes();
        let aux: &[u8] = &[b'X', b'i', b'i', val[0], val[1], val[2], val[3]];
        assert_eq!(find_string_tag(aux, b"Xi"), None);
    }

    #[test]
    fn test_find_int_tag_finds_i_type() {
        let val = 100_000_i32.to_le_bytes();
        let aux: &[u8] = &[b'X', b'i', b'i', val[0], val[1], val[2], val[3]];
        assert_eq!(find_int_tag(aux, b"Xi"), Some(100_000));
    }

    #[test]
    fn test_find_tag_type_finds_i_type() {
        let val = 100_000_i32.to_le_bytes();
        let aux: &[u8] = &[b'X', b'i', b'i', val[0], val[1], val[2], val[3]];
        assert_eq!(find_tag_type(aux, b"Xi"), Some(b'i'));
    }

    // --- I (uint32) type ---

    #[test]
    fn test_find_string_tag_cannot_find_upper_i_type() {
        let val = 3_000_000_000_u32.to_le_bytes();
        let aux: &[u8] = &[b'X', b'I', b'I', val[0], val[1], val[2], val[3]];
        assert_eq!(find_string_tag(aux, b"XI"), None);
    }

    #[test]
    fn test_find_int_tag_finds_upper_i_type() {
        let val = 3_000_000_000_u32.to_le_bytes();
        let aux: &[u8] = &[b'X', b'I', b'I', val[0], val[1], val[2], val[3]];
        // find_int_tag returns i64, so this should work for large uint32 values
        assert_eq!(find_int_tag(aux, b"XI"), Some(3_000_000_000));
    }

    #[test]
    fn test_find_tag_type_finds_upper_i_type() {
        let val = 3_000_000_000_u32.to_le_bytes();
        let aux: &[u8] = &[b'X', b'I', b'I', val[0], val[1], val[2], val[3]];
        assert_eq!(find_tag_type(aux, b"XI"), Some(b'I'));
    }

    // --- A (single char) type ---

    #[test]
    fn test_find_string_tag_cannot_find_a_type() {
        let aux = [b'X', b'A', b'A', b'G']; // XA:A:G
        assert_eq!(find_string_tag(&aux, b"XA"), None);
    }

    #[test]
    fn test_find_int_tag_cannot_find_a_type() {
        let aux = [b'X', b'A', b'A', b'G'];
        assert_eq!(find_int_tag(&aux, b"XA"), None);
    }

    #[test]
    fn test_find_tag_type_finds_a_type() {
        let aux = [b'X', b'A', b'A', b'G'];
        assert_eq!(find_tag_type(&aux, b"XA"), Some(b'A'));
    }

    // --- f (float) type ---

    #[test]
    fn test_find_string_tag_cannot_find_f_type() {
        let mut aux = vec![b'X', b'F', b'f'];
        aux.extend_from_slice(&1.5f32.to_le_bytes());
        assert_eq!(find_string_tag(&aux, b"XF"), None);
    }

    #[test]
    fn test_find_int_tag_cannot_find_f_type() {
        let mut aux = vec![b'X', b'F', b'f'];
        aux.extend_from_slice(&1.5f32.to_le_bytes());
        assert_eq!(find_int_tag(&aux, b"XF"), None);
    }

    #[test]
    fn test_find_tag_type_finds_f_type() {
        let mut aux = vec![b'X', b'F', b'f'];
        aux.extend_from_slice(&1.5f32.to_le_bytes());
        assert_eq!(find_tag_type(&aux, b"XF"), Some(b'f'));
    }

    // --- Tag after a B-array tag (verifies traversal works) ---

    #[test]
    fn test_find_string_tag_after_b_array() {
        let mut aux = make_b_int_array_tag(*b"pa", &[1, 2, 3]);
        aux.extend_from_slice(b"RXZhello\x00");
        // Can find the Z-tag after the B-array
        assert_eq!(find_string_tag(&aux, b"RX"), Some(b"hello".as_ref()));
    }

    #[test]
    fn test_find_int_tag_after_b_array() {
        let mut aux = make_b_int_array_tag(*b"pa", &[1, 2, 3]);
        aux.extend_from_slice(&[b'N', b'M', b'C', 5]); // NM:C:5
        assert_eq!(find_int_tag(&aux, b"NM"), Some(5));
    }

    #[test]
    fn test_find_tag_type_after_b_array() {
        let mut aux = make_b_int_array_tag(*b"pa", &[1, 2, 3]);
        aux.extend_from_slice(b"RXZhello\x00");
        assert_eq!(find_tag_type(&aux, b"RX"), Some(b'Z'));
    }

    // --- Dedup pa-tag validation bug reproduction ---
    // This test mirrors the exact check in dedup.rs lines 932-935

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

    // ========================================================================
    // pack_sequence_into
    // ========================================================================

    #[test]
    fn test_pack_sequence_into_empty() {
        let mut dst = Vec::new();
        pack_sequence_into(&mut dst, b"");
        assert!(dst.is_empty());
    }

    #[test]
    fn test_pack_sequence_into_even() {
        let mut dst = Vec::new();
        pack_sequence_into(&mut dst, b"ACGT");
        // A=1, C=2, G=4, T=8
        assert_eq!(dst, [0x12, 0x48]);
    }

    #[test]
    fn test_pack_sequence_into_odd() {
        let mut dst = Vec::new();
        pack_sequence_into(&mut dst, b"ACG");
        // A=1, C=2, G=4, pad=0
        assert_eq!(dst, [0x12, 0x40]);
    }

    #[test]
    fn test_pack_sequence_into_single_base() {
        let mut dst = Vec::new();
        pack_sequence_into(&mut dst, b"T");
        // T=8, pad=0
        assert_eq!(dst, [0x80]);
    }

    #[test]
    fn test_pack_sequence_into_17_bases() {
        // 17 bases: exercises 16-base chunked path + 1-base remainder
        let mut dst = Vec::new();
        pack_sequence_into(&mut dst, b"ACGTACGTACGTACGTA");
        assert_eq!(dst.len(), 9); // ceil(17/2)
        // First 8 bytes (16 bases chunked)
        assert_eq!(dst[0], 0x12); // AC
        assert_eq!(dst[1], 0x48); // GT
        assert_eq!(dst[2], 0x12); // AC
        assert_eq!(dst[3], 0x48); // GT
        assert_eq!(dst[4], 0x12); // AC
        assert_eq!(dst[5], 0x48); // GT
        assert_eq!(dst[6], 0x12); // AC
        assert_eq!(dst[7], 0x48); // GT
        // Last byte (1 base + padding)
        assert_eq!(dst[8], 0x10); // A=1, pad=0
    }

    #[test]
    fn test_pack_sequence_into_n_bases() {
        let mut dst = Vec::new();
        pack_sequence_into(&mut dst, b"NN");
        assert_eq!(dst, [0xFF]); // N=15, N=15
    }

    // ========================================================================
    // append_int_tag
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
        assert_eq!(rec, [b'c', b'M', b'c', (-5i8) as u8]);
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
    // append_float_tag
    // ========================================================================

    #[test]
    fn test_append_float_tag() {
        let mut rec = Vec::new();
        append_float_tag(&mut rec, b"cE", 0.05);
        let v = 0.05f32.to_le_bytes();
        assert_eq!(rec, [b'c', b'E', b'f', v[0], v[1], v[2], v[3]]);
    }

    // ========================================================================
    // append_i16_array_tag
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
    // UnmappedBamRecordBuilder
    // ========================================================================

    #[test]
    fn test_builder_basic_unmapped_record() {
        let mut builder = UnmappedBamRecordBuilder::new();
        builder.build_record(b"read1", flags::UNMAPPED, b"ACGT", &[30, 25, 35, 40]);

        let bam = builder.as_bytes();

        // Verify fixed header fields via read primitives
        assert_eq!(ref_id(bam), -1);
        assert_eq!(pos(bam), -1);
        assert_eq!(l_read_name(bam), 6); // "read1" + NUL
        assert_eq!(mapq(bam), 0);
        assert_eq!(n_cigar_op(bam), 0);
        assert_eq!(flags(bam), flags::UNMAPPED);
        assert_eq!(l_seq(bam), 4);
        assert_eq!(mate_ref_id(bam), -1);
        assert_eq!(mate_pos(bam), -1);
        assert_eq!(read_name(bam), b"read1");

        // Verify packed sequence via read primitives
        let so = seq_offset(bam);
        assert_eq!(get_base(bam, so, 0), 1); // A
        assert_eq!(get_base(bam, so, 1), 2); // C
        assert_eq!(get_base(bam, so, 2), 4); // G
        assert_eq!(get_base(bam, so, 3), 8); // T

        // Verify quality scores via read primitives
        let qo = qual_offset(bam);
        assert_eq!(get_qual(bam, qo, 0), 30);
        assert_eq!(get_qual(bam, qo, 1), 25);
        assert_eq!(get_qual(bam, qo, 2), 35);
        assert_eq!(get_qual(bam, qo, 3), 40);

        // Verify no aux data
        assert!(aux_data_slice(bam).is_empty());
    }

    #[test]
    fn test_builder_with_tags() {
        let mut builder = UnmappedBamRecordBuilder::new();
        builder.build_record(b"cons:1:ACG", flags::UNMAPPED, b"ACGT", &[30, 30, 30, 30]);
        builder.append_string_tag(b"RG", b"sample1");
        builder.append_string_tag(b"MI", b"42");
        builder.append_int_tag(b"cD", 10);
        builder.append_int_tag(b"cM", 3);
        builder.append_float_tag(b"cE", 0.05);
        builder.append_i16_array_tag(b"cd", &[10, 8, 12, 9]);

        let bam = builder.as_bytes();
        let aux = aux_data_slice(bam);

        // Verify string tags
        assert_eq!(find_string_tag(aux, b"RG"), Some(b"sample1" as &[u8]));
        assert_eq!(find_string_tag(aux, b"MI"), Some(b"42" as &[u8]));

        // Verify integer tags
        assert_eq!(find_int_tag(aux, b"cD"), Some(10));
        assert_eq!(find_int_tag(aux, b"cM"), Some(3));

        // Verify float tag
        let ce = find_float_tag(aux, b"cE").unwrap();
        assert!((ce - 0.05).abs() < 1e-7);
    }

    #[test]
    fn test_builder_reuse() {
        let mut builder = UnmappedBamRecordBuilder::new();

        // Build first record
        builder.build_record(b"r1", flags::UNMAPPED, b"ACGT", &[30, 30, 30, 30]);
        builder.append_string_tag(b"MI", b"1");
        let len1 = builder.as_bytes().len();
        assert_eq!(l_seq(builder.as_bytes()), 4);
        assert_eq!(read_name(builder.as_bytes()), b"r1");

        // Reuse for different-length record
        builder.clear();
        builder.build_record(b"r2", flags::UNMAPPED, b"AC", &[30, 30]);
        builder.append_string_tag(b"MI", b"2");

        let bam = builder.as_bytes();
        assert_ne!(bam.len(), len1);
        assert_eq!(l_seq(bam), 2);
        assert_eq!(read_name(bam), b"r2");
        assert_eq!(find_string_tag(aux_data_slice(bam), b"MI"), Some(b"2" as &[u8]));
    }

    #[test]
    fn test_builder_write_with_block_size() {
        let mut builder = UnmappedBamRecordBuilder::new();
        builder.build_record(b"r1", flags::UNMAPPED, b"ACGT", &[30, 30, 30, 30]);

        let mut output = Vec::new();
        builder.write_with_block_size(&mut output);

        // First 4 bytes = block_size
        let block_size = u32::from_le_bytes([output[0], output[1], output[2], output[3]]);
        assert_eq!(block_size as usize, builder.as_bytes().len());
        assert_eq!(&output[4..], builder.as_bytes());
    }

    #[test]
    fn test_builder_paired_consensus_flags() {
        let mut builder = UnmappedBamRecordBuilder::new();

        // R1
        let r1_flags =
            flags::PAIRED | flags::UNMAPPED | flags::MATE_UNMAPPED | flags::FIRST_SEGMENT;
        builder.build_record(b"cons:1:ACG", r1_flags, b"ACGT", &[30, 30, 30, 30]);
        assert_eq!(flags(builder.as_bytes()), r1_flags);

        // R2
        builder.clear();
        let r2_flags = flags::PAIRED | flags::UNMAPPED | flags::MATE_UNMAPPED | flags::LAST_SEGMENT;
        builder.build_record(b"cons:1:ACG", r2_flags, b"TGCA", &[25, 35, 30, 40]);
        assert_eq!(flags(builder.as_bytes()), r2_flags);
    }

    #[test]
    fn test_builder_empty_sequence() {
        let mut builder = UnmappedBamRecordBuilder::new();
        builder.build_record(b"empty", flags::UNMAPPED, b"", &[]);
        builder.append_string_tag(b"RG", b"rg0");

        let bam = builder.as_bytes();
        assert_eq!(l_seq(bam), 0);
        assert_eq!(find_string_tag(aux_data_slice(bam), b"RG"), Some(b"rg0" as &[u8]));
    }

    #[test]
    fn test_builder_missing_quality() {
        let mut builder = UnmappedBamRecordBuilder::new();
        // Empty quals with non-empty bases -> fills with 0xFF
        builder.build_record(b"r1", flags::UNMAPPED, b"ACGT", &[]);

        let bam = builder.as_bytes();
        let qo = qual_offset(bam);
        for i in 0..4 {
            assert_eq!(get_qual(bam, qo, i), 0xFF);
        }
    }

    #[test]
    fn test_builder_header_matches_vendored_default() {
        // Cross-validate header layout against known bytes from
        // vendored encoder test_encode_with_default_fields.
        // Default RecordBuf: name="*", flags=UNMAPPED, mapq=255, seq=empty
        // Our builder: custom name, flags, mapq=0 (hardcoded for consensus)
        //
        // Verify shared unmapped constants: ref_id=-1, pos=-1, bin=4680,
        // n_cigar_op=0, next_ref_id=-1, next_pos=-1, tlen=0
        let mut builder = UnmappedBamRecordBuilder::new();
        builder.build_record(b"*", flags::UNMAPPED, b"", &[]);

        let bam = builder.as_bytes();
        // ref_id = -1
        assert_eq!(&bam[0..4], &0xFF_FF_FF_FFu32.to_le_bytes());
        // pos = -1
        assert_eq!(&bam[4..8], &0xFF_FF_FF_FFu32.to_le_bytes());
        // l_read_name = 2 ("*" + NUL)
        assert_eq!(bam[8], 0x02);
        // mapq = 0 (consensus-specific, vendored default is 255)
        assert_eq!(bam[9], 0x00);
        // bin = 4680 = 0x1248
        assert_eq!(&bam[10..12], &[0x48, 0x12]);
        // n_cigar_op = 0
        assert_eq!(&bam[12..14], &[0x00, 0x00]);
        // flags = UNMAPPED = 4
        assert_eq!(&bam[14..16], &[0x04, 0x00]);
        // l_seq = 0
        assert_eq!(&bam[16..20], &[0x00, 0x00, 0x00, 0x00]);
        // next_ref_id = -1
        assert_eq!(&bam[20..24], &0xFF_FF_FF_FFu32.to_le_bytes());
        // next_pos = -1
        assert_eq!(&bam[24..28], &0xFF_FF_FF_FFu32.to_le_bytes());
        // tlen = 0
        assert_eq!(&bam[28..32], &[0x00, 0x00, 0x00, 0x00]);
        // name = "*\0"
        assert_eq!(&bam[32..34], &[b'*', 0x00]);
    }

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

    #[test]
    fn test_builder_phred33_string_tag() {
        let mut builder = UnmappedBamRecordBuilder::new();
        builder.build_record(b"test", flags::UNMAPPED, b"ACGT", &[30, 30, 30, 30]);
        builder.append_phred33_string_tag(b"aq", &[0, 10, 30, 40]);
        let bytes = builder.as_bytes();
        let parsed = ParsedBamRecord::from_bytes(bytes);
        let aq = parsed.get_string_tag(b"aq").unwrap();
        assert_eq!(aq, b"!+?I");
    }

    #[test]
    fn test_parsed_bam_record_roundtrip() {
        let mut builder = UnmappedBamRecordBuilder::new();
        builder.build_record(
            b"myread",
            flags::UNMAPPED | flags::PAIRED | flags::FIRST_SEGMENT | flags::MATE_UNMAPPED,
            b"ACGTNN",
            &[30, 25, 20, 15, 10, 5],
        );
        builder.append_string_tag(b"RG", b"sample1");
        builder.append_int_tag(b"cD", 5);
        builder.append_float_tag(b"cE", 0.01);
        builder.append_i16_array_tag(b"cd", &[3, 4, 5, 3, 2, 1]);

        let mut output = Vec::new();
        builder.write_with_block_size(&mut output);

        let records = ParsedBamRecord::parse_all(&output);
        assert_eq!(records.len(), 1);
        let rec = &records[0];
        assert_eq!(rec.name, b"myread");
        assert_eq!(
            rec.flag,
            flags::UNMAPPED | flags::PAIRED | flags::FIRST_SEGMENT | flags::MATE_UNMAPPED
        );
        assert_eq!(rec.bases, b"ACGTNN");
        assert_eq!(rec.quals, vec![30, 25, 20, 15, 10, 5]);
        assert_eq!(rec.get_string_tag(b"RG").unwrap(), b"sample1");
        assert_eq!(rec.get_int_tag(b"cD").unwrap(), 5);
        let ce = rec.get_float_tag(b"cE").unwrap();
        assert!((ce - 0.01).abs() < 0.001);
        assert_eq!(rec.get_i16_array_tag(b"cd").unwrap(), vec![3, 4, 5, 3, 2, 1]);
    }

    // ========================================================================
    // Helper: encode a CIGAR op as raw u32 (op_type in low 4 bits, length in upper 28)
    // ========================================================================

    /// Encode a single CIGAR op.  `op_type`: M=0, I=1, D=2, N=3, S=4, H=5, P=6, `=7`, X=8.
    fn encode_op(op_type: u32, len: usize) -> u32 {
        ((len as u32) << 4) | op_type
    }

    /// Build a `make_bam_bytes` record with a custom template length (tlen).
    #[allow(clippy::too_many_arguments)]
    fn make_bam_bytes_with_tlen(
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

    // ========================================================================
    // upgrade_clipping_raw tests
    // ========================================================================

    #[test]
    fn test_upgrade_clipping_raw_from_start_soft_only() {
        // 5S10M: clip_amount=3, should upgrade 3S to H
        // existing: H=0, S=5, clip_amount=3 => upgrade 3 => new_hard=3, remaining_soft=2
        let cigar = &[encode_op(4, 5), encode_op(0, 10)]; // 5S10M
        let (result, ref_consumed) = upgrade_clipping_raw(cigar, 3, true, encode_op);
        assert_eq!(ref_consumed, 0);
        assert_eq!(result.len(), 3); // 3H, 2S, 10M
        assert_eq!(result[0], encode_op(5, 3)); // 3H
        assert_eq!(result[1], encode_op(4, 2)); // 2S
        assert_eq!(result[2], encode_op(0, 10)); // 10M
    }

    #[test]
    fn test_upgrade_clipping_raw_from_start_hard_and_soft() {
        // 2H5S10M: clip_amount=4, existing_hard=2, existing_soft=5
        // length_to_upgrade = min(5, 4-2) = 2
        // new_hard = 2+2 = 4, remaining_soft = 5-2 = 3
        let cigar = &[encode_op(5, 2), encode_op(4, 5), encode_op(0, 10)]; // 2H5S10M
        let (result, ref_consumed) = upgrade_clipping_raw(cigar, 4, true, encode_op);
        assert_eq!(ref_consumed, 0);
        assert_eq!(result.len(), 3); // 4H, 3S, 10M
        assert_eq!(result[0], encode_op(5, 4)); // 4H
        assert_eq!(result[1], encode_op(4, 3)); // 3S
        assert_eq!(result[2], encode_op(0, 10)); // 10M
    }

    #[test]
    fn test_upgrade_clipping_raw_from_start_all_soft_to_hard() {
        // 5S10M: clip_amount=5, upgrade all 5S to H
        let cigar = &[encode_op(4, 5), encode_op(0, 10)]; // 5S10M
        let (result, ref_consumed) = upgrade_clipping_raw(cigar, 5, true, encode_op);
        assert_eq!(ref_consumed, 0);
        assert_eq!(result.len(), 2); // 5H, 10M (no remaining S)
        assert_eq!(result[0], encode_op(5, 5)); // 5H
        assert_eq!(result[1], encode_op(0, 10)); // 10M
    }

    #[test]
    fn test_upgrade_clipping_raw_from_start_clip_amount_equals_hard() {
        // 3H5S10M: clip_amount=3, existing_hard=3
        // length_to_upgrade = min(5, 3-3) = 0
        // new_hard=3, remaining_soft=5
        let cigar = &[encode_op(5, 3), encode_op(4, 5), encode_op(0, 10)]; // 3H5S10M
        let (result, ref_consumed) = upgrade_clipping_raw(cigar, 3, true, encode_op);
        assert_eq!(ref_consumed, 0);
        assert_eq!(result.len(), 3); // 3H, 5S, 10M
        assert_eq!(result[0], encode_op(5, 3));
        assert_eq!(result[1], encode_op(4, 5));
        assert_eq!(result[2], encode_op(0, 10));
    }

    #[test]
    fn test_upgrade_clipping_raw_from_end_soft_only() {
        // 10M5S: clip_amount=3, upgrade 3S at end to H
        let cigar = &[encode_op(0, 10), encode_op(4, 5)]; // 10M5S
        let (result, ref_consumed) = upgrade_clipping_raw(cigar, 3, false, encode_op);
        assert_eq!(ref_consumed, 0);
        assert_eq!(result.len(), 3); // 10M, 2S, 3H
        assert_eq!(result[0], encode_op(0, 10)); // 10M
        assert_eq!(result[1], encode_op(4, 2)); // 2S
        assert_eq!(result[2], encode_op(5, 3)); // 3H
    }

    #[test]
    fn test_upgrade_clipping_raw_from_end_hard_and_soft() {
        // 10M5S2H: clip_amount=5, existing_hard=2, existing_soft=5
        // length_to_upgrade = min(5, 5-2) = 3
        // new_hard = 2+3 = 5, remaining_soft = 5-3 = 2
        let cigar = &[encode_op(0, 10), encode_op(4, 5), encode_op(5, 2)]; // 10M5S2H
        let (result, ref_consumed) = upgrade_clipping_raw(cigar, 5, false, encode_op);
        assert_eq!(ref_consumed, 0);
        assert_eq!(result.len(), 3); // 10M, 2S, 5H
        assert_eq!(result[0], encode_op(0, 10));
        assert_eq!(result[1], encode_op(4, 2));
        assert_eq!(result[2], encode_op(5, 5));
    }

    #[test]
    fn test_upgrade_clipping_raw_from_end_all_soft_to_hard() {
        // 10M5S: clip_amount=5, upgrade all soft
        let cigar = &[encode_op(0, 10), encode_op(4, 5)]; // 10M5S
        let (result, ref_consumed) = upgrade_clipping_raw(cigar, 5, false, encode_op);
        assert_eq!(ref_consumed, 0);
        assert_eq!(result.len(), 2); // 10M, 5H (no remaining S)
        assert_eq!(result[0], encode_op(0, 10));
        assert_eq!(result[1], encode_op(5, 5));
    }

    // ========================================================================
    // clip_cigar_ops_raw tests
    // ========================================================================

    #[test]
    fn test_clip_cigar_ops_raw_zero_clip() {
        let cigar = &[encode_op(0, 10)]; // 10M
        let (result, ref_consumed) = clip_cigar_ops_raw(cigar, 0, true);
        assert_eq!(ref_consumed, 0);
        assert_eq!(result, cigar);
    }

    #[test]
    fn test_clip_cigar_ops_raw_empty_cigar() {
        let (result, ref_consumed) = clip_cigar_ops_raw(&[], 5, true);
        assert_eq!(ref_consumed, 0);
        assert!(result.is_empty());
    }

    #[test]
    fn test_clip_cigar_ops_raw_upgrade_path_from_start() {
        // 5S10M: clip_amount=3, existing_clip=5, 3<=5 => upgrade path
        let cigar = &[encode_op(4, 5), encode_op(0, 10)]; // 5S10M
        let (result, ref_consumed) = clip_cigar_ops_raw(cigar, 3, true);
        assert_eq!(ref_consumed, 0);
        assert_eq!(result[0], encode_op(5, 3)); // 3H
        assert_eq!(result[1], encode_op(4, 2)); // 2S
        assert_eq!(result[2], encode_op(0, 10)); // 10M
    }

    #[test]
    fn test_clip_cigar_ops_raw_upgrade_path_from_end() {
        // 10M5S: clip_amount=3, existing_clip=5, 3<=5 => upgrade path
        let cigar = &[encode_op(0, 10), encode_op(4, 5)]; // 10M5S
        let (result, ref_consumed) = clip_cigar_ops_raw(cigar, 3, false);
        assert_eq!(ref_consumed, 0);
        assert_eq!(result[0], encode_op(0, 10)); // 10M
        assert_eq!(result[1], encode_op(4, 2)); // 2S
        assert_eq!(result[2], encode_op(5, 3)); // 3H
    }

    #[test]
    fn test_clip_cigar_ops_raw_alignment_clip_from_start() {
        // 10M: clip_amount=3, no existing clips => clips 3 bases from alignment
        let cigar = &[encode_op(0, 10)]; // 10M
        let (result, ref_consumed) = clip_cigar_ops_raw(cigar, 3, true);
        assert_eq!(ref_consumed, 3);
        assert_eq!(result.len(), 2); // 3H, 7M
        assert_eq!(result[0], encode_op(5, 3)); // 3H
        assert_eq!(result[1], encode_op(0, 7)); // 7M
    }

    #[test]
    fn test_clip_cigar_ops_raw_alignment_clip_from_end() {
        // 10M: clip_amount=3, no existing clips => clips 3 bases from alignment end
        let cigar = &[encode_op(0, 10)]; // 10M
        let (result, ref_consumed) = clip_cigar_ops_raw(cigar, 3, false);
        assert_eq!(ref_consumed, 0); // end-clipping doesn't consume ref for position adjustment
        assert_eq!(result.len(), 2); // 7M, 3H
        assert_eq!(result[0], encode_op(0, 7)); // 7M
        assert_eq!(result[1], encode_op(5, 3)); // 3H
    }

    #[test]
    fn test_clip_cigar_ops_raw_clip_past_existing_from_start() {
        // 2S10M: clip_amount=5, existing_clip=2, alignment_clip=3
        let cigar = &[encode_op(4, 2), encode_op(0, 10)]; // 2S10M
        let (result, ref_consumed) = clip_cigar_ops_raw(cigar, 5, true);
        assert_eq!(ref_consumed, 3);
        assert_eq!(result.len(), 2); // 5H, 7M
        assert_eq!(result[0], encode_op(5, 5)); // 5H (2S + 3 from alignment)
        assert_eq!(result[1], encode_op(0, 7)); // 7M
    }

    #[test]
    fn test_clip_cigar_ops_raw_clip_past_existing_from_end() {
        // 10M2S: clip_amount=5, existing_clip=2, alignment_clip=3
        let cigar = &[encode_op(0, 10), encode_op(4, 2)]; // 10M2S
        let (result, ref_consumed) = clip_cigar_ops_raw(cigar, 5, false);
        assert_eq!(ref_consumed, 0);
        assert_eq!(result.len(), 2); // 7M, 5H
        assert_eq!(result[0], encode_op(0, 7)); // 7M
        assert_eq!(result[1], encode_op(5, 5)); // 5H (2S + 3 from alignment)
    }

    #[test]
    fn test_clip_cigar_ops_raw_with_insertion_at_boundary_start() {
        // 10M3I5M: clip 10 from start should clip through 10M and consume entire 3I
        let cigar = &[encode_op(0, 10), encode_op(1, 3), encode_op(0, 5)]; // 10M3I5M
        let (result, ref_consumed) = clip_cigar_ops_raw(cigar, 10, true);
        // Clips 10M (10 read bases, 10 ref bases), then hits I:
        // After 10M, read_bases_clipped=10 == clip_amount=10, done
        assert_eq!(ref_consumed, 10);
        assert_eq!(result.len(), 3); // 10H, 3I, 5M
        assert_eq!(result[0], encode_op(5, 10)); // 10H
        assert_eq!(result[1], encode_op(1, 3)); // 3I preserved
        assert_eq!(result[2], encode_op(0, 5)); // 5M preserved
    }

    #[test]
    fn test_clip_cigar_ops_raw_with_deletion_at_boundary_start() {
        // 5M2D10M: clip 5 from start should clip 5M and skip the 2D
        let cigar = &[encode_op(0, 5), encode_op(2, 2), encode_op(0, 10)]; // 5M2D10M
        let (result, ref_consumed) = clip_cigar_ops_raw(cigar, 5, true);
        // Clips 5M (5 read + 5 ref), then hits deletion at boundary (read_bases_clipped==5),
        // skip 2D (ref_consumed += 2)
        assert_eq!(ref_consumed, 7); // 5 from M + 2 from D
        assert_eq!(result.len(), 2); // 5H, 10M
        assert_eq!(result[0], encode_op(5, 5)); // 5H
        assert_eq!(result[1], encode_op(0, 10)); // 10M
    }

    #[test]
    fn test_clip_cigar_ops_raw_with_deletion_at_boundary_end() {
        // 10M2D5M: clip 5 from end should clip 5M and skip the 2D
        let cigar = &[encode_op(0, 10), encode_op(2, 2), encode_op(0, 5)]; // 10M2D5M
        let (result, ref_consumed) = clip_cigar_ops_raw(cigar, 5, false);
        // Clips 5M from end (5 read bases), then skip 2D adjacent
        assert_eq!(ref_consumed, 0);
        assert_eq!(result.len(), 2); // 10M, 5H
        assert_eq!(result[0], encode_op(0, 10)); // 10M
        assert_eq!(result[1], encode_op(5, 5)); // 5H
    }

    #[test]
    fn test_clip_cigar_ops_raw_split_match_from_start() {
        // 10M: clip 4, splits the M
        let cigar = &[encode_op(0, 10)]; // 10M
        let (result, ref_consumed) = clip_cigar_ops_raw(cigar, 4, true);
        assert_eq!(ref_consumed, 4);
        assert_eq!(result.len(), 2); // 4H, 6M
        assert_eq!(result[0], encode_op(5, 4));
        assert_eq!(result[1], encode_op(0, 6));
    }

    #[test]
    fn test_clip_cigar_ops_raw_split_match_from_end() {
        // 10M: clip 4, splits the M from end
        let cigar = &[encode_op(0, 10)]; // 10M
        let (result, ref_consumed) = clip_cigar_ops_raw(cigar, 4, false);
        assert_eq!(ref_consumed, 0);
        assert_eq!(result.len(), 2); // 6M, 4H
        assert_eq!(result[0], encode_op(0, 6));
        assert_eq!(result[1], encode_op(5, 4));
    }

    #[test]
    fn test_clip_cigar_ops_raw_insertion_consumed_at_boundary_start() {
        // 5M3I5M: clip 6 from start
        // First 5M: 5 read bases clipped, 5 ref bases consumed
        // Then 3I: need 1 more, but insertion is consumed entirely
        // read_bases_clipped = 5 + 3 = 8 > 6 (entire I consumed at boundary)
        let cigar = &[encode_op(0, 5), encode_op(1, 3), encode_op(0, 5)]; // 5M3I5M
        let (result, ref_consumed) = clip_cigar_ops_raw(cigar, 6, true);
        // 5M consumed (5 read, 5 ref), then 3I consumed entirely (3 read)
        // total read bases = 8, ref_consumed = 5
        assert_eq!(ref_consumed, 5);
        assert_eq!(result.len(), 2); // 8H, 5M
        assert_eq!(result[0], encode_op(5, 8)); // 8H
        assert_eq!(result[1], encode_op(0, 5)); // 5M
    }

    #[test]
    fn test_clip_cigar_ops_raw_with_eq_and_x_ops() {
        // 5=3X: clip 4 from start splits the = op, consuming 4 ref bases
        let cigar = &[encode_op(7, 5), encode_op(8, 3)]; // 5= 3X
        let (result, ref_consumed) = clip_cigar_ops_raw(cigar, 4, true);
        assert_eq!(ref_consumed, 4);
        assert_eq!(result.len(), 3); // 4H, 1=, 3X
        assert_eq!(result[0], encode_op(5, 4));
        assert_eq!(result[1], encode_op(7, 1));
        assert_eq!(result[2], encode_op(8, 3));
    }

    #[test]
    fn test_clip_cigar_ops_raw_clip_entire_alignment_from_start() {
        // 10M: clip all 10 bases from start
        let cigar = &[encode_op(0, 10)]; // 10M
        let (result, ref_consumed) = clip_cigar_ops_raw(cigar, 10, true);
        assert_eq!(ref_consumed, 10);
        assert_eq!(result.len(), 1); // 10H only
        assert_eq!(result[0], encode_op(5, 10));
    }

    #[test]
    fn test_clip_cigar_ops_raw_clip_entire_alignment_from_end() {
        // 10M: clip all 10 bases from end
        let cigar = &[encode_op(0, 10)]; // 10M
        let (result, ref_consumed) = clip_cigar_ops_raw(cigar, 10, false);
        assert_eq!(ref_consumed, 0);
        assert_eq!(result.len(), 1); // 10H only
        assert_eq!(result[0], encode_op(5, 10));
    }

    #[test]
    fn test_clip_cigar_ops_raw_complex_cigar_start() {
        // 3S10M2I5M4S: clip_amount=8
        // existing_clip = 3 (3S), alignment_clip = 8 - 3 = 5
        // After stripping S: 10M2I5M4S
        // Clip 5 from 10M: split into 5 clipped + 5M remaining, ref=5
        let cigar = &[
            encode_op(4, 3),  // 3S
            encode_op(0, 10), // 10M
            encode_op(1, 2),  // 2I
            encode_op(0, 5),  // 5M
            encode_op(4, 4),  // 4S
        ];
        let (result, ref_consumed) = clip_cigar_ops_raw(cigar, 8, true);
        assert_eq!(ref_consumed, 5);
        // Result: 8H (3S+5 from alignment), 5M, 2I, 5M, 4S
        assert_eq!(result[0], encode_op(5, 8)); // 8H
        assert_eq!(result[1], encode_op(0, 5)); // 5M (remaining from 10M)
        assert_eq!(result[2], encode_op(1, 2)); // 2I
        assert_eq!(result[3], encode_op(0, 5)); // 5M
        assert_eq!(result[4], encode_op(4, 4)); // 4S (trailing)
    }

    #[test]
    fn test_clip_cigar_ops_raw_complex_cigar_end() {
        // 3S10M2I5M4S: clip_amount=8
        // existing_clip at end = 4 (4S), alignment_clip = 8 - 4 = 4
        // After stripping trailing S: 3S10M2I5M
        // Clip 4 from end of 5M: split into 1M remaining + 4 clipped
        let cigar = &[
            encode_op(4, 3),  // 3S
            encode_op(0, 10), // 10M
            encode_op(1, 2),  // 2I
            encode_op(0, 5),  // 5M
            encode_op(4, 4),  // 4S
        ];
        let (result, ref_consumed) = clip_cigar_ops_raw(cigar, 8, false);
        assert_eq!(ref_consumed, 0);
        // Result: 3S, 10M, 2I, 1M, 8H (4S+4 from alignment)
        assert_eq!(result[0], encode_op(4, 3)); // 3S (leading)
        assert_eq!(result[1], encode_op(0, 10)); // 10M
        assert_eq!(result[2], encode_op(1, 2)); // 2I
        assert_eq!(result[3], encode_op(0, 1)); // 1M (remaining from 5M)
        assert_eq!(result[4], encode_op(5, 8)); // 8H (4S+4)
    }

    // ========================================================================
    // is_fr_pair_raw tests
    // ========================================================================

    #[test]
    fn test_is_fr_pair_raw_not_paired() {
        // Not paired => false
        let rec = make_bam_bytes(0, 100, 0, b"rea", &[encode_op(0, 10)], 10, 0, 200, &[]);
        assert!(!is_fr_pair_raw(&rec));
    }

    #[test]
    fn test_is_fr_pair_raw_unmapped() {
        // Paired but unmapped => false
        let rec = make_bam_bytes(
            0,
            100,
            flags::PAIRED | flags::UNMAPPED,
            b"rea",
            &[encode_op(0, 10)],
            10,
            0,
            200,
            &[],
        );
        assert!(!is_fr_pair_raw(&rec));
    }

    #[test]
    fn test_is_fr_pair_raw_mate_unmapped() {
        // Paired, mapped but mate unmapped => false
        let rec = make_bam_bytes(
            0,
            100,
            flags::PAIRED | flags::MATE_UNMAPPED,
            b"rea",
            &[encode_op(0, 10)],
            10,
            -1,
            -1,
            &[],
        );
        assert!(!is_fr_pair_raw(&rec));
    }

    #[test]
    fn test_is_fr_pair_raw_different_references() {
        // Paired, both mapped, but different references => false
        let rec = make_bam_bytes(
            0,
            100,
            flags::PAIRED | flags::MATE_REVERSE,
            b"rea",
            &[encode_op(0, 10)],
            10,
            1,
            200,
            &[],
        );
        assert!(!is_fr_pair_raw(&rec));
    }

    #[test]
    fn test_is_fr_pair_raw_same_strand_ff() {
        // Paired, same reference, but both forward (FF) => false
        let rec =
            make_bam_bytes(0, 100, flags::PAIRED, b"rea", &[encode_op(0, 10)], 10, 0, 200, &[]);
        assert!(!is_fr_pair_raw(&rec));
    }

    #[test]
    fn test_is_fr_pair_raw_same_strand_rr() {
        // Paired, same reference, both reverse (RR) => false
        let rec = make_bam_bytes(
            0,
            100,
            flags::PAIRED | flags::REVERSE | flags::MATE_REVERSE,
            b"rea",
            &[encode_op(0, 10)],
            10,
            0,
            200,
            &[],
        );
        assert!(!is_fr_pair_raw(&rec));
    }

    #[test]
    fn test_is_fr_pair_raw_fr_positive_strand_read() {
        // FR pair: this read is forward, mate is reverse, on same reference
        // positive_five_prime = alignment_start = 101
        // negative_five_prime = alignment_start + insert_size = 101 + 200 = 301
        // 101 < 301 => FR => true
        let rec = make_bam_bytes_with_tlen(
            0,
            100,
            flags::PAIRED | flags::MATE_REVERSE,
            b"rea",
            &[encode_op(0, 10)],
            10,
            0,
            200,
            200,
            &[],
        );
        assert!(is_fr_pair_raw(&rec));
    }

    #[test]
    fn test_is_fr_pair_raw_fr_negative_strand_read() {
        // FR pair: this read is reverse, mate is forward
        // positive_five_prime = mate_start = 101
        // negative_five_prime = alignment_end = 101 + 10 - 1 = 110
        // Since mate at 101 < end at 110, this is FR => true
        let rec = make_bam_bytes_with_tlen(
            0,
            100,
            flags::PAIRED | flags::REVERSE,
            b"rea",
            &[encode_op(0, 10)],
            10,
            0,
            100,
            -10,
            &[],
        );
        assert!(is_fr_pair_raw(&rec));
    }

    #[test]
    fn test_is_fr_pair_raw_rf_orientation() {
        // RF pair: this read is forward, mate is reverse, but mate is upstream
        // Read at pos 200, mate at pos 100
        // positive_five_prime = alignment_start = 201
        // negative_five_prime = alignment_start + insert_size = 201 + (-100) = 101
        // 201 > 101 => NOT FR (it's RF) => false
        let rec = make_bam_bytes_with_tlen(
            0,
            200,
            flags::PAIRED | flags::MATE_REVERSE,
            b"rea",
            &[encode_op(0, 10)],
            10,
            0,
            100,
            -100,
            &[],
        );
        assert!(!is_fr_pair_raw(&rec));
    }

    // ========================================================================
    // compute_bases_past_ref_pos tests
    // ========================================================================

    #[test]
    fn test_compute_bases_past_ref_pos_simple_match() {
        // 10M starting at ref pos 100 (1-based)
        // target_ref_pos = 105: should find read_pos at offset 5
        let cigar = &[encode_op(0, 10)]; // 10M
        let result = compute_bases_past_ref_pos(&[], cigar, 100, 105);
        assert_eq!(result, 6); // 1-based: read pos 6 at ref pos 105
    }

    #[test]
    fn test_compute_bases_past_ref_pos_at_start() {
        // 10M starting at ref pos 100
        // target_ref_pos = 100: first position
        let cigar = &[encode_op(0, 10)];
        let result = compute_bases_past_ref_pos(&[], cigar, 100, 100);
        assert_eq!(result, 1);
    }

    #[test]
    fn test_compute_bases_past_ref_pos_at_end() {
        // 10M starting at ref pos 100
        // target_ref_pos = 109: last position
        let cigar = &[encode_op(0, 10)];
        let result = compute_bases_past_ref_pos(&[], cigar, 100, 109);
        assert_eq!(result, 10);
    }

    #[test]
    fn test_compute_bases_past_ref_pos_past_alignment() {
        // 10M starting at ref pos 100
        // target_ref_pos = 110: beyond alignment
        let cigar = &[encode_op(0, 10)];
        let result = compute_bases_past_ref_pos(&[], cigar, 100, 110);
        assert_eq!(result, 0); // outside alignment
    }

    #[test]
    fn test_compute_bases_past_ref_pos_with_insertion() {
        // 5M3I5M: insertion adds 3 query bases without consuming reference
        // At ref 100: 5M covers ref 100-104, 3I adds 3 query bases,
        // then 5M covers ref 105-109
        // target=107: in second 5M, offset 2 from ref 105
        // query pos = 5 (from first M) + 3 (from I) + 3 (offset in second M) = 11
        let cigar = &[encode_op(0, 5), encode_op(1, 3), encode_op(0, 5)]; // 5M3I5M
        let result = compute_bases_past_ref_pos(&[], cigar, 100, 107);
        assert_eq!(result, 11);
    }

    #[test]
    fn test_compute_bases_past_ref_pos_in_deletion() {
        // 5M3D5M: deletion spans ref 105-107 without consuming query
        // target=106: falls in the deletion
        let cigar = &[encode_op(0, 5), encode_op(2, 3), encode_op(0, 5)]; // 5M3D5M
        let result = compute_bases_past_ref_pos(&[], cigar, 100, 106);
        assert_eq!(result, 0); // position is in a deletion
    }

    #[test]
    fn test_compute_bases_past_ref_pos_with_soft_clip() {
        // 3S10M: soft clip consumes 3 query bases but no reference
        // Alignment starts at ref 100, so 10M covers ref 100-109
        // target=102: offset 2 in the M, query pos = 3 (from S) + 3 = 6
        let cigar = &[encode_op(4, 3), encode_op(0, 10)]; // 3S10M
        let result = compute_bases_past_ref_pos(&[], cigar, 100, 102);
        assert_eq!(result, 6);
    }

    // ========================================================================
    // compute_bases_before_ref_pos tests
    // ========================================================================

    #[test]
    fn test_compute_bases_before_ref_pos_simple_match() {
        // 10M starting at ref pos 100
        // target_ref_pos = 105: read_pos increments to 6, but returns 6-1=5
        let cigar = &[encode_op(0, 10)];
        let result = compute_bases_before_ref_pos(&[], cigar, 100, 105);
        assert_eq!(result, 5);
    }

    #[test]
    fn test_compute_bases_before_ref_pos_at_start() {
        // 10M starting at ref pos 100
        // target_ref_pos = 100: first position, read_pos=1, returns 1-1=0
        let cigar = &[encode_op(0, 10)];
        let result = compute_bases_before_ref_pos(&[], cigar, 100, 100);
        assert_eq!(result, 0);
    }

    #[test]
    fn test_compute_bases_before_ref_pos_past_alignment() {
        // 10M starting at ref pos 100
        // target_ref_pos = 110: beyond alignment
        let cigar = &[encode_op(0, 10)];
        let result = compute_bases_before_ref_pos(&[], cigar, 100, 110);
        assert_eq!(result, 0);
    }

    #[test]
    fn test_compute_bases_before_ref_pos_with_insertion() {
        // 5M3I5M: at ref 107 (in second M block)
        // query consumed: 5(M) + 3(I) + 3(into second M) = 11, returns 11-1=10
        let cigar = &[encode_op(0, 5), encode_op(1, 3), encode_op(0, 5)]; // 5M3I5M
        let result = compute_bases_before_ref_pos(&[], cigar, 100, 107);
        assert_eq!(result, 10);
    }

    #[test]
    fn test_compute_bases_before_ref_pos_in_deletion() {
        // 5M3D5M: deletion at ref 105-107
        // target=106: falls in deletion
        let cigar = &[encode_op(0, 5), encode_op(2, 3), encode_op(0, 5)]; // 5M3D5M
        let result = compute_bases_before_ref_pos(&[], cigar, 100, 106);
        assert_eq!(result, 0);
    }

    #[test]
    fn test_compute_bases_before_ref_pos_with_soft_clip() {
        // 3S10M: soft clip consumes 3 query bases
        // At ref 102: query = 3(S) + 3(into M) = 6, returns 6-1=5
        let cigar = &[encode_op(4, 3), encode_op(0, 10)]; // 3S10M
        let result = compute_bases_before_ref_pos(&[], cigar, 100, 102);
        assert_eq!(result, 5);
    }

    // ========================================================================
    // num_bases_extending_past_mate_raw tests
    // ========================================================================

    #[test]
    fn test_num_bases_extending_past_mate_raw_not_paired() {
        // Not paired => 0
        let rec = make_bam_bytes(0, 100, 0, b"rea", &[encode_op(0, 10)], 10, 0, 200, &[]);
        assert_eq!(num_bases_extending_past_mate_raw(&rec), 0);
    }

    #[test]
    fn test_num_bases_extending_past_mate_raw_unmapped() {
        // Paired but unmapped => 0
        let rec = make_bam_bytes(
            0,
            100,
            flags::PAIRED | flags::UNMAPPED | flags::MATE_REVERSE,
            b"rea",
            &[encode_op(0, 10)],
            10,
            0,
            200,
            &[],
        );
        assert_eq!(num_bases_extending_past_mate_raw(&rec), 0);
    }

    #[test]
    fn test_num_bases_extending_past_mate_raw_mate_unmapped() {
        // Paired, mapped but mate unmapped => 0
        let rec = make_bam_bytes(
            0,
            100,
            flags::PAIRED | flags::MATE_UNMAPPED | flags::MATE_REVERSE,
            b"rea",
            &[encode_op(0, 10)],
            10,
            -1,
            -1,
            &[],
        );
        assert_eq!(num_bases_extending_past_mate_raw(&rec), 0);
    }

    #[test]
    fn test_num_bases_extending_past_mate_raw_same_strand() {
        // Both same strand => 0
        let mut aux = Vec::new();
        aux.extend_from_slice(b"MCZ10M\x00");
        let rec = make_bam_bytes(
            0,
            100,
            flags::PAIRED, // both forward
            b"rea",
            &[encode_op(0, 10)],
            10,
            0,
            200,
            &aux,
        );
        assert_eq!(num_bases_extending_past_mate_raw(&rec), 0);
    }

    #[test]
    fn test_num_bases_extending_past_mate_raw_different_references() {
        // Different references => 0
        let mut aux = Vec::new();
        aux.extend_from_slice(b"MCZ10M\x00");
        let rec = make_bam_bytes(
            0,
            100,
            flags::PAIRED | flags::MATE_REVERSE,
            b"rea",
            &[encode_op(0, 10)],
            10,
            1,
            200,
            &aux,
        );
        assert_eq!(num_bases_extending_past_mate_raw(&rec), 0);
    }

    #[test]
    fn test_num_bases_extending_past_mate_raw_no_mc_tag() {
        // Paired FR but no MC tag => 0
        let rec = make_bam_bytes(
            0,
            100,
            flags::PAIRED | flags::MATE_REVERSE,
            b"rea",
            &[encode_op(0, 10)],
            10,
            0,
            200,
            &[],
        );
        assert_eq!(num_bases_extending_past_mate_raw(&rec), 0);
    }

    #[test]
    fn test_num_bases_extending_past_mate_raw_positive_strand_overlap() {
        // Positive strand read extending past reverse mate's unclipped end
        // Read: forward at pos 100, 20M (alignment_end = 101 + 20 - 1 = 120, 1-based)
        // Mate: reverse at pos 105, MC=10M (mate unclipped_end = 106 + 10 - 1 = 115, 1-based)
        // alignment_end (120) >= mate_ue (115), so need to compute bases past ref 115
        // At ref 115 (1-based), starting from alignment_start=101:
        // offset in 20M = 115 - 101 + 1 = 15 -> read_pos 15
        // read_length = 20
        // result = 20 - 15 = 5
        let mut aux = Vec::new();
        aux.extend_from_slice(b"MCZ10M\x00");
        let rec = make_bam_bytes(
            0,
            100, // 0-based pos
            flags::PAIRED | flags::MATE_REVERSE,
            b"rea",
            &[encode_op(0, 20)],
            20,
            0,
            105, // 0-based mate pos
            &aux,
        );
        assert_eq!(num_bases_extending_past_mate_raw(&rec), 5);
    }

    #[test]
    fn test_num_bases_extending_past_mate_raw_positive_strand_no_overlap() {
        // Positive strand read NOT extending past reverse mate's unclipped end
        // Read: forward at pos 100, 10M (alignment_end = 101 + 10 - 1 = 110)
        // Mate: reverse at pos 200, MC=10M (mate unclipped_end = 201 + 10 - 1 = 210)
        // alignment_end (110) < mate_ue (210), check trailing soft clips
        // No trailing soft clips => trailing_sc.saturating_sub(gap) = 0
        let mut aux = Vec::new();
        aux.extend_from_slice(b"MCZ10M\x00");
        let rec = make_bam_bytes(
            0,
            100,
            flags::PAIRED | flags::MATE_REVERSE,
            b"rea",
            &[encode_op(0, 10)],
            10,
            0,
            200,
            &aux,
        );
        assert_eq!(num_bases_extending_past_mate_raw(&rec), 0);
    }

    #[test]
    fn test_num_bases_extending_past_mate_raw_negative_strand_overlap() {
        // Negative strand read extending before forward mate's unclipped start
        // Read: reverse at pos 100, 20M
        // Mate: forward at pos 105, MC=10M (mate unclipped_start = 106)
        // this_pos (101) <= mate_us (106), so compute bases before ref pos 106
        // 20M from pos 101: at ref 106, read_pos = 6, returns 6-1 = 5
        let mut aux = Vec::new();
        aux.extend_from_slice(b"MCZ10M\x00");
        let rec = make_bam_bytes(
            0,
            100, // 0-based pos
            flags::PAIRED | flags::REVERSE,
            b"rea",
            &[encode_op(0, 20)],
            20,
            0,
            105, // 0-based mate pos
            &aux,
        );
        assert_eq!(num_bases_extending_past_mate_raw(&rec), 5);
    }

    #[test]
    fn test_num_bases_extending_past_mate_raw_negative_strand_no_overlap() {
        // Negative strand read NOT extending before forward mate
        // Read: reverse at pos 200, 10M
        // Mate: forward at pos 100, MC=10M (mate unclipped_start = 101)
        // this_pos (201) > mate_us (101), check leading soft clips
        // No leading soft clips, gap = 201 - 101 = 100, 0.saturating_sub(100) = 0
        let mut aux = Vec::new();
        aux.extend_from_slice(b"MCZ10M\x00");
        let rec = make_bam_bytes(
            0,
            200,
            flags::PAIRED | flags::REVERSE,
            b"rea",
            &[encode_op(0, 10)],
            10,
            0,
            100,
            &aux,
        );
        assert_eq!(num_bases_extending_past_mate_raw(&rec), 0);
    }

    #[test]
    fn test_num_bases_extending_past_mate_raw_negative_strand_gap_with_soft_clip() {
        // Negative strand read with soft clip, this_pos > mate_us
        // Read: reverse at pos 110 (0-based), 3S10M (query_len=13)
        // Mate: forward at pos 105 (0-based), MC=10M (mate unclipped_start = 106)
        // this_pos = 111, mate_us = 106
        // this_pos > mate_us, so gap = 111 - 106 = 5
        // leading_soft_clip = 3, 3.saturating_sub(5) = 0
        let mut aux = Vec::new();
        aux.extend_from_slice(b"MCZ10M\x00");
        let rec = make_bam_bytes(
            0,
            110,
            flags::PAIRED | flags::REVERSE,
            b"rea",
            &[encode_op(4, 3), encode_op(0, 10)],
            13,
            0,
            105,
            &aux,
        );
        assert_eq!(num_bases_extending_past_mate_raw(&rec), 0);
    }

    #[test]
    fn test_num_bases_extending_past_mate_raw_positive_strand_gap_with_soft_clip() {
        // Positive strand read with trailing soft clip, alignment_end < mate_ue
        // Read: forward at pos 100 (0-based), 10M3S (query_len=13)
        // Mate: reverse at pos 200 (0-based), MC=10M (mate unclipped_end = 201+10-1=210)
        // alignment_end = 101+10-1 = 110, mate_ue = 210
        // 110 < 210, gap = 210 - 110 = 100
        // trailing_soft_clip = 3, 3.saturating_sub(100) = 0
        let mut aux = Vec::new();
        aux.extend_from_slice(b"MCZ10M\x00");
        let rec = make_bam_bytes(
            0,
            100,
            flags::PAIRED | flags::MATE_REVERSE,
            b"rea",
            &[encode_op(0, 10), encode_op(4, 3)],
            13,
            0,
            200,
            &aux,
        );
        assert_eq!(num_bases_extending_past_mate_raw(&rec), 0);
    }

    // ========================================================================
    // read_pos_at_ref_pos_raw tests
    // ========================================================================

    #[test]
    fn test_read_pos_at_ref_pos_raw_simple() {
        // 10M starting at position 100: ref pos 102 => query pos 3
        let cigar = &[encode_op(0, 10)];
        assert_eq!(read_pos_at_ref_pos_raw(cigar, 100, 102, false), Some(3));
    }

    #[test]
    fn test_read_pos_at_ref_pos_raw_before_alignment() {
        // ref_pos < alignment_start => None
        let cigar = &[encode_op(0, 10)];
        assert_eq!(read_pos_at_ref_pos_raw(cigar, 100, 99, false), None);
    }

    #[test]
    fn test_read_pos_at_ref_pos_raw_after_alignment() {
        // ref_pos after alignment => None
        let cigar = &[encode_op(0, 10)];
        assert_eq!(read_pos_at_ref_pos_raw(cigar, 100, 110, false), None);
    }

    #[test]
    fn test_read_pos_at_ref_pos_raw_in_deletion() {
        // 5M3D5M: deletion at ref 105-107
        // ref_pos=106 is in deletion => None (without return_last_base_if_deleted)
        let cigar = &[encode_op(0, 5), encode_op(2, 3), encode_op(0, 5)];
        assert_eq!(read_pos_at_ref_pos_raw(cigar, 100, 106, false), None);
    }

    #[test]
    fn test_read_pos_at_ref_pos_raw_in_deletion_return_last() {
        // 5M3D5M: deletion at ref 105-107
        // ref_pos=106 is in deletion, return_last_base_if_deleted=true => returns 5
        let cigar = &[encode_op(0, 5), encode_op(2, 3), encode_op(0, 5)];
        assert_eq!(read_pos_at_ref_pos_raw(cigar, 100, 106, true), Some(5));
    }

    #[test]
    fn test_read_pos_at_ref_pos_raw_with_insertion() {
        // 5M3I5M: insertion doesn't consume reference
        // ref 105 is first base of second 5M, query pos = 5 (from first M) + 3 (from I) + 1 = 9
        let cigar = &[encode_op(0, 5), encode_op(1, 3), encode_op(0, 5)];
        assert_eq!(read_pos_at_ref_pos_raw(cigar, 100, 105, false), Some(9));
    }

    #[test]
    fn test_read_pos_at_ref_pos_raw_with_soft_clip() {
        // 3S10M: soft clip consumes query but not reference
        // ref 100 = first ref base in 10M, query pos = 3 (from S) + 1 = 4
        let cigar = &[encode_op(4, 3), encode_op(0, 10)];
        assert_eq!(read_pos_at_ref_pos_raw(cigar, 100, 100, false), Some(4));
    }

    #[test]
    fn test_read_pos_at_ref_pos_raw_at_exact_start() {
        // 10M at position 100: query at ref 100 = 1
        let cigar = &[encode_op(0, 10)];
        assert_eq!(read_pos_at_ref_pos_raw(cigar, 100, 100, false), Some(1));
    }

    #[test]
    fn test_read_pos_at_ref_pos_raw_at_exact_end() {
        // 10M at position 100: query at ref 109 = 10
        let cigar = &[encode_op(0, 10)];
        assert_eq!(read_pos_at_ref_pos_raw(cigar, 100, 109, false), Some(10));
    }

    // ========================================================================
    // query_length_from_cigar tests
    // ========================================================================

    #[test]
    fn test_query_length_from_cigar_simple() {
        let cigar = &[encode_op(0, 10)]; // 10M
        assert_eq!(query_length_from_cigar(cigar), 10);
    }

    #[test]
    fn test_query_length_from_cigar_with_clips_and_insertion() {
        // 3S10M2I5M4S: query consuming = 3+10+2+5+4 = 24
        let cigar = &[
            encode_op(4, 3),  // 3S
            encode_op(0, 10), // 10M
            encode_op(1, 2),  // 2I
            encode_op(0, 5),  // 5M
            encode_op(4, 4),  // 4S
        ];
        assert_eq!(query_length_from_cigar(cigar), 24);
    }

    #[test]
    fn test_query_length_from_cigar_deletion_not_counted() {
        // 10M3D5M: deletion doesn't consume query
        let cigar = &[encode_op(0, 10), encode_op(2, 3), encode_op(0, 5)];
        assert_eq!(query_length_from_cigar(cigar), 15);
    }

    #[test]
    fn test_query_length_from_cigar_hard_clip_not_counted() {
        // 3H10M2H: hard clips don't consume query
        let cigar = &[encode_op(5, 3), encode_op(0, 10), encode_op(5, 2)];
        assert_eq!(query_length_from_cigar(cigar), 10);
    }

    #[test]
    fn test_query_length_from_cigar_eq_and_x() {
        // 5=3X: both consume query, total 8
        let cigar = &[encode_op(7, 5), encode_op(8, 3)];
        assert_eq!(query_length_from_cigar(cigar), 8);
    }

    #[test]
    fn test_query_length_from_cigar_empty() {
        assert_eq!(query_length_from_cigar(&[]), 0);
    }

    // ========================================================================
    // trailing_soft_clip_from_ops and leading_soft_clip_from_ops tests
    // ========================================================================

    #[test]
    fn test_trailing_soft_clip_from_ops_none() {
        let cigar = &[encode_op(0, 10)]; // 10M
        assert_eq!(trailing_soft_clip_from_ops(cigar), 0);
    }

    #[test]
    fn test_trailing_soft_clip_from_ops_with_soft() {
        let cigar = &[encode_op(0, 10), encode_op(4, 5)]; // 10M5S
        assert_eq!(trailing_soft_clip_from_ops(cigar), 5);
    }

    #[test]
    fn test_trailing_soft_clip_from_ops_with_hard_after_soft() {
        let cigar = &[encode_op(0, 10), encode_op(4, 5), encode_op(5, 3)]; // 10M5S3H
        assert_eq!(trailing_soft_clip_from_ops(cigar), 5);
    }

    #[test]
    fn test_leading_soft_clip_from_ops_none() {
        let cigar = &[encode_op(0, 10)]; // 10M
        assert_eq!(leading_soft_clip_from_ops(cigar), 0);
    }

    #[test]
    fn test_leading_soft_clip_from_ops_with_soft() {
        let cigar = &[encode_op(4, 5), encode_op(0, 10)]; // 5S10M
        assert_eq!(leading_soft_clip_from_ops(cigar), 5);
    }

    #[test]
    fn test_leading_soft_clip_from_ops_with_hard_before_soft() {
        let cigar = &[encode_op(5, 3), encode_op(4, 5), encode_op(0, 10)]; // 3H5S10M
        assert_eq!(leading_soft_clip_from_ops(cigar), 5);
    }

    // ========================================================================
    // UnmappedBamRecordBuilder additional tests
    // ========================================================================

    #[test]
    fn test_builder_default() {
        let builder = UnmappedBamRecordBuilder::default();
        // Default builder has empty buffer and is not sealed
        assert!(builder.buf.is_empty());
        assert!(!builder.sealed);
    }

    #[test]
    fn test_builder_with_capacity() {
        let builder = UnmappedBamRecordBuilder::with_capacity(1024);
        assert!(builder.buf.capacity() >= 1024);
    }

    #[test]
    fn test_builder_odd_length_sequence() {
        // Odd-length sequence to exercise the packing edge case
        let mut builder = UnmappedBamRecordBuilder::new();
        builder.build_record(b"read1", flags::UNMAPPED, b"ACG", &[30, 25, 20]);

        let bam = builder.as_bytes();
        assert_eq!(l_seq(bam), 3);
        let so = seq_offset(bam);
        assert_eq!(get_base(bam, so, 0), 1); // A
        assert_eq!(get_base(bam, so, 1), 2); // C
        assert_eq!(get_base(bam, so, 2), 4); // G
    }

    #[test]
    fn test_builder_long_sequence() {
        // Longer sequence to exercise multi-byte packing
        let mut builder = UnmappedBamRecordBuilder::new();
        let bases = b"ACGTACGTACGTACGT"; // 16 bases
        let quals = vec![30u8; 16];
        builder.build_record(b"r1", flags::UNMAPPED, bases, &quals);

        let bam = builder.as_bytes();
        assert_eq!(l_seq(bam), 16);
        let so = seq_offset(bam);
        for (i, &base) in bases.iter().enumerate() {
            let expected = match base {
                b'A' => 1,
                b'C' => 2,
                b'G' => 4,
                b'T' => 8,
                _ => 0,
            };
            assert_eq!(get_base(bam, so, i), expected, "base mismatch at position {i}");
        }
    }

    #[test]
    fn test_builder_multiple_records_via_write_with_block_size() {
        // Build two records into the same output buffer
        let mut builder = UnmappedBamRecordBuilder::new();
        let mut output = Vec::new();

        builder.build_record(b"r1", flags::UNMAPPED, b"ACGT", &[30, 30, 30, 30]);
        builder.append_string_tag(b"MI", b"1");
        builder.write_with_block_size(&mut output);

        builder.clear();
        builder.build_record(b"r2", flags::UNMAPPED, b"TGCA", &[25, 25, 25, 25]);
        builder.append_string_tag(b"MI", b"2");
        builder.write_with_block_size(&mut output);

        let records = ParsedBamRecord::parse_all(&output);
        assert_eq!(records.len(), 2);
        assert_eq!(records[0].name, b"r1");
        assert_eq!(records[0].bases, b"ACGT");
        assert_eq!(records[0].get_string_tag(b"MI").unwrap(), b"1");
        assert_eq!(records[1].name, b"r2");
        assert_eq!(records[1].bases, b"TGCA");
        assert_eq!(records[1].get_string_tag(b"MI").unwrap(), b"2");
    }

    // ========================================================================
    // clip_cigar_start_raw / clip_cigar_end_raw additional edge cases
    // ========================================================================

    #[test]
    fn test_clip_cigar_start_raw_with_existing_hard_clip() {
        // 3H10M: clip 5 from start
        // existing_hard=3, existing_soft=0, post_clip_ops=[10M]
        // alignment_clip = 5 (called from clip_cigar_ops_raw: 5 - 3 = 2)
        // But directly calling clip_cigar_start_raw with clip_amount=2:
        // clips 2 from 10M => 8M remains, ref_consumed=2
        // total_hard = 3 + 0 + 2 = 5
        let cigar = &[encode_op(5, 3), encode_op(0, 10)]; // 3H10M
        let (result, ref_consumed) = clip_cigar_start_raw(cigar, 2, encode_op);
        assert_eq!(ref_consumed, 2);
        assert_eq!(result.len(), 2); // 5H, 8M
        assert_eq!(result[0], encode_op(5, 5)); // 3 existing + 2 new
        assert_eq!(result[1], encode_op(0, 8));
    }

    #[test]
    fn test_clip_cigar_end_raw_with_existing_hard_clip() {
        // 10M3H: clip 2 from end
        // existing_hard=3, existing_soft=0, post_clip_ops=[10M]
        // clips 2 from end of 10M => 8M remains
        // total_hard = 3 + 0 + 2 = 5
        let cigar = &[encode_op(0, 10), encode_op(5, 3)]; // 10M3H
        let (result, ref_consumed) = clip_cigar_end_raw(cigar, 2, encode_op);
        assert_eq!(ref_consumed, 0);
        assert_eq!(result.len(), 2); // 8M, 5H
        assert_eq!(result[0], encode_op(0, 8));
        assert_eq!(result[1], encode_op(5, 5)); // 3 existing + 2 new
    }

    #[test]
    fn test_clip_cigar_start_raw_insertion_at_exact_boundary() {
        // 3I10M: clip 1 from start
        // The insertion has 3 read bases. clip_amount=1, but I consumes entire (3 bases)
        // total read bases clipped = 3, ref_consumed = 0
        let cigar = &[encode_op(1, 3), encode_op(0, 10)]; // 3I10M
        let (result, ref_consumed) = clip_cigar_start_raw(cigar, 1, encode_op);
        assert_eq!(ref_consumed, 0);
        // 3I consumed entirely => 3H, 10M
        assert_eq!(result.len(), 2);
        assert_eq!(result[0], encode_op(5, 3));
        assert_eq!(result[1], encode_op(0, 10));
    }

    #[test]
    fn test_clip_cigar_end_raw_insertion_at_exact_boundary() {
        // 10M3I: clip 1 from end
        // The insertion has 3 read bases. clip_amount=1, but I consumes entire (3 bases)
        let cigar = &[encode_op(0, 10), encode_op(1, 3)]; // 10M3I
        let (result, ref_consumed) = clip_cigar_end_raw(cigar, 1, encode_op);
        assert_eq!(ref_consumed, 0);
        // 3I consumed entirely => 10M, 3H
        assert_eq!(result.len(), 2);
        assert_eq!(result[0], encode_op(0, 10));
        assert_eq!(result[1], encode_op(5, 3));
    }

    // ========================================================================
    // alignment_start_from_raw / alignment_end_from_raw tests
    // ========================================================================

    #[test]
    fn test_alignment_start_from_raw_mapped() {
        let rec = make_bam_bytes(0, 100, 0, b"rea", &[encode_op(0, 10)], 10, -1, -1, &[]);
        assert_eq!(alignment_start_from_raw(&rec), Some(101));
    }

    #[test]
    fn test_alignment_start_from_raw_unmapped() {
        let rec = make_bam_bytes(-1, -1, flags::UNMAPPED, b"rea", &[], 0, -1, -1, &[]);
        assert_eq!(alignment_start_from_raw(&rec), None);
    }

    #[test]
    fn test_alignment_end_from_raw_mapped() {
        // 10M at pos 100 (0-based): end = 100 + 10 = 110 (exclusive)
        let rec = make_bam_bytes(0, 100, 0, b"rea", &[encode_op(0, 10)], 10, -1, -1, &[]);
        assert_eq!(alignment_end_from_raw(&rec), Some(110));
    }

    #[test]
    fn test_alignment_end_from_raw_unmapped() {
        let rec = make_bam_bytes(-1, -1, flags::UNMAPPED, b"rea", &[], 0, -1, -1, &[]);
        assert_eq!(alignment_end_from_raw(&rec), None);
    }

    #[test]
    fn test_alignment_end_from_raw_no_cigar() {
        // Mapped but no CIGAR => ref_len=0 => None
        let rec = make_bam_bytes(0, 100, 0, b"rea", &[], 0, -1, -1, &[]);
        assert_eq!(alignment_end_from_raw(&rec), None);
    }

    // ========================================================================
    // template_length tests
    // ========================================================================

    #[test]
    fn test_template_length_zero() {
        let rec = make_bam_bytes(0, 100, 0, b"rea", &[], 0, -1, -1, &[]);
        assert_eq!(template_length(&rec), 0);
    }

    #[test]
    fn test_template_length_positive() {
        let rec = make_bam_bytes_with_tlen(0, 100, 0, b"rea", &[], 0, 0, 200, 150, &[]);
        assert_eq!(template_length(&rec), 150);
    }

    #[test]
    fn test_template_length_negative() {
        let rec = make_bam_bytes_with_tlen(0, 200, 0, b"rea", &[], 0, 0, 100, -150, &[]);
        assert_eq!(template_length(&rec), -150);
    }

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

    // ========================================================================
    // extract_sequence / quality_scores_slice tests
    // ========================================================================

    #[test]
    fn test_extract_sequence() {
        let mut rec = make_bam_bytes(0, 0, 0, b"rea", &[], 4, -1, -1, &[]);
        let so = seq_offset(&rec);
        rec[so] = 0x12; // A, C
        rec[so + 1] = 0x48; // G, T
        let seq = extract_sequence(&rec);
        assert_eq!(seq, b"ACGT");
    }

    #[test]
    fn test_quality_scores_slice() {
        let mut rec = make_bam_bytes(0, 0, 0, b"rea", &[], 4, -1, -1, &[]);
        let qo = qual_offset(&rec);
        rec[qo] = 30;
        rec[qo + 1] = 25;
        rec[qo + 2] = 20;
        rec[qo + 3] = 15;
        let quals = quality_scores_slice(&rec);
        assert_eq!(quals, &[30, 25, 20, 15]);
    }

    #[test]
    fn test_quality_scores_slice_mut() {
        let mut rec = make_bam_bytes(0, 0, 0, b"rea", &[], 4, -1, -1, &[]);
        let quals = quality_scores_slice_mut(&mut rec);
        quals[0] = 30;
        quals[1] = 25;
        assert_eq!(quality_scores_slice(&rec), &[30, 25, 0, 0]);
    }

    // ========================================================================
    // set_base tests
    // ========================================================================

    #[test]
    fn test_set_base() {
        let mut rec = make_bam_bytes(0, 0, 0, b"rea", &[], 4, -1, -1, &[]);
        let so = seq_offset(&rec);
        // Start with AAAA (all 0x11 pairs)
        rec[so] = 0x11; // A=1, A=1
        rec[so + 1] = 0x11;

        // set_base takes ASCII bases, not codes
        // Set position 0 to T (code 8)
        set_base(&mut rec, so, 0, b'T');
        assert_eq!(get_base(&rec, so, 0), 8); // T
        assert_eq!(get_base(&rec, so, 1), 1); // A unchanged

        // Set position 1 to C (code 2)
        set_base(&mut rec, so, 1, b'C');
        assert_eq!(get_base(&rec, so, 0), 8); // T unchanged
        assert_eq!(get_base(&rec, so, 1), 2); // C

        // Set position 3 to G (code 4) - odd position
        set_base(&mut rec, so, 3, b'G');
        assert_eq!(get_base(&rec, so, 2), 1); // A unchanged
        assert_eq!(get_base(&rec, so, 3), 4); // G
    }

    // ========================================================================
    // Field accessor tests — documents the BAM binary format layout
    // ========================================================================

    #[test]
    fn test_field_accessors_roundtrip() {
        // Build a record with known values and verify every accessor
        let mut rec = make_bam_bytes_with_tlen(
            3,                                                     // tid
            200,                                                   // pos
            flags::PAIRED | flags::REVERSE | flags::FIRST_SEGMENT, // flag
            b"read1",                                              // name
            &[encode_op(0, 10)],                                   // 10M cigar
            6,                                                     // seq_len
            5,                                                     // mate_tid
            400,                                                   // mate_pos
            150,                                                   // tlen
            &[],                                                   // aux
        );
        rec[9] = 42; // mapq

        assert_eq!(ref_id(&rec), 3);
        assert_eq!(pos(&rec), 200);
        assert_eq!(l_read_name(&rec), 6); // "read1" + null
        assert_eq!(mapq(&rec), 42);
        assert_eq!(n_cigar_op(&rec), 1);
        assert_eq!(flags(&rec), flags::PAIRED | flags::REVERSE | flags::FIRST_SEGMENT);
        assert_eq!(l_seq(&rec), 6);
        assert_eq!(mate_ref_id(&rec), 5);
        assert_eq!(mate_pos(&rec), 400);
        assert_eq!(template_length(&rec), 150);
        assert_eq!(read_name(&rec), b"read1");
    }

    #[test]
    fn test_flags_unmapped_record() {
        let rec = make_bam_bytes(-1, -1, flags::UNMAPPED, b"r", &[], 0, -1, -1, &[]);
        assert_eq!(flags(&rec), flags::UNMAPPED);
        assert_eq!(ref_id(&rec), -1);
        assert_eq!(pos(&rec), -1);
    }

    // ========================================================================
    // aux_data_offset_from_record / aux_data_slice / seq_offset / qual_offset
    // ========================================================================

    #[test]
    fn test_aux_data_offset_from_record_basic() {
        let rec = make_bam_bytes(0, 0, 0, b"rd", &[encode_op(0, 4)], 4, -1, -1, &[]);
        let offset = aux_data_offset_from_record(&rec).unwrap();
        // 32 + l_read_name(3) + cigar(4) + seq_bytes(2) + qual(4) = 45
        assert_eq!(offset, 32 + 3 + 4 + 2 + 4);
    }

    #[test]
    fn test_aux_data_offset_from_record_too_short() {
        let short = vec![0u8; 10];
        assert!(aux_data_offset_from_record(&short).is_none());
    }

    #[test]
    fn test_aux_data_slice_with_tags() {
        let aux = b"MIZtest\0";
        let rec = make_bam_bytes(0, 0, 0, b"rd", &[], 4, -1, -1, aux);
        let slice = aux_data_slice(&rec);
        assert_eq!(slice, aux.as_slice());
    }

    #[test]
    fn test_aux_data_slice_empty() {
        let rec = make_bam_bytes(0, 0, 0, b"rd", &[], 4, -1, -1, &[]);
        let slice = aux_data_slice(&rec);
        assert!(slice.is_empty());
    }

    #[test]
    fn test_seq_offset_and_qual_offset() {
        let rec = make_bam_bytes(0, 0, 0, b"rd", &[encode_op(0, 4)], 4, -1, -1, &[]);
        let so = seq_offset(&rec);
        let qo = qual_offset(&rec);
        // seq_offset = 32 + l_read_name(3) + cigar(4) = 39
        assert_eq!(so, 32 + 3 + 4);
        // qual_offset = seq_offset + seq_bytes(2) = 41
        assert_eq!(qo, so + 2);
        // quality should be 4 bytes at that offset
        assert_eq!(&rec[qo..qo + 4], &[0, 0, 0, 0]);
    }

    // ========================================================================
    // alignment_end_from_raw — additional cases
    // ========================================================================

    #[test]
    fn test_alignment_end_from_raw_with_deletion() {
        // 5M2D5M at pos=0 → ref_len = 5+2+5 = 12
        let cigar = &[encode_op(0, 5), encode_op(2, 2), encode_op(0, 5)];
        let rec = make_bam_bytes(0, 0, 0, b"rd", cigar, 10, -1, -1, &[]);
        assert_eq!(alignment_end_from_raw(&rec), Some(12));
    }

    // ========================================================================
    // simplify_cigar_from_raw — additional complex patterns
    // ========================================================================

    #[test]
    fn test_simplify_cigar_from_raw_hard_clips() {
        use noodles::sam::alignment::record::cigar::op::Kind;
        // 5H10M5H → H becomes M, coalesced: 20M
        let cigar = &[encode_op(5, 5), encode_op(0, 10), encode_op(5, 5)];
        let result = simplify_cigar_from_raw(cigar);
        assert_eq!(result, vec![(Kind::Match, 20)]);
    }

    #[test]
    fn test_simplify_cigar_from_raw_mixed_operations() {
        use noodles::sam::alignment::record::cigar::op::Kind;
        // 2H3S5M2I3M1D4M4S1H → H,S,M coalesce to M; I,D stay
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
        // 5M3N5M → N (skip) preserved, not converted to M
        let cigar = &[encode_op(0, 5), encode_op(3, 3), encode_op(0, 5)];
        let result = simplify_cigar_from_raw(cigar);
        assert_eq!(result, vec![(Kind::Match, 5), (Kind::Skip, 3), (Kind::Match, 5)]);
    }

    #[test]
    fn test_simplify_cigar_from_raw_all_soft_clips() {
        use noodles::sam::alignment::record::cigar::op::Kind;
        // 10S → becomes 10M
        let cigar = &[encode_op(4, 10)];
        let result = simplify_cigar_from_raw(cigar);
        assert_eq!(result, vec![(Kind::Match, 10)]);
    }

    // ========================================================================
    // extract_sequence — additional edge cases
    // ========================================================================

    #[test]
    fn test_extract_sequence_odd_length() {
        // Odd-length seq: the last nibble in the packed byte is in the high nibble
        let mut rec = make_bam_bytes(0, 0, 0, b"rd", &[], 3, -1, -1, &[]);
        let so = seq_offset(&rec);
        rec[so] = 0x12; // A(1), C(2)
        rec[so + 1] = 0x40; // G(4), padding(0)
        let seq = extract_sequence(&rec);
        assert_eq!(seq, b"ACG");
    }

    #[test]
    fn test_extract_sequence_all_n() {
        let mut rec = make_bam_bytes(0, 0, 0, b"rd", &[], 4, -1, -1, &[]);
        let so = seq_offset(&rec);
        rec[so] = 0xFF; // N(15), N(15)
        rec[so + 1] = 0xFF; // N(15), N(15)
        let seq = extract_sequence(&rec);
        assert_eq!(seq, b"NNNN");
    }

    #[test]
    fn test_extract_sequence_single_base() {
        let mut rec = make_bam_bytes(0, 0, 0, b"rd", &[], 1, -1, -1, &[]);
        let so = seq_offset(&rec);
        rec[so] = 0x80; // T(8)
        let seq = extract_sequence(&rec);
        assert_eq!(seq, b"T");
    }

    // ========================================================================
    // quality_scores_slice — additional edge cases
    // ========================================================================

    #[test]
    fn test_quality_scores_slice_high_qualities() {
        let mut rec = make_bam_bytes(0, 0, 0, b"rd", &[], 3, -1, -1, &[]);
        let qo = qual_offset(&rec);
        rec[qo] = 40;
        rec[qo + 1] = 93; // max Phred
        rec[qo + 2] = 0; // min Phred
        let quals = quality_scores_slice(&rec);
        assert_eq!(quals, &[40, 93, 0]);
    }

    #[test]
    fn test_quality_scores_slice_mut_modify_all() {
        let mut rec = make_bam_bytes(0, 0, 0, b"rd", &[], 4, -1, -1, &[]);
        let quals = quality_scores_slice_mut(&mut rec);
        for (i, q) in quals.iter_mut().enumerate() {
            *q = (i * 10) as u8;
        }
        assert_eq!(quality_scores_slice(&rec), &[0, 10, 20, 30]);
    }

    // ========================================================================
    // set_base — additional edge cases
    // ========================================================================

    #[test]
    fn test_set_base_all_bases() {
        // Verify set_base/get_base roundtrip for all standard bases
        let mut rec = make_bam_bytes(0, 0, 0, b"rd", &[], 8, -1, -1, &[]);
        let so = seq_offset(&rec);
        let bases = [b'A', b'C', b'G', b'T', b'N', b'R', b'Y', b'M'];
        let codes: [u8; 8] = [1, 2, 4, 8, 15, 5, 10, 3];
        for (i, &base) in bases.iter().enumerate() {
            set_base(&mut rec, so, i, base);
        }
        for (i, &expected_code) in codes.iter().enumerate() {
            assert_eq!(
                get_base(&rec, so, i),
                expected_code,
                "base {} at pos {i}",
                bases[i] as char
            );
        }
    }

    #[test]
    fn test_set_base_preserves_neighbor() {
        // Ensure setting an even position doesn't corrupt the odd neighbor and vice versa
        let mut rec = make_bam_bytes(0, 0, 0, b"rd", &[], 2, -1, -1, &[]);
        let so = seq_offset(&rec);
        set_base(&mut rec, so, 0, b'A');
        set_base(&mut rec, so, 1, b'T');
        assert_eq!(get_base(&rec, so, 0), 1); // A
        assert_eq!(get_base(&rec, so, 1), 8); // T
        // Now change only position 0
        set_base(&mut rec, so, 0, b'G');
        assert_eq!(get_base(&rec, so, 0), 4); // G
        assert_eq!(get_base(&rec, so, 1), 8); // T still intact
    }

    // ========================================================================
    // set_flags
    // ========================================================================

    #[test]
    fn test_set_flags_roundtrip() {
        let mut rec = make_bam_bytes(0, 0, 0, b"rd", &[], 0, -1, -1, &[]);
        assert_eq!(flags(&rec), 0);
        set_flags(&mut rec, flags::PAIRED | flags::REVERSE | flags::DUPLICATE);
        assert_eq!(flags(&rec), flags::PAIRED | flags::REVERSE | flags::DUPLICATE);
        set_flags(&mut rec, 0);
        assert_eq!(flags(&rec), 0);
    }

    // ========================================================================
    // raw_records_to_record_bufs
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

    // ========================================================================
    // set_qual — additional edge cases
    // ========================================================================

    #[test]
    fn test_set_qual_boundary_values() {
        let mut rec = make_bam_bytes(0, 0, 0, b"rd", &[], 3, -1, -1, &[]);
        let qo = qual_offset(&rec);
        set_qual(&mut rec, qo, 0, 0); // min
        set_qual(&mut rec, qo, 1, 93); // max typical
        set_qual(&mut rec, qo, 2, 255); // max possible
        assert_eq!(get_qual(&rec, qo, 0), 0);
        assert_eq!(get_qual(&rec, qo, 1), 93);
        assert_eq!(get_qual(&rec, qo, 2), 255);
    }

    // ========================================================================
    // mask_base
    // ========================================================================

    #[test]
    fn test_mask_base_sets_n() {
        let mut rec = make_bam_bytes(0, 0, 0, b"rd", &[], 4, -1, -1, &[]);
        let so = seq_offset(&rec);
        // Set all bases to ACGT
        set_base(&mut rec, so, 0, b'A');
        set_base(&mut rec, so, 1, b'C');
        set_base(&mut rec, so, 2, b'G');
        set_base(&mut rec, so, 3, b'T');
        // Mask positions 1 and 3
        mask_base(&mut rec, so, 1);
        mask_base(&mut rec, so, 3);
        assert_eq!(get_base(&rec, so, 0), 1); // A preserved
        assert_eq!(get_base(&rec, so, 1), 15); // N
        assert_eq!(get_base(&rec, so, 2), 4); // G preserved
        assert_eq!(get_base(&rec, so, 3), 15); // N
    }
}
