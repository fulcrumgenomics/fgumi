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
}
