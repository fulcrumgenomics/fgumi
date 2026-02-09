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
                // Integer types
                b'c' => Some((i8::from_le_bytes([aux_data[pos + 3]]) as u64, true)),
                b'C' => Some((u64::from(aux_data[pos + 3]), true)),
                b's' => {
                    Some((i16::from_le_bytes([aux_data[pos + 3], aux_data[pos + 4]]) as u64, true))
                }
                b'S' => Some((
                    u64::from(u16::from_le_bytes([aux_data[pos + 3], aux_data[pos + 4]])),
                    true,
                )),
                b'i' => Some((
                    i32::from_le_bytes([
                        aux_data[pos + 3],
                        aux_data[pos + 4],
                        aux_data[pos + 5],
                        aux_data[pos + 6],
                    ]) as u64,
                    true,
                )),
                b'I' => Some((
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

    // Parse the numeric part
    let mut value = 0u64;
    for &c in num_part {
        if c.is_ascii_digit() {
            value = value.saturating_mul(10).saturating_add(u64::from(c - b'0'));
        } else {
            return Some((0, true)); // Invalid format, return default
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

/// Calculate the size of a tag value based on its type.
fn tag_value_size(val_type: u8, data: &[u8]) -> Option<usize> {
    Some(match val_type {
        b'A' | b'c' | b'C' => 1,
        b's' | b'S' => 2,
        b'i' | b'I' | b'f' => 4,
        b'Z' | b'H' => data.iter().position(|&b| b == 0)? + 1,
        b'B' => {
            if data.len() < 5 {
                return None;
            }
            let elem_type = data[0];
            let count = u32::from_le_bytes([data[1], data[2], data[3], data[4]]) as usize;
            let elem_size = match elem_type {
                b'c' | b'C' => 1,
                b's' | b'S' => 2,
                b'i' | b'I' | b'f' => 4,
                _ => return None,
            };
            5 + count * elem_size
        }
        _ => return None,
    })
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

    // Get CIGAR ops for unclipped position calculation
    let a_cigar = get_cigar_ops(a);
    let b_cigar = get_cigar_ops(b);

    // Get MC tags for mate's unclipped position
    let a_mc = find_mc_tag_in_record(a);
    let b_mc = find_mc_tag_in_record(b);

    // Calculate unclipped 5' positions
    let a_unclipped_pos = unclipped_5prime(a_pos, a_reverse, a_cigar);
    let b_unclipped_pos = unclipped_5prime(b_pos, b_reverse, b_cigar);

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
        let num: i32 = cigar[num_start..i].parse().unwrap_or(1);

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

        let num: i32 = cigar[num_start..i].parse().unwrap_or(1);
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
#[allow(unsafe_code)]
pub fn get_cigar_ops(bam: &[u8]) -> &[u32] {
    let l_read_name = bam[8] as usize;
    let n_cigar_op = u16::from_le_bytes([bam[12], bam[13]]) as usize;

    if n_cigar_op == 0 {
        return &[];
    }

    let cigar_start = 32 + l_read_name;
    let cigar_end = cigar_start + n_cigar_op * 4;

    if cigar_end > bam.len() {
        return &[];
    }

    // SAFETY: CIGAR operations are u32 aligned in BAM format per the BAM spec.
    // The BAM file format guarantees that the CIGAR data starts at offset 32 + l_read_name,
    // which is always aligned to 4 bytes since the fixed header is 32 bytes and l_read_name
    // includes the null terminator and any padding needed for alignment.
    #[allow(clippy::cast_ptr_alignment)]
    unsafe {
        std::slice::from_raw_parts(bam[cigar_start..].as_ptr().cast::<u32>(), n_cigar_op)
    }
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

    #[test]
    fn test_find_mi_tag_unsigned_int_types() {
        // MI:C:200 (unsigned byte)
        let aux = [b'M', b'I', b'C', 200];
        assert_eq!(find_mi_tag(&aux), Some((200, true)));

        // MI:S:5000 (unsigned short)
        let mut aux = [b'M', b'I', b'S', 0, 0];
        aux[3..5].copy_from_slice(&5000u16.to_le_bytes());
        assert_eq!(find_mi_tag(&aux), Some((5000, true)));

        // MI:I:100000 (unsigned int)
        let mut aux = [b'M', b'I', b'I', 0, 0, 0, 0];
        aux[3..7].copy_from_slice(&100_000u32.to_le_bytes());
        assert_eq!(find_mi_tag(&aux), Some((100_000, true)));
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
        // MI:Z:abc\0 -> non-numeric chars return default (0, true)
        let aux = b"MIZabc\x00";
        assert_eq!(find_mi_tag(aux), Some((0, true)));
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
}
