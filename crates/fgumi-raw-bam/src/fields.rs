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

/// BAM file magic bytes (`BAM\1`).
pub const BAM_MAGIC: &[u8; 4] = b"BAM\x01";

/// Minimum length of a valid BAM record (the 32-byte fixed header).
/// All read-primitive functions that access fixed-offset fields (`flags`, `mapq`,
/// `ref_id`, `pos`, etc.) require at least this many bytes.
pub const MIN_BAM_RECORD_LEN: usize = 32;

/// Borrowed read-only view over a complete BAM record's bytes.
///
/// Wraps `&'a [u8]`. All BAM-spec accessors (`flags`, `pos`, `cigar`, `sequence`,
/// `tags`, …) are inherent methods on this type.
///
/// All accessors assume `bytes.len() >= MIN_BAM_RECORD_LEN`. Use [`RawRecordView::new`]
/// when the caller can vouch for that (hot path); use [`RawRecordView::try_new`] for
/// untrusted input.
#[derive(Copy, Clone, Debug)]
pub struct RawRecordView<'a>(&'a [u8]);

impl<'a> RawRecordView<'a> {
    /// Wrap raw BAM bytes without validation.
    ///
    /// In debug builds, asserts `bytes.len() >= MIN_BAM_RECORD_LEN`. In release
    /// builds, accessor methods on a too-short slice will panic on out-of-bounds
    /// indexing.
    #[inline]
    #[must_use]
    pub const fn new(bytes: &'a [u8]) -> Self {
        debug_assert!(bytes.len() >= MIN_BAM_RECORD_LEN);
        Self(bytes)
    }

    /// Wrap raw BAM bytes, returning `None` if too short for the fixed header.
    #[inline]
    #[must_use]
    pub const fn try_new(bytes: &'a [u8]) -> Option<Self> {
        if bytes.len() >= MIN_BAM_RECORD_LEN { Some(Self(bytes)) } else { None }
    }

    /// Returns the underlying byte slice.
    #[inline]
    #[must_use]
    pub const fn as_bytes(&self) -> &'a [u8] {
        self.0
    }
}

impl<'a> RawRecordView<'a> {
    // -- header (fixed-offset) --

    /// Extract reference sequence ID from this record.
    #[inline]
    #[must_use]
    pub fn ref_id(&self) -> i32 {
        ref_id(self.0)
    }

    /// Extract 0-based leftmost position from this record.
    #[inline]
    #[must_use]
    pub fn pos(&self) -> i32 {
        pos(self.0)
    }

    /// Extract mapping quality from this record.
    #[inline]
    #[must_use]
    pub fn mapq(&self) -> u8 {
        mapq(self.0)
    }

    /// Extract bitwise flags from this record.
    #[inline]
    #[must_use]
    pub fn flags(&self) -> u16 {
        flags(self.0)
    }

    /// Extract `l_read_name` (read-name length including NUL) from this record.
    #[inline]
    #[must_use]
    pub fn l_read_name(&self) -> u8 {
        l_read_name(self.0)
    }

    /// Extract number of CIGAR operations from this record.
    #[inline]
    #[must_use]
    pub fn n_cigar_op(&self) -> u16 {
        n_cigar_op(self.0)
    }

    /// Extract sequence length from this record.
    #[inline]
    #[must_use]
    pub fn l_seq(&self) -> u32 {
        l_seq(self.0)
    }

    /// Extract mate reference sequence ID from this record.
    #[inline]
    #[must_use]
    pub fn mate_ref_id(&self) -> i32 {
        mate_ref_id(self.0)
    }

    /// Extract mate 0-based position from this record.
    #[inline]
    #[must_use]
    pub fn mate_pos(&self) -> i32 {
        mate_pos(self.0)
    }

    /// Extract template length (tlen) from this record.
    #[inline]
    #[must_use]
    pub fn template_length(&self) -> i32 {
        template_length(self.0)
    }

    /// Extract read name (without null terminator) from this record.
    #[inline]
    #[must_use]
    pub fn read_name(&self) -> &'a [u8] {
        read_name(self.0)
    }

    /// Extract the BAM bin field (computed by the upstream encoder from `pos` and CIGAR).
    #[inline]
    #[must_use]
    pub fn bin(&self) -> u16 {
        u16::from_le_bytes([self.0[10], self.0[11]])
    }

    // -- flag-bit convenience --

    /// Returns `true` if the read is paired in sequencing.
    #[inline]
    #[must_use]
    pub fn is_paired(&self) -> bool {
        self.flags() & flags::PAIRED != 0
    }

    /// Returns `true` if the read is unmapped.
    #[inline]
    #[must_use]
    pub fn is_unmapped(&self) -> bool {
        self.flags() & flags::UNMAPPED != 0
    }

    /// Returns `true` if the mate is unmapped.
    #[inline]
    #[must_use]
    pub fn is_mate_unmapped(&self) -> bool {
        self.flags() & flags::MATE_UNMAPPED != 0
    }

    /// Returns `true` if the read is reverse complemented.
    #[inline]
    #[must_use]
    pub fn is_reverse(&self) -> bool {
        self.flags() & flags::REVERSE != 0
    }

    /// Returns `true` if the mate is reverse complemented.
    #[inline]
    #[must_use]
    pub fn is_mate_reverse(&self) -> bool {
        self.flags() & flags::MATE_REVERSE != 0
    }

    /// Returns `true` if the read is the first segment (R1).
    #[inline]
    #[must_use]
    pub fn is_first_segment(&self) -> bool {
        self.flags() & flags::FIRST_SEGMENT != 0
    }

    /// Returns `true` if the read is the last segment (R2).
    #[inline]
    #[must_use]
    pub fn is_last_segment(&self) -> bool {
        self.flags() & flags::LAST_SEGMENT != 0
    }

    /// Returns `true` if the read is a secondary alignment.
    #[inline]
    #[must_use]
    pub fn is_secondary(&self) -> bool {
        self.flags() & flags::SECONDARY != 0
    }

    /// Returns `true` if the read fails quality controls.
    #[inline]
    #[must_use]
    pub fn is_qc_fail(&self) -> bool {
        self.flags() & flags::QC_FAIL != 0
    }

    /// Returns `true` if the read is a PCR or optical duplicate.
    #[inline]
    #[must_use]
    pub fn is_duplicate(&self) -> bool {
        self.flags() & flags::DUPLICATE != 0
    }

    /// Returns `true` if the read is a supplementary alignment.
    #[inline]
    #[must_use]
    pub fn is_supplementary(&self) -> bool {
        self.flags() & flags::SUPPLEMENTARY != 0
    }
}

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

/// Fixed-size tag value sizes indexed by type byte.
/// 0 = variable/unknown (needs special handling).
pub(crate) const TAG_FIXED_SIZES: [u8; 256] = {
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

/// Extract template length (tlen) from a BAM record.
#[inline]
#[must_use]
pub fn template_length(bam: &[u8]) -> i32 {
    i32::from_le_bytes([bam[28], bam[29], bam[30], bam[31]])
}

/// Extract read name (without null terminator) from a BAM record.
#[inline]
#[must_use]
pub fn read_name(bam: &[u8]) -> &[u8] {
    let l = bam[8] as usize;
    if l > 1 { &bam[32..32 + l - 1] } else { &[] }
}

// ============================================================================
// Write Primitives
// ============================================================================

/// Set flags (u16) in a BAM record.
#[inline]
pub fn set_flags(bam: &mut [u8], new_flags: u16) {
    bam[14..16].copy_from_slice(&new_flags.to_le_bytes());
}

/// Set the mapping quality in a BAM record.
#[inline]
pub fn set_mapq(bam: &mut [u8], new_mapq: u8) {
    bam[9] = new_mapq;
}

/// Set the mate reference sequence ID in a BAM record.
#[inline]
pub fn set_mate_ref_id(bam: &mut [u8], val: i32) {
    bam[20..24].copy_from_slice(&val.to_le_bytes());
}

/// Set the mate 0-based position in a BAM record.
#[inline]
pub fn set_mate_pos(bam: &mut [u8], val: i32) {
    bam[24..28].copy_from_slice(&val.to_le_bytes());
}

/// Set the template length (TLEN) in a BAM record.
#[inline]
pub fn set_template_length(bam: &mut [u8], val: i32) {
    bam[28..32].copy_from_slice(&val.to_le_bytes());
}

/// Set the reference sequence ID in a BAM record.
#[inline]
pub fn set_ref_id(bam: &mut [u8], val: i32) {
    bam[0..4].copy_from_slice(&val.to_le_bytes());
}

/// Set the 0-based leftmost position in a BAM record.
#[inline]
pub fn set_pos(bam: &mut [u8], val: i32) {
    bam[4..8].copy_from_slice(&val.to_le_bytes());
}

/// Calculate the offset to auxiliary data in a BAM record.
///
/// `aux_offset = 32 + l_read_name + n_cigar_op*4 + (l_seq+1)/2 + l_seq`
#[inline]
#[must_use]
pub fn aux_data_offset(l_read_name: usize, n_cigar_op: usize, l_seq: usize) -> usize {
    32 + l_read_name + n_cigar_op * 4 + l_seq.div_ceil(2) + l_seq
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

/// Extract coordinate key fields directly from BAM bytes.
///
/// Returns (tid, pos, reverse, name).
#[inline]
#[must_use]
pub fn extract_coordinate_fields(bam_bytes: &[u8]) -> (i32, i32, bool, &[u8]) {
    debug_assert!(
        bam_bytes.len() >= MIN_BAM_RECORD_LEN,
        "BAM record too short ({} < {MIN_BAM_RECORD_LEN})",
        bam_bytes.len()
    );

    let tid = ref_id(bam_bytes);
    let p = pos(bam_bytes);
    let flag = flags(bam_bytes);
    let reverse = (flag & flags::REVERSE) != 0;
    let unmapped = (flag & flags::UNMAPPED) != 0;
    let name = read_name(bam_bytes);

    if unmapped { (i32::MAX, i32::MAX, false, name) } else { (tid, p, reverse, name) }
}

/// Packed boolean flags for [`TemplateCoordFields`].
///
/// Stores five boolean fields in a single `u8` to keep the parent struct
/// under the threshold for `clippy::struct_excessive_bools`.
#[derive(Clone, Copy, Debug, Default, PartialEq, Eq)]
pub struct TemplateCoordFlags(u8);

impl TemplateCoordFlags {
    const REVERSE: u8 = 1;
    const MATE_REVERSE: u8 = 1 << 1;
    const UNMAPPED: u8 = 1 << 2;
    const MATE_UNMAPPED: u8 = 1 << 3;
    const PAIRED: u8 = 1 << 4;

    /// Build flags from a BAM flag word, extracting the relevant bits.
    #[must_use]
    pub fn from_flag(flag: u16) -> Self {
        let mut bits = 0u8;
        if (flag & flags::REVERSE) != 0 {
            bits |= Self::REVERSE;
        }
        if (flag & flags::MATE_REVERSE) != 0 {
            bits |= Self::MATE_REVERSE;
        }
        if (flag & flags::UNMAPPED) != 0 {
            bits |= Self::UNMAPPED;
        }
        if (flag & flags::MATE_UNMAPPED) != 0 {
            bits |= Self::MATE_UNMAPPED;
        }
        if (flag & flags::PAIRED) != 0 {
            bits |= Self::PAIRED;
        }
        Self(bits)
    }

    /// Whether the read is on the reverse strand.
    #[inline]
    #[must_use]
    pub fn reverse(self) -> bool {
        self.0 & Self::REVERSE != 0
    }

    /// Whether the mate is on the reverse strand.
    #[inline]
    #[must_use]
    pub fn mate_reverse(self) -> bool {
        self.0 & Self::MATE_REVERSE != 0
    }

    /// Whether the read is unmapped.
    #[inline]
    #[must_use]
    pub fn unmapped(self) -> bool {
        self.0 & Self::UNMAPPED != 0
    }

    /// Whether the mate is unmapped.
    #[inline]
    #[must_use]
    pub fn mate_unmapped(self) -> bool {
        self.0 & Self::MATE_UNMAPPED != 0
    }

    /// Whether the read is paired.
    #[inline]
    #[must_use]
    pub fn paired(self) -> bool {
        self.0 & Self::PAIRED != 0
    }
}

/// Fields extracted for template-coordinate sorting.
pub struct TemplateCoordFields<'a> {
    pub tid: i32,
    pub pos: i32,
    pub mate_tid: i32,
    pub mate_pos: i32,
    pub flags: TemplateCoordFlags,
    pub name: &'a [u8],
    pub l_read_name: usize,
    pub n_cigar_op: usize,
    pub l_seq: usize,
}

/// Extract template-coordinate key fields directly from BAM bytes.
#[inline]
#[must_use]
pub fn extract_template_coordinate_fields(bam_bytes: &[u8]) -> TemplateCoordFields<'_> {
    debug_assert!(
        bam_bytes.len() >= MIN_BAM_RECORD_LEN,
        "BAM record too short ({} < {MIN_BAM_RECORD_LEN})",
        bam_bytes.len()
    );

    let tid = ref_id(bam_bytes);
    let p = pos(bam_bytes);
    let l_rn = l_read_name(bam_bytes) as usize;
    let n_co = n_cigar_op(bam_bytes) as usize;
    let flag = flags(bam_bytes);
    let l_s = l_seq(bam_bytes) as usize;
    let m_tid = mate_ref_id(bam_bytes);
    let m_pos = mate_pos(bam_bytes);
    let name = read_name(bam_bytes);
    let tc_flags = TemplateCoordFlags::from_flag(flag);

    TemplateCoordFields {
        tid,
        pos: p,
        mate_tid: m_tid,
        mate_pos: m_pos,
        flags: tc_flags,
        name,
        l_read_name: l_rn,
        n_cigar_op: n_co,
        l_seq: l_s,
    }
}

impl<'a> RawRecordView<'a> {
    /// Extract `(tid, pos, reverse, name)` for coordinate-key sorting.
    /// Unmapped reads return `(i32::MAX, i32::MAX, false, name)`.
    #[inline]
    #[must_use]
    pub fn coordinate_fields(&self) -> (i32, i32, bool, &'a [u8]) {
        extract_coordinate_fields(self.0)
    }

    /// Extract template-coordinate sort fields in a single pass.
    #[inline]
    #[must_use]
    pub fn template_coordinate_fields(&self) -> TemplateCoordFields<'a> {
        extract_template_coordinate_fields(self.0)
    }
}

/// Borrowed fixed-length mutable view over a complete BAM record's bytes.
///
/// Wraps `&'a mut [u8]`. Provides every fixed-offset setter plus per-base /
/// per-quality writes, but cannot change the record's byte length.
///
/// Read methods are accessed via [`RawRecordMut::view`].
#[derive(Debug)]
pub struct RawRecordMut<'a>(&'a mut [u8]);

impl<'a> RawRecordMut<'a> {
    /// Wrap raw BAM bytes mutably without validation.
    ///
    /// In debug builds, asserts `bytes.len() >= MIN_BAM_RECORD_LEN`. In release
    /// builds, setter methods on a too-short slice will panic on out-of-bounds
    /// indexing.
    #[inline]
    #[must_use]
    pub fn new(bytes: &'a mut [u8]) -> Self {
        debug_assert!(bytes.len() >= MIN_BAM_RECORD_LEN);
        Self(bytes)
    }

    /// Wrap raw BAM bytes mutably, returning `None` if too short for the fixed header.
    #[inline]
    #[must_use]
    pub fn try_new(bytes: &'a mut [u8]) -> Option<Self> {
        if bytes.len() >= MIN_BAM_RECORD_LEN { Some(Self(bytes)) } else { None }
    }

    /// Borrow as a read-only view (lifetime tied to `&self`).
    #[inline]
    #[must_use]
    pub fn view(&self) -> RawRecordView<'_> {
        RawRecordView(self.0)
    }

    /// Returns the underlying byte slice as a shared reference.
    #[inline]
    #[must_use]
    pub fn as_bytes(&self) -> &[u8] {
        self.0
    }

    /// Returns the underlying byte slice as a mutable reference.
    #[inline]
    #[must_use]
    pub fn as_bytes_mut(&mut self) -> &mut [u8] {
        self.0
    }

    // -- header (fixed-offset) setters --

    /// Set the reference sequence ID.
    #[inline]
    pub fn set_ref_id(&mut self, v: i32) {
        set_ref_id(self.0, v);
    }

    /// Set the 0-based leftmost position.
    #[inline]
    pub fn set_pos(&mut self, v: i32) {
        set_pos(self.0, v);
    }

    /// Set the mapping quality.
    #[inline]
    pub fn set_mapq(&mut self, v: u8) {
        set_mapq(self.0, v);
    }

    /// Set the bitwise flags.
    #[inline]
    pub fn set_flags(&mut self, v: u16) {
        set_flags(self.0, v);
    }

    /// Set the mate reference sequence ID.
    #[inline]
    pub fn set_mate_ref_id(&mut self, v: i32) {
        set_mate_ref_id(self.0, v);
    }

    /// Set the mate 0-based position.
    #[inline]
    pub fn set_mate_pos(&mut self, v: i32) {
        set_mate_pos(self.0, v);
    }

    /// Set the template length (TLEN).
    #[inline]
    pub fn set_template_length(&mut self, v: i32) {
        set_template_length(self.0, v);
    }

    /// Set the BAM bin field (bytes 10–11, little-endian u16).
    #[inline]
    pub fn set_bin(&mut self, v: u16) {
        self.0[10..12].copy_from_slice(&v.to_le_bytes());
    }
}

#[cfg(test)]
#[allow(clippy::identity_op)]
mod tests {
    use super::*;
    use crate::testutil::*;

    // ========================================================================
    // extract_coordinate_fields tests
    // ========================================================================

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

    // ========================================================================
    // aux_data_offset tests
    // ========================================================================

    #[test]
    fn test_aux_data_offset() {
        // l_read_name=5, n_cigar_op=2, l_seq=50
        // offset = 32 + 5 + 2*4 + (50+1)/2 + 50 = 32 + 5 + 8 + 25 + 50 = 120
        assert_eq!(aux_data_offset(5, 2, 50), 120);
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

    #[test]
    fn test_tag_value_size_unknown_type() {
        // Unknown type byte -> None
        assert_eq!(tag_value_size(b'?', &[0; 10]), None);
    }

    #[test]
    fn test_tag_value_size_b_array_invalid_elem_type() {
        // B-array with elem_type 'Z' (not a fixed-size type) -> None
        let mut data = vec![b'Z']; // invalid elem_type
        data.extend_from_slice(&3u32.to_le_bytes()); // count
        data.extend_from_slice(&[0; 12]); // data
        assert_eq!(tag_value_size(b'B', &data), None);
    }

    #[test]
    fn test_tag_value_size_b_array_truncated_header() {
        // B-array with less than 5 bytes of header -> None
        let data = [b'i', 0, 0]; // only 3 bytes, need 5
        assert_eq!(tag_value_size(b'B', &data), None);
    }

    #[test]
    fn test_tag_value_size_string_no_null() {
        // Z-type string without null terminator -> None
        let data = b"hello";
        assert_eq!(tag_value_size(b'Z', data), None);
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
    fn test_set_flags_roundtrip() {
        let mut rec = make_bam_bytes(0, 0, 0, b"rd", &[], 0, -1, -1, &[]);
        assert_eq!(flags(&rec), 0);
        set_flags(&mut rec, flags::PAIRED | flags::REVERSE | flags::DUPLICATE);
        assert_eq!(flags(&rec), flags::PAIRED | flags::REVERSE | flags::DUPLICATE);
        set_flags(&mut rec, 0);
        assert_eq!(flags(&rec), 0);
    }

    // ========================================================================
    // aux_data_offset_from_record tests
    // ========================================================================

    #[test]
    fn test_aux_data_offset_from_record_truncated_no_panic() {
        // aux_data_offset_from_record returns None for records too short to read header
        assert_eq!(aux_data_offset_from_record(&[0u8; 19]), None);
        assert!(aux_data_offset_from_record(&[0u8; 20]).is_some());
    }

    #[test]
    fn test_aux_data_offset_from_record_basic() {
        let rec = make_bam_bytes(0, 0, 0, b"rd", &[encode_op(0, 4)], 4, -1, -1, &[]);
        let offset = aux_data_offset_from_record(&rec)
            .expect("record should have valid header for aux offset");
        // 32 + l_read_name(3) + cigar(4) + seq_bytes(2) + qual(4) = 45
        assert_eq!(offset, 32 + 3 + 4 + 2 + 4);
    }

    #[test]
    fn test_aux_data_offset_from_record_too_short() {
        let short = vec![0u8; 10];
        assert!(aux_data_offset_from_record(&short).is_none());
    }

    // ========================================================================
    // aux_data_slice tests
    // ========================================================================

    #[test]
    fn test_aux_data_slice_out_of_bounds_returns_empty() {
        // When aux_data_offset exceeds record length, aux_data_slice returns empty
        let mut bam = make_bam_bytes(0, 0, 0, b"rea", &[], 10, -1, -1, &[]);

        // Modify l_seq to unreasonably large value -- offset will exceed record length
        bam[16..20].copy_from_slice(&1000u32.to_le_bytes());

        let offset = aux_data_offset_from_record(&bam)
            .expect("record should have valid header for aux offset");
        assert!(offset > bam.len());

        // aux_data_slice returns empty instead of panicking
        assert!(aux_data_slice(&bam).is_empty());
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

    // ========================================================================
    // seq_offset / qual_offset tests
    // ========================================================================

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
    // read_name edge cases
    // ========================================================================

    #[test]
    fn test_read_name_empty() {
        // l_read_name = 1 means just the null terminator -> empty name
        let mut rec = vec![0u8; 34];
        rec[8] = 1; // l_read_name = 1
        rec[32] = 0; // null terminator only
        assert_eq!(read_name(&rec), b"");
    }

    #[test]
    fn test_read_name_zero_length() {
        // l_read_name = 0 -> empty name (edge case)
        let mut rec = vec![0u8; 33];
        rec[8] = 0;
        assert_eq!(read_name(&rec), b"");
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
    // Field accessor roundtrip tests
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
    // TemplateCoordFlags tests
    // ========================================================================

    #[test]
    fn test_template_coord_flags_all_set() {
        let flag = flags::REVERSE
            | flags::MATE_REVERSE
            | flags::UNMAPPED
            | flags::MATE_UNMAPPED
            | flags::PAIRED;
        let tcf = TemplateCoordFlags::from_flag(flag);
        assert!(tcf.reverse());
        assert!(tcf.mate_reverse());
        assert!(tcf.unmapped());
        assert!(tcf.mate_unmapped());
        assert!(tcf.paired());
    }

    #[test]
    fn test_template_coord_flags_none_set() {
        let tcf = TemplateCoordFlags::from_flag(0);
        assert!(!tcf.reverse());
        assert!(!tcf.mate_reverse());
        assert!(!tcf.unmapped());
        assert!(!tcf.mate_unmapped());
        assert!(!tcf.paired());
    }

    #[test]
    fn test_template_coord_flags_individual_bits() {
        let tcf = TemplateCoordFlags::from_flag(flags::REVERSE | flags::PAIRED);
        assert!(tcf.reverse());
        assert!(!tcf.mate_reverse());
        assert!(!tcf.unmapped());
        assert!(!tcf.mate_unmapped());
        assert!(tcf.paired());
    }

    #[test]
    fn test_template_coord_flags_default() {
        let tcf = TemplateCoordFlags::default();
        assert!(!tcf.reverse());
        assert!(!tcf.mate_reverse());
        assert!(!tcf.unmapped());
        assert!(!tcf.mate_unmapped());
        assert!(!tcf.paired());
    }

    // ========================================================================
    // extract_template_coordinate_fields tests
    // ========================================================================

    #[test]
    fn test_extract_template_coordinate_fields_basic() {
        let rec = make_bam_bytes_with_tlen(
            3,                                                    // tid
            200,                                                  // pos
            flags::PAIRED | flags::REVERSE | flags::MATE_REVERSE, // flag
            b"read1",                                             // name
            &[encode_op(0, 10)],                                  // 10M cigar
            6,                                                    // seq_len
            5,                                                    // mate_tid
            400,                                                  // mate_pos
            150,                                                  // tlen
            &[],                                                  // aux
        );
        let fields = extract_template_coordinate_fields(&rec);
        assert_eq!(fields.tid, 3);
        assert_eq!(fields.pos, 200);
        assert_eq!(fields.mate_tid, 5);
        assert_eq!(fields.mate_pos, 400);
        assert!(fields.flags.reverse());
        assert!(fields.flags.mate_reverse());
        assert!(!fields.flags.unmapped());
        assert!(!fields.flags.mate_unmapped());
        assert!(fields.flags.paired());
        assert_eq!(fields.name, b"read1");
        assert_eq!(fields.l_read_name, 6); // "read1" + NUL
        assert_eq!(fields.n_cigar_op, 1);
        assert_eq!(fields.l_seq, 6);
    }

    #[test]
    fn test_extract_template_coordinate_fields_unmapped() {
        let rec = make_bam_bytes(
            -1,
            -1,
            flags::UNMAPPED | flags::PAIRED | flags::MATE_UNMAPPED,
            b"rd",
            &[],
            0,
            -1,
            -1,
            &[],
        );
        let fields = extract_template_coordinate_fields(&rec);
        assert_eq!(fields.tid, -1);
        assert_eq!(fields.pos, -1);
        assert!(fields.flags.unmapped());
        assert!(fields.flags.mate_unmapped());
        assert!(fields.flags.paired());
        assert_eq!(fields.name, b"rd");
    }

    // ========================================================================
    // Field setter tests
    // ========================================================================

    #[test]
    fn test_set_mapq_roundtrip() {
        let mut rec = make_bam_bytes(0, 0, 0, b"r1", &[], 4, -1, -1, &[]);
        assert_eq!(mapq(&rec), 0);
        set_mapq(&mut rec, 42);
        assert_eq!(mapq(&rec), 42);
    }

    #[test]
    fn test_set_mate_ref_id_roundtrip() {
        let mut rec = make_bam_bytes(0, 0, 0, b"r1", &[], 4, -1, -1, &[]);
        assert_eq!(mate_ref_id(&rec), -1);
        set_mate_ref_id(&mut rec, 5);
        assert_eq!(mate_ref_id(&rec), 5);
    }

    #[test]
    fn test_set_mate_pos_roundtrip() {
        let mut rec = make_bam_bytes(0, 0, 0, b"r1", &[], 4, -1, -1, &[]);
        assert_eq!(mate_pos(&rec), -1);
        set_mate_pos(&mut rec, 12345);
        assert_eq!(mate_pos(&rec), 12345);
    }

    #[test]
    fn test_set_template_length_roundtrip() {
        let mut rec = make_bam_bytes(0, 0, 0, b"r1", &[], 4, -1, -1, &[]);
        assert_eq!(template_length(&rec), 0);
        set_template_length(&mut rec, -300);
        assert_eq!(template_length(&rec), -300);
    }

    #[test]
    fn test_set_ref_id_roundtrip() {
        let mut rec = make_bam_bytes(0, 0, 0, b"r1", &[], 4, -1, -1, &[]);
        assert_eq!(ref_id(&rec), 0);
        set_ref_id(&mut rec, 3);
        assert_eq!(ref_id(&rec), 3);
    }

    #[test]
    fn test_set_pos_roundtrip() {
        let mut rec = make_bam_bytes(0, 0, 0, b"r1", &[], 4, -1, -1, &[]);
        assert_eq!(pos(&rec), 0);
        set_pos(&mut rec, 9999);
        assert_eq!(pos(&rec), 9999);
    }

    #[test]
    fn test_raw_record_view_new_and_as_bytes() {
        let bytes = vec![0u8; 32];
        let view = RawRecordView::new(&bytes);
        assert_eq!(view.as_bytes().len(), 32);
    }

    #[test]
    fn test_raw_record_view_try_new_too_short_returns_none() {
        let bytes = vec![0u8; 31];
        assert!(RawRecordView::try_new(&bytes).is_none());
    }

    #[test]
    fn test_raw_record_view_try_new_min_length_succeeds() {
        let bytes = vec![0u8; MIN_BAM_RECORD_LEN];
        assert!(RawRecordView::try_new(&bytes).is_some());
    }

    #[test]
    fn test_raw_record_view_header_accessors() {
        use crate::testutil::*;
        let mut rec = make_bam_bytes(
            3,
            200,
            flags::PAIRED | flags::REVERSE,
            b"read1",
            &[encode_op(0, 10)],
            6,
            5,
            400,
            &[],
        );
        rec[9] = 42; // mapq
        // bin
        rec[10..12].copy_from_slice(&4680u16.to_le_bytes());

        let v = RawRecordView::new(&rec);
        assert_eq!(v.ref_id(), 3);
        assert_eq!(v.pos(), 200);
        assert_eq!(v.mapq(), 42);
        assert_eq!(v.flags(), flags::PAIRED | flags::REVERSE);
        assert_eq!(v.bin(), 4680);
        assert_eq!(v.l_read_name(), 6);
        assert_eq!(v.n_cigar_op(), 1);
        assert_eq!(v.l_seq(), 6);
        assert_eq!(v.mate_ref_id(), 5);
        assert_eq!(v.mate_pos(), 400);
        assert_eq!(v.read_name(), b"read1");
        assert!(v.is_paired());
        assert!(v.is_reverse());
        assert!(!v.is_unmapped());
    }

    #[test]
    fn test_view_batched_coord_helpers() {
        use crate::testutil::*;
        let rec = make_bam_bytes(3, 200, flags::REVERSE, b"r", &[], 0, -1, -1, &[]);
        let v = RawRecordView::new(&rec);
        let (tid, pos, reverse, name) = v.coordinate_fields();
        assert_eq!((tid, pos, reverse), (3, 200, true));
        assert_eq!(name, b"r");

        let f = v.template_coordinate_fields();
        assert_eq!(f.tid, 3);
        assert!(f.flags.reverse());
    }

    // ========================================================================
    // RawRecordMut tests
    // ========================================================================

    #[test]
    fn test_raw_record_mut_constructors_and_view() {
        let mut bytes = vec![0u8; 32];
        let m = RawRecordMut::new(&mut bytes);
        assert_eq!(m.view().as_bytes().len(), 32);
    }

    #[test]
    fn test_raw_record_mut_setters() {
        use crate::testutil::*;
        let mut rec = make_bam_bytes(0, 0, 0, b"r", &[], 0, -1, -1, &[]);
        let mut m = RawRecordMut::new(&mut rec);
        m.set_ref_id(7);
        m.set_pos(123);
        m.set_mapq(40);
        m.set_flags(flags::PAIRED | flags::REVERSE);
        m.set_bin(4680);
        m.set_mate_ref_id(8);
        m.set_mate_pos(456);
        m.set_template_length(300);
        let v = m.view();
        assert_eq!(v.ref_id(), 7);
        assert_eq!(v.pos(), 123);
        assert_eq!(v.mapq(), 40);
        assert_eq!(v.flags(), flags::PAIRED | flags::REVERSE);
        assert_eq!(v.bin(), 4680);
        assert_eq!(v.mate_ref_id(), 8);
        assert_eq!(v.mate_pos(), 456);
        assert_eq!(v.template_length(), 300);
    }
}
