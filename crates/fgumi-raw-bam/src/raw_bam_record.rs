//! Raw BAM record reader that bypasses noodles' Record type.
//!
//! This module provides [`RawRecord`], a zero-overhead wrapper around raw BAM bytes
//! that enables high-performance sorting without full record decoding.
//!
//! # Why Not Use `noodles::bam::Record`?
//!
//! The standard `noodles::bam::Record` type wraps raw bytes but doesn't expose them
//! via `AsRef<[u8]>`. This module provides direct access to raw bytes for:
//! - Zero-copy field extraction using fixed byte offsets
//! - Direct writes to output without re-encoding
//! - 3-4x lower memory usage than decoded `RecordBuf`
//!
//! # BAM Record Binary Layout
//!
//! ```text
//! Offset  Size  Field
//! ------  ----  -----
//! 0-3     4     refID (i32) - reference sequence ID
//! 4-7     4     pos (i32) - 0-based leftmost position
//! 8       1     l_read_name (u8) - length including null
//! 9       1     mapq (u8) - mapping quality
//! 10-11   2     bin (u16) - BAM bin
//! 12-13   2     n_cigar_op (u16) - CIGAR operation count
//! 14-15   2     flag (u16) - bitwise flags
//! 16-19   4     l_seq (u32) - sequence length
//! 20-23   4     next_refID (i32) - mate reference ID
//! 24-27   4     next_pos (i32) - mate position
//! 28-31   4     tlen (i32) - template length
//! 32+     var   read_name, CIGAR, sequence, quality, aux data
//! ```
//!
//! This module will be removed once noodles exposes raw bytes from Record
//! (<https://github.com/zaeleus/noodles/pull/373>).

use std::io::{self, Read, Write};

/// A raw BAM record stored as bytes.
///
/// This is a zero-overhead wrapper that provides `AsRef<[u8]>` access to the
/// underlying BAM record bytes, enabling high-performance field extraction
/// and direct output writes.
#[derive(Clone, Debug, Default, Eq, PartialEq)]
pub struct RawRecord(Vec<u8>);

impl RawRecord {
    /// Creates a new empty raw record.
    #[inline]
    #[must_use]
    pub fn new() -> Self {
        Self(Vec::new())
    }

    /// Creates a raw record with the given capacity.
    #[inline]
    #[must_use]
    pub fn with_capacity(capacity: usize) -> Self {
        Self(Vec::with_capacity(capacity))
    }

    /// Returns the length of the record in bytes.
    #[inline]
    #[must_use]
    pub fn len(&self) -> usize {
        self.0.len()
    }

    /// Returns the allocated capacity of the underlying byte buffer.
    ///
    /// Useful for memory-estimation purposes (e.g., [`MemoryEstimate`] implementations).
    #[inline]
    #[must_use]
    pub fn capacity(&self) -> usize {
        self.0.capacity()
    }

    /// Returns true if the record is empty.
    #[inline]
    #[must_use]
    pub fn is_empty(&self) -> bool {
        self.0.is_empty()
    }

    /// Clears the record, removing all bytes.
    #[inline]
    pub fn clear(&mut self) {
        self.0.clear();
    }

    /// Returns the inner bytes, consuming the record.
    #[inline]
    #[must_use]
    pub fn into_inner(self) -> Vec<u8> {
        self.0
    }

    /// Returns a mutable reference to the inner byte vector.
    ///
    /// Exposes the raw `Vec<u8>` for callers that need resizing operations
    /// (e.g., tag splice/replace) that require `&mut Vec<u8>`.
    ///
    /// **Invariants are the caller's responsibility.** [`RawRecord::as_mut_vec`]
    /// (and [`DerefMut`](std::ops::DerefMut) below) bypass every record-shape
    /// invariant: callers must ensure that after any mutation `l_read_name`,
    /// `n_cigar_op`, and `l_seq` (and the packed sequence / quality byte
    /// lengths derived from them) remain consistent with the underlying buffer.
    #[inline]
    pub fn as_mut_vec(&mut self) -> &mut Vec<u8> {
        &mut self.0
    }
}

impl AsRef<[u8]> for RawRecord {
    #[inline]
    fn as_ref(&self) -> &[u8] {
        &self.0
    }
}

impl std::ops::Deref for RawRecord {
    type Target = [u8];

    #[inline]
    fn deref(&self) -> &Self::Target {
        &self.0
    }
}

/// Direct mutable byte access. **Bypasses all record invariants** —
/// callers must keep `l_read_name`, `n_cigar_op`, and `l_seq` consistent with
/// the buffer after any mutation. See [`RawRecord::as_mut_vec`] for context.
impl std::ops::DerefMut for RawRecord {
    #[inline]
    fn deref_mut(&mut self) -> &mut Self::Target {
        &mut self.0
    }
}

impl From<Vec<u8>> for RawRecord {
    #[inline]
    fn from(buf: Vec<u8>) -> Self {
        Self(buf)
    }
}

/// Reads a single raw BAM record from the given reader.
///
/// This function reads:
/// 1. The 4-byte `block_size` prefix (little-endian u32)
/// 2. The record data (`block_size` bytes)
///
/// Returns the number of bytes read (excluding the 4-byte prefix), or 0 at EOF.
///
/// # Errors
///
/// Returns an error if:
/// - The reader returns an error
/// - EOF is reached in the middle of a record
pub fn read_raw_record<R>(reader: &mut R, record: &mut RawRecord) -> io::Result<usize>
where
    R: Read,
{
    let block_size = match read_block_size(reader)? {
        0 => return Ok(0),
        n => n,
    };

    record.0.resize(block_size, 0);
    reader.read_exact(&mut record.0)?;

    Ok(block_size)
}

/// Reads the 4-byte block size prefix.
///
/// Returns 0 at EOF (no bytes available).
fn read_block_size<R>(reader: &mut R) -> io::Result<usize>
where
    R: Read,
{
    let mut buf = [0u8; 4];

    // Try to read the first byte to detect EOF
    match reader.read(&mut buf[..1]) {
        Ok(0) => return Ok(0), // EOF
        Ok(1) => {}
        Ok(_) => unreachable!(),
        Err(ref e) if e.kind() == io::ErrorKind::Interrupted => {
            return read_block_size(reader);
        }
        Err(e) => return Err(e),
    }

    // Read the remaining 3 bytes
    reader.read_exact(&mut buf[1..])?;

    let n = u32::from_le_bytes(buf);
    usize::try_from(n).map_err(|e| io::Error::new(io::ErrorKind::InvalidData, e))
}

use crate::fields::{RawRecordMut, RawRecordView};
use crate::tags::{RawTagsEditor, RawTagsMut, RawTagsView};

impl RawRecord {
    /// Borrow as a read-only view.
    #[inline]
    #[must_use]
    pub fn view(&self) -> RawRecordView<'_> {
        RawRecordView::new(&self.0)
    }

    /// Borrow as a fixed-length mutable view.
    #[inline]
    #[must_use]
    pub fn view_mut(&mut self) -> RawRecordMut<'_> {
        RawRecordMut::new(&mut self.0)
    }

    /// Read-only view over this record's auxiliary tag section.
    #[inline]
    #[must_use]
    pub fn tags(&self) -> RawTagsView<'_> {
        self.view().tags()
    }

    /// Fixed-length mutable view over this record's auxiliary tag section.
    #[inline]
    #[must_use]
    pub fn tags_mut(&mut self) -> RawTagsMut<'_> {
        let len = self.0.len();
        let off = crate::fields::aux_data_offset_from_record(&self.0).unwrap_or(len);
        RawTagsMut::new(&mut self.0[off..])
    }

    /// Length-changing tag editor with cached aux offset.
    #[inline]
    #[must_use]
    pub fn tags_editor(&mut self) -> RawTagsEditor<'_> {
        RawTagsEditor::from_vec(&mut self.0)
    }

    // -- Header reads (return primitives — no lifetime borrow from self) --

    /// Extract reference sequence ID from this record.
    #[inline]
    #[must_use]
    pub fn ref_id(&self) -> i32 {
        RawRecordView::new(&self.0).ref_id()
    }

    /// Extract 0-based leftmost position from this record.
    #[inline]
    #[must_use]
    pub fn pos(&self) -> i32 {
        RawRecordView::new(&self.0).pos()
    }

    /// Extract mapping quality from this record.
    #[inline]
    #[must_use]
    pub fn mapq(&self) -> u8 {
        RawRecordView::new(&self.0).mapq()
    }

    /// Extract bitwise flags from this record.
    #[inline]
    #[must_use]
    pub fn flags(&self) -> u16 {
        RawRecordView::new(&self.0).flags()
    }

    /// Extract BAM bin field from this record.
    #[inline]
    #[must_use]
    pub fn bin(&self) -> u16 {
        RawRecordView::new(&self.0).bin()
    }

    /// Extract `l_read_name` (read-name length including NUL) from this record.
    #[inline]
    #[must_use]
    pub fn l_read_name(&self) -> u8 {
        RawRecordView::new(&self.0).l_read_name()
    }

    /// Extract number of CIGAR operations from this record.
    #[inline]
    #[must_use]
    pub fn n_cigar_op(&self) -> u16 {
        RawRecordView::new(&self.0).n_cigar_op()
    }

    /// Extract sequence length from this record.
    #[inline]
    #[must_use]
    pub fn l_seq(&self) -> u32 {
        RawRecordView::new(&self.0).l_seq()
    }

    /// Extract mate reference sequence ID from this record.
    #[inline]
    #[must_use]
    pub fn mate_ref_id(&self) -> i32 {
        RawRecordView::new(&self.0).mate_ref_id()
    }

    /// Extract mate 0-based position from this record.
    #[inline]
    #[must_use]
    pub fn mate_pos(&self) -> i32 {
        RawRecordView::new(&self.0).mate_pos()
    }

    /// Extract template length from this record.
    #[inline]
    #[must_use]
    pub fn template_length(&self) -> i32 {
        RawRecordView::new(&self.0).template_length()
    }

    // -- Flag-bit convenience reads --

    /// Returns `true` if the read is paired in sequencing.
    #[inline]
    #[must_use]
    pub fn is_paired(&self) -> bool {
        RawRecordView::new(&self.0).is_paired()
    }

    /// Returns `true` if the read is unmapped.
    #[inline]
    #[must_use]
    pub fn is_unmapped(&self) -> bool {
        RawRecordView::new(&self.0).is_unmapped()
    }

    /// Returns `true` if the mate is unmapped.
    #[inline]
    #[must_use]
    pub fn is_mate_unmapped(&self) -> bool {
        RawRecordView::new(&self.0).is_mate_unmapped()
    }

    /// Returns `true` if the read is reverse complemented.
    #[inline]
    #[must_use]
    pub fn is_reverse(&self) -> bool {
        RawRecordView::new(&self.0).is_reverse()
    }

    /// Returns `true` if the mate is reverse complemented.
    #[inline]
    #[must_use]
    pub fn is_mate_reverse(&self) -> bool {
        RawRecordView::new(&self.0).is_mate_reverse()
    }

    /// Returns `true` if the read is the first segment (R1).
    #[inline]
    #[must_use]
    pub fn is_first_segment(&self) -> bool {
        RawRecordView::new(&self.0).is_first_segment()
    }

    /// Returns `true` if the read is the last segment (R2).
    #[inline]
    #[must_use]
    pub fn is_last_segment(&self) -> bool {
        RawRecordView::new(&self.0).is_last_segment()
    }

    /// Returns `true` if the read is a secondary alignment.
    #[inline]
    #[must_use]
    pub fn is_secondary(&self) -> bool {
        RawRecordView::new(&self.0).is_secondary()
    }

    /// Returns `true` if the read fails quality controls.
    #[inline]
    #[must_use]
    pub fn is_qc_fail(&self) -> bool {
        RawRecordView::new(&self.0).is_qc_fail()
    }

    /// Returns `true` if the read is a PCR or optical duplicate.
    #[inline]
    #[must_use]
    pub fn is_duplicate(&self) -> bool {
        RawRecordView::new(&self.0).is_duplicate()
    }

    /// Returns `true` if the read is a supplementary alignment.
    #[inline]
    #[must_use]
    pub fn is_supplementary(&self) -> bool {
        RawRecordView::new(&self.0).is_supplementary()
    }

    // -- CIGAR reads --

    /// Compute the reference length consumed by this record's CIGAR.
    #[inline]
    #[must_use]
    pub fn reference_length(&self) -> i32 {
        RawRecordView::new(&self.0).reference_length()
    }

    /// Compute the query (read) length from this record's CIGAR.
    #[inline]
    #[must_use]
    pub fn query_length(&self) -> usize {
        RawRecordView::new(&self.0).query_length()
    }

    /// Compute reference-consuming and query-consuming lengths in one pass.
    ///
    /// Returns `(ref_length, query_length)`. See [`RawRecordView::cigar_lengths`].
    #[inline]
    #[must_use]
    pub fn cigar_lengths(&self) -> (i32, usize) {
        self.view().cigar_lengths()
    }

    /// Return the 1-based alignment end position, or `None` if unmapped.
    #[inline]
    #[must_use]
    pub fn alignment_end_1based(&self) -> Option<usize> {
        RawRecordView::new(&self.0).alignment_end_1based()
    }

    /// Return the 1-based alignment start position, or `None` if unmapped.
    #[inline]
    #[must_use]
    pub fn alignment_start_1based(&self) -> Option<usize> {
        RawRecordView::new(&self.0).alignment_start_1based()
    }

    /// Return the 1-based unclipped 5′ position.
    #[inline]
    #[must_use]
    pub fn unclipped_5prime_1based(&self) -> i32 {
        RawRecordView::new(&self.0).unclipped_5prime_1based()
    }

    /// Convert this record's CIGAR to a human-readable string.
    #[inline]
    #[must_use]
    pub fn cigar_to_string(&self) -> String {
        RawRecordView::new(&self.0).cigar_to_string()
    }

    /// Return the raw CIGAR ops as a `Vec<u32>`.
    #[inline]
    #[must_use]
    pub fn cigar_ops_vec(&self) -> Vec<u32> {
        RawRecordView::new(&self.0).cigar_ops_vec()
    }

    /// Zero-allocation iterator yielding raw CIGAR ops as `u32`.
    #[inline]
    pub fn cigar_ops_iter(&self) -> impl Iterator<Item = u32> + '_ {
        // Mirror `cigar_ops_typed` so the borrow is tied to &self.0 directly.
        use crate::fields::{l_read_name, n_cigar_op};
        let bam: &[u8] = &self.0;
        let lrn = l_read_name(bam) as usize;
        let n = n_cigar_op(bam) as usize;
        let start = 32 + lrn;
        let end = start + n * 4;
        let bytes = if end <= bam.len() { &bam[start..end] } else { &bam[..0] };
        bytes.chunks_exact(4).map(|c| u32::from_le_bytes([c[0], c[1], c[2], c[3]]))
    }

    /// Typed iteration over CIGAR ops. Each [`crate::cigar::CigarOp`] is a `Copy` newtype over
    /// the raw u32 — zero-cost ergonomic equivalent of noodles `cigar().iter()`.
    #[inline]
    pub fn cigar_ops_typed(&self) -> impl Iterator<Item = crate::cigar::CigarOp> + '_ {
        // Delegate via view() but keep lifetime tied to self's bytes, not the
        // temporary RawRecordView. We re-implement the same two-step here so
        // that the compiler can see that the borrow is from &self.0 directly.
        use crate::cigar::CigarOp;
        use crate::fields::{l_read_name, n_cigar_op};
        let bam: &[u8] = &self.0;
        let lrn = l_read_name(bam) as usize;
        let n = n_cigar_op(bam) as usize;
        let start = 32 + lrn;
        let end = start + n * 4;
        let bytes = if end <= bam.len() { &bam[start..end] } else { &bam[..0] };
        bytes
            .chunks_exact(4)
            .map(|c| CigarOp::from_raw(u32::from_le_bytes([c[0], c[1], c[2], c[3]])))
    }

    /// Decode the sequence bases to ASCII.
    #[inline]
    #[must_use]
    pub fn sequence_vec(&self) -> Vec<u8> {
        RawRecordView::new(&self.0).sequence_vec()
    }

    // -- Methods that return references borrowed from self --
    //    Hand-written so lifetime elision pins the returned reference to &self.

    /// Return the read name (without null terminator).
    #[inline]
    #[must_use]
    pub fn read_name(&self) -> &[u8] {
        self.view().read_name()
    }

    /// Return the raw CIGAR bytes.
    #[inline]
    #[must_use]
    pub fn cigar_raw_bytes(&self) -> &[u8] {
        self.view().cigar_raw_bytes()
    }

    /// Return the quality score bytes.
    #[inline]
    #[must_use]
    pub fn quality_scores(&self) -> &[u8] {
        self.view().quality_scores()
    }

    /// Return the packed (2-bit-encoded) sequence bytes.
    #[inline]
    #[must_use]
    pub fn sequence_packed(&self) -> &[u8] {
        self.view().sequence_packed()
    }

    /// Return the ASCII base at `position`.
    #[inline]
    #[must_use]
    pub fn get_base(&self, position: usize) -> u8 {
        self.view().get_base(position)
    }

    /// Return `true` if the base at `position` is N (missing/unknown).
    #[inline]
    #[must_use]
    pub fn is_base_n(&self, position: usize) -> bool {
        self.view().is_base_n(position)
    }

    /// Return the quality score at `position`.
    #[inline]
    #[must_use]
    pub fn get_qual(&self, position: usize) -> u8 {
        self.view().get_qual(position)
    }

    /// Return a mutable slice of quality scores.
    #[inline]
    pub fn quality_scores_mut(&mut self) -> &mut [u8] {
        crate::sequence::quality_scores_slice_mut(&mut self.0)
    }

    /// Set the ASCII base at `position`.
    #[inline]
    pub fn set_base(&mut self, position: usize, base: u8) {
        let seq_off = crate::fields::seq_offset(&self.0);
        crate::sequence::set_base(&mut self.0, seq_off, position, base);
    }

    /// Mask the base at `position` to N.
    #[inline]
    pub fn mask_base(&mut self, position: usize) {
        let seq_off = crate::fields::seq_offset(&self.0);
        crate::sequence::mask_base(&mut self.0, seq_off, position);
    }

    /// Set the quality score at `position`.
    #[inline]
    pub fn set_qual(&mut self, position: usize, value: u8) {
        let qual_off = crate::fields::qual_offset(&self.0);
        crate::sequence::set_qual(&mut self.0, qual_off, position, value);
    }

    // -- Header writes (fixed-length, no return) --

    /// Set the reference sequence ID.
    #[inline]
    pub fn set_ref_id(&mut self, v: i32) {
        RawRecordMut::new(&mut self.0).set_ref_id(v);
    }

    /// Set the 0-based leftmost position.
    #[inline]
    pub fn set_pos(&mut self, v: i32) {
        RawRecordMut::new(&mut self.0).set_pos(v);
    }

    /// Set the mapping quality.
    #[inline]
    pub fn set_mapq(&mut self, v: u8) {
        RawRecordMut::new(&mut self.0).set_mapq(v);
    }

    /// Set the bitwise flags.
    #[inline]
    pub fn set_flags(&mut self, v: u16) {
        RawRecordMut::new(&mut self.0).set_flags(v);
    }

    /// Set the BAM bin field.
    #[inline]
    pub fn set_bin(&mut self, v: u16) {
        RawRecordMut::new(&mut self.0).set_bin(v);
    }

    /// Set the mate reference sequence ID.
    #[inline]
    pub fn set_mate_ref_id(&mut self, v: i32) {
        RawRecordMut::new(&mut self.0).set_mate_ref_id(v);
    }

    /// Set the mate 0-based position.
    #[inline]
    pub fn set_mate_pos(&mut self, v: i32) {
        RawRecordMut::new(&mut self.0).set_mate_pos(v);
    }

    /// Set the template length (TLEN).
    #[inline]
    pub fn set_template_length(&mut self, v: i32) {
        RawRecordMut::new(&mut self.0).set_template_length(v);
    }

    // -- Length-changing edits --

    /// Replace the CIGAR operation array.
    ///
    /// Updates `n_cigar_op` and splices the new ops in. Each op is a `u32` in
    /// the standard BAM packing (`(op_len << 4) | op_type`).
    ///
    /// All other fields (sequence, quality, aux tags) are preserved.
    ///
    /// # Panics
    /// Panics if `new_ops.len() > u16::MAX as usize`.
    pub fn set_cigar_ops(&mut self, new_ops: &[u32]) {
        let new_n = u16::try_from(new_ops.len())
            .unwrap_or_else(|_| panic!("CIGAR op count {} overflows u16", new_ops.len()));
        let l_rn = self.0[8] as usize;
        let old_n = u16::from_le_bytes([self.0[12], self.0[13]]) as usize;
        let start = 32 + l_rn;
        let end = start + old_n * 4;
        let mut bytes = Vec::with_capacity(new_ops.len() * 4);
        for &op in new_ops {
            bytes.extend_from_slice(&op.to_le_bytes());
        }
        self.0.splice(start..end, bytes);
        self.0[12..14].copy_from_slice(&new_n.to_le_bytes());
    }

    /// Replace both the sequence and the quality scores.
    ///
    /// `seq` is ASCII (`A`/`C`/`G`/`T`/`N`/IUPAC); `qual` is raw Phred (not +33).
    /// `seq.len()` and `qual.len()` must match (BAM spec).
    ///
    /// Updates `l_seq`; read name, CIGAR, and aux tags are preserved.
    ///
    /// # Panics
    ///
    /// Panics if `seq.len() != qual.len()` or if `seq.len() > u32::MAX as usize`.
    pub fn set_sequence_and_qualities(&mut self, seq: &[u8], qual: &[u8]) {
        assert_eq!(seq.len(), qual.len(), "sequence and quality must have equal length");
        let new_l = u32::try_from(seq.len()).expect("seq length overflows u32");
        let l_rn = self.0[8] as usize;
        let n_co = u16::from_le_bytes([self.0[12], self.0[13]]) as usize;
        let old_l = u32::from_le_bytes([self.0[16], self.0[17], self.0[18], self.0[19]]) as usize;

        let body_start = 32 + l_rn + n_co * 4;
        let old_packed = old_l.div_ceil(2);
        let body_end = body_start + old_packed + old_l;
        let new_packed_len = seq.len().div_ceil(2);
        let new_body_len = new_packed_len + qual.len();
        let old_body_len = body_end - body_start;

        // Resize in place by shifting aux data to its new position, then write
        // packed seq and qual directly into the record buffer. This avoids the
        // intermediate Vec allocation and extra memcpy that the splice approach
        // required.
        if new_body_len > old_body_len {
            let grow = new_body_len - old_body_len;
            let old_len = self.0.len();
            self.0.resize(old_len + grow, 0);
            self.0.copy_within(body_end..old_len, body_end + grow);
        } else if new_body_len < old_body_len {
            let shrink = old_body_len - new_body_len;
            self.0.copy_within(body_end.., body_end - shrink);
            self.0.truncate(self.0.len() - shrink);
        }

        let (seq_slot, qual_slot) =
            self.0[body_start..body_start + new_body_len].split_at_mut(new_packed_len);
        crate::sequence::pack_sequence_into_slice(seq_slot, seq);
        qual_slot.copy_from_slice(qual);

        self.0[16..20].copy_from_slice(&new_l.to_le_bytes());
    }

    /// Convenience: replace quality scores when their length matches `l_seq`.
    ///
    /// For length-changing replacement, use [`RawRecord::set_sequence_and_qualities`].
    ///
    /// # Panics
    ///
    /// Panics if `qual.len() != l_seq`.
    pub fn set_qualities(&mut self, qual: &[u8]) {
        let l = self.l_seq() as usize;
        assert_eq!(
            l,
            qual.len(),
            "qualities length must match l_seq; use set_sequence_and_qualities to resize"
        );
        self.quality_scores_mut().copy_from_slice(qual);
    }

    /// Replace the read name. Updates `l_read_name` and splices the name +
    /// NUL terminator into the record.
    ///
    /// All other fields (CIGAR, sequence, quality, aux tags) are preserved
    /// because they live after the name and shift accordingly.
    ///
    /// # Panics
    ///
    /// Panics if `new_name.len() + 1 > u8::MAX as usize` (the BAM spec encodes
    /// `l_read_name` as a single byte), or if `new_name` contains an embedded
    /// NUL byte (the BAM spec terminates the read name with a single NUL, so
    /// embedded NULs would produce a truncated QNAME for downstream readers).
    pub fn set_read_name(&mut self, new_name: &[u8]) {
        assert!(!new_name.contains(&0), "read name must not contain embedded NUL bytes");
        let new_l = new_name.len() + 1; // include NUL
        let new_l_u8 = u8::try_from(new_l)
            .unwrap_or_else(|_| panic!("read name too long: {} bytes", new_name.len()));
        let old_l = self.0[8] as usize;
        let start = 32;
        let end = start + old_l;
        // Splice: replace old name+NUL with new name+NUL
        let mut replacement = Vec::with_capacity(new_l);
        replacement.extend_from_slice(new_name);
        replacement.push(0);
        self.0.splice(start..end, replacement);
        // Update l_read_name byte at offset 8
        self.0[8] = new_l_u8;
    }
}

/// A reader for raw BAM records.
///
/// This wraps any `Read` type (typically a BGZF reader) and provides
/// methods for reading raw BAM records as bytes.
pub struct RawBamReader<R> {
    inner: R,
}

impl<R: Read> RawBamReader<R> {
    /// Creates a new raw BAM reader wrapping the given reader.
    ///
    /// Note: The caller is responsible for reading/skipping the BAM header
    /// before reading records. Use `skip_header` or `read_header_with_noodles`
    /// to handle the header.
    #[inline]
    pub fn new(inner: R) -> Self {
        Self { inner }
    }

    /// Returns a reference to the inner reader.
    #[inline]
    pub fn get_ref(&self) -> &R {
        &self.inner
    }

    /// Returns a mutable reference to the inner reader.
    #[inline]
    pub fn get_mut(&mut self) -> &mut R {
        &mut self.inner
    }

    /// Consumes the reader and returns the inner reader.
    #[inline]
    pub fn into_inner(self) -> R {
        self.inner
    }

    /// Reads a single raw BAM record.
    ///
    /// Returns the number of bytes read, or 0 at EOF.
    ///
    /// # Errors
    ///
    /// Returns an error if reading fails or if EOF is encountered mid-record.
    #[inline]
    pub fn read_record(&mut self, record: &mut RawRecord) -> io::Result<usize> {
        read_raw_record(&mut self.inner, record)
    }
}

/// Write a single raw BAM record to `w`: 4-byte `block_size` (u32 LE) prefix
/// followed by the record bytes.
///
/// This is the inverse of [`read_raw_record`]. The caller is responsible for
/// writing the BAM header before any records.
///
/// # Errors
///
/// Returns an error if writing to `w` fails.
///
/// # Panics
///
/// Panics if `rec.len()` exceeds `u32::MAX`.
#[inline]
pub fn write_raw_record<W: Write>(w: &mut W, rec: &RawRecord) -> io::Result<()> {
    let block_size =
        u32::try_from(rec.len()).expect("record length exceeds u32 BAM block_size limit");
    w.write_all(&block_size.to_le_bytes())?;
    w.write_all(rec.as_ref())
}

/// Write a slice of raw BAM records to `w`. Equivalent to calling
/// [`write_raw_record`] for each record, but written as one helper.
///
/// # Errors
///
/// Returns on the first write error.
#[inline]
pub fn write_raw_records<W: Write>(w: &mut W, recs: &[RawRecord]) -> io::Result<()> {
    for rec in recs {
        write_raw_record(w, rec)?;
    }
    Ok(())
}

/// Thin writer wrapper that mirrors [`RawBamReader`]: handles the
/// `block_size` framing on writes.
///
/// Does NOT write BAM magic/header/ref-list — those are the caller's
/// responsibility (typically already produced by noodles' writer or by
/// copying from an input BAM).
pub struct RawBamWriter<W> {
    inner: W,
}

impl<W: Write> RawBamWriter<W> {
    /// Creates a new `RawBamWriter` wrapping the given writer.
    #[inline]
    pub fn new(inner: W) -> Self {
        Self { inner }
    }

    /// Returns a reference to the inner writer.
    #[inline]
    pub fn get_ref(&self) -> &W {
        &self.inner
    }

    /// Returns a mutable reference to the inner writer.
    #[inline]
    pub fn get_mut(&mut self) -> &mut W {
        &mut self.inner
    }

    /// Consumes the writer and returns the inner writer.
    #[inline]
    pub fn into_inner(self) -> W {
        self.inner
    }

    /// Write a single record with its `block_size` prefix.
    ///
    /// # Errors
    ///
    /// Returns an error if writing fails.
    #[inline]
    pub fn write_record(&mut self, rec: &RawRecord) -> io::Result<()> {
        write_raw_record(&mut self.inner, rec)
    }

    /// Write many records.
    ///
    /// # Errors
    ///
    /// Returns on the first write error.
    #[inline]
    pub fn write_records(&mut self, recs: &[RawRecord]) -> io::Result<()> {
        write_raw_records(&mut self.inner, recs)
    }

    /// Flush the underlying writer.
    ///
    /// # Errors
    ///
    /// Returns an error if flushing fails.
    #[inline]
    pub fn flush(&mut self) -> io::Result<()> {
        self.inner.flush()
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use fgumi_tag::SamTag;

    #[test]
    fn test_raw_record_inherent_delegators() {
        use crate::flags;
        use crate::testutil::*;
        let bytes =
            make_bam_bytes(3, 200, flags::PAIRED, b"rd", &[encode_op(0, 4)], 4, -1, -1, b"NMc\x05");
        let mut rec = RawRecord::from(bytes);
        // Reads via inherent delegators
        assert_eq!(rec.ref_id(), 3);
        assert_eq!(rec.pos(), 200);
        assert!(rec.is_paired());
        assert_eq!(rec.tags().find_int(SamTag::NM), Some(5));

        // Fixed-length writes via inherent delegators
        rec.set_pos(999);
        rec.set_mapq(40);
        assert_eq!(rec.pos(), 999);
        assert_eq!(rec.mapq(), 40);

        // Length-changing edits via tags_editor
        {
            let mut ed = rec.tags_editor();
            ed.update_int(SamTag::NM, 8);
        }
        assert_eq!(rec.tags().find_int(SamTag::NM), Some(8));
    }

    #[test]
    fn test_raw_record_new() {
        let record = RawRecord::new();
        assert!(record.is_empty());
        assert_eq!(record.len(), 0);
    }

    #[test]
    fn test_raw_record_as_ref() {
        let record = RawRecord::from(vec![1, 2, 3, 4]);
        assert_eq!(record.as_ref(), &[1, 2, 3, 4]);
    }

    #[test]
    fn test_raw_record_deref() {
        let record = RawRecord::from(vec![1, 2, 3, 4]);
        assert_eq!(&*record, &[1, 2, 3, 4]);
    }

    #[test]
    fn test_read_raw_record_success() {
        // block_size = 8, followed by 8 bytes of data
        let data = [
            0x08, 0x00, 0x00, 0x00, // block_size = 8
            0x01, 0x02, 0x03, 0x04, 0x05, 0x06, 0x07, 0x08, // record data
        ];
        let mut reader = &data[..];
        let mut record = RawRecord::new();

        let n = read_raw_record(&mut reader, &mut record)
            .expect("reading valid 8-byte record should succeed");
        assert_eq!(n, 8);
        assert_eq!(record.as_ref(), &[1, 2, 3, 4, 5, 6, 7, 8]);
    }

    #[test]
    fn test_read_raw_record_eof() {
        let data: &[u8] = &[];
        let mut reader = data;
        let mut record = RawRecord::new();

        let n = read_raw_record(&mut reader, &mut record)
            .expect("reading from empty input should return Ok(0)");
        assert_eq!(n, 0);
    }

    #[test]
    fn test_read_raw_record_truncated_size() {
        let data = [0x08, 0x00]; // Incomplete block_size
        let mut reader = &data[..];
        let mut record = RawRecord::new();

        let result = read_raw_record(&mut reader, &mut record);
        assert!(result.is_err());
    }

    #[test]
    fn test_read_raw_record_truncated_data() {
        let data = [
            0x08, 0x00, 0x00, 0x00, // block_size = 8
            0x01, 0x02, 0x03, // Only 3 bytes of data
        ];
        let mut reader = &data[..];
        let mut record = RawRecord::new();

        let result = read_raw_record(&mut reader, &mut record);
        assert!(result.is_err());
    }

    #[test]
    fn test_set_read_name_resize() {
        use crate::testutil::*;
        let bytes = make_bam_bytes(0, 0, 0, b"oldname", &[encode_op(0, 4)], 4, -1, -1, b"NMc\x05");
        let mut rec = RawRecord::from(bytes);
        let pre_cigar = rec.cigar_ops_vec();
        let pre_seq = rec.sequence_vec();
        let pre_qual = rec.quality_scores().to_vec();
        let pre_nm = rec.tags().find_int(SamTag::NM);

        rec.set_read_name(b"new");
        assert_eq!(rec.read_name(), b"new");
        assert_eq!(rec.l_read_name(), 4); // 3 bytes + NUL
        // All other fields preserved
        assert_eq!(rec.cigar_ops_vec(), pre_cigar);
        assert_eq!(rec.sequence_vec(), pre_seq);
        assert_eq!(rec.quality_scores(), pre_qual.as_slice());
        assert_eq!(rec.tags().find_int(SamTag::NM), pre_nm);

        // Empty name (still NUL-terminated)
        rec.set_read_name(b"");
        assert_eq!(rec.read_name(), b"");
        assert_eq!(rec.l_read_name(), 1);
    }

    #[test]
    fn test_set_cigar_ops_resize() {
        use crate::testutil::*;
        let bytes = make_bam_bytes(0, 0, 0, b"r", &[encode_op(0, 10)], 10, -1, -1, b"NMc\x05");
        let mut rec = RawRecord::from(bytes);
        let pre_name = rec.read_name().to_vec();
        let pre_seq = rec.sequence_vec();
        let pre_qual = rec.quality_scores().to_vec();
        let pre_nm = rec.tags().find_int(SamTag::NM);

        let new_ops = vec![encode_op(0, 4), encode_op(2, 1), encode_op(0, 6)];
        rec.set_cigar_ops(&new_ops);
        assert_eq!(rec.cigar_ops_vec(), new_ops);
        assert_eq!(rec.n_cigar_op() as usize, new_ops.len());

        assert_eq!(rec.read_name(), pre_name.as_slice());
        assert_eq!(rec.sequence_vec(), pre_seq);
        assert_eq!(rec.quality_scores(), pre_qual.as_slice());
        assert_eq!(rec.tags().find_int(SamTag::NM), pre_nm);
    }

    #[test]
    fn test_set_sequence_and_qualities_grow() {
        use crate::testutil::*;
        let bytes = make_bam_bytes(0, 0, 0, b"r", &[encode_op(0, 4)], 4, -1, -1, b"NMc\x05");
        let mut rec = RawRecord::from(bytes);
        let pre_name = rec.read_name().to_vec();
        let pre_cigar = rec.cigar_ops_vec();
        let pre_nm = rec.tags().find_int(SamTag::NM);

        let new_seq = b"ACGTACGTACGT";
        let new_qual = vec![30u8; 12];
        rec.set_sequence_and_qualities(new_seq, &new_qual);

        assert_eq!(rec.l_seq() as usize, 12);
        assert_eq!(rec.sequence_vec(), new_seq.to_vec());
        assert_eq!(rec.quality_scores(), new_qual.as_slice());
        assert_eq!(rec.read_name(), pre_name.as_slice());
        assert_eq!(rec.cigar_ops_vec(), pre_cigar);
        assert_eq!(rec.tags().find_int(SamTag::NM), pre_nm);
    }

    #[test]
    fn test_set_sequence_and_qualities_shrink_with_odd_length() {
        use crate::testutil::*;
        let bytes = make_bam_bytes(0, 0, 0, b"r", &[encode_op(0, 4)], 4, -1, -1, &[]);
        let mut rec = RawRecord::from(bytes);
        rec.set_sequence_and_qualities(b"ACG", &[20, 21, 22]);
        assert_eq!(rec.l_seq(), 3);
        assert_eq!(rec.sequence_vec(), b"ACG".to_vec());
        assert_eq!(rec.quality_scores(), &[20u8, 21, 22]);
    }

    #[test]
    #[should_panic(expected = "must have equal length")]
    fn test_set_sequence_and_qualities_mismatch_panics() {
        use crate::testutil::*;
        let mut rec = RawRecord::from(make_bam_bytes(0, 0, 0, b"r", &[], 0, -1, -1, &[]));
        rec.set_sequence_and_qualities(b"AC", &[20]);
    }

    #[test]
    fn test_set_qualities_fixed_length() {
        use crate::testutil::*;
        let bytes = make_bam_bytes(0, 0, 0, b"r", &[encode_op(0, 4)], 4, -1, -1, &[]);
        let mut rec = RawRecord::from(bytes);
        rec.set_qualities(&[10, 20, 30, 40]);
        assert_eq!(rec.quality_scores(), &[10u8, 20, 30, 40]);
    }

    #[test]
    #[should_panic(expected = "qualities length must match l_seq")]
    fn test_set_qualities_wrong_length_panics() {
        use crate::testutil::*;
        let bytes = make_bam_bytes(0, 0, 0, b"r", &[encode_op(0, 4)], 4, -1, -1, &[]);
        let mut rec = RawRecord::from(bytes);
        rec.set_qualities(&[10, 20]);
    }

    #[test]
    fn test_read_multiple_records() {
        let data = [
            0x04, 0x00, 0x00, 0x00, // block_size = 4
            0x01, 0x02, 0x03, 0x04, // record 1
            0x02, 0x00, 0x00, 0x00, // block_size = 2
            0x05, 0x06, // record 2
        ];
        let mut reader = &data[..];
        let mut record = RawRecord::new();

        let n =
            read_raw_record(&mut reader, &mut record).expect("reading first record should succeed");
        assert_eq!(n, 4);
        assert_eq!(record.as_ref(), &[1, 2, 3, 4]);

        let n = read_raw_record(&mut reader, &mut record)
            .expect("reading second record should succeed");
        assert_eq!(n, 2);
        assert_eq!(record.as_ref(), &[5, 6]);

        let n =
            read_raw_record(&mut reader, &mut record).expect("reading at EOF should return Ok(0)");
        assert_eq!(n, 0); // EOF
    }

    // ── Edge-case tests for length-changing edits ──────────────────────────

    /// A 254-byte name (+ NUL = 255 = `u8::MAX`) is the longest name the BAM
    /// spec allows; `set_read_name` must accept it without panic.
    #[test]
    fn test_set_read_name_max_length_254_bytes_succeeds() {
        use crate::testutil::*;
        let bytes = make_bam_bytes(0, 0, 0, b"r", &[encode_op(0, 4)], 4, -1, -1, &[]);
        let mut rec = RawRecord::from(bytes);
        let max_name = vec![b'A'; 254]; // + NUL = 255 = u8::MAX
        rec.set_read_name(&max_name);
        assert_eq!(rec.l_read_name(), 255);
        assert_eq!(rec.read_name().len(), 254);
        assert_eq!(rec.read_name(), max_name.as_slice());
    }

    /// A 255-byte name (+ NUL = 256) overflows `u8`; `set_read_name` must panic.
    #[test]
    #[should_panic(expected = "read name too long")]
    fn test_set_read_name_over_limit_panics() {
        use crate::testutil::*;
        let bytes = make_bam_bytes(0, 0, 0, b"r", &[encode_op(0, 4)], 4, -1, -1, &[]);
        let mut rec = RawRecord::from(bytes);
        let too_long = vec![b'A'; 255]; // + NUL = 256 > u8::MAX
        rec.set_read_name(&too_long);
    }

    /// An empty name is valid: `l_read_name` should be 1 (just the NUL
    /// terminator) and `read_name()` should return an empty slice.
    #[test]
    fn test_set_read_name_empty() {
        use crate::testutil::*;
        let bytes = make_bam_bytes(0, 0, 0, b"somename", &[encode_op(0, 4)], 4, -1, -1, &[]);
        let mut rec = RawRecord::from(bytes);
        rec.set_read_name(b"");
        assert_eq!(rec.l_read_name(), 1); // just the NUL
        assert_eq!(rec.read_name(), b"");
    }

    /// Embedded NULs would give downstream C-string-style readers a truncated
    /// QNAME while `l_read_name` still spans the full field; reject them.
    #[test]
    #[should_panic(expected = "embedded NUL")]
    fn test_set_read_name_rejects_embedded_nul() {
        use crate::testutil::*;
        let bytes = make_bam_bytes(0, 0, 0, b"r", &[encode_op(0, 4)], 4, -1, -1, &[]);
        let mut rec = RawRecord::from(bytes);
        rec.set_read_name(b"a\0b");
    }

    /// `set_cigar_ops(&[])` clears the CIGAR entirely: `n_cigar_op` must be 0
    /// and aux tags must survive the splice unchanged.
    #[test]
    fn test_set_cigar_ops_empty_clears_cigar() {
        use crate::testutil::*;
        let bytes = make_bam_bytes(
            0,
            0,
            0,
            b"r",
            &[encode_op(0, 10), encode_op(2, 1), encode_op(0, 5)],
            15,
            -1,
            -1,
            b"NMc\x05",
        );
        let mut rec = RawRecord::from(bytes);
        let pre_nm = rec.tags().find_int(SamTag::NM);

        rec.set_cigar_ops(&[]);

        assert_eq!(rec.n_cigar_op(), 0);
        assert_eq!(rec.cigar_ops_vec(), Vec::<u32>::new());
        // Aux tags survive the clear.
        assert_eq!(rec.tags().find_int(SamTag::NM), pre_nm);
    }

    /// `set_sequence_and_qualities(b"", &[])` is legal: `l_seq` becomes 0, the
    /// packed-sequence and quality slices are empty, and aux tags are preserved.
    #[test]
    fn test_set_sequence_and_qualities_zero_length() {
        use crate::testutil::*;
        let bytes = make_bam_bytes(0, 0, 0, b"r", &[encode_op(0, 4)], 4, -1, -1, b"NMc\x05");
        let mut rec = RawRecord::from(bytes);
        let pre_nm = rec.tags().find_int(SamTag::NM);

        rec.set_sequence_and_qualities(b"", &[]);

        assert_eq!(rec.l_seq(), 0);
        assert_eq!(rec.sequence_vec(), Vec::<u8>::new());
        assert_eq!(rec.quality_scores(), &[] as &[u8]);
        // Aux tags survive the removal of all seq+qual bytes.
        assert_eq!(rec.tags().find_int(SamTag::NM), pre_nm);
    }

    /// A `B:S` (u16) array tag followed by an int tag must be byte-for-byte
    /// identical after a seq+qual resize (grow from 10 to 16 bases).
    #[test]
    fn test_set_sequence_and_qualities_preserves_array_aux_tags() {
        use crate::testutil::*;
        // Build aux: B:S array of 4 u16 values, then NM:i:7
        let mut aux = Vec::new();
        aux.extend_from_slice(b"bqBS");
        aux.extend_from_slice(&4u32.to_le_bytes());
        for v in [100u16, 200, 300, 400] {
            aux.extend_from_slice(&v.to_le_bytes());
        }
        aux.extend_from_slice(b"NMc\x07");

        let bytes = make_bam_bytes(0, 0, 0, b"r", &[encode_op(0, 10)], 10, -1, -1, &aux);
        let mut rec = RawRecord::from(bytes);
        let pre_aux = rec.tags().as_bytes().to_vec();

        // Grow seq from 10 to 16 bases.
        rec.set_sequence_and_qualities(b"ACGTACGTACGTACGT", &[30u8; 16]);

        // Array tag must decode identically.
        let post_arr =
            rec.tags().find_array(SamTag::BQ).expect("B:S array tag should be preserved");
        assert_eq!(post_arr.count, 4);
        let decoded: Vec<u16> = (0..4)
            .map(|i| u16::from_le_bytes([post_arr.data[i * 2], post_arr.data[i * 2 + 1]]))
            .collect();
        assert_eq!(decoded, vec![100u16, 200, 300, 400]);

        assert_eq!(rec.tags().find_int(SamTag::NM), Some(7));
        // Byte-for-byte preservation of the entire aux section.
        assert_eq!(rec.tags().as_bytes(), pre_aux.as_slice());
    }

    /// A `B:i` (i32) array tag must survive a read-name grow (length-changing
    /// edit in the variable-length section before seq/qual/aux).
    #[test]
    fn test_set_read_name_with_array_aux_tags() {
        use crate::testutil::*;
        // Build aux: B:i array of 3 i32 values.
        let mut aux = Vec::new();
        aux.extend_from_slice(b"bqBi");
        aux.extend_from_slice(&3u32.to_le_bytes());
        for v in [1_000_000i32, -2_000_000, 3_000_000] {
            aux.extend_from_slice(&v.to_le_bytes());
        }

        let bytes = make_bam_bytes(0, 0, 0, b"oldname", &[encode_op(0, 4)], 4, -1, -1, &aux);
        let mut rec = RawRecord::from(bytes);
        let pre_aux = rec.tags().as_bytes().to_vec();

        // Grow the read name so the splice shifts aux.
        rec.set_read_name(b"new_longer_name_xyz");

        assert_eq!(rec.read_name(), b"new_longer_name_xyz");
        // Aux bytes must be byte-for-byte identical after the shift.
        assert_eq!(rec.tags().as_bytes(), pre_aux.as_slice());
    }

    // ── write_raw_record / write_raw_records / RawBamWriter tests ─────────

    /// Write a single record with `write_raw_record`, read it back with
    /// `read_raw_record`, and assert byte-identical round-trip.
    #[test]
    fn test_write_raw_record_round_trips_via_read_raw_record() {
        let original = RawRecord::from(vec![0x01, 0x02, 0x03, 0x04, 0x05]);
        let mut buf: Vec<u8> = Vec::new();
        write_raw_record(&mut buf, &original).expect("write should succeed");

        // Expected bytes: 4-byte LE block_size (5) + the 5 data bytes
        assert_eq!(&buf[..4], &5u32.to_le_bytes());
        assert_eq!(&buf[4..], original.as_ref());

        // Round-trip via read_raw_record
        let mut cursor = std::io::Cursor::new(&buf);
        let mut recovered = RawRecord::new();
        let n = read_raw_record(&mut cursor, &mut recovered).expect("read should succeed");
        assert_eq!(n, 5);
        assert_eq!(recovered, original);
    }

    /// Write 3 records with `write_raw_records`, read them back in order.
    #[test]
    fn test_write_raw_records_round_trips() {
        let recs = vec![
            RawRecord::from(vec![0xAA, 0xBB]),
            RawRecord::from(vec![0x11, 0x22, 0x33]),
            RawRecord::from(vec![0xFF]),
        ];
        let mut buf: Vec<u8> = Vec::new();
        write_raw_records(&mut buf, &recs).expect("write_raw_records should succeed");

        let mut cursor = std::io::Cursor::new(&buf);
        for expected in &recs {
            let mut got = RawRecord::new();
            let n = read_raw_record(&mut cursor, &mut got).expect("read should succeed");
            assert_eq!(n, expected.len());
            assert_eq!(&got, expected);
        }
        // No more records.
        let mut got = RawRecord::new();
        let n = read_raw_record(&mut cursor, &mut got).expect("EOF read should return Ok(0)");
        assert_eq!(n, 0);
    }

    /// `RawBamWriter` round-trips a single record the same way the free
    /// function does.
    #[test]
    fn test_raw_bam_writer_round_trips() {
        let original = RawRecord::from(vec![0xDE, 0xAD, 0xBE, 0xEF]);
        let mut writer = RawBamWriter::new(Vec::<u8>::new());
        writer.write_record(&original).expect("write_record should succeed");
        writer.flush().expect("flush should succeed");
        let buf = writer.into_inner();

        assert_eq!(&buf[..4], &4u32.to_le_bytes());
        assert_eq!(&buf[4..], original.as_ref());

        let mut cursor = std::io::Cursor::new(&buf);
        let mut recovered = RawRecord::new();
        let n = read_raw_record(&mut cursor, &mut recovered).expect("read should succeed");
        assert_eq!(n, 4);
        assert_eq!(recovered, original);
    }
}
