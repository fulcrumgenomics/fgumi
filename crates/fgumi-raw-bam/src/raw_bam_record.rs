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

use std::io::{self, Read};

/// A raw BAM record stored as bytes.
///
/// This is a zero-overhead wrapper that provides `AsRef<[u8]>` access to the
/// underlying BAM record bytes, enabling high-performance field extraction
/// and direct output writes.
#[derive(Clone, Default, Eq, PartialEq)]
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
        // Clamp against truncated/corrupt records where the parsed offset can
        // exceed the buffer length; falls back to an empty aux view rather than panicking.
        let off = crate::fields::aux_data_offset_from_record(&self.0)
            .filter(|&off| off <= len)
            .unwrap_or(len);
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

    /// Replace the read name. Updates `l_read_name` and splices the name +
    /// NUL terminator into the record.
    ///
    /// All other fields (CIGAR, sequence, quality, aux tags) are preserved
    /// because they live after the name and shift accordingly.
    ///
    /// # Panics
    ///
    /// Panics if `new_name.len() + 1 > u8::MAX as usize` (the BAM spec encodes
    /// `l_read_name` as a single byte) or if `new_name` contains an embedded NUL
    /// (BAM read names are NUL-terminated; an interior NUL would produce an
    /// invalid record that downstream readers misparse).
    pub fn set_read_name(&mut self, new_name: &[u8]) {
        assert!(!new_name.contains(&0), "read name must not contain NUL bytes");
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

#[cfg(test)]
mod tests {
    use super::*;

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
        assert_eq!(rec.tags().find_int(b"NM"), Some(5));

        // Fixed-length writes via inherent delegators
        rec.set_pos(999);
        rec.set_mapq(40);
        assert_eq!(rec.pos(), 999);
        assert_eq!(rec.mapq(), 40);

        // Length-changing edits via tags_editor
        {
            let mut ed = rec.tags_editor();
            ed.update_int(b"NM", 8);
        }
        assert_eq!(rec.tags().find_int(b"NM"), Some(8));
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
    fn test_tags_mut_clamps_truncated_aux_offset() {
        // Build a record then truncate below the parsed aux offset. The header
        // fields still imply aux starts past the buffer end, so tags_mut() must
        // clamp rather than panic on the out-of-bounds slice.
        use crate::testutil::*;
        let bytes = make_bam_bytes(0, 0, 0, b"r", &[encode_op(0, 4)], 4, -1, -1, b"NMc\x05");
        let mut truncated: Vec<u8> = bytes;
        // Drop the aux bytes plus one byte of the seq/qual region so the parsed
        // aux offset exceeds the buffer length.
        truncated.truncate(truncated.len() - 5);
        let mut rec = RawRecord::from(truncated);
        let tags = rec.tags_mut(); // must not panic
        assert!(tags.view().is_empty());
    }

    #[test]
    fn test_set_read_name_resize() {
        use crate::testutil::*;
        let bytes = make_bam_bytes(0, 0, 0, b"oldname", &[encode_op(0, 4)], 4, -1, -1, b"NMc\x05");
        let mut rec = RawRecord::from(bytes);
        let pre_cigar = rec.cigar_ops_vec();
        let pre_seq = rec.sequence_vec();
        let pre_qual = rec.quality_scores().to_vec();
        let pre_nm = rec.tags().find_int(b"NM");

        rec.set_read_name(b"new");
        assert_eq!(rec.read_name(), b"new");
        assert_eq!(rec.l_read_name(), 4); // 3 bytes + NUL
        // All other fields preserved
        assert_eq!(rec.cigar_ops_vec(), pre_cigar);
        assert_eq!(rec.sequence_vec(), pre_seq);
        assert_eq!(rec.quality_scores(), pre_qual.as_slice());
        assert_eq!(rec.tags().find_int(b"NM"), pre_nm);

        // Empty name (still NUL-terminated)
        rec.set_read_name(b"");
        assert_eq!(rec.read_name(), b"");
        assert_eq!(rec.l_read_name(), 1);
    }

    #[test]
    #[should_panic(expected = "read name must not contain NUL bytes")]
    fn test_set_read_name_rejects_embedded_nul() {
        use crate::testutil::*;
        let bytes = make_bam_bytes(0, 0, 0, b"r", &[encode_op(0, 4)], 4, -1, -1, b"");
        let mut rec = RawRecord::from(bytes);
        rec.set_read_name(b"bad\0name");
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
}
