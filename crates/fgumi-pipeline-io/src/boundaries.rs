//! `BoundaryState` / `BoundaryBatch`: the BAM record-boundary state machine.
//!
//! Lives in `fgumi-pipeline-io` so both `FindBamBoundaries` (in the `fgumi`
//! crate, via a re-export shim) and the fused `SortBuffer` ingest can share it.
//!
//! Scans decompressed BGZF block data for BAM record boundaries (header skip,
//! cross-block record carryover, EOF validation) without decoding records.
//! Driven by the `FindBamBoundaries` step (`super::bam`).
//!
//! Relocated from the legacy `bam.rs` (deleted in the issue #330 migration);
//! the boundary-finding logic is reused verbatim so the new framework's
//! boundary semantics match the legacy pipeline's exactly.

use std::io;

/// Output of `FindBoundaries` step: buffer + record offsets for parallel decoding.
///
/// This struct enables parallel BAM record decoding by pre-computing where
/// each record starts in the decompressed data. The actual parsing/decoding
/// can then be parallelized across multiple threads.
#[derive(Debug, Clone)]
pub struct BoundaryBatch {
    /// The decompressed bytes (with leftover prepended, suffix removed).
    pub buffer: Vec<u8>,
    /// Byte offsets where each record starts (offsets into buffer).
    /// Length = `num_records` + 1 (last entry is `buffer.len()` for easy slicing).
    pub offsets: Vec<usize>,
}

/// State for the `FindBoundaries` step (sequential).
///
/// This state maintains leftover bytes from incomplete records that span
/// across BGZF block boundaries. The boundary finding is very fast (~0.1μs
/// per block) since it only reads 4-byte integers without decoding records.
///
/// Uses a reusable work buffer to minimize allocations on the hot path.
pub struct BoundaryState {
    /// Leftover bytes from previous block (incomplete record at end).
    leftover: Vec<u8>,
    /// Reusable working buffer to avoid per-call allocations.
    work_buffer: Vec<u8>,
    /// Whether the BAM header has been skipped.
    header_skipped: bool,
    /// Length of the previous call's `offsets` Vec, used to pre-size the next
    /// one. Adjacent BGZF blocks hold near-identical record counts, so this
    /// collapses the per-block push-regrowth (~8 reallocations) to ~1. The
    /// returned `offsets` Vec is moved into `BoundaryBatch`, so it cannot be a
    /// reused buffer; pre-sizing is the cheap, correctness-neutral alternative.
    prev_offsets_len: usize,
}

/// Return the byte length of the BAM header at the start of `data`.
///
/// The BAM header consists of:
/// - 4-byte magic (`BAM\x01`)
/// - 4-byte `l_text` (little-endian u32)
/// - `l_text` bytes of plain-text header
/// - 4-byte `n_ref` (little-endian u32)
/// - for each of the `n_ref` references: 4-byte `l_name` + `l_name` bytes of name + 4-byte `l_ref`
///
/// Returns:
/// - `Ok(Some(offset))` — `offset` bytes consume the complete header; the first record starts there.
/// - `Ok(None)` — `data` is too short to contain a complete header; more bytes are needed.
/// - `Err(InvalidData)` — `data` does not begin with the BAM magic. This path skips a BAM header,
///   so a wrong magic means the stream is not what the caller declared; fail closed rather than
///   silently treating arbitrary bytes as headerless records. Genuinely headerless streams must use
///   [`BoundaryState::new_no_header`], which never calls this.
///
/// # Errors
///
/// Returns `InvalidData` when `data` is long enough to check the magic but does not start with it.
pub fn bam_header_len(data: &[u8]) -> io::Result<Option<usize>> {
    // BAM header structure:
    // - magic: 4 bytes ("BAM\1")
    // - l_text: 4 bytes (header text length)
    // - text: l_text bytes
    // - n_ref: 4 bytes (number of references)
    // - for each reference:
    //   - l_name: 4 bytes
    //   - name: l_name bytes
    //   - l_ref: 4 bytes

    if data.len() < 8 {
        return Ok(None);
    }

    // Check magic. This function is only reached on the header-skipping path
    // (`header_skipped == false`); a wrong magic there means the stream is not
    // the BAM the caller declared, so fail closed instead of misinterpreting the
    // bytes as headerless records.
    if &data[0..4] != fgumi_raw_bam::BAM_MAGIC {
        return Err(io::Error::new(io::ErrorKind::InvalidData, "invalid BAM magic"));
    }

    let l_text = u32::from_le_bytes([data[4], data[5], data[6], data[7]]) as usize;
    let mut offset = 8 + l_text;

    if data.len() < offset + 4 {
        return Ok(None);
    }

    let n_ref =
        u32::from_le_bytes([data[offset], data[offset + 1], data[offset + 2], data[offset + 3]])
            as usize;
    offset += 4;

    // Parse each reference
    for _ in 0..n_ref {
        if data.len() < offset + 4 {
            return Ok(None);
        }
        let l_name = u32::from_le_bytes([
            data[offset],
            data[offset + 1],
            data[offset + 2],
            data[offset + 3],
        ]) as usize;
        offset += 4 + l_name + 4; // l_name + name + l_ref

        if data.len() < offset {
            return Ok(None);
        }
    }

    Ok(Some(offset))
}

impl BoundaryState {
    /// Create a new boundary state.
    #[must_use]
    pub fn new() -> Self {
        Self {
            leftover: Vec::new(),
            work_buffer: Vec::new(),
            header_skipped: false,
            prev_offsets_len: 0,
        }
    }

    /// Create a new boundary state that doesn't skip the header.
    /// Use this when the input stream is already positioned past the header.
    #[must_use]
    pub fn new_no_header() -> Self {
        Self {
            leftover: Vec::new(),
            work_buffer: Vec::new(),
            header_skipped: true,
            prev_offsets_len: 0,
        }
    }

    /// Find record boundaries in decompressed data.
    ///
    /// This is FAST (~0.1μs per block) because it only scans 4-byte integers
    /// to find where records start - no actual record decoding is performed.
    ///
    /// # Arguments
    ///
    /// * `decompressed` - Decompressed bytes from one or more BGZF blocks
    ///
    /// # Returns
    ///
    /// A `BoundaryBatch` containing the complete records and their offsets.
    /// Any incomplete record at the end is saved as leftover for the next call.
    ///
    /// # Errors
    ///
    /// Returns an I/O error if the BAM header is malformed.
    ///
    /// # Record-level validation
    ///
    /// This function does NOT validate individual record `block_size` values
    /// against a malformed (but self-consistent) BAM stream. The per-record
    /// cross-check below (offset delta vs. the stored prefix) is a
    /// `debug_assertions`-only regression tripwire for this scanner's own
    /// arithmetic — it re-reads the same `block_size` bytes the scan already
    /// trusted, so it can only catch an internal bookkeeping bug, never input
    /// corruption. Authoritative release-build validation of record structure
    /// (out-of-bounds record end, trailing partial record) is performed
    /// downstream by `parse_records` / `parse_record_ranges` on the same bytes,
    /// which hard-error in all build modes. The `offsets` vector this returns
    /// is not consumed in release builds (`FindBamBoundaries` forwards only
    /// `buffer`), so promoting the cross-check to release would re-validate a
    /// tautology at a per-record cost for no correctness benefit.
    pub fn find_boundaries(&mut self, decompressed: &[u8]) -> io::Result<BoundaryBatch> {
        let (offsets, range) = self.scan(decompressed)?;
        // Owning copy of the complete records — for callers that need an owned
        // buffer (`FindBamBoundaries`). The zero-copy ingest path
        // ([`scan`](Self::scan) + [`records_bytes`](Self::records_bytes)) skips
        // this allocation + copy entirely.
        let buffer = self.work_buffer[range].to_vec();

        // Debug-only regression tripwire (NOT input validation): cross-check
        // each record's stored block_size prefix against the offset delta this
        // scan just computed. Both derive from the same bytes with no
        // intervening mutation, so this only catches an internal arithmetic /
        // indexing bug in the scan above — a corrupt-but-self-consistent
        // block_size passes trivially. Authoritative release validation lives
        // in parse_records / parse_record_ranges downstream (see the
        // `find_boundaries` doc comment).
        #[cfg(debug_assertions)]
        for i in 0..offsets.len().saturating_sub(1) {
            let start = offsets[i];
            let end = offsets[i + 1];
            if end > start + 4 {
                let stored = u32::from_le_bytes([
                    buffer[start],
                    buffer[start + 1],
                    buffer[start + 2],
                    buffer[start + 3],
                ]) as usize;
                let expected = end - start - 4;
                debug_assert_eq!(
                    stored, expected,
                    "find_boundaries: block_size mismatch at record {i}: stored={stored}, expected={expected}"
                );
            }
        }

        Ok(BoundaryBatch { buffer, offsets })
    }

    /// Zero-copy core of [`find_boundaries`](Self::find_boundaries): combine
    /// leftover + `decompressed` into the reusable `work_buffer`, skip the header
    /// (first call), scan record boundaries, and stash the trailing partial
    /// record as leftover — **without** copying the complete records out.
    ///
    /// Returns `(offsets, range)` where the complete records live in
    /// `self.work_buffer[range]` (valid until the next `scan`/`find_boundaries`
    /// call) and `offsets[i] .. offsets[i+1]` slices record `i` *relative to
    /// `range.start`* (so record `i`'s bytes are
    /// `records_bytes()[offsets[i] .. offsets[i+1]]`). The caller must consume
    /// the records before the next call. A header-only / incomplete-header block
    /// yields `(vec![0], 0..0)`.
    ///
    /// # Errors
    ///
    /// Returns an I/O error if the BAM header is malformed.
    pub fn scan(
        &mut self,
        decompressed: &[u8],
    ) -> io::Result<(Vec<usize>, std::ops::Range<usize>)> {
        // Step 1: Combine leftover with new data into reusable work_buffer.
        self.work_buffer.clear();
        if !self.leftover.is_empty() {
            self.work_buffer.append(&mut self.leftover);
        }
        self.work_buffer.extend_from_slice(decompressed);

        // Step 2: Skip header if not already done.
        let mut cursor = 0usize;
        if !self.header_skipped {
            let Some(header_size) = bam_header_len(&self.work_buffer)? else {
                // Not enough data to parse header; save as leftover, empty range.
                std::mem::swap(&mut self.leftover, &mut self.work_buffer);
                return Ok((vec![0], 0..0));
            };
            cursor = header_size;
            self.header_skipped = true;
        }

        // Step 3: Scan for record boundaries (FAST - just read integers).
        let start_cursor = cursor;
        let mut offsets = Vec::with_capacity(self.prev_offsets_len.max(1));
        offsets.push(0usize);
        while cursor + 4 <= self.work_buffer.len() {
            let block_size = u32::from_le_bytes([
                self.work_buffer[cursor],
                self.work_buffer[cursor + 1],
                self.work_buffer[cursor + 2],
                self.work_buffer[cursor + 3],
            ]) as usize;
            let record_end = cursor + 4 + block_size;
            if record_end > self.work_buffer.len() {
                break; // Incomplete record - becomes leftover.
            }
            cursor = record_end;
            offsets.push(cursor - start_cursor);
        }
        self.prev_offsets_len = offsets.len();

        // Step 4: Save trailing partial record as leftover (small — at most one
        // record). The complete records stay in `work_buffer[start_cursor..cursor]`
        // for the caller to read borrowed, with no copy.
        self.leftover.clear();
        self.leftover.extend_from_slice(&self.work_buffer[cursor..]);

        Ok((offsets, start_cursor..cursor))
    }

    /// Borrow the complete records produced by the most recent [`scan`](Self::scan),
    /// given the `range` it returned. Valid until the next `scan`/`find_boundaries`.
    #[must_use]
    pub fn records_bytes(&self, range: std::ops::Range<usize>) -> &[u8] {
        &self.work_buffer[range]
    }

    /// Call at EOF to get any remaining leftover.
    ///
    /// This validates that any remaining bytes form complete records.
    /// If there are incomplete bytes at EOF, an error is returned.
    ///
    /// # Errors
    ///
    /// Returns an I/O error if there are incomplete BAM records at EOF.
    pub fn finish(&mut self) -> io::Result<Option<BoundaryBatch>> {
        if self.leftover.is_empty() {
            return Ok(None);
        }

        // Try to parse remaining leftover
        let mut offsets = vec![0usize];
        let mut cursor = 0usize;

        while cursor + 4 <= self.leftover.len() {
            let block_size = u32::from_le_bytes([
                self.leftover[cursor],
                self.leftover[cursor + 1],
                self.leftover[cursor + 2],
                self.leftover[cursor + 3],
            ]) as usize;

            let record_end = cursor + 4 + block_size;
            if record_end > self.leftover.len() {
                return Err(io::Error::new(
                    io::ErrorKind::UnexpectedEof,
                    format!(
                        "Incomplete BAM record at EOF: need {} bytes, have {}",
                        record_end - cursor,
                        self.leftover.len() - cursor
                    ),
                ));
            }

            cursor = record_end;
            offsets.push(cursor);
        }

        // The loop only advances `cursor` by whole records. If it stops with
        // bytes still unconsumed (`cursor < leftover.len()`), those 1-3 trailing
        // bytes are too short to even hold a 4-byte block-size prefix — i.e. a
        // truncated BAM record. Surface it as an error rather than dropping the
        // bytes and masking corruption.
        if cursor < self.leftover.len() {
            return Err(io::Error::new(
                io::ErrorKind::UnexpectedEof,
                format!(
                    "Incomplete BAM record at EOF: {} trailing byte(s) cannot form a complete record",
                    self.leftover.len() - cursor
                ),
            ));
        }

        Ok(Some(BoundaryBatch { buffer: std::mem::take(&mut self.leftover), offsets }))
    }
}

impl Default for BoundaryState {
    fn default() -> Self {
        Self::new()
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn bam_header_len_matches_constructed_header() {
        // Minimal BAM header: magic + l_text=0 + n_ref=0 → 12 bytes.
        let mut h = Vec::new();
        h.extend_from_slice(b"BAM\x01");
        h.extend_from_slice(&0u32.to_le_bytes()); // l_text
        h.extend_from_slice(&0u32.to_le_bytes()); // n_ref
        assert_eq!(bam_header_len(&h).unwrap(), Some(12));

        // One reference: magic + l_text=0 + n_ref=1 + l_name=3 + name "r1\0" + l_ref=100
        // → 4 + 4 + 4 + 4 + 3 + 4 = 23 bytes, i.e. 12 + 4 + 3 + 4.
        let mut h2 = Vec::new();
        h2.extend_from_slice(b"BAM\x01");
        h2.extend_from_slice(&0u32.to_le_bytes()); // l_text
        h2.extend_from_slice(&1u32.to_le_bytes()); // n_ref
        h2.extend_from_slice(&3u32.to_le_bytes()); // l_name
        h2.extend_from_slice(b"r1\0"); // name
        h2.extend_from_slice(&100u32.to_le_bytes()); // l_ref
        assert_eq!(bam_header_len(&h2).unwrap(), Some(12 + 4 + 3 + 4));

        // Non-BAM magic → InvalidData error (fail closed; genuinely headerless
        // streams must use `BoundaryState::new_no_header`).
        let non_bam = b"NOT\x01\x00\x00\x00\x00\x00\x00\x00\x00";
        let err = bam_header_len(non_bam).expect_err("non-BAM magic must error");
        assert_eq!(err.kind(), io::ErrorKind::InvalidData);

        // Too short to hold a complete header → Ok(None).
        assert_eq!(bam_header_len(b"BAM").unwrap(), None);
        assert_eq!(bam_header_len(&[]).unwrap(), None);

        // Exactly 8 bytes with BAM magic but no room for n_ref → None.
        let mut short = Vec::new();
        short.extend_from_slice(b"BAM\x01");
        short.extend_from_slice(&0u32.to_le_bytes()); // l_text=0, offset would be 8
        // data.len() == 8, need offset+4 == 12 → Ok(None)
        assert_eq!(bam_header_len(&short).unwrap(), None);
    }
}
