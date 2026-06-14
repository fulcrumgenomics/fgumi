//! `BoundaryState` / `BoundaryBatch`: the BAM record-boundary state machine.
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

    /// Parse BAM header and return the number of bytes consumed.
    /// Returns None if more data is needed.
    fn parse_header_size(data: &[u8]) -> Option<usize> {
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
            return None;
        }

        // Check magic
        if &data[0..4] != fgumi_raw_bam::BAM_MAGIC {
            // Not a valid BAM file, but let's not error here
            // Just return 0 so records start immediately
            return Some(0);
        }

        let l_text = u32::from_le_bytes([data[4], data[5], data[6], data[7]]) as usize;
        let mut offset = 8 + l_text;

        if data.len() < offset + 4 {
            return None;
        }

        let n_ref = u32::from_le_bytes([
            data[offset],
            data[offset + 1],
            data[offset + 2],
            data[offset + 3],
        ]) as usize;
        offset += 4;

        // Parse each reference
        for _ in 0..n_ref {
            if data.len() < offset + 4 {
                return None;
            }
            let l_name = u32::from_le_bytes([
                data[offset],
                data[offset + 1],
                data[offset + 2],
                data[offset + 3],
            ]) as usize;
            offset += 4 + l_name + 4; // l_name + name + l_ref

            if data.len() < offset {
                return None;
            }
        }

        Some(offset)
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
    pub fn find_boundaries(&mut self, decompressed: &[u8]) -> io::Result<BoundaryBatch> {
        // Step 1: Combine leftover with new data into reusable work_buffer
        // This avoids allocating a new Vec on every call
        self.work_buffer.clear();
        if !self.leftover.is_empty() {
            self.work_buffer.append(&mut self.leftover);
        }
        self.work_buffer.extend_from_slice(decompressed);

        // Step 2: Skip header if not already done
        let mut cursor = 0usize;
        if !self.header_skipped {
            if let Some(header_size) = Self::parse_header_size(&self.work_buffer) {
                cursor = header_size;
                self.header_skipped = true;
            } else {
                // Not enough data to parse header, save as leftover and return empty batch
                std::mem::swap(&mut self.leftover, &mut self.work_buffer);
                return Ok(BoundaryBatch { buffer: Vec::new(), offsets: vec![0] });
            }
        }

        // Step 3: Scan for record boundaries (FAST - just read integers)
        let start_cursor = cursor;
        // Pre-size from the previous block's record count so the per-record
        // pushes below don't trigger repeated Vec regrowth (adjacent blocks
        // hold near-identical record counts). First offset is 0 (relative to
        // the start of records).
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
                break; // Incomplete record - becomes leftover
            }

            cursor = record_end;
            // Offset is relative to start of records (after header)
            offsets.push(cursor - start_cursor);
        }

        // Remember this block's offset count to pre-size the next call.
        self.prev_offsets_len = offsets.len();

        // Step 4: Save leftover for next block (reuse allocation)
        // Split work_buffer: [0..start_cursor | start_cursor..cursor | cursor..]
        //                     header (discard) | records (output)    | leftover
        self.leftover.clear();
        self.leftover.extend_from_slice(&self.work_buffer[cursor..]);

        // Extract the records buffer - this allocation is unavoidable as we return ownership
        let buffer = self.work_buffer[start_cursor..cursor].to_vec();

        // Validate: verify each record's block_size matches the offset difference
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
