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

/// Upper bound on a single BAM record's `block_size`. The `block_size` prefix is
/// read straight from the input bytes, so a corrupt or hostile value (e.g.
/// `0xFFFF_FFFF`) would otherwise make the scanner buffer the whole unconsumed
/// tail into `leftover` — up to ~4 GiB from a 4-byte corruption — and only fail
/// at [`BoundaryState::finish`]. A record larger than this is corruption, not
/// data; reject it at the point the prefix is read.
const MAX_RECORD_BLOCK_SIZE: usize = 64 * 1024 * 1024;

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
            if block_size > MAX_RECORD_BLOCK_SIZE {
                return Err(io::Error::new(
                    io::ErrorKind::InvalidData,
                    format!(
                        "BAM record block_size {block_size} exceeds the maximum \
                         {MAX_RECORD_BLOCK_SIZE} (corrupt or truncated input)"
                    ),
                ));
            }
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

            if block_size > MAX_RECORD_BLOCK_SIZE {
                return Err(io::Error::new(
                    io::ErrorKind::InvalidData,
                    format!(
                        "BAM record block_size {block_size} exceeds the maximum \
                         {MAX_RECORD_BLOCK_SIZE} (corrupt or truncated input)"
                    ),
                ));
            }
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
    use rstest::rstest;

    /// Build a BAM header: magic, `l_text` + text, `n_ref` + one entry per `(name, l_ref)`.
    ///
    /// Every test header is constructed here rather than committed as a fixture, so a
    /// header's declared lengths and its actual bytes cannot drift apart.
    fn header(text: &str, refs: &[(&str, u32)]) -> Vec<u8> {
        let mut h = Vec::new();
        h.extend_from_slice(fgumi_raw_bam::BAM_MAGIC);
        h.extend_from_slice(&u32::try_from(text.len()).unwrap().to_le_bytes());
        h.extend_from_slice(text.as_bytes());
        h.extend_from_slice(&u32::try_from(refs.len()).unwrap().to_le_bytes());
        for (name, l_ref) in refs {
            // l_name counts the trailing NUL, matching the BAM spec.
            let name_bytes = format!("{name}\0");
            h.extend_from_slice(&u32::try_from(name_bytes.len()).unwrap().to_le_bytes());
            h.extend_from_slice(name_bytes.as_bytes());
            h.extend_from_slice(&l_ref.to_le_bytes());
        }
        h
    }

    /// Build one BAM record: a 4-byte little-endian `block_size` followed by
    /// `payload_len` bytes of `fill`. The scanner only reads the length prefix, so
    /// the payload just has to be the declared size and be identifiable.
    fn record(payload_len: usize, fill: u8) -> Vec<u8> {
        let mut r = u32::try_from(payload_len).unwrap().to_le_bytes().to_vec();
        r.extend(std::iter::repeat_n(fill, payload_len));
        r
    }

    /// The expected total on-disk size of a record with `payload_len` bytes.
    fn record_len(payload_len: usize) -> usize {
        payload_len + 4
    }

    // ---------------------------------------------------------------------
    // bam_header_len
    // ---------------------------------------------------------------------

    /// What a `bam_header_len` case expects: a parsed length (or `None` for
    /// "need more bytes"), or a hard `InvalidData` rejection.
    #[derive(Debug, Clone, Copy)]
    enum Expect {
        Header(Option<usize>),
        InvalidMagic,
    }

    #[rstest]
    #[case::empty(vec![], Expect::Header(None))]
    #[case::shorter_than_magic(b"BAM".to_vec(), Expect::Header(None))]
    #[case::magic_only_no_room_for_n_ref(
        // 8 bytes: magic + l_text=0. Needs offset+4 == 12 to read n_ref.
        [&fgumi_raw_bam::BAM_MAGIC[..], &0u32.to_le_bytes()[..]].concat(),
        Expect::Header(None)
    )]
    #[case::no_refs(header("", &[]), Expect::Header(Some(12)))]
    #[case::with_header_text(header("@HD\tVN:1.6\n", &[]), Expect::Header(Some(12 + 11)))]
    // 12 + (4 l_name + 3 name + 4 l_ref) == 23
    #[case::one_ref(header("", &[("r1", 100)]), Expect::Header(Some(23)))]
    #[case::two_refs(header("", &[("r1", 100), ("chr2", 200)]), Expect::Header(Some(23 + 4 + 5 + 4)))]
    #[case::truncated_mid_ref_name(header("", &[("r1", 100)])[..20].to_vec(), Expect::Header(None))]
    #[case::truncated_before_l_name(header("", &[("r1", 100)])[..14].to_vec(), Expect::Header(None))]
    #[case::bad_magic(b"NOT\x01\x00\x00\x00\x00\x00\x00\x00\x00".to_vec(), Expect::InvalidMagic)]
    fn bam_header_len_cases(#[case] data: Vec<u8>, #[case] expected: Expect) {
        match expected {
            Expect::Header(len) => assert_eq!(bam_header_len(&data).unwrap(), len),
            Expect::InvalidMagic => {
                let err = bam_header_len(&data).expect_err("expected an error");
                assert_eq!(err.kind(), io::ErrorKind::InvalidData);
            }
        }
    }

    // ---------------------------------------------------------------------
    // Constructors
    // ---------------------------------------------------------------------

    #[test]
    fn new_skips_the_header_and_new_no_header_does_not() {
        let hdr = header("", &[]);
        let rec = record(8, 0xAB);

        // `new` consumes the header, so the first record starts after it.
        let mut with_header = BoundaryState::new();
        let (offsets, range) = with_header.scan(&[hdr.clone(), rec.clone()].concat()).unwrap();
        assert_eq!(range, hdr.len()..hdr.len() + record_len(8));
        assert_eq!(offsets, vec![0, record_len(8)]);
        assert_eq!(with_header.records_bytes(range), rec.as_slice());

        // `new_no_header` treats byte 0 as the first record.
        let mut headerless = BoundaryState::new_no_header();
        let (offsets, range) = headerless.scan(&rec).unwrap();
        assert_eq!(range, 0..record_len(8));
        assert_eq!(offsets, vec![0, record_len(8)]);
    }

    #[test]
    fn default_matches_new_and_still_skips_the_header() {
        let data = [header("", &[]), record(4, 0x11)].concat();
        let mut from_default = BoundaryState::default();
        let mut from_new = BoundaryState::new();
        assert_eq!(from_default.scan(&data).unwrap(), from_new.scan(&data).unwrap());
    }

    // ---------------------------------------------------------------------
    // scan
    // ---------------------------------------------------------------------

    #[rstest]
    #[case::single(vec![8])]
    #[case::several_same_size(vec![4, 4, 4])]
    #[case::mixed_sizes(vec![1, 16, 3, 32])]
    #[case::zero_length_payload(vec![0, 0])]
    fn scan_finds_every_complete_record(#[case] payloads: Vec<usize>) {
        let mut data = header("", &[]);
        for (i, len) in payloads.iter().enumerate() {
            data.extend(record(*len, u8::try_from(i).unwrap()));
        }

        let mut state = BoundaryState::new();
        let (offsets, range) = state.scan(&data).unwrap();

        // One offset per record plus the terminating end offset.
        assert_eq!(offsets.len(), payloads.len() + 1);
        let mut running = 0usize;
        for (i, len) in payloads.iter().enumerate() {
            assert_eq!(offsets[i], running, "record {i} start");
            running += record_len(*len);
        }
        assert_eq!(*offsets.last().unwrap(), running);
        assert_eq!(range.len(), running);
        assert!(state.leftover.is_empty(), "no leftover when every record is complete");
    }

    #[test]
    fn scan_holds_back_a_trailing_partial_record_as_leftover() {
        let complete = record(8, 0x01);
        let partial = &record(64, 0x02)[..10]; // declares 64 bytes, supplies 6
        let data = [header("", &[]), complete.clone(), partial.to_vec()].concat();

        let mut state = BoundaryState::new();
        let (offsets, range) = state.scan(&data).unwrap();

        // Only the complete record is emitted.
        assert_eq!(offsets, vec![0, record_len(8)]);
        assert_eq!(state.records_bytes(range), complete.as_slice());
        assert_eq!(state.leftover, partial, "the partial record is carried forward verbatim");
    }

    #[test]
    fn scan_reassembles_a_record_split_across_two_blocks() {
        let rec = record(32, 0x7E);
        let data = [header("", &[]), rec.clone()].concat();
        let split = data.len() - 20; // cut mid-record

        let mut state = BoundaryState::new();

        // First block ends mid-record: nothing complete yet.
        let (offsets, range) = state.scan(&data[..split]).unwrap();
        assert_eq!(offsets, vec![0], "no complete record in the first block");
        assert_eq!(range.len(), 0);
        assert!(!state.leftover.is_empty());

        // Second block completes it, and the bytes match the original record.
        let (offsets, range) = state.scan(&data[split..]).unwrap();
        assert_eq!(offsets, vec![0, record_len(32)]);
        assert_eq!(state.records_bytes(range), rec.as_slice());
        assert!(state.leftover.is_empty());
    }

    #[test]
    fn scan_defers_when_the_header_itself_is_incomplete() {
        let hdr = header("@HD\tVN:1.6\n", &[("r1", 100)]);
        let mut state = BoundaryState::new();

        // A prefix too short to hold the whole header yields an empty batch and
        // stashes everything for the next call.
        let (offsets, range) = state.scan(&hdr[..10]).unwrap();
        assert_eq!(offsets, vec![0]);
        assert_eq!(range, 0..0);
        assert_eq!(state.leftover, hdr[..10]);
        assert!(!state.header_skipped, "header must not be marked skipped yet");

        // The rest of the header plus a record then parses normally.
        let rec = record(8, 0x5A);
        let (offsets, range) = state.scan(&[&hdr[10..], rec.as_slice()].concat()).unwrap();
        assert!(state.header_skipped);
        assert_eq!(offsets, vec![0, record_len(8)]);
        assert_eq!(state.records_bytes(range), rec.as_slice());
    }

    #[test]
    fn scan_propagates_a_bad_magic_as_invalid_data() {
        let mut state = BoundaryState::new();
        let err = state.scan(b"NOPE\x00\x00\x00\x00\x00\x00\x00\x00").expect_err("bad magic");
        assert_eq!(err.kind(), io::ErrorKind::InvalidData);
    }

    // ---------------------------------------------------------------------
    // find_boundaries (owned-buffer wrapper over scan)
    // ---------------------------------------------------------------------

    #[test]
    fn find_boundaries_returns_the_same_bytes_scan_would_borrow() {
        let recs = [record(8, 0x01), record(16, 0x02), record(2, 0x03)].concat();
        let data = [header("", &[("r1", 10)]), recs.clone()].concat();

        let mut owned = BoundaryState::new();
        let batch = owned.find_boundaries(&data).unwrap();

        let mut borrowed = BoundaryState::new();
        let (offsets, range) = borrowed.scan(&data).unwrap();

        assert_eq!(batch.offsets, offsets);
        assert_eq!(batch.buffer, borrowed.records_bytes(range));
        assert_eq!(batch.buffer, recs, "the owned copy is the record bytes, header excluded");

        // Offsets slice the buffer into the original records.
        for i in 0..batch.offsets.len() - 1 {
            let rec = &batch.buffer[batch.offsets[i]..batch.offsets[i + 1]];
            let declared = u32::from_le_bytes(rec[..4].try_into().unwrap()) as usize;
            assert_eq!(declared, rec.len() - 4, "record {i} length prefix matches its slice");
        }
    }

    #[test]
    fn find_boundaries_on_a_header_only_block_yields_an_empty_batch() {
        let mut state = BoundaryState::new();
        let batch = state.find_boundaries(&header("", &[])).unwrap();
        assert!(batch.buffer.is_empty());
        assert_eq!(batch.offsets, vec![0]);
    }

    // ---------------------------------------------------------------------
    // finish
    // ---------------------------------------------------------------------

    #[test]
    fn finish_returns_none_when_nothing_is_pending() {
        let mut state = BoundaryState::new_no_header();
        assert!(state.finish().unwrap().is_none());

        // Also none after a scan that consumed every record.
        let (_, _) = state.scan(&record(4, 0x09)).unwrap();
        assert!(state.finish().unwrap().is_none());
    }

    #[test]
    fn finish_emits_leftover_that_forms_complete_records() {
        // `scan` never leaves a *complete* record behind, so the pending buffer is
        // seeded directly to exercise the success branch.
        let mut state = BoundaryState::new_no_header();
        state.leftover = [record(4, 0xA1), record(8, 0xA2)].concat();

        let batch = state.finish().unwrap().expect("complete records must be emitted");
        assert_eq!(batch.offsets, vec![0, record_len(4), record_len(4) + record_len(8)]);
        assert_eq!(batch.buffer.len(), record_len(4) + record_len(8));
        assert!(state.leftover.is_empty(), "finish takes the pending bytes");
    }

    #[rstest]
    // Declares a 64-byte payload but only supplies part of it.
    #[case::truncated_payload(record(64, 0x02)[..10].to_vec())]
    // Fewer than 4 bytes cannot even hold a block-size prefix.
    #[case::one_trailing_byte(vec![0x00])]
    #[case::three_trailing_bytes(vec![0x00, 0x01, 0x02])]
    // A whole record followed by an unusable tail.
    #[case::complete_then_stray_bytes([record(4, 0x03), vec![0xFF, 0xFF]].concat())]
    fn finish_rejects_incomplete_trailing_bytes(#[case] pending: Vec<u8>) {
        let mut state = BoundaryState::new_no_header();
        state.leftover = pending;
        let err = state.finish().expect_err("truncated input must fail closed");
        assert_eq!(err.kind(), io::ErrorKind::UnexpectedEof);
    }
}
