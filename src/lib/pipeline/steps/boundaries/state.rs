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

/// Upper bound on bytes the scanner will carry forward waiting for one record
/// (or the BAM header) to complete.
///
/// This is a **corruption backstop, not a biological limit**. The scanner is a
/// push scanner: it is handed decompressed blocks and cannot ask for "the rest
/// of this record" the way a pull reader (htslib's `bam_read1`) can, so a
/// corrupt length field makes every block look like "not enough data yet" and
/// the carry grows for the rest of the stream. That still terminates —
/// [`BoundaryState::finish`] reports the shortfall at EOF — but only after
/// buffering the remaining input, so a large corrupt file exhausts memory
/// before it can print the diagnostic it was about to print.
///
/// The value sits in the gap between what BAM permits and what BAM files
/// contain. `block_size` is a `u32`, so the format allows a ~4 GiB record, while
/// the largest records seen in practice — ONT ultra-long reads carrying
/// methylation tags — run a few MB. 256 MiB is roughly 25-50x above that and
/// still bounds memory hard.
///
/// It is deliberately NOT tied to the pipeline's queue budget: byte-bounded
/// queues always admit at least one item regardless of size, so a legitimate
/// multi-MB record flows through a 4 MiB per-step budget, and bounding the carry
/// there would reject valid BAMs.
const MAX_CARRY_BYTES: usize = 256 * 1024 * 1024;

/// Ceiling on the bytes carried for a single *open queryname run* under
/// [`BatchCut::Queryname`], enforced independently of `min_emit_bytes` and of
/// [`MAX_CARRY_BYTES`].
///
/// The two existing bounds do not cover this source of growth. `min_emit_bytes`
/// gates only the *closed* prefix, so it never bounds a run that has not closed.
/// `MAX_CARRY_BYTES` bounds one *incomplete* trailing record; a queryname run is
/// made of individually *complete* records, so it never trips that backstop.
/// Without a dedicated cap a stream that never closes a run — every record
/// sharing one read name, or a run of empty BAM names (`l_read_name == 1`) —
/// grows `qn_buffer` until allocation fails.
///
/// A queryname run is the records sharing one name: a template's primary pair
/// plus its secondary/supplementary alignments. Even a pathologically deep
/// template is a few MB; 256 MiB is orders of magnitude above any real run and
/// still bounds memory hard. It is a *separate* limit from `MAX_CARRY_BYTES`
/// (whose value is justified for one incomplete record, not for a complete run)
/// even though the two happen to share a value today.
const MAX_QUERYNAME_RUN_BYTES: usize = 256 * 1024 * 1024;

/// Where the boundary scanner is allowed to cut a batch.
///
/// The scanner always cuts on *record* boundaries (never mid-record). This
/// selects whether it may also cut mid-run of a queryname, or must keep every
/// queryname run whole within one emitted batch.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Default)]
pub enum BatchCut {
    /// Cut at any record boundary (today's behaviour). A queryname run may span
    /// two emitted batches. Every non-grouping chain uses this.
    #[default]
    Record,
    /// Cut only at queryname boundaries: the trailing partial queryname run is
    /// carried to the next batch, so every emitted batch is *closed under
    /// queryname* (no run straddles a batch). Enables the parallel
    /// `AssembleTemplates` grouper. The trailing carry is at most one whole
    /// queryname run rather than one partial record.
    Queryname,
}

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
    /// Where this scanner may cut a batch (record vs queryname boundary).
    cut: BatchCut,
    /// Coalescing target for `BatchCut::Queryname`: keep absorbing input into
    /// the carry until the *emittable* (closed) prefix reaches this many bytes
    /// before emitting, so downstream batches are ~1-4 MiB rather than one
    /// ~64 KiB BGZF block. `0` disables coalescing (emit as soon as anything is
    /// closed). Ignored for `BatchCut::Record`.
    min_emit_bytes: usize,
    /// Read name of the last record in the batch this scanner most recently
    /// emitted, under `BatchCut::Queryname`. The next emitted batch's first
    /// record must have a *different* name, or the closed-under-queryname
    /// invariant is broken; the scanner asserts this once per emit as an
    /// always-on self-check (a cut-logic regression fails on the first
    /// straddling run rather than silently corrupting output). Empty when no
    /// batch has been emitted yet.
    last_emitted_name: Vec<u8>,
    /// Coalescing accumulator for [`BatchCut::Queryname`]: the bytes of all
    /// complete, already-scanned records not yet emitted. New records are
    /// appended here as they are scanned, so they are never re-scanned across
    /// calls (the fix for the O(n^2) carry re-walk). Emitting splits this at a
    /// queryname boundary and retains the trailing run. Unused in
    /// [`BatchCut::Record`] mode.
    qn_buffer: Vec<u8>,
    /// Record-end offsets into [`Self::qn_buffer`], `num_records + 1` entries
    /// with a leading `0`. Grows as records are scanned; split on emit. Unused
    /// in `Record` mode.
    qn_offsets: Vec<usize>,
    /// Ceiling on the bytes of a single open queryname run carried in
    /// `qn_buffer` (see [`MAX_QUERYNAME_RUN_BYTES`]). A field rather than a bare
    /// constant so tests can drive the rejection path with a small cap instead
    /// of buffering 256 MiB; production callers get [`MAX_QUERYNAME_RUN_BYTES`].
    max_queryname_run_bytes: usize,
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
            cut: BatchCut::Record,
            min_emit_bytes: 0,
            last_emitted_name: Vec::new(),
            qn_buffer: Vec::new(),
            qn_offsets: vec![0],
            max_queryname_run_bytes: MAX_QUERYNAME_RUN_BYTES,
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
            cut: BatchCut::Record,
            min_emit_bytes: 0,
            last_emitted_name: Vec::new(),
            qn_buffer: Vec::new(),
            qn_offsets: vec![0],
            max_queryname_run_bytes: MAX_QUERYNAME_RUN_BYTES,
        }
    }

    /// Select the batch cut mode (default [`BatchCut::Record`]). For
    /// [`BatchCut::Queryname`], `min_emit_bytes` sets the coalescing target
    /// (0 = emit as soon as any run is closed). Chainable on either constructor.
    #[must_use]
    pub fn with_cut(mut self, cut: BatchCut, min_emit_bytes: usize) -> Self {
        self.cut = cut;
        self.min_emit_bytes = min_emit_bytes;
        self
    }

    /// Override the open-queryname-run carry ceiling. Test-facing: it makes the
    /// [`MAX_QUERYNAME_RUN_BYTES`] rejection reachable without buffering 256 MiB.
    #[cfg(test)]
    #[must_use]
    fn with_max_queryname_run_bytes(mut self, max_queryname_run_bytes: usize) -> Self {
        self.max_queryname_run_bytes = max_queryname_run_bytes;
        self
    }

    /// Parse BAM header and return the number of bytes consumed.
    ///
    /// `Ok(None)` means more data is needed. `Err` means the header's own length
    /// fields are corrupt — their sum overflows `usize`, so no amount of further
    /// data could satisfy them.
    ///
    /// # Errors
    ///
    /// Returns `InvalidData` when `l_text` or a reference's `l_name` overflows
    /// the running offset. Unreachable on 64-bit targets (both are `u32` widened
    /// to `usize`); present because a length field taken from the stream must be
    /// folded in with checked arithmetic rather than trusted to fit.
    fn parse_header_size(data: &[u8]) -> io::Result<Option<usize>> {
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

        // Check magic
        if &data[0..4] != fgumi_raw_bam::BAM_MAGIC {
            // Not a valid BAM file, but let's not error here
            // Just return 0 so records start immediately
            return Ok(Some(0));
        }

        let l_text = u32::from_le_bytes([data[4], data[5], data[6], data[7]]) as usize;
        let Some(mut offset) = 8usize.checked_add(l_text) else {
            return Err(io::Error::new(
                io::ErrorKind::InvalidData,
                format!(
                    "FindBamBoundaries: BAM header l_text={l_text} overflows the header offset"
                ),
            ));
        };

        if data.len() < offset + 4 {
            return Ok(None);
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
                return Ok(None);
            }
            let l_name = u32::from_le_bytes([
                data[offset],
                data[offset + 1],
                data[offset + 2],
                data[offset + 3],
            ]) as usize;
            // l_name + name + l_ref, folded in with checked arithmetic for the
            // same reason as `l_text` above.
            let Some(next) = offset.checked_add(8).and_then(|o| o.checked_add(l_name)) else {
                return Err(io::Error::new(
                    io::ErrorKind::InvalidData,
                    format!(
                        "FindBamBoundaries: BAM header l_name={l_name} overflows \
                         the header offset at {offset}"
                    ),
                ));
            };
            offset = next;

            if data.len() < offset {
                return Ok(None);
            }
        }

        Ok(Some(offset))
    }

    /// Read the queryname of record `i` (0-based) framed by `offsets` relative
    /// to `records_start` within `buf`. Record `i`'s body is
    /// `buf[records_start + offsets[i] + 4 .. records_start + offsets[i+1]]`
    /// (the 4-byte `block_size` prefix skipped). Fails closed (`InvalidData`)
    /// rather than panicking on a body too short to hold its own read name.
    fn record_name_in<'a>(
        buf: &'a [u8],
        records_start: usize,
        offsets: &[usize],
        i: usize,
    ) -> io::Result<&'a [u8]> {
        let body = &buf[records_start + offsets[i] + 4..records_start + offsets[i + 1]];
        if body.len() < fgumi_raw_bam::MIN_BAM_RECORD_LEN {
            return Err(io::Error::new(
                io::ErrorKind::InvalidData,
                format!(
                    "FindBamBoundaries: BAM record too short for a queryname cut ({} < {})",
                    body.len(),
                    fgumi_raw_bam::MIN_BAM_RECORD_LEN,
                ),
            ));
        }
        let l_read_name = body[8] as usize;
        if l_read_name == 0 || body.len() < 32 + l_read_name {
            return Err(io::Error::new(
                io::ErrorKind::InvalidData,
                format!(
                    "FindBamBoundaries: BAM record truncated: l_read_name={l_read_name} \
                     but only {} bytes after the fixed header",
                    body.len().saturating_sub(32),
                ),
            ));
        }
        Ok(fgumi_raw_bam::read_name(body))
    }

    /// Read the queryname of record `i` in the coalescing accumulator
    /// `qn_buffer`, framed by `offsets` (absolute into `qn_buffer`). Thin
    /// wrapper over [`Self::record_name_in`] with `records_start = 0`.
    ///
    /// # Errors
    ///
    /// `InvalidData` if the record body is too short to hold its read name.
    fn qn_record_name<'a>(buf: &'a [u8], offsets: &[usize], i: usize) -> io::Result<&'a [u8]> {
        Self::record_name_in(buf, 0, offsets, i)
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
    /// Returns `InvalidData` if the bytes carried forward waiting for the header
    /// or for one record to complete exceed this module's `MAX_CARRY_BYTES`
    /// corruption backstop. A well-formed stream always completes a record
    /// within a bounded carry, so exceeding it means a length field is corrupt.
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
            if let Some(header_size) = Self::parse_header_size(&self.work_buffer)? {
                cursor = header_size;
                self.header_skipped = true;
            } else {
                // Not enough data to parse header yet. Bound the carry: a corrupt
                // `l_text` / `l_name` makes this branch unreachable-to-satisfy, so
                // without the check every remaining block accumulates here.
                if self.work_buffer.len() > MAX_CARRY_BYTES {
                    return Err(io::Error::new(
                        io::ErrorKind::InvalidData,
                        format!(
                            "FindBamBoundaries: BAM header still unresolved after {} byte(s) \
                             (limit {MAX_CARRY_BYTES}); the header length fields are corrupt",
                            self.work_buffer.len(),
                        ),
                    ));
                }
                // Save as leftover and return an empty batch.
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

            // Checked: `block_size` comes from the stream. Unreachable on
            // 64-bit, but the framing rule is checked-before-use.
            let Some(record_end) = cursor.checked_add(4).and_then(|c| c.checked_add(block_size))
            else {
                return Err(io::Error::new(
                    io::ErrorKind::InvalidData,
                    format!(
                        "FindBamBoundaries: record end overflows \
                         (cursor={cursor}, block_size={block_size})"
                    ),
                ));
            };
            if record_end > self.work_buffer.len() {
                break; // Incomplete record - becomes leftover
            }

            cursor = record_end;
            // Offset is relative to start of records (after header)
            offsets.push(cursor - start_cursor);
        }

        // Step 4a: Split the incomplete trailing record into `leftover` for the
        // next call. `offsets.last()` is the byte length of the complete-records
        // region (relative to `start_cursor`); everything after `cursor` is the
        // partial tail. Bound the carry before retaining it: a corrupt
        // `block_size` keeps the scan breaking at the same record forever, so
        // without this the leftover grows to the size of the remaining input.
        // `finish` would still catch the shortfall at EOF, but only after the
        // memory was already spent. (This bound is on the genuinely-incomplete
        // tail only; the queryname coalescing carry lives in `qn_buffer`, and its
        // trailing open run is bounded by `MAX_QUERYNAME_RUN_BYTES` — not this
        // backstop, and not `min_emit_bytes`, which is only the coalescing target.)
        let carry_len = self.work_buffer.len() - cursor;
        if carry_len > MAX_CARRY_BYTES {
            return Err(io::Error::new(
                io::ErrorKind::InvalidData,
                format!(
                    "FindBamBoundaries: incomplete BAM record has carried {carry_len} byte(s) \
                     (limit {MAX_CARRY_BYTES}); its block_size prefix is corrupt",
                ),
            ));
        }
        self.leftover.clear();
        self.leftover.extend_from_slice(&self.work_buffer[cursor..]);

        // The complete-records region is `work_buffer[start_cursor..cursor]`,
        // described by `offsets` (relative to `start_cursor`, leading `0`).
        let records = &self.work_buffer[start_cursor..cursor];

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
                    records[start],
                    records[start + 1],
                    records[start + 2],
                    records[start + 3],
                ]) as usize;
                let expected = end - start - 4;
                debug_assert_eq!(
                    stored, expected,
                    "find_boundaries: block_size mismatch at record {i}: stored={stored}, expected={expected}"
                );
            }
        }

        // Step 4b: emit.
        if self.cut != BatchCut::Queryname {
            // Record mode (default, unchanged): return this call's complete
            // records directly. No cross-call accumulation.
            self.prev_offsets_len = offsets.len();
            let buffer = records.to_vec();
            return Ok(BoundaryBatch { buffer, offsets });
        }

        // Queryname mode: append this call's complete records to the coalescing
        // accumulator (they are never re-scanned), then decide the emit split on
        // the accumulator without re-scanning already-scanned records. This is
        // the O(1)-amortized-per-record path that replaces the old whole-carry
        // re-walk.
        let base = self.qn_buffer.len();
        self.qn_buffer.extend_from_slice(records);
        // `offsets[0]` is the leading 0 already present as the last entry of
        // `qn_offsets` (a record-end sentinel); append the rest, shifted into
        // `qn_buffer` coordinates.
        for w in offsets.windows(2) {
            self.qn_offsets.push(base + w[1]);
        }
        self.prev_offsets_len = offsets.len();

        self.emit_queryname_batch()
    }

    /// Decide and perform the queryname emit split on the coalescing
    /// accumulator (`qn_buffer` / `qn_offsets`). Walks back from the last record
    /// to the start of the trailing queryname run (bounded by the run length,
    /// not the accumulator size), applies the `min_emit_bytes` coalescing gate,
    /// and — when a non-empty closed prefix is ready — returns it as a
    /// `BoundaryBatch`, retaining the trailing run (and its offsets) for the
    /// next emit. Returns an empty batch when nothing is closed yet (keep
    /// accumulating). `force` (set by `finish` at EOF) emits the entire
    /// accumulator regardless of the coalescing gate.
    ///
    /// # Errors
    ///
    /// `InvalidData` if a record is too short to hold its read name, or if the
    /// closed-under-queryname self-check trips (a run straddled a batch).
    fn emit_queryname_batch(&mut self) -> io::Result<BoundaryBatch> {
        self.emit_queryname_batch_inner(false)
    }

    fn emit_queryname_batch_inner(&mut self, force: bool) -> io::Result<BoundaryBatch> {
        let n_records = self.qn_offsets.len() - 1;
        if n_records == 0 {
            return Ok(BoundaryBatch { buffer: Vec::new(), offsets: vec![0] });
        }

        // `j` = index of the first record of the trailing queryname run (the run
        // to carry). 0 = the whole accumulator is one run (or coalescing says
        // "not enough closed yet"): emit nothing unless forced.
        let mut j = n_records - 1;
        let mut cur_name = Self::qn_record_name(&self.qn_buffer, &self.qn_offsets, j)?;
        while j > 0 {
            let prev_name = Self::qn_record_name(&self.qn_buffer, &self.qn_offsets, j - 1)?;
            if prev_name != cur_name {
                break;
            }
            j -= 1;
            cur_name = prev_name;
        }
        // First record of the trailing (open) queryname run, captured *before*
        // the coalescing gate below rewrites `j` to 0. The cap must measure this
        // run alone — measuring from the post-coalescing `j` would fold the
        // closed, coalesced prefix into the count and could reject valid input.
        let run_start = j;

        if force {
            // EOF flush: emit everything (the trailing run is complete too).
            j = n_records;
        } else if self.min_emit_bytes > 0 && self.qn_offsets[j] < self.min_emit_bytes {
            // Closed prefix below the coalescing target: keep accumulating.
            j = 0;
        }

        // Bound the trailing *open* queryname run (`qn_buffer[qn_offsets[run_start]..]`)
        // — not the coalesced closed prefix before it, which `min_emit_bytes`
        // already bounds. Neither `min_emit_bytes` nor the incomplete-tail
        // `MAX_CARRY_BYTES` backstop bounds the open run, so a stream that never
        // closes a run — every record one name, or empty read names — would grow
        // `qn_buffer` without limit. Reject it as corrupt / not-queryname-grouped.
        // Skipped at EOF: `force` emits the complete final run rather than
        // carrying it, so there is nothing to bound.
        if !force {
            let run_len = self.qn_buffer.len() - self.qn_offsets[run_start];
            if run_len > self.max_queryname_run_bytes {
                let limit = self.max_queryname_run_bytes;
                return Err(io::Error::new(
                    io::ErrorKind::InvalidData,
                    format!(
                        "FindBamBoundaries: an open queryname run has carried {run_len} \
                         byte(s) (limit {limit}); the stream is not queryname-grouped or \
                         its read names are corrupt",
                    ),
                ));
            }
        }

        if j == 0 {
            // Nothing to emit yet.
            return Ok(BoundaryBatch { buffer: Vec::new(), offsets: vec![0] });
        }

        // Emit records 0..j; retain j.. (the trailing run) for the next emit.
        let split = self.qn_offsets[j];
        let emit_offsets: Vec<usize> = self.qn_offsets[..=j].to_vec();

        // Self-check + carry the boundary name. The first emitted record's name
        // must differ from the previous emitted batch's last name, else a
        // queryname run straddled a batch — the invariant `AssembleTemplates`
        // relies on and cannot itself detect.
        let first_name = Self::qn_record_name(&self.qn_buffer, &emit_offsets, 0)?;
        if !self.last_emitted_name.is_empty() && self.last_emitted_name == first_name {
            return Err(io::Error::new(
                io::ErrorKind::InvalidData,
                format!(
                    "FindBamBoundaries: queryname cut invariant violated \
                         (a run of read name {:?} straddled a batch boundary)",
                    String::from_utf8_lossy(first_name),
                ),
            ));
        }
        let last_name = Self::qn_record_name(&self.qn_buffer, &emit_offsets, j - 1)?;
        self.last_emitted_name.clear();
        self.last_emitted_name.extend_from_slice(last_name);

        let buffer = self.qn_buffer[..split].to_vec();

        // Retain the trailing run: shift its bytes to the front of `qn_buffer`
        // and rebase its offsets to start at 0.
        self.qn_buffer.drain(..split);
        let retained: Vec<usize> = self.qn_offsets[j..].iter().map(|&o| o - split).collect();
        self.qn_offsets = retained;

        Ok(BoundaryBatch { buffer, offsets: emit_offsets })
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
        // Queryname mode: complete records were scanned into the coalescing
        // accumulator as they arrived, so only a genuinely-incomplete tail can
        // remain in `leftover` — that is a truncated record at EOF (error).
        // Flush whatever the coalescing gate held back in `qn_buffer` as the
        // final batch.
        if self.cut == BatchCut::Queryname {
            if !self.leftover.is_empty() {
                return Err(io::Error::new(
                    io::ErrorKind::UnexpectedEof,
                    format!(
                        "Incomplete BAM record at EOF: {} trailing byte(s) cannot form a \
                         complete record",
                        self.leftover.len(),
                    ),
                ));
            }
            let batch = self.emit_queryname_batch_inner(true)?;
            if batch.buffer.is_empty() {
                return Ok(None);
            }
            return Ok(Some(batch));
        }

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

            // Checked for the same reason as the scan loop above.
            let Some(record_end) = cursor.checked_add(4).and_then(|c| c.checked_add(block_size))
            else {
                return Err(io::Error::new(
                    io::ErrorKind::InvalidData,
                    format!(
                        "FindBamBoundaries: record end overflows at EOF \
                         (cursor={cursor}, block_size={block_size})"
                    ),
                ));
            };
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

    /// Frame `payload` as one BAM record: `block_size: u32 LE` + body. The
    /// scanner never inspects the body, so an opaque payload is sufficient —
    /// and keeps each case's record lengths readable in the table.
    fn record(payload: &[u8]) -> Vec<u8> {
        let mut framed = Vec::with_capacity(4 + payload.len());
        framed.extend_from_slice(
            &u32::try_from(payload.len()).expect("payload fits u32").to_le_bytes(),
        );
        framed.extend_from_slice(payload);
        framed
    }

    /// A well-formed binary BAM header: magic, `l_text` + text, `n_ref`, then
    /// `(l_name, name\0, l_ref)` per reference. This is exactly the prefix
    /// `parse_header_size` walks, so building it honestly is what makes the
    /// header-skip assertions meaningful.
    fn bam_header(text: &str, refs: &[(&str, u32)]) -> Vec<u8> {
        let mut header = Vec::new();
        header.extend_from_slice(fgumi_raw_bam::BAM_MAGIC);
        header.extend_from_slice(&u32::try_from(text.len()).expect("text fits u32").to_le_bytes());
        header.extend_from_slice(text.as_bytes());
        header.extend_from_slice(&u32::try_from(refs.len()).expect("n_ref fits u32").to_le_bytes());
        for (name, ref_len) in refs {
            let mut name_bytes = name.as_bytes().to_vec();
            name_bytes.push(0); // names are NUL-terminated, and l_name counts the NUL
            header.extend_from_slice(
                &u32::try_from(name_bytes.len()).expect("l_name fits u32").to_le_bytes(),
            );
            header.extend_from_slice(&name_bytes);
            header.extend_from_slice(&ref_len.to_le_bytes());
        }
        header
    }

    /// Concatenate framed records, as one decompressed block would hold them.
    fn records(payloads: &[&[u8]]) -> Vec<u8> {
        payloads.iter().flat_map(|p| record(p)).collect()
    }

    // ========================================================================
    // Header skipping
    // ========================================================================

    /// The header must be consumed and excluded from the emitted buffer, with
    /// only record bytes surviving. Parameterized over header shapes because
    /// `parse_header_size`'s arithmetic differs per section: text length and
    /// per-reference `(l_name, name, l_ref)` walking are separate loops, and a
    /// bug in either would pass the other's case.
    #[rstest]
    #[case::no_text_no_refs("", &[])]
    #[case::text_only("@HD\tVN:1.6\n", &[])]
    #[case::one_ref("@HD\tVN:1.6\n", &[("chr1", 1000)])]
    #[case::several_refs("@HD\tVN:1.6\n", &[("chr1", 1000), ("chr2", 2000), ("chrM", 16569)])]
    #[case::empty_text_with_refs("", &[("chr1", 248_956_422)])]
    fn find_boundaries_skips_the_header_and_emits_only_record_bytes(
        #[case] text: &str,
        #[case] refs: &[(&str, u32)],
    ) {
        let payloads: [&[u8]; 2] = [&[1u8; 8], &[2u8; 16]];
        let mut block = bam_header(text, refs);
        block.extend_from_slice(&records(&payloads));

        let mut state = BoundaryState::new();
        let batch = state.find_boundaries(&block).expect("well-formed block scans");

        assert_eq!(
            batch.buffer,
            records(&payloads),
            "header bytes must not reach the output buffer"
        );
        assert_eq!(batch.offsets, vec![0, 12, 32], "offsets are relative to the first record");
    }

    /// A `BoundaryState::new_no_header()` treats byte 0 as the first record —
    /// this is the runall-spliced sub-pipeline case, where an upstream stage
    /// already consumed the header.
    #[test]
    fn new_no_header_treats_the_first_byte_as_a_record_boundary() {
        let payloads: [&[u8]; 2] = [&[7u8; 4], &[9u8; 4]];
        let block = records(&payloads);

        let mut state = BoundaryState::new_no_header();
        let batch = state.find_boundaries(&block).expect("record-aligned block scans");

        assert_eq!(batch.buffer, block, "every byte is a record byte when the header is skipped");
        assert_eq!(batch.offsets, vec![0, 8, 16]);
    }

    /// `Default` must agree with `new()` — i.e. it expects a header. Asserted
    /// behaviorally (the header is consumed) rather than by comparing private
    /// fields, so the test survives a field-level refactor.
    #[test]
    fn default_expects_a_header_like_new() {
        let payloads: [&[u8]; 1] = [&[1u8; 8]];
        let mut block = bam_header("", &[]);
        block.extend_from_slice(&records(&payloads));

        let batch = BoundaryState::default().find_boundaries(&block).expect("scans");
        assert_eq!(batch.buffer, records(&payloads));
    }

    /// A header split across blocks must be buffered, not mis-parsed: the first
    /// call cannot determine the header size, so it emits an empty batch and
    /// carries every byte forward. Parameterized over cut points that land in
    /// each distinct section of `parse_header_size`'s walk, since each has its
    /// own "need more data" early return.
    #[rstest]
    #[case::mid_magic(2)]
    #[case::after_magic_before_l_text(4)]
    #[case::mid_text(10)]
    #[case::before_n_ref(15)]
    #[case::mid_ref_name(24)]
    fn a_header_split_across_blocks_is_buffered_until_complete(#[case] cut: usize) {
        let payloads: [&[u8]; 1] = [&[3u8; 8]];
        let header = bam_header("@HD\tVN:1.6\n", &[("chr1", 1000)]);
        assert!(
            cut < header.len(),
            "cut must fall inside the header for this case to mean anything"
        );
        let mut full = header.clone();
        full.extend_from_slice(&records(&payloads));

        let mut state = BoundaryState::new();
        let first = state.find_boundaries(&full[..cut]).expect("partial header is not an error");
        assert!(first.buffer.is_empty(), "no records can be emitted before the header is parsed");
        assert_eq!(first.offsets, vec![0]);

        let second = state.find_boundaries(&full[cut..]).expect("rest of the header completes it");
        assert_eq!(second.buffer, records(&payloads), "records resume once the header is consumed");
        assert_eq!(second.offsets, vec![0, 12]);
    }

    /// Data that does not start with `BAM\1` is treated as record bytes from
    /// offset 0 rather than erroring — the documented "not a valid BAM file,
    /// but let's not error here" behavior. Worth pinning because it is
    /// surprising: a caller who feeds the wrong stream gets garbage records,
    /// not a diagnostic.
    #[test]
    fn a_non_bam_magic_prefix_is_scanned_as_records_from_offset_zero() {
        let payloads: [&[u8]; 1] = [&[5u8; 8]];
        let block = records(&payloads);
        assert_ne!(&block[0..4], fgumi_raw_bam::BAM_MAGIC, "fixture must not look like a header");

        let mut state = BoundaryState::new();
        let batch = state.find_boundaries(&block).expect("scans");

        assert_eq!(batch.buffer, block);
        assert_eq!(batch.offsets, vec![0, 12]);
    }

    /// Fewer than 8 bytes cannot even hold magic + `l_text`, so the scanner
    /// must ask for more data instead of reading past the end.
    #[test]
    fn fewer_than_eight_bytes_cannot_resolve_the_header() {
        let mut state = BoundaryState::new();
        let batch = state.find_boundaries(b"BAM\x01").expect("short input is not an error");
        assert!(batch.buffer.is_empty());
        assert_eq!(batch.offsets, vec![0]);
    }

    // ========================================================================
    // Record scanning and cross-block carryover
    // ========================================================================

    /// The core carryover contract: a record straddling two blocks is held back
    /// on the first call and emitted whole on the second, never split. The
    /// concatenated output across both calls must equal the input records
    /// exactly — that is the invariant downstream parsing depends on.
    #[rstest]
    #[case::split_mid_body(18)]
    #[case::split_on_a_record_start(12)]
    #[case::split_mid_length_prefix(14)]
    #[case::split_one_byte_in(1)]
    fn a_record_spanning_two_blocks_is_emitted_whole(#[case] cut: usize) {
        let payloads: [&[u8]; 3] = [&[1u8; 8], &[2u8; 8], &[3u8; 8]];
        let all = records(&payloads);
        assert!(cut < all.len());

        let mut state = BoundaryState::new_no_header();
        let first = state.find_boundaries(&all[..cut]).expect("first half scans");
        let second = state.find_boundaries(&all[cut..]).expect("second half scans");

        let mut seen = first.buffer.clone();
        seen.extend_from_slice(&second.buffer);
        assert_eq!(seen, all, "no byte may be dropped or duplicated across the split");
        assert!(state.finish().expect("nothing left over").is_none());
    }

    /// Offsets must always start at 0, end at `buffer.len()`, and step by each
    /// record's framed size — this is the contract that lets a consumer slice
    /// records with `offsets.windows(2)`.
    #[test]
    fn offsets_bracket_every_record_in_the_emitted_buffer() {
        let payloads: [&[u8]; 3] = [&[1u8; 4], &[2u8; 20], &[3u8; 0]];
        let block = records(&payloads);

        let mut state = BoundaryState::new_no_header();
        let batch = state.find_boundaries(&block).expect("scans");

        assert_eq!(batch.offsets, vec![0, 8, 32, 36]);
        assert_eq!(*batch.offsets.last().expect("non-empty"), batch.buffer.len());
        for (i, window) in batch.offsets.windows(2).enumerate() {
            let body = &batch.buffer[window[0] + 4..window[1]];
            assert_eq!(body, payloads[i], "record {i} must round-trip through its offsets");
        }
    }

    /// An empty block is a legal no-op: nothing to scan, nothing carried.
    #[test]
    fn an_empty_block_yields_an_empty_batch() {
        let mut state = BoundaryState::new_no_header();
        let batch = state.find_boundaries(&[]).expect("empty input is not an error");
        assert!(batch.buffer.is_empty());
        assert_eq!(batch.offsets, vec![0]);
    }

    /// Feeding many blocks in sequence must give the same records as feeding
    /// one concatenated block. This is what exercises the `prev_offsets_len`
    /// pre-sizing path across calls, and it is the property that actually
    /// matters: block framing must not be observable in the output.
    #[test]
    fn streaming_many_blocks_matches_scanning_one_concatenated_block() {
        let payloads: Vec<Vec<u8>> = (0..32u8).map(|i| vec![i; usize::from(i) + 1]).collect();
        let refs: Vec<&[u8]> = payloads.iter().map(Vec::as_slice).collect();
        let all = records(&refs);

        let one_shot = {
            let mut state = BoundaryState::new_no_header();
            let batch = state.find_boundaries(&all).expect("scans");
            assert!(state.finish().expect("nothing left over").is_none());
            batch.buffer
        };

        let streamed = {
            let mut state = BoundaryState::new_no_header();
            let mut out = Vec::new();
            for chunk in all.chunks(7) {
                out.extend_from_slice(&state.find_boundaries(chunk).expect("scans").buffer);
            }
            assert!(state.finish().expect("nothing left over").is_none());
            out
        };

        assert_eq!(streamed, one_shot, "block framing must not be observable downstream");
        assert_eq!(one_shot, all);
    }

    // ========================================================================
    // finish(): EOF validation
    // ========================================================================

    /// With no carryover there is nothing to flush.
    #[test]
    fn finish_returns_none_when_no_bytes_were_carried_over() {
        let mut state = BoundaryState::new_no_header();
        state.find_boundaries(&records(&[&[1u8; 8]])).expect("scans");
        assert!(state.finish().expect("no leftover").is_none());
    }

    /// `finish` flushes carried-over bytes when they happen to form whole
    /// records. That is narrower than it sounds, and the narrowness is the
    /// point: `find_boundaries` consumes *every* complete record it can see, so
    /// after a successful scan the carryover is by construction an INCOMPLETE
    /// record — making `None` or an error the normal outcomes at EOF.
    ///
    /// The one reachable path to a flushed batch is a stream that ended before
    /// the header could even be resolved (< 8 bytes, so `parse_header_size`
    /// never ran) whose bytes nonetheless frame cleanly. Pinning it keeps the
    /// flush path honest and documents why it is nearly dead code.
    #[test]
    fn finish_flushes_carried_over_bytes_that_form_whole_records() {
        // Seven bytes: too short for magic + `l_text`, so the whole block is
        // carried over unparsed — but a clean `block_size = 3` + 3-byte body.
        let block = record(&[7u8; 3]);
        assert_eq!(block.len(), 7);

        let mut state = BoundaryState::new();
        let batch = state.find_boundaries(&block).expect("short input is not an error");
        assert!(batch.buffer.is_empty(), "the header was never resolved, so nothing is emitted");

        let flushed =
            state.finish().expect("carryover frames cleanly").expect("a batch is flushed");
        assert_eq!(flushed.buffer, block);
        assert_eq!(flushed.offsets, vec![0, 7]);
    }

    /// After a successful scan the carryover is always a partial record, so the
    /// normal EOF outcome is `None` — the complete records were already emitted
    /// by `find_boundaries` itself.
    #[test]
    fn finish_returns_none_after_a_scan_that_consumed_every_record() {
        let payloads: [&[u8]; 2] = [&[1u8; 8], &[2u8; 8]];
        let all = records(&payloads);

        let mut state = BoundaryState::new_no_header();
        let batch = state.find_boundaries(&all[..12]).expect("scans");
        assert_eq!(batch.buffer, record(&[1u8; 8]), "the first record is emitted immediately");

        assert!(
            state.finish().expect("no carryover").is_none(),
            "a scan that ended on a record boundary leaves nothing to flush"
        );
    }

    /// A truncated record at EOF is corruption and must surface as
    /// `UnexpectedEof` rather than silently dropping bytes. Two distinct
    /// truncations reach two different error sites: a declared `block_size`
    /// that outruns the buffer, and trailing bytes too short to even hold the
    /// 4-byte length prefix.
    #[rstest]
    #[case::body_shorter_than_declared_block_size(&[0x20, 0, 0, 0, 1, 2, 3])]
    #[case::three_trailing_bytes(&[1, 2, 3])]
    #[case::one_trailing_byte(&[9])]
    fn finish_rejects_an_incomplete_record_at_eof(#[case] tail: &[u8]) {
        let mut block = record(&[1u8; 8]);
        block.extend_from_slice(tail);

        let mut state = BoundaryState::new_no_header();
        let batch = state.find_boundaries(&block).expect("the complete record scans");
        assert_eq!(batch.buffer, record(&[1u8; 8]));

        let err = state.finish().expect_err("a truncated tail must not be dropped");
        assert_eq!(err.kind(), io::ErrorKind::UnexpectedEof);
    }

    /// `finish` takes the leftover, so a second call reports nothing remaining
    /// rather than re-emitting the same bytes. `FindBamBoundaries` guards this
    /// with its own `finalized` flag, but the state machine must not depend on
    /// that guard for correctness of the buffer it returns.
    #[test]
    fn finish_does_not_re_emit_the_same_leftover_twice() {
        // Same short-input setup as the flush test above — the only shape that
        // leaves whole records parked in carryover.
        let block = record(&[7u8; 3]);
        let mut state = BoundaryState::new();
        state.find_boundaries(&block).expect("short input buffers");

        assert!(state.finish().expect("first flush").is_some());
        assert!(state.finish().expect("second flush").is_none(), "leftover must be consumed once");
    }

    // ========================================================================
    // Bounded carryover on corrupt framing
    // ========================================================================

    /// A corrupt `block_size` must be rejected once the carry passes
    /// [`MAX_CARRY_BYTES`], rather than buffering the rest of the stream.
    ///
    /// Without the bound this terminates *correctly* — `finish` catches the
    /// shortfall at EOF — but only after carrying every remaining byte, so peak
    /// memory tracks input size and a large corrupt file OOMs before it can
    /// print the diagnostic it was about to print.
    #[test]
    fn a_corrupt_block_size_is_rejected_once_the_carry_passes_the_bound() {
        let mut state = BoundaryState::new_no_header();
        // Declare a body far larger than any real record, and never supply it.
        let mut first =
            u32::try_from(MAX_CARRY_BYTES * 4).expect("fits u32").to_le_bytes().to_vec();
        first.extend_from_slice(&[7u8; 64]);
        state.find_boundaries(&first).expect("first block just starts the carry");

        // Feed blocks until the bound trips. It must trip well before the
        // declared size is reached.
        let block = vec![9u8; 1024 * 1024];
        let mut fed = first.len();
        let err = loop {
            match state.find_boundaries(&block) {
                Ok(_) => {
                    fed += block.len();
                    assert!(
                        fed <= MAX_CARRY_BYTES + block.len(),
                        "carry ran past the bound without erroring ({fed} bytes fed)",
                    );
                }
                Err(e) => break e,
            }
        };
        assert_eq!(err.kind(), io::ErrorKind::InvalidData);
    }

    /// A header whose declared `l_text` never resolves must be bounded the same
    /// way — it is the other path that carries the whole stream forward.
    #[test]
    fn a_header_that_never_resolves_is_rejected_once_the_carry_passes_the_bound() {
        let mut header = fgumi_raw_bam::BAM_MAGIC.to_vec();
        header.extend_from_slice(&u32::MAX.to_le_bytes()); // absurd l_text
        let mut state = BoundaryState::new();
        state.find_boundaries(&header).expect("first block just starts the carry");

        let block = vec![0u8; 1024 * 1024];
        let mut fed = header.len();
        let err = loop {
            match state.find_boundaries(&block) {
                Ok(_) => {
                    fed += block.len();
                    assert!(
                        fed <= MAX_CARRY_BYTES + block.len(),
                        "header carry ran past the bound without erroring ({fed} bytes fed)",
                    );
                }
                Err(e) => break e,
            }
        };
        assert_eq!(err.kind(), io::ErrorKind::InvalidData);
    }

    /// The bound must not reject legitimately large records. An ONT ultra-long
    /// read with methylation tags can reach a few MB, and such a record spans
    /// many BGZF blocks — so it is carried across dozens of calls before it
    /// completes. That is exactly the shape the bound must let through, which is
    /// why it sits far above any real record rather than at the queue budget.
    #[test]
    fn a_multi_megabyte_record_spanning_many_blocks_is_still_accepted() {
        const BODY: usize = 8 * 1024 * 1024; // ~8 MB — above any real BAM record
        // Compile-time: if the bound is ever lowered under this fixture, the
        // test must fail to build rather than silently stop proving anything.
        const { assert!(BODY < MAX_CARRY_BYTES) };

        let mut stream = u32::try_from(BODY).expect("fits u32").to_le_bytes().to_vec();
        stream.extend_from_slice(&vec![3u8; BODY]);

        let mut state = BoundaryState::new_no_header();
        let mut emitted = Vec::new();
        for chunk in stream.chunks(64 * 1024) {
            emitted.extend_from_slice(
                &state.find_boundaries(chunk).expect("a large but legal record must scan").buffer,
            );
        }
        assert!(state.finish().expect("nothing left over").is_none());
        assert_eq!(emitted, stream, "the whole record must survive the carry intact");
    }

    // -- Queryname cut (BatchCut::Queryname) --

    fn named_record(name: &[u8]) -> Vec<u8> {
        let l_read_name = name.len() + 1;
        let mut body = vec![0u8; 32];
        body[8] = u8::try_from(l_read_name).expect("name < 255");
        body.extend_from_slice(name);
        body.push(0);
        record(&body)
    }

    fn qn_block(records: &[Vec<u8>]) -> Vec<u8> {
        records.iter().flatten().copied().collect()
    }

    fn queryname_cut_names(
        input_blocks: &[Vec<u8>],
        min_emit_bytes: usize,
    ) -> io::Result<Vec<Vec<Vec<u8>>>> {
        let mut state =
            BoundaryState::new_no_header().with_cut(BatchCut::Queryname, min_emit_bytes);
        let mut out: Vec<Vec<Vec<u8>>> = Vec::new();
        let collect = |b: &BoundaryBatch, out: &mut Vec<Vec<Vec<u8>>>| {
            if b.offsets.len() <= 1 {
                return;
            }
            let mut names = Vec::new();
            for i in 0..b.offsets.len() - 1 {
                let body = &b.buffer[b.offsets[i] + 4..b.offsets[i + 1]];
                names.push(fgumi_raw_bam::read_name(body).to_vec());
            }
            out.push(names);
        };
        for blk in input_blocks {
            let batch = state.find_boundaries(blk)?;
            collect(&batch, &mut out);
        }
        if let Some(batch) = state.finish()? {
            collect(&batch, &mut out);
        }
        Ok(out)
    }

    fn qn_flat(batches: &[Vec<Vec<u8>>]) -> Vec<Vec<u8>> {
        batches.iter().flatten().cloned().collect()
    }

    fn qn_no_straddle(batches: &[Vec<Vec<u8>>]) {
        for w in batches.windows(2) {
            assert_ne!(w[0].last(), w[1].first(), "a run straddled a batch boundary");
        }
    }

    #[test]
    fn queryname_cut_carries_a_run_split_across_two_blocks() {
        let b0 = qn_block(&[named_record(b"a"), named_record(b"a"), named_record(b"b")]);
        let b1 = qn_block(&[named_record(b"b"), named_record(b"c")]);
        let batches = queryname_cut_names(&[b0, b1], 0).expect("cut ok");
        qn_no_straddle(&batches);
        assert_eq!(
            qn_flat(&batches),
            vec![b"a".to_vec(), b"a".to_vec(), b"b".to_vec(), b"b".to_vec(), b"c".to_vec()],
        );
    }

    #[test]
    fn queryname_cut_whole_block_is_one_run_carries_until_it_ends() {
        let b0 = qn_block(&[named_record(b"x"), named_record(b"x")]);
        let b1 = qn_block(&[named_record(b"x")]);
        let b2 = qn_block(&[named_record(b"y")]);
        let batches = queryname_cut_names(&[b0, b1, b2], 0).expect("cut ok");
        assert_eq!(
            qn_flat(&batches),
            vec![b"x".to_vec(), b"x".to_vec(), b"x".to_vec(), b"y".to_vec()],
        );
        qn_no_straddle(&batches);
    }

    #[test]
    fn queryname_cut_run_boundary_exactly_at_block_end() {
        let b0 = qn_block(&[named_record(b"a"), named_record(b"a")]);
        let b1 = qn_block(&[named_record(b"b"), named_record(b"b")]);
        let batches = queryname_cut_names(&[b0, b1], 0).expect("cut ok");
        assert_eq!(
            qn_flat(&batches),
            vec![b"a".to_vec(), b"a".to_vec(), b"b".to_vec(), b"b".to_vec()],
        );
        qn_no_straddle(&batches);
    }

    #[test]
    fn queryname_cut_single_record_templates() {
        let b0 = qn_block(&[named_record(b"r0"), named_record(b"r1"), named_record(b"r2")]);
        let batches = queryname_cut_names(&[b0], 0).expect("cut ok");
        assert_eq!(qn_flat(&batches), vec![b"r0".to_vec(), b"r1".to_vec(), b"r2".to_vec()]);
    }

    #[test]
    fn queryname_cut_coalesces_to_min_emit_bytes() {
        let blocks: Vec<Vec<u8>> =
            (0..8).map(|i| qn_block(&[named_record(format!("q{i}").as_bytes())])).collect();
        let batches = queryname_cut_names(&blocks, 1 << 20).expect("cut ok");
        assert_eq!(batches.len(), 1, "coalescing should hold until finish()");
        assert_eq!(batches[0].len(), 8);
    }

    #[test]
    fn queryname_cut_preserves_all_bytes_across_boundaries() {
        let b0 = qn_block(&[named_record(b"long"), named_record(b"long"), named_record(b"long")]);
        let b1 = qn_block(&[named_record(b"long"), named_record(b"solo")]);
        let batches = queryname_cut_names(&[b0, b1], 0).expect("cut ok");
        assert_eq!(
            qn_flat(&batches),
            vec![
                b"long".to_vec(),
                b"long".to_vec(),
                b"long".to_vec(),
                b"long".to_vec(),
                b"solo".to_vec(),
            ],
        );
        qn_no_straddle(&batches);
    }

    #[test]
    fn queryname_cut_rejects_a_record_too_short_for_its_name() {
        let short_body = record(&[0u8; 16]);
        let mut state = BoundaryState::new_no_header().with_cut(BatchCut::Queryname, 0);
        let err = state.find_boundaries(&short_body).expect_err("must reject");
        assert_eq!(err.kind(), io::ErrorKind::InvalidData);
    }

    #[test]
    fn queryname_cut_rejects_an_unbounded_open_run() {
        // Every record shares one name, so the run never closes: `qn_buffer`
        // would grow for the whole stream. With a small ceiling the cutter
        // rejects it as `InvalidData` instead of buffering without bound.
        let mut state = BoundaryState::new_no_header()
            .with_cut(BatchCut::Queryname, 0)
            .with_max_queryname_run_bytes(64);
        let mut err = None;
        for _ in 0..64 {
            if let Err(e) = state.find_boundaries(&qn_block(&[named_record(b"same")])) {
                err = Some(e);
                break;
            }
        }
        let err = err.expect("an unbounded open run must be rejected");
        assert_eq!(err.kind(), io::ErrorKind::InvalidData);
        assert!(err.to_string().contains("queryname run"), "got: {err}");
    }

    #[test]
    fn queryname_cut_cap_measures_only_the_open_run_not_the_coalesced_prefix() {
        // A large closed prefix (distinct names, held back by a big coalescing
        // target) plus a small trailing open run: the whole accumulator exceeds
        // the cap, but the open run alone is well under it. The cap must measure
        // only the open run, so this valid queryname-grouped input is accepted.
        let cap = 200;
        let mut state = BoundaryState::new_no_header()
            .with_cut(BatchCut::Queryname, 10_000) // coalescing holds the closed prefix
            .with_max_queryname_run_bytes(cap);
        // 8 distinct-name records (~39 bytes each => ~312 bytes) exceed `cap`,
        // while every open run is a single ~39-byte record well under it.
        let names: Vec<Vec<u8>> =
            (0..8).map(|i| named_record(format!("q{i}").as_bytes())).collect();
        let batch = state
            .find_boundaries(&qn_block(&names))
            .expect("must not reject: the open run is under the cap");
        assert!(batch.offsets.len() <= 1, "coalescing holds everything: nothing emitted yet");
        // EOF force-flush emits the whole run without applying the cap.
        let flushed = state.finish().expect("finish ok").expect("a final batch");
        assert_eq!(flushed.offsets.len() - 1, 8, "all 8 records flushed at EOF");
    }

    #[test]
    fn queryname_cut_rejects_a_run_of_empty_read_names() {
        // A run of empty BAM names (`l_read_name == 1`, just the NUL) also never
        // closes; the same ceiling must bound it.
        let mut state = BoundaryState::new_no_header()
            .with_cut(BatchCut::Queryname, 0)
            .with_max_queryname_run_bytes(64);
        let mut err = None;
        for _ in 0..64 {
            if let Err(e) = state.find_boundaries(&qn_block(&[named_record(b"")])) {
                err = Some(e);
                break;
            }
        }
        let err = err.expect("a run of empty read names must be rejected");
        assert_eq!(err.kind(), io::ErrorKind::InvalidData);
    }

    /// Concatenate the emitted record bytes across all batches (including the
    /// `finish` flush) so a test can assert byte-identity against the input.
    fn queryname_cut_bytes(input_blocks: &[Vec<u8>], min_emit_bytes: usize) -> Vec<u8> {
        let mut state =
            BoundaryState::new_no_header().with_cut(BatchCut::Queryname, min_emit_bytes);
        let mut out = Vec::new();
        for blk in input_blocks {
            let batch = state.find_boundaries(blk).expect("cut ok");
            out.extend_from_slice(&batch.buffer);
        }
        if let Some(batch) = state.finish().expect("finish ok") {
            out.extend_from_slice(&batch.buffer);
        }
        out
    }

    /// Regression for the O(n^2) carry re-scan: feed many small blocks under a
    /// coalescing target large enough that several blocks accumulate before an
    /// emit, forcing repeated drain/rebase of the trailing run in the
    /// accumulator. Asserts the emitted bytes are exactly the input bytes, in
    /// order, and no queryname run straddles a batch. The old implementation
    /// re-scanned the entire carry on every call (correct output but quadratic);
    /// this proves the accumulator path stays correct across many emits.
    #[test]
    fn queryname_cut_many_blocks_coalesce_without_losing_bytes() {
        // 60 templates, each a 2-record run (paired reads), one template per
        // block — so runs never straddle a single block but coalescing must
        // hold several blocks before each emit.
        let mut blocks = Vec::new();
        let mut expected = Vec::new();
        for i in 0..60u32 {
            let name = format!("read{i:04}");
            let r0 = named_record(name.as_bytes());
            let r1 = named_record(name.as_bytes());
            expected.extend_from_slice(&r0);
            expected.extend_from_slice(&r1);
            blocks.push(qn_block(&[r0, r1]));
        }
        // Coalesce ~4 blocks per emit.
        let one_block_len = blocks[0].len();
        let bytes = queryname_cut_bytes(&blocks, one_block_len * 4);
        assert_eq!(bytes, expected, "coalescing dropped or reordered bytes");

        // And the batch structure is closed under queryname.
        let batches = queryname_cut_names(&blocks, one_block_len * 4).expect("cut ok");
        qn_no_straddle(&batches);
        let flat = qn_flat(&batches);
        assert_eq!(flat.len(), 120, "every record emitted exactly once");
    }

    /// A run that spans several coalesced blocks is still carried whole: the
    /// accumulator must not emit a partial run at a coalescing boundary. Ten
    /// blocks all one run, then a block that ends it.
    #[test]
    fn queryname_cut_long_run_spanning_many_blocks_carried_whole() {
        let mut blocks: Vec<Vec<u8>> =
            (0..10).map(|_| qn_block(&[named_record(b"same")])).collect();
        blocks.push(qn_block(&[named_record(b"next")]));
        // A small coalescing target that a single block already exceeds — the
        // run must still be held until it ends, never split mid-run.
        let batches = queryname_cut_names(&blocks, 1).expect("cut ok");
        qn_no_straddle(&batches);
        let flat = qn_flat(&batches);
        assert_eq!(flat.iter().filter(|n| n.as_slice() == b"same").count(), 10);
        assert_eq!(flat.iter().filter(|n| n.as_slice() == b"next").count(), 1);
    }

    /// A header that never completes before EOF leaves the partial header in
    /// `leftover`; `finish` then reports it as corruption rather than emitting
    /// header bytes as if they were records.
    #[test]
    fn finish_reports_a_header_that_never_completed_as_corruption() {
        let header = bam_header("@HD\tVN:1.6\n", &[("chr1", 1000)]);
        let mut state = BoundaryState::new();
        let batch = state.find_boundaries(&header[..12]).expect("partial header buffers");
        assert!(batch.buffer.is_empty());

        let err = state.finish().expect_err("a partial header is not a valid record stream");
        assert_eq!(err.kind(), io::ErrorKind::UnexpectedEof);
    }
}
