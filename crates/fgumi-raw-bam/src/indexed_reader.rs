//! Indexed BAM reader that yields raw bytes without decoding to `RecordBuf`.
//!
//! [`IndexedRawBamReader`] delegates index parsing and BGZF seek to noodles,
//! then yields [`RawRecord`] directly — bypassing the noodles `Record` type
//! whose inner bytes are not publicly accessible (noodles#373).
//!
//! # Feature gate
//!
//! This module is only compiled when the `noodles` feature is enabled, because
//! it depends on noodles for BGZF seek, BAI index parsing, and SAM header
//! decoding.

use std::io::{self, Read, Seek};

use anyhow::Context as _;
use noodles::bam::bai;
use noodles::bgzf;
use noodles::bgzf::io::Seek as BgzfSeek;
use noodles::core::region::Interval;
use noodles::csi::BinningIndex as _;
use noodles::sam;

use crate::raw_bam_record::{RawRecord, read_raw_record};

// ─── public type ────────────────────────────────────────────────────────────

/// An indexed BAM reader that yields [`RawRecord`] bytes.
///
/// Uses a BAI index to compute the BGZF virtual-offset chunks that overlap a
/// query region, then seeks the BGZF stream and reads raw BAM records from
/// those chunks.  Per-record overlap is verified via a cheap raw-byte
/// intersection test before yielding.
///
/// # Constructing
///
/// ```rust,ignore
/// use fgumi_raw_bam::IndexedRawBamReader;
/// use noodles::bam::bai;
///
/// let index = bai::fs::read("sample.bam.bai")?;
/// let mut reader = IndexedRawBamReader::from_path("sample.bam", index)?;
/// let header = reader.read_header()?;
/// let region = "chr1:100-200".parse()?;
/// for result in reader.query(&header, &region)? {
///     let raw = result?;
///     // use raw bytes …
/// }
/// ```
pub struct IndexedRawBamReader<R> {
    /// Noodles BAM reader wrapping the BGZF decompressor.
    inner: noodles::bam::io::Reader<bgzf::io::Reader<R>>,
    /// BAI (or CSI) index held as a boxed `BinningIndex`.
    index: Box<dyn noodles::csi::BinningIndex>,
    /// BGZF virtual position of the first record (captured after `read_header`).
    /// Used as the `query_unmapped` fallback when the index lacks a
    /// `last_first_record_start_position` marker, so the iterator rewinds to
    /// the first record rather than resuming at a mid-stream cursor.
    first_record_position: Option<bgzf::VirtualPosition>,
}

// ─── convenience constructors ───────────────────────────────────────────────

impl IndexedRawBamReader<std::fs::File> {
    /// Open a BAM file from a path and pair it with a BAI index.
    ///
    /// The BGZF layer is created automatically from the file handle.
    ///
    /// # Errors
    ///
    /// Returns an error if the file cannot be opened.
    pub fn from_path(path: impl AsRef<std::path::Path>, index: bai::Index) -> io::Result<Self> {
        let file = std::fs::File::open(path)?;
        Ok(Self::new(file, index))
    }
}

// ─── generic impl ───────────────────────────────────────────────────────────

impl<R: Read + Seek> IndexedRawBamReader<R> {
    /// Create from any seekable readable stream and a BAI index.
    #[must_use]
    pub fn new(reader: R, index: bai::Index) -> Self {
        Self {
            inner: noodles::bam::io::Reader::new(reader),
            index: Box::new(index),
            first_record_position: None,
        }
    }
}

impl<R: Read> IndexedRawBamReader<R> {
    /// Read (and consume) the BAM header.
    ///
    /// Must be called exactly once before [`query`], because the header
    /// encodes the reference-sequence dictionary needed to resolve region
    /// names to reference-sequence IDs.
    ///
    /// # Errors
    ///
    /// Returns an error if the header cannot be read or parsed.
    pub fn read_header(&mut self) -> anyhow::Result<sam::Header> {
        let header = self.inner.read_header().context("reading BAM header")?;
        // Remember the post-header BGZF position so `query_unmapped` can rewind
        // to the first record when the index lacks a location marker.
        self.first_record_position = Some(self.inner.get_ref().virtual_position());
        Ok(header)
    }
}

impl<R: Read + Seek> IndexedRawBamReader<R> {
    /// Return an iterator over unmapped [`RawRecord`]s.
    ///
    /// Seeks to the position of the last placed read's start (or the first
    /// record if the index has no such marker), then yields only records
    /// flagged as unmapped (`FLAG & 4`).
    ///
    /// # Errors
    ///
    /// - Seek or read failures.
    /// - The index has no `last_first_record_start_position` AND
    ///   [`read_header`](Self::read_header) has not been called — without
    ///   either, the first-record position is unknown and falling back to
    ///   BGZF offset 0 would attempt to parse the BAM header as a record.
    pub fn query_unmapped(&mut self) -> anyhow::Result<UnmappedIter<'_, R>> {
        // Prefer the index marker; otherwise use the first-record position
        // captured during `read_header`. Refusing to guess prevents
        // UnmappedIter from parsing header bytes as a record.
        let virtual_pos = self
            .index
            .last_first_record_start_position()
            .or(self.first_record_position)
            .ok_or_else(|| {
                anyhow::anyhow!(
                    "query_unmapped: index has no unmapped start marker and \
                     read_header() has not been called",
                )
            })?;
        self.inner
            .get_mut()
            .seek_to_virtual_position(virtual_pos)
            .context("seeking to unmapped position")?;
        Ok(UnmappedIter { inner: &mut self.inner, record: RawRecord::new() })
    }

    /// Return an iterator over [`RawRecord`]s that overlap `region`.
    ///
    /// The iterator yields only records whose alignment interval intersects
    /// the query region; records that fall entirely outside the region (but
    /// in the same BGZF chunk) are silently skipped.
    ///
    /// # Arguments
    ///
    /// * `header` — SAM header obtained from [`read_header`]; used to
    ///   resolve the region name to a reference-sequence ID.
    /// * `region` — Genomic region to query (e.g. `"chr1:100-200".parse()`).
    ///
    /// # Errors
    ///
    /// Returns an error if the region name is not found in the header or if
    /// the index cannot be queried.
    pub fn query<'a>(
        &'a mut self,
        header: &'a sam::Header,
        region: &noodles::core::Region,
    ) -> anyhow::Result<RawQueryIter<'a, R>> {
        // Resolve the region name to a 0-based reference-sequence index.
        let ref_id = header.reference_sequences().get_index_of(region.name()).ok_or_else(|| {
            anyhow::anyhow!("region reference sequence {:?} not found in BAM header", region.name())
        })?;

        let interval = region.interval();

        // Ask the index for the BGZF chunks covering this interval.
        let chunks = self.index.query(ref_id, interval).context("querying BAI index")?;

        // Hand `&mut bgzf_reader` to `csi::io::Query` so it can seek between
        // chunks and feed us decompressed data in order.
        let bgzf_reader = self.inner.get_mut();
        let csi_query = noodles::csi::io::Query::new(bgzf_reader, chunks);

        Ok(RawQueryIter { reader: csi_query, ref_id, interval, record: RawRecord::new() })
    }

    /// Return an iterator over [`RawRecord`]s overlapping **any** of `regions`,
    /// in a single coordinate-ordered pass, with each record yielded exactly
    /// once.
    ///
    /// This is the multi-interval analog of [`query`](Self::query) and mirrors
    /// htsjdk's `SamReader.query(QueryInterval[])` (and therefore fgbio's
    /// `SamSource.query(loci)`): the BAI chunks for every region are gathered
    /// and [`merge_chunks`](noodles::csi::binning_index::merge_chunks)-ed into a
    /// minimal, sorted, non-overlapping set, then scanned once. Because the
    /// merged chunks are disjoint and read in virtual-position order (which, for
    /// a coordinate-sorted BAM, is genomic-coordinate order), a record
    /// overlapping several regions is visited a **single** time and records are
    /// emitted in coordinate order. The caller therefore needs no external
    /// buffering, sort, or de-duplication.
    ///
    /// `regions` may overlap, be unsorted, and span multiple reference
    /// sequences; the output ordering and single-visit guarantee hold
    /// regardless.
    ///
    /// # Arguments
    ///
    /// * `header` — SAM header obtained from [`read_header`](Self::read_header);
    ///   resolves each region name to a reference-sequence ID.
    /// * `regions` — genomic regions to query.
    ///
    /// # Errors
    ///
    /// Returns an error if a region name is not found in the header or if the
    /// index cannot be queried.
    pub fn query_intervals<'a>(
        &'a mut self,
        header: &'a sam::Header,
        regions: &[noodles::core::Region],
    ) -> anyhow::Result<RawMultiQueryIter<'a, R>> {
        use noodles::csi::binning_index::index::reference_sequence::bin::Chunk;
        use std::collections::BTreeMap;

        // Resolve every region to a reference ID + interval, grouping the
        // intervals per reference.
        let mut intervals_by_ref: BTreeMap<usize, Vec<Interval>> = BTreeMap::new();
        for region in regions {
            let ref_id =
                header.reference_sequences().get_index_of(region.name()).ok_or_else(|| {
                    anyhow::anyhow!(
                        "region reference sequence {:?} not found in BAM header",
                        region.name()
                    )
                })?;
            intervals_by_ref.entry(ref_id).or_default().push(region.interval());
        }

        // Optimize each reference's intervals into a sorted, non-overlapping set
        // (htsjdk `QueryInterval.optimizeIntervals`), then gather the BAI chunks
        // for every optimized interval. Optimizing first shrinks both the index
        // queries and the per-record overlap filter below.
        let mut chunks: Vec<Chunk> = Vec::new();
        for (ref_id, intervals) in &mut intervals_by_ref {
            *intervals = merge_intervals(std::mem::take(intervals));
            for interval in intervals.iter() {
                chunks.extend(self.index.query(*ref_id, *interval).context("querying BAI index")?);
            }
        }

        // Merge into a minimal, sorted, non-overlapping chunk set so the scan
        // reads each covered BGZF region exactly once. Overlapping chunks from
        // adjacent regions collapse here, which is what makes a read spanning
        // several regions surface only once (htsjdk `optimizeIntervals` analog).
        let merged = noodles::csi::binning_index::merge_chunks(&chunks);

        let bgzf_reader = self.inner.get_mut();
        let csi_query = noodles::csi::io::Query::new(bgzf_reader, merged);

        Ok(RawMultiQueryIter { reader: csi_query, intervals_by_ref, record: RawRecord::new() })
    }
}

// ─── iterator ───────────────────────────────────────────────────────────────

/// Iterator over raw BAM records overlapping a queried region.
///
/// Created by [`IndexedRawBamReader::query`].  Reuses a single
/// [`RawRecord`] buffer internally; each `next()` call clones the bytes
/// into a new owned [`RawRecord`] that is yielded to the caller.
pub struct RawQueryIter<'r, R> {
    reader: noodles::csi::io::Query<'r, bgzf::io::Reader<R>>,
    ref_id: usize,
    interval: Interval,
    record: RawRecord,
}

impl<R: Read + Seek> Iterator for RawQueryIter<'_, R> {
    type Item = anyhow::Result<RawRecord>;

    fn next(&mut self) -> Option<Self::Item> {
        loop {
            match read_raw_record(&mut self.reader, &mut self.record) {
                Ok(0) => return None, // EOF / all chunks consumed
                Ok(_) => {
                    match intersects_region(&self.record, self.ref_id, self.interval) {
                        Ok(true) => {
                            // Yield a clone so the caller owns the bytes.
                            return Some(Ok(self.record.clone()));
                        }
                        Ok(false) => {} // record outside query region; skip
                        Err(e) => return Some(Err(e)),
                    }
                }
                Err(e) => return Some(Err(e.into())),
            }
        }
    }
}

// ─── multi-interval iterator ─────────────────────────────────────────────────

/// Iterator over raw BAM records overlapping any of several queried regions.
///
/// Created by [`IndexedRawBamReader::query_intervals`]. Scans the merged BAI
/// chunk set for all queried regions in a single virtual-position-ordered pass,
/// yielding each overlapping record exactly once and in coordinate order. Reuses
/// a single [`RawRecord`] buffer internally; each `next()` clones the yielded
/// bytes into a new owned [`RawRecord`].
pub struct RawMultiQueryIter<'r, R> {
    reader: noodles::csi::io::Query<'r, bgzf::io::Reader<R>>,
    /// Queried intervals grouped by reference-sequence ID. A record is yielded
    /// only when it overlaps one of the intervals for its own reference ID; the
    /// chunk scan can surface records outside every queried region (chunk
    /// granularity) or on a reference with no queried region, which are skipped.
    intervals_by_ref: std::collections::BTreeMap<usize, Vec<Interval>>,
    record: RawRecord,
}

impl<R: Read + Seek> Iterator for RawMultiQueryIter<'_, R> {
    type Item = anyhow::Result<RawRecord>;

    fn next(&mut self) -> Option<Self::Item> {
        loop {
            match read_raw_record(&mut self.reader, &mut self.record) {
                Ok(0) => return None, // EOF / all chunks consumed
                Ok(_) => {
                    // Validate the record's structure *before* touching any
                    // fixed-offset field: `record_span_1based` fails closed on a
                    // record shorter than the 32-byte core, so a corrupt BGZF block
                    // declaring a tiny nonzero size cannot reach the `ref_id`/`pos`
                    // slice reads below and panic. Returns `Ok(None)` for unmapped
                    // records (no reference span, cannot overlap).
                    let (rec_start, rec_end) = match record_span_1based(&self.record) {
                        Ok(Some(span)) => span,
                        Ok(None) => continue, // no reference span (unmapped) cannot overlap
                        Err(e) => return Some(Err(e)), // malformed record: fail closed
                    };

                    // The fixed core is now known to be present, so reading `ref_id`
                    // can no longer index past the buffer. Keep the record only if it
                    // overlaps one of the intervals queried on that reference.
                    let record_ref_id = crate::fields::ref_id(&self.record);
                    let Ok(ref_id) = usize::try_from(record_ref_id) else {
                        continue; // unmapped / negative ref id cannot overlap a region
                    };
                    let Some(intervals) = self.intervals_by_ref.get(&ref_id) else {
                        continue; // record on a reference with no queried region
                    };

                    // `intervals` is sorted by start and non-overlapping (`merge_intervals`),
                    // so its ends are strictly increasing. Binary-search for the first
                    // interval that could reach the record (its end `>= rec_start`); the
                    // record overlaps iff that interval's start is `<= rec_end`. This scans
                    // O(log intervals) rather than every queried interval per record.
                    let idx = intervals.partition_point(|iv| interval_end(iv) < rec_start);
                    let overlaps =
                        intervals.get(idx).is_some_and(|iv| interval_start(iv) <= rec_end);
                    if overlaps {
                        return Some(Ok(self.record.clone()));
                    }
                }
                Err(e) => return Some(Err(e.into())),
            }
        }
    }
}

// ─── unmapped iterator ──────────────────────────────────────────────────────

/// Iterator over unmapped raw BAM records.
///
/// Created by [`IndexedRawBamReader::query_unmapped`].
pub struct UnmappedIter<'r, R> {
    inner: &'r mut noodles::bam::io::Reader<bgzf::io::Reader<R>>,
    record: RawRecord,
}

impl<R: Read> Iterator for UnmappedIter<'_, R> {
    type Item = anyhow::Result<RawRecord>;

    fn next(&mut self) -> Option<Self::Item> {
        loop {
            match read_raw_record(self.inner.get_mut(), &mut self.record) {
                Ok(0) => return None,
                Ok(_) => {
                    if self.record.is_unmapped() {
                        return Some(Ok(self.record.clone()));
                    }
                    // skip mapped records
                }
                Err(e) => return Some(Err(e.into())),
            }
        }
    }
}

// ─── interval optimization ───────────────────────────────────────────────────

/// Inclusive 1-based start of `interval` as a plain `usize` (unbounded → 1).
fn interval_start(interval: &Interval) -> usize {
    interval.start().map_or(1, usize::from)
}

/// Inclusive 1-based end of `interval` as a plain `usize` (unbounded → `usize::MAX`).
fn interval_end(interval: &Interval) -> usize {
    interval.end().map_or(usize::MAX, usize::from)
}

/// Build an [`Interval`] from 1-based `start`/`end` positions.
///
/// `end == usize::MAX` is the sentinel [`interval_end`] returns for an unbounded
/// end (whole-contig `chr1` / right-open `chr1:100-` queries). It is preserved as
/// a right-open [`Interval`] rather than round-tripped through `Position`, which
/// would pin a finite end and stop those queries from covering the tail of a
/// reference.
fn make_interval(start: usize, end: usize) -> Interval {
    use noodles::core::Position;
    let start = Position::try_from(start.max(1)).expect("start >= 1 is a valid Position");
    if end == usize::MAX {
        return Interval::from(start..);
    }
    let end = Position::try_from(end.max(1)).expect("end >= 1 is a valid Position");
    Interval::from(start..=end)
}

/// Sort and merge a reference's intervals into a minimal, non-overlapping set,
/// mirroring htsjdk's `QueryInterval.optimizeIntervals`.
///
/// Two intervals merge when they overlap **or abut** — i.e. the next interval's
/// start is `<= end + 1` of the running interval (1-based inclusive coordinates),
/// so `[1520,1521]` and `[1522,1525]` collapse to `[1520,1525]` while
/// `[1520,1521]` and `[1523,1525]` (a gap at 1522) stay separate. Optimizing
/// keeps the per-record overlap filter small and matches the disjoint-interval
/// precondition htsjdk's multi-interval query assumes.
fn merge_intervals(mut intervals: Vec<Interval>) -> Vec<Interval> {
    intervals.sort_by_key(|iv| (interval_start(iv), interval_end(iv)));

    let mut merged: Vec<Interval> = Vec::with_capacity(intervals.len());
    for interval in intervals {
        let (start, end) = (interval_start(&interval), interval_end(&interval));
        if let Some(last) = merged.last_mut() {
            let last_end = interval_end(last);
            if start <= last_end.saturating_add(1) {
                // Overlapping or abutting: extend the running interval if needed.
                if end > last_end {
                    *last = make_interval(interval_start(last), end);
                }
                continue;
            }
        }
        merged.push(make_interval(start, end));
    }
    merged
}

// ─── overlap test ───────────────────────────────────────────────────────────

/// Compute a mapped record's 1-based inclusive reference span `[start, end]`.
///
/// Returns:
/// - `Ok(Some((start, end)))` for a mapped record with a valid reference span. A
///   record with a zero-length reference span (unmapped mate placed at a position,
///   or no CIGAR) is treated as covering its single start position, so its span
///   equals `intersects_region`'s single-position fallback.
/// - `Ok(None)` when the record has no reference position (unmapped, or a negative
///   `pos`) and therefore cannot overlap any queried interval.
/// - `Err(..)` when the record is malformed — shorter than the fixed core or its
///   declared CIGAR extends past the record buffer. Rather than silently treating a
///   truncated record as a point alignment at its start, the span cannot be trusted,
///   so callers fail closed and propagate the error.
fn record_span_1based(record: &RawRecord) -> anyhow::Result<Option<(usize, usize)>> {
    use crate::cigar::reference_length_from_raw_bam_checked;
    use crate::fields::pos;

    let bytes: &[u8] = record;
    // Validate the record's structure (fixed core + read-name + CIGAR all within
    // bounds) before trusting any field; this also guarantees `pos` is in bounds.
    let ref_len_raw = reference_length_from_raw_bam_checked(bytes).context(
        "malformed BAM record: truncated, or declared CIGAR/read-name length \
             exceeds record buffer",
    )?;

    let p = pos(bytes);
    if p < 0 {
        return Ok(None); // unmapped / no position cannot overlap a region
    }
    #[allow(clippy::cast_sign_loss)] // guarded by `p >= 0` check above
    let start_1based = p as usize + 1;

    #[allow(clippy::cast_sign_loss)] // reference_length_from_raw_bam_checked returns non-negative
    let ref_len = ref_len_raw as usize;
    let end_1based = if ref_len == 0 { start_1based } else { start_1based + ref_len - 1 };
    Ok(Some((start_1based, end_1based)))
}

/// Returns `true` if `record` overlaps `(ref_id, interval)`.
///
/// Mirrors noodles' `query::intersects` but operates on raw BAM bytes so
/// no decoding to `Record` is needed.
fn intersects_region(
    record: &RawRecord,
    ref_id: usize,
    interval: Interval,
) -> anyhow::Result<bool> {
    use crate::cigar::reference_length_from_raw_bam_checked;
    use crate::fields::{pos, ref_id as raw_ref_id};
    use noodles::core::Position;

    let bytes: &[u8] = record;
    // Validate the record's structure (fixed core + read-name + CIGAR all within
    // bounds) before trusting any field; a truncated/malformed record fails closed
    // instead of being mistaken for a point alignment at its start.
    let ref_len_raw = reference_length_from_raw_bam_checked(bytes).context(
        "malformed BAM record: truncated, or declared CIGAR/read-name length \
             exceeds record buffer",
    )?;

    // Reference-sequence ID check.
    let record_ref_id = raw_ref_id(bytes);
    if record_ref_id < 0 {
        return Ok(false); // unmapped
    }
    if usize::try_from(record_ref_id).ok() != Some(ref_id) {
        return Ok(false);
    }

    // Position: 0-based in BAM, convert to 1-based for noodles Interval.
    let p = pos(bytes);
    if p < 0 {
        return Ok(false);
    }
    #[allow(clippy::cast_sign_loss)] // guarded by `p >= 0` check above
    let start_1based = p as usize + 1;

    // Alignment end: start_1based + ref_length - 1
    #[allow(clippy::cast_sign_loss)] // reference_length_from_raw_bam_checked returns non-negative
    let ref_len = ref_len_raw as usize;
    if ref_len == 0 {
        // Unmapped or no CIGAR; treat as single-position
        let start_pos =
            Position::try_from(start_1based).context("converting alignment start to Position")?;
        let record_interval = Interval::from(start_pos..=start_pos);
        return Ok(interval.intersects(record_interval));
    }
    let end_1based = start_1based + ref_len - 1;

    let start_pos =
        Position::try_from(start_1based).context("converting alignment start to Position")?;
    let end_pos = Position::try_from(end_1based).context("converting alignment end to Position")?;

    let record_interval = Interval::from(start_pos..=end_pos);
    Ok(interval.intersects(record_interval))
}

// ─── tests ───────────────────────────────────────────────────────────────────

#[cfg(test)]
mod tests {
    use super::*;
    use crate::testutil::*;
    use noodles::bam::{self, bai};
    use noodles::sam::alignment::io::Write as _;
    use rstest::rstest;
    use std::io::Cursor;

    // ── helpers ──────────────────────────────────────────────────────────────

    /// Build a coordinate-sorted SAM header with `ref_names` as reference sequences.
    ///
    /// The header includes `SO:coordinate` so that `bam::fs::index` will accept
    /// the BAM file without error.
    fn make_header(ref_names: &[(&str, usize)]) -> sam::Header {
        use bstr::BString;
        use noodles::sam::header::record::value::{Map, map::ReferenceSequence};
        use std::num::NonZeroUsize;

        // Parse a header line that includes SO:coordinate so that bam::fs::index
        // does not reject the file with "invalid sort order".
        let hd_line: sam::Header = "@HD\tVN:1.6\tSO:coordinate\n".parse().expect("parse @HD line");

        let mut builder = sam::Header::builder();
        if let Some(hd) = hd_line.header() {
            builder = builder.set_header(hd.clone());
        }
        for (name, len) in ref_names {
            builder = builder.add_reference_sequence(
                BString::from(*name),
                Map::<ReferenceSequence>::new(
                    NonZeroUsize::try_from(*len).expect("non-zero ref len"),
                ),
            );
        }
        builder.build()
    }

    /// Encode a `RecordBuf` to a BGZF-compressed BAM `Vec<u8>`, including
    /// the BAM magic + header preamble.
    fn encode_bam(header: &sam::Header, records: &[sam::alignment::RecordBuf]) -> Vec<u8> {
        let mut writer = bam::io::Writer::new(Vec::new());
        writer.write_header(header).expect("write header");
        for rec in records {
            writer.write_alignment_record(header, rec).expect("write record");
        }
        writer.into_inner().finish().expect("bgzf finish")
    }

    /// Build and index a BAM from `records` using noodles' `bam::fs::index`.
    /// Returns `(bam_bytes, bai_index)`.
    fn build_and_index(
        header: &sam::Header,
        records: &[sam::alignment::RecordBuf],
    ) -> (Vec<u8>, bai::Index) {
        use std::io::Write;
        use tempfile::NamedTempFile;

        let bam_bytes = encode_bam(header, records);

        // Write to a temp file so noodles can index it.
        let mut tmp = NamedTempFile::new().expect("tmp file");
        tmp.write_all(&bam_bytes).expect("write bam");
        tmp.flush().expect("flush");

        let index = bam::fs::index(tmp.path()).expect("index bam");
        (bam_bytes, index)
    }

    /// Create a simple mapped `RecordBuf` at `(ref_id, pos_1based)` with
    /// `read_len` M bases.  The `RecordBuf` is decoded from a `SamBuilder` raw
    /// record, which is sufficient for writing into a BAM file.
    #[allow(clippy::cast_possible_truncation, clippy::cast_possible_wrap)]
    fn mapped_record(
        ref_id: usize,
        pos_1based: usize,
        read_len: usize,
    ) -> sam::alignment::RecordBuf {
        use crate::SamBuilder;

        let mut b = SamBuilder::new();
        b.ref_id(ref_id as i32)
            .pos((pos_1based as i32) - 1)
            .mapq(60)
            .cigar_ops(&[encode_op(0, read_len)])
            .sequence(&vec![b'A'; read_len])
            .qualities(&vec![30u8; read_len])
            .flags(crate::flags::PAIRED | crate::flags::FIRST_SEGMENT);
        crate::noodles_compat::raw_records_to_record_bufs(&[b.build().as_ref().to_vec()])
            .expect("decode mapped record")
            .pop()
            .expect("one record")
    }

    /// Create an unmapped `RecordBuf` (`ref_id` = -1, flags UNMAPPED).
    fn unmapped_record() -> sam::alignment::RecordBuf {
        use crate::SamBuilder;

        let mut b = SamBuilder::new();
        b.ref_id(-1)
            .pos(-1)
            .mapq(0)
            .sequence(b"ACGT")
            .qualities(&[30u8; 4])
            .flags(crate::flags::UNMAPPED);
        crate::noodles_compat::raw_records_to_record_bufs(&[b.build().as_ref().to_vec()])
            .expect("decode unmapped record")
            .pop()
            .expect("one record")
    }

    // ── tests ─────────────────────────────────────────────────────────────────

    /// A well-formed mapped record yields its 1-based inclusive reference span.
    #[test]
    fn record_span_1based_accepts_well_formed_record() {
        // pos=100 (0-based) -> 101 (1-based), 10M -> [101, 110].
        let bytes = make_bam_bytes(0, 100, 0, b"rea", &[(10 << 4)], 10, -1, -1, &[]);
        let record = RawRecord::from(bytes);
        assert_eq!(record_span_1based(&record).expect("span"), Some((101, 110)));
    }

    /// A record whose declared CIGAR overruns its buffer must fail closed rather
    /// than be mistaken for a point alignment at its start position.
    #[test]
    fn record_span_1based_rejects_truncated_cigar() {
        // Well-formed 10M record, then corrupt n_cigar_op to claim 100 ops.
        let mut bytes = make_bam_bytes(0, 100, 0, b"rea", &[(10 << 4)], 10, -1, -1, &[]);
        bytes[12..14].copy_from_slice(&100u16.to_le_bytes());
        let record = RawRecord::from(bytes);
        assert!(record_span_1based(&record).is_err());
    }

    /// A record shorter than the 32-byte fixed core must fail closed rather than
    /// let a fixed-offset field read index past the buffer and panic. This is the
    /// guard `RawMultiQueryIter::next` runs *before* reading `ref_id`, so a corrupt
    /// BGZF block declaring a tiny nonzero size fails closed instead of crashing
    /// the scan (a bare `fields::ref_id` would slice `bam[0..4]` and panic here).
    #[test]
    fn record_span_1based_rejects_sub_core_length_record() {
        // Three bytes: shorter than the fixed core, and shorter than the 4 bytes a
        // bare `ref_id` read would index.
        let record = RawRecord::from(vec![0u8; 3]);
        assert!(record_span_1based(&record).is_err());
    }

    /// The single-interval overlap test shares the same malformed-CIGAR guard and
    /// also fails closed on a truncated record.
    #[test]
    fn intersects_region_rejects_truncated_cigar() {
        use noodles::core::Position;
        let mut bytes = make_bam_bytes(0, 100, 0, b"rea", &[(10 << 4)], 10, -1, -1, &[]);
        bytes[12..14].copy_from_slice(&100u16.to_le_bytes());
        let record = RawRecord::from(bytes);
        let interval =
            Interval::from(Position::try_from(1).unwrap()..=Position::try_from(1000).unwrap());
        assert!(intersects_region(&record, 0, interval).is_err());
    }

    /// Query a region; verify the right records are returned.
    #[test]
    fn query_returns_overlapping_records() {
        let header = make_header(&[("chr1", 1000)]);

        // Three records: pos 10-19, pos 50-59, pos 100-109 (all 10bp M).
        let r1 = mapped_record(0, 10, 10);
        let r2 = mapped_record(0, 50, 10);
        let r3 = mapped_record(0, 100, 10);

        let records = [r1, r2, r3];
        let (bam_bytes, index) = build_and_index(&header, &records);

        let mut reader = IndexedRawBamReader::new(Cursor::new(bam_bytes), index);
        let hdr = reader.read_header().expect("header");

        // Query chr1:45-55 — should hit r2 (pos 50-59) only.
        let region = "chr1:45-55".parse().expect("region");
        let results: Vec<_> = reader.query(&hdr, &region).expect("query").collect();
        assert_eq!(results.len(), 1, "expected exactly 1 overlapping record");
        let raw = results[0].as_ref().expect("record ok");
        // Verify position (0-based pos = 49 → 1-based start = 50).
        assert_eq!(raw.pos(), 49);
    }

    /// Ported from htsjdk `QueryIntervalTest.testOptimizeIntervals`: overlapping
    /// and abutting intervals merge; a one-base gap keeps them separate; unsorted
    /// input is sorted. Compared as `(start, end)` pairs so the assertion is
    /// independent of `Interval`'s own equality. Each case names the scenario so a
    /// failure reports the exact interval set and expected merged output.
    #[rstest]
    #[case::overlapping(vec![(1520, 1521), (1521, 1525)], vec![(1520, 1525)])]
    #[case::abutting_zero_gap(vec![(1520, 1521), (1522, 1525)], vec![(1520, 1525)])]
    #[case::one_base_gap_kept_separate(
        vec![(1520, 1521), (1523, 1525)],
        vec![(1520, 1521), (1523, 1525)]
    )]
    #[case::unsorted_overlapping(vec![(20, 30), (10, 20)], vec![(10, 30)])]
    fn merge_intervals_matches_htsjdk_optimize_intervals(
        #[case] input: Vec<(usize, usize)>,
        #[case] expected: Vec<(usize, usize)>,
    ) {
        use noodles::core::Position;
        let iv = |s: usize, e: usize| {
            Interval::from(Position::try_from(s).unwrap()..=Position::try_from(e).unwrap())
        };
        let pairs = |ivs: Vec<Interval>| -> Vec<(usize, usize)> {
            ivs.iter().map(|i| (interval_start(i), interval_end(i))).collect()
        };

        let intervals: Vec<Interval> = input.iter().map(|&(s, e)| iv(s, e)).collect();
        assert_eq!(pairs(merge_intervals(intervals)), expected);
    }

    /// Ported from htsjdk `BAMFileIndexTest.testMultiIntervalQuery`: a
    /// multi-interval query returns exactly the union of the per-interval
    /// queries, with each record yielded once and in coordinate order. Records
    /// are keyed by `(ref_id, pos)`, which is unique here by construction.
    #[test]
    fn query_intervals_equals_union_of_single_queries() {
        use noodles::core::Region;
        use std::collections::BTreeSet;

        let header = make_header(&[("chr1", 1000), ("chr2", 1000)]);
        // SPAN starts before every region and covers them all — a naive
        // per-region loop would surface it once per region; the merged scan must
        // surface it exactly once.
        let records = [
            mapped_record(0, 5, 300),  // chr1:5-304 (spans every chr1 region)
            mapped_record(0, 50, 10),  // chr1:50-59
            mapped_record(0, 100, 10), // chr1:100-109 (inside two overlapping regions)
            mapped_record(0, 150, 10), // chr1:150-159
            mapped_record(0, 400, 10), // chr1:400-409 (not queried)
            mapped_record(1, 30, 10),  // chr2:30-39
        ];
        let (bam, index) = build_and_index(&header, &records);

        let regions: Vec<Region> = ["chr1:45-105", "chr1:100-160", "chr2:25-35"]
            .iter()
            .map(|s| s.parse().expect("region"))
            .collect();

        let mut reader = IndexedRawBamReader::new(Cursor::new(bam), index);
        let hdr = reader.read_header().expect("header");

        // Union of the per-interval single-region queries.
        let mut union: BTreeSet<(i32, i32)> = BTreeSet::new();
        for region in &regions {
            for result in reader.query(&hdr, region).expect("query") {
                let rec = result.expect("record ok");
                union.insert((rec.ref_id(), rec.pos()));
            }
        }

        // Multi-interval query over the same reader.
        let multi: Vec<_> = reader
            .query_intervals(&hdr, &regions)
            .expect("query_intervals")
            .map(|r| r.expect("record ok"))
            .collect();
        let multi_keys: Vec<(i32, i32)> = multi.iter().map(|r| (r.ref_id(), r.pos())).collect();

        // Each record appears exactly once.
        let multi_set: BTreeSet<(i32, i32)> = multi_keys.iter().copied().collect();
        assert_eq!(
            multi_set.len(),
            multi_keys.len(),
            "multi-interval query yielded a record more than once: {multi_keys:?}"
        );
        // Same set as the union of the single-interval queries.
        assert_eq!(multi_set, union, "multi-interval set differs from single-query union");
        // Emitted in coordinate order.
        let mut sorted = multi_keys.clone();
        sorted.sort_unstable();
        assert_eq!(multi_keys, sorted, "records not in coordinate order: {multi_keys:?}");
    }

    /// `make_interval` must keep an unbounded end (`usize::MAX` sentinel) open
    /// rather than pinning it to a finite `Position`, while a bounded end still
    /// round-trips inclusively.
    #[test]
    fn make_interval_preserves_unbounded_end() {
        let unbounded = make_interval(100, usize::MAX);
        assert_eq!(interval_start(&unbounded), 100);
        assert!(unbounded.end().is_none(), "unbounded end must not be pinned to a finite Position");
        assert_eq!(interval_end(&unbounded), usize::MAX);

        let bounded = make_interval(100, 200);
        assert_eq!((interval_start(&bounded), interval_end(&bounded)), (100, 200));
        assert_eq!(bounded.end().map(usize::from), Some(200));
    }

    /// A whole-contig (`chr1`) query carries an unbounded interval; every record on
    /// the contig — including ones deep in its tail — must be returned, proving the
    /// unbounded end is not pinned to a finite position.
    #[test]
    fn query_intervals_whole_contig_covers_tail() {
        use noodles::core::Region;
        use std::collections::BTreeSet;

        let header = make_header(&[("chr1", 100_000)]);
        let records = [
            mapped_record(0, 10, 10),     // chr1:10-19
            mapped_record(0, 100, 10),    // chr1:100-109
            mapped_record(0, 50_000, 10), // chr1:50000-50009 (deep in the tail)
        ];
        let (bam, index) = build_and_index(&header, &records);
        let mut reader = IndexedRawBamReader::new(Cursor::new(bam), index);
        let hdr = reader.read_header().expect("header");

        let regions: Vec<Region> = vec!["chr1".parse().expect("region")];
        let got: BTreeSet<i32> = reader
            .query_intervals(&hdr, &regions)
            .expect("query_intervals")
            .map(|r| r.expect("record ok").pos())
            .collect();
        // 0-based positions of all three records (10,100,50000 → -1).
        let expected: BTreeSet<i32> = [9, 99, 49_999].into_iter().collect();
        assert_eq!(got, expected, "whole-contig query must cover the entire reference");
    }

    /// Overlap boundaries for the merged-interval binary search: a record whose
    /// span touches an interval at either edge overlaps, while one that falls in
    /// the gap between two intervals — or entirely before/after them — does not.
    /// Guards the `partition_point` bounds (`end >= rec_start` / `start <= rec_end`)
    /// against off-by-one regressions.
    #[test]
    fn query_intervals_overlap_boundaries() {
        use noodles::core::Region;
        use std::collections::BTreeSet;

        let header = make_header(&[("chr1", 1000)]);
        // Two point regions at chr1:100 and chr1:200 → merged intervals
        // [(100,100), (200,200)]. Records chosen to sit on each boundary and in
        // the gaps; only those whose 1-based span includes 100 or 200 overlap.
        let records = [
            mapped_record(0, 50, 10),  // [50,59]    before first  → excluded
            mapped_record(0, 91, 10),  // [91,100]   ends at 100    → included
            mapped_record(0, 100, 1),  // [100,100]  exact point    → included
            mapped_record(0, 101, 10), // [101,110]  gap            → excluded
            mapped_record(0, 196, 5),  // [196,200]  ends at 200    → included
            mapped_record(0, 200, 5),  // [200,204]  starts at 200  → included
            mapped_record(0, 205, 5),  // [205,209]  after last     → excluded
        ];
        let (bam, index) = build_and_index(&header, &records);

        let regions: Vec<Region> =
            ["chr1:100-100", "chr1:200-200"].iter().map(|s| s.parse().expect("region")).collect();

        let mut reader = IndexedRawBamReader::new(Cursor::new(bam), index);
        let hdr = reader.read_header().expect("header");

        let got: BTreeSet<i32> = reader
            .query_intervals(&hdr, &regions)
            .expect("query_intervals")
            .map(|r| r.expect("record ok").pos())
            .collect();

        // 0-based positions of the four overlapping records (91,100,196,200 → -1).
        let expected: BTreeSet<i32> = [90, 99, 195, 199].into_iter().collect();
        assert_eq!(got, expected, "boundary overlap set mismatch");
    }

    /// `query_intervals` with no regions yields nothing.
    #[test]
    fn query_intervals_no_regions_is_empty() {
        let header = make_header(&[("chr1", 1000)]);
        let (bam, index) = build_and_index(&header, &[mapped_record(0, 10, 10)]);
        let mut reader = IndexedRawBamReader::new(Cursor::new(bam), index);
        let hdr = reader.read_header().expect("header");
        let count = reader.query_intervals(&hdr, &[]).expect("query_intervals").count();
        assert_eq!(count, 0, "no regions must yield no records");
    }

    /// Query a region with no overlapping records — iterator must be empty.
    #[test]
    fn query_empty_region_returns_no_records() {
        let header = make_header(&[("chr1", 1000)]);
        let r1 = mapped_record(0, 10, 10);
        let (bam_bytes, index) = build_and_index(&header, &[r1]);

        let mut reader = IndexedRawBamReader::new(Cursor::new(bam_bytes), index);
        let hdr = reader.read_header().expect("header");

        // Query chr1:500-600 — no records there.
        let region = "chr1:500-600".parse().expect("region");
        let results: Vec<_> = reader.query(&hdr, &region).expect("query").collect();
        assert_eq!(results.len(), 0, "expected no records in empty region");
    }

    /// Multiple records in the same region; verify all are returned in order.
    #[test]
    fn query_returns_multiple_overlapping_records() {
        let header = make_header(&[("chr1", 1000)]);
        let r1 = mapped_record(0, 100, 50);
        let r2 = mapped_record(0, 120, 50);
        let r3 = mapped_record(0, 200, 50); // outside query window

        let records = [r1, r2, r3];
        let (bam_bytes, index) = build_and_index(&header, &records);

        let mut reader = IndexedRawBamReader::new(Cursor::new(bam_bytes), index);
        let hdr = reader.read_header().expect("header");

        // Query chr1:100-169 — hits r1 and r2, not r3.
        let region = "chr1:100-169".parse().expect("region");
        let results: Vec<_> = reader.query(&hdr, &region).expect("query").collect();
        assert_eq!(results.len(), 2, "expected 2 overlapping records");

        let pos0 = results[0].as_ref().expect("ok").pos();
        let pos1 = results[1].as_ref().expect("ok").pos();
        assert_eq!(pos0, 99); // 1-based 100 → 0-based 99
        assert_eq!(pos1, 119); // 1-based 120 → 0-based 119
    }

    /// Two queries on the same reader; verify seeking works correctly.
    #[test]
    fn multiple_queries_on_same_reader() {
        let header = make_header(&[("chr1", 1000)]);
        let r1 = mapped_record(0, 10, 10);
        let r2 = mapped_record(0, 200, 10);

        let records = [r1, r2];
        let (bam_bytes, index) = build_and_index(&header, &records);

        let mut reader = IndexedRawBamReader::new(Cursor::new(bam_bytes), index);
        let hdr = reader.read_header().expect("header");

        // First query: chr1:1-20 → r1 only.
        let region1 = "chr1:1-20".parse().expect("region");
        let res1: Vec<_> = reader.query(&hdr, &region1).expect("query").collect();
        assert_eq!(res1.len(), 1);
        assert_eq!(res1[0].as_ref().expect("ok").pos(), 9);

        // Second query: chr1:190-210 → r2 only.
        let region2 = "chr1:190-210".parse().expect("region");
        let res2: Vec<_> = reader.query(&hdr, &region2).expect("query").collect();
        assert_eq!(res2.len(), 1);
        assert_eq!(res2[0].as_ref().expect("ok").pos(), 199);
    }

    /// `query_unmapped` returns only unmapped records.
    #[test]
    fn query_unmapped_returns_unmapped_records() {
        let header = make_header(&[("chr1", 1000)]);

        // One mapped record and one unmapped record.
        let r_mapped = mapped_record(0, 10, 10);
        let r_unmap = unmapped_record();

        let records = [r_mapped, r_unmap];
        let (bam_bytes, index) = build_and_index(&header, &records);

        let mut reader = IndexedRawBamReader::new(Cursor::new(bam_bytes), index);
        let _hdr = reader.read_header().expect("header");

        let results: Vec<_> = reader.query_unmapped().expect("query_unmapped").collect();
        assert_eq!(results.len(), 1, "expected exactly 1 unmapped record");
        assert!(results[0].as_ref().expect("ok").is_unmapped());
    }

    /// If neither the index nor `read_header` provides a first-record
    /// position, `query_unmapped` must return an explicit error rather than
    /// silently fall back to BGZF offset 0 (which would parse the BAM header
    /// as a record).
    #[test]
    fn query_unmapped_errors_when_no_marker_and_header_unread() {
        let header = make_header(&[("chr1", 1000)]);

        // Empty BAM → noodles does not set last_first_record_start_position.
        let records: [sam::alignment::RecordBuf; 0] = [];
        let (bam_bytes, index) = build_and_index(&header, &records);
        assert!(index.last_first_record_start_position().is_none());

        let mut reader = IndexedRawBamReader::new(Cursor::new(bam_bytes), index);
        // Intentionally NOT calling read_header — the error path depends on
        // first_record_position being None.
        let Err(err) = reader.query_unmapped() else {
            panic!("expected query_unmapped to error without a position");
        };
        let msg = format!("{err}");
        assert!(
            msg.contains("read_header"),
            "error message should mention read_header(); got: {msg}",
        );
    }

    /// When the index lacks a `last_first_record_start_position` marker,
    /// `query_unmapped` must rewind to the first record (captured in
    /// `read_header`) rather than resume wherever the BGZF cursor was last
    /// left. Exercises the empty-BAM case where noodles leaves the marker
    /// unset; verifies the fallback seek succeeds and two consecutive calls
    /// yield the same (empty) result.
    #[test]
    fn query_unmapped_falls_back_to_first_record_when_no_marker() {
        let header = make_header(&[("chr1", 1000)]);

        // Empty BAM → noodles does not set last_first_record_start_position.
        let records: [sam::alignment::RecordBuf; 0] = [];
        let (bam_bytes, index) = build_and_index(&header, &records);
        assert!(
            index.last_first_record_start_position().is_none(),
            "test precondition: no last_first_record_start_position",
        );

        let mut reader = IndexedRawBamReader::new(Cursor::new(bam_bytes), index);
        let _hdr = reader.read_header().expect("header");

        // First call consumes the (empty) unmapped section and advances the
        // BGZF cursor. Second call without the fix would silently return an
        // iterator rooted at that advanced cursor; with the fix, it rewinds
        // to the captured first-record position.
        let a: Vec<_> = reader.query_unmapped().expect("1st query_unmapped").collect();
        assert_eq!(a.len(), 0);
        let b: Vec<_> = reader.query_unmapped().expect("2nd query_unmapped").collect();
        assert_eq!(b.len(), 0);
    }
}
