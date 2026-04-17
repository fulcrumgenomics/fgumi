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
        Self { inner: noodles::bam::io::Reader::new(reader), index: Box::new(index) }
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
        self.inner.read_header().context("reading BAM header")
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
    /// Returns an error if seeking or reading fails.
    pub fn query_unmapped(&mut self) -> anyhow::Result<UnmappedIter<'_, R>> {
        if let Some(pos) = self.index.last_first_record_start_position() {
            self.inner
                .get_mut()
                .seek_to_virtual_position(pos)
                .context("seeking to unmapped position")?;
        }
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

// ─── overlap test ───────────────────────────────────────────────────────────

/// Returns `true` if `record` overlaps `(ref_id, interval)`.
///
/// Mirrors noodles' `query::intersects` but operates on raw BAM bytes so
/// no decoding to `Record` is needed.
fn intersects_region(
    record: &RawRecord,
    ref_id: usize,
    interval: Interval,
) -> anyhow::Result<bool> {
    use crate::cigar::reference_length_from_raw_bam;
    use crate::fields::{pos, ref_id as raw_ref_id};
    use noodles::core::Position;

    let bytes: &[u8] = record;
    if bytes.len() < 8 {
        return Ok(false);
    }

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
    let ref_len_raw = reference_length_from_raw_bam(bytes);
    #[allow(clippy::cast_sign_loss)] // reference_length_from_raw_bam returns non-negative
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
}
