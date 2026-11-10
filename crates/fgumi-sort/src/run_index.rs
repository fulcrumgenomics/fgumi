//! In-memory index of a spill run.
//!
//! A zstd spill run is a sequence of independently decodable frames, each
//! holding a whole number of `[key?][u32 len][body]` records. [`RunIndex`]
//! records, for every entry (one frame, or `stride` consecutive frames once
//! sparsified), its compressed and raw sizes, its record count, and the
//! serialized sort key of its first record. It is a data type, the helpers
//! that derive one frame's metadata, and the run-file footer format: its
//! encoding ([`RunIndex::write_footer`]), which the spill writers append after
//! a run's frames, and its location and decoding (`locate_footer`,
//! `read_footer`), positional reads through a
//! [`PositionalSource`] and the only reads of a run file here.
//!
//! # Footer
//!
//! An indexed zstd run carries its index at the end of its own file, so no
//! closed run's index is resident while phase 1 fills the sort buffer:
//!
//! ```text
//! index block : b"FGIX" | u32 version | u16 key_width | u16 pad | u32 stride |
//!               u64 frames | u32 max_frame_raw | u64 n | offsets (n+1) x u64 |
//!               raw_len n x u32 | records n x u32 | keys n x key_width
//! file        : [magic][frame]...[frame][index block][u64 block length][b"FGUMIIDX"]
//! ```
//!
//! Little-endian throughout. Whether a run has a footer is its writer's
//! word: the opener passes where the writer ended the frames, and the bytes
//! after them must be exactly one footer, whose 16-byte trailer names the
//! block's length. The block's magic, its exact length, its entry count
//! against its frames and stride, and its offsets tiling the frame bytes
//! validate it.

use std::io::{self, Cursor, Write};

use fgumi_bam_io::pread::{PositionalSource, read_at_exact};

use crate::keys::{RawCoordinateKey, RawSortKey};
use crate::spill_block::{RECORD_LEN_BYTES, SpillKeyKind};

/// Serialized key width, in bytes, that a [`RunIndex`] stores per entry for
/// `kind`, or `None` for kinds that are not indexed (the queryname orders).
///
/// A template kind's width is the key its spill records carry
/// ([`SpillKeyKind::serialized_key_bytes`]); a coordinate run's records carry
/// no key, so its entries store the [`RawCoordinateKey`] derived from the
/// first record's body. The `const` assertions below keep the table in step
/// with both.
#[must_use]
pub const fn spill_key_width(kind: SpillKeyKind) -> Option<u16> {
    match kind {
        SpillKeyKind::Coordinate => Some(8),
        SpillKeyKind::TemplateK24 => Some(24),
        SpillKeyKind::TemplateCb32 | SpillKeyKind::TemplateTert32 => Some(32),
        SpillKeyKind::TemplateK40 => Some(40),
        SpillKeyKind::QuerynameLex | SpillKeyKind::QuerynameNatural => None,
    }
}

/// Whether `kind`'s entry width is `width` bytes.
const fn width_is(kind: SpillKeyKind, width: usize) -> bool {
    match spill_key_width(kind) {
        Some(w) => w as usize == width,
        None => false,
    }
}

const _: () = {
    assert!(width_is(SpillKeyKind::Coordinate, RawCoordinateKey::SIZE));
    assert!(width_is(SpillKeyKind::TemplateK24, SpillKeyKind::TemplateK24.serialized_key_bytes()));
    assert!(width_is(
        SpillKeyKind::TemplateCb32,
        SpillKeyKind::TemplateCb32.serialized_key_bytes()
    ));
    assert!(width_is(
        SpillKeyKind::TemplateTert32,
        SpillKeyKind::TemplateTert32.serialized_key_bytes()
    ));
    assert!(width_is(SpillKeyKind::TemplateK40, SpillKeyKind::TemplateK40.serialized_key_bytes()));
};

/// Per-entry index of one spill run (struct of arrays; one entry per frame, or
/// per `stride` consecutive frames once [`RunIndex::sparsify`]d).
///
/// `offsets` holds `len() + 1` absolute file offsets: entry `i` spans
/// `[offsets[i], offsets[i + 1])`, so the last element is the end of the last
/// entry. A freshly built index starts at offset 0; the spill writer
/// [`rebase`](RunIndex::rebase)s it by the file position its frames start at.
#[derive(Clone, Debug, PartialEq, Eq)]
pub struct RunIndex {
    key_width: u16,
    /// Frames each entry covers (1 until sparsified).
    stride: u32,
    /// Frames covered by all entries together (the last entry of a sparsified
    /// index may cover fewer than `stride`).
    frames: u64,
    /// The largest raw length of any one frame the index covers. Kept when the
    /// index is sparsified, so it bounds one frame's decode buffer at any
    /// stride.
    max_frame_raw: u32,
    offsets: Vec<u64>,
    raw_len: Vec<u32>,
    records: Vec<u32>,
    /// First-record keys, `key_width` bytes per entry, back to back.
    keys: Vec<u8>,
}

impl Default for RunIndex {
    fn default() -> Self {
        Self::new(0)
    }
}

impl RunIndex {
    /// An empty index whose frame keys are `key_width` bytes wide.
    #[must_use]
    pub fn new(key_width: u16) -> Self {
        Self {
            key_width,
            stride: 1,
            frames: 0,
            max_frame_raw: 0,
            offsets: vec![0],
            raw_len: Vec::new(),
            records: Vec::new(),
            keys: Vec::new(),
        }
    }

    /// An empty index whose entries will each cover `stride` frames, built
    /// with [`Self::append_folded`] so no denser index exists on the way.
    ///
    /// # Panics
    ///
    /// If `stride` is 0.
    #[must_use]
    pub fn new_strided(key_width: u16, stride: u32) -> Self {
        assert!(stride > 0, "a run index's stride is at least 1");
        Self { stride, ..Self::new(key_width) }
    }

    /// Heap bytes one entry costs (its `u64` offset, `u32` raw length and
    /// record count, and its key) for keys `key_width` bytes wide.
    #[must_use]
    pub const fn entry_bytes(key_width: u16) -> u64 {
        (size_of::<u64>() + 2 * size_of::<u32>()) as u64 + key_width as u64
    }

    /// Append one frame's metadata. `key` is the first record's serialized key.
    /// `comp_len` is the frame's file span (its `u32` length prefix plus the
    /// compressed frame).
    ///
    /// # Panics
    ///
    /// If the index has been sparsified (it is append-only at stride 1).
    pub fn push_frame(&mut self, comp_len: u32, raw_len: u32, records: u32, key: &[u8]) {
        debug_assert_eq!(key.len(), usize::from(self.key_width), "frame key width mismatch");
        assert_eq!(self.stride, 1, "cannot push a frame onto a sparsified index");
        self.offsets.push(self.end_of_last() + u64::from(comp_len));
        self.raw_len.push(raw_len);
        self.records.push(records);
        self.keys.extend_from_slice(key);
        self.frames += 1;
        self.max_frame_raw = self.max_frame_raw.max(raw_len);
    }

    /// End offset of the last entry (the first entry's start when empty).
    fn end_of_last(&self) -> u64 {
        *self.offsets.last().expect("offsets always hold a trailing end")
    }

    /// Append every entry of `other`, which must have the same key width and
    /// stride, shifting its offsets by `base` (`other`'s offsets are relative
    /// to the start of its own frames). `base` should be this index's end.
    ///
    /// # Panics
    ///
    /// If the key widths or strides differ, or `other` does not continue this
    /// index (this index's end is not `base` plus `other`'s first offset).
    pub fn append_rebased(&mut self, other: &RunIndex, base: u64) {
        assert_eq!(self.key_width, other.key_width, "appending indexes of different widths");
        assert_eq!(self.stride, other.stride, "appending indexes of different strides");
        assert_eq!(
            self.end_of_last(),
            base + other.offsets[0],
            "appended entries must continue the index"
        );
        // Drop our trailing end: `other`'s first offset replaces it.
        self.offsets.pop();
        self.offsets.extend(other.offsets.iter().map(|&o| o + base));
        self.raw_len.extend_from_slice(&other.raw_len);
        self.records.extend_from_slice(&other.records);
        self.keys.extend_from_slice(&other.keys);
        self.frames += other.frames;
        self.max_frame_raw = self.max_frame_raw.max(other.max_frame_raw);
    }

    /// Append the frames of `frames`, a dense (stride 1) index of the bytes
    /// that continue this index's file, shifting their offsets by `base`,
    /// folded at this index's own stride: they fill this index's short last
    /// entry first, then start a new entry every `stride` frames. The result
    /// is the index [`Self::sparsify`] would make of all the frames dense,
    /// built without the dense intermediate.
    ///
    /// # Errors
    ///
    /// [`io::ErrorKind::InvalidData`], leaving the index unchanged, when an
    /// entry's raw length or record count would not fit a `u32` (as
    /// [`Self::sparsify`]).
    ///
    /// # Panics
    ///
    /// If the key widths differ, `frames` is sparsified, or it does not
    /// continue this index (this index's end is not `base` plus its first
    /// offset).
    pub fn append_folded(&mut self, frames: &RunIndex, base: u64) -> io::Result<()> {
        assert_eq!(self.key_width, frames.key_width, "appending indexes of different widths");
        assert_eq!(frames.stride, 1, "only dense frames are folded");
        assert_eq!(
            self.end_of_last(),
            base + frames.offsets[0],
            "appended entries must continue the index"
        );
        let stride = u64::from(self.stride);
        let overflow = |what: &str| {
            io::Error::new(
                io::ErrorKind::InvalidData,
                format!("folding a run index at stride {stride}: {what} overflows u32"),
            )
        };
        // Check every entry's sums first, so an overflow changes nothing.
        let (mut raw, mut records) = if self.frames.is_multiple_of(stride) {
            (0, 0)
        } else {
            (
                u64::from(*self.raw_len.last().expect("a short entry")),
                u64::from(*self.records.last().expect("a short entry")),
            )
        };
        for (i, (&r, &n)) in frames.raw_len.iter().zip(&frames.records).enumerate() {
            if (self.frames + i as u64).is_multiple_of(stride) {
                (raw, records) = (0, 0);
            }
            raw += u64::from(r);
            records += u64::from(n);
            if u32::try_from(raw).is_err() {
                return Err(overflow("an entry's raw length"));
            }
            if u32::try_from(records).is_err() {
                return Err(overflow("an entry's record count"));
            }
        }
        // Drop our trailing end: each new entry's start, then `frames`' end,
        // replace it.
        self.offsets.pop();
        for i in 0..frames.len() {
            if self.frames.is_multiple_of(stride) {
                self.offsets.push(frames.offsets[i] + base);
                self.raw_len.push(frames.raw_len[i]);
                self.records.push(frames.records[i]);
                self.keys.extend_from_slice(frames.key_bytes(i));
            } else {
                // Checked above: every entry's sums fit.
                *self.raw_len.last_mut().expect("a short entry") += frames.raw_len[i];
                *self.records.last_mut().expect("a short entry") += frames.records[i];
            }
            self.frames += 1;
        }
        self.offsets.push(frames.end_of_last() + base);
        self.max_frame_raw = self.max_frame_raw.max(frames.max_frame_raw);
        Ok(())
    }

    /// Append every entry of `other`, which continues this index's file
    /// (`other`'s offsets are relative to this index's end).
    #[cfg(test)]
    pub(crate) fn append(&mut self, other: &RunIndex) {
        let base = self.end_of_last();
        self.append_rebased(other, base);
    }

    /// Shift every offset by `base`.
    pub fn rebase(&mut self, base: u64) {
        for o in &mut self.offsets {
            *o += base;
        }
    }

    /// Number of entries.
    #[must_use]
    pub fn len(&self) -> usize {
        self.raw_len.len()
    }

    /// Whether the index has no entries.
    #[must_use]
    pub fn is_empty(&self) -> bool {
        self.raw_len.is_empty()
    }

    /// Width in bytes of each entry's serialized key.
    #[must_use]
    pub fn key_width(&self) -> u16 {
        self.key_width
    }

    /// Frames each entry covers.
    #[must_use]
    pub fn stride(&self) -> u32 {
        self.stride
    }

    /// Frames covered by the whole index: exact even after [`Self::sparsify`]
    /// leaves a short last entry.
    #[must_use]
    pub fn covered_frames(&self) -> u64 {
        self.frames
    }

    /// The largest raw length of any one frame the index covers (0 when
    /// empty), at any stride.
    #[must_use]
    pub fn max_frame_raw(&self) -> u32 {
        self.max_frame_raw
    }

    /// The largest raw length of any one entry (0 when empty).
    #[must_use]
    pub fn max_entry_raw(&self) -> u32 {
        self.raw_len.iter().copied().max().unwrap_or(0)
    }

    /// File offset where entry `i` starts.
    #[must_use]
    pub fn offset(&self, i: usize) -> u64 {
        self.offsets[i]
    }

    /// File offset where entry `i` ends.
    #[must_use]
    pub fn end(&self, i: usize) -> u64 {
        self.offsets[i + 1]
    }

    /// Bytes entry `i` spans in the run file: for each frame it covers, the
    /// `u32` length prefix plus the compressed frame.
    #[must_use]
    pub fn comp_len(&self, i: usize) -> u64 {
        self.offsets[i + 1] - self.offsets[i]
    }

    /// Decompressed length of entry `i`.
    #[must_use]
    pub fn raw_len(&self, i: usize) -> u32 {
        self.raw_len[i]
    }

    /// Record count of entry `i`.
    #[must_use]
    pub fn records(&self, i: usize) -> u32 {
        self.records[i]
    }

    /// Serialized first-record key of entry `i`.
    #[must_use]
    pub fn key_bytes(&self, i: usize) -> &[u8] {
        let w = usize::from(self.key_width);
        &self.keys[i * w..(i + 1) * w]
    }

    /// Decode entry `i`'s first-record key as `K`.
    ///
    /// # Errors
    ///
    /// Returns an error if `K`'s serialized size is not this index's key width,
    /// or the key bytes do not decode.
    pub fn key_at<K: RawSortKey>(&self, i: usize) -> io::Result<K> {
        if K::SERIALIZED_SIZE != Some(usize::from(self.key_width)) {
            return Err(io::Error::new(
                io::ErrorKind::InvalidInput,
                format!(
                    "run index key width {} does not match requested key type ({:?})",
                    self.key_width,
                    K::SERIALIZED_SIZE
                ),
            ));
        }
        K::read_from(&mut Cursor::new(self.key_bytes(i)))
    }

    /// Total decompressed bytes across all entries.
    #[must_use]
    pub fn total_raw_bytes(&self) -> u64 {
        self.raw_len.iter().map(|&n| u64::from(n)).sum()
    }

    /// Total record count across all entries.
    #[must_use]
    pub fn total_records(&self) -> u64 {
        self.records.iter().map(|&n| u64::from(n)).sum()
    }

    /// Keep every `factor`-th entry (0, `factor`, 2 x `factor`, ...), in place.
    /// A kept entry spans its whole group: its raw length and record count are
    /// the group's sums, its key the group's first key, and the last entry
    /// keeps the run's end. `stride` is multiplied by `factor`. A `factor` of
    /// 0 or 1 changes nothing.
    ///
    /// # Errors
    ///
    /// Returns [`io::ErrorKind::InvalidData`], leaving the index unchanged,
    /// when a group's raw length or record count, or the new stride, does
    /// not fit a `u32` (a run of very long records, frames of up to 64 MiB
    /// each). A caller drops the index then, so the merge falls back to
    /// serial rather than reading a corrupt index.
    pub fn sparsify(&mut self, factor: u32) -> io::Result<()> {
        if factor <= 1 {
            return Ok(());
        }
        let f = factor as usize;
        let overflow = |what: &str| {
            io::Error::new(
                io::ErrorKind::InvalidData,
                format!("sparsifying a run index by {factor}: {what} overflows u32"),
            )
        };
        let stride = self.stride.checked_mul(factor).ok_or_else(|| overflow("the stride"))?;
        for (raw, records) in self.raw_len.chunks(f).zip(self.records.chunks(f)) {
            let raw: u64 = raw.iter().map(|&n| u64::from(n)).sum();
            let records: u64 = records.iter().map(|&n| u64::from(n)).sum();
            if u32::try_from(raw).is_err() {
                return Err(overflow("a group's raw length"));
            }
            if u32::try_from(records).is_err() {
                return Err(overflow("a group's record count"));
            }
        }
        let kw = usize::from(self.key_width);
        let n = self.len();
        let kept = n.div_ceil(f);
        for i in 0..n {
            let w = i / f;
            if i % f == 0 {
                self.offsets[w] = self.offsets[i];
                self.raw_len[w] = self.raw_len[i];
                self.records[w] = self.records[i];
                if w != i {
                    self.keys.copy_within(i * kw..(i + 1) * kw, w * kw);
                }
            } else {
                // Checked above: every group's sums fit.
                self.raw_len[w] += self.raw_len[i];
                self.records[w] += self.records[i];
            }
        }
        self.offsets[kept] = self.offsets[n];
        self.offsets.truncate(kept + 1);
        self.raw_len.truncate(kept);
        self.records.truncate(kept);
        self.keys.truncate(kept * kw);
        self.stride = stride;
        Ok(())
    }

    /// Release spare capacity, so [`Self::heap_bytes`] is exact.
    pub fn shrink_to_fit(&mut self) {
        self.offsets.shrink_to_fit();
        self.raw_len.shrink_to_fit();
        self.records.shrink_to_fit();
        self.keys.shrink_to_fit();
    }

    /// Heap bytes held by this index (capacity based).
    #[must_use]
    pub fn heap_bytes(&self) -> usize {
        self.offsets.capacity() * size_of::<u64>()
            + (self.raw_len.capacity() + self.records.capacity()) * size_of::<u32>()
            + self.keys.capacity()
    }
}

/// Magic that starts a run-index footer block.
const FOOTER_BLOCK_MAGIC: [u8; 4] = *b"FGIX";

/// Footer block format version.
const FOOTER_VERSION: u32 = 1;

/// Bytes of the footer block's fixed header (magic through `n`).
const FOOTER_HEADER_BYTES: usize = 4 + 4 + 2 + 2 + 4 + 8 + 4 + 8;

/// Magic that ends a run file carrying an index footer.
const FOOTER_TRAILER_MAGIC: [u8; 8] = *b"FGUMIIDX";

/// Bytes of the trailer after the index block: its `u64` length and
/// [`FOOTER_TRAILER_MAGIC`].
const FOOTER_TRAILER_BYTES: usize = 16;

/// The index block length a run file's last [`FOOTER_TRAILER_BYTES`] bytes
/// name, or `None` when they are not a footer trailer.
fn parse_footer_trailer(tail: &[u8; FOOTER_TRAILER_BYTES]) -> Option<u64> {
    let (len, magic) = tail.split_at(8);
    (magic == FOOTER_TRAILER_MAGIC)
        .then(|| u64::from_le_bytes(len.try_into().expect("8-byte length")))
}

fn invalid_footer(what: impl std::fmt::Display) -> io::Error {
    io::Error::new(io::ErrorKind::InvalidData, format!("corrupt run-index footer: {what}"))
}

/// `(offset, length)` of the run-index footer block of a run file
/// `file_len` bytes long whose frames span `[body_start, frames_end)`, as
/// its writer recorded them, or `None` when the frames run to the file's
/// end (the writer wrote no footer).
///
/// Whether a footer exists is the writer's word (`frames_end`), never a
/// guess from the file's last bytes: an unindexed zstd run can end in a
/// frame stored raw, whose last bytes are a record's own and may look like a
/// trailer. The bytes after the frames must be exactly one footer: an index
/// block starting with its magic, then a trailer naming the block's length.
///
/// # Errors
///
/// A failed read; [`io::ErrorKind::InvalidData`] when `frames_end` is not
/// within `[body_start, file_len]` or the bytes after the frames are not
/// one footer (a damaged run file, or a writer's figure that disagrees with
/// the file, which must not be read as frames).
pub(crate) fn locate_footer<S: PositionalSource + ?Sized>(
    src: &S,
    body_start: u64,
    frames_end: u64,
    file_len: u64,
) -> io::Result<Option<(u64, u64)>> {
    if frames_end < body_start || frames_end > file_len {
        return Err(io::Error::new(
            io::ErrorKind::InvalidData,
            format!(
                "the run's frames were written to end at byte {frames_end}, but the file's \
                 frames start at byte {body_start} and it holds {file_len} bytes"
            ),
        ));
    }
    let footer_len = file_len - frames_end;
    if footer_len == 0 {
        return Ok(None);
    }
    let not_a_footer = |what: String| {
        invalid_footer(format!(
            "the {footer_len} bytes after the frames at byte {frames_end} {what}"
        ))
    };
    let Some(block_len) = footer_len.checked_sub(FOOTER_TRAILER_BYTES as u64) else {
        return Err(not_a_footer("are too few for a footer".to_owned()));
    };
    let mut tail = [0u8; FOOTER_TRAILER_BYTES];
    read_at_exact(src, &mut tail, file_len - FOOTER_TRAILER_BYTES as u64)?;
    match parse_footer_trailer(&tail) {
        Some(named) if named == block_len => {}
        Some(named) => {
            return Err(not_a_footer(format!(
                "end in a trailer naming a {named}-byte block, not {block_len}"
            )));
        }
        None => return Err(not_a_footer("do not end in a footer trailer".to_owned())),
    }
    let mut head = [0u8; FOOTER_BLOCK_MAGIC.len()];
    if block_len < head.len() as u64 {
        return Err(not_a_footer("hold no block magic".to_owned()));
    }
    read_at_exact(src, &mut head, frames_end)?;
    if head != FOOTER_BLOCK_MAGIC {
        return Err(not_a_footer("do not start with the block magic".to_owned()));
    }
    Ok(Some((frames_end, block_len)))
}

/// A footer block's fixed header (magic through `n`), checked against
/// itself: magic, version, a non-zero stride and `n` entries covering
/// `frames` frames at that stride.
struct FooterHeader {
    key_width: u16,
    stride: u32,
    frames: u64,
    max_frame_raw: u32,
    n: usize,
}

impl FooterHeader {
    fn parse(head: &[u8; FOOTER_HEADER_BYTES]) -> io::Result<Self> {
        let le_u16 = |at: usize| u16::from_le_bytes(head[at..at + 2].try_into().expect("2 bytes"));
        let le_u32 = |at: usize| u32::from_le_bytes(head[at..at + 4].try_into().expect("4 bytes"));
        let le_u64 = |at: usize| u64::from_le_bytes(head[at..at + 8].try_into().expect("8 bytes"));
        // magic 0..4, version 4..8, key width 8..10, pad 10..12, stride
        // 12..16, frames 16..24, max frame 24..28, n 28..36.
        if head[..4] != FOOTER_BLOCK_MAGIC {
            return Err(invalid_footer("bad block magic"));
        }
        let version = le_u32(4);
        if version != FOOTER_VERSION {
            return Err(invalid_footer(format!("unsupported version {version}")));
        }
        let (key_width, stride, frames, max_frame_raw, n) =
            (le_u16(8), le_u32(12), le_u64(16), le_u32(24), le_u64(28));
        if stride == 0 {
            return Err(invalid_footer(format!("stride 0 for {frames} frames in {n} entries")));
        }
        if frames.div_ceil(u64::from(stride)) != n {
            return Err(invalid_footer(format!(
                "{n} entries, but {frames} frames at stride {stride} make {}",
                frames.div_ceil(u64::from(stride))
            )));
        }
        let n =
            usize::try_from(n).map_err(|_| invalid_footer(format!("entry count {n} overflows")))?;
        Ok(Self { key_width, stride, frames, max_frame_raw, n })
    }

    /// Bytes of the `offsets`, `raw_len` and `records` arrays, and of the
    /// keys, that follow the header.
    fn payload_lens(&self) -> io::Result<(usize, usize)> {
        let n = self.n;
        let overflows = || invalid_footer(format!("entry count {n} overflows"));
        let arrays =
            n.checked_mul(8 + 4 + 4).and_then(|b| b.checked_add(8)).ok_or_else(overflows)?;
        let keys = n.checked_mul(usize::from(self.key_width)).ok_or_else(overflows)?;
        arrays.checked_add(keys).ok_or_else(overflows)?;
        Ok((arrays, keys))
    }

    /// The index of this header's arrays (`arrays`, exactly
    /// [`Self::payload_lens`] long) and `keys`, checked against the run's
    /// frames `[body_start, body_end)`.
    fn decode(
        self,
        arrays: &[u8],
        keys: Vec<u8>,
        body_start: u64,
        body_end: u64,
    ) -> io::Result<RunIndex> {
        let n = self.n;
        let (offsets, rest) = arrays.split_at((n + 1) * 8);
        let (raw_len, records) = rest.split_at(n * 4);
        let u64s = |b: &[u8]| -> Vec<u64> {
            b.chunks_exact(8).map(|c| u64::from_le_bytes(c.try_into().expect("8 bytes"))).collect()
        };
        let u32s = |b: &[u8]| -> Vec<u32> {
            b.chunks_exact(4).map(|c| u32::from_le_bytes(c.try_into().expect("4 bytes"))).collect()
        };
        let offsets = u64s(offsets);
        if offsets.windows(2).any(|w| w[0] > w[1]) {
            return Err(invalid_footer("offsets decrease"));
        }
        if offsets[0] != body_start || offsets[n] != body_end {
            return Err(invalid_footer(format!(
                "entries span [{}, {}), the run's frames [{body_start}, {body_end})",
                offsets[0], offsets[n]
            )));
        }
        let Self { key_width, stride, frames, max_frame_raw, .. } = self;
        Ok(RunIndex {
            key_width,
            stride,
            frames,
            max_frame_raw,
            offsets,
            raw_len: u32s(raw_len),
            records: u32s(records),
            keys,
        })
    }
}

/// Decode the footer index block at `offset`, `block_len` bytes long, of a
/// run in `src` whose frames span `[body_start, body_end)`. The header is
/// read first and checked, so a corrupt entry count cannot ask for more
/// memory than the block holds; the arrays are then read and decoded into
/// exactly sized vectors, and the keys are read straight into the index's
/// own buffer. The block must be exactly as long as its header says, its
/// offsets must not decrease, and its entries must tile the frames exactly
/// (first offset `body_start`, last `body_end`).
///
/// # Errors
///
/// A failed read; [`io::ErrorKind::InvalidData`] for a block that is
/// truncated, too long, of the wrong magic or version, or inconsistent with
/// itself or the body.
pub(crate) fn read_footer<S: PositionalSource + ?Sized>(
    src: &S,
    offset: u64,
    block_len: u64,
    body_start: u64,
    body_end: u64,
) -> io::Result<RunIndex> {
    let mut head = [0u8; FOOTER_HEADER_BYTES];
    if block_len < head.len() as u64 {
        return Err(invalid_footer("block is truncated"));
    }
    read_at_exact(src, &mut head, offset)?;
    let header = FooterHeader::parse(&head)?;
    let (arrays_len, keys_len) = header.payload_lens()?;
    let payload = block_len - head.len() as u64;
    if payload != (arrays_len + keys_len) as u64 {
        return Err(invalid_footer(format!(
            "{} entries of {}-byte keys need {} bytes after the header, the block holds \
             {payload}",
            header.n,
            header.key_width,
            arrays_len + keys_len
        )));
    }
    let mut arrays = vec![0u8; arrays_len];
    read_at_exact(src, &mut arrays, offset + head.len() as u64)?;
    let mut keys = vec![0u8; keys_len];
    read_at_exact(src, &mut keys, offset + (head.len() + arrays_len) as u64)?;
    header.decode(&arrays, keys, body_start, body_end)
}

/// A footer block held in memory, read positionally (tests).
#[cfg(test)]
struct BlockSource<'a>(&'a [u8]);

#[cfg(test)]
impl PositionalSource for BlockSource<'_> {
    fn read_at(&self, buf: &mut [u8], offset: u64) -> io::Result<usize> {
        let at = usize::try_from(offset).unwrap_or(usize::MAX).min(self.0.len());
        let n = buf.len().min(self.0.len() - at);
        buf[..n].copy_from_slice(&self.0[at..at + n]);
        Ok(n)
    }

    fn byte_len(&self) -> io::Result<u64> {
        Ok(self.0.len() as u64)
    }
}

impl RunIndex {
    /// Bytes [`Self::write_footer`] writes for this index: the block plus the
    /// trailer.
    #[must_use]
    pub fn footer_bytes(&self) -> u64 {
        let n = self.len() as u64;
        FOOTER_HEADER_BYTES as u64
            + (n + 1) * 8
            + n * 8
            + self.keys.len() as u64
            + FOOTER_TRAILER_BYTES as u64
    }

    /// Append this index as a run-file footer (index block, then trailer; see
    /// the module docs) to `w`. Returns the bytes written
    /// ([`Self::footer_bytes`]).
    ///
    /// # Errors
    ///
    /// Any error `w` returns.
    pub fn write_footer<W: Write>(&self, w: &mut W) -> io::Result<u64> {
        let block_len = self.footer_bytes() - FOOTER_TRAILER_BYTES as u64;
        w.write_all(&FOOTER_BLOCK_MAGIC)?;
        w.write_all(&FOOTER_VERSION.to_le_bytes())?;
        w.write_all(&self.key_width.to_le_bytes())?;
        w.write_all(&0u16.to_le_bytes())?;
        w.write_all(&self.stride.to_le_bytes())?;
        w.write_all(&self.frames.to_le_bytes())?;
        w.write_all(&self.max_frame_raw.to_le_bytes())?;
        w.write_all(&(self.len() as u64).to_le_bytes())?;
        for o in &self.offsets {
            w.write_all(&o.to_le_bytes())?;
        }
        for n in &self.raw_len {
            w.write_all(&n.to_le_bytes())?;
        }
        for n in &self.records {
            w.write_all(&n.to_le_bytes())?;
        }
        w.write_all(&self.keys)?;
        w.write_all(&block_len.to_le_bytes())?;
        w.write_all(&FOOTER_TRAILER_MAGIC)?;
        Ok(block_len + FOOTER_TRAILER_BYTES as u64)
    }

    /// Decode a footer index block (without its trailer) held in memory, of
    /// a run whose frames span `[body_start, body_end)` (see
    /// [`read_footer`]).
    ///
    /// # Errors
    ///
    /// As [`read_footer`].
    #[cfg(test)]
    pub(crate) fn read_from(block: &[u8], body_start: u64, body_end: u64) -> io::Result<RunIndex> {
        read_footer(&BlockSource(block), 0, block.len() as u64, body_start, body_end)
    }
}

fn truncated(what: &str) -> io::Error {
    io::Error::new(io::ErrorKind::UnexpectedEof, format!("truncated spill frame: {what}"))
}

/// Walk one raw (decompressed) frame of `[key][u32 len][body]` records.
///
/// Returns the record count and writes the first record's serialized key into
/// `key_out` (cleared first): the leading key bytes for template kinds, the key
/// derived from the BAM body for coordinate.
///
/// # Errors
///
/// Returns an error for an empty frame, a truncated record, or a kind that is
/// not indexed (queryname).
pub fn frame_meta(kind: SpillKeyKind, raw: &[u8], key_out: &mut Vec<u8>) -> io::Result<u32> {
    key_out.clear();
    if spill_key_width(kind).is_none() {
        return Err(io::Error::new(
            io::ErrorKind::InvalidInput,
            format!("spill key kind {kind:?} is not indexed"),
        ));
    }
    if raw.is_empty() {
        return Err(io::Error::new(io::ErrorKind::InvalidData, "empty spill frame"));
    }
    let key_bytes = kind.serialized_key_bytes();
    let mut pos = 0usize;
    let mut count = 0u32;
    while pos < raw.len() {
        let rec_start = pos;
        let len_at = rec_start + key_bytes;
        let body_at = len_at + RECORD_LEN_BYTES;
        let len_field = raw.get(len_at..body_at).ok_or_else(|| truncated("record header"))?;
        let body_len = u32::from_le_bytes([len_field[0], len_field[1], len_field[2], len_field[3]]);
        let body_end = usize::try_from(body_len)
            .ok()
            .and_then(|n| body_at.checked_add(n))
            .filter(|&end| end <= raw.len())
            .ok_or_else(|| truncated("record body"))?;
        if count == 0 {
            if key_bytes > 0 {
                key_out.extend_from_slice(&raw[rec_start..len_at]);
            } else {
                if body_end - body_at < fgumi_raw_bam::MIN_BAM_RECORD_LEN {
                    return Err(io::Error::new(
                        io::ErrorKind::InvalidData,
                        "coordinate spill record too short to hold its key",
                    ));
                }
                RawCoordinateKey::extract_from_record(&raw[body_at..body_end]).write_to(key_out)?;
            }
        }
        count = count.checked_add(1).ok_or_else(|| truncated("record count overflow"))?;
        pos = body_end;
    }
    Ok(count)
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::{
        CbKey32, TemplateKey, TemplateKey24, TertKey32, frame_keyed_record_into,
        inline::TemplateKey32,
    };
    use rstest::rstest;
    use std::fmt::Debug;

    /// A minimal BAM record body: tid, pos, flag set; the rest zero.
    fn bam_body(tid: i32, pos: i32, flag: u16) -> Vec<u8> {
        let mut b = vec![0u8; 40];
        b[0..4].copy_from_slice(&tid.to_le_bytes());
        b[4..8].copy_from_slice(&pos.to_le_bytes());
        b[14..16].copy_from_slice(&flag.to_le_bytes());
        b
    }

    /// Frame three records with `keys`; assert count and first-key decode.
    fn check_template<K: RawSortKey + Debug>(kind: SpillKeyKind, keys: &[K; 3]) {
        let mut raw = Vec::new();
        for (i, k) in keys.iter().enumerate() {
            frame_keyed_record_into(&mut raw, k, &bam_body(0, i32::try_from(i).unwrap(), 0))
                .unwrap();
        }
        let mut out = vec![0xFF; 3];
        assert_eq!(frame_meta(kind, &raw, &mut out).unwrap(), 3);
        assert_eq!(K::read_from(&mut Cursor::new(&out)).unwrap(), keys[0]);
        assert_eq!(out.len(), kind.serialized_key_bytes());
    }

    /// Every kind's frame of three records reports a count of 3 and its
    /// first record's key.
    #[rstest]
    #[case::coordinate(SpillKeyKind::Coordinate)]
    #[case::k24(SpillKeyKind::TemplateK24)]
    #[case::cb32(SpillKeyKind::TemplateCb32)]
    #[case::tert32(SpillKeyKind::TemplateTert32)]
    #[case::k40(SpillKeyKind::TemplateK40)]
    fn frame_meta_reports_first_key_and_count(#[case] kind: SpillKeyKind) {
        let k32 = |n: u64| TemplateKey32 {
            primary: n,
            secondary: n + 1,
            opt: n + 2,
            name_hash_upper: n + 3,
        };
        match kind {
            SpillKeyKind::Coordinate => {
                let mut raw = Vec::new();
                let bodies = [bam_body(2, 100, 16), bam_body(2, 200, 0), bam_body(3, 5, 0)];
                for b in &bodies {
                    frame_keyed_record_into(&mut raw, &RawCoordinateKey::default(), b).unwrap();
                }
                let mut out = Vec::new();
                assert_eq!(frame_meta(kind, &raw, &mut out).unwrap(), 3);
                let got = RawCoordinateKey::read_from(&mut Cursor::new(&out)).unwrap();
                assert_eq!(got, RawCoordinateKey::new(2, 100, true, 3));
            }
            SpillKeyKind::TemplateK24 => {
                let k =
                    |n: u64| TemplateKey24 { primary: n, secondary: n + 1, name_hash_upper: n + 2 };
                check_template(kind, &[k(1), k(10), k(20)]);
            }
            SpillKeyKind::TemplateCb32 => {
                check_template(kind, &[CbKey32(k32(1)), CbKey32(k32(10)), CbKey32(k32(20))]);
            }
            SpillKeyKind::TemplateTert32 => {
                check_template(kind, &[TertKey32(k32(1)), TertKey32(k32(10)), TertKey32(k32(20))]);
            }
            SpillKeyKind::TemplateK40 => {
                let k = |n: u64| TemplateKey {
                    primary: n,
                    secondary: n + 1,
                    cb_hash: n + 2,
                    tertiary: n + 3,
                    name_hash_upper: n + 4,
                };
                check_template(kind, &[k(1), k(10), k(20)]);
            }
            _ => unreachable!("queryname kinds are covered by frame_meta_rejects"),
        }
    }

    #[rstest]
    #[case::empty(SpillKeyKind::Coordinate, vec![], "empty")]
    #[case::truncated_len(SpillKeyKind::Coordinate, vec![1, 0], "truncated spill frame: record header")]
    #[case::truncated_body(
        SpillKeyKind::Coordinate,
        [100u32.to_le_bytes().to_vec(), vec![0u8; 10]].concat(),
        "truncated spill frame: record body"
    )]
    #[case::short_coordinate_body(
        SpillKeyKind::Coordinate,
        [4u32.to_le_bytes().to_vec(), vec![0u8; 4]].concat(),
        "too short"
    )]
    #[case::truncated_template_key(
        SpillKeyKind::TemplateK24,
        vec![0u8; 10],
        "truncated spill frame: record header"
    )]
    #[case::queryname(SpillKeyKind::QuerynameLex, vec![0u8; 8], "not indexed")]
    fn frame_meta_rejects(#[case] kind: SpillKeyKind, #[case] raw: Vec<u8>, #[case] needle: &str) {
        let err = frame_meta(kind, &raw, &mut Vec::new()).unwrap_err();
        assert!(err.to_string().contains(needle), "{err}");
    }

    #[rstest]
    #[case::coordinate(SpillKeyKind::Coordinate, Some(8))]
    #[case::k24(SpillKeyKind::TemplateK24, Some(24))]
    #[case::cb32(SpillKeyKind::TemplateCb32, Some(32))]
    #[case::tert32(SpillKeyKind::TemplateTert32, Some(32))]
    #[case::k40(SpillKeyKind::TemplateK40, Some(40))]
    #[case::qn_lex(SpillKeyKind::QuerynameLex, None)]
    #[case::qn_natural(SpillKeyKind::QuerynameNatural, None)]
    fn spill_key_width_by_kind(#[case] kind: SpillKeyKind, #[case] want: Option<u16>) {
        assert_eq!(spill_key_width(kind), want);
    }

    fn ckey(n: u64) -> [u8; 8] {
        RawCoordinateKey { sort_key: n }.sort_key.to_le_bytes()
    }

    #[test]
    fn run_index_round_trips() {
        let mut a = RunIndex::new(8);
        assert!(a.is_empty());
        a.push_frame(10, 100, 4, &ckey(1));
        a.push_frame(20, 200, 5, &ckey(2));
        a.push_frame(30, 300, 6, &ckey(3));
        let mut b = RunIndex::new(8);
        b.push_frame(40, 400, 7, &ckey(4));
        b.push_frame(50, 500, 8, &ckey(5));
        a.append(&b);

        assert_eq!(a.len(), 5);
        assert!(!a.is_empty());
        assert_eq!(a.key_width(), 8);
        assert_eq!(a.comp_len(3), 40);
        assert_eq!(a.raw_len(1), 200);
        assert_eq!(a.records(2), 6);
        assert_eq!(a.key_bytes(0), ckey(1));
        assert_eq!(a.total_raw_bytes(), 1500);
        assert_eq!(a.total_records(), 30);
        assert_eq!(a.key_at::<RawCoordinateKey>(4).unwrap(), RawCoordinateKey { sort_key: 5 });
        assert!(a.heap_bytes() >= 6 * 8 + 5 * 8 + 40);
        assert_eq!(a.offset(3), 60);
        assert_eq!(a.end(4), 150);
    }

    fn pushed(frames: &[(u32, u32, u32)]) -> RunIndex {
        let mut idx = RunIndex::new(8);
        for (i, &(comp, raw, records)) in frames.iter().enumerate() {
            idx.push_frame(comp, raw, records, &ckey(i as u64));
        }
        idx
    }

    #[test]
    fn offsets_are_cumulative_and_rebased() {
        let mut idx = pushed(&[(10, 100, 3), (20, 200, 4), (30, 300, 5)]);
        idx.rebase(7);
        assert_eq!(idx.offset(0), 7);
        assert_eq!(idx.end(0), 17);
        assert_eq!(idx.offset(2), 37);
        assert_eq!(idx.end(2), 67);
        assert_eq!(idx.comp_len(1), 20);
        assert_eq!(idx.covered_frames(), 3);
        assert_eq!(idx.stride(), 1);
    }

    #[test]
    fn append_rebased_continues_the_file() {
        let mut a = pushed(&[(10, 100, 1), (20, 200, 1)]);
        a.rebase(4);
        let b = pushed(&[(5, 50, 1)]);
        a.append_rebased(&b, 34);
        assert_eq!(a.offset(2), 34);
        assert_eq!(a.end(2), 39);
        assert_eq!(a.len(), 3);
        assert_eq!(a.covered_frames(), 3);
    }

    #[rstest]
    #[case::by_two(2, vec![(0, 30, 300, 7), (30, 70, 700, 11), (100, 50, 500, 5)])]
    #[case::by_four(4, vec![(0, 100, 1000, 18), (100, 50, 500, 5)])]
    #[case::by_one_is_identity(1, vec![
        (0, 10, 100, 3), (10, 20, 200, 4), (30, 30, 300, 5), (60, 40, 400, 6), (100, 50, 500, 5),
    ])]
    fn sparsify_table(#[case] factor: u32, #[case] expected: Vec<(u64, u64, u32, u32)>) {
        let mut idx =
            pushed(&[(10, 100, 3), (20, 200, 4), (30, 300, 5), (40, 400, 6), (50, 500, 5)]);
        idx.sparsify(factor).unwrap();
        assert_eq!(idx.len(), expected.len());
        for (i, &(offset, comp, raw, records)) in expected.iter().enumerate() {
            assert_eq!(idx.offset(i), offset, "offset {i}");
            assert_eq!(idx.comp_len(i), comp, "comp_len {i}");
            assert_eq!(idx.raw_len(i), raw, "raw_len {i}");
            assert_eq!(idx.records(i), records, "records {i}");
            assert_eq!(idx.key_bytes(i), ckey((i * factor as usize) as u64), "key {i}");
        }
        assert_eq!(idx.covered_frames(), 5);
        assert_eq!(idx.max_frame_raw(), 500);
        assert_eq!(idx.stride(), factor);
        assert_eq!(idx.total_raw_bytes(), 1500);
        assert_eq!(idx.total_records(), 23);
        assert_eq!(idx.end(idx.len() - 1), 150);
    }

    /// A group whose sums overflow `u32` fails and leaves the index as it
    /// was. Columns: frames (raw, records), factor, error needle.
    #[rstest]
    #[case::raw_overflow(vec![(u32::MAX - 10, 1), (11, 1)], 2, "raw length")]
    #[case::records_overflow(vec![(10, u32::MAX), (10, 1)], 2, "record count")]
    #[case::fits_exactly(vec![(u32::MAX - 10, 1), (10, 1)], 2, "")]
    fn sparsify_overflow_leaves_the_index_unchanged(
        #[case] frames: Vec<(u32, u32)>,
        #[case] factor: u32,
        #[case] needle: &str,
    ) {
        let mut idx = RunIndex::new(8);
        for (i, &(raw, records)) in frames.iter().enumerate() {
            idx.push_frame(10, raw, records, &ckey(i as u64));
        }
        let before = idx.clone();
        match idx.sparsify(factor) {
            Ok(()) => {
                assert!(needle.is_empty(), "expected an error naming {needle:?}");
                assert_eq!(idx.len(), 1);
                assert_eq!(idx.max_frame_raw(), before.max_frame_raw());
            }
            Err(e) => {
                assert!(!needle.is_empty() && e.to_string().contains(needle), "{e}");
                assert_eq!(idx, before);
            }
        }
    }

    /// Folding a run's blocks at a stride, block by block, makes the index
    /// sparsifying all its frames dense would, whatever the block cuts (a
    /// short entry carried across blocks, a block shorter than the stride).
    /// Columns: frames per block, stride.
    #[rstest]
    #[case::dense(vec![3, 4], 1)]
    #[case::short_entry_carried(vec![3, 4, 2], 2)]
    #[case::blocks_shorter_than_the_stride(vec![1, 1, 1, 1, 1], 4)]
    #[case::one_block(vec![9], 4)]
    #[case::exact_groups(vec![4, 8], 4)]
    fn append_folded_matches_sparsifying_the_dense_index(
        #[case] blocks: Vec<u32>,
        #[case] stride: u32,
    ) {
        let base = 4u64;
        let mut dense = RunIndex::new(8);
        let mut folded = RunIndex::new_strided(8, stride);
        folded.rebase(base);
        let mut at = base;
        let mut frame = 0u64;
        for len in blocks {
            let mut block = RunIndex::new(8);
            for _ in 0..len {
                let (comp, raw, records) = (10 + frame % 7, 1000 + frame, 1 + frame % 3);
                block.push_frame(
                    u32::try_from(comp).unwrap(),
                    u32::try_from(raw).unwrap(),
                    u32::try_from(records).unwrap(),
                    &ckey(frame),
                );
                dense.push_frame(
                    u32::try_from(comp).unwrap(),
                    u32::try_from(raw).unwrap(),
                    u32::try_from(records).unwrap(),
                    &ckey(frame),
                );
                frame += 1;
            }
            folded.append_folded(&block, at).unwrap();
            at = folded.end_of_last();
        }
        dense.rebase(base);
        dense.sparsify(stride).unwrap();
        assert_eq!(folded, dense);
    }

    /// An entry whose sums would overflow leaves the folded index as it was.
    /// Columns: frames already folded (raw, records), frames appended, error
    /// needle (empty: fits).
    #[rstest]
    #[case::raw_overflow(vec![(u32::MAX - 10, 1)], vec![(11, 1)], "raw length")]
    #[case::records_overflow(vec![(10, u32::MAX)], vec![(10, 1)], "record count")]
    #[case::fits_exactly(vec![(u32::MAX - 10, 1)], vec![(10, 1)], "")]
    #[case::next_entry_starts_fresh(vec![(u32::MAX - 5, 1), (5, 1)], vec![(u32::MAX, 1)], "")]
    fn append_folded_overflow_leaves_the_index_unchanged(
        #[case] first: Vec<(u32, u32)>,
        #[case] then: Vec<(u32, u32)>,
        #[case] needle: &str,
    ) {
        let block = |frames: &[(u32, u32)], from: usize| {
            let mut b = RunIndex::new(8);
            for (i, &(raw, records)) in frames.iter().enumerate() {
                b.push_frame(10, raw, records, &ckey((from + i) as u64));
            }
            b
        };
        let mut idx = RunIndex::new_strided(8, 2);
        idx.append_folded(&block(&first, 0), 0).unwrap();
        let before = idx.clone();
        let end = idx.end_of_last();
        match idx.append_folded(&block(&then, first.len()), end) {
            Ok(()) => {
                assert!(needle.is_empty(), "expected an error naming {needle:?}");
                assert_eq!(idx.covered_frames(), (first.len() + then.len()) as u64);
            }
            Err(e) => {
                assert!(!needle.is_empty() && e.to_string().contains(needle), "{e}");
                assert_eq!(idx, before);
            }
        }
    }

    /// The largest frame survives appends and sparsifying. Columns: first
    /// run's raw lengths, second run's, factor, expected max frame, max entry.
    #[rstest]
    #[case::dense(vec![100, 300], vec![200], 1, 300, 300)]
    #[case::appended_max(vec![100, 300], vec![700], 1, 700, 700)]
    #[case::sparse(vec![100, 300], vec![700, 10], 2, 700, 710)]
    fn max_frame_raw_table(
        #[case] a: Vec<u32>,
        #[case] b: Vec<u32>,
        #[case] factor: u32,
        #[case] max_frame: u32,
        #[case] max_entry: u32,
    ) {
        let build = |raws: &[u32]| {
            let mut idx = RunIndex::new(8);
            for (i, &raw) in raws.iter().enumerate() {
                idx.push_frame(10, raw, 1, &ckey(i as u64));
            }
            idx
        };
        let mut idx = build(&a);
        idx.append(&build(&b));
        idx.sparsify(factor).unwrap();
        assert_eq!(idx.max_frame_raw(), max_frame);
        assert_eq!(idx.max_entry_raw(), max_entry);
    }

    #[rstest]
    #[case::narrow(8, 24)]
    #[case::wide(40, 56)]
    fn entry_bytes_by_width(#[case] width: u16, #[case] expected: u64) {
        assert_eq!(RunIndex::entry_bytes(width), expected);
    }

    #[test]
    fn shrink_to_fit_makes_heap_bytes_exact() {
        let mut idx = RunIndex::new(8);
        for i in 0..1000u64 {
            idx.push_frame(10, 100, 1, &ckey(i));
        }
        idx.shrink_to_fit();
        assert_eq!(idx.heap_bytes(), 1001 * 8 + 1000 * (4 + 4 + 8));
    }

    /// An index of `frames` frames with `key_width`-byte keys (frame `i`:
    /// `10 + i % 7` compressed bytes, `1000 + i` raw, `1 + i % 3` records),
    /// starting at byte 4 and sparsified to `stride`.
    fn footer_index(key_width: u16, frames: u32, stride: u32) -> RunIndex {
        let mut idx = RunIndex::new(key_width);
        for i in 0..frames {
            let key: Vec<u8> =
                (0..key_width).map(|b| u8::try_from((i + u32::from(b)) % 256).unwrap()).collect();
            idx.push_frame(10 + i % 7, 1000 + i, 1 + i % 3, &key);
        }
        idx.rebase(4);
        idx.sparsify(stride).unwrap();
        idx
    }

    /// The footer and its block decode to the index written, at every key
    /// width, for an empty, a one-frame and a many-frame run, dense and
    /// sparsified; the trailer names the block's length.
    #[rstest]
    fn footer_round_trips(
        #[values(8, 24, 32, 40)] key_width: u16,
        #[values(0, 1, 1_000)] frames: u32,
        #[values(1, 4)] stride: u32,
    ) {
        let idx = footer_index(key_width, frames, stride);
        let mut footer = Vec::new();
        let written = idx.write_footer(&mut footer).unwrap();
        assert_eq!(written, footer.len() as u64);
        assert_eq!(written, idx.footer_bytes());
        let (block, trailer) = footer.split_at(footer.len() - FOOTER_TRAILER_BYTES);
        assert_eq!(parse_footer_trailer(trailer.try_into().unwrap()), Some(block.len() as u64));
        assert_eq!(block[..4], FOOTER_BLOCK_MAGIC);
        let body_end = idx.end_of_last();
        assert_eq!(RunIndex::read_from(block, 4, body_end).unwrap(), idx);
    }

    /// A run file's last 16 bytes are a trailer only with the footer magic.
    #[rstest]
    #[case::footer(*b"\x07\0\0\0\0\0\0\0FGUMIIDX", Some(7))]
    #[case::frame_bytes([0xAB; 16], None)]
    #[case::bgzf_eof_tail(*b"\0\0\0\0\0\0\0\0\0\0\0\0\0\0\0\0", None)]
    fn footer_trailer_needs_its_magic(#[case] tail: [u8; 16], #[case] expected: Option<u64>) {
        assert_eq!(parse_footer_trailer(&tail), expected);
    }

    /// How a test corrupts a footer block of a three-frame, 8-byte-key index
    /// spanning `[4, 34)`.
    #[derive(Clone, Copy, Debug)]
    enum Corrupt {
        /// Drop the block's last byte.
        Truncated,
        /// Append a byte.
        TooLong,
        /// Overwrite the magic.
        Magic,
        /// Set the version to 2.
        Version,
        /// Claim one more entry than the frames and stride make.
        EntryCount,
        /// Claim more frames than the entries cover at the stride.
        FrameCount,
        /// Claim a stride the entries do not have.
        Stride,
        /// Claim one more entry and frame than the arrays hold, consistently.
        EntryAndFrameCount,
        /// Read the block against a body ending one byte later.
        BodyEnd,
        /// Read the block against a body starting one byte earlier.
        BodyStart,
        /// Swap two offsets so they decrease.
        Decreasing,
        /// Set the stride to 0.
        ZeroStride,
    }

    /// A block that is truncated, too long, of the wrong magic or version,
    /// inconsistent with itself (its entry count is not its frames at its
    /// stride, rounded up), or not tiling the run's frames exactly is
    /// `InvalidData`. Columns: the corruption, a needle of the error.
    #[rstest]
    #[case::truncated(Corrupt::Truncated, "the block holds")]
    #[case::too_long(Corrupt::TooLong, "the block holds")]
    #[case::magic(Corrupt::Magic, "bad block magic")]
    #[case::version(Corrupt::Version, "unsupported version 2")]
    #[case::entry_count(Corrupt::EntryCount, "4 entries, but 3 frames at stride 1 make 3")]
    #[case::frame_count(Corrupt::FrameCount, "3 entries, but 5 frames at stride 1 make 5")]
    #[case::stride(Corrupt::Stride, "3 entries, but 3 frames at stride 2 make 2")]
    #[case::entry_and_frame_count(Corrupt::EntryAndFrameCount, "the block holds")]
    #[case::body_end(Corrupt::BodyEnd, "the run's frames")]
    #[case::body_start(Corrupt::BodyStart, "the run's frames")]
    #[case::decreasing(Corrupt::Decreasing, "offsets decrease")]
    #[case::zero_stride(Corrupt::ZeroStride, "stride 0")]
    fn footer_rejects(#[case] corrupt: Corrupt, #[case] needle: &str) {
        let idx = footer_index(8, 3, 1);
        assert_eq!((idx.offset(0), idx.end_of_last()), (4, 4 + 10 + 11 + 12));
        let mut block = Vec::new();
        idx.write_footer(&mut block).unwrap();
        block.truncate(block.len() - FOOTER_TRAILER_BYTES);
        let (mut start, mut end) = (4, idx.end_of_last());
        // Header: magic 0..4, version 4..8, key width 8..10, pad 10..12,
        // stride 12..16, frames 16..24, max frame 24..28, n 28..36, offsets.
        match corrupt {
            Corrupt::Truncated => {
                block.pop();
            }
            Corrupt::TooLong => block.push(0),
            Corrupt::Magic => block[0] = b'X',
            Corrupt::Version => block[4..8].copy_from_slice(&2u32.to_le_bytes()),
            Corrupt::EntryCount => block[28..36].copy_from_slice(&4u64.to_le_bytes()),
            Corrupt::FrameCount => block[16..24].copy_from_slice(&5u64.to_le_bytes()),
            Corrupt::Stride => block[12..16].copy_from_slice(&2u32.to_le_bytes()),
            Corrupt::EntryAndFrameCount => {
                block[16..24].copy_from_slice(&4u64.to_le_bytes());
                block[28..36].copy_from_slice(&4u64.to_le_bytes());
            }
            Corrupt::BodyEnd => end += 1,
            Corrupt::BodyStart => start -= 1,
            Corrupt::Decreasing => {
                let (a, b) = (block[36..44].to_vec(), block[44..52].to_vec());
                block[36..44].copy_from_slice(&b);
                block[44..52].copy_from_slice(&a);
            }
            Corrupt::ZeroStride => block[12..16].copy_from_slice(&0u32.to_le_bytes()),
        }
        let err = RunIndex::read_from(&block, start, end).unwrap_err();
        assert_eq!(err.kind(), io::ErrorKind::InvalidData, "{err}");
        assert!(err.to_string().contains(needle), "{err}");
    }

    #[test]
    fn key_at_rejects_width_mismatch() {
        let mut idx = RunIndex::new(24);
        idx.push_frame(1, 1, 1, &[0u8; 24]);
        assert!(idx.key_at::<TemplateKey>(0).is_err());
    }
}
