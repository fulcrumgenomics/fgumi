//! Cooperative consolidation of spill runs for the arena sort path.
//!
//! The arena spill path (`SpillGather` → `SpillBlockCompress` → `SpillWrite`)
//! writes one file per sorted run. `--max-temp-files` bounds how many of those
//! runs may be live at once: when the limit would be exceeded, a contiguous
//! range of runs is merged into one. This module is the merge kernel.
//!
//! It is built from the same pieces the arena path already uses, rather than
//! the `RawExternalSorter` consolidation (which needs a `SortWorkerPool` and
//! spawns a reader thread per input):
//!
//! - inputs are read with a plain [`BufReader`] per file and one shared
//!   [`SpillBlockDecompressor`], so any codec the spill writer produced reads
//!   back;
//! - records are parsed with the spill framing `[key?][u32 LE len][body]`,
//!   embedded keys being re-extracted from the body exactly as the Phase-2 slot
//!   parser does;
//! - the output is framed with [`frame_keyed_record_into`] into raw blocks,
//!   compressed with a [`SpillBlockCompressor`] and bracketed by
//!   [`spill_magic`] / [`spill_trailer`], so it is indistinguishable from a
//!   spill file the writer produced directly.
//!
//! The merge is **resumable**: [`RunMergerDyn::step`] merges at most a bounded
//! number of records and returns, so a pipeline step can drive it cooperatively
//! without blocking inside `try_run`.
//!
//! # Output identity
//!
//! Ties between equal keys are broken by input position (the loser tree's leaf
//! index), and the caller passes inputs in run order. A stable merge of a
//! contiguous range of runs, placed at that range's position, is therefore a
//! re-association of the final merge and leaves the sorted output unchanged.

use std::fs::{File, OpenOptions};
use std::io::{self, BufReader, BufWriter, Seek, SeekFrom, Write};
use std::path::{Path, PathBuf};

use crate::codec::SpillCodec;
use crate::inline::{CbKey32, TemplateKey, TemplateKey24, TertKey32};
use crate::keys::{RawCoordinateKey, RawQuerynameKey, RawQuerynameLexKey, RawSortKey};
use crate::loser_tree::LoserTree;
use crate::spill_block::{
    SpillBlockCompressor, frame_keyed_record_into, spill_magic, spill_trailer,
};
use crate::spill_block_reader::SpillBlockDecompressor;

/// Read-buffer capacity per consolidation input. Small on purpose: a
/// consolidation may hold dozens of inputs open at once and the buffers are not
/// charged against `--max-memory`.
const INPUT_BUFFER_BYTES: usize = 64 * 1024;

/// The sort key type a spill run was written with.
///
/// Consolidation must decode keys at the width they were written, and for
/// template-coordinate that width is only chosen at runtime (from the first
/// record), so the kind travels with the spilled data instead of being fixed
/// when the chain is built.
#[derive(Clone, Copy, Debug, PartialEq, Eq, Hash)]
pub enum SpillKeyKind {
    /// Coordinate order ([`RawCoordinateKey`], embedded in the record).
    Coordinate,
    /// Queryname order, lexicographic ([`RawQuerynameLexKey`]).
    QuerynameLex,
    /// Queryname order, natural ([`RawQuerynameKey`]).
    QuerynameNatural,
    /// Template-coordinate, 24-byte core-only lane ([`TemplateKey24`]).
    TemplateK24,
    /// Template-coordinate, 32-byte cell-barcode lane ([`CbKey32`]).
    TemplateCb32,
    /// Template-coordinate, 32-byte tertiary lane ([`TertKey32`]).
    TemplateTert32,
    /// Template-coordinate, full 40-byte key ([`TemplateKey`]).
    TemplateK40,
}

/// Outcome of one bounded [`RunMergerDyn::step`].
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum RunMergeProgress {
    /// Records were merged and more remain.
    Working,
    /// The merge is complete and the output file is finished and closed.
    Done {
        /// Records written to the output run.
        records: u64,
    },
}

/// Object-safe view over a key-typed run merger, so a key-agnostic pipeline
/// step can drive whichever [`SpillKeyKind`] its runs were written with.
pub trait RunMergerDyn: Send {
    /// Merge records into the output until `max_records` records or `max_bytes`
    /// bytes of record data have been written by this call, whichever is first
    /// (at least one record is merged per call while any remain).
    ///
    /// Returns [`RunMergeProgress::Done`] once every input is exhausted and the
    /// output has been finalized; further calls return `Done` again.
    ///
    /// # Errors
    ///
    /// Propagates read, decompression, framing, compression and write errors.
    /// A spill run that ends in the middle of a record is reported as
    /// [`io::ErrorKind::UnexpectedEof`].
    fn step(&mut self, max_records: usize, max_bytes: u64) -> io::Result<RunMergeProgress>;

    /// Records written to the output so far.
    fn records_written(&self) -> u64;
}

/// What to merge and where to write it.
#[derive(Debug)]
pub struct RunMergeSpec<'a> {
    /// Input spill runs, in run order (oldest first). Ties between equal keys
    /// are resolved in this order.
    pub inputs: &'a [PathBuf],
    /// Output path. Must not exist.
    pub output: &'a Path,
    /// Codec to write the output with. Inputs are sniffed individually.
    pub codec: SpillCodec,
    /// Compression level for the output (bgzf level, or zstd level ≥ 1).
    pub compression: u32,
    /// Target raw (uncompressed) block size for the output.
    pub block_size: usize,
}

/// Build a merger for runs written with `kind`, opening every input and the
/// output.
///
/// Opening primes each input with its first record, so the returned merger
/// holds one open descriptor per input plus one for the output until it
/// reports [`RunMergeProgress::Done`] or is dropped.
///
/// # Errors
///
/// Returns an error if an input cannot be opened or its first record cannot be
/// read, or if the output cannot be created.
pub fn new_run_merger(
    kind: SpillKeyKind,
    spec: &RunMergeSpec<'_>,
) -> io::Result<Box<dyn RunMergerDyn>> {
    Ok(match kind {
        SpillKeyKind::Coordinate => Box::new(RunMerger::<RawCoordinateKey>::open(spec)?),
        SpillKeyKind::QuerynameLex => Box::new(RunMerger::<RawQuerynameLexKey>::open(spec)?),
        SpillKeyKind::QuerynameNatural => Box::new(RunMerger::<RawQuerynameKey>::open(spec)?),
        SpillKeyKind::TemplateK24 => Box::new(RunMerger::<TemplateKey24>::open(spec)?),
        SpillKeyKind::TemplateCb32 => Box::new(RunMerger::<CbKey32>::open(spec)?),
        SpillKeyKind::TemplateTert32 => Box::new(RunMerger::<TertKey32>::open(spec)?),
        SpillKeyKind::TemplateK40 => Box::new(RunMerger::<TemplateKey>::open(spec)?),
    })
}

// ============================================================================
// Input: one spill run, read record by record
// ============================================================================

/// Sequential record reader over one spill run file.
struct SpillRunReader {
    file: BufReader<File>,
    codec: SpillCodec,
    /// Current decompressed block and the read position within it.
    block: Vec<u8>,
    pos: usize,
    /// `true` once the file has no further blocks.
    eof: bool,
    path: PathBuf,
}

impl SpillRunReader {
    /// Open `path`, detecting its codec from the file magic the same way
    /// `open_spill_slot` does: a `ZSP1` prologue selects zstd (and is consumed);
    /// anything else is a BGZF stream read from offset 0.
    fn open(path: &Path) -> io::Result<Self> {
        let mut file = File::open(path).map_err(|e| {
            io::Error::new(e.kind(), format!("failed to open spill run {}: {e}", path.display()))
        })?;
        // Same sniff as `open_spill_slot`: a short read of the prologue is a
        // truncated file, not a BGZF stream, so it must not fall through to
        // the BGZF reader as a partial magic would.
        let mut magic = [0u8; 4];
        let filled = crate::external::read_exact_or_eof(&mut file, &mut magic)?;
        let codec = if filled {
            SpillCodec::from_magic(&magic).unwrap_or(SpillCodec::Bgzf)
        } else {
            SpillCodec::Bgzf
        };
        if codec == SpillCodec::Bgzf {
            file.seek(SeekFrom::Start(0))?;
        }
        Ok(Self {
            file: BufReader::with_capacity(INPUT_BUFFER_BYTES, file),
            codec,
            block: Vec::new(),
            pos: 0,
            eof: false,
            path: path.to_path_buf(),
        })
    }

    /// Ensure at least one unread byte is buffered. Returns `false` at the end
    /// of the file.
    fn fill(&mut self, dec: &mut SpillBlockDecompressor) -> io::Result<bool> {
        while self.pos == self.block.len() {
            if self.eof {
                return Ok(false);
            }
            // Empty blocks (e.g. a BGZF EOF marker) are skipped by looping.
            if let Some(block) = dec.read_blocks(&mut self.file, self.codec, 1)?.pop() {
                self.block = block;
                self.pos = 0;
            } else {
                self.eof = true;
                self.block.clear();
                self.pos = 0;
            }
        }
        Ok(true)
    }

    /// Append exactly `n` bytes to `out`, crossing block boundaries as needed.
    fn take_into(
        &mut self,
        dec: &mut SpillBlockDecompressor,
        n: usize,
        out: &mut Vec<u8>,
    ) -> io::Result<()> {
        let mut remaining = n;
        while remaining > 0 {
            if !self.fill(dec)? {
                return Err(io::Error::new(
                    io::ErrorKind::UnexpectedEof,
                    format!("spill run {} ends in the middle of a record", self.path.display()),
                ));
            }
            let take = remaining.min(self.block.len() - self.pos);
            out.extend_from_slice(&self.block[self.pos..self.pos + take]);
            self.pos += take;
            remaining -= take;
        }
        Ok(())
    }

    /// Read the next record's body into `body` and return its key, or `None`
    /// at a clean end of file. `key_buf` is scratch for a serialized key prefix.
    fn next_record<K: RawSortKey + Default + 'static>(
        &mut self,
        dec: &mut SpillBlockDecompressor,
        key_size: usize,
        key_buf: &mut Vec<u8>,
        body: &mut Vec<u8>,
    ) -> io::Result<Option<K>> {
        if !self.fill(dec)? {
            return Ok(None);
        }
        key_buf.clear();
        self.take_into(dec, key_size, key_buf)?;
        let key_len = key_buf.len();
        // The length prefix lands after the key in the same reused scratch, so
        // reading a record allocates nothing once the buffers have grown.
        self.take_into(dec, 4, key_buf)?;
        let len = u32::from_le_bytes(key_buf[key_len..].try_into().expect("4-byte length prefix"));
        key_buf.truncate(key_len);
        body.clear();
        self.take_into(dec, len as usize, body)?;
        parse_key::<K>(key_buf, body).map(Some).map_err(|e| {
            io::Error::new(e.kind(), format!("spill run {}: {e}", self.path.display()))
        })
    }
}

/// Serialized key prefix size for `K`. The final merge's slot parser owns the
/// framing contract; consolidation reads the same files, so it shares the helper
/// rather than keeping a second copy that could drift from it.
fn key_prefix_size<K: RawSortKey>() -> io::Result<usize> {
    crate::external::slot_key_size::<K>()
        .map_err(|e| io::Error::new(io::ErrorKind::InvalidInput, format!("{e:#}")))
}

/// Recover a record's sort key exactly as the final merge's slot parser does.
fn parse_key<K: RawSortKey + Default + 'static>(key_bytes: &[u8], body: &[u8]) -> io::Result<K> {
    crate::external::slot_parse_key::<K>(key_bytes, body)
        .map_err(|e| io::Error::new(io::ErrorKind::InvalidData, format!("{e:#}")))
}

// ============================================================================
// Output: one spill run, written block by block
// ============================================================================

/// Writes framed records to a new spill file in the spill writer's format.
struct SpillRunWriter {
    file: BufWriter<File>,
    compressor: SpillBlockCompressor,
    codec: SpillCodec,
    raw: Vec<u8>,
    block_size: usize,
}

impl SpillRunWriter {
    fn create(
        path: &Path,
        codec: SpillCodec,
        compression: u32,
        block_size: usize,
    ) -> io::Result<Self> {
        // `create_new` fails closed rather than overwrite a run still in use.
        let file = OpenOptions::new().write(true).create_new(true).open(path).map_err(|e| {
            io::Error::new(e.kind(), format!("failed to create spill run {}: {e}", path.display()))
        })?;
        let mut file = BufWriter::new(file);
        file.write_all(spill_magic(codec))?;
        let block_size = block_size.max(1);
        Ok(Self {
            file,
            compressor: SpillBlockCompressor::new(codec, compression)?,
            codec,
            raw: Vec::with_capacity(block_size),
            block_size,
        })
    }

    fn write_record<K: RawSortKey>(&mut self, key: &K, body: &[u8]) -> io::Result<()> {
        frame_keyed_record_into(&mut self.raw, key, body)?;
        if self.raw.len() >= self.block_size {
            self.flush_block()?;
        }
        Ok(())
    }

    fn flush_block(&mut self) -> io::Result<()> {
        if !self.raw.is_empty() {
            let compressed = self.compressor.compress_block(&self.raw)?;
            self.file.write_all(&compressed)?;
            self.raw.clear();
        }
        Ok(())
    }

    fn finish(mut self) -> io::Result<()> {
        self.flush_block()?;
        self.file.write_all(spill_trailer(self.codec))?;
        self.file.flush()
    }
}

// ============================================================================
// The merger
// ============================================================================

/// Resumable k-way merge of spill runs keyed by `K`.
struct RunMerger<K: RawSortKey + Default + 'static> {
    inputs: Vec<SpillRunReader>,
    decompressor: SpillBlockDecompressor,
    key_size: usize,
    key_buf: Vec<u8>,
    /// Current record body per active leaf (parallel to `source_map`).
    records: Vec<Vec<u8>>,
    /// Loser-tree leaf → index into `inputs`.
    source_map: Vec<usize>,
    /// `None` when every input was empty from the start.
    tree: Option<LoserTree<K>>,
    /// `None` once the output has been finalized.
    writer: Option<SpillRunWriter>,
    records_written: u64,
}

impl<K: RawSortKey + Default + 'static> RunMerger<K> {
    fn open(spec: &RunMergeSpec<'_>) -> io::Result<Self> {
        let key_size = key_prefix_size::<K>()?;
        let mut decompressor = SpillBlockDecompressor::new();
        let mut key_buf = Vec::with_capacity(key_size);
        let mut inputs = Vec::with_capacity(spec.inputs.len());
        let mut keys = Vec::with_capacity(spec.inputs.len());
        let mut records = Vec::with_capacity(spec.inputs.len());
        let mut source_map = Vec::with_capacity(spec.inputs.len());
        for path in spec.inputs {
            let mut reader = SpillRunReader::open(path)?;
            let mut body = Vec::new();
            // Leaves are pushed in input order, skipping empty inputs, so the
            // tree's tie-break (leaf index) follows run order.
            if let Some(key) =
                reader.next_record::<K>(&mut decompressor, key_size, &mut key_buf, &mut body)?
            {
                keys.push(key);
                records.push(body);
                source_map.push(inputs.len());
            }
            inputs.push(reader);
        }
        let writer =
            SpillRunWriter::create(spec.output, spec.codec, spec.compression, spec.block_size)?;
        let tree = if keys.is_empty() { None } else { Some(LoserTree::new(keys)) };
        Ok(Self {
            inputs,
            decompressor,
            key_size,
            key_buf,
            records,
            source_map,
            tree,
            writer: Some(writer),
            records_written: 0,
        })
    }
}

impl<K: RawSortKey + Default + 'static> RunMergerDyn for RunMerger<K> {
    fn step(&mut self, max_records: usize, max_bytes: u64) -> io::Result<RunMergeProgress> {
        let Some(writer) = self.writer.as_mut() else {
            return Ok(RunMergeProgress::Done { records: self.records_written });
        };
        if let Some(tree) = self.tree.as_mut() {
            let (mut merged, mut bytes) = (0usize, 0u64);
            while merged < max_records && bytes < max_bytes && tree.winner_is_active() {
                let leaf = tree.winner();
                writer.write_record(tree.winner_key(), &self.records[leaf])?;
                self.records_written += 1;
                merged += 1;
                bytes += self.records[leaf].len() as u64;
                let input = &mut self.inputs[self.source_map[leaf]];
                match input.next_record::<K>(
                    &mut self.decompressor,
                    self.key_size,
                    &mut self.key_buf,
                    &mut self.records[leaf],
                )? {
                    Some(key) => tree.replace_winner(key),
                    None => tree.remove_winner(),
                }
            }
            if tree.winner_is_active() {
                return Ok(RunMergeProgress::Working);
            }
        }
        let writer = self.writer.take().expect("writer present until finished");
        writer.finish()?;
        // Release the input descriptors as soon as the merge is done.
        self.inputs.clear();
        Ok(RunMergeProgress::Done { records: self.records_written })
    }

    fn records_written(&self) -> u64 {
        self.records_written
    }
}

#[cfg(test)]
mod tests;
