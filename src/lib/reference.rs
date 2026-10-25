//! Reference genome FASTA reading with all sequences loaded into memory.
//!
//! This module provides thread-safe access to reference genome sequences, which is needed
//! for tasks like NM/UQ/MD tag calculation and variant calling.
//!
//! Following fgbio's approach, the entire reference is loaded into memory at startup
//! to ensure O(1) lookup performance for each read during tag regeneration.
//!
//! Uses FAI index for fast raw-byte reading (htsjdk-style) instead of line-by-line parsing.
//!
//! # Memory Usage
//!
//! Reference sequences are stored as raw bytes, matching htsjdk's approach.
//! For a typical human reference (~3GB), this uses approximately 3GB of memory
//! but provides the fastest possible load times (~9s vs ~22s for compressed storage).
//!
//! Where fetches are rare, `IndexedReferenceReader` instead holds only the `.fai`
//! index and reads each region from the FASTA when it is fetched, returning the
//! same bytes; `ReferenceSource` lets a consumer take either (both crate-internal).
//!
//! # Future improvement
//!
//! The custom FAI-based raw-byte reading (`read_sequence_raw`) could be replaced with
//! noodles' built-in indexed reader once <https://github.com/zaeleus/noodles/pull/365>
//! is merged and released, which adds the same optimization to noodles.
use crate::errors::FgumiError;
use anyhow::{Context, Result};
use log::debug;
use noodles::core::Position;
use noodles::fasta::fai;
use std::collections::HashMap;
use std::fs::File;
use std::io::{BufReader, Read, Seek, SeekFrom};
use std::path::{Path, PathBuf};
use std::sync::Arc;

/// Read a sequence from a FASTA file using FAI index metadata.
/// Uses raw byte reading with mathematical newline handling (htsjdk-style).
///
/// Optimized: Reads entire sequence bytes in one syscall, then strips newlines in memory.
/// Fast path: If sequence fits in a single line, skip newline stripping entirely.
#[allow(clippy::cast_possible_truncation)]
fn read_sequence_raw(file: &mut File, record: &fai::Record) -> Result<Vec<u8>> {
    let line_bases = record.line_base_count().get() as usize;
    let line_width = record.line_width().get() as usize;
    let seq_len = record.length() as usize;
    let offset = record.position();

    // Fast path: if sequence fits in a single line, no newlines to strip
    if seq_len <= line_bases {
        file.seek(SeekFrom::Start(offset))?;
        let mut sequence = vec![0u8; seq_len];
        file.read_exact(&mut sequence)?;
        return Ok(sequence);
    }

    // Calculate number of complete lines and remaining bases
    let complete_lines = seq_len / line_bases;
    let remaining_bases = seq_len % line_bases;

    // Total bytes = (complete_lines * line_width) + remaining_bases
    // But last line might not have a terminator, so we calculate conservatively
    let total_bytes = if remaining_bases > 0 {
        complete_lines * line_width + remaining_bases
    } else if complete_lines > 0 {
        // All bases fit exactly in complete lines, last line has no terminator at end of seq
        (complete_lines - 1) * line_width + line_bases
    } else {
        0
    };

    // Seek and read all bytes at once
    file.seek(SeekFrom::Start(offset))?;
    let mut raw_bytes = vec![0u8; total_bytes];
    file.read_exact(&mut raw_bytes)?;

    // Strip newlines in memory (much faster than seeking)
    let mut sequence = Vec::with_capacity(seq_len);
    let terminator_len = line_width - line_bases;

    let mut pos = 0;
    while sequence.len() < seq_len && pos < raw_bytes.len() {
        // Read up to line_bases bytes
        let bases_to_copy = (seq_len - sequence.len()).min(line_bases).min(raw_bytes.len() - pos);
        sequence.extend_from_slice(&raw_bytes[pos..pos + bases_to_copy]);
        pos += bases_to_copy;

        // Skip line terminator if present and we need more bases
        if sequence.len() < seq_len && pos < raw_bytes.len() {
            pos += terminator_len;
        }
    }

    Ok(sequence)
}

/// Whether `path` is a regular file, which can be opened more than once and
/// read at any offset. A FIFO or `/dev/fd/N` (`--ref <(cat ref.fa)`) is not: it
/// can be read only once, from the start.
fn is_regular_file(path: &Path) -> bool {
    std::fs::metadata(path).is_ok_and(|m| m.is_file())
}

/// Whether the regular file at `path` starts with the gzip magic (`1f 8b`), as a
/// gzip or BGZF file does. A file that is not regular gives `false` without being
/// opened, since reading its first bytes would consume them; a file that cannot
/// be read gives `false`, so whoever opens it reports the error.
fn is_gzip(path: &Path) -> bool {
    if !is_regular_file(path) {
        return false;
    }
    let Ok(mut file) = File::open(path) else {
        return false;
    };
    let mut magic = [0u8; 2];
    file.read_exact(&mut magic).is_ok() && magic == [0x1f, 0x8b]
}

/// Find a sibling file for a FASTA file by trying two naming conventions:
/// 1. Replace extension (e.g., `ref.fasta` → `ref.<replace_ext>`)
/// 2. Append extension to full path (e.g., `ref.fa` → `ref.fa.<append_ext>`)
fn find_sibling_file(fasta_path: &Path, replace_ext: &str, append_ext: &str) -> Option<PathBuf> {
    let replaced = fasta_path.with_extension(replace_ext);
    if replaced.exists() {
        return Some(replaced);
    }

    let appended = PathBuf::from(format!("{}.{append_ext}", fasta_path.display()));
    if appended.exists() {
        return Some(appended);
    }

    None
}

/// Find FAI index path for a FASTA file.
///
/// Tries multiple naming conventions:
/// 1. Replace extension with `.fa.fai` (e.g., `ref.fasta` → `ref.fa.fai`)
/// 2. Append `.fai` to full path (e.g., `ref.fa` → `ref.fa.fai`, `ref.fasta` → `ref.fasta.fai`)
fn find_fai_path(fasta_path: &Path) -> Option<PathBuf> {
    find_sibling_file(fasta_path, "fa.fai", "fai")
}

/// Find sequence dictionary path for a FASTA file.
///
/// Tries multiple naming conventions used by different tools:
/// 1. Replace extension with `.dict` (fgbio/HTSJDK/Picard convention: `ref.fa` → `ref.dict`)
/// 2. Append `.dict` to full path (GATK convention: `ref.fa` → `ref.fa.dict`)
///
/// # Arguments
/// * `fasta_path` - Path to the FASTA file
///
/// # Returns
/// The path to the dictionary file if found, or `None` if not found.
///
/// # Examples
/// ```no_run
/// use std::path::Path;
/// use fgumi_lib::reference::find_dict_path;
///
/// // Will find either "ref.dict" or "ref.fa.dict"
/// if let Some(dict_path) = find_dict_path(Path::new("ref.fa")) {
///     println!("Found dictionary: {}", dict_path.display());
/// }
/// ```
#[must_use]
pub fn find_dict_path(fasta_path: &Path) -> Option<PathBuf> {
    find_sibling_file(fasta_path, "dict", "dict")
}

/// The 0-based `[start, end)` indices of the 1-based inclusive region
/// `chrom:start-end` on a sequence of `len` bases, or the error both reference
/// readers report for a region past the end of the sequence (or an empty one).
fn region_indices(
    chrom: &str,
    start: Position,
    end: Position,
    len: usize,
) -> Result<(usize, usize)> {
    let start_idx = usize::from(start) - 1;
    let end_idx = usize::from(end);
    if end_idx > len || start_idx >= end_idx {
        return Err(FgumiError::InvalidParameter {
            parameter: "region".to_string(),
            reason: format!("Requested region {chrom}:{start}-{end} exceeds sequence length {len}"),
        }
        .into());
    }
    Ok((start_idx, end_idx))
}

/// A thread-safe reference genome reader with all sequences preloaded into memory.
///
/// This reader loads the entire FASTA file into memory at construction time,
/// providing O(1) lookup performance for sequence fetches. This approach matches
/// fgbio's `nmUqMdTagRegeneratingWriter` which reads all contigs into a Map upfront.
///
/// For a typical human reference (e.g., hs38DH at ~3GB), this uses approximately
/// 3GB of memory (raw byte storage like htsjdk) and provides the fastest load times.
#[derive(Clone)]
pub struct ReferenceReader {
    /// All sequences loaded into memory as raw bytes, keyed by sequence name
    sequences: Arc<HashMap<String, Vec<u8>>>,
}

/// The reference FASTAs a chain has loaded, by path, so stages that read the
/// same FASTA (a methylation-aware consensus caller, realign, clip) share one
/// in-memory copy instead of each loading their own.
#[derive(Default)]
pub(crate) struct ReferenceCache(Vec<(std::path::PathBuf, ReferenceReader)>);

impl ReferenceCache {
    /// The reader for `path` if a stage already loaded it, without loading it.
    pub(crate) fn get(&self, path: &Path) -> Option<ReferenceReader> {
        self.0.iter().find(|(p, _)| p == path).map(|(_, reader)| reader.clone())
    }

    /// The reader for `path`, loading it on first use.
    ///
    /// # Errors
    ///
    /// Returns an error if the FASTA cannot be loaded ([`ReferenceReader::new`]).
    pub(crate) fn load(&mut self, path: &Path) -> Result<ReferenceReader> {
        if let Some((_, reader)) = self.0.iter().find(|(p, _)| p == path) {
            return Ok(reader.clone());
        }
        let reader = ReferenceReader::new(path)?;
        self.0.push((path.to_path_buf(), reader.clone()));
        Ok(reader)
    }
}

impl ReferenceReader {
    /// Creates a new reference reader, loading all sequences into memory.
    ///
    /// This reads the entire FASTA file into memory at construction time.
    /// For a typical human reference (~3GB), this takes a few seconds but
    /// provides O(1) lookup performance for all subsequent fetches.
    ///
    /// A gzip or BGZF FASTA (detected by its magic bytes, whatever its name)
    /// is decompressed as it is read, from the start: its `.fai` (as `samtools
    /// faidx` writes for a BGZF FASTA) holds uncompressed offsets, which cannot
    /// be used to seek in the compressed file.
    ///
    /// # Arguments
    /// * `path` - Path to the reference FASTA file (may be gzip or BGZF compressed)
    ///
    /// # Errors
    /// Returns an error if:
    /// - The file does not exist
    /// - The file cannot be read or parsed as FASTA
    ///
    /// # Examples
    /// ```no_run
    /// use fgumi_lib::reference::ReferenceReader;
    ///
    /// let reader = ReferenceReader::new("reference.fasta")?;
    /// # Ok::<(), anyhow::Error>(())
    /// ```
    pub fn new<P: AsRef<Path>>(path: P) -> Result<Self> {
        let path = path.as_ref();

        // Verify the file exists and is readable
        if !path.exists() {
            return Err(FgumiError::InvalidFileFormat {
                file_type: "Reference FASTA".to_string(),
                path: path.display().to_string(),
                reason: "File does not exist".to_string(),
            }
            .into());
        }

        debug!("Reading reference FASTA into memory: {}", path.display());

        // A FIFO can be opened and read only once, so it is neither probed for
        // compression nor read through an index beside it.
        if !is_regular_file(path) {
            debug!("Not a regular file, using sequential reading");
            return Self::new_sequential(path);
        }

        if is_gzip(path) {
            debug!("Decompressing a gzip/BGZF FASTA sequentially");
            return Self::new_compressed(path);
        }

        // Try FAI-based fast reading first
        if let Some(fai_path) = find_fai_path(path) {
            debug!("Using FAI index for fast loading: {}", fai_path.display());
            return Self::new_with_fai(path, &fai_path);
        }

        // Fall back to noodles for non-indexed files
        debug!("No FAI index found, using sequential reading");
        Self::new_sequential(path)
    }

    /// Load sequences using FAI index for fast raw-byte reading (htsjdk-style).
    fn new_with_fai(fasta_path: &Path, fai_path: &Path) -> Result<Self> {
        let index = fai::fs::read(fai_path)
            .with_context(|| format!("Failed to read FAI index: {}", fai_path.display()))?;
        let records: &[fai::Record] = index.as_ref();
        let mut file = File::open(fasta_path)
            .with_context(|| format!("Failed to open FASTA: {}", fasta_path.display()))?;

        let mut sequences = HashMap::with_capacity(records.len());

        for record in records {
            let raw_sequence = read_sequence_raw(&mut file, record)?;
            let name = String::from_utf8_lossy(record.name().as_ref()).into_owned();
            sequences.insert(name, raw_sequence);
        }

        debug!("Loaded {} contigs into memory (FAI-indexed)", sequences.len());
        Ok(Self { sequences: Arc::new(sequences) })
    }

    /// Load sequences using noodles sequential reading (fallback for non-indexed files).
    fn new_sequential(path: &Path) -> Result<Self> {
        let reader = noodles::fasta::io::reader::Builder.build_from_path(path)?;
        Self::read_all(reader)
    }

    /// Load the sequences of a gzip or BGZF FASTA, decompressing it from the
    /// start (`MultiGzDecoder` reads both: BGZF is a series of gzip members).
    fn new_compressed(path: &Path) -> Result<Self> {
        let file = File::open(path)
            .with_context(|| format!("Failed to open FASTA: {}", path.display()))?;
        let decoder = flate2::read::MultiGzDecoder::new(BufReader::new(file));
        let reader = noodles::fasta::io::Reader::new(BufReader::new(decoder));
        Self::read_all(reader)
            .with_context(|| format!("Failed to read compressed FASTA: {}", path.display()))
    }

    /// Reads every record of `reader` into memory.
    fn read_all<R: std::io::BufRead>(mut reader: noodles::fasta::io::Reader<R>) -> Result<Self> {
        let mut sequences = HashMap::new();
        for result in reader.records() {
            let record = result?;
            let name = std::str::from_utf8(record.name())?.to_string();
            let raw_sequence: Vec<u8> = record.sequence().as_ref().to_vec();
            sequences.insert(name, raw_sequence);
        }

        debug!("Loaded {} contigs into memory (sequential)", sequences.len());
        Ok(Self { sequences: Arc::new(sequences) })
    }

    /// Retrieves a subsequence from the reference genome.
    ///
    /// Since all sequences are preloaded into memory, this is an O(1) lookup
    /// followed by a slice copy.
    ///
    /// # Arguments
    /// * `chrom` - Chromosome/sequence name (e.g., "chr1", "1")
    /// * `start` - Start position (1-based, inclusive)
    /// * `end` - End position (1-based, inclusive)
    ///
    /// # Returns
    /// The requested subsequence as a vector of bytes (preserving original case)
    ///
    /// # Errors
    /// Returns an error if:
    /// - The chromosome is not found in the reference
    /// - The requested region exceeds the chromosome length
    ///
    /// # Examples
    /// ```no_run
    /// use fgumi_lib::reference::ReferenceReader;
    /// use noodles::core::Position;
    ///
    /// let reader = ReferenceReader::new("reference.fasta")?;
    ///
    /// // Fetch first 100 bases of chr1
    /// let seq = reader.fetch("chr1", Position::try_from(1)?, Position::try_from(100)?)?;
    /// assert_eq!(seq.len(), 100);
    /// # Ok::<(), anyhow::Error>(())
    /// ```
    pub fn fetch(&self, chrom: &str, start: Position, end: Position) -> Result<Vec<u8>> {
        Ok(self.fetch_slice(chrom, start, end)?.to_vec())
    }

    /// Borrowed-slice variant of [`Self::fetch`] that returns a reference into the
    /// in-memory sequence — same lookup semantics, no allocation.
    ///
    /// Prefer this over [`Self::fetch`] in hot loops (e.g., per-record reference
    /// access) where allocating a fresh `Vec<u8>` per call adds up to gigabytes of churn.
    ///
    /// # Errors
    ///
    /// Same as [`Self::fetch`]: returns an error if the chromosome is missing or the
    /// requested region exceeds the sequence length.
    pub fn fetch_slice(&self, chrom: &str, start: Position, end: Position) -> Result<&[u8]> {
        let sequence = self
            .sequences
            .get(chrom)
            .ok_or_else(|| FgumiError::ReferenceNotFound { ref_name: chrom.to_string() })?;
        let (start_idx, end_idx) = region_indices(chrom, start, end, sequence.len())?;
        Ok(&sequence[start_idx..end_idx])
    }

    /// Gets a single base from the reference at the specified position.
    ///
    /// This is a convenience method that delegates to `fetch()` with a single-base region.
    ///
    /// # Arguments
    /// * `chrom` - Chromosome/sequence name (e.g., "chr1", "1")
    /// * `pos` - Position (1-based)
    ///
    /// # Returns
    /// The base at the specified position (preserving original case from FASTA)
    ///
    /// # Errors
    /// Returns an error if:
    /// - The chromosome is not found in the reference
    /// - The position exceeds the chromosome length
    ///
    /// # Examples
    /// ```no_run
    /// use fgumi_lib::reference::ReferenceReader;
    /// use noodles::core::Position;
    ///
    /// let reader = ReferenceReader::new("reference.fasta")?;
    ///
    /// // Get the base at position 1000 of chr1
    /// let base = reader.base_at("chr1", Position::try_from(1000)?)?;
    /// assert!(matches!(base, b'A' | b'C' | b'G' | b'T' | b'N'));
    /// # Ok::<(), anyhow::Error>(())
    /// ```
    pub fn base_at(&self, chrom: &str, pos: Position) -> Result<u8> {
        let sequence = self
            .sequences
            .get(chrom)
            .ok_or_else(|| FgumiError::ReferenceNotFound { ref_name: chrom.to_string() })?;

        // Convert from 1-based to 0-based indexing
        let pos_idx = usize::from(pos) - 1;

        sequence.get(pos_idx).copied().ok_or_else(|| {
            FgumiError::InvalidParameter {
                parameter: "position".to_string(),
                reason: format!(
                    "Position {}:{} exceeds sequence length {}",
                    chrom,
                    pos,
                    sequence.len()
                ),
            }
            .into()
        })
    }
}

/// One contig's layout in the FASTA, from its `.fai` line.
#[derive(Debug, Clone, Copy)]
struct FaiLayout {
    /// Bases in the contig.
    length: usize,
    /// Byte offset of the contig's first base.
    offset: u64,
    /// Bases per full line.
    line_bases: usize,
    /// Bytes per full line, terminator included.
    line_width: usize,
}

impl FaiLayout {
    /// The file offset of the 0-based base `pos`.
    fn byte_offset(&self, pos: usize) -> u64 {
        let line = pos / self.line_bases;
        let column = pos % self.line_bases;
        self.offset + (line * self.line_width + column) as u64
    }
}

/// Why a FASTA cannot be read on demand by [`IndexedReferenceReader`].
#[derive(Debug, Clone, PartialEq, Eq)]
pub(crate) enum Unavailable {
    /// The FASTA is not a regular file (a FIFO, such as `--ref <(cat ref.fa)`),
    /// so it cannot be read at an offset.
    NotRegular(PathBuf),
    /// The FASTA is gzip or BGZF compressed.
    Compressed(PathBuf),
    /// No `.fai` index was found beside the FASTA.
    NoIndex(PathBuf),
}

impl Unavailable {
    /// What would make the FASTA readable on demand.
    pub(crate) fn advice(&self) -> &'static str {
        match self {
            Self::NotRegular(_) => "pass a regular, indexed FASTA file",
            Self::Compressed(_) => "decompress it and index it with `samtools faidx`",
            Self::NoIndex(_) => "index it with `samtools faidx`",
        }
    }
}

impl std::fmt::Display for Unavailable {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        match self {
            Self::NotRegular(path) => {
                write!(f, "{} is not a regular file but a pipe or other stream", path.display())
            }
            Self::Compressed(path) => write!(f, "{} is compressed (gzip or BGZF)", path.display()),
            Self::NoIndex(path) => write!(f, "no .fai index found at {}.fai", path.display()),
        }
    }
}

/// A reference reader that keeps only the `.fai` index in memory and reads each
/// requested region from the FASTA on demand.
///
/// The on-demand counterpart of [`ReferenceReader`]: it returns the same bytes
/// (case preserved, as written in the FASTA) and the same errors for an unknown
/// contig or an out-of-range region, but each fetch is a positional read
/// (`pread`) of just that region, so threads share one open file without a lock.
/// Use it where fetches are rare, so the disk reads cost less than holding the
/// whole FASTA (about 3 GB for a human reference).
///
/// The FASTA must be uncompressed and indexed (`samtools faidx`). Positional
/// reads use `pread`, so it is Unix-only, like the sort engine.
#[derive(Clone)]
pub(crate) struct IndexedReferenceReader {
    /// The FASTA, shared by every clone and read only positionally.
    file: Arc<File>,
    /// Each contig's layout, by name.
    contigs: Arc<HashMap<String, FaiLayout>>,
}

impl IndexedReferenceReader {
    /// Opens `path` and reads its `.fai` index, unless the FASTA cannot be read
    /// on demand: then the inner result says why ([`Unavailable`]), so a caller
    /// can fall back to holding it in memory. A regular FASTA is probed once; one
    /// that is not regular (a FIFO) is declined without being opened, so the
    /// fallback can still read it. No sequence is read until fetched.
    ///
    /// # Errors
    ///
    /// Returns an error if the FASTA does not exist or cannot be opened, or its
    /// index cannot be read.
    pub(crate) fn open(path: &Path) -> Result<std::result::Result<Self, Unavailable>> {
        if !path.exists() {
            return Err(FgumiError::InvalidFileFormat {
                file_type: "Reference FASTA".to_string(),
                path: path.display().to_string(),
                reason: "File does not exist".to_string(),
            }
            .into());
        }
        if !is_regular_file(path) {
            return Ok(Err(Unavailable::NotRegular(path.to_path_buf())));
        }
        if is_gzip(path) {
            return Ok(Err(Unavailable::Compressed(path.to_path_buf())));
        }
        let Some(fai_path) = find_fai_path(path) else {
            return Ok(Err(Unavailable::NoIndex(path.to_path_buf())));
        };
        let file = File::open(path)
            .with_context(|| format!("Failed to open FASTA: {}", path.display()))?;
        let index = fai::fs::read(&fai_path)
            .with_context(|| format!("Failed to read FAI index: {}", fai_path.display()))?;
        let records: &[fai::Record] = index.as_ref();
        let mut contigs = HashMap::with_capacity(records.len());
        for record in records {
            let name = String::from_utf8_lossy(record.name().as_ref()).into_owned();
            let layout = FaiLayout {
                length: usize::try_from(record.length())?,
                offset: record.position(),
                line_bases: usize::try_from(record.line_base_count().get())?,
                line_width: usize::try_from(record.line_width().get())?,
            };
            anyhow::ensure!(
                layout.line_width >= layout.line_bases,
                "FAI index {}: contig {name} has a line width ({}) below its bases per line ({})",
                fai_path.display(),
                layout.line_width,
                layout.line_bases
            );
            contigs.insert(name, layout);
        }
        debug!("Indexed {} contigs for on-demand reads: {}", contigs.len(), fai_path.display());
        Ok(Ok(Self { file: Arc::new(file), contigs: Arc::new(contigs) }))
    }

    /// [`Self::open`], failing with the advice of [`Unavailable`] when the FASTA
    /// cannot be read on demand.
    #[cfg(test)]
    fn new<P: AsRef<Path>>(path: P) -> Result<Self> {
        Self::open(path.as_ref())?.map_err(|why| {
            anyhow::anyhow!(
                "cannot read the reference FASTA on demand: {why}; {}, or load it into memory",
                why.advice()
            )
        })
    }

    /// Retrieves `chrom:start-end` (1-based, inclusive) from the FASTA, exactly as
    /// [`ReferenceReader::fetch`] does from memory.
    ///
    /// # Errors
    ///
    /// Returns the same errors as [`ReferenceReader::fetch`] for an unknown
    /// contig or an out-of-range region, and an I/O error if the read fails.
    pub(crate) fn fetch(&self, chrom: &str, start: Position, end: Position) -> Result<Vec<u8>> {
        let layout = self
            .contigs
            .get(chrom)
            .ok_or_else(|| FgumiError::ReferenceNotFound { ref_name: chrom.to_string() })?;
        let (start_idx, end_idx) = region_indices(chrom, start, end, layout.length)?;
        let first = layout.byte_offset(start_idx);
        let raw_len = usize::try_from(layout.byte_offset(end_idx - 1) + 1 - first)?;
        let mut raw = vec![0u8; raw_len];
        fgumi_bam_io::pread::read_at_exact(&*self.file, &mut raw, first).with_context(|| {
            format!("Failed to read {chrom}:{start}-{end} from the reference FASTA")
        })?;
        let wanted = end_idx - start_idx;
        if raw_len == wanted {
            return Ok(raw);
        }
        // Drop the line terminators in place: move the rest of the first line,
        // then whole lines, down over them. The write index never passes the
        // read index, and `copy_within` handles the overlap.
        let terminator = layout.line_width - layout.line_bases;
        let (mut read, mut written) = (0, 0);
        let mut line_left = layout.line_bases - start_idx % layout.line_bases;
        while written < wanted {
            let take = line_left.min(wanted - written);
            raw.copy_within(read..read + take, written);
            written += take;
            read += take + terminator;
            line_left = layout.line_bases;
        }
        raw.truncate(wanted);
        Ok(raw)
    }
}

impl fgumi_sam::ReferenceProvider for IndexedReferenceReader {
    fn fetch(&self, chrom: &str, start: Position, end: Position) -> Result<Vec<u8>> {
        self.fetch(chrom, start, end)
    }
}

/// A reference FASTA held in memory ([`ReferenceReader`]) or read on demand
/// through its index ([`IndexedReferenceReader`]), for a consumer that can use
/// either: realign, whose fetches are rare enough that it reads them from disk
/// by default.
#[derive(Clone)]
pub(crate) enum ReferenceSource {
    /// The whole FASTA in memory.
    InMemory(ReferenceReader),
    /// Regions read from the FASTA as they are fetched.
    Indexed(IndexedReferenceReader),
}

impl fgumi_sam::ReferenceProvider for ReferenceSource {
    fn fetch(&self, chrom: &str, start: Position, end: Position) -> Result<Vec<u8>> {
        match self {
            Self::InMemory(r) => r.fetch(chrom, start, end),
            Self::Indexed(r) => r.fetch(chrom, start, end),
        }
    }

    fn fetch_borrowed(
        &self,
        chrom: &str,
        start: Position,
        end: Position,
    ) -> Result<std::borrow::Cow<'_, [u8]>> {
        match self {
            Self::InMemory(r) => Ok(std::borrow::Cow::Borrowed(r.fetch_slice(chrom, start, end)?)),
            Self::Indexed(r) => Ok(std::borrow::Cow::Owned(r.fetch(chrom, start, end)?)),
        }
    }
}

impl fgumi_sam::ReferenceProvider for ReferenceReader {
    fn fetch(
        &self,
        chrom: &str,
        start: noodles::core::Position,
        end: noodles::core::Position,
    ) -> anyhow::Result<Vec<u8>> {
        self.fetch(chrom, start, end)
    }

    /// Returns a borrow into the in-memory sequence, avoiding the per-call `Vec`
    /// allocation `fetch` performs.
    fn fetch_borrowed(
        &self,
        chrom: &str,
        start: noodles::core::Position,
        end: noodles::core::Position,
    ) -> anyhow::Result<std::borrow::Cow<'_, [u8]>> {
        Ok(std::borrow::Cow::Borrowed(self.fetch_slice(chrom, start, end)?))
    }
}

#[cfg(feature = "consensus")]
impl fgumi_consensus::methylation::RefBaseProvider for ReferenceReader {
    fn base_at_0based(&self, chrom: &str, pos: u64) -> Option<u8> {
        let sequence = self.sequences.get(chrom)?;
        sequence.get(usize::try_from(pos).ok()?).copied()
    }

    fn sequence_for(&self, chrom: &str) -> Option<&[u8]> {
        self.sequences.get(chrom).map(Vec::as_slice)
    }
}

#[cfg(test)]
mod tests {
    /// Stages asking for the same FASTA share one loaded copy; another path
    /// loads its own.
    #[test]
    fn reference_cache_shares_one_copy_per_path() -> anyhow::Result<()> {
        use crate::sam::builder::create_test_fasta;
        let a = create_test_fasta(&[("chr1", "ACGT")])?;
        let b = create_test_fasta(&[("chr1", "ACGT")])?;
        let mut cache = ReferenceCache::default();
        let first = cache.load(a.path())?;
        let again = cache.load(a.path())?;
        let other = cache.load(b.path())?;
        assert!(Arc::ptr_eq(&first.sequences, &again.sequences), "same path: one copy");
        assert!(!Arc::ptr_eq(&first.sequences, &other.sequences), "another path loads anew");
        Ok(())
    }

    use super::*;
    use crate::sam::builder::create_default_test_fasta;

    /// `fetch_borrowed` must hand back a borrow, not an owned copy — that is the
    /// entire point of the override. A `Cow::Owned` here means the default trait
    /// implementation is being used and the per-record allocation is still happening.
    #[test]
    fn test_fetch_borrowed_borrows_and_matches_fetch() -> Result<()> {
        use fgumi_sam::ReferenceProvider;

        let fasta = create_default_test_fasta()?;
        let reader = ReferenceReader::new(fasta.path())?;
        let (start, end) = (Position::try_from(1)?, Position::try_from(4)?);

        let borrowed = ReferenceProvider::fetch_borrowed(&reader, "chr1", start, end)?;
        assert!(
            matches!(borrowed, std::borrow::Cow::Borrowed(_)),
            "expected a borrowed slice; an owned value means the allocation is still happening"
        );

        // Same bytes as the allocating path.
        let owned = ReferenceProvider::fetch(&reader, "chr1", start, end)?;
        assert_eq!(borrowed.as_ref(), owned.as_slice());

        // Errors propagate identically for a missing chromosome.
        assert!(ReferenceProvider::fetch_borrowed(&reader, "nope", start, end).is_err());
        Ok(())
    }

    /// The production call sites (`clip.rs`, `filter.rs`) pass `&ReferenceReader`
    /// into a `R: ReferenceProvider` generic, so they resolve through the `T: Deref`
    /// blanket impl rather than `ReferenceReader`'s own. If that blanket stops
    /// forwarding `fetch_borrowed`, it silently falls back to the trait default —
    /// which returns `Cow::Owned`, reinstating the per-record allocation with no
    /// output change and no other test failing. This pins the forwarding.
    #[test]
    fn test_fetch_borrowed_forwards_through_reference() -> Result<()> {
        use fgumi_sam::ReferenceProvider;

        /// Mirrors how `regenerate_alignment_tags_raw` takes its provider.
        fn borrows_via_generic<R: ReferenceProvider>(
            provider: R,
            chrom: &str,
            start: Position,
            end: Position,
        ) -> Result<bool> {
            Ok(matches!(provider.fetch_borrowed(chrom, start, end)?, std::borrow::Cow::Borrowed(_)))
        }

        let fasta = create_default_test_fasta()?;
        let reader = ReferenceReader::new(fasta.path())?;
        let (start, end) = (Position::try_from(1)?, Position::try_from(4)?);

        assert!(
            borrows_via_generic(&reader, "chr1", start, end)?,
            "&ReferenceReader must still borrow; an owned value means the Deref blanket \
             impl stopped forwarding fetch_borrowed and the allocation is back"
        );
        Ok(())
    }

    #[test]
    fn test_fetch_subsequence() -> Result<()> {
        let fasta = create_default_test_fasta()?;
        let reader = ReferenceReader::new(fasta.path())?;

        // Fetch from chr1
        let seq = reader.fetch("chr1", Position::try_from(1)?, Position::try_from(4)?)?;
        assert_eq!(seq, b"ACGT");

        // Fetch from chr2
        let seq = reader.fetch("chr2", Position::try_from(5)?, Position::try_from(8)?)?;
        assert_eq!(seq, b"CCCC");

        Ok(())
    }

    #[test]
    fn test_fetch_slice_returns_borrowed_bytes() -> Result<()> {
        let fasta = create_default_test_fasta()?;
        let reader = ReferenceReader::new(fasta.path())?;

        // Borrowed alternative to `fetch` — same bytes, no allocation.
        let slice: &[u8] =
            reader.fetch_slice("chr1", Position::try_from(1)?, Position::try_from(4)?)?;
        assert_eq!(slice, b"ACGT");

        // Same bytes as `fetch` would return.
        let owned = reader.fetch("chr1", Position::try_from(1)?, Position::try_from(4)?)?;
        assert_eq!(slice, owned.as_slice());

        // Bounds errors propagate just like `fetch`.
        assert!(
            reader
                .fetch_slice("chr1", Position::try_from(1)?, Position::try_from(10_000)?)
                .is_err()
        );
        assert!(
            reader.fetch_slice("nope", Position::try_from(1)?, Position::try_from(2)?).is_err()
        );
        // Inverted interval (start > end) must also error.
        assert!(
            reader.fetch_slice("chr1", Position::try_from(5)?, Position::try_from(4)?).is_err()
        );

        Ok(())
    }

    #[test]
    fn test_base_at() -> Result<()> {
        let fasta = create_default_test_fasta()?;
        let reader = ReferenceReader::new(fasta.path())?;

        assert_eq!(reader.base_at("chr1", Position::try_from(1)?)?, b'A');
        assert_eq!(reader.base_at("chr1", Position::try_from(2)?)?, b'C');
        assert_eq!(reader.base_at("chr2", Position::try_from(1)?)?, b'G');

        Ok(())
    }

    #[test]
    fn test_all_sequences_loaded() -> Result<()> {
        let fasta = create_default_test_fasta()?;
        let reader = ReferenceReader::new(fasta.path())?;

        // All sequences should be available immediately after construction
        let seq1 = reader.fetch("chr1", Position::try_from(1)?, Position::try_from(4)?)?;
        assert_eq!(seq1, b"ACGT");

        let seq2 = reader.fetch("chr2", Position::try_from(1)?, Position::try_from(4)?)?;
        assert_eq!(seq2, b"GGGG");

        // Fetching chr1 again should still work (all in memory)
        let seq1_again = reader.fetch("chr1", Position::try_from(1)?, Position::try_from(4)?)?;
        assert_eq!(seq1_again, b"ACGT");

        Ok(())
    }

    #[test]
    fn test_nonexistent_sequence() {
        let fasta = create_default_test_fasta().expect("creating test FASTA should succeed");
        let reader =
            ReferenceReader::new(fasta.path()).expect("creating reference reader should succeed");

        let result = reader.fetch(
            "chr999",
            Position::try_from(1).expect("position conversion should succeed"),
            Position::try_from(4).expect("position conversion should succeed"),
        );
        assert!(result.is_err());
    }

    #[test]
    fn test_out_of_bounds() {
        let fasta = create_default_test_fasta().expect("creating test FASTA should succeed");
        let reader =
            ReferenceReader::new(fasta.path()).expect("creating reference reader should succeed");

        // chr1 is only 12 bases long
        let result = reader.fetch(
            "chr1",
            Position::try_from(1).expect("position conversion should succeed"),
            Position::try_from(100).expect("position conversion should succeed"),
        );
        assert!(result.is_err());
    }

    #[test]
    fn test_reference_case_preserved() -> Result<()> {
        // Test that FASTA case is preserved (important for MD tag generation)
        // fgbio preserves case in reference sequences for proper MD tag output
        use crate::sam::builder::create_test_fasta;

        let file = create_test_fasta(&[("chr1", "AcGtNnAaCcGgTt")])?; // Mixed case sequence

        let reader = ReferenceReader::new(file.path())?;

        // Fetch the full sequence and verify case is preserved
        let seq = reader.fetch("chr1", Position::try_from(1)?, Position::try_from(14)?)?;
        assert_eq!(seq, b"AcGtNnAaCcGgTt");

        // Verify individual bases preserve case
        assert_eq!(reader.base_at("chr1", Position::try_from(1)?)?, b'A'); // uppercase A
        assert_eq!(reader.base_at("chr1", Position::try_from(2)?)?, b'c'); // lowercase c
        assert_eq!(reader.base_at("chr1", Position::try_from(3)?)?, b'G'); // uppercase G
        assert_eq!(reader.base_at("chr1", Position::try_from(4)?)?, b't'); // lowercase t
        assert_eq!(reader.base_at("chr1", Position::try_from(5)?)?, b'N'); // uppercase N
        assert_eq!(reader.base_at("chr1", Position::try_from(6)?)?, b'n'); // lowercase n

        Ok(())
    }

    #[test]
    fn test_n_bases_at_various_positions() -> Result<()> {
        use crate::sam::builder::create_test_fasta;

        // N at start, middle, and end
        let file = create_test_fasta(&[("chr1", "NACGTNACGTN")])?;
        let reader = ReferenceReader::new(file.path())?;

        assert_eq!(reader.base_at("chr1", Position::try_from(1)?)?, b'N');
        assert_eq!(reader.base_at("chr1", Position::try_from(6)?)?, b'N');
        assert_eq!(reader.base_at("chr1", Position::try_from(11)?)?, b'N');

        // Fetch range including N
        let seq = reader.fetch("chr1", Position::try_from(1)?, Position::try_from(6)?)?;
        assert_eq!(seq, b"NACGTN");

        Ok(())
    }

    #[test]
    fn test_all_n_sequence() -> Result<()> {
        use crate::sam::builder::create_test_fasta;

        let file = create_test_fasta(&[("chrN", "NNNNNNNNNN")])?;
        let reader = ReferenceReader::new(file.path())?;

        let seq = reader.fetch("chrN", Position::try_from(1)?, Position::try_from(10)?)?;
        assert_eq!(seq, b"NNNNNNNNNN");

        for i in 1..=10 {
            assert_eq!(reader.base_at("chrN", Position::try_from(i)?)?, b'N');
        }

        Ok(())
    }

    #[test]
    fn test_long_sequence_boundaries() -> Result<()> {
        use crate::sam::builder::create_test_fasta;

        // Create a 100-base sequence with markers at specific positions
        // We want to test fetching across the 32-base boundary (for bit-packed storage)
        let mut seq = String::new();

        // Positions 1-28: ACGT repeated (7 times)
        for _ in 0..7 {
            seq.push_str("ACGT");
        }
        // Positions 29-32: NNNN (straddles u64 boundary at position 32)
        seq.push_str("NNNN");
        // Positions 33-60: ACGT repeated (7 times)
        for _ in 0..7 {
            seq.push_str("ACGT");
        }
        // Positions 61-64: lowercase acgt (straddles second u64 boundary)
        seq.push_str("acgt");
        // Positions 65-100: ACGT repeated (9 times)
        for _ in 0..9 {
            seq.push_str("ACGT");
        }

        assert_eq!(seq.len(), 100);

        let file = create_test_fasta(&[("chr1", &seq)])?;
        let reader = ReferenceReader::new(file.path())?;

        // Verify positions around first marker (NNNN at 29-32)
        assert_eq!(reader.base_at("chr1", Position::try_from(28)?)?, b'T');
        assert_eq!(reader.base_at("chr1", Position::try_from(29)?)?, b'N');
        assert_eq!(reader.base_at("chr1", Position::try_from(32)?)?, b'N');
        assert_eq!(reader.base_at("chr1", Position::try_from(33)?)?, b'A');

        // Verify positions around second marker (acgt at 61-64)
        assert_eq!(reader.base_at("chr1", Position::try_from(60)?)?, b'T');
        assert_eq!(reader.base_at("chr1", Position::try_from(61)?)?, b'a');
        assert_eq!(reader.base_at("chr1", Position::try_from(64)?)?, b't');
        assert_eq!(reader.base_at("chr1", Position::try_from(65)?)?, b'A');

        // Fetch across boundaries
        let cross_first = reader.fetch("chr1", Position::try_from(27)?, Position::try_from(34)?)?;
        assert_eq!(cross_first, b"GTNNNNAC");

        let cross_second =
            reader.fetch("chr1", Position::try_from(59)?, Position::try_from(66)?)?;
        assert_eq!(cross_second, b"GTacgtAC");

        Ok(())
    }

    #[test]
    fn test_single_base_sequence() -> Result<()> {
        use crate::sam::builder::create_test_fasta;

        let file = create_test_fasta(&[("chr1", "A"), ("chr2", "N"), ("chr3", "g")])?;
        let reader = ReferenceReader::new(file.path())?;

        assert_eq!(reader.base_at("chr1", Position::try_from(1)?)?, b'A');
        assert_eq!(reader.base_at("chr2", Position::try_from(1)?)?, b'N');
        assert_eq!(reader.base_at("chr3", Position::try_from(1)?)?, b'g');

        let seq = reader.fetch("chr1", Position::try_from(1)?, Position::try_from(1)?)?;
        assert_eq!(seq, b"A");

        Ok(())
    }

    #[test]
    fn test_fetch_full_sequence() -> Result<()> {
        use crate::sam::builder::create_test_fasta;

        let original = "ACGTNacgtn";
        let file = create_test_fasta(&[("chr1", original)])?;
        let reader = ReferenceReader::new(file.path())?;

        let seq = reader.fetch("chr1", Position::try_from(1)?, Position::try_from(10)?)?;
        assert_eq!(seq, original.as_bytes());

        Ok(())
    }

    #[test]
    fn test_fetch_last_base() -> Result<()> {
        use crate::sam::builder::create_test_fasta;

        let file = create_test_fasta(&[("chr1", "ACGTN")])?;
        let reader = ReferenceReader::new(file.path())?;

        // Fetch last base
        let seq = reader.fetch("chr1", Position::try_from(5)?, Position::try_from(5)?)?;
        assert_eq!(seq, b"N");

        assert_eq!(reader.base_at("chr1", Position::try_from(5)?)?, b'N');

        Ok(())
    }

    #[test]
    fn test_multiple_chromosomes_isolation() -> Result<()> {
        use crate::sam::builder::create_test_fasta;

        let file = create_test_fasta(&[
            ("chr1", "AAAA"),
            ("chr2", "CCCC"),
            ("chr3", "GGGG"),
            ("chr4", "TTTT"),
        ])?;
        let reader = ReferenceReader::new(file.path())?;

        // Verify each chromosome has its own sequence
        assert_eq!(reader.fetch("chr1", Position::try_from(1)?, Position::try_from(4)?)?, b"AAAA");
        assert_eq!(reader.fetch("chr2", Position::try_from(1)?, Position::try_from(4)?)?, b"CCCC");
        assert_eq!(reader.fetch("chr3", Position::try_from(1)?, Position::try_from(4)?)?, b"GGGG");
        assert_eq!(reader.fetch("chr4", Position::try_from(1)?, Position::try_from(4)?)?, b"TTTT");

        Ok(())
    }

    #[test]
    fn test_mixed_case_all_bases() -> Result<()> {
        use crate::sam::builder::create_test_fasta;

        // Test all 10 possible base values: A, C, G, T, N (upper and lower)
        let file = create_test_fasta(&[("chr1", "ACGTNacgtn")])?;
        let reader = ReferenceReader::new(file.path())?;

        let expected = [b'A', b'C', b'G', b'T', b'N', b'a', b'c', b'g', b't', b'n'];
        for (i, &expected_base) in expected.iter().enumerate() {
            let pos = Position::try_from(i + 1)?;
            assert_eq!(
                reader.base_at("chr1", pos)?,
                expected_base,
                "Mismatch at position {}",
                i + 1
            );
        }

        Ok(())
    }

    #[test]
    fn test_runs_of_n_bases() -> Result<()> {
        use crate::sam::builder::create_test_fasta;

        // Simulate masked regions with runs of N
        let file = create_test_fasta(&[("chr1", "ACGTNNNNNNNNACGT")])?;
        let reader = ReferenceReader::new(file.path())?;

        // Fetch run of N's
        let n_run = reader.fetch("chr1", Position::try_from(5)?, Position::try_from(12)?)?;
        assert_eq!(n_run, b"NNNNNNNN");

        // Fetch across N boundary
        let across = reader.fetch("chr1", Position::try_from(3)?, Position::try_from(14)?)?;
        assert_eq!(across, b"GTNNNNNNNNAC");

        Ok(())
    }

    #[test]
    fn test_position_one_based() -> Result<()> {
        use crate::sam::builder::create_test_fasta;

        let file = create_test_fasta(&[("chr1", "ACGTN")])?;
        let reader = ReferenceReader::new(file.path())?;

        // Position 1 should be first base 'A', not second
        assert_eq!(reader.base_at("chr1", Position::try_from(1)?)?, b'A');
        assert_eq!(reader.base_at("chr1", Position::try_from(2)?)?, b'C');

        // fetch(1, 1) should return single base 'A'
        let seq = reader.fetch("chr1", Position::try_from(1)?, Position::try_from(1)?)?;
        assert_eq!(seq, b"A");

        // fetch(1, 2) should return "AC"
        let seq = reader.fetch("chr1", Position::try_from(1)?, Position::try_from(2)?)?;
        assert_eq!(seq, b"AC");

        Ok(())
    }

    #[test]
    fn test_find_dict_path_replacing_convention() -> Result<()> {
        // Test the fgbio/HTSJDK/Picard convention: ref.fa -> ref.dict
        let temp_dir = tempfile::tempdir()?;
        let fasta_path = temp_dir.path().join("ref.fa");
        let dict_path = temp_dir.path().join("ref.dict");

        // Create empty files
        std::fs::write(&fasta_path, "")?;
        std::fs::write(&dict_path, "")?;

        let found = find_dict_path(&fasta_path);
        assert!(found.is_some());
        assert_eq!(found.expect("should find dictionary path"), dict_path);

        Ok(())
    }

    #[test]
    fn test_find_dict_path_appending_convention() -> Result<()> {
        // Test the GATK convention: ref.fa -> ref.fa.dict
        let temp_dir = tempfile::tempdir()?;
        let fasta_path = temp_dir.path().join("ref.fa");
        let dict_path = temp_dir.path().join("ref.fa.dict");

        // Create empty files
        std::fs::write(&fasta_path, "")?;
        std::fs::write(&dict_path, "")?;

        let found = find_dict_path(&fasta_path);
        assert!(found.is_some());
        assert_eq!(found.expect("should find dictionary path"), dict_path);

        Ok(())
    }

    #[test]
    fn test_find_dict_path_prefers_replacing_convention() -> Result<()> {
        // When both exist, prefer the fgbio/HTSJDK convention (ref.dict)
        let temp_dir = tempfile::tempdir()?;
        let fasta_path = temp_dir.path().join("ref.fa");
        let dict_replacing = temp_dir.path().join("ref.dict");
        let dict_appending = temp_dir.path().join("ref.fa.dict");

        // Create all files
        std::fs::write(&fasta_path, "")?;
        std::fs::write(&dict_replacing, "")?;
        std::fs::write(&dict_appending, "")?;

        let found = find_dict_path(&fasta_path);
        assert!(found.is_some());
        // Should find ref.dict first (fgbio/HTSJDK convention)
        assert_eq!(found.expect("should find dictionary path"), dict_replacing);

        Ok(())
    }

    #[test]
    fn test_find_dict_path_not_found() {
        // When no dict file exists, return None
        let temp_dir = tempfile::tempdir().expect("creating temp file/dir should succeed");
        let fasta_path = temp_dir.path().join("ref.fa");

        // Create only the FASTA file, no dict
        std::fs::write(&fasta_path, "").expect("writing file should succeed");

        let found = find_dict_path(&fasta_path);
        assert!(found.is_none());
    }

    #[test]
    fn test_find_dict_path_fasta_extension() -> Result<()> {
        // Test with .fasta extension: ref.fasta -> ref.dict
        let temp_dir = tempfile::tempdir()?;
        let fasta_path = temp_dir.path().join("ref.fasta");
        let dict_path = temp_dir.path().join("ref.dict");

        std::fs::write(&fasta_path, "")?;
        std::fs::write(&dict_path, "")?;

        let found = find_dict_path(&fasta_path);
        assert!(found.is_some());
        assert_eq!(found.expect("should find dictionary path"), dict_path);

        Ok(())
    }

    #[test]
    fn test_find_dict_path_fasta_appending_convention() -> Result<()> {
        // Test with .fasta extension: ref.fasta -> ref.fasta.dict
        let temp_dir = tempfile::tempdir()?;
        let fasta_path = temp_dir.path().join("ref.fasta");
        let dict_path = temp_dir.path().join("ref.fasta.dict");

        std::fs::write(&fasta_path, "")?;
        std::fs::write(&dict_path, "")?;

        let found = find_dict_path(&fasta_path);
        assert!(found.is_some());
        assert_eq!(found.expect("should find dictionary path"), dict_path);

        Ok(())
    }

    /// [`IndexedReferenceReader`] against the in-memory [`ReferenceReader`] on a
    /// wrapped, soft-masked, multi-contig FASTA indexed by noodles.
    mod indexed {
        use super::super::*;
        use fgumi_sam::ReferenceProvider;
        use rstest::{fixture, rstest};
        use std::fmt::Write as _;
        use tempfile::TempDir;

        /// Bases per FASTA line.
        const LINE: usize = 10;

        /// The contigs: lowercase (soft-masked) runs, a contig an exact multiple
        /// of the line length, one shorter than a line, and one a base past a
        /// full line.
        const CONTIGS: [(&str, &str); 4] = [
            ("chrA", "ACGTacgtNNACGTTGCAacgtACGTGGCCttaaGGCCAAAAAcccccGGGGGtttttACGTA"),
            ("chrB", "acgtacgtacGTACGTACGTaaaaaccccc"),
            ("chrC", "ggCCa"),
            ("chrD", "ACGTACGTACg"),
        ];

        #[derive(Debug, Clone, Copy)]
        enum Eol {
            Lf,
            CrLf,
        }

        struct Fasta {
            _dir: TempDir,
            path: PathBuf,
        }

        /// Writes [`CONTIGS`] wrapped at [`LINE`] bases with `eol` line endings
        /// and indexes it with noodles (`ref.fa.fai`).
        fn write_fasta(eol: Eol) -> Fasta {
            let eol = match eol {
                Eol::Lf => "\n",
                Eol::CrLf => "\r\n",
            };
            let dir = tempfile::tempdir().expect("tempdir");
            let path = dir.path().join("ref.fa");
            let mut text = String::new();
            for (name, seq) in CONTIGS {
                write!(text, ">{name} description{eol}").expect("write");
                for line in seq.as_bytes().chunks(LINE) {
                    text.push_str(std::str::from_utf8(line).expect("ascii"));
                    text.push_str(eol);
                }
            }
            std::fs::write(&path, text).expect("write fasta");
            let index = noodles::fasta::fs::index(&path).expect("index");
            fai::fs::write(dir.path().join("ref.fa.fai"), &index).expect("write fai");
            Fasta { _dir: dir, path }
        }

        #[fixture]
        fn lf() -> Fasta {
            write_fasta(Eol::Lf)
        }

        fn pos(p: usize) -> Position {
            Position::try_from(p).expect("1-based position")
        }

        fn truth(chrom: &str) -> &'static str {
            CONTIGS.iter().find(|(n, _)| *n == chrom).expect("contig").1
        }

        /// Every region of every contig, under both line endings: the indexed
        /// reader returns the in-memory reader's bytes, which are the FASTA's own
        /// (lowercase kept). Covers every line-boundary straddle and contig end.
        #[rstest]
        fn every_region_matches_in_memory(#[values(Eol::Lf, Eol::CrLf)] eol: Eol) {
            let fasta = write_fasta(eol);
            let memory = ReferenceReader::new(&fasta.path).expect("in-memory");
            let indexed = IndexedReferenceReader::new(&fasta.path).expect("indexed");
            for (name, seq) in CONTIGS {
                for start in 1..=seq.len() {
                    for end in start..=seq.len() {
                        let want = &seq.as_bytes()[start - 1..end];
                        let got = indexed.fetch(name, pos(start), pos(end)).expect("fetch");
                        assert_eq!(got, want, "{name}:{start}-{end} ({eol:?})");
                        let mem = memory.fetch(name, pos(start), pos(end)).expect("fetch");
                        assert_eq!(got, mem, "{name}:{start}-{end} ({eol:?})");
                    }
                }
            }
        }

        /// Named regions, read through the `ReferenceProvider` both readers serve
        /// realign with.
        #[rstest]
        #[case::first_base("chrA", 1, 1)]
        #[case::whole_first_line("chrA", 1, 10)]
        #[case::last_base_of_a_line("chrA", 10, 10)]
        #[case::first_base_of_a_line("chrA", 11, 11)]
        #[case::straddles_one_boundary("chrA", 8, 13)]
        #[case::straddles_several_boundaries("chrA", 9, 42)]
        #[case::soft_masked_run("chrA", 5, 8)]
        #[case::contig_end("chrA", 60, 63)]
        #[case::whole_contig("chrA", 1, 63)]
        #[case::exact_multiple_end("chrB", 21, 30)]
        #[case::exact_multiple_whole("chrB", 1, 30)]
        #[case::shorter_than_a_line("chrC", 2, 5)]
        #[case::one_past_a_full_line("chrD", 10, 11)]
        #[case::last_contig_last_base("chrD", 11, 11)]
        fn region_matches(
            lf: Fasta,
            #[case] chrom: &str,
            #[case] start: usize,
            #[case] end: usize,
        ) {
            let memory = ReferenceReader::new(&lf.path).expect("in-memory");
            let indexed = IndexedReferenceReader::new(&lf.path).expect("indexed");
            let want = &truth(chrom).as_bytes()[start - 1..end];
            let got = ReferenceProvider::fetch(&indexed, chrom, pos(start), pos(end)).expect("ok");
            assert_eq!(got, want);
            let mem = ReferenceProvider::fetch(&memory, chrom, pos(start), pos(end)).expect("ok");
            assert_eq!(got, mem);
        }

        /// A few thousand pseudo-random regions, from threads sharing one reader
        /// (positional reads need no lock).
        #[rstest]
        fn random_regions_from_threads_match(lf: Fasta) {
            let memory = ReferenceReader::new(&lf.path).expect("in-memory");
            let indexed = IndexedReferenceReader::new(&lf.path).expect("indexed");
            std::thread::scope(|scope| {
                for seed in 1..=4u64 {
                    let (memory, indexed) = (&memory, &indexed);
                    scope.spawn(move || {
                        let mut state = seed;
                        let mut next = |n: usize| {
                            state ^= state << 13;
                            state ^= state >> 7;
                            state ^= state << 17;
                            usize::try_from(state % n as u64).expect("fits")
                        };
                        for _ in 0..1000 {
                            let (name, seq) = CONTIGS[next(CONTIGS.len())];
                            let start = 1 + next(seq.len());
                            let end = start + next(seq.len() - start + 1);
                            let got = indexed.fetch(name, pos(start), pos(end)).expect("fetch");
                            let mem = memory.fetch(name, pos(start), pos(end)).expect("fetch");
                            assert_eq!(got, mem, "{name}:{start}-{end}");
                        }
                    });
                }
            });
        }

        /// Unknown contigs and out-of-range regions fail with the in-memory
        /// reader's message.
        #[rstest]
        #[case::unknown_contig("chrZ", 1, 1, "chrZ")]
        #[case::end_past_contig("chrA", 60, 64, "exceeds sequence length 63")]
        #[case::start_past_contig("chrC", 6, 6, "exceeds sequence length 5")]
        #[case::start_after_end("chrA", 5, 4, "exceeds sequence length 63")]
        fn errors_match_in_memory(
            lf: Fasta,
            #[case] chrom: &str,
            #[case] start: usize,
            #[case] end: usize,
            #[case] message: &str,
        ) {
            let memory = ReferenceReader::new(&lf.path).expect("in-memory");
            let indexed = IndexedReferenceReader::new(&lf.path).expect("indexed");
            let got = indexed.fetch(chrom, pos(start), pos(end)).expect_err("indexed fails");
            let mem = memory.fetch(chrom, pos(start), pos(end)).expect_err("in-memory fails");
            assert_eq!(format!("{got:#}"), format!("{mem:#}"));
            assert!(format!("{got:#}").contains(message), "{got:#}");
        }

        /// Opening needs an existing, uncompressed, `.fai`-indexed FASTA, and
        /// says which of these is missing, advising only what that one needs.
        #[rstest]
        #[case::missing_fasta("missing", "File does not exist", None)]
        #[case::missing_fai(
            "no-fai",
            "no .fai index found at",
            Some("index it with `samtools faidx`")
        )]
        #[case::gzip("gzip", "is compressed (gzip or BGZF)", Some("decompress it and index it"))]
        fn open_errors(#[case] form: &str, #[case] message: &str, #[case] advice: Option<&str>) {
            use std::io::Write;
            let dir = tempfile::tempdir().expect("tempdir");
            let path = dir.path().join("ref.fa");
            match form {
                "missing" => {}
                "no-fai" => std::fs::write(&path, ">chr1\nACGT\n").expect("write"),
                "gzip" => {
                    let file = std::fs::File::create(&path).expect("create");
                    let mut gz = flate2::write::GzEncoder::new(file, flate2::Compression::fast());
                    gz.write_all(b">chr1\nACGT\n").expect("write");
                    gz.finish().expect("finish");
                    std::fs::write(dir.path().join("ref.fa.fai"), "chr1\t4\t6\t4\t5\n")
                        .expect("fai");
                }
                other => unreachable!("{other}"),
            }
            let err = IndexedReferenceReader::new(&path).err().expect("open fails");
            let rendered = format!("{err:#}");
            assert!(rendered.contains(message), "{rendered}");
            if let Some(advice) = advice {
                assert!(rendered.contains(advice), "{rendered}");
            }
            if form == "no-fai" {
                assert!(rendered.contains(&format!("{}.fai", path.display())), "{rendered}");
                assert!(!rendered.contains("decompress"), "nothing to decompress: {rendered}");
            }
            // The reason realign logs when it falls back to holding the FASTA in
            // memory; a missing FASTA is an error instead.
            match IndexedReferenceReader::open(&path) {
                Err(e) => {
                    assert_eq!(form, "missing", "{e:#}");
                    assert!(format!("{e:#}").contains(message), "{e:#}");
                }
                Ok(Err(why)) => assert!(why.to_string().contains(message), "{why}"),
                Ok(Ok(_)) => panic!("{form} must not open"),
            }
        }

        /// An uncompressed, indexed FASTA can be read on demand.
        #[rstest]
        fn indexed_fasta_opens(lf: Fasta) {
            assert!(matches!(IndexedReferenceReader::open(&lf.path), Ok(Ok(_))));
        }

        #[derive(Debug, Clone, Copy)]
        enum Compression {
            Bgzf,
            Gzip,
        }

        /// A gzip or BGZF FASTA is decompressed into memory, never read as raw
        /// bytes: not at the uncompressed offsets of a `.fai` beside it (which
        /// `samtools faidx` writes for a BGZF FASTA), and not by its name.
        #[rstest]
        #[case::bgzf_with_fai(Compression::Bgzf, "ref.fa.gz", true)]
        #[case::bgzf_without_fai(Compression::Bgzf, "ref.fa.gz", false)]
        #[case::gzip_with_fai(Compression::Gzip, "ref.fa.gz", true)]
        #[case::gzip_named_like_plain_text(Compression::Gzip, "ref.fa", false)]
        fn compressed_fasta_loads_into_memory(
            lf: Fasta,
            #[case] compression: Compression,
            #[case] name: &str,
            #[case] with_fai: bool,
        ) {
            use std::io::Write;
            let text = std::fs::read(&lf.path).expect("plain FASTA");
            let dir = tempfile::tempdir().expect("tempdir");
            let path = dir.path().join(name);
            let file = std::fs::File::create(&path).expect("create");
            match compression {
                Compression::Bgzf => {
                    let mut writer = noodles::bgzf::io::Writer::new(file);
                    writer.write_all(&text).expect("write");
                    writer.finish().expect("finish");
                }
                Compression::Gzip => {
                    let mut gz = flate2::write::GzEncoder::new(file, flate2::Compression::fast());
                    gz.write_all(&text).expect("write");
                    gz.finish().expect("finish");
                }
            }
            if with_fai {
                // The plain FASTA's index: offsets into the uncompressed text.
                let fai = format!("{}.fai", path.display());
                std::fs::copy(format!("{}.fai", lf.path.display()), fai).expect("fai");
            }
            let reader = ReferenceReader::new(&path).expect("compressed FASTA loads");
            for (contig, bases) in CONTIGS {
                let got = reader.fetch(contig, pos(1), pos(bases.len())).expect("fetch");
                assert_eq!(got, bases.as_bytes(), "{contig} ({compression:?}, {name})");
            }
        }

        /// A FASTA given as a FIFO (`--ref <(cat ref.fa)`) is read once, from the
        /// start: nothing probes it first, and an index beside it is not used to
        /// seek in it. The on-demand reader declines it without opening it, so
        /// realign's fallback to memory reads the whole stream.
        #[cfg(unix)]
        #[rstest]
        #[case::without_fai(false)]
        #[case::with_fai(true)]
        fn fifo_fasta_loads_into_memory(lf: Fasta, #[case] with_fai: bool) {
            use std::io::Write;
            let text = std::fs::read(&lf.path).expect("plain FASTA");
            let dir = tempfile::tempdir().expect("tempdir");
            let path = dir.path().join("ref.fa");
            let status =
                std::process::Command::new("mkfifo").arg(&path).status().expect("run mkfifo");
            assert!(status.success(), "mkfifo {} failed: {status}", path.display());
            if with_fai {
                std::fs::copy(format!("{}.fai", lf.path.display()), dir.path().join("ref.fa.fai"))
                    .expect("fai");
            }
            // One writer, as from a process substitution: the FASTA can be opened
            // and read only once.
            let writer = {
                let path = path.clone();
                std::thread::spawn(move || {
                    if let Ok(mut fifo) = std::fs::OpenOptions::new().write(true).open(&path) {
                        let _ = fifo.write_all(&text);
                    }
                })
            };
            // A second open of the FIFO would block forever, so load on a thread
            // and fail on a timeout instead of hanging.
            let (tx, rx) = std::sync::mpsc::channel();
            {
                let path = path.clone();
                std::thread::spawn(move || {
                    let declined = IndexedReferenceReader::open(&path).map(Result::err);
                    let loaded = ReferenceReader::new(&path);
                    let _ = tx.send((declined, loaded));
                });
            }
            let (declined, loaded) = rx
                .recv_timeout(std::time::Duration::from_secs(30))
                .expect("the FIFO is opened once and read to its end");
            let why = declined.expect("open").expect("the on-demand reader declines a FIFO");
            assert!(why.to_string().contains("is not a regular file"), "{why}");
            let reader = loaded.expect("the FIFO loads into memory");
            for (contig, bases) in CONTIGS {
                let got = reader.fetch(contig, pos(1), pos(bases.len())).expect("fetch");
                assert_eq!(got, bases.as_bytes(), "{contig}");
            }
            writer.join().expect("writer");
        }

        /// `ReferenceSource` borrows from the in-memory copy and reads the
        /// indexed one, with the same bytes.
        #[rstest]
        fn reference_source_dispatches(lf: Fasta) {
            let memory =
                ReferenceSource::InMemory(ReferenceReader::new(&lf.path).expect("in-memory"));
            let indexed =
                ReferenceSource::Indexed(IndexedReferenceReader::new(&lf.path).expect("indexed"));
            let (start, end) = (pos(8), pos(13));
            let borrowed = memory.fetch_borrowed("chrA", start, end).expect("ok");
            assert!(matches!(borrowed, std::borrow::Cow::Borrowed(_)));
            let read = indexed.fetch_borrowed("chrA", start, end).expect("ok");
            assert_eq!(borrowed, read);
            assert_eq!(read.as_ref(), b"tNNACG");
        }
    }
}
