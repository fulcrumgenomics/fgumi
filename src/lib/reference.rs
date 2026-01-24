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
//! # Memory Efficiency
//!
//! Reference sequences are stored using 2-bit encoding with separate tracking for N bases
//! and case information. This reduces memory usage by approximately 2.7x compared to
//! storing raw bytes (from ~3GB to ~1.1GB for a human genome).
//!
//! # Future improvement
//!
//! The custom FAI-based raw-byte reading (`read_sequence_raw`) could be replaced with
//! noodles' built-in indexed reader once <https://github.com/zaeleus/noodles/pull/365>
//! is merged and released, which adds the same optimization to noodles.

use crate::compressed_seq::CompressedSequence;
use crate::errors::FgumiError;
use anyhow::{Context, Result};
use log::debug;
use noodles::core::Position;
use noodles::fasta::fai;
use std::collections::HashMap;
use std::fs::File;
use std::io::{Read, Seek, SeekFrom};
use std::path::{Path, PathBuf};
use std::sync::Arc;

/// Read a sequence from a FASTA file using FAI index metadata.
/// Uses raw byte reading with mathematical newline handling (htsjdk-style).
///
/// Optimized: Reads entire sequence bytes in one syscall, then strips newlines in memory.
/// Fast path: If sequence fits in a single line, skip newline stripping entirely.
fn read_sequence_raw(file: &mut File, record: &fai::Record) -> Result<Vec<u8>> {
    let line_bases = record.line_bases() as usize;
    let line_width = record.line_width() as usize;
    let seq_len = record.length() as usize;
    let offset = record.offset();

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

/// Find FAI index path for a FASTA file.
fn find_fai_path(fasta_path: &Path) -> Option<PathBuf> {
    // Try .fai extension
    let fai_path = fasta_path.with_extension("fa.fai");
    if fai_path.exists() {
        return Some(fai_path);
    }

    // Try appending .fai to full path
    let fai_path = PathBuf::from(format!("{}.fai", fasta_path.display()));
    if fai_path.exists() {
        return Some(fai_path);
    }

    None
}

/// A thread-safe reference genome reader with all sequences preloaded into memory.
///
/// This reader loads the entire FASTA file into memory at construction time,
/// providing O(1) lookup performance for sequence fetches. This approach matches
/// fgbio's `nmUqMdTagRegeneratingWriter` which reads all contigs into a Map upfront.
///
/// For a typical human reference (e.g., hs38DH at ~3GB), this uses approximately
/// 1.1GB of memory (using compressed 2-bit encoding) and provides dramatic speedup
/// for operations that access many reads across different chromosomes.
#[derive(Clone)]
pub struct ReferenceReader {
    /// All sequences loaded into memory using compressed storage, keyed by sequence name
    sequences: Arc<HashMap<String, CompressedSequence>>,
}

impl ReferenceReader {
    /// Creates a new reference reader, loading all sequences into memory.
    ///
    /// This reads the entire FASTA file into memory at construction time.
    /// For a typical human reference (~3GB), this takes a few seconds but
    /// provides O(1) lookup performance for all subsequent fetches.
    ///
    /// # Arguments
    /// * `path` - Path to the reference FASTA file (may be gzipped)
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
            let compressed = CompressedSequence::from_bytes(&raw_sequence);
            let name = String::from_utf8_lossy(record.name().as_ref()).into_owned();
            sequences.insert(name, compressed);
        }

        debug!("Loaded {} contigs into memory (FAI-indexed, compressed)", sequences.len());
        Ok(Self { sequences: Arc::new(sequences) })
    }

    /// Load sequences using noodles sequential reading (fallback for non-indexed files).
    fn new_sequential(path: &Path) -> Result<Self> {
        use noodles::fasta;

        let mut sequences = HashMap::new();
        let mut reader = fasta::io::reader::Builder.build_from_path(path)?;

        for result in reader.records() {
            let record = result?;
            let name = std::str::from_utf8(record.name())?.to_string();
            let raw_sequence: &[u8] = record.sequence().as_ref();
            let compressed = CompressedSequence::from_bytes(raw_sequence);
            sequences.insert(name, compressed);
        }

        debug!("Loaded {} contigs into memory (sequential, compressed)", sequences.len());
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
        let sequence = self
            .sequences
            .get(chrom)
            .ok_or_else(|| FgumiError::ReferenceNotFound { ref_name: chrom.to_string() })?;

        // Convert from 1-based inclusive to 0-based [start, end) indexing
        let start_idx = usize::from(start) - 1;
        let end_idx = usize::from(end);

        sequence.fetch(start_idx, end_idx).ok_or_else(|| {
            FgumiError::InvalidParameter {
                parameter: "region".to_string(),
                reason: format!(
                    "Requested region {}:{}-{} exceeds sequence length {}",
                    chrom,
                    start,
                    end,
                    sequence.len()
                ),
            }
            .into()
        })
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

        sequence.base_at(pos_idx).ok_or_else(|| {
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

#[cfg(test)]
mod tests {
    use super::*;
    use crate::sam::builder::create_default_test_fasta;

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
        let fasta = create_default_test_fasta().unwrap();
        let reader = ReferenceReader::new(fasta.path()).unwrap();

        let result =
            reader.fetch("chr999", Position::try_from(1).unwrap(), Position::try_from(4).unwrap());
        assert!(result.is_err());
    }

    #[test]
    fn test_out_of_bounds() {
        let fasta = create_default_test_fasta().unwrap();
        let reader = ReferenceReader::new(fasta.path()).unwrap();

        // chr1 is only 12 bases long
        let result =
            reader.fetch("chr1", Position::try_from(1).unwrap(), Position::try_from(100).unwrap());
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
}
