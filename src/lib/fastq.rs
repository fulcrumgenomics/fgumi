//! FASTQ file parsing and read structure handling.
//!
//! This module provides functionality for reading and parsing FASTQ files with support for
//! complex read structures. It can handle reads that contain multiple segments of different
//! types (template, UMI, sample barcodes, cell barcodes) and provides utilities for
//! extracting and manipulating these segments.
//!
//! # Read Structures
//!
//! Read structures describe the layout of bases within a FASTQ read. For example, a read
//! might be structured as "8M143T" meaning 8 bases of molecular barcode (UMI) followed by
//! 143 bases of template sequence.
//!
//! # Example
//!
//! ```rust,ignore
//! use read_structure::ReadStructure;
//! use std::fs::File;
//! use std::io::BufReader;
//! use seq_io::fastq::Reader;
//!
//! let rs = ReadStructure::from_str("8M143T").unwrap();
//! let file = File::open("reads.fq").unwrap();
//! let reader = Reader::new(BufReader::new(file));
//! let mut iter = ReadSetIterator::new(rs, reader, vec![]);
//!
//! for read_set in iter {
//!     // Process each read set...
//! }
//! ```

use anyhow::{Result, anyhow};
use read_structure::ReadStructure;
use read_structure::SegmentType;
use seq_io::fastq::Reader as FastqReader;
use seq_io::fastq::Record;
use std::fmt::Display;
use std::io::BufRead;
use std::iter::Filter;
use std::slice::Iter;
use std::str::FromStr;

/// Type alias for segment type iter functions, which iterate over the segments of a `FastqSet`
/// filtering for a specific type.
type SegmentIter<'a> = Filter<Iter<'a, FastqSegment>, fn(&&FastqSegment) -> bool>;

/// A segment of a FASTQ record representing bases and qualities of a specific type.
///
/// FASTQ records can be divided into segments based on their biological function
/// (template, UMI, barcode, etc.). Each segment contains the bases, quality scores,
/// and the segment type.
#[derive(PartialEq, Eq, Debug, Clone)]
pub struct FastqSegment {
    /// The nucleotide bases of this segment (A, C, G, T, N)
    pub seq: Vec<u8>,
    /// Phred-scaled quality scores (typically ASCII-encoded as Phred+33)
    pub quals: Vec<u8>,
    /// The biological/functional type of this segment
    pub segment_type: SegmentType,
}

////////////////////////////////////////////////////////////////////////////////
// FastqSet and its impls
////////////////////////////////////////////////////////////////////////////////

/// Reasons why a read might be skipped during processing.
///
/// When processing FASTQ reads with read structures, some reads may not meet
/// the minimum requirements and need to be skipped. This enum categorizes the
/// reasons for skipping.
#[derive(Eq, Hash, PartialEq, Debug, Clone, Copy)]
pub enum SkipReason {
    /// The read had too few bases for the required segments in the read structure
    TooFewBases,
}

impl Display for SkipReason {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        match self {
            SkipReason::TooFewBases => write!(f, "Too few bases"),
        }
    }
}

impl FromStr for SkipReason {
    type Err = anyhow::Error;
    fn from_str(string: &str) -> Result<Self, Self::Err> {
        match string {
            "too few bases" | "too-few-bases" | "toofewbases" => Ok(SkipReason::TooFewBases),
            _ => Err(anyhow!(
                "Invalid skip reason: '{string}' (valid values: 'too-few-bases', 'toofewbases')"
            )),
        }
    }
}

/// A FASTQ record parsed into its component segments according to a read structure.
///
/// When reading FASTQ files with complex structures (e.g., containing UMIs, barcodes,
/// and template sequences), this structure holds the parsed segments along with the
/// original header. Each segment is classified by type (template, molecular barcode,
/// sample barcode, etc.).
///
/// # Fields
///
/// * `header` - The original FASTQ header line (without the '@' prefix)
/// * `segments` - Vector of parsed segments in order
/// * `skip_reason` - If Some, indicates why this read should be skipped during processing
#[derive(PartialEq, Debug, Clone)]
pub struct FastqSet {
    /// The FASTQ header line (without '@' prefix)
    pub header: Vec<u8>,
    /// Ordered list of parsed segments from this read
    pub segments: Vec<FastqSegment>,
    /// Reason for skipping this read, or None if the read should be processed
    pub skip_reason: Option<SkipReason>,
}

impl FastqSet {
    /// Returns an iterator over template segments.
    ///
    /// Template segments contain the actual biological sequence to be analyzed
    /// (as opposed to barcodes or UMIs).
    ///
    /// # Returns
    ///
    /// Iterator over references to template segments in this read set
    pub fn template_segments(&self) -> SegmentIter<'_> {
        self.segments.iter().filter(|s| s.segment_type == SegmentType::Template)
    }

    /// Returns an iterator over sample barcode segments.
    ///
    /// Sample barcodes identify which sample/library a read belongs to in
    /// multiplexed sequencing runs.
    ///
    /// # Returns
    ///
    /// Iterator over references to sample barcode segments
    pub fn sample_barcode_segments(&self) -> SegmentIter<'_> {
        self.segments.iter().filter(|s| s.segment_type == SegmentType::SampleBarcode)
    }

    /// Returns an iterator over molecular barcode (UMI) segments.
    ///
    /// Molecular barcodes (UMIs) uniquely tag individual molecules for
    /// duplicate detection and quantification.
    ///
    /// # Returns
    ///
    /// Iterator over references to molecular barcode segments
    pub fn molecular_barcode_segments(&self) -> SegmentIter<'_> {
        self.segments.iter().filter(|s| s.segment_type == SegmentType::MolecularBarcode)
    }

    /// Returns an iterator over cell barcode segments.
    ///
    /// Cell barcodes identify which cell a read originated from in
    /// single-cell sequencing experiments.
    ///
    /// # Returns
    ///
    /// Iterator over references to cell barcode segments
    pub fn cell_barcode_segments(&self) -> SegmentIter<'_> {
        self.segments.iter().filter(|s| s.segment_type == SegmentType::CellularBarcode)
    }

    /// Combines multiple `FastqSet` instances into a single set.
    ///
    /// Merges the segments from multiple read sets, preserving the header from
    /// the first set. Useful for combining segments from multiple FASTQ files
    /// that represent different parts of the same read (e.g., R1, R2, I1, I2).
    ///
    /// # Arguments
    ///
    /// * `readsets` - Vector of read sets to combine
    ///
    /// # Returns
    ///
    /// A new `FastqSet` with all segments combined
    ///
    /// # Panics
    ///
    /// Panics if the input vector is empty or contains no segments
    #[must_use]
    pub fn combine_readsets(readsets: Vec<Self>) -> Self {
        let total_segments: usize = readsets.iter().map(|s| s.segments.len()).sum();
        assert!(total_segments > 0, "Cannot call combine readsets on an empty vec!");

        let mut readset_iter = readsets.into_iter();
        let mut first = readset_iter
            .next()
            .expect("readset_iter is non-empty because total_segments > 0 was checked");
        first.segments.reserve_exact(total_segments - first.segments.len());

        for next_readset in readset_iter {
            first.segments.extend(next_readset.segments);
        }

        first
    }

    /// Creates a `FastqSet` from raw FASTQ record data and a read structure.
    ///
    /// This method applies the read structure to segment the sequence and quality
    /// data, producing a `FastqSet` with appropriately typed segments.
    ///
    /// # Arguments
    ///
    /// * `header` - The FASTQ header/name (without '@' prefix)
    /// * `sequence` - The full sequence bases
    /// * `quality` - The full quality scores
    /// * `read_structure` - The read structure describing segment layout
    /// * `skip_reasons` - List of reasons to skip reads gracefully instead of failing
    ///
    /// # Returns
    ///
    /// A `FastqSet` with segments parsed according to the read structure.
    /// If the read doesn't have enough bases and `TooFewBases` is in `skip_reasons`,
    /// returns a `FastqSet` with `skip_reason` set.
    ///
    /// # Panics
    ///
    /// Panics if the read has too few bases and `TooFewBases` is not in `skip_reasons`.
    #[must_use]
    pub fn from_record_with_structure(
        header: &[u8],
        sequence: &[u8],
        quality: &[u8],
        read_structure: &ReadStructure,
        skip_reasons: &[SkipReason],
    ) -> Self {
        let mut segments = Vec::with_capacity(read_structure.number_of_segments());

        // Check minimum length required by the read structure
        let min_len: usize = read_structure.iter().map(|s| s.length().unwrap_or(0)).sum();
        if sequence.len() < min_len {
            if skip_reasons.contains(&SkipReason::TooFewBases) {
                return Self {
                    header: header.to_vec(),
                    segments: vec![],
                    skip_reason: Some(SkipReason::TooFewBases),
                };
            }
            let read_name_str = String::from_utf8_lossy(header);
            panic!(
                "Read {} had too few bases to demux {} vs. {} needed in read structure {}.",
                read_name_str,
                sequence.len(),
                min_len,
                read_structure
            );
        }

        // Extract segments according to the read structure
        for (segment_index, read_segment) in read_structure.iter().enumerate() {
            let (seq, quals) =
                read_segment.extract_bases_and_quals(sequence, quality).unwrap_or_else(|e| {
                    let read_name_str = String::from_utf8_lossy(header);
                    panic!(
                        "Error extracting bases (len: {}) or quals (len: {}) for the {}th \
                         segment ({}) in read structure ({}) from record {}; {}",
                        sequence.len(),
                        quality.len(),
                        segment_index,
                        read_segment,
                        read_structure,
                        read_name_str,
                        e
                    );
                });

            segments.push(FastqSegment {
                seq: seq.to_vec(),
                quals: quals.to_vec(),
                segment_type: read_segment.kind,
            });
        }

        Self { header: header.to_vec(), segments, skip_reason: None }
    }
}

////////////////////////////////////////////////////////////////////////////////
// ReadSetIterator and its impls
////////////////////////////////////////////////////////////////////////////////

/// Iterator for parsing FASTQ files according to a read structure.
///
/// This iterator reads FASTQ records and parses them into `FastqSet` objects based on
/// the provided read structure. It handles validation and can optionally skip reads that
/// don't meet requirements (e.g., too few bases) rather than panicking.
///
/// # Example
///
/// ```rust,ignore
/// let read_structure = ReadStructure::from_str("8M143T")?;
/// let file = File::open("reads.fq")?;
/// let reader = Reader::new(BufReader::new(file));
/// let mut iterator = ReadSetIterator::new(read_structure, reader, vec![]);
///
/// for read_set in iterator {
///     // Process each parsed read set
/// }
/// ```
pub struct ReadSetIterator {
    /// Read structure describing the layout of bases in each read
    read_structure: ReadStructure,
    /// FASTQ file reader
    source: FastqReader<Box<dyn BufRead + Send>>,
    /// Reasons to skip reads instead of panicking (e.g., too few bases)
    skip_reasons: Vec<SkipReason>,
}

impl Iterator for ReadSetIterator {
    type Item = FastqSet;

    fn next(&mut self) -> Option<Self::Item> {
        let rec = self.source.next()?;
        let record = match rec {
            Ok(record) => record,
            Err(e) => {
                panic!("Error parsing FASTQ record: {e}");
            }
        };
        Some(FastqSet::from_record_with_structure(
            record.head(),
            record.seq(),
            record.qual(),
            &self.read_structure,
            &self.skip_reasons,
        ))
    }
}

impl ReadSetIterator {
    /// Creates a new iterator for parsing FASTQ records with a read structure.
    ///
    /// # Arguments
    ///
    /// * `read_structure` - Read structure describing the segment layout
    /// * `source` - FASTQ reader for input
    /// * `skip_reasons` - List of reasons to skip reads gracefully instead of panicking.
    ///   If a read fails for a reason not in this list, the iterator will panic.
    ///
    /// # Returns
    ///
    /// A new `ReadSetIterator` instance
    ///
    /// # Example
    ///
    /// ```rust,ignore
    /// // Skip reads that are too short instead of panicking
    /// let iterator = ReadSetIterator::new(
    ///     read_structure,
    ///     reader,
    ///     vec![SkipReason::TooFewBases]
    /// );
    /// ```
    #[must_use]
    pub fn new(
        read_structure: ReadStructure,
        source: FastqReader<Box<dyn BufRead + Send>>,
        skip_reasons: Vec<SkipReason>,
    ) -> Self {
        Self { read_structure, source, skip_reasons }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    fn read_set(segments: Vec<FastqSegment>) -> FastqSet {
        FastqSet { header: "NOT_IMPORTANT".as_bytes().to_owned(), segments, skip_reason: None }
    }

    fn seg(bases: &[u8], segment_type: SegmentType) -> FastqSegment {
        let quals = vec![b'#'; bases.len()];
        FastqSegment { seq: bases.to_vec(), quals, segment_type }
    }

    #[test]
    fn test_template_segments() {
        let segments = vec![
            seg("AGCT".as_bytes(), SegmentType::SampleBarcode),
            seg("GATA".as_bytes(), SegmentType::Template),
            seg("CAC".as_bytes(), SegmentType::Template),
            seg("GACCCC".as_bytes(), SegmentType::MolecularBarcode),
        ];
        let expected = vec![
            seg("GATA".as_bytes(), SegmentType::Template),
            seg("CAC".as_bytes(), SegmentType::Template),
        ];

        let read_set = read_set(segments);

        assert_eq!(expected, read_set.template_segments().cloned().collect::<Vec<_>>());
    }
    #[test]
    fn test_sample_barcode_segments() {
        let segments = vec![
            seg("AGCT".as_bytes(), SegmentType::Template),
            seg("GATA".as_bytes(), SegmentType::SampleBarcode),
            seg("CAC".as_bytes(), SegmentType::SampleBarcode),
            seg("GACCCC".as_bytes(), SegmentType::MolecularBarcode),
        ];
        let expected = vec![
            seg("GATA".as_bytes(), SegmentType::SampleBarcode),
            seg("CAC".as_bytes(), SegmentType::SampleBarcode),
        ];

        let read_set = read_set(segments);

        assert_eq!(expected, read_set.sample_barcode_segments().cloned().collect::<Vec<_>>());
    }
    #[test]
    fn test_molecular_barcode_segments() {
        let segments = vec![
            seg("AGCT".as_bytes(), SegmentType::Template),
            seg("GATA".as_bytes(), SegmentType::MolecularBarcode),
            seg("CAC".as_bytes(), SegmentType::MolecularBarcode),
            seg("GACCCC".as_bytes(), SegmentType::SampleBarcode),
        ];
        let expected = vec![
            seg("GATA".as_bytes(), SegmentType::MolecularBarcode),
            seg("CAC".as_bytes(), SegmentType::MolecularBarcode),
        ];

        let read_set = read_set(segments);

        assert_eq!(expected, read_set.molecular_barcode_segments().cloned().collect::<Vec<_>>());
    }

    #[test]
    fn test_combine_readsets() {
        let segments1 = vec![
            seg("A".as_bytes(), SegmentType::Template),
            seg("G".as_bytes(), SegmentType::Template),
            seg("C".as_bytes(), SegmentType::MolecularBarcode),
            seg("T".as_bytes(), SegmentType::SampleBarcode),
        ];
        let read_set1 = read_set(segments1.clone());
        let segments2 = vec![
            seg("AA".as_bytes(), SegmentType::Template),
            seg("AG".as_bytes(), SegmentType::Template),
            seg("AC".as_bytes(), SegmentType::MolecularBarcode),
            seg("AT".as_bytes(), SegmentType::SampleBarcode),
        ];
        let read_set2 = read_set(segments2.clone());
        let segments3 = vec![
            seg("AAA".as_bytes(), SegmentType::Template),
            seg("AAG".as_bytes(), SegmentType::Template),
            seg("AAC".as_bytes(), SegmentType::MolecularBarcode),
            seg("AAT".as_bytes(), SegmentType::SampleBarcode),
        ];
        let read_set3 = read_set(segments3.clone());

        let mut expected_segments = Vec::new();
        expected_segments.extend(segments1);
        expected_segments.extend(segments2);
        expected_segments.extend(segments3);
        let expected = read_set(expected_segments);

        let result = FastqSet::combine_readsets(vec![read_set1, read_set2, read_set3]);

        assert_eq!(result, expected);
    }

    #[test]
    #[should_panic(expected = "Cannot call combine readsets on an empty vec!")]
    fn test_combine_readsets_fails_on_empty_vector() {
        let _result = FastqSet::combine_readsets(Vec::new());
    }

    // =====================================================================
    // Cell barcode segment tests
    // =====================================================================

    #[test]
    fn test_cell_barcode_segments() {
        let segments = vec![
            seg("AGCT".as_bytes(), SegmentType::Template),
            seg("GATA".as_bytes(), SegmentType::CellularBarcode),
            seg("CAC".as_bytes(), SegmentType::CellularBarcode),
            seg("GACCCC".as_bytes(), SegmentType::MolecularBarcode),
        ];
        let expected = vec![
            seg("GATA".as_bytes(), SegmentType::CellularBarcode),
            seg("CAC".as_bytes(), SegmentType::CellularBarcode),
        ];

        let read_set = read_set(segments);

        assert_eq!(expected, read_set.cell_barcode_segments().cloned().collect::<Vec<_>>());
    }

    #[test]
    fn test_cell_barcode_segments_empty() {
        let segments = vec![
            seg("AGCT".as_bytes(), SegmentType::Template),
            seg("GATA".as_bytes(), SegmentType::MolecularBarcode),
        ];

        let read_set = read_set(segments);

        assert_eq!(0, read_set.cell_barcode_segments().count());
    }

    // =====================================================================
    // SkipReason Display and FromStr tests
    // =====================================================================

    #[test]
    fn test_skip_reason_display() {
        assert_eq!(SkipReason::TooFewBases.to_string(), "Too few bases");
    }

    #[test]
    fn test_skip_reason_from_str_valid() {
        // Test various valid forms
        assert_eq!(SkipReason::from_str("too few bases").unwrap(), SkipReason::TooFewBases);
        assert_eq!(SkipReason::from_str("too-few-bases").unwrap(), SkipReason::TooFewBases);
        assert_eq!(SkipReason::from_str("toofewbases").unwrap(), SkipReason::TooFewBases);
    }

    #[test]
    fn test_skip_reason_from_str_invalid() {
        let result = SkipReason::from_str("invalid-reason");
        assert!(result.is_err());
        let err = result.unwrap_err();
        assert!(err.to_string().contains("Invalid skip reason"));
        assert!(err.to_string().contains("invalid-reason"));
    }

    // =====================================================================
    // ReadSetIterator tests
    // =====================================================================

    #[test]
    fn test_read_set_iterator_basic() {
        use std::io::Cursor;

        // Create a simple FASTQ record
        let fastq_data = b"@read1\nACGTACGT\n+\nIIIIIIII\n";
        let cursor = Cursor::new(fastq_data.to_vec());
        let reader: FastqReader<Box<dyn BufRead + Send>> = FastqReader::new(Box::new(cursor));

        // Simple read structure: 4M (molecular barcode) + 4T (template)
        let read_structure = ReadStructure::from_str("4M4T").unwrap();
        let mut iterator = ReadSetIterator::new(read_structure, reader, vec![]);

        let result = iterator.next();
        assert!(result.is_some());

        let fastq_set = result.unwrap();
        assert_eq!(fastq_set.header, b"read1");
        assert_eq!(fastq_set.segments.len(), 2);
        assert!(fastq_set.skip_reason.is_none());

        // First segment should be molecular barcode
        assert_eq!(fastq_set.segments[0].seq, b"ACGT");
        assert_eq!(fastq_set.segments[0].segment_type, SegmentType::MolecularBarcode);

        // Second segment should be template
        assert_eq!(fastq_set.segments[1].seq, b"ACGT");
        assert_eq!(fastq_set.segments[1].segment_type, SegmentType::Template);

        // No more records
        assert!(iterator.next().is_none());
    }

    #[test]
    fn test_read_set_iterator_skip_too_few_bases() {
        use std::io::Cursor;

        // Create a FASTQ record that's too short for the read structure
        let fastq_data = b"@read1\nACGT\n+\nIIII\n";
        let cursor = Cursor::new(fastq_data.to_vec());
        let reader: FastqReader<Box<dyn BufRead + Send>> = FastqReader::new(Box::new(cursor));

        // Read structure requires 10 bases total
        let read_structure = ReadStructure::from_str("4M6T").unwrap();
        let mut iterator =
            ReadSetIterator::new(read_structure, reader, vec![SkipReason::TooFewBases]);

        let result = iterator.next();
        assert!(result.is_some());

        let fastq_set = result.unwrap();
        assert_eq!(fastq_set.header, b"read1");
        assert!(fastq_set.segments.is_empty());
        assert_eq!(fastq_set.skip_reason, Some(SkipReason::TooFewBases));
    }

    #[test]
    #[should_panic(expected = "too few bases")]
    fn test_read_set_iterator_panic_on_too_few_bases_without_skip() {
        use std::io::Cursor;

        // Create a FASTQ record that's too short
        let fastq_data = b"@read1\nACGT\n+\nIIII\n";
        let cursor = Cursor::new(fastq_data.to_vec());
        let reader: FastqReader<Box<dyn BufRead + Send>> = FastqReader::new(Box::new(cursor));

        // Read structure requires 10 bases total
        let read_structure = ReadStructure::from_str("4M6T").unwrap();
        let mut iterator = ReadSetIterator::new(read_structure, reader, vec![]);

        // This should panic because TooFewBases is not in skip_reasons
        let _ = iterator.next();
    }

    #[test]
    fn test_read_set_iterator_multiple_records() {
        use std::io::Cursor;

        let fastq_data = b"@read1\nACGTAAAA\n+\nIIIIIIII\n@read2\nTGCATTTT\n+\nIIIIIIII\n";
        let cursor = Cursor::new(fastq_data.to_vec());
        let reader: FastqReader<Box<dyn BufRead + Send>> = FastqReader::new(Box::new(cursor));

        let read_structure = ReadStructure::from_str("4M4T").unwrap();
        let mut iterator = ReadSetIterator::new(read_structure, reader, vec![]);

        // First record
        let first = iterator.next().unwrap();
        assert_eq!(first.header, b"read1");
        assert_eq!(first.segments[0].seq, b"ACGT");
        assert_eq!(first.segments[1].seq, b"AAAA");

        // Second record
        let second = iterator.next().unwrap();
        assert_eq!(second.header, b"read2");
        assert_eq!(second.segments[0].seq, b"TGCA");
        assert_eq!(second.segments[1].seq, b"TTTT");

        // No more records
        assert!(iterator.next().is_none());
    }

    #[test]
    fn test_read_set_iterator_variable_length_segment() {
        use std::io::Cursor;

        // Test with variable-length template segment (last segment gets remaining bases)
        let fastq_data = b"@read1\nACGTTTTTTTTT\n+\nIIIIIIIIIIII\n";
        let cursor = Cursor::new(fastq_data.to_vec());
        let reader: FastqReader<Box<dyn BufRead + Send>> = FastqReader::new(Box::new(cursor));

        // 4M + variable T (remaining bases go to template)
        let read_structure = ReadStructure::from_str("4M+T").unwrap();
        let mut iterator = ReadSetIterator::new(read_structure, reader, vec![]);

        let result = iterator.next().unwrap();
        assert_eq!(result.segments.len(), 2);

        // Fixed molecular barcode segment
        assert_eq!(result.segments[0].seq, b"ACGT");
        assert_eq!(result.segments[0].segment_type, SegmentType::MolecularBarcode);

        // Variable template gets remaining 8 bases
        assert_eq!(result.segments[1].seq, b"TTTTTTTT");
        assert_eq!(result.segments[1].segment_type, SegmentType::Template);
    }

    #[test]
    fn test_fastq_set_with_skip_reason_only() {
        // Test creating a FastqSet with skip_reason and empty segments
        let fastq_set = FastqSet {
            header: b"skipped_read".to_vec(),
            segments: vec![],
            skip_reason: Some(SkipReason::TooFewBases),
        };

        assert_eq!(fastq_set.skip_reason, Some(SkipReason::TooFewBases));
        assert!(fastq_set.segments.is_empty());
        assert_eq!(fastq_set.template_segments().count(), 0);
        assert_eq!(fastq_set.molecular_barcode_segments().count(), 0);
        assert_eq!(fastq_set.sample_barcode_segments().count(), 0);
        assert_eq!(fastq_set.cell_barcode_segments().count(), 0);
    }
}
