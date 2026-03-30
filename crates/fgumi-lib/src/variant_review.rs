//! Support for reviewing consensus variants
//!
//! This module provides structures and utilities for reviewing variant calls
//! from consensus reads, including tracking base counts and pileup information.

use crate::phred::NO_CALL_BASE;
use serde::Serialize;

/// Detailed information about a consensus read carrying a variant.
///
/// Each row contains information about the variant site (repeated for each
/// consensus read) and specific information about that consensus read.
#[derive(Debug, Clone, Serialize)]
pub struct ConsensusVariantReviewInfo {
    // Variant site information (same for all reads at this position)
    /// Chromosome/contig name
    pub chrom: String,
    /// Position of the variant (1-based)
    pub pos: i32,
    /// Reference allele
    pub ref_allele: String,
    /// Genotype string (if available from VCF)
    pub genotype: String,
    /// Comma-separated list of filters (or "PASS")
    pub filters: String,

    // Consensus-level counts (same for all reads at this position)
    /// Count of A observations at this position across all consensus reads
    #[serde(rename = "A")]
    pub consensus_a: usize,
    /// Count of C observations
    #[serde(rename = "C")]
    pub consensus_c: usize,
    /// Count of G observations
    #[serde(rename = "G")]
    pub consensus_g: usize,
    /// Count of T observations
    #[serde(rename = "T")]
    pub consensus_t: usize,
    /// Count of N observations
    #[serde(rename = "N")]
    pub consensus_n: usize,

    // Per-consensus-read information
    /// Consensus read name
    pub consensus_read: String,
    /// Insert description (chr:start-end | orientation)
    pub consensus_insert: String,
    /// Base call from the consensus read
    pub consensus_call: char,
    /// Quality score from the consensus read
    pub consensus_qual: u8,

    // Raw read counts supporting this consensus base
    /// Number of As in raw reads supporting this consensus base
    pub a: usize,
    /// Number of Cs in raw reads
    pub c: usize,
    /// Number of Gs in raw reads
    pub g: usize,
    /// Number of Ts in raw reads
    pub t: usize,
    /// Number of Ns in raw reads
    pub n: usize,
}

/// Represents a variant position
#[derive(Debug, Clone)]
pub struct Variant {
    pub chrom: String,
    pub pos: i32,
    pub ref_base: char,
    pub genotype: Option<String>,
    pub filters: Option<String>,
}

impl Variant {
    #[must_use]
    pub fn new(chrom: String, pos: i32, ref_base: char) -> Self {
        Self { chrom, pos, ref_base, genotype: None, filters: None }
    }

    #[must_use]
    pub fn with_genotype(mut self, genotype: String) -> Self {
        self.genotype = Some(genotype);
        self
    }

    #[must_use]
    pub fn with_filters(mut self, filters: String) -> Self {
        self.filters = Some(filters);
        self
    }
}

/// Base counts at a position
#[derive(Debug, Default, Clone)]
pub struct BaseCounts {
    pub a: usize,
    pub c: usize,
    pub g: usize,
    pub t: usize,
    pub n: usize,
}

impl BaseCounts {
    pub fn add_base(&mut self, base: u8) {
        match base.to_ascii_uppercase() {
            b'A' => self.a += 1,
            b'C' => self.c += 1,
            b'G' => self.g += 1,
            b'T' => self.t += 1,
            NO_CALL_BASE => self.n += 1,
            _ => {}
        }
    }

    #[must_use]
    pub fn total(&self) -> usize {
        self.a + self.c + self.g + self.t + self.n
    }
}

/// Generates a /1 or /2 suffix based on read pair information
#[must_use]
pub fn read_number_suffix(is_first_of_pair: bool) -> &'static str {
    if is_first_of_pair { "/1" } else { "/2" }
}

/// Generates an insert string for a consensus read
///
/// Returns a string like "chr1:100-200 | F1R2" for FR pairs, or "NA" for non-FR pairs
#[must_use]
#[allow(clippy::cast_possible_truncation, clippy::cast_possible_wrap)]
pub fn format_insert_string(
    record: &noodles::sam::alignment::RecordBuf,
    header: &noodles::sam::Header,
) -> String {
    let flags = record.flags();

    // Only generate strings for FR pairs, all other reads get "NA"
    // Must be paired, both reads mapped, and on same reference with FR orientation
    if !flags.is_segmented() {
        return "NA".to_string();
    }

    if flags.is_unmapped() || flags.is_mate_unmapped() {
        return "NA".to_string();
    }

    // Check that mate is on same reference
    let Some(ref_idx) = record.reference_sequence_id() else {
        return "NA".to_string();
    };

    let Some(mate_ref_idx) = record.mate_reference_sequence_id() else {
        return "NA".to_string();
    };

    if ref_idx != mate_ref_idx {
        return "NA".to_string();
    }

    // Check for FR orientation: one forward, one reverse
    let is_reverse = flags.is_reverse_complemented();
    let mate_is_reverse = flags.is_mate_reverse_complemented();

    if is_reverse == mate_is_reverse {
        return "NA".to_string(); // Both same orientation (FF or RR)
    }

    // Get reference sequence name
    let ref_name = match record.reference_sequence(header) {
        Some(Ok((name, _))) => String::from_utf8_lossy(name).to_string(),
        _ => return "NA".to_string(),
    };

    // Calculate outer position based on Scala logic:
    // val outer = if (r.negativeStrand) r.end else r.start
    let outer = if is_reverse {
        // If negative strand, use alignment end
        match record.alignment_end() {
            Some(pos) => usize::from(pos) as i32,
            None => return "NA".to_string(),
        }
    } else {
        // If positive strand, use alignment start
        match record.alignment_start() {
            Some(pos) => usize::from(pos) as i32,
            None => return "NA".to_string(),
        }
    };

    // Get insert size (signed)
    let isize = record.template_length();

    // Calculate start and end coordinates using Scala logic:
    // val Seq(start, end) = Seq(outer, outer + isize + (if (isize < 0) 1 else -1)).sorted
    let other_coord = outer + isize + if isize < 0 { 1 } else { -1 };
    let (start, end) =
        if outer < other_coord { (outer, other_coord) } else { (other_coord, outer) };

    // Determine pairing using Scala logic:
    // val pairing = (r.firstOfPair, start == outer) match {
    //   case (true, true)   => "F1R2"
    //   case (true, false)  => "F2R1"
    //   case (false, true)  => "F2R1"
    //   case (false, false) => "F1R2"
    // }
    let is_first = flags.is_first_segment();
    let start_equals_outer = start == outer;

    let pairing = match (is_first, start_equals_outer) {
        (true, true) | (false, false) => "F1R2",
        (true, false) | (false, true) => "F2R1",
    };

    format!("{ref_name}:{start}-{end} | {pairing}")
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::sam::builder::RecordBuilder;
    use noodles::sam::Header;
    use noodles::sam::alignment::RecordBuf;

    #[test]
    fn test_read_number_suffix() {
        assert_eq!(read_number_suffix(true), "/1");
        assert_eq!(read_number_suffix(false), "/2");
    }

    #[test]
    fn test_format_insert_string_unpaired() {
        // Create a simple SAM header
        let header = Header::builder().build();

        // Create an unpaired read (no segmented flag)
        let record = RecordBuilder::new().sequence("ACGT").build();

        let result = format_insert_string(&record, &header);
        assert_eq!(result, "NA");
    }

    #[test]
    fn test_base_counts_new() {
        let counts = BaseCounts::default();
        assert_eq!(counts.a, 0);
        assert_eq!(counts.c, 0);
        assert_eq!(counts.g, 0);
        assert_eq!(counts.t, 0);
        assert_eq!(counts.n, 0);
        assert_eq!(counts.total(), 0);
    }

    #[test]
    fn test_base_counts_add_base() {
        let mut counts = BaseCounts::default();

        counts.add_base(b'A');
        assert_eq!(counts.a, 1);
        assert_eq!(counts.total(), 1);

        counts.add_base(b'a'); // Test lowercase
        assert_eq!(counts.a, 2);

        counts.add_base(b'C');
        counts.add_base(b'G');
        counts.add_base(b'T');
        counts.add_base(b'N');

        assert_eq!(counts.c, 1);
        assert_eq!(counts.g, 1);
        assert_eq!(counts.t, 1);
        assert_eq!(counts.n, 1);
        assert_eq!(counts.total(), 6);
    }

    #[test]
    fn test_base_counts_ignore_non_bases() {
        let mut counts = BaseCounts::default();

        counts.add_base(b'X');
        counts.add_base(b'-');
        counts.add_base(b'.');

        assert_eq!(counts.total(), 0);
    }

    #[test]
    fn test_variant_new() {
        let variant = Variant::new("chr1".to_string(), 100, 'A');
        assert_eq!(variant.chrom, "chr1");
        assert_eq!(variant.pos, 100);
        assert_eq!(variant.ref_base, 'A');
        assert!(variant.genotype.is_none());
        assert!(variant.filters.is_none());
    }

    #[test]
    fn test_variant_with_genotype_and_filters() {
        let variant = Variant::new("chr1".to_string(), 100, 'A')
            .with_genotype("0/1".to_string())
            .with_filters("PASS".to_string());

        assert_eq!(variant.genotype, Some("0/1".to_string()));
        assert_eq!(variant.filters, Some("PASS".to_string()));
    }

    // Helper to create a test header with chr1
    fn create_test_header() -> Header {
        use bstr::BString;
        use noodles::sam::header::record::value::{Map, map::ReferenceSequence};
        use std::num::NonZeroUsize;

        Header::builder()
            .add_reference_sequence(
                BString::from("chr1"),
                Map::<ReferenceSequence>::new(NonZeroUsize::try_from(1000).unwrap()),
            )
            .build()
    }

    // Helper to create a paired read with specific parameters
    fn create_paired_read(
        start: usize,
        mate_start: usize,
        is_reverse: bool,
        mate_is_reverse: bool,
        is_first: bool,
        template_len: i32,
    ) -> RecordBuf {
        RecordBuilder::new()
            .sequence(&"A".repeat(10))
            .cigar("10M")
            .reference_sequence_id(0)
            .alignment_start(start)
            .mate_reference_sequence_id(0)
            .mate_alignment_start(mate_start)
            .template_length(template_len)
            .first_segment(is_first)
            .reverse_complement(is_reverse)
            .mate_reverse_complement(mate_is_reverse)
            .build()
    }

    #[test]
    fn test_format_insert_string_unmapped_returns_na() {
        use noodles::sam::alignment::record::Flags;

        let header = create_test_header();
        let mut record =
            RecordBuilder::new().sequence("ACGT").first_segment(true).unmapped(true).build();

        // Ensure SEGMENTED | UNMAPPED flags are set
        *record.flags_mut() = Flags::SEGMENTED | Flags::UNMAPPED;

        let result = format_insert_string(&record, &header);
        assert_eq!(result, "NA");
    }

    #[test]
    fn test_format_insert_string_mate_unmapped_returns_na() {
        let header = create_test_header();
        let mut record = create_paired_read(100, 200, false, true, true, 101);

        // Set mate as unmapped
        *record.flags_mut() |= noodles::sam::alignment::record::Flags::MATE_UNMAPPED;

        let result = format_insert_string(&record, &header);
        assert_eq!(result, "NA");
    }

    #[test]
    fn test_format_insert_string_same_orientation_returns_na() {
        let header = create_test_header();

        // Both forward (FF)
        let record = create_paired_read(100, 200, false, false, true, 101);
        let result = format_insert_string(&record, &header);
        assert_eq!(result, "NA");

        // Both reverse (RR)
        let record = create_paired_read(100, 200, true, true, true, -101);
        let result = format_insert_string(&record, &header);
        assert_eq!(result, "NA");
    }

    #[test]
    fn test_format_insert_string_cross_contig_returns_na() {
        let header = create_test_header();
        let mut record = create_paired_read(100, 200, false, true, true, 101);

        // Set mate on different reference
        *record.mate_reference_sequence_id_mut() = Some(1);

        let result = format_insert_string(&record, &header);
        assert_eq!(result, "NA");
    }

    #[test]
    fn test_format_insert_string_fr_pair_f1r2() {
        let header = create_test_header();

        // First read forward (start=100), second read reverse (start=191)
        // Template length = 101 (positive for forward read)
        // Expected: chr1:100-200 | F1R2
        let record = create_paired_read(100, 191, false, true, true, 101);
        let result = format_insert_string(&record, &header);
        assert_eq!(result, "chr1:100-200 | F1R2");
    }

    #[test]
    fn test_format_insert_string_fr_pair_f2r1() {
        let header = create_test_header();

        // First read reverse (start=191), second read forward (start=100)
        // Template length = -101 (negative for reverse read)
        // Expected: chr1:100-200 | F2R1
        let record = create_paired_read(191, 100, true, false, true, -101);
        let result = format_insert_string(&record, &header);
        assert_eq!(result, "chr1:100-200 | F2R1");
    }

    #[test]
    fn test_format_insert_string_second_of_pair_forward() {
        let header = create_test_header();

        // Second read forward (start=100), first read reverse (start=191)
        // is_first=false, forward orientation
        // Expected: chr1:100-200 | F2R1
        let record = create_paired_read(100, 191, false, true, false, 101);
        let result = format_insert_string(&record, &header);
        assert_eq!(result, "chr1:100-200 | F2R1");
    }

    #[test]
    fn test_format_insert_string_second_of_pair_reverse() {
        let header = create_test_header();

        // Second read reverse (start=191), first read forward (start=100)
        // is_first=false, reverse orientation
        // Expected: chr1:100-200 | F1R2
        let record = create_paired_read(191, 100, true, false, false, -101);
        let result = format_insert_string(&record, &header);
        assert_eq!(result, "chr1:100-200 | F1R2");
    }
}
