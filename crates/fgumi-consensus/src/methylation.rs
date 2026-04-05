//! Methylation-aware consensus calling for EM-Seq and TAPs data.
//!
//! EM-Seq enzymatically converts unmethylated cytosines to thymine before PCR.
//! At a reference C position: methylated → reads show C; unmethylated → reads show T.
//!
//! TAPs/Illumina 5-base converts methylated cytosines to thymine.
//! At a reference C position: unmethylated → reads show C; methylated → reads show T.
//!
//! The base counting (unconverted C vs converted T) is identical in both chemistries.
//! Only the MM/ML probability interpretation differs: EM-Seq uses unconverted/total,
//! TAPs uses converted/total.
//!
//! This module provides:
//! - Reference position mapping from consensus coordinates
//! - Per-base methylation evidence tracking
//! - SAM-spec MM/ML tag generation for methylation calls
//! - Dense cu/ct count tags for unconverted/converted evidence

use crate::vanilla_caller::SourceRead;
use fgumi_raw_bam::flags;
use noodles::sam::alignment::record::cigar::op::Kind;

/// Per-base methylation evidence at a single consensus position.
#[derive(Debug, Clone, Default)]
pub struct MethylationEvidence {
    /// Whether this position is a reference cytosine (eligible for methylation call).
    /// For top-strand reads: ref=C. For bottom-strand reads (after RC): ref=G.
    pub is_ref_c: bool,
    /// Number of reads showing C (unconverted) at this ref-C position.
    /// For bottom strand (after RC): number showing G (complement of unconverted C).
    /// Stored as u32 to avoid overflow at high coverage; clamped to i16 on output.
    pub unconverted_count: u32,
    /// Number of reads showing T (converted) at this ref-C position.
    /// For bottom strand (after RC): number showing A (complement of converted T).
    /// Stored as u32 to avoid overflow at high coverage; clamped to i16 on output.
    pub converted_count: u32,
}

/// Methylation annotation for an entire consensus read.
#[derive(Debug, Clone)]
pub struct MethylationAnnotation {
    /// Per-base methylation evidence (same length as consensus read).
    pub evidence: Vec<MethylationEvidence>,
}

impl MethylationAnnotation {
    /// Returns the unconverted counts as an i16 array for the `cu` tag.
    /// Values are clamped to `i16::MAX` to fit the BAM tag format.
    #[must_use]
    pub fn unconverted_counts(&self) -> Vec<i16> {
        self.evidence
            .iter()
            .map(|e| i16::try_from(e.unconverted_count).unwrap_or(i16::MAX))
            .collect()
    }

    /// Returns the converted counts as an i16 array for the `ct` tag.
    /// Values are clamped to `i16::MAX` to fit the BAM tag format.
    #[must_use]
    pub fn converted_counts(&self) -> Vec<i16> {
        self.evidence.iter().map(|e| i16::try_from(e.converted_count).unwrap_or(i16::MAX)).collect()
    }

    /// Returns a truncated copy of this annotation with only the first `len` positions.
    #[must_use]
    pub fn truncate(&self, len: usize) -> Self {
        Self { evidence: self.evidence[..len.min(self.evidence.len())].to_vec() }
    }
}

/// Maps each query position to a reference position using a simplified CIGAR.
///
/// For forward-strand reads, walks from `alignment_start` forward.
/// For reverse-strand reads (where CIGAR and bases have been reversed in
/// `create_source_read`), we reconstruct the original reference span and
/// map positions from the rightmost position backward.
///
/// Returns `Vec<Option<i64>>` where `None` = insertion (no ref base).
#[must_use]
#[expect(
    clippy::cast_possible_wrap,
    reason = "CIGAR lengths are small enough that usize→i64 won't wrap"
)]
pub fn query_to_ref_positions(
    simplified_cigar: &[(Kind, usize)],
    alignment_start: i64,
    is_reverse: bool,
    original_cigar: &[(Kind, usize)],
) -> Vec<Option<i64>> {
    // Calculate total query length from the (possibly reversed) cigar
    let query_len: usize =
        simplified_cigar.iter().filter(|(k, _)| k.consumes_read()).map(|(_, len)| *len).sum();

    let mut positions = Vec::with_capacity(query_len);

    if is_reverse {
        // For reverse strand: the CIGAR has been reversed in create_source_read.
        // We need to compute the alignment end from the *original* cigar, then
        // walk the reversed cigar mapping positions from right to left in reference space.
        let ref_span: i64 = original_cigar
            .iter()
            .filter(|(k, _)| k.consumes_reference())
            .map(|(_, len)| *len as i64)
            .sum();
        let alignment_end = alignment_start + ref_span - 1; // 0-based inclusive end

        let mut ref_pos = alignment_end;
        for &(kind, len) in simplified_cigar {
            match kind {
                Kind::Match | Kind::SequenceMatch | Kind::SequenceMismatch => {
                    for _ in 0..len {
                        positions.push(Some(ref_pos));
                        ref_pos -= 1;
                    }
                }
                Kind::Insertion | Kind::SoftClip => {
                    for _ in 0..len {
                        positions.push(None);
                    }
                }
                Kind::Deletion | Kind::Skip => {
                    ref_pos -= len as i64;
                }
                _ => {}
            }
        }
    } else {
        // Forward strand: walk from alignment_start forward
        let mut ref_pos = alignment_start;
        for &(kind, len) in simplified_cigar {
            match kind {
                Kind::Match | Kind::SequenceMatch | Kind::SequenceMismatch => {
                    for _ in 0..len {
                        positions.push(Some(ref_pos));
                        ref_pos += 1;
                    }
                }
                Kind::Insertion | Kind::SoftClip => {
                    for _ in 0..len {
                        positions.push(None);
                    }
                }
                Kind::Deletion | Kind::Skip => {
                    ref_pos += len as i64;
                }
                _ => {}
            }
        }
    }

    positions
}

/// Annotates simplex consensus methylation from source reads and reference.
///
/// For each consensus position that aligns to a reference cytosine (top strand)
/// or reference guanine (bottom strand, after RC), counts source reads showing
/// unconverted vs converted bases. Does not modify the consensus bases.
///
/// # Arguments
/// * `consensus_bases` - Consensus bases (used for length/position mapping)
/// * `source_reads` - Source reads used to build this consensus
/// * `ref_bases_at_positions` - Reference bases at aligned positions (`None` for insertions)
/// * `is_top_strand` - Whether this consensus represents the top (forward) strand
pub(crate) fn annotate_simplex_methylation(
    consensus_bases: &[u8],
    source_reads: &[SourceRead],
    ref_bases_at_positions: &[Option<u8>],
    is_top_strand: bool,
) -> MethylationAnnotation {
    let len = consensus_bases.len();
    let mut evidence = vec![MethylationEvidence::default(); len];

    // Determine which reference base indicates a C position and what conversion looks like
    // Top strand: ref=C, unconverted=C, converted=T
    // Bottom strand (after RC): ref=G (complement of C), unconverted=G, converted=A
    let (ref_target, unconverted_base, converted_base) =
        if is_top_strand { (b'C', b'C', b'T') } else { (b'G', b'G', b'A') };

    for (i, ev) in evidence.iter_mut().enumerate() {
        // Check if this position aligns to a reference C/G
        let ref_base = ref_bases_at_positions.get(i).and_then(|b| *b);
        let Some(rb) = ref_base else { continue };
        let rb_upper = rb.to_ascii_uppercase();
        if rb_upper != ref_target {
            continue;
        }

        ev.is_ref_c = true;

        // Count unconverted vs converted in source reads
        for sr in source_reads {
            if i >= sr.bases.len() {
                continue;
            }
            let base = sr.bases[i].to_ascii_uppercase();
            if base == unconverted_base {
                ev.unconverted_count = ev.unconverted_count.saturating_add(1);
            } else if base == converted_base {
                ev.converted_count = ev.converted_count.saturating_add(1);
            }
        }

        // Base normalization (T→C / A→G at ref-C positions) is NOT done here — it is
        // handled by the caller (annotate_and_normalize) which normalizes source reads
        // after annotation so that conversion events don't inflate consensus error counts.
    }

    MethylationAnnotation { evidence }
}

/// Builds SAM-spec MM:Z and ML:B:C tags from methylation annotation.
///
/// MM format: `C+m,skip1,skip2,...;` listing skip counts between modified C bases.
/// ML companion array: one probability [0-255] per modification listed in MM.
///
/// For top-strand reads, we track `C+m` modifications (5mC on same strand as SEQ).
/// For bottom-strand reads (after RC), the consensus has G bases where the original
/// bottom-strand had C. Per the SAM spec, opposite-strand 5mC is encoded as `G-m`
/// (minus marker indicates the modification is on the opposite strand from SEQ).
///
/// The `methylation_mode` parameter controls the probability calculation:
/// - EM-Seq: prob = unconverted/total (C stayed as C because it was methylated)
/// - TAPs: prob = converted/total (C was converted to T because it was methylated)
///
/// Returns `(mm_string, ml_array)`. Returns `None` if no ref-C positions exist.
///
/// # Panics
///
/// Panics if `consensus_bases` and `annotation.evidence` have different lengths.
#[must_use]
pub fn build_mm_ml_tags(
    consensus_bases: &[u8],
    annotation: &MethylationAnnotation,
    is_top_strand: bool,
    methylation_mode: crate::MethylationMode,
) -> Option<(String, Vec<u8>)> {
    assert_eq!(
        consensus_bases.len(),
        annotation.evidence.len(),
        "consensus_bases and annotation.evidence must have the same length"
    );

    // The base we track in MM depends on strand
    let track_base = if is_top_strand { b'C' } else { b'G' };

    let mut skips = Vec::new();
    let mut probs = Vec::new();
    let mut skip_count: usize = 0;

    for (i, ev) in annotation.evidence.iter().enumerate() {
        let base_upper = consensus_bases[i].to_ascii_uppercase();
        if base_upper != track_base {
            continue;
        }

        if ev.is_ref_c {
            // This is a ref-C position with a C/G in consensus
            let total = u64::from(ev.unconverted_count) + u64::from(ev.converted_count);
            if total > 0 {
                // EM-Seq: methylation prob = unconverted/total (C = methylated, stayed as C)
                // TAPs:   methylation prob = converted/total  (T = methylated, converted from C)
                let numerator = match methylation_mode {
                    crate::MethylationMode::EmSeq => u64::from(ev.unconverted_count),
                    crate::MethylationMode::Taps => u64::from(ev.converted_count),
                    crate::MethylationMode::Disabled => return None,
                };
                let prob = (numerator * 255 / total).min(255) as u8;
                skips.push(skip_count);
                probs.push(prob);
                skip_count = 0;
            } else {
                skip_count += 1;
            }
        } else {
            // C/G in consensus but not at a ref-C position — just skip it
            skip_count += 1;
        }
    }

    if skips.is_empty() {
        return None;
    }

    // Build MM string: "C+m,skip1,skip2,...;" (top) or "G-m,skip1,...;" (bottom)
    let (base_char, strand_marker) = if is_top_strand { ('C', '+') } else { ('G', '-') };
    let mut mm = format!("{base_char}{strand_marker}m");
    for s in &skips {
        use std::fmt::Write;
        write!(mm, ",{s}").unwrap();
    }
    mm.push(';');

    Some((mm, probs))
}

/// Builds an MM-format tag string without ML companion (for per-strand am/bm tags).
///
/// Same format as `build_mm_ml_tags` but returns only the MM:Z string.
#[must_use]
pub fn build_mm_tag_no_ml(
    consensus_bases: &[u8],
    annotation: &MethylationAnnotation,
    is_top_strand: bool,
    methylation_mode: crate::MethylationMode,
) -> Option<String> {
    build_mm_ml_tags(consensus_bases, annotation, is_top_strand, methylation_mode).map(|(mm, _)| mm)
}

/// Fetches reference bases for aligned positions.
///
/// Given a set of query-to-ref position mappings and a reference fetch function,
/// returns the reference base at each position (or `None` for insertions/unmapped).
pub fn fetch_ref_bases_at_positions(
    ref_positions: &[Option<i64>],
    ref_name: &str,
    reference: &dyn RefBaseProvider,
) -> Vec<Option<u8>> {
    // Use sequence_for for O(1) per-base access when available, falling back to
    // per-base HashMap lookups only when the implementor doesn't provide it.
    if let Some(seq) = reference.sequence_for(ref_name) {
        ref_positions
            .iter()
            .map(|pos| pos.and_then(|p| usize::try_from(p).ok().and_then(|i| seq.get(i).copied())))
            .collect()
    } else {
        ref_positions
            .iter()
            .map(|pos| {
                pos.and_then(|p| {
                    u64::try_from(p).ok().and_then(|pos| reference.base_at_0based(ref_name, pos))
                })
            })
            .collect()
    }
}

/// Trait for providing reference bases (allows testing without full `ReferenceReader`).
pub trait RefBaseProvider {
    /// Returns the base at a 0-based position, or None if out of bounds.
    fn base_at_0based(&self, chrom: &str, pos: u64) -> Option<u8>;

    /// Returns the full sequence for a chromosome, or None if not found.
    ///
    /// Default implementation returns None (forcing per-base lookup fallback).
    /// Implementors with in-memory sequences should override for O(1) access.
    fn sequence_for(&self, _chrom: &str) -> Option<&[u8]> {
        None
    }
}

/// Determines whether a `SourceRead` was originally on the top (forward) strand.
///
/// In EM-Seq, the "top strand" is the forward strand of the original molecule.
/// For paired-end reads:
/// - R1 forward (not reverse) = top strand
/// - R1 reverse = bottom strand
/// - R2 follows mate orientation (opposite of R1)
#[must_use]
pub fn is_top_strand(source_read_flags: u16) -> bool {
    let is_reverse = source_read_flags & flags::REVERSE != 0;
    let is_r2 = source_read_flags & flags::LAST_SEGMENT != 0;
    // Top strand: R1 forward or R2 reverse
    // Bottom strand: R1 reverse or R2 forward
    is_reverse == is_r2
}

/// Combines two strand methylation annotations into a duplex annotation.
///
/// Sums unconverted and converted counts from both strands at each position.
#[must_use]
pub fn combine_methylation_annotations(
    ab: &MethylationAnnotation,
    ba: &MethylationAnnotation,
    len: usize,
) -> MethylationAnnotation {
    let mut evidence = Vec::with_capacity(len);
    for i in 0..len {
        let ab_ev = ab.evidence.get(i);
        let ba_ev = ba.evidence.get(i);
        let is_ref_c = ab_ev.is_some_and(|e| e.is_ref_c) || ba_ev.is_some_and(|e| e.is_ref_c);
        let unconverted = ab_ev
            .map_or(0, |e| e.unconverted_count)
            .saturating_add(ba_ev.map_or(0, |e| e.unconverted_count));
        let converted = ab_ev
            .map_or(0, |e| e.converted_count)
            .saturating_add(ba_ev.map_or(0, |e| e.converted_count));
        evidence.push(MethylationEvidence {
            is_ref_c,
            unconverted_count: unconverted,
            converted_count: converted,
        });
    }
    MethylationAnnotation { evidence }
}

#[cfg(test)]
#[allow(clippy::cast_possible_truncation, clippy::cast_possible_wrap, clippy::needless_range_loop)]
pub(crate) mod tests {
    use super::*;
    use std::collections::HashMap;

    /// Simple in-memory reference for testing. Reused across crate tests.
    pub(crate) struct TestRef {
        sequences: HashMap<String, Vec<u8>>,
    }

    impl TestRef {
        pub(crate) fn new(seqs: &[(&str, &[u8])]) -> Self {
            let mut sequences = HashMap::new();
            for (name, seq) in seqs {
                sequences.insert((*name).to_string(), seq.to_vec());
            }
            Self { sequences }
        }
    }

    impl RefBaseProvider for TestRef {
        fn base_at_0based(&self, chrom: &str, pos: u64) -> Option<u8> {
            self.sequences.get(chrom).and_then(|seq| seq.get(pos as usize).copied())
        }

        fn sequence_for(&self, chrom: &str) -> Option<&[u8]> {
            self.sequences.get(chrom).map(Vec::as_slice)
        }
    }

    #[test]
    fn test_query_to_ref_positions_all_matches() {
        // 10M cigar, forward strand
        let cigar = vec![(Kind::Match, 10)];
        let positions = query_to_ref_positions(&cigar, 100, false, &cigar);
        assert_eq!(positions.len(), 10);
        for (i, pos) in positions.iter().enumerate() {
            assert_eq!(*pos, Some(100 + i as i64));
        }
    }

    #[test]
    fn test_query_to_ref_positions_with_insertion() {
        // 5M2I3M
        let cigar = vec![(Kind::Match, 5), (Kind::Insertion, 2), (Kind::Match, 3)];
        let positions = query_to_ref_positions(&cigar, 100, false, &cigar);
        assert_eq!(positions.len(), 10);
        // First 5: ref 100-104
        for i in 0..5 {
            assert_eq!(positions[i], Some(100 + i as i64));
        }
        // Insertions: None
        assert_eq!(positions[5], None);
        assert_eq!(positions[6], None);
        // Next 3: ref 105-107
        for i in 0..3 {
            assert_eq!(positions[7 + i], Some(105 + i as i64));
        }
    }

    #[test]
    fn test_query_to_ref_positions_with_deletion() {
        // 5M2D5M
        let cigar = vec![(Kind::Match, 5), (Kind::Deletion, 2), (Kind::Match, 5)];
        let positions = query_to_ref_positions(&cigar, 100, false, &cigar);
        assert_eq!(positions.len(), 10);
        for i in 0..5 {
            assert_eq!(positions[i], Some(100 + i as i64));
        }
        // After 2D, ref skips to 107
        for i in 0..5 {
            assert_eq!(positions[5 + i], Some(107 + i as i64));
        }
    }

    #[test]
    fn test_query_to_ref_positions_reverse_strand() {
        // Original cigar: 10M at position 100
        // After reversal in create_source_read: cigar is still 10M (symmetric)
        // Reverse strand should map positions from alignment_end backward
        let original_cigar = vec![(Kind::Match, 10)];
        let reversed_cigar = vec![(Kind::Match, 10)]; // Same since 10M reversed is 10M
        let positions = query_to_ref_positions(&reversed_cigar, 100, true, &original_cigar);
        assert_eq!(positions.len(), 10);
        // Should map from 109 down to 100
        for i in 0..10 {
            assert_eq!(positions[i], Some(109 - i as i64));
        }
    }

    #[test]
    fn test_annotate_simplex_all_methylated() {
        // All source reads show C at ref-C positions → methylated
        let consensus = b"ACGT".to_vec();
        let sr1 = make_test_source_read(b"ACGT", 0);
        let sr2 = make_test_source_read(b"ACGT", 0);
        let ref_bases = vec![Some(b'A'), Some(b'C'), Some(b'G'), Some(b'T')];

        let annot = annotate_simplex_methylation(
            &consensus,
            &[sr1, sr2],
            &ref_bases,
            true, // top strand
        );

        // Position 1 (ref=C): both reads show C → methylated, count=2
        assert!(annot.evidence[1].is_ref_c);
        assert_eq!(annot.evidence[1].unconverted_count, 2);
        assert_eq!(annot.evidence[1].converted_count, 0);
        // Consensus base should remain C
        assert_eq!(consensus[1], b'C');
    }

    #[test]
    fn test_annotate_simplex_all_unmethylated() {
        // All source reads show T at ref-C position → unmethylated
        let consensus = b"ATGT".to_vec(); // consensus called T at pos 1
        let sr1 = make_test_source_read(b"ATGT", 0);
        let sr2 = make_test_source_read(b"ATGT", 0);
        let ref_bases = vec![Some(b'A'), Some(b'C'), Some(b'G'), Some(b'T')];

        let annot = annotate_simplex_methylation(&consensus, &[sr1, sr2], &ref_bases, true);

        assert!(annot.evidence[1].is_ref_c);
        assert_eq!(annot.evidence[1].unconverted_count, 0);
        assert_eq!(annot.evidence[1].converted_count, 2);
        // Consensus base is NOT replaced — methylation state is tracked in cu/ct tags and MM/ML,
        // and base replacement would interfere with bwameth re-alignment
        assert_eq!(consensus[1], b'T');
    }

    #[test]
    fn test_annotate_simplex_mixed() {
        // Some C, some T at ref-C position
        let consensus = b"ACGT".to_vec();
        let sr1 = make_test_source_read(b"ACGT", 0); // C at pos 1
        let sr2 = make_test_source_read(b"ATGT", 0); // T at pos 1
        let ref_bases = vec![Some(b'A'), Some(b'C'), Some(b'G'), Some(b'T')];

        let annot = annotate_simplex_methylation(&consensus, &[sr1, sr2], &ref_bases, true);

        assert!(annot.evidence[1].is_ref_c);
        assert_eq!(annot.evidence[1].unconverted_count, 1);
        assert_eq!(annot.evidence[1].converted_count, 1);
    }

    #[test]
    fn test_annotate_simplex_non_c_positions() {
        // Non-ref-C positions should not be annotated
        let consensus = b"AGGT".to_vec();
        let sr1 = make_test_source_read(b"AGGT", 0);
        let ref_bases = vec![Some(b'A'), Some(b'G'), Some(b'G'), Some(b'T')];

        let annot = annotate_simplex_methylation(&consensus, &[sr1], &ref_bases, true);

        for ev in &annot.evidence {
            assert!(!ev.is_ref_c);
        }
    }

    #[test]
    fn test_annotate_simplex_reverse_strand() {
        // Bottom strand: ref=G (complement of C), unconverted=G, converted=A
        let consensus = b"CAGT".to_vec(); // A at pos 1 = converted on bottom strand
        let sr1 = make_test_source_read(b"CAGT", flags::REVERSE);
        let ref_bases = vec![Some(b'T'), Some(b'G'), Some(b'C'), Some(b'A')]; // at reversed positions

        let annot = annotate_simplex_methylation(
            &consensus,
            &[sr1],
            &ref_bases,
            false, // bottom strand
        );

        // Position 1: ref=G → eligible for bottom-strand methylation
        assert!(annot.evidence[1].is_ref_c);
        assert_eq!(annot.evidence[1].unconverted_count, 0); // A, not G
        assert_eq!(annot.evidence[1].converted_count, 1); // A = converted on bottom strand
        // Consensus base is NOT replaced — methylation state is tracked in cu/ct tags and MM/ML
        assert_eq!(consensus[1], b'A');
    }

    #[test]
    fn test_build_mm_ml_tags_basic() {
        let consensus = b"ACGCAC".to_vec();
        let annotation = MethylationAnnotation {
            evidence: vec![
                MethylationEvidence { is_ref_c: false, unconverted_count: 0, converted_count: 0 }, // A
                MethylationEvidence { is_ref_c: true, unconverted_count: 3, converted_count: 0 }, // C - methylated
                MethylationEvidence { is_ref_c: false, unconverted_count: 0, converted_count: 0 }, // G
                MethylationEvidence { is_ref_c: true, unconverted_count: 0, converted_count: 3 }, // C - unmethylated
                MethylationEvidence { is_ref_c: false, unconverted_count: 0, converted_count: 0 }, // A
                MethylationEvidence { is_ref_c: false, unconverted_count: 0, converted_count: 0 }, // C - not ref-C
            ],
        };

        let result = build_mm_ml_tags(&consensus, &annotation, true, crate::MethylationMode::EmSeq);
        assert!(result.is_some());
        let (mm, ml) = result.unwrap();
        // Two C bases that are ref-C: skip 0 to first, skip 0 to second
        // Third C is not ref-C
        assert_eq!(mm, "C+m,0,0;");
        assert_eq!(ml.len(), 2);
        assert_eq!(ml[0], 255); // fully methylated
        assert_eq!(ml[1], 0); // fully unmethylated
    }

    #[test]
    fn test_build_mm_ml_tags_no_modifications() {
        // No ref-C positions at all
        let consensus = b"AGGT".to_vec();
        let annotation = MethylationAnnotation {
            evidence: vec![
                MethylationEvidence::default(),
                MethylationEvidence::default(),
                MethylationEvidence::default(),
                MethylationEvidence::default(),
            ],
        };

        let result = build_mm_ml_tags(&consensus, &annotation, true, crate::MethylationMode::EmSeq);
        assert!(result.is_none());
    }

    #[test]
    fn test_build_mm_tag_no_ml() {
        let consensus = b"ACGT".to_vec();
        let annotation = MethylationAnnotation {
            evidence: vec![
                MethylationEvidence { is_ref_c: false, unconverted_count: 0, converted_count: 0 },
                MethylationEvidence { is_ref_c: true, unconverted_count: 2, converted_count: 1 },
                MethylationEvidence { is_ref_c: false, unconverted_count: 0, converted_count: 0 },
                MethylationEvidence { is_ref_c: false, unconverted_count: 0, converted_count: 0 },
            ],
        };

        let result =
            build_mm_tag_no_ml(&consensus, &annotation, true, crate::MethylationMode::EmSeq);
        assert!(result.is_some());
        assert_eq!(result.unwrap(), "C+m,0;");
    }

    #[test]
    fn test_is_top_strand() {
        // R1 forward = top strand
        assert!(is_top_strand(flags::PAIRED | flags::FIRST_SEGMENT));
        // R1 reverse = bottom strand
        assert!(!is_top_strand(flags::PAIRED | flags::FIRST_SEGMENT | flags::REVERSE));
        // R2 forward = bottom strand
        assert!(!is_top_strand(flags::PAIRED | flags::LAST_SEGMENT));
        // R2 reverse = top strand
        assert!(is_top_strand(flags::PAIRED | flags::LAST_SEGMENT | flags::REVERSE));
    }

    #[test]
    fn test_combine_methylation_annotations() {
        let ab = MethylationAnnotation {
            evidence: vec![
                MethylationEvidence { is_ref_c: true, unconverted_count: 2, converted_count: 1 },
                MethylationEvidence { is_ref_c: false, unconverted_count: 0, converted_count: 0 },
            ],
        };
        let ba = MethylationAnnotation {
            evidence: vec![
                MethylationEvidence { is_ref_c: true, unconverted_count: 1, converted_count: 2 },
                MethylationEvidence { is_ref_c: false, unconverted_count: 0, converted_count: 0 },
            ],
        };
        let combined = combine_methylation_annotations(&ab, &ba, 2);
        assert!(combined.evidence[0].is_ref_c);
        assert_eq!(combined.evidence[0].unconverted_count, 3);
        assert_eq!(combined.evidence[0].converted_count, 3);
        assert!(!combined.evidence[1].is_ref_c);
    }

    #[test]
    fn test_fetch_ref_bases_at_positions() {
        // Uses the sequence_for fast path (TestRef implements it)
        let reference = TestRef::new(&[("chr1", b"ACGTACGT")]);
        let positions = vec![Some(0), Some(1), None, Some(3), Some(7)];
        let bases = fetch_ref_bases_at_positions(&positions, "chr1", &reference);
        assert_eq!(bases, vec![Some(b'A'), Some(b'C'), None, Some(b'T'), Some(b'T')]);
    }

    /// A `RefBaseProvider` that only supports per-base lookups (no `sequence_for`).
    struct PerBaseLookupRef {
        sequences: HashMap<String, Vec<u8>>,
    }

    impl RefBaseProvider for PerBaseLookupRef {
        fn base_at_0based(&self, chrom: &str, pos: u64) -> Option<u8> {
            self.sequences.get(chrom).and_then(|seq| seq.get(pos as usize).copied())
        }
        // sequence_for not overridden → returns None → fallback path
    }

    #[test]
    fn test_fetch_ref_bases_fallback_path() {
        let mut sequences = HashMap::new();
        sequences.insert("chr1".to_string(), b"ACGTACGT".to_vec());
        let reference = PerBaseLookupRef { sequences };
        let positions = vec![Some(0), Some(1), None, Some(3), Some(7)];
        let bases = fetch_ref_bases_at_positions(&positions, "chr1", &reference);
        assert_eq!(bases, vec![Some(b'A'), Some(b'C'), None, Some(b'T'), Some(b'T')]);
    }

    #[test]
    fn test_build_mm_ml_with_skips() {
        // Consensus: CCACC — three C bases, but only positions 0 and 3 are ref-C
        let consensus = b"CCACC".to_vec();
        let annotation = MethylationAnnotation {
            evidence: vec![
                MethylationEvidence { is_ref_c: true, unconverted_count: 5, converted_count: 0 }, // C at 0: ref-C, methylated
                MethylationEvidence { is_ref_c: false, unconverted_count: 0, converted_count: 0 }, // C at 1: NOT ref-C
                MethylationEvidence { is_ref_c: false, unconverted_count: 0, converted_count: 0 }, // A at 2
                MethylationEvidence { is_ref_c: true, unconverted_count: 0, converted_count: 5 }, // C at 3: ref-C, unmethylated
                MethylationEvidence { is_ref_c: false, unconverted_count: 0, converted_count: 0 }, // C at 4: NOT ref-C
            ],
        };

        let (mm, ml) =
            build_mm_ml_tags(&consensus, &annotation, true, crate::MethylationMode::EmSeq).unwrap();
        // First ref-C is the 1st C (skip 0), second ref-C is the 4th C (skip 1 non-ref C + skip 1 more)
        // Walking: C at 0 (ref-C, skip=0), C at 1 (not ref-C, skip++), C at 3 (ref-C, skip=1), C at 4 (not ref-C)
        assert_eq!(mm, "C+m,0,1;");
        assert_eq!(ml, vec![255, 0]);
    }

    #[test]
    fn test_build_mm_ml_tags_bottom_strand() {
        // Bottom-strand consensus has G bases at methylation sites
        let consensus = b"AGCGAG".to_vec();
        let annotation = MethylationAnnotation {
            evidence: vec![
                MethylationEvidence { is_ref_c: false, unconverted_count: 0, converted_count: 0 }, // A
                MethylationEvidence { is_ref_c: true, unconverted_count: 3, converted_count: 0 }, // G - methylated
                MethylationEvidence { is_ref_c: false, unconverted_count: 0, converted_count: 0 }, // C
                MethylationEvidence { is_ref_c: true, unconverted_count: 0, converted_count: 3 }, // G - unmethylated
                MethylationEvidence { is_ref_c: false, unconverted_count: 0, converted_count: 0 }, // A
                MethylationEvidence { is_ref_c: false, unconverted_count: 0, converted_count: 0 }, // G - not ref-C
            ],
        };

        let (mm, ml) =
            build_mm_ml_tags(&consensus, &annotation, false, crate::MethylationMode::EmSeq)
                .expect("should have tags");
        // Per SAM spec: opposite-strand 5mC uses G-m (minus = opposite strand of SEQ)
        assert_eq!(mm, "G-m,0,0;");
        assert_eq!(ml.len(), 2);
        assert_eq!(ml[0], 255); // fully methylated
        assert_eq!(ml[1], 0); // fully unmethylated
    }

    #[test]
    fn test_methylation_counters_saturate() {
        let annotation = MethylationAnnotation {
            evidence: vec![MethylationEvidence {
                is_ref_c: true,
                unconverted_count: u32::MAX,
                converted_count: u32::MAX,
            }],
        };
        let ab = &annotation;
        let ba = &annotation;
        let combined = combine_methylation_annotations(ab, ba, 1);
        // Should saturate at u32::MAX, not wrap or panic
        assert_eq!(combined.evidence[0].unconverted_count, u32::MAX);
        assert_eq!(combined.evidence[0].converted_count, u32::MAX);
    }

    #[test]
    fn test_build_mm_ml_tags_taps_all_methylated() {
        // All converted (T) at ref-C → TAPs prob = converted/total = 255
        let bases = vec![b'C'; 5]; // consensus restored to C
        let annotation = MethylationAnnotation {
            evidence: vec![
                MethylationEvidence {
                    is_ref_c: true,
                    unconverted_count: 0,
                    converted_count: 3
                };
                5
            ],
        };
        let result = build_mm_ml_tags(&bases, &annotation, true, crate::MethylationMode::Taps);
        let (mm, ml) = result.unwrap();
        assert!(mm.starts_with("C+m"));
        assert_eq!(ml, vec![255u8; 5]);
    }

    #[test]
    fn test_build_mm_ml_tags_taps_all_unmethylated() {
        // All unconverted (C) at ref-C → TAPs prob = converted/total = 0/3 = 0
        let bases = vec![b'C'; 5];
        let annotation = MethylationAnnotation {
            evidence: vec![
                MethylationEvidence {
                    is_ref_c: true,
                    unconverted_count: 3,
                    converted_count: 0
                };
                5
            ],
        };
        let result = build_mm_ml_tags(&bases, &annotation, true, crate::MethylationMode::Taps);
        let (_, ml) = result.unwrap();
        assert_eq!(ml, vec![0u8; 5]);
    }

    #[test]
    fn test_build_mm_ml_tags_emseq_unchanged() {
        // Verify EM-seq behavior is unchanged: unconverted/total
        let bases = vec![b'C'; 5];
        let annotation = MethylationAnnotation {
            evidence: vec![
                MethylationEvidence {
                    is_ref_c: true,
                    unconverted_count: 3,
                    converted_count: 0
                };
                5
            ],
        };
        let result = build_mm_ml_tags(&bases, &annotation, true, crate::MethylationMode::EmSeq);
        let (_, ml) = result.unwrap();
        assert_eq!(ml, vec![255u8; 5]); // EM-seq: 3/3 unconverted = 255
    }

    /// Helper to create a `SourceRead` for testing.
    fn make_test_source_read(bases: &[u8], flg: u16) -> SourceRead {
        SourceRead {
            original_idx: 0,
            bases: bases.to_vec(),
            quals: vec![30; bases.len()],
            simplified_cigar: vec![(Kind::Match, bases.len())],
            flags: flg,
            ref_id: 0,
            alignment_start: 0,
            original_cigar: vec![(Kind::Match, bases.len())],
        }
    }
}
