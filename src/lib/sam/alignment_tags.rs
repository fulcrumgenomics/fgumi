//! Alignment tag regeneration (NM, UQ, MD) after base masking.
//!
//! When consensus reads have bases masked to 'N', alignment tags need to be recalculated
//! to accurately reflect the new sequence. This module provides functions to regenerate:
//!
//! - **NM**: Edit distance to the reference (mismatches + indels)
//! - **UQ**: Phred likelihood of the segment (sum of mismatch qualities)
//! - **MD**: Mismatched and deleted reference bases

use anyhow::{Context, Result};
use noodles::core::Position;
use noodles::sam::Header;
use noodles::sam::alignment::record::cigar::Cigar as CigarTrait;
use noodles::sam::alignment::record::cigar::op::Kind;
use noodles::sam::alignment::record::data::field::Tag;
use noodles::sam::alignment::record_buf::RecordBuf;
use noodles::sam::alignment::record_buf::data::field::Value;

use crate::reference::ReferenceReader;

/// Tags used for alignment information
#[must_use]
pub fn nm_tag() -> Tag {
    Tag::from([b'N', b'M'])
}

#[must_use]
pub fn md_tag() -> Tag {
    Tag::from([b'M', b'D'])
}

#[must_use]
pub fn uq_tag() -> Tag {
    Tag::from([b'U', b'Q'])
}

/// Regenerates NM, UQ, and MD tags for a record after base masking
///
/// For unmapped reads, the tags are removed (set to null) to match fgbio behavior.
/// For mapped reads, the tags are recalculated based on the alignment and reference.
///
/// # Arguments
/// * `record` - The record to regenerate tags for (modified in place)
/// * `header` - SAM header (needed to resolve reference sequence names)
/// * `reference` - Reference genome reader
///
/// # Returns
/// True if tags were regenerated, false if read is unmapped (tags are nulled)
pub fn regenerate_alignment_tags(
    record: &mut RecordBuf,
    header: &Header,
    reference: &ReferenceReader,
) -> Result<bool> {
    // For unmapped reads, null out the tags (matching fgbio behavior)
    if record.flags().is_unmapped() {
        record.data_mut().remove(&nm_tag());
        record.data_mut().remove(&uq_tag());
        record.data_mut().remove(&md_tag());
        return Ok(false);
    }

    // Get reference sequence ID and look up name in header
    let ref_seq_id = record.reference_sequence_id().context("Missing reference sequence ID")?;
    let ref_seqs = header.reference_sequences();
    let (ref_name_bytes, _) =
        ref_seqs.get_index(ref_seq_id).context("Reference sequence ID not found in header")?;
    let ref_name = std::str::from_utf8(ref_name_bytes.as_ref())?;

    let ref_start = record.alignment_start().context("Missing alignment start")?;
    let cigar = record.cigar();
    let seq = record.sequence();
    let qual = record.quality_scores();

    // Calculate total reference span from CIGAR and fetch entire alignment span once
    // This is a key optimization - instead of fetching per CIGAR operation, we fetch once
    // Note: Skip (N in CIGAR) is NOT included - it doesn't affect NM/MD/UQ calculation
    let ref_span: usize = cigar
        .iter()
        .filter_map(|op| op.ok())
        .map(|op| match op.kind() {
            Kind::Match | Kind::SequenceMatch | Kind::SequenceMismatch | Kind::Deletion => op.len(),
            _ => 0,
        })
        .sum();

    // Handle edge case: CIGAR with no reference-consuming operations (e.g., pure insertion "4I")
    // Set tags to sensible defaults and return early
    if ref_span == 0 {
        record.data_mut().insert(nm_tag(), Value::from(0u32));
        record.data_mut().insert(md_tag(), Value::String("0".to_owned().into()));
        record.data_mut().insert(uq_tag(), Value::from(0u32));
        return Ok(false);
    }

    // Fetch entire reference span in one call
    let ref_end = Position::new(usize::from(ref_start) + ref_span - 1)
        .context("Invalid reference end position")?;
    let all_ref_bases = reference.fetch(ref_name, ref_start, ref_end)?;

    // Calculate edit distance and mismatch quality
    let mut nm = 0; // Edit distance
    let mut uq = 0u32; // Mismatch quality sum
    let mut md_string = String::new();

    let mut ref_offset = 0; // Offset into all_ref_bases
    let mut seq_pos = 0;
    let mut match_count = 0;

    for result in cigar.iter() {
        let op = result?;
        let kind = op.kind();
        let len = op.len();

        match kind {
            Kind::Match | Kind::SequenceMatch | Kind::SequenceMismatch => {
                // Use slice of pre-fetched reference bases
                let ref_bases = &all_ref_bases[ref_offset..ref_offset + len];

                // Compare each base
                for &ref_base in ref_bases.iter() {
                    let seq_base = seq
                        .as_ref()
                        .get(seq_pos)
                        .copied()
                        .context("Sequence index out of bounds")?;
                    let qual_score = qual
                        .as_ref()
                        .get(seq_pos)
                        .copied()
                        .context("Quality index out of bounds")?;

                    if seq_base == b'N' {
                        // Masked base: count as mismatch and add to MD
                        nm += 1;
                        uq += u32::from(qual_score);

                        // Always push match count (even 0) before mismatch per SAM spec
                        md_string.push_str(&match_count.to_string());
                        match_count = 0;
                        // Preserve reference case (matches fgbio behavior)
                        md_string.push(ref_base as char);
                    } else if !seq_base.eq_ignore_ascii_case(&ref_base) {
                        // Mismatch
                        nm += 1;
                        uq += u32::from(qual_score);

                        // Always push match count (even 0) before mismatch per SAM spec
                        md_string.push_str(&match_count.to_string());
                        match_count = 0;
                        // Preserve reference case (matches fgbio behavior)
                        md_string.push(ref_base as char);
                    } else {
                        // Match
                        match_count += 1;
                    }

                    seq_pos += 1;
                }

                ref_offset += len;
            }
            Kind::Insertion => {
                // Insertion: counts toward NM, but NOT UQ (UQ only counts mismatches)
                // MD only tracks reference bases (mismatches and deletions), not insertions
                nm += len;

                // Don't sum qualities for inserted bases - UQ only counts mismatches
                // Don't push to MD - insertions are invisible in MD string
                // Just advance sequence position
                seq_pos += len;
            }
            Kind::Deletion => {
                // Deletion: counts toward NM, add to MD
                nm += len;

                // Always push match count (even 0) before deletion per SAM spec
                md_string.push_str(&match_count.to_string());
                match_count = 0;

                md_string.push('^');

                // Use slice of pre-fetched reference bases
                let ref_bases = &all_ref_bases[ref_offset..ref_offset + len];

                // Preserve reference case (matches fgbio behavior)
                for &base in ref_bases {
                    md_string.push(base as char);
                }

                ref_offset += len;
            }
            Kind::SoftClip => {
                // Soft clip: advance sequence position only
                seq_pos += len;
            }
            Kind::HardClip | Kind::Pad | Kind::Skip => {
                // These don't consume sequence or reference for tag calculation
                // Note: Skip (N in CIGAR) represents spliced alignments - the skipped
                // reference region doesn't affect NM/MD/UQ calculation
            }
        }
    }

    // Add final match count to MD (always, even if 0, per SAM spec)
    md_string.push_str(&match_count.to_string());

    // Update tags
    record.data_mut().insert(nm_tag(), Value::from(nm as i32));
    record.data_mut().insert(uq_tag(), Value::from(uq as i32));
    record.data_mut().insert(md_tag(), Value::from(md_string));

    Ok(true)
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::sam::builder::RecordBuilder;
    use noodles::sam::alignment::record::Flags;
    use noodles::sam::header::record::value::map::ReferenceSequence;
    use std::io::Write;
    use tempfile::NamedTempFile;

    fn create_test_reference() -> Result<NamedTempFile> {
        let mut file = NamedTempFile::new()?;
        writeln!(file, ">chr1")?;
        writeln!(file, "ACGTACGTACGTACGT")?;
        file.flush()?;
        Ok(file)
    }

    fn create_test_header() -> Header {
        use noodles::sam::header::record::value::Map;
        use std::num::NonZeroUsize;

        let mut header_builder = Header::builder();
        let ref_seq = Map::<ReferenceSequence>::new(NonZeroUsize::new(16).unwrap());
        header_builder = header_builder.add_reference_sequence(b"chr1", ref_seq);
        header_builder.build()
    }

    /// Helper to create a mapped record for alignment tag tests
    fn create_mapped_record(seq: &str, quals: &[u8], cigar: &str, start: usize) -> RecordBuf {
        RecordBuilder::new()
            .sequence(seq)
            .qualities(quals)
            .cigar(cigar)
            .reference_sequence_id(0)
            .alignment_start(start)
            .build()
    }

    #[test]
    fn test_perfect_match() -> Result<()> {
        let fasta = create_test_reference()?;
        let reference = ReferenceReader::new(fasta.path())?;
        let header = create_test_header();

        // Perfect match: ACGT
        let mut record = create_mapped_record("ACGT", &[30, 30, 30, 30], "4M", 1);

        let result = regenerate_alignment_tags(&mut record, &header, &reference)?;
        assert!(result); // Should successfully regenerate

        // Should have NM=0, UQ=0, MD=4
        assert_eq!(record.data().get(&nm_tag()), Some(&Value::from(0)));
        assert_eq!(record.data().get(&uq_tag()), Some(&Value::from(0)));
        assert_eq!(record.data().get(&md_tag()), Some(&Value::from("4".to_string())));

        Ok(())
    }

    #[test]
    fn test_one_mismatch() -> Result<()> {
        let fasta = create_test_reference()?;
        let reference = ReferenceReader::new(fasta.path())?;
        let header = create_test_header();

        // Mismatch at position 2: ATGT vs ACGT
        let mut record = create_mapped_record("ATGT", &[30, 30, 30, 30], "4M", 1);

        regenerate_alignment_tags(&mut record, &header, &reference)?;

        // Should have NM=1, UQ=30, MD=1C2
        assert_eq!(record.data().get(&nm_tag()), Some(&Value::from(1)));
        assert_eq!(record.data().get(&uq_tag()), Some(&Value::from(30)));
        assert_eq!(record.data().get(&md_tag()), Some(&Value::from("1C2".to_string())));

        Ok(())
    }

    #[test]
    fn test_masked_base() -> Result<()> {
        let fasta = create_test_reference()?;
        let reference = ReferenceReader::new(fasta.path())?;
        let header = create_test_header();

        // Masked base at position 2: ANGT vs ACGT
        let mut record = create_mapped_record("ANGT", &[30, 0, 30, 30], "4M", 1);

        regenerate_alignment_tags(&mut record, &header, &reference)?;

        // Masked base counts as mismatch: NM=1, UQ=0, MD=1C2
        assert_eq!(record.data().get(&nm_tag()), Some(&Value::from(1)));
        assert_eq!(record.data().get(&uq_tag()), Some(&Value::from(0)));
        assert_eq!(record.data().get(&md_tag()), Some(&Value::from("1C2".to_string())));

        Ok(())
    }

    #[test]
    fn test_unmapped_read() -> Result<()> {
        let fasta = create_test_reference()?;
        let reference = ReferenceReader::new(fasta.path())?;
        let header = create_test_header();

        let mut record = RecordBuilder::new()
            .sequence("ACGT")
            .flags(Flags::UNMAPPED)
            .tag("NM", 7i32)
            .tag("MD", "6A7C8T9G")
            .tag("UQ", 237i32)
            .build();

        let result = regenerate_alignment_tags(&mut record, &header, &reference)?;
        assert!(!result); // Should return false for unmapped

        // Verify tags are nulled out (matching fgbio behavior)
        assert!(record.data().get(&nm_tag()).is_none());
        assert!(record.data().get(&md_tag()).is_none());
        assert!(record.data().get(&uq_tag()).is_none());

        Ok(())
    }

    /// Test case matching fgbio's `BamsTest` "regenerate tags on mapped reads"
    /// Uses an all-A reference and a sequence with mismatches
    #[test]
    fn test_regenerate_tags_fgbio_equivalent() -> Result<()> {
        // Create a reference with all A's (like fgbio's DummyRefWalker)
        let mut file = NamedTempFile::new()?;
        writeln!(file, ">chr1")?;
        writeln!(file, "AAAAAAAAAAAAAAAAAAAA")?; // 20 A's
        file.flush()?;
        let reference = ReferenceReader::new(file.path())?;
        let header = create_test_header();

        // Sequence "AAACAAAATA" - mismatches at positions 4 (C vs A) and 9 (T vs A)
        // Pre-populate with wrong values (like fgbio test)
        let mut record = RecordBuilder::new()
            .sequence("AAACAAAATA")
            .qualities(&[20u8; 10])
            .cigar("10M")
            .reference_sequence_id(0)
            .alignment_start(1)
            .tag("NM", 7i32)
            .tag("MD", "6A7C8T9G")
            .tag("UQ", 237i32)
            .build();

        regenerate_alignment_tags(&mut record, &header, &reference)?;

        // Should match fgbio: NM=2, MD="3A4A1", UQ=40
        assert_eq!(record.data().get(&nm_tag()), Some(&Value::from(2)));
        assert_eq!(record.data().get(&md_tag()), Some(&Value::from("3A4A1".to_string())));
        assert_eq!(record.data().get(&uq_tag()), Some(&Value::from(40)));

        Ok(())
    }

    #[test]
    fn test_insertion() -> Result<()> {
        let fasta = create_test_reference()?;
        let reference = ReferenceReader::new(fasta.path())?;
        let header = create_test_header();

        // Alignment with insertion: 2M2I2M vs AC--GT
        let mut record = create_mapped_record("ACTTGT", &[30, 30, 25, 25, 30, 30], "2M2I2M", 1);

        regenerate_alignment_tags(&mut record, &header, &reference)?;

        // NM=2 (insertion), UQ=0 (insertions don't contribute to UQ), MD=4 (insertions invisible, 4 consecutive matches)
        assert_eq!(record.data().get(&nm_tag()), Some(&Value::from(2)));
        assert_eq!(record.data().get(&uq_tag()), Some(&Value::from(0)));
        assert_eq!(record.data().get(&md_tag()), Some(&Value::from("4".to_string())));

        Ok(())
    }

    #[test]
    fn test_deletion() -> Result<()> {
        let fasta = create_test_reference()?;
        let reference = ReferenceReader::new(fasta.path())?;
        let header = create_test_header();

        // Alignment with deletion: 2M2D2M vs ACGTAC but read is AC--AC
        let mut record = create_mapped_record("ACAC", &[30, 30, 30, 30], "2M2D2M", 1);

        regenerate_alignment_tags(&mut record, &header, &reference)?;

        // NM=2 (deletion), UQ=0, MD=2^GT2
        assert_eq!(record.data().get(&nm_tag()), Some(&Value::from(2)));
        assert_eq!(record.data().get(&uq_tag()), Some(&Value::from(0)));
        assert_eq!(record.data().get(&md_tag()), Some(&Value::from("2^GT2".to_string())));

        Ok(())
    }

    #[test]
    fn test_soft_clip() -> Result<()> {
        let fasta = create_test_reference()?;
        let reference = ReferenceReader::new(fasta.path())?;
        let header = create_test_header();

        // 2S4M2S: soft clips don't affect tags
        let mut record =
            create_mapped_record("TTACGTGG", &[20, 20, 30, 30, 30, 30, 20, 20], "2S4M2S", 1);

        regenerate_alignment_tags(&mut record, &header, &reference)?;

        // Soft clips ignored: NM=0, UQ=0, MD=4
        assert_eq!(record.data().get(&nm_tag()), Some(&Value::from(0)));
        assert_eq!(record.data().get(&uq_tag()), Some(&Value::from(0)));
        assert_eq!(record.data().get(&md_tag()), Some(&Value::from("4".to_string())));

        Ok(())
    }

    #[test]
    fn test_hard_clip() -> Result<()> {
        let fasta = create_test_reference()?;
        let reference = ReferenceReader::new(fasta.path())?;
        let header = create_test_header();

        // 2H4M2H: hard clips don't affect tags
        let mut record = create_mapped_record("ACGT", &[30, 30, 30, 30], "2H4M2H", 1);

        regenerate_alignment_tags(&mut record, &header, &reference)?;

        // Hard clips ignored: NM=0, UQ=0, MD=4
        assert_eq!(record.data().get(&nm_tag()), Some(&Value::from(0)));
        assert_eq!(record.data().get(&uq_tag()), Some(&Value::from(0)));
        assert_eq!(record.data().get(&md_tag()), Some(&Value::from("4".to_string())));

        Ok(())
    }

    #[test]
    fn test_multiple_mismatches() -> Result<()> {
        let fasta = create_test_reference()?;
        let reference = ReferenceReader::new(fasta.path())?;
        let header = create_test_header();

        // Multiple mismatches: AATT vs ACGT
        let mut record = create_mapped_record("AATT", &[30, 25, 20, 35], "4M", 1);

        regenerate_alignment_tags(&mut record, &header, &reference)?;

        // NM=2, UQ=45, MD=1C0G1 (0 between consecutive mismatches)
        assert_eq!(record.data().get(&nm_tag()), Some(&Value::from(2)));
        assert_eq!(record.data().get(&uq_tag()), Some(&Value::from(45)));
        assert_eq!(record.data().get(&md_tag()), Some(&Value::from("1C0G1".to_string())));

        Ok(())
    }

    #[test]
    fn test_multiple_masked_bases() -> Result<()> {
        let fasta = create_test_reference()?;
        let reference = ReferenceReader::new(fasta.path())?;
        let header = create_test_header();

        // Multiple masked bases: ANNN vs ACGT
        let mut record = create_mapped_record("ANNN", &[30, 0, 0, 0], "4M", 1);

        regenerate_alignment_tags(&mut record, &header, &reference)?;

        // NM=3, UQ=0, MD=1C0G0T0 (0 between consecutive mismatches, ends with 0)
        assert_eq!(record.data().get(&nm_tag()), Some(&Value::from(3)));
        assert_eq!(record.data().get(&uq_tag()), Some(&Value::from(0)));
        assert_eq!(record.data().get(&md_tag()), Some(&Value::from("1C0G0T0".to_string())));

        Ok(())
    }

    #[test]
    fn test_complex_cigar() -> Result<()> {
        let fasta = create_test_reference()?;
        let reference = ReferenceReader::new(fasta.path())?;
        let header = create_test_header();

        // Complex: 2M1I1M1D2M
        // Ref: ACGTACGT
        // Read: ACTCAC with 1 insertion (T) and 1 deletion (G->)
        let mut record = create_mapped_record("ACTCAC", &[30, 30, 25, 30, 30, 30], "2M1I1M1D2M", 1);

        regenerate_alignment_tags(&mut record, &header, &reference)?;

        // NM=3 (1 insertion + 1 mismatch + 1 deletion), UQ=30 (only mismatch qual, insertions excluded), MD=2G0^T2 (0 between mismatch and deletion)
        assert_eq!(record.data().get(&nm_tag()), Some(&Value::from(3)));
        assert_eq!(record.data().get(&uq_tag()), Some(&Value::from(30)));
        assert_eq!(record.data().get(&md_tag()), Some(&Value::from("2G0^T2".to_string())));

        Ok(())
    }

    #[test]
    fn test_sequence_match_and_mismatch_ops() -> Result<()> {
        let fasta = create_test_reference()?;
        let reference = ReferenceReader::new(fasta.path())?;
        let header = create_test_header();

        // Use explicit sequence match/mismatch ops: 2=1X1=
        let mut record = create_mapped_record("ACTT", &[30, 30, 25, 30], "2=1X1=", 1);

        regenerate_alignment_tags(&mut record, &header, &reference)?;

        // NM=1, UQ=25, MD=2G1
        assert_eq!(record.data().get(&nm_tag()), Some(&Value::from(1)));
        assert_eq!(record.data().get(&uq_tag()), Some(&Value::from(25)));
        assert_eq!(record.data().get(&md_tag()), Some(&Value::from("2G1".to_string())));

        Ok(())
    }

    #[test]
    fn test_pad_operation() -> Result<()> {
        let fasta = create_test_reference()?;
        let reference = ReferenceReader::new(fasta.path())?;
        let header = create_test_header();

        // CIGAR with Pad operation: 2M2P2M (Pad doesn't affect sequence or reference)
        let mut record = create_mapped_record("ACGT", &[30, 30, 30, 30], "2M2P2M", 1);

        regenerate_alignment_tags(&mut record, &header, &reference)?;

        // Pad ignored: NM=0, UQ=0, MD=4
        assert_eq!(record.data().get(&nm_tag()), Some(&Value::from(0)));
        assert_eq!(record.data().get(&uq_tag()), Some(&Value::from(0)));
        assert_eq!(record.data().get(&md_tag()), Some(&Value::from("4".to_string())));

        Ok(())
    }

    #[test]
    fn test_skip_operation() -> Result<()> {
        let fasta = create_test_reference()?;
        let reference = ReferenceReader::new(fasta.path())?;
        let header = create_test_header();

        // CIGAR with Skip operation (N): 2M2N2M (spliced alignment)
        let mut record = create_mapped_record("ACGT", &[30, 30, 30, 30], "2M2N2M", 1);

        regenerate_alignment_tags(&mut record, &header, &reference)?;

        // Skip ignored: NM=0, UQ=0, MD=4
        assert_eq!(record.data().get(&nm_tag()), Some(&Value::from(0)));
        assert_eq!(record.data().get(&uq_tag()), Some(&Value::from(0)));
        assert_eq!(record.data().get(&md_tag()), Some(&Value::from("4".to_string())));

        Ok(())
    }

    #[test]
    fn test_case_insensitive_matching() -> Result<()> {
        let fasta = create_test_reference()?;
        let reference = ReferenceReader::new(fasta.path())?;
        let header = create_test_header();

        // Lowercase sequence should match uppercase reference
        let mut record = create_mapped_record("acgt", &[30, 30, 30, 30], "4M", 1);

        regenerate_alignment_tags(&mut record, &header, &reference)?;

        // Should match: NM=0, UQ=0, MD=4
        assert_eq!(record.data().get(&nm_tag()), Some(&Value::from(0)));
        assert_eq!(record.data().get(&uq_tag()), Some(&Value::from(0)));
        assert_eq!(record.data().get(&md_tag()), Some(&Value::from("4".to_string())));

        Ok(())
    }

    #[test]
    fn test_insertion_at_end() -> Result<()> {
        let fasta = create_test_reference()?;
        let reference = ReferenceReader::new(fasta.path())?;
        let header = create_test_header();

        // Insertion at end: 4M2I
        let mut record = create_mapped_record("ACGTTT", &[30, 30, 30, 30, 20, 20], "4M2I", 1);

        regenerate_alignment_tags(&mut record, &header, &reference)?;

        // NM=2, UQ=0 (insertions don't contribute to UQ), MD=4
        assert_eq!(record.data().get(&nm_tag()), Some(&Value::from(2)));
        assert_eq!(record.data().get(&uq_tag()), Some(&Value::from(0)));
        assert_eq!(record.data().get(&md_tag()), Some(&Value::from("4".to_string())));

        Ok(())
    }

    #[test]
    fn test_deletion_at_end() -> Result<()> {
        let fasta = create_test_reference()?;
        let reference = ReferenceReader::new(fasta.path())?;
        let header = create_test_header();

        // Deletion at end: 2M2D
        let mut record = create_mapped_record("AC", &[30, 30], "2M2D", 1);

        regenerate_alignment_tags(&mut record, &header, &reference)?;

        // NM=2, UQ=0, MD=2^GT0 (always ends with match count)
        assert_eq!(record.data().get(&nm_tag()), Some(&Value::from(2)));
        assert_eq!(record.data().get(&uq_tag()), Some(&Value::from(0)));
        assert_eq!(record.data().get(&md_tag()), Some(&Value::from("2^GT0".to_string())));

        Ok(())
    }

    #[test]
    fn test_mixed_matches_and_masks() -> Result<()> {
        let fasta = create_test_reference()?;
        let reference = ReferenceReader::new(fasta.path())?;
        let header = create_test_header();

        // Mixed: ACNTTC vs ACGTAC (N at pos 3, A->T at pos 4)
        let mut record = create_mapped_record("ACNTTC", &[30, 30, 0, 30, 25, 30], "6M", 1);

        regenerate_alignment_tags(&mut record, &header, &reference)?;

        // NM=2 (masked + 1 mismatch), UQ=25, MD=2G1A1
        assert_eq!(record.data().get(&nm_tag()), Some(&Value::from(2)));
        assert_eq!(record.data().get(&uq_tag()), Some(&Value::from(25)));
        assert_eq!(record.data().get(&md_tag()), Some(&Value::from("2G1A1".to_string())));

        Ok(())
    }
}
