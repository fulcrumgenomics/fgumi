//! Alignment tag regeneration (NM, UQ, MD) after base masking.
//!
//! When consensus reads have bases masked to 'N', alignment tags need to be recalculated
//! to accurately reflect the new sequence. This module provides functions to regenerate:
//!
//! - **NM**: Edit distance to the reference (mismatches + indels)
//! - **UQ**: Phred likelihood of the segment (sum of mismatch qualities)
//! - **MD**: Mismatched and deleted reference bases

use std::fmt::Write;

use anyhow::{Context, Result};
use noodles::core::Position;
use noodles::sam::Header;
use noodles::sam::alignment::record::cigar::Cigar as CigarTrait;
use noodles::sam::alignment::record::cigar::op::Kind;
use noodles::sam::alignment::record::data::field::Tag;
use noodles::sam::alignment::record_buf::RecordBuf;
use noodles::sam::alignment::record_buf::data::field::Value;

use crate::ReferenceProvider;
use crate::SamTag;

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
/// For unmapped reads, the tags are removed (set to null) to match fgbio behavior, as they are
/// for a read flagged mapped but with no reference sequence id, which has nothing to recompute
/// against. For mapped reads, the tags are recalculated based on the alignment and reference.
/// A mapped read with no bases (`SEQ` `*`) has nothing to recompute from, so its tags are
/// left as they are; fgbio's `Bams.regenerateNmUqMdTags` fails on such a read instead. A read
/// with bases but no qualities (`QUAL` `*`) gets NM and MD, and keeps its UQ unchanged, as in
/// fgbio. The checks run in the same order as in [`regenerate_alignment_tags_raw_with_scoring`],
/// so both forms give the same result for the same record.
///
/// Every base differing from the reference counts ([`ConversionScoring::Literal`]); this
/// record-level form has no methylation-aware scoring. Use
/// [`regenerate_alignment_tags_raw_with_scoring`] to hide conversions from NM/UQ.
///
/// # Arguments
/// * `record` - The record to regenerate tags for (modified in place)
/// * `header` - SAM header (needed to resolve reference sequence names)
/// * `reference` - Reference genome provider
///
/// # Returns
/// True if tags were regenerated, false if the read is unmapped or has no reference sequence
/// id (tags are nulled) or has no bases (tags are left unchanged)
///
/// # Errors
///
/// Returns an error if the reference sequence ID is not found in the header, the alignment
/// start is missing, or the reference bases cannot be fetched. These checks run before the
/// no-bases check, so a read with no bases still fails on them.
#[allow(clippy::too_many_lines, clippy::cast_possible_truncation, clippy::cast_possible_wrap)]
pub fn regenerate_alignment_tags(
    record: &mut RecordBuf,
    header: &Header,
    reference: &impl ReferenceProvider,
) -> Result<bool> {
    // For unmapped reads, null out the tags (matching fgbio behavior)
    if record.flags().is_unmapped() {
        record.data_mut().remove(&nm_tag());
        record.data_mut().remove(&uq_tag());
        record.data_mut().remove(&md_tag());
        return Ok(false);
    }

    // A read flagged mapped but with no reference is malformed: strip its tags, as the raw form
    // does for a negative reference id.
    let Some(ref_seq_id) = record.reference_sequence_id() else {
        record.data_mut().remove(&nm_tag());
        record.data_mut().remove(&uq_tag());
        record.data_mut().remove(&md_tag());
        return Ok(false);
    };
    let ref_seqs = header.reference_sequences();
    let (ref_name_bytes, _) =
        ref_seqs.get_index(ref_seq_id).context("Reference sequence ID not found in header")?;
    let ref_name = std::str::from_utf8(ref_name_bytes.as_ref())?;

    let ref_start = record.alignment_start().context("Missing alignment start")?;

    // A mapped read with no bases: leave its alignment tags as they are. This sits after the
    // header and position checks so it never hides those errors.
    let seq = record.sequence();
    if seq.is_empty() {
        return Ok(false);
    }

    let cigar = record.cigar();
    let qual = record.quality_scores();
    // A read with bases but no qualities (`QUAL` `*`) has no UQ to compute; keep the one it has.
    let has_quals = !qual.as_ref().is_empty();

    // Calculate total reference span from CIGAR and fetch entire alignment span once
    // This is a key optimization - instead of fetching per CIGAR operation, we fetch once
    // Skip (N) is included because it advances the reference position
    let ref_span: usize = cigar
        .iter()
        .filter_map(Result::ok)
        .map(|op| match op.kind() {
            Kind::Match
            | Kind::SequenceMatch
            | Kind::SequenceMismatch
            | Kind::Deletion
            | Kind::Skip => op.len(),
            _ => 0,
        })
        .sum();

    // Handle edge case: CIGAR with no reference-consuming operations (e.g., pure insertion "4I")
    // Set tags to sensible defaults and return early
    if ref_span == 0 {
        record.data_mut().insert(nm_tag(), Value::from(0u32));
        record.data_mut().insert(md_tag(), Value::String("0".to_owned().into()));
        if has_quals {
            record.data_mut().insert(uq_tag(), Value::from(0u32));
        }
        return Ok(true);
    }

    // Fetch entire reference span in one call. `fetch_borrowed` lets an in-memory
    // provider hand back a slice instead of allocating a fresh `Vec` per record;
    // this runs once per record for every `fgumi clip` / `fgumi filter` invocation.
    let ref_end = Position::new(usize::from(ref_start) + ref_span - 1)
        .context("Invalid reference end position")?;
    let all_ref_bases = reference.fetch_borrowed(ref_name, ref_start, ref_end)?;

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
                for &ref_base in ref_bases {
                    let seq_base = seq
                        .as_ref()
                        .get(seq_pos)
                        .copied()
                        .context("Sequence index out of bounds")?;
                    let qual_score = if has_quals {
                        qual.as_ref()
                            .get(seq_pos)
                            .copied()
                            .context("Quality index out of bounds")?
                    } else {
                        0
                    };

                    if seq_base == b'N' {
                        // Masked base: count as mismatch and add to MD
                        nm += 1;
                        uq += u32::from(qual_score);

                        // Always push match count (even 0) before mismatch per SAM spec
                        write!(md_string, "{match_count}").expect("write to String is infallible");
                        match_count = 0;
                        // Preserve reference case (matches fgbio behavior)
                        md_string.push(ref_base as char);
                    } else if !seq_base.eq_ignore_ascii_case(&ref_base) {
                        // Mismatch
                        nm += 1;
                        uq += u32::from(qual_score);

                        // Always push match count (even 0) before mismatch per SAM spec
                        write!(md_string, "{match_count}").expect("write to String is infallible");
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
                write!(md_string, "{match_count}").expect("write to String is infallible");
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
            Kind::Skip => {
                // Skip (N) advances reference position but doesn't affect NM/MD/UQ
                ref_offset += len;
            }
            Kind::HardClip | Kind::Pad => {
                // These don't consume sequence or reference
            }
        }
    }

    // Add final match count to MD (always, even if 0, per SAM spec)
    write!(md_string, "{match_count}").expect("write to String is infallible");

    // Update tags
    record.data_mut().insert(nm_tag(), Value::from(nm as i32));
    if has_quals {
        record.data_mut().insert(uq_tag(), Value::from(uq.min(i32::MAX as u32) as i32));
    }
    record.data_mut().insert(md_tag(), Value::from(md_string));

    Ok(true)
}

// ============================================================================
// Raw-byte alignment tag regeneration
// ============================================================================

use fgumi_raw_bam::{self, RawRecordView, RawTagsEditor};

/// How methylation conversions (EM-seq / bisulfite / TAPs) are scored when regenerating
/// NM and UQ. MD is always SAM-literal.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum ConversionScoring {
    /// Every base differing from the reference is a mismatch (SAM-literal).
    Literal,
    /// A base showing the conversion of the read's original strand is not counted in NM/UQ,
    /// so a well-converted read does not look divergent. Assumes a directional library:
    /// OT-derived reads (R1 or fragment forward, R2 reverse) hide `ref C × read T`; OB-derived
    /// reads (R1 or fragment reverse, R2 forward) hide `ref G × read A`, all in genomic
    /// orientation. The opposite direction, other substitutions, `N` and indels still count.
    /// MD still lists every hidden base, so CIGAR + SEQ + MD reconstruct the reference.
    Hidden,
}

/// Whether a read of a directional bisulfite/EM-seq/TAPs library derives from the original top
/// strand (OT): R1 (or a fragment) aligned forward, or R2 aligned reverse.
///
/// R1/fragment reads carry the original strand's own sequence and R2 reads its complement, so a
/// read is OT-derived exactly when it is R2 aligned reverse or R1/fragment aligned forward. Its
/// conversions then show as `C`→`T` in genomic orientation; OB-derived reads show `G`→`A`.
#[must_use]
pub fn is_top_strand(flags: u16) -> bool {
    let is_last_segment = flags & fgumi_raw_bam::flags::LAST_SEGMENT != 0;
    let is_reverse = flags & fgumi_raw_bam::flags::REVERSE != 0;
    is_last_segment == is_reverse
}

/// Returns the `(reference, read)` base pair hidden by `scoring` for a record with `flags`, in
/// genomic orientation, or `None` when every difference counts.
fn hidden_conversion(scoring: ConversionScoring, flags: u16) -> Option<(u8, u8)> {
    match scoring {
        ConversionScoring::Literal => None,
        ConversionScoring::Hidden => {
            if is_top_strand(flags) {
                Some((b'C', b'T'))
            } else {
                Some((b'G', b'A'))
            }
        }
    }
}

/// Regenerates NM, UQ, and MD tags for a raw BAM record after base masking.
///
/// For unmapped reads, the tags are removed. For mapped reads, the tags are
/// recalculated based on the alignment and reference. A record flagged mapped
/// but carrying a negative reference id is malformed and has nothing to
/// recompute against, so its tags are removed as well rather than left stale.
/// A mapped record with no bases (`SEQ` `*`, e.g. a secondary alignment as
/// `bwa mem -a` writes it) has nothing to recompute from, so its tags are left
/// as they are: recomputing them would walk a CIGAR that consumes bases the
/// record does not carry. fgbio's `Bams.regenerateNmUqMdTags` fails on such a
/// record instead. This check runs after the header and position checks, so it
/// never hides their errors. A record with bases but no qualities (`QUAL` `*`,
/// stored as `0xFF` filler) gets NM and MD and keeps its UQ unchanged, as in
/// fgbio, rather than summing the filler bytes. The typed
/// [`regenerate_alignment_tags`] runs the same checks in the same order.
///
/// Returns `Ok(true)` if tags were regenerated, `Ok(false)` if they were removed
/// (an unmapped read, or a mapped read with no reference id) or left unchanged
/// (a mapped read with no bases).
///
/// Every base differing from the reference counts ([`ConversionScoring::Literal`]); see
/// [`regenerate_alignment_tags_raw_with_scoring`] to hide methylation conversions from NM/UQ.
///
/// # Errors
///
/// Returns an error if the record is too short, the reference sequence ID is not found
/// in the header, the alignment start is invalid, or the CIGAR operations reference
/// beyond the available sequence or reference data.
pub fn regenerate_alignment_tags_raw(
    record: &mut Vec<u8>,
    header: &Header,
    reference: &impl ReferenceProvider,
) -> Result<bool> {
    regenerate_alignment_tags_raw_with_scoring(
        record,
        header,
        reference,
        ConversionScoring::Literal,
    )
}

/// Regenerates NM, UQ, and MD tags for a raw BAM record, scoring methylation conversions in
/// NM/UQ per `scoring` (MD is always SAM-literal). Otherwise identical to [`regenerate_alignment_tags_raw`].
///
/// # Errors
///
/// See [`regenerate_alignment_tags_raw`].
#[allow(
    clippy::too_many_lines,
    clippy::cast_sign_loss,
    clippy::cast_possible_truncation,
    clippy::cast_possible_wrap
)]
pub fn regenerate_alignment_tags_raw_with_scoring(
    record: &mut Vec<u8>,
    header: &Header,
    reference: &impl ReferenceProvider,
    scoring: ConversionScoring,
) -> Result<bool> {
    if record.len() < fgumi_raw_bam::MIN_BAM_RECORD_LEN {
        anyhow::bail!(
            "BAM record too short ({} bytes, minimum {})",
            record.len(),
            fgumi_raw_bam::MIN_BAM_RECORD_LEN
        );
    }
    // For unmapped reads, remove alignment tags
    if RawRecordView::new(record).is_unmapped() {
        // Strip all three alignment tags in a single aux pass.
        RawTagsEditor::from_vec(record)
            .rebuild_with(&[SamTag::NM.into(), SamTag::UQ.into(), SamTag::MD.into()], &[]);
        return Ok(false);
    }

    // Get reference sequence ID and look up name in header
    let ref_seq_id = RawRecordView::new(record).ref_id();
    if ref_seq_id < 0 {
        // A record flagged mapped but carrying no reference is malformed: there
        // is nothing to recompute against. Strip the alignment tags rather than
        // leaving whatever stale NM/UQ/MD the record arrived with, matching the
        // unmapped branch above and fgbio's `regenerateNmUqMdTags`. Failing
        // closed here means a downstream consumer sees no tag instead of a
        // wrong one.
        // Strip all three alignment tags in a single aux pass.
        RawTagsEditor::from_vec(record)
            .rebuild_with(&[SamTag::NM.into(), SamTag::UQ.into(), SamTag::MD.into()], &[]);
        return Ok(false);
    }
    let ref_seqs = header.reference_sequences();
    let (ref_name_bytes, _) = ref_seqs
        .get_index(ref_seq_id as usize)
        .context("Reference sequence ID not found in header")?;
    let ref_name = std::str::from_utf8(ref_name_bytes.as_ref())?;

    let alignment_start_0based = RawRecordView::new(record).pos();
    if alignment_start_0based < 0 {
        anyhow::bail!("Invalid alignment start position: {alignment_start_0based}");
    }
    let ref_start = Position::new((alignment_start_0based + 1) as usize)
        .context("Invalid alignment start position")?;

    // A mapped record with no bases: leave its alignment tags as they are. This sits after the
    // header and position checks so it never hides those errors.
    let l_seq = fgumi_raw_bam::l_seq(record) as usize;
    if l_seq == 0 {
        return Ok(false);
    }
    // BAM stores a missing QUAL (`*`) as 0xFF filler. htsjdk treats a record whose first quality
    // byte is 0xFF as having no qualities (`BAMRecord.decodeBaseQualities`), and fgbio then leaves
    // UQ unchanged (`Bams.regenerateNmUqMdTags`), so do the same instead of summing the filler.
    let qual_off = fgumi_raw_bam::qual_offset(record);
    let has_quals = record.get(qual_off).is_some_and(|&q| q != 0xFF);

    // Calculate reference span directly from the raw CIGAR bytes (zero allocation).
    let ref_span = usize::try_from(fgumi_raw_bam::reference_length_from_raw_bam(record))
        .context("CIGAR-derived reference span is negative")?;

    // Handle edge case: CIGAR with no reference-consuming operations
    if ref_span == 0 {
        let mut editor = RawTagsEditor::from_vec(record);
        editor.update_int(SamTag::NM, 0);
        if has_quals {
            editor.update_int(SamTag::UQ, 0);
        }
        editor.update_string(SamTag::MD, b"0");
        return Ok(true);
    }

    // Fetch entire reference span in one call. `fetch_borrowed` lets an in-memory
    // provider hand back a slice instead of allocating a fresh `Vec` per record;
    // this runs once per record for every `fgumi clip` / `fgumi filter` invocation.
    let ref_end = Position::new(usize::from(ref_start) + ref_span - 1)
        .context("Invalid reference end position")?;
    let all_ref_bases = reference.fetch_borrowed(ref_name, ref_start, ref_end)?;

    // Get seq/qual offsets and validate bounds
    let seq_off = fgumi_raw_bam::seq_offset(record);
    let seq_bytes = l_seq.div_ceil(2);
    if seq_off + seq_bytes > record.len() || qual_off + l_seq > record.len() {
        anyhow::bail!("Truncated BAM record: seq/qual extends past record end");
    }

    let hidden = hidden_conversion(scoring, RawRecordView::new(record).flags());

    // Calculate NM, UQ, MD
    let mut nm: i32 = 0;
    let mut uq: u32 = 0;
    let mut md_string = String::new();
    let mut ref_offset = 0;
    let mut seq_pos = 0;
    let mut match_count: usize = 0;

    for op in RawRecordView::new(record).cigar_ops_iter() {
        let op_type = op & 0xF;
        let op_len = (op >> 4) as usize;

        match op_type {
            0 | 7 | 8 => {
                // M (0), = (7), X (8) — alignment match/mismatch
                if ref_offset + op_len > all_ref_bases.len() {
                    anyhow::bail!("CIGAR references beyond fetched reference span");
                }
                if seq_pos + op_len > l_seq {
                    anyhow::bail!("CIGAR consumes more bases than sequence length");
                }
                let ref_bases = &all_ref_bases[ref_offset..ref_offset + op_len];
                for &ref_base in ref_bases {
                    let seq_base = fgumi_raw_bam::BAM_BASE_TO_ASCII
                        [fgumi_raw_bam::get_base(record, seq_off, seq_pos) as usize];
                    let qual_score = fgumi_raw_bam::get_qual(record, qual_off, seq_pos);

                    if seq_base == b'N' {
                        // Masked base: count as mismatch
                        nm += 1;
                        uq += u32::from(qual_score);
                        write!(md_string, "{match_count}").expect("write to String is infallible");
                        match_count = 0;
                        md_string.push(ref_base as char);
                    } else if !seq_base.eq_ignore_ascii_case(&ref_base) {
                        // Mismatch. MD always lists it so that SEQ + CIGAR + MD reconstruct the
                        // reference; NM/UQ skip a hidden conversion.
                        if hidden != Some((ref_base.to_ascii_uppercase(), seq_base)) {
                            nm += 1;
                            uq += u32::from(qual_score);
                        }
                        write!(md_string, "{match_count}").expect("write to String is infallible");
                        match_count = 0;
                        md_string.push(ref_base as char);
                    } else {
                        match_count += 1;
                    }
                    seq_pos += 1;
                }
                ref_offset += op_len;
            }
            1 => {
                // I — insertion
                if seq_pos + op_len > l_seq {
                    anyhow::bail!("CIGAR insertion consumes more bases than sequence length");
                }
                nm += op_len as i32;
                seq_pos += op_len;
            }
            2 => {
                // D — deletion
                if ref_offset + op_len > all_ref_bases.len() {
                    anyhow::bail!("CIGAR deletion references beyond fetched reference span");
                }
                nm += op_len as i32;
                write!(md_string, "{match_count}").expect("write to String is infallible");
                match_count = 0;
                md_string.push('^');
                let ref_bases = &all_ref_bases[ref_offset..ref_offset + op_len];
                for &base in ref_bases {
                    md_string.push(base as char);
                }
                ref_offset += op_len;
            }
            4 => {
                // S — soft clip
                if seq_pos + op_len > l_seq {
                    anyhow::bail!("CIGAR soft clip consumes more bases than sequence length");
                }
                seq_pos += op_len;
            }
            3 => {
                // N — skip (spliced alignment), advances reference but not NM/MD/UQ
                ref_offset += op_len;
            }
            // H (5), P (6), and others — no sequence or ref consumed for tag calc
            _ => {}
        }
    }

    // Add final match count
    write!(md_string, "{match_count}").expect("write to String is infallible");

    // Update tags
    let mut editor = RawTagsEditor::from_vec(record);
    editor.update_int(SamTag::NM, nm);
    if has_quals {
        editor.update_int(SamTag::UQ, uq.min(i32::MAX as u32) as i32);
    }
    editor.update_string(SamTag::MD, md_string.as_bytes());

    Ok(true)
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::builder::RecordBuilder;
    use noodles::sam::alignment::record::Flags;
    use noodles::sam::alignment::record_buf::{QualityScores, Sequence};
    use noodles::sam::header::record::value::map::ReferenceSequence;
    use std::collections::HashMap;
    use std::io::Write;
    use tempfile::NamedTempFile;

    /// Mock reference provider for tests
    struct MockReference {
        sequences: HashMap<String, Vec<u8>>,
    }

    impl MockReference {
        fn from_fasta(path: &std::path::Path) -> Result<Self> {
            let content = std::fs::read_to_string(path)?;
            let mut sequences = HashMap::new();
            let mut current_name = String::new();
            let mut current_seq = Vec::new();

            for line in content.lines() {
                if let Some(name) = line.strip_prefix('>') {
                    if !current_name.is_empty() {
                        sequences.insert(current_name.clone(), current_seq.clone());
                        current_seq.clear();
                    }
                    current_name = name.to_string();
                } else {
                    current_seq.extend_from_slice(line.as_bytes());
                }
            }
            if !current_name.is_empty() {
                sequences.insert(current_name, current_seq);
            }

            Ok(Self { sequences })
        }
    }

    impl ReferenceProvider for MockReference {
        fn fetch(&self, chrom: &str, start: Position, end: Position) -> Result<Vec<u8>> {
            let sequence = self
                .sequences
                .get(chrom)
                .ok_or_else(|| anyhow::anyhow!("Reference not found: {chrom}"))?;
            let start_idx = usize::from(start) - 1;
            let end_idx = usize::from(end);
            Ok(sequence[start_idx..end_idx].to_vec())
        }
    }

    fn create_test_reference() -> Result<(NamedTempFile, MockReference)> {
        let mut file = NamedTempFile::new()?;
        writeln!(file, ">chr1")?;
        writeln!(file, "ACGTACGTACGTACGT")?;
        file.flush()?;
        let reference = MockReference::from_fasta(file.path())?;
        Ok((file, reference))
    }

    fn create_test_header() -> Header {
        use noodles::sam::header::record::value::Map;
        use std::num::NonZeroUsize;

        let mut header_builder = Header::builder();
        let ref_seq = Map::<ReferenceSequence>::new(
            NonZeroUsize::new(16).expect("reference sequence length must be non-zero"),
        );
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

    /// Encode a `RecordBuf` to raw BAM bytes for testing raw tag regeneration.
    #[allow(clippy::cast_sign_loss)]
    fn encode_record_buf_to_raw(header: &Header, record: &RecordBuf) -> Result<Vec<u8>> {
        use noodles::sam::alignment::io::Write as AlignmentWrite;
        use std::io::{Cursor, Read};

        // Write a complete BAM file in memory
        let mut bam_data = Vec::new();
        {
            let mut writer = noodles::bam::io::Writer::new(&mut bam_data);
            writer.write_header(header)?;
            writer.write_alignment_record(header, record)?;
            writer.get_mut().try_finish()?;
        }

        // Read back: skip the header, then read the raw record bytes
        let mut reader = noodles::bam::io::Reader::new(Cursor::new(&bam_data));
        let _header = reader.read_header()?;

        // Each BAM record is prefixed with block_size:i32
        let mut size_buf = [0u8; 4];
        reader.get_mut().read_exact(&mut size_buf)?;
        let block_size = i32::from_le_bytes(size_buf) as usize;
        let mut record_bytes = vec![0u8; block_size];
        reader.get_mut().read_exact(&mut record_bytes)?;

        Ok(record_bytes)
    }

    #[test]
    fn test_perfect_match() -> Result<()> {
        let (_fasta, reference) = create_test_reference()?;
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
        let (_fasta, reference) = create_test_reference()?;
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
        let (_fasta, reference) = create_test_reference()?;
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
        let (_fasta, reference) = create_test_reference()?;
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
        let reference = MockReference::from_fasta(file.path())?;
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
        let (_fasta, reference) = create_test_reference()?;
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
        let (_fasta, reference) = create_test_reference()?;
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
        let (_fasta, reference) = create_test_reference()?;
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
        let (_fasta, reference) = create_test_reference()?;
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
        let (_fasta, reference) = create_test_reference()?;
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
        let (_fasta, reference) = create_test_reference()?;
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
        let (_fasta, reference) = create_test_reference()?;
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
        let (_fasta, reference) = create_test_reference()?;
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
        let (_fasta, reference) = create_test_reference()?;
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
        let (_fasta, reference) = create_test_reference()?;
        let header = create_test_header();

        // CIGAR with Skip operation (N): 2M2N2M (spliced alignment)
        // Ref: ACGTACGTACGTACGT
        // Read: AC|GT aligned as 2M at pos 1-2, skip 2 ref bases, 2M at pos 5-6
        // Pos 1-2: AC vs AC (match), Pos 5-6: GT vs AC (2 mismatches)
        let mut record = create_mapped_record("ACGT", &[30, 30, 30, 30], "2M2N2M", 1);

        regenerate_alignment_tags(&mut record, &header, &reference)?;

        // Skip advances ref but doesn't contribute to NM/MD/UQ directly;
        // the bases after the skip align to different ref positions.
        assert_eq!(record.data().get(&nm_tag()), Some(&Value::from(2)));
        assert_eq!(record.data().get(&uq_tag()), Some(&Value::from(60)));
        assert_eq!(record.data().get(&md_tag()), Some(&Value::from("2A0C0".to_string())));

        Ok(())
    }

    #[test]
    fn test_case_insensitive_matching() -> Result<()> {
        let (_fasta, reference) = create_test_reference()?;
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
        let (_fasta, reference) = create_test_reference()?;
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
        let (_fasta, reference) = create_test_reference()?;
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
        let (_fasta, reference) = create_test_reference()?;
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

    /// A record flagged mapped but carrying `ref_id < 0` is malformed: there is
    /// no reference to recompute against. It must not keep whatever stale
    /// NM/UQ/MD it arrived with, since a downstream consumer cannot tell a
    /// stale tag from a correct one. Failing closed matches both the unmapped
    /// branch and fgbio's `regenerateNmUqMdTags`.
    #[test]
    fn test_regenerate_alignment_tags_raw_strips_stale_tags_on_negative_ref_id() -> Result<()> {
        let (_fasta, reference) = create_test_reference()?;
        let header = create_test_header();
        let record = create_mapped_record("ACGTACGT", &[30, 30, 30, 30, 30, 30, 30, 30], "8M", 1);
        let mut raw = encode_record_buf_to_raw(&header, &record)?;

        // Plant stale alignment tags, then blank the reference id while leaving
        // the record flagged as mapped.
        {
            let mut editor = RawTagsEditor::from_vec(&mut raw);
            editor.update_int(SamTag::NM, 99);
            editor.update_int(SamTag::UQ, 12_345);
            editor.update_string(SamTag::MD, b"8");
        }
        assert!(
            fgumi_raw_bam::find_int_tag(fgumi_raw_bam::aux_data_slice(&raw), SamTag::NM).is_some(),
            "precondition: the stale NM tag should be present before the call",
        );

        // These bytes exclude the block_size prefix, so refID sits at 0..4.
        raw[0..4].copy_from_slice(&(-1i32).to_le_bytes());
        assert!(
            !RawRecordView::new(&raw).is_unmapped(),
            "precondition: the record must still be flagged mapped, or the \
             unmapped branch would handle it instead",
        );

        let changed = regenerate_alignment_tags_raw(&mut raw, &header, &reference)?;
        assert!(!changed, "no tags can be recomputed without a reference");
        assert!(
            !RawRecordView::new(&raw).is_unmapped(),
            "the record must remain flagged mapped after tag cleanup",
        );

        let aux = fgumi_raw_bam::aux_data_slice(&raw);
        for tag in [SamTag::NM, SamTag::UQ, SamTag::MD] {
            assert!(
                fgumi_raw_bam::find_tag_type(aux, tag).is_none(),
                "{tag:?} must be stripped rather than left stale on a negative ref_id",
            );
        }

        Ok(())
    }

    #[test]
    fn test_regenerate_alignment_tags_raw_validates_bounds() -> Result<()> {
        // Create a valid record, encode to raw, then truncate it
        let (_fasta, reference) = create_test_reference()?;
        let header = create_test_header();
        let record = create_mapped_record("ACGTACGT", &[30, 30, 30, 30, 30, 30, 30, 30], "8M", 1);

        let raw = encode_record_buf_to_raw(&header, &record)?;

        // Truncate to remove quality scores (but keep seq)
        let qual_off = fgumi_raw_bam::qual_offset(&raw);
        let mut truncated = raw[..qual_off].to_vec();

        let result = regenerate_alignment_tags_raw(&mut truncated, &header, &reference);
        assert!(result.is_err());
        assert!(result.unwrap_err().to_string().contains("Truncated"));

        Ok(())
    }

    #[test]
    fn test_regenerate_alignment_tags_raw_rejects_short_record() -> Result<()> {
        let (_fasta, reference) = create_test_reference()?;
        let header = create_test_header();

        // Record shorter than MIN_BAM_RECORD_LEN (36 bytes)
        let mut too_short = vec![0u8; 10];
        let result = regenerate_alignment_tags_raw(&mut too_short, &header, &reference);
        assert!(result.is_err());
        assert!(result.unwrap_err().to_string().contains("too short"));

        Ok(())
    }

    /// Round-trip test: encode `RecordBuf` → regenerate raw tags → verify NM/UQ/MD match `RecordBuf` path.
    #[test]
    fn test_regenerate_alignment_tags_raw_happy_path() -> Result<()> {
        let (_fasta, reference) = create_test_reference()?;
        let header = create_test_header();

        // Record with mismatches: ATGT vs ref ACGT (mismatch at pos 2)
        let mut record_buf = create_mapped_record("ATGT", &[30, 30, 25, 30], "4M", 1);
        regenerate_alignment_tags(&mut record_buf, &header, &reference)?;

        // Encode to raw bytes and regenerate tags via raw path
        let mut raw = encode_record_buf_to_raw(&header, &record_buf)?;
        regenerate_alignment_tags_raw(&mut raw, &header, &reference)?;

        // Read tags back from raw bytes
        let aux_off = fgumi_raw_bam::aux_data_offset_from_record(&raw).unwrap_or(raw.len());
        let aux = &raw[aux_off..];

        let nm = fgumi_raw_bam::find_int_tag(aux, SamTag::NM);
        let uq = fgumi_raw_bam::find_int_tag(aux, SamTag::UQ);
        let md = fgumi_raw_bam::find_string_tag(aux, SamTag::MD);

        // Verify raw path produces same results as RecordBuf path
        assert_eq!(nm, Some(1i64), "NM should be 1 (one mismatch)");
        assert_eq!(uq, Some(30i64), "UQ should be 30 (quality at mismatch position)");
        assert_eq!(
            md.map(|s| std::str::from_utf8(s).expect("MD tag should be valid UTF-8")),
            Some("1C2"),
            "MD should be 1C2"
        );

        Ok(())
    }

    /// Round-trip test with masked bases (N): verify raw path matches `RecordBuf`.
    #[test]
    fn test_regenerate_alignment_tags_raw_with_masked_bases() -> Result<()> {
        let (_fasta, reference) = create_test_reference()?;
        let header = create_test_header();

        // Record with masked base: ANGT vs ref ACGT
        let mut record_buf = create_mapped_record("ANGT", &[30, 0, 30, 30], "4M", 1);
        regenerate_alignment_tags(&mut record_buf, &header, &reference)?;

        let mut raw = encode_record_buf_to_raw(&header, &record_buf)?;
        regenerate_alignment_tags_raw(&mut raw, &header, &reference)?;

        let aux_off = fgumi_raw_bam::aux_data_offset_from_record(&raw).unwrap_or(raw.len());
        let aux = &raw[aux_off..];

        let nm = fgumi_raw_bam::find_int_tag(aux, SamTag::NM);
        let uq = fgumi_raw_bam::find_int_tag(aux, SamTag::UQ);
        let md = fgumi_raw_bam::find_string_tag(aux, SamTag::MD);

        assert_eq!(nm, Some(1i64), "NM should be 1 (masked base = mismatch)");
        assert_eq!(uq, Some(0i64), "UQ should be 0 (masked base quality is 0)");
        assert_eq!(
            md.map(|s| std::str::from_utf8(s).expect("MD tag should be valid UTF-8")),
            Some("1C2")
        );

        Ok(())
    }

    /// How a test record's segment and strand flags are set.
    #[derive(Clone, Copy, Debug)]
    enum ReadType {
        Fragment,
        R1,
        R2,
    }

    /// Hidden conversions follow the original strand a read derives from: OT-type reads
    /// (R1-fwd, R2-rev, fragment-fwd) hide `ref C × seq T`, OB-type reads (R1-rev, R2-fwd,
    /// fragment-rev) hide `ref G × seq A`. Everything else (the opposite direction, other
    /// substitutions, `N`) stays a mismatch. MD is SAM-literal in every case: it lists every
    /// difference, hidden or not, so SEQ + CIGAR + MD reconstruct the reference. Reference
    /// `ACGTACGT`, read `ATATATGA`: C→T at 2 and 6, G→A at 3, T→A at 8. Qualities are 10..=17
    /// so UQ identifies which bases counted.
    #[rstest::rstest]
    #[case::literal_ignores_read_type(
        ReadType::R1,
        false,
        ConversionScoring::Literal,
        "ATATATGA",
        4,
        55,
        "1C0G2C1T0"
    )]
    #[case::fragment_fwd_hides_c_to_t(
        ReadType::Fragment,
        false,
        ConversionScoring::Hidden,
        "ATATATGA",
        2,
        29,
        "1C0G2C1T0"
    )]
    #[case::fragment_rev_hides_g_to_a(
        ReadType::Fragment,
        true,
        ConversionScoring::Hidden,
        "ATATATGA",
        3,
        43,
        "1C0G2C1T0"
    )]
    #[case::r1_fwd_hides_c_to_t(
        ReadType::R1,
        false,
        ConversionScoring::Hidden,
        "ATATATGA",
        2,
        29,
        "1C0G2C1T0"
    )]
    #[case::r1_rev_hides_g_to_a(
        ReadType::R1,
        true,
        ConversionScoring::Hidden,
        "ATATATGA",
        3,
        43,
        "1C0G2C1T0"
    )]
    #[case::r2_fwd_hides_g_to_a(
        ReadType::R2,
        false,
        ConversionScoring::Hidden,
        "ATATATGA",
        3,
        43,
        "1C0G2C1T0"
    )]
    #[case::r2_rev_hides_c_to_t(
        ReadType::R2,
        true,
        ConversionScoring::Hidden,
        "ATATATGA",
        2,
        29,
        "1C0G2C1T0"
    )]
    #[case::masked_base_at_ref_c_still_counts(
        ReadType::R1,
        false,
        ConversionScoring::Hidden,
        "ANGTACGT",
        1,
        11,
        "1C6"
    )]
    fn test_regenerate_alignment_tags_raw_conversion_scoring(
        #[case] read_type: ReadType,
        #[case] reverse: bool,
        #[case] scoring: ConversionScoring,
        #[case] seq: &str,
        #[case] expected_nm: i64,
        #[case] expected_uq: i64,
        #[case] expected_md: &str,
    ) -> Result<()> {
        let (_fasta, reference) = create_test_reference()?;
        let header = create_test_header();
        let builder = RecordBuilder::new()
            .sequence(seq)
            .qualities(&[10, 11, 12, 13, 14, 15, 16, 17])
            .cigar("8M")
            .reference_sequence_id(0)
            .alignment_start(1)
            .reverse_complement(reverse);
        let builder = match read_type {
            ReadType::Fragment => builder,
            ReadType::R1 => builder.first_segment(true),
            ReadType::R2 => builder.first_segment(false),
        };
        let mut raw = encode_record_buf_to_raw(&header, &builder.build())?;

        regenerate_alignment_tags_raw_with_scoring(&mut raw, &header, &reference, scoring)?;

        let aux = fgumi_raw_bam::aux_data_slice(&raw);
        assert_eq!(fgumi_raw_bam::find_int_tag(aux, SamTag::NM), Some(expected_nm), "NM");
        assert_eq!(fgumi_raw_bam::find_int_tag(aux, SamTag::UQ), Some(expected_uq), "UQ");
        assert_eq!(
            fgumi_raw_bam::find_string_tag(aux, SamTag::MD)
                .map(|s| std::str::from_utf8(s).expect("MD tag should be valid UTF-8")),
            Some(expected_md),
            "MD"
        );
        Ok(())
    }

    /// A mapped `8M` record at position 1 with no bases (`SEQ` and `QUAL` `*`).
    fn empty_sequence_record() -> RecordBuf {
        // The builder fills in bases for a CIGAR-only record, so clear them afterwards.
        let mut record = create_mapped_record("ACGTACGT", &[30; 8], "8M", 1);
        *record.sequence_mut() = Sequence::default();
        *record.quality_scores_mut() = QualityScores::default();
        record
    }

    /// A mapped record with no bases (`SEQ` `*`, as `bwa mem -a` writes secondaries) keeps its
    /// NM/UQ/MD unchanged and returns `Ok(false)`: there are no bases to compare against the
    /// reference, and walking its `8M` CIGAR would run past the empty sequence.
    #[test]
    fn test_regenerate_alignment_tags_leaves_empty_sequence_tags_unchanged() -> Result<()> {
        let (_fasta, reference) = create_test_reference()?;
        let header = create_test_header();

        let mut record = empty_sequence_record();
        assert!(record.sequence().is_empty(), "precondition: the record has no bases");
        record.data_mut().insert(nm_tag(), Value::from(7u32));
        record.data_mut().insert(uq_tag(), Value::from(70u32));
        record.data_mut().insert(md_tag(), Value::from("8".to_string()));

        let regenerated = regenerate_alignment_tags(&mut record, &header, &reference)?;

        assert!(!regenerated, "nothing was regenerated");
        assert_eq!(record.data().get(&nm_tag()), Some(&Value::from(7u32)));
        assert_eq!(record.data().get(&uq_tag()), Some(&Value::from(70u32)));
        assert_eq!(record.data().get(&md_tag()), Some(&Value::from("8".to_string())));
        Ok(())
    }

    /// Raw-path counterpart of
    /// `test_regenerate_alignment_tags_leaves_empty_sequence_tags_unchanged`, for both scorings.
    #[rstest::rstest]
    fn test_regenerate_alignment_tags_raw_leaves_empty_sequence_tags_unchanged(
        #[values(ConversionScoring::Literal, ConversionScoring::Hidden)] scoring: ConversionScoring,
    ) -> Result<()> {
        let (_fasta, reference) = create_test_reference()?;
        let header = create_test_header();

        let mut raw = encode_record_buf_to_raw(&header, &empty_sequence_record())?;
        {
            let mut editor = RawTagsEditor::from_vec(&mut raw);
            editor.update_int(SamTag::NM, 7);
            editor.update_int(SamTag::UQ, 70);
            editor.update_string(SamTag::MD, b"8");
        }
        assert_eq!(fgumi_raw_bam::l_seq(&raw), 0, "precondition: the record has no bases");
        assert!(!RawRecordView::new(&raw).is_unmapped(), "precondition: the record is mapped");
        let before = raw.clone();

        let regenerated =
            regenerate_alignment_tags_raw_with_scoring(&mut raw, &header, &reference, scoring)?;

        assert!(!regenerated, "nothing was regenerated");
        assert_eq!(raw, before, "the record is left byte-for-byte unchanged");
        Ok(())
    }

    /// The bases a parity-test record carries.
    #[derive(Clone, Copy, Debug)]
    enum ParitySeq {
        /// `ACGTACGA`, one mismatch against the reference, all qualities 30.
        Bases,
        /// The same bases with no qualities (`QUAL` `*`, 0xFF filler in BAM).
        NoQuals,
        /// No bases (`SEQ` and `QUAL` `*`).
        Empty,
    }

    /// The reference id a parity-test record carries.
    #[derive(Clone, Copy, Debug)]
    enum ParityRef {
        /// `chr1`, present in the header.
        Valid,
        /// No reference id (typed `None`, raw `-1`) on a record flagged mapped.
        Missing,
        /// A reference id the header does not have.
        OutOfRange,
    }

    /// NM, UQ and MD of a record, `None` where the tag is absent.
    type AlignmentTags = (Option<i64>, Option<i64>, Option<String>);

    /// The typed and raw regenerators run their checks in the same order, so they agree on every
    /// record, including the malformed ones where several checks apply at once: no bases with a
    /// missing or out-of-range reference id. Each case pins the shared outcome: `Some((returned,
    /// tags))` on success, `None` on error. Every record starts with NM 7, UQ 70 and MD `8`.
    #[rstest::rstest]
    #[case::bases(ParitySeq::Bases, ParityRef::Valid, false,
        Some((true, (Some(1), Some(30), Some("7T0".to_string())))))]
    #[case::no_quals_keeps_uq(ParitySeq::NoQuals, ParityRef::Valid, false,
        Some((true, (Some(1), Some(70), Some("7T0".to_string())))))]
    #[case::empty_keeps_tags(ParitySeq::Empty, ParityRef::Valid, false,
        Some((false, (Some(7), Some(70), Some("8".to_string())))))]
    #[case::empty_missing_ref_strips(ParitySeq::Empty, ParityRef::Missing, false,
        Some((false, (None, None, None))))]
    #[case::bases_missing_ref_strips(ParitySeq::Bases, ParityRef::Missing, false,
        Some((false, (None, None, None))))]
    #[case::empty_out_of_range_ref_errors(ParitySeq::Empty, ParityRef::OutOfRange, false, None)]
    #[case::bases_out_of_range_ref_errors(ParitySeq::Bases, ParityRef::OutOfRange, false, None)]
    #[case::empty_unmapped_strips(ParitySeq::Empty, ParityRef::Valid, true,
        Some((false, (None, None, None))))]
    fn test_regenerate_alignment_tags_typed_and_raw_agree(
        #[case] seq: ParitySeq,
        #[case] ref_id: ParityRef,
        #[case] unmapped: bool,
        #[case] expected: Option<(bool, AlignmentTags)>,
    ) -> Result<()> {
        let (_fasta, reference) = create_test_reference()?;
        let header = create_test_header();

        let mut typed = create_mapped_record("ACGTACGA", &[30; 8], "8M", 1);
        match seq {
            ParitySeq::Bases => {}
            ParitySeq::NoQuals => *typed.quality_scores_mut() = QualityScores::default(),
            ParitySeq::Empty => {
                *typed.sequence_mut() = Sequence::default();
                *typed.quality_scores_mut() = QualityScores::default();
            }
        }
        if unmapped {
            *typed.flags_mut() = Flags::UNMAPPED;
        }
        typed.data_mut().insert(nm_tag(), Value::from(7u32));
        typed.data_mut().insert(uq_tag(), Value::from(70u32));
        typed.data_mut().insert(md_tag(), Value::from("8".to_string()));
        // Encode before breaking the reference id: the BAM writer validates it.
        let mut raw = encode_record_buf_to_raw(&header, &typed)?;
        match ref_id {
            ParityRef::Valid => {}
            ParityRef::Missing => {
                *typed.reference_sequence_id_mut() = None;
                fgumi_raw_bam::set_ref_id(&mut raw, -1);
            }
            ParityRef::OutOfRange => {
                *typed.reference_sequence_id_mut() = Some(5);
                fgumi_raw_bam::set_ref_id(&mut raw, 5);
            }
        }

        let typed_result = regenerate_alignment_tags(&mut typed, &header, &reference);
        let raw_result = regenerate_alignment_tags_raw(&mut raw, &header, &reference);

        let typed_tags: AlignmentTags = (
            typed.data().get(&nm_tag()).and_then(Value::as_int),
            typed.data().get(&uq_tag()).and_then(Value::as_int),
            match typed.data().get(&md_tag()) {
                Some(Value::String(md)) => Some(md.to_string()),
                _ => None,
            },
        );
        let aux = fgumi_raw_bam::aux_data_slice(&raw);
        let raw_tags: AlignmentTags = (
            fgumi_raw_bam::find_int_tag(aux, SamTag::NM),
            fgumi_raw_bam::find_int_tag(aux, SamTag::UQ),
            fgumi_raw_bam::find_string_tag(aux, SamTag::MD)
                .map(|md| String::from_utf8_lossy(md).into_owned()),
        );

        let typed_outcome = typed_result.ok().map(|regenerated| (regenerated, typed_tags));
        let raw_outcome = raw_result.ok().map(|regenerated| (regenerated, raw_tags));
        assert_eq!(typed_outcome, raw_outcome, "typed and raw forms agree");
        assert_eq!(typed_outcome, expected, "outcome");
        Ok(())
    }

    /// Zero-ref-span CIGAR (pure insertion) should return `Ok(true)` since NM/MD/UQ are written.
    #[test]
    fn test_regenerate_alignment_tags_zero_ref_span_returns_true() -> Result<()> {
        let (_fasta, reference) = create_test_reference()?;
        let header = create_test_header();

        // Pure-insertion CIGAR: 4I. Mapped record but ref_span = 0.
        let mut record = create_mapped_record("ACGT", &[30, 30, 30, 30], "4I", 1);
        let regenerated = regenerate_alignment_tags(&mut record, &header, &reference)?;

        assert!(regenerated, "tags were written for zero-ref-span record; return must be true");
        assert_eq!(record.data().get(&nm_tag()), Some(&Value::from(0u32)));
        assert_eq!(record.data().get(&uq_tag()), Some(&Value::from(0u32)));
        assert_eq!(record.data().get(&md_tag()), Some(&Value::from("0".to_string())));

        Ok(())
    }

    /// Raw-path zero-ref-span CIGAR (pure insertion) should return `Ok(true)` since NM/MD/UQ
    /// are written.
    #[test]
    fn test_regenerate_alignment_tags_raw_zero_ref_span_returns_true() -> Result<()> {
        let (_fasta, reference) = create_test_reference()?;
        let header = create_test_header();

        let record_buf = create_mapped_record("ACGT", &[30, 30, 30, 30], "4I", 1);
        let mut raw = encode_record_buf_to_raw(&header, &record_buf)?;
        let regenerated = regenerate_alignment_tags_raw(&mut raw, &header, &reference)?;

        assert!(regenerated, "tags were written for zero-ref-span raw record; return must be true");
        let aux_off = fgumi_raw_bam::aux_data_offset_from_record(&raw).unwrap_or(raw.len());
        let aux = &raw[aux_off..];
        assert_eq!(fgumi_raw_bam::find_int_tag(aux, SamTag::NM), Some(0i64));
        assert_eq!(fgumi_raw_bam::find_int_tag(aux, SamTag::UQ), Some(0i64));
        assert_eq!(
            fgumi_raw_bam::find_string_tag(aux, SamTag::MD)
                .map(|s| std::str::from_utf8(s).expect("MD tag should be valid UTF-8")),
            Some("0")
        );

        Ok(())
    }
}
