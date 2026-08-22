//! Read clipping utilities for BAM/SAM records.
//!
//! This module provides functionality for clipping reads in various ways (soft, hard, etc.)
//! and is used by tools like `clip` and consensus calling tools.

use noodles::sam::alignment::record::cigar::op::Kind;

/// Base-oriented aux tags that are reverse-complemented (with `SEQ`) when a read is
/// reverse-complemented, mirroring htsjdk `SAMRecord.TAGS_TO_REVERSE_COMPLEMENT`.
const TAGS_TO_REVERSE_COMPLEMENT: [fgumi_raw_bam::SamTag; 2] =
    [fgumi_raw_bam::SamTag::E2, fgumi_raw_bam::SamTag::SQ];

/// Quality-oriented aux tags that are reversed (with `QUAL`) when a read is
/// reverse-complemented, mirroring htsjdk `SAMRecord.TAGS_TO_REVERSE`.
const TAGS_TO_REVERSE: [fgumi_raw_bam::SamTag; 2] =
    [fgumi_raw_bam::SamTag::OQ, fgumi_raw_bam::SamTag::U2];

/// Modes of clipping that can be applied to reads
#[derive(Debug, Clone, Copy, PartialEq, Eq, clap::ValueEnum)]
pub enum ClippingMode {
    /// Soft clip: convert bases to S operators in CIGAR, keep bases and qualities
    Soft,
    /// Soft clip with masking: convert to S operators and mask bases to N, qualities to min
    #[value(name = "soft-with-mask")]
    SoftWithMask,
    /// Hard clip: remove bases, qualities, and convert to H operators in CIGAR
    Hard,
}

impl std::fmt::Display for ClippingMode {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        match self {
            Self::Soft => write!(f, "soft"),
            Self::SoftWithMask => write!(f, "soft-with-mask"),
            Self::Hard => write!(f, "hard"),
        }
    }
}

/// A raw-byte clipper that operates directly on [`fgumi_raw_bam::RawRecord`].
///
/// This is the clipper every production caller uses. It operates on `RawRecord` bytes
/// instead of noodles `RecordBuf` objects, following the same fgbio-derived semantics (same
/// methods, same `ClippingMode`) as the former `RecordBuf`-based `SamRecordClipper` it replaced.
///
/// All positions in the raw-byte world are 0-based (BAM spec). The algorithms convert
/// to/from 1-based where required for CIGAR reference-position arithmetic.
pub struct RawRecordClipper {
    mode: ClippingMode,
    auto_clip_attributes: bool,
}

impl RawRecordClipper {
    /// Creates a new raw-byte clipper with the specified mode.
    #[must_use]
    pub fn new(mode: ClippingMode) -> Self {
        Self { mode, auto_clip_attributes: false }
    }

    /// Creates a new raw-byte clipper with auto-clip attributes enabled.
    ///
    /// When enabled with hard clipping mode, any string or array tags that are the same
    /// length as the read's sequence will be automatically clipped to match.
    #[must_use]
    pub fn with_auto_clip(mode: ClippingMode, auto_clip_attributes: bool) -> Self {
        Self { mode, auto_clip_attributes }
    }

    /// Clip per-base tags whose length equals `old_length` when in hard-clip mode.
    ///
    /// Collects tags that must change, then re-applies them through `RawTagsEditor`.
    fn clip_extended_attributes_raw(
        &self,
        record: &mut fgumi_raw_bam::RawRecord,
        remove: usize,
        from_start: bool,
    ) {
        use fgumi_raw_bam::TagValue;

        if !matches!(self.mode, ClippingMode::Hard) || remove == 0 || !self.auto_clip_attributes {
            return;
        }

        let new_length = record.l_seq() as usize;
        let old_length = new_length + remove;

        // Collect tags to update: must release the immutable borrow before mutating.
        let aux = fgumi_raw_bam::aux_data_slice(record.as_ref()).to_vec();
        let view = fgumi_raw_bam::RawTagsView::new(&aux);

        // Tag update instructions: (tag, type_byte, new_bytes)
        let mut string_updates: Vec<([u8; 2], Vec<u8>)> = Vec::new();
        let mut array_updates: Vec<([u8; 2], u8, Vec<u8>)> = Vec::new();

        for entry in view.iter_typed() {
            let (tag, value) = entry;
            match value {
                TagValue::String(s) => {
                    if s.len() == old_length {
                        let (start, end) =
                            if from_start { (remove, old_length) } else { (0, new_length) };
                        string_updates.push((tag, s[start..end].to_vec()));
                    }
                }
                TagValue::Array(arr) => {
                    if arr.count == old_length {
                        let (start, end) =
                            if from_start { (remove, old_length) } else { (0, new_length) };
                        let start_byte = start * arr.elem_size;
                        let end_byte = end * arr.elem_size;
                        array_updates.push((
                            tag,
                            arr.elem_type,
                            arr.data[start_byte..end_byte].to_vec(),
                        ));
                    }
                }
                _ => {}
            }
        }

        // Apply updates via RawTagsEditor
        let mut editor = record.tags_editor();
        for (tag, value) in &string_updates {
            editor.update_string(tag, value);
        }
        for (tag, elem_type, data) in &array_updates {
            raw_clip_update_array_tag(&mut editor, *tag, *elem_type, data);
        }
    }

    /// Number of bases available to be clipped from the alignment.
    ///
    /// The sum of the alignment (`M`/`=`/`X`) and insertion (`I`) CIGAR op lengths.
    /// Raw-byte equivalent of the former `SamRecordClipper::number_of_clippable_bases`.
    fn number_of_clippable_bases_raw(ops: &[u32]) -> usize {
        ops.iter()
            .filter(|&&op| matches!(op & 0xF, 0 | 1 | 7 | 8)) // M, I, =, X
            .map(|&op| (op >> 4) as usize)
            .sum()
    }

    /// Unmaps a read the way htsjdk `SAMUtils.makeReadUnmapped` does.
    ///
    /// Raw-byte equivalent of the former `SamRecordClipper::make_read_unmapped`: clears the
    /// reference, position, mapping quality, CIGAR, and template length, clears the
    /// duplicate/secondary/supplementary/proper-pair flags, and sets the unmapped flag.
    /// Reverse-strand reads have their bases reverse-complemented, qualities reversed, and
    /// strand-sensitive aux tags reoriented (the base-oriented [`TAGS_TO_REVERSE_COMPLEMENT`]
    /// reverse-complemented, the quality-oriented [`TAGS_TO_REVERSE`] reversed) before the
    /// reverse flag is cleared. Bases and qualities are otherwise preserved.
    fn make_read_unmapped_raw(record: &mut fgumi_raw_bam::RawRecord) {
        use fgumi_raw_bam::flags as rflags;

        let flags = record.flags();
        if flags & rflags::REVERSE != 0 {
            let rc = fgumi_dna::reverse_complement(&record.sequence_vec());
            let mut quals = record.quality_scores().to_vec();
            quals.reverse();
            record.set_sequence_and_qualities(&rc, &quals);
            // Reorient the strand-sensitive aux tags alongside SEQ/QUAL, matching htsjdk
            // `SAMRecord.reverseComplement`. Only the `Z`-string encoding is handled -- the sole
            // form these tags take in practice; htsjdk also reverses array-encoded (`B:...`)
            // values, but none of these tags is emitted as an array by real tooling.
            let mut tags = record.tags_mut();
            for tag in TAGS_TO_REVERSE_COMPLEMENT {
                tags.reverse_complement_string(tag);
            }
            for tag in TAGS_TO_REVERSE {
                tags.reverse_string(tag);
            }
        }

        let new_flags = (flags
            & !(rflags::REVERSE
                | rflags::DUPLICATE
                | rflags::SECONDARY
                | rflags::SUPPLEMENTARY
                | rflags::PROPER_PAIR))
            | rflags::UNMAPPED;
        record.set_flags(new_flags);
        record.set_ref_id(-1);
        record.set_pos(-1);
        record.set_mapq(0);
        record.set_template_length(0);
        record.set_cigar_ops(&[]);
        // Now unmapped (POS = -1): the bin must become the SAM unmapped bin (4680).
        record.recompute_bin();
    }

    /// Clips a specified number of bases from the start (left side) of the alignment.
    ///
    /// Returns the number of bases actually clipped.
    ///
    /// # Panics
    ///
    /// Panics if the reference-position delta from clipping does not fit in `i32`
    /// or would advance `POS` past `i32::MAX`. CIGAR lengths are bounded by the
    /// BAM read length (well below `i32::MAX`), so this should not happen for
    /// well-formed records.
    #[expect(
        clippy::too_many_lines,
        reason = "mirrors SamRecordClipper::clip_start_of_alignment with raw-byte API"
    )]
    #[expect(
        clippy::cast_possible_truncation,
        reason = "CIGAR lengths are bounded by BAM read length which fits in u32"
    )]
    pub fn clip_start_of_alignment(
        &self,
        record: &mut fgumi_raw_bam::RawRecord,
        bases_to_clip: usize,
    ) -> usize {
        if bases_to_clip == 0 {
            return 0;
        }

        // Don't clip unmapped reads
        if record.flags() & fgumi_raw_bam::flags::UNMAPPED != 0 {
            return 0;
        }

        if record.l_seq() == 0 {
            return 0;
        }

        let old_ops: Vec<u32> = record.cigar_ops_vec();

        // If the read has no more clippable (aligned/inserted) bases than requested, it
        // cannot retain an alignment; unmap it (fgbio SamRecordClipper.scala:111-114).
        let num_clippable = Self::number_of_clippable_bases_raw(&old_ops);
        if num_clippable <= bases_to_clip {
            Self::make_read_unmapped_raw(record);
            return num_clippable;
        }

        // Extract existing hard and soft clips from the start
        let existing_hard_clip: usize = old_ops
            .iter()
            .take_while(|&&op| (op & 0xF) == 5) // HardClip
            .map(|&op| (op >> 4) as usize)
            .sum();

        let existing_soft_clip: usize = old_ops
            .iter()
            .skip_while(|&&op| (op & 0xF) == 5) // skip HardClip
            .take_while(|&&op| (op & 0xF) == 4) // SoftClip
            .map(|&op| (op >> 4) as usize)
            .sum();

        // Skip to operations after existing clips
        let post_clip_ops: Vec<u32> =
            old_ops.iter().copied().skip_while(|&op| matches!(op & 0xF, 4 | 5)).collect();

        let mut read_bases_clipped: usize = 0;
        let mut ref_bases_clipped: usize = 0;
        let mut new_ops: Vec<u32> = Vec::new();
        let mut iter = post_clip_ops.iter().peekable();

        while read_bases_clipped < bases_to_clip
            || (read_bases_clipped == bases_to_clip
                && new_ops.is_empty()
                && iter.peek().map(|&op| op & 0xF) == Some(2))
        // Deletion
        {
            let Some(&op) = iter.next() else { break };
            let kind = op & 0xF;
            let len = (op >> 4) as usize;

            // Match SamRecordClipper parity: count M/I/=/X as read-consuming and M/D/=/X
            // as ref-consuming. N is deliberately excluded from ref here so the raw and
            // typed clippers agree on bases_to_clip and alignment_start for spliced reads.
            let consumes_read = matches!(kind, 0 | 1 | 7 | 8); // M, I, =, X
            let consumes_ref = matches!(kind, 0 | 2 | 7 | 8); // M, D, =, X

            if consumes_read && len > (bases_to_clip - read_bases_clipped) {
                if kind == 1 {
                    // Insertion: consume entire op at clip boundary
                    read_bases_clipped += len;
                } else {
                    let remaining_clip = bases_to_clip - read_bases_clipped;
                    let remaining_length = len - remaining_clip;
                    read_bases_clipped += remaining_clip;
                    ref_bases_clipped += remaining_clip;
                    new_ops.push((remaining_length as u32) << 4 | kind);
                }
            } else {
                if consumes_read {
                    read_bases_clipped += len;
                }
                if consumes_ref {
                    ref_bases_clipped += len;
                }
            }
        }

        // Add remaining operations
        new_ops.extend(iter.copied());

        // Prepend appropriate clipping operators
        let (final_ops, bases_to_remove) = match self.mode {
            ClippingMode::Hard => {
                let added_hard_clip = existing_soft_clip + read_bases_clipped;
                let total_hard_clip = existing_hard_clip + added_hard_clip;
                let mut result = Vec::with_capacity(new_ops.len() + 1);
                result.push((total_hard_clip as u32) << 4 | 5); // HardClip
                result.extend_from_slice(&new_ops);
                (result, added_hard_clip)
            }
            ClippingMode::Soft | ClippingMode::SoftWithMask => {
                let total_soft_clip = existing_soft_clip + read_bases_clipped;
                let mut result = Vec::new();
                if existing_hard_clip > 0 {
                    result.push((existing_hard_clip as u32) << 4 | 5); // HardClip
                }
                result.push((total_soft_clip as u32) << 4 | 4); // SoftClip
                result.extend_from_slice(&new_ops);
                (result, 0)
            }
        };

        // Update CIGAR
        record.set_cigar_ops(&final_ops);

        // Update alignment start position (0-based in BAM)
        if ref_bases_clipped > 0 {
            let old_pos = record.pos();
            if old_pos >= 0 {
                let delta =
                    i32::try_from(ref_bases_clipped).expect("reference clip delta must fit in i32");
                let new_pos =
                    old_pos.checked_add(delta).expect("clipping would move BAM POS past i32::MAX");
                record.set_pos(new_pos);
            }
        }

        // Handle sequence and quality updates
        match self.mode {
            ClippingMode::Soft => {
                // Keep sequence and qualities as-is
            }
            ClippingMode::SoftWithMask => {
                let seq_len = record.l_seq() as usize;
                let mut new_seq = record.sequence_vec();
                let mut new_qual = record.quality_scores().to_vec();
                let total_soft_clip = existing_soft_clip + read_bases_clipped;
                for i in 0..total_soft_clip.min(seq_len) {
                    new_seq[i] = fgumi_dna::NO_CALL_BASE;
                    new_qual[i] = fgumi_dna::MIN_PHRED;
                }
                record.set_sequence_and_qualities(&new_seq, &new_qual);
            }
            ClippingMode::Hard => {
                if bases_to_remove > 0 {
                    let seq = record.sequence_vec();
                    let qual = record.quality_scores().to_vec();
                    let new_seq = seq[bases_to_remove..].to_vec();
                    let new_qual = qual[bases_to_remove..].to_vec();
                    record.set_sequence_and_qualities(&new_seq, &new_qual);
                    self.clip_extended_attributes_raw(record, bases_to_remove, true);
                }
            }
        }

        // POS and/or the CIGAR changed above; refresh the BAM bin so the verbatim
        // raw write emits a correct index bin (the encoder step that would normally
        // recompute it is bypassed by the raw pipeline).
        record.recompute_bin();
        read_bases_clipped
    }

    /// Clips a specified number of bases from the end (right side) of the alignment.
    ///
    /// Returns the number of bases actually clipped.
    #[expect(
        clippy::cast_possible_truncation,
        reason = "CIGAR lengths are bounded by BAM read length which fits in u32"
    )]
    #[expect(
        clippy::too_many_lines,
        reason = "mirrors clip_start_of_alignment with symmetric end-clipping logic"
    )]
    pub fn clip_end_of_alignment(
        &self,
        record: &mut fgumi_raw_bam::RawRecord,
        bases_to_clip: usize,
    ) -> usize {
        if bases_to_clip == 0 {
            return 0;
        }

        if record.flags() & fgumi_raw_bam::flags::UNMAPPED != 0 {
            return 0;
        }

        if record.l_seq() == 0 {
            return 0;
        }

        let old_ops: Vec<u32> = record.cigar_ops_vec();

        // If the read has no more clippable (aligned/inserted) bases than requested, it
        // cannot retain an alignment; unmap it (fgbio SamRecordClipper.scala:157-160).
        let num_clippable = Self::number_of_clippable_bases_raw(&old_ops);
        if num_clippable <= bases_to_clip {
            Self::make_read_unmapped_raw(record);
            return num_clippable;
        }

        // Extract existing hard and soft clips from the end (reverse order)
        let existing_hard_clip: usize = old_ops
            .iter()
            .rev()
            .take_while(|&&op| (op & 0xF) == 5)
            .map(|&op| (op >> 4) as usize)
            .sum();

        let existing_soft_clip: usize = old_ops
            .iter()
            .rev()
            .skip_while(|&&op| (op & 0xF) == 5)
            .take_while(|&&op| (op & 0xF) == 4)
            .map(|&op| (op >> 4) as usize)
            .sum();

        // Strip trailing clips
        let mut post_clip_ops: Vec<u32> =
            old_ops.iter().copied().rev().skip_while(|&op| matches!(op & 0xF, 4 | 5)).collect();
        post_clip_ops.reverse();

        let mut read_bases_clipped: usize = 0;
        let mut new_ops: Vec<u32> = Vec::new();
        let mut iter = post_clip_ops.iter().rev().peekable();

        while read_bases_clipped < bases_to_clip
            || (read_bases_clipped == bases_to_clip
                && new_ops.is_empty()
                && iter.peek().map(|&op| op & 0xF) == Some(2))
        // Deletion
        {
            let Some(&op) = iter.next() else { break };
            let kind = op & 0xF;
            let len = (op >> 4) as usize;

            let consumes_read = matches!(kind, 0 | 1 | 7 | 8);

            if consumes_read && len > (bases_to_clip - read_bases_clipped) {
                if kind == 1 {
                    // Insertion: consume entire op
                    read_bases_clipped += len;
                } else {
                    let remaining_clip = bases_to_clip - read_bases_clipped;
                    let remaining_length = len - remaining_clip;
                    read_bases_clipped += remaining_clip;
                    new_ops.push((remaining_length as u32) << 4 | kind);
                }
            } else if consumes_read {
                read_bases_clipped += len;
            }
        }

        // Collect remaining (iter is reversed, so just collect and re-reverse)
        let remaining: Vec<u32> = iter.copied().collect();
        new_ops.extend(remaining.iter());
        new_ops.reverse();

        // Append appropriate clipping operators
        let (final_ops, bases_to_remove) = match self.mode {
            ClippingMode::Hard => {
                let added_hard_clip = existing_soft_clip + read_bases_clipped;
                let total_hard_clip = existing_hard_clip + added_hard_clip;
                let mut result = new_ops;
                result.push((total_hard_clip as u32) << 4 | 5); // HardClip
                (result, added_hard_clip)
            }
            ClippingMode::Soft | ClippingMode::SoftWithMask => {
                let total_soft_clip = existing_soft_clip + read_bases_clipped;
                let mut result = new_ops;
                result.push((total_soft_clip as u32) << 4 | 4); // SoftClip
                if existing_hard_clip > 0 {
                    result.push((existing_hard_clip as u32) << 4 | 5); // HardClip
                }
                (result, 0)
            }
        };

        // Update CIGAR
        record.set_cigar_ops(&final_ops);

        // Handle sequence and quality updates
        let seq_len = record.l_seq() as usize;
        match self.mode {
            ClippingMode::Soft => {}
            ClippingMode::SoftWithMask => {
                let mut new_seq = record.sequence_vec();
                let mut new_qual = record.quality_scores().to_vec();
                let total_soft_clip = existing_soft_clip + read_bases_clipped;
                let start_mask = seq_len.saturating_sub(total_soft_clip);
                for i in start_mask..seq_len {
                    new_seq[i] = fgumi_dna::NO_CALL_BASE;
                    new_qual[i] = fgumi_dna::MIN_PHRED;
                }
                record.set_sequence_and_qualities(&new_seq, &new_qual);
            }
            ClippingMode::Hard => {
                if bases_to_remove > 0 {
                    let seq = record.sequence_vec();
                    let qual = record.quality_scores().to_vec();
                    let keep_len = seq_len.saturating_sub(bases_to_remove);
                    let new_seq = seq[..keep_len].to_vec();
                    let new_qual = qual[..keep_len].to_vec();
                    record.set_sequence_and_qualities(&new_seq, &new_qual);
                    self.clip_extended_attributes_raw(record, bases_to_remove, false);
                }
            }
        }

        // End-clipping leaves POS unchanged but shrinks the aligned span, which can
        // move the read into a finer bin; refresh it for the verbatim raw write.
        record.recompute_bin();
        read_bases_clipped
    }

    /// Clips bases from the 5' end of the read (strand-aware).
    ///
    /// Returns the number of bases actually clipped.
    pub fn clip_5_prime_end_of_alignment(
        &self,
        record: &mut fgumi_raw_bam::RawRecord,
        bases_to_clip: usize,
    ) -> usize {
        if record.flags() & fgumi_raw_bam::flags::REVERSE != 0 {
            self.clip_end_of_alignment(record, bases_to_clip)
        } else {
            self.clip_start_of_alignment(record, bases_to_clip)
        }
    }

    /// Clips bases from the 3' end of the read (strand-aware).
    ///
    /// Returns the number of bases actually clipped.
    pub fn clip_3_prime_end_of_alignment(
        &self,
        record: &mut fgumi_raw_bam::RawRecord,
        bases_to_clip: usize,
    ) -> usize {
        if record.flags() & fgumi_raw_bam::flags::REVERSE != 0 {
            self.clip_start_of_alignment(record, bases_to_clip)
        } else {
            self.clip_end_of_alignment(record, bases_to_clip)
        }
    }

    /// Helper: count query bases in a raw CIGAR corresponding to a reference region.
    fn calculate_query_bases_for_ref_region_raw(
        ops: &[u32],
        ref_bases: usize,
        from_start: bool,
    ) -> usize {
        let mut remaining_ref = ref_bases;
        let mut query_bases: usize = 0;

        let iter: Box<dyn Iterator<Item = &u32>> =
            if from_start { Box::new(ops.iter()) } else { Box::new(ops.iter().rev()) };

        for &op in iter {
            if remaining_ref == 0 {
                break;
            }
            let kind = op & 0xF;
            let len = (op >> 4) as usize;

            // Match SamRecordClipper parity (see calculate_query_bases_for_ref_region):
            // M/D/=/X consume ref; N is excluded so the raw and typed paths agree.
            let consumes_ref = matches!(kind, 0 | 2 | 7 | 8);
            let consumes_query = matches!(kind, 0 | 1 | 7 | 8); // M, I, =, X (not S)

            if consumes_ref {
                let ref_consumed = len.min(remaining_ref);
                remaining_ref -= ref_consumed;
                if consumes_query {
                    query_bases += ref_consumed;
                }
            } else if consumes_query && remaining_ref > 0 {
                // Insertion before we've consumed all the ref we need
                query_bases += len;
            }
        }

        query_bases
    }

    /// Clips overlapping portions of an FR read pair.
    ///
    /// Returns `(bases_clipped_r1, bases_clipped_r2)`.
    pub fn clip_overlapping_reads(
        &self,
        r1: &mut fgumi_raw_bam::RawRecord,
        r2: &mut fgumi_raw_bam::RawRecord,
    ) -> (usize, usize) {
        // Gate on the symmetric, order-independent per-*pair* FR classifier. The per-record
        // `is_fr_pair_raw` derives the mate 5' from TLEN on its forward-strand arm, which
        // misclassifies dovetail FR pairs (htsjdk/samtools#1771) and would skip overlap clipping
        // on a valid pair. Both records are in hand here, so `is_primary_fr_pair_raw` (the same
        // gate the sibling `clip_extending_past_mate_ends` uses) is the correct check.
        if !fgumi_raw_bam::is_primary_fr_pair_raw(r1.as_ref(), r2.as_ref()) {
            return (0, 0);
        }

        // Normalize by strand so the positive-strand read is treated as r1 and the
        // negative-strand read as r2 (see the typed `clip_overlapping_reads` for the
        // rationale; mirrors fgbio `SamRecordClipper.clipOverlappingReads:305`).
        let swapped = r1.flags() & fgumi_raw_bam::flags::REVERSE != 0;
        let (r1, r2) = if swapped { (r2, r1) } else { (r1, r2) };

        let Some(r1_start) = r1.alignment_start_1based() else {
            return (0, 0);
        };
        let Some(r2_start) = r2.alignment_start_1based() else {
            return (0, 0);
        };

        let r1_ops = r1.cigar_ops_vec();
        let r2_ops = r2.cigar_ops_vec();
        let r1_ref_len = crate::record_utils::cigar_reference_length_raw(&r1_ops);
        let r2_ref_len = crate::record_utils::cigar_reference_length_raw(&r2_ops);
        let r1_end = r1_start + r1_ref_len.saturating_sub(1);
        let r2_end = r2_start + r2_ref_len.saturating_sub(1);

        let overlap_start = r1_start.max(r2_start);
        let overlap_end = r1_end.min(r2_end);

        if overlap_start > overlap_end {
            return (0, 0);
        }

        let mut midpoint = usize::midpoint(r1_start, r2_end);

        if midpoint > r1_end {
            midpoint = r1_end;
        } else if midpoint < r2_start {
            midpoint = r2_start.saturating_sub(1);
        }

        let r1_bases_to_clip = if r1_end > midpoint {
            let ref_bases_to_clip = r1_end - midpoint;
            Self::calculate_query_bases_for_ref_region_raw(&r1_ops, ref_bases_to_clip, false)
        } else {
            0
        };

        let r2_bases_to_clip = if midpoint + 1 > r2_start {
            let ref_bases_to_clip = midpoint + 1 - r2_start;
            Self::calculate_query_bases_for_ref_region_raw(&r2_ops, ref_bases_to_clip, true)
        } else {
            0
        };

        let clipped_r1 =
            if r1_bases_to_clip > 0 { self.clip_end_of_alignment(r1, r1_bases_to_clip) } else { 0 };
        let clipped_r2 = if r2_bases_to_clip > 0 {
            self.clip_start_of_alignment(r2, r2_bases_to_clip)
        } else {
            0
        };

        if matches!(self.mode, ClippingMode::Hard) {
            let _ = self.upgrade_all_clipping_raw(r1);
            let _ = self.upgrade_all_clipping_raw(r2);
        }

        // Map clip counts back to the caller's original (r1, r2) argument order.
        if swapped { (clipped_r2, clipped_r1) } else { (clipped_r1, clipped_r2) }
    }

    /// Clips reads that extend beyond their mate's alignment ends.
    ///
    /// Returns `(bases_clipped_r1, bases_clipped_r2)`.
    ///
    /// Delegates the past-mate distance to
    /// [`fgumi_raw_bam::num_bases_extending_past_mate_vs_mate_raw`], the query-space, mate-in-hand
    /// count fgbio#1172 lands (fixing issue #760). This raw-byte path backs the LIVE `fgumi clip`
    /// command; the reference-space `RecordBuf` past-mate methods on the former
    /// `SamRecordClipper`, which `fgumi clip` never constructed, were removed as redundant. That function
    /// gates on the symmetric [`fgumi_raw_bam::is_primary_fr_pair_raw`] internally, so no separate
    /// FR-pair check is needed here.
    pub fn clip_extending_past_mate_ends(
        &self,
        r1: &mut fgumi_raw_bam::RawRecord,
        r2: &mut fgumi_raw_bam::RawRecord,
    ) -> (usize, usize) {
        // NB: compute BOTH counts before clipping either read (fgbio#1172): otherwise the second
        // read is measured against an alignment the first clip already shortened. The FR-pair
        // gate is inside num_bases_extending_past_mate_vs_mate_raw (returns 0 for non-FR).
        let n1 = fgumi_raw_bam::num_bases_extending_past_mate_vs_mate_raw(r1.as_ref(), r2.as_ref());
        let n2 = fgumi_raw_bam::num_bases_extending_past_mate_vs_mate_raw(r2.as_ref(), r1.as_ref());
        let hard1 = Self::existing_hard_clip_3_prime_raw(r1);
        let hard2 = Self::existing_hard_clip_3_prime_raw(r2);
        let c1 = if n1 > 0 { self.clip_3_prime_end_of_read_raw(r1, n1 + hard1) } else { 0 };
        let c2 = if n2 > 0 { self.clip_3_prime_end_of_read_raw(r2, n2 + hard2) } else { 0 };
        (c1, c2)
    }

    /// Existing hard-clipping at the record's 3' end (trailing for +strand, leading for
    /// -strand). `clip_3_prime_end_of_read_raw` takes TOTAL desired clipping incl. existing; the
    /// query-space count excludes hard clips, so they must be added back (fgbio#1172). Zero when
    /// none present, so applying it unconditionally is a no-op in that case.
    fn existing_hard_clip_3_prime_raw(rec: &fgumi_raw_bam::RawRecord) -> usize {
        let ops = rec.cigar_ops_vec();
        if rec.flags() & fgumi_raw_bam::flags::REVERSE != 0 {
            ops.iter().take_while(|&&op| (op & 0xF) == 5).map(|&op| (op >> 4) as usize).sum()
        } else {
            ops.iter().rev().take_while(|&&op| (op & 0xF) == 5).map(|&op| (op >> 4) as usize).sum()
        }
    }

    /// Ensures at least `clip_length` bases are clipped at the start of the read.
    pub fn clip_start_of_read_raw(
        &self,
        record: &mut fgumi_raw_bam::RawRecord,
        clip_length: usize,
    ) -> usize {
        let ops = record.cigar_ops_vec();
        let existing_clipping: usize = ops
            .iter()
            .take_while(|&&op| matches!(op & 0xF, 4 | 5))
            .map(|&op| (op >> 4) as usize)
            .sum();

        if clip_length > existing_clipping {
            self.clip_start_of_alignment(record, clip_length - existing_clipping)
        } else {
            self.upgrade_clipping_raw(record, clip_length, true);
            0
        }
    }

    /// Ensures at least `clip_length` bases are clipped at the end of the read.
    pub fn clip_end_of_read_raw(
        &self,
        record: &mut fgumi_raw_bam::RawRecord,
        clip_length: usize,
    ) -> usize {
        let ops = record.cigar_ops_vec();
        let existing_clipping: usize = ops
            .iter()
            .rev()
            .take_while(|&&op| matches!(op & 0xF, 4 | 5))
            .map(|&op| (op >> 4) as usize)
            .sum();

        if clip_length > existing_clipping {
            self.clip_end_of_alignment(record, clip_length - existing_clipping)
        } else {
            self.upgrade_clipping_raw(record, clip_length, false);
            0
        }
    }

    /// Ensures at least `clip_length` bases are clipped at the 5' end (strand-aware).
    pub fn clip_5_prime_end_of_read_raw(
        &self,
        record: &mut fgumi_raw_bam::RawRecord,
        clip_length: usize,
    ) -> usize {
        if record.flags() & fgumi_raw_bam::flags::REVERSE != 0 {
            self.clip_end_of_read_raw(record, clip_length)
        } else {
            self.clip_start_of_read_raw(record, clip_length)
        }
    }

    /// Ensures at least `clip_length` bases are clipped at the 3' end (strand-aware).
    pub fn clip_3_prime_end_of_read_raw(
        &self,
        record: &mut fgumi_raw_bam::RawRecord,
        clip_length: usize,
    ) -> usize {
        if record.flags() & fgumi_raw_bam::flags::REVERSE != 0 {
            self.clip_start_of_read_raw(record, clip_length)
        } else {
            self.clip_end_of_read_raw(record, clip_length)
        }
    }

    /// Upgrades soft clipping to hard clipping in a raw record.
    #[expect(
        clippy::cast_possible_truncation,
        reason = "CIGAR lengths are bounded by BAM read length which fits in u32"
    )]
    fn upgrade_clipping_raw(
        &self,
        record: &mut fgumi_raw_bam::RawRecord,
        length: usize,
        from_start: bool,
    ) {
        if self.mode == ClippingMode::Soft || length == 0 {
            return;
        }

        let old_ops: Vec<u32> = record.cigar_ops_vec();

        let (hard_clipped, soft_clipped) = if from_start {
            let hard: usize = old_ops
                .iter()
                .take_while(|&&op| (op & 0xF) == 5)
                .map(|&op| (op >> 4) as usize)
                .sum();
            let soft: usize = old_ops
                .iter()
                .skip_while(|&&op| (op & 0xF) == 5)
                .take_while(|&&op| (op & 0xF) == 4)
                .map(|&op| (op >> 4) as usize)
                .sum();
            (hard, soft)
        } else {
            let hard: usize = old_ops
                .iter()
                .rev()
                .take_while(|&&op| (op & 0xF) == 5)
                .map(|&op| (op >> 4) as usize)
                .sum();
            let soft: usize = old_ops
                .iter()
                .rev()
                .skip_while(|&&op| (op & 0xF) == 5)
                .take_while(|&&op| (op & 0xF) == 4)
                .map(|&op| (op >> 4) as usize)
                .sum();
            (hard, soft)
        };

        if hard_clipped >= length || soft_clipped == 0 {
            return;
        }

        let length_to_upgrade = soft_clipped.min(length - hard_clipped);

        if self.mode == ClippingMode::Hard {
            let ops_to_process: Vec<u32> =
                if from_start { old_ops.clone() } else { old_ops.iter().copied().rev().collect() };

            let mut new_ops: Vec<u32> = Vec::new();
            let mut i = 0;
            let mut existing_hard: usize = 0;
            let mut existing_soft: usize = 0;

            // Count leading hard clips
            while i < ops_to_process.len() && (ops_to_process[i] & 0xF) == 5 {
                existing_hard += (ops_to_process[i] >> 4) as usize;
                i += 1;
            }
            // Count leading soft clips
            while i < ops_to_process.len() && (ops_to_process[i] & 0xF) == 4 {
                existing_soft += (ops_to_process[i] >> 4) as usize;
                i += 1;
            }

            let new_hard_count = existing_hard + length_to_upgrade;
            new_ops.push((new_hard_count as u32) << 4 | 5);
            if existing_soft > length_to_upgrade {
                new_ops.push(((existing_soft - length_to_upgrade) as u32) << 4 | 4);
            }
            new_ops.extend_from_slice(&ops_to_process[i..]);

            let final_ops: Vec<u32> =
                if from_start { new_ops } else { new_ops.iter().copied().rev().collect() };
            record.set_cigar_ops(&final_ops);

            // Update sequence and quals
            let seq = record.sequence_vec();
            let qual = record.quality_scores().to_vec();
            if from_start {
                let new_seq = seq[length_to_upgrade..].to_vec();
                let new_qual = qual[length_to_upgrade..].to_vec();
                record.set_sequence_and_qualities(&new_seq, &new_qual);
            } else {
                let new_len = seq.len().saturating_sub(length_to_upgrade);
                let new_seq = seq[..new_len].to_vec();
                let new_qual = qual[..new_len].to_vec();
                record.set_sequence_and_qualities(&new_seq, &new_qual);
            }
            self.clip_extended_attributes_raw(record, length_to_upgrade, from_start);
        } else if self.mode == ClippingMode::SoftWithMask {
            // SoftWithMask: mask the upgraded soft-clipped bases to N/min-quality without
            // altering the CIGAR (see the typed `upgrade_clipping` for rationale).
            let mut new_seq = record.sequence_vec();
            let mut new_qual = record.quality_scores().to_vec();
            let seq_len = new_seq.len();
            let (mask_start, mask_end) = if from_start {
                (0, length_to_upgrade)
            } else {
                (seq_len - length_to_upgrade, seq_len)
            };
            for i in mask_start..mask_end {
                new_seq[i] = fgumi_dna::NO_CALL_BASE;
                new_qual[i] = fgumi_dna::MIN_PHRED;
            }
            record.set_sequence_and_qualities(&new_seq, &new_qual);
        }
    }

    /// Upgrades all existing clipping in a raw record to the current clipping mode.
    ///
    /// In `Hard` mode soft clips are converted to hard clips (bases removed); in
    /// `SoftWithMask` mode the CIGAR is left intact and the soft-clipped bases are
    /// masked to N with minimum quality. `Soft` mode is a no-op.
    ///
    /// Returns `(leading_soft, trailing_soft)` — the soft-clip counts upgraded at
    /// each end (converted in `Hard` mode, masked in `SoftWithMask` mode).
    ///
    /// # Errors
    ///
    /// This function does not currently return errors but uses `anyhow::Result` for
    /// API symmetry with `upgrade_all_clipping`.
    #[expect(
        clippy::too_many_lines,
        reason = "mirrors SamRecordClipper::upgrade_all_clipping with raw-byte CIGAR surgery"
    )]
    #[expect(
        clippy::cast_possible_truncation,
        reason = "CIGAR lengths are bounded by BAM read length which fits in u32"
    )]
    pub fn upgrade_all_clipping_raw(
        &self,
        record: &mut fgumi_raw_bam::RawRecord,
    ) -> anyhow::Result<(usize, usize)> {
        if matches!(self.mode, ClippingMode::Soft) {
            return Ok((0, 0));
        }
        if record.flags() & fgumi_raw_bam::flags::UNMAPPED != 0 {
            return Ok((0, 0));
        }

        let ops: Vec<u32> = record.cigar_ops_vec();
        let has_soft_clips = ops.iter().any(|&op| (op & 0xF) == 4);
        if !has_soft_clips {
            return Ok((0, 0));
        }

        let mut leading_hard: usize = 0;
        let mut leading_soft: usize = 0;
        let mut trailing_soft: usize = 0;
        let mut trailing_hard: usize = 0;

        for &op in &ops {
            match op & 0xF {
                5 => leading_hard += (op >> 4) as usize, // HardClip
                4 => {
                    // SoftClip
                    leading_soft += (op >> 4) as usize;
                    break;
                }
                _ => break,
            }
        }
        for &op in ops.iter().rev() {
            match op & 0xF {
                5 => trailing_hard += (op >> 4) as usize,
                4 => {
                    trailing_soft += (op >> 4) as usize;
                    break;
                }
                _ => break,
            }
        }

        // SoftWithMask: mask existing soft-clipped bases at both ends via upgrade_clipping_raw,
        // leaving the CIGAR intact (see the typed upgrade_all_clipping for rationale).
        if self.mode == ClippingMode::SoftWithMask {
            if leading_soft > 0 {
                self.upgrade_clipping_raw(record, leading_hard + leading_soft, true);
            }
            if trailing_soft > 0 {
                self.upgrade_clipping_raw(record, trailing_hard + trailing_soft, false);
            }
            return Ok((leading_soft, trailing_soft));
        }

        let old_seq_len = record.l_seq() as usize;
        let mut new_cigar_ops: Vec<u32> = Vec::new();
        let mut seq_pos: usize = 0;
        let mut new_sequence: Vec<u8> = Vec::new();
        let mut new_qualities: Vec<u8> = Vec::new();
        let mut is_leading = true;
        let seq = record.sequence_vec();
        let qualities = record.quality_scores().to_vec();

        for &op in &ops {
            let kind = op & 0xF;
            let len = (op >> 4) as usize;

            match kind {
                4 => {
                    // SoftClip -> convert to HardClip
                    if is_leading && new_cigar_ops.is_empty() && leading_hard > 0 {
                        new_cigar_ops.push(((leading_hard + len) as u32) << 4 | 5);
                    } else if new_cigar_ops.last().map(|o| o & 0xF) == Some(5) {
                        let last_len = (new_cigar_ops.last().copied().unwrap_or(0) >> 4) as usize;
                        new_cigar_ops.pop();
                        new_cigar_ops.push(((last_len + len) as u32) << 4 | 5);
                    } else {
                        new_cigar_ops.push((len as u32) << 4 | 5);
                    }
                    seq_pos += len;
                }
                5 => {
                    // HardClip — merge with adjacent if needed
                    if new_cigar_ops.last().map(|o| o & 0xF) == Some(5) {
                        let last_len = (new_cigar_ops.last().copied().unwrap_or(0) >> 4) as usize;
                        new_cigar_ops.pop();
                        new_cigar_ops.push(((last_len + len) as u32) << 4 | 5);
                    } else if !is_leading || new_cigar_ops.is_empty() {
                        new_cigar_ops.push(op);
                    }
                }
                _ => {
                    is_leading = false;
                    new_cigar_ops.push(op);
                    // Copy bases/quals for query-consuming operations (M, I, =, X).
                    // Ref-only ops (D, N, P) must not advance seq_pos.
                    let consumes_query = matches!(kind, 0 | 1 | 7 | 8);
                    if consumes_query {
                        for j in 0..len {
                            if seq_pos + j < old_seq_len {
                                new_sequence.push(seq[seq_pos + j]);
                                new_qualities.push(qualities[seq_pos + j]);
                            }
                        }
                        seq_pos += len;
                    }
                }
            }
        }

        record.set_cigar_ops(&new_cigar_ops);

        if self.auto_clip_attributes && (leading_soft > 0 || trailing_soft > 0) {
            // Collect tags matching old length, then update them
            let aux = fgumi_raw_bam::aux_data_slice(record.as_ref()).to_vec();
            let view = fgumi_raw_bam::RawTagsView::new(&aux);
            let mut string_updates: Vec<([u8; 2], Vec<u8>)> = Vec::new();
            let mut array_updates: Vec<([u8; 2], u8, Vec<u8>)> = Vec::new();

            for (tag, value) in view.iter_typed() {
                use fgumi_raw_bam::TagValue;
                match value {
                    TagValue::String(s) => {
                        if s.len() == old_seq_len {
                            let start = leading_soft;
                            let end = old_seq_len - trailing_soft;
                            string_updates.push((tag, s[start..end].to_vec()));
                        }
                    }
                    TagValue::Array(arr) => {
                        if arr.count == old_seq_len {
                            let start = leading_soft * arr.elem_size;
                            let end = (old_seq_len - trailing_soft) * arr.elem_size;
                            array_updates.push((tag, arr.elem_type, arr.data[start..end].to_vec()));
                        }
                    }
                    _ => {}
                }
            }

            let mut editor = record.tags_editor();
            for (tag, value) in &string_updates {
                editor.update_string(tag, value);
            }
            for (tag, elem_type, data) in &array_updates {
                raw_clip_update_array_tag(&mut editor, *tag, *elem_type, data);
            }
        }

        record.set_sequence_and_qualities(&new_sequence, &new_qualities);

        Ok((leading_soft, trailing_soft))
    }
}

/// Dispatch a raw-byte array tag update by element type.
///
/// Routes to the appropriate `RawTagsEditor::update_array_*` variant based on `elem_type`.
/// Unrecognised element types are silently ignored.
fn raw_clip_update_array_tag(
    editor: &mut fgumi_raw_bam::RawTagsEditor<'_>,
    tag: [u8; 2],
    elem_type: u8,
    raw_bytes: &[u8],
) {
    match elem_type {
        b'c' => {
            let vals: Vec<i8> = raw_bytes.iter().map(|&b| b.cast_signed()).collect();
            editor.update_array_i8(tag, &vals);
        }
        b'C' => editor.update_array_u8(tag, raw_bytes),
        b's' => {
            let vals: Vec<i16> =
                raw_bytes.chunks_exact(2).map(|c| i16::from_le_bytes([c[0], c[1]])).collect();
            editor.update_array_i16(tag, &vals);
        }
        b'S' => {
            let vals: Vec<u16> =
                raw_bytes.chunks_exact(2).map(|c| u16::from_le_bytes([c[0], c[1]])).collect();
            editor.update_array_u16(tag, &vals);
        }
        b'i' => {
            let vals: Vec<i32> = raw_bytes
                .chunks_exact(4)
                .map(|c| i32::from_le_bytes([c[0], c[1], c[2], c[3]]))
                .collect();
            editor.update_array_i32(tag, &vals);
        }
        b'I' => {
            let vals: Vec<u32> = raw_bytes
                .chunks_exact(4)
                .map(|c| u32::from_le_bytes([c[0], c[1], c[2], c[3]]))
                .collect();
            editor.update_array_u32(tag, &vals);
        }
        b'f' => {
            let vals: Vec<f32> = raw_bytes
                .chunks_exact(4)
                .map(|c| f32::from_le_bytes([c[0], c[1], c[2], c[3]]))
                .collect();
            editor.update_array_f32(tag, &vals);
        }
        _ => {} // Unknown BAM array subtype — leave untouched
    }
}

/// Helper functions for CIGAR manipulation
pub mod cigar_utils {
    use super::Kind;
    use noodles::sam::alignment::record::Cigar as CigarTrait;

    /// Counts the number of aligned bases in a CIGAR string
    #[must_use]
    #[expect(
        clippy::redundant_closure_for_method_calls,
        reason = "Op::len is not a method on the trait"
    )]
    pub fn aligned_bases(cigar: &impl CigarTrait) -> usize {
        cigar
            .iter()
            .filter_map(Result::ok)
            .filter(|op| {
                matches!(op.kind(), Kind::Match | Kind::SequenceMatch | Kind::SequenceMismatch)
            })
            .map(|op| op.len())
            .sum()
    }

    /// Counts the number of clipped bases in a CIGAR string
    #[must_use]
    #[expect(
        clippy::redundant_closure_for_method_calls,
        reason = "Op::len is not a method on the trait"
    )]
    pub fn clipped_bases(cigar: &impl CigarTrait) -> usize {
        cigar
            .iter()
            .filter_map(Result::ok)
            .filter(|op| matches!(op.kind(), Kind::SoftClip | Kind::HardClip))
            .map(|op| op.len())
            .sum()
    }

    /// Counts reference-consuming operations.
    ///
    /// Re-exported from [`crate::record_utils::reference_length`] to avoid duplication.
    pub use crate::record_utils::reference_length;

    /// Simplified CIGAR representation: a vector of (operation kind, length) pairs.
    pub type SimplifiedCigar = Vec<(Kind, usize)>;

    /// Simplifies a CIGAR string by:
    /// 1. Converting S, EQ, X, H operations to M (match)
    /// 2. Coalescing adjacent operations of the same type
    ///
    /// This is useful for comparing read alignments where only indel positions matter.
    #[must_use]
    pub fn simplify_cigar(cigar: &noodles::sam::alignment::record_buf::Cigar) -> SimplifiedCigar {
        let mut simplified: SimplifiedCigar = Vec::new();

        for op in cigar.as_ref() {
            let (kind, len) = (op.kind(), op.len());

            // Convert S, EQ, X, H to M; keep I, D, and M as-is
            let new_kind = match kind {
                Kind::SoftClip | Kind::SequenceMatch | Kind::SequenceMismatch | Kind::HardClip => {
                    Kind::Match
                }
                _ => kind,
            };

            // Coalesce adjacent operations of the same type
            if let Some((last_kind, last_len)) = simplified.last_mut()
                && *last_kind == new_kind
            {
                *last_len += len;
                continue;
            }

            simplified.push((new_kind, len));
        }

        simplified
    }

    /// Checks if `cigar_a` is a prefix of `cigar_b`.
    ///
    /// This allows for reads of different lengths to be grouped together if they share
    /// the same alignment pattern (same indel positions and lengths).
    ///
    /// For example, 10M is a prefix of 20M, but 10M1I is not a prefix of 10M1D.
    #[must_use]
    pub fn is_cigar_prefix(a: &[(Kind, usize)], b: &[(Kind, usize)]) -> bool {
        if a.len() > b.len() {
            return false;
        }

        let last_index = a.len().saturating_sub(1);

        for (i, &(op_a, len_a)) in a.iter().enumerate() {
            let (op_b, len_b) = b[i];
            if op_a != op_b {
                return false;
            }
            // For the last element, a's length can be <= b's length (prefix)
            // For other elements, lengths must match exactly
            if i == last_index {
                if len_a > len_b {
                    return false;
                }
            } else if len_a != len_b {
                return false;
            }
        }
        true
    }
}

/// Test-only adapter that drives [`RawRecordClipper`] through the `RecordBuf` API by
/// round-tripping each record to raw BAM bytes and back.
///
/// The behavioral clipper tests were originally written against the `RecordBuf`-based
/// `SamRecordClipper`, which has been removed (superseded by the raw-byte [`RawRecordClipper`]
/// that every production caller uses). Rather than rewrite each test's assertions into the raw
/// API, this adapter re-encodes the input `RecordBuf` to a `RawRecord`, runs the real
/// [`RawRecordClipper`] method on it, and decodes the result back to a `RecordBuf`, so the
/// tests' setup and assertions run verbatim while exercising the raw clipper's actual logic.
/// The buf→raw→buf round-trip is lossless for the fields the adapter-driven tests assert on
/// (CIGAR, POS, sequence, quality, flags, and *UTF-8* tag values), so a mismatch reflects a real
/// `RawRecordClipper` discrepancy, not a round-trip artifact. Tag values round-trip only when
/// UTF-8: a non-UTF-8 `Z` tag cannot be represented in a noodles `RecordBuf`, so a test needing
/// non-UTF-8 tag bytes must build the `RawRecord` directly and assert on the raw bytes (as
/// `test_auto_clip_attributes_string_non_utf8_preserves_bytes` and
/// `test_upgrade_all_clipping_string_non_utf8_preserves_bytes` do) rather than drive this adapter.
#[cfg(test)]
mod clip_test_adapter {
    use super::{ClippingMode, RawRecordClipper};
    use fgumi_raw_bam::{RawRecord, encode_record_buf_to_raw, raw_record_to_record_buf};
    use noodles::sam::alignment::RecordBuf;

    /// Build a SAM header with a single 100 000-base reference, matching the clipper tests'
    /// coordinate space.
    pub(super) fn test_header() -> noodles::sam::Header {
        use noodles::sam::header::record::value::Map;
        use noodles::sam::header::record::value::map::ReferenceSequence;
        use std::num::NonZeroUsize;
        let ref_seq = Map::<ReferenceSequence>::new(
            NonZeroUsize::new(100_000).expect("ref length must be nonzero"),
        );
        noodles::sam::Header::builder().add_reference_sequence(b"chr1", ref_seq).build()
    }

    /// Encode a `RecordBuf` to a `RawRecord` through the shared test header.
    pub(super) fn to_raw(record: &RecordBuf) -> RawRecord {
        encode_record_buf_to_raw(record, &test_header())
            .expect("encode_record_buf_to_raw should succeed")
    }

    /// Decode a `RawRecord` back to a `RecordBuf` through the shared test header.
    pub(super) fn to_buf(raw: &RawRecord) -> RecordBuf {
        raw_record_to_record_buf(raw, &test_header())
            .expect("raw_record_to_record_buf should succeed")
    }

    /// `RecordBuf`-API façade over [`RawRecordClipper`]. `mode` / `auto_clip_attributes` are
    /// mirrored as fields so the construction-sanity tests can read them like the old struct.
    pub(super) struct RawClipperOnBuf {
        inner: RawRecordClipper,
        pub(super) mode: ClippingMode,
        pub(super) auto_clip_attributes: bool,
    }

    impl RawClipperOnBuf {
        pub(super) fn new(mode: ClippingMode) -> Self {
            Self { inner: RawRecordClipper::new(mode), mode, auto_clip_attributes: false }
        }

        pub(super) fn with_auto_clip(mode: ClippingMode, auto_clip_attributes: bool) -> Self {
            Self {
                inner: RawRecordClipper::with_auto_clip(mode, auto_clip_attributes),
                mode,
                auto_clip_attributes,
            }
        }

        /// Round-trip a single record through the raw clipper: `RecordBuf` → raw → clip → back.
        fn on_one(
            &self,
            record: &mut RecordBuf,
            run: impl FnOnce(&RawRecordClipper, &mut RawRecord) -> usize,
        ) -> usize {
            let mut raw = to_raw(record);
            let ret = run(&self.inner, &mut raw);
            *record = to_buf(&raw);
            ret
        }

        pub(super) fn clip_start_of_alignment(&self, r: &mut RecordBuf, n: usize) -> usize {
            self.on_one(r, |c, raw| c.clip_start_of_alignment(raw, n))
        }

        pub(super) fn clip_end_of_alignment(&self, r: &mut RecordBuf, n: usize) -> usize {
            self.on_one(r, |c, raw| c.clip_end_of_alignment(raw, n))
        }

        pub(super) fn clip_5_prime_end_of_alignment(&self, r: &mut RecordBuf, n: usize) -> usize {
            self.on_one(r, |c, raw| c.clip_5_prime_end_of_alignment(raw, n))
        }

        pub(super) fn clip_3_prime_end_of_alignment(&self, r: &mut RecordBuf, n: usize) -> usize {
            self.on_one(r, |c, raw| c.clip_3_prime_end_of_alignment(raw, n))
        }

        pub(super) fn clip_start_of_read(&self, r: &mut RecordBuf, n: usize) -> usize {
            self.on_one(r, |c, raw| c.clip_start_of_read_raw(raw, n))
        }

        pub(super) fn clip_end_of_read(&self, r: &mut RecordBuf, n: usize) -> usize {
            self.on_one(r, |c, raw| c.clip_end_of_read_raw(raw, n))
        }

        pub(super) fn clip_overlapping_reads(
            &self,
            r1: &mut RecordBuf,
            r2: &mut RecordBuf,
        ) -> (usize, usize) {
            let mut raw1 = to_raw(r1);
            let mut raw2 = to_raw(r2);
            let ret = self.inner.clip_overlapping_reads(&mut raw1, &mut raw2);
            *r1 = to_buf(&raw1);
            *r2 = to_buf(&raw2);
            ret
        }

        pub(super) fn upgrade_all_clipping(
            &self,
            r: &mut RecordBuf,
        ) -> anyhow::Result<(usize, usize)> {
            let mut raw = to_raw(r);
            let ret = self.inner.upgrade_all_clipping_raw(&mut raw)?;
            *r = to_buf(&raw);
            Ok(ret)
        }
    }
}

#[cfg(test)]
#[expect(clippy::similar_names, reason = "test code uses clipper/clipped naming pattern")]
mod tests {
    use super::clip_test_adapter::{RawClipperOnBuf, to_raw};
    use super::*;
    use crate::builder::RecordBuilder;
    use crate::record_utils;
    use fgumi_dna::{MIN_PHRED, NO_CALL_BASE};
    use noodles::sam::alignment::record::cigar::Cigar;
    use noodles::sam::alignment::record::cigar::Op as CigarOp;
    use noodles::sam::alignment::record_buf::Cigar as CigarBuf;
    use noodles::sam::alignment::record_buf::data::field::Value;
    use proptest::prelude::*;
    use rstest::rstest;

    #[test]
    fn test_clipping_mode() {
        let clipper = RawClipperOnBuf::new(ClippingMode::Soft);
        assert_eq!(clipper.mode, ClippingMode::Soft);

        let clipper = RawClipperOnBuf::new(ClippingMode::Hard);
        assert_eq!(clipper.mode, ClippingMode::Hard);
    }

    #[test]
    fn test_auto_clip() {
        let clipper = RawClipperOnBuf::with_auto_clip(ClippingMode::Soft, true);
        assert!(clipper.auto_clip_attributes);

        let clipper_disabled = RawClipperOnBuf::with_auto_clip(ClippingMode::Hard, false);
        assert!(!clipper_disabled.auto_clip_attributes);
    }

    /// `raw_clip_update_array_tag` must update all seven BAM array subtypes (c/C/s/S/i/I/f)
    /// so hard-clipping a per-base `B:c` or `B:I` tag trims it to the new read length,
    /// matching typed-path behavior.
    #[test]
    fn test_raw_clip_update_array_tag_handles_all_subtypes() {
        let mut rec = fgumi_raw_bam::testutil::make_bam_bytes(0, 0, 0, b"r", &[], 0, -1, -1, &[]);

        // Append a B:c tag "XC" with 4 elements [-5, -4, -3, -2].
        {
            let mut ed = fgumi_raw_bam::RawTagsEditor::from_vec(&mut rec);
            ed.append_array_i8(b"XC", &[-5i8, -4, -3, -2]);
        }

        // Trim to 2 elements (bytes [-5, -4] interpreted as i8-wide bytes).
        {
            let mut ed = fgumi_raw_bam::RawTagsEditor::from_vec(&mut rec);
            let truncated: [u8; 2] = [(-5i8).cast_unsigned(), (-4i8).cast_unsigned()];
            raw_clip_update_array_tag(&mut ed, *b"XC", b'c', &truncated);
        }
        let aux_off = fgumi_raw_bam::aux_data_offset_from_record(&rec).expect("aux present");
        let arr = fgumi_raw_bam::find_array_tag(&rec[aux_off..], b"XC").expect("XC present");
        assert_eq!(arr.elem_type, b'c', "B:c subtype preserved");
        assert_eq!(arr.count, 2, "B:c count trimmed");
        let decoded: Vec<i8> = arr.data.iter().map(|&b| b.cast_signed()).collect();
        assert_eq!(decoded, vec![-5i8, -4]);

        // Same drill for B:I: append 4 u32 values then trim to 3.
        {
            let mut ed = fgumi_raw_bam::RawTagsEditor::from_vec(&mut rec);
            ed.append_array_u32(b"XI", &[10u32, 20, 30, 40]);
        }
        {
            let mut ed = fgumi_raw_bam::RawTagsEditor::from_vec(&mut rec);
            let mut truncated = Vec::with_capacity(12);
            for v in [10u32, 20, 30] {
                truncated.extend_from_slice(&v.to_le_bytes());
            }
            raw_clip_update_array_tag(&mut ed, *b"XI", b'I', &truncated);
        }
        let aux_off = fgumi_raw_bam::aux_data_offset_from_record(&rec).expect("aux present");
        let arr = fgumi_raw_bam::find_array_tag(&rec[aux_off..], b"XI").expect("XI present");
        assert_eq!(arr.elem_type, b'I', "B:I subtype preserved");
        assert_eq!(arr.count, 3, "B:I count trimmed");
        let decoded: Vec<u32> = arr
            .data
            .chunks_exact(4)
            .map(|c| u32::from_le_bytes([c[0], c[1], c[2], c[3]]))
            .collect();
        assert_eq!(decoded, vec![10u32, 20, 30]);
    }

    use noodles::core::Position;
    use noodles::sam::alignment::record::Cigar as CigarTrait;
    use noodles::sam::alignment::record::Flags;
    use noodles::sam::alignment::record_buf::RecordBuf;

    /// Helper to format a CIGAR string for comparison
    fn format_cigar(cigar: &impl CigarTrait) -> String {
        use std::fmt::Write;
        cigar.iter().filter_map(Result::ok).fold(String::new(), |mut acc, op| {
            let kind_char = match op.kind() {
                Kind::Match => 'M',
                Kind::Insertion => 'I',
                Kind::Deletion => 'D',
                Kind::Skip => 'N',
                Kind::SoftClip => 'S',
                Kind::HardClip => 'H',
                Kind::Pad => 'P',
                Kind::SequenceMatch => '=',
                Kind::SequenceMismatch => 'X',
            };
            let _ = write!(acc, "{}{}", op.len(), kind_char);
            acc
        })
    }

    /// Helper to create a simple test record with given CIGAR, sequence, and position
    /// Query-consuming (SEQ) length of a CIGAR string: the sum of `M`/`I`/`S`/`=`/`X` run
    /// lengths. Used to pad placeholder sequences to a length the raw BAM encoder accepts.
    fn cigar_query_len(cigar_str: &str) -> usize {
        let mut total = 0usize;
        let mut num = 0usize;
        for ch in cigar_str.chars() {
            if let Some(d) = ch.to_digit(10) {
                num = num * 10 + d as usize;
            } else {
                if matches!(ch, 'M' | 'I' | 'S' | '=' | 'X') {
                    total += num;
                }
                num = 0;
            }
        }
        total
    }

    fn create_test_record(cigar_str: &str, seq: &str, start_pos: usize) -> RecordBuf {
        // The raw BAM path (and the spec) require SEQ length to equal the CIGAR's query length.
        // Many clipper tests pass a placeholder sequence (too short or too long) and only assert
        // on the resulting CIGAR/POS, so normalize the sequence to exactly the query length --
        // keeping any provided prefix and padding the rest with `A` -- to make the record a valid
        // BAM record without changing the behavior under test. A content-asserting test already
        // supplies a query-length sequence, so it is left untouched.
        let query_len = cigar_query_len(cigar_str);
        let mut seq: String = seq.chars().take(query_len).collect();
        if seq.len() < query_len {
            seq.push_str(&"A".repeat(query_len - seq.len()));
        }
        RecordBuilder::mapped_read()
            .sequence(&seq)
            .cigar(cigar_str)
            .alignment_start(start_pos)
            .build()
    }

    // ===================================================================
    // clip_5_prime (clipStartOfAlignment) tests - Basic functionality
    // ===================================================================

    #[test]
    fn test_clip_start_of_alignment_soft_matched_bases() {
        let clipper = RawClipperOnBuf::new(ClippingMode::Soft);
        let mut record = create_test_record("50M", "ACGTACGTACGT", 10);

        let clipped = clipper.clip_start_of_alignment(&mut record, 10);
        assert_eq!(clipped, 10);
        assert_eq!(record.alignment_start(), Position::new(20));

        let cigar_str = format_cigar(&record.cigar());
        assert_eq!(cigar_str, "10S40M");
    }

    #[test]
    fn test_clip_start_of_alignment_soft_with_insertion() {
        let clipper = RawClipperOnBuf::new(ClippingMode::Soft);
        let mut record = create_test_record("4M2I44M", "ACGTACGTACGT", 10);

        let clipped = clipper.clip_start_of_alignment(&mut record, 10);
        assert_eq!(clipped, 10);
        assert_eq!(record.alignment_start(), Position::new(18)); // 10 + (4+2 ref consumed)

        let cigar_str = format_cigar(&record.cigar());
        assert_eq!(cigar_str, "10S40M");
    }

    #[test]
    fn test_clip_start_of_alignment_soft_with_deletion() {
        let clipper = RawClipperOnBuf::new(ClippingMode::Soft);
        let mut record = create_test_record("6M2D44M", "ACGTACGTACGT", 10);

        let clipped = clipper.clip_start_of_alignment(&mut record, 10);
        assert_eq!(clipped, 10);
        assert_eq!(record.alignment_start(), Position::new(22)); // 10 + 6 + 2 (deletion) + 4

        let cigar_str = format_cigar(&record.cigar());
        assert_eq!(cigar_str, "10S40M");
    }

    #[test]
    fn test_clip_start_of_alignment_soft_additional_bases() {
        let clipper = RawClipperOnBuf::new(ClippingMode::Soft);
        let mut record = create_test_record("10S40M", "ACGTACGTACGT", 10);

        let clipped = clipper.clip_start_of_alignment(&mut record, 10);
        assert_eq!(clipped, 10);
        assert_eq!(record.alignment_start(), Position::new(20));

        let cigar_str = format_cigar(&record.cigar());
        assert_eq!(cigar_str, "20S30M");
    }

    #[test]
    fn test_clip_start_of_alignment_preserve_hard_clip() {
        let clipper = RawClipperOnBuf::new(ClippingMode::Soft);
        let mut record = create_test_record("10H40M", "ACGTACGTACGT", 10);

        let clipped = clipper.clip_start_of_alignment(&mut record, 10);
        assert_eq!(clipped, 10);
        assert_eq!(record.alignment_start(), Position::new(20));

        let cigar_str = format_cigar(&record.cigar());
        assert_eq!(cigar_str, "10H10S30M");
    }

    #[test]
    fn test_clip_start_of_alignment_complex_cigar() {
        let clipper = RawClipperOnBuf::new(ClippingMode::Soft);
        // 2H4S16M10I5M5I10M
        let mut record = create_test_record("2H4S16M10I5M5I10M", "ACGTACGTACGT", 10);

        let clipped = clipper.clip_start_of_alignment(&mut record, 10);
        assert_eq!(clipped, 10);
        assert_eq!(record.alignment_start(), Position::new(20));

        let cigar_str = format_cigar(&record.cigar());
        assert_eq!(cigar_str, "2H14S6M10I5M5I10M");
    }

    #[test]
    fn test_clip_start_of_alignment_consumes_trailing_insertion() {
        let clipper = RawClipperOnBuf::new(ClippingMode::Soft);
        let mut record = create_test_record("8M4I38M", "ACGTACGTACGT", 10);

        // Ask to clip 10 bases, but should consume 12 bases (8M + 4I) because
        // the insertion at the clip boundary is consumed entirely
        let clipped = clipper.clip_start_of_alignment(&mut record, 10);
        assert_eq!(clipped, 12);
        assert_eq!(record.alignment_start(), Position::new(18));

        let cigar_str = format_cigar(&record.cigar());
        assert_eq!(cigar_str, "12S38M");
    }

    #[test]
    fn test_clip_start_of_alignment_preserve_insertion_after_clip() {
        let clipper = RawClipperOnBuf::new(ClippingMode::Soft);
        let mut record = create_test_record("10M4I36M", "ACGTACGTACGT", 10);

        let clipped = clipper.clip_start_of_alignment(&mut record, 10);
        assert_eq!(clipped, 10);
        assert_eq!(record.alignment_start(), Position::new(20));

        let cigar_str = format_cigar(&record.cigar());
        assert_eq!(cigar_str, "10S4I36M");
    }

    #[test]
    fn test_clip_start_of_alignment_remove_deletion_after_clip() {
        let clipper = RawClipperOnBuf::new(ClippingMode::Soft);
        let mut record = create_test_record("10M4D40M", "ACGTACGTACGT", 10);

        let clipped = clipper.clip_start_of_alignment(&mut record, 10);
        assert_eq!(clipped, 10);
        assert_eq!(record.alignment_start(), Position::new(24)); // 10 + 10 + 4 (deletion)

        let cigar_str = format_cigar(&record.cigar());
        assert_eq!(cigar_str, "10S40M");
    }

    #[test]
    fn test_clip_start_of_alignment_preserve_distant_deletion() {
        let clipper = RawClipperOnBuf::new(ClippingMode::Soft);
        let mut record = create_test_record("25M4D25M", "ACGTACGTACGT", 10);

        let clipped = clipper.clip_start_of_alignment(&mut record, 10);
        assert_eq!(clipped, 10);
        assert_eq!(record.alignment_start(), Position::new(20));

        let cigar_str = format_cigar(&record.cigar());
        assert_eq!(cigar_str, "10S15M4D25M");
    }

    #[test]
    fn test_clip_start_of_alignment_soft_with_mask() {
        let clipper = RawClipperOnBuf::new(ClippingMode::SoftWithMask);
        let seq = "ACGTACGTAC"; // 10 bases
        let mut record = create_test_record("10M", seq, 10);

        let clipped = clipper.clip_start_of_alignment(&mut record, 5);
        assert_eq!(clipped, 5);

        let cigar_str = format_cigar(&record.cigar());
        assert_eq!(cigar_str, "5S5M");

        // Check first 5 bases are masked to N
        let bases = record.sequence();
        for i in 0..5 {
            assert_eq!(bases.as_ref()[i], b'N', "Base at position {i} should be N");
        }

        // Check first 5 quals are set to min
        let quals = record.quality_scores();
        for i in 0..5 {
            assert_eq!(quals.as_ref()[i], 2, "Quality at position {i} should be 2");
        }
    }

    #[test]
    fn test_clip_start_of_alignment_soft_with_mask_existing_soft_clip() {
        let clipper = RawClipperOnBuf::new(ClippingMode::SoftWithMask);
        let seq = "ACGTACGTACGTACGTACGT"; // 20 bases
        let mut record = create_test_record("10S10M", seq, 10);

        // Only 10 clippable (aligned) bases; requesting 10 consumes the whole alignment,
        // so fgbio unmaps the read (SamRecordClipper.scala:111-114) rather than masking.
        let clipped = clipper.clip_start_of_alignment(&mut record, 10);
        assert_eq!(clipped, 10);

        assert!(record.flags().is_unmapped(), "read should be unmapped");
        assert_eq!(record.cigar().iter().count(), 0, "CIGAR should be empty");

        // Bases are preserved (NOT masked) on unmap — makeReadUnmapped keeps seq/qual.
        assert_eq!(
            record.sequence().as_ref(),
            seq.as_bytes(),
            "bases should be preserved, not masked"
        );
    }

    #[test]
    fn test_clip_start_of_alignment_unmapped_read_does_nothing() {
        let clipper = RawClipperOnBuf::new(ClippingMode::Soft);
        let mut record = create_test_record("50M", "ACGTACGTACGTACGTACGT", 1000);

        // Set unmapped flag
        *record.flags_mut() = Flags::UNMAPPED;

        let clipped = clipper.clip_start_of_alignment(&mut record, 10);

        // Should not clip unmapped reads
        assert_eq!(clipped, 0);
        assert_eq!(record.alignment_start(), Position::new(1000));

        // CIGAR should remain unchanged
        let cigar_str = format_cigar(&record.cigar());
        assert_eq!(cigar_str, "50M");
    }

    #[test]
    fn test_clip_start_of_alignment_soft() {
        let clipper = RawClipperOnBuf::new(ClippingMode::Soft);
        let mut record = create_test_record("100M", "ACGTACGTACGT", 1000);

        let clipped = clipper.clip_start_of_alignment(&mut record, 5);
        assert_eq!(clipped, 5);

        // Check that CIGAR starts with soft clip
        let cigar = record.cigar();
        let first_op = cigar
            .iter()
            .next()
            .expect("CIGAR should have at least one op")
            .expect("failed to parse CIGAR op");
        assert_eq!(first_op.kind(), Kind::SoftClip);
        assert_eq!(first_op.len(), 5);
    }

    #[test]
    fn test_clip_end_of_alignment_soft() {
        let clipper = RawClipperOnBuf::new(ClippingMode::Soft);
        let mut record = create_test_record("100M", "ACGTACGTACGT", 1000);

        let clipped = clipper.clip_end_of_alignment(&mut record, 5);
        assert_eq!(clipped, 5);

        // Check that CIGAR ends with soft clip
        let cigar = record.cigar();
        let ops: Vec<_> = cigar.iter().filter_map(Result::ok).collect();
        let last_op = ops.last().expect("CIGAR should have at least one op");
        assert_eq!(last_op.kind(), Kind::SoftClip);
        assert_eq!(last_op.len(), 5);
    }

    // ===================================================================
    // clip_3_prime (clipEndOfAlignment) tests - Extended coverage
    // ===================================================================

    #[test]
    fn test_clip_end_of_alignment_soft_matched_bases() {
        let clipper = RawClipperOnBuf::new(ClippingMode::Soft);
        // 50 bases to match 50M CIGAR
        let seq = "ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTAC";
        assert_eq!(seq.len(), 50, "Sequence length should be 50");
        let mut record = create_test_record("50M", seq, 10);

        let clipped = clipper.clip_end_of_alignment(&mut record, 10);
        assert_eq!(clipped, 10);
        assert_eq!(record.alignment_start(), Position::new(10));

        let cigar_str = format_cigar(&record.cigar());
        assert_eq!(cigar_str, "40M10S");
    }

    #[test]
    fn test_clip_end_of_alignment_hard_existing_soft_clip() {
        let clipper = RawClipperOnBuf::new(ClippingMode::Hard);
        let seq = "ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTAC"; // 50 bases
        assert_eq!(seq.len(), 50, "Sequence length should be 50");
        let mut record = create_test_record("40M10S", seq, 10);

        let orig_seq: Vec<u8> = record.sequence().as_ref().to_vec();
        let orig_qual: Vec<u8> = record.quality_scores().as_ref().to_vec();

        let clipped = clipper.clip_end_of_alignment(&mut record, 10);
        assert_eq!(clipped, 10);
        assert_eq!(record.alignment_start(), Position::new(10));

        let cigar_str = format_cigar(&record.cigar());
        assert_eq!(cigar_str, "30M20H");

        // Check sequence and qualities are truncated
        assert_eq!(record.sequence().as_ref(), &orig_seq[..30]);
        assert_eq!(record.quality_scores().as_ref(), &orig_qual[..30]);
    }

    #[test]
    fn test_clip_start_of_alignment_hard() {
        let clipper = RawClipperOnBuf::new(ClippingMode::Hard);
        let mut record = create_test_record("12M", "ACGTACGTACGT", 1000);

        let original_len = record.sequence().len();
        let clipped = clipper.clip_start_of_alignment(&mut record, 5);
        assert_eq!(clipped, 5);

        // Check that sequence is shortened
        assert_eq!(record.sequence().len(), original_len - 5);

        // Check that CIGAR starts with hard clip
        let cigar = record.cigar();
        let first_op = cigar
            .iter()
            .next()
            .expect("CIGAR should have at least one op")
            .expect("failed to parse CIGAR op");
        assert_eq!(first_op.kind(), Kind::HardClip);
        assert_eq!(first_op.len(), 5);
    }

    #[test]
    fn test_clip_overlapping_reads_no_overlap() {
        let clipper = RawClipperOnBuf::new(ClippingMode::Soft);

        // R1: 1000-1100, R2: 1200-1300 (no overlap)
        let mut r1 = create_test_record("100M", "ACGTACGTACGT", 1000);
        let mut r2 = create_test_record("100M", "TGCATGCATGCA", 1200);

        let (clipped_r1, clipped_r2) = clipper.clip_overlapping_reads(&mut r1, &mut r2);

        // No clipping should occur
        assert_eq!(clipped_r1, 0);
        assert_eq!(clipped_r2, 0);
    }

    #[test]
    fn test_clip_overlapping_reads_with_overlap() {
        let clipper = RawClipperOnBuf::new(ClippingMode::Soft);

        // R1: 1000-1099 (forward), R2: 1050-1149 (reverse) - 50bp overlap
        // Create proper FR pair
        let mut r1 = create_paired_record("100M", "ACGTACGTACGT", 1000, false, true, 1050, "100M");
        let mut r2 = create_paired_record("100M", "TGCATGCATGCA", 1050, true, false, 1000, "100M");

        let (clipped_r1, clipped_r2) = clipper.clip_overlapping_reads(&mut r1, &mut r2);

        // With midpoint algorithm, both reads should be clipped
        // Midpoint = (1000 + 1149) / 2 = 1074
        // R1 clips from 1075 onwards (approximately 25 bases)
        // R2 clips from start to 1074 (approximately 25 bases)
        assert!(clipped_r1 > 0, "R1 should be clipped");
        assert!(clipped_r2 > 0, "R2 should be clipped");

        // Both should clip roughly equal amounts (within a few bases due to rounding)
        let diff = (i32::try_from(clipped_r1).expect("clipped_r1 should fit in i32")
            - i32::try_from(clipped_r2).expect("clipped_r2 should fit in i32"))
        .abs();
        assert!(
            diff <= 2,
            "Clipping should be roughly equal, but got {clipped_r1} vs {clipped_r2}"
        );
    }

    #[test]
    fn test_cigar_utils_reference_length() {
        let record = create_test_record("50M10I40M", "ACGTACGTACGT", 1000);
        let ref_len = cigar_utils::reference_length(&record.cigar());

        // 50M + 40M = 90 (insertion doesn't consume reference)
        assert_eq!(ref_len, 90);
    }

    #[test]
    fn test_cigar_utils_aligned_bases() {
        let record = create_test_record("50M10I40M", "ACGTACGTACGT", 1000);
        let aligned = cigar_utils::aligned_bases(&record.cigar());

        // 50M + 40M = 90 (insertion is not aligned)
        assert_eq!(aligned, 90);
    }

    #[test]
    fn test_cigar_utils_clipped_bases() {
        let record = create_test_record("5S90M5S", "ACGTACGTACGT", 1000);
        let clipped = cigar_utils::clipped_bases(&record.cigar());

        // 5S + 5S = 10
        assert_eq!(clipped, 10);
    }

    #[test]
    fn test_clip_start_of_alignment_with_existing_soft_clip() {
        let clipper = RawClipperOnBuf::new(ClippingMode::Soft);
        let mut record = create_test_record("5S95M", "ACGTACGTACGT", 1000);

        // Clip additional 3 bases from 5' end
        let clipped = clipper.clip_start_of_alignment(&mut record, 3);
        assert_eq!(clipped, 3);

        // Total soft clip should be 5 + 3 = 8
        let cigar = record.cigar();
        let first_op = cigar
            .iter()
            .next()
            .expect("CIGAR should have at least one op")
            .expect("failed to parse CIGAR op");
        assert_eq!(first_op.kind(), Kind::SoftClip);
    }

    #[test]
    fn test_clip_entire_read() {
        let clipper = RawClipperOnBuf::new(ClippingMode::Soft);
        let seq = "ACGTACGTACGT";
        let mut record = create_test_record("12M", seq, 1000);

        // Try to clip more than the read length
        let clipped = clipper.clip_start_of_alignment(&mut record, 20);

        // Should only clip up to sequence length
        assert_eq!(clipped, seq.len());
    }

    #[test]
    fn test_clip_zero_bases() {
        let clipper = RawClipperOnBuf::new(ClippingMode::Soft);
        let mut record = create_test_record("100M", "ACGTACGTACGT", 1000);

        let clipped = clipper.clip_start_of_alignment(&mut record, 0);
        assert_eq!(clipped, 0);

        // CIGAR should remain unchanged
        let cigar = record.cigar();
        let first_op = cigar
            .iter()
            .next()
            .expect("CIGAR should have at least one op")
            .expect("failed to parse CIGAR op");
        assert_eq!(first_op.kind(), Kind::Match);
    }

    /// Helper to create a paired read with specific flags and positions
    fn create_paired_record(
        cigar_str: &str,
        seq: &str,
        start_pos: usize,
        is_reverse: bool,
        mate_reverse: bool,
        mate_pos: usize,
        mate_cigar_str: &str,
    ) -> RecordBuf {
        let mut record = create_test_record(cigar_str, seq, start_pos);

        // Set paired flags
        let mut flags = Flags::SEGMENTED;
        if is_reverse {
            flags |= Flags::REVERSE_COMPLEMENTED;
        }
        if mate_reverse {
            flags |= Flags::MATE_REVERSE_COMPLEMENTED;
        }
        *record.flags_mut() = flags;

        // Set mate position
        *record.mate_alignment_start_mut() = Position::new(mate_pos);

        // Set same reference for both reads
        *record.mate_reference_sequence_id_mut() = record.reference_sequence_id();

        // Calculate and set template_length for FR pair detection
        // Parse mate CIGAR to calculate its reference length
        let mate_cigar_ops: Vec<CigarOp> = mate_cigar_str
            .split(|c: char| !c.is_numeric())
            .filter(|s| !s.is_empty())
            .zip(mate_cigar_str.chars().filter(|c| c.is_alphabetic()))
            .map(|(len_str, kind_char)| {
                let len: usize = len_str.parse().expect("failed to parse CIGAR operation length");
                let kind = match kind_char {
                    'M' => Kind::Match,
                    'I' => Kind::Insertion,
                    'D' => Kind::Deletion,
                    'S' => Kind::SoftClip,
                    'H' => Kind::HardClip,
                    'N' => Kind::Skip,
                    _ => panic!("Invalid CIGAR operation: {kind_char}"),
                };
                CigarOp::new(kind, len)
            })
            .collect();
        let mate_cigar = CigarBuf::from(mate_cigar_ops);
        let mate_ref_len = cigar_utils::reference_length(&mate_cigar);
        let ref_len = cigar_utils::reference_length(&record.cigar());

        let tlen = if is_reverse {
            // This read is reverse, mate is forward
            // tlen = -(this_end - mate_start + 1)
            let this_end = start_pos + ref_len - 1;
            if this_end >= mate_pos {
                -i32::try_from(this_end - mate_pos + 1).expect("template length should fit in i32")
            } else {
                i32::try_from(mate_pos - this_end - 1).expect("template length should fit in i32")
            }
        } else {
            // This read is forward, mate is reverse
            // tlen = mate_end - this_start + 1
            let mate_end = mate_pos + mate_ref_len - 1;
            if mate_end >= start_pos {
                i32::try_from(mate_end - start_pos + 1).expect("template length should fit in i32")
            } else {
                -i32::try_from(start_pos - mate_end - 1).expect("template length should fit in i32")
            }
        };

        *record.template_length_mut() = tlen;

        record
    }

    #[test]
    fn test_is_fr_pair_valid() {
        // Create valid FR pair: R1 forward, R2 reverse
        let r1 = create_paired_record("100M", "ACGTACGTACGT", 1000, false, true, 1100, "100M");
        let r2 = create_paired_record("100M", "TGCATGCATGCA", 1100, true, false, 1000, "100M");

        assert!(record_utils::is_fr_pair(&r1, &r2));
    }

    #[test]
    fn test_is_fr_pair_both_forward() {
        // FF orientation - should not be considered FR
        let r1 = create_paired_record("100M", "ACGTACGTACGT", 1000, false, false, 1100, "100M");
        let r2 = create_paired_record("100M", "TGCATGCATGCA", 1100, false, false, 1000, "100M");

        assert!(!record_utils::is_fr_pair(&r1, &r2));
    }

    #[test]
    fn test_is_fr_pair_both_reverse() {
        // RR orientation - should not be considered FR
        let r1 = create_paired_record("100M", "ACGTACGTACGT", 1000, true, true, 1100, "100M");
        let r2 = create_paired_record("100M", "TGCATGCATGCA", 1100, true, true, 1000, "100M");

        assert!(!record_utils::is_fr_pair(&r1, &r2));
    }

    #[test]
    fn test_is_fr_pair_rf_orientation() {
        // RF orientation (reverse of FR) - should not be considered FR
        let r1 = create_paired_record("100M", "ACGTACGTACGT", 1000, true, false, 1100, "100M");
        let r2 = create_paired_record("100M", "TGCATGCATGCA", 1100, false, true, 1000, "100M");

        assert!(!record_utils::is_fr_pair(&r1, &r2));
    }

    #[test]
    fn test_is_fr_pair_unmapped() {
        let mut r1 = create_paired_record("100M", "ACGTACGTACGT", 1000, false, true, 1100, "100M");
        let r2 = create_paired_record("100M", "TGCATGCATGCA", 1100, true, false, 1000, "150M");

        // Set R1 as unmapped
        *r1.flags_mut() |= Flags::UNMAPPED;

        assert!(!record_utils::is_fr_pair(&r1, &r2));
    }

    #[test]
    fn test_is_fr_pair_mate_unmapped() {
        let mut r1 = create_paired_record("100M", "ACGTACGTACGT", 1000, false, true, 1100, "100M");
        let r2 = create_paired_record("100M", "TGCATGCATGCA", 1100, true, false, 1000, "150M");

        // Set R1's mate as unmapped
        *r1.flags_mut() |= Flags::MATE_UNMAPPED;

        assert!(!record_utils::is_fr_pair(&r1, &r2));
    }

    #[test]
    fn test_is_fr_pair_different_chromosomes() {
        let r1 = create_paired_record("100M", "ACGTACGTACGT", 1000, false, true, 1100, "100M");
        let mut r2 = create_paired_record("100M", "TGCATGCATGCA", 1100, true, false, 1000, "150M");

        // Set R2 to different reference sequence
        *r2.reference_sequence_id_mut() = Some(1);

        assert!(!record_utils::is_fr_pair(&r1, &r2));
    }

    #[test]
    fn test_is_fr_pair_not_paired() {
        let mut r1 = create_paired_record("100M", "ACGTACGTACGT", 1000, false, true, 1100, "100M");
        let r2 = create_paired_record("100M", "TGCATGCATGCA", 1100, true, false, 1000, "150M");

        // Remove paired flag from R1
        *r1.flags_mut() = Flags::empty();

        assert!(!record_utils::is_fr_pair(&r1, &r2));
    }

    #[test]
    fn test_clip_overlapping_reads_non_fr_pair() {
        let clipper = RawClipperOnBuf::new(ClippingMode::Soft);

        // Create FF pair (both forward) with overlap
        // R1: 1000-1100, R2: 1050-1150 (50bp overlap)
        let mut r1 = create_paired_record("100M", "ACGTACGTACGT", 1000, false, false, 1050, "100M");
        let mut r2 = create_paired_record("100M", "TGCATGCATGCA", 1050, false, false, 1000, "100M");

        let (clipped_r1, clipped_r2) = clipper.clip_overlapping_reads(&mut r1, &mut r2);

        // Should NOT clip because this is not an FR pair
        assert_eq!(clipped_r1, 0);
        assert_eq!(clipped_r2, 0);
    }

    /// Inconsistent pair — r1 flags look FR but r2 itself is unmapped — must not be clipped.
    /// Raw path must validate both records, not just r1.
    #[test]
    fn test_raw_clip_overlapping_rejects_inconsistent_pair_when_r2_is_unmapped() {
        use fgumi_raw_bam::encode_record_buf_to_raw;
        use noodles::sam::header::record::value::Map;
        use noodles::sam::header::record::value::map::ReferenceSequence;
        use std::num::NonZeroUsize;

        // Build a valid FR pair, then mark r2 as self-unmapped. r1's flags still claim the mate
        // is mapped, but the pair is not a valid primary FR pair because r2 itself is unmapped,
        // so the symmetric gate must reject it (validating both records, not just r1).
        let r1_buf =
            create_paired_record("100M", &"A".repeat(100), 1000, false, true, 1100, "100M");
        let mut r2_buf =
            create_paired_record("100M", &"T".repeat(100), 1100, true, false, 1000, "100M");
        *r2_buf.flags_mut() |= Flags::UNMAPPED;

        let ref_seq = Map::<ReferenceSequence>::new(
            NonZeroUsize::new(100_000).expect("ref length must be nonzero"),
        );
        let header =
            noodles::sam::Header::builder().add_reference_sequence(b"chr1", ref_seq).build();
        let mut r1 = encode_record_buf_to_raw(&r1_buf, &header).expect("encode r1");
        let mut r2 = encode_record_buf_to_raw(&r2_buf, &header).expect("encode r2");

        // Not a valid primary FR pair: r2 is self-unmapped.
        assert!(!fgumi_raw_bam::is_primary_fr_pair_raw(r1.as_ref(), r2.as_ref()));

        let clipper = RawRecordClipper::new(ClippingMode::Soft);
        assert_eq!(clipper.clip_overlapping_reads(&mut r1, &mut r2), (0, 0));
        assert_eq!(clipper.clip_extending_past_mate_ends(&mut r1, &mut r2), (0, 0));
    }

    /// Retro sibling of #839: `clip_overlapping_reads` must clip a valid dovetail FR pair whose
    /// FORWARD read carries a leftmost-to-rightmost `TLEN`. `is_fr_pair_raw`'s per-record
    /// forward arm derives the mate 5' from `TLEN`, so it mis-reports the forward read as non-FR
    /// on such dovetails (htsjdk/samtools#1771); the old `is_fr_pair_raw(r1) || is_fr_pair_raw(r2)`
    /// gate then short-circuited to `(0, 0)` and left the read-through/adapter overlap unclipped.
    /// The symmetric `is_primary_fr_pair_raw(r1, r2)` gate — both records in hand, the same fix
    /// the sibling `clip_extending_past_mate_ends` already uses — classifies the pair correctly.
    #[test]
    fn test_raw_clip_overlapping_reads_dovetail_fr_forward_tlen() {
        use fgumi_raw_bam::encode_record_buf_to_raw;
        use noodles::sam::header::record::value::Map;
        use noodles::sam::header::record::value::map::ReferenceSequence;
        use std::num::NonZeroUsize;

        // Dovetail FR pair: forward 100M @ 101 (ref 101..200), reverse 100M @ 61 (ref 61..160).
        // Overlap is ref 101..160; the reference midpoint between the forward 5' (101) and the
        // reverse 3' (160) is 130, so each read clips 70 query bases back to the midpoint.
        let seq = "A".repeat(100);
        let mut fwd_buf = create_paired_record("100M", &seq, 101, false, true, 61, "100M");
        // An aligner writing TLEN under SAM v1 §1.4's leftmost-to-rightmost convention emits a
        // negative TLEN on the rightmost (here forward) read, which is what mis-drives
        // is_fr_pair_raw's per-record forward arm. create_paired_record's own 5'-to-5' TLEN would
        // sidestep the bug, so override it to the leftmost-to-rightmost value (span 61..200).
        *fwd_buf.template_length_mut() = -140;
        let rev_buf = create_paired_record("100M", &seq, 61, true, false, 101, "100M");

        let ref_seq =
            Map::<ReferenceSequence>::new(NonZeroUsize::new(100_000).expect("ref length nonzero"));
        let header =
            noodles::sam::Header::builder().add_reference_sequence(b"chr1", ref_seq).build();
        let mut fwd = encode_record_buf_to_raw(&fwd_buf, &header).expect("encode fwd");
        let mut rev = encode_record_buf_to_raw(&rev_buf, &header).expect("encode rev");

        // The pair is a valid primary FR pair (classified on the reverse record's CIGAR arm),
        // even though the forward read's per-record TLEN arm misclassifies it.
        assert!(fgumi_raw_bam::is_primary_fr_pair_raw(fwd.as_ref(), rev.as_ref()));

        // Before the fix this returned (0, 0); the dovetail overlap must be clipped to the
        // reference midpoint on both reads.
        let clipper = RawRecordClipper::new(ClippingMode::Soft);
        assert_eq!(clipper.clip_overlapping_reads(&mut fwd, &mut rev), (70, 70));
    }

    #[test]
    fn test_auto_clip_attributes_string_5_prime() {
        use noodles::sam::alignment::record::data::field::Tag;

        let clipper = RawClipperOnBuf::with_auto_clip(ClippingMode::Hard, true);
        let mut record = create_test_record("10M", "ACGTACGTAC", 1000);

        // Add a string attribute that matches the read length
        let tag = Tag::from([b'X', b'S']);
        let value = Value::from("0123456789");
        record.data_mut().insert(tag, value);

        // Clip 3 bases from 5' end
        let clipped = clipper.clip_start_of_alignment(&mut record, 3);
        assert_eq!(clipped, 3);

        // Check that the attribute was clipped
        if let Some(Value::String(s)) = record.data().get(&tag) {
            let bytes: &[u8] = s.as_ref();
            assert_eq!(bytes, b"3456789");
        } else {
            panic!("Tag XS not found or wrong type");
        }
    }

    #[test]
    fn test_auto_clip_attributes_string_3_prime() {
        use noodles::sam::alignment::record::data::field::Tag;

        let clipper = RawClipperOnBuf::with_auto_clip(ClippingMode::Hard, true);
        let mut record = create_test_record("10M", "ACGTACGTAC", 1000);

        // Add a string attribute that matches the read length
        let tag = Tag::from([b'X', b'S']);
        let value = Value::from("0123456789");
        record.data_mut().insert(tag, value);

        // Clip 3 bases from 3' end
        let clipped = clipper.clip_end_of_alignment(&mut record, 3);
        assert_eq!(clipped, 3);

        // Check that the attribute was clipped
        if let Some(Value::String(s)) = record.data().get(&tag) {
            let bytes: &[u8] = s.as_ref();
            assert_eq!(bytes, b"0123456");
        } else {
            panic!("Tag XS not found or wrong type");
        }
    }

    /// A non-UTF-8 string tag must be clipped byte-for-byte, not silently
    /// replaced with an empty string (or corrupted with U+FFFD). Regression test
    /// for the `from_utf8(..).unwrap_or("")` / `from_utf8_lossy` handling.
    #[test]
    fn test_auto_clip_attributes_string_non_utf8_preserves_bytes() {
        // A non-UTF-8 `Z` tag cannot be represented in a noodles `RecordBuf`, so this exercises
        // the raw clipper directly: encode a clean 10M read, append the non-UTF-8 XX:Z tag to the
        // raw bytes, then hard-clip.
        let buf = create_test_record("10M", "ACGTACGTAC", 1000);
        let mut rec = to_raw(&buf).as_ref().to_vec();
        // A 10-byte tag value that is NOT valid UTF-8 (0xFF/0xFE never appear in well-formed
        // UTF-8), matching the 10-base read length so it is clipped.
        let raw: Vec<u8> = vec![b'0', b'1', 0xFF, 0xFE, b'4', b'5', b'6', b'7', b'8', b'9'];
        {
            let mut ed = fgumi_raw_bam::RawTagsEditor::from_vec(&mut rec);
            ed.append_string(b"XX", &raw);
        }
        let mut record = fgumi_raw_bam::RawRecord::from(rec);

        // Clip 3 bases from the 5' end (== start of the alignment for a forward read).
        assert_eq!(
            RawRecordClipper::with_auto_clip(ClippingMode::Hard, true)
                .clip_start_of_alignment(&mut record, 3),
            3
        );

        // The tag must hold exactly bytes 3..10 of the original — not "" and not a
        // U+FFFD-mangled string.
        let bytes = fgumi_raw_bam::find_string_tag_in_record(record.as_ref(), b"XX")
            .expect("XX tag present");
        assert_eq!(bytes, &raw[3..], "non-UTF-8 tag bytes must be preserved");
    }

    #[test]
    fn test_auto_clip_attributes_array_5_prime() {
        use noodles::sam::alignment::record::data::field::Tag;
        use noodles::sam::alignment::record_buf::data::field::value::Array;

        let clipper = RawClipperOnBuf::with_auto_clip(ClippingMode::Hard, true);
        let mut record = create_test_record("10M", "ACGTACGTAC", 1000);

        // Add a UInt8Array attribute that matches the read length
        let tag = Tag::from([b'X', b'A']);
        let array: Vec<u8> = vec![0, 1, 2, 3, 4, 5, 6, 7, 8, 9];
        let value = Value::from(array);
        record.data_mut().insert(tag, value);

        // Clip 3 bases from 5' end
        let clipped = clipper.clip_start_of_alignment(&mut record, 3);
        assert_eq!(clipped, 3);

        // Check that the attribute was clipped
        if let Some(Value::Array(Array::UInt8(arr))) = record.data().get(&tag) {
            let vec_data: Vec<u8> = arr.clone();
            assert_eq!(vec_data, vec![3, 4, 5, 6, 7, 8, 9]);
        } else {
            panic!("Tag XA not found or wrong type");
        }
    }

    #[test]
    fn test_auto_clip_attributes_array_3_prime() {
        use noodles::sam::alignment::record::data::field::Tag;
        use noodles::sam::alignment::record_buf::data::field::value::Array;

        let clipper = RawClipperOnBuf::with_auto_clip(ClippingMode::Hard, true);
        let mut record = create_test_record("10M", "ACGTACGTAC", 1000);

        // Add a UInt8Array attribute that matches the read length
        let tag = Tag::from([b'X', b'A']);
        let array: Vec<u8> = vec![0, 1, 2, 3, 4, 5, 6, 7, 8, 9];
        let value = Value::from(array);
        record.data_mut().insert(tag, value);

        // Clip 3 bases from 3' end
        let clipped = clipper.clip_end_of_alignment(&mut record, 3);
        assert_eq!(clipped, 3);

        // Check that the attribute was clipped
        if let Some(Value::Array(Array::UInt8(arr))) = record.data().get(&tag) {
            let vec_data: Vec<u8> = arr.clone();
            assert_eq!(vec_data, vec![0, 1, 2, 3, 4, 5, 6]);
        } else {
            panic!("Tag XA not found or wrong type");
        }
    }

    #[test]
    fn test_auto_clip_attributes_only_in_hard_mode() {
        use noodles::sam::alignment::record::data::field::Tag;

        // Test with Soft mode - attributes should NOT be clipped
        let clipper = RawClipperOnBuf::with_auto_clip(ClippingMode::Soft, true);
        let mut record = create_test_record("10M", "ACGTACGTAC", 1000);

        let tag = Tag::from([b'X', b'S']);
        let value = Value::from("0123456789");
        record.data_mut().insert(tag, value);

        clipper.clip_start_of_alignment(&mut record, 3);

        // Attribute should remain unchanged in Soft mode
        if let Some(Value::String(s)) = record.data().get(&tag) {
            let bytes: &[u8] = s.as_ref();
            assert_eq!(bytes, b"0123456789");
        }
    }

    #[test]
    fn test_auto_clip_attributes_only_when_enabled() {
        use noodles::sam::alignment::record::data::field::Tag;

        // Test with auto_clip_attributes disabled - attributes should NOT be clipped
        let clipper = RawClipperOnBuf::with_auto_clip(ClippingMode::Hard, false);
        let mut record = create_test_record("10M", "ACGTACGTAC", 1000);

        let tag = Tag::from([b'X', b'S']);
        let value = Value::from("0123456789");
        record.data_mut().insert(tag, value);

        clipper.clip_start_of_alignment(&mut record, 3);

        // Attribute should remain unchanged when auto-clip is disabled
        if let Some(Value::String(s)) = record.data().get(&tag) {
            let bytes: &[u8] = s.as_ref();
            assert_eq!(bytes, b"0123456789");
        }
    }

    #[test]
    fn test_auto_clip_attributes_only_matching_length() {
        use noodles::sam::alignment::record::data::field::Tag;

        let clipper = RawClipperOnBuf::with_auto_clip(ClippingMode::Hard, true);
        let mut record = create_test_record("10M", "ACGTACGTAC", 1000);

        // Add two attributes: one matching length, one not
        let tag1 = Tag::from([b'X', b'1']);
        let value1 = Value::from("0123456789"); // Matches length (10)
        record.data_mut().insert(tag1, value1);

        let tag2 = Tag::from([b'X', b'2']);
        let value2 = Value::from("01234"); // Does not match length (5)
        record.data_mut().insert(tag2, value2);

        clipper.clip_start_of_alignment(&mut record, 3);

        // Check tag1 was clipped
        if let Some(Value::String(s)) = record.data().get(&tag1) {
            let bytes: &[u8] = s.as_ref();
            assert_eq!(bytes, b"3456789");
        }

        // Check tag2 was NOT clipped
        if let Some(Value::String(s)) = record.data().get(&tag2) {
            let bytes: &[u8] = s.as_ref();
            assert_eq!(bytes, b"01234");
        }
    }

    // ===================================================================
    // clip_3_prime tests - Additional coverage matching Scala
    // ===================================================================

    #[test]
    fn test_clip_end_of_alignment_soft_with_insertion() {
        let clipper = RawClipperOnBuf::new(ClippingMode::Soft);
        // 44M + 2I + 4M = 50 bases
        let seq = "ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTAC"; // 50 bases
        assert_eq!(seq.len(), 50);
        let mut record = create_test_record("44M2I4M", seq, 10);

        let clipped = clipper.clip_end_of_alignment(&mut record, 10);
        assert_eq!(clipped, 10);
        assert_eq!(record.alignment_start(), Position::new(10));

        let cigar_str = format_cigar(&record.cigar());
        assert_eq!(cigar_str, "40M10S");
    }

    #[test]
    fn test_clip_end_of_alignment_soft_with_deletion() {
        let clipper = RawClipperOnBuf::new(ClippingMode::Soft);
        // 44M + 2D + 6M = 50 bases in query
        let seq = "ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTAC"; // 50 bases
        assert_eq!(seq.len(), 50);
        let mut record = create_test_record("44M2D6M", seq, 10);

        let clipped = clipper.clip_end_of_alignment(&mut record, 10);
        assert_eq!(clipped, 10);
        assert_eq!(record.alignment_start(), Position::new(10));

        let cigar_str = format_cigar(&record.cigar());
        assert_eq!(cigar_str, "40M10S");
    }

    #[test]
    fn test_clip_end_of_alignment_soft_additional_bases() {
        let clipper = RawClipperOnBuf::new(ClippingMode::Soft);
        let seq = "ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTAC"; // 50 bases
        let mut record = create_test_record("40M10S", seq, 10);

        let clipped = clipper.clip_end_of_alignment(&mut record, 10);
        assert_eq!(clipped, 10);
        assert_eq!(record.alignment_start(), Position::new(10));

        let cigar_str = format_cigar(&record.cigar());
        assert_eq!(cigar_str, "30M20S");
    }

    #[test]
    fn test_clip_end_of_alignment_preserve_hard_clip() {
        let clipper = RawClipperOnBuf::new(ClippingMode::Soft);
        let seq = "ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTAC"; // 40 bases (no hard clip in sequence)
        let mut record = create_test_record("40M10H", &seq[..40], 10);

        let clipped = clipper.clip_end_of_alignment(&mut record, 10);
        assert_eq!(clipped, 10);
        assert_eq!(record.alignment_start(), Position::new(10));

        let cigar_str = format_cigar(&record.cigar());
        assert_eq!(cigar_str, "30M10S10H");
    }

    #[test]
    fn test_clip_end_of_alignment_complex_cigar() {
        let clipper = RawClipperOnBuf::new(ClippingMode::Soft);
        // 10M + 5I + 5M + 10I + 16M + 4S + 2H = 50 query bases
        let seq = "ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTAC"; // 50 bases
        let mut record = create_test_record("10M5I5M10I16M4S2H", seq, 10);

        let clipped = clipper.clip_end_of_alignment(&mut record, 10);
        assert_eq!(clipped, 10);
        assert_eq!(record.alignment_start(), Position::new(10));

        let cigar_str = format_cigar(&record.cigar());
        assert_eq!(cigar_str, "10M5I5M10I6M14S2H");
    }

    #[test]
    fn test_clip_end_of_alignment_consumes_leading_insertion() {
        let clipper = RawClipperOnBuf::new(ClippingMode::Soft);
        // 38M + 4I + 8M = 50 bases
        let seq = "ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTAC"; // 50 bases
        let mut record = create_test_record("38M4I8M", seq, 10);

        // Ask to clip 10 bases, but should consume 12 bases (4I + 8M) because
        // the insertion at the clip boundary is consumed entirely
        let clipped = clipper.clip_end_of_alignment(&mut record, 10);
        assert_eq!(clipped, 12);
        assert_eq!(record.alignment_start(), Position::new(10));

        let cigar_str = format_cigar(&record.cigar());
        assert_eq!(cigar_str, "38M12S");
    }

    #[test]
    fn test_clip_end_of_alignment_preserve_insertion_after_clip() {
        let clipper = RawClipperOnBuf::new(ClippingMode::Soft);
        // 36M + 4I + 10M = 50 bases
        let seq = "ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTAC"; // 50 bases
        let mut record = create_test_record("36M4I10M", seq, 10);

        let clipped = clipper.clip_end_of_alignment(&mut record, 10);
        assert_eq!(clipped, 10);
        assert_eq!(record.alignment_start(), Position::new(10));

        let cigar_str = format_cigar(&record.cigar());
        assert_eq!(cigar_str, "36M4I10S");
    }

    #[test]
    fn test_clip_end_of_alignment_remove_deletion_before_clip() {
        let clipper = RawClipperOnBuf::new(ClippingMode::Soft);
        // 40M + 4D + 10M = 50 query bases
        let seq = "ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTAC"; // 50 bases
        let mut record = create_test_record("40M4D10M", seq, 10);

        let clipped = clipper.clip_end_of_alignment(&mut record, 10);
        assert_eq!(clipped, 10);
        assert_eq!(record.alignment_start(), Position::new(10));

        let cigar_str = format_cigar(&record.cigar());
        assert_eq!(cigar_str, "40M10S");
    }

    #[test]
    fn test_clip_end_of_alignment_preserve_distant_deletion() {
        let clipper = RawClipperOnBuf::new(ClippingMode::Soft);
        // 25M + 4D + 25M = 50 query bases
        let seq = "ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTAC"; // 50 bases
        let mut record = create_test_record("25M4D25M", seq, 10);

        let clipped = clipper.clip_end_of_alignment(&mut record, 10);
        assert_eq!(clipped, 10);
        assert_eq!(record.alignment_start(), Position::new(10));

        let cigar_str = format_cigar(&record.cigar());
        assert_eq!(cigar_str, "25M4D15M10S");
    }

    #[test]
    fn test_clip_end_of_alignment_unmapped_read_does_nothing() {
        let clipper = RawClipperOnBuf::new(ClippingMode::Soft);
        let seq = "ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTAC"; // 50 bases
        let mut record = create_test_record("50M", seq, 10);

        // Set unmapped flag
        *record.flags_mut() = Flags::UNMAPPED;

        let clipped = clipper.clip_end_of_alignment(&mut record, 10);

        // Should not clip unmapped reads
        assert_eq!(clipped, 0);
        assert_eq!(record.alignment_start(), Position::new(10));

        // CIGAR should remain unchanged
        let cigar_str = format_cigar(&record.cigar());
        assert_eq!(cigar_str, "50M");
    }

    #[test]
    fn test_clip_end_of_alignment_soft_with_mask() {
        let clipper = RawClipperOnBuf::new(ClippingMode::SoftWithMask);
        let seq = "ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTAC"; // 50 bases
        let mut record = create_test_record("50M", seq, 10);

        let clipped = clipper.clip_end_of_alignment(&mut record, 10);
        assert_eq!(clipped, 10);
        assert_eq!(record.alignment_start(), Position::new(10));

        let cigar_str = format_cigar(&record.cigar());
        assert_eq!(cigar_str, "40M10S");

        // Check last 10 bases are masked to N
        let bases = record.sequence();
        for i in 40..50 {
            assert_eq!(bases.as_ref()[i], b'N', "Base at position {i} should be N");
        }

        // Check last 10 quals are set to min
        let quals = record.quality_scores();
        for i in 40..50 {
            assert_eq!(quals.as_ref()[i], 2, "Quality at position {i} should be 2");
        }
    }

    #[test]
    fn test_clip_end_of_alignment_soft_with_mask_existing_soft_clip() {
        let clipper = RawClipperOnBuf::new(ClippingMode::SoftWithMask);
        let seq = "ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTAC"; // 50 bases
        let mut record = create_test_record("40M10S", seq, 10);

        let clipped = clipper.clip_end_of_alignment(&mut record, 10);
        assert_eq!(clipped, 10);
        assert_eq!(record.alignment_start(), Position::new(10));

        let cigar_str = format_cigar(&record.cigar());
        assert_eq!(cigar_str, "30M20S");

        // Check last 20 bases are masked to N
        let bases = record.sequence();
        for i in 30..50 {
            assert_eq!(bases.as_ref()[i], b'N', "Base at position {i} should be N");
        }

        // Check last 20 quals are set to min
        let quals = record.quality_scores();
        for i in 30..50 {
            assert_eq!(quals.as_ref()[i], 2, "Quality at position {i} should be 2");
        }
    }

    #[test]
    fn test_clip_end_of_alignment_hard() {
        let clipper = RawClipperOnBuf::new(ClippingMode::Hard);
        let seq = "ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTAC"; // 50 bases
        let mut record = create_test_record("50M", seq, 10);

        let orig_seq: Vec<u8> = record.sequence().as_ref().to_vec();
        let orig_qual: Vec<u8> = record.quality_scores().as_ref().to_vec();

        let clipped = clipper.clip_end_of_alignment(&mut record, 10);
        assert_eq!(clipped, 10);
        assert_eq!(record.alignment_start(), Position::new(10));

        let cigar_str = format_cigar(&record.cigar());
        assert_eq!(cigar_str, "40M10H");

        // Check sequence and qualities are truncated
        assert_eq!(record.sequence().as_ref(), &orig_seq[..40]);
        assert_eq!(record.quality_scores().as_ref(), &orig_qual[..40]);
    }

    #[test]
    fn test_clip_start_of_alignment_hard_existing_soft_clip() {
        let clipper = RawClipperOnBuf::new(ClippingMode::Hard);
        let seq = "ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTAC"; // 50 bases
        assert_eq!(seq.len(), 50);
        let mut record = create_test_record("10S40M", seq, 10);

        let orig_seq: Vec<u8> = record.sequence().as_ref().to_vec();
        let orig_qual: Vec<u8> = record.quality_scores().as_ref().to_vec();

        let clipped = clipper.clip_start_of_alignment(&mut record, 10);
        assert_eq!(clipped, 10);
        assert_eq!(record.alignment_start(), Position::new(20));

        let cigar_str = format_cigar(&record.cigar());
        assert_eq!(cigar_str, "20H30M");

        // Check sequence and qualities are truncated (removed 10S + 10M = 20 bases)
        assert_eq!(record.sequence().as_ref(), &orig_seq[20..]);
        assert_eq!(record.quality_scores().as_ref(), &orig_qual[20..]);
    }

    // ===================================================================
    // Auto-trim attribute tests - all mode combinations
    // ===================================================================

    #[test]
    fn test_clip_start_of_alignment_auto_trim_soft_mode_false() {
        use noodles::sam::alignment::record::data::field::Tag;

        let clipper = RawClipperOnBuf::with_auto_clip(ClippingMode::Soft, false);
        let seq = "ACGTACGTACGTACGTACGT"; // 20 bases
        let mut record = create_test_record("20M", seq, 10);

        let a1_tag = Tag::from([b'A', b'1']);
        record.data_mut().insert(a1_tag, Value::from("AB".repeat(10)));

        let clipped = clipper.clip_start_of_alignment(&mut record, 5);
        assert_eq!(clipped, 5);

        // In Soft mode with auto=false, attributes should NOT be modified
        if let Some(Value::String(s)) = record.data().get(&a1_tag) {
            let bytes: &[u8] = s.as_ref();
            assert_eq!(bytes, "AB".repeat(10).as_bytes());
        }
    }

    #[test]
    fn test_clip_start_of_alignment_auto_trim_soft_mode_true() {
        use noodles::sam::alignment::record::data::field::Tag;

        let clipper = RawClipperOnBuf::with_auto_clip(ClippingMode::Soft, true);
        let seq = "ACGTACGTACGTACGTACGT"; // 20 bases
        let mut record = create_test_record("20M", seq, 10);

        let a1_tag = Tag::from([b'A', b'1']);
        record.data_mut().insert(a1_tag, Value::from("AB".repeat(10)));

        let clipped = clipper.clip_start_of_alignment(&mut record, 5);
        assert_eq!(clipped, 5);

        // In Soft mode, even with auto=true, attributes should NOT be modified
        if let Some(Value::String(s)) = record.data().get(&a1_tag) {
            let bytes: &[u8] = s.as_ref();
            assert_eq!(bytes, "AB".repeat(10).as_bytes());
        }
    }

    #[test]
    fn test_clip_start_of_alignment_auto_trim_soft_with_mask_mode_false() {
        use noodles::sam::alignment::record::data::field::Tag;

        let clipper = RawClipperOnBuf::with_auto_clip(ClippingMode::SoftWithMask, false);
        let seq = "ACGTACGTACGTACGTACGT"; // 20 bases
        let mut record = create_test_record("20M", seq, 10);

        let a1_tag = Tag::from([b'A', b'1']);
        record.data_mut().insert(a1_tag, Value::from("AB".repeat(10)));

        let clipped = clipper.clip_start_of_alignment(&mut record, 5);
        assert_eq!(clipped, 5);

        // In SoftWithMask mode with auto=false, attributes should NOT be modified
        if let Some(Value::String(s)) = record.data().get(&a1_tag) {
            let bytes: &[u8] = s.as_ref();
            assert_eq!(bytes, "AB".repeat(10).as_bytes());
        }
    }

    #[test]
    fn test_clip_start_of_alignment_auto_trim_soft_with_mask_mode_true() {
        use noodles::sam::alignment::record::data::field::Tag;

        let clipper = RawClipperOnBuf::with_auto_clip(ClippingMode::SoftWithMask, true);
        let seq = "ACGTACGTACGTACGTACGT"; // 20 bases
        let mut record = create_test_record("20M", seq, 10);

        let a1_tag = Tag::from([b'A', b'1']);
        record.data_mut().insert(a1_tag, Value::from("AB".repeat(10)));

        let clipped = clipper.clip_start_of_alignment(&mut record, 5);
        assert_eq!(clipped, 5);

        // In SoftWithMask mode, even with auto=true, attributes should NOT be modified
        if let Some(Value::String(s)) = record.data().get(&a1_tag) {
            let bytes: &[u8] = s.as_ref();
            assert_eq!(bytes, "AB".repeat(10).as_bytes());
        }
    }

    #[test]
    fn test_clip_start_of_alignment_auto_trim_hard_mode_false() {
        use noodles::sam::alignment::record::data::field::Tag;
        use noodles::sam::alignment::record_buf::data::field::value::Array;

        let clipper = RawClipperOnBuf::with_auto_clip(ClippingMode::Hard, false);
        let seq = "ACGTACGTACGTACGTACGT"; // 20 bases
        let mut record = create_test_record("20M", seq, 10);

        let a1_tag = Tag::from([b'A', b'1']);
        let a2_tag = Tag::from([b'A', b'2']);

        record.data_mut().insert(a1_tag, Value::from("AB".repeat(10)));
        record.data_mut().insert(a2_tag, Value::from((1..=20).collect::<Vec<i32>>()));

        let clipped = clipper.clip_start_of_alignment(&mut record, 5);
        assert_eq!(clipped, 5);

        // In Hard mode with auto=false, attributes should NOT be modified
        if let Some(Value::String(s)) = record.data().get(&a1_tag) {
            let bytes: &[u8] = s.as_ref();
            assert_eq!(bytes, "AB".repeat(10).as_bytes());
        }
        if let Some(Value::Array(Array::Int32(arr))) = record.data().get(&a2_tag) {
            let vec: Vec<i32> = arr.clone();
            assert_eq!(vec, (1..=20).collect::<Vec<i32>>());
        }
    }

    #[test]
    fn test_clip_start_of_alignment_auto_trim_hard_mode_true() {
        use noodles::sam::alignment::record::data::field::Tag;
        use noodles::sam::alignment::record_buf::data::field::value::Array;

        let clipper = RawClipperOnBuf::with_auto_clip(ClippingMode::Hard, true);
        let seq = "ACGTACGTACGTACGTACGT"; // 20 bases
        let mut record = create_test_record("20M", seq, 10);

        let a1_tag = Tag::from([b'A', b'1']);
        let a2_tag = Tag::from([b'A', b'2']);
        let b1_tag = Tag::from([b'B', b'1']);
        let b2_tag = Tag::from([b'B', b'2']);

        record.data_mut().insert(a1_tag, Value::from("AB".repeat(10)));
        record.data_mut().insert(a2_tag, Value::from((1..=20).collect::<Vec<i32>>()));
        record.data_mut().insert(b1_tag, Value::from("A".repeat(10)));
        record.data_mut().insert(b2_tag, Value::from((1..=10).collect::<Vec<i32>>()));

        let clipped = clipper.clip_start_of_alignment(&mut record, 5);
        assert_eq!(clipped, 5);

        // In Hard mode with auto=true, attributes matching read length should be clipped
        if let Some(Value::String(s)) = record.data().get(&a1_tag) {
            // "ABABABABABABABABABAB" -> remove first 5 -> "BABABABABABABAB"
            let bytes: &[u8] = s.as_ref();
            assert_eq!(bytes, "BABABABABABABAB".as_bytes());
        }
        if let Some(Value::Array(Array::Int32(arr))) = record.data().get(&a2_tag) {
            let vec: Vec<i32> = arr.clone();
            assert_eq!(vec, (6..=20).collect::<Vec<i32>>());
        }
        // B1 and B2 should NOT be modified (length doesn't match)
        if let Some(Value::String(s)) = record.data().get(&b1_tag) {
            let bytes: &[u8] = s.as_ref();
            assert_eq!(bytes, "A".repeat(10).as_bytes());
        }
        if let Some(Value::Array(Array::Int32(arr))) = record.data().get(&b2_tag) {
            let vec: Vec<i32> = arr.clone();
            assert_eq!(vec, (1..=10).collect::<Vec<i32>>());
        }
    }

    #[test]
    fn test_clip_end_of_alignment_auto_trim_soft_mode_false() {
        use noodles::sam::alignment::record::data::field::Tag;

        let clipper = RawClipperOnBuf::with_auto_clip(ClippingMode::Soft, false);
        let seq = "ACGTACGTACGTACGTACGT"; // 20 bases
        let mut record = create_test_record("20M", seq, 10);

        let a1_tag = Tag::from([b'A', b'1']);
        record.data_mut().insert(a1_tag, Value::from("AB".repeat(10)));

        let clipped = clipper.clip_end_of_alignment(&mut record, 5);
        assert_eq!(clipped, 5);

        // In Soft mode with auto=false, attributes should NOT be modified
        if let Some(Value::String(s)) = record.data().get(&a1_tag) {
            let bytes: &[u8] = s.as_ref();
            assert_eq!(bytes, "AB".repeat(10).as_bytes());
        }
    }

    #[test]
    fn test_clip_end_of_alignment_auto_trim_soft_mode_true() {
        use noodles::sam::alignment::record::data::field::Tag;

        let clipper = RawClipperOnBuf::with_auto_clip(ClippingMode::Soft, true);
        let seq = "ACGTACGTACGTACGTACGTACGT"; // 20 bases
        let mut record = create_test_record("20M", seq, 10);

        let a1_tag = Tag::from([b'A', b'1']);
        record.data_mut().insert(a1_tag, Value::from("AB".repeat(10)));

        let clipped = clipper.clip_end_of_alignment(&mut record, 5);
        assert_eq!(clipped, 5);

        // In Soft mode, even with auto=true, attributes should NOT be modified
        if let Some(Value::String(s)) = record.data().get(&a1_tag) {
            let bytes: &[u8] = s.as_ref();
            assert_eq!(bytes, "AB".repeat(10).as_bytes());
        }
    }

    #[test]
    fn test_clip_end_of_alignment_auto_trim_soft_with_mask_mode_false() {
        use noodles::sam::alignment::record::data::field::Tag;

        let clipper = RawClipperOnBuf::with_auto_clip(ClippingMode::SoftWithMask, false);
        let seq = "ACGTACGTACGTACGTACGT"; // 20 bases
        let mut record = create_test_record("20M", seq, 10);

        let a1_tag = Tag::from([b'A', b'1']);
        record.data_mut().insert(a1_tag, Value::from("AB".repeat(10)));

        let clipped = clipper.clip_end_of_alignment(&mut record, 5);
        assert_eq!(clipped, 5);

        // In SoftWithMask mode with auto=false, attributes should NOT be modified
        if let Some(Value::String(s)) = record.data().get(&a1_tag) {
            let bytes: &[u8] = s.as_ref();
            assert_eq!(bytes, "AB".repeat(10).as_bytes());
        }
    }

    #[test]
    fn test_clip_end_of_alignment_auto_trim_soft_with_mask_mode_true() {
        use noodles::sam::alignment::record::data::field::Tag;

        let clipper = RawClipperOnBuf::with_auto_clip(ClippingMode::SoftWithMask, true);
        let seq = "ACGTACGTACGTACGTACGT"; // 20 bases
        let mut record = create_test_record("20M", seq, 10);

        let a1_tag = Tag::from([b'A', b'1']);
        record.data_mut().insert(a1_tag, Value::from("AB".repeat(10)));

        let clipped = clipper.clip_end_of_alignment(&mut record, 5);
        assert_eq!(clipped, 5);

        // In SoftWithMask mode, even with auto=true, attributes should NOT be modified
        if let Some(Value::String(s)) = record.data().get(&a1_tag) {
            let bytes: &[u8] = s.as_ref();
            assert_eq!(bytes, "AB".repeat(10).as_bytes());
        }
    }

    #[test]
    fn test_clip_end_of_alignment_auto_trim_hard_mode_false() {
        use noodles::sam::alignment::record::data::field::Tag;
        use noodles::sam::alignment::record_buf::data::field::value::Array;

        let clipper = RawClipperOnBuf::with_auto_clip(ClippingMode::Hard, false);
        let seq = "ACGTACGTACGTACGTACGT"; // 20 bases
        let mut record = create_test_record("20M", seq, 10);

        let a1_tag = Tag::from([b'A', b'1']);
        let a2_tag = Tag::from([b'A', b'2']);

        record.data_mut().insert(a1_tag, Value::from("AB".repeat(10)));
        record.data_mut().insert(a2_tag, Value::from((1..=20).collect::<Vec<i32>>()));

        let clipped = clipper.clip_end_of_alignment(&mut record, 5);
        assert_eq!(clipped, 5);

        // In Hard mode with auto=false, attributes should NOT be modified
        if let Some(Value::String(s)) = record.data().get(&a1_tag) {
            let bytes: &[u8] = s.as_ref();
            assert_eq!(bytes, "AB".repeat(10).as_bytes());
        }
        if let Some(Value::Array(Array::Int32(arr))) = record.data().get(&a2_tag) {
            let vec: Vec<i32> = arr.clone();
            assert_eq!(vec, (1..=20).collect::<Vec<i32>>());
        }
    }

    #[test]
    fn test_clip_end_of_alignment_auto_trim_hard_mode_true() {
        use noodles::sam::alignment::record::data::field::Tag;
        use noodles::sam::alignment::record_buf::data::field::value::Array;

        let clipper = RawClipperOnBuf::with_auto_clip(ClippingMode::Hard, true);
        let seq = "ACGTACGTACGTACGTACGT"; // 20 bases
        let mut record = create_test_record("20M", seq, 10);

        let a1_tag = Tag::from([b'A', b'1']);
        let a2_tag = Tag::from([b'A', b'2']);
        let b1_tag = Tag::from([b'B', b'1']);
        let b2_tag = Tag::from([b'B', b'2']);

        record.data_mut().insert(a1_tag, Value::from("AB".repeat(10)));
        record.data_mut().insert(a2_tag, Value::from((1..=20).collect::<Vec<i32>>()));
        record.data_mut().insert(b1_tag, Value::from("A".repeat(10)));
        record.data_mut().insert(b2_tag, Value::from((1..=10).collect::<Vec<i32>>()));

        let clipped = clipper.clip_end_of_alignment(&mut record, 5);
        assert_eq!(clipped, 5);

        // In Hard mode with auto=true, attributes matching read length should be clipped
        if let Some(Value::String(s)) = record.data().get(&a1_tag) {
            // "ABABABABABABABABABAB" -> remove last 5 -> "ABABABABABABABA"
            let bytes: &[u8] = s.as_ref();
            assert_eq!(bytes, "ABABABABABABABA".as_bytes());
        }
        if let Some(Value::Array(Array::Int32(arr))) = record.data().get(&a2_tag) {
            let vec: Vec<i32> = arr.clone();
            assert_eq!(vec, (1..=15).collect::<Vec<i32>>());
        }
        // B1 and B2 should NOT be modified (length doesn't match)
        if let Some(Value::String(s)) = record.data().get(&b1_tag) {
            let bytes: &[u8] = s.as_ref();
            assert_eq!(bytes, "A".repeat(10).as_bytes());
        }
        if let Some(Value::Array(Array::Int32(arr))) = record.data().get(&b2_tag) {
            let vec: Vec<i32> = arr.clone();
            assert_eq!(vec, (1..=10).collect::<Vec<i32>>());
        }
    }

    // ===================================================================
    // upgrade_clipping tests
    // ===================================================================

    #[test]
    fn test_upgrade_all_clipping_convert_leading_trailing() {
        use noodles::sam::alignment::record::data::field::Tag;

        // Test without auto-clip
        let clipper = RawClipperOnBuf::new(ClippingMode::Hard);
        let seq = "12345678901234567890123456789012345678901234567890"; // 50 bases
        let mut no_auto = create_test_record("5S35M10S", seq, 10);
        let az_tag = Tag::from([b'a', b'z']);
        no_auto
            .data_mut()
            .insert(az_tag, Value::from("12345678901234567890123456789012345678901234567890"));

        let result = clipper
            .upgrade_all_clipping(&mut no_auto)
            .expect("upgrade_all_clipping should succeed");
        assert_eq!(result, (5, 10));

        let cigar_str = format_cigar(&no_auto.cigar());
        assert_eq!(cigar_str, "5H35M10H");
        assert_eq!(no_auto.sequence().len(), 35);

        // Attributes should NOT be modified without auto-clip
        if let Some(Value::String(s)) = no_auto.data().get(&az_tag) {
            let bytes: &[u8] = s.as_ref();
            assert_eq!(bytes, "12345678901234567890123456789012345678901234567890".as_bytes());
        }

        // Test with auto-clip
        let clipper_auto = RawClipperOnBuf::with_auto_clip(ClippingMode::Hard, true);
        let mut with_auto = create_test_record("5S35M10S", seq, 10);
        with_auto
            .data_mut()
            .insert(az_tag, Value::from("12345678901234567890123456789012345678901234567890"));

        let result = clipper_auto
            .upgrade_all_clipping(&mut with_auto)
            .expect("upgrade_all_clipping should succeed");
        assert_eq!(result, (5, 10));

        let cigar_str = format_cigar(&with_auto.cigar());
        assert_eq!(cigar_str, "5H35M10H");
        assert_eq!(with_auto.sequence().len(), 35);

        // Attributes SHOULD be modified with auto-clip (remove first 5 and last 10)
        if let Some(Value::String(s)) = with_auto.data().get(&az_tag) {
            let bytes: &[u8] = s.as_ref();
            assert_eq!(bytes, "67890123456789012345678901234567890".as_bytes());
        }
    }

    /// The soft→hard upgrade path (`upgrade_all_clipping`) must also clip a
    /// non-UTF-8 string tag byte-for-byte, not mangle it with `from_utf8_lossy`.
    #[test]
    fn test_upgrade_all_clipping_string_non_utf8_preserves_bytes() {
        // A non-UTF-8 `Z` tag cannot round-trip through a noodles `RecordBuf`, so exercise the
        // raw clipper directly: encode a clean 5S35M10S read, append the non-UTF-8 aa:Z tag to
        // the raw bytes, then upgrade all clipping to hard with auto-clip on.
        let seq = "12345678901234567890123456789012345678901234567890"; // 50 bases
        let buf = create_test_record("5S35M10S", seq, 10);
        let mut rec = to_raw(&buf).as_ref().to_vec();
        // A 50-byte tag value that is not valid UTF-8, matching the read length.
        let mut raw: Vec<u8> = (0..50u8).map(|i| b'0' + (i % 10)).collect();
        raw[7] = 0xFF;
        raw[42] = 0xFE;
        {
            let mut ed = fgumi_raw_bam::RawTagsEditor::from_vec(&mut rec);
            ed.append_string(b"aa", &raw);
        }
        let mut record = fgumi_raw_bam::RawRecord::from(rec);

        let result = RawRecordClipper::with_auto_clip(ClippingMode::Hard, true)
            .upgrade_all_clipping_raw(&mut record)
            .expect("upgrade should succeed");
        assert_eq!(result, (5, 10));

        // Auto-clip removes the first 5 and last 10 bytes → raw[5..40], verbatim.
        let bytes = fgumi_raw_bam::find_string_tag_in_record(record.as_ref(), b"aa")
            .expect("aa tag present");
        assert_eq!(bytes, &raw[5..40], "non-UTF-8 tag bytes must be preserved on upgrade");
    }

    #[test]
    fn test_upgrade_all_clipping_soft_clips_after_hard_clips() {
        use noodles::sam::alignment::record::data::field::Tag;

        // Test without auto-clip
        let clipper = RawClipperOnBuf::new(ClippingMode::Hard);
        let seq = "12345678901234567890123456789012345678901234567890"; // 50 bases
        let mut no_auto = create_test_record("5H5S35M10S5H", seq, 10);
        let az_tag = Tag::from([b'a', b'z']);
        no_auto
            .data_mut()
            .insert(az_tag, Value::from("12345678901234567890123456789012345678901234567890"));

        let result = clipper
            .upgrade_all_clipping(&mut no_auto)
            .expect("upgrade_all_clipping should succeed");
        assert_eq!(result, (5, 10));

        let cigar_str = format_cigar(&no_auto.cigar());
        assert_eq!(cigar_str, "10H35M15H");
        assert_eq!(no_auto.sequence().len(), 35);

        // Attributes should NOT be modified without auto-clip
        if let Some(Value::String(s)) = no_auto.data().get(&az_tag) {
            let bytes: &[u8] = s.as_ref();
            assert_eq!(bytes, "12345678901234567890123456789012345678901234567890".as_bytes());
        }

        // Test with auto-clip
        let clipper_auto = RawClipperOnBuf::with_auto_clip(ClippingMode::Hard, true);
        let mut with_auto = create_test_record("5H5S35M10S5H", seq, 10);
        with_auto
            .data_mut()
            .insert(az_tag, Value::from("12345678901234567890123456789012345678901234567890"));

        let result = clipper_auto
            .upgrade_all_clipping(&mut with_auto)
            .expect("upgrade_all_clipping should succeed");
        assert_eq!(result, (5, 10));

        let cigar_str = format_cigar(&with_auto.cigar());
        assert_eq!(cigar_str, "10H35M15H");
        assert_eq!(with_auto.sequence().len(), 35);

        // Attributes SHOULD be modified with auto-clip
        if let Some(Value::String(s)) = with_auto.data().get(&az_tag) {
            let bytes: &[u8] = s.as_ref();
            assert_eq!(bytes, "67890123456789012345678901234567890".as_bytes());
        }
    }

    #[test]
    fn test_upgrade_all_clipping_no_soft_clipping() {
        use noodles::sam::alignment::record::data::field::Tag;

        let clipper = RawClipperOnBuf::new(ClippingMode::Hard);
        let az_tag = Tag::from([b'a', b'z']);

        // Test with no soft clips
        let seq1 = "1234567890123456789012345678901234567890123456789012345"; // 55 bases
        let mut no_soft = create_test_record("55M", seq1, 10);
        no_soft
            .data_mut()
            .insert(az_tag, Value::from("12345678901234567890123456789012345678901234567890"));

        let result = clipper
            .upgrade_all_clipping(&mut no_soft)
            .expect("upgrade_all_clipping should succeed");
        assert_eq!(result, (0, 0));

        let cigar_str = format_cigar(&no_soft.cigar());
        assert_eq!(cigar_str, "55M");
        assert_eq!(no_soft.sequence().len(), 55);

        // Test with only hard clips (no soft clips)
        let mut hard_only = create_test_record("5H55M10H", seq1, 10);
        hard_only
            .data_mut()
            .insert(az_tag, Value::from("12345678901234567890123456789012345678901234567890"));

        let result = clipper
            .upgrade_all_clipping(&mut hard_only)
            .expect("upgrade_all_clipping should succeed");
        assert_eq!(result, (0, 0));

        let cigar_str = format_cigar(&hard_only.cigar());
        assert_eq!(cigar_str, "5H55M10H");
        assert_eq!(hard_only.sequence().len(), 55);
    }

    #[test]
    fn test_upgrade_all_clipping_unmapped_or_wrong_mode() {
        use noodles::sam::alignment::record::data::field::Tag;

        let az_tag = Tag::from([b'a', b'z']);
        let seq = "1234567890123456789012345678901234567890123456789012345"; // 55 bases

        // Test Soft mode (should not convert)
        let clipper_soft = RawClipperOnBuf::new(ClippingMode::Soft);
        let mut mapped = create_test_record("55M", seq, 10);
        mapped
            .data_mut()
            .insert(az_tag, Value::from("12345678901234567890123456789012345678901234567890"));

        let result = clipper_soft
            .upgrade_all_clipping(&mut mapped)
            .expect("upgrade_all_clipping should succeed");
        assert_eq!(result, (0, 0));

        // Test SoftWithMask mode (should not convert)
        let clipper_mask = RawClipperOnBuf::new(ClippingMode::SoftWithMask);
        let mut mapped2 = create_test_record("55M", seq, 10);
        mapped2
            .data_mut()
            .insert(az_tag, Value::from("12345678901234567890123456789012345678901234567890"));

        let result = clipper_mask
            .upgrade_all_clipping(&mut mapped2)
            .expect("upgrade_all_clipping should succeed");
        assert_eq!(result, (0, 0));

        // Test unmapped read (should not convert)
        let clipper_hard = RawClipperOnBuf::new(ClippingMode::Hard);
        let mut unmapped = create_test_record("55M", seq, 10);
        *unmapped.flags_mut() = Flags::UNMAPPED;
        unmapped
            .data_mut()
            .insert(az_tag, Value::from("12345678901234567890123456789012345678901234567890"));

        let result = clipper_hard
            .upgrade_all_clipping(&mut unmapped)
            .expect("upgrade_all_clipping should succeed");
        assert_eq!(result, (0, 0));
    }

    // ===================================================================
    // clipOverlappingReads tests - complex overlap scenarios
    // ===================================================================

    fn create_pair(
        r1_start: usize,
        r1_cigar: &str,
        r1_seq: &str,
        r2_start: usize,
        r2_cigar: &str,
        r2_seq: &str,
    ) -> (RecordBuf, RecordBuf) {
        let mut r1 = create_test_record(r1_cigar, r1_seq, r1_start);
        let mut r2 = create_test_record(r2_cigar, r2_seq, r2_start);

        // Calculate template length (insert size)
        // For FR pair: template length = r2_end - r1_start + 1
        let r2_len = cigar_utils::reference_length(&r2.cigar());
        let r2_end = r2_start + r2_len - 1;
        let tlen = i32::try_from(r2_end).expect("r2_end should fit in i32")
            - i32::try_from(r1_start).expect("r1_start should fit in i32")
            + 1;

        // Set up FR pair flags
        // R1 is forward strand, R2 is reverse strand (typical FR pair)
        *r1.flags_mut() = Flags::SEGMENTED | Flags::MATE_REVERSE_COMPLEMENTED;
        *r2.flags_mut() = Flags::SEGMENTED | Flags::REVERSE_COMPLEMENTED;

        // Set mate information
        *r1.mate_reference_sequence_id_mut() = Some(0);
        *r1.mate_alignment_start_mut() = Position::new(r2_start);
        *r1.template_length_mut() = tlen;
        *r2.mate_reference_sequence_id_mut() = Some(0);
        *r2.mate_alignment_start_mut() = Position::new(r1_start);
        *r2.template_length_mut() = -tlen;

        (r1, r2)
    }

    /// CLIP3-01 (typed): `clip_overlapping_reads` must normalize by *strand*, not by
    /// argument order. Passing the negative-strand read first must produce the same
    /// per-record clipping as passing the positive-strand read first, with the
    /// returned tuple swapped. Mirrors fgbio `SamRecordClipper.clipOverlappingReads`
    /// (`if (rec.negativeStrand) clipOverlappingReads(rec=mate, mate=rec).swap`).
    #[rstest]
    #[case::one_base_overlap(1, 100)]
    #[case::large_overlap(1, 50)]
    #[case::most_overlap(1, 20)]
    fn test_clip_overlapping_reads_normalizes_by_strand_typed(
        #[case] fwd_start: usize,
        #[case] rev_start: usize,
        #[values(ClippingMode::Soft, ClippingMode::Hard)] mode: ClippingMode,
    ) {
        let clipper = RawClipperOnBuf::new(mode);
        let seq = "A".repeat(100);

        // Forward-first: positive-strand read passed first (canonical, already handled).
        let mut fwd_a =
            create_paired_record("100M", &seq, fwd_start, false, true, rev_start, "100M");
        let mut rev_a =
            create_paired_record("100M", &seq, rev_start, true, false, fwd_start, "100M");
        let (clip_fwd_a, clip_rev_a) = clipper.clip_overlapping_reads(&mut fwd_a, &mut rev_a);

        // Reverse-first: same biology, negative-strand read passed first.
        let mut rev_b =
            create_paired_record("100M", &seq, rev_start, true, false, fwd_start, "100M");
        let mut fwd_b =
            create_paired_record("100M", &seq, fwd_start, false, true, rev_start, "100M");
        let (clip_rev_b, clip_fwd_b) = clipper.clip_overlapping_reads(&mut rev_b, &mut fwd_b);

        // Same per-record clip counts regardless of argument order.
        assert_eq!(clip_fwd_a, clip_fwd_b, "forward-read clip count differs by arg order");
        assert_eq!(clip_rev_a, clip_rev_b, "reverse-read clip count differs by arg order");

        // Same resulting record state (CIGAR + alignment start) for each read.
        assert_eq!(format_cigar(&fwd_a.cigar()), format_cigar(&fwd_b.cigar()), "forward CIGAR");
        assert_eq!(fwd_a.alignment_start(), fwd_b.alignment_start(), "forward start");
        assert_eq!(format_cigar(&rev_a.cigar()), format_cigar(&rev_b.cigar()), "reverse CIGAR");
        assert_eq!(rev_a.alignment_start(), rev_b.alignment_start(), "reverse start");
    }

    /// CLIP3-01 (raw): same strand-normalization requirement for the raw-byte
    /// `clip_overlapping_reads` used by the `clip` command's hot path.
    #[rstest]
    #[case::one_base_overlap(1, 100)]
    #[case::large_overlap(1, 50)]
    #[case::most_overlap(1, 20)]
    fn test_clip_overlapping_reads_normalizes_by_strand_raw(
        #[case] fwd_start: usize,
        #[case] rev_start: usize,
        #[values(ClippingMode::Soft, ClippingMode::Hard)] mode: ClippingMode,
    ) {
        use fgumi_raw_bam::encode_record_buf_to_raw;
        use noodles::sam::header::record::value::Map;
        use noodles::sam::header::record::value::map::ReferenceSequence;
        use std::num::NonZeroUsize;

        let seq = "A".repeat(100);
        let ref_seq =
            Map::<ReferenceSequence>::new(NonZeroUsize::new(100_000).expect("ref length nonzero"));
        let header =
            noodles::sam::Header::builder().add_reference_sequence(b"chr1", ref_seq).build();
        let clipper = RawRecordClipper::new(mode);
        let enc = |buf: &RecordBuf| encode_record_buf_to_raw(buf, &header).expect("encode raw");

        // Forward-first.
        let mut fwd_a =
            enc(&create_paired_record("100M", &seq, fwd_start, false, true, rev_start, "100M"));
        let mut rev_a =
            enc(&create_paired_record("100M", &seq, rev_start, true, false, fwd_start, "100M"));
        let (clip_fwd_a, clip_rev_a) = clipper.clip_overlapping_reads(&mut fwd_a, &mut rev_a);

        // Reverse-first.
        let mut rev_b =
            enc(&create_paired_record("100M", &seq, rev_start, true, false, fwd_start, "100M"));
        let mut fwd_b =
            enc(&create_paired_record("100M", &seq, fwd_start, false, true, rev_start, "100M"));
        let (clip_rev_b, clip_fwd_b) = clipper.clip_overlapping_reads(&mut rev_b, &mut fwd_b);

        assert_eq!(clip_fwd_a, clip_fwd_b, "forward-read clip count differs by arg order");
        assert_eq!(clip_rev_a, clip_rev_b, "reverse-read clip count differs by arg order");
        assert_eq!(fwd_a.cigar_ops_vec(), fwd_b.cigar_ops_vec(), "forward CIGAR");
        assert_eq!(fwd_a.alignment_start_1based(), fwd_b.alignment_start_1based(), "forward start");
        assert_eq!(rev_a.cigar_ops_vec(), rev_b.cigar_ops_vec(), "reverse CIGAR");
        assert_eq!(rev_a.alignment_start_1based(), rev_b.alignment_start_1based(), "reverse start");
    }

    #[test]
    fn test_clip_overlapping_reads_one_base_overlap() {
        let clipper = RawClipperOnBuf::new(ClippingMode::Soft);
        let seq = "A".repeat(100);
        let (mut r1, mut r2) = create_pair(1, "100M", &seq, 100, "100M", &seq);

        let (clipped_r1, clipped_r2) = clipper.clip_overlapping_reads(&mut r1, &mut r2);
        assert_eq!(clipped_r1, 0);
        assert_eq!(clipped_r2, 1);

        assert_eq!(r1.alignment_start(), Position::new(1));
        let cigar_r1 = format_cigar(&r1.cigar());
        assert_eq!(cigar_r1, "100M");

        assert_eq!(r2.alignment_start(), Position::new(101));
        let cigar_r2 = format_cigar(&r2.cigar());
        assert_eq!(cigar_r2, "1S99M");
    }

    #[test]
    fn test_clip_overlapping_reads_two_base_overlap() {
        let clipper = RawClipperOnBuf::new(ClippingMode::Soft);
        let seq = "A".repeat(100);
        let (mut r1, mut r2) = create_pair(2, "100M", &seq, 100, "100M", &seq);

        let (clipped_r1, clipped_r2) = clipper.clip_overlapping_reads(&mut r1, &mut r2);
        assert_eq!(clipped_r1, 1);
        assert_eq!(clipped_r2, 1);

        assert_eq!(r1.alignment_start(), Position::new(2));
        let cigar_r1 = format_cigar(&r1.cigar());
        assert_eq!(cigar_r1, "99M1S");

        assert_eq!(r2.alignment_start(), Position::new(101));
        let cigar_r2 = format_cigar(&r2.cigar());
        assert_eq!(cigar_r2, "1S99M");
    }

    #[test]
    fn test_clip_overlapping_reads_with_deletion() {
        let clipper = RawClipperOnBuf::new(ClippingMode::Soft);
        let seq = "A".repeat(90);
        let (mut r1, mut r2) = create_pair(2, "80M10D10M", &seq, 70, "10M10D80M", &seq);

        let (clipped_r1, clipped_r2) = clipper.clip_overlapping_reads(&mut r1, &mut r2);
        assert_eq!(clipped_r1, 10);
        assert_eq!(clipped_r2, 10);

        assert_eq!(r1.alignment_start(), Position::new(2));
        let cigar_r1 = format_cigar(&r1.cigar());
        assert_eq!(cigar_r1, "80M10S");

        assert_eq!(r2.alignment_start(), Position::new(90));
        let cigar_r2 = format_cigar(&r2.cigar());
        assert_eq!(cigar_r2, "10S80M");
    }

    #[test]
    fn test_clip_overlapping_reads_full_overlap() {
        let clipper = RawClipperOnBuf::new(ClippingMode::Soft);
        let seq = "A".repeat(100);
        let (mut r1, mut r2) = create_pair(1, "100M", &seq, 1, "100M", &seq);

        let (clipped_r1, clipped_r2) = clipper.clip_overlapping_reads(&mut r1, &mut r2);
        assert_eq!(clipped_r1, 50);
        assert_eq!(clipped_r2, 50);

        assert_eq!(r1.alignment_start(), Position::new(1));
        let cigar_r1 = format_cigar(&r1.cigar());
        assert_eq!(cigar_r1, "50M50S");

        assert_eq!(r2.alignment_start(), Position::new(51));
        let cigar_r2 = format_cigar(&r2.cigar());
        assert_eq!(cigar_r2, "50S50M");
    }

    #[test]
    fn test_clip_overlapping_reads_extend_past_each_other() {
        let clipper = RawClipperOnBuf::new(ClippingMode::Soft);
        let seq = "A".repeat(100);
        let (mut r1, mut r2) = create_pair(50, "100M", &seq, 1, "100M", &seq);

        // Midpoint should be (50+100)/2 = 75
        let (clipped_r1, clipped_r2) = clipper.clip_overlapping_reads(&mut r1, &mut r2);
        assert_eq!(clipped_r1, 74);
        assert_eq!(clipped_r2, 75);

        assert_eq!(r1.alignment_start(), Position::new(50));
        let cigar_r1 = format_cigar(&r1.cigar());
        assert_eq!(cigar_r1, "26M74S");

        assert_eq!(r2.alignment_start(), Position::new(76));
        let cigar_r2 = format_cigar(&r2.cigar());
        assert_eq!(cigar_r2, "75S25M");
    }

    #[test]
    fn test_clip_overlapping_reads_forward_much_longer() {
        let clipper = RawClipperOnBuf::with_auto_clip(ClippingMode::Hard, false);
        let r1_seq = "A".repeat(100);
        let r2_seq = "A".repeat(100);
        let (mut r1, mut r2) = create_pair(1, "100M", &r1_seq, 30, "80S20M", &r2_seq);

        let (clipped_r1, clipped_r2) = clipper.clip_overlapping_reads(&mut r1, &mut r2);
        assert_eq!(clipped_r1, 71);
        assert_eq!(clipped_r2, 0);

        assert_eq!(r1.alignment_start(), Position::new(1));
        let cigar_r1 = format_cigar(&r1.cigar());
        assert_eq!(cigar_r1, "29M71H");

        assert_eq!(r2.alignment_start(), Position::new(30));
        let cigar_r2 = format_cigar(&r2.cigar());
        assert_eq!(cigar_r2, "80H20M");
    }

    #[test]
    fn test_clip_overlapping_reads_reverse_much_longer() {
        let clipper = RawClipperOnBuf::with_auto_clip(ClippingMode::Hard, false);
        let r1_seq = "A".repeat(100);
        let r2_seq = "A".repeat(100);
        let (mut r1, mut r2) = create_pair(50, "20M80S", &r1_seq, 1, "100M", &r2_seq);

        let (clipped_r1, clipped_r2) = clipper.clip_overlapping_reads(&mut r1, &mut r2);
        assert_eq!(clipped_r1, 0);
        assert_eq!(clipped_r2, 69);

        assert_eq!(r1.alignment_start(), Position::new(50));
        let cigar_r1 = format_cigar(&r1.cigar());
        assert_eq!(cigar_r1, "20M80H");

        assert_eq!(r2.alignment_start(), Position::new(70));
        let cigar_r2 = format_cigar(&r2.cigar());
        assert_eq!(cigar_r2, "69H31M");
    }

    #[test]
    fn test_clip_overlapping_reads_with_one_end_deletion() {
        let clipper = RawClipperOnBuf::new(ClippingMode::Soft);
        let r1_seq = "A".repeat(100);
        let r2_seq = "A".repeat(110);
        let (mut r1, mut r2) = create_pair(1, "60M10D40M", &r1_seq, 50, "10M10D80M10D10M", &r2_seq);

        let (clipped_r1, clipped_r2) = clipper.clip_overlapping_reads(&mut r1, &mut r2);
        assert_eq!(clipped_r1, 25);
        assert_eq!(clipped_r2, 26);

        assert_eq!(r1.alignment_start(), Position::new(1));
        let cigar_r1 = format_cigar(&r1.cigar());
        assert_eq!(cigar_r1, "60M10D15M25S");

        assert_eq!(r2.alignment_start(), Position::new(86));
        let cigar_r2 = format_cigar(&r2.cigar());
        assert_eq!(cigar_r2, "26S64M10D10M");
    }

    #[test]
    fn test_clip_overlapping_reads_both_ends_deletion() {
        let clipper = RawClipperOnBuf::new(ClippingMode::Soft);
        let seq = "A".repeat(100);
        let (mut r1, mut r2) = create_pair(1, "50M10D50M", &seq, 3, "47M10D53M", &seq);

        let (clipped_r1, clipped_r2) = clipper.clip_overlapping_reads(&mut r1, &mut r2);
        assert_eq!(clipped_r1, 50);
        assert_eq!(clipped_r2, 47);

        assert_eq!(r1.alignment_start(), Position::new(1));
        let cigar_r1 = format_cigar(&r1.cigar());
        assert_eq!(cigar_r1, "50M50S");

        // R2 clips through the 47M and 10D, so the new start is at position 60 (where 53M begins)
        assert_eq!(r2.alignment_start(), Position::new(60));
        let cigar_r2 = format_cigar(&r2.cigar());
        assert_eq!(cigar_r2, "47S53M");
    }

    #[test]
    fn test_clip_overlapping_reads_no_overlap_far_apart() {
        let clipper = RawClipperOnBuf::new(ClippingMode::Soft);
        let seq = "A".repeat(100);
        let (mut r1, mut r2) = create_pair(1000, "100M", &seq, 1, "100M", &seq);

        let (clipped_r1, clipped_r2) = clipper.clip_overlapping_reads(&mut r1, &mut r2);
        assert_eq!(clipped_r1, 0);
        assert_eq!(clipped_r2, 0);

        assert_eq!(r1.alignment_start(), Position::new(1000));
        let cigar_r1 = format_cigar(&r1.cigar());
        assert_eq!(cigar_r1, "100M");

        assert_eq!(r2.alignment_start(), Position::new(1));
        let cigar_r2 = format_cigar(&r2.cigar());
        assert_eq!(cigar_r2, "100M");
    }

    // Tests for clip_extending_past_mate

    // Tests for overlapping reads with insertions

    #[test]
    fn test_clip_overlapping_reads_with_insertions() {
        // Based on Scala test: "handle reads that contain insertions"
        let clipper = RawClipperOnBuf::new(ClippingMode::Soft);
        let r1_seq = "A".repeat(50);
        let r2_seq = "A".repeat(50);
        // R1: 100-147 with 2bp insertion (40M2I8M)
        // R2: 130-169 with 2bp insertion (10M2I38M)
        let (mut r1, mut r2) = create_pair(100, "40M2I8M", &r1_seq, 130, "10M2I38M", &r2_seq);

        // Before clipping, r1.end >= r2.start (they overlap)
        let r1_end_before = 100 + 48; // 40M + 8M = 48 reference bases
        let r2_start_before = 130;
        assert!(r1_end_before >= r2_start_before);

        let (clipped_r1, clipped_r2) = clipper.clip_overlapping_reads(&mut r1, &mut r2);

        // Both should be clipped
        assert!(clipped_r1 > 0 || clipped_r2 > 0);

        // After clipping, they should not overlap (abutting is OK)
        let r1_end_after = usize::from(
            r1.alignment_start().expect("r1 should have alignment start after clipping"),
        ) + cigar_utils::reference_length(&r1.cigar());
        let r2_start_after = usize::from(
            r2.alignment_start().expect("r2 should have alignment start after clipping"),
        );
        assert!(
            r1_end_after <= r2_start_after,
            "After clipping, r1 end ({r1_end_after}) should be at or before r2 start ({r2_start_after})"
        );
    }

    #[test]
    fn test_clip_overlapping_reads_with_multiple_insertions() {
        let clipper = RawClipperOnBuf::new(ClippingMode::Soft);
        let seq = "A".repeat(60);
        // R1: Multiple insertions
        // R2: Multiple insertions
        let (mut r1, mut r2) = create_pair(100, "20M2I20M3I10M", &seq, 120, "15M2I25M3I10M", &seq);

        let (clipped_r1, clipped_r2) = clipper.clip_overlapping_reads(&mut r1, &mut r2);

        // Should clip without errors
        assert!(clipped_r1 > 0 || clipped_r2 > 0);
    }

    #[test]
    fn test_clip_overlapping_reads_insertion_at_overlap_boundary() {
        let clipper = RawClipperOnBuf::new(ClippingMode::Soft);
        let seq = "A".repeat(55);
        // R1: 100-149, with insertion near the end
        // R2: 145-194
        let (mut r1, mut r2) = create_pair(100, "48M2I5M", &seq, 145, "50M", &seq);

        let (clipped_r1, clipped_r2) = clipper.clip_overlapping_reads(&mut r1, &mut r2);

        // Should handle insertion at overlap boundary
        assert!(clipped_r1 > 0 || clipped_r2 > 0);
    }

    // Tests for SoftWithMask mode

    #[test]
    fn test_clip_end_of_alignment_soft_with_mask_masking() {
        let clipper = RawClipperOnBuf::new(ClippingMode::SoftWithMask);
        let mut record = create_test_record("50M", &"A".repeat(50), 1000);

        let clipped = clipper.clip_end_of_alignment(&mut record, 10);
        assert_eq!(clipped, 10);

        let cigar = format_cigar(&record.cigar());
        assert_eq!(cigar, "40M10S");

        // Check that the last 10 bases are masked to N
        let seq = record.sequence();
        let seq_bytes: Vec<u8> = seq.as_ref().to_vec();
        for (i, &base) in seq_bytes.iter().enumerate().skip(40).take(10) {
            assert_eq!(base, b'N', "Base at position {i} should be N");
        }

        // Check that qualities are set to min (2)
        let quals = record.quality_scores();
        let qual_bytes: Vec<u8> = quals.as_ref().to_vec();
        for (i, &qual) in qual_bytes.iter().enumerate().skip(40).take(10) {
            assert_eq!(qual, 2, "Quality at position {i} should be 2");
        }
    }

    #[test]
    fn test_clip_start_of_alignment_soft_with_mask_masking() {
        let clipper = RawClipperOnBuf::new(ClippingMode::SoftWithMask);
        let mut record = create_test_record("50M", &"A".repeat(50), 1000);

        let clipped = clipper.clip_start_of_alignment(&mut record, 10);
        assert_eq!(clipped, 10);

        let cigar = format_cigar(&record.cigar());
        assert_eq!(cigar, "10S40M");

        // Check that the first 10 bases are masked to N
        let seq = record.sequence();
        let seq_bytes: Vec<u8> = seq.as_ref().to_vec();
        for (i, &base) in seq_bytes.iter().enumerate().take(10) {
            assert_eq!(base, b'N', "Base at position {i} should be N");
        }

        // Check that qualities are set to min (2)
        let quals = record.quality_scores();
        let qual_bytes: Vec<u8> = quals.as_ref().to_vec();
        for (i, &qual) in qual_bytes.iter().enumerate().take(10) {
            assert_eq!(qual, 2, "Quality at position {i} should be 2");
        }
    }

    #[test]
    fn test_clip_overlapping_reads_soft_with_mask() {
        let clipper = RawClipperOnBuf::new(ClippingMode::SoftWithMask);
        let seq = "A".repeat(100);
        let (mut r1, mut r2) = create_pair(1, "100M", &seq, 50, "100M", &seq);

        let (clipped_r1, clipped_r2) = clipper.clip_overlapping_reads(&mut r1, &mut r2);

        assert!(clipped_r1 > 0 || clipped_r2 > 0);

        // Check that clipped bases are masked
        if clipped_r1 > 0 {
            let seq = r1.sequence();
            let seq_bytes: Vec<u8> = seq.as_ref().to_vec();
            let start_mask = 100 - clipped_r1;
            for (i, &base) in seq_bytes.iter().enumerate().skip(start_mask).take(100 - start_mask) {
                assert_eq!(base, b'N', "R1 base at position {i} should be N");
            }
        }
    }

    // Edge cases

    #[test]
    fn test_clip_very_short_read() {
        let clipper = RawClipperOnBuf::new(ClippingMode::Soft);
        let mut record = create_test_record("5M", "ACGTA", 1000);

        // Try to clip more than the read length
        let clipped = clipper.clip_end_of_alignment(&mut record, 3);
        assert_eq!(clipped, 3);

        let cigar = format_cigar(&record.cigar());
        assert_eq!(cigar, "2M3S");
    }

    #[test]
    fn test_clip_entire_read_3_prime() {
        let clipper = RawClipperOnBuf::new(ClippingMode::Soft);
        let mut record = create_test_record("10M", &"A".repeat(10), 1000);

        // Clipping all 10 aligned bases leaves no alignment, so fgbio unmaps the read
        // (SamRecordClipper.scala:157-160) instead of emitting an all-soft-clip CIGAR.
        let clipped = clipper.clip_end_of_alignment(&mut record, 10);
        assert_eq!(clipped, 10);

        assert!(record.flags().is_unmapped(), "read should be unmapped");
        assert_eq!(record.cigar().iter().count(), 0, "CIGAR should be empty");
    }

    #[test]
    fn test_clip_overlapping_reads_complex_cigar() {
        let clipper = RawClipperOnBuf::new(ClippingMode::Soft);
        let r1_seq = "A".repeat(70);
        let r2_seq = "A".repeat(70);
        // Complex CIGAR with multiple operation types
        let (mut r1, mut r2) =
            create_pair(100, "20M5I10M3D15M2I10M", &r1_seq, 130, "15M2I20M5D10M3I10M", &r2_seq);

        let (_clipped_r1, _clipped_r2) = clipper.clip_overlapping_reads(&mut r1, &mut r2);

        // Should handle complex CIGAR without panicking (test passes if no panic occurs)
    }

    #[test]
    fn test_clip_with_leading_and_trailing_soft_clips() {
        let clipper = RawClipperOnBuf::new(ClippingMode::Soft);
        let mut record = create_test_record("10S30M10S", &"A".repeat(50), 1000);

        // Clip additional 5 bases from 3' end
        let clipped = clipper.clip_end_of_alignment(&mut record, 5);
        assert_eq!(clipped, 5);

        let cigar = format_cigar(&record.cigar());
        assert_eq!(cigar, "10S25M15S");
    }

    #[test]
    fn test_clip_overlapping_reads_one_base_different_positions() {
        // Test overlapping by exactly one base at different reference positions
        let clipper = RawClipperOnBuf::new(ClippingMode::Soft);
        let seq = "A".repeat(50);

        // Test at low position
        let (mut r1, mut r2) = create_pair(10, "50M", &seq, 59, "50M", &seq);
        let (c1, c2) = clipper.clip_overlapping_reads(&mut r1, &mut r2);
        assert!(c1 > 0 || c2 > 0);

        // Test at high position
        let (mut r1, mut r2) = create_pair(10000, "50M", &seq, 10049, "50M", &seq);
        let (c1, c2) = clipper.clip_overlapping_reads(&mut r1, &mut r2);
        assert!(c1 > 0 || c2 > 0);
    }

    // Tests for NOT upgrading clipping (same mode or going backwards)

    #[test]
    fn test_not_upgrade_soft_to_soft() {
        // Same mode - should not upgrade
        let clipper = RawClipperOnBuf::new(ClippingMode::Soft);
        let seq = "12345678901234567890123456789012345678901234567890";
        let mut record = create_test_record("5S35M10S", seq, 10);

        let result =
            clipper.upgrade_all_clipping(&mut record).expect("upgrade_all_clipping should succeed");
        assert_eq!(result, (0, 0));

        let cigar_str = format_cigar(&record.cigar());
        assert_eq!(cigar_str, "5S35M10S"); // No change
        assert_eq!(record.sequence().len(), 50); // No sequence change
    }

    #[test]
    fn test_not_upgrade_soft_with_mask_to_soft() {
        // Going backwards - should not upgrade
        let clipper = RawClipperOnBuf::new(ClippingMode::Soft);
        let seq = "12345678901234567890123456789012345678901234567890";
        let mut record = create_test_record("5S35M10S", seq, 10);

        let result =
            clipper.upgrade_all_clipping(&mut record).expect("upgrade_all_clipping should succeed");
        assert_eq!(result, (0, 0));

        let cigar_str = format_cigar(&record.cigar());
        assert_eq!(cigar_str, "5S35M10S");
    }

    #[test]
    fn test_upgrade_all_clipping_soft_with_mask_masks_existing_soft_clips() {
        // CLIP3-03: SoftWithMask must mask existing soft-clipped bases to N/min-quality
        // (CIGAR unchanged), matching fgbio upgradeAllClipping. Previously a no-op.
        let clipper = RawClipperOnBuf::new(ClippingMode::SoftWithMask);
        let seq = "ACGTGACGTGACGTGACGTGACGTGACGTGACGTGACGTGACGTGACGTG"; // 50 bases, no N
        let mut record = create_test_record("5S35M10S", seq, 10);

        let result =
            clipper.upgrade_all_clipping(&mut record).expect("upgrade_all_clipping should succeed");
        assert_eq!(result, (5, 10), "returns the (leading, trailing) soft counts masked");

        // CIGAR is unchanged; only the bases/quals of the soft-clipped regions change.
        assert_eq!(format_cigar(&record.cigar()), "5S35M10S");

        let bases = record.sequence().as_ref().to_vec();
        assert!(bases[..5].iter().all(|&b| b == NO_CALL_BASE), "leading soft bases masked to N");
        assert!(bases[40..].iter().all(|&b| b == NO_CALL_BASE), "trailing soft bases masked to N");
        assert!(bases[5..40].iter().all(|&b| b != NO_CALL_BASE), "aligned bases untouched");

        let quals = record.quality_scores().as_ref().to_vec();
        assert!(quals[..5].iter().all(|&q| q == MIN_PHRED), "leading soft quals masked");
        assert!(quals[40..].iter().all(|&q| q == MIN_PHRED), "trailing soft quals masked");
    }

    #[test]
    fn test_upgrade_all_clipping_raw_soft_with_mask_masks_existing_soft_clips() {
        // CLIP3-03 (raw): same masking behavior for the raw-byte upgrade path.
        use fgumi_raw_bam::encode_record_buf_to_raw;
        use noodles::sam::header::record::value::Map;
        use noodles::sam::header::record::value::map::ReferenceSequence;
        use std::num::NonZeroUsize;

        let seq = "ACGTGACGTGACGTGACGTGACGTGACGTGACGTGACGTGACGTGACGTG"; // 50 bases, no N
        let buf = create_test_record("5S35M10S", seq, 10);
        let ref_seq =
            Map::<ReferenceSequence>::new(NonZeroUsize::new(100_000).expect("ref length nonzero"));
        let header =
            noodles::sam::Header::builder().add_reference_sequence(b"chr1", ref_seq).build();
        let mut raw = encode_record_buf_to_raw(&buf, &header).expect("encode raw");

        let clipper = RawRecordClipper::new(ClippingMode::SoftWithMask);
        let result =
            clipper.upgrade_all_clipping_raw(&mut raw).expect("upgrade_all_clipping_raw succeeds");
        assert_eq!(result, (5, 10));

        // CIGAR unchanged: 5S 35M 10S.
        assert_eq!(raw.cigar_ops_vec(), vec![(5u32 << 4) | 4, 35u32 << 4, (10u32 << 4) | 4]);

        let bases = raw.sequence_vec();
        assert!(bases[..5].iter().all(|&b| b == fgumi_dna::NO_CALL_BASE), "leading masked to N");
        assert!(bases[40..].iter().all(|&b| b == fgumi_dna::NO_CALL_BASE), "trailing masked to N");
        assert!(bases[5..40].iter().all(|&b| b != fgumi_dna::NO_CALL_BASE), "aligned untouched");

        let quals = raw.quality_scores().to_vec();
        assert!(quals[..5].iter().all(|&q| q == fgumi_dna::MIN_PHRED), "leading quals masked");
        assert!(quals[40..].iter().all(|&q| q == fgumi_dna::MIN_PHRED), "trailing quals masked");
    }

    #[test]
    fn test_not_upgrade_hard_to_hard() {
        // Same mode - should not upgrade
        let clipper = RawClipperOnBuf::new(ClippingMode::Hard);
        let seq = "12345678901234567890123456789012345"; // 35 bases
        let mut record = create_test_record("5H35M10H", seq, 10);

        let result =
            clipper.upgrade_all_clipping(&mut record).expect("upgrade_all_clipping should succeed");
        assert_eq!(result, (0, 0));

        let cigar_str = format_cigar(&record.cigar());
        assert_eq!(cigar_str, "5H35M10H");
        assert_eq!(record.sequence().len(), 35);
    }

    #[test]
    fn test_not_upgrade_hard_to_soft_with_mask() {
        // Going backwards - should not upgrade (Hard is already most restrictive)
        let clipper = RawClipperOnBuf::new(ClippingMode::SoftWithMask);
        let seq = "12345678901234567890123456789012345";
        let mut record = create_test_record("5H35M10H", seq, 10);

        let result =
            clipper.upgrade_all_clipping(&mut record).expect("upgrade_all_clipping should succeed");
        assert_eq!(result, (0, 0));

        let cigar_str = format_cigar(&record.cigar());
        assert_eq!(cigar_str, "5H35M10H");
    }

    #[test]
    fn test_not_upgrade_hard_to_soft() {
        // Going backwards - should not upgrade
        let clipper = RawClipperOnBuf::new(ClippingMode::Soft);
        let seq = "12345678901234567890123456789012345";
        let mut record = create_test_record("5H35M10H", seq, 10);

        let result =
            clipper.upgrade_all_clipping(&mut record).expect("upgrade_all_clipping should succeed");
        assert_eq!(result, (0, 0));

        let cigar_str = format_cigar(&record.cigar());
        assert_eq!(cigar_str, "5H35M10H");
    }

    #[test]
    fn test_not_upgrade_unmapped_read() {
        // Unmapped reads should not be upgraded
        let clipper = RawClipperOnBuf::new(ClippingMode::Hard);
        let seq = "12345678901234567890123456789012345678901234567890";
        let mut record = create_test_record("5S35M10S", seq, 10);

        // Mark as unmapped
        *record.flags_mut() = Flags::UNMAPPED;

        let result =
            clipper.upgrade_all_clipping(&mut record).expect("upgrade_all_clipping should succeed");
        assert_eq!(result, (0, 0));

        let cigar_str = format_cigar(&record.cigar());
        assert_eq!(cigar_str, "5S35M10S");
    }

    #[test]
    fn test_upgrade_soft_to_hard_with_existing_hard_clips() {
        // Soft with existing hard clips -> Hard mode
        let clipper = RawClipperOnBuf::new(ClippingMode::Hard);
        let seq = "1234567890123456789012345678901234567890"; // 40 bases (2H + 5S + 30M + 3S)
        let mut record = create_test_record("2H5S30M3S", seq, 10);

        let result =
            clipper.upgrade_all_clipping(&mut record).expect("upgrade_all_clipping should succeed");
        assert_eq!(result, (5, 3));

        let cigar_str = format_cigar(&record.cigar());
        assert_eq!(cigar_str, "7H30M3H"); // 2H + 5S -> 7H, 3S -> 3H
        assert_eq!(record.sequence().len(), 30);
    }

    #[test]
    fn test_upgrade_soft_with_mask_to_hard() {
        // SoftWithMask -> Hard should upgrade
        let clipper = RawClipperOnBuf::new(ClippingMode::Hard);
        let seq = "12345678901234567890123456789012345678901234567890";
        let mut record = create_test_record("5S35M10S", seq, 10);

        let result =
            clipper.upgrade_all_clipping(&mut record).expect("upgrade_all_clipping should succeed");
        assert_eq!(result, (5, 10));

        let cigar_str = format_cigar(&record.cigar());
        assert_eq!(cigar_str, "5H35M10H");
        assert_eq!(record.sequence().len(), 35);
    }

    #[test]
    fn test_upgrade_clipping_with_only_leading_soft_clip() {
        let clipper = RawClipperOnBuf::new(ClippingMode::Hard);
        let seq = "12345678901234567890123456789012345678901234567890";
        let mut record = create_test_record("10S40M", seq, 10);

        let result =
            clipper.upgrade_all_clipping(&mut record).expect("upgrade_all_clipping should succeed");
        assert_eq!(result, (10, 0));

        let cigar_str = format_cigar(&record.cigar());
        assert_eq!(cigar_str, "10H40M");
        assert_eq!(record.sequence().len(), 40);
    }

    #[test]
    fn test_upgrade_clipping_with_only_trailing_soft_clip() {
        let clipper = RawClipperOnBuf::new(ClippingMode::Hard);
        let seq = "12345678901234567890123456789012345678901234567890";
        let mut record = create_test_record("40M10S", seq, 10);

        let result =
            clipper.upgrade_all_clipping(&mut record).expect("upgrade_all_clipping should succeed");
        assert_eq!(result, (0, 10));

        let cigar_str = format_cigar(&record.cigar());
        assert_eq!(cigar_str, "40M10H");
        assert_eq!(record.sequence().len(), 40);
    }

    #[test]
    fn test_upgrade_clipping_complex_cigar_with_indels() {
        let clipper = RawClipperOnBuf::new(ClippingMode::Hard);
        let seq = "123456789012345678901234567890123456789012345"; // 45 bases (5S + 20M + 5I + 10M + 5S)
        let mut record = create_test_record("5S20M5I10M5S", seq, 10);

        let result =
            clipper.upgrade_all_clipping(&mut record).expect("upgrade_all_clipping should succeed");
        assert_eq!(result, (5, 5));

        let cigar_str = format_cigar(&record.cigar());
        assert_eq!(cigar_str, "5H20M5I10M5H");
        assert_eq!(record.sequence().len(), 35); // 20 + 5 + 10
    }

    #[test]
    fn test_upgrade_clipping_no_soft_clips_present() {
        let clipper = RawClipperOnBuf::new(ClippingMode::Hard);
        let seq = "12345678901234567890123456789012345678901234567890";
        let mut record = create_test_record("50M", seq, 10);

        let result =
            clipper.upgrade_all_clipping(&mut record).expect("upgrade_all_clipping should succeed");
        assert_eq!(result, (0, 0));

        let cigar_str = format_cigar(&record.cigar());
        assert_eq!(cigar_str, "50M");
        assert_eq!(record.sequence().len(), 50);
    }

    // ===== Additional tests ported from fgbio to match clipExtendingPastMateEnds coverage =====

    #[test]
    fn test_clip_5_prime_end_of_alignment_positive_strand() {
        // fgbio test: "add more clipping to the 5' end" - positive strand
        let clipper = RawClipperOnBuf::new(ClippingMode::Soft);

        // Positive strand, no existing clipping
        let mut rec1 = create_test_record("50M", &"A".repeat(50), 10);
        *rec1.flags_mut() = Flags::empty(); // positive strand
        let clipped = clipper.clip_5_prime_end_of_alignment(&mut rec1, 10);
        assert_eq!(clipped, 10);
        assert_eq!(format_cigar(&rec1.cigar()), "10S40M");

        // Positive strand, existing soft clipping
        let mut rec2 = create_test_record("10S40M", &"A".repeat(50), 10);
        *rec2.flags_mut() = Flags::empty(); // positive strand
        let clipped = clipper.clip_5_prime_end_of_alignment(&mut rec2, 10);
        assert_eq!(clipped, 10);
        assert_eq!(format_cigar(&rec2.cigar()), "20S30M");
    }

    #[test]
    fn test_clip_5_prime_end_of_alignment_negative_strand() {
        // fgbio test: "add more clipping to the 5' end" - negative strand
        let clipper = RawClipperOnBuf::new(ClippingMode::Soft);

        // Negative strand, no existing clipping
        let mut rec3 = create_test_record("50M", &"A".repeat(50), 10);
        *rec3.flags_mut() = Flags::REVERSE_COMPLEMENTED;
        let clipped = clipper.clip_5_prime_end_of_alignment(&mut rec3, 10);
        assert_eq!(clipped, 10);
        assert_eq!(format_cigar(&rec3.cigar()), "40M10S");

        // Negative strand, existing soft clipping
        let mut rec4 = create_test_record("40M10S", &"A".repeat(50), 10);
        *rec4.flags_mut() = Flags::REVERSE_COMPLEMENTED;
        let clipped = clipper.clip_5_prime_end_of_alignment(&mut rec4, 10);
        assert_eq!(clipped, 10);
        assert_eq!(format_cigar(&rec4.cigar()), "30M20S");
    }

    #[test]
    fn test_clip_3_prime_end_of_alignment_negative_strand() {
        // fgbio test: "add more clipping to the 3' end" - negative strand
        let clipper = RawClipperOnBuf::new(ClippingMode::Soft);

        // Negative strand, no existing clipping
        let mut rec1 = create_test_record("50M", &"A".repeat(50), 10);
        *rec1.flags_mut() = Flags::REVERSE_COMPLEMENTED;
        let clipped = clipper.clip_3_prime_end_of_alignment(&mut rec1, 10);
        assert_eq!(clipped, 10);
        assert_eq!(format_cigar(&rec1.cigar()), "10S40M");

        // Negative strand, existing soft clipping
        let mut rec2 = create_test_record("10S40M", &"A".repeat(50), 10);
        *rec2.flags_mut() = Flags::REVERSE_COMPLEMENTED;
        let clipped = clipper.clip_3_prime_end_of_alignment(&mut rec2, 10);
        assert_eq!(clipped, 10);
        assert_eq!(format_cigar(&rec2.cigar()), "20S30M");
    }

    #[test]
    fn test_clip_3_prime_end_of_alignment_positive_strand() {
        // fgbio test: "add more clipping to the 3' end" - positive strand
        let clipper = RawClipperOnBuf::new(ClippingMode::Soft);

        // Positive strand, no existing clipping
        let mut rec3 = create_test_record("50M", &"A".repeat(50), 10);
        *rec3.flags_mut() = Flags::empty(); // positive strand
        let clipped = clipper.clip_3_prime_end_of_alignment(&mut rec3, 10);
        assert_eq!(clipped, 10);
        assert_eq!(format_cigar(&rec3.cigar()), "40M10S");

        // Positive strand, existing soft clipping
        let mut rec4 = create_test_record("40M10S", &"A".repeat(50), 10);
        *rec4.flags_mut() = Flags::empty(); // positive strand
        let clipped = clipper.clip_3_prime_end_of_alignment(&mut rec4, 10);
        assert_eq!(clipped, 10);
        assert_eq!(format_cigar(&rec4.cigar()), "30M20S");
    }

    // ===================================================================
    // fgbio#1172 query-space regression for the RawRecordClipper past-mate path
    // (clip_extending_past_mate_ends / num_bases_extending_past_mate_vs_mate_raw),
    // ported from fgbio's SamRecordClipperTest.scala. See fgumi-raw-bam's
    // bases_extending_past_mate_ops (Task 1) for the shared query-space core.
    // ===================================================================

    /// Local strand marker for the fgbio#1172 ported cases, mapped to the
    /// `REVERSE_COMPLEMENTED` / `MATE_REVERSE_COMPLEMENTED` flags.
    #[derive(Clone, Copy, Debug, PartialEq, Eq)]
    enum Strand {
        Plus,
        Minus,
    }

    impl Strand {
        fn is_reverse(self) -> bool {
            self == Strand::Minus
        }
    }

    /// Signed template length (insert size) for a record given its own strand/position/
    /// reference length and its mate's position/reference length; matches the sign
    /// convention `record_utils::is_fr_pair` (via `get_pair_orientation`) relies on.
    fn template_length_for(
        is_reverse: bool,
        own_pos: usize,
        own_ref_len: usize,
        mate_pos: usize,
        mate_ref_len: usize,
    ) -> i32 {
        if is_reverse {
            let own_end = i32::try_from(own_pos + own_ref_len.saturating_sub(1))
                .expect("own_end fits in i32");
            let mate_pos = i32::try_from(mate_pos).expect("mate_pos fits in i32");
            -(own_end - mate_pos + 1)
        } else {
            let mate_end = i32::try_from(mate_pos + mate_ref_len.saturating_sub(1))
                .expect("mate_end fits in i32");
            let own_pos = i32::try_from(own_pos).expect("own_pos fits in i32");
            mate_end - own_pos + 1
        }
    }

    /// Builds an FR-oriented `RecordBuf` pair for the fgbio#1172 ported cases: sets
    /// paired/first-second/reverse flags, mate ref+pos+strand, `MC` tags, and a
    /// `template_length` consistent with the given positions/CIGARs (needed for
    /// `record_utils::is_fr_pair`). Sequence/qualities are auto-generated by
    /// `RecordBuilder` to match each CIGAR's query length.
    fn fr_pair(rec1: (&str, usize, Strand), rec2: (&str, usize, Strand)) -> (RecordBuf, RecordBuf) {
        let (cigar1, pos1, strand1) = rec1;
        let (cigar2, pos2, strand2) = rec2;

        let mut r1 = RecordBuilder::mapped_read()
            .cigar(cigar1)
            .alignment_start(pos1)
            .paired(true)
            .first_segment(true)
            .reverse_complement(strand1.is_reverse())
            .mate_reverse_complement(strand2.is_reverse())
            .mate_reference_sequence_id(0)
            .mate_alignment_start(pos2)
            .tag("MC", cigar2)
            .build();
        let mut r2 = RecordBuilder::mapped_read()
            .cigar(cigar2)
            .alignment_start(pos2)
            .paired(true)
            .first_segment(false)
            .reverse_complement(strand2.is_reverse())
            .mate_reverse_complement(strand1.is_reverse())
            .mate_reference_sequence_id(0)
            .mate_alignment_start(pos1)
            .tag("MC", cigar1)
            .build();

        let ref_len1 = cigar_utils::reference_length(&r1.cigar());
        let ref_len2 = cigar_utils::reference_length(&r2.cigar());
        *r1.template_length_mut() =
            template_length_for(strand1.is_reverse(), pos1, ref_len1, pos2, ref_len2);
        *r2.template_length_mut() =
            template_length_for(strand2.is_reverse(), pos2, ref_len2, pos1, ref_len1);

        (r1, r2)
    }

    // ===================================================================
    // fgbio#1172 query-space regression, RawRecordClipper: this is the LIVE `fgumi clip`
    // path (RawRecordClipper::clip_extending_past_mate_ends). SamRecordClipper never had
    // this fixed and its RecordBuf-only past-mate methods have been removed as redundant
    // (issue #760) — this is the single remaining home for the query-space count.
    // ===================================================================

    /// Build a SAM header with a single reference sequence long enough for the fgbio#1172
    /// ported cases (`disjoint_no_shared_pos` uses position 1020).
    fn raw_test_header() -> noodles::sam::Header {
        use noodles::sam::header::record::value::Map;
        use noodles::sam::header::record::value::map::ReferenceSequence;
        use std::num::NonZeroUsize;
        let ref_seq = Map::<ReferenceSequence>::new(
            NonZeroUsize::new(100_000).expect("ref length must be nonzero"),
        );
        noodles::sam::Header::builder().add_reference_sequence(b"chr1", ref_seq).build()
    }

    /// Builds an FR-oriented raw pair for the fgbio#1172 ported cases: reuses [`fr_pair`]'s
    /// `RecordBuf` construction (paired/first-second/reverse flags, mate ref+pos+strand, `MC`
    /// tags, consistent `template_length`) and encodes both records to `RawRecord` bytes, since
    /// `RawRecordClipper::clip_extending_past_mate_ends` is the LIVE `fgumi clip` path this
    /// issue fixes.
    fn raw_fr_pair(
        rec1: (&str, usize, Strand),
        rec2: (&str, usize, Strand),
    ) -> (fgumi_raw_bam::RawRecord, fgumi_raw_bam::RawRecord) {
        let (r1_buf, r2_buf) = fr_pair(rec1, rec2);
        let header = raw_test_header();
        let r1 = fgumi_raw_bam::encode_record_buf_to_raw(&r1_buf, &header).expect("encode r1");
        let r2 = fgumi_raw_bam::encode_record_buf_to_raw(&r2_buf, &header).expect("encode r2");
        (r1, r2)
    }

    /// Format a raw CIGAR (BAM `(len<<4)|op` u32s) as a SAM CIGAR string, for asserting the
    /// exact post-clip shapes in the fgbio#1172 ported cases.
    fn raw_format_cigar(ops: &[u32]) -> String {
        use std::fmt::Write as _;
        let mut out = String::new();
        for &op in ops {
            let len = op >> 4;
            let kind_char = match op & 0xF {
                0 => 'M',
                1 => 'I',
                2 => 'D',
                3 => 'N',
                4 => 'S',
                5 => 'H',
                6 => 'P',
                7 => '=',
                8 => 'X',
                _ => '?',
            };
            let _ = write!(out, "{len}{kind_char}");
        }
        out
    }

    /// Table A (fgbio#1172, ported verbatim, `RawRecordClipper`): `clip_extending_past_mate_ends`
    /// returns additional `(clip_r1, clip_r2)`. Confirms the LIVE `fgumi clip` path (#760).
    #[rstest]
    // Control, Soft: no indel. rec 100M@100/+, mate 90M10S@60/-. Expect (40,40); rec->60M40S, mate->40S50M10S.
    #[case::no_indel_control(ClippingMode::Soft, ("100M", 100, Strand::Plus), ("90M10S", 60, Strand::Minus), (40, 40))]
    // Over-clip regression, Soft: deletion at mate's un-soft-clipped end. OLD unmapped rec. NEW (2,0).
    #[case::deletion_at_mate_end(ClippingMode::Soft, ("2S124M1D3M", 101, Strand::Plus), ("3S124M2S", 100, Strand::Minus), (2, 0))]
    // Knock-on, Soft: rec's deletion at mate end; mate must still be clipped (order-independent). (2,2).
    #[case::deletion_knock_on(ClippingMode::Soft, ("2S124M1D3M", 101, Strand::Plus), ("115M14S", 97, Strand::Minus), (2, 2))]
    // Under-clip regression, Hard: insertion before mate end (the #1090 example). (3,0); rec->70M10I20M50H.
    #[case::insertion_before_mate_end(ClippingMode::Hard, ("70M10I23M47S", 100, Strand::Plus), ("50S70M30S", 100, Strand::Minus), (3, 0))]
    // Disjoint continuity, Hard: share NO ref position (mate at 1020). Bounded by soft clip -> 0 aligned clipped.
    #[case::disjoint_no_shared_pos(ClippingMode::Hard, ("20M80S", 1000, Strand::Plus), ("80S20M", 1020, Strand::Minus), (0, 0))]
    fn raw_clip_extending_past_mate_ends_query_space(
        #[case] mode: ClippingMode,
        #[case] rec1: (&str, usize, Strand),
        #[case] rec2: (&str, usize, Strand),
        #[case] expected: (usize, usize),
    ) {
        let (mut r1, mut r2) = raw_fr_pair(rec1, rec2);
        let clipper = RawRecordClipper::new(mode);
        assert_eq!(clipper.clip_extending_past_mate_ends(&mut r1, &mut r2), expected);
        // Regression: neither read may be unmapped by a past-mate clip.
        assert!(r1.alignment_start_1based().is_some());
        assert!(r2.alignment_start_1based().is_some());
    }

    /// Hard-clip add-back case (fgbio#1172), `RawRecordClipper`: pre-clip each 3' end by 10
    /// (upgrading to hard clip via `clip_3_prime_end_of_read_raw`), THEN measure/clip past-mate.
    /// `clip_3_prime_end_of_read_raw`/`clip_end_of_read_raw` take TOTAL desired clipping and
    /// subtract existing hard+soft, so the second call must add the already-hard-clipped bases
    /// back before re-requesting, or it under-clips by that amount (see
    /// `existing_hard_clip_3_prime_raw`, applied unconditionally in `clip_extending_past_mate_ends`).
    #[test]
    fn raw_clip_extending_past_mate_ends_query_space_already_hard_clipped_3prime() {
        let clipper = RawRecordClipper::new(ClippingMode::Hard);
        let (mut r1, mut r2) =
            raw_fr_pair(("100M", 100, Strand::Plus), ("90M10S", 60, Strand::Minus));

        clipper.clip_3_prime_end_of_read_raw(&mut r1, 10);
        clipper.clip_3_prime_end_of_read_raw(&mut r2, 10);
        assert_eq!(raw_format_cigar(&r1.cigar_ops_vec()), "90M10H");
        assert_eq!(raw_format_cigar(&r2.cigar_ops_vec()), "10H80M10S");

        let result = clipper.clip_extending_past_mate_ends(&mut r1, &mut r2);
        assert_eq!(result, (30, 30));
        assert_eq!(raw_format_cigar(&r1.cigar_ops_vec()), "60M40H");
        assert_eq!(raw_format_cigar(&r2.cigar_ops_vec()), "40H50M10S");
    }

    /// Table B (fgbio#1172, ported verbatim, `RawRecordClipper`):
    /// `num_bases_extending_past_mate_vs_mate_raw(rec, mate)` count, checked in both directions.
    #[rstest]
    // Indel near end: deletion (r1) and its mate. Shared through ref 223 / first-shared 101. Expect 2 & 2.
    #[case::count_deletion(("2S124M1D3M", 101, Strand::Plus), ("3S124M2S", 100, Strand::Minus), 2, 2)]
    // Indel near end: insertion. Shared through 169 / first-shared 100. Expect 50 & 50.
    #[case::count_insertion(("70M10I23M47S", 100, Strand::Plus), ("50S70M30S", 100, Strand::Minus), 50, 50)]
    // Disjoint sharing exactly one ref position (mate at 1019): 61 & 61.
    #[case::count_disjoint_one_shared(("20M80S", 1000, Strand::Plus), ("80S20M", 1019, Strand::Minus), 61, 61)]
    // Disjoint sharing none (mate at 1020): must be 60 & 60, not 0 (continuity).
    #[case::count_disjoint_none(("20M80S", 1000, Strand::Plus), ("80S20M", 1020, Strand::Minus), 60, 60)]
    fn raw_num_bases_extending_past_mate_counts(
        #[case] rec1: (&str, usize, Strand),
        #[case] rec2: (&str, usize, Strand),
        #[case] expect_r1: usize,
        #[case] expect_r2: usize,
    ) {
        let (r1, r2) = raw_fr_pair(rec1, rec2);
        assert_eq!(
            fgumi_raw_bam::num_bases_extending_past_mate_vs_mate_raw(r1.as_ref(), r2.as_ref()),
            expect_r1
        );
        assert_eq!(
            fgumi_raw_bam::num_bases_extending_past_mate_vs_mate_raw(r2.as_ref(), r1.as_ref()),
            expect_r2
        );
    }

    // ===================================================================
    // fgbio#1172 query-space regression: property tests over a generated FR-pair
    // strategy, modeled on the fgbio#1172 sweep geometry (overlapping / touching /
    // disjoint mate blocks, ~half carrying a single indel near an alignment end).
    // Built on the `fr_pair`/`RecordBuilder` machinery above so every generated
    // pair clears the same FR-pair gate and carries the same MC-tag/mate-field
    // wiring as the ported fgbio#1172 cases.
    // ===================================================================

    /// One 100bp read's CIGAR shape: leading/trailing soft-clip lengths (each 0-80, summing to
    /// at most 90 so at least 10 aligned bases always remain) plus an optional single indel
    /// (insertion or deletion, 1-6bp) placed adjacent to either the leading or trailing aligned
    /// edge. Mirrors the shapes exercised by the fgbio#1172 ported cases above
    /// (`no_indel_control` = no indel; `deletion_at_mate_end`/`insertion_before_mate_end` = an
    /// indel adjacent to the mate-facing end).
    #[derive(Clone, Copy, Debug)]
    struct CigarShape {
        leading_clip: usize,
        trailing_clip: usize,
        /// `(near_start, is_insertion, len)`, or `None` for a plain `{aligned}M` interior.
        indel: Option<(bool, bool, usize)>,
    }

    impl CigarShape {
        const READ_LEN: usize = 100;

        /// Renders this shape as a CIGAR string whose query length is exactly `READ_LEN`
        /// (`RecordBuilder::cigar` auto-generates SEQ/QUAL to match, so no explicit sequence is
        /// needed).
        fn to_cigar(self) -> String {
            use std::fmt::Write as _;

            let avail = Self::READ_LEN - self.leading_clip - self.trailing_clip;
            let mut cigar = String::new();
            if self.leading_clip > 0 {
                write!(cigar, "{}S", self.leading_clip).expect("write! to String never fails");
            }
            match self.indel {
                None => write!(cigar, "{avail}M").expect("write! to String never fails"),
                Some((near_start, true, len)) => {
                    // Insertion consumes query: the aligned (M) total shrinks by `len`.
                    // `avail >= 10` and `len <= 6` by construction, so `m_total >= 4`.
                    let m_total = avail - len;
                    let (m1, m2) = if near_start { (1, m_total - 1) } else { (m_total - 1, 1) };
                    write!(cigar, "{m1}M{len}I{m2}M").expect("write! to String never fails");
                }
                Some((near_start, false, len)) => {
                    // Deletion doesn't consume query: the aligned (M) total is `avail`.
                    let (m1, m2) = if near_start { (1, avail - 1) } else { (avail - 1, 1) };
                    write!(cigar, "{m1}M{len}D{m2}M").expect("write! to String never fails");
                }
            }
            if self.trailing_clip > 0 {
                write!(cigar, "{}S", self.trailing_clip).expect("write! to String never fails");
            }
            cigar
        }
    }

    /// Strategy for [`CigarShape`]: soft clips independently in `0..=80`, constrained so their
    /// sum never exceeds 90 (guaranteeing `avail >= 10`), plus a single indel roughly half the
    /// time.
    fn arb_cigar_shape() -> impl Strategy<Value = CigarShape> {
        (0..=80usize)
            .prop_flat_map(|leading_clip| {
                let max_trailing = (90usize.saturating_sub(leading_clip)).min(80);
                (Just(leading_clip), 0..=max_trailing)
            })
            .prop_flat_map(|(leading_clip, trailing_clip)| {
                let indel_strategy = prop_oneof![
                    1 => Just(None),
                    1 => (any::<bool>(), any::<bool>(), 1..=6usize)
                        .prop_map(|(near_start, is_insertion, len)| Some((
                            near_start,
                            is_insertion,
                            len
                        ))),
                ];
                (Just(leading_clip), Just(trailing_clip), indel_strategy)
            })
            .prop_map(|(leading_clip, trailing_clip, indel)| CigarShape {
                leading_clip,
                trailing_clip,
                indel,
            })
    }

    /// How r1's and r2's aligned reference blocks relate: their reference blocks share bases
    /// (`Overlapping`), are immediately adjacent and share none (`Touching`), or are separated
    /// by a gap (`Disjoint`). Mirrors the geometries exercised by the fgbio#1172 ported cases
    /// above (`deletion_knock_on` overlaps, `disjoint_no_shared_pos` touches,
    /// `count_disjoint_none` is disjoint).
    #[derive(Clone, Copy, Debug)]
    enum PairGeometry {
        /// Overlap magnitude in bases, clamped below `min(ref_len1, ref_len2)` so some
        /// reference span always separates r1's start from r2's end (required for
        /// `record_utils::is_fr_pair`'s forward-5'-before-reverse-5' check). Drawn up to the
        /// 100bp read length so the clamp binds and the extreme "one alignment almost entirely
        /// swallows the other" geometry — which historically produced the full-read-unmap bug
        /// the `never_unmaps` property guards against — is exercised, not just shallow overlaps.
        Overlapping(usize),
        Touching,
        /// Gap in bases strictly separating the two reference blocks.
        Disjoint(usize),
    }

    fn arb_pair_geometry() -> impl Strategy<Value = PairGeometry> {
        prop_oneof![
            (0..100usize).prop_map(PairGeometry::Overlapping),
            Just(PairGeometry::Touching),
            (0..60usize).prop_map(PairGeometry::Disjoint),
        ]
    }

    /// Strategy for a valid FR-oriented `(RecordBuf, RecordBuf)` pair: two independently
    /// generated 100bp `CigarShape`s (r1 forward, r2 reverse) positioned per a generated
    /// [`PairGeometry`], then built via [`fr_pair`] so mate ref/pos/strand fields, `MC` tags,
    /// and `template_length` all come from the same helper the fgbio#1172 ported cases above
    /// use (keeping this generator's pairs indistinguishable, from the clipper's point of view,
    /// from a real FR pair).
    fn arb_fr_pair() -> impl Strategy<Value = (RecordBuf, RecordBuf)> {
        // Large enough that `pos1 + ref_len1 + offset` never underflows for any generated
        // overlap/gap magnitude.
        const BASE_POS: usize = 10_000;

        (arb_cigar_shape(), arb_cigar_shape(), arb_pair_geometry()).prop_map(
            |(shape1, shape2, geometry)| {
                let cigar1 = shape1.to_cigar();
                let cigar2 = shape2.to_cigar();
                let ref_len1 = crate::builder::cigar_ref_len(&cigar1);
                let ref_len2 = crate::builder::cigar_ref_len(&cigar2);

                let pos1 = BASE_POS;
                let offset: isize = match geometry {
                    PairGeometry::Overlapping(magnitude) => {
                        let cap = ref_len1.min(ref_len2).saturating_sub(1);
                        -isize::try_from(magnitude.min(cap)).expect("small overlap fits isize")
                    }
                    PairGeometry::Touching => 0,
                    PairGeometry::Disjoint(gap) => {
                        isize::try_from(gap + 1).expect("small gap fits isize")
                    }
                };
                let pos1_end = isize::try_from(pos1 + ref_len1).expect("fits isize");
                let pos2 =
                    usize::try_from(pos1_end + offset).expect("pos2 stays positive for BASE_POS");

                fr_pair((&cigar1, pos1, Strand::Plus), (&cigar2, pos2, Strand::Minus))
            },
        )
    }

    /// Strategy for a dovetail read-through `(RecordBuf, RecordBuf)` pair: each 100bp read aligns
    /// a short block and soft-clips the rest, with the two aligned blocks meeting near one shared
    /// reference position (`shift` jitters that meeting point). Both reads therefore read through
    /// past the far end of the mate — the short-insert / adapter-read-through geometry that
    /// `clip_extending_past_mate_ends` actually clips, and that [`arb_fr_pair`]'s
    /// forward-left/reverse-right layout never produces. Mirrors the fgbio#1172
    /// `disjoint_no_shared_pos` / `count_disjoint_*` cases parametrically.
    fn arb_readthrough_fr_pair() -> impl Strategy<Value = (RecordBuf, RecordBuf)> {
        const BASE_POS: usize = 10_000;
        (10..=40usize, 10..=40usize, -3isize..=3isize).prop_map(|(a, b, shift)| {
            let cigar1 = format!("{a}M{}S", CigarShape::READ_LEN - a); // fwd: short block + read-through
            let cigar2 = format!("{}S{b}M", CigarShape::READ_LEN - b); // rev: read-through + short block
            let pos1 = BASE_POS;
            // Reverse aligned block starts at (roughly) the forward aligned block's end, so the
            // two share ~one reference position and their soft-clipped tails cross.
            let pos2 = usize::try_from(isize::try_from(pos1 + a).expect("fits isize") - 1 + shift)
                .expect("pos2 stays positive for BASE_POS");
            fr_pair((&cigar1, pos1, Strand::Plus), (&cigar2, pos2, Strand::Minus))
        })
    }

    /// Strategy for an "overhang" read-through pair that clips *aligned* bases: the forward read
    /// is fully aligned (`100M`) and extends past the far end of a reverse mate that starts to its
    /// left and soft-clips its trailing (5') bases. The forward read's aligned 3' therefore runs
    /// past the mate's un-soft-clipped end, so `clip_extending_past_mate_ends` removes real aligned
    /// bases (a nonzero clip) rather than merely re-clipping soft-clipped tails. Mirrors the
    /// fgbio#1172 `no_indel_control` case (`100M` + `90M10S` -> `(40, 40)`), which
    /// [`arb_fr_pair`]'s cap — which keeps the forward block from ever reaching past the reverse
    /// block — cannot reach.
    fn arb_overhang_fr_pair() -> impl Strategy<Value = (RecordBuf, RecordBuf)> {
        const BASE_POS: usize = 10_000;
        // rev_aligned: length of the reverse read's aligned (M) block; back_shift: how far left of
        // the forward read the reverse block starts. Bounds keep the forward 5' strictly left of
        // the reverse 5' (FR-valid) while the forward 3' overhangs the reverse's un-soft-clipped end.
        (20..=90usize, 1..=39usize).prop_map(|(rev_aligned, back_shift)| {
            let cigar1 = format!("{}M", CigarShape::READ_LEN); // fully-aligned forward, overhangs right
            let cigar2 = format!("{rev_aligned}M{}S", CigarShape::READ_LEN - rev_aligned);
            let pos1 = BASE_POS;
            let pos2 = pos1 - back_shift; // reverse block starts to the left of the forward block
            fr_pair((&cigar1, pos1, Strand::Plus), (&cigar2, pos2, Strand::Minus))
        })
    }

    /// Strategy for a valid FR-oriented `(RawRecord, RawRecord)` pair. Mixes three geometries and
    /// encodes both records to `RawRecord` bytes, exercising the LIVE `fgumi clip` path
    /// (`RawRecordClipper`): [`arb_fr_pair`]'s aligned-block sweep (mostly no-clip),
    /// [`arb_readthrough_fr_pair`]'s dovetail (nonzero read-through *count*, clip absorbed by soft
    /// clipping), and [`arb_overhang_fr_pair`]'s overhang (nonzero *aligned-base* clip). See
    /// `arb_raw_fr_pair_exercises_nonzero_clips`.
    fn arb_raw_fr_pair()
    -> impl Strategy<Value = (fgumi_raw_bam::RawRecord, fgumi_raw_bam::RawRecord)> {
        prop_oneof![arb_fr_pair(), arb_readthrough_fr_pair(), arb_overhang_fr_pair()].prop_map(
            |(r1_buf, r2_buf)| {
                let header = raw_test_header();
                let r1 = fgumi_raw_bam::encode_record_buf_to_raw(&r1_buf, &header)
                    .expect("encode r1 to raw");
                let r2 = fgumi_raw_bam::encode_record_buf_to_raw(&r2_buf, &header)
                    .expect("encode r2 to raw");
                (r1, r2)
            },
        )
    }

    proptest! {
        /// A past-mate clip must never unmap an already-mapped FR read (RawRecordClipper).
        #[test]
        fn raw_past_mate_clip_never_unmaps(pair in arb_raw_fr_pair()) {
            let (mut r1, mut r2) = pair;
            let clipper = RawRecordClipper::new(ClippingMode::Soft);
            let _ = clipper.clip_extending_past_mate_ends(&mut r1, &mut r2);
            prop_assert!(r1.alignment_start_1based().is_some());
            prop_assert!(r2.alignment_start_1based().is_some());
        }

        /// A second pass over an already-clipped pair must clip nothing further
        /// (RawRecordClipper).
        #[test]
        fn raw_past_mate_clip_is_idempotent(pair in arb_raw_fr_pair()) {
            let (mut r1, mut r2) = pair;
            let clipper = RawRecordClipper::new(ClippingMode::Soft);
            let _ = clipper.clip_extending_past_mate_ends(&mut r1, &mut r2);
            prop_assert_eq!(clipper.clip_extending_past_mate_ends(&mut r1, &mut r2), (0, 0));
        }

        /// Swapping the two records passed in swaps the returned clip amounts
        /// (RawRecordClipper).
        #[test]
        fn raw_past_mate_clip_is_symmetric(pair in arb_raw_fr_pair()) {
            let (mut a1, mut a2) = pair.clone();
            let (mut b1, mut b2) = pair;
            let clipper = RawRecordClipper::new(ClippingMode::Soft);
            let (c1, c2) = clipper.clip_extending_past_mate_ends(&mut a1, &mut a2);
            let (d2, d1) = clipper.clip_extending_past_mate_ends(&mut b2, &mut b1);
            prop_assert_eq!((c1, c2), (d1, d2));
        }
    }

    /// The `never_unmaps` / `idempotent` / `symmetric` properties above would pass vacuously if
    /// the generator rarely produced read-through past the mate. Pin that it does not: over a
    /// deterministic sample, a meaningful fraction of generated pairs must clip a nonzero number
    /// of bases from at least one read, so the properties are exercised against real clipping
    /// rather than a stream of no-op pairs.
    #[test]
    fn arb_raw_fr_pair_exercises_nonzero_clips() {
        use proptest::strategy::{Strategy, ValueTree};
        use proptest::test_runner::TestRunner;

        let mut runner = TestRunner::deterministic();
        let clipper = RawRecordClipper::new(ClippingMode::Soft);
        let strategy = arb_raw_fr_pair();
        let total = 512;
        let mut nonzero = 0;
        for _ in 0..total {
            let (mut r1, mut r2) =
                strategy.new_tree(&mut runner).expect("generate a pair").current();
            let (c1, c2) = clipper.clip_extending_past_mate_ends(&mut r1, &mut r2);
            if c1 > 0 || c2 > 0 {
                nonzero += 1;
            }
        }
        assert!(
            nonzero >= total / 10,
            "generator exercises read-through too rarely: only {nonzero}/{total} pairs clipped \
             nonzero — the past-mate properties risk passing vacuously"
        );
    }
}

// ============================================================================
// Round-trip fidelity tests: RawRecordClipper direct == via RecordBuf adapter
// ============================================================================
//
// These tests once cross-checked `RawRecordClipper` against the removed
// `SamRecordClipper`. They now compare the raw clipper run directly on a
// `RawRecord` against the same clip run through the `RawClipperOnBuf` adapter
// (`RecordBuf` -> raw -> clip -> `RecordBuf`), confirming the buf<->raw round-trip
// preserves CIGAR, alignment start, sequence, and quality across all three
// `ClippingMode` values for a variety of representative inputs.
//
// Because both compared paths now run the SAME `RawRecordClipper`, this module guards only
// buf<->raw round-trip fidelity -- it cannot catch a clip-logic regression (a bug would produce
// identical wrong output on both sides and still pass here). Clip-arithmetic correctness rests
// on the concrete-expectation assertions in `mod tests` (hard-coded fgbio-derived CIGAR/SEQ
// values), not on this module.
// ============================================================================

#[cfg(test)]
mod crosscheck_tests {
    use super::clip_test_adapter::RawClipperOnBuf;
    use super::*;
    use crate::builder::RecordBuilder;
    use fgumi_raw_bam::{RawRecord, encode_record_buf_to_raw};
    use noodles::sam::alignment::RecordBuf;
    use noodles::sam::alignment::record::Sequence as SequenceTrait;
    use noodles::sam::alignment::record::cigar::Cigar;
    use noodles::sam::alignment::record::cigar::op::Kind;
    use rstest::rstest;

    /// Build a SAM header with a single reference sequence of length 100 000.
    fn test_header() -> noodles::sam::Header {
        use noodles::sam::header::record::value::Map;
        use noodles::sam::header::record::value::map::ReferenceSequence;
        use std::num::NonZeroUsize;
        let ref_seq = Map::<ReferenceSequence>::new(
            NonZeroUsize::new(100_000).expect("ref length must be nonzero"),
        );
        noodles::sam::Header::builder().add_reference_sequence(b"chr1", ref_seq).build()
    }

    /// Encode a `RecordBuf` to a `RawRecord` using the shared test header.
    fn to_raw(record: &RecordBuf) -> RawRecord {
        let header = test_header();
        encode_record_buf_to_raw(record, &header).expect("encode_record_buf_to_raw should succeed")
    }

    /// Format a raw CIGAR (vec of u32) as a human-readable SAM CIGAR string.
    fn raw_cigar_str(ops: &[u32]) -> String {
        use std::fmt::Write;
        let mut out = String::new();
        for &op in ops {
            let len = op >> 4;
            let kind_char = match op & 0xF {
                0 => 'M',
                1 => 'I',
                2 => 'D',
                3 => 'N',
                4 => 'S',
                5 => 'H',
                6 => 'P',
                7 => '=',
                8 => 'X',
                _ => '?',
            };
            let _ = write!(out, "{len}{kind_char}");
        }
        out
    }

    /// Format a noodles CIGAR as a SAM string.
    fn buf_cigar_str(cigar: &noodles::sam::alignment::record_buf::Cigar) -> String {
        use std::fmt::Write;
        cigar.as_ref().iter().fold(String::new(), |mut acc, op| {
            let kind_char = match op.kind() {
                Kind::Match => 'M',
                Kind::Insertion => 'I',
                Kind::Deletion => 'D',
                Kind::Skip => 'N',
                Kind::SoftClip => 'S',
                Kind::HardClip => 'H',
                Kind::Pad => 'P',
                Kind::SequenceMatch => '=',
                Kind::SequenceMismatch => 'X',
            };
            let _ = write!(acc, "{}{}", op.len(), kind_char);
            acc
        })
    }

    /// Assert that a raw record has the same CIGAR, pos, seq, and qual as a `RecordBuf`.
    fn assert_raw_matches_buf(raw: &RawRecord, buf: &RecordBuf, context: &str) {
        let raw_cigar = raw_cigar_str(&raw.cigar_ops_vec());
        let buf_cigar = buf_cigar_str(buf.cigar());
        assert_eq!(raw_cigar, buf_cigar, "{context}: CIGAR mismatch");

        let raw_pos = raw.pos(); // 0-based
        // buf.alignment_start() is 1-based; convert to 0-based for comparison with raw.pos()
        let buf_pos = i32::try_from(buf.alignment_start().map_or(0, usize::from))
            .expect("alignment start fits in i32")
            - 1;
        assert_eq!(raw_pos, buf_pos, "{context}: pos mismatch");

        let raw_seq = raw.sequence_vec();
        let buf_seq: Vec<u8> = buf.sequence().iter().collect();
        assert_eq!(raw_seq, buf_seq, "{context}: sequence mismatch");

        let raw_qual = raw.quality_scores().to_vec();
        let buf_qual: Vec<u8> = buf.quality_scores().as_ref().to_vec();
        assert_eq!(raw_qual, buf_qual, "{context}: quality mismatch");
    }

    // =========================================================================
    // BAM bin field is recomputed after raw clip mutations (fgbio parity)
    //
    // The raw pipeline emits record bytes verbatim, so a clip that moves POS or
    // rewrites the CIGAR must refresh the bin (bytes 10-11) itself. A read placed
    // at 16300 with 200M straddles the 16 kb bin boundary (level-4 bin 585);
    // clipping that changes its span moves it into a level-5 (16 kb) bin.
    // =========================================================================

    /// Build a raw mapped read at `pos` (0-based) with `cigar`, encoded through
    /// noodles so the starting bin is correct.
    fn raw_mapped(pos: i32, cigar: &str, seq_len: usize) -> RawRecord {
        let start_1based = usize::try_from(pos + 1).expect("1-based start must be non-negative");
        let buf = RecordBuilder::mapped_read()
            .sequence(&"A".repeat(seq_len))
            .cigar(cigar)
            .alignment_start(start_1based)
            .build();
        to_raw(&buf)
    }

    #[test]
    fn raw_clip_start_of_alignment_recomputes_bin_after_pos_shift() {
        let mut raw = raw_mapped(16300, "200M", 200);
        assert_eq!(raw.bin(), fgumi_raw_bam::reg2bin(16300, 16500), "precondition: bin 585");
        // Soft-clip 116 from the start: pos -> 16416, CIGAR 116S84M, span [16416, 16500).
        RawRecordClipper::new(ClippingMode::Soft).clip_start_of_alignment(&mut raw, 116);
        assert_eq!(raw.pos(), 16416, "precondition: pos shifted");
        assert_eq!(raw.bin(), fgumi_raw_bam::reg2bin(16416, 16416 + 84));
    }

    #[test]
    fn raw_clip_end_of_alignment_recomputes_bin_after_span_shrinks() {
        let mut raw = raw_mapped(16300, "200M", 200);
        // Soft-clip 116 from the end: pos stays 16300, span shrinks to [16300, 16384),
        // which no longer straddles the 16 kb boundary -> level-5 bin.
        RawRecordClipper::new(ClippingMode::Soft).clip_end_of_alignment(&mut raw, 116);
        assert_eq!(raw.pos(), 16300, "precondition: end-clip leaves pos unchanged");
        assert_eq!(raw.bin(), fgumi_raw_bam::reg2bin(16300, 16300 + 84));
    }

    #[test]
    fn raw_clip_that_unmaps_read_sets_unmapped_bin() {
        let mut raw = raw_mapped(16300, "200M", 200);
        // Clipping every clippable base unmaps the read (fgbio SamRecordClipper).
        RawRecordClipper::new(ClippingMode::Soft).clip_start_of_alignment(&mut raw, 200);
        assert!(raw.flags() & fgumi_raw_bam::flags::UNMAPPED != 0, "precondition: unmapped");
        assert_eq!(raw.bin(), fgumi_raw_bam::UNMAPPED_BIN);
    }

    // =========================================================================
    // clip_start_of_alignment cross-checks
    // =========================================================================

    #[test]
    fn crosscheck_clip_start_soft_simple_match() {
        for bases in [5usize, 10, 15] {
            for mode in [ClippingMode::Soft, ClippingMode::SoftWithMask, ClippingMode::Hard] {
                let mut buf = RecordBuilder::mapped_read()
                    .sequence(&"A".repeat(50))
                    .cigar("50M")
                    .alignment_start(100)
                    .build();
                let mut raw = to_raw(&buf);

                let buf_clipped =
                    RawClipperOnBuf::new(mode).clip_start_of_alignment(&mut buf, bases);
                let raw_clipped =
                    RawRecordClipper::new(mode).clip_start_of_alignment(&mut raw, bases);

                assert_eq!(
                    buf_clipped, raw_clipped,
                    "mode={mode:?} bases={bases}: clip count mismatch"
                );
                assert_raw_matches_buf(
                    &raw,
                    &buf,
                    &format!("clip_start 50M mode={mode:?} bases={bases}"),
                );
            }
        }
    }

    #[test]
    fn crosscheck_clip_start_soft_with_insertion() {
        for mode in [ClippingMode::Soft, ClippingMode::SoftWithMask, ClippingMode::Hard] {
            let mut buf = RecordBuilder::mapped_read()
                .sequence(&"A".repeat(50))
                .cigar("4M2I44M")
                .alignment_start(100)
                .build();
            let mut raw = to_raw(&buf);

            let buf_clipped = RawClipperOnBuf::new(mode).clip_start_of_alignment(&mut buf, 10);
            let raw_clipped = RawRecordClipper::new(mode).clip_start_of_alignment(&mut raw, 10);

            assert_eq!(buf_clipped, raw_clipped, "mode={mode:?}: clip count");
            assert_raw_matches_buf(&raw, &buf, &format!("clip_start 4M2I44M mode={mode:?}"));
        }
    }

    /// N (ref skip) handling must agree between raw and typed clippers so spliced reads
    /// get the same clip count and alignment start on both code paths.
    #[test]
    fn crosscheck_clip_start_with_skip_n() {
        for mode in [ClippingMode::Soft, ClippingMode::SoftWithMask, ClippingMode::Hard] {
            let mut buf = RecordBuilder::mapped_read()
                .sequence(&"A".repeat(50))
                .cigar("10M100N40M")
                .alignment_start(100)
                .build();
            let mut raw = to_raw(&buf);

            let buf_clipped = RawClipperOnBuf::new(mode).clip_start_of_alignment(&mut buf, 5);
            let raw_clipped = RawRecordClipper::new(mode).clip_start_of_alignment(&mut raw, 5);

            assert_eq!(buf_clipped, raw_clipped, "mode={mode:?}: clip count");
            assert_raw_matches_buf(&raw, &buf, &format!("clip_start 10M100N40M mode={mode:?}"));
        }
    }

    #[test]
    fn crosscheck_clip_end_with_skip_n() {
        for mode in [ClippingMode::Soft, ClippingMode::SoftWithMask, ClippingMode::Hard] {
            let mut buf = RecordBuilder::mapped_read()
                .sequence(&"A".repeat(50))
                .cigar("40M100N10M")
                .alignment_start(100)
                .build();
            let mut raw = to_raw(&buf);

            let buf_clipped = RawClipperOnBuf::new(mode).clip_end_of_alignment(&mut buf, 5);
            let raw_clipped = RawRecordClipper::new(mode).clip_end_of_alignment(&mut raw, 5);

            assert_eq!(buf_clipped, raw_clipped, "mode={mode:?}: clip count");
            assert_raw_matches_buf(&raw, &buf, &format!("clip_end 40M100N10M mode={mode:?}"));
        }
    }

    #[test]
    fn crosscheck_clip_start_soft_with_deletion() {
        for mode in [ClippingMode::Soft, ClippingMode::SoftWithMask, ClippingMode::Hard] {
            let mut buf = RecordBuilder::mapped_read()
                .sequence(&"A".repeat(50))
                .cigar("6M2D44M")
                .alignment_start(100)
                .build();
            let mut raw = to_raw(&buf);

            let buf_clipped = RawClipperOnBuf::new(mode).clip_start_of_alignment(&mut buf, 10);
            let raw_clipped = RawRecordClipper::new(mode).clip_start_of_alignment(&mut raw, 10);

            assert_eq!(buf_clipped, raw_clipped, "mode={mode:?}: clip count");
            assert_raw_matches_buf(&raw, &buf, &format!("clip_start 6M2D44M mode={mode:?}"));
        }
    }

    #[test]
    fn crosscheck_clip_start_existing_clips() {
        for mode in [ClippingMode::Soft, ClippingMode::SoftWithMask, ClippingMode::Hard] {
            // Existing soft clip
            let mut buf = RecordBuilder::mapped_read()
                .sequence(&"A".repeat(50))
                .cigar("10S40M")
                .alignment_start(100)
                .build();
            let mut raw = to_raw(&buf);

            let buf_clipped = RawClipperOnBuf::new(mode).clip_start_of_alignment(&mut buf, 10);
            let raw_clipped = RawRecordClipper::new(mode).clip_start_of_alignment(&mut raw, 10);

            assert_eq!(buf_clipped, raw_clipped, "10S40M mode={mode:?}");
            assert_raw_matches_buf(&raw, &buf, &format!("clip_start 10S40M mode={mode:?}"));
        }
    }

    #[test]
    fn crosscheck_clip_start_existing_hard_clip() {
        for mode in [ClippingMode::Soft, ClippingMode::Hard] {
            let mut buf = RecordBuilder::mapped_read()
                .sequence(&"A".repeat(40))
                .cigar("10H40M")
                .alignment_start(100)
                .build();
            let mut raw = to_raw(&buf);

            let buf_clipped = RawClipperOnBuf::new(mode).clip_start_of_alignment(&mut buf, 10);
            let raw_clipped = RawRecordClipper::new(mode).clip_start_of_alignment(&mut raw, 10);

            assert_eq!(buf_clipped, raw_clipped, "10H40M mode={mode:?}");
            assert_raw_matches_buf(&raw, &buf, &format!("clip_start 10H40M mode={mode:?}"));
        }
    }

    // =========================================================================
    // clip_end_of_alignment cross-checks
    // =========================================================================

    #[test]
    fn crosscheck_clip_end_soft_simple_match() {
        for bases in [5usize, 10, 15] {
            for mode in [ClippingMode::Soft, ClippingMode::SoftWithMask, ClippingMode::Hard] {
                let mut buf = RecordBuilder::mapped_read()
                    .sequence(&"A".repeat(50))
                    .cigar("50M")
                    .alignment_start(100)
                    .build();
                let mut raw = to_raw(&buf);

                let buf_clipped = RawClipperOnBuf::new(mode).clip_end_of_alignment(&mut buf, bases);
                let raw_clipped =
                    RawRecordClipper::new(mode).clip_end_of_alignment(&mut raw, bases);

                assert_eq!(buf_clipped, raw_clipped, "mode={mode:?} bases={bases}");
                assert_raw_matches_buf(
                    &raw,
                    &buf,
                    &format!("clip_end 50M mode={mode:?} bases={bases}"),
                );
            }
        }
    }

    #[test]
    fn crosscheck_clip_end_with_deletion() {
        for mode in [ClippingMode::Soft, ClippingMode::SoftWithMask, ClippingMode::Hard] {
            // 44M2D4M: query length = 44+4 = 48
            let mut buf = RecordBuilder::mapped_read()
                .sequence(&"A".repeat(48))
                .cigar("44M2D4M")
                .alignment_start(100)
                .build();
            let mut raw = to_raw(&buf);

            let buf_clipped = RawClipperOnBuf::new(mode).clip_end_of_alignment(&mut buf, 10);
            let raw_clipped = RawRecordClipper::new(mode).clip_end_of_alignment(&mut raw, 10);

            assert_eq!(buf_clipped, raw_clipped, "mode={mode:?}");
            assert_raw_matches_buf(&raw, &buf, &format!("clip_end 44M2D4M mode={mode:?}"));
        }
    }

    #[test]
    fn crosscheck_clip_end_trailing_insertion() {
        for mode in [ClippingMode::Soft, ClippingMode::SoftWithMask, ClippingMode::Hard] {
            let mut buf = RecordBuilder::mapped_read()
                .sequence(&"A".repeat(50))
                .cigar("38M4I8M")
                .alignment_start(100)
                .build();
            let mut raw = to_raw(&buf);

            let buf_clipped = RawClipperOnBuf::new(mode).clip_end_of_alignment(&mut buf, 10);
            let raw_clipped = RawRecordClipper::new(mode).clip_end_of_alignment(&mut raw, 10);

            assert_eq!(buf_clipped, raw_clipped, "mode={mode:?}");
            assert_raw_matches_buf(&raw, &buf, &format!("clip_end 38M4I8M mode={mode:?}"));
        }
    }

    // =========================================================================
    // clip_start_of_read / clip_end_of_read cross-checks
    // =========================================================================

    #[test]
    fn crosscheck_clip_start_of_read_upgrade() {
        // clip_length <= existing clipping => upgrade_clipping is triggered
        // 10S40M: query length = 50
        let mode = ClippingMode::Hard;
        let mut buf = RecordBuilder::mapped_read()
            .sequence(&"A".repeat(50))
            .cigar("10S40M")
            .alignment_start(100)
            .build();
        let mut raw = to_raw(&buf);

        let buf_clipped = RawClipperOnBuf::new(mode).clip_start_of_read(&mut buf, 5);
        let raw_clipped = RawRecordClipper::new(mode).clip_start_of_read_raw(&mut raw, 5);

        assert_eq!(buf_clipped, raw_clipped, "10S40M Hard upgrade start");
        assert_raw_matches_buf(&raw, &buf, "clip_start_of_read 10S40M Hard upgrade");
    }

    #[test]
    fn crosscheck_clip_end_of_read_upgrade() {
        // 40M10S: query length = 50
        let mode = ClippingMode::Hard;
        let mut buf = RecordBuilder::mapped_read()
            .sequence(&"A".repeat(50))
            .cigar("40M10S")
            .alignment_start(100)
            .build();
        let mut raw = to_raw(&buf);

        let buf_clipped = RawClipperOnBuf::new(mode).clip_end_of_read(&mut buf, 5);
        let raw_clipped = RawRecordClipper::new(mode).clip_end_of_read_raw(&mut raw, 5);

        assert_eq!(buf_clipped, raw_clipped, "40M10S Hard upgrade end");
        assert_raw_matches_buf(&raw, &buf, "clip_end_of_read 40M10S Hard upgrade");
    }

    /// Raw upgrade-only path must trim per-base tags to match the shrunken
    /// sequence, matching `clip_extended_attributes` on the typed path.
    #[test]
    fn clip_start_of_read_raw_upgrade_trims_per_base_tags() {
        // 10S40M: query length = 50, upgrade 5 of the existing soft clips to hard.
        let buf = RecordBuilder::mapped_read()
            .sequence(&"A".repeat(50))
            .cigar("10S40M")
            .alignment_start(100)
            .build();
        let mut raw = to_raw(&buf);

        // Attach a per-base B:c tag of length 50 to the raw record.
        {
            let mut ed = raw.tags_editor();
            let values: Vec<i8> = (0i8..50i8).collect();
            ed.append_array_i8(b"XC", &values);
        }

        let clipper = RawRecordClipper::with_auto_clip(ClippingMode::Hard, true);
        clipper.clip_start_of_read_raw(&mut raw, 5);

        // Sequence shrank by 5; the per-base B:c must also shrink by 5 (from the start).
        let aux = fgumi_raw_bam::aux_data_slice(raw.as_ref());
        let arr = fgumi_raw_bam::find_array_tag(aux, b"XC").expect("XC present");
        assert_eq!(arr.elem_type, b'c', "B:c subtype preserved");
        assert_eq!(arr.count, 45, "B:c count matches new sequence length");
        let decoded: Vec<i8> = arr.data.iter().map(|&b| b.cast_signed()).collect();
        let expected: Vec<i8> = (5i8..50i8).collect();
        assert_eq!(decoded, expected, "B:c trimmed from the start");
    }

    // =========================================================================
    // upgrade_all_clipping cross-checks
    // =========================================================================

    // Hard converts the soft clips to hard clips (removing bases); SoftWithMask masks the
    // soft-clipped bases to N/min-quality while leaving the CIGAR intact. Both are checked for
    // raw/typed parity here. Soft is a documented no-op for upgrade and is omitted.
    #[rstest]
    #[case::hard(ClippingMode::Hard)]
    #[case::soft_with_mask(ClippingMode::SoftWithMask)]
    fn crosscheck_upgrade_all_clipping(#[case] mode: ClippingMode) {
        // 10S30M10S: query length = 50
        let mut buf = RecordBuilder::mapped_read()
            .sequence(&"A".repeat(50))
            .cigar("10S30M10S")
            .alignment_start(100)
            .build();
        let mut raw = to_raw(&buf);

        let buf_result = RawClipperOnBuf::new(mode).upgrade_all_clipping(&mut buf);
        let raw_result = RawRecordClipper::new(mode).upgrade_all_clipping_raw(&mut raw);

        assert_eq!(
            buf_result.unwrap(),
            raw_result.unwrap(),
            "upgrade_all_clipping return value mode={mode:?}",
        );
        assert_raw_matches_buf(
            &raw,
            &buf,
            &format!("upgrade_all_clipping mode={mode:?} 10S30M10S"),
        );
    }

    // With no soft clips both modes early-return (0, 0) and leave the record untouched; check
    // raw/typed parity for each.
    #[rstest]
    #[case::hard(ClippingMode::Hard)]
    #[case::soft_with_mask(ClippingMode::SoftWithMask)]
    fn crosscheck_upgrade_all_clipping_no_soft_clips(#[case] mode: ClippingMode) {
        let mut buf = RecordBuilder::mapped_read()
            .sequence(&"A".repeat(50))
            .cigar("50M")
            .alignment_start(100)
            .build();
        let mut raw = to_raw(&buf);

        let buf_result = RawClipperOnBuf::new(mode).upgrade_all_clipping(&mut buf);
        let raw_result = RawRecordClipper::new(mode).upgrade_all_clipping_raw(&mut raw);

        assert_eq!(
            buf_result.unwrap(),
            raw_result.unwrap(),
            "no soft clips return value mode={mode:?}",
        );
        assert_raw_matches_buf(
            &raw,
            &buf,
            &format!("upgrade_all_clipping mode={mode:?} 50M (no-op)"),
        );
    }

    /// Ref-only CIGAR ops (D/N/P) must not advance `seq_pos` in the raw path,
    /// or trailing query bases get shifted vs. the typed path. Uses a varied
    /// sequence and varied qualities so any mis-aligned copy (Hard) or mis-aligned
    /// mask (`SoftWithMask`) is visible. Both modes are checked for raw/typed parity.
    #[rstest]
    #[case::hard(ClippingMode::Hard)]
    #[case::soft_with_mask(ClippingMode::SoftWithMask)]
    fn crosscheck_upgrade_all_clipping_with_deletion(#[case] mode: ClippingMode) {
        // 10S20M5D20M10S: query length = 60 (D does not consume query).
        let bases = b"ACGT";
        let seq: String = (0..60).map(|i| bases[i % bases.len()] as char).collect();
        let quals: Vec<u8> = (0..60u8).map(|i| 33 + (i % 40)).collect();
        let mut buf = RecordBuilder::mapped_read()
            .sequence(&seq)
            .qualities(&quals)
            .cigar("10S20M5D20M10S")
            .alignment_start(100)
            .build();
        let mut raw = to_raw(&buf);

        let buf_result = RawClipperOnBuf::new(mode).upgrade_all_clipping(&mut buf);
        let raw_result = RawRecordClipper::new(mode).upgrade_all_clipping_raw(&mut raw);

        assert_eq!(
            buf_result.unwrap(),
            raw_result.unwrap(),
            "upgrade_all_clipping return value with deletion mode={mode:?}",
        );
        assert_raw_matches_buf(
            &raw,
            &buf,
            &format!("upgrade_all_clipping mode={mode:?} 10S20M5D20M10S"),
        );
    }

    // =========================================================================
    // clip_5_prime / clip_3_prime strand-aware cross-checks
    // =========================================================================

    #[test]
    fn crosscheck_clip_5prime_alignment_positive_strand() {
        for mode in [ClippingMode::Soft, ClippingMode::Hard] {
            let mut buf = RecordBuilder::mapped_read()
                .sequence(&"A".repeat(50))
                .cigar("50M")
                .alignment_start(100)
                .build();
            let mut raw = to_raw(&buf);

            let buf_c = RawClipperOnBuf::new(mode).clip_5_prime_end_of_alignment(&mut buf, 10);
            let raw_c = RawRecordClipper::new(mode).clip_5_prime_end_of_alignment(&mut raw, 10);

            assert_eq!(buf_c, raw_c, "5' positive mode={mode:?}");
            assert_raw_matches_buf(&raw, &buf, &format!("clip_5prime pos mode={mode:?}"));
        }
    }

    #[test]
    fn crosscheck_clip_3prime_alignment_positive_strand() {
        for mode in [ClippingMode::Soft, ClippingMode::Hard] {
            let mut buf = RecordBuilder::mapped_read()
                .sequence(&"A".repeat(50))
                .cigar("50M")
                .alignment_start(100)
                .build();
            let mut raw = to_raw(&buf);

            let buf_c = RawClipperOnBuf::new(mode).clip_3_prime_end_of_alignment(&mut buf, 10);
            let raw_c = RawRecordClipper::new(mode).clip_3_prime_end_of_alignment(&mut raw, 10);

            assert_eq!(buf_c, raw_c, "3' positive mode={mode:?}");
            assert_raw_matches_buf(&raw, &buf, &format!("clip_3prime pos mode={mode:?}"));
        }
    }

    // =========================================================================
    // Unmap-when-fully-clipped cross-checks (CLIP-01)
    //
    // fgbio `SamRecordClipper.clip{Start,End}OfAlignment` unmaps the read (via
    // htsjdk `SAMUtils.makeReadUnmapped`) when the requested clip is >= the number
    // of clippable (alignment + insertion) bases, keeping bases/quals intact.
    // Refs: SamRecordClipper.scala:106-128,152-173; SamRecordClipperTest.scala:168,311.
    // =========================================================================

    /// Assert both clippers put the record into the fgbio `makeReadUnmapped` state
    /// and that the raw and typed records agree on sequence and qualities.
    fn assert_both_unmapped(raw: &RawRecord, buf: &RecordBuf, context: &str) {
        use fgumi_raw_bam::flags as rflags;

        // --- RecordBuf (typed) ---
        assert!(buf.flags().is_unmapped(), "{context}: buf not unmapped");
        assert!(!buf.flags().is_reverse_complemented(), "{context}: buf reverse not cleared");
        assert!(!buf.flags().is_duplicate(), "{context}: buf duplicate not cleared");
        assert!(!buf.flags().is_secondary(), "{context}: buf secondary not cleared");
        assert!(!buf.flags().is_supplementary(), "{context}: buf supplementary not cleared");
        assert!(!buf.flags().is_properly_segmented(), "{context}: buf proper-pair not cleared");
        assert_eq!(buf.reference_sequence_id(), None, "{context}: buf ref id");
        assert_eq!(buf.alignment_start(), None, "{context}: buf alignment start");
        // MAPQ must be 0 (htsjdk NO_MAPPING_QUALITY), NOT the "unavailable" sentinel
        // that a typed `None` serializes to (255). This must equal the raw path's 0
        // below — the two clippers agree on the same on-wire MAPQ.
        assert_eq!(
            buf.mapping_quality().map(u8::from),
            Some(0u8),
            "{context}: buf mapq should be 0 (htsjdk makeReadUnmapped), matching raw path"
        );
        assert_eq!(buf.template_length(), 0, "{context}: buf tlen");
        assert_eq!(buf.cigar().iter().count(), 0, "{context}: buf cigar not empty");

        // --- RawRecord ---
        assert!(raw.flags() & rflags::UNMAPPED != 0, "{context}: raw not unmapped");
        assert!(raw.flags() & rflags::REVERSE == 0, "{context}: raw reverse not cleared");
        assert!(raw.flags() & rflags::DUPLICATE == 0, "{context}: raw duplicate not cleared");
        assert!(raw.flags() & rflags::SECONDARY == 0, "{context}: raw secondary not cleared");
        assert!(
            raw.flags() & rflags::SUPPLEMENTARY == 0,
            "{context}: raw supplementary not cleared"
        );
        assert!(raw.flags() & rflags::PROPER_PAIR == 0, "{context}: raw proper-pair not cleared");
        assert_eq!(raw.ref_id(), -1, "{context}: raw ref id");
        assert_eq!(raw.pos(), -1, "{context}: raw pos");
        assert_eq!(raw.mapq(), 0, "{context}: raw mapq");
        assert_eq!(raw.template_length(), 0, "{context}: raw tlen");
        assert!(raw.cigar_ops_vec().is_empty(), "{context}: raw cigar not empty");

        // --- Cross-check seq & qual (bases/quals are preserved on unmap) ---
        let buf_seq: Vec<u8> = buf.sequence().iter().collect();
        assert_eq!(raw.sequence_vec(), buf_seq, "{context}: seq mismatch");
        assert_eq!(
            raw.quality_scores().to_vec(),
            buf.quality_scores().as_ref().to_vec(),
            "{context}: qual mismatch"
        );
    }

    #[rstest]
    fn crosscheck_unmap_clip_start_exceeds_alignment(
        #[values(ClippingMode::Soft, ClippingMode::SoftWithMask, ClippingMode::Hard)]
        mode: ClippingMode,
    ) {
        // fgbio SamRecordClipperTest.scala:168 — 10S40M, clip 50 from the start.
        // Only 40 clippable bases (the 40M); 40 <= 50 => unmap, return 40.
        let mut buf = RecordBuilder::mapped_read()
            .sequence(&"A".repeat(50))
            .cigar("10S40M")
            .alignment_start(100)
            .mapping_quality(60)
            .build();
        let mut raw = to_raw(&buf);

        let buf_c = RawClipperOnBuf::new(mode).clip_start_of_alignment(&mut buf, 50);
        let raw_c = RawRecordClipper::new(mode).clip_start_of_alignment(&mut raw, 50);

        assert_eq!(buf_c, 40, "buf return count mode={mode:?}");
        assert_eq!(raw_c, 40, "raw return count mode={mode:?}");
        assert_eq!(buf.sequence().iter().count(), 50, "buf seq len preserved mode={mode:?}");
        assert_eq!(raw.l_seq(), 50, "raw seq len preserved mode={mode:?}");
        assert_both_unmapped(&raw, &buf, &format!("clip_start 10S40M clip=50 mode={mode:?}"));
    }

    #[rstest]
    fn crosscheck_unmap_clip_end_exceeds_alignment(
        #[values(ClippingMode::Soft, ClippingMode::SoftWithMask, ClippingMode::Hard)]
        mode: ClippingMode,
    ) {
        // fgbio SamRecordClipperTest.scala:311 — 40M10S, clip 50 from the end.
        let mut buf = RecordBuilder::mapped_read()
            .sequence(&"A".repeat(50))
            .cigar("40M10S")
            .alignment_start(100)
            .mapping_quality(60)
            .build();
        let mut raw = to_raw(&buf);

        let buf_c = RawClipperOnBuf::new(mode).clip_end_of_alignment(&mut buf, 50);
        let raw_c = RawRecordClipper::new(mode).clip_end_of_alignment(&mut raw, 50);

        assert_eq!(buf_c, 40, "buf return count mode={mode:?}");
        assert_eq!(raw_c, 40, "raw return count mode={mode:?}");
        assert_both_unmapped(&raw, &buf, &format!("clip_end 40M10S clip=50 mode={mode:?}"));
    }

    #[rstest]
    fn crosscheck_unmap_boundary_clip_equals_clippable(
        #[values(ClippingMode::Soft, ClippingMode::Hard)] mode: ClippingMode,
    ) {
        // Boundary: requesting exactly the clippable-base count also unmaps
        // (fgbio uses `numClippable <= numberOfBasesToClip`).
        let mut buf = RecordBuilder::mapped_read()
            .sequence(&"A".repeat(50))
            .cigar("50M")
            .alignment_start(100)
            .mapping_quality(60)
            .build();
        let mut raw = to_raw(&buf);

        let buf_c = RawClipperOnBuf::new(mode).clip_start_of_alignment(&mut buf, 50);
        let raw_c = RawRecordClipper::new(mode).clip_start_of_alignment(&mut raw, 50);

        assert_eq!(buf_c, 50, "buf return count mode={mode:?}");
        assert_eq!(raw_c, 50, "raw return count mode={mode:?}");
        assert_both_unmapped(&raw, &buf, &format!("clip_start 50M clip=50 mode={mode:?}"));
    }

    #[rstest]
    fn crosscheck_no_unmap_when_clip_below_clippable(
        #[values(ClippingMode::Soft, ClippingMode::Hard)] mode: ClippingMode,
    ) {
        // Regression guard: clipping fewer than the clippable bases must NOT unmap.
        let mut buf = RecordBuilder::mapped_read()
            .sequence(&"A".repeat(50))
            .cigar("50M")
            .alignment_start(100)
            .mapping_quality(60)
            .build();
        let mut raw = to_raw(&buf);

        let buf_c = RawClipperOnBuf::new(mode).clip_start_of_alignment(&mut buf, 49);
        let raw_c = RawRecordClipper::new(mode).clip_start_of_alignment(&mut raw, 49);

        assert_eq!(buf_c, 49, "buf return count mode={mode:?}");
        assert_eq!(raw_c, 49, "raw return count mode={mode:?}");
        assert!(!buf.flags().is_unmapped(), "buf must stay mapped mode={mode:?}");
        assert!(raw.flags() & fgumi_raw_bam::flags::UNMAPPED == 0, "raw must stay mapped");
        assert_raw_matches_buf(&raw, &buf, &format!("clip_start 50M clip=49 mode={mode:?}"));
    }

    #[rstest]
    fn crosscheck_unmap_negative_strand_revcomps_bases(
        #[values(ClippingMode::Soft, ClippingMode::Hard)] mode: ClippingMode,
    ) {
        // A negative-strand read that gets fully clipped is unmapped with its stored
        // (reference-forward) SEQ reverse-complemented and QUAL reversed, mirroring
        // htsjdk `SAMUtils.makeReadUnmapped` (reverseComplement + clear reverse flag).
        let seq = format!("{}{}", "A".repeat(25), "C".repeat(25));
        let expected_rc = fgumi_dna::reverse_complement(seq.as_bytes());
        let mut buf = RecordBuilder::mapped_read()
            .sequence(&seq)
            .cigar("50M")
            .alignment_start(100)
            .mapping_quality(60)
            .reverse_complement(true)
            .build();
        let mut raw = to_raw(&buf);

        RawClipperOnBuf::new(mode).clip_start_of_alignment(&mut buf, 60);
        RawRecordClipper::new(mode).clip_start_of_alignment(&mut raw, 60);

        assert_both_unmapped(&raw, &buf, &format!("neg-strand unmap mode={mode:?}"));
        assert_eq!(raw.sequence_vec(), expected_rc, "neg-strand seq revcomp mode={mode:?}");
    }

    /// Reads a `Z`-type string tag back from a typed [`RecordBuf`].
    fn typed_string_tag(buf: &RecordBuf, tag: [u8; 2]) -> Vec<u8> {
        use noodles::sam::alignment::record::data::field::Tag;
        use noodles::sam::alignment::record_buf::data::field::Value;
        match buf.data().get(&Tag::new(tag[0], tag[1])) {
            Some(Value::String(s)) => s.to_vec(),
            other => {
                panic!("expected string tag {:?}, got {other:?}", String::from_utf8_lossy(&tag))
            }
        }
    }

    /// Reads a `Z`-type string tag back from a [`RawRecord`].
    fn raw_string_tag(raw: &RawRecord, tag: [u8; 2]) -> Vec<u8> {
        fgumi_raw_bam::find_string_tag(fgumi_raw_bam::aux_data_slice(raw.as_ref()), tag)
            .unwrap_or_else(|| panic!("expected string tag {:?}", String::from_utf8_lossy(&tag)))
            .to_vec()
    }

    #[rstest]
    fn crosscheck_unmap_negative_strand_reorients_tags(
        #[values(ClippingMode::Soft, ClippingMode::Hard)] mode: ClippingMode,
    ) {
        use fgumi_raw_bam::SamTag;
        // htsjdk `SAMUtils.makeReadUnmapped` calls `reverseComplement(true)`, which — for
        // reverse-strand reads — also reverse-complements the base-oriented tags
        // (E2, SQ) and reverses the quality-oriented tags (OQ, U2), matching
        // `SAMRecord.TAGS_TO_REVERSE_COMPLEMENT` / `TAGS_TO_REVERSE`.
        let oq = format!("{}{}", "5".repeat(25), "F".repeat(25)); // reversed -> F*25 5*25
        let u2 = format!("{}{}", "#".repeat(25), "2".repeat(25)); // reversed -> 2*25 #*25
        let e2 = format!("{}{}", "A".repeat(25), "C".repeat(25)); // revcomp  -> G*25 T*25
        let expected_oq = format!("{}{}", "F".repeat(25), "5".repeat(25));
        let expected_u2 = format!("{}{}", "2".repeat(25), "#".repeat(25));
        let expected_e2 = format!("{}{}", "G".repeat(25), "T".repeat(25));

        let mut buf = RecordBuilder::mapped_read()
            .sequence(&"A".repeat(50))
            .cigar("50M")
            .alignment_start(100)
            .mapping_quality(60)
            .reverse_complement(true)
            .tag("OQ", oq.as_str())
            .tag("U2", u2.as_str())
            .tag("E2", e2.as_str())
            .build();
        let mut raw = to_raw(&buf);

        RawClipperOnBuf::new(mode).clip_start_of_alignment(&mut buf, 60);
        RawRecordClipper::new(mode).clip_start_of_alignment(&mut raw, 60);

        assert_both_unmapped(&raw, &buf, &format!("neg-strand tag reorient mode={mode:?}"));
        for (tag, want) in
            [(*SamTag::OQ, &expected_oq), (*SamTag::U2, &expected_u2), (*SamTag::E2, &expected_e2)]
        {
            let label = String::from_utf8_lossy(&tag).into_owned();
            assert_eq!(typed_string_tag(&buf, tag), want.as_bytes(), "typed {label} mode={mode:?}");
            assert_eq!(raw_string_tag(&raw, tag), want.as_bytes(), "raw {label} mode={mode:?}");
        }
    }

    #[rstest]
    fn crosscheck_unmap_positive_strand_leaves_tags_untouched(
        #[values(ClippingMode::Soft, ClippingMode::Hard)] mode: ClippingMode,
    ) {
        use fgumi_raw_bam::SamTag;
        // For forward-strand reads, htsjdk `makeReadUnmapped` does NOT call
        // `reverseComplement`, so the base/quality-oriented tags are left as-is.
        let oq = format!("{}{}", "5".repeat(25), "F".repeat(25));
        let e2 = format!("{}{}", "A".repeat(25), "C".repeat(25));

        let mut buf = RecordBuilder::mapped_read()
            .sequence(&"A".repeat(50))
            .cigar("50M")
            .alignment_start(100)
            .mapping_quality(60)
            .tag("OQ", oq.as_str())
            .tag("E2", e2.as_str())
            .build();
        let mut raw = to_raw(&buf);

        RawClipperOnBuf::new(mode).clip_start_of_alignment(&mut buf, 60);
        RawRecordClipper::new(mode).clip_start_of_alignment(&mut raw, 60);

        assert_both_unmapped(&raw, &buf, &format!("fwd-strand unmap mode={mode:?}"));
        for (tag, want) in [(*SamTag::OQ, &oq), (*SamTag::E2, &e2)] {
            let label = String::from_utf8_lossy(&tag).into_owned();
            assert_eq!(typed_string_tag(&buf, tag), want.as_bytes(), "typed {label} unchanged");
            assert_eq!(raw_string_tag(&raw, tag), want.as_bytes(), "raw {label} unchanged");
        }
    }
}
