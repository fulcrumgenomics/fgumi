//! Editing SAM-spec base modification tags (`MM`/`ML`) when SEQ changes.
//!
//! Independent of the consensus callers, so filtering can keep `MM`/`ML` consistent with SEQ
//! whichever callers are built.

use fgumi_raw_bam as bam_fields;
use fgumi_raw_bam::{RawRecordView, SamTag};

/// Rewrites MM/ML after bases were masked to `N`.
///
/// MM skip counts index the bases of each group's type in SEQ, so turning one of them into
/// `N` shifts every later entry. For each group this drops the entries (and their ML values,
/// one per modification code) whose base was masked and recomputes the skips over the bases
/// that remain. `before` and `after` are SEQ before and after masking, both in MM (original
/// read) orientation; masking must only have replaced bases with `N`. Groups on base `N`
/// count every base and are left unchanged. An `ml` of the wrong length for `mm` is an error
/// (an `ml` that is present but empty included); pass `None` for an MM-only tag (such as
/// `am`/`bm`) or an MM without ML, and ignore the returned ML.
///
/// Returns `None` if `mm` is malformed or does not fit `before`.
#[must_use]
pub fn drop_masked_modifications(
    mm: &str,
    ml: Option<&[u8]>,
    before: &[u8],
    after: &[u8],
) -> Option<(String, Vec<u8>)> {
    if before.len() != after.len() {
        return None;
    }
    rewrite_modifications(mm, ml, before, |base, pos| base == b'N' || after[pos] == base, |_| true)
}

/// Drops the MM/ML entries of the bases at `drop` (positions in MM, i.e. original read,
/// orientation) while SEQ is unchanged.
///
/// The dropped bases stay in SEQ, so they still count in the skips of later entries: with the
/// `?` flag an unlisted base has unknown status, which is what a dropped call is. Under the
/// `.` flag (or none) an unlisted base is unmodified, so a call cannot be dropped from such a
/// group without turning it into a confident "unmodified"; fgumi always writes `?`, and such
/// an MM is reported as not editable. `seq` is SEQ in MM orientation. As for
/// [`drop_masked_modifications`], pass `None` for an MM-only tag and ignore the returned ML.
///
/// Returns `None` if `mm` is malformed, does not fit `seq`, or a call must be dropped from a
/// group without the `?` flag.
#[must_use]
pub fn drop_modifications_at(
    mm: &str,
    ml: Option<&[u8]>,
    seq: &[u8],
    drop: &[usize],
) -> Option<(String, Vec<u8>)> {
    let mut dropped = vec![false; seq.len()];
    for &pos in drop {
        if let Some(d) = dropped.get_mut(pos) {
            *d = true;
        }
    }
    rewrite_modifications(mm, ml, seq, |_, _| true, |pos| !dropped[pos])
}

/// Rewrites each MM group, keeping the entry at a listed position when `keep(pos)` and its base
/// is `still_tracked(base, pos)`, and recomputing the skips over the bases still tracked. An
/// entry dropped by `keep` stays tracked (it counts in later skips); one no longer tracked does
/// not. Dropping an entry by `keep` from a group without the `?` flag returns `None`, as its
/// base would then read as unmodified.
fn rewrite_modifications(
    mm: &str,
    ml: Option<&[u8]>,
    before: &[u8],
    still_tracked: impl Fn(u8, usize) -> bool,
    keep: impl Fn(usize) -> bool,
) -> Option<(String, Vec<u8>)> {
    use std::fmt::Write;

    let check_ml = ml.is_some();
    let ml = ml.unwrap_or_default();
    let mut edited_tag = String::with_capacity(mm.len());
    let mut edited_probabilities = Vec::with_capacity(ml.len());
    let mut ml_offset: usize = 0;

    let groups = mm.strip_suffix(';')?;
    for group in groups.split(';') {
        let mut fields = group.split(',');
        let header = fields.next()?;
        let base = header.as_bytes().first()?.to_ascii_uppercase();
        let unlisted_is_unknown = header.ends_with('?');
        let codes = header.get(2..)?.trim_end_matches(['.', '?']);
        let codes_per_base =
            if codes.bytes().all(|b| b.is_ascii_digit()) { 1 } else { codes.len() };
        edited_tag.push_str(header);

        // Positions (in `before`) of the bases of this group's type; masking removes some.
        let tracked: Vec<usize> =
            (0..before.len()).filter(|&i| base == b'N' || before[i] == base).collect();
        let is_tracked = |pos: usize| still_tracked(base, pos);

        let mut ordinal: usize = 0; // index into `tracked` of the next listed entry
        let mut cursor = 0; // entries of `tracked` scanned so far
        let mut pending_skip = 0; // still-tracked bases since the last kept entry
        for skip in fields {
            // A skip is a plain decimal count (`parse` would also accept a leading `+`).
            if skip.is_empty() || !skip.bytes().all(|b| b.is_ascii_digit()) {
                return None;
            }
            ordinal = ordinal.checked_add(skip.parse::<usize>().ok()?)?;
            let &pos = tracked.get(ordinal)?;
            let values = if check_ml {
                let values = ml.get(ml_offset..ml_offset.checked_add(codes_per_base)?)?;
                ml_offset += codes_per_base;
                values
            } else {
                &[]
            };
            ordinal += 1;

            for &p in &tracked[cursor..ordinal - 1] {
                if is_tracked(p) {
                    pending_skip += 1;
                }
            }
            cursor = ordinal;
            if !is_tracked(pos) {
                continue;
            }
            if keep(pos) {
                write!(edited_tag, ",{pending_skip}").expect("write to String is infallible");
                edited_probabilities.extend_from_slice(values);
                pending_skip = 0;
            } else if unlisted_is_unknown {
                pending_skip += 1;
            } else {
                return None;
            }
        }
        edited_tag.push(';');
    }
    if check_ml && ml_offset != ml.len() {
        return None;
    }
    Some((edited_tag, edited_probabilities))
}

/// Keeps `MM`/`ML` and the per-strand `am`/`bm` consistent with SEQ after masking.
///
/// `pre_mask_seq` is SEQ (genomic orientation, as stored) before any masking. The tags index
/// SEQ in original read orientation, so both sequences are reverse-complemented for
/// reverse-mapped records before each tag is edited: the entries of masked bases are dropped
/// and the skips recomputed. A tag that no longer fits SEQ is removed (with `ML`/`MN` for `MM`) rather than
/// left pointing at the wrong bases, as is an `MM` whose `ML` is not `B:C` and so cannot be
/// edited in step. An `MM` without `ML` is edited alone. `MN` (SEQ length) is unaffected by
/// masking; an `MN` that no longer matches SEQ (SEQ was hard-clipped after the tags were
/// written) means no tag fits, and all are removed. Returns whether any tag was removed.
pub fn drop_masked_modifications_raw(record: &mut Vec<u8>, pre_mask_seq: &[u8]) -> bool {
    if !has_modification_tags(record) {
        return false;
    }
    if stale_modification_length(record) {
        return remove_modification_tags(record);
    }
    let view = RawRecordView::new(record);
    let post_mask_seq = view.sequence_vec();
    if post_mask_seq == pre_mask_seq {
        return false;
    }
    let (before, after) = if view.is_reverse() {
        (
            fgumi_dna::dna::reverse_complement(pre_mask_seq),
            fgumi_dna::dna::reverse_complement(&post_mask_seq),
        )
    } else {
        (pre_mask_seq.to_vec(), post_mask_seq)
    };
    edit_modification_tags_raw(record, |tag, ml| {
        drop_masked_modifications(tag, ml, &before, &after)
    })
}

/// Keeps `MM`/`ML`, `MN` and the per-strand `am`/`bm` consistent with SEQ after clipping.
///
/// Clipping can mask bases to `N` (`soft-with-mask`), remove bases from the ends of SEQ (hard
/// clipping), and, when it unmaps a reverse-mapped read, reverse-complement SEQ. The tags index
/// SEQ in original read orientation, so both sequences are first put in that orientation:
/// `pre_clip_seq` (SEQ as stored before clipping) by `pre_clip_reverse` (the record's reverse
/// flag then), and SEQ now by the record's reverse flag now. `removed_start` is the number of
/// bases removed from the start of the stored pre-clip SEQ, or `None` when it is not known, in
/// which case a shorter SEQ cannot be placed.
///
/// SEQ must be a window of the pre-clip SEQ, apart from bases masked to `N`. The calls outside
/// the window or on a masked base are dropped and the skips recomputed over the bases that
/// remain (a group on base `N` counts every base of the window), and `MN`, when present, is set
/// to the new length. All the tags are removed instead when they did not fit SEQ before
/// clipping (`MN` differs from its length) or SEQ is not such a window; a tag that cannot be
/// edited (see [`drop_masked_modifications_raw`]) is removed alone. A record whose SEQ is `*`
/// before and after is left alone. Returns whether any tag was removed.
pub fn trim_clipped_modifications_raw(
    record: &mut Vec<u8>,
    pre_clip_seq: &[u8],
    pre_clip_reverse: bool,
    removed_start: Option<usize>,
) -> bool {
    if !has_modification_tags(record) {
        return false;
    }
    let view = RawRecordView::new(record);
    let post_reverse = view.is_reverse();
    let post_clip_seq = view.sequence_vec();
    if pre_clip_seq.is_empty() && post_clip_seq.is_empty() {
        return false;
    }
    let mn = bam_fields::find_int_tag(bam_fields::aux_data_slice(record), SamTag::MN);
    if mn.is_some_and(|mn| usize::try_from(mn).ok() != Some(pre_clip_seq.len())) {
        return remove_modification_tags(record);
    }
    if pre_clip_reverse == post_reverse && post_clip_seq == pre_clip_seq {
        return false;
    }
    let in_read_orientation = |seq: &[u8], reverse: bool| {
        if reverse { fgumi_dna::dna::reverse_complement(seq) } else { seq.to_vec() }
    };
    let before = in_read_orientation(pre_clip_seq, pre_clip_reverse);
    let after = in_read_orientation(&post_clip_seq, post_reverse);
    // Where the window starts in `before`. Bases removed from the start of a reverse-mapped
    // read's stored SEQ are removed from the end of the read.
    let lo = if after.len() == before.len() {
        Some(0)
    } else {
        removed_start.filter(|_| pre_clip_reverse == post_reverse).and_then(|start| {
            if pre_clip_reverse {
                before.len().checked_sub(start)?.checked_sub(after.len())
            } else {
                Some(start)
            }
        })
    };
    let window = lo.and_then(|lo| Some(lo..lo.checked_add(after.len())?));
    let fits = window.as_ref().and_then(|w| before.get(w.clone())).is_some_and(|kept| {
        kept.iter().zip(&after).all(|(&b, &a)| a == b || a.eq_ignore_ascii_case(&b'N'))
    });
    let Some(window) = window.filter(|_| fits) else {
        return remove_modification_tags(record);
    };
    if after == before {
        return false;
    }
    let removed = edit_modification_tags_raw(record, |tag, ml| {
        rewrite_modifications(
            tag,
            ml,
            &before,
            |base, pos| {
                window.contains(&pos) && (base == b'N' || after[pos - window.start] == base)
            },
            |_| true,
        )
    });
    if after.len() != before.len()
        && bam_fields::find_int_tag(bam_fields::aux_data_slice(record), SamTag::MN).is_some()
    {
        // `MN` is a 32-bit signed tag; a SEQ too long for it cannot be described.
        let Ok(l_seq) = i32::try_from(after.len()) else {
            return remove_modification_tags(record);
        };
        bam_fields::RawTagsEditor::from_vec(record).update_int(SamTag::MN, l_seq);
    }
    removed
}

/// Drops the methylation calls at `positions` (genomic orientation) from `MM`/`ML` and the
/// per-strand `am`/`bm`, leaving SEQ unchanged.
///
/// Used for duplex consensus records, whose SEQ is the molecule's sequence and whose methylation
/// is carried by these tags (see [`mask_methylation_sites_raw`](crate::filter::mask_methylation_sites_raw)). The tags index SEQ in original
/// read orientation, so positions are mirrored for reverse-mapped records. A tag that cannot be
/// edited is removed as in [`drop_masked_modifications_raw`]. Returns whether any tag was
/// removed.
pub fn drop_modifications_at_raw(record: &mut Vec<u8>, positions: &[usize]) -> bool {
    if positions.is_empty() || !has_modification_tags(record) {
        return false;
    }
    if stale_modification_length(record) {
        return remove_modification_tags(record);
    }
    let view = RawRecordView::new(record);
    let seq = view.sequence_vec();
    let (seq, drop) = if view.is_reverse() {
        let last = seq.len().saturating_sub(1);
        let drop: Vec<usize> = positions.iter().map(|&p| last.saturating_sub(p)).collect();
        (fgumi_dna::dna::reverse_complement(&seq), drop)
    } else {
        (seq, positions.to_vec())
    };
    edit_modification_tags_raw(record, |tag, ml| drop_modifications_at(tag, ml, &seq, &drop))
}

/// Whether the record's `MN` (the SEQ length the modification tags were written against)
/// differs from its current SEQ length, as after hard clipping: then no tag fits SEQ.
fn stale_modification_length(record: &[u8]) -> bool {
    let aux = bam_fields::aux_data_slice(record);
    bam_fields::find_int_tag(aux, SamTag::MN)
        .is_some_and(|mn| mn != i64::from(RawRecordView::new(record).l_seq()))
}

/// Removes `MM`, `ML`, `MN`, `am` and `bm`. Returns whether any was present.
fn remove_modification_tags(record: &mut Vec<u8>) -> bool {
    let present = has_modification_tags(record);
    let mut editor = bam_fields::RawTagsEditor::from_vec(record);
    for tag in [SamTag::MM, SamTag::ML, SamTag::MN, SamTag::AM_BASES, SamTag::BM_BASES] {
        editor.remove(tag);
    }
    present
}

/// Whether the record carries any of `MM`, `am` or `bm`.
#[must_use]
pub fn has_modification_tags(record: &[u8]) -> bool {
    let aux = bam_fields::aux_data_slice(record);
    [SamTag::MM, SamTag::AM_BASES, SamTag::BM_BASES]
        .into_iter()
        .any(|tag| bam_fields::find_tag_type(aux, tag).is_some())
}

/// Applies `edit` to `MM` (with its `ML`) and to `am`/`bm`. `edit` receives the tag and the `ML`
/// bytes (empty for `am`/`bm` and for an `MM` without `ML`) and returns the rewritten tag and ML,
/// or `None` when the tag cannot be edited, in which case it is removed (with `ML`/`MN` for
/// `MM`). An `MM` whose `ML` is not `B:C` cannot be edited in step and is removed. `MN` (SEQ
/// length) is unaffected. Returns whether any tag was removed.
fn edit_modification_tags_raw(
    record: &mut Vec<u8>,
    edit: impl Fn(&str, Option<&[u8]>) -> Option<(String, Vec<u8>)>,
) -> bool {
    // `None`: no ML. `Some(None)`: an ML that is not `B:C`. `Some(Some(_))`: the ML bytes.
    let (mm, ml, am, bm) = {
        let aux = bam_fields::aux_data_slice(record);
        let string = |tag| bam_fields::find_string_tag(aux, tag).map(<[u8]>::to_vec);
        let ml = bam_fields::find_tag_type(aux, SamTag::ML).map(|_| {
            bam_fields::find_array_tag(aux, SamTag::ML)
                .filter(|a| a.elem_type == b'C')
                .map(|a| a.data.to_vec())
        });
        (string(SamTag::MM), ml, string(SamTag::AM_BASES), string(SamTag::BM_BASES))
    };
    if mm.is_none() && am.is_none() && bm.is_none() {
        return false;
    }
    let mut removed = false;
    let edit = |tag: &[u8], ml: Option<&[u8]>| edit(std::str::from_utf8(tag).ok()?, ml);

    let mut editor = bam_fields::RawTagsEditor::from_vec(record);
    if let Some(mm) = mm {
        let edited = match &ml {
            None => edit(&mm, None),
            Some(Some(ml)) => edit(&mm, Some(ml)),
            Some(None) => None,
        };
        if let Some((new_mm, new_ml)) = edited {
            editor.update_string(SamTag::MM, new_mm.as_bytes());
            if ml.is_some() {
                editor.update_array_u8(SamTag::ML, &new_ml);
            }
        } else {
            editor.remove(SamTag::MM);
            editor.remove(SamTag::ML);
            editor.remove(SamTag::MN);
            removed = true;
        }
    }
    for (tag, value) in [(SamTag::AM_BASES, am), (SamTag::BM_BASES, bm)] {
        let Some(value) = value else { continue };
        if let Some((new_value, _)) = edit(&value, None) {
            editor.update_string(tag, new_value.as_bytes());
        } else {
            editor.remove(tag);
            removed = true;
        }
    }
    removed
}

#[cfg(test)]
mod tests {
    use super::*;

    /// Masking a base to `N` removes it from the tracked bases MM counts over: its entry (and
    /// ML values) are dropped and the skip counts are recomputed over the bases that remain.
    /// Sequences are in MM (original read) orientation. `CACACC` tracks C at 0, 2, 4, 5; the
    /// single-group MM lists the Cs at 0 and 4.
    #[rstest::rstest]
    #[case::nothing_masked(b"CACACC", "C+m?,0,1;", &[10, 20], "C+m?,0,1;", &[10, 20])]
    #[case::listed_base_masked(b"NACACC", "C+m?,0,1;", &[10, 20], "C+m?,1;", &[20])]
    #[case::unlisted_base_between_masked(b"CANACC", "C+m?,0,1;", &[10, 20], "C+m?,0,0;", &[10, 20])]
    #[case::last_listed_base_masked(b"CACANC", "C+m?,0,1;", &[10, 20], "C+m?,0;", &[10])]
    #[case::trailing_unlisted_base_masked(b"CACACN", "C+m?,0,1;", &[10, 20], "C+m?,0,1;", &[10, 20])]
    #[case::untracked_base_masked(b"CNCACC", "C+m?,0,1;", &[10, 20], "C+m?,0,1;", &[10, 20])]
    #[case::every_listed_base_masked(b"NANANC", "C+m?,0,1;", &[10, 20], "C+m?;", &[])]
    #[case::two_codes_per_base(b"NACACC", "C+mh?,0,1;", &[10, 11, 20, 21], "C+mh?,1;", &[20, 21])]
    #[case::second_group_edited_independently(b"CANACC", "C+m?,0;G-m?,0;", &[10, 30], "C+m?,0;G-m?;", &[10])]
    fn test_drop_masked_modifications(
        #[case] after: &[u8],
        #[case] mm: &str,
        #[case] ml: &[u8],
        #[case] expected_mm: &str,
        #[case] expected_ml: &[u8],
    ) {
        // The G in the two-group case sits where `after` holds N at 2: `CAGACC` before masking.
        let before: &[u8] = if mm.contains('G') { b"CAGACC" } else { b"CACACC" };
        let (new_mm, new_ml) =
            drop_masked_modifications(mm, Some(ml), before, after).expect("well-formed MM");
        assert_eq!(new_mm, expected_mm);
        assert_eq!(new_ml, expected_ml);
    }

    /// Dropping a call keeps its base in SEQ, so it still counts in the skips of later entries.
    /// `CACACC` tracks C at 0, 2, 4, 5; the MM lists the Cs at 0 and 4.
    #[rstest::rstest]
    #[case::nothing_dropped(&[], "C+m?,0,1;", &[10, 20], "C+m?,0,1;", &[10, 20])]
    #[case::first_listed_dropped(&[0], "C+m?,0,1;", &[10, 20], "C+m?,2;", &[20])]
    #[case::last_listed_dropped(&[4], "C+m?,0,1;", &[10, 20], "C+m?,0;", &[10])]
    #[case::unlisted_position_ignored(&[2], "C+m?,0,1;", &[10, 20], "C+m?,0,1;", &[10, 20])]
    #[case::all_dropped(&[0, 4], "C+m?,0,1;", &[10, 20], "C+m?;", &[])]
    #[case::second_group_independent(&[2], "C+m?,0;G-m?,0;", &[10, 30], "C+m?,0;G-m?;", &[10])]
    fn test_drop_modifications_at(
        #[case] drop: &[usize],
        #[case] mm: &str,
        #[case] ml: &[u8],
        #[case] expected_mm: &str,
        #[case] expected_ml: &[u8],
    ) {
        let seq: &[u8] = if mm.contains('G') { b"CAGACC" } else { b"CACACC" };
        let (new_mm, new_ml) =
            drop_modifications_at(mm, Some(ml), seq, drop).expect("well-formed MM");
        assert_eq!(new_mm, expected_mm);
        assert_eq!(new_ml, expected_ml);
    }

    /// An MM whose skips run past the tracked bases, or whose ML is short, cannot be edited.
    #[rstest::rstest]
    #[case::skips_past_tracked_bases("C+m?,0,5;", &[10, 20])]
    #[case::ml_too_short("C+m?,0,1;", &[10])]
    #[case::missing_terminator("C+m?,0,1", &[10, 20])]
    #[case::skip_overflows("C+m?,0,18446744073709551615;", &[10, 20])]
    #[case::signed_skip("C+m?,+0,1;", &[10, 20])]
    #[case::empty_skip("C+m?,,1;", &[10, 20])]
    #[case::present_but_empty_ml("C+m?,0,1;", &[])]
    fn test_drop_masked_modifications_rejects_malformed(#[case] mm: &str, #[case] ml: &[u8]) {
        assert!(drop_masked_modifications(mm, Some(ml), b"CACACC", b"NACACC").is_none());
    }

    /// The parser's other branches: an `N` group counts every base, a `ChEBI` code is one
    /// modification per base, and an MM without ML (or an MM-only tag) is edited alone.
    /// `CACACC` with C at 0, 2, 4, 5; masking the C at 0.
    #[rstest::rstest]
    #[case::n_group_counts_every_base("N+n?,1,2;", Some(&[10u8, 20][..]), "N+n?,1,2;", vec![10, 20])]
    #[case::chebi_code("C+27551?,0,1;", Some(&[10u8, 20][..]), "C+27551?,1;", vec![20])]
    #[case::no_ml("C+m?,0,1;", None, "C+m?,1;", vec![])]
    fn test_drop_masked_modifications_parser_branches(
        #[case] mm: &str,
        #[case] ml: Option<&[u8]>,
        #[case] expected_mm: &str,
        #[case] expected_ml: Vec<u8>,
    ) {
        let (new_mm, new_ml) =
            drop_masked_modifications(mm, ml, b"CACACC", b"NACACC").expect("well-formed MM");
        assert_eq!(new_mm, expected_mm);
        assert_eq!(new_ml, expected_ml);
    }

    /// Under the `.` flag, or with no flag, an unlisted base is unmodified, so a call cannot be
    /// dropped from the group without asserting that; such an MM is not editable. Masking (the
    /// base leaves the group's bases) is still fine.
    #[rstest::rstest]
    #[case::explicit_dot("C+m.,0,1;")]
    #[case::implicit("C+m,0,1;")]
    fn test_drop_modifications_at_needs_unknown_flag(#[case] mm: &str) {
        assert!(drop_modifications_at(mm, Some(&[10, 20]), b"CACACC", &[0]).is_none());
        assert_eq!(
            drop_modifications_at(mm, Some(&[10, 20]), b"CACACC", &[2]),
            Some((mm.to_string(), vec![10, 20])),
            "dropping an unlisted position changes nothing"
        );
        assert!(drop_masked_modifications(mm, Some(&[10, 20]), b"CACACC", b"NACACC").is_some());
    }
}
