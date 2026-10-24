//! Per-end primary selection on re-aligned consensus templates.
//!
//! Input is one template as emitted by the align backend's zipper merge: the
//! primary record of each end, possibly supplementary records, the consensus tags
//! and the `rw` window copied over from the unmapped consensus reads. For each end
//! (R1, R2, or the single fragment read) [`select_primary`] keeps bwa's primary if
//! it overlaps the end's window, else promotes the first `XA` alternative that does,
//! else swaps the primary with an in-window split (`SA`) piece of the end, else
//! leaves the primary as it is. `rw` is stripped from every record.

use anyhow::{Context, Result, anyhow};
use fgumi_raw_bam::{CigarKind, RawRecord, RawRecordView, SamTag, flags};
use fgumi_sam::alignment_tags::regenerate_alignment_tags_raw;
use fgumi_sam::clipper::{ClippingMode, RawRecordClipper};
use fgumi_sam::realign_window::RealignWindow;
use fgumi_sam::xa::{XaHit, parse_xa, write_xa};
use noodles::sam::Header;

use crate::aligner::BwaScoring;
use crate::tag_reversal::reverse_per_base_tags_raw;
use crate::template::{PairOrientation, Template};
use fgumi_sam::ReferenceProvider;

/// Counts of per-end selection outcomes.
#[derive(Debug, Default, Clone, PartialEq, Eq)]
pub struct RealignCounts {
    /// Ends seen (every outcome except [`Selection::Absent`]).
    pub reads: u64,
    /// Ends whose bwa primary was already inside the window.
    pub in_window: u64,
    /// Ends whose primary was replaced by an in-window `XA` alternative.
    pub swapped: u64,
    /// Ends whose primary swapped roles with an in-window supplementary (split,
    /// `SA`) piece of the same end.
    pub swapped_supplementary: u64,
    /// Ends outside their window with an in-window `XA` hit or split piece that
    /// were not swapped because bwa had flagged the pair proper and the mate's
    /// primary is inside the mate's own window.
    pub kept_proper_pair: u64,
    /// Mapped ends with neither the primary, an `XA` hit nor a split piece inside
    /// the window.
    pub outside_window: u64,
    /// Ends that did not map.
    pub unmapped: u64,
    /// Ends without an `rw` window.
    pub no_window: u64,
    /// Ends whose `XA` tag did not parse; they are treated as having no `XA`.
    pub malformed_xa: u64,
    /// `XA` hits skipped because their CIGAR is not a valid mapped CIGAR.
    pub skipped_xa_hits: u64,
}

impl RealignCounts {
    /// Adds `other`'s counts to these. Destructures `other` so a new field must be
    /// added here before it compiles.
    pub(crate) fn merge(&mut self, other: &Self) {
        let Self {
            reads,
            in_window,
            swapped,
            swapped_supplementary,
            kept_proper_pair,
            outside_window,
            unmapped,
            no_window,
            malformed_xa,
            skipped_xa_hits,
        } = other;
        self.reads += reads;
        self.in_window += in_window;
        self.swapped += swapped;
        self.swapped_supplementary += swapped_supplementary;
        self.kept_proper_pair += kept_proper_pair;
        self.outside_window += outside_window;
        self.unmapped += unmapped;
        self.no_window += no_window;
        self.malformed_xa += malformed_xa;
        self.skipped_xa_hits += skipped_xa_hits;
    }

    /// Adds the outcomes of one template's ends.
    pub(crate) fn tally(&mut self, outcome: &SelectOutcome) {
        self.malformed_xa += u64::from(outcome.malformed_xa);
        self.skipped_xa_hits += u64::from(outcome.skipped_xa_hits);
        for &sel in &outcome.ends {
            let counter = match sel {
                Selection::Absent => continue,
                Selection::Kept => &mut self.in_window,
                Selection::Swapped => &mut self.swapped,
                Selection::SwappedSupplementary => &mut self.swapped_supplementary,
                Selection::KeptProperPair => &mut self.kept_proper_pair,
                Selection::OutsideWindow => &mut self.outside_window,
                Selection::Unmapped => &mut self.unmapped,
                Selection::NoWindow => &mut self.no_window,
            };
            *counter += 1;
            self.reads += 1;
        }
    }
}

/// What [`select_primary`] did to one template.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub(crate) struct SelectOutcome {
    /// The outcome for `[R1, R2]` (a fragment's read is R1).
    pub(crate) ends: [Selection; 2],
    /// Ends whose `XA` tag did not parse (treated as having no `XA`).
    pub(crate) malformed_xa: u32,
    /// `XA` hits skipped for a CIGAR that is not a valid mapped CIGAR.
    pub(crate) skipped_xa_hits: u32,
}

/// The outcome of primary selection for one end.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub(crate) enum Selection {
    /// The primary overlaps the window and is kept.
    Kept,
    /// An in-window `XA` alternative replaced the primary.
    Swapped,
    /// No `XA` hit overlaps the window, but a split (`SA`) piece of the end
    /// does: the piece and the primary swapped roles.
    SwappedSupplementary,
    /// The primary is outside the window and an `XA` hit or split piece is
    /// inside it, but bwa had flagged the pair proper and the mate's primary is
    /// inside the mate's own window, so the primary is kept (the proper-pair
    /// guard).
    KeptProperPair,
    /// Neither the primary, an `XA` hit nor a split piece overlaps the window;
    /// the primary is kept.
    OutsideWindow,
    /// The end is unmapped.
    Unmapped,
    /// The end has no `rw` window.
    NoWindow,
    /// The template has no record for this end (the second end of a fragment).
    Absent,
}

/// Chooses the primary alignment of each end of `merged` using its `rw` window.
///
/// Returns the outcome for `[R1, R2]`; a fragment's read is reported as R1 and its
/// second slot is [`Selection::Absent`]. A malformed `XA` never fails the run: an
/// `XA` tag that does not parse is treated as absent, and a hit whose CIGAR is not
/// a valid mapped CIGAR is skipped; both are counted in the returned
/// [`SelectOutcome`].
///
/// Per end, in this order:
///
/// 1. bwa's primary overlaps the window: it is kept ([`Selection::Kept`]).
/// 2. else the first `XA` hit that overlaps the window (and whose CIGAR covers
///    the read) is promoted ([`Selection::Swapped`]; see below).
/// 3. else, if a split piece of the end overlaps the window (on its contig and
///    strand, by the same test as an `XA` hit), the piece and the primary swap
///    roles ([`Selection::SwappedSupplementary`]; see below). A piece is a
///    supplementary record of the end, or a secondary one carrying `SA` (an
///    aligner run with `-M`, as `bwa-mem3 mem --meth` is, flags the pieces of a
///    split secondary). Of several in-window pieces the one with the longest
///    reference span is taken, then the first.
/// 4. else the primary is kept ([`Selection::OutsideWindow`]).
///
/// The proper-pair guard: steps 2 and 3 are skipped, and the primary kept
/// ([`Selection::KeptProperPair`]), when bwa flagged the pair proper (`0x2` on
/// both primaries) and the mate's primary overlaps the mate's own window. The
/// mate then agrees with bwa's placement of this end, so this end's window is
/// the suspect: for example, the raw reads were placed on a poly-G artefact
/// that the consensus read, with a cleaner tail, no longer prefers. The guard is
/// decided from both ends' state before either is changed, so it does not
/// depend on which end is processed first; the mate's window is read only for
/// an end that has a candidate to swap in. Fragments, and pairs whose mate is
/// unmapped, outside its window or without one, are never guarded.
///
/// An `XA` promotion:
///
/// - moves the primary to the hit with MAPQ 0. Its `AS` is the score bwa lists
///   for the hit in a six-field `XA` (`bwa-mem3 mem -u`). With the four-field
///   `XA` it is recomputed with `scoring` from its CIGAR and `XA`'s edit distance
///   (mismatches = `NM` minus indel bases), so bisulfite conversions bwa did not
///   count are not penalised. That approximates bwa's score: bwa scores a base
///   against an `N` (in the read or the reference) as -1, while here every
///   mismatch counted in `NM` costs `-B`, and bwa's local-alignment clipping is
///   taken as given by the `XA` CIGAR. `AS` and `XS` are stored in the smallest
///   signed integer type, as the zipper merge stores them.
/// - removes the end's split records: its supplementary records, and its
///   secondary records that carry `SA`. Secondary records without `SA` (`-a`
///   alternatives) are kept.
///
/// A split-piece swap:
///
/// - makes the piece the primary: its supplementary/secondary flag is cleared,
///   hard clips become soft clips, and it gets the full SEQ/QUAL from the old
///   primary in its own strand's orientation. It keeps the aux tags bwa and the
///   zipper merge gave it, which include its own MAPQ, `NM`, `MD`, `AS` and `XS`
///   (the MAPQ is not zeroed: bwa scored this piece's placement itself) and the
///   consensus tags, which the merge copies onto every record of the end with
///   the per-base ones oriented for that record's strand. `XA` is removed (on a
///   piece it described that piece's alternatives, not the read's), and so are
///   `tc` and `pa`. The bisulfite tags bwa computed for the piece's own
///   placement stay: `XG` and `YD`, and `XM` when it has one call per base of
///   the full read (a hard-clipped piece's `XM` covers only its own bases, so
///   it is removed).
/// - makes the old primary the piece, flagged as the piece was (supplementary,
///   or secondary under `-M`). If the piece was hard-clipped, so is it: its SEQ
///   and QUAL are trimmed to the soft-clipped CIGAR's aligned bases, and so is
///   its `XM`. Its per-base consensus tags stay full-length, as the merge leaves
///   them on a hard-clipped supplementary. It keeps its own MAPQ and aux tags.
/// - rewrites `SA` on every split piece of the end, the new primary first, in
///   bwa's format `rname,pos,strand,CIGAR,mapQ,NM;` (hard clips written as soft
///   clips, as bwa writes them). The end's other pieces stay as they are.
///
/// After any swap:
///
/// - mate fields are recomputed with [`Template::fix_mate_info`], and the
///   remaining secondary records are pointed at the other end's moved or kept
///   primary as well; a lone mate (its other end absent) gets
///   [`Template::set_lone_mate_info`] instead.
/// - the proper-pair flag is decided once for the template and applied to
///   every record of it, split pieces and secondaries included. It is kept only
///   if bwa had flagged the pair proper and the pair is still FR on one contig
///   with a `|TLEN|` no larger than before the swap, and cleared otherwise. The
///   windows say where each end belongs, not that the molecule's raw reads were
///   a concordant pair (nothing upstream filters discordant molecules), so they
///   cannot stand in for bwa's insert-size check.
/// - the `tc` (template-coordinate) tags of secondary and supplementary records,
///   which the merge computed from the old primaries, are recomputed (and a
///   promoted piece's own `tc` removed).
///
/// `rw` is removed from every record whatever the outcome.
///
/// # Errors
///
/// Returns an error if an `rw` value is malformed, a record names a contig
/// missing from `header`, or the promoted alignment cannot be scored against
/// `reference`.
pub(crate) fn select_primary(
    merged: &mut Template,
    reference: &impl ReferenceProvider,
    header: &Header,
    scoring: &BwaScoring,
) -> Result<SelectOutcome> {
    let mut outcome =
        SelectOutcome { ends: [Selection::Absent; 2], malformed_xa: 0, skipped_xa_hits: 0 };
    let before = PairState::of(merged);
    // The guard needs both ends as bwa left them. R1's looks at R2's window,
    // which is unchanged until R2 is selected. R2's looks at R1's, which may
    // have changed by then; but R1 was in its window, before any change,
    // exactly when R1 was kept.
    let proper = PairState::both_proper(merged);
    if let Some((idx, _)) = merged.r1 {
        let guard = match merged.r2 {
            Some((mate, _)) if proper => Guard::MateInWindow(mate),
            _ => Guard::Off,
        };
        outcome.ends[0] =
            select_end(merged, idx, true, guard, reference, header, scoring, &mut outcome)?;
    }
    if let Some((idx, _)) = merged.r2 {
        let guard =
            if proper && outcome.ends[0] == Selection::Kept { Guard::On } else { Guard::Off };
        outcome.ends[1] =
            select_end(merged, idx, false, guard, reference, header, scoring, &mut outcome)?;
    }
    let out = outcome.ends;

    if out.iter().any(|s| matches!(s, Selection::Swapped | Selection::SwappedSupplementary)) {
        // An XA-swapped end's split records described the old primary's split.
        for (slot, first) in [(0, true), (1, false)] {
            if out[slot] == Selection::Swapped {
                merged.remove_split_records(first)?;
            }
        }
        if merged.lone_mate_segment().is_some() {
            merged.set_lone_mate_info();
        } else {
            merged.fix_mate_info()?;
            merged.fix_secondary_mate_info();
            set_proper_pair(merged, before);
        }
        merged.add_template_coordinate_tags();
    }

    for record in &mut merged.records {
        fgumi_raw_bam::remove_tag(record.as_mut_vec(), SamTag::RW);
    }
    // Tag edits shifted aux bytes, so a cached UMI offset would be stale.
    merged.cached_umi_position = None;
    Ok(outcome)
}

/// The proper-pair guard of one end (see [`select_primary`]), resolved only
/// when the end has a candidate to swap in.
#[derive(Debug, Clone, Copy)]
enum Guard {
    /// Not guarded: a fragment or lone mate, or bwa did not flag the pair proper.
    Off,
    /// Guarded: the pair is proper and the mate's primary was in its window.
    On,
    /// The pair is proper: guarded if the mate's primary, the unchanged record
    /// at this index, overlaps its window.
    MateInWindow(usize),
}

impl Guard {
    /// Whether the end is guarded, reading the mate's window if needed.
    fn holds(self, merged: &Template, header: &Header) -> Result<bool> {
        match self {
            Self::Off => Ok(false),
            Self::On => Ok(true),
            Self::MateInWindow(mate) => primary_in_window(&merged.records[mate], header),
        }
    }
}

/// Whether the mapped `record` overlaps its `rw` window; `false` when it is
/// unmapped or has no window.
fn primary_in_window(record: &RawRecord, header: &Header) -> Result<bool> {
    let Some(rw) =
        fgumi_raw_bam::find_string_tag(fgumi_raw_bam::aux_data_slice(record), SamTag::RW)
    else {
        return Ok(false);
    };
    if record.flags() & flags::UNMAPPED != 0 {
        return Ok(false);
    }
    let window = RealignWindow::parse(rw).with_context(|| invalid_tag(record, "rw"))?;
    let (contig, start, end, reverse) = placement(record, header)?;
    Ok(window.overlaps(contig, start, end, reverse))
}

/// Selects the primary of one end (`first`: R1 or a fragment's read), whose
/// primary record is `merged.records[idx]`: keeps it, promotes an in-window `XA`
/// hit, or swaps it with an in-window split piece. When `guard` holds (the
/// proper-pair guard of [`select_primary`]) an end that would be swapped is
/// kept instead. A malformed `XA` tag or hit CIGAR is counted in `outcome`, not
/// an error.
#[allow(clippy::too_many_arguments)]
fn select_end(
    merged: &mut Template,
    idx: usize,
    first: bool,
    guard: Guard,
    reference: &impl ReferenceProvider,
    header: &Header,
    scoring: &BwaScoring,
    outcome: &mut SelectOutcome,
) -> Result<Selection> {
    // The candidates are found while the window borrows the record; the record
    // is changed only after its last use.
    let record = &merged.records[idx];
    let Some(rw) =
        fgumi_raw_bam::find_string_tag(fgumi_raw_bam::aux_data_slice(record), SamTag::RW)
    else {
        return Ok(Selection::NoWindow);
    };
    let window = RealignWindow::parse(rw).with_context(|| invalid_tag(record, "rw"))?;

    if RawRecordView::new(record).flags() & flags::UNMAPPED != 0 {
        return Ok(Selection::Unmapped);
    }
    let (contig, start, end, reverse) = placement(record, header)?;
    if window.overlaps(contig, start, end, reverse) {
        return Ok(Selection::Kept);
    }

    let candidate = match in_window_xa_hit(record, &window, outcome) {
        Some((hits, i, ops)) => Candidate::XaHit(hits, i, ops),
        None => match in_window_split_piece(merged, idx, first, &window, header)? {
            Some(piece) => Candidate::SplitPiece(piece),
            None => return Ok(Selection::OutsideWindow),
        },
    };
    if guard.holds(merged, header)? {
        return Ok(Selection::KeptProperPair);
    }
    match candidate {
        Candidate::XaHit(hits, i, ops) => {
            promote(&mut merged.records[idx], &hits, i, &ops, reference, header, scoring)?;
            Ok(Selection::Swapped)
        }
        Candidate::SplitPiece(piece) => {
            swap_with_split_piece(merged, idx, first, piece, reference, header)?;
            Ok(Selection::SwappedSupplementary)
        }
    }
}

/// What an end outside its window could be swapped to.
enum Candidate {
    /// An in-window `XA` hit: the parsed hits, the hit's index and its CIGAR.
    XaHit(Vec<XaHit>, usize, Vec<u32>),
    /// The index of an in-window split piece of the end.
    SplitPiece(usize),
}

/// `(contig, 1-based start, 1-based inclusive end, reverse)` of a mapped record.
fn placement<'h>(record: &[u8], header: &'h Header) -> Result<(&'h str, u32, u32, bool)> {
    let view = RawRecordView::new(record);
    let contig = contig_name(header, view.ref_id())?;
    let start = u32::try_from(view.pos() + 1).context("record has a negative position")?;
    let span = u32::try_from(fgumi_raw_bam::reference_length_from_raw_bam(record))
        .context("record has a negative reference span")?;
    let end = start.saturating_add(span).saturating_sub(1).max(start);
    Ok((contig, start, end, view.flags() & flags::REVERSE != 0))
}

/// The first in-window `XA` hit of `record` whose CIGAR covers the read, as
/// `(hits, index of the hit, its CIGAR ops)`, or `None` without one. A
/// malformed `XA` tag or hit CIGAR is counted in `outcome` and treated as
/// absent.
fn in_window_xa_hit(
    record: &RawRecord,
    window: &RealignWindow<'_>,
    outcome: &mut SelectOutcome,
) -> Option<(Vec<XaHit>, usize, Vec<u32>)> {
    let xa = fgumi_raw_bam::find_string_tag(fgumi_raw_bam::aux_data_slice(record), SamTag::XA)?;
    let hits = match parse_xa(xa) {
        Ok(hits) => hits,
        Err(e) => {
            log::debug!("{}: {e:#}; treating it as absent", invalid_tag(record, "XA"));
            outcome.malformed_xa += 1;
            return None;
        }
    };
    let l_seq = fgumi_raw_bam::l_seq(record) as usize;
    let mut candidate = None;
    for (i, hit) in hits.iter().enumerate() {
        let Some(ops) = hit.cigar_ops() else {
            log::debug!("{}: skipping hit with CIGAR '{}'", invalid_tag(record, "XA"), hit.cigar);
            outcome.skipped_xa_hits += 1;
            continue;
        };
        // A hit whose CIGAR does not cover the read's bases cannot replace the primary.
        if window.overlaps(&hit.contig, hit.start, hit.end(&ops), hit.reverse)
            && fgumi_raw_bam::query_length_from_cigar(&ops) == l_seq
        {
            candidate = Some((i, ops));
            break;
        }
    }
    candidate.map(|(i, ops)| (hits, i, ops))
}

/// Replaces the primary alignment of `record` with `hits[i]` (whose CIGAR is `ops`).
fn promote(
    record: &mut RawRecord,
    hits: &[XaHit],
    i: usize,
    ops: &[u32],
    reference: &impl ReferenceProvider,
    header: &Header,
    scoring: &BwaScoring,
) -> Result<()> {
    let hit = &hits[i];
    let new_ref_id = header
        .reference_sequences()
        .get_index_of(hit.contig.as_bytes())
        .ok_or_else(|| anyhow!("XA contig '{}' is not in the header", hit.contig))?;
    let new_ref_id = i32::try_from(new_ref_id).context("reference id overflows i32")?;
    let new_pos = i32::try_from(hit.start).context("XA position overflows i32")? - 1;

    // The old primary, which becomes an XA alternative.
    let mut old = old_primary_hit(record, header, reference)?;
    let aux = fgumi_raw_bam::aux_data_slice(record);
    let old_score = fgumi_raw_bam::find_int_tag(aux, SamTag::AS);
    // In the six-field `XA` of `bwa-mem3 mem -u`, the old primary lists its own
    // score and MAPQ like the other hits.
    if hits.iter().any(|h| h.score_mapq.is_some()) {
        let score = old_score.and_then(|v| i32::try_from(v).ok()).unwrap_or(0);
        old.score_mapq = Some((score, u32::from(RawRecordView::new(record).mapq())));
    }
    let had_uq = fgumi_raw_bam::find_int_tag(aux, SamTag::UQ).is_some();

    if hit.reverse != old.reverse {
        flip_strand(record)?;
    }
    record.set_ref_id(new_ref_id);
    record.set_pos(new_pos);
    record.set_cigar_ops(ops);
    record.set_mapq(0);
    record.recompute_bin();
    let record = record.as_mut_vec();

    // MD from the reference (under bisulfite alignment it lists conversions, as
    // bwa-mem3 --meth's MD does). NM and the mismatch count AS is scored with
    // come from the XA hit: bwa's edit distance, which under --meth does not
    // count conversions.
    regenerate_alignment_tags_raw(record, header, reference).with_context(|| {
        format!("regenerating MD of the promoted alignment of {}", read_name(record))
    })?;
    let indel: i64 =
        ops.iter().filter(|&&op| matches!(op & 0xF, 1 | 2)).map(|&op| i64::from(op >> 4)).sum();
    let mismatches = u32::try_from((i64::from(hit.nm) - indel).max(0)).unwrap_or(u32::MAX);
    if !had_uq {
        fgumi_raw_bam::remove_tag(record, SamTag::UQ);
    }
    fgumi_raw_bam::update_int_tag(record, SamTag::NM, i32::try_from(hit.nm).unwrap_or(i32::MAX));
    let score = hit.score_mapq.map_or_else(|| scoring.score(ops, mismatches), |(score, _)| score);
    fgumi_raw_bam::update_int_tag(record, SamTag::AS, score);
    match old_score.and_then(|v| i32::try_from(v).ok()) {
        Some(v) => fgumi_raw_bam::update_int_tag(record, SamTag::XS, v),
        None => fgumi_raw_bam::remove_tag(record, SamTag::XS),
    }
    // In the smallest signed type, as the zipper merge stores every other `AS`/`XS`.
    fgumi_raw_bam::normalize_int_tag_to_smallest_signed(record, SamTag::AS);
    fgumi_raw_bam::normalize_int_tag_to_smallest_signed(record, SamTag::XS);

    let mut alternatives: Vec<XaHit> =
        hits.iter().enumerate().filter(|&(j, _)| j != i).map(|(_, h)| h.clone()).collect();
    alternatives.push(old);
    let mut xa = Vec::new();
    write_xa(&alternatives, &mut xa);
    fgumi_raw_bam::update_string_tag(record, SamTag::XA, &xa);
    // Tags the aligner computed for the old placement, which the XA hit carries
    // nothing to recompute: the split (`SA`) and the per-alignment tags of
    // [`PLACEMENT_TAGS`].
    fgumi_raw_bam::remove_tag(record, SamTag::SA);
    for tag in PLACEMENT_TAGS {
        fgumi_raw_bam::remove_tag(record, tag);
    }
    Ok(())
}

/// Per-alignment tags removed from a record promoted to primary: bwa-mem3's
/// score ratio (`pa`) and the bisulfite tags of `bwa-mem3 mem --meth` / bwameth
/// that depend on the aligned strand or locus (`XG`, `XM`, `YD`). `XR` (the
/// read's conversion) describes the read, not its alignment, and is kept.
const PLACEMENT_TAGS: [SamTag; 4] = [SamTag::PA, SamTag::XG, SamTag::XM, SamTag::YD];

/// The edit distance of `record`: its `NM`, or recomputed against the reference
/// if it has none.
fn edit_distance(
    record: &[u8],
    header: &Header,
    reference: &impl ReferenceProvider,
) -> Result<i64> {
    if let Some(nm) = fgumi_raw_bam::find_int_tag(fgumi_raw_bam::aux_data_slice(record), SamTag::NM)
    {
        return Ok(nm);
    }
    let mut copy = record.to_vec();
    regenerate_alignment_tags_raw(&mut copy, header, reference)?;
    Ok(fgumi_raw_bam::find_int_tag(fgumi_raw_bam::aux_data_slice(&copy), SamTag::NM).unwrap_or(0))
}

/// The current primary of `record` as an `XA` hit. `NM` is taken from the record,
/// or recomputed against the reference if the record has none.
fn old_primary_hit(
    record: &[u8],
    header: &Header,
    reference: &impl ReferenceProvider,
) -> Result<XaHit> {
    let view = RawRecordView::new(record);
    let nm = edit_distance(record, header, reference)?;
    Ok(XaHit {
        contig: contig_name(header, view.ref_id())?.to_string(),
        start: u32::try_from(view.pos() + 1).context("primary has a negative position")?,
        reverse: view.flags() & flags::REVERSE != 0,
        cigar: fgumi_raw_bam::cigar_to_string_from_raw(record),
        nm: u32::try_from(nm).unwrap_or(0),
        score_mapq: None,
    })
}

/// `ops` with every hard clip turned into a soft clip of the same length.
fn hard_to_soft(ops: &[u32]) -> Vec<u32> {
    let (hard, soft) = (CigarKind::HardClip as u32, CigarKind::SoftClip as u32);
    ops.iter().map(|&op| if op & 0xF == hard { (op & !0xF) | soft } else { op }).collect()
}

/// The split piece of one end to swap with its out-of-window primary
/// (`merged.records[primary]`): of the pieces that overlap `window` and whose
/// CIGAR, with hard clips as soft clips, covers the read, the one with the
/// longest reference span, then the first.
fn in_window_split_piece(
    merged: &Template,
    primary: usize,
    first: bool,
    window: &RealignWindow<'_>,
    header: &Header,
) -> Result<Option<usize>> {
    let l_seq = fgumi_raw_bam::l_seq(&merged.records[primary]) as usize;
    let mut best: Option<(usize, u32)> = None;
    for i in merged.split_piece_indices(first) {
        let record = &merged.records[i];
        if record.flags() & flags::UNMAPPED != 0 {
            continue;
        }
        let (contig, start, end, reverse) = placement(record, header)?;
        let ops = hard_to_soft(&record.cigar_ops_vec());
        if !window.overlaps(contig, start, end, reverse)
            || fgumi_raw_bam::query_length_from_cigar(&ops) != l_seq
        {
            continue;
        }
        let span = end - start + 1;
        if best.is_none_or(|(_, longest)| span > longest) {
            best = Some((i, span));
        }
    }
    Ok(best.map(|(i, _)| i))
}

/// Swaps the roles of one end's primary (`merged.records[primary]`) and its split
/// piece `merged.records[piece]`; see [`select_primary`].
fn swap_with_split_piece(
    merged: &mut Template,
    primary: usize,
    first: bool,
    piece: usize,
    reference: &impl ReferenceProvider,
    header: &Header,
) -> Result<()> {
    let old = &merged.records[primary];
    let old_reverse = old.flags() & flags::REVERSE != 0;
    let piece_reverse = merged.records[piece].flags() & flags::REVERSE != 0;
    let hard_clipped = merged.records[piece]
        .cigar_ops_vec()
        .iter()
        .any(|&op| op & 0xF == CigarKind::HardClip as u32);

    // The full read in the piece's orientation.
    let mut seq = fgumi_raw_bam::extract_sequence(old);
    let mut qual = fgumi_raw_bam::quality_scores_slice(old).to_vec();
    if old_reverse != piece_reverse {
        seq = fgumi_dna::dna::reverse_complement(&seq);
        qual.reverse();
    }

    // Swaps the two records and their flags (see the Template method).
    merged.swap_primary_with_split_piece(first, piece)?;

    // The piece, now in the primary slot, becomes the primary.
    let promoted = &mut merged.records[primary];
    let ops = hard_to_soft(&promoted.cigar_ops_vec());
    promoted.set_cigar_ops(&ops);
    promoted.set_sequence_and_qualities(&seq, &qual);
    let aux = promoted.as_mut_vec();
    for tag in [SamTag::XA, SamTag::TC, SamTag::PA] {
        fgumi_raw_bam::remove_tag(aux, tag);
    }
    let covers_read =
        fgumi_raw_bam::find_string_tag(fgumi_raw_bam::aux_data_slice(aux), SamTag::XM)
            .is_none_or(|xm| xm.len() == seq.len());
    if !covers_read {
        fgumi_raw_bam::remove_tag(aux, SamTag::XM);
    }

    // The old primary, now in the piece's slot, becomes the piece, clipped as
    // the piece was. The default clipper leaves aux tags alone: the per-base
    // consensus tags stay full-length, as the merge leaves them on a
    // hard-clipped supplementary.
    if hard_clipped {
        let demoted = &mut merged.records[piece];
        let (lead, trail) =
            RawRecordClipper::new(ClippingMode::Hard).upgrade_all_clipping_raw(demoted)?;
        trim_xm(demoted, lead, trail);
    }

    rewrite_sa(merged, first, reference, header)
}

/// Trims `record`'s `XM` (one Bismark methylation call per base of SEQ) after
/// `lead` and `trail` bases were clipped off the ends of SEQ, so it matches
/// SEQ again; an `XM` of any other length than SEQ had is removed.
fn trim_xm(record: &mut RawRecord, lead: usize, trail: usize) {
    let l_seq = fgumi_raw_bam::l_seq(record) as usize;
    let Some(xm) =
        fgumi_raw_bam::find_string_tag(fgumi_raw_bam::aux_data_slice(record), SamTag::XM)
    else {
        return;
    };
    if xm.len() == lead + l_seq + trail {
        let trimmed = xm[lead..lead + l_seq].to_vec();
        fgumi_raw_bam::update_string_tag(record.as_mut_vec(), SamTag::XM, &trimmed);
    } else {
        fgumi_raw_bam::remove_tag(record.as_mut_vec(), SamTag::XM);
    }
}

/// Rewrites `SA` on every split piece of one end (its primary and
/// [`Template::split_piece_indices`]): each lists all the others, the primary
/// first, as `rname,pos,strand,CIGAR,mapQ,NM;` with hard clips written as soft
/// clips (bwa's format). A piece without others loses `SA`.
fn rewrite_sa(
    merged: &mut Template,
    first: bool,
    reference: &impl ReferenceProvider,
    header: &Header,
) -> Result<()> {
    let primary = if first { merged.r1 } else { merged.r2 };
    let Some((primary, _)) = primary else {
        return Ok(());
    };
    let pieces: Vec<usize> =
        std::iter::once(primary).chain(merged.split_piece_indices(first)).collect();
    let mut entries = Vec::with_capacity(pieces.len());
    for &i in &pieces {
        let record = &merged.records[i];
        let (contig, start, _, reverse) = placement(record, header)?;
        let cigar = record.cigar_to_string().replace('H', "S");
        let nm = edit_distance(record, header, reference)?;
        let strand = if reverse { '-' } else { '+' };
        entries.push(format!("{contig},{start},{strand},{cigar},{},{nm};", record.mapq()));
    }
    for (k, &i) in pieces.iter().enumerate() {
        let sa: String =
            entries.iter().enumerate().filter(|&(j, _)| j != k).map(|(_, e)| e.as_str()).collect();
        let aux = merged.records[i].as_mut_vec();
        if sa.is_empty() {
            fgumi_raw_bam::remove_tag(aux, SamTag::SA);
        } else {
            fgumi_raw_bam::update_string_tag(aux, SamTag::SA, sa.as_bytes());
        }
    }
    Ok(())
}

/// Moves `record` to the opposite strand: reverse-complements SEQ, reverses QUAL,
/// re-orients the per-base consensus tags and toggles the reverse flag.
fn flip_strand(record: &mut RawRecord) -> Result<()> {
    // `reverse_per_base_tags_raw` acts only on a reverse-strand record, so call it
    // while the record is on the reverse strand: before the flip when leaving it,
    // after when entering it. Either way the tags are reversed exactly once.
    let was_reverse = RawRecordView::new(record).flags() & flags::REVERSE != 0;
    if was_reverse {
        reverse_per_base_tags_raw(record)?;
    }
    let seq = fgumi_dna::dna::reverse_complement(&fgumi_raw_bam::extract_sequence(record));
    let mut qual = fgumi_raw_bam::quality_scores_slice(record).to_vec();
    qual.reverse();
    record.set_sequence_and_qualities(&seq, &qual);
    record.set_flags(record.flags() ^ flags::REVERSE);
    if !was_reverse {
        reverse_per_base_tags_raw(record)?;
    }
    Ok(())
}

/// A pair's proper-pair flag and `|TLEN|` as bwa left them, before any swap.
#[derive(Debug, Clone, Copy)]
struct PairState {
    proper: bool,
    tlen: u32,
}

impl PairState {
    /// Whether bwa flagged the pair proper: `0x2` on both primaries (`false`
    /// without both).
    fn both_proper(template: &Template) -> bool {
        match (template.r1, template.r2) {
            (Some((r1, _)), Some((r2, _))) => {
                template.records[r1].flags() & template.records[r2].flags() & flags::PROPER_PAIR
                    != 0
            }
            _ => false,
        }
    }

    /// The state of `template`'s R1 primary, or `None` without both ends.
    fn of(template: &Template) -> Option<Self> {
        let ((r1, _), _) = template.r1.zip(template.r2)?;
        let record = &template.records[r1];
        Some(Self {
            proper: record.flags() & flags::PROPER_PAIR != 0,
            tlen: fgumi_raw_bam::template_length(record).unsigned_abs(),
        })
    }
}

/// Decides a pair's proper-pair flag after a swap and applies it to every
/// record of the template (primaries, split pieces and secondaries alike, as
/// bwa flags every record of a proper pair). It is kept only if bwa had flagged
/// the pair proper (`before`) and the pair is still FR on one contig with a
/// `|TLEN|` no larger than before. Otherwise it is cleared: bwa set it for a
/// placement the swap replaced, and the windows cannot stand in for bwa's
/// insert-size check (they come from every mapped raw read of the molecule,
/// discordant or not).
fn set_proper_pair(template: &mut Template, before: Option<PairState>) {
    let (Some(before), Some((r1, _))) = (before, template.r1) else {
        return;
    };
    let fr = template.pair_orientation() == Some(PairOrientation::FR);
    let tlen = fgumi_raw_bam::template_length(&template.records[r1]).unsigned_abs();
    let proper = before.proper && fr && tlen <= before.tlen;
    for record in &mut template.records {
        let flag = record.flags();
        let flag = if proper { flag | flags::PROPER_PAIR } else { flag & !flags::PROPER_PAIR };
        record.set_flags(flag);
    }
}

/// The name of reference `ref_id` in `header`.
fn contig_name(header: &Header, ref_id: i32) -> Result<&str> {
    let (name, _) = usize::try_from(ref_id)
        .ok()
        .and_then(|i| header.reference_sequences().get_index(i))
        .ok_or_else(|| anyhow!("reference id {ref_id} is not in the header"))?;
    std::str::from_utf8(name.as_ref()).context("reference name is not UTF-8")
}

fn read_name(record: &[u8]) -> String {
    String::from_utf8_lossy(fgumi_raw_bam::read_name(record)).into_owned()
}

fn invalid_tag(record: &[u8], tag: &str) -> String {
    format!("invalid {tag} tag on read {}", read_name(record))
}

#[cfg(test)]
#[allow(clippy::too_many_arguments)] // rstest case tables and their `expect` helper
mod tests {
    use super::*;
    use crate::reference::ReferenceReader;
    use fgumi_raw_bam::SamBuilder;
    use fgumi_raw_bam::flags::{
        FIRST_SEGMENT, LAST_SEGMENT, MATE_UNMAPPED, PAIRED, PROPER_PAIR, REVERSE, SECONDARY,
        SUPPLEMENTARY, UNMAPPED,
    };
    use fgumi_sam::builder::create_test_fasta;
    use noodles::sam::header::record::value::{Map, map::ReferenceSequence};
    use rstest::{fixture, rstest};
    use std::num::NonZeroUsize;
    use tempfile::NamedTempFile;

    const CHR1: &str = "GCTAAAGACAATTACATAACATACACGTCAGCACGAAACTTGTTGGCCCAGTGTGAATCGCTTAAGGGTTAAGTAAGTGTGATGCATACGCCTTTACTTGCTGTGTCCACCCCATCGGACTGGCATTTTTATTACACTCAGAAACAGAACTCGGGTAATTTTGACAGGTCACGCAGAGGCGCGCCCTCCTGAAGTGCGTG";
    const CHR2: &str = "GACACTCGCTATGAATCTCTGATTTACCCACTCTGCCAAACTCCAGCGCGGTCAGTTCCATCACCCTAAGTAACCGAATAATGCGTTCGCTCTATTGACTACGACGCGCTCATTCCCTTGTCGGAGAGTTATGGAACAAGGACGCTGTCTGAGACTAGAAGACAGATAGTGCACACGACCGGCGTCGGAGAAACTCTATT";

    /// The window both R1 (and a fragment) and R2 truly come from.
    const RW1: &str = "chr1:91-130:+";
    const RW2: &str = "chr1:150-190:-";
    /// An R2 window its `chr1:161-180:-` placement does not overlap.
    const RW2_ELSEWHERE: &str = "chr1:185-199:-";

    struct Fixture {
        _fasta: NamedTempFile,
        reference: ReferenceReader,
        header: Header,
    }

    #[fixture]
    fn fx() -> Fixture {
        let fasta = create_test_fasta(&[("chr1", CHR1), ("chr2", CHR2)]).expect("fasta");
        let reference = ReferenceReader::new(fasta.path()).expect("reference");
        let len = NonZeroUsize::new(200).expect("non-zero");
        let header = Header::builder()
            .add_reference_sequence(b"chr1", Map::<ReferenceSequence>::new(len))
            .add_reference_sequence(b"chr2", Map::<ReferenceSequence>::new(len))
            .build();
        Fixture { _fasta: fasta, reference, header }
    }

    /// One end of the merged input. Its bases are those of its true in-window
    /// placement (R1: `chr1:101-120:+`, R2: `chr1:161-180:-`), stored in the
    /// orientation of `placement`.
    #[derive(Debug, Clone, Copy)]
    struct EndSpec {
        /// `(tid, 1-based pos, reverse, CIGAR)`; `None` = unmapped.
        placement: Option<(i32, i32, bool, &'static str)>,
        nm: i32,
        score: i32,
        xa: Option<&'static str>,
        rw: Option<&'static str>,
        supplementary: bool,
        /// The primary carries no `NM` (recomputed from the reference on demotion).
        without_nm: bool,
        /// Read index (genome orientation of the true placement) changed to `A`.
        mismatch_at: Option<usize>,
        /// Read indices (genome orientation) of reference `C`s read as `T`, as a
        /// bisulfite/EM-seq conversion reads them.
        conversions: &'static [usize],
    }

    const fn end(
        placement: Option<(i32, i32, bool, &'static str)>,
        nm: i32,
        score: i32,
        xa: Option<&'static str>,
        rw: Option<&'static str>,
    ) -> EndSpec {
        EndSpec {
            placement,
            nm,
            score,
            xa,
            rw,
            supplementary: false,
            without_nm: false,
            mismatch_at: None,
            conversions: &[],
        }
    }

    impl EndSpec {
        const fn with_supplementary(mut self) -> Self {
            self.supplementary = true;
            self
        }
        const fn without_nm(mut self) -> Self {
            self.without_nm = true;
            self
        }
        const fn with_mismatch(mut self, at: usize) -> Self {
            self.mismatch_at = Some(at);
            self
        }
        const fn with_conversions(mut self, at: &'static [usize]) -> Self {
            self.conversions = at;
            self
        }
    }

    /// Expected primary fields for one end after selection.
    #[derive(Debug, Clone, Copy)]
    struct Expect {
        tid: i32,
        pos1: i32,
        reverse: bool,
        cigar: &'static str,
        mapq: u8,
        nm: i64,
        score: i64,
        xs: Option<i64>,
        xa: Option<&'static str>,
        md: &'static str,
        tlen: i32,
    }

    const fn expect(
        tid: i32,
        pos1: i32,
        reverse: bool,
        cigar: &'static str,
        mapq: u8,
        nm: i64,
        score: i64,
        xs: Option<i64>,
        xa: Option<&'static str>,
        md: &'static str,
        tlen: i32,
    ) -> Expect {
        Expect { tid, pos1, reverse, cigar, mapq, nm, score, xs, xa, md, tlen }
    }

    #[derive(Debug, Clone, Copy, PartialEq, Eq)]
    enum Seg {
        Fragment,
        R1,
        R2,
    }

    /// `(genome bases, reverse)` of the end's true placement.
    fn truth(seg: Seg) -> (Vec<u8>, bool) {
        match seg {
            Seg::Fragment | Seg::R1 => (CHR1.as_bytes()[100..120].to_vec(), false),
            Seg::R2 => (CHR1.as_bytes()[160..180].to_vec(), true),
        }
    }

    fn quals() -> Vec<u8> {
        (10..30).collect()
    }
    fn ad() -> Vec<i16> {
        (0..20).collect()
    }
    fn ae() -> Vec<i16> {
        (100..120).collect()
    }

    fn seg_flags(seg: Seg) -> u16 {
        match seg {
            Seg::Fragment => 0,
            Seg::R1 => PAIRED | FIRST_SEGMENT,
            Seg::R2 => PAIRED | LAST_SEGMENT,
        }
    }

    /// The stored (BAM-orientation) bases of an end.
    fn stored_seq(spec: &EndSpec, seg: Seg) -> Vec<u8> {
        let (mut genome, truth_rev) = truth(seg);
        if let Some(i) = spec.mismatch_at {
            genome[i] = b'A';
        }
        for &i in spec.conversions {
            assert_eq!(genome[i], b'C', "a conversion needs a reference C");
            genome[i] = b'T';
        }
        let rev = spec.placement.map_or(truth_rev, |p| p.2);
        if rev == truth_rev { genome } else { fgumi_dna::dna::reverse_complement(&genome) }
    }

    fn build_records(spec: &EndSpec, seg: Seg, proper: bool) -> Vec<RawRecord> {
        let seq = stored_seq(spec, seg);
        let mut b = SamBuilder::new();
        let mut flag = seg_flags(seg);
        if proper {
            flag |= PROPER_PAIR;
        }
        b.read_name(b"q1").sequence(&seq).qualities(&quals());
        match spec.placement {
            Some((tid, pos1, rev, cigar)) => {
                if rev {
                    flag |= REVERSE;
                }
                b.ref_id(tid)
                    .pos(pos1 - 1)
                    .mapq(60)
                    .cigar_ops(&fgumi_raw_bam::parse_cigar_string_ops(cigar).expect("cigar"));
                if !spec.without_nm {
                    b.add_int_tag(SamTag::NM, spec.nm);
                }
                b.add_int_tag(SamTag::AS, spec.score).add_string_tag(SamTag::MD, b"20");
            }
            None => {
                flag |= UNMAPPED;
            }
        }
        b.flags(flag);
        if let Some(xa) = spec.xa {
            b.add_string_tag(SamTag::XA, xa.as_bytes());
        }
        if spec.supplementary {
            b.add_string_tag(SamTag::SA, b"chr2,150,+,10S10M,0,0;");
        }
        if let Some(rw) = spec.rw {
            b.add_string_tag(SamTag::RW, rw.as_bytes());
        }
        b.add_array_i16(SamTag::AD_BASES, &ad())
            .add_array_i16(SamTag::AE_BASES, &ae())
            .add_int_tag(SamTag::CD, 7)
            .add_float_tag(SamTag::CE, 0.01);
        let mut records = vec![b.build()];
        if spec.supplementary {
            let mut s = SamBuilder::new();
            s.read_name(b"q1")
                .flags(seg_flags(seg) | SUPPLEMENTARY)
                .ref_id(1)
                .pos(149)
                .mapq(60)
                .cigar_ops(&fgumi_raw_bam::parse_cigar_string_ops("10H10M").expect("cigar"))
                .sequence(&seq[10..])
                .qualities(&quals()[10..]);
            s.add_int_tag(SamTag::NM, 0).add_int_tag(SamTag::AS, 10);
            if let Some(rw) = spec.rw {
                s.add_string_tag(SamTag::RW, rw.as_bytes());
            }
            records.push(s.build());
        }
        records
    }

    fn build_template(r1: &EndSpec, r2: Option<&EndSpec>, proper: bool) -> Template {
        let mut records = Vec::new();
        match r2 {
            Some(r2) => {
                records.extend(build_records(r1, Seg::R1, proper));
                records.extend(build_records(r2, Seg::R2, proper));
            }
            None => records.extend(build_records(r1, Seg::Fragment, false)),
        }
        let mut template = Template::from_records(records).expect("template");
        if r2.is_some() {
            template.fix_mate_info().expect("fix_mate_info");
        }
        template
    }

    fn strip_rw(record: &RawRecord) -> Vec<u8> {
        let mut bytes = record.as_ref().to_vec();
        fgumi_raw_bam::remove_tag(&mut bytes, SamTag::RW);
        bytes
    }

    fn aux(record: &RawRecord) -> &[u8] {
        fgumi_raw_bam::aux_data_slice(record)
    }

    fn i16_array(record: &RawRecord, tag: SamTag) -> Vec<i16> {
        let arr = fgumi_raw_bam::find_array_tag(aux(record), tag).expect("array tag");
        (0..arr.count)
            .map(|i| fgumi_raw_bam::array_tag_element_u16(&arr, i).cast_signed())
            .collect()
    }

    fn assert_end(record: &RawRecord, exp: &Expect) {
        let view = RawRecordView::new(record);
        assert_eq!(view.ref_id(), exp.tid, "tid");
        assert_eq!(view.pos() + 1, exp.pos1, "pos");
        assert_eq!(view.flags() & REVERSE != 0, exp.reverse, "strand");
        assert_eq!(fgumi_raw_bam::cigar_to_string_from_raw(record), exp.cigar, "cigar");
        assert_eq!(view.mapq(), exp.mapq, "mapq");
        assert_eq!(fgumi_raw_bam::find_int_tag(aux(record), SamTag::NM), Some(exp.nm), "NM");
        assert_eq!(fgumi_raw_bam::find_int_tag(aux(record), SamTag::AS), Some(exp.score), "AS");
        assert_eq!(fgumi_raw_bam::find_int_tag(aux(record), SamTag::XS), exp.xs, "XS");
        assert_eq!(
            fgumi_raw_bam::find_string_tag(aux(record), SamTag::XA),
            exp.xa.map(str::as_bytes),
            "XA"
        );
        assert_eq!(
            fgumi_raw_bam::find_string_tag(aux(record), SamTag::MD),
            Some(exp.md.as_bytes()),
            "MD"
        );
        assert_eq!(fgumi_raw_bam::template_length(record), exp.tlen, "TLEN");
    }

    fn supplementary_count(template: &Template, first: bool) -> usize {
        template
            .records
            .iter()
            .filter(|r| {
                let f = RawRecordView::new(r).flags();
                f & SUPPLEMENTARY != 0 && (f & LAST_SEGMENT == 0) == first
            })
            .count()
    }

    // R1/fragment truth: chr1:101-120:+. R2 truth: chr1:161-180:-. TLEN of the
    // in-window pair is 80.
    #[rstest]
    #[case::primary_in_window(
        end(Some((0, 101, false, "20M")), 0, 20, Some("chr2,+51,20M,3;"), Some(RW1)),
        Some(end(Some((0, 161, true, "20M")), 0, 20, None, Some(RW2))),
        true,
        [Selection::Kept, Selection::Kept],
        expect(0, 101, false, "20M", 60, 0, 20, None, Some("chr2,+51,20M,3;"), "20", 80),
        Some(expect(0, 161, true, "20M", 60, 0, 20, None, None, "20", -80)),
        Some(true),
        [false, false],
    )]
    #[case::primary_out_xa_in(
        end(
            Some((1, 51, false, "20M")), 2, 12,
            Some("chr2,+181,20M,4;chr1,+101,20M,1;"), Some(RW1),
        ).with_supplementary().with_mismatch(5),
        Some(end(Some((0, 161, true, "20M")), 0, 20, None, Some(RW2))),
        false,
        [Selection::Swapped, Selection::Kept],
        expect(
            0, 101, false, "20M", 0, 1, 15, Some(12),
            Some("chr2,+181,20M,4;chr2,+51,20M,2;"), "5T14", 80,
        ),
        Some(expect(0, 161, true, "20M", 60, 0, 20, None, None, "20", -80)),
        Some(false),
        [false, false],
    )]
    // Methylation: bwa's XA NM (0) does not count the two C->T conversions, so
    // neither does AS; MD, regenerated from the reference, lists them.
    #[case::methylation_conversions_not_scored(
        end(Some((1, 51, false, "20M")), 0, 20, Some("chr1,+101,20M,0;"), Some(RW1))
            .with_conversions(&[6, 9]),
        Some(end(Some((0, 161, true, "20M")), 0, 20, None, Some(RW2))),
        false,
        [Selection::Swapped, Selection::Kept],
        expect(0, 101, false, "20M", 0, 0, 20, Some(20), Some("chr2,+51,20M,0;"), "6C2C10", 80),
        Some(expect(0, 161, true, "20M", 60, 0, 20, None, None, "20", -80)),
        Some(false),
        [false, false],
    )]
    #[case::swap_flips_strand(
        end(Some((1, 51, true, "20M")), 0, 20, Some("chr1,+103,2S18M,0;"), Some(RW1)),
        Some(end(Some((0, 161, true, "20M")), 0, 20, None, Some(RW2))),
        false,
        [Selection::Swapped, Selection::Kept],
        expect(0, 103, false, "2S18M", 0, 0, 18, Some(20), Some("chr2,-51,20M,0;"), "18", 78),
        Some(expect(0, 161, true, "20M", 60, 0, 20, None, None, "20", -78)),
        Some(false),
        [true, false],
    )]
    #[case::primary_out_xa_wrong_strand(
        end(Some((1, 51, false, "20M")), 0, 20, Some("chr1,-101,20M,0;"), Some(RW1)),
        Some(end(Some((0, 161, true, "20M")), 0, 20, None, Some(RW2))),
        false,
        [Selection::OutsideWindow, Selection::Kept],
        expect(1, 51, false, "20M", 60, 0, 20, None, Some("chr1,-101,20M,0;"), "20", 0),
        Some(expect(0, 161, true, "20M", 60, 0, 20, None, None, "20", 0)),
        Some(false),
        [false, false],
    )]
    #[case::no_xa(
        end(Some((1, 51, false, "20M")), 0, 20, None, Some(RW1)),
        Some(end(Some((0, 161, true, "20M")), 0, 20, None, Some(RW2))),
        false,
        [Selection::OutsideWindow, Selection::Kept],
        expect(1, 51, false, "20M", 60, 0, 20, None, None, "20", 0),
        Some(expect(0, 161, true, "20M", 60, 0, 20, None, None, "20", 0)),
        Some(false),
        [false, false],
    )]
    // Both ends end up in their windows, FR on one contig, but bwa had not
    // flagged the pair proper, so it stays improper: a swap never sets the flag.
    #[case::r1_swapped_r2_kept(
        end(Some((0, 11, false, "20M")), 0, 20, Some("chr1,+101,20M,0;"), Some(RW1)),
        Some(end(Some((0, 161, true, "20M")), 0, 20, None, Some(RW2))),
        false,
        [Selection::Swapped, Selection::Kept],
        expect(0, 101, false, "20M", 0, 0, 20, Some(20), Some("chr1,+11,20M,0;"), "20", 80),
        Some(expect(0, 161, true, "20M", 60, 0, 20, None, None, "20", -80)),
        Some(false),
        [false, false],
    )]
    // The proper-pair guard: bwa flagged the pair proper and R2 is in its
    // window, so R1's in-window XA hit is not promoted and nothing changes.
    #[case::proper_pair_guard_blocks_xa_swap(
        end(Some((0, 11, false, "20M")), 0, 20, Some("chr1,+101,20M,0;"), Some(RW1)),
        Some(end(Some((0, 161, true, "20M")), 0, 20, None, Some(RW2))),
        true,
        [Selection::KeptProperPair, Selection::Kept],
        expect(0, 11, false, "20M", 60, 0, 20, None, Some("chr1,+101,20M,0;"), "20", 170),
        Some(expect(0, 161, true, "20M", 60, 0, 20, None, None, "20", -170)),
        Some(true),
        [false, false],
    )]
    // The same guard for R2 (R1 in its window), so it does not depend on which
    // end is processed first.
    #[case::proper_pair_guard_blocks_r2_xa_swap(
        end(Some((0, 101, false, "20M")), 0, 20, None, Some(RW1)),
        Some(end(Some((0, 191, true, "20M")), 0, 20, Some("chr1,-161,20M,0;"), Some(RW2))),
        true,
        [Selection::Kept, Selection::KeptProperPair],
        expect(0, 101, false, "20M", 60, 0, 20, None, None, "20", 110),
        Some(expect(0, 191, true, "20M", 60, 0, 20, None, Some("chr1,-161,20M,0;"), "20", -110)),
        Some(true),
        [false, false],
    )]
    // The guard is decided before either end moves: both ends start outside
    // their windows, so neither is guarded, though after R1's swap R1 would be
    // in its window. bwa had flagged the pair proper and the swaps leave it FR on
    // one contig with a smaller insert (120 -> 80), so the flag is kept.
    #[case::proper_pair_guard_uses_pre_swap_state(
        end(Some((1, 51, false, "20M")), 0, 20, Some("chr1,+101,20M,0;"), Some(RW1)),
        Some(end(Some((1, 151, false, "20M")), 0, 20, Some("chr1,-161,20M,0;"), Some(RW2))),
        true,
        [Selection::Swapped, Selection::Swapped],
        expect(0, 101, false, "20M", 0, 0, 20, Some(20), Some("chr2,+51,20M,0;"), "20", 80),
        Some(expect(0, 161, true, "20M", 0, 0, 20, Some(20), Some("chr2,+151,20M,0;"), "20", -80)),
        Some(true),
        [false, true],
    )]
    // A proper pair whose mate has no window is not guarded; the swap leaves it
    // FR on one contig with a smaller insert (170 -> 80), so the flag is kept.
    #[case::proper_pair_mate_without_window_not_guarded(
        end(Some((0, 11, false, "20M")), 0, 20, Some("chr1,+101,20M,0;"), Some(RW1)),
        Some(end(Some((0, 161, true, "20M")), 0, 20, None, None)),
        true,
        [Selection::Swapped, Selection::NoWindow],
        expect(0, 101, false, "20M", 0, 0, 20, Some(20), Some("chr1,+11,20M,0;"), "20", 80),
        Some(expect(0, 161, true, "20M", 60, 0, 20, None, None, "20", -80)),
        Some(true),
        [false, false],
    )]
    // The windows are no insert-size check: a swap that widens the insert
    // (40 -> 80) clears bwa's flag. (R2 is outside its window, so the
    // proper-pair guard does not apply.)
    #[case::r1_swapped_insert_grows_clears_proper(
        end(Some((0, 141, false, "20M")), 0, 20, Some("chr1,+101,20M,0;"), Some(RW1)),
        Some(end(Some((0, 161, true, "20M")), 0, 20, None, Some(RW2_ELSEWHERE))),
        true,
        [Selection::Swapped, Selection::OutsideWindow],
        expect(0, 101, false, "20M", 0, 0, 20, Some(20), Some("chr1,+141,20M,0;"), "20", 80),
        Some(expect(0, 161, true, "20M", 60, 0, 20, None, None, "20", -80)),
        Some(false),
        [false, false],
    )]
    // Both ends swapped (R2 onto the other strand), each with a supplementary:
    // both supplementaries are dropped and the mate fields describe the two
    // promoted placements. bwa had not paired the old placements, so the pair
    // stays improper.
    #[case::both_swapped(
        end(Some((1, 51, false, "20M")), 0, 20, Some("chr1,+101,20M,0;"), Some(RW1))
            .with_supplementary(),
        Some(
            end(Some((1, 151, false, "20M")), 0, 20, Some("chr1,-161,20M,0;"), Some(RW2))
                .with_supplementary(),
        ),
        false,
        [Selection::Swapped, Selection::Swapped],
        expect(0, 101, false, "20M", 0, 0, 20, Some(20), Some("chr2,+51,20M,0;"), "20", 80),
        Some(expect(0, 161, true, "20M", 0, 0, 20, Some(20), Some("chr2,+151,20M,0;"), "20", -80)),
        Some(false),
        [false, true],
    )]
    // Window edges: a primary ending one base before its window starts is
    // outside it...
    #[case::primary_ends_before_window(
        end(Some((0, 71, false, "20M")), 0, 20, None, Some(RW1)),
        Some(end(Some((0, 161, true, "20M")), 0, 20, None, Some(RW2))),
        false,
        [Selection::OutsideWindow, Selection::Kept],
        expect(0, 71, false, "20M", 60, 0, 20, None, None, "20", 110),
        Some(expect(0, 161, true, "20M", 60, 0, 20, None, None, "20", -110)),
        Some(false),
        [false, false],
    )]
    // ...and one whose last base is the window's first base is inside it.
    #[case::primary_last_base_on_window_start(
        end(Some((0, 72, false, "20M")), 0, 20, None, Some(RW1)),
        Some(end(Some((0, 161, true, "20M")), 0, 20, None, Some(RW2))),
        false,
        [Selection::Kept, Selection::Kept],
        expect(0, 72, false, "20M", 60, 0, 20, None, None, "20", 109),
        Some(expect(0, 161, true, "20M", 60, 0, 20, None, None, "20", -109)),
        Some(false),
        [false, false],
    )]
    // R2 is outside its window, so the proper-pair guard does not apply; bwa
    // had flagged the pair proper and the swap leaves it FR on one contig with a
    // smaller insert, so the flag is kept...
    #[case::r1_swapped_r2_outside_keeps_bwa_proper(
        end(Some((0, 11, false, "20M")), 0, 20, Some("chr1,+101,20M,0;"), Some(RW1)),
        Some(end(Some((0, 161, true, "20M")), 0, 20, None, Some(RW2_ELSEWHERE))),
        true,
        [Selection::Swapped, Selection::OutsideWindow],
        expect(0, 101, false, "20M", 0, 0, 20, Some(20), Some("chr1,+11,20M,0;"), "20", 80),
        Some(expect(0, 161, true, "20M", 60, 0, 20, None, None, "20", -80)),
        Some(true),
        [false, false],
    )]
    // ...and an improper pair is left improper.
    #[case::r1_swapped_r2_outside_stays_improper(
        end(Some((0, 11, false, "20M")), 0, 20, Some("chr1,+101,20M,0;"), Some(RW1)),
        Some(end(Some((0, 161, true, "20M")), 0, 20, None, Some(RW2_ELSEWHERE))),
        false,
        [Selection::Swapped, Selection::OutsideWindow],
        expect(0, 101, false, "20M", 0, 0, 20, Some(20), Some("chr1,+11,20M,0;"), "20", 80),
        Some(expect(0, 161, true, "20M", 60, 0, 20, None, None, "20", -80)),
        Some(false),
        [false, false],
    )]
    // (R2 is outside its window, so the proper-pair guard does not apply.)
    #[case::r1_swapped_moves_pair_off_contig(
        end(Some((0, 11, false, "20M")), 0, 20, Some("chr2,+101,20M,16;"), Some("chr2:91-130:+")),
        Some(end(Some((0, 161, true, "20M")), 0, 20, None, Some(RW2_ELSEWHERE))),
        true,
        [Selection::Swapped, Selection::OutsideWindow],
        expect(
            1, 101, false, "20M", 0, 16, 0, Some(20), Some("chr1,+11,20M,0;"),
            "0A0C1A0C0G1G0C0T1A0T0T0C1C0T0T0G0", 0,
        ),
        Some(expect(0, 161, true, "20M", 60, 0, 20, None, None, "20", 0)),
        Some(false),
        [false, false],
    )]
    #[case::fragment(
        end(Some((1, 51, false, "20M")), 0, 20, Some("chr1,+101,20M,0;"), Some(RW1))
            .with_supplementary(),
        None,
        false,
        [Selection::Swapped, Selection::Absent],
        expect(0, 101, false, "20M", 0, 0, 20, Some(20), Some("chr2,+51,20M,0;"), "20", 0),
        None,
        None,
        [false, false],
    )]
    #[case::unmapped_end(
        end(Some((0, 101, false, "20M")), 0, 20, None, Some(RW1)),
        Some(end(None, 0, 0, None, Some(RW2))),
        false,
        [Selection::Kept, Selection::Unmapped],
        expect(0, 101, false, "20M", 60, 0, 20, None, None, "20", 0),
        None,
        None,
        [false, false],
    )]
    #[case::no_window(
        end(Some((1, 51, false, "20M")), 0, 20, Some("chr1,+101,20M,0;"), None),
        Some(end(Some((0, 161, true, "20M")), 0, 20, None, Some(RW2))),
        false,
        [Selection::NoWindow, Selection::Kept],
        expect(1, 51, false, "20M", 60, 0, 20, None, Some("chr1,+101,20M,0;"), "20", 0),
        Some(expect(0, 161, true, "20M", 60, 0, 20, None, None, "20", 0)),
        Some(false),
        [false, false],
    )]
    // An R2 swap from a forward primary onto its reverse-strand window: the
    // per-base tags are reversed after entering the reverse strand.
    #[case::r2_swapped_forward_to_reverse(
        end(Some((0, 101, false, "20M")), 0, 20, None, Some(RW1)),
        Some(end(Some((1, 51, false, "20M")), 0, 20, Some("chr1,-161,20M,0;"), Some(RW2))),
        false,
        [Selection::Kept, Selection::Swapped],
        expect(0, 101, false, "20M", 60, 0, 20, None, None, "20", 80),
        Some(expect(0, 161, true, "20M", 0, 0, 20, Some(20), Some("chr2,+51,20M,0;"), "20", -80)),
        Some(false),
        [false, true],
    )]
    // R1's supplementary is dropped and R2's survives with its range shifted.
    #[case::r1_swapped_both_have_supplementaries(
        end(Some((1, 51, false, "20M")), 0, 20, Some("chr1,+101,20M,0;"), Some(RW1))
            .with_supplementary(),
        Some(end(Some((0, 161, true, "20M")), 0, 20, None, Some(RW2)).with_supplementary()),
        false,
        [Selection::Swapped, Selection::Kept],
        expect(0, 101, false, "20M", 0, 0, 20, Some(20), Some("chr2,+51,20M,0;"), "20", 80),
        Some(expect(0, 161, true, "20M", 60, 0, 20, None, None, "20", -80)),
        Some(false),
        [false, false],
    )]
    // Of two in-window hits the first is promoted.
    #[case::first_in_window_hit_wins(
        end(
            Some((1, 51, false, "20M")), 0, 20,
            Some("chr1,+101,20M,1;chr1,+102,1S19M,0;"), Some(RW1),
        ),
        Some(end(Some((0, 161, true, "20M")), 0, 20, None, Some(RW2))),
        false,
        [Selection::Swapped, Selection::Kept],
        expect(
            0, 101, false, "20M", 0, 1, 15, Some(20),
            Some("chr1,+102,1S19M,0;chr2,+51,20M,0;"), "20", 80,
        ),
        Some(expect(0, 161, true, "20M", 60, 0, 20, None, None, "20", -80)),
        Some(false),
        [false, false],
    )]
    // An in-window hit whose CIGAR does not cover the read cannot be promoted.
    #[case::in_window_hit_short_cigar(
        end(Some((1, 51, false, "20M")), 0, 20, Some("chr1,+101,19M,0;"), Some(RW1)),
        Some(end(Some((0, 161, true, "20M")), 0, 20, None, Some(RW2))),
        false,
        [Selection::OutsideWindow, Selection::Kept],
        expect(1, 51, false, "20M", 60, 0, 20, None, Some("chr1,+101,19M,0;"), "20", 0),
        Some(expect(0, 161, true, "20M", 60, 0, 20, None, None, "20", 0)),
        Some(false),
        [false, false],
    )]
    // A demoted primary without NM gets its edit distance from the reference.
    #[case::demoted_primary_without_nm(
        end(Some((1, 51, false, "20M")), 0, 20, Some("chr1,+101,20M,0;"), Some(RW1)).without_nm(),
        Some(end(Some((0, 161, true, "20M")), 0, 20, None, Some(RW2))),
        false,
        [Selection::Swapped, Selection::Kept],
        expect(0, 101, false, "20M", 0, 0, 20, Some(20), Some("chr2,+51,20M,13;"), "20", 80),
        Some(expect(0, 161, true, "20M", 60, 0, 20, None, None, "20", -80)),
        Some(false),
        [false, false],
    )]
    // The proper-pair guard needs a candidate: an end outside its window with
    // neither an in-window XA hit nor a split piece is counted outside it, even
    // in a proper pair whose mate is in its window, and nothing changes.
    #[case::proper_pair_without_candidate_is_outside(
        end(Some((1, 51, false, "20M")), 0, 20, None, Some(RW1)),
        Some(end(Some((0, 161, true, "20M")), 0, 20, None, Some(RW2))),
        true,
        [Selection::OutsideWindow, Selection::Kept],
        expect(1, 51, false, "20M", 60, 0, 20, None, None, "20", 0),
        Some(expect(0, 161, true, "20M", 60, 0, 20, None, None, "20", 0)),
        Some(true),
        [false, false],
    )]
    // A swap that leaves |TLEN| as it was keeps bwa's flag: the bound is
    // inclusive. R1 moves from chr1:82-101 reverse (5' end 101) onto the other
    // strand at 101, so TLEN stays 80. (R2 is outside its window, so the
    // proper-pair guard does not apply.)
    #[case::r1_swapped_equal_insert_keeps_proper(
        end(Some((0, 82, true, "20M")), 0, 20, Some("chr1,+101,20M,0;"), Some(RW1)),
        Some(end(Some((0, 161, true, "20M")), 0, 20, None, Some(RW2_ELSEWHERE))),
        true,
        [Selection::Swapped, Selection::OutsideWindow],
        expect(0, 101, false, "20M", 0, 0, 20, Some(20), Some("chr1,-82,20M,0;"), "20", 80),
        Some(expect(0, 161, true, "20M", 60, 0, 20, None, None, "20", -80)),
        Some(true),
        [true, false],
    )]
    fn select_primary_cases(
        fx: Fixture,
        #[case] r1: EndSpec,
        #[case] r2: Option<EndSpec>,
        #[case] proper_in: bool,
        #[case] expected: [Selection; 2],
        #[case] exp_r1: Expect,
        #[case] exp_r2: Option<Expect>,
        #[case] expected_proper: Option<bool>,
        #[case] flipped: [bool; 2],
    ) {
        let input = build_template(&r1, r2.as_ref(), proper_in);
        let mut merged = input.clone();
        let got = select_primary(&mut merged, &fx.reference, &fx.header, &BwaScoring::default())
            .expect("select_primary")
            .ends;
        assert_eq!(got, expected);

        // `rw` is gone from every record, whatever the outcome.
        for rec in &merged.records {
            assert!(fgumi_raw_bam::find_string_tag(aux(rec), SamTag::RW).is_none());
        }

        let out_r1 = merged.r1().expect("r1");
        assert_end(out_r1, &exp_r1);
        if let Some(exp_r2) = exp_r2 {
            assert_end(merged.r2().expect("r2"), &exp_r2);
        }

        if !expected.contains(&Selection::Swapped) {
            // Nothing moved: every record is byte-identical apart from `rw`.
            let want: Vec<Vec<u8>> = input.records.iter().map(strip_rw).collect();
            let have: Vec<Vec<u8>> = merged.records.iter().map(|r| r.as_ref().to_vec()).collect();
            assert_eq!(have, want);
        }

        // Supplementary ranges stay valid and point at their own end's records.
        for (range, first) in [(merged.r1_supplementals, true), (merged.r2_supplementals, false)] {
            if let Some((start, stop)) = range {
                for rec in &merged.records[start..stop] {
                    let f = RawRecordView::new(rec).flags();
                    assert!(f & SUPPLEMENTARY != 0, "supplementary range");
                    assert_eq!(f & LAST_SEGMENT == 0, first, "supplementary range end");
                }
            }
        }

        let ends = [(merged.r1(), input.r1()), (merged.r2(), input.r2())];
        for (slot, (out_rec, in_rec)) in ends.into_iter().enumerate() {
            if expected[slot] != Selection::Swapped {
                continue;
            }
            let (out_rec, in_rec) = (out_rec.expect("output end"), in_rec.expect("input end"));
            let first = slot == 0;
            assert!(fgumi_raw_bam::find_string_tag(aux(out_rec), SamTag::SA).is_none(), "SA");
            assert_eq!(supplementary_count(&merged, first), 0, "swapped end's supplementaries");
            let in_seq = fgumi_raw_bam::extract_sequence(in_rec);
            let out_seq = fgumi_raw_bam::extract_sequence(out_rec);
            let in_qual = fgumi_raw_bam::quality_scores_slice(in_rec).to_vec();
            let out_qual = fgumi_raw_bam::quality_scores_slice(out_rec).to_vec();
            let (mut want_ad, mut want_ae) = (ad(), ae());
            let (want_seq, want_qual) = if flipped[slot] {
                want_ad.reverse();
                want_ae.reverse();
                let mut q = in_qual.clone();
                q.reverse();
                (fgumi_dna::dna::reverse_complement(&in_seq), q)
            } else {
                (in_seq, in_qual)
            };
            assert_eq!(out_seq, want_seq, "SEQ");
            assert_eq!(out_qual, want_qual, "QUAL");
            assert_eq!(i16_array(out_rec, SamTag::AD_BASES), want_ad, "ad");
            assert_eq!(i16_array(out_rec, SamTag::AE_BASES), want_ae, "ae");
            assert_eq!(fgumi_raw_bam::find_int_tag(aux(out_rec), SamTag::CD), Some(7), "cD");
            assert_eq!(fgumi_raw_bam::find_float_tag(aux(out_rec), SamTag::CE), Some(0.01), "cE");
        }

        match merged.r2() {
            Some(out_r2) => {
                // Mate fields agree in both directions.
                for (rec, mate) in [(out_r1, out_r2), (out_r2, out_r1)] {
                    let (rv, mv) = (RawRecordView::new(rec), RawRecordView::new(mate));
                    if mv.flags() & UNMAPPED == 0 && rv.flags() & UNMAPPED == 0 {
                        assert_eq!(rv.mate_ref_id(), mv.ref_id(), "mate ref");
                        assert_eq!(rv.mate_pos(), mv.pos(), "mate pos");
                        assert_eq!(rv.is_mate_reverse(), mv.flags() & REVERSE != 0, "mate strand");
                    }
                }
                assert_eq!(
                    fgumi_raw_bam::template_length(out_r1),
                    -fgumi_raw_bam::template_length(out_r2),
                    "TLEN symmetry"
                );
                if let Some(proper) = expected_proper {
                    for rec in [out_r1, out_r2] {
                        assert_eq!(
                            RawRecordView::new(rec).flags() & PROPER_PAIR != 0,
                            proper,
                            "proper pair"
                        );
                    }
                }
            }
            None => assert_eq!(merged.records.len(), 1, "fragment keeps one record"),
        }
    }

    /// A swap drops the aligner's tags for the old placement (`SA`, `pa` and the
    /// bisulfite `XG`/`XM`/`YD`) and keeps the read-level `XR`; a kept primary
    /// keeps them all.
    #[rstest]
    #[case::swapped((1, 51, false, "20M"), Selection::Swapped, false)]
    #[case::kept((0, 101, false, "20M"), Selection::Kept, true)]
    fn swap_drops_placement_tags(
        fx: Fixture,
        #[case] placement: (i32, i32, bool, &'static str),
        #[case] expected: Selection,
        #[case] placement_tags_kept: bool,
    ) {
        let r1 = end(Some(placement), 0, 20, Some("chr1,+101,20M,0;"), Some(RW1));
        let mut template = build_template(&r1, None, false);
        let primary = template.r1.expect("r1").0;
        let record = template.records[primary].as_mut_vec();
        // Bismark conversion values (`CT`), written as `str` bytes.
        let ct = "CT".as_bytes();
        fgumi_raw_bam::update_string_tag(record, SamTag::XG, ct);
        fgumi_raw_bam::update_string_tag(record, SamTag::XM, b"....................");
        fgumi_raw_bam::update_string_tag(record, SamTag::XR, ct);
        fgumi_raw_bam::update_string_tag(record, SamTag::YD, b"f");
        fgumi_raw_bam::append_raw_tag(record, SamTag::PA, b'f', &1.0f32.to_le_bytes());
        let got = select_primary(&mut template, &fx.reference, &fx.header, &BwaScoring::default())
            .expect("select_primary");
        assert_eq!(got.ends[0], expected);
        let out = template.r1().expect("r1");
        for tag in [SamTag::XG, SamTag::XM, SamTag::YD, SamTag::PA] {
            assert_eq!(
                fgumi_raw_bam::find_tag_type(aux(out), tag).is_some(),
                placement_tags_kept,
                "{tag}"
            );
        }
        assert_eq!(
            fgumi_raw_bam::find_string_tag(aux(out), SamTag::XR),
            Some("CT".as_bytes()),
            "XR"
        );
    }

    /// A secondary record of `seg` at chr2:101, with `SA` when `split`.
    fn secondary(seg: Seg, split: bool) -> RawRecord {
        let mut b = SamBuilder::new();
        b.read_name(b"q1")
            .flags(seg_flags(seg) | SECONDARY)
            .ref_id(1)
            .pos(100)
            .mapq(0)
            .cigar_ops(&fgumi_raw_bam::parse_cigar_string_ops("20M").expect("cigar"))
            .sequence(&truth(seg).0)
            .qualities(&quals());
        b.add_int_tag(SamTag::AS, 15);
        if split {
            b.add_string_tag(SamTag::SA, b"chr2,150,+,10S10M,0,0;");
        }
        b.build()
    }

    /// After R1 is swapped: its split pieces go, whether flagged supplementary
    /// or, as an aligner run with `-M` (`bwa-mem3 mem --meth`) flags them,
    /// secondary with `SA`; its `-a` alternative (secondary without `SA`) and
    /// R2's split piece stay. Every remaining secondary and supplementary points
    /// at the other end's current primary as its mate, and its `tc` is
    /// recomputed from the moved primaries.
    #[rstest]
    fn swap_drops_split_pieces_and_repoints_the_rest(fx: Fixture) {
        let r1 = end(Some((1, 51, false, "20M")), 0, 20, Some("chr1,+101,20M,0;"), Some(RW1))
            .with_supplementary();
        let r2 = end(Some((0, 161, true, "20M")), 0, 20, None, Some(RW2)).with_supplementary();
        let mut records = build_template(&r1, Some(&r2), false).into_records();
        records.extend([
            secondary(Seg::R1, true),
            secondary(Seg::R1, false),
            secondary(Seg::R2, true),
        ]);
        let mut template = Template::from_records(records).expect("template");
        template.add_template_coordinate_tags();
        let stale_tc = fgumi_raw_bam::find_array_tag(aux(&template.records[2]), SamTag::TC)
            .map(|a| a.data.to_vec());

        let got = select_primary(&mut template, &fx.reference, &fx.header, &BwaScoring::default())
            .expect("select_primary");
        assert_eq!(got.ends, [Selection::Swapped, Selection::Kept]);

        // R1 keeps only its SA-less secondary; R2 keeps its supplementary and its
        // SA-carrying secondary. (`SA` on a primary marks it split.)
        let kinds: Vec<(bool, u16, bool)> = template
            .records
            .iter()
            .map(|r| {
                let f = RawRecordView::new(r).flags();
                let sa = fgumi_raw_bam::find_string_tag(aux(r), SamTag::SA).is_some();
                (f & FIRST_SEGMENT != 0, f & (SECONDARY | SUPPLEMENTARY), sa)
            })
            .collect();
        assert_eq!(
            kinds,
            vec![
                (true, 0, false),
                (false, 0, true),
                (false, SUPPLEMENTARY, false),
                (true, SECONDARY, false),
                (false, SECONDARY, true),
            ]
        );

        let (r1_i, r2_i) = (template.r1.expect("r1").0, template.r2.expect("r2").0);
        let pos_of = |i: usize| RawRecordView::new(&template.records[i]).pos();
        let mut fresh = template.clone();
        fresh.add_template_coordinate_tags();
        for (i, record) in template.records.iter().enumerate() {
            let view = RawRecordView::new(record);
            if view.flags() & (SECONDARY | SUPPLEMENTARY) == 0 {
                continue;
            }
            let mate = if view.flags() & FIRST_SEGMENT != 0 { r2_i } else { r1_i };
            assert_eq!(view.mate_pos(), pos_of(mate), "record {i}: mate pos");
            let tc =
                fgumi_raw_bam::find_array_tag(aux(record), SamTag::TC).map(|a| a.data.to_vec());
            let want = fgumi_raw_bam::find_array_tag(aux(&fresh.records[i]), SamTag::TC)
                .map(|a| a.data.to_vec());
            assert!(tc.is_some(), "record {i}: tc");
            assert_eq!(tc, want, "record {i}: tc from the moved primaries");
            assert_ne!(tc, stale_tc, "record {i}: tc must not describe the old R1");
        }
    }

    /// A promoted six-field (`bwa-mem3 mem -u`) hit's `AS` is the score bwa
    /// lists for it, not the rescored value; four-field hits are rescored. Either
    /// way `AS` and `XS` are stored in the smallest signed type, as the zipper
    /// merge stores them.
    #[rstest]
    #[case::six_field_score_used("chr1,+101,20M,0,17,3;", 17)]
    #[case::four_field_rescored("chr1,+101,20M,0;", 20)]
    fn promoted_as_and_tag_types(fx: Fixture, #[case] xa: &'static str, #[case] expected_as: i64) {
        let r1 = end(Some((1, 51, false, "20M")), 0, 150, Some(xa), Some(RW1));
        let mut template = build_template(&r1, None, false);
        let got = select_primary(&mut template, &fx.reference, &fx.header, &BwaScoring::default())
            .expect("select_primary");
        assert_eq!(got.ends[0], Selection::Swapped);
        let out = template.r1().expect("r1");
        assert_eq!(fgumi_raw_bam::find_int_tag(aux(out), SamTag::AS), Some(expected_as));
        assert_eq!(fgumi_raw_bam::find_int_tag(aux(out), SamTag::XS), Some(150));
        assert_eq!(fgumi_raw_bam::find_tag_type(aux(out), SamTag::AS), Some(b'c'), "AS type");
        assert_eq!(fgumi_raw_bam::find_tag_type(aux(out), SamTag::XS), Some(b's'), "XS type");
    }

    /// A lone mate (its other end filtered away) that is swapped stays one end
    /// of a pair whose mate is unmapped: its segment bit, `RNEXT`/`PNEXT` at its
    /// new placement, `TLEN` 0, no `MC`/`MQ`, not proper.
    #[rstest]
    fn swapped_lone_mate_keeps_its_mate_fields(fx: Fixture) {
        let r2 = end(Some((1, 51, false, "20M")), 0, 20, Some("chr1,-161,20M,0;"), Some(RW2));
        let mut template =
            Template::from_records(build_records(&r2, Seg::R2, false)).expect("template");
        template.set_lone_mate_info();
        let got = select_primary(&mut template, &fx.reference, &fx.header, &BwaScoring::default())
            .expect("select_primary");
        assert_eq!(got.ends, [Selection::Absent, Selection::Swapped]);
        let out = template.r2().expect("lone R2");
        let view = RawRecordView::new(out);
        assert_eq!(view.flags() & (PAIRED | LAST_SEGMENT), PAIRED | LAST_SEGMENT);
        assert_ne!(view.flags() & MATE_UNMAPPED, 0, "mate unmapped");
        assert_eq!(view.flags() & PROPER_PAIR, 0, "never proper");
        assert_eq!((view.ref_id(), view.pos() + 1), (0, 161), "promoted placement");
        assert_eq!((view.mate_ref_id(), view.mate_pos()), (view.ref_id(), view.pos()));
        assert_eq!(fgumi_raw_bam::template_length(out), 0, "TLEN");
        for tag in [SamTag::MC, SamTag::MQ] {
            assert!(fgumi_raw_bam::find_tag_type(aux(out), tag).is_none(), "{tag}");
        }
    }

    /// A malformed `rw` is fgumi's own tag gone wrong, so it fails naming the read.
    #[rstest]
    fn select_primary_rejects_malformed_rw(fx: Fixture) {
        let r1 = end(Some((0, 101, false, "20M")), 0, 20, None, Some("chr1:130-91:+"));
        let mut template = build_template(&r1, None, false);
        let err = select_primary(&mut template, &fx.reference, &fx.header, &BwaScoring::default())
            .expect_err("malformed rw");
        assert!(format!("{err:#}").contains("invalid rw tag on read q1"), "{err:#}");
    }

    /// A malformed `XA` never fails selection: a tag that does not parse is
    /// treated as absent, and a hit with an invalid CIGAR is skipped; each is
    /// counted. A six-field (`bwa-mem3 mem -u`) `XA` is promoted from and
    /// rewritten in the six-field form.
    #[rstest]
    #[case::unparsable_tag("chr1,+101,20M;", Selection::OutsideWindow, 1, 0, "chr1,+101,20M;")]
    #[case::invalid_cigar_skipped(
        "chr1,+101,268435476M,0;chr1,+101,20M,0;",
        Selection::Swapped,
        0,
        1,
        "chr1,+101,268435476M,0;chr2,+51,20M,0;"
    )]
    #[case::only_hit_invalid(
        "chr1,+101,20M5,0;",
        Selection::OutsideWindow,
        0,
        1,
        "chr1,+101,20M5,0;"
    )]
    #[case::six_field(
        "chr2,+181,20M,4,4,0;chr1,+101,20M,0,20,0;",
        Selection::Swapped,
        0,
        0,
        "chr2,+181,20M,4,4,0;chr2,+51,20M,0,20,60;"
    )]
    fn select_primary_xa_edge_cases(
        fx: Fixture,
        #[case] xa: &'static str,
        #[case] expected: Selection,
        #[case] malformed_xa: u32,
        #[case] skipped_xa_hits: u32,
        #[case] expected_xa: &str,
    ) {
        let r1 = end(Some((1, 51, false, "20M")), 0, 20, Some(xa), Some(RW1));
        let r2 = end(Some((0, 161, true, "20M")), 0, 20, None, Some(RW2));
        let mut merged = build_template(&r1, Some(&r2), false);
        let got = select_primary(&mut merged, &fx.reference, &fx.header, &BwaScoring::default())
            .expect("select_primary");
        assert_eq!(
            got,
            SelectOutcome { ends: [expected, Selection::Kept], malformed_xa, skipped_xa_hits }
        );
        let r1 = merged.r1().expect("r1");
        assert_eq!(
            fgumi_raw_bam::find_string_tag(aux(r1), SamTag::XA),
            Some(expected_xa.as_bytes())
        );
    }

    /// The bases of the chimeric test read in read (sequencing) orientation.
    const CHIMERA: &[u8; 20] = b"ACGTTGCAAGGCTTACCGTA";

    /// One alignment of a chimeric end: placement, MAPQ, `NM`, `AS` and its
    /// record kind (0 for the primary, else `SUPPLEMENTARY` or `SECONDARY`).
    #[derive(Debug, Clone, Copy)]
    struct Aln {
        tid: i32,
        pos1: i32,
        reverse: bool,
        cigar: &'static str,
        mapq: u8,
        nm: i32,
        kind: u16,
    }

    const fn aln(tid: i32, pos1: i32, reverse: bool, cigar: &'static str, mapq: u8) -> Aln {
        Aln { tid, pos1, reverse, cigar, mapq, nm: 0, kind: 0 }
    }

    impl Aln {
        const fn supp(mut self) -> Self {
            self.kind = SUPPLEMENTARY;
            self
        }
        const fn sec(mut self) -> Self {
            self.kind = SECONDARY;
            self
        }
        /// `rname,pos,strand,CIGAR,mapQ,NM;` with hard clips as soft clips.
        fn sa_entry(&self) -> String {
            let contig = if self.tid == 0 { "chr1" } else { "chr2" };
            let strand = if self.reverse { '-' } else { '+' };
            let cigar = self.cigar.replace('H', "S");
            format!("{contig},{},{strand},{cigar},{},{};", self.pos1, self.mapq, self.nm)
        }
    }

    /// The read's full SEQ, QUAL, `ad` and `ae` stored on a record of strand
    /// `reverse` (reversed, and SEQ reverse-complemented, on the reverse strand).
    fn stored_chimera(reverse: bool) -> (Vec<u8>, Vec<u8>, Vec<i16>, Vec<i16>) {
        let (mut seq, mut qual, mut ad, mut ae) = (CHIMERA.to_vec(), quals(), ad(), ae());
        if reverse {
            seq = fgumi_dna::dna::reverse_complement(&seq);
            qual.reverse();
            ad.reverse();
            ae.reverse();
        }
        (seq, qual, ad, ae)
    }

    /// The leading and trailing hard-clip lengths of `ops`.
    fn hard_clips(ops: &[u32]) -> (usize, usize) {
        let hard = CigarKind::HardClip as u32;
        let h =
            |op: Option<&u32>| op.filter(|&&o| o & 0xF == hard).map_or(0, |&o| (o >> 4) as usize);
        (h(ops.first()), if ops.len() > 1 { h(ops.last()) } else { 0 })
    }

    /// One chimeric end as the zipper merge leaves it: `alns[0]` the primary,
    /// the rest its split pieces. Every record carries the consensus tags (the
    /// per-base ones full-length, oriented for its strand), `rw`, and `SA`
    /// listing the others, primary first; a hard-clipped record's SEQ/QUAL hold
    /// only its unclipped bases. The primary carries `xa`, if any.
    fn chimeric_records(seg: Seg, alns: &[Aln], rw: &str, xa: Option<&str>) -> Vec<RawRecord> {
        alns.iter()
            .enumerate()
            .map(|(i, a)| {
                let (seq, qual, ad, ae) = stored_chimera(a.reverse);
                let ops = fgumi_raw_bam::parse_cigar_string_ops(a.cigar).expect("cigar");
                let (lead, trail) = hard_clips(&ops);
                let keep = lead..seq.len() - trail;
                let mut b = SamBuilder::new();
                b.read_name(b"q1")
                    .flags(seg_flags(seg) | a.kind | if a.reverse { REVERSE } else { 0 })
                    .ref_id(a.tid)
                    .pos(a.pos1 - 1)
                    .mapq(a.mapq)
                    .cigar_ops(&ops)
                    .sequence(&seq[keep.clone()])
                    .qualities(&qual[keep]);
                b.add_int_tag(SamTag::NM, a.nm)
                    .add_int_tag(SamTag::AS, 8)
                    .add_string_tag(SamTag::MD, b"8");
                if i == 0
                    && let Some(xa) = xa
                {
                    b.add_string_tag(SamTag::XA, xa.as_bytes());
                }
                let sa: String = alns
                    .iter()
                    .enumerate()
                    .filter(|&(j, _)| j != i)
                    .map(|(_, other)| other.sa_entry())
                    .collect();
                b.add_string_tag(SamTag::SA, sa.as_bytes())
                    .add_string_tag(SamTag::RW, rw.as_bytes())
                    .add_array_i16(SamTag::AD_BASES, &ad)
                    .add_array_i16(SamTag::AE_BASES, &ae)
                    .add_int_tag(SamTag::CD, 7)
                    .add_float_tag(SamTag::CE, 0.01)
                    .add_string_tag(SamTag::MI, b"7");
                b.build()
            })
            .collect()
    }

    /// One end of a split-piece test: a plain end, or a chimeric one.
    #[derive(Debug, Clone)]
    enum EndInput {
        Plain(EndSpec),
        Chimeric { alns: Vec<Aln>, rw: &'static str, xa: Option<&'static str> },
    }

    fn chimeric(alns: &[Aln], rw: &'static str) -> EndInput {
        EndInput::Chimeric { alns: alns.to_vec(), rw, xa: None }
    }

    impl EndInput {
        fn records(&self, seg: Seg, proper: bool) -> Vec<RawRecord> {
            match self {
                Self::Plain(spec) => build_records(spec, seg, proper),
                Self::Chimeric { alns, rw, xa } => {
                    let mut records = chimeric_records(seg, alns, rw, *xa);
                    if proper {
                        for r in &mut records {
                            r.set_flags(r.flags() | PROPER_PAIR);
                        }
                    }
                    records
                }
            }
        }
    }

    /// An output record of a chimeric end: `rname,pos,strand,CIGAR` (as written
    /// in `SA`, but with its own clip style), MAPQ, its kind (0, `SUPPLEMENTARY`
    /// or `SECONDARY`) and its `SA`.
    type Rec = (&'static str, u8, u16, &'static str);

    fn rec_of(record: &RawRecord) -> (String, u8, u16, String) {
        let view = RawRecordView::new(record);
        let contig = if view.ref_id() == 0 { "chr1" } else { "chr2" };
        let strand = if view.flags() & REVERSE != 0 { '-' } else { '+' };
        let placement = format!(
            "{contig},{},{strand},{}",
            view.pos() + 1,
            fgumi_raw_bam::cigar_to_string_from_raw(record)
        );
        let sa = fgumi_raw_bam::find_string_tag(aux(record), SamTag::SA)
            .map(|v| String::from_utf8_lossy(v).into_owned())
            .unwrap_or_default();
        (placement, view.mapq(), view.flags() & (SECONDARY | SUPPLEMENTARY), sa)
    }

    fn owned(rec: &Rec) -> (String, u8, u16, String) {
        (rec.0.to_string(), rec.1, rec.2, rec.3.to_string())
    }

    /// Checks every record of one chimeric end after selection: the primary and
    /// the pieces (`expected`, primary first then the pieces in record order)
    /// and, on each, SEQ/QUAL (the full read on the primary and on a soft-clipped
    /// piece, the unclipped bases on a hard-clipped one, in its strand's
    /// orientation) and full-length per-base tags oriented for its strand.
    fn assert_chimeric_end(merged: &Template, first: bool, expected: &[Rec]) {
        let primary = if first { merged.r1 } else { merged.r2 }.expect("primary").0;
        let indices: Vec<usize> =
            std::iter::once(primary).chain(merged.split_piece_indices(first)).collect();
        let got: Vec<_> = indices.iter().map(|&i| rec_of(&merged.records[i])).collect();
        let want: Vec<_> = expected.iter().map(owned).collect();
        assert_eq!(got, want, "records of the end");
        for (k, &i) in indices.iter().enumerate() {
            let record = &merged.records[i];
            let reverse = record.flags() & REVERSE != 0;
            let (seq, qual, per_base_depths, per_base_errors) = stored_chimera(reverse);
            let (lead, trail) = hard_clips(&record.cigar_ops_vec());
            let keep = lead..seq.len() - trail;
            assert_eq!(fgumi_raw_bam::extract_sequence(record), seq[keep.clone()], "SEQ {k}");
            assert_eq!(fgumi_raw_bam::quality_scores_slice(record), &qual[keep], "QUAL {k}");
            assert_eq!(i16_array(record, SamTag::AD_BASES), per_base_depths, "ad {k}");
            assert_eq!(i16_array(record, SamTag::AE_BASES), per_base_errors, "ae {k}");
            assert_eq!(fgumi_raw_bam::find_int_tag(aux(record), SamTag::CD), Some(7), "cD {k}");
            assert_eq!(fgumi_raw_bam::find_float_tag(aux(record), SamTag::CE), Some(0.01), "cE");
            assert_eq!(
                fgumi_raw_bam::find_string_tag(aux(record), SamTag::MI),
                Some(&b"7"[..]),
                "MI {k}"
            );
            let tc = fgumi_raw_bam::find_array_tag(aux(record), SamTag::TC).is_some();
            assert_eq!(tc, k > 0, "tc only on the pieces ({k})");
        }
        let promoted = &merged.records[primary];
        assert!(fgumi_raw_bam::find_string_tag(aux(promoted), SamTag::XA).is_none(), "XA");
    }

    const PIECE_RW_FWD: &str = "chr1:91-130:+";
    const PIECE_RW_REV: &str = "chr1:91-130:-";

    /// R2 of the split-piece cases: kept in its window, chr1:161-180 reverse.
    fn plain_r2() -> EndInput {
        EndInput::Plain(end(Some((0, 161, true, "20M")), 0, 20, None, Some(RW2)))
    }

    /// R2 at chr1:161-180 reverse, outside its window and with no candidate.
    fn r2_outside_its_window() -> EndInput {
        EndInput::Plain(end(Some((0, 161, true, "20M")), 0, 20, None, Some(RW2_ELSEWHERE)))
    }

    // R1 is chimeric (primary on chr2, a piece in R1's window on chr1) unless
    // stated otherwise; R2 is kept in its window at chr1:161-180 reverse.
    #[rstest]
    #[case::hard_clipped_same_strand(
        chimeric(&[aln(1, 51, false, "12M8S", 60), aln(0, 101, false, "12H8M", 30).supp()], PIECE_RW_FWD),
        plain_r2(),
        false,
        [Selection::SwappedSupplementary, Selection::Kept],
        &[
            ("chr1,101,+,12S8M", 30, 0, "chr2,51,+,12M8S,60,0;"),
            ("chr2,51,+,12M8H", 60, SUPPLEMENTARY, "chr1,101,+,12S8M,30,0;"),
        ],
        &[],
        false,
    )]
    #[case::hard_clipped_strand_flip(
        chimeric(&[aln(1, 51, false, "12M8S", 60), aln(0, 101, true, "8M12H", 30).supp()], PIECE_RW_REV),
        plain_r2(),
        false,
        [Selection::SwappedSupplementary, Selection::Kept],
        &[
            ("chr1,101,-,8M12S", 30, 0, "chr2,51,+,12M8S,60,0;"),
            ("chr2,51,+,12M8H", 60, SUPPLEMENTARY, "chr1,101,-,8M12S,30,0;"),
        ],
        &[],
        false,
    )]
    // `-Y`: the piece is soft-clipped, so the demoted primary stays soft-clipped.
    #[case::soft_clipped(
        chimeric(&[aln(1, 51, true, "8S12M", 60), aln(0, 101, false, "12S8M", 30).supp()], PIECE_RW_FWD),
        plain_r2(),
        false,
        [Selection::SwappedSupplementary, Selection::Kept],
        &[
            ("chr1,101,+,12S8M", 30, 0, "chr2,51,-,8S12M,60,0;"),
            ("chr2,51,-,8S12M", 60, SUPPLEMENTARY, "chr1,101,+,12S8M,30,0;"),
        ],
        &[],
        false,
    )]
    // `-M` (as under `--meth`): the piece is flagged secondary, and so is the
    // demoted primary.
    #[case::meth_secondary_piece(
        chimeric(&[aln(1, 51, false, "12M8S", 60), aln(0, 101, false, "12H8M", 30).sec()], PIECE_RW_FWD),
        plain_r2(),
        false,
        [Selection::SwappedSupplementary, Selection::Kept],
        &[
            ("chr1,101,+,12S8M", 30, 0, "chr2,51,+,12M8S,60,0;"),
            ("chr2,51,+,12M8H", 60, SECONDARY, "chr1,101,+,12S8M,30,0;"),
        ],
        &[],
        false,
    )]
    // Of two pieces only the chr1 one is in the window; the other stays where it
    // is, with its `SA` rewritten. (Template order lists the pieces in reverse.)
    #[case::one_of_two_pieces_in_window(
        chimeric(
            &[
                aln(1, 51, false, "10M10S", 60),
                aln(1, 151, false, "10H5M5H", 20).supp(),
                aln(0, 101, false, "15H5M", 30).supp(),
            ],
            PIECE_RW_FWD,
        ),
        plain_r2(),
        false,
        [Selection::SwappedSupplementary, Selection::Kept],
        &[
            ("chr1,101,+,15S5M", 30, 0, "chr2,51,+,10M10S,60,0;chr2,151,+,10S5M5S,20,0;"),
            ("chr2,51,+,10M10H", 60, SUPPLEMENTARY, "chr1,101,+,15S5M,30,0;chr2,151,+,10S5M5S,20,0;"),
            ("chr2,151,+,10H5M5H", 20, SUPPLEMENTARY, "chr1,101,+,15S5M,30,0;chr2,51,+,10M10S,60,0;"),
        ],
        &[],
        false,
    )]
    // Two in-window pieces: the longer reference span wins.
    #[case::longest_in_window_piece(
        chimeric(
            &[
                aln(1, 51, false, "10M10S", 60),
                aln(0, 101, false, "12H8M", 30).supp(),
                aln(0, 111, false, "10H5M5H", 40).supp(),
            ],
            PIECE_RW_FWD,
        ),
        plain_r2(),
        false,
        [Selection::SwappedSupplementary, Selection::Kept],
        &[
            ("chr1,101,+,12S8M", 30, 0, "chr1,111,+,10S5M5S,40,0;chr2,51,+,10M10S,60,0;"),
            ("chr1,111,+,10H5M5H", 40, SUPPLEMENTARY, "chr1,101,+,12S8M,30,0;chr2,51,+,10M10S,60,0;"),
            ("chr2,51,+,10M10H", 60, SUPPLEMENTARY, "chr1,101,+,12S8M,30,0;chr1,111,+,10S5M5S,40,0;"),
        ],
        &[],
        false,
    )]
    // Two in-window pieces of equal span: the first in record order wins (the
    // template lists the pieces in reverse of their input order).
    #[case::equal_spans_first_wins(
        chimeric(
            &[
                aln(1, 51, false, "10M10S", 60),
                aln(0, 101, false, "10H5M5H", 30).supp(),
                aln(0, 111, false, "15H5M", 40).supp(),
            ],
            PIECE_RW_FWD,
        ),
        plain_r2(),
        false,
        [Selection::SwappedSupplementary, Selection::Kept],
        &[
            ("chr1,111,+,15S5M", 40, 0, "chr2,51,+,10M10S,60,0;chr1,101,+,10S5M5S,30,0;"),
            ("chr2,51,+,10M10H", 60, SUPPLEMENTARY, "chr1,111,+,15S5M,40,0;chr1,101,+,10S5M5S,30,0;"),
            ("chr1,101,+,10H5M5H", 30, SUPPLEMENTARY, "chr1,111,+,15S5M,40,0;chr2,51,+,10M10S,60,0;"),
        ],
        &[],
        false,
    )]
    // A primary on chr1 outside the window, proper with R2 (insert 160), and R2
    // outside its own window (so the proper-pair guard does not apply): the swap
    // leaves the pair FR on chr1 with a smaller insert (80), so bwa's flag stays...
    #[case::proper_pair_kept(
        chimeric(&[aln(0, 21, false, "12M8S", 60), aln(0, 101, false, "12H8M", 30).supp()], PIECE_RW_FWD),
        r2_outside_its_window(),
        true,
        [Selection::SwappedSupplementary, Selection::OutsideWindow],
        &[
            ("chr1,101,+,12S8M", 30, 0, "chr1,21,+,12M8S,60,0;"),
            ("chr1,21,+,12M8H", 60, SUPPLEMENTARY, "chr1,101,+,12S8M,30,0;"),
        ],
        &[],
        true,
    )]
    // ...and so it does under `-M`, where the piece and then the demoted
    // primary are secondary: every record of the pair keeps the flag...
    #[case::proper_pair_kept_meth_secondary(
        chimeric(&[aln(0, 21, false, "12M8S", 60), aln(0, 101, false, "12H8M", 30).sec()], PIECE_RW_FWD),
        r2_outside_its_window(),
        true,
        [Selection::SwappedSupplementary, Selection::OutsideWindow],
        &[
            ("chr1,101,+,12S8M", 30, 0, "chr1,21,+,12M8S,60,0;"),
            ("chr1,21,+,12M8H", 60, SECONDARY, "chr1,101,+,12S8M,30,0;"),
        ],
        &[],
        true,
    )]
    // ...and one that widens the insert (40 -> 80) clears it on every record,
    // the demoted primary included, whether it is supplementary...
    #[case::proper_pair_cleared(
        chimeric(&[aln(0, 141, false, "12M8S", 60), aln(0, 101, false, "12H8M", 30).supp()], PIECE_RW_FWD),
        r2_outside_its_window(),
        true,
        [Selection::SwappedSupplementary, Selection::OutsideWindow],
        &[
            ("chr1,101,+,12S8M", 30, 0, "chr1,141,+,12M8S,60,0;"),
            ("chr1,141,+,12M8H", 60, SUPPLEMENTARY, "chr1,101,+,12S8M,30,0;"),
        ],
        &[],
        false,
    )]
    // ...or, under `-M`, secondary.
    #[case::proper_pair_cleared_meth_secondary(
        chimeric(&[aln(0, 141, false, "12M8S", 60), aln(0, 101, false, "12H8M", 30).sec()], PIECE_RW_FWD),
        r2_outside_its_window(),
        true,
        [Selection::SwappedSupplementary, Selection::OutsideWindow],
        &[
            ("chr1,101,+,12S8M", 30, 0, "chr1,141,+,12M8S,60,0;"),
            ("chr1,141,+,12M8H", 60, SECONDARY, "chr1,101,+,12S8M,30,0;"),
        ],
        &[],
        false,
    )]
    // The proper-pair guard: bwa flagged the pair proper and R2 is in its
    // window, so R1's in-window piece is not swapped in and every record of R1
    // stays as it was.
    #[case::proper_pair_guard_blocks_split_piece_swap(
        chimeric(&[aln(0, 21, false, "12M8S", 60), aln(0, 101, false, "12H8M", 30).supp()], PIECE_RW_FWD),
        plain_r2(),
        true,
        [Selection::KeptProperPair, Selection::Kept],
        &[
            ("chr1,21,+,12M8S", 60, 0, "chr1,101,+,12S8M,30,0;"),
            ("chr1,101,+,12H8M", 30, SUPPLEMENTARY, "chr1,21,+,12M8S,60,0;"),
        ],
        &[],
        true,
    )]
    // An improper pair is not guarded: the piece is swapped in as usual.
    #[case::improper_pair_not_guarded(
        chimeric(&[aln(0, 21, false, "12M8S", 60), aln(0, 101, false, "12H8M", 30).supp()], PIECE_RW_FWD),
        plain_r2(),
        false,
        [Selection::SwappedSupplementary, Selection::Kept],
        &[
            ("chr1,101,+,12S8M", 30, 0, "chr1,21,+,12M8S,60,0;"),
            ("chr1,21,+,12M8H", 60, SUPPLEMENTARY, "chr1,101,+,12S8M,30,0;"),
        ],
        &[],
        false,
    )]
    // R1 promotes an XA hit; R2's in-window piece swaps with its primary.
    #[case::r1_by_xa_r2_by_split_piece(
        EndInput::Plain(end(Some((1, 51, false, "20M")), 0, 20, Some("chr1,+101,20M,0;"), Some(RW1))),
        chimeric(&[aln(1, 51, false, "12M8S", 60), aln(0, 161, true, "8M12H", 30).supp()], RW2),
        false,
        [Selection::Swapped, Selection::SwappedSupplementary],
        &[],
        &[
            ("chr1,161,-,8M12S", 30, 0, "chr2,51,+,12M8S,60,0;"),
            ("chr2,51,+,12M8H", 60, SUPPLEMENTARY, "chr1,161,-,8M12S,30,0;"),
        ],
        false,
    )]
    fn split_piece_swap_cases(
        fx: Fixture,
        #[case] r1: EndInput,
        #[case] r2: EndInput,
        #[case] proper_in: bool,
        #[case] expected: [Selection; 2],
        #[case] exp_r1: &[Rec],
        #[case] exp_r2: &[Rec],
        #[case] expected_proper: bool,
    ) {
        let mut records = r1.records(Seg::R1, proper_in);
        records.extend(r2.records(Seg::R2, proper_in));
        let mut merged = Template::from_records(records).expect("template");
        merged.fix_mate_info().expect("fix_mate_info");
        merged.add_template_coordinate_tags();
        let got = select_primary(&mut merged, &fx.reference, &fx.header, &BwaScoring::default())
            .expect("select_primary");
        assert_eq!(got.ends, expected);

        for (first, exp) in [(true, exp_r1), (false, exp_r2)] {
            if !exp.is_empty() {
                assert_chimeric_end(&merged, first, exp);
            }
        }
        // An XA-swapped end keeps no split pieces.
        for (slot, first) in [(0, true), (1, false)] {
            if expected[slot] == Selection::Swapped {
                assert!(merged.split_piece_indices(first).is_empty(), "split pieces dropped");
            }
        }

        let (r1_i, r2_i) = (merged.r1.expect("r1").0, merged.r2.expect("r2").0);
        for (i, record) in merged.records.iter().enumerate() {
            let view = RawRecordView::new(record);
            assert!(fgumi_raw_bam::find_string_tag(aux(record), SamTag::RW).is_none(), "rw");
            let mate = if view.flags() & FIRST_SEGMENT != 0 { r2_i } else { r1_i };
            let mate_view = RawRecordView::new(&merged.records[mate]);
            assert_eq!(
                (view.mate_ref_id(), view.mate_pos()),
                (mate_view.ref_id(), mate_view.pos()),
                "record {i}: mate"
            );
            assert_eq!(
                view.is_mate_reverse(),
                mate_view.flags() & REVERSE != 0,
                "record {i}: mate strand"
            );
        }
        for (i, record) in merged.records.iter().enumerate() {
            assert_eq!(record.flags() & PROPER_PAIR != 0, expected_proper, "record {i}: proper");
        }
        assert_eq!(
            fgumi_raw_bam::template_length(&merged.records[r1_i]),
            -fgumi_raw_bam::template_length(&merged.records[r2_i]),
            "TLEN symmetry"
        );
    }

    /// An in-window `XA` hit wins over an in-window split piece: the hit is
    /// promoted and the end's pieces are dropped.
    #[rstest]
    fn in_window_xa_hit_wins_over_split_piece(fx: Fixture) {
        let r1 = EndInput::Chimeric {
            alns: vec![aln(1, 51, false, "12M8S", 60), aln(0, 101, false, "12H8M", 30).supp()],
            rw: PIECE_RW_FWD,
            xa: Some("chr1,+103,20M,2;"),
        };
        let mut records = r1.records(Seg::R1, false);
        records.extend(plain_r2().records(Seg::R2, false));
        let mut merged = Template::from_records(records).expect("template");
        merged.fix_mate_info().expect("fix_mate_info");
        let got = select_primary(&mut merged, &fx.reference, &fx.header, &BwaScoring::default())
            .expect("select_primary");
        assert_eq!(got.ends, [Selection::Swapped, Selection::Kept]);
        let r1 = merged.r1().expect("r1");
        assert_eq!(rec_of(r1).0, "chr1,103,+,20M", "promoted XA hit");
        assert_eq!(RawRecordView::new(r1).mapq(), 0, "XA promotion zeroes MAPQ");
        assert!(merged.split_piece_indices(true).is_empty(), "pieces dropped");
        assert!(fgumi_raw_bam::find_string_tag(aux(r1), SamTag::SA).is_none(), "SA dropped");
    }

    /// A malformed `XA` is counted and treated as absent, so an in-window split
    /// piece is still swapped in; a piece whose CIGAR (hard clips as soft) does
    /// not cover the read is not a candidate.
    #[rstest]
    #[case::malformed_xa_then_piece(
        Some("chr1,+101,20M;"),
        "12H8M",
        Selection::SwappedSupplementary,
        1
    )]
    #[case::piece_cigar_short(None, "11H8M", Selection::OutsideWindow, 0)]
    fn split_piece_edge_cases(
        fx: Fixture,
        #[case] xa: Option<&'static str>,
        #[case] piece_cigar: &'static str,
        #[case] expected: Selection,
        #[case] malformed_xa: u32,
    ) {
        let r1 = EndInput::Chimeric {
            alns: vec![aln(1, 51, false, "12M8S", 60), aln(0, 101, false, piece_cigar, 30).supp()],
            rw: PIECE_RW_FWD,
            xa,
        };
        let mut merged =
            Template::from_records(r1.records(Seg::Fragment, false)).expect("template");
        let got = select_primary(&mut merged, &fx.reference, &fx.header, &BwaScoring::default())
            .expect("select_primary");
        assert_eq!(
            got,
            SelectOutcome { ends: [expected, Selection::Absent], malformed_xa, skipped_xa_hits: 0 }
        );
    }

    /// `swap_primary_with_split_piece` only accepts one of the end's own pieces.
    #[test]
    fn swap_rejects_a_record_that_is_not_a_piece() {
        let r1 = chimeric(
            &[aln(1, 51, false, "12M8S", 60), aln(0, 101, false, "12H8M", 30).supp()],
            PIECE_RW_FWD,
        );
        let mut records = r1.records(Seg::R1, false);
        records.extend(plain_r2().records(Seg::R2, false));
        let mut merged = Template::from_records(records).expect("template");
        let r2 = merged.r2.expect("r2").0;
        assert!(merged.swap_primary_with_split_piece(true, r2).is_err(), "R2 primary");
        let piece = merged.split_piece_indices(true)[0];
        assert!(merged.swap_primary_with_split_piece(false, piece).is_err(), "R1 piece for R2");
        merged.swap_primary_with_split_piece(true, piece).expect("R1 piece");
    }

    /// The bisulfite tags of a split-piece swap under `bwa-mem3 mem --meth`:
    /// the promoted piece keeps the `XG`/`YD` bwa computed for its own
    /// placement, and its `XM` only if it has a call for every base of the full
    /// read; the demoted primary's `XM` is trimmed with its SEQ when it is
    /// hard-clipped. On every record `XM` has one call per base of SEQ.
    #[rstest]
    #[case::hard_clipped("12H8M", None, Some("hhhhhhhhhhhh"))]
    #[case::soft_clipped("12S8M", Some("zzzzzzzzzzzzZZZZZZZZ"), Some("hhhhhhhhhhhhxxxxxxxx"))]
    fn split_piece_swap_keeps_bisulfite_tags_consistent(
        fx: Fixture,
        #[case] piece_cigar: &'static str,
        #[case] promoted_xm: Option<&str>,
        #[case] demoted_xm: Option<&str>,
    ) {
        let alns = [aln(1, 51, false, "12M8S", 60), aln(0, 101, false, piece_cigar, 30).supp()];
        let mut records = chimeric_records(Seg::R1, &alns, PIECE_RW_FWD, None);
        // The old primary's calls: `h` for its aligned bases, `x` for its clip.
        let tags: [(&[u8], &[u8], &[u8]); 2] = [
            (b"hhhhhhhhhhhhxxxxxxxx", "CT".as_bytes(), b"f"),
            (
                if piece_cigar.contains('H') { b"ZZZZZZZZ" } else { b"zzzzzzzzzzzzZZZZZZZZ" },
                "GA".as_bytes(),
                b"r",
            ),
        ];
        for (record, (xm, xg, yd)) in records.iter_mut().zip(tags) {
            let aux = record.as_mut_vec();
            fgumi_raw_bam::update_string_tag(aux, SamTag::XM, xm);
            fgumi_raw_bam::update_string_tag(aux, SamTag::XG, xg);
            fgumi_raw_bam::update_string_tag(aux, SamTag::YD, yd);
        }
        records.extend(plain_r2().records(Seg::R2, false));
        let mut merged = Template::from_records(records).expect("template");
        merged.fix_mate_info().expect("fix_mate_info");
        let got = select_primary(&mut merged, &fx.reference, &fx.header, &BwaScoring::default())
            .expect("select_primary");
        assert_eq!(got.ends, [Selection::SwappedSupplementary, Selection::Kept]);

        let string = |record: &RawRecord, tag| {
            fgumi_raw_bam::find_string_tag(aux(record), tag)
                .map(|v| String::from_utf8_lossy(v).into_owned())
        };
        let promoted = merged.r1().expect("r1");
        let demoted = &merged.records[merged.split_piece_indices(true)[0]];
        assert_eq!(string(promoted, SamTag::XG).as_deref(), Some("GA"), "promoted XG");
        assert_eq!(string(promoted, SamTag::YD).as_deref(), Some("r"), "promoted YD");
        assert_eq!(string(promoted, SamTag::XM).as_deref(), promoted_xm, "promoted XM");
        assert_eq!(string(demoted, SamTag::XG).as_deref(), Some("CT"), "demoted XG");
        assert_eq!(string(demoted, SamTag::XM).as_deref(), demoted_xm, "demoted XM");
        for record in &merged.records {
            if let Some(xm) = string(record, SamTag::XM) {
                assert_eq!(xm.len(), fgumi_raw_bam::l_seq(record) as usize, "XM length is l_seq");
            }
        }
    }

    /// `merge` adds every counter: two counts with distinct values in every
    /// field sum field by field, so no two counters can be crossed.
    #[test]
    fn realign_counts_merge_adds_every_field() {
        let a = RealignCounts {
            reads: 1,
            in_window: 2,
            swapped: 3,
            swapped_supplementary: 4,
            kept_proper_pair: 5,
            outside_window: 6,
            unmapped: 7,
            no_window: 8,
            malformed_xa: 9,
            skipped_xa_hits: 10,
        };
        let b = RealignCounts {
            reads: 100,
            in_window: 200,
            swapped: 300,
            swapped_supplementary: 400,
            kept_proper_pair: 500,
            outside_window: 600,
            unmapped: 700,
            no_window: 800,
            malformed_xa: 900,
            skipped_xa_hits: 1000,
        };
        let mut merged = a;
        merged.merge(&b);
        assert_eq!(
            merged,
            RealignCounts {
                reads: 101,
                in_window: 202,
                swapped: 303,
                swapped_supplementary: 404,
                kept_proper_pair: 505,
                outside_window: 606,
                unmapped: 707,
                no_window: 808,
                malformed_xa: 909,
                skipped_xa_hits: 1010,
            }
        );
    }

    #[rstest]
    #[case::paired(
        [Selection::Kept, Selection::Swapped], 0, 0,
        RealignCounts { reads: 2, in_window: 1, swapped: 1, ..RealignCounts::default() },
    )]
    #[case::proper_pair_guard(
        [Selection::KeptProperPair, Selection::Kept], 0, 0,
        RealignCounts { reads: 2, kept_proper_pair: 1, in_window: 1, ..RealignCounts::default() },
    )]
    #[case::xa_and_split_piece(
        [Selection::Swapped, Selection::SwappedSupplementary], 0, 0,
        RealignCounts { reads: 2, swapped: 1, swapped_supplementary: 1, ..RealignCounts::default() },
    )]
    #[case::fragment(
        [Selection::OutsideWindow, Selection::Absent], 0, 0,
        RealignCounts { reads: 1, outside_window: 1, ..RealignCounts::default() },
    )]
    #[case::unmapped_and_no_window(
        [Selection::Unmapped, Selection::NoWindow], 0, 0,
        RealignCounts { reads: 2, unmapped: 1, no_window: 1, ..RealignCounts::default() },
    )]
    #[case::xa_issues(
        [Selection::OutsideWindow, Selection::Absent], 1, 2,
        RealignCounts {
            reads: 1,
            outside_window: 1,
            malformed_xa: 1,
            skipped_xa_hits: 2,
            ..RealignCounts::default()
        },
    )]
    fn realign_counts_tally(
        #[case] ends: [Selection; 2],
        #[case] malformed_xa: u32,
        #[case] skipped_xa_hits: u32,
        #[case] expected: RealignCounts,
    ) {
        let mut counts = RealignCounts::default();
        counts.tally(&SelectOutcome { ends, malformed_xa, skipped_xa_hits });
        assert_eq!(counts, expected);
    }
}
