//! Per-molecule realign windows for `fgumi group --realign-window`.
//!
//! A *segment* is one read position of one molecule: the raw primary mapped
//! records that share `(molecule id, tid, unclipped 5' position, strand)`. Under
//! template-coordinate grouping this pairs a duplex molecule's `A` R1 with its
//! `B` R2 (they start at the same 5' position on the same strand), so no `/A`
//! `/B` special-casing is needed: the molecule id excludes the strand suffix.
//!
//! The window of a segment spans the minimum unclipped start to the maximum
//! unclipped end over all of its raw records. It is stamped on each of those
//! records as the `rw` tag so a later re-alignment of the consensus read can
//! prefer the alignment that falls inside it.

use std::sync::Arc;

use ahash::AHashMap;
use anyhow::{Result, anyhow};
use fgumi_raw_bam::RawRecordView;
use fgumi_raw_bam::flags;
use fgumi_sam::record_utils::{unclipped_end_raw, unclipped_start_raw};
use fgumi_umi::MoleculeId;
use noodles::sam::Header;

use crate::sam::SamTag;
use crate::template::Template;

/// `(molecule id, tid, unclipped 5' position, reverse)`.
type SegmentKey = (u64, i32, i32, bool);

/// The reference names of a header, indexed by reference id, built once per
/// run and shared by every batch's [`SegmentWindows`] (so a window holds no
/// contig name of its own).
#[derive(Debug)]
pub(crate) struct ContigNames(Vec<Vec<u8>>);

impl ContigNames {
    /// The reference names of `header`, in reference-id order.
    pub(crate) fn from_header(header: &Header) -> Self {
        Self(header.reference_sequences().keys().map(|name| name.to_vec()).collect())
    }

    fn get(&self, tid: i32) -> Option<&[u8]> {
        usize::try_from(tid).ok().and_then(|i| self.0.get(i)).map(Vec::as_slice)
    }
}

/// Windows for every segment of a batch of templates: each segment's unclipped
/// `(start, end)` span, 1-based inclusive, keyed by its [`SegmentKey`]; the key
/// carries the contig (by id, named through [`ContigNames`]) and the strand.
#[derive(Debug, Clone)]
pub(crate) struct SegmentWindows {
    map: AHashMap<SegmentKey, (u32, u32)>,
    contigs: Arc<ContigNames>,
    /// The global molecule-id offset applied to `Template::mi` after the windows
    /// were computed with local ids ([`Self::rebase`]); subtracted on lookup.
    base: u64,
}

/// The segment key of a raw primary mapped record and its unclipped 5' end, or
/// `None` if the record has no molecule id, is unmapped, or is not primary.
/// `five_prime` computes the unclipped 5' position from the record and strand.
fn segment_key_with(
    id: Option<u64>,
    raw: &[u8],
    five_prime: impl FnOnce(&[u8], bool) -> Option<usize>,
) -> Option<SegmentKey> {
    let id = id?;
    let view = RawRecordView::new(raw);
    let flag = view.flags();
    if flag & (flags::UNMAPPED | flags::SECONDARY | flags::SUPPLEMENTARY) != 0 {
        return None;
    }
    let tid = view.ref_id();
    if tid < 0 {
        return None;
    }
    let reverse = flag & flags::REVERSE != 0;
    Some((id, tid, i32::try_from(five_prime(raw, reverse)?).ok()?, reverse))
}

/// The segment key of a raw primary mapped record (one CIGAR walk).
fn segment_key(id: Option<u64>, raw: &[u8]) -> Option<SegmentKey> {
    segment_key_with(id, raw, |raw, reverse| {
        if reverse { unclipped_end_raw(raw) } else { unclipped_start_raw(raw) }
    })
}

/// Compute the window of every segment among `templates` (which carry their
/// molecule ids in `Template::mi`), naming contigs through `contigs`.
///
/// # Errors
/// Returns an error if a mapped record's reference id is not in `contigs` or a
/// coordinate does not fit a `u32`.
pub(crate) fn compute_segment_windows(
    templates: &[Template],
    contigs: &Arc<ContigNames>,
) -> Result<SegmentWindows> {
    let mut map: AHashMap<SegmentKey, (u32, u32)> = AHashMap::new();
    for template in templates {
        for raw in [template.r1(), template.r2()].into_iter().flatten() {
            let raw: &[u8] = raw;
            // One CIGAR walk per end: the key's 5' position is one of the two.
            let (Some(start), Some(end)) = (unclipped_start_raw(raw), unclipped_end_raw(raw))
            else {
                continue;
            };
            let five_prime = |_: &[u8], reverse: bool| Some(if reverse { end } else { start });
            let Some(key) = segment_key_with(template.mi.id(), raw, five_prime) else { continue };
            if contigs.get(key.1).is_none() {
                return Err(anyhow!("reference id {} is not in the header", key.1));
            }
            // A soft-clip running off the contig start saturates the unclipped start at 0;
            // the window is 1-based, so clamp it to the first base.
            let start =
                u32::try_from(start.max(1)).map_err(|_| anyhow!("unclipped start too large"))?;
            let end = u32::try_from(end).map_err(|_| anyhow!("unclipped end too large"))?;
            map.entry(key)
                .and_modify(|(s, e)| {
                    *s = (*s).min(start);
                    *e = (*e).max(end);
                })
                .or_insert((start, end));
        }
    }
    Ok(SegmentWindows { map, contigs: Arc::clone(contigs), base: 0 })
}

impl SegmentWindows {
    /// The `rw` value of the segment `raw` belongs to within molecule `mi`, if
    /// any, appended to `out`; returns whether one was written.
    ///
    /// `mi` is the template's molecule id after any [`rebase`](Self::rebase).
    pub(crate) fn write_for_record(&self, mi: MoleculeId, raw: &[u8], out: &mut Vec<u8>) -> bool {
        let id = mi.id().and_then(|id| id.checked_sub(self.base));
        let Some(key) = segment_key(id, raw) else {
            return false;
        };
        let (Some(&(start, end)), Some(contig)) = (self.map.get(&key), self.contigs.get(key.1))
        else {
            return false;
        };
        fgumi_sam::realign_window::write_value(contig, start, end, key.3, out);
        true
    }

    /// Stamps the `rw` of `raw`'s segment (if any) onto `scratch`, a copy of
    /// `raw` being serialized, using `rw_buf` as scratch space. The one place
    /// both group emitters (the terminal `SerializeGroups` and the fused
    /// group → consensus bridge) stamp `rw`, so they cannot drift apart.
    pub(crate) fn stamp_rw(
        &self,
        mi: MoleculeId,
        raw: &[u8],
        scratch: &mut Vec<u8>,
        rw_buf: &mut Vec<u8>,
    ) {
        rw_buf.clear();
        if self.write_for_record(mi, raw, rw_buf) {
            fgumi_raw_bam::update_string_tag(scratch, SamTag::RW, rw_buf);
        }
    }

    /// Record the global molecule-id offset `base` applied to `Template::mi`
    /// after the windows were computed with local ids; lookups subtract it.
    pub(crate) fn rebase(&mut self, base: u64) {
        self.base = base;
    }

    /// Heap bytes held by the window map (the shared contig names are not
    /// counted: one copy serves every batch).
    pub(crate) fn estimate_heap_size(&self) -> usize {
        self.map.capacity()
            * (std::mem::size_of::<SegmentKey>() + std::mem::size_of::<(u32, u32)>() + 1)
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use fgumi_raw_bam::SamBuilder;
    use fgumi_raw_bam::flags::{FIRST_SEGMENT, LAST_SEGMENT, PAIRED, REVERSE, UNMAPPED};
    use noodles::sam::header::record::value::{Map, map::ReferenceSequence};
    use rstest::rstest;
    use std::num::NonZeroUsize;

    /// One raw read: `(tid, 0-based pos, cigar, reverse)`; `None` = unmapped.
    type Rd = Option<(i32, i32, &'static str, bool)>;
    /// One template: `(mi, r1, r2)`.
    type Tpl = (MoleculeId, Rd, Rd);

    fn header() -> Header {
        let mut b = Header::builder();
        for name in ["chr1", "chr2"] {
            b = b.add_reference_sequence(
                name,
                Map::<ReferenceSequence>::new(NonZeroUsize::new(1_000_000).unwrap()),
            );
        }
        b.build()
    }

    fn cigar_ops(cigar: &str) -> Vec<u32> {
        fgumi_raw_bam::parse_cigar_string_ops(cigar).expect("cigar")
    }

    fn record(name: &[u8], seg: u16, rd: Rd) -> fgumi_raw_bam::RawRecord {
        let mut b = SamBuilder::new();
        b.read_name(name).sequence(b"ACGT").qualities(&[30; 4]);
        match rd {
            Some((tid, pos, cigar, rev)) => {
                b.flags(PAIRED | seg | if rev { REVERSE } else { 0 })
                    .ref_id(tid)
                    .pos(pos)
                    .cigar_ops(&cigar_ops(cigar));
            }
            None => {
                b.flags(PAIRED | seg | UNMAPPED).ref_id(-1).pos(-1);
            }
        }
        b.build()
    }

    fn templates(specs: &[Tpl]) -> Vec<Template> {
        specs
            .iter()
            .enumerate()
            .map(|(i, (mi, r1, r2))| {
                let name = format!("q{i}");
                let mut t = Template::from_records(vec![
                    record(name.as_bytes(), FIRST_SEGMENT, *r1),
                    record(name.as_bytes(), LAST_SEGMENT, *r2),
                ])
                .expect("template");
                t.mi = *mi;
                t
            })
            .collect()
    }

    fn contigs() -> Arc<ContigNames> {
        Arc::new(ContigNames::from_header(&header()))
    }

    fn window_str(windows: &SegmentWindows, mi: MoleculeId, raw: &[u8]) -> Option<String> {
        let mut out = Vec::new();
        windows.write_for_record(mi, raw, &mut out).then(|| String::from_utf8(out).unwrap())
    }

    const S0: MoleculeId = MoleculeId::Single(0);
    const S1: MoleculeId = MoleculeId::Single(1);

    #[rstest]
    #[case::single_pair(
        &[(S0, Some((0, 100, "10M", false)), Some((0, 200, "10M", true)))],
        &[(0, 1, Some("chr1:101-110:+")), (0, 2, Some("chr1:201-210:-"))],
    )]
    #[case::family_varied_3p_ends(
        &[
            (S0, Some((0, 100, "10M", false)), Some((0, 200, "10M", true))),
            (S0, Some((0, 100, "12M", false)), Some((0, 198, "12M", true))),
        ],
        &[
            (0, 1, Some("chr1:101-112:+")), (1, 1, Some("chr1:101-112:+")),
            (0, 2, Some("chr1:199-210:-")), (1, 2, Some("chr1:199-210:-")),
        ],
    )]
    #[case::soft_clipped_member(
        &[
            (S0, Some((0, 100, "10M", false)), Some((1, 200, "10M", true))),
            (S0, Some((0, 103, "3S5M4S", false)), Some((1, 200, "10M", true))),
        ],
        &[(0, 1, Some("chr1:101-112:+")), (1, 1, Some("chr1:101-112:+")), (1, 2, Some("chr2:201-210:-"))],
    )]
    #[case::clipped_off_contig_start(
        &[(S0, Some((0, 0, "5S10M", false)), Some((0, 200, "10M", true)))],
        &[(0, 1, Some("chr1:1-10:+")), (0, 2, Some("chr1:201-210:-"))],
    )]
    #[case::duplex_ab(
        &[
            (MoleculeId::PairedA(0), Some((0, 100, "10M", false)), Some((0, 200, "10M", true))),
            (MoleculeId::PairedB(0), Some((0, 198, "12M", true)), Some((0, 100, "12M", false))),
        ],
        &[
            (0, 1, Some("chr1:101-112:+")), (1, 2, Some("chr1:101-112:+")),
            (0, 2, Some("chr1:199-210:-")), (1, 1, Some("chr1:199-210:-")),
        ],
    )]
    #[case::one_end_unmapped_in_all(
        &[
            (S0, Some((0, 100, "10M", false)), None),
            (S0, Some((0, 100, "11M", false)), None),
        ],
        &[(0, 1, Some("chr1:101-111:+")), (0, 2, None), (1, 2, None)],
    )]
    #[case::two_molecules_same_position(
        &[
            (S0, Some((0, 100, "10M", false)), Some((0, 200, "10M", true))),
            (S1, Some((0, 100, "15M", false)), Some((0, 195, "15M", true))),
        ],
        &[
            (0, 1, Some("chr1:101-110:+")), (1, 1, Some("chr1:101-115:+")),
            (0, 2, Some("chr1:201-210:-")), (1, 2, Some("chr1:196-210:-")),
        ],
    )]
    fn compute_segment_windows_cases(
        #[case] specs: &[Tpl],
        #[case] expected: &[(usize, u8, Option<&str>)],
    ) {
        let ts = templates(specs);
        let windows = compute_segment_windows(&ts, &contigs()).expect("windows");
        for &(idx, which, want) in expected {
            let t = &ts[idx];
            let raw: &[u8] = if which == 1 { t.r1().unwrap() } else { t.r2().unwrap() };
            let got = window_str(&windows, t.mi, raw);
            assert_eq!(got.as_deref(), want, "template {idx} R{which}");
        }
    }

    /// After `rebase`, lookups take the offset molecule ids the templates then
    /// carry; an id below the offset has no window.
    #[rstest]
    #[case::local_id_before_rebase(0, S0, Some("chr1:101-110:+"))]
    #[case::offset_id_after_rebase(7, MoleculeId::Single(7), Some("chr1:101-110:+"))]
    #[case::local_id_after_rebase(7, S0, None)]
    #[case::other_offset_id(7, MoleculeId::Single(8), None)]
    fn rebase_shifts_lookup_ids(
        #[case] base: u64,
        #[case] lookup: MoleculeId,
        #[case] expected: Option<&str>,
    ) {
        let ts = templates(&[(S0, Some((0, 100, "10M", false)), None)]);
        let mut windows = compute_segment_windows(&ts, &contigs()).expect("windows");
        windows.rebase(base);
        let raw: &[u8] = ts[0].r1().unwrap();
        assert_eq!(window_str(&windows, lookup, raw).as_deref(), expected);
    }

    /// `stamp_rw` adds the segment's `rw` to the record being serialized, and
    /// leaves a record without a window unchanged.
    #[rstest]
    #[case::windowed(Some((0, 100, "10M", false)), Some("chr1:101-110:+"))]
    #[case::unmapped(None, None)]
    fn stamp_rw_cases(#[case] rd: Rd, #[case] expected: Option<&str>) {
        let ts = templates(&[(S0, rd, None)]);
        let windows = compute_segment_windows(&ts, &contigs()).expect("windows");
        let raw: &[u8] = ts[0].r1().unwrap();
        let mut scratch = raw.to_vec();
        let mut rw_buf = Vec::new();
        windows.stamp_rw(S0, raw, &mut scratch, &mut rw_buf);
        let got =
            fgumi_raw_bam::find_string_tag(fgumi_raw_bam::aux_data_slice(&scratch), SamTag::RW);
        assert_eq!(got, expected.map(str::as_bytes));
    }

    /// A tid missing from the contig names is an error.
    #[test]
    fn unknown_reference_id_errors() {
        let ts = templates(&[(S0, Some((5, 100, "10M", false)), None)]);
        let err = compute_segment_windows(&ts, &contigs()).expect_err("unknown tid");
        assert!(err.to_string().contains("reference id 5 is not in the header"), "{err}");
    }
}
