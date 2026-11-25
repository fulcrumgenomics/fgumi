//! Half-mapped templates that are copies of a fully mapped pair.
//!
//! A pair with one unmapped read has only one position, so it is in another position group than the
//! pairs it is a copy of. Their ends are kept in [`pair_ends`] until it is reached. `fgumi group`
//! gives it their molecule and puts it there ([`copies`]). It is then in their position group by its
//! `tc` tag if its output is sorted and grouped again, where [`attach`] finds its molecule.

use std::borrow::Cow;

use fgumi_bam_io::GroupKey;
use fgumi_raw_bam::RawRecord;
use fgumi_umi::MoleculeId;

use crate::template::Template;

pub(crate) mod copies;
pub(crate) mod pair_ends;

/// Longest UMI that fits in a [`UmiCode`]; longer UMIs are kept as text.
const MAX_PACKED_BASES: usize = 15;

/// A UMI's bases, two bits each, below a marker bit that gives their number.
type UmiCode = u32;

/// The code of the empty UMI, which stands for every UMI under `--no-umi`.
const NO_UMI: UmiCode = 1;

/// Packs `umi`, or returns `None` if it is too long or has a character other than A, C, G, T.
fn pack(umi: &[u8]) -> Option<UmiCode> {
    if umi.len() > MAX_PACKED_BASES {
        return None;
    }
    let (packed, codes) = umi.iter().fold((1, 0), |(packed, codes), &base| {
        let code = BASE_CODES[usize::from(base)];
        (packed << 2 | UmiCode::from(code & 3), codes | code)
    });
    (codes & NOT_A_BASE == 0).then_some(packed)
}

/// Set in the code of a character other than A, C, G, T.
const NOT_A_BASE: u8 = 4;

/// The code of each character.
const BASE_CODES: [u8; 256] = {
    let mut codes = [NOT_A_BASE; 256];
    let mut code = 0;
    while code < 4 {
        let base = b"ACGT"[code as usize];
        codes[base as usize] = code;
        codes[base.to_ascii_lowercase() as usize] = code;
        code += 1;
    }
    codes
};

fn packed_len(packed: UmiCode) -> u32 {
    packed.ilog2() / 2
}

fn unpack(packed: UmiCode) -> ([u8; MAX_PACKED_BASES], usize) {
    let mut bases = [0u8; MAX_PACKED_BASES];
    let len = packed_len(packed) as usize;
    for (i, base) in bases[..len].iter_mut().enumerate() {
        *base = b"ACGT"[((packed >> (2 * (len - 1 - i))) & 3) as usize];
    }
    (bases, len)
}

/// Which side of a position an end is on: bit 0 for the reverse strand, and bit 1 for R2 (always
/// clear unless ends are matched by segment).
fn side(reverse: bool, is_r2: bool) -> u32 {
    u32::from(reverse) | (u32::from(is_r2) << 1)
}

/// The UMI read by one end of a template.
#[derive(Clone)]
struct EndUmi<'a> {
    bases: Cow<'a, [u8]>,
    packed: Option<UmiCode>,
    /// The template's whole UMI rather than one half of a pair of UMIs.
    whole: bool,
}

impl<'a> EndUmi<'a> {
    fn new(bases: impl Into<Cow<'a, [u8]>>, whole: bool) -> Option<Self> {
        let bases = bases.into();
        let packed = pack(&bases);
        (!bases.is_empty()).then_some(Self { bases, packed, whole })
    }

    fn text(&self) -> Box<[u8]> {
        self.bases.to_ascii_uppercase().into_boxed_slice()
    }

    fn upper(self) -> Cow<'a, [u8]> {
        if self.bases.iter().any(u8::is_ascii_lowercase) {
            Cow::Owned(self.bases.to_ascii_uppercase())
        } else {
            self.bases
        }
    }
}

/// The mapped end of a half-mapped template.
#[derive(Debug)]
pub(crate) struct HalfMapped {
    /// The template's index in its position group.
    pub(crate) index: usize,
    lane: Lane,
    contig: i32,
    pos: i32,
    side: u32,
    /// The code, or the text of a UMI that has none.
    umi: Result<UmiCode, Box<[u8]>>,
    truncate: bool,
}

impl HalfMapped {
    fn umi_text(&self) -> Cow<'_, [u8]> {
        match &self.umi {
            Ok(packed) => {
                let (bases, len) = unpack(*packed);
                Cow::Owned(bases[..len].to_vec())
            }
            Err(text) => Cow::Borrowed(text),
        }
    }
}

/// Reads the mapped ends and their UMIs off templates.
#[derive(Debug, Clone, Copy)]
pub(crate) struct EndReader {
    pub(crate) umi_tag: [u8; 2],
    pub(crate) no_umi: bool,
    /// Match an R1 end only with R1 ends and an R2 end only with R2 ends.
    pub(crate) by_segment: bool,
    /// Compare single UMIs over the length of the shorter one (`--min-umi-length`).
    pub(crate) truncate_umis: bool,
}

impl EndReader {
    pub(crate) fn new(
        umi_tag: [u8; 2],
        no_umi: bool,
        by_segment: bool,
        min_umi_length: Option<usize>,
    ) -> Self {
        Self { umi_tag, no_umi, by_segment, truncate_umis: min_umi_length.is_some() }
    }

    /// The mapped end of the template at `index` of the group keyed `key`, if it is a pair with
    /// one mapped primary read and that end has a UMI.
    pub(crate) fn read_half(
        self,
        template: &Template,
        index: usize,
        key: &GroupKey,
    ) -> Option<HalfMapped> {
        let mapped = |read: Option<&RawRecord>| read.is_some_and(|read| !read.is_unmapped());
        let is_r1 = mapped(template.r1());
        if is_r1 == mapped(template.r2()) {
            return None;
        }
        let read = if is_r1 { template.r1() } else { template.r2() }?;
        if !read.is_paired() || !read.is_mate_unmapped() {
            return None;
        }
        let (r1_umi, r2_umi) = self.end_umis(template)?;
        let umi = if is_r1 { r1_umi } else { r2_umi }?;
        let contig = fgumi_raw_bam::ref_id(read);
        let pos = fgumi_raw_bam::unclipped_5prime_from_raw_bam(read);
        (contig >= 0 && pos != i32::MAX).then(|| HalfMapped {
            index,
            lane: (key.library_idx, key.cell_hash),
            contig,
            pos,
            side: side(read.is_reverse(), self.by_segment && !is_r1),
            truncate: umi.whole && self.truncate_umis,
            umi: umi.packed.ok_or_else(|| umi.text()),
        })
    }

    /// The ends of `template`, of the group keyed `key`, if both its primary reads are mapped: the
    /// UMI, whether it is the key's second position, and whether it is R2 (if matched by segment).
    fn pair_ends<'a>(
        self,
        template: &'a Template,
        key: &GroupKey,
    ) -> [Option<(EndUmi<'a>, bool, bool)>; 2] {
        let (Some(r1), Some(r2)) = (template.r1(), template.r2()) else { return [None, None] };
        if r1.is_unmapped() || r2.is_unmapped() || r1.is_mate_unmapped() || r2.is_mate_unmapped() {
            return [None, None];
        }
        let Some((r1_umi, r2_umi)) = self.end_umis(template) else { return [None, None] };
        let r1_is_second = if key.strand1 == key.strand2 {
            let r1_pos = fgumi_raw_bam::unclipped_5prime_from_raw_bam(r1);
            (fgumi_raw_bam::ref_id(r1), r1_pos) != (key.ref_id1, key.pos1)
        } else {
            u8::from(r1.is_reverse()) != key.strand1
        };
        [
            r1_umi.map(|umi| (umi, r1_is_second, false)),
            r2_umi.map(|umi| (umi, !r1_is_second, self.by_segment)),
        ]
    }

    /// The UMIs of the R1 and R2 ends: the halves of a two-segment UMI, or else the whole UMI.
    /// An empty one is `None`.
    fn end_umis(self, template: &Template) -> Option<(Option<EndUmi<'_>>, Option<EndUmi<'_>>)> {
        if self.no_umi {
            let umi = EndUmi { bases: Cow::Borrowed(b""), packed: Some(NO_UMI), whole: true };
            return Some((Some(umi.clone()), Some(umi)));
        }
        let umi = template.cached_umi().or_else(|| {
            let read = template.r1().or_else(|| template.r2())?;
            fgumi_raw_bam::find_string_tag(fgumi_raw_bam::aux_data_slice(read), self.umi_tag)
        })?;
        if let Some(dash) = umi.iter().position(|&b| b == b'-') {
            let (r1_umi, r2_umi) = (&umi[..dash], &umi[dash + 1..]);
            if !r2_umi.contains(&b'-') {
                return Some((EndUmi::new(r1_umi, false), EndUmi::new(r2_umi, false)));
            }
        }
        let whole = if self.truncate_umis {
            let is_base = |b: &u8| matches!(b.to_ascii_uppercase(), b'A' | b'C' | b'G' | b'T');
            EndUmi::new(umi.iter().copied().filter(is_base).collect::<Vec<u8>>(), true)
        } else {
            EndUmi::new(umi, true)
        };
        Some((whole.clone(), whole))
    }
}

/// Whether exactly one primary read of a pair is mapped, as its own flags say.
pub(crate) fn is_half_mapped(template: &Template) -> bool {
    let mapped = |read: Option<&RawRecord>| {
        read.is_some_and(|read| !read.is_unmapped() && read.is_paired() && read.is_mate_unmapped())
    };
    let unmapped = |read: Option<&RawRecord>| read.is_none_or(RawRecord::is_unmapped);
    (mapped(template.r1()) && unmapped(template.r2()))
        || (mapped(template.r2()) && unmapped(template.r1()))
}

/// Removes and returns the half-mapped templates that are in a group of pairs by their `tc` tag.
pub(crate) fn split_off_placed(templates: &mut Vec<Template>, key: &GroupKey) -> Vec<Template> {
    if !key.has_mate_position() || !templates.iter().any(is_half_mapped) {
        return Vec::new();
    }
    let (placed, rest) = std::mem::take(templates).into_iter().partition(is_half_mapped);
    *templates = rest;
    placed
}

/// The fully mapped pairs of its position group that a half-mapped template matches.
#[derive(Debug, Clone, Copy, Default)]
pub(crate) struct Placement {
    /// Whether the template has been given their molecule, which is the only one as close.
    pub(crate) attached: bool,
    /// Whether one of them has the same end as the same segment, so that the two would be in one
    /// family if the template's other read were mapped.
    pub(crate) same_strand: bool,
}

/// Gives each of `placed` the molecule of the fully mapped pairs in `templates` whose UMI at the
/// same end is closest. A template as close to two molecules as to any gets neither. The molecules
/// of `templates` do not change.
pub(crate) fn attach(
    templates: &[Template],
    placed: &mut [Template],
    key: &GroupKey,
    reader: EndReader,
    edits: u32,
) -> Vec<Placement> {
    let paired = templates.iter().any(|template| template.mi.is_paired());
    let segments = EndReader { by_segment: true, ..reader };
    let mut ends: Vec<(bool, bool, Cow<'_, [u8]>, u64)> = Vec::new();
    for template in templates {
        let Some(id) = template.mi.id() else { continue };
        for (umi, second, is_r2) in segments.pair_ends(template, key).into_iter().flatten() {
            ends.push((second, is_r2, umi.upper(), id));
        }
    }
    ends.sort_unstable();
    ends.dedup();

    let place = |template: &mut Template| {
        let mut placement = Placement::default();
        let Some(half) = segments.read_half(template, 0, key) else { return placement };
        let end = (half.contig, half.pos, half.side & 1);
        let second = end != (key.ref_id1, key.pos1, u32::from(key.strand1));
        if second && end != (key.ref_id2, key.pos2, u32::from(key.strand2)) {
            return placement;
        }
        let is_r2 = half.side & 2 != 0;
        let umi = half.umi_text();
        let mut closest: Option<(usize, u64)> = None;
        let mut tied = false;
        for (known_second, known_is_r2, known, id) in &ends {
            if *known_second != second || (known.len() != umi.len() && !half.truncate) {
                continue;
            }
            let distance = known.iter().zip(umi.iter()).filter(|(a, b)| a != b).count();
            if distance > edits as usize {
                continue;
            }
            placement.same_strand |= *known_is_r2 == is_r2;
            if reader.by_segment && *known_is_r2 != is_r2 {
                continue;
            }
            match closest {
                Some((best, _)) if distance > best => {}
                Some((best, closest_id)) if distance == best => tied |= closest_id != *id,
                _ => {
                    closest = Some((distance, *id));
                    tied = false;
                }
            }
        }
        if let Some((_, id)) = closest.filter(|_| !tied) {
            template.mi = match (paired, is_r2 == second) {
                (false, _) => MoleculeId::Single(id),
                (true, true) => MoleculeId::PairedA(id),
                (true, false) => MoleculeId::PairedB(id),
            };
            placement.attached = true;
        }
        placement
    };
    placed.iter_mut().map(place).collect()
}

/// A library and cell barcode hash, as in [`GroupKey`].
type Lane = (u16, u64);
