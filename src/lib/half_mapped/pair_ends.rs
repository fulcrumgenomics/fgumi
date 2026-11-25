//! Finds the fully mapped pairs with the same end as a half-mapped template.
//!
//! Template-coordinate order compares `tid1`, `tid2`, then `pos1`, and a half-mapped template
//! has the largest `tid2`. It therefore comes after every pair with an end on its contig, so one
//! pass over the input is enough: [`PairEnds`] collects the ends of the pairs as they go by, and
//! releases a contig's ends once the input has moved past that contig. An end can carry its
//! molecule ([`MoleculeEnd`]) or nothing (`()`).

use std::ops::RangeInclusive;

use ahash::AHashMap;
use fgumi_bam_io::GroupKey;

use super::{EndReader, HalfMapped, Lane, UmiCode, packed_len, side, unpack};
use crate::template::Template;

/// Fewest entries allocated for one chunk of a contig's ends.
const MIN_CHUNK_ENTRIES: usize = 1 << 6;

/// Most entries held in one chunk of a contig's ends.
const MAX_CHUNK_ENTRIES: usize = 1 << 16;

/// Most times a chunk is searched before it is sorted.
const MAX_UNSORTED_SEARCHES: u32 = 64;

/// Set in a [`UmiCode`] when the other bits index a list of texts.
const TEXT_UMI: UmiCode = 1 << 31;

/// How two UMIs are compared.
#[derive(Debug, Clone, Copy)]
struct UmiTolerance {
    edits: u32,
    /// Compare only as many bases as the shorter UMI has.
    truncate: bool,
}

impl UmiTolerance {
    /// The distance of two UMIs, if it is tolerated.
    fn packed_distance(self, a: UmiCode, b: UmiCode) -> Option<u32> {
        let (len_a, len_b) = (packed_len(a), packed_len(b));
        if len_a != len_b && !self.truncate {
            return None;
        }
        let len = len_a.min(len_b);
        let diff = (a >> (2 * (len_a - len))) ^ (b >> (2 * (len_b - len)));
        let distance = ((diff | (diff >> 1)) & 0x5555_5555).count_ones();
        (distance <= self.edits).then_some(distance)
    }

    fn text_distance(self, a: &[u8], b: &[u8]) -> Option<u32> {
        if a.len() != b.len() && !self.truncate {
            return None;
        }
        let distance = u32::try_from(a.iter().zip(b).filter(|(x, y)| x != y).count()).ok()?;
        (distance <= self.edits).then_some(distance)
    }
}

/// An end of a fully mapped pair, placed by its group's key.
#[derive(Debug, Clone, Copy, PartialEq, Eq, PartialOrd, Ord)]
struct PairEnd<M> {
    /// The second position of the group key rather than the first.
    second: bool,
    is_r2: bool,
    umi: UmiCode,
    /// The molecule's id in its position group, or `()`.
    molecule: M,
}

/// What an end's entry keeps of its molecule.
pub(crate) trait EndMolecule: Copy + Ord + std::fmt::Debug {
    /// The molecule of `template`, by its id in its position group, if it has one.
    fn in_group(template: &Template) -> Option<Self>;

    /// The molecule among all those of a contig, whose first was `contig_base`, when the ids of
    /// its position group start at `group_base`. `None` if it does not fit.
    fn in_contig(self, group_base: u64, contig_base: u64, second: bool) -> Option<Self>;
}

impl EndMolecule for () {
    fn in_group(_template: &Template) -> Option<Self> {
        Some(())
    }

    fn in_contig(self, _group_base: u64, _contig_base: u64, _second: bool) -> Option<Self> {
        Some(())
    }
}

/// A molecule's id, less the contig's first one, and whether the end is its second.
#[derive(Debug, Clone, Copy, PartialEq, Eq, PartialOrd, Ord)]
pub(crate) struct MoleculeEnd(u32);

impl EndMolecule for MoleculeEnd {
    fn in_group(template: &Template) -> Option<Self> {
        u32::try_from(template.mi.id()?).ok().map(Self)
    }

    fn in_contig(self, group_base: u64, contig_base: u64, second: bool) -> Option<Self> {
        let id = u32::try_from(group_base + u64::from(self.0) - contig_base).ok()?;
        (id < 1 << 31).then_some(Self(id << 1 | u32::from(second)))
    }
}

/// The ends of the fully mapped pairs in a batch of position groups.
#[derive(Debug)]
pub(crate) struct PairEndBuffer<M> {
    /// The distinct ends of each of `groups` in turn.
    ends: Vec<PairEnd<M>>,
    /// The keys of the groups that have pair ends, their indices, and how many.
    groups: Vec<(GroupKey, usize, usize)>,
    texts: Vec<Box<[u8]>>,
}

impl<M> Default for PairEndBuffer<M> {
    fn default() -> Self {
        Self { ends: Vec::new(), groups: Vec::new(), texts: Vec::new() }
    }
}

impl<M: EndMolecule> PairEndBuffer<M> {
    /// Adds the ends of the pairs among `templates`, of the group at `index` of its batch, that
    /// have both primary reads mapped.
    pub(crate) fn add_group(
        &mut self,
        reader: EndReader,
        key: &GroupKey,
        index: usize,
        templates: &[Template],
    ) {
        if !key.has_mate_position() || key.ref_id1 < 0 {
            return;
        }
        let start = self.ends.len();
        for template in templates {
            self.add_pair(reader, key, template);
        }
        if self.ends.len() - start > 2 {
            self.ends[start..].sort_unstable();
            let mut unique = start + 1;
            for next in start + 1..self.ends.len() {
                if self.ends[next] != self.ends[unique - 1] {
                    self.ends[unique] = self.ends[next];
                    unique += 1;
                }
            }
            self.ends.truncate(unique);
        }
        if self.ends.len() > start {
            self.groups.push((*key, index, self.ends.len() - start));
        }
    }

    fn add_pair(&mut self, reader: EndReader, key: &GroupKey, template: &Template) {
        let Some(molecule) = M::in_group(template) else { return };
        for (umi, second, is_r2) in reader.pair_ends(template, key).into_iter().flatten() {
            let umi = umi.packed.unwrap_or_else(|| {
                self.texts.push(umi.text());
                TEXT_UMI | u32::try_from(self.texts.len() - 1).expect("too many UMIs in a batch")
            });
            self.ends.push(PairEnd { second, is_r2, umi, molecule });
        }
    }

    /// Reads the ends from the first group on.
    pub(crate) fn reader(&self) -> PairEndReader<'_, M> {
        PairEndReader { ends: &self.ends, groups: &self.groups, texts: &self.texts }
    }

    pub(crate) fn heap_size(&self) -> usize {
        self.ends.capacity() * std::mem::size_of::<PairEnd<M>>()
            + self.groups.capacity() * std::mem::size_of::<(GroupKey, usize, usize)>()
            + self.texts.iter().map(|text| text.len()).sum::<usize>()
            + self.texts.capacity() * std::mem::size_of::<Box<[u8]>>()
    }
}

/// The ends in a [`PairEndBuffer`] that have not been added to a [`PairEnds`].
pub(crate) struct PairEndReader<'a, M> {
    ends: &'a [PairEnd<M>],
    groups: &'a [(GroupKey, usize, usize)],
    texts: &'a [Box<[u8]>],
}

/// One pair end in 12 bytes, or 16 with its molecule.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
struct Entry<M> {
    pos: u32,
    lane_and_side: u32,
    umi: UmiCode,
    molecule: M,
}

impl<M> Entry<M> {
    fn place(&self) -> u64 {
        place(self.pos, self.lane_and_side)
    }
}

fn place(pos: u32, lane_and_side: u32) -> u64 {
    u64::from(pos) << 32 | u64::from(lane_and_side)
}

/// Ends that were added one after another, so most are near each other.
#[derive(Debug)]
struct Chunk<M> {
    entries: Vec<Entry<M>>,
    /// The lowest and highest [`position_bits`].
    range: RangeInclusive<u32>,
    is_sorted: bool,
    /// How many times the chunk was searched.
    searches: u32,
}

impl<M: EndMolecule> Chunk<M> {
    fn new(capacity: usize) -> Self {
        let range = RangeInclusive::new(u32::MAX, 0);
        Self { entries: Vec::with_capacity(capacity), range, is_sorted: false, searches: 0 }
    }

    fn push(&mut self, entry: Entry<M>) {
        let (low, high) = (*self.range.start(), *self.range.end());
        self.range = low.min(entry.pos)..=high.max(entry.pos);
        self.entries.push(entry);
    }

    /// Counts a search, and sorts the entries once there have been enough for that to pay.
    fn prepare(&mut self) {
        self.searches += 1;
        if !self.is_sorted && self.searches > MAX_UNSORTED_SEARCHES {
            self.entries.sort_unstable_by_key(|e| (e.place(), e.umi, e.molecule));
            self.entries.dedup();
            self.is_sorted = true;
        }
    }

    /// The entries at `place` with the UMI `exact`, if any, or else with any UMI.
    fn at(&self, place: u64, exact: Option<UmiCode>) -> impl Iterator<Item = &Entry<M>> {
        let from = (place, exact.unwrap_or(0));
        let start = if self.is_sorted {
            self.entries.partition_point(|entry| (entry.place(), entry.umi) < from)
        } else {
            0
        };
        let is_match = move |entry: &&Entry<M>| {
            entry.place() == place && exact.is_none_or(|umi| umi == entry.umi)
        };
        let is_sorted = self.is_sorted;
        self.entries[start..]
            .iter()
            .take_while(move |entry| !is_sorted || is_match(entry))
            .filter(is_match)
    }
}

/// The distinct ends of the fully mapped pairs on one contig.
#[derive(Debug)]
struct ContigEnds<M> {
    chunks: Vec<Chunk<M>>,
    len: usize,
    /// Set once the contig's half-mapped templates have been reached, after which no ends are added.
    sealed: bool,
    /// The id of the first molecule added.
    base: u64,
    lanes: AHashMap<Lane, u32>,
    last_lane: Option<(Lane, u32)>,
    texts: Vec<Box<[u8]>>,
    text_codes: AHashMap<Box<[u8]>, UmiCode>,
}

impl<M: EndMolecule> ContigEnds<M> {
    fn new(base: u64) -> Self {
        Self {
            chunks: Vec::new(),
            len: 0,
            sealed: false,
            base,
            lanes: AHashMap::new(),
            last_lane: None,
            texts: Vec::new(),
            text_codes: AHashMap::new(),
        }
    }

    fn insert(&mut self, pos: i32, lane_and_side: u32, umi: UmiCode, molecule: M) {
        let is_full = |chunk: &Chunk<M>| chunk.entries.len() == chunk.entries.capacity();
        if self.chunks.last().is_none_or(is_full) {
            let capacity = self.len.clamp(MIN_CHUNK_ENTRIES, MAX_CHUNK_ENTRIES);
            self.chunks.push(Chunk::new(capacity));
        }
        let entry = Entry { pos: position_bits(pos), lane_and_side, umi, molecule };
        self.chunks.last_mut().expect("a chunk with room was just ensured").push(entry);
        self.len += 1;
    }

    /// The number of `lane`, shifted to leave room for an end's side.
    fn lane_number(&mut self, lane: Lane) -> Option<u32> {
        if let Some((last, number)) = self.last_lane
            && last == lane
        {
            return Some(number);
        }
        let next = u32::try_from(self.lanes.len()).ok().filter(|&n| n < 1 << 30)?;
        let number = *self.lanes.entry(lane).or_insert(next << 2);
        self.last_lane = Some((lane, number));
        Some(number)
    }

    fn text_code(&mut self, text: &[u8]) -> Option<UmiCode> {
        if let Some(&code) = self.text_codes.get(text) {
            return Some(code);
        }
        let code = TEXT_UMI | u32::try_from(self.texts.len()).ok().filter(|&n| n < TEXT_UMI)?;
        self.texts.push(text.into());
        self.text_codes.insert(text.into(), code);
        Some(code)
    }

    /// The ends that are the same as that of `half`, with a UMI within `edits` of it, and the
    /// number of bases that differ.
    fn matches<'a>(
        &'a mut self,
        half: &'a HalfMapped,
        edits: u32,
    ) -> impl Iterator<Item = (u32, M)> + 'a {
        self.sealed = true;
        let pos = position_bits(half.pos);
        for chunk in self.chunks.iter_mut().filter(|chunk| chunk.range.contains(&pos)) {
            chunk.prepare();
        }
        let place = self.lanes.get(&half.lane).map(|lane| place(pos, lane | half.side));
        let tolerance = UmiTolerance { edits, truncate: half.truncate };
        let exact = half.umi.as_ref().ok().copied().filter(|_| edits == 0 && !half.truncate);
        let ends: &'a Self = self;
        place.into_iter().flat_map(move |place| {
            let chunks = ends.chunks.iter().filter(move |chunk| chunk.range.contains(&pos));
            chunks.flat_map(move |chunk| chunk.at(place, exact)).filter_map(move |entry| {
                let distance = ends.umi_distance(&half.umi, entry.umi, tolerance)?;
                Some((distance, entry.molecule))
            })
        })
    }

    /// How many bases of two UMIs differ, if they match.
    fn umi_distance(
        &self,
        umi: &Result<UmiCode, Box<[u8]>>,
        known: UmiCode,
        tolerance: UmiTolerance,
    ) -> Option<u32> {
        let known_text =
            (known & TEXT_UMI != 0).then(|| self.texts[(known & !TEXT_UMI) as usize].as_ref());
        match (umi, known_text) {
            (Ok(packed), None) => tolerance.packed_distance(*packed, known),
            (Ok(packed), Some(known)) => {
                let (bases, len) = unpack(*packed);
                tolerance.text_distance(&bases[..len], known)
            }
            (Err(text), None) => {
                let (bases, len) = unpack(known);
                tolerance.text_distance(text, &bases[..len])
            }
            (Err(text), Some(known)) => tolerance.text_distance(text, known),
        }
    }
}

#[allow(clippy::cast_sign_loss)]
fn position_bits(pos: i32) -> u32 {
    pos as u32
}

/// The ends of the fully mapped pairs seen so far, by contig.
///
/// Must see every batch in input order.
#[derive(Debug)]
pub(crate) struct PairEnds<M> {
    contigs: Vec<(i32, ContigEnds<M>)>,
    /// The index in `contigs` of each contig.
    index: AHashMap<i32, usize>,
    /// The contig looked up last.
    last: Option<(i32, usize)>,
    /// The [`contig_rank`] the input has reached.
    rank: u16,
}

impl<M> Default for PairEnds<M> {
    fn default() -> Self {
        Self { contigs: Vec::new(), index: AHashMap::new(), last: None, rank: 0 }
    }
}

impl<M: EndMolecule> PairEnds<M> {
    /// The ends on `contig`, whose first molecule is `base` if there are none yet.
    fn contig(&mut self, contig: i32, base: u64) -> &mut ContigEnds<M> {
        let index = match self.last {
            Some((last, index)) if last == contig => index,
            _ => *self.index.entry(contig).or_insert_with(|| {
                self.contigs.push((contig, ContigEnds::new(base)));
                self.contigs.len() - 1
            }),
        };
        self.last = Some((contig, index));
        &mut self.contigs[index].1
    }

    /// Adds the ends of a batch's groups before `index`. `base` is the first molecule id of the
    /// group at an index.
    pub(crate) fn add(
        &mut self,
        unread: &mut PairEndReader<'_, M>,
        index: usize,
        base: impl Fn(usize) -> u64,
    ) {
        while let Some(((key, group, count), groups)) =
            unread.groups.split_first().filter(|((_, group, _), _)| *group < index)
        {
            let (group_ends, ends) = unread.ends.split_at(*count);
            (unread.groups, unread.ends) = (groups, ends);
            let base = base(*group);
            for (second, contig, pos, strand) in [
                (false, key.ref_id1, key.pos1, key.strand1),
                (true, key.ref_id2, key.pos2, key.strand2),
            ] {
                let known = self.contig(contig, base);
                if known.sealed {
                    continue;
                }
                let lane = (key.library_idx, key.cell_hash);
                let Some(lane) = known.lane_number(lane) else { continue };
                for end in group_ends.iter().filter(|end| end.second == second) {
                    let umi = if end.umi & TEXT_UMI == 0 {
                        Some(end.umi)
                    } else {
                        known.text_code(&unread.texts[(end.umi & !TEXT_UMI) as usize])
                    };
                    let molecule = end.molecule.in_contig(base, known.base, second);
                    if let (Some(umi), Some(molecule)) = (umi, molecule) {
                        known.insert(pos, lane | side(strand != 0, end.is_r2), umi, molecule);
                    }
                }
            }
        }
    }

    /// Releases the contigs before the last of the batch, once its half-mapped templates have
    /// been looked up.
    pub(crate) fn release(&mut self, buffer: &PairEndBuffer<M>) {
        let Some((key, ..)) = buffer.groups.last() else { return };
        let rank = contig_rank(key.ref_id1);
        if rank != self.rank {
            self.rank = rank;
            self.contigs.retain(|(contig, _)| contig_rank(*contig) >= rank);
            let indices = self.contigs.iter().enumerate();
            self.index = indices.map(|(index, (contig, _))| (*contig, index)).collect();
            self.last = None;
        }
    }

    /// Whether a fully mapped pair has the same end as `half`, with a UMI within `edits` of it.
    pub(crate) fn contains(&mut self, half: &HalfMapped, edits: u32) -> bool {
        self.contig(half.contig, 0).matches(half, edits).next().is_some()
    }
}

impl PairEnds<MoleculeEnd> {
    /// The molecule with the same end as `half` and the closest UMI, within `edits` of it, and
    /// whether that end is the molecule's second. `None` if several molecules are as close.
    pub(crate) fn closest(&mut self, half: &HalfMapped, edits: u32) -> Option<(u64, bool)> {
        let known = self.contig(half.contig, 0);
        let base = known.base;
        let (mut closest, mut tied) = (None, false);
        for (distance, MoleculeEnd(end)) in known.matches(half, edits) {
            match closest {
                Some((best, _)) if distance > best => {}
                Some((best, first)) if distance == best => tied |= first >> 1 != end >> 1,
                _ => (closest, tied) = (Some((distance, end)), false),
            }
        }
        let (_, end) = closest.filter(|_| !tied)?;
        Some((base + u64::from(end >> 1), end & 1 == 1))
    }
}

/// The contig's place in template-coordinate order, whose key keeps the low 16 bits.
#[allow(clippy::cast_sign_loss, clippy::cast_possible_truncation)]
fn contig_rank(contig: i32) -> u16 {
    contig as u16
}

#[cfg(test)]
mod tests {
    use super::super::pack;
    use super::*;
    use rstest::rstest;

    fn pair_key(ref_id1: i32, pos1: i32, ref_id2: i32, pos2: i32) -> GroupKey {
        GroupKey::paired(ref_id1, pos1, 0, ref_id2, pos2, 1, 0, 0, 0)
    }

    fn end(second: bool, umi: &[u8]) -> PairEnd<()> {
        PairEnd { second, is_r2: false, umi: pack(umi).expect("a packable UMI"), molecule: () }
    }

    fn half(contig: i32, pos: i32, reverse: bool, umi: &[u8]) -> HalfMapped {
        HalfMapped {
            index: 0,
            lane: (0, 0),
            contig,
            pos,
            side: side(reverse, false),
            umi: pack(umi).ok_or_else(|| umi.into()),
            truncate: false,
        }
    }

    /// A batch of one group whose pairs have the UMI `ACGT` at the first end of `key` and `TTTT`
    /// at the second.
    fn batch(key: GroupKey) -> PairEndBuffer<()> {
        let ends = vec![end(false, b"ACGT"), end(true, b"TTTT")];
        PairEndBuffer { ends, groups: vec![(key, 0, 2)], texts: Vec::new() }
    }

    fn add(known: &mut PairEnds<()>, buffer: &PairEndBuffer<()>) {
        known.add(&mut buffer.reader(), usize::MAX, |_| 0);
    }

    #[rstest]
    #[case::empty(b"")]
    #[case::all_bases(b"ACGT")]
    #[case::leading_a(b"AAC")]
    #[case::longest(b"TGCATGCATGCATGC")]
    fn packed_umis_unpack(#[case] umi: &[u8]) {
        let (bases, len) = unpack(pack(umi).unwrap());
        assert_eq!(&bases[..len], umi);
    }

    #[rstest]
    #[case::too_long(b"ACGTACGTACGTACGT")]
    #[case::not_a_base(b"ACRT")]
    fn other_umis_do_not_pack(#[case] umi: &[u8]) {
        assert_eq!(pack(umi), None);
    }

    /// Packed and text comparisons agree, whichever of the two a UMI is held as.
    #[rstest]
    #[case::same(b"ACGT", b"ACGT", 0, false, true)]
    #[case::every_base_differs(b"AAAA", b"TTTT", 3, false, false)]
    #[case::within_edits(b"AAAA", b"CATA", 2, false, true)]
    #[case::beyond_edits(b"AAAA", b"CATA", 1, false, false)]
    #[case::longer(b"ACG", b"ACGA", 1, false, false)]
    #[case::truncated(b"ACG", b"ACGA", 0, true, true)]
    #[case::truncated_mismatch(b"ACT", b"ACGA", 0, true, false)]
    #[case::leading_a_counts(b"AACG", b"ACG", 0, false, false)]
    fn umis_match(
        #[case] a: &[u8],
        #[case] b: &[u8],
        #[case] edits: u32,
        #[case] truncate: bool,
        #[case] expected: bool,
    ) {
        let tolerance = UmiTolerance { edits, truncate };
        let (packed_a, packed_b) = (pack(a).unwrap(), pack(b).unwrap());
        assert_eq!(tolerance.packed_distance(packed_a, packed_b).is_some(), expected);
        assert_eq!(tolerance.packed_distance(packed_b, packed_a).is_some(), expected);
        assert_eq!(tolerance.text_distance(a, b).is_some(), expected);
        assert_eq!(tolerance.text_distance(b, a).is_some(), expected);
    }

    #[test]
    fn ends_are_found_in_every_chunk() {
        let mut known = PairEnds::default();
        let count = i32::try_from(3 * MAX_CHUNK_ENTRIES / 2).unwrap();
        for pos in (0..count).rev() {
            add(&mut known, &batch(pair_key(0, pos, 0, pos + 1_000_000)));
        }
        assert!(known.contig(0, 0).chunks.len() > 2);
        for edits in [0, 1] {
            for pos in 0..count {
                assert!(known.contains(&half(0, pos, false, b"ACGT"), edits));
                assert!(known.contains(&half(0, pos + 1_000_000, true, b"TTTT"), edits));
                assert!(!known.contains(&half(0, pos, true, b"ACGT"), edits));
                assert!(!known.contains(&half(0, pos, false, b"TTTT"), edits));
                for lane in [(1, 0), (0, 1)] {
                    let half = HalfMapped { lane, ..half(0, pos, false, b"ACGT") };
                    assert!(!known.contains(&half, edits));
                }
            }
            assert!(!known.contains(&half(0, count, false, b"ACGT"), edits));
        }
    }

    #[test]
    fn a_contig_is_released_once_the_input_has_moved_past_it() {
        let mut known = PairEnds::default();
        for key in [pair_key(0, 100, 0, 300), pair_key(0, 100, 2, 300)] {
            add(&mut known, &batch(key));
            known.release(&batch(key));
        }
        assert!(known.contains(&half(0, 100, false, b"ACGT"), 0));

        let next = batch(pair_key(1, 100, 1, 300));
        add(&mut known, &next);
        assert!(known.contains(&half(0, 100, false, b"ACGT"), 0), "until the batch is released");
        known.release(&next);
        let mut contigs: Vec<_> = known.contigs.iter().map(|(contig, _)| *contig).collect();
        contigs.sort_unstable();
        assert_eq!(contigs, [1, 2]);
        assert!(known.contains(&half(2, 300, true, b"TTTT"), 0));
    }

    /// The sort key keeps 16 bits of a contig index, so contigs that share them are interleaved.
    #[test]
    fn contigs_that_sort_together_are_released_together() {
        let mut known = PairEnds::default();
        let halves = || [half(0, 100, false, b"ACGT"), half(1 << 16, 200, false, b"ACGT")];
        for half in halves() {
            let pairs = batch(pair_key(half.contig, half.pos, half.contig, 300));
            add(&mut known, &pairs);
            known.release(&pairs);
        }
        for half in halves() {
            assert!(known.contains(&half, 0));
            let elsewhere = HalfMapped { contig: (1 << 16) - half.contig, ..half };
            assert!(!known.contains(&elsewhere, 0));
        }

        known.release(&batch(pair_key(1, 100, 1, 300)));
        assert!(known.contigs.iter().all(|(contig, _)| *contig == 1));
    }

    /// Input in template-coordinate order has no pairs on a contig after its half-mapped templates.
    #[test]
    fn ends_after_a_contigs_half_mapped_templates_are_ignored() {
        let mut known = PairEnds::default();
        assert!(!known.contains(&half(0, 100, false, b"ACGT"), 0));
        add(&mut known, &batch(pair_key(0, 100, 0, 300)));
        assert!(!known.contains(&half(0, 100, false, b"ACGT"), 0));
        assert_eq!(known.contig(0, 0).len, 0);
    }

    /// Each end is at position 100 of contig 1 and the first end of its molecule, or the second.
    fn molecules(ends: &[(u64, bool, &str)]) -> PairEnds<MoleculeEnd> {
        let mut known = PairEnds::default();
        for (index, &(id, second, umi)) in ends.iter().enumerate() {
            let key = if second {
                GroupKey::paired(0, 50, 0, 1, 100, 0, 0, 0, 0)
            } else {
                GroupKey::paired(1, 100, 0, 2, 50, 0, 0, 0, 0)
            };
            // The group's ids start below the molecule's.
            let molecule = MoleculeEnd(3);
            let umi = pack(umi.as_bytes()).unwrap();
            let end = PairEnd { second, is_r2: false, umi, molecule };
            let buffer =
                PairEndBuffer { ends: vec![end], groups: vec![(key, index, 1)], texts: vec![] };
            known.add(&mut buffer.reader(), usize::MAX, |group| {
                assert_eq!(group, index);
                id - 3
            });
        }
        known
    }

    #[rstest]
    #[case::none(&[], None)]
    #[case::one(&[(7, false, "ACGT")], Some((7, false)))]
    #[case::second_end(&[(7, true, "ACGT")], Some((7, true)))]
    #[case::far_ids(&[(7, true, "AAAA"), ((1 << 31) + 3, false, "ACGT")], Some(((1 << 31) + 3, false)))]
    #[case::too_far_ids(&[(7, true, "AAAA"), ((1 << 31) + 4, false, "ACGT")], None)]
    #[case::closer_first(&[(7, false, "ACGT"), (8, false, "ACGA")], Some((7, false)))]
    #[case::closer_last(&[(7, false, "ACGA"), (8, false, "ACGT")], Some((8, false)))]
    #[case::one_mismatch(&[(7, false, "ACGA")], Some((7, false)))]
    #[case::two_mismatches(&[(7, false, "ACAA")], None)]
    #[case::tie(&[(7, false, "ACGA"), (8, true, "ACGC")], None)]
    #[case::same_molecule_twice(&[(7, false, "ACGA"), (7, false, "ACGC")], Some((7, false)))]
    #[case::tie_beyond_the_closest(
        &[(7, false, "ACGA"), (8, false, "ACGC"), (9, false, "ACGT")],
        Some((9, false))
    )]
    fn the_closest_molecule_is_found(
        #[case] ends: &[(u64, bool, &str)],
        #[case] expected: Option<(u64, bool)>,
    ) {
        assert_eq!(molecules(ends).closest(&half(1, 100, false, b"ACGT"), 1), expected);
    }

    /// The ends of a group are added once the groups before have their ids.
    #[test]
    fn groups_are_added_up_to_an_index() {
        let key = |pos| GroupKey::paired(0, pos, 0, 0, pos + 500, 1, 0, 0, 0);
        let end = |molecule| PairEnd {
            second: false,
            is_r2: false,
            umi: pack(b"ACGT").unwrap(),
            molecule: MoleculeEnd(molecule),
        };
        let buffer = PairEndBuffer {
            ends: vec![end(0), end(1), end(0)],
            groups: vec![(key(100), 0, 2), (key(200), 5, 1)],
            texts: vec![],
        };
        let mut known = PairEnds::default();
        let mut reader = buffer.reader();
        known.add(&mut reader, 5, |group| [10][group]);
        assert_eq!(known.contig(0, 0).len, 2);
        known.add(&mut reader, usize::MAX, |group| [10, 0, 0, 0, 0, 20][group]);
        assert_eq!(known.closest(&half(0, 200, false, b"ACGT"), 0), Some((20, false)));
        assert_eq!(known.closest(&half(0, 100, false, b"ACGT"), 0), None, "two molecules");
    }
}
