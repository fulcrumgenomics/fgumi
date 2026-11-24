//! Finds the fully mapped pairs of half-mapped templates that `fgumi sort` did not place.
//!
//! Template-coordinate order compares `tid1`, `tid2`, then `pos1`, and a half-mapped template
//! has the largest `tid2`. It therefore comes after every pair with an end on its contig, so one
//! pass over the input is enough: [`PairEnds`] collects the ends of the pairs as they go by, and
//! releases a contig's ends once the input has moved past that contig.

use ahash::AHashMap;
use fgumi_bam_io::GroupKey;

use super::{EndReader, HalfMapped, Lane, UmiCode, packed_len, side, unpack};
use crate::template::Template;

/// Fewest entries allocated for one chunk of a contig's ends.
const MIN_CHUNK_ENTRIES: usize = 1 << 6;

/// Most entries held in one chunk of a contig's ends (48 MiB).
const MAX_CHUNK_ENTRIES: usize = 1 << 22;

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
    fn packed_match(self, a: UmiCode, b: UmiCode) -> bool {
        let (len_a, len_b) = (packed_len(a), packed_len(b));
        if len_a != len_b && !self.truncate {
            return false;
        }
        let len = len_a.min(len_b);
        let diff = (a >> (2 * (len_a - len))) ^ (b >> (2 * (len_b - len)));
        ((diff | (diff >> 1)) & 0x5555_5555).count_ones() <= self.edits
    }

    fn text_match(self, a: &[u8], b: &[u8]) -> bool {
        if a.len() != b.len() && !self.truncate {
            return false;
        }
        a.iter().zip(b).filter(|(x, y)| x != y).count() <= self.edits as usize
    }
}

/// An end of a fully mapped pair, placed by its group's key.
#[derive(Debug, Clone, Copy, PartialEq, Eq, PartialOrd, Ord)]
struct PairEnd {
    /// The second position of the group key rather than the first.
    second: bool,
    is_r2: bool,
    umi: UmiCode,
}

/// The ends of the fully mapped pairs in a batch of position groups.
#[derive(Debug, Default)]
pub(crate) struct PairEndBuffer {
    /// The distinct ends of each of `groups` in turn.
    ends: Vec<PairEnd>,
    /// The keys of the groups that have pair ends, and how many.
    groups: Vec<(GroupKey, usize)>,
    texts: Vec<Box<[u8]>>,
}

impl PairEndBuffer {
    /// Adds the ends of the pairs among `templates` that have both primary reads mapped.
    pub(crate) fn add_group(&mut self, reader: EndReader, key: &GroupKey, templates: &[Template]) {
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
            self.groups.push((*key, self.ends.len() - start));
        }
    }

    fn add_pair(&mut self, reader: EndReader, key: &GroupKey, template: &Template) {
        for (umi, second, is_r2) in reader.pair_ends(template, key).into_iter().flatten() {
            let umi = umi.packed.unwrap_or_else(|| {
                self.texts.push(umi.text());
                TEXT_UMI | u32::try_from(self.texts.len() - 1).expect("too many UMIs in a batch")
            });
            self.ends.push(PairEnd { second, is_r2, umi });
        }
    }

    pub(crate) fn heap_size(&self) -> usize {
        self.ends.capacity() * std::mem::size_of::<PairEnd>()
            + self.groups.capacity() * std::mem::size_of::<(GroupKey, usize)>()
            + self.texts.iter().map(|text| text.len()).sum::<usize>()
            + self.texts.capacity() * std::mem::size_of::<Box<[u8]>>()
    }
}

/// One pair end in 12 bytes; field order is the search order.
#[derive(Debug, Clone, Copy, PartialEq, Eq, PartialOrd, Ord)]
struct Entry {
    pos: u32,
    lane_and_side: u32,
    umi: UmiCode,
}

/// The distinct ends of the fully mapped pairs on one contig.
#[derive(Debug, Default)]
struct ContigEnds {
    /// Each chunk is sorted once [`Self::seal`] has run.
    chunks: Vec<Vec<Entry>>,
    len: usize,
    lanes: AHashMap<Lane, u32>,
    last_lane: Option<(Lane, u32)>,
    texts: Vec<Box<[u8]>>,
    text_codes: AHashMap<Box<[u8]>, UmiCode>,
}

impl ContigEnds {
    fn insert(&mut self, pos: i32, lane_and_side: u32, umi: UmiCode) {
        if self.chunks.last().is_none_or(|chunk| chunk.len() == chunk.capacity()) {
            let capacity = self.len.clamp(MIN_CHUNK_ENTRIES, MAX_CHUNK_ENTRIES);
            self.chunks.push(Vec::with_capacity(capacity));
        }
        let entry = Entry { pos: position_bits(pos), lane_and_side, umi };
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

    fn seal(mut self) -> Self {
        for chunk in &mut self.chunks {
            chunk.sort_unstable();
            chunk.dedup();
        }
        self.text_codes = AHashMap::new();
        self
    }

    fn contains(&self, half: &HalfMapped, edits: u32) -> bool {
        let Some(&lane) = self.lanes.get(&half.lane) else { return false };
        let pos = position_bits(half.pos);
        let lane_and_side = lane | half.side;
        let tolerance = UmiTolerance { edits, truncate: half.truncate };

        if let Ok(umi) = half.umi
            && edits == 0
            && !half.truncate
        {
            let entry = Entry { pos, lane_and_side, umi };
            return self.chunks.iter().any(|chunk| chunk.binary_search(&entry).is_ok());
        }

        let first = Entry { pos, lane_and_side, umi: 0 };
        self.chunks.iter().any(|chunk| {
            chunk[chunk.partition_point(|entry| *entry < first)..]
                .iter()
                .take_while(|entry| entry.pos == pos && entry.lane_and_side == lane_and_side)
                .any(|entry| self.umi_matches(&half.umi, entry.umi, tolerance))
        })
    }

    fn umi_matches(
        &self,
        umi: &Result<UmiCode, Box<[u8]>>,
        known: UmiCode,
        tolerance: UmiTolerance,
    ) -> bool {
        let known_text =
            (known & TEXT_UMI != 0).then(|| self.texts[(known & !TEXT_UMI) as usize].as_ref());
        match (umi, known_text) {
            (Ok(packed), None) => tolerance.packed_match(*packed, known),
            (Ok(packed), Some(known)) => {
                let (bases, len) = unpack(*packed);
                tolerance.text_match(&bases[..len], known)
            }
            (Err(text), None) => {
                let (bases, len) = unpack(known);
                tolerance.text_match(text, &bases[..len])
            }
            (Err(text), Some(known)) => tolerance.text_match(text, known),
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
#[derive(Debug, Default)]
pub(crate) struct PairEnds {
    open: AHashMap<i32, ContigEnds>,
    /// Contigs whose half-mapped templates have been reached; no more ends are added to these.
    sealed: AHashMap<i32, ContigEnds>,
    /// The [`contig_rank`] the input has reached.
    rank: u16,
}

impl PairEnds {
    /// Adds a batch's ends.
    pub(crate) fn add(&mut self, buffer: &PairEndBuffer) {
        let mut ends = buffer.ends.as_slice();
        for (key, count) in &buffer.groups {
            let (group_ends, rest) = ends.split_at(*count);
            ends = rest;
            for (second, contig, pos, strand) in [
                (false, key.ref_id1, key.pos1, key.strand1),
                (true, key.ref_id2, key.pos2, key.strand2),
            ] {
                if self.sealed.contains_key(&contig) {
                    continue;
                }
                let known = self.open.entry(contig).or_default();
                let lane = (key.library_idx, key.cell_hash);
                let Some(lane) = known.lane_number(lane) else { continue };
                for end in group_ends.iter().filter(|end| end.second == second) {
                    let umi = if end.umi & TEXT_UMI == 0 {
                        Some(end.umi)
                    } else {
                        known.text_code(&buffer.texts[(end.umi & !TEXT_UMI) as usize])
                    };
                    if let Some(umi) = umi {
                        known.insert(pos, lane | side(strand != 0, end.is_r2), umi);
                    }
                }
            }
        }
    }

    /// Releases the contigs before the last of the batch, once its half-mapped templates have
    /// been looked up.
    pub(crate) fn release(&mut self, buffer: &PairEndBuffer) {
        let Some((key, _)) = buffer.groups.last() else { return };
        let rank = contig_rank(key.ref_id1);
        if rank != self.rank {
            self.rank = rank;
            self.open.retain(|&contig, _| contig_rank(contig) >= rank);
            self.sealed.retain(|&contig, _| contig_rank(contig) >= rank);
        }
    }

    /// Whether a fully mapped pair has the same end as `half`, with a UMI within `edits` of it.
    pub(crate) fn contains(&mut self, half: &HalfMapped, edits: u32) -> bool {
        let open = &mut self.open;
        self.sealed
            .entry(half.contig)
            .or_insert_with(|| open.remove(&half.contig).unwrap_or_default().seal())
            .contains(half, edits)
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

    fn end(second: bool, umi: &[u8]) -> PairEnd {
        PairEnd { second, is_r2: false, umi: pack(umi).expect("a packable UMI") }
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
    fn batch(key: GroupKey) -> PairEndBuffer {
        let ends = vec![end(false, b"ACGT"), end(true, b"TTTT")];
        PairEndBuffer { ends, groups: vec![(key, 2)], texts: Vec::new() }
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
        assert_eq!(tolerance.packed_match(packed_a, packed_b), expected);
        assert_eq!(tolerance.packed_match(packed_b, packed_a), expected);
        assert_eq!(tolerance.text_match(a, b), expected);
        assert_eq!(tolerance.text_match(b, a), expected);
    }

    #[test]
    fn ends_are_found_in_every_chunk() {
        let mut known = PairEnds::default();
        let count = 10 * i32::try_from(MIN_CHUNK_ENTRIES).unwrap();
        for pos in (0..count).rev() {
            known.add(&batch(pair_key(0, pos, 0, pos + 1_000_000)));
        }
        assert!(known.open[&0].chunks.len() > 2);
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
            known.add(&batch(key));
            known.release(&batch(key));
        }
        assert!(known.contains(&half(0, 100, false, b"ACGT"), 0));

        let next = batch(pair_key(1, 100, 1, 300));
        known.add(&next);
        assert!(known.contains(&half(0, 100, false, b"ACGT"), 0), "until the batch is released");
        known.release(&next);
        assert!(known.sealed.is_empty());
        assert_eq!(known.open.len(), 2, "contigs 1 and 2");
        assert!(known.contains(&half(2, 300, true, b"TTTT"), 0));
    }

    /// The sort key keeps 16 bits of a contig index, so contigs that share them are interleaved.
    #[test]
    fn contigs_that_sort_together_are_released_together() {
        let mut known = PairEnds::default();
        let halves = || [half(0, 100, false, b"ACGT"), half(1 << 16, 200, false, b"ACGT")];
        for half in halves() {
            let pairs = batch(pair_key(half.contig, half.pos, half.contig, 300));
            known.add(&pairs);
            known.release(&pairs);
        }
        for half in halves() {
            assert!(known.contains(&half, 0));
            let elsewhere = HalfMapped { contig: (1 << 16) - half.contig, ..half };
            assert!(!known.contains(&elsewhere, 0));
        }

        known.release(&batch(pair_key(1, 100, 1, 300)));
        assert!(known.open.is_empty() && known.sealed.is_empty());
    }

    /// Input in template-coordinate order has no pairs on a contig after its half-mapped templates.
    #[test]
    fn ends_after_a_contigs_half_mapped_templates_are_ignored() {
        let mut known = PairEnds::default();
        assert!(!known.contains(&half(0, 100, false, b"ACGT"), 0));
        known.add(&batch(pair_key(0, 100, 0, 300)));
        assert!(!known.contains(&half(0, 100, false, b"ACGT"), 0));
        assert!(known.open.is_empty());
    }
}
