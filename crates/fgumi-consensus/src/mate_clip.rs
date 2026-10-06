//! Read-through clipping for the simplex and duplex callers: the bases at a read's 3' end that
//! extend past its mate, with the mate CIGAR taken from the read's `MC` tag or, when that tag is
//! absent, from the mate record in the same group.

use ahash::AHashMap;
use anyhow::anyhow;
use fgumi_raw_bam::{RawRecordView, flags};

/// Computes the read-through clip — the bases a read extends past its mate — for the reads of
/// one group, taking the mate CIGAR from the read's `MC` tag or, when that tag is absent, from
/// the mate record in the group.
///
/// fgbio's `UmiConsensusCaller.updateMateCigars` (`UmiConsensusCaller.scala:328-343`, called at
/// `:356`) fills in `MC` from the in-group mate before any read is trimmed, so input without
/// `MC` is still clipped at the mate's end; without the backfill the clip would be 0 and
/// read-through adapter would reach the consensus. Deliberate differences from fgbio:
/// - a read that carries `MC` is measured against that tag, even when it disagrees with the
///   mate's actual CIGAR (fgbio overwrites both reads' `MC` when either read lacks it);
/// - a read lacking `MC` whose mate is not in the group gets `None` from [`Self::clip`], which
///   the callers turn into an error naming the read ([`missing_mate_error`]), and a name seen
///   more than twice uses the first mate found (fgbio fails both with a `MatchError`);
/// - a mate is needed by every primary read mapped to the same reference as its mate on the
///   opposite strand, RF (outie) pairs included, where fgbio only backfills, and so only fails
///   on, pairs it judges FR from `TLEN` — which misclassifies dovetailed FR pairs;
/// - the callers fail the run on a missing mate only for a read that can still reach a
///   consensus — simplex for a read in an emitted consensus, duplex for any read of a group that
///   passes its minimum-reads and strand-orientation checks — whereas fgbio backfills every
///   group before any filtering (`UmiConsensusCaller.scala:356`), so it also fails on groups it
///   then rejects;
/// - an in-group mate that is unmapped, or whose CIGAR is not a valid CIGAR for a mapped
///   record, gives a clip of 0, as an invalid `MC` tag does.
///
/// Only primary reads are backfilled or indexed as mates; a secondary or supplementary read
/// keeps its own `MC` tag (or no clip). fgbio never backfills them either: its
/// `ConsensusCallingIterator` drops secondary and supplementary reads before any caller sees
/// the group (`ConsensusCallingIterator.scala:57`).
///
/// The mate index is built on the first read that needs it, so input that carries `MC` — the
/// norm — pays only the one `MC` lookup per read that clipping always needed.
pub(crate) struct MateClipper<'a, I> {
    /// Every read of the group: where `mates` is built from.
    records: I,
    /// Primary paired reads keyed by `(read name, is first segment)`; built on first use.
    mates: Option<MatesByName<'a>>,
    /// Reused buffer for the decoded CIGAR of the mate being clipped against.
    mate_ops: Vec<u32>,
}

/// Primary paired reads of a group keyed by `(read name, is first segment)`.
type MatesByName<'a> = AHashMap<(&'a [u8], bool), &'a [u8]>;

impl<'a, I> MateClipper<'a, I>
where
    I: Iterator<Item = &'a [u8]> + Clone,
{
    /// Creates a clipper over `records`, all reads of one group. Nothing is indexed yet.
    pub(crate) fn new(records: I) -> Self {
        Self { records, mates: None, mate_ops: Vec::new() }
    }

    /// Returns the number of bases `rec` extends past its mate (0 unless an FR pair), or `None`
    /// when `rec` has no `MC` tag, needs its mate's CIGAR, and that mate is not in the group.
    pub(crate) fn clip(&mut self, rec: &[u8]) -> Option<usize> {
        if let Some(clip) = fgumi_raw_bam::num_bases_extending_past_mate_from_mc_raw(rec) {
            return Some(clip);
        }
        if !needs_mate(rec) {
            return Some(0);
        }
        let records = &self.records;
        let mates = self.mates.get_or_insert_with(|| index_mates(records.clone()));
        let mate = RawRecordView::new(mates.get(&mate_key(rec))?);
        if mate.flags() & flags::UNMAPPED != 0 {
            return Some(0);
        }
        self.mate_ops.clear();
        self.mate_ops.extend(mate.cigar_ops_iter());
        Some(fgumi_raw_bam::num_bases_extending_past_mate_cigar_raw(rec, &self.mate_ops))
    }
}

/// The error for a read whose [`MateClipper::clip`] is `None`: it needs its mate's CIGAR and has
/// neither an `MC` tag nor its mate in the group.
pub(crate) fn missing_mate_error(rec: &[u8]) -> anyhow::Error {
    let name = String::from_utf8_lossy(RawRecordView::new(rec).read_name());
    anyhow!(
        "Mate cigar (MC SAM tag) needed for read '{name}': the read has no MC tag and its primary \
         mate is not in the same group. Add MC tags (e.g. with fgumi zipper or samtools fixmate), \
         or keep both reads of each template in the same group."
    )
}

/// Indexes the primary paired reads of `records` by name and segment, keeping the first read
/// seen for each key.
fn index_mates<'a>(records: impl Iterator<Item = &'a [u8]>) -> MatesByName<'a> {
    let mut mates = MatesByName::new();
    for rec in records {
        if RawRecordView::new(rec).flags() & flags::PAIRED != 0 && is_primary(rec) {
            mates.entry(key(rec)).or_insert(rec);
        }
    }
    mates
}

fn is_primary(rec: &[u8]) -> bool {
    RawRecordView::new(rec).flags() & (flags::SECONDARY | flags::SUPPLEMENTARY) == 0
}

/// True when `rec`, a read lacking `MC`, could still extend past its mate: a primary read and its
/// mate both mapped to the same reference on opposite strands. FR orientation is not required
/// here because, without the mate CIGAR, it can only be inferred from `TLEN`, which
/// misclassifies dovetailed pairs (htsjdk/samtools#1771).
fn needs_mate(rec: &[u8]) -> bool {
    is_primary(rec) && fgumi_raw_bam::is_mapped_opposite_strand_pair_raw(rec)
}

fn key(rec: &[u8]) -> (&[u8], bool) {
    let view = RawRecordView::new(rec);
    (view.read_name(), view.flags() & flags::FIRST_SEGMENT != 0)
}

fn mate_key(rec: &[u8]) -> (&[u8], bool) {
    let (name, is_first) = key(rec);
    (name, !is_first)
}

/// Read-through pair fixture shared by the mate-clip, simplex and duplex tests.
#[cfg(test)]
pub(crate) mod read_through_fixture {
    use fgumi_dna::dna::reverse_complement;
    use fgumi_raw_bam::{RawRecord, SamBuilder, encode_op, flags};
    use fgumi_sam::SamTag;

    /// 80 bp genomic insert of the read-through pair.
    pub(crate) const READ_THROUGH_INSERT: &[u8] =
        b"ACCTCTCCATCTGACCCAAGATTGTGCTTGTTCAATTCTTCTTAACGTGATAACAGAATCAAACCTGCCAGGCGGTCGTC";
    /// 20 bp adapter sequenced past the end of the insert.
    pub(crate) const READ_THROUGH_ADAPTER: &[u8] = b"AGATCGGAAGAGCACACGTC";

    /// Which reads of a read-through pair carry an `MC` tag.
    #[derive(Clone, Copy, Debug)]
    pub(crate) enum McTags {
        /// Both reads carry the correct mate CIGAR.
        Both,
        /// Neither read carries an `MC` tag.
        Neither,
        /// Only the reverse-strand read carries `MC`; the forward-strand read lacks it.
        ForwardMissing,
        /// Both reads carry a (stale) `100M`, which is not the mate's real CIGAR.
        StaleBoth,
    }

    /// Builds a read-through FR pair `(R1, R2)`: a forward `80M20S` read and a reverse `20S80M`
    /// read, both 100 bp starting at 1-based position 1001, so each reads 20 bp of adapter past
    /// its mate. When `r1_forward` is set R1 is the forward read (a simplex or AB-strand
    /// template), otherwise R2 is (a BA-strand template).
    pub(crate) fn read_through_pair(
        name: &[u8],
        mi: &[u8],
        mc: McTags,
        r1_forward: bool,
    ) -> (RawRecord, RawRecord) {
        let mut fwd_seq = READ_THROUGH_INSERT.to_vec();
        fwd_seq.extend_from_slice(READ_THROUGH_ADAPTER);
        let mut rev_seq = reverse_complement(READ_THROUGH_ADAPTER);
        rev_seq.extend_from_slice(READ_THROUGH_INSERT);
        // (forward MC, reverse MC): each read's MC names its mate's CIGAR.
        let (fwd_mc, rev_mc): (Option<&[u8]>, Option<&[u8]>) = match mc {
            McTags::Both => (Some(b"20S80M"), Some(b"80M20S")),
            McTags::Neither => (None, None),
            McTags::ForwardMissing => (None, Some(b"80M20S")),
            McTags::StaleBoth => (Some(b"100M"), Some(b"100M")),
        };
        let build = |segment: u16, reverse: bool, mc: Option<&[u8]>| -> RawRecord {
            let (strand, tlen, seq, cigar) = if reverse {
                (flags::REVERSE, -80, &rev_seq, [encode_op(4, 20), encode_op(0, 80)])
            } else {
                (flags::MATE_REVERSE, 80, &fwd_seq, [encode_op(0, 80), encode_op(4, 20)])
            };
            let mut b = SamBuilder::new();
            b.read_name(name)
                .ref_id(0)
                .pos(1000)
                .flags(flags::PAIRED | flags::PROPER_PAIR | segment | strand)
                .mate_ref_id(0)
                .mate_pos(1000)
                .template_length(tlen)
                .sequence(seq)
                .qualities(&[40; 100])
                .cigar_ops(&cigar)
                .add_string_tag(SamTag::MI, mi);
            if let Some(mc) = mc {
                b.add_string_tag(SamTag::MC, mc);
            }
            b.build()
        };
        let (r1_mc, r2_mc) = if r1_forward { (fwd_mc, rev_mc) } else { (rev_mc, fwd_mc) };
        (
            build(flags::FIRST_SEGMENT, !r1_forward, r1_mc),
            build(flags::LAST_SEGMENT, r1_forward, r2_mc),
        )
    }
}

#[cfg(test)]
mod tests {
    use super::read_through_fixture::{McTags, read_through_pair};
    use super::*;
    use fgumi_raw_bam::{RawRecord, SamBuilder, encode_op};
    use fgumi_sam::SamTag;
    use rstest::rstest;

    /// Clips every read of `records` with one [`MateClipper`] over the group.
    fn clips(records: &[RawRecord]) -> Vec<Option<usize>> {
        let mut clipper = MateClipper::new(records.iter().map(AsRef::as_ref));
        records.iter().map(|r| clipper.clip(r.as_ref())).collect()
    }

    /// One read of a same-reference pair for the [`MateClipper`] tests, with `cigar` at 0-based
    /// `pos` and its mate at 0-based `mate_pos` on the opposite strand (the read is forward unless
    /// `reverse`), with `extra_flags` OR-ed in.
    #[expect(clippy::too_many_arguments, reason = "test fixture builder")]
    fn clipper_read(
        name: &[u8],
        first: bool,
        reverse: bool,
        pos: i32,
        mate_pos: i32,
        cigar: &[u32],
        mc: Option<&[u8]>,
        extra_flags: u16,
    ) -> RawRecord {
        let segment = if first { flags::FIRST_SEGMENT } else { flags::LAST_SEGMENT };
        let strand = if reverse { flags::REVERSE } else { flags::MATE_REVERSE };
        let len: usize = cigar
            .iter()
            .filter(|&&op| matches!(op & 0xf, 0 | 1 | 4))
            .map(|&op| (op >> 4) as usize)
            .sum();
        let mut b = SamBuilder::new();
        b.read_name(name)
            .ref_id(0)
            .pos(pos)
            .flags(flags::PAIRED | segment | strand | extra_flags)
            .mate_ref_id(0)
            .mate_pos(mate_pos)
            .sequence(&vec![b'A'; len])
            .qualities(&vec![40; len])
            .cigar_ops(cigar);
        if let Some(mc) = mc {
            b.add_string_tag(SamTag::MC, mc);
        }
        b.build()
    }

    /// `80M20S`, `20S80M` and `100M`, packed as `(len << 4) | op_code`.
    const FWD: [u32; 2] = [80 << 4, (20 << 4) | 4];
    const REV: [u32; 2] = [(20 << 4) | 4, 80 << 4];
    const FULL: [u32; 1] = [100 << 4];

    /// The groups the [`MateClipper`] clip tests run over.
    #[derive(Clone, Copy, Debug)]
    enum ClipperGroup {
        /// A read-through pair with the given `MC` tags.
        Pair(McTags),
        /// A forward R1 lacking `MC`, then two R2s of the same name: the read-through mate
        /// (`20S80M`) first and a `100M` one second.
        DuplicateNameReadThroughMateFirst,
        /// As [`ClipperGroup::DuplicateNameReadThroughMateFirst`], with the `100M` R2 first.
        DuplicateNameFullMateFirst,
        /// A read-through pair lacking `MC` whose reads carry no `SEQ`/`QUAL` (`*`).
        EmptySeqWithoutMc,
        /// An RF (outie) pair lacking `MC`: the forward read starts after its reverse mate ends.
        RfPairWithoutMc,
        /// A read-through pair lacking `MC`, plus a supplementary R1 of the same name that also
        /// lacks `MC`.
        SupplementaryWithoutMc,
    }

    impl ClipperGroup {
        fn records(self) -> Vec<RawRecord> {
            let r1 = || clipper_read(b"d", true, false, 1000, 1000, &FWD, None, 0);
            match self {
                Self::Pair(mc) => {
                    let (r1, r2) = read_through_pair(b"p", b"1", mc, true);
                    vec![r1, r2]
                }
                Self::DuplicateNameReadThroughMateFirst | Self::DuplicateNameFullMateFirst => {
                    let rt = clipper_read(b"d", false, true, 1000, 1000, &REV, Some(b"80M20S"), 0);
                    let fm = clipper_read(b"d", false, true, 1000, 1000, &FULL, Some(b"80M20S"), 0);
                    if matches!(self, Self::DuplicateNameReadThroughMateFirst) {
                        vec![r1(), rt, fm]
                    } else {
                        vec![r1(), fm, rt]
                    }
                }
                Self::EmptySeqWithoutMc => [(true, false, &FWD), (false, true, &REV)]
                    .into_iter()
                    .map(|(first, reverse, cigar)| {
                        let segment =
                            if first { flags::FIRST_SEGMENT } else { flags::LAST_SEGMENT };
                        let strand = if reverse { flags::REVERSE } else { flags::MATE_REVERSE };
                        let mut b = SamBuilder::new();
                        b.read_name(b"e")
                            .ref_id(0)
                            .pos(1000)
                            .flags(flags::PAIRED | segment | strand)
                            .mate_ref_id(0)
                            .mate_pos(1000)
                            .cigar_ops(cigar);
                        b.build()
                    })
                    .collect(),
                Self::RfPairWithoutMc => vec![
                    clipper_read(b"outie", true, false, 1200, 1000, &FULL, None, 0),
                    clipper_read(b"outie", false, true, 1000, 1200, &FULL, None, 0),
                ],
                Self::SupplementaryWithoutMc => {
                    let (r1, r2) = read_through_pair(b"p", b"1", McTags::Neither, true);
                    let supplementary = clipper_read(
                        b"p",
                        true,
                        false,
                        1000,
                        1000,
                        &FWD,
                        None,
                        flags::SUPPLEMENTARY,
                    );
                    vec![r1, r2, supplementary]
                }
            }
        }
    }

    /// [`MateClipper::clip`] measures a read with `MC` against the tag, and a read without `MC`
    /// against the CIGAR of its in-group mate. The read-through pair's reads each extend 20 bp
    /// past their mate; a stale `100M` on the forward read hides its overhang. When a name has
    /// more than two reads, the first mate in group order is used. The clip is measured on the
    /// CIGARs alone, so a read with no `SEQ` is clipped like any other. An RF pair never
    /// extends past its mate. A supplementary read is never backfilled, matching fgbio, which
    /// drops it before backfilling.
    #[rstest]
    #[case::both_mc(ClipperGroup::Pair(McTags::Both), &[20, 20])]
    #[case::neither_mc_backfilled(ClipperGroup::Pair(McTags::Neither), &[20, 20])]
    #[case::forward_mc_backfilled(ClipperGroup::Pair(McTags::ForwardMissing), &[20, 20])]
    #[case::stale_mc_used(ClipperGroup::Pair(McTags::StaleBoth), &[0, 20])]
    #[case::duplicate_name_first_mate_used(ClipperGroup::DuplicateNameReadThroughMateFirst, &[20, 20, 0])]
    #[case::duplicate_name_first_mate_used_swapped(ClipperGroup::DuplicateNameFullMateFirst, &[0, 0, 20])]
    #[case::empty_seq_measured_on_cigars(ClipperGroup::EmptySeqWithoutMc, &[20, 20])]
    #[case::rf_pair_not_clipped(ClipperGroup::RfPairWithoutMc, &[0, 0])]
    #[case::supplementary_not_backfilled(ClipperGroup::SupplementaryWithoutMc, &[20, 20, 0])]
    fn test_mate_clipper_clips_from_mc_or_in_group_mate(
        #[case] group: ClipperGroup,
        #[case] expected: &[usize],
    ) {
        let expected: Vec<Option<usize>> = expected.iter().copied().map(Some).collect();
        assert_eq!(clips(&group.records()), expected);
    }

    /// The in-group mate, at 1-based 1051, that the forward `80M20S` R1 at 1001 (lacking `MC`)
    /// is clipped against.
    #[derive(Clone, Copy, Debug)]
    enum MateRecord {
        /// A `30M` mate ending at 1080, where the R1's alignment ends: the usable control, past
        /// which the R1's 20 bp soft clip extends.
        Valid,
        /// Unmapped with no CIGAR, while the R1's own `MATE_UNMAPPED` flag is (inconsistently)
        /// clear.
        Unmapped,
        /// Unmapped but still carrying the `30M`, so only its flag rules it out.
        UnmappedWithCigar,
        /// Flagged mapped, with no CIGAR.
        MappedWithoutCigar,
        /// Fully soft-clipped (`100S`): no reference-consuming operation.
        SoftClipOnly,
        /// The `30M` mate truncated after its name, so its CIGAR is cut off.
        Truncated,
    }

    /// An in-group mate that is unmapped or whose CIGAR is not a valid CIGAR for a mapped record
    /// is not clipped against, exactly as an `MC` tag of `*` or `100S` is not. Before that check
    /// such a mate's empty reference span passed the FR gate and clipped the forward read at the
    /// mate's start, taking 30 aligned bases with its soft clip (a clip of 50).
    #[rstest]
    #[case::valid_control(MateRecord::Valid, Some(20))]
    #[case::unmapped(MateRecord::Unmapped, Some(0))]
    #[case::unmapped_with_cigar(MateRecord::UnmappedWithCigar, Some(0))]
    #[case::mapped_without_cigar(MateRecord::MappedWithoutCigar, Some(0))]
    #[case::soft_clip_only(MateRecord::SoftClipOnly, Some(0))]
    #[case::truncated(MateRecord::Truncated, Some(0))]
    fn test_mate_clipper_does_not_clip_against_an_unusable_mate(
        #[case] mate: MateRecord,
        #[case] expected: Option<usize>,
    ) {
        let r1 = clipper_read(b"d", true, false, 1000, 1050, &FWD, None, 0);
        let thirty_m = [encode_op(0, 30)];
        let (cigar, extra_flags): (&[u32], u16) = match mate {
            MateRecord::Valid | MateRecord::Truncated => (&thirty_m, 0),
            MateRecord::Unmapped => (&[], flags::UNMAPPED),
            MateRecord::UnmappedWithCigar => (&thirty_m, flags::UNMAPPED),
            MateRecord::MappedWithoutCigar => (&[], 0),
            MateRecord::SoftClipOnly => (&[encode_op(4, 100)], 0),
        };
        let r2 = clipper_read(b"d", false, true, 1050, 1000, cigar, None, extra_flags);
        // Keep the 32 fixed bytes and the name (`d\0`), but cut off the CIGAR and the rest.
        let r2_bytes =
            if matches!(mate, MateRecord::Truncated) { &r2.as_ref()[..34] } else { r2.as_ref() };
        let mut clipper = MateClipper::new([r1.as_ref(), r2_bytes].into_iter());
        assert_eq!(clipper.clip(r1.as_ref()), expected);
    }

    /// `clip` is `None` — the mate is required — for a primary read lacking `MC` whose mate is
    /// mapped to its reference on the opposite strand but is absent from the group. Every other
    /// case differs from `needs_mate` in exactly one guard, so it is that guard which makes the
    /// absent mate unnecessary; `clip` then returns 0 without consulting the group.
    #[rstest]
    #[case::needs_mate(flags::PAIRED | flags::MATE_REVERSE, 0, None, None)]
    #[case::read_carries_mc(flags::PAIRED | flags::MATE_REVERSE, 0, Some(&b"100M"[..]), Some(0))]
    #[case::mate_unmapped(flags::PAIRED | flags::MATE_REVERSE | flags::MATE_UNMAPPED, 0, None, Some(0))]
    #[case::read_unmapped(flags::PAIRED | flags::MATE_REVERSE | flags::UNMAPPED, 0, None, Some(0))]
    #[case::mate_on_other_reference(flags::PAIRED | flags::MATE_REVERSE, 1, None, Some(0))]
    #[case::mate_on_same_strand(flags::PAIRED, 0, None, Some(0))]
    #[case::unpaired(flags::MATE_REVERSE, 0, None, Some(0))]
    #[case::secondary(flags::PAIRED | flags::MATE_REVERSE | flags::SECONDARY, 0, None, Some(0))]
    #[case::supplementary(flags::PAIRED | flags::MATE_REVERSE | flags::SUPPLEMENTARY, 0, None, Some(0))]
    fn test_mate_clipper_requires_the_mate_only_when_its_cigar_is_needed(
        #[case] flag: u16,
        #[case] mate_ref_id: i32,
        #[case] mc: Option<&[u8]>,
        #[case] expected: Option<usize>,
    ) {
        let mut b = SamBuilder::new();
        b.read_name(b"lonely")
            .ref_id(0)
            .pos(1000)
            .flags(flag | flags::FIRST_SEGMENT)
            .mate_ref_id(mate_ref_id)
            .mate_pos(1100)
            .sequence(&[b'A'; 100])
            .qualities(&[40; 100])
            .cigar_ops(&FULL);
        if let Some(mc) = mc {
            b.add_string_tag(SamTag::MC, mc);
        }
        let read = b.build();
        assert_eq!(clips(std::slice::from_ref(&read)), vec![expected]);
    }

    /// The mate must be a primary read of the group: an FR or RF pair whose mate is absent, or
    /// present only as a secondary alignment, cannot be backfilled. An RF pair is included,
    /// which fgbio would not backfill or fail on: without the mate CIGAR, orientation can only
    /// be judged from `TLEN`, which misclassifies dovetailed FR pairs.
    #[rstest]
    #[case::fr_mate_absent(false, 1000, None)]
    #[case::rf_mate_absent(false, 1200, None)]
    #[case::mate_only_secondary(true, 1000, Some(flags::SECONDARY))]
    fn test_mate_clipper_needs_a_primary_mate_in_the_group(
        #[case] include_mate: bool,
        #[case] r1_pos: i32,
        #[case] mate_extra_flags: Option<u16>,
    ) {
        // R1 forward at `r1_pos`, its reverse mate at 1100: FR when R1 starts first, RF after.
        let r1 = clipper_read(b"needy", true, false, r1_pos, 1100, &FULL, None, 0);
        let mut records = vec![r1];
        if include_mate {
            let flags = mate_extra_flags.unwrap_or(0);
            records.push(clipper_read(b"needy", false, true, 1100, r1_pos, &FULL, None, flags));
        }
        assert_eq!(clips(&records)[0], None);
    }

    /// The missing-mate error names the read and says how to fix the input.
    #[test]
    fn test_missing_mate_error_names_the_read_and_the_remedy() {
        let (r1, _) = read_through_pair(b"orphan_read", b"1", McTags::Neither, true);
        assert_eq!(
            format!("{:#}", missing_mate_error(r1.as_ref())),
            "Mate cigar (MC SAM tag) needed for read 'orphan_read': the read has no MC tag and \
             its primary mate is not in the same group. Add MC tags (e.g. with fgumi zipper or \
             samtools fixmate), or keep both reads of each template in the same group."
        );
    }
}
