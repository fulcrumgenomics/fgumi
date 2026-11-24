//! Places half-mapped pairs beside the fully mapped pair they are copies of.
//!
//! A pair with one unmapped read has one known end, so the template-coordinate sort would put it
//! after every fully mapped pair on its contig. [`HalfMappedPlacer`] keeps the records of such
//! templates out of the runs and, as keys are extracted, writes each fully mapped pair's
//! coordinate and the hashes of its two UMIs to a side list (24 bytes a pair, in a temporary file
//! once it is large). When the input ends it reads the side list once and looks for each
//! half-mapped template's end and UMI. A template found with exactly one second end gets that
//! pair's coordinate in its `tc` tag, which both the sort and the grouper key it by. The records
//! then join the merge as one more sorted run.

use std::borrow::Cow;
use std::fs::File;
use std::io::{self, Read, Seek, Write};
use std::path::PathBuf;
use std::sync::{Arc, OnceLock};

use ahash::AHashMap;
use anyhow::Result;
use bytemuck::{Pod, Zeroable};
use fgumi_raw_bam::{RawRecordView, SamTag, TemplateCoordinateTag, flags};
use parking_lot::Mutex;

use crate::arena_pool::PooledSegmentedBuf;
use crate::inline::TemplateKey;
use crate::segmented_buf::SegmentedBuf;
use crate::template_arena::{TemplateArenaRefs, TemplateMemChunk, seal_template_refs};
use crate::template_key::{TemplateKeyContext, dropped_lane_error, violation_read_name};

/// Bytes of the side list kept in memory before it goes to a temporary file.
const SPOOL_BUFFER_BYTES: usize = 16 << 20;

/// The most reference sequences a key tells apart, as it keeps 16 bits of their index.
pub const MAX_REFERENCES: usize = 0xFFFF;

const PAIR_BYTES: usize = std::mem::size_of::<PairEnds>();

/// A fully mapped pair's coordinate, as the sort key packs it, and the hash of each end's UMI.
#[derive(Debug, Clone, Copy, Pod, Zeroable)]
#[repr(C)]
struct PairEnds {
    primary: u64,
    secondary: u64,
    umis: [u32; 2],
}

impl PairEnds {
    fn ends(&self) -> [End; 2] {
        let TemplateKeyEnds { lower, upper } = TemplateKeyEnds::new(self.primary, self.secondary);
        [lower, upper]
    }
}

/// Contig, 0-based unclipped 5' position, and whether on the reverse strand.
type End = (i32, i32, bool);

struct TemplateKeyEnds {
    lower: End,
    upper: End,
}

impl TemplateKeyEnds {
    /// Unpacks [`TemplateKey::new`]'s first two lanes.
    #[allow(clippy::cast_possible_truncation, clippy::cast_possible_wrap)]
    fn new(primary: u64, secondary: u64) -> Self {
        let tid = |shift: u32| ((primary >> shift) & 0xFFFF) as i32;
        let pos = |packed: u64| (packed as u32 ^ 0x8000_0000) as i32;
        Self {
            lower: (tid(48), pos(primary), secondary & 2 == 0),
            upper: (tid(32), pos(secondary >> 32), secondary & 1 == 0),
        }
    }
}

/// Hashes ends and UMIs; the seeds are fixed so that runs are reproducible.
fn hasher() -> ahash::RandomState {
    ahash::RandomState::with_seeds(
        0x9e37_79b9_7f4a_7c15,
        0xbf58_476d_1ce4_e5b9,
        0x94d0_49bb_1331_11eb,
        0x2545_f491_4f6c_dd1d,
    )
}

/// What one batch of records adds to a [`HalfMappedPlacer`].
#[derive(Debug, Default)]
pub struct PlacerBatch {
    pairs: Vec<PairEnds>,
    held: Vec<Vec<u8>>,
}

impl PlacerBatch {
    pub(crate) fn len(&self) -> usize {
        self.pairs.len() + self.held.len()
    }
}

/// The side list, in memory until it is large.
#[derive(Debug)]
struct Spool {
    buffer: Vec<u8>,
    file: Option<File>,
    dir: PathBuf,
}

impl Spool {
    fn extend(&mut self, bytes: &[u8]) -> io::Result<()> {
        self.buffer.extend_from_slice(bytes);
        if self.buffer.len() >= SPOOL_BUFFER_BYTES {
            if self.file.is_none() {
                self.file = Some(tempfile::tempfile_in(&self.dir)?);
            }
            if let Some(file) = &mut self.file {
                file.write_all(&self.buffer)?;
                self.buffer.clear();
            }
        }
        Ok(())
    }

    /// Gives every pair to `f`, and empties the list.
    fn drain(&mut self, mut f: impl FnMut(&PairEnds)) -> io::Result<()> {
        let mut each = |bytes: &[u8]| {
            for pair in bytes.chunks_exact(PAIR_BYTES) {
                f(&bytemuck::pod_read_unaligned(pair));
            }
        };
        if let Some(mut file) = self.file.take() {
            file.rewind()?;
            let mut block = vec![0; 1 << 20];
            block.truncate(block.len() / PAIR_BYTES * PAIR_BYTES);
            let mut filled = 0;
            loop {
                let n = file.read(&mut block[filled..])?;
                filled += n;
                if n == 0 || filled == block.len() {
                    each(&block[..filled]);
                    if n == 0 {
                        break;
                    }
                    filled = 0;
                }
            }
        }
        each(&self.buffer);
        self.buffer = Vec::new();
        Ok(())
    }
}

/// The records of half-mapped templates, each with its place in the input.
type Held = Vec<(u64, Vec<Vec<u8>>)>;

/// See the module documentation.
#[derive(Debug)]
pub struct HalfMappedPlacer {
    hasher: ahash::RandomState,
    pairs: Mutex<Spool>,
    held: Mutex<Held>,
    context: OnceLock<Arc<TemplateKeyContext>>,
    error: Mutex<Option<io::Error>>,
}

impl HalfMappedPlacer {
    /// A placer whose side list goes to a temporary file in `dir` once it is large.
    #[must_use]
    pub fn new(dir: PathBuf) -> Arc<Self> {
        Arc::new(Self {
            hasher: hasher(),
            pairs: Mutex::new(Spool { buffer: Vec::new(), file: None, dir }),
            held: Mutex::new(Vec::new()),
            context: OnceLock::new(),
            error: Mutex::new(None),
        })
    }

    /// Keeps the context that keys the placed records; the first call wins.
    pub(crate) fn set_context(&self, context: impl FnOnce() -> Arc<TemplateKeyContext>) {
        self.context.get_or_init(context);
    }

    /// Reads the record `body`, whose full key is `key`, into `batch`. Returns whether the record
    /// is one of a half-mapped template, which must then be left out of the run.
    pub(crate) fn observe(&self, body: &[u8], key: &TemplateKey, batch: &mut PlacerBatch) -> bool {
        let flag = RawRecordView::new(body).flags();
        if flag & flags::PAIRED == 0 {
            return false;
        }
        let unmapped = flag & flags::UNMAPPED != 0;
        let mate_unmapped = flag & flags::MATE_UNMAPPED != 0;
        let is_primary = flag & (flags::SECONDARY | flags::SUPPLEMENTARY) == 0;
        if unmapped != mate_unmapped || (!is_primary && has_one_primary(body, flag)) {
            batch.held.push(body.to_vec());
            return true;
        }
        let is_pair_r1 = is_primary && flag & flags::FIRST_SEGMENT != 0;
        if is_pair_r1
            && !unmapped
            && let (Some(r1), Some(r2)) = (self.umi(body, key, true), self.umi(body, key, false))
        {
            let r1_is_upper = key.name_hash_upper & 1 == 1;
            let umis = if r1_is_upper { [r2, r1] } else { [r1, r2] };
            batch.pairs.push(PairEnds { primary: key.primary, secondary: key.secondary, umis });
        }
        false
    }

    /// Adds `batch`, the `order`th of the input. [`Self::take_chunk`] reports an error in this.
    pub(crate) fn commit(&self, batch: PlacerBatch, order: u64) {
        if !batch.pairs.is_empty()
            && let Err(error) = self.pairs.lock().extend(bytemuck::cast_slice(&batch.pairs))
        {
            *self.error.lock() = Some(error);
        }
        if !batch.held.is_empty() {
            self.held.lock().push((order, batch.held));
        }
    }

    /// Hashes the library, cell and the UMI of a read: the half of a two-segment UMI that the read
    /// carries, or the whole of any other. `None` if that is empty.
    fn umi(&self, body: &[u8], key: &TemplateKey, is_r1: bool) -> Option<u32> {
        let umi = fgumi_raw_bam::find_string_tag(fgumi_raw_bam::aux_data_slice(body), SamTag::RX);
        let umi = match umi {
            None => None,
            Some(umi) => Some(match umi.iter().position(|&b| b == b'-') {
                Some(dash) if !umi[dash + 1..].contains(&b'-') => {
                    if is_r1 {
                        &umi[..dash]
                    } else {
                        &umi[dash + 1..]
                    }
                }
                _ => umi,
            }),
        };
        if umi.is_some_and(<[u8]>::is_empty) {
            return None;
        }
        let library = key.tertiary >> 48;
        let umi = umi.map(|umi| {
            if umi.iter().any(u8::is_ascii_lowercase) {
                Cow::Owned(umi.to_ascii_uppercase())
            } else {
                Cow::Borrowed(umi)
            }
        });
        #[allow(clippy::cast_possible_truncation)]
        Some(self.hasher.hash_one((library, key.cb_hash, umi)) as u32)
    }

    fn slot(&self, end: End, umi: u32) -> u64 {
        self.hasher.hash_one((end, umi))
    }

    /// The half-mapped templates' records as a sorted chunk, each template placed if it can be.
    /// `None` if there are none.
    ///
    /// # Errors
    ///
    /// If the side list could not be written or read, or a record has a value in a key lane that `--key-types`
    /// dropped.
    pub fn take_chunk(&self) -> Result<Option<TemplateMemChunk>> {
        if let Some(error) = self.error.lock().take() {
            return Err(error.into());
        }
        let mut held = std::mem::take(&mut *self.held.lock());
        let Some(context) = self.context.get().filter(|_| !held.is_empty()) else {
            return Ok(None);
        };
        held.sort_by_key(|(order, _)| *order);
        let mut records: Vec<Vec<u8>> = held.into_iter().flat_map(|(_, r)| r).collect();

        let mut templates: AHashMap<Vec<u8>, Vec<usize>> = AHashMap::new();
        for (i, record) in records.iter().enumerate() {
            let name = fgumi_raw_bam::read_name(record);
            match templates.get_mut(name) {
                Some(indices) => indices.push(i),
                None => {
                    templates.insert(name.to_vec(), vec![i]);
                }
            }
        }

        let mut found: AHashMap<u64, Candidate> = AHashMap::new();
        for indices in templates.values() {
            if let Some((end, umi)) = self.mapped_end(context, &records, indices) {
                found.entry(self.slot(end, umi)).or_default();
            }
        }
        self.pairs.lock().drain(|pair| {
            let [lower, upper] = pair.ends();
            for (end, umi) in [(lower, pair.umis[0]), (upper, pair.umis[1])] {
                if let Some(candidate) = found.get_mut(&self.slot(end, umi)) {
                    candidate.add((pair.primary, pair.secondary));
                }
            }
        })?;

        let mut placed = 0u64;
        for indices in templates.values() {
            let Some((end, umi)) = self.mapped_end(context, &records, indices) else { continue };
            let tc = match found[&self.slot(end, umi)] {
                Candidate::One(primary, secondary) => {
                    let TemplateKeyEnds { lower, upper } = TemplateKeyEnds::new(primary, secondary);
                    let one_based = |(tid, pos, neg): End| (tid, pos + 1, neg);
                    let tc =
                        TemplateCoordinateTag::from_primaries(one_based(lower), one_based(upper));
                    (lower == end || upper == end).then_some(tc)
                }
                Candidate::None | Candidate::Many => None,
            };
            let one_primary = TemplateCoordinateTag::one_primary((end.0, end.1 + 1, end.2));
            placed += u64::from(place(&mut records, indices, tc, one_primary));
        }
        log::info!(
            "Placed {placed} of {} half-mapped templates beside a fully mapped pair",
            templates.len()
        );

        let mut arena = SegmentedBuf::new();
        let mut refs = TemplateArenaRefs::for_variant(context.variant);
        for record in records {
            let (key, violation) = context.extract_one(&record);
            if let Some(violation) = violation {
                return Err(dropped_lane_error(&violation_read_name(&record), violation));
            }
            let offset = arena.extend_from_slice(&record) as u64;
            refs.push_full(&key, offset, u32::try_from(record.len())?);
        }
        let pool = rayon::ThreadPoolBuilder::new().num_threads(1).build()?;
        let arena = Arc::new(PooledSegmentedBuf::unpooled(arena));
        Ok(Some(seal_template_refs(&mut refs, arena, 1, &pool)))
    }

    /// The mapped end of a half-mapped template and the hash of its UMI.
    fn mapped_end(
        &self,
        context: &TemplateKeyContext,
        records: &[Vec<u8>],
        indices: &[usize],
    ) -> Option<(End, u32)> {
        let mut mapped = indices.iter().map(|&i| &records[i]).filter(|record| {
            let flag = RawRecordView::new(record).flags();
            flag & (flags::UNMAPPED | flags::SECONDARY | flags::SUPPLEMENTARY) == 0
        });
        let record = mapped.next().filter(|_| mapped.next().is_none())?;
        let view = RawRecordView::new(record);
        let flag = view.flags();
        let pos = fgumi_raw_bam::unclipped_5prime_from_raw_bam(record);
        if view.ref_id() < 0 || pos == i32::MAX {
            return None;
        }
        let (key, _) = context.extract_one(record);
        let umi = self.umi(record, &key, flag & flags::FIRST_SEGMENT != 0)?;
        Some(((view.ref_id(), pos - 1, flag & flags::REVERSE != 0), umi))
    }
}

/// Whether the `tc` tag of a secondary or supplementary record says that its template has one
/// mapped primary read.
fn has_one_primary(body: &[u8], flag: u16) -> bool {
    TemplateCoordinateTag::from_record_aux(fgumi_raw_bam::aux_data_slice(body), flag)
        .is_some_and(|tc| tc.has_one_primary())
}

/// The fully mapped pairs seen with one end and UMI, by where their other end is.
#[derive(Debug, Clone, Copy, Default)]
enum Candidate {
    #[default]
    None,
    One(u64, u64),
    Many,
}

impl Candidate {
    fn add(&mut self, pair: (u64, u64)) {
        *self = match *self {
            Self::None => Self::One(pair.0, pair.1),
            Self::One(primary, secondary) if (primary, secondary) == pair => *self,
            Self::One(..) | Self::Many => Self::Many,
        };
    }
}

/// Gives the records of a template the coordinate `tc`, and returns whether that worked. Without
/// one, undoes an earlier sort's placement, which left `one_primary` on no record but secondary
/// and supplementary ones.
fn place(
    records: &mut [Vec<u8>],
    indices: &[usize],
    tc: Option<TemplateCoordinateTag>,
    one_primary: TemplateCoordinateTag,
) -> bool {
    let is_primary = |record: &[u8]| {
        RawRecordView::new(record).flags() & (flags::SECONDARY | flags::SUPPLEMENTARY) == 0
    };
    let set = |record: &mut Vec<u8>, tc: &TemplateCoordinateTag| {
        fgumi_raw_bam::remove_tag(record, SamTag::TC);
        fgumi_raw_bam::append_i32_array_tag(record, SamTag::TC, &tc.to_values());
    };
    let was_placed =
        indices.iter().any(|&i| TemplateCoordinateTag::from_half_mapped(&records[i]).is_some());

    if let Some(tc) = tc {
        let updated: Vec<Vec<u8>> = indices
            .iter()
            .map(|&i| {
                let mut record = records[i].clone();
                set(&mut record, &tc);
                record
            })
            .collect();
        let keyed_alike = updated.iter().all(|record| {
            !is_primary(record) || TemplateCoordinateTag::from_half_mapped(record) == Some(tc)
        });
        if keyed_alike {
            for (&i, record) in indices.iter().zip(updated) {
                records[i] = record;
            }
            return true;
        }
    }
    if was_placed {
        for &i in indices {
            if is_primary(&records[i]) {
                fgumi_raw_bam::remove_tag(&mut records[i], SamTag::TC);
            } else {
                set(&mut records[i], &one_primary);
            }
        }
    }
    false
}

#[cfg(test)]
mod tests {
    use super::*;
    use rstest::rstest;

    #[rstest]
    #[case::one_contig((0, 100, false), (0, 300, true))]
    #[case::two_contigs((3, 100, true), (7, 50, false))]
    #[case::before_the_contig((1, -5, false), (1, 0, true))]
    #[case::last_contig((0, 0, false), (0xFFFE, i32::MAX - 1, true))]
    fn ends_are_unpacked_from_the_key(#[case] lower: End, #[case] upper: End) {
        let (tid1, pos1, neg1) = lower;
        let (tid2, pos2, neg2) = upper;
        let key = TemplateKey::new(tid1, pos1, neg1, tid2, pos2, neg2, 0, 0, (0, false), 0, false);
        let pair = PairEnds { primary: key.primary, secondary: key.secondary, umis: [0, 0] };
        assert_eq!(pair.ends(), [lower, upper]);
    }

    #[test]
    fn the_side_list_is_read_back_whole() {
        let dir = tempfile::tempdir().unwrap();
        let mut spool = Spool { buffer: Vec::new(), file: None, dir: dir.path().into() };
        let count = 3 * SPOOL_BUFFER_BYTES / PAIR_BYTES / 2;
        let pair = |i: usize| PairEnds { primary: i as u64, secondary: 7, umis: [1, 2] };
        for i in 0..count {
            spool.extend(bytemuck::bytes_of(&pair(i))).unwrap();
        }
        assert!(spool.file.is_some(), "the list is large enough for a file");
        let mut next = 0;
        spool
            .drain(|read| {
                assert_eq!(bytemuck::bytes_of(read), bytemuck::bytes_of(&pair(next)));
                next += 1;
            })
            .unwrap();
        assert_eq!(next, count);
    }

    #[test]
    fn a_template_has_a_place_only_if_every_pair_agrees() {
        let mut candidate = Candidate::default();
        candidate.add((1, 2));
        candidate.add((1, 2));
        assert!(matches!(candidate, Candidate::One(1, 2)));
        candidate.add((1, 3));
        candidate.add((1, 2));
        assert!(matches!(candidate, Candidate::Many));
    }
}
