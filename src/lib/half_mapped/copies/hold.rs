//! Keeps the batches of a stream back until the half-mapped templates that go into them are known.
//!
//! A stream cannot be rewritten as a file can, so its batches go to a temporary file until the end
//! of the input, and then follow with the templates in them.

use std::fs::File;
use std::io::{self, BufReader, BufWriter, Read, Seek, Write};
use std::sync::Arc;
use std::sync::atomic::{AtomicBool, Ordering};

use fgumi_pipeline_core::Unpushed;
use fgumi_pipeline_core::held::HeldSlot;
use fgumi_pipeline_core::item::{HeapSize, Ordered};
use fgumi_pipeline_core::outputs::OrderedBytesSingle;
use fgumi_pipeline_core::queues::QueueSpec;
use fgumi_pipeline_core::reorder::BranchOrdering;
use fgumi_pipeline_core::step::{Step, StepCtx, StepKind, StepOutcome, StepProfile};

use fgumi_raw_bam::{RawRecord, write_framed_record};

use rayon::prelude::*;

use super::splice::rewrite;
use super::{BatchCopies, CopyRef, CopyStore, insert, insert_framed, split_records};
use crate::mi_group::MiGroup;
use crate::pipeline::steps::types::{BgzfBlock, DecompressedBlock};

/// Most batches taken from the input at a time.
const MAX_BATCHES_PER_RUN: usize = 8;

/// Batches put back together at a time, in parallel.
const RELEASE_WINDOW: usize = 256;

/// Zstandard level of the batches kept back.
const HOLD_LEVEL: i32 = 1;

/// Compresses records to keep them back, after their length.
fn pack(bytes: &[u8]) -> io::Result<Vec<u8>> {
    let len = u32::try_from(bytes.len()).map_err(io::Error::other)?;
    let mut packed = len.to_le_bytes().to_vec();
    packed.extend(zstd::bulk::compress(bytes, HOLD_LEVEL)?);
    Ok(packed)
}

fn unpack(packed: &[u8]) -> io::Result<Vec<u8>> {
    let (len, frame) = packed
        .split_first_chunk::<4>()
        .ok_or_else(|| io::Error::new(io::ErrorKind::UnexpectedEof, "a kept batch is truncated"))?;
    zstd::bulk::decompress(frame, u32::from_le_bytes(*len) as usize)
}

/// Whether batches have to be kept back, and the templates that go into them.
#[derive(Debug, Clone)]
pub(crate) struct CopyHold {
    /// Set unless the input is known to have no half-mapped templates.
    pub(crate) needed: Arc<AtomicBool>,
    pub(crate) store: Arc<CopyStore>,
}

impl CopyHold {
    pub(crate) fn is_needed(&self) -> bool {
        self.needed.load(Ordering::Relaxed)
    }
}

/// A batch of position groups that can be kept back as compressed records.
pub(crate) trait Holdable: HeapSize + Ordered + Send + Sized + 'static {
    /// The batch's records as BGZF blocks.
    fn into_blocks(self) -> io::Result<Vec<u8>>;

    /// The batch `serial` with `copies` put into it, as [`insert`] does.
    fn from_blocks(
        serial: u64,
        blocks: Vec<u8>,
        copies: BatchCopies,
        restore: &Restore,
    ) -> io::Result<Self>;
}

/// What batches need to be put back together.
pub(crate) struct Restore {
    store: Arc<CopyStore>,
    /// Level of the blocks that are written out.
    compression_level: u32,
    rules: MiGroupRules,
}

/// Serial step that passes batches on, or keeps them back until the end of the input.
pub(crate) struct HoldForCopies<T: Holdable> {
    hold: CopyHold,
    state: State<T>,
    held: HeldSlot<Unpushed<T>>,
    restore: Restore,
    limit_bytes: u64,
}

enum State<T> {
    Passing,
    Holding { file: BufWriter<File>, batches: Vec<(u64, u64)> },
    Releasing(Release<T>),
}

struct Release<T> {
    file: BufReader<File>,
    batches: std::vec::IntoIter<(u64, u64)>,
    copies: std::iter::Peekable<std::vec::IntoIter<(usize, Vec<CopyRef>)>>,
    ready: std::vec::IntoIter<T>,
}

impl<T: Holdable> HoldForCopies<T> {
    pub(crate) fn new(hold: CopyHold, compression_level: u32, limit_bytes: u64) -> Self {
        let store = Arc::clone(&hold.store);
        let restore = Restore { store, compression_level, rules: MiGroupRules::default() };
        Self { hold, state: State::Passing, held: HeldSlot::new(), restore, limit_bytes }
    }

    pub(crate) fn with_mi_group_rules(mut self, rules: MiGroupRules) -> Self {
        self.restore.rules = rules;
        self
    }

    fn hold_back(&mut self, batch: T) -> io::Result<()> {
        if matches!(self.state, State::Passing) {
            let file = BufWriter::with_capacity(1 << 20, tempfile::tempfile()?);
            self.state = State::Holding { file, batches: Vec::new() };
        }
        let State::Holding { file, batches } = &mut self.state else {
            return Err(io::Error::other("a batch came after the end of the input"));
        };
        let serial = batch.ordinal();
        let blocks = batch.into_blocks()?;
        file.write_all(&blocks)?;
        batches.push((serial, blocks.len() as u64));
        Ok(())
    }

    /// The next batch that was kept back, with its templates.
    fn release(&mut self) -> io::Result<Option<T>> {
        if let State::Holding { .. } = self.state
            && let State::Holding { file, batches } =
                std::mem::replace(&mut self.state, State::Passing)
        {
            let mut file = file.into_inner().map_err(io::IntoInnerError::into_error)?;
            file.rewind()?;
            self.state = State::Releasing(Release {
                file: BufReader::with_capacity(1 << 20, file),
                batches: batches.into_iter(),
                copies: self.hold.store.by_batch().into_iter().peekable(),
                ready: Vec::new().into_iter(),
            });
        }
        let State::Releasing(Release { file, batches, copies, ready }) = &mut self.state else {
            return Ok(None);
        };
        if let Some(batch) = ready.next() {
            return Ok(Some(batch));
        }
        let mut window = Vec::with_capacity(RELEASE_WINDOW);
        for (serial, len) in batches.take(RELEASE_WINDOW) {
            let mut blocks = vec![0; usize::try_from(len).map_err(io::Error::other)?];
            file.read_exact(&mut blocks)?;
            let copies = match copies.next_if(|(batch, _)| *batch as u64 == serial) {
                Some((_, refs)) => self.hold.store.read(&refs)?,
                None => Vec::new(),
            };
            window.push((serial, blocks, copies));
        }
        if window.is_empty() && copies.next().is_some() {
            return Err(io::Error::other("a template's molecule was not found"));
        }
        let restore = &self.restore;
        let restored: Vec<T> = window
            .into_par_iter()
            .map(|(serial, blocks, copies)| T::from_blocks(serial, blocks, copies, restore))
            .collect::<io::Result<_>>()?;
        *ready = restored.into_iter();
        Ok(ready.next())
    }
}

impl<T: Holdable> Step for HoldForCopies<T> {
    type Input = T;
    type Outputs = OrderedBytesSingle<T>;

    fn profile(&self) -> StepProfile {
        StepProfile {
            name: "HoldForCopies",
            kind: StepKind::Serial,
            sticky: false,
            output_queues: vec![QueueSpec::ByteBounded { limit_bytes: self.limit_bytes }],
            branch_ordering: vec![BranchOrdering::ByItemOrdinal],
        }
    }

    fn try_run(&mut self, ctx: &mut StepCtx<'_, Self>) -> io::Result<StepOutcome> {
        if let Some(unpushed) = self.held.take()
            && let Err(again) = ctx.outputs.retry(unpushed)
        {
            self.held.put(again);
            return Ok(StepOutcome::Contention);
        }

        let mut progress = false;
        for _ in 0..MAX_BATCHES_PER_RUN {
            let Some(batch) = ctx.input.pop() else { break };
            progress = true;
            if self.hold.is_needed() {
                self.hold_back(batch)?;
            } else if let Err(unpushed) = ctx.outputs.push(batch) {
                self.held.put(unpushed);
                break;
            }
        }
        if progress {
            return Ok(StepOutcome::Progress);
        }
        if !ctx.input.is_drained() {
            return Ok(StepOutcome::NoProgress);
        }
        for _ in 0..MAX_BATCHES_PER_RUN {
            let Some(batch) = self.release()? else {
                return Ok(if progress { StepOutcome::Progress } else { StepOutcome::Finished });
            };
            progress = true;
            if let Err(unpushed) = ctx.outputs.push(batch) {
                self.held.put(unpushed);
                break;
            }
        }
        Ok(StepOutcome::Progress)
    }
}

impl Holdable for BgzfBlock {
    fn into_blocks(self) -> io::Result<Vec<u8>> {
        Ok(self.bytes.to_vec())
    }

    fn from_blocks(
        serial: u64,
        blocks: Vec<u8>,
        copies: BatchCopies,
        restore: &Restore,
    ) -> io::Result<Self> {
        let bytes = if copies.is_empty() {
            blocks
        } else {
            let (bytes, changes) = rewrite(&blocks, copies, restore.compression_level)?;
            restore.store.add_size_changes(&changes);
            bytes
        };
        Ok(Self { batch_serial: serial, bytes: bytes.into(), uncompressed_size: 0, index: None })
    }
}

impl Holdable for DecompressedBlock {
    fn into_blocks(self) -> io::Result<Vec<u8>> {
        pack(&self.bytes)
    }

    fn from_blocks(
        serial: u64,
        blocks: Vec<u8>,
        copies: BatchCopies,
        restore: &Restore,
    ) -> io::Result<Self> {
        let mut bytes = unpack(&blocks)?;
        if !copies.is_empty() {
            let (merged, changes) = insert_framed(&bytes, copies)?;
            restore.store.add_size_changes(&changes);
            bytes = merged;
        }
        Ok(Self { batch_serial: serial, bytes })
    }
}

/// The records of a batch of MI groups, kept back as BGZF blocks.
#[derive(Debug)]
pub(crate) struct HeldMiGroups {
    blocks: Vec<u8>,
    copies: BatchCopies,
    store: Option<Arc<CopyStore>>,
    rules: MiGroupRules,
}

/// How records are put into MI groups.
#[derive(Debug, Clone, Copy, Default)]
pub(crate) struct MiGroupRules {
    pub(crate) strip_strand: bool,
    pub(crate) allow_unmapped: bool,
}

impl HeldMiGroups {
    pub(crate) fn pack(groups: &[MiGroup], rules: MiGroupRules) -> io::Result<Self> {
        let mut bytes = Vec::new();
        for record in groups.iter().flat_map(|group| &group.records) {
            write_framed_record(&mut bytes, record)?;
        }
        let blocks = pack(&bytes)?;
        Ok(Self { blocks, copies: Vec::new(), store: None, rules })
    }

    pub(crate) fn unpack(self) -> io::Result<Vec<MiGroup>> {
        let Self { blocks, mut copies, store, rules } = self;
        let mut records = split_records(&unpack(&blocks)?)?;
        if let Some(store) = store.filter(|_| !copies.is_empty()) {
            let keep = |record: &_| {
                crate::commands::common::consensus_pregroup_keep_raw(record, rules.allow_unmapped)
            };
            copies.retain(|copy| keep(copy));
            let (merged, changes) = insert(records, copies);
            store.add_size_changes(&changes);
            records = merged;
        }
        let mut groups: Vec<MiGroup> = Vec::new();
        for record in records {
            let aux = fgumi_raw_bam::aux_data_slice(&record);
            let mi = fgumi_raw_bam::find_string_tag(aux, crate::sam::SamTag::MI)
                .ok_or_else(|| io::Error::other("a grouped record has no MI"))?;
            let mi = String::from_utf8_lossy(mi);
            let mi = if rules.strip_strand { fgumi_umi::extract_mi_base(&mi) } else { &mi };
            match groups.last_mut() {
                Some(group) if group.mi == mi => group.records.push(record),
                _ => groups.push(MiGroup::new(mi.to_owned(), vec![record])),
            }
        }
        Ok(groups)
    }

    pub(crate) fn heap_size(&self) -> usize {
        self.blocks.capacity() + self.copies.iter().map(RawRecord::capacity).sum::<usize>()
    }
}

impl Holdable for crate::pipeline::steps::group::mi::BatchedMiGroups {
    fn into_blocks(self) -> io::Result<Vec<u8>> {
        self.held.map(|held| held.blocks).ok_or_else(|| io::Error::other("a batch is not packed"))
    }

    fn from_blocks(
        serial: u64,
        blocks: Vec<u8>,
        copies: BatchCopies,
        restore: &Restore,
    ) -> io::Result<Self> {
        let store = Some(Arc::clone(&restore.store));
        let rules = restore.rules;
        Ok(Self::held(serial, HeldMiGroups { blocks, copies, store, rules }))
    }
}

#[cfg(test)]
mod tests {
    use super::super::tests::template;
    use super::*;
    use rstest::rstest;

    /// Groups that were kept back are those that were not, with the copies in them.
    #[rstest]
    #[case::simplex(false, false, &[("3/A", 3), ("3/B", 2), ("4/A", 2)])]
    #[case::duplex(true, false, &[("3", 5), ("4", 2)])]
    #[case::allow_unmapped(false, true, &[("3/A", 4), ("3/B", 2), ("4/A", 2)])]
    fn held_groups_are_grouped_again(
        #[case] strip_strand: bool,
        #[case] allow_unmapped: bool,
        #[case] expected: &[(&str, usize)],
    ) {
        let groups: Vec<_> = [("a", "3/A"), ("b", "3/B"), ("c", "4/A")]
            .into_iter()
            .map(|(name, mi)| MiGroup::new(mi.to_string(), template(name, Some(mi), false)))
            .collect();
        let rules = MiGroupRules { strip_strand, allow_unmapped };
        let mut held = HeldMiGroups::pack(&groups, rules).unwrap();
        let store = Arc::new(CopyStore::default());
        held.copies = template("d", Some("3/A"), true);
        held.store = Some(Arc::clone(&store));
        let groups = held.unpack().unwrap();
        let sizes: Vec<_> = groups.iter().map(|g| (g.mi.as_str(), g.records.len())).collect();
        assert_eq!(sizes, expected);
        assert_eq!(store.take_size_changes().len(), 1);
    }
}
