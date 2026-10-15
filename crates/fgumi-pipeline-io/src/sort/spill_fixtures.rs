//! Spill-run fixtures shared by the sort's unit tests (`spill_write`,
//! `spill_read_planner`, `spill_decompress`, `consolidate_runs`): keyed spill
//! runs, run files written and read back at the frame level, loaded merge
//! slots, the stable-merge oracle, and `KChain`, the consolidation pipeline
//! driven step by step.
#![allow(clippy::cast_possible_truncation, clippy::cast_possible_wrap, clippy::type_complexity)]

use std::io::{self, Write as _};
use std::ops::Range;
use std::path::{Path, PathBuf};
use std::sync::Arc;
use std::sync::atomic::AtomicBool;

use fgumi_bam_io::pread::SliceBufferPool;
use fgumi_pipeline_core::handles::BranchInputHandle;
use fgumi_pipeline_core::step::{InputHandle, Step, StepOutcome};
use fgumi_pipeline_core::testing::StepProbe;
use fgumi_pipeline_core::{PhasedCap, StepCap};
use fgumi_sort::{
    InMemoryChunk, MergeDemand, RawCoordinateKey, RawSortKey, SortMergeSlot, SpillBlockCompressor,
    SpillCodec, SpillKeyKind, TemplateKey, TemplateMemChunk, TmpDirAllocator,
    frame_keyed_record_into,
};
use parking_lot::Mutex;
use tempfile::TempDir;

use crate::pread::{PreadSlices, ReadRequest, ReadSlice, ReadTarget};
use crate::sort::consolidate_runs::ConsolidateRuns;
use crate::sort::consolidation::ConsolidationRegistry;
use crate::sort::protocol::{
    ConsolidationInput, ConsolidationJob, MemoryChunkErased, SortPhase1Event, SortPhase2Event,
    SpillBlockEvent,
};
use crate::sort::spill_block_compress::SpillBlockCompress;
use crate::sort::spill_decompress::SortSpillDecompress;
use crate::sort::spill_read_planner::SpillReadPlanner;
use crate::sort::spill_write::SpillWrite;
use crate::sort::supply_ledger::SpillSupply;
use crate::sort::write_consolidated_run::WriteConsolidatedRun;

/// A template key whose order is `primary`'s.
pub(crate) fn template_key(primary: u64) -> TemplateKey {
    TemplateKey::new(
        0,
        primary as i32,
        false,
        i32::MAX,
        i32::MAX,
        false,
        0,
        0,
        (0, false),
        0,
        false,
    )
}

/// Kernel-compress one raw block for `codec`, mirroring `SpillBlockCompress`.
pub(crate) fn compress(codec: SpillCodec, raw: &[u8]) -> Vec<u8> {
    SpillBlockCompressor::new(codec, 1).unwrap().compress_block(raw).unwrap()
}

/// A single-block spill run of real keyed records (template-coordinate, so the
/// key is a prefix and any body will do), as `SpillGather` + `SpillBlockCompress`
/// would deliver it. Each body is [`run_body`]`(file_id, key)`, unique across
/// runs even where keys repeat.
pub(crate) fn keyed_run(codec: SpillCodec, file_id: u32, keys: Range<u64>) -> SpillBlockEvent {
    SpillBlockEvent::Block {
        ordinal: 0,
        file_id,
        key_kind: SpillKeyKind::TemplateK40,
        is_last_in_file: true,
        records_ingested_so_far: u64::from(file_id) + 1,
        bytes: compress(codec, &framed_template_records_of(file_id, keys)),
    }
}

/// The body [`keyed_run`] gives the record keyed `key` in run `file_id`.
pub(crate) fn run_body(file_id: u32, key: u64) -> Vec<u8> {
    [file_id.to_le_bytes().as_slice(), key.to_le_bytes().as_slice()].concat()
}

/// `template_key(k)` + `run_body(body_id, k)` framed for every `k`, in one
/// raw buffer.
pub(crate) fn framed_template_records_of(body_id: u32, keys: Range<u64>) -> Vec<u8> {
    let mut raw = Vec::new();
    for k in keys {
        frame_keyed_record_into(&mut raw, &template_key(k), &run_body(body_id, k)).unwrap();
    }
    raw
}

/// [`framed_template_records_of`] with body id 0.
pub(crate) fn framed_template_records(keys: Range<u64>) -> Vec<u8> {
    framed_template_records_of(0, keys)
}

/// Write a complete spill run file (magic, one block of `keyed_run`'s records,
/// trailer) at `path`; returns its length. The body of each record is
/// `run_body(body_id, key)`.
pub(crate) fn write_run_file(
    path: &Path,
    codec: SpillCodec,
    body_id: u32,
    keys: Range<u64>,
) -> u64 {
    let SpillBlockEvent::Block { bytes, .. } = keyed_run(codec, body_id, keys) else {
        unreachable!()
    };
    let mut f = std::fs::OpenOptions::new().write(true).create_new(true).open(path).unwrap();
    f.write_all(fgumi_sort::spill_magic(codec)).unwrap();
    f.write_all(&bytes).unwrap();
    f.write_all(fgumi_sort::spill_trailer(codec)).unwrap();
    drop(f);
    std::fs::metadata(path).unwrap().len()
}

/// Every `(key bytes, body)` in a run file, in file order (40-byte `TemplateK40`
/// key, u32 length, body — the layout `keyed_run` frames).
pub(crate) fn read_run_records(path: &Path, codec: SpillCodec) -> Vec<(Vec<u8>, Vec<u8>)> {
    let bytes = std::fs::read(path).unwrap();
    let body = &bytes[fgumi_sort::spill_magic(codec).len()..];
    let raw: Vec<u8> = fgumi_sort::SpillBlockDecompressor::new()
        .read_blocks(&mut &body[..], codec, 4096)
        .unwrap()
        .concat();
    let mut out = Vec::new();
    let mut at = 0;
    while at < raw.len() {
        let len = u32::from_le_bytes(raw[at + 40..at + 44].try_into().unwrap()) as usize;
        out.push((raw[at..at + 40].to_vec(), raw[at + 44..at + 44 + len].to_vec()));
        at += 44 + len;
    }
    out
}

/// The stable merge of `runs` (in run order): sort by parsed key, ties in run
/// order then file order. This is `reference_sort`'s contract specialised to
/// keyed run records; it never calls the code under test.
pub(crate) fn merged_reference(runs: &[(u32, Range<u64>)]) -> Vec<(Vec<u8>, Vec<u8>)> {
    let mut all: Vec<(Vec<u8>, Vec<u8>)> = runs
        .iter()
        .flat_map(|(id, keys)| {
            keys.clone().map(move |k| {
                let mut key = Vec::new();
                template_key(k).write_to(&mut key).unwrap();
                (key, run_body(*id, k))
            })
        })
        .collect();
    // `sort_by` is stable.
    all.sort_by(|a, b| {
        TemplateKey::read_from(&mut &a.0[..])
            .unwrap()
            .cmp(&TemplateKey::read_from(&mut &b.0[..]).unwrap())
    });
    all
}

// ── Slots, frames and run files at the frame level ──────────────────────────

/// A slot already holding `raw` (one decompressed block of framed records) and
/// at EOF: the merge driver never stalls on it. Mirrors PR 6's `empty_slot` +
/// feeder (`bp_commit_read` then `bp_insert_drain_finalize`).
pub(crate) fn loaded_slot(file_id: u32, slot_id: u32, raw: Vec<u8>) -> Arc<SortMergeSlot> {
    let slot = Arc::new(
        SortMergeSlot::new(
            file_id,
            Arc::new(tempfile::tempfile().unwrap()),
            0,
            0,
            SpillCodec::Zstd,
        )
        .with_slot_id(slot_id),
    );
    if raw.is_empty() {
        slot.bp_commit_read(0, true);
        slot.bp_drain_and_finalize();
    } else {
        slot.bp_commit_read(1, true);
        slot.bp_insert_drain_finalize(0, vec![raw], 1);
    }
    slot
}

/// The minimal coordinate body: tid/pos where the key is extracted from, plus
/// a distinguishing tag, zero-padded to `len` (at least 36) bytes.
pub(crate) fn coordinate_body_len(tid: i32, pos: i32, tag: u32, len: usize) -> Vec<u8> {
    let mut body = vec![0u8; len.max(36)];
    body[0..4].copy_from_slice(&tid.to_le_bytes());
    body[4..8].copy_from_slice(&pos.to_le_bytes());
    body[32..36].copy_from_slice(&tag.to_le_bytes());
    body
}

/// [`coordinate_body_len`] at 36 bytes (the retired tests' body).
pub(crate) fn coordinate_body(tid: i32, pos: i32, tag: u32) -> Vec<u8> {
    coordinate_body_len(tid, pos, tag, 36)
}

/// Records for the framing test, already sorted with distinct keys, body
/// lengths varied 8..=300 bytes (36..=300 for coordinate, whose key is read
/// from the body), so 100-byte blocks see records that overflow alone.
/// `Coordinate` frames the embedded key; `TemplateK40` the 40-byte prefix: the
/// two framing paths.
#[derive(Clone, Copy, Debug)]
pub(crate) enum FrameShape {
    Coordinate,
    TemplateK40,
}

impl FrameShape {
    pub(crate) fn key_kind(self) -> SpillKeyKind {
        match self {
            Self::Coordinate => SpillKeyKind::Coordinate,
            Self::TemplateK40 => SpillKeyKind::TemplateK40,
        }
    }

    /// `n` sorted records as the gather's chunk.
    pub(crate) fn sorted_chunk(self, n: u64) -> MemoryChunkErased {
        let len = |i: u64| 8 + (i * 37 % 293) as usize;
        match self {
            Self::Coordinate => {
                let recs = (0..n)
                    .map(|i| {
                        let body = coordinate_body_len(0, i as i32, i as u32, len(i));
                        (RawCoordinateKey::extract_from_record(&body), body)
                    })
                    .collect();
                MemoryChunkErased::Coordinate(InMemoryChunk::from_owned_records(recs))
            }
            Self::TemplateK40 => {
                let recs =
                    (0..n).map(|i| (template_key(i), vec![(i % 251) as u8; len(i)])).collect();
                MemoryChunkErased::TemplateCoordinate(TemplateMemChunk::K40(
                    InMemoryChunk::from_owned_records(recs),
                ))
            }
        }
    }
}

/// Frames → a spill run file, cut every `block_size` raw bytes regardless of
/// record boundaries (records straddle blocks, as the retired writer's odd
/// block sizes made them). Returns the file length.
pub(crate) fn write_frames_run(
    path: &Path,
    codec: SpillCodec,
    frames: &[Vec<u8>],
    block_size: usize,
) -> u64 {
    let mut comp = SpillBlockCompressor::new(codec, 1).unwrap();
    let mut f = std::fs::OpenOptions::new().write(true).create_new(true).open(path).unwrap();
    f.write_all(fgumi_sort::spill_magic(codec)).unwrap();
    let raw: Vec<u8> = frames.concat();
    for block in raw.chunks(block_size) {
        f.write_all(&comp.compress_block(block).unwrap()).unwrap();
    }
    f.write_all(fgumi_sort::spill_trailer(codec)).unwrap();
    drop(f);
    std::fs::metadata(path).unwrap().len()
}

/// One keyed frame: `[key][u32 len][body]`.
pub(crate) fn keyed_frame(key: &TemplateKey, body: &[u8]) -> Vec<u8> {
    let mut v = Vec::new();
    frame_keyed_record_into(&mut v, key, body).unwrap();
    v
}

/// One embedded-key frame: `[u32 len][body]`.
pub(crate) fn embedded_frame(body: &[u8]) -> Vec<u8> {
    let mut v = Vec::new();
    fgumi_sort::frame_embedded_record_into(&mut v, body).unwrap();
    v
}

/// A run file's frames, re-cut from the decompressed stream: `key_len` key
/// bytes (0 for embedded kinds), a u32 length, the body.
pub(crate) fn read_frames(path: &Path, codec: SpillCodec, key_len: usize) -> Vec<Vec<u8>> {
    let bytes = std::fs::read(path).unwrap();
    let raw: Vec<u8> = fgumi_sort::SpillBlockDecompressor::new()
        .read_blocks(&mut &bytes[fgumi_sort::spill_magic(codec).len()..], codec, 4096)
        .unwrap()
        .concat();
    let (mut out, mut at) = (Vec::new(), 0);
    while at < raw.len() {
        let len =
            u32::from_le_bytes(raw[at + key_len..at + key_len + 4].try_into().unwrap()) as usize;
        out.push(raw[at..at + key_len + 4 + len].to_vec());
        at += key_len + 4 + len;
    }
    out
}

// ── KChain: the consolidation loop, pumped by hand ───────────────────────────

/// [`KChain`] options.
#[derive(Clone)]
pub(crate) struct KOpts {
    pub limit: usize,
    pub total_memory: u64,
    pub phased: Option<Arc<PhasedCap>>,
    pub codec: SpillCodec,
    pub block_size: usize,
    pub per_call: (usize, usize),
}

impl Default for KOpts {
    fn default() -> Self {
        Self {
            limit: 2,
            total_memory: 64 << 20,
            phased: None,
            codec: SpillCodec::Zstd,
            block_size: 64 << 10,
            per_call: (usize::MAX, usize::MAX),
        }
    }
}

/// A recorded job: `(job_id, slot ids, input paths, a clone of each slot's file)`.
pub(crate) type RecordedJob =
    (u32, Vec<u32>, Vec<PathBuf>, Vec<Arc<dyn fgumi_bam_io::pread::PositionalSource>>);

/// `SpillWrite → SpillReadPlanner ⇒ (.0 PreadSpillSlices → 2 × SortSpillDecompress)
/// (.1 collected phase-2 events) (.2 ConsolidateRuns → SpillBlockCompress →
/// WriteConsolidatedRun)`, real steps through `StepProbe`s pumped on the test
/// thread, so every interleaving is deterministic.
pub(crate) struct KChain {
    pub dir: TempDir,
    sw: (SpillWrite, StepProbe<SpillWrite>),
    sw_out: BranchInputHandle<SortPhase1Event>,
    pl: (SpillReadPlanner, StepProbe<SpillReadPlanner>),
    pl_reads: BranchInputHandle<ReadRequest>,
    pl_events: BranchInputHandle<SortPhase2Event>,
    pl_jobs: BranchInputHandle<ConsolidationJob>,
    pr: (PreadSlices, StepProbe<PreadSlices>),
    pr_out: BranchInputHandle<ReadSlice>,
    dc: Vec<(SortSpillDecompress, StepProbe<SortSpillDecompress>)>,
    next_dc: usize,
    cr: (ConsolidateRuns, StepProbe<ConsolidateRuns>),
    cr_out: BranchInputHandle<SpillBlockEvent>,
    sb: (SpillBlockCompress, StepProbe<SpillBlockCompress>),
    sb_out: BranchInputHandle<SpillBlockEvent>,
    wr: (WriteConsolidatedRun, StepProbe<WriteConsolidatedRun>),
    pub registry: Arc<ConsolidationRegistry>,
    pub cons_demand: Arc<MergeDemand>,
    pub merge_demand: Arc<MergeDemand>,
    /// Every read request in issue order: `(slot_id, len)`.
    pub requests: Vec<(u32, u32)>,
    /// Every job as forwarded on branch 2.
    pub jobs: Vec<RecordedJob>,
    /// Index into `requests` at which each job's retire passed the planner.
    pub retired_at: Vec<(u32, usize)>,
    pub phase2: Vec<SortPhase2Event>,
    pub codec: SpillCodec,
    /// Called once per job when it is forwarded on branch 2 (before any of its
    /// reads land).
    pub on_job: Option<Box<dyn FnMut(&ConsolidationJob)>>,
}

impl KChain {
    pub(crate) fn new(opts: KOpts) -> Self {
        let KOpts { limit, total_memory, phased, codec, block_size, per_call } = opts;
        let dir = TempDir::new().unwrap();
        let alloc = TmpDirAllocator::new(vec![dir.path().to_path_buf()]).unwrap();
        let registry = ConsolidationRegistry::new();
        let supply = SpillSupply::new();
        let cons_demand = Arc::new(MergeDemand::new());
        let merge_demand = Arc::clone(&supply.demand);
        let all_announced = Arc::new(AtomicBool::new(false));
        let cap = phased.clone().map(StepCap::Phased);
        let slices = SliceBufferPool::new(8);
        let sw = SpillWrite::new(Arc::new(Mutex::new(alloc)), codec, 1 << 20, Arc::new(Vec::new()))
            .with_consolidation(limit, Arc::clone(&registry));
        let pl = SpillReadPlanner::new(total_memory, 2, 2, 1 << 20, &supply, Arc::clone(&slices))
            .with_consolidation(
                Arc::clone(&cons_demand),
                Arc::clone(&all_announced),
                phased.clone(),
            );
        let pr = PreadSlices::spill(slices, 1 << 20).with_step_cap(cap.clone());
        // Two worker copies of one step: they share its registry, as the
        // pool's clones do.
        let dc0 = SortSpillDecompress::new(1 << 20, &supply)
            .with_consolidation_demand(Arc::clone(&cons_demand), all_announced)
            .with_step_cap(cap.clone());
        let dc = vec![dc0.new_worker_copy(), dc0];
        let cr = ConsolidateRuns::new(block_size, 1 << 20)
            .with_merge_demand(Arc::clone(&cons_demand))
            .with_supply_ledger(Arc::clone(&supply.ledger))
            .with_per_call(per_call.0, per_call.1);
        let sb = SpillBlockCompress::new(codec, 1, 1 << 20)
            .named("SpillBlockCompressConsolidation")
            .with_step_cap(cap);
        let wr = WriteConsolidatedRun::new(codec, Arc::clone(&registry));
        let mut swp = StepProbe::new(&sw);
        let mut plp = StepProbe::new(&pl);
        let mut prp = StepProbe::new(&pr);
        let mut crp = StepProbe::new(&cr);
        let mut sbp = StepProbe::new(&sb);
        let wrp = StepProbe::new(&wr);
        Self {
            sw_out: swp.take_output(0),
            pl_reads: plp.take_output(0),
            pl_events: plp.take_output(1),
            pl_jobs: plp.take_output(2),
            pr_out: prp.take_output(0),
            cr_out: crp.take_output(0),
            sb_out: sbp.take_output(0),
            dc: dc
                .into_iter()
                .map(|d| {
                    let p = StepProbe::new(&d);
                    (d, p)
                })
                .collect(),
            next_dc: 0,
            sw: (sw, swp),
            pl: (pl, plp),
            pr: (pr, prp),
            cr: (cr, crp),
            sb: (sb, sbp),
            wr: (wr, wrp),
            dir,
            registry,
            cons_demand,
            merge_demand,
            requests: Vec::new(),
            jobs: Vec::new(),
            retired_at: Vec::new(),
            phase2: Vec::new(),
            codec,
            on_job: None,
        }
    }

    pub(crate) fn feed(&self, ev: SpillBlockEvent) {
        self.sw.1.push_input(ev);
    }

    pub(crate) fn close(&self) {
        self.sw.1.close_input();
    }

    /// Start a job directly at the planner (no `SpillWrite` run stack): what
    /// `SpillWrite` would emit for `inputs` (file ids 0..n), with `output`
    /// registered.
    pub(crate) fn start_job(
        &self,
        job_id: u32,
        key_kind: SpillKeyKind,
        inputs: &[PathBuf],
        output: &Path,
    ) {
        self.registry.begin(job_id, output.to_path_buf());
        let inputs = inputs
            .iter()
            .enumerate()
            .map(|(i, p)| ConsolidationInput {
                file_id: i as u32,
                path: p.clone(),
                len: std::fs::metadata(p).unwrap().len(),
            })
            .collect();
        self.pl.1.push_input(SortPhase1Event::ConsolidationStart { job_id, inputs, key_kind });
    }

    /// Live runs in `SpillWrite`'s stack.
    pub(crate) fn live_runs(&self) -> usize {
        self.sw.0.live_runs_for_test()
    }

    /// One pass: run every step once, then forward every output. Returns
    /// whether any step progressed or any item moved.
    pub(crate) fn pump(&mut self) -> io::Result<bool> {
        let p = StepOutcome::Progress;
        let mut moved = self.sw.1.try_run(&mut self.sw.0)? == p;
        moved |= self.pl.1.try_run(&mut self.pl.0)? == p;
        moved |= self.pr.1.try_run(&mut self.pr.0)? == p;
        for (d, probe) in &mut self.dc {
            moved |= probe.try_run(d)? == p;
        }
        moved |= self.cr.1.try_run(&mut self.cr.0)? == p;
        moved |= self.sb.1.try_run(&mut self.sb.0)? == p;
        moved |= self.wr.1.try_run(&mut self.wr.0)? == p;
        while let Some(ev) = self.sw_out.pop() {
            if let SortPhase1Event::ConsolidationRetire { job_id } = &ev {
                self.retired_at.push((*job_id, self.requests.len()));
            }
            self.pl.1.push_input(ev);
            moved = true;
        }
        while let Some(r) = self.pl_reads.pop() {
            let ReadTarget::Spill { slot, .. } = &r.target else { panic!("spill reads only") };
            self.requests.push((slot.slot_id(), r.len));
            self.pr.1.push_input(r);
            moved = true;
        }
        while let Some(ev) = self.pl_events.pop() {
            self.phase2.push(ev);
            moved = true;
        }
        while let Some(job) = self.pl_jobs.pop() {
            self.jobs.push((
                job.job_id,
                job.slots.iter().map(|s| s.slot_id()).collect(),
                job.input_paths.clone(),
                job.slots.iter().map(|s| Arc::clone(s.source())).collect(),
            ));
            if let Some(f) = &mut self.on_job {
                f(&job);
            }
            self.cr.1.push_input(job);
            moved = true;
        }
        while let Some(s) = self.pr_out.pop() {
            self.dc[self.next_dc].1.push_input(s); // alternate the copies
            self.next_dc = (self.next_dc + 1) % self.dc.len();
            moved = true;
        }
        while let Some(b) = self.cr_out.pop() {
            self.sb.1.push_input(b);
            moved = true;
        }
        while let Some(b) = self.sb_out.pop() {
            self.wr.1.push_input(b);
            moved = true;
        }
        Ok(moved)
    }

    /// Pump until two consecutive passes move nothing (the second pass is the
    /// "one more pass" every retired entry needs to be dropped).
    ///
    /// # Panics
    ///
    /// After 1 000 000 passes.
    pub(crate) fn settle(&mut self) -> io::Result<()> {
        let mut idle = 0;
        for _ in 0..1_000_000 {
            idle = if self.pump()? { 0 } else { idle + 1 };
            if idle == 2 {
                return Ok(());
            }
        }
        panic!("KChain did not settle");
    }

    /// The decompress registry's slot ids, the planner's planned slot ids, and
    /// both demands' published ids.
    pub(crate) fn live_slot_ids(&self) -> Vec<u32> {
        let mut ids: Vec<u32> =
            self.dc.iter().flat_map(|(d, _)| d.registered_slot_ids_for_test()).collect();
        ids.extend(self.pl.0.planned_slot_ids_for_test());
        for d in [&self.cons_demand, &self.merge_demand] {
            ids.extend([d.awaited(), d.predicted(), d.frontier()].into_iter().flatten());
        }
        ids
    }

    /// Records of every surviving run named by a phase-2 `SpillReady`, per run.
    pub(crate) fn surviving_runs(&self) -> Vec<Vec<(Vec<u8>, Vec<u8>)>> {
        self.phase2
            .iter()
            .filter_map(|e| match e {
                SortPhase2Event::SpillReady { path, .. } => {
                    Some(read_run_records(path, self.codec))
                }
                _ => None,
            })
            .collect()
    }
}
