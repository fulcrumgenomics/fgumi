//! Unit tests for the split arena front's steps.

use std::sync::mpsc;
use std::time::Duration;

use fgumi_pipeline_core::builder::{Pipeline, PipelineConfig};
use fgumi_pipeline_core::step::InputHandle as _;
use fgumi_pipeline_core::testing::StepProbe;
use fgumi_raw_bam::testutil::make_bam_bytes;
use fgumi_sort::{DroppedLaneViolation, KeyTypesSpec, SegmentedBuf};
use noodles::sam::Header;

use super::*;
use crate::sort::{ArenaSortStrategy as _, QuerynameStrategy};

/// A mapped, paired template record body named `rd_{i}`.
fn body(i: usize) -> Vec<u8> {
    let pos = i32::try_from(1000 + (i * 7919) % 50_000).unwrap();
    make_bam_bytes(0, pos, 0x3, format!("rd_{i}").as_bytes(), &[], 20, 0, pos + 50, &[])
}

/// `[u32 len][body]` frames.
fn frames(bodies: &[Vec<u8>]) -> Vec<u8> {
    let mut out = Vec::new();
    for b in bodies {
        out.extend_from_slice(&u32::try_from(b.len()).unwrap().to_le_bytes());
        out.extend_from_slice(b);
    }
    out
}

/// A minimal binary BAM header with one reference.
fn bam_header() -> Vec<u8> {
    let mut h = b"BAM\x01".to_vec();
    h.extend_from_slice(&0u32.to_le_bytes());
    h.extend_from_slice(&1u32.to_le_bytes());
    h.extend_from_slice(&2u32.to_le_bytes());
    h.extend_from_slice(b"r\0");
    h.extend_from_slice(&1_000_000u32.to_le_bytes());
    h
}

/// One arena holding `prefix` zero bytes then `data`; returns the arena and the
/// offset `data` starts at.
fn arena_with(prefix: usize, data: &[u8]) -> (Arc<PooledSegmentedBuf>, u64) {
    let mut buf = SegmentedBuf::with_capacity(0, prefix + data.len() + 64);
    if prefix > 0 {
        buf.extend_from_slice(&vec![0u8; prefix]);
    }
    let at = buf.extend_from_slice(data) as u64;
    (Arc::new(PooledSegmentedBuf::unpooled(buf)), at)
}

fn block(
    arena: &Arc<PooledSegmentedBuf>,
    ordinal: u64,
    offset: u64,
    len: usize,
    run_seq: u32,
    seals_to_spill: bool,
) -> InflatedBlock {
    InflatedBlock {
        arena: Arc::clone(arena),
        ordinal,
        offset,
        len: u32::try_from(len).unwrap(),
        is_last_of_run: true,
        run_seq,
        seals_to_spill,
    }
}

fn seed() -> TemplateKeySeed {
    TemplateKeySeed::from_header(&Header::default(), None, KeyTypesSpec::Auto)
}

/// Ordinals are dense from 0 across runs, every run ends in exactly one batch
/// carrying its `RunEnd` — even an empty one when the run's records fill whole
/// batches — and record indices are run-relative.
#[test]
fn ordinals_are_dense_with_every_terminal_batch() {
    let mut step = FindTemplateBoundaries::new(seed(), 1 << 20).with_batch_records(64);
    let run0: Vec<Vec<u8>> = (0..128).map(body).collect();
    let data0 = [bam_header(), frames(&run0)].concat();
    let (a0, at0) = arena_with(0, &data0);
    assert_eq!(step.ingest_block(&block(&a0, 0, at0, data0.len(), 0, true)).unwrap(), 128);
    let run1: Vec<Vec<u8>> = (128..198).map(body).collect();
    let data1 = frames(&run1);
    let (a1, at1) = arena_with(FRONT_REGION, &data1);
    assert_eq!(step.ingest_block(&block(&a1, 1, at1, data1.len(), 1, false)).unwrap(), 70);
    let batches: Vec<_> = step.pending.drain(..).collect();
    let got: Vec<(u64, u32, u64, bool)> = batches
        .iter()
        .map(|b| (b.ordinal, b.n_records, b.first_record_index, b.run_end.is_some()))
        .collect();
    assert_eq!(
        got,
        vec![
            (0, 64, 0, false),
            (1, 64, 64, false),
            (2, 0, 128, true),
            (3, 64, 0, false),
            (4, 6, 64, true),
        ]
    );
    let ends: Vec<RunEnd> = batches.iter().filter_map(|b| b.run_end).collect();
    assert!(ends[0].seals_to_spill && !ends[1].seals_to_spill);
    assert_eq!((ends[0].records_ingested_so_far, ends[1].records_ingested_so_far), (128, 198));
    assert!(step.finalized, "the residual run finishes the scan");
}

/// A run whose bytes end mid-record carries that partial record into the next
/// run's front region, where it becomes record 0 of the next run's first batch.
#[test]
fn a_straddler_is_record_zero_of_the_next_run() {
    let mut step = FindTemplateBoundaries::new(seed(), 1 << 20).with_batch_records(64);
    let bodies: Vec<Vec<u8>> = (0..10).map(body).collect();
    let all = frames(&bodies);
    let cut = all.len() - 20; // the last record straddles
    let data0 = [bam_header(), all[..cut].to_vec()].concat();
    let (a0, at0) = arena_with(0, &data0);
    assert_eq!(step.ingest_block(&block(&a0, 0, at0, data0.len(), 0, true)).unwrap(), 9);
    let (a1, at1) = arena_with(FRONT_REGION, &all[cut..]);
    assert_eq!(step.ingest_block(&block(&a1, 1, at1, all.len() - cut, 1, false)).unwrap(), 1);
    let batches: Vec<_> = step.pending.drain(..).collect();
    let last = batches.last().expect("run 1's batch");
    assert_eq!((last.n_records, last.first_record_index), (1, 0));
    let len = usize::try_from(last.byte_end - last.byte_start).unwrap();
    let start = usize::try_from(last.byte_start).unwrap();
    assert_eq!(&last.arena.slice(start, len)[4..], &bodies[9][..], "the carried record");
}

/// A run's published bytes are its complete records' frames only: run 0's
/// BAM header and the partial record carried into run 1 are not part of run
/// 0's mean (they would inflate the mean that sizes run 1's cap), and the
/// carried record counts toward run 1.
#[test]
fn run_bytes_count_only_complete_records() {
    let mut step = FindTemplateBoundaries::new(seed(), 1 << 20);
    let bodies: Vec<Vec<u8>> = (0..10).map(body).collect();
    let all = frames(&bodies);
    let cut = all.len() - 20; // the last record straddles
    let text = vec![b'@'; 10_000];
    let mut header = b"BAM\x01".to_vec();
    header.extend_from_slice(&u32::try_from(text.len()).unwrap().to_le_bytes());
    header.extend_from_slice(&text);
    header.extend_from_slice(&bam_header()[8..]);
    let data0 = [header, all[..cut].to_vec()].concat();
    let (a0, at0) = arena_with(0, &data0);
    assert_eq!(step.ingest_block(&block(&a0, 0, at0, data0.len(), 0, true)).unwrap(), 9);
    let (a1, at1) = arena_with(FRONT_REGION, &all[cut..]);
    assert_eq!(step.ingest_block(&block(&a1, 1, at1, all.len() - cut, 1, false)).unwrap(), 1);
    let ends: Vec<RunEnd> = step.pending.drain(..).filter_map(|b| b.run_end).collect();
    let last = frames(&bodies[9..]).len() as u64;
    assert_eq!(ends[0].run_bytes, all.len() as u64 - last, "run 0: nine whole records");
    assert_eq!(ends[1].run_bytes, last, "run 1: the carried record");
}

/// A run sealed to spill before the BAM header ends — the header outgrew one
/// run's arena cap — names that cause and the fix, not "input too short";
/// the same short span as the final run is genuinely truncated input.
#[test]
fn a_header_larger_than_a_run_says_so() {
    let mut header = b"BAM\x01".to_vec();
    header.extend_from_slice(&50_000u32.to_le_bytes());
    header.extend_from_slice(&[b'@'; 1_000]);
    let (a, at) = arena_with(0, &header);
    let mut spill = FindTemplateBoundaries::new(seed(), 1 << 20);
    let err = spill.ingest_block(&block(&a, 0, at, header.len(), 0, true)).unwrap_err();
    assert_eq!(err.kind(), io::ErrorKind::InvalidData);
    assert_eq!(
        err.to_string(),
        "FindTemplateBoundaries: the BAM header does not fit in one run's memory budget (1008 \
         bytes read without reaching its end) — raise --max-memory"
    );
    let mut last = FindTemplateBoundaries::new(seed(), 1 << 20);
    let err = last.ingest_block(&block(&a, 0, at, header.len(), 0, false)).unwrap_err();
    assert_eq!(
        err.to_string(),
        "FindTemplateBoundaries: input too short to contain a valid BAM header"
    );
}

/// The boundary step and the run sort are `Detached`: the framework never
/// copies them, and a copy would silently drop the run ledger (and the run
/// sort's reserve budget), so asking for one panics instead.
#[rstest::rstest]
#[should_panic(expected = "only Parallel steps need to override this")]
fn detached_split_steps_refuse_a_worker_copy(#[values(true, false)] boundaries: bool) {
    if boundaries {
        let _ = FindTemplateBoundaries::new(seed(), 1 << 20).new_worker_copy();
    } else {
        let _ = SortRun::<TemplateKeys>::new(1, 1 << 20).new_worker_copy();
    }
}

/// A final record cut short at the end of input is the truncation error the
/// record-input front raises too.
#[test]
fn a_truncated_final_record_is_an_error() {
    let mut step = FindTemplateBoundaries::new(seed(), 1 << 20);
    let all = frames(&(0..3).map(body).collect::<Vec<_>>());
    let data = [bam_header(), all[..all.len() - 5].to_vec()].concat();
    let (a, at) = arena_with(0, &data);
    let err = step.ingest_block(&block(&a, 0, at, data.len(), 0, false)).unwrap_err();
    assert_eq!(err.to_string(), "truncated BAM: final record incomplete at end of input");
}

/// A frame whose `block_size` is below the 32 fixed BAM bytes cannot be a
/// record: the scan rejects it with `InvalidData` before provisioning or any
/// key extractor indexes its fixed fields.
#[test]
fn scan_rejects_a_frame_shorter_than_a_bam_record() {
    let mut step = FindTemplateBoundaries::new(seed(), 1 << 20);
    let data = [bam_header(), frames(&[body(0), vec![0u8; 12]])].concat();
    let (a, at) = arena_with(0, &data);
    let err = step.ingest_block(&block(&a, 0, at, data.len(), 0, false)).unwrap_err();
    assert_eq!(err.kind(), io::ErrorKind::InvalidData);
    assert_eq!(
        err.to_string(),
        "FindTemplateBoundaries: block_size 12 is shorter than a BAM record (32 bytes)"
    );
}

/// One `RefBatch` for `run_end`-less tests.
fn ref_batch(
    ctx: &Arc<TemplateKeyContext>,
    arena: &Arc<PooledSegmentedBuf>,
    ordinal: u64,
    violation_at: Option<u64>,
    run_end: Option<RunEnd>,
) -> TemplateRefBatch {
    let now = Instant::now();
    TemplateRefBatch {
        arena: Arc::clone(arena),
        ctx: Some(Arc::clone(ctx)),
        ordinal,
        run_seq: 0,
        first_record_index: ordinal * 64,
        refs: Some(ExtractedRefs::K24(Box::new([]))),
        violation: violation_at.map(|i| KeyViolation {
            record_index: i,
            violation: DroppedLaneViolation::Mi,
            name: format!("rd_{i}"),
        }),
        keyed_at: now,
        run_end,
    }
}

/// Emits its items, one per call, in the order given — whatever their
/// ordinals.
struct OutOfOrderSource {
    items: Vec<TemplateRefBatch>,
    held: HeldSlot<Unpushed<TemplateRefBatch>>,
}

impl Step for OutOfOrderSource {
    type Input = ();
    type Outputs = OrderedBytesSingle<TemplateRefBatch>;

    fn profile(&self) -> StepProfile {
        StepProfile {
            name: "OutOfOrderSource",
            kind: StepKind::Exclusive,
            sticky: true,
            output_queues: vec![QueueSpec::ByteBounded { limit_bytes: 1 << 20 }],
            branch_ordering: vec![BranchOrdering::None],
        }
    }

    fn try_run(&mut self, ctx: &mut StepCtx<'_, Self>) -> io::Result<StepOutcome> {
        if let Some(unpushed) = self.held.take()
            && let Err(again) = ctx.outputs.retry(unpushed)
        {
            self.held.put(again);
            return Ok(StepOutcome::Contention);
        }
        if self.items.is_empty() {
            return Ok(StepOutcome::Finished);
        }
        let item = self.items.remove(0);
        if let Err(unpushed) = ctx.outputs.push(item) {
            self.held.put(unpushed);
        }
        Ok(StepOutcome::Progress)
    }
}

/// A `Parallel` pass-through — the extraction step's place in the chain — so
/// its output edge carries the production `ByItemOrdinal` reorder stage.
#[derive(Default)]
struct ParallelPassThrough {
    held: HeldSlot<Unpushed<TemplateRefBatch>>,
}

impl Step for ParallelPassThrough {
    type Input = TemplateRefBatch;
    type Outputs = OrderedBytesSingle<TemplateRefBatch>;

    fn profile(&self) -> StepProfile {
        StepProfile {
            name: "ParallelPassThrough",
            kind: StepKind::Parallel,
            sticky: false,
            output_queues: vec![QueueSpec::ByteBounded { limit_bytes: 1 << 20 }],
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
        let Some(item) = ctx.input.pop() else {
            return Ok(if ctx.input.is_drained() {
                StepOutcome::Finished
            } else {
                StepOutcome::NoProgress
            });
        };
        if let Err(unpushed) = ctx.outputs.push(item) {
            self.held.put(unpushed);
        }
        Ok(StepOutcome::Progress)
    }

    fn new_worker_copy(&self) -> Self {
        Self::default()
    }
}

/// Drops every chunk event.
struct NullSink;

impl Step for NullSink {
    type Input = SortChunkEvent;
    type Outputs = ();

    fn profile(&self) -> StepProfile {
        StepProfile {
            name: "NullSink",
            kind: StepKind::Serial,
            sticky: false,
            output_queues: vec![],
            branch_ordering: vec![],
        }
    }

    fn try_run(&mut self, ctx: &mut StepCtx<'_, Self>) -> io::Result<StepOutcome> {
        Ok(match ctx.input.pop() {
            Some(_) => StepOutcome::Progress,
            None if ctx.input.is_drained() => StepOutcome::Finished,
            None => StepOutcome::NoProgress,
        })
    }
}

/// `SortRun` behind the production `ByItemOrdinal` reorder sees batches in
/// input order whatever order they were extracted in, so it reports the
/// lowest-indexed violation — record 7, in batch 0, though batch 5 (record
/// 300) arrived first — and only that one.
#[test]
fn sort_run_reports_the_lowest_indexed_violation_whatever_the_arrival_order() {
    let first = body(0);
    let ctx = seed().provision(&first);
    let (arena, _) = arena_with(0, b"x");
    let end = RunEnd {
        seals_to_spill: false,
        records_ingested_so_far: 384,
        run_bytes: 0,
        scan_end: Instant::now(),
    };
    let mut items = vec![ref_batch(&ctx, &arena, 5, Some(300), Some(end))];
    items.extend((1..5).map(|o| ref_batch(&ctx, &arena, o, None, None)));
    items.push(ref_batch(&ctx, &arena, 0, Some(7), None));
    let builder = Pipeline::builder();
    builder
        .chain(OutOfOrderSource { items, held: HeldSlot::new() })
        .chain(ParallelPassThrough::default())
        .chain(SortTemplateRun::new(1, 1 << 20))
        .chain(NullSink)
        .into_sink_marker();
    let err = builder
        .build()
        .unwrap()
        .run(PipelineConfig { threads: 2, ..Default::default() })
        .expect_err("a violation fails the sort");
    let msg = format!("{err:#}");
    assert!(msg.contains("rd_7"), "names the first violation in input order: {msg}");
    assert!(!msg.contains("rd_300"), "and only it: {msg}");
}

/// A key batch over three real records.
fn key_batch(ordinal: u64) -> TemplateKeyBatch {
    key_batch_with(ordinal, 3, 0)
}

/// A key batch over `n` real records, in run `run_seq`.
fn key_batch_with(ordinal: u64, n: u32, run_seq: u32) -> TemplateKeyBatch {
    let bodies: Vec<Vec<u8>> = (0..n as usize).map(body).collect();
    let data = frames(&bodies);
    let (arena, at) = arena_with(16, &data);
    TemplateKeyBatch {
        arena,
        ctx: Some(seed().provision(&bodies[0])),
        ordinal,
        run_seq,
        first_record_index: 0,
        n_records: n,
        byte_start: at,
        byte_end: at + data.len() as u64,
        run_end: None,
    }
}

/// A full cap refuses before the pop: the batch stays queued, and once the
/// permit frees the same batch is extracted.
#[test]
fn refused_extract_leaves_the_batch_queued() {
    let cap = PhaseCap::new("sort-phase1", 1);
    let mut step = ExtractTemplateKeys::new(1 << 20).with_phase_cap(Some(Arc::clone(&cap)));
    let mut probe = StepProbe::new(&step);
    probe.push_input(key_batch(0));
    let held = cap.try_acquire().expect("the only permit");
    assert_eq!(probe.try_run(&mut step).unwrap(), StepOutcome::Capped);
    assert!(!probe.input_is_empty(), "the refused batch is still queued");
    drop(held);
    assert_eq!(probe.try_run(&mut step).unwrap(), StepOutcome::Progress);
    let out = probe.take_output::<TemplateRefBatch>(0);
    let batch = out.pop().expect("the extracted batch");
    assert_eq!((batch.ordinal, batch.len()), (0, 3));
}

/// An empty poll takes no permit: with the only permit held elsewhere, an open
/// empty input is `NoProgress` and a drained one `Finished`, not `Capped`, and
/// neither counts as a refusal.
#[test]
fn extract_does_not_take_a_permit_without_input() {
    let cap = PhaseCap::new("sort-phase1", 1);
    let mut step = ExtractTemplateKeys::new(1 << 20).with_phase_cap(Some(Arc::clone(&cap)));
    let probe = StepProbe::new(&step);
    let _held = cap.try_acquire().expect("the only permit");
    assert_eq!(probe.try_run(&mut step).unwrap(), StepOutcome::NoProgress);
    probe.close_input();
    assert_eq!(probe.try_run(&mut step).unwrap(), StepOutcome::Finished);
    assert_eq!(cap.refused(), 0, "no refusal was counted");
}

/// The extraction step owes every phase-1 pool step's admission contract.
#[test]
fn extract_follows_the_admission_contract() {
    fgumi_pipeline_core::testing::assert_admission_contract(
        "sort-phase1",
        |cap| ExtractTemplateKeys::new(1 << 20).with_phase_cap(cap),
        key_batch(0),
    );
}

/// What a probe extractor's seal saw.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
struct SealSeen {
    cap_active: usize,
    pool_threads: usize,
}

/// The probe extractor's shared state: the cap it reports on and the channel
/// its seal reports to.
struct ProbeShared {
    cap: Arc<PhaseCap>,
    entered: parking_lot::Mutex<mpsc::Sender<SealSeen>>,
}

/// An extractor whose seal reports the cap's active permits and the pool width.
struct ProbeKeys;

impl ArenaKeyExtractor for ProbeKeys {
    const STEP_NAMES: [&'static str; 3] = ["ProbeBoundaries", "ProbeExtract", "ProbeRun"];
    type Seed = ();
    type Ctx = Arc<ProbeShared>;
    type Refs = u32;
    type Accum = u32;
    type Scratch = ();

    fn provision((): &(), _first_body: &[u8]) -> Self::Ctx {
        unreachable!("tests build the context")
    }

    fn ref_width(_ctx: &Self::Ctx) -> usize {
        8
    }

    fn extract_batch(
        _ctx: &Self::Ctx,
        (): &mut (),
        _arena: &PooledSegmentedBuf,
        _byte_start: u64,
        _byte_end: u64,
        n: u32,
        _first_index: u64,
    ) -> io::Result<(u32, Option<KeyViolation>)> {
        Ok((n, None))
    }

    fn refs_len(refs: &u32) -> usize {
        *refs as usize
    }

    fn refs_byte_len(refs: &u32) -> usize {
        *refs as usize * 8
    }

    fn accum_for(_ctx: &Self::Ctx) -> u32 {
        0
    }

    fn reserve(_acc: &mut u32, _n: usize) {}

    fn append(acc: &mut u32, refs: u32) {
        *acc += refs;
    }

    fn sorts_in_parallel(_acc: &u32, sort_threads: usize) -> bool {
        sort_threads > 1
    }

    fn seal(
        _acc: &mut u32,
        ctx: &Self::Ctx,
        _arena: Arc<PooledSegmentedBuf>,
        sort_threads: usize,
        pool: &BoundedSortPool,
    ) -> MemoryChunkErased {
        let seen = pool.bounded(sort_threads).install(|| SealSeen {
            cap_active: ctx.cap.active(),
            pool_threads: rayon::current_num_threads(),
        });
        ctx.entered.lock().send(seen).unwrap();
        Self::empty_chunk()
    }

    fn empty_chunk() -> MemoryChunkErased {
        MemoryChunkErased::Coordinate(InMemoryChunk::default())
    }

    fn pool_name() -> &'static str {
        "probe-sort"
    }
}

/// A probe run's only batch, ending the run.
fn probe_batch(shared: &Arc<ProbeShared>) -> RefBatch<ProbeKeys> {
    let (arena, _) = arena_with(0, b"x");
    let now = Instant::now();
    RefBatch {
        arena,
        ctx: Some(Arc::clone(shared)),
        ordinal: 0,
        run_seq: 0,
        first_record_index: 0,
        refs: Some(4),
        violation: None,
        keyed_at: now,
        run_end: Some(RunEnd {
            seals_to_spill: false,
            records_ingested_so_far: 4,
            run_bytes: 0,
            scan_end: now,
        }),
    }
}

/// With a phase-1 cap of 2, the run sort first reserves the whole cap: while a
/// pool step holds one permit the reservation is pending — new admissions are
/// refused, though a permit is free — and the sort has not started; once the
/// pool step releases, the sort runs holding both permits on a 2-thread pool.
#[test]
fn run_sort_reserves_the_whole_cap() {
    let cap = PhaseCap::new("sort-phase1", 2);
    let signal = fgumi_pipeline_core::signal::PipelineSignal::new();
    cap.bind_signal_for_test(&signal);
    let (tx, entered) = mpsc::channel();
    let shared =
        Arc::new(ProbeShared { cap: Arc::clone(&cap), entered: parking_lot::Mutex::new(tx) });
    let pool_step = cap.try_acquire().expect("a pool step holds one permit");

    let run_cap = Arc::clone(&cap);
    let batch = probe_batch(&shared);
    let sorter = std::thread::spawn(move || {
        let mut step = SortRun::<ProbeKeys>::new(2, 1 << 20).with_phase_cap(Some(run_cap));
        let probe = StepProbe::new(&step);
        probe.push_input(batch);
        probe.try_run(&mut step).map(|_| ())
    });
    let watchdog = Instant::now();
    while cap.pending_reservations() == 0 {
        assert!(watchdog.elapsed() < Duration::from_secs(5), "watchdog: no reservation");
        std::thread::yield_now();
    }
    assert_eq!(cap.active(), 1);
    assert!(cap.try_acquire().is_none(), "a pending reservation refuses new admissions");
    assert!(entered.try_recv().is_err(), "the sort waits for the whole cap");
    drop(pool_step);
    let seen = entered.recv_timeout(Duration::from_secs(5)).expect("the sort ran");
    assert_eq!(seen, SealSeen { cap_active: 2, pool_threads: 2 });
    sorter.join().unwrap().unwrap();
    assert_eq!(cap.active(), 0, "released after the sort");
}

/// A run sort waiting for the cap gives up when the pipeline is cancelled.
#[test]
fn run_sort_reservation_observes_cancel() {
    let cap = PhaseCap::new("sort-phase1", 2);
    let signal = fgumi_pipeline_core::signal::PipelineSignal::new();
    cap.bind_signal_for_test(&signal);
    let (tx, _entered) = mpsc::channel();
    let shared =
        Arc::new(ProbeShared { cap: Arc::clone(&cap), entered: parking_lot::Mutex::new(tx) });
    let _forever = cap.try_acquire().expect("held for the whole test");
    let run_cap = Arc::clone(&cap);
    let batch = probe_batch(&shared);
    let (done_tx, done) = mpsc::channel();
    std::thread::spawn(move || {
        let mut step = SortRun::<ProbeKeys>::new(2, 1 << 20).with_phase_cap(Some(run_cap));
        let probe = StepProbe::new(&step);
        probe.push_input(batch);
        let _ = done_tx.send(probe.try_run(&mut step).map_err(|e| e.kind()));
    });
    while cap.pending_reservations() == 0 {
        std::thread::yield_now();
    }
    signal.cancel();
    let outcome = done.recv_timeout(Duration::from_secs(5)).expect("woken by the cancel");
    assert_eq!(outcome, Err(io::ErrorKind::Interrupted));
}

/// The three steps' profiles: boundaries and run sort on the `sort-coord`
/// driver, extraction a byte-bounded, ordinal-reordered pool step off the
/// reader worker.
#[test]
fn split_step_profiles() {
    let boundaries = FindTemplateBoundaries::new(seed(), 4096);
    let extract = ExtractTemplateKeys::new(4096);
    let run = SortTemplateRun::new(4, 4096);
    for (name, kind, group) in [
        (boundaries.profile().name, boundaries.profile().kind, boundaries.detached_group()),
        (run.profile().name, run.profile().kind, run.detached_group()),
    ] {
        assert_eq!(kind, StepKind::Detached, "{name}");
        assert_eq!(group, DetachedGroup::Shared(crate::sort::SORT_COORD_GROUP), "{name}");
    }
    let p = extract.profile();
    assert_eq!(p.name, "ExtractTemplateKeys");
    assert_eq!(p.kind, StepKind::Parallel);
    assert!(matches!(p.output_queues[..], [QueueSpec::ByteBounded { limit_bytes: 4096 }]));
    assert!(matches!(p.branch_ordering[..], [BranchOrdering::ByItemOrdinal]));
    assert_eq!(extract.pool_placement(), PoolPlacement::ExcludeReader);
}

/// A full batch of 24-byte-key records commits 4096 × 40 bytes of refs.
#[test]
fn a_full_k24_batch_commits_its_ref_bytes() {
    let mut batch = key_batch(0);
    batch.n_records = u32::try_from(KEY_BATCH_RECORDS).unwrap();
    assert_eq!(TemplateKeys::ref_width(batch.ctx.as_ref().unwrap()), 40);
    assert_eq!(batch.heap_size(), KEY_BATCH_RECORDS * 40);
}

/// The census classifies each batch against its own run's scan end: a batch
/// keyed before the scan cut the run's last batch overlapped the scan, one
/// keyed after it is tail. Run 0 (5 records) is keyed before its scan end,
/// run 1 (3 records) after.
#[test]
fn census_classifies_batches_against_their_runs_scan_end() {
    let mut extract = ExtractTemplateKeys::new(1 << 20);
    let mut ex = StepProbe::new(&extract);
    ex.push_input(key_batch_with(0, 5, 0));
    ex.push_input(key_batch_with(1, 3, 1));
    for _ in 0..2 {
        assert_eq!(ex.try_run(&mut extract).unwrap(), StepOutcome::Progress);
    }
    let out = ex.take_output::<TemplateRefBatch>(0);
    let (mut a, mut b) = (out.pop().expect("run 0"), out.pop().expect("run 1"));
    let t0 = Instant::now();
    let end = |ingested, scan_end| RunEnd {
        seals_to_spill: true,
        records_ingested_so_far: ingested,
        run_bytes: 0,
        scan_end,
    };
    (a.keyed_at, a.run_end) = (t0, Some(end(5, t0 + Duration::from_secs(1))));
    (b.keyed_at, b.run_end) =
        (t0 + Duration::from_secs(2), Some(end(8, t0 + Duration::from_secs(1))));
    let mut sort = SortTemplateRun::new(1, 1 << 20);
    let sp = StepProbe::new(&sort);
    sp.push_input(a);
    sp.push_input(b);
    for _ in 0..4 {
        sp.try_run(&mut sort).unwrap();
    }
    let c = sort.census();
    assert_eq!((c.runs, c.records, c.overlapped_records, c.tail_records), (2, 8, 5, 3));
}

/// The record path (`QuerynameStrategy::push_record`, batch-closed every 4096
/// records) and the split path (`QuerynameKeys::extract_batch`, one call per
/// 4096-record batch, one builder carried across them) produce identical refs
/// — including the batch that re-learns its prefix class when the instrument
/// changes, and a short final batch.
#[test]
fn extract_queryname_keys_shares_the_one_extraction() {
    use fgumi_sort::{QN_CLASS_BATCH_RECORDS, QnLex, QnRefBuilder, QnSortState};
    let n = 2 * QN_CLASS_BATCH_RECORDS + 100;
    let bodies: Vec<Vec<u8>> = (0..n)
        .map(|i| {
            let name = if i < QN_CLASS_BATCH_RECORDS {
                format!("A00132:53:HFHJ:1:1101:{i}:7")
            } else {
                format!("M04567:12:000000000-ABCDE:1:1101:{i}:9")
            };
            make_bam_bytes(-1, -1, 0x4, name.as_bytes(), &[], 0, -1, -1, &[])
        })
        .collect();
    let (arena, at) = arena_with(16, &frames(&bodies));

    let mut strategy = QuerynameStrategy::<QnLex>::new(QnSortState::new(false));
    let mut off = at;
    for b in &bodies {
        let len = u32::try_from(b.len()).unwrap();
        strategy.push_record(b, off + 4, len).unwrap();
        off += 4 + u64::from(len);
    }
    let want = strategy.refs_for_test(&arena);

    let ctx = QuerynameKeys::<QnLex>::provision(&QnSortState::new(false), &bodies[0]);
    let mut builder = QnRefBuilder::default();
    let mut got = Vec::with_capacity(n);
    let mut start = at;
    for (i, batch) in bodies.chunks(QN_CLASS_BATCH_RECORDS).enumerate() {
        let bytes: u64 = batch.iter().map(|b| 4 + b.len() as u64).sum();
        let (refs, violation) = QuerynameKeys::<QnLex>::extract_batch(
            &ctx,
            &mut builder,
            &arena,
            start,
            start + bytes,
            u32::try_from(batch.len()).unwrap(),
            (i * QN_CLASS_BATCH_RECORDS) as u64,
        )
        .unwrap();
        assert!(violation.is_none(), "queryname keys have no lanes to violate");
        got.extend(refs);
        start += bytes;
    }
    // Every field of every ref (QnKey has no `PartialEq`: the comparator is the
    // only order on it).
    let fields = |refs: &[fgumi_sort::QnRef]| -> Vec<_> {
        refs.iter()
            .map(|(k, off, len)| (k.window, k.flags, k.tail_len, k.class, *off, *len))
            .collect()
    };
    assert_eq!(got.len(), want.len());
    assert!(fields(&got) == fields(&want), "the split path's refs differ from the record path's");
    assert!(want.iter().any(|r| r.0.class != 0), "the fixture learns a prefix class");
}

/// A queryname batch whose frames do not match the scan's cut is an error,
/// not a short batch.
#[test]
fn queryname_extract_rejects_mismatched_frames() {
    use fgumi_sort::{QnLex, QnRefBuilder, QnSortState};
    let bodies: Vec<Vec<u8>> = (0..3).map(body).collect();
    let data = frames(&bodies);
    let (arena, at) = arena_with(16, &data);
    let ctx = QuerynameKeys::<QnLex>::provision(&QnSortState::new(false), &bodies[0]);
    let end = at + data.len() as u64;
    for (n, end) in [(4, end), (2, end), (3, end - 1)] {
        let err = QuerynameKeys::<QnLex>::extract_batch(
            &ctx,
            &mut QnRefBuilder::default(),
            &arena,
            at,
            end,
            n,
            0,
        )
        .expect_err("mismatched frames");
        assert_eq!(err.kind(), io::ErrorKind::InvalidData, "n={n}: {err}");
    }
}

/// The queryname split steps carry their own names, and the run sort's pool
/// is the queryname pool.
#[test]
fn queryname_split_step_names() {
    use fgumi_sort::{QnLex, QnNatural};
    for names in [QuerynameKeys::<QnLex>::STEP_NAMES, QuerynameKeys::<QnNatural>::STEP_NAMES] {
        assert_eq!(names, ["FindQuerynameBoundaries", "ExtractQuerynameKeys", "SortQuerynameRun"]);
    }
    assert_eq!(QuerynameKeys::<QnLex>::pool_name(), "qname-sort");
}

/// A coordinate record body: `tid`/`pos` (`-1` for unmapped), reverse strand on
/// odd `i`, named `c{i}`.
fn coord_body(i: usize, tid: i32, pos: i32) -> Vec<u8> {
    let flags = if i % 2 == 1 { 0x10 } else { 0 } | if tid < 0 { 0x4 } else { 0 };
    make_bam_bytes(tid, pos, flags, format!("c{i}").as_bytes(), &[], 20, -1, -1, &[])
}

/// A coordinate run extracted in 64-record batches and sealed is the chunk the
/// copy-based `CoordinateChunkSorter` builds over the same records — ties
/// (equal coordinates) in input order, unmapped records last. (Twin of the
/// deleted `find_boundaries_and_sort_single_run_matches_oracle`.)
#[rstest::rstest]
fn coordinate_run_matches_the_copy_based_sorter(#[values(1usize, 3)] sort_threads: usize) {
    use fgumi_sort::CoordinateChunkSorter;
    let n_ref = 4u32;
    let bodies: Vec<Vec<u8>> = (0..300usize)
        .map(|i| {
            let tid = if i % 11 == 0 { -1 } else { i32::try_from(i * 7 % 4).unwrap() };
            let pos = if tid < 0 { -1 } else { i32::try_from(i * 31 % 50).unwrap() };
            coord_body(i, tid, pos)
        })
        .collect();
    let (arena, at) = arena_with(16, &frames(&bodies));
    let ctx = CoordinateKeys::provision(&n_ref, &bodies[0]);
    let mut acc = CoordinateKeys::accum_for(&ctx);
    let mut start = at;
    for (b, batch) in bodies.chunks(64).enumerate() {
        let bytes: u64 = batch.iter().map(|r| 4 + r.len() as u64).sum();
        let (refs, violation) = CoordinateKeys::extract_batch(
            &ctx,
            &mut (),
            &arena,
            start,
            start + bytes,
            u32::try_from(batch.len()).unwrap(),
            b as u64 * 64,
        )
        .unwrap();
        assert!(violation.is_none());
        assert_eq!(CoordinateKeys::refs_byte_len(&refs), batch.len() * 24);
        CoordinateKeys::append(&mut acc, refs);
        start += bytes;
    }
    let pool = BoundedSortPool::new(CoordinateKeys::pool_name());
    let MemoryChunkErased::Coordinate(chunk) =
        CoordinateKeys::seal(&mut acc, &ctx, arena, sort_threads, &pool)
    else {
        panic!("a coordinate run seals a coordinate chunk");
    };
    assert!(acc.is_empty(), "the seal leaves the accumulator empty for the next run");
    let mut oracle = CoordinateChunkSorter::for_test(usize::MAX, n_ref);
    for r in &bodies {
        let _ = oracle.push(r).unwrap();
    }
    let oc = oracle.take_sorted_chunk();
    let got: Vec<&[u8]> = (0..chunk.len()).map(|i| chunk.record_bytes(i)).collect();
    let want: Vec<&[u8]> = (0..oc.len()).map(|i| oc.record_bytes(i)).collect();
    assert!(got == want, "t{sort_threads}: the sealed run differs from the copy-based sorter");
}

/// A run whose input drains before `ReadBlocks` marked its last block is sealed
/// by `finalize` as the residual: one terminal batch carrying a non-spilling
/// `RunEnd`, and a second `finalize` cuts nothing. (Twin of the deleted
/// `finalize_stages_residual_before_all_announced_each_once`, boundary half.)
#[test]
fn boundaries_finalize_cuts_the_open_run_as_the_residual() {
    let mut step = FindCoordinateBoundaries::new(1, 1 << 20).with_batch_records(64);
    let bodies: Vec<Vec<u8>> = (0..10).map(|i| coord_body(i, 0, 10)).collect();
    let data = [bam_header(), frames(&bodies)].concat();
    let (a, at) = arena_with(0, &data);
    let mut open = block(&a, 0, at, data.len(), 0, false);
    open.is_last_of_run = false;
    assert_eq!(step.ingest_block(&open).unwrap(), 10);
    assert!(step.pending.iter().all(|b| b.run_end.is_none()), "no run end before finalize");
    step.finalize().unwrap();
    let ends: Vec<RunEnd> = step.pending.iter().filter_map(|b| b.run_end).collect();
    assert_eq!(ends.len(), 1);
    assert!(!ends[0].seals_to_spill && ends[0].records_ingested_so_far == 10);
    let n = step.pending.len();
    step.finalize().unwrap();
    assert_eq!(step.pending.len(), n, "a second finalize cuts nothing");
}

/// `SortRun` emits a residual run as `Residual` then `AllAnnounced`, each once
/// and in that order, then finishes. (Twin of the deleted
/// `finalize_stages_residual_before_all_announced_each_once`, sort half.)
#[test]
fn sort_run_stages_residual_before_all_announced_each_once() {
    let bodies: Vec<Vec<u8>> =
        (0..5).map(|i| coord_body(i, 1, 50 - i32::try_from(i).unwrap())).collect();
    let data = frames(&bodies);
    let (arena, at) = arena_with(16, &data);
    let ctx = CoordinateKeys::provision(&2, &bodies[0]);
    let (refs, _) =
        CoordinateKeys::extract_batch(&ctx, &mut (), &arena, at, at + data.len() as u64, 5, 0)
            .unwrap();
    let now = Instant::now();
    let batch = RefBatch::<CoordinateKeys> {
        arena,
        ctx: Some(ctx),
        ordinal: 0,
        run_seq: 0,
        first_record_index: 0,
        refs: Some(refs),
        violation: None,
        keyed_at: now,
        run_end: Some(RunEnd {
            seals_to_spill: false,
            records_ingested_so_far: 5,
            run_bytes: data.len() as u64,
            scan_end: now,
        }),
    };
    let mut sort = SortCoordinateRun::new(1, 1 << 20);
    let mut sp = StepProbe::new(&sort);
    sp.push_input(batch);
    sp.close_input();
    let mut last = StepOutcome::Progress;
    for _ in 0..8 {
        last = sp.try_run(&mut sort).unwrap();
        if last == StepOutcome::Finished {
            break;
        }
    }
    assert_eq!(last, StepOutcome::Finished);
    let out = sp.take_output::<SortChunkEvent>(0);
    let events: Vec<SortChunkEvent> = std::iter::from_fn(|| out.pop()).collect();
    assert_eq!(events.len(), 2, "Residual and AllAnnounced, once each");
    assert!(matches!(&events[0], SortChunkEvent::Residual { chunk, .. } if chunk.len() == 5));
    assert!(matches!(
        events[1],
        SortChunkEvent::AllAnnounced { slot_count: 0, memory_chunk_count: 1, total_records: 5 }
    ));
}

/// The run sort sizes the order's bounded pool to the phase-1 width it was
/// built with (`--sort-threads`), observed on the queryname pool, which every
/// queryname seal installs. (Twin of the deleted
/// `phase1_sort_threads_sizes_the_queryname_worker_pool`.)
#[rstest::rstest]
fn sort_run_sizes_the_queryname_pool_to_sort_threads(#[values(1usize, 3)] sort_threads: usize) {
    use fgumi_sort::{QnLex, QnRefBuilder, QnSortState};
    let bodies: Vec<Vec<u8>> = (0..4).map(body).collect();
    let data = frames(&bodies);
    let (arena, at) = arena_with(16, &data);
    let ctx = QuerynameKeys::<QnLex>::provision(&QnSortState::new(false), &bodies[0]);
    let (refs, _) = QuerynameKeys::<QnLex>::extract_batch(
        &ctx,
        &mut QnRefBuilder::default(),
        &arena,
        at,
        at + data.len() as u64,
        4,
        0,
    )
    .unwrap();
    let now = Instant::now();
    let batch = RefBatch::<QuerynameKeys<QnLex>> {
        arena,
        ctx: Some(ctx),
        ordinal: 0,
        run_seq: 0,
        first_record_index: 0,
        refs: Some(refs),
        violation: None,
        keyed_at: now,
        run_end: Some(RunEnd {
            seals_to_spill: false,
            records_ingested_so_far: 4,
            run_bytes: data.len() as u64,
            scan_end: now,
        }),
    };
    let mut sort = SortQuerynameRun::<QnLex>::new(sort_threads, 1 << 20);
    assert_eq!(sort.sort_pool_threads(), None, "the pool is built on the first seal");
    let sp = StepProbe::new(&sort);
    sp.push_input(batch);
    sp.try_run(&mut sort).unwrap();
    assert_eq!(sort.sort_pool_threads(), Some(sort_threads));
}

/// The coordinate split steps carry their own names and pool name.
#[test]
fn coordinate_split_step_names() {
    assert_eq!(
        CoordinateKeys::STEP_NAMES,
        ["FindCoordinateBoundaries", "ExtractCoordinateKeys", "SortCoordinateRun"]
    );
    assert_eq!(CoordinateKeys::pool_name(), "coord-sort");
    assert_eq!(CoordinateKeys::ref_width(&CoordinateKeys::provision(&1, &[])), 24);
}

/// A coordinate run sorts in parallel — and so reserves the whole phase-1 cap —
/// exactly when the radix would spawn more than one thread for it.
#[test]
fn coordinate_sorts_in_parallel_exactly_when_the_radix_does() {
    let mut saw_parallel = false;
    for len in [0usize, 1, 1 << 10, 1 << 16, 1 << 20] {
        let acc = vec![fgumi_sort::RecordRef::new(0, 0, 0); len];
        for threads in [1usize, 2, 8] {
            let want = fgumi_sort::radix_sorts_in_parallel(len, threads);
            assert_eq!(CoordinateKeys::sorts_in_parallel(&acc, threads), want, "{len} t{threads}");
            saw_parallel |= want;
        }
    }
    assert!(saw_parallel, "some run is large enough to sort in parallel");
}

/// A run accumulator whose capacity a test can read.
trait AccCapacity {
    fn cap(&self) -> usize;
}

impl AccCapacity for Vec<fgumi_sort::RecordRef> {
    fn cap(&self) -> usize {
        self.capacity()
    }
}

impl AccCapacity for fgumi_sort::TemplateArenaRefs {
    fn cap(&self) -> usize {
        self.capacity()
    }
}

impl AccCapacity for Vec<fgumi_sort::QnRef> {
    fn cap(&self) -> usize {
        self.capacity()
    }
}

/// One run's single ref batch over `bodies`, extracted through `X`.
fn run_batch<X: ArenaKeyExtractor>(
    ctx: &X::Ctx,
    bodies: &[Vec<u8>],
    ordinal: u64,
    run_seq: u32,
    seals_to_spill: bool,
) -> RefBatch<X> {
    let data = frames(bodies);
    let (arena, at) = arena_with(16, &data);
    let n = u32::try_from(bodies.len()).unwrap();
    let (refs, violation) =
        X::extract_batch(ctx, &mut X::Scratch::default(), &arena, at, at + data.len() as u64, n, 0)
            .unwrap();
    let now = Instant::now();
    RefBatch {
        arena,
        ctx: Some(ctx.clone()),
        ordinal,
        run_seq,
        first_record_index: 0,
        refs: Some(refs),
        violation,
        keyed_at: now,
        run_end: Some(RunEnd {
            seals_to_spill,
            records_ingested_so_far: u64::from(n),
            run_bytes: data.len() as u64,
            scan_end: now,
        }),
    }
}

/// Feed `batches` to `sort` one at a time, returning each sealed chunk's
/// record bodies and the accumulator's capacity right after each seal.
fn seal_runs<X: ArenaKeyExtractor>(
    sort: &mut SortRun<X>,
    batches: Vec<RefBatch<X>>,
) -> Vec<(Vec<Vec<u8>>, usize)>
where
    X::Accum: AccCapacity,
{
    let mut sp = StepProbe::new(&*sort);
    let out = sp.take_output::<SortChunkEvent>(0);
    let mut sealed = Vec::new();
    for batch in batches {
        sp.push_input(batch);
        sp.try_run(sort).unwrap();
        let cap = sort.acc.as_ref().map_or(0, AccCapacity::cap);
        while let Some(ev) = out.pop() {
            if let SortChunkEvent::Spill { chunk, .. } | SortChunkEvent::Residual { chunk, .. } = ev
            {
                sealed.push((
                    (0..chunk.len()).map(|i| chunk.record_bytes(i).to_vec()).collect(),
                    cap,
                ));
            }
        }
    }
    sealed
}

/// Every order's seal moves the run's refs into its chunk, so the run sort
/// keeps no ref buffer beside a chunk that is spilling (or, after the final
/// run, through the merge): right after each seal — the spilled run 1 and the
/// residual run 2 — the accumulator holds no capacity, and run 2's chunk is
/// exactly the chunk a fresh run sort seals for the same records.
#[rstest::rstest]
fn seals_keep_no_ref_buffer(#[values("template", "coordinate", "queryname")] order: &str) {
    fn check<X: ArenaKeyExtractor>(ctx: &X::Ctx, run1: &[Vec<u8>], run2: &[Vec<u8>])
    where
        X::Accum: AccCapacity,
    {
        let mut sort = SortRun::<X>::new(1, 1 << 20);
        let sealed = seal_runs(
            &mut sort,
            vec![run_batch::<X>(ctx, run1, 0, 0, true), run_batch::<X>(ctx, run2, 1, 1, false)],
        );
        assert_eq!(sealed.len(), 2);
        assert_eq!(sealed[0].1, 0, "the spilled run's refs moved into its chunk");
        assert_eq!(sealed[1].1, 0, "the final run's refs moved into its chunk");
        let mut fresh = SortRun::<X>::new(1, 1 << 20);
        let alone = seal_runs(&mut fresh, vec![run_batch::<X>(ctx, run2, 0, 0, false)]);
        assert!(sealed[1].0 == alone[0].0, "run 2 seals as on a fresh run sort");
    }
    let run1: Vec<Vec<u8>> = (0..1500).map(body).collect();
    let run2: Vec<Vec<u8>> = (1500..2300).map(body).collect();
    match order {
        "template" => check::<TemplateKeys>(&seed().provision(&run1[0]), &run1, &run2),
        "queryname" => {
            let ctx = QuerynameKeys::<fgumi_sort::QnLex>::provision(
                &fgumi_sort::QnSortState::new(false),
                &run1[0],
            );
            check::<QuerynameKeys<fgumi_sort::QnLex>>(&ctx, &run1, &run2);
        }
        _ => {
            let run1: Vec<Vec<u8>> =
                (0..1500).map(|i| coord_body(i, 1, i32::try_from(i * 37 % 900).unwrap())).collect();
            let run2: Vec<Vec<u8>> =
                (0..800).map(|i| coord_body(i, 0, i32::try_from(i * 11 % 300).unwrap())).collect();
            check::<CoordinateKeys>(&CoordinateKeys::provision(&2, &run1[0]), &run1, &run2);
        }
    }
}

/// A coordinate or template seal converts the refs into the chunk's tuples in
/// place: the chunk's records live where the accumulator's refs did.
#[rstest::rstest]
fn coordinate_and_template_seals_move_the_refs(#[values(true, false)] template: bool) {
    fn check<X: ArenaKeyExtractor>(ctx: &X::Ctx, bodies: &[Vec<u8>], at: fn(&X::Accum) -> *const u8)
    where
        X::Accum: AccCapacity,
    {
        let batch = run_batch::<X>(ctx, bodies, 0, 0, false);
        let mut acc = X::accum_for(ctx);
        X::append(&mut acc, batch.refs.unwrap());
        let ptr = at(&acc);
        let chunk = X::seal(&mut acc, ctx, batch.arena, 1, &BoundedSortPool::new("test-sort"));
        let records = match &chunk {
            MemoryChunkErased::Coordinate(c) => c.records_ptr().cast::<u8>(),
            MemoryChunkErased::TemplateCoordinate(fgumi_sort::TemplateMemChunk::K24(c)) => {
                c.records_ptr().cast::<u8>()
            }
            _ => panic!("an unexpected chunk"),
        };
        assert_eq!(chunk.len(), bodies.len());
        assert!(std::ptr::eq(records, ptr), "the chunk's tuples are the refs' allocation");
        assert_eq!(acc.cap(), 0, "the accumulator gave its buffer to the chunk");
    }
    if template {
        let bodies: Vec<Vec<u8>> = (0..500).map(body).collect();
        check::<TemplateKeys>(&seed().provision(&bodies[0]), &bodies, |acc| match acc {
            fgumi_sort::TemplateArenaRefs::K24(v) => v.as_ptr().cast(),
            _ => panic!("the fixture provisions the 24-byte lane"),
        });
    } else {
        let bodies: Vec<Vec<u8>> =
            (0..500).map(|i| coord_body(i, 1, i32::try_from(i * 37 % 900).unwrap())).collect();
        check::<CoordinateKeys>(&CoordinateKeys::provision(&2, &bodies[0]), &bodies, |acc| {
            acc.as_ptr().cast()
        });
    }
}

/// A queryname seal moves the run's refs into the chunk: the chunk's record
/// tuples live where the accumulator's refs did. (Its accumulator therefore
/// starts every run empty.)
#[test]
fn queryname_seal_moves_the_refs_without_copying() {
    use fgumi_sort::{QnLex, QnSortState};
    let bodies: Vec<Vec<u8>> = (0..500).map(body).collect();
    let ctx = QuerynameKeys::<QnLex>::provision(&QnSortState::new(false), &bodies[0]);
    let batch = run_batch::<QuerynameKeys<QnLex>>(&ctx, &bodies, 0, 0, false);
    let mut acc = QuerynameKeys::<QnLex>::accum_for(&ctx);
    QuerynameKeys::<QnLex>::append(&mut acc, batch.refs.unwrap());
    let at = acc.as_ptr();
    let pool = BoundedSortPool::new("qname-sort");
    let MemoryChunkErased::QuerynameLex(chunk) =
        QuerynameKeys::<QnLex>::seal(&mut acc, &ctx, batch.arena, 2, &pool)
    else {
        panic!("a lexicographic queryname chunk");
    };
    assert_eq!(chunk.len(), 500);
    assert!(std::ptr::eq(chunk.records_ptr().cast::<u8>(), at.cast::<u8>()), "refs moved");
    assert_eq!(acc.capacity(), 0, "the accumulator gave its buffer to the chunk");
}
