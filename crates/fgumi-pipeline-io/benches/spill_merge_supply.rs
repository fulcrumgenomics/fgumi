#![deny(unsafe_code)]
//! The merge supply in isolation: `SpillReadPlanner → PreadSpillSlices →
//! SortSpillDecompress` over pre-written zstd spills, drained by a synthetic
//! consumer that pops the slots' FIFOs in a loser-tree-like interleaving
//! (runs of consecutive blocks per source with p50 = 1 and p99 = 1024) and
//! publishes the awaited and predicted slots as the merge does.
//!
//! Per `k ∈ {27, 256, 1024}` it reports blocks/s (criterion), and once per
//! case on stderr: the consumer's wait for a block it stalled on (p50/p99),
//! the self-serve share of claims, the stash peak, and allocations per block
//! (a separate pass under a counting allocator, so it does not distort the
//! timing). A second group, `pread_slices/{1,4,8}-streams`, measures
//! `PreadSlices` alone reading a 1 GiB file in 4 MiB fills.
//!
//! Spills go to a tempdir, or to `FGUMI_BENCH_DIR` (a real disk);
//! `FGUMI_BENCH_PREAD_BYTES` overrides the `pread_slices` file size.
//!
//!     cargo bench -p fgumi-pipeline-io --bench spill_merge_supply
use std::collections::VecDeque;
use std::io::{self, Write};
use std::path::{Path, PathBuf};
use std::sync::Arc;
use std::sync::atomic::{AtomicU64, Ordering};
use std::time::{Duration, Instant};

use criterion::{Criterion, Throughput, criterion_group, criterion_main};
use fgumi_bam_io::pread::{PositionalSource, ReadStreamsPolicy, SliceBufferPool, slice_ranges};
use fgumi_pipeline_core::outputs::Single;
use fgumi_pipeline_core::queues::QueueSpec;
use fgumi_pipeline_core::reorder::BranchOrdering;
use fgumi_pipeline_core::step::{Step, StepCtx, StepKind, StepOutcome, StepProfile};
use fgumi_pipeline_core::{HeldSlot, Pipeline, PipelineConfig, Unpushed};
use fgumi_pipeline_io::pread::{PreadSlices, ReadRequest, ReadSlice, ReadTarget};
use fgumi_pipeline_io::sort::protocol::{SortPhase1Event, SortPhase2Event};
use fgumi_pipeline_io::sort::{SortSpillDecompress, SpillReadPlanner, SpillSupply, SupplyLedger};
use fgumi_sort::{MergeDemand, SortMergeSlot, SpillCodec};

#[global_allocator]
static ALLOC: dhat::Alloc = dhat::Alloc;

/// Blocks across all spills of one case (split evenly over `k`).
const TOTAL_BLOCKS: usize = 16_384;
/// Records per block (~8 KiB of keyed records per block).
const RECORDS_PER_BLOCK: usize = 64;
const THREADS: usize = 8;
const LIMIT: u64 = 4 << 20;

fn bench_dir() -> tempfile::TempDir {
    match std::env::var_os("FGUMI_BENCH_DIR") {
        Some(d) => tempfile::tempdir_in(d).expect("FGUMI_BENCH_DIR tempdir"),
        None => tempfile::tempdir().expect("tempdir"),
    }
}

/// splitmix64: a seeded, dependency-free generator.
fn mix(mut x: u64) -> u64 {
    x = x.wrapping_add(0x9E37_79B9_7F4A_7C15);
    x = (x ^ (x >> 30)).wrapping_mul(0xBF58_476D_1CE4_E5B9);
    x = (x ^ (x >> 27)).wrapping_mul(0x94D0_49BB_1331_11EB);
    x ^ (x >> 31)
}

/// Write `k` zstd spills of `blocks` blocks each; returns their paths.
fn write_spills(dir: &Path, k: usize, blocks: usize) -> Vec<PathBuf> {
    (0..k)
        .map(|f| {
            let path = dir.join(format!("bench{f}.spill"));
            let mut file = std::io::BufWriter::new(std::fs::File::create(&path).unwrap());
            file.write_all(fgumi_sort::spill_magic(SpillCodec::Zstd)).unwrap();
            let mut c = fgumi_sort::SpillBlockCompressor::new(SpillCodec::Zstd, 1).unwrap();
            let mut raw = Vec::new();
            for b in 0..blocks {
                raw.clear();
                for i in 0..RECORDS_PER_BLOCK {
                    let n = ((b * RECORDS_PER_BLOCK + i) * k + f) as u64;
                    let key = fgumi_sort::RawCoordinateKey { sort_key: n };
                    let body: Vec<u8> = (0..112).map(|j| (mix(n ^ j) & 0x3) as u8 + b'A').collect();
                    fgumi_sort::frame_keyed_record_into(&mut raw, &key, &body).unwrap();
                }
                file.write_all(&c.compress_block(&raw).unwrap()).unwrap();
            }
            file.write_all(fgumi_sort::spill_trailer(SpillCodec::Zstd)).unwrap();
            drop(file);
            path
        })
        .collect()
}

/// Open the spills as fresh slots (a slot carries one run's state).
fn open_slots(paths: &[PathBuf]) -> Vec<Arc<SortMergeSlot>> {
    paths
        .iter()
        .enumerate()
        .map(|(f, p)| fgumi_sort::open_spill_slot(p, u32::try_from(f).unwrap()).unwrap())
        .collect()
}

/// `Exclusive` source of the phase events.
struct Events {
    events: VecDeque<SortPhase1Event>,
    held: HeldSlot<Unpushed<SortPhase1Event>>,
}

impl Step for Events {
    type Input = ();
    type Outputs = Single<SortPhase1Event>;

    fn profile(&self) -> StepProfile {
        StepProfile {
            name: "Events",
            kind: StepKind::Exclusive,
            sticky: true,
            output_queues: vec![QueueSpec::ByteBounded { limit_bytes: LIMIT }],
            branch_ordering: vec![BranchOrdering::None],
        }
    }

    fn try_run(&mut self, ctx: &mut StepCtx<'_, Self>) -> io::Result<StepOutcome> {
        if let Some(u) = self.held.take()
            && let Err(again) = ctx.outputs.retry(u)
        {
            self.held.put(again);
            return Ok(StepOutcome::Contention);
        }
        let Some(e) = self.events.pop_front() else { return Ok(StepOutcome::Finished) };
        if let Err(u) = ctx.outputs.push(e) {
            self.held.put(u);
        }
        Ok(StepOutcome::Progress)
    }
}

/// What the consumer measured.
#[derive(Default)]
struct ConsumerReport {
    blocks: u64,
    waits: Vec<Duration>,
}

/// `Detached` consumer: registers the slots, then pops runs of consecutive
/// blocks per source, publishing awaited (on an empty FIFO) and predicted (the
/// next source) through the merge demand.
struct Consumer {
    demand: Arc<MergeDemand>,
    slots: Vec<Arc<SortMergeSlot>>,
    announced: bool,
    live: Vec<usize>,
    current: usize,
    run_left: u64,
    rng: u64,
    waiting_since: Option<Instant>,
    report: Arc<parking_lot::Mutex<ConsumerReport>>,
}

impl Consumer {
    /// Run length: p50 = 1, p99 = 1024 consecutive blocks (log-uniform tail).
    fn run_length(&mut self) -> u64 {
        self.rng = mix(self.rng);
        let u = self.rng % 1000;
        if u < 500 { 1 } else { 1 << (1 + (u - 500) * 10 / 495).min(10) }
    }

    fn next_source(&mut self) {
        self.rng = mix(self.rng);
        self.current = self.live[usize::try_from(self.rng % self.live.len() as u64).unwrap()];
        self.run_left = self.run_length();
        self.rng = mix(self.rng);
        let predicted = self.live[usize::try_from(self.rng % self.live.len() as u64).unwrap()];
        self.demand.set_predicted(Some(self.slots[predicted].file_id));
    }
}

impl Step for Consumer {
    type Input = SortPhase2Event;
    type Outputs = ();

    fn profile(&self) -> StepProfile {
        StepProfile {
            name: "Consumer",
            kind: StepKind::Detached,
            sticky: true,
            output_queues: vec![],
            branch_ordering: vec![],
        }
    }

    fn try_run(&mut self, ctx: &mut StepCtx<'_, Self>) -> io::Result<StepOutcome> {
        while let Some(e) = ctx.input.pop() {
            match e {
                SortPhase2Event::SpillReady { slot, .. } => self.slots.push(slot),
                SortPhase2Event::AllAnnounced { .. } => {
                    self.announced = true;
                    self.live = (0..self.slots.len()).collect();
                    self.next_source();
                }
                SortPhase2Event::MemoryChunk { .. } => {}
            }
        }
        if !self.announced {
            return Ok(if ctx.input.is_drained() {
                StepOutcome::Finished
            } else {
                StepOutcome::NoProgress
            });
        }
        let mut progressed = false;
        for _ in 0..256 {
            if self.live.is_empty() {
                return Ok(StepOutcome::Finished);
            }
            let slot = Arc::clone(&self.slots[self.current]);
            if let Some(_block) = slot.pop_decompressed() {
                if let Some(t0) = self.waiting_since.take() {
                    self.report.lock().waits.push(t0.elapsed());
                }
                self.report.lock().blocks += 1;
                progressed = true;
                self.run_left -= 1;
                if self.run_left == 0 {
                    self.next_source();
                }
                continue;
            }
            if slot.queue_eof() {
                assert!(!slot.has_error(), "a spill failed");
                self.live.retain(|&i| i != self.current);
                if !self.live.is_empty() {
                    self.next_source();
                }
                continue;
            }
            if self.demand.await_slot(&slot) {
                continue;
            }
            self.waiting_since.get_or_insert_with(Instant::now);
            break;
        }
        Ok(if progressed { StepOutcome::Progress } else { StepOutcome::NoProgress })
    }
}

/// What one run measured.
struct SupplyRun {
    elapsed: Duration,
    blocks: u64,
    waits: Vec<Duration>,
    ledger: Arc<SupplyLedger>,
}

fn run_supply(slots: &[Arc<SortMergeSlot>]) -> SupplyRun {
    let demand = Arc::new(MergeDemand::new());
    let ledger = Arc::new(SupplyLedger::default());
    let mut events: VecDeque<_> = slots
        .iter()
        .map(|s| SortPhase1Event::SpillReady {
            slot: Arc::clone(s),
            path: PathBuf::new(),
            records_ingested_so_far: 0,
        })
        .collect();
    events.push_back(SortPhase1Event::AllAnnounced {
        slot_count: u32::try_from(slots.len()).unwrap(),
        memory_chunk_count: 0,
        total_records: 0,
    });
    let supply = SpillSupply { demand: Arc::clone(&demand), ledger: Arc::clone(&ledger) };
    let slices = SliceBufferPool::new(2 * THREADS + 8);
    let planner =
        SpillReadPlanner::new(12 << 30, THREADS, THREADS, LIMIT, &supply, Arc::clone(&slices))
            .with_read_streams(ReadStreamsPolicy::fixed(1));
    let pread = PreadSlices::spill(slices, LIMIT);
    let decompress = SortSpillDecompress::new(LIMIT, &supply);
    let report = Arc::new(parking_lot::Mutex::new(ConsumerReport::default()));
    let consumer = Consumer {
        demand,
        slots: Vec::new(),
        announced: false,
        live: Vec::new(),
        current: 0,
        run_left: 0,
        rng: 0x5EED,
        waiting_since: None,
        report: Arc::clone(&report),
    };
    let builder = Pipeline::builder();
    let fan = builder.chain(Events { events, held: HeldSlot::new() }).chain(planner).into_multi();
    fan.b0.chain(pread).chain(decompress).into_sink_marker();
    fan.b1.chain(consumer).into_sink_marker();
    let pipeline = builder.build().unwrap();
    let t0 = Instant::now();
    pipeline.run(PipelineConfig { threads: THREADS, ..Default::default() }).unwrap();
    let elapsed = t0.elapsed();
    let r = std::mem::take(&mut *report.lock());
    SupplyRun { elapsed, blocks: r.blocks, waits: r.waits, ledger }
}

fn percentile(sorted: &[Duration], p: f64) -> Duration {
    if sorted.is_empty() {
        return Duration::ZERO;
    }
    #[allow(clippy::cast_possible_truncation, clippy::cast_sign_loss, clippy::cast_precision_loss)]
    let i = ((sorted.len() - 1) as f64 * p).round() as usize;
    sorted[i]
}

fn bench_spill_merge_supply(c: &mut Criterion) {
    let mut group = c.benchmark_group("spill_merge_supply");
    group.sample_size(10);
    for k in [27usize, 256, 1024] {
        let dir = bench_dir();
        let paths = write_spills(dir.path(), k, TOTAL_BLOCKS / k);
        let blocks = ((TOTAL_BLOCKS / k) * k) as u64;
        group.throughput(Throughput::Elements(blocks));
        let mut last = None;
        group.bench_function(format!("k{k}"), |b| {
            b.iter_custom(|iters| {
                let mut total = Duration::ZERO;
                for _ in 0..iters {
                    let slots = open_slots(&paths);
                    let run = run_supply(&slots);
                    assert_eq!(run.blocks, blocks, "every block consumed");
                    total += run.elapsed;
                    last = Some(run);
                }
                total
            });
        });
        if let Some(mut run) = last {
            run.waits.sort_unstable();
            let (workers, consumer) = run.ledger.claims();
            let slots = open_slots(&paths);
            let profiler = dhat::Profiler::builder().testing().build();
            let before = dhat::HeapStats::get().total_blocks;
            let counted = run_supply(&slots);
            let allocs = dhat::HeapStats::get().total_blocks - before;
            drop(profiler);
            #[allow(clippy::cast_precision_loss)]
            let share = consumer as f64 / (workers + consumer).max(1) as f64;
            #[allow(clippy::cast_precision_loss)]
            let per_block = allocs as f64 / counted.blocks.max(1) as f64;
            eprintln!(
                "spill_merge_supply/k{k}: consumer wait p50 {:?} p99 {:?} over {} waits; \
                 self-serve share {:.1}%; stash peak {} KiB; {per_block:.1} allocations/block",
                percentile(&run.waits, 0.50),
                percentile(&run.waits, 0.99),
                run.waits.len(),
                share * 100.0,
                run.ledger.stash_peak_bytes() >> 10,
            );
        }
    }
    group.finish();
}

/// Serial source of `ReadRequest`s tiling a file in 4 MiB fills of `streams`
/// slices each.
struct Requests {
    source: Arc<dyn PositionalSource>,
    len: u64,
    next: u64,
    ordinal: u64,
    streams: usize,
    pending: VecDeque<ReadRequest>,
    held: HeldSlot<Unpushed<ReadRequest>>,
}

impl Step for Requests {
    type Input = ();
    type Outputs = Single<ReadRequest>;

    fn profile(&self) -> StepProfile {
        StepProfile {
            name: "Requests",
            kind: StepKind::Exclusive,
            sticky: true,
            output_queues: vec![QueueSpec::ByteBounded { limit_bytes: 32 << 20 }],
            branch_ordering: vec![BranchOrdering::None],
        }
    }

    fn try_run(&mut self, ctx: &mut StepCtx<'_, Self>) -> io::Result<StepOutcome> {
        if let Some(u) = self.held.take()
            && let Err(again) = ctx.outputs.retry(u)
        {
            self.held.put(again);
            return Ok(StepOutcome::Contention);
        }
        if self.pending.is_empty() && self.next < self.len {
            let want = (4u64 << 20).min(self.len - self.next);
            for (offset, len) in
                slice_ranges(self.next, usize::try_from(want).unwrap(), self.streams)
            {
                self.pending.push_back(ReadRequest {
                    ordinal: self.ordinal,
                    stream: 0,
                    seq: u32::try_from(self.ordinal).unwrap(),
                    source: Arc::clone(&self.source),
                    offset,
                    len,
                    last: offset + u64::from(len) == self.len,
                    target: ReadTarget::Input,
                });
                self.ordinal += 1;
            }
            self.next += want;
        }
        let Some(r) = self.pending.pop_front() else { return Ok(StepOutcome::Finished) };
        if let Err(u) = ctx.outputs.push(r) {
            self.held.put(u);
        }
        Ok(StepOutcome::Progress)
    }
}

/// Pool sink dropping the slices (their buffers return to the pool).
struct Drop1(Arc<AtomicU64>);

impl Step for Drop1 {
    type Input = ReadSlice;
    type Outputs = ();

    fn profile(&self) -> StepProfile {
        StepProfile {
            name: "Drop",
            kind: StepKind::Parallel,
            sticky: false,
            output_queues: vec![],
            branch_ordering: vec![],
        }
    }

    fn try_run(&mut self, ctx: &mut StepCtx<'_, Self>) -> io::Result<StepOutcome> {
        match ctx.input.pop() {
            Some(s) => {
                self.0.fetch_add(s.bytes.len() as u64, Ordering::Relaxed);
                Ok(StepOutcome::Progress)
            }
            None if ctx.input.is_drained() => Ok(StepOutcome::Finished),
            None => Ok(StepOutcome::NoProgress),
        }
    }

    fn new_worker_copy(&self) -> Self {
        Drop1(Arc::clone(&self.0))
    }
}

fn bench_pread_slices(c: &mut Criterion) {
    let len: u64 = std::env::var("FGUMI_BENCH_PREAD_BYTES")
        .ok()
        .and_then(|v| v.parse().ok())
        .unwrap_or(1 << 30);
    let dir = bench_dir();
    let path = dir.path().join("pread.bin");
    {
        let mut f = std::io::BufWriter::new(std::fs::File::create(&path).unwrap());
        let mut buf = vec![0u8; 1 << 20];
        let mut x = 1u64;
        let mut written = 0u64;
        while written < len {
            for chunk in buf.chunks_mut(8) {
                x = mix(x);
                chunk.copy_from_slice(&x.to_le_bytes()[..chunk.len()]);
            }
            let n = usize::try_from((len - written).min(buf.len() as u64)).unwrap();
            f.write_all(&buf[..n]).unwrap();
            written += n as u64;
        }
    }
    let file: Arc<dyn PositionalSource> = Arc::new(std::fs::File::open(&path).unwrap());
    let mut group = c.benchmark_group("pread_slices");
    group.throughput(Throughput::Bytes(len));
    group.sample_size(10);
    for streams in [1usize, 4, 8] {
        group.bench_function(format!("{streams}-streams"), |b| {
            b.iter(|| {
                let got = Arc::new(AtomicU64::new(0));
                let builder = Pipeline::builder();
                builder
                    .chain(Requests {
                        source: Arc::clone(&file),
                        len,
                        next: 0,
                        ordinal: 0,
                        streams,
                        pending: VecDeque::new(),
                        held: HeldSlot::new(),
                    })
                    .chain(PreadSlices::spill(SliceBufferPool::new(2 * THREADS + 8), 32 << 20))
                    .chain(Drop1(Arc::clone(&got)))
                    .into_sink_marker();
                builder
                    .build()
                    .unwrap()
                    .run(PipelineConfig { threads: THREADS, ..Default::default() })
                    .unwrap();
                assert_eq!(got.load(Ordering::Relaxed), len);
            });
        });
    }
    group.finish();
}

criterion_group!(benches, bench_spill_merge_supply, bench_pread_slices);
criterion_main!(benches);
