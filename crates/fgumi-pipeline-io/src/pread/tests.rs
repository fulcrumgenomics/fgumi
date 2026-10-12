use std::collections::VecDeque;
use std::sync::Arc;
use std::sync::atomic::Ordering;

use fgumi_bam_io::pread::test_sources::{FaultySource, MemSource, ShortSource};
use fgumi_bam_io::pread::{PositionalSource, SliceBufferPool};
use fgumi_pipeline_core::builder::{Pipeline, PipelineConfig};
use fgumi_pipeline_core::testing::{StepProbe, assert_admission_contract};
use fgumi_pipeline_core::{
    BranchOrdering, HeldSlot, PhaseCap, PipelineError, PoolPlacement, Unpushed,
    outputs::OrderedBytesSingle,
    queues::QueueSpec,
    step::{Step, StepCtx, StepKind, StepOutcome, StepProfile},
};
use rstest::rstest;

use super::*;

/// `Serial` source that emits prepared requests in order.
struct RequestSource {
    reqs: VecDeque<ReadRequest>,
    held: HeldSlot<Unpushed<ReadRequest>>,
}

impl Step for RequestSource {
    type Input = ();
    type Outputs = OrderedBytesSingle<ReadRequest>;

    fn profile(&self) -> StepProfile {
        StepProfile {
            name: "RequestSource",
            kind: StepKind::Serial,
            sticky: true,
            output_queues: vec![QueueSpec::ByteBounded { limit_bytes: 1 << 20 }],
            branch_ordering: vec![BranchOrdering::ByItemOrdinal],
        }
    }

    fn try_run(&mut self, ctx: &mut StepCtx<'_, Self>) -> io::Result<StepOutcome> {
        if matches!(
            ctx.outputs.retry_held(&mut self.held),
            fgumi_pipeline_core::HeldRetry::StillHeld
        ) {
            return Ok(StepOutcome::Contention);
        }
        let Some(req) = self.reqs.pop_front() else { return Ok(StepOutcome::Finished) };
        if let Err(unpushed) = ctx.outputs.push(req) {
            self.held.put(unpushed);
        }
        Ok(StepOutcome::Progress)
    }
}

/// `Serial` sink recording slices in arrival order.
struct SliceSink {
    got: Arc<parking_lot::Mutex<Vec<ReadSlice>>>,
}

impl Step for SliceSink {
    type Input = ReadSlice;
    type Outputs = ();

    fn profile(&self) -> StepProfile {
        StepProfile {
            name: "SliceSink",
            kind: StepKind::Serial,
            sticky: false,
            output_queues: vec![],
            branch_ordering: vec![],
        }
    }

    fn try_run(&mut self, ctx: &mut StepCtx<'_, Self>) -> io::Result<StepOutcome> {
        match ctx.input.pop() {
            Some(s) => {
                self.got.lock().push(s);
                Ok(StepOutcome::Progress)
            }
            None if ctx.input.is_drained() => Ok(StepOutcome::Finished),
            None => Ok(StepOutcome::NoProgress),
        }
    }
}

/// Requests covering `[0, len)` of `src` in `slice`-byte reads: dense ordinals
/// and seqs, `last` on the final one.
fn requests_over(src: &Arc<dyn PositionalSource>, len: usize, slice: usize) -> Vec<ReadRequest> {
    let n = len.div_ceil(slice);
    (0..n)
        .map(|i| {
            let off = i * slice;
            let this = slice.min(len - off);
            ReadRequest {
                ordinal: i as u64,
                stream: 0,
                seq: u32::try_from(i).unwrap(),
                source: Arc::clone(src),
                offset: off as u64,
                len: u32::try_from(this).unwrap(),
                last: i + 1 == n,
                target: ReadTarget::Input,
            }
        })
        .collect()
}

fn run_step(
    reqs: Vec<ReadRequest>,
    step: PreadSlices,
    threads: usize,
) -> Result<Vec<ReadSlice>, PipelineError> {
    let got = Arc::new(parking_lot::Mutex::new(Vec::new()));
    let builder = Pipeline::builder();
    builder
        .chain(RequestSource { reqs: reqs.into(), held: HeldSlot::new() })
        .chain(step)
        .chain(SliceSink { got: Arc::clone(&got) })
        .into_sink_marker();
    let pipeline = builder.build().expect("build");
    pipeline.run(PipelineConfig { threads, ..Default::default() })?;
    Ok(std::mem::take(&mut *got.lock()))
}

fn instance(input: bool, pool: Arc<SliceBufferPool>) -> PreadSlices {
    if input {
        PreadSlices::input(pool, Arc::new(InputLedger::default()), 4 << 20)
    } else {
        PreadSlices::spill(pool, 4 << 20)
    }
}

fn run_pread(reqs: Vec<ReadRequest>, input: bool, threads: usize) -> Vec<ReadSlice> {
    run_step(reqs, instance(input, SliceBufferPool::new(8)), threads).expect("run")
}

fn run_pread_expect_err(reqs: Vec<ReadRequest>, input: bool, threads: usize) -> PipelineError {
    match run_step(reqs, instance(input, SliceBufferPool::new(8)), threads) {
        Ok(_) => panic!("the run must fail"),
        Err(e) => e,
    }
}

#[rstest]
#[case::input_ordered(true)]
#[case::spill_unordered(false)]
fn every_slice_is_delivered_once_with_its_bytes(#[case] input: bool) {
    let data: Vec<u8> = (0..(3 << 20)).map(|i| u8::try_from(i % 251).unwrap()).collect();
    let src: Arc<dyn PositionalSource> = Arc::new(MemSource::new(data.clone()));
    let reqs = requests_over(&src, data.len(), 512 << 10);
    let want_offsets: Vec<u64> = reqs.iter().map(|r| r.offset).collect();
    let slices = run_pread(reqs, input, 4);
    let mut got_offsets: Vec<u64> = slices.iter().map(|s| s.offset).collect();
    got_offsets.sort_unstable();
    assert_eq!(got_offsets, want_offsets, "every request delivered exactly once");
    for s in &slices {
        let off = usize::try_from(s.offset).unwrap();
        assert_eq!(&*s.bytes, &data[off..off + s.bytes.len()]);
    }
    assert_eq!(slices.iter().filter(|s| s.last).count(), 1);
    if input {
        assert!(
            slices.windows(2).all(|w| w[0].ordinal < w[1].ordinal),
            "the input edge is ordinal-ordered"
        );
    }
}

/// A source shorter than its reported length fails the run.
#[test]
fn pread_slices_short_read_fails_closed() {
    let src: Arc<dyn PositionalSource> = Arc::new(ShortSource::new(vec![0u8; 100_000], 1 << 20));
    let err = run_pread_expect_err(requests_over(&src, 1 << 20, 256 << 10), true, 2);
    assert!(err.to_string().contains("EOF"), "{err}");
}

#[test]
fn pread_error_surfaces() {
    let src: Arc<dyn PositionalSource> = Arc::new(FaultySource::new(vec![0u8; 1 << 20], 300 << 10));
    let err = run_pread_expect_err(requests_over(&src, 1 << 20, 256 << 10), false, 2);
    assert!(err.to_string().contains("injected"), "{err}");
}

#[test]
fn input_instance_bumps_landed_per_slice() {
    let ledger = Arc::new(InputLedger::default());
    let src: Arc<dyn PositionalSource> = Arc::new(MemSource::new(vec![1u8; 2 << 20]));
    let step = PreadSlices::input(SliceBufferPool::new(4), Arc::clone(&ledger), 4 << 20);
    let _ = run_step(requests_over(&src, 2 << 20, 1 << 20), step, 2).expect("run");
    assert_eq!(ledger.landed.load(Ordering::Relaxed), 2 << 20);
}

#[test]
fn buffers_are_recycled_through_the_pool() {
    let pool = SliceBufferPool::new(8);
    let src: Arc<dyn PositionalSource> = Arc::new(MemSource::new(vec![0u8; 8 << 20]));
    let slices =
        run_step(requests_over(&src, 8 << 20, 1 << 20), instance(false, Arc::clone(&pool)), 2)
            .expect("run");
    assert_eq!(pool.resident_bytes(), 8 << 20, "every collected slice is a live lease");
    drop(slices);
    assert_eq!(pool.resident_bytes(), 0);
    assert!(pool.free_len() >= 1);
}

#[test]
fn profile_and_placement() {
    let pool = SliceBufferPool::new(1);
    let i = PreadSlices::input(Arc::clone(&pool), Arc::new(InputLedger::default()), 1 << 20);
    let s = PreadSlices::spill(pool, 1 << 20);
    assert_eq!(i.profile().name, "PreadInputSlices");
    assert_eq!(s.profile().name, "PreadSpillSlices");
    assert_eq!(i.profile().kind, StepKind::Parallel);
    assert_eq!(i.pool_placement(), PoolPlacement::ExcludeReader);
    assert_eq!(s.pool_placement(), PoolPlacement::AllWorkers);
    assert_eq!(i.profile().branch_ordering, vec![BranchOrdering::ByItemOrdinal]);
    assert_eq!(s.profile().branch_ordering, vec![BranchOrdering::None]);
    assert!(i.phase_cap().is_none() && s.phase_cap().is_none(), "uncapped unless given a cap");
}

/// A refused permit leaves the request queued (the oracle is the input queue,
/// not the cap's counters).
#[test]
fn refused_permit_leaves_the_request_queued() {
    let cap = PhaseCap::new("sort-phase1", 1);
    let held = cap.try_acquire().expect("the test holds the only permit");
    let src: Arc<dyn PositionalSource> = Arc::new(MemSource::new(vec![7u8; 1 << 20]));
    let mut step =
        PreadSlices::spill(SliceBufferPool::new(2), 1 << 20).with_phase_cap(Some(Arc::clone(&cap)));
    let mut probe = StepProbe::new(&step);
    probe.push_input(requests_over(&src, 1 << 20, 1 << 20).remove(0));
    assert_eq!(probe.try_run(&mut step).unwrap(), StepOutcome::Capped);
    assert!(!probe.input_is_empty(), "the refused step must not pop its request");
    drop(held);
    assert_eq!(probe.try_run(&mut step).unwrap(), StepOutcome::Progress);
    let out = probe.take_output::<ReadSlice>(0);
    let slice = fgumi_pipeline_core::InputHandle::pop(&out).expect("one slice");
    assert!(slice.bytes.len() == 1 << 20 && slice.bytes.iter().all(|&b| b == 7));
}

/// No permit is taken without input: with the only permit held elsewhere and
/// an empty, open input the step reports `NoProgress`, never `Capped`.
#[test]
fn pread_does_not_take_a_permit_without_input() {
    let cap = PhaseCap::new("sort-phase1", 1);
    let _held = cap.try_acquire().unwrap();
    let mut step =
        PreadSlices::spill(SliceBufferPool::new(2), 1 << 20).with_phase_cap(Some(Arc::clone(&cap)));
    let probe = StepProbe::new(&step);
    let refused_before = cap.refused();
    assert_eq!(probe.try_run(&mut step).unwrap(), StepOutcome::NoProgress);
    assert_eq!(cap.refused(), refused_before, "an idle poll is not a refusal");
}

/// The full admission contract every capped pool step owes.
#[rstest]
#[case::input(true)]
#[case::spill(false)]
fn pread_slices_meets_the_admission_contract(#[case] input: bool) {
    let src: Arc<dyn PositionalSource> = Arc::new(MemSource::new(vec![3u8; 4096]));
    let item = requests_over(&src, 4096, 4096).remove(0);
    assert_admission_contract(
        "sort-phase1",
        |cap| instance(input, SliceBufferPool::new(2)).with_phase_cap(cap),
        item,
    );
}

/// The histogram is optional and, when present, sees every request.
#[test]
fn hist_records_every_request() {
    let hist = Arc::new(RequestSizeHist::default());
    let src: Arc<dyn PositionalSource> = Arc::new(MemSource::new(vec![0u8; 3 << 20]));
    let step = PreadSlices::spill(SliceBufferPool::new(4), 4 << 20).with_hist(Arc::clone(&hist));
    let _ = run_step(requests_over(&src, 3 << 20, 1 << 20), step, 2).expect("run");
    assert_eq!(hist.count(), 3);
    assert_eq!(hist.p50_p90_min(), (1 << 20, 1 << 20, 1 << 20));
    let (mean, max) = hist.inflight_mean_max();
    assert!(mean >= 1.0 && max >= 1, "{mean} {max}");
}

#[rstest]
#[case::empty(&[], (0, 0, 0))]
#[case::one_small(&[100], (1 << 10, 1 << 10, 100))]
#[case::mixed(&[512 << 10, 512 << 10, 1 << 20, 1 << 20, 4 << 20], (1 << 20, 4 << 20, 512 << 10))]
fn hist_percentiles_are_bucket_lower_bounds(#[case] lens: &[u32], #[case] want: (u64, u64, u64)) {
    let h = RequestSizeHist::default();
    for &l in lens {
        h.record(l);
    }
    assert_eq!(h.p50_p90_min(), want);
}
